//
// BAGEL - Parallel electron correlation program.
// Filename: fci_compute.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <src/fci/harrison.h>
#include <src/fci/space.h>
#include <src/util/davidson.h>
#include <src/util/constants.h>
#include <vector>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace bagel;

/* Implementing the method as described by Harrison and Zarrabian */
shared_ptr<Dvec> HarrisonZarrabian::form_sigma(shared_ptr<const Dvec> ccvec, shared_ptr<const MOFile> jop,
                     const vector<int>& conv) const { // d and e are scratch area for D and E intermediates 
  const int ij = norb_*norb_; 

  const int nstate = ccvec->ij();

  auto sigmavec = make_shared<Dvec>(ccvec->det(), nstate);
  sigmavec->zero();

  shared_ptr<Determinants> base_det = space_->finddet(0,0);
  shared_ptr<Determinants> int_det = space_->finddet(-1,-1);

  /* d and e are only used in the alpha-beta case and exist in the (nalpha-1)(nbeta-1) spaces */
  auto d = make_shared<Dvec>(int_det, ij);
  auto e = make_shared<Dvec>(int_det, ij);

  for (int istate = 0; istate != nstate; ++istate) {
    Timer pdebug(2);
    if (conv[istate]) continue;
    shared_ptr<const Civec> cc = ccvec->data(istate);  
    shared_ptr<Civec> sigma = sigmavec->data(istate);  

    // (task1) one-electron alpha: sigma(Psib, Psi'a) += sign h'(ij) C(Psib, Psia) 
    sigma_1(cc, sigma, jop);

    // (task2) two electron contributions
    // (2aa) alpha-alpha contributions
    sigma_2aa(cc,sigma,jop);
    pdebug.tick_print("task2aa");
    
    // (2bb) beta-beta contributions
    /* Mostly the same as the alpha-alpha, except for data storage */
    sigma_2bb(cc, sigma, jop);
    pdebug.tick_print("task2bb");

    // (2ab) alpha-beta contributions
    /* Resembles more the Knowles & Handy FCI terms */
    d->zero();

    sigma_2ab_1(cc, d);
    pdebug.tick_print("task2ab-1");

    sigma_2ab_2(d, e, jop);
    pdebug.tick_print("task2ab-2");

    sigma_2ab_3(sigma, e);
    pdebug.tick_print("task2ab-3");
    
    // (task3) one-electron beta: sigma(Psib', Psia) += sign h'(ij) C(Psib, Psia)
    sigma_3(cc, sigma, jop);
    pdebug.tick_print("task3");

  }

  return sigmavec;
}

void HarrisonZarrabian::sigma_1(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, shared_ptr<const MOFile> jop) const {
  assert(cc->det() == sigma->det());
  const int ij = nij(); 
  const int lb = cc->lenb();
  for (int ip = 0; ip != ij; ++ip) {
    const double h = jop->mo1e(ip);
    for (auto& iter : cc->det()->phia(ip)) {
      const double hc = h * iter.sign;
      daxpy_(lb, hc, cc->element_ptr(0, iter.source), 1, sigma->element_ptr(0, iter.target), 1); 
    }
  }
}

void HarrisonZarrabian::sigma_3(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, shared_ptr<const MOFile> jop) const {
  const int la = cc->lena();
  const int ij = nij();

  for (int i = 0; i < la; ++i) {
    double* const target_array0 = sigma->element_ptr(0, i);
    const double* const source_array0 = cc->element_ptr(0, i);
    for (int ip = 0; ip != ij; ++ip) {
      const double h = jop->mo1e(ip);
      for (auto& iter : cc->det()->phib(ip)) {
        const double hc = h * iter.sign;
        target_array0[iter.target] += hc * source_array0[iter.source];
      }
    }
  }
}

void HarrisonZarrabian::sigma_2aa(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, shared_ptr<const MOFile> jop) const {
  shared_ptr<Determinants> base_det = space_->finddet(0,0);

  const int norb = norb_;

  const double* const source_base = cc->data();
  double* target_base = sigma->data();
  const int lb = sigma->lenb();
  
  for (auto aiter = base_det->stringa().begin(); aiter != base_det->stringa().end(); ++aiter, target_base+=lb) {
    bitset<nbit__> nstring = *aiter;
    for (int i = 0; i != norb; ++i) {
      if (!nstring[i]) continue;
      for (int j = 0; j < i; ++j) {
        if(!nstring[j]) continue;
        const int ij_phase = base_det->sign(nstring,i,j);
        bitset<nbit__> string_ij = nstring; 
        string_ij.reset(i); string_ij.reset(j);
        for (int l = 0; l != norb; ++l) {
          if (string_ij[l]) continue;
          for (int k = 0; k < l; ++k) {
            if (string_ij[k]) continue;
            const int kl_phase = base_det->sign(string_ij,l,k);
            const double phase = -static_cast<double>(ij_phase*kl_phase);
            bitset<nbit__> string_ijkl = string_ij;
            string_ijkl.set(k); string_ijkl.set(l);
            const double temp = phase * ( jop->mo2e_hz(i,j,k,l) - jop->mo2e_hz(i,j,l,k) );
            const double* source = source_base + base_det->lexical<0>(string_ijkl)*lb;
            daxpy_(lb, temp, source, 1, target_base, 1);
          }
        }
      }
    }
  }
}

void HarrisonZarrabian::sigma_2bb(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, shared_ptr<const MOFile> jop) const {
  const shared_ptr<Determinants> base_det = space_->finddet(0,0);

  const double* const source_base = cc->data();
  double* target_base = sigma->data();
  const int norb = norb_;

  const int la = sigma->lena();
  const int lb = sigma->lenb();
  for (auto biter = base_det->stringb().begin(); biter != base_det->stringb().end(); ++biter,++target_base) {
    bitset<nbit__> nstring = *biter;
    for (int i = 0; i != norb; ++i) {
      if (!nstring[i]) continue;
      for (int j = 0; j < i; ++j) {
        if(!nstring[j]) continue;
        const int ij_phase = base_det->sign(nstring,i,j);
        bitset<nbit__> string_ij = nstring;
        string_ij.reset(i); string_ij.reset(j);
        for (int l = 0; l != norb; ++l) {
          if (string_ij[l]) continue;
          for (int k = 0; k < l; ++k) {
            if (string_ij[k]) continue;
            const int kl_phase = base_det->sign(string_ij,l,k);
            const double phase = -static_cast<double>(ij_phase*kl_phase);
            bitset<nbit__> string_ijkl = string_ij;
            string_ijkl.set(k); string_ijkl.set(l);
            const double temp = phase * ( jop->mo2e_hz(i,j,k,l) - jop->mo2e_hz(i,j,l,k) );
            const double* source = source_base + base_det->lexical<1>(string_ijkl);
            daxpy_(la, temp, source, lb, target_base, lb);
          }
        }
      }
    }
  }
}

void HarrisonZarrabian::sigma_2ab_1(shared_ptr<const Civec> cc, shared_ptr<Dvec> d) const {
  const int norb = norb_;

  shared_ptr<Determinants> int_det = space_->finddet(-1,-1);
  shared_ptr<Determinants> base_det = space_->finddet(0,0);

  const int lbt = int_det->lenb();
  const int lbs = base_det->lenb();
  const double* source_base = cc->data();

  for (int k = 0; k < norb; ++k) {
    for (int l = 0; l < norb; ++l) {
      double* target_base = d->data(k*norb + l)->data();
      for (auto& aiter : int_det->phiupa(k)) {
        double *target = target_base + aiter.source*lbt;
        const double *source = source_base + aiter.target*lbs;
        for (auto& biter : int_det->phiupb(l)) {
          const double sign = aiter.sign * biter.sign;
          target[biter.source] += sign * source[biter.target];
        }
      }
    }
  }
}

void HarrisonZarrabian::sigma_2ab_2(shared_ptr<Dvec> d, shared_ptr<Dvec> e, shared_ptr<const MOFile> jop) const {
  const int la = d->lena();
  const int lb = d->lenb();
  const int ij = d->ij();
  const int lenab = la*lb;
  dgemm_("n", "n", lenab, ij, ij, 1.0, d->data(), lenab, jop->mo2e_ptr(), ij, 0.0, e->data(), lenab);
}

void HarrisonZarrabian::sigma_2ab_3(shared_ptr<Civec> sigma, shared_ptr<Dvec> e) const {
  const shared_ptr<Determinants> base_det = space_->finddet(0,0);
  const shared_ptr<Determinants> int_det = space_->finddet(-1,-1);

  const int norb = norb_;
  const int lbt = base_det->lenb();
  const int lbs = int_det->lenb();
  double* target_base = sigma->data();

  for (int i = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      const double* source_base = e->data(i*norb + j)->data();
      for (auto& aiter : int_det->phiupa(i)) {
        double *target = target_base + aiter.target*lbt;
        const double *source = source_base + aiter.source*lbs;
        for (auto& biter : int_det->phiupb(j)) {
          const double sign = aiter.sign * biter.sign;
          target[biter.target] += sign * source[biter.source];
        }
      } 
    }
  }
}
