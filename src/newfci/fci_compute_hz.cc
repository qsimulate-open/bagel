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
#include <src/newfci/fci_hz.h>
#include <src/newfci/space.h>
#include <src/util/davidson.h>
#include <src/util/constants.h>
#include <vector>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace std::chrono;
using namespace bagel;

/* Implementing the method as described by Harrison and Zarrabian */
shared_ptr<NewDvec> NewFCI_HZ::form_sigma(shared_ptr<const NewDvec> ccvec, shared_ptr<const NewMOFile> jop,
                     const vector<int>& conv) const { // d and e are scratch area for D and E intermediates 
  const int ij = norb_*norb_; 

  const int nstate = ccvec->ij();

  shared_ptr<NewDvec> sigmavec(new NewDvec(ccvec->det(), nstate));
  sigmavec->zero();

  shared_ptr<NewDeterminants> base_det = space_->finddet(0,0);
  shared_ptr<NewDeterminants> int_det = space_->finddet(-1,-1);

  /* d and e are only used in the alpha-beta case and exist in the (nalpha-1)(nbeta-1) spaces */
  shared_ptr<NewDvec> d(new NewDvec(int_det, ij));
  shared_ptr<NewDvec> e(new NewDvec(int_det, ij));

  for (int istate = 0; istate != nstate; ++istate) {
    if (conv[istate]) continue;
    shared_ptr<const NewCivec> cc = ccvec->data(istate);  
    shared_ptr<NewCivec> sigma = sigmavec->data(istate);  

    vector<pair<string, double> > timing;
    int start = ::clock();

    // (task1) one-electron alpha: sigma(Psib, Psi'a) += sign h'(ij) C(Psib, Psia) 
    sigma_1(cc, sigma, jop);
    if (tprint) print_timing_("task1", start, timing);

    // (task2) two electron contributions

    // (2aa) alpha-alpha contributions
    sigma_2aa(cc,sigma,jop);
    if (tprint) print_timing_("task2aa", start, timing);
    
    // (2bb) beta-beta contributions
    /* Mostly the same as the alpha-alpha, except for data storage */
    sigma_2bb(cc, sigma, jop);
    if (tprint) print_timing_("task2bb", start, timing);

    // (2ab) alpha-beta contributions
    /* Resembles more the Knowles & Handy FCI terms */
    d->zero();

    sigma_2ab_1(cc, d);
    if (tprint) print_timing_("task2ab-1", start, timing);

    sigma_2ab_2(d, e, jop);
    if (tprint) print_timing_("task2ab-2", start, timing);

    sigma_2ab_3(sigma, e);
    if (tprint) print_timing_("task2ab-3", start, timing);
    
    // (task3) one-electron beta: sigma(Psib', Psia) += sign h'(ij) C(Psib, Psia)
    sigma_3(cc, sigma, jop);

    if (tprint) {
      print_timing_("task3", start, timing);
      cout << "     timing info" << endl;
      for (auto iter = timing.begin(); iter != timing.end(); ++iter)
        cout << "    " << setw(10) << iter->first << setw(10) << setprecision(2) << iter ->second << endl;
    }
  }

  return sigmavec;
}

void NewFCI_HZ::sigma_1(shared_ptr<const NewCivec> cc, shared_ptr<NewCivec> sigma, shared_ptr<const NewMOFile> jop) const {
  assert(cc->det() == sigma->det());
  const int ij = nij(); 
  const int lb = cc->lenb();
  for (int ip = 0; ip != ij; ++ip) {
    const double h = jop->mo1e(ip);
    for (auto iter = cc->det()->phia(ip).begin();  iter != cc->det()->phia(ip).end(); ++iter) {
      const double hc = h * get<1>(*iter);
      daxpy_(lb, hc, cc->element_ptr(0, get<2>(*iter)), 1, sigma->element_ptr(0, get<0>(*iter)), 1); 
    }
  }
}

void NewFCI_HZ::sigma_3(shared_ptr<const NewCivec> cc, shared_ptr<NewCivec> sigma, shared_ptr<const NewMOFile> jop) const {
  const int la = cc->lena();
  const int ij = nij();

  for (int i = 0; i < la; ++i) {
    double* const target_array0 = sigma->element_ptr(0, i);
    const double* const source_array0 = cc->element_ptr(0, i);
    for (int ip = 0; ip != ij; ++ip) {
      const double h = jop->mo1e(ip);
      for (auto iter = cc->det()->phib(ip).begin();  iter != cc->det()->phib(ip).end(); ++iter) {
        const double hc = h * get<1>(*iter);
        target_array0[get<0>(*iter)] += hc * source_array0[get<2>(*iter)];
      }
    }
  }
}

void NewFCI_HZ::sigma_2aa(shared_ptr<const NewCivec> cc, shared_ptr<NewCivec> sigma, shared_ptr<const NewMOFile> jop) const {
  shared_ptr<NewDeterminants> base_det = space_->finddet(0,0);

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

void NewFCI_HZ::sigma_2bb(shared_ptr<const NewCivec> cc, shared_ptr<NewCivec> sigma, shared_ptr<const NewMOFile> jop) const {
  const shared_ptr<NewDeterminants> base_det = space_->finddet(0,0);

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

void NewFCI_HZ::sigma_2ab_1(shared_ptr<const NewCivec> cc, shared_ptr<NewDvec> d) const {
  const int norb = norb_;

  shared_ptr<NewDeterminants> int_det = space_->finddet(-1,-1);
  shared_ptr<NewDeterminants> base_det = space_->finddet(0,0);

  const int lbt = int_det->lenb();
  const int lbs = base_det->lenb();
  const double* source_base = cc->data();

  for (int k = 0; k < norb; ++k) {
    for (int l = 0; l < norb; ++l) {
      double* target_base = d->data(k*norb + l)->data();
      for (auto aiter = int_det->phiupa(k).begin(); aiter != int_det->phiupa(k).end(); ++aiter) {
        for (auto biter = int_det->phiupb(l).begin(); biter != int_det->phiupb(l).end(); ++biter) {
          const double sign = static_cast<double>(get<1>(*aiter)*get<1>(*biter));
          target_base[get<2>(*aiter)*lbt + get<2>(*biter)] += sign * source_base[get<0>(*aiter)*lbs + get<0>(*biter)];
        }
      }
    }
  }
}

void NewFCI_HZ::sigma_2ab_2(shared_ptr<NewDvec> d, shared_ptr<NewDvec> e, shared_ptr<const NewMOFile> jop) const {
  const int la = d->lena();
  const int lb = d->lenb();
  const int ij = d->ij();
  const int lenab = la*lb;
  dgemm_("n", "n", lenab, ij, ij, 1.0, d->data(), lenab, jop->mo2e_ptr(), ij, 0.0, e->data(), lenab);
}

void NewFCI_HZ::sigma_2ab_3(shared_ptr<NewCivec> sigma, shared_ptr<NewDvec> e) const {
  const shared_ptr<NewDeterminants> base_det = space_->finddet(0,0);
  const shared_ptr<NewDeterminants> int_det = space_->finddet(-1,-1);

  const int norb = norb_;
  const int lbt = base_det->lenb();
  const int lbs = int_det->lenb();
  double* target_base = sigma->data();

  for (int i = 0; i < norb; ++i) {
    for (int j = 0; j < norb; ++j) {
      const double* source_base = e->data(i*norb + j)->data();
      for (auto aiter = int_det->phiupa(i).begin(); aiter != int_det->phiupa(i).end(); ++aiter) {
        for (auto biter = int_det->phiupb(j).begin(); biter != int_det->phiupb(j).end(); ++biter) {
          const double sign = static_cast<double>(get<1>(*aiter)*get<1>(*biter));
          target_base[get<0>(*aiter)*lbt + get<0>(*biter)] += sign * source_base[get<2>(*aiter)*lbs + get<2>(*biter)];
        }
      } 
    }
  }
}
