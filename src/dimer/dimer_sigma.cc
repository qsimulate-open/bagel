//
// BAGEL - Parallel electron correlation program.
// Filename: dimer_sigma.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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
#include <vector>
#include <stdexcept>

#include <src/dimer/dimer.h>
#include <src/fci/harrison.h>
#include <src/fci/space.h>
#include <src/util/matrix.h>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace std::chrono;
using namespace bagel;

shared_ptr<Dvec> Dimer::form_sigma_1e(shared_ptr<const Dvec> ccvec, double* hdata, const int ij) const {
  const int nstates = ccvec->ij();

  shared_ptr<Dvec> sigmavec(new Dvec(ccvec->det(), nstates));
  shared_ptr<Determinants> det = space_->finddet(0,0);

  for (int istate = 0; istate < nstates; ++istate) {
    shared_ptr<const Civec> cc = ccvec->data(istate);
    shared_ptr<Civec> sigma = sigmavec->data(istate);

    const int lb = cc->lenb();
    const int la = cc->lena();
    for (int ip = 0; ip != ij; ++ip) {
      const double h = hdata[ip];
      for (auto& iter : det->phia(ip)) {
        const double hc = h * get<1>(iter);
        daxpy_(lb, hc, cc->element_ptr(0, get<2>(iter)), 1, sigma->element_ptr(0, get<0>(iter)), 1); 
      }
    }

    for (int i = 0; i < la; ++i) {
      double* const target_array0 = sigma->element_ptr(0, i);
      const double* const source_array0 = cc->element_ptr(0, i);
      for (int ip = 0; ip != ij; ++ip) {
        const double h = hdata[ip];
        for (auto& iter : det->phib(ip)) {
          const double hc = h * get<1>(iter);
          target_array0[get<0>(iter)] += hc * source_array0[get<2>(iter)];
        }
      }
    }
  }

  return sigmavec;
}

/* Implementing the method as described by Harrison and Zarrabian */
shared_ptr<Dvec> Dimer::form_sigma_2e(shared_ptr<const Dvec> ccvec, double* mo2e_ptr, const int nact) const {
  const int nstate = ccvec->ij();
  const int ij = nact*nact;

  shared_ptr<Dvec> sigmavec(new Dvec(ccvec->det(), nstate));
  sigmavec->zero();

  shared_ptr<Determinants> base_det = space_->finddet(0,0);
  shared_ptr<Determinants> int_det = space_->finddet(-1,-1);

  /* d and e are only used in the alpha-beta case and exist in the (nalpha-1)(nbeta-1) spaces */
  shared_ptr<Dvec> d(new Dvec(int_det, ij));
  shared_ptr<Dvec> e(new Dvec(int_det, ij));

  for (int istate = 0; istate != nstate; ++istate) { 
    shared_ptr<const Civec> cc = ccvec->data(istate);  
    shared_ptr<Civec> sigma = sigmavec->data(istate);  

    vector<pair<string, double> > timing;
    int start = ::clock();

    // (2aa) alpha-alpha contributions
    sigma_2aa(cc, sigma, mo2e_ptr, nact);
    // (2bb) beta-beta contributions
    sigma_2bb(cc, sigma, mo2e_ptr, nact);

    // (2ab) alpha-beta contributions
    d->zero();

    sigma_2ab_1(cc, d, nact);
    sigma_2ab_2(d, e, mo2e_ptr);
    sigma_2ab_3(sigma, e, nact);
  }

  return sigmavec;
}

void Dimer::sigma_2aa(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, double* mo2e_ptr, const int nact) const {
  shared_ptr<Determinants> base_det = space_->finddet(0,0);

  const double* const source_base = cc->data();
  double* target_base = sigma->data();
  const int lb = sigma->lenb();
  
  for (auto aiter = base_det->stringa().begin(); aiter != base_det->stringa().end(); ++aiter, target_base+=lb) {
    bitset<nbit__> nstring = *aiter;
    for (int i = 0; i != nact; ++i) {
      if (!nstring[i]) continue;
      for (int j = 0; j < i; ++j) {
        if(!nstring[j]) continue;
        double* mo2e_ij = mo2e_ptr + nact*nact*j + nact*nact*nact*i;
        const int ij_phase = base_det->sign(nstring,i,j);
        bitset<nbit__> string_ij = nstring; 
        string_ij.reset(i); string_ij.reset(j);
        for (int l = 0; l != nact; ++l) {
          if (string_ij[l]) continue;
          for (int k = 0; k < l; ++k) {
            if (string_ij[k]) continue;
            const int kl_phase = base_det->sign(string_ij,l,k);
            const double phase = -static_cast<double>(ij_phase*kl_phase);
            bitset<nbit__> string_ijkl = string_ij;
            string_ijkl.set(k); string_ijkl.set(l);
            const double temp = phase * ( mo2e_ij[l + nact*k] - mo2e_ij[k + nact*l] );
            const double* source = source_base + base_det->lexical<0>(string_ijkl)*lb;
            daxpy_(lb, temp, source, 1, target_base, 1);
          }
        }
      }
    }
  }
}

void Dimer::sigma_2bb(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, double* mo2e_ptr, const int nact) const {
  const shared_ptr<Determinants> base_det = space_->finddet(0,0);

  const double* const source_base = cc->data();
  double* target_base = sigma->data();

  const int la = sigma->lena();
  const int lb = sigma->lenb();
  for (auto biter = base_det->stringb().begin(); biter != base_det->stringb().end(); ++biter,++target_base) {
    bitset<nbit__> nstring = *biter;
    for (int i = 0; i != nact; ++i) {
      if (!nstring[i]) continue;
      for (int j = 0; j < i; ++j) {
        if(!nstring[j]) continue;
        const int ij_phase = base_det->sign(nstring,i,j);
        double* mo2e_ij = mo2e_ptr + nact*nact*j + nact*nact*nact*i;
        bitset<nbit__> string_ij = nstring;
        string_ij.reset(i); string_ij.reset(j);
        for (int l = 0; l != nact; ++l) {
          if (string_ij[l]) continue;
          for (int k = 0; k < l; ++k) {
            if (string_ij[k]) continue;
            const int kl_phase = base_det->sign(string_ij,l,k);
            const double phase = -static_cast<double>(ij_phase*kl_phase);
            bitset<nbit__> string_ijkl = string_ij;
            string_ijkl.set(k); string_ijkl.set(l);
            const double temp = phase * ( mo2e_ij[l + nact*k] - mo2e_ij[k + nact*l] );
            const double* source = source_base + base_det->lexical<1>(string_ijkl);
            daxpy_(la, temp, source, lb, target_base, lb);
          }
        }
      }
    }
  }
}

void Dimer::sigma_2ab_1(shared_ptr<const Civec> cc, shared_ptr<Dvec> d, const int nact) const {

  shared_ptr<Determinants> int_det = space_->finddet(-1,-1);
  shared_ptr<Determinants> base_det = space_->finddet(0,0);

  const int lbt = int_det->lenb();
  const int lbs = base_det->lenb();
  const double* source_base = cc->data();

  for (int k = 0; k < nact; ++k) {
    for (int l = 0; l < nact; ++l) {
      double* target_base = d->data(k*nact + l)->data();
      for (auto aiter = int_det->phiupa(k).begin(); aiter != int_det->phiupa(k).end(); ++aiter) {
        double *target = target_base + get<2>(*aiter)*lbt;
        const double *source = source_base + get<0>(*aiter)*lbs;
        for (auto biter = int_det->phiupb(l).begin(); biter != int_det->phiupb(l).end(); ++biter) {
          const double sign = static_cast<double>(get<1>(*aiter)*get<1>(*biter));
          target[get<2>(*biter)] += sign * source[get<0>(*biter)];
        }
      }
    }
  }
}

void Dimer::sigma_2ab_2(shared_ptr<Dvec> d, shared_ptr<Dvec> e, double* mo2e_ptr) const {
  const int la = d->lena();
  const int lb = d->lenb();
  const int ij = d->ij();
  const int lenab = la*lb;
  dgemm_("n", "n", lenab, ij, ij, 1.0, d->data(), lenab, mo2e_ptr, ij, 0.0, e->data(), lenab);
}

void Dimer::sigma_2ab_3(shared_ptr<Civec> sigma, shared_ptr<Dvec> e, const int nact) const {
  const shared_ptr<Determinants> base_det = space_->finddet(0,0);
  const shared_ptr<Determinants> int_det = space_->finddet(-1,-1);

  const int lbt = base_det->lenb();
  const int lbs = int_det->lenb();
  double* target_base = sigma->data();

  for (int i = 0; i < nact; ++i) {
    for (int j = 0; j < nact; ++j) {
      const double* source_base = e->data(i*nact + j)->data();
      for (auto aiter = int_det->phiupa(i).begin(); aiter != int_det->phiupa(i).end(); ++aiter) {
        double *target = target_base + get<0>(*aiter)*lbt;
        const double *source = source_base + get<2>(*aiter)*lbs;
        for (auto biter = int_det->phiupb(j).begin(); biter != int_det->phiupb(j).end(); ++biter) {
          const double sign = static_cast<double>(get<1>(*aiter)*get<1>(*biter));
          target[get<0>(*biter)] += sign * source[get<2>(*biter)];
        }
      } 
    }
  }
}
