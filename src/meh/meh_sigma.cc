//
// BAGEL - Parallel electron correlation program.
// Filename: meh_sigma.cc
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

#include <src/meh/meh.h>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace bagel;

shared_ptr<Dvec> MultiExcitonHamiltonian::form_sigma_1e(shared_ptr<const Dvec> ccvec, const double* hdata) const {
  const int nstates = ccvec->ij();

  shared_ptr<const Determinants> det = ccvec->det();
  auto sigmavec = make_shared<Dvec>(det, nstates);

  const int ij = det->norb() * det->norb();
  const int lb = det->lenb();
  const int la = det->lena();

  for (int istate = 0; istate < nstates; ++istate) {
    shared_ptr<const Civec> cc = ccvec->data(istate);
    shared_ptr<Civec> sigma = sigmavec->data(istate);

    for (int ip = 0; ip != ij; ++ip) {
      const double h = hdata[ip];
      for (auto& iter : det->phia(ip)) {
        const double hc = h * iter.sign;
        daxpy_(lb, hc, cc->element_ptr(0, iter.source), 1, sigma->element_ptr(0, iter.target), 1);
      }
    }
  }

  for (int istate = 0; istate < nstates; ++istate) {
    shared_ptr<const Civec> cc = ccvec->data(istate);
    shared_ptr<Civec> sigma = sigmavec->data(istate);

    for (int i = 0; i < la; ++i) {
      double* const target_array0 = sigma->element_ptr(0, i);
      const double* const source_array0 = cc->element_ptr(0, i);
      for (int ip = 0; ip != ij; ++ip) {
        const double h = hdata[ip];
        for (auto& iter : det->phib(ip)) {
          const double hc = h * iter.sign;
          target_array0[iter.target] += hc * source_array0[iter.source];
        }
      }
    }
  }

  return sigmavec;
}

/* Implementing the method as described by Harrison and Zarrabian */
shared_ptr<Dvec> MultiExcitonHamiltonian::form_sigma_2e(shared_ptr<const Dvec> ccvec, const double* mo2e_ptr) const {
  const int nstate = ccvec->ij();
  shared_ptr<const Determinants> base_det = ccvec->det();
  shared_ptr<const Determinants> int_det = base_det->remalpha()->rembeta();

  const int norb = base_det->norb();
  const int ij = norb*norb;

  auto sigmavec = make_shared<Dvec>(base_det, nstate);

  /* d and e are only used in the alpha-beta case and exist in the (nalpha-1)(nbeta-1) spaces */
  auto d = make_shared<Dvec>(int_det, ij);
  auto e = make_shared<Dvec>(int_det, ij);

  for (int istate = 0; istate != nstate; ++istate) {
    shared_ptr<const Civec> cc = ccvec->data(istate);
    shared_ptr<Civec> sigma = sigmavec->data(istate);

    vector<pair<string, double>> timing;
    int start = ::clock();

    // (2aa) alpha-alpha contributions
    sigma_2aa(cc, sigma, mo2e_ptr, norb);
    // (2bb) beta-beta contributions
    sigma_2bb(cc, sigma, mo2e_ptr, norb);

    // (2ab) alpha-beta contributions
    d->zero();

    sigma_2ab_1(cc, d, norb);
    sigma_2ab_2(d, e, mo2e_ptr);
    sigma_2ab_3(sigma, e, norb);
  }

  return sigmavec;
}

void MultiExcitonHamiltonian::sigma_2aa(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, const double* mo2e_ptr, const int norb) const {
  shared_ptr<const Determinants> base_det = cc->det();

  const double* source_base = cc->data();
  double* target_base = sigma->data();
  const int lb = sigma->lenb();

  for (auto aiter = base_det->stringa().begin(); aiter != base_det->stringa().end(); ++aiter, target_base+=lb) {
    bitset<nbit__> nstring = *aiter;
    for (int i = 0; i != norb; ++i) {
      if (!nstring[i]) continue;
      for (int j = 0; j < i; ++j) {
        if (!nstring[j]) continue;
        const double* mo2e_ij = mo2e_ptr + norb*norb*j + norb*norb*norb*i;
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
            const double temp = phase * ( mo2e_ij[l + norb*k] - mo2e_ij[k + norb*l] );
            const double* source = source_base + base_det->lexical<0>(string_ijkl)*lb;
            daxpy_(lb, temp, source, 1, target_base, 1);
          }
        }
      }
    }
  }
}

void MultiExcitonHamiltonian::sigma_2bb(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, const double* mo2e_ptr, const int norb) const {
  shared_ptr<const Determinants> base_det = cc->det();

  const double* source_base = cc->data();
  double* target_base = sigma->data();

  const int la = sigma->lena();
  const int lb = sigma->lenb();
  for (auto biter = base_det->stringb().begin(); biter != base_det->stringb().end(); ++biter,++target_base) {
    bitset<nbit__> nstring = *biter;
    for (int i = 0; i != norb; ++i) {
      if (!nstring[i]) continue;
      for (int j = 0; j < i; ++j) {
        if (!nstring[j]) continue;
        const int ij_phase = base_det->sign(nstring,i,j);
        const double* mo2e_ij = mo2e_ptr + norb*norb*j + norb*norb*norb*i;
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
            const double temp = phase * ( mo2e_ij[l + norb*k] - mo2e_ij[k + norb*l] );
            const double* source = source_base + base_det->lexical<1>(string_ijkl);
            daxpy_(la, temp, source, lb, target_base, lb);
          }
        }
      }
    }
  }
}

void MultiExcitonHamiltonian::sigma_2ab_1(shared_ptr<const Civec> cc, shared_ptr<Dvec> d, const int norb) const {

  shared_ptr<const Determinants> base_det = cc->det();
  shared_ptr<const Determinants> int_det = base_det->remalpha()->rembeta();

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

void MultiExcitonHamiltonian::sigma_2ab_2(shared_ptr<Dvec> d, shared_ptr<Dvec> e, const double* mo2e_ptr) const {
  const int lenab = d->lena() * d->lenb();
  const int ij = d->ij();

  dgemm_("n", "n", lenab, ij, ij, 1.0, d->data(), lenab, mo2e_ptr, ij, 0.0, e->data(), lenab);
}

void MultiExcitonHamiltonian::sigma_2ab_3(shared_ptr<Civec> sigma, shared_ptr<Dvec> e, const int norb) const {
  shared_ptr<const Determinants> base_det = sigma->det();
  shared_ptr<const Determinants> int_det = base_det->remalpha()->rembeta();

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
          const double sign = static_cast<double>(aiter.sign * biter.sign);
          target[biter.target] += sign * source[biter.source];
        }
      }
    }
  }
}
