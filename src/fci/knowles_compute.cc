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
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

#include <src/fci/knowles.h>
#include <src/math/davidson.h>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace bagel;

shared_ptr<Dvec> KnowlesHandy::form_sigma(shared_ptr<const Dvec> ccvec, shared_ptr<const MOFile> jop,
                     const vector<int>& conv) const { // d and e are scratch area for D and E intermediates

  const int ij = nij();
  const size_t lenab = ccvec->det()->lena() * ccvec->det()->lenb();

  const int nstate = ccvec->ij();
  auto sigmavec = make_shared<Dvec>(ccvec->det(), nstate);
  sigmavec->zero();

  // we need two vectors for intermediate quantities
  auto d = make_shared<Matrix>(lenab, ij);
  auto e = make_shared<Matrix>(lenab, ij);


  for (int istate = 0; istate != nstate; ++istate) {
    Timer pdebug(2);
    if (conv[istate]) continue;
    shared_ptr<const Civec> cc = ccvec->data(istate);
    shared_ptr<Civec> sigma = sigmavec->data(istate);

    // (task1) one-electron alpha: sigma(Psib, Psi'a) += sign h'(ij) C(Psib, Psia)
    sigma_1(cc, sigma, jop);
    pdebug.tick_print("task1");

    // (task2) two electron contributions
    d->zero();

    // step (c) (task2a-1) D(Phib, Phia, ij) += sign C(Psib, Phi'a)
    sigma_2a1(cc, d);
    pdebug.tick_print("task2a-1");

    // step (d) (task2a-2) D(Phib, Phia, ij) += sign C(Psib', Phia)
    sigma_2a2(cc, d);
    pdebug.tick_print("task2a-2");

    // step (e) (task2b) E(Phib, Phia, kl) = D(Psib, Phia, ij) (ij|kl)
    sigma_2b(d, e, jop);
    pdebug.tick_print("task2b");

    // step (f) (task2c-1) sigma(Phib, Phia') += sign E(Psib, Phia, kl)
    sigma_2c1(sigma, e);
    pdebug.tick_print("task2c-1");

    // step (g) (task2c-2) sigma(Phib', Phia) += sign E(Psib, Phia, kl)
    sigma_2c2(sigma, e);
    pdebug.tick_print("task2c-2");

    // (task3) one-electron beta: sigma(Psib', Psia) += sign h'(ij) C(Psib, Psia)
    sigma_3(cc, sigma, jop);
    pdebug.tick_print("task3");

  }

  return sigmavec;
}

void KnowlesHandy::sigma_2a1(shared_ptr<const Civec> cc, shared_ptr<Matrix> d) const {
  assert(d->ndim() == cc->size());
  const int lb = cc->lenb();
  const int ij = nij();
  const double* const source_base = cc->data();
  for (int ip = 0; ip != ij; ++ip) {
    double* const target_base = d->element_ptr(0, ip);
    for (auto& iter : cc->det()->phia(ip)) {
      const double sign = static_cast<double>(iter.sign);
      double* const target_array = target_base + iter.source*lb;
      daxpy_(lb, sign, source_base + iter.target*lb, 1, target_array, 1);
    }
  }
}

void KnowlesHandy::sigma_2a2(shared_ptr<const Civec> cc, shared_ptr<Matrix> d) const {
  assert(d->ndim() == cc->size());
  const int la = cc->lena();
  const int lb = cc->lenb();
  const int ij = nij();
  for (int i = 0; i < la; ++i) {
    const double* const source_array0 = cc->element_ptr(0, i);
    for (int ip = 0; ip != ij; ++ip) {
      double* const target_array0 = d->element_ptr(i*lb, ip);
      for (auto& iter : cc->det()->phib(ip)) {
        const double sign = static_cast<double>(iter.sign);
        target_array0[iter.source] += sign * source_array0[iter.target];
      }
    }
  }
}

void KnowlesHandy::sigma_1(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, shared_ptr<const MOFile> jop) const {
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

void KnowlesHandy::sigma_2c1(shared_ptr<Civec> sigma, shared_ptr<const Matrix> e) const {
  const int lb = sigma->lenb();
  const int ij = nij();
  for (int ip = 0; ip != ij; ++ip) {
    const double* const source_base = e->element_ptr(0, ip);
    for (auto& iter : sigma->det()->phia(ip)) {
      const double sign = static_cast<double>(iter.sign);
      double* const target_array = sigma->element_ptr(0, iter.target);
      daxpy_(lb, sign, source_base + lb*iter.source, 1, target_array, 1);
    }
  }
}

void KnowlesHandy::sigma_2c2(shared_ptr<Civec> sigma, shared_ptr<const Matrix> e) const {
  const int la = sigma->lena();
  const int lb = sigma->lenb();
  const int ij = nij();
  for (int i = 0; i < la; ++i) {
    double* const target_array0 = sigma->element_ptr(0, i);
    for (int ip = 0; ip != ij; ++ip) {
      const double* const source_array0 = e->element_ptr(i*lb, ip);
      for (auto& iter : sigma->det()->phib(ip)) {
        const double sign = static_cast<double>(iter.sign);
        target_array0[iter.target] += sign * source_array0[iter.source];
      }
    }
  }
}


void KnowlesHandy::sigma_3(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, shared_ptr<const MOFile> jop) const {
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

void KnowlesHandy::sigma_2b(shared_ptr<Matrix> d, shared_ptr<Matrix> e, shared_ptr<const MOFile> jop) const {
  const int ij = d->mdim();
  const int lenab = d->ndim();
  dgemm_("n", "n", lenab, ij, ij, 0.5, d->data(), lenab, jop->mo2e_ptr(), ij,
                                  0.0, e->data(), lenab);
}
