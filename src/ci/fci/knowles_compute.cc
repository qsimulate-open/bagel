//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_compute.cc
// Copyright (C) 2011 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/ci/fci/knowles.h>

// toggle for timing print out.
static const bool tprint = false;

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::KnowlesHandy)

using namespace std;
using namespace bagel;

shared_ptr<Dvec> KnowlesHandy::form_sigma(shared_ptr<const Dvec> ccvec, shared_ptr<const MOFile> jop,
                     const vector<int>& conv) const { // d and e are scratch area for D and E intermediates

  const int ij = (norb_*(norb_+1))/2;

  const int nstate = ccvec->ij();
  auto sigmavec = make_shared<Dvec>(ccvec->det(), nstate);
  sigmavec->zero();

  // we need two vectors for intermediate quantities
  auto d = make_shared<Dvec>(ccvec->det(), ij);
  auto e = make_shared<Dvec>(ccvec->det(), ij);


  for (int istate = 0; istate != nstate; ++istate) {
    Timer pdebug(3);
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

// The first two are a part of Base because they are needed in the RDM parts
void FCI::sigma_2a1(shared_ptr<const Civec> cc, shared_ptr<Dvec> d) const {
  assert(d->det() == cc->det());
  const int lb = d->lenb();
  const int ij = d->ij();
  const double* const source_base = cc->data();
  for (int ip = 0; ip != ij; ++ip) {
    double* const target_base = d->data(ip)->data();
    for (auto& iter : cc->det()->phia(ip)) {
      const double sign = static_cast<double>(iter.sign);
      double* const target_array = target_base + iter.source*lb;
      blas::ax_plus_y_n(sign, source_base + iter.target*lb, lb, target_array);
    }
  }
}

void FCI::sigma_2a2(shared_ptr<const Civec> cc, shared_ptr<Dvec> d) const {
  assert(d->det() == cc->det());
  const int la = d->lena();
  const int ij = d->ij();
  for (int i = 0; i < la; ++i) {
    const double* const source_array0 = cc->element_ptr(0, i);
    for (int ip = 0; ip != ij; ++ip) {
      double* const target_array0 = d->data(ip)->element_ptr(0, i);
      for (auto& iter : cc->det()->phib(ip)) {
        const double sign = static_cast<double>(iter.sign);
        target_array0[iter.source] += sign * source_array0[iter.target];
      }
    }
  }
}

void KnowlesHandy::sigma_1(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, shared_ptr<const MOFile> jop) const {
  assert(cc->det() == sigma->det());
  const int ij = (norb_*(norb_+1))/2;
  const int lb = cc->lenb();
  for (int ip = 0; ip != ij; ++ip) {
    const double h = jop->mo1e(ip);
    for (auto& iter : cc->det()->phia(ip)) {
      const double hc = h * iter.sign;
      blas::ax_plus_y_n(hc, cc->element_ptr(0, iter.source), lb, sigma->element_ptr(0, iter.target));
    }
  }
}

void KnowlesHandy::sigma_2c1(shared_ptr<Civec> sigma, shared_ptr<const Dvec> e) const {
  const int lb = e->lenb();
  const int ij = e->ij();
  for (int ip = 0; ip != ij; ++ip) {
    const double* const source_base = e->data(ip)->data();
    for (auto& iter : e->det()->phia(ip)) {
      const double sign = static_cast<double>(iter.sign);
      double* const target_array = sigma->element_ptr(0, iter.target);
      blas::ax_plus_y_n(sign, source_base + lb*iter.source, lb, target_array);
    }
  }
}

void KnowlesHandy::sigma_2c2(shared_ptr<Civec> sigma, shared_ptr<const Dvec> e) const {
  const int la = e->lena();
  const int ij = e->ij();
  for (int i = 0; i < la; ++i) {
    double* const target_array0 = sigma->element_ptr(0, i);
    for (int ip = 0; ip != ij; ++ip) {
      const double* const source_array0 = e->data(ip)->element_ptr(0, i);
      for (auto& iter : e->det()->phib(ip)) {
        const double sign = static_cast<double>(iter.sign);
        target_array0[iter.target] += sign * source_array0[iter.source];
      }
    }
  }
}


void KnowlesHandy::sigma_3(shared_ptr<const Civec> cc, shared_ptr<Civec> sigma, shared_ptr<const MOFile> jop) const {
  const int la = cc->lena();
  const int ij = (norb_*(norb_+1))/2;

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

void KnowlesHandy::sigma_2b(shared_ptr<Dvec> d, shared_ptr<Dvec> e, shared_ptr<const MOFile> jop) const {
  const int la = d->lena();
  const int lb = d->lenb();
  const int ij = d->ij();
  const int lenab = la*lb;
  dgemm_("n", "n", lenab, ij, ij, 0.5, d->data(), lenab, jop->mo2e_ptr(), ij,
                                  0.0, e->data(), lenab);
}

