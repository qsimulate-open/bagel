//
// BAGEL - Parallel electron correlation program.
// Filename: zknowles_compute.cc
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell <michaelcaldwell2013@u.northwestern.edu>
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
#include <src/zfci/zfci.h>
#include <src/zfci/zknowles.h>
#include <src/math/zdavidson.h>
#include <src/util/constants.h>
#include <vector>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace bagel;

shared_ptr<ZDvec> ZKnowlesHandy::form_sigma(shared_ptr<const ZDvec> ccvec, shared_ptr<const MOFile> jop,
                     const vector<int>& conv) const { // d and e are scratch area for D and E intermediates

  const int ij = nij();
  const int nstate = ccvec->ij();
  auto sigmavec = make_shared<ZDvec>(ccvec->det(), nstate);
  sigmavec->zero();
  // we need two vectors for intermediate quantities
  auto d = make_shared<ZDvec>(ccvec->det(), ij);
  auto e = make_shared<ZDvec>(ccvec->det(), ij);


  for (int istate = 0; istate != nstate; ++istate) {
    Timer pdebug(2);
    if (conv[istate]) continue;
    shared_ptr<const ZCivec> cc = ccvec->data(istate);
    shared_ptr<ZCivec> sigma = sigmavec->data(istate);

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
void ZFCI::sigma_2a1(shared_ptr<const ZCivec> cc, shared_ptr<ZDvec> d) const {
  assert(d->det() == cc->det());
  const int lb = d->lenb();
  const int ij = d->ij();
  const complex<double>* const source_base = cc->data();
  for (int ip = 0; ip != ij; ++ip) {
    complex<double>* const target_base = d->data(ip)->data();
    for (auto& iter : cc->det()->phia(ip)) {
      const double sign = static_cast<double>(iter.sign);
      complex<double>* const target_array = target_base + iter.source*lb;
      zaxpy_(lb, sign, source_base + iter.target*lb, 1, target_array, 1);
    }
  }
}

void ZFCI::sigma_2a2(shared_ptr<const ZCivec> cc, shared_ptr<ZDvec> d) const {
  assert(d->det() == cc->det());
  const int la = d->lena();
  const int ij = d->ij();
  for (int i = 0; i < la; ++i) {
    const complex<double>* const source_array0 = cc->element_ptr(0, i);
    for (int ip = 0; ip != ij; ++ip) {
      complex<double>* const target_array0 = d->data(ip)->element_ptr(0, i);
      for (auto& iter : cc->det()->phib(ip)) {
        const double sign = static_cast<double>(iter.sign);
        target_array0[iter.source] += sign * source_array0[iter.target];
      }
    }
  }
}

void ZKnowlesHandy::sigma_1(shared_ptr<const ZCivec> cc, shared_ptr<ZCivec> sigma, shared_ptr<const MOFile> jop) const {
  assert(cc->det() == sigma->det());
  const int ij = nij();
  const int lb = cc->lenb();
  for (int ip = 0; ip != ij; ++ip) {
    const double h = jop->mo1e(ip);
    for (auto& iter : cc->det()->phia(ip)) {
      const double hc = h * iter.sign;
      zaxpy_(lb, hc, cc->element_ptr(0, iter.source), 1, sigma->element_ptr(0, iter.target), 1);
    }
  }
}

void ZKnowlesHandy::sigma_2c1(shared_ptr<ZCivec> sigma, shared_ptr<const ZDvec> e) const {
  const int lb = e->lenb();
  const int ij = e->ij();
  for (int ip = 0; ip != ij; ++ip) {
    const complex<double>* const source_base = e->data(ip)->data();
    for (auto& iter : (e->det())->phia(ip)) {
      const double sign = static_cast<double>(iter.sign);
      complex<double>* const target_array = sigma->element_ptr(0, iter.target);
      zaxpy_(lb, sign, source_base + lb*iter.source, 1, target_array, 1);
    }
  }
}

void ZKnowlesHandy::sigma_2c2(shared_ptr<ZCivec> sigma, shared_ptr<const ZDvec> e) const {
  const int la = e->lena();
  const int ij = e->ij();
  for (int i = 0; i < la; ++i) {
    complex<double>* const target_array0 = sigma->element_ptr(0, i);
    for (int ip = 0; ip != ij; ++ip) {
      const complex<double>* const source_array0 = e->data(ip)->element_ptr(0, i);
      for (auto& iter : e->det()->phib(ip)) {
        const double sign = static_cast<double>(iter.sign);
        target_array0[iter.target] += sign * source_array0[iter.source];
      }
    }
  }
}

void ZKnowlesHandy::sigma_3(shared_ptr<const ZCivec> cc, shared_ptr<ZCivec> sigma, shared_ptr<const MOFile> jop) const {
  const int la = cc->lena();
  const int ij = nij();

  for (int i = 0; i < la; ++i) {
    complex<double>* const target_array0 = sigma->element_ptr(0, i);
    const complex<double>* const source_array0 = cc->element_ptr(0, i);
    for (int ip = 0; ip != ij; ++ip) {
      const double h = jop->mo1e(ip);
      for (auto& iter : cc->det()->phib(ip)) {
        const double hc = h * iter.sign;
        target_array0[iter.target] += hc * source_array0[iter.source];
      }
    }
  }
}

void ZKnowlesHandy::sigma_2b(shared_ptr<ZDvec> d, shared_ptr<ZDvec> e, shared_ptr<const MOFile> jop) const {
  const int la = d->lena();
  const int lb = d->lenb();
  const int ij = d->ij();
  const int lenab = la*lb;
  //temporary fix to allow for non-complex mo2e.
  //is lenab the correct dimension here?
  complex<double> complex_mo2e[ij*ij];
  for (int j=0;j<(ij);j++) {
    for (int i=0; i<(ij); i++) {
      complex_mo2e[i+j*ij] = complex<double>(jop->mo2e(i,j), 0.0);
    }
  }
//argument 9 was originally jop->mo2e_ptr() for dgemm
  zgemm3m_("n", "n", lenab, ij, ij, 0.5, d->data(), lenab, complex_mo2e, ij, 0.0, e->data(), lenab);
}
