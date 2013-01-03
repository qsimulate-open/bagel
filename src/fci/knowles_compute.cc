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
#include <src/fci/fci.h>
#include <src/fci/knowles.h>
#include <src/util/davidson.h>
#include <src/util/constants.h>
#include <vector>

// toggle for timing print out.
static const bool tprint = false;

using namespace std;
using namespace std::chrono;
using namespace bagel;

shared_ptr<Dvec> KnowlesHandy::form_sigma(shared_ptr<const Dvec> ccvec, shared_ptr<const MOFile> jop,
                     const vector<int>& conv) const { // d and e are scratch area for D and E intermediates 

  const int ij = nij(); 

  const int nstate = ccvec->ij();
  shared_ptr<Dvec> sigmavec(new Dvec(ccvec->det(), nstate));
  sigmavec->zero();

  // we need two vectors for intermediate quantities
  shared_ptr<Dvec> d(new Dvec(ccvec->det(), ij));
  shared_ptr<Dvec> e(new Dvec(ccvec->det(), ij));


  for (int istate = 0; istate != nstate; ++istate) {
    if (conv[istate]) continue;
    shared_ptr<const Civec> cc = ccvec->data(istate);  
    shared_ptr<Civec> sigma = sigmavec->data(istate);  

    vector<pair<string, double> > timing;
    auto start = high_resolution_clock::now();

    // (task1) one-electron alpha: sigma(Psib, Psi'a) += sign h'(ij) C(Psib, Psia) 
    sigma_1(cc, sigma, jop);
    if (tprint) print_timing_("task1", start, timing);

    // (task2) two electron contributions
    d->zero();

    // step (c) (task2a-1) D(Phib, Phia, ij) += sign C(Psib, Phi'a)
    sigma_2a1(cc, d);
    if (tprint) print_timing_("task2a-1", start, timing);

    // step (d) (task2a-2) D(Phib, Phia, ij) += sign C(Psib', Phia)
    sigma_2a2(cc, d);
    if (tprint) print_timing_("task2a-2", start, timing);

    // step (e) (task2b) E(Phib, Phia, kl) = D(Psib, Phia, ij) (ij|kl)
    sigma_2b(d, e, jop);
    if (tprint) print_timing_("task2b", start, timing);

    // step (f) (task2c-1) sigma(Phib, Phia') += sign E(Psib, Phia, kl)
    sigma_2c1(sigma, e);
    if (tprint) print_timing_("task2c-1", start, timing);

    // step (g) (task2c-2) sigma(Phib', Phia) += sign E(Psib, Phia, kl)
    sigma_2c2(sigma, e);
    if (tprint) print_timing_("task2c-2", start, timing);

    // (task3) one-electron beta: sigma(Psib', Psia) += sign h'(ij) C(Psib, Psia)
    sigma_3(cc, sigma, jop);

    if (tprint) {
      print_timing_("task3", start, timing);
      cout << "     timing info" << endl;
      for (auto& iter : timing)
        cout << "    " << setw(10) << iter.first << setw(10) << setprecision(2) << iter.second << endl;
    }
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
      daxpy_(lb, sign, source_base + iter.target*lb, 1, target_array, 1);
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

void KnowlesHandy::sigma_2c1(shared_ptr<Civec> sigma, shared_ptr<const Dvec> e) const {
  const int lb = e->lenb();
  const int ij = e->ij();
  for (int ip = 0; ip != ij; ++ip) { 
    const double* const source_base = e->data(ip)->data();
    for (auto& iter : e->det()->phia(ip)) {
      const double sign = static_cast<double>(iter.sign); 
      double* const target_array = sigma->element_ptr(0, iter.target);
      daxpy_(lb, sign, source_base + lb*iter.source, 1, target_array, 1);
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

void KnowlesHandy::sigma_2b(shared_ptr<Dvec> d, shared_ptr<Dvec> e, shared_ptr<const MOFile> jop) const {
  const int la = d->lena();
  const int lb = d->lenb();
  const int ij = d->ij();
  const int lenab = la*lb;
  dgemm_("n", "n", lenab, ij, ij, 0.5, d->data(), lenab, jop->mo2e_ptr(), ij,
                                  0.0, e->data(), lenab);
}

