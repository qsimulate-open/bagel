//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ras/rasci_denom.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/ci/ras/rasci.h>
#include <src/ci/ras/denomtask.h>

using namespace std;
using namespace bagel;

void RAS::DenomTask::compute() {
  const int nspin = abit_.count() - stringb_->nele();
  const int nspin2 = nspin*nspin;
  const int norb = stringb_->norb();

  double* iter = data_;
  const bitset<nbit__> ia = abit_;
  for (auto& ib : *stringb_) {
    const int nopen = (ia^ib).count();
    const double F = (nopen >> 1) ? (static_cast<double>(nspin2 - nopen)/(nopen*(nopen-1))) : 0.0;
    *iter = 0.0;
    for (int i = 0; i < norb; ++i) {
      const int nia = ia[i];
      const int nib = ib[i];
      const int niab = (nia + nib);
      const int Ni = (nia ^ nib);
      for (int j = 0; j < i; ++j) {
        const int nja = ia[j];
        const int njb = ib[j];
        const int Nj = nja ^ njb;
        const int addj = niab * (nja + njb);

        *iter += jop_[j+norb*i] * 2.0 * addj - kop_[j+norb*i] * (F*Ni*Nj + addj);
      }
      *iter += h_[i] * niab - kop_[i+norb*i] * 0.5 * (Ni - niab*niab);
    }
    ++iter;
  }
}

void RASCI::const_denom() {
  Timer denom_t;
  unique_ptr<double[]> h(new double[norb_]);
  unique_ptr<double[]> jop(new double[norb_*norb_]);
  unique_ptr<double[]> kop(new double[norb_*norb_]);

  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop[i*norb_+j] = jop[j*norb_+i] = 0.5*jop_->mo2e_hz(j, i, j, i);
      kop[i*norb_+j] = kop[j*norb_+i] = 0.5*jop_->mo2e_hz(j, i, i, j);
    }
    h[i] = jop_->mo1e(i,i);
  }
  denom_t.tick_print("jop, kop");

  denom_ = make_shared<RASCivec>(det());

  size_t tasksize = 0;
  for (auto& iblock : denom_->blocks()) { if (iblock) tasksize += iblock->lena(); }
  TaskQueue<RAS::DenomTask> tasks(tasksize);

  for (auto& iblock : denom_->blocks()) {
    if ( !iblock ) continue;
    double* iter = iblock->data();
    for (auto& ia : *iblock->stringsa()) {
      tasks.emplace_back(iter, ia, iblock->stringsb(), jop.get(), kop.get(), h.get());
      iter += iblock->lenb();
    }
  }

  tasks.compute();
  denom_t.tick_print("denom");
}


void RASCI::update(shared_ptr<const Coeff> c) {
  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  Timer timer;
  // Same Jop as used in FCI
  jop_ = make_shared<Jop>(ref_, ncore_, ncore_+norb_, c, /*store*/false, "HZ");

  // right now full basis is used.
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;

  const_denom();
}
