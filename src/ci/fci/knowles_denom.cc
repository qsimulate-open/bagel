//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fci_init.cc
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


#include <iomanip>
#include <stdexcept>
#include <src/ci/fci/knowles.h>
#include <src/ci/fci/mofile.h>
#include <src/util/combination.hpp>
#include <iostream>

using namespace std;
using namespace bagel;

//
// averaged diagonal elements as defined in Knowles & Handy (1989) Compt. Phys. Comm.
//
void KnowlesHandy::const_denom() {
  vector<double> jop, kop, fk;
  jop.resize(norb_*norb_);
  kop.resize(norb_*norb_);
  fk.resize(norb_);
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop[i*norb_+j] = jop[j*norb_+i] = 0.5*jop_->mo2e_kh(j, j, i, i);
    }
  }
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      kop[i*norb_+j] = kop[j*norb_+i] = 0.5*jop_->mo2e_kh(j, i, j, i);
    }
  }
  for (int i = 0; i != norb_; ++i) {
    fk[i] = 0.0;
    for (int j = 0; j != norb_; ++j) {
      fk[i] += kop[i*norb_+j];
    }
  }

  denom_ = make_shared<Civec>(det());
  const int nspin = det()->nspin();
  const int nspin2 = nspin*nspin;

  double* iter = denom_->data();
  for (auto& ia : det()->string_bits_a()) {
    for (auto& ib : det()->string_bits_b()) {
      const int nopen = (ia^ib).count();
      const double F = (nopen >> 1) ? (static_cast<double>(nspin2 - nopen)/(nopen*(nopen-1))) : 0.0;
      *iter = 0.0;
      for (int i = 0; i != norb_; ++i) {
        const int nia = ia[i];
        const int nib = ib[i];
        const int niab = nia + nib;
        const int Ni = (nia ^ nib);
        for (int j = 0; j != i; ++j) {
          const int nja = ia[j];
          const int njb = ib[j];
          const int Nj = (nja ^ njb);
          const int addj = niab * (nja + njb);
          *iter += jop[j+norb_*i] * 2.0 * addj - kop[j+norb_*i] * (F*Ni*Nj + addj);
        }
        *iter += (jop_->mo1e(i,i) + fk[i]) * niab - kop[i+norb_*i] * 0.5 * (Ni - niab*niab);
      }
      ++iter;
    }
  }
}

void KnowlesHandy::update(shared_ptr<const Matrix> c) {
  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  Timer timer;
  coeff_ = c;
  jop_ = make_shared<Jop>(ref_, ncore_, ncore_+norb_, coeff_, store_half_ints_, "KH");

  // right now full basis is used.
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;

  const_denom();
}
