//
// BAGEL - Parallel electron correlation program.
// Filename: zknowles_denom.cc
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>
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


#include <iomanip>
#include <stdexcept>
#include <src/zfci/zknowles.h>
#include <src/zfci/relmofile.h>
#include <src/zfci/zmofile.h>
#include <src/integral/rys/eribatch.h>
#include <src/util/combination.hpp>
#include <iostream>

using namespace std;
using namespace bagel;

//
// averaged diagonal elements as defined in Knowles & Handy (1989) Compt. Phys. Comm.
//
void ZKnowlesHandy::const_denom() {
  vector<complex<double>> jop, kop, fk;
  jop.resize(norb_*norb_);
  kop.resize(norb_*norb_);
  fk.resize(norb_);
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      jop[i*norb_+j] = 0.5*jop_->mo2e_kh(j, j, i, i);
      jop[j*norb_+i] = 0.5*jop_->mo2e_kh(i, i, j, j);
    }
  }
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j <= i; ++j) {
      kop[i*norb_+j] = 0.5*jop_->mo2e_kh(j, i, i, j);
      kop[j*norb_+i] = 0.5*jop_->mo2e_kh(i, j, j, i);
    }
  }
  for (int i = 0; i != norb_; ++i) {
    fk[i] = 0.0;
    for (int j = 0; j != norb_; ++j) {
      fk[i] += kop[j*norb_+i];
    }
  }
  denom_ = make_shared<Civec>(det());
  const int nspin = det()->nspin();
  const int nspin2 = nspin*nspin;
  complex<double> med = 0.0;
  double* iter = denom_->data();
  for (auto& ia : det()->stringa()) {
    for (auto& ib : det()->stringb()) {
      const int nopen = (ia^ib).count();
      const double F = (nopen >> 1) ? (static_cast<double>(nspin2 - nopen)/(nopen*(nopen-1))) : 0.0;
      *iter = 0.0;
      med = 0.0;
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

          med += (jop[j+norb_*i]+jop[i+norb_*j]) * static_cast<double>(addj) - (0.5 * (kop[j+norb_*i]+kop[i+norb_*j]) * (F*Ni*Nj + addj));
        }
        med += (jop_->mo1e(i,i) + fk[i]) * static_cast<double>(niab) - kop[i+norb_*i] * 0.5 * static_cast<double>(Ni - niab*niab);
      }
      assert(fabs(med.imag())<1e-8);
      *iter = med.real();
      ++iter;
    }
  }
}
void ZKnowlesHandy::update(shared_ptr<const Coeff> c) {
  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  Timer timer;
  jop_ = make_shared<ZJop>(ref_, ncore_, ncore_+norb_, c, "KH");

  // right now full basis is used.
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;
  const_denom();
  mult_phase_factor();
}

void ZKnowlesHandy::relupdate() {
  // iiii file to be created (MO transformation).
  // now jop_->mo1e() and jop_->mo2e() contains one and two body part of Hamiltonian
  Timer timer;
  auto jop_ = make_shared<RelJop>(ref_, ncore_, ncore_+norb_, "KH");

  // right now full basis is used.
  cout << "    * Integral transformation done. Elapsed time: " << setprecision(2) << timer.tick() << endl << endl;

  throw logic_error("testing...");
  const_denom();
}
