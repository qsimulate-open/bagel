//
// BAGEL - Parallel electron correlation program.
// Filename: zknowles_phase.cc
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>>
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

#include <cmath>
#include <random>
#include <src/zfci/zknowles.h>

using namespace std;
using namespace bagel;

void ZKnowlesHandy::mult_phase_factor() {

  // (1) implement ZMOFile and replace MOFile in zfci.h with ZMOFile (and check)
  // (2) these two lines should be in complex
  // (3) think about compression
  // (4) make a matrix like below


  ZMatrix phase(norb_, norb_);
  for (int i = 0; i != norb_; ++i) {
    const double ran = rand();
    const complex<double> fac(cos(ran), sin(ran));
    phase(i,i) = fac;
  }
  const int ij = nij();
  auto mo1e = make_shared<ZMatrix>(norb_, norb_);
  unique_ptr<complex<double>[]> mo2e(new complex<double>[ij*ij]);

//transform 1e integrals.
//TODO check indices
  for (int i = 0; i != norb_; i++) {
    for (int j = 0; j != norb_; j++) {
      mo1e->element(i, j) = jop_->mo1e(j, i);
    }
  }
  *mo1e = phase % *mo1e * phase;

  mo1e->norm();

  //iterating over subordinate norb_, norb_ matrices of mo2e
  for (int i = 0; i != norb_; ++i) {
    for (int j = 0; j != norb_; ++j) {
      unique_ptr<complex<double>[]> in(new complex<double>[ij]);
      unique_ptr<complex<double>[]> tmp(new complex<double>[ij]);
      unique_ptr<complex<double>[]> out(new complex<double>[ij]);

      //assigning array 'in' by iterating over k and l to produce and norb_*norb_ array
      for (int k = 0; k != norb_; ++k) {
        for (int l = 0; l != norb_; ++l) {
          in[l+norb_*k] = jop_->mo2e_kh(i, j, k, l);
        }
      }

      // in*U
      zgemm3m_("n","n", norb_, norb_, norb_, 1.0, in.get(), norb_, phase.data(), norb_, 1.0, tmp.get(), norb_);
      // UH * (in*U)
      zgemm3m_("c","n", norb_, norb_, norb_, 1.0, phase.data(), norb_, tmp.get(), norb_, 1.0, out.get(), norb_);

      // assigning transformed arrays to positions in mo2e identical to the original positions in jop_->mo2e
      for (int k = 0; k!= norb_; ++k) {
        for (int l = 0; l != norb_; ++l) {
          mo2e[l + norb_*k + norb_*norb_*j + norb_*norb_*norb_*i] = out[l+norb_*k];
        }
      }

    }
  }

  complex<double> norm = zdotc_(ij*ij, mo2e, 1, mo2e, 1);
  assert(abs(norm.imag())<1e-8);

//apply transformation to hamiltonian
  jop_ =  make_shared<ZHtilde>(ref_, 0, 0, mo1e, move(mo2e));

}
