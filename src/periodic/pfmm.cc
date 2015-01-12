//
// BAGEL - Parallel electron correlation program.
// Filename: pfmm.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/periodic/pfmm.h>

using namespace std;
using namespace bagel;

PFMM::PFMM(shared_ptr<const Lattice> lattice, const int lmax) : lattice_(lattice), lmax_(lmax) {

  assert(lmax_ <= ANG_HRR_END);
  num_multipoles_ = (lmax_ + 1) * (lmax_ + 1);
}


vector<vector<complex<double>>> PFMM::multipoles(const shared_ptr<const PData> density) {

  vector<vector<complex<double>>> out(num_multipoles_);
  const int ncell = density->nblock();
  const int nbasis = density->blocksize();

  auto pm = make_shared<const PMultipole>(lattice_, lmax_);
  for (int n = 0; n != num_multipoles_; ++n) {
    PData olm = (*pm)[n];
    vector<complex<double>> olm_n(ncell);
    for (int i = 0; i != ncell; ++i) {
      for (int j = 0; j != nbasis; ++j) {
        for (int k = 0; k != nbasis; ++k)
          olm_n[i] += *(olm(i)->data()+k+j*nbasis) * *((*density)(i)->data()+k+j*nbasis);
      }
    }
    out[n] = olm_n;
  }

  return out;
}


void PFMM::print_multipoles(vector<complex<double>> multipoles) const {
  cout << "LMAX = " << lmax_ << endl;
  int cnt = 0;
  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m, ++cnt)
      cout << setprecision(9) << multipoles[cnt] << "   ";
    cout << endl;
  }
}
