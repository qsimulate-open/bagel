//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cap.cc
// Copyright (C) 2018 Toru Shiozaki
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

#include <src/mat1e/cap.h>
#include <src/scf/ks/dftgrid.h>

using namespace std;
using namespace bagel;

CAP::CAP(shared_ptr<const Molecule> mol) : Matrix1e(mol), mol_(mol) {

}


void CAP::compute() {
#if 1
  // Original Becke grid. nrad is abitrary, nang should be the value in src/scf/ks/lebedevlist.h 
  BLGrid grid(/*nrad*/ 100, /*nang*/ 302, mol_);
#else
  // standard DFT grid
  DefaultGrid grid(mol_);
#endif

  shared_ptr<const Matrix> ao = grid.grid()->basis();
  shared_ptr<const Matrix> gp = grid.grid()->data(); // x,y,z,weight
  shared_ptr<Matrix> aow = ao->copy();

  const array<double,3> origin = mol_->charge_center();
  const array<double,3> size = mol_->cap();

  auto r = [](const double r0, const double size) {
    const double ra = fabs(r0);
    const double onset = size * 0.5;
    return ra < onset ? 0.0 : ra - onset;
  };

  for (size_t i = 0; i != ao->mdim(); ++i) {
    const double x = gp->element(0, i) - origin[0];
    const double y = gp->element(1, i) - origin[1];
    const double z = gp->element(2, i) - origin[2];
    const double w = gp->element(3, i);

    // here we calculate the CAP amplitude at this grid point
    // omega is set to one 
    const double cap = pow(r(x, size[0]), 2) + pow(r(y, size[1]), 2) + pow(r(z, size[2]), 2);
    blas::scale_n(cap * w, aow->element_ptr(0, i), aow->ndim());
  }

  *static_cast<Matrix*>(this) = *ao ^ *aow;
}

