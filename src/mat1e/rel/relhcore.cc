//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relhcore.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/util/constants.h>
#include <src/mat1e/rel/relhcore.h>
#include <src/mat1e/rel/reldipole.h>

using namespace std;
using namespace bagel;

void RelHcore::compute_() {
  const int n = geom_->nbasis();

  Matrix nai(2*n, 2*n);
  nai.copy_block(0, 0, n, n, nai_);
  nai.copy_block(n, n, n, n, nai_);

  Matrix kinetic(2*n, 2*n);
  kinetic.copy_block(0, 0, n, n, kinetic_);
  kinetic.copy_block(n, n, n, n, kinetic_);

  const complex<double> w(0.25/(c__*c__));
  const complex<double> wi(0.0, w.real());
  ZMatrix zsnai(2*n, 2*n);
  zsnai.add_real_block(  w, 0, 0, n, n, (*smallnai_)[0]);
  zsnai.add_real_block(  w, n, n, n, n, (*smallnai_)[0]);
  zsnai.add_real_block( wi, 0, 0, n, n, (*smallnai_)[1]);
  zsnai.add_real_block(-wi, n, n, n, n, (*smallnai_)[1]);
  zsnai.add_real_block( wi, 0, n, n, n, (*smallnai_)[2]);
  zsnai.add_real_block( wi, n, 0, n, n, (*smallnai_)[2]);
  zsnai.add_real_block(  w, 0, n, n, n, (*smallnai_)[3]);
  zsnai.add_real_block( -w, n, 0, n, n, (*smallnai_)[3]);

  // RKB hcore: T is off-diagonal block matrices, V is first main diagonal, and 1/4m^2c^2W-T is second main diagonal
  copy_real_block(1.0,  0,   0, 2*n, 2*n, nai);
  copy_real_block(1.0,  0, 2*n, 2*n, 2*n, kinetic);
  copy_real_block(1.0,2*n,   0, 2*n, 2*n, kinetic);
  copy_block(2*n, 2*n, 2*n, 2*n, zsnai);
  add_real_block(-1.0, 2*n, 2*n, 2*n, 2*n, kinetic);

  if (geom_->external()) {
    RelDipole dipole(geom_);
    array<shared_ptr<const ZMatrix>,3> mat = dipole.compute_matrices();
    ax_plus_y(geom_->external(0), mat[0]);
    ax_plus_y(geom_->external(1), mat[1]);
    ax_plus_y(geom_->external(2), mat[2]);
  }
}

