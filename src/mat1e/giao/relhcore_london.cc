//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: relhcore_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#include <src/mat1e/giao/relhcore_london.h>

using namespace std;
using namespace bagel;

void RelHcore_London::compute_() {
  const int n = geom_->nbasis();

  ZMatrix nai(2*n, 2*n);
  nai.copy_block(0, 0, n, n, nai_);
  nai.copy_block(n, n, n, n, nai_);

  ZMatrix kinetic(2*n, 2*n);
  kinetic.copy_block(0, 0, n, n, kinetic_);
  kinetic.copy_block(n, n, n, n, kinetic_);

  const complex<double> w(0.25/(c__*c__));
  const complex<double> wi(0.0, w.real());
  ZMatrix zsnai(2*n, 2*n);
  zsnai.add_block(  w, 0, 0, n, n, (*smallnai_)[0]);
  zsnai.add_block(  w, n, n, n, n, (*smallnai_)[0]);
  zsnai.add_block( wi, 0, 0, n, n, (*smallnai_)[1]);
  zsnai.add_block(-wi, n, n, n, n, (*smallnai_)[1]);
  zsnai.add_block( wi, 0, n, n, n, (*smallnai_)[2]);
  zsnai.add_block( wi, n, 0, n, n, (*smallnai_)[2]);
  zsnai.add_block(  w, 0, n, n, n, (*smallnai_)[3]);
  zsnai.add_block( -w, n, 0, n, n, (*smallnai_)[3]);

  const complex<double> rh (-0.5);
  const complex<double> ih (0.0, rh.real());
  ZMatrix zeeman(2*n, 2*n);

  zeeman.add_block( rh*geom_->magnetic_field(2), 0, 0, n, n, overlap_);
  zeeman.add_block(-rh*geom_->magnetic_field(2), n, n, n, n, overlap_);
  zeeman.add_block( rh*geom_->magnetic_field(0), 0, n, n, n, overlap_);
  zeeman.add_block( rh*geom_->magnetic_field(0), n, 0, n, n, overlap_);
  zeeman.add_block(-ih*geom_->magnetic_field(1), 0, n, n, n, overlap_);
  zeeman.add_block( ih*geom_->magnetic_field(1), n, 0, n, n, overlap_);

  // RMB hcore: all 1-electron integrals over magnetically balanced 2-spinors
  zero();
  copy_block(  0,   0, 2*n, 2*n, nai);
  copy_block(  0, 2*n, 2*n, 2*n, kinetic);
  copy_block(2*n,   0, 2*n, 2*n, kinetic);
  copy_block(2*n, 2*n, 2*n, 2*n, zsnai);
  add_block(-1.0, 2*n, 2*n, 2*n, 2*n, kinetic);
  add_block(-1.0,   0, 2*n, 2*n, 2*n, zeeman);
  add_block(-1.0, 2*n,   0, 2*n, 2*n, zeeman);
  add_block( 1.0, 2*n, 2*n, 2*n, 2*n, zeeman);

}

