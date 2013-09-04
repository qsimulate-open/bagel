//
// BAGEL - Parallel electron correlation program.
// Filename: relhcore.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/util/constants.h>
#include <src/rel/relhcore.h>

using namespace std;
using namespace bagel;

void RelHcore::compute_() {
  const int n = geom_->nbasis();

  auto nai     = make_shared<Matrix>(2*n, 2*n);
  nai->copy_block(0, 0, n, n, nai_);
  nai->copy_block(n, n, n, n, nai_);

  auto kinetic = make_shared<Matrix>(2*n, 2*n);
  kinetic->copy_block(0, 0, n, n, kinetic_);
  kinetic->copy_block(n, n, n, n, kinetic_);

  const complex<double> w(0.25/(c__*c__));
  const complex<double> wi(0.0, w.real());
  auto zsnai = make_shared<ZMatrix>(2*n, 2*n);
  zsnai->add_real_block(  w, 0, 0, n, n, (*smallnai_)[0]);
  zsnai->add_real_block(  w, n, n, n, n, (*smallnai_)[0]);
  zsnai->add_real_block( wi, 0, 0, n, n, (*smallnai_)[1]);
  zsnai->add_real_block(-wi, n, n, n, n, (*smallnai_)[1]);
  zsnai->add_real_block( wi, 0, n, n, n, (*smallnai_)[2]);
  zsnai->add_real_block( wi, n, 0, n, n, (*smallnai_)[2]);
  zsnai->add_real_block(  w, 0, n, n, n, (*smallnai_)[3]);
  zsnai->add_real_block( -w, n, 0, n, n, (*smallnai_)[3]);

  // RKB hcore: T is off diagonal block matrices, V is first main diagonal, and 1/4m^2c^2W-T is second main diagonal
  zero();
  copy_real_block(1.0,  0,   0, 2*n, 2*n, nai);
  copy_real_block(1.0,  0, 2*n, 2*n, 2*n, kinetic);
  copy_real_block(1.0,2*n,   0, 2*n, 2*n, kinetic);
  copy_block(2*n, 2*n, 2*n, 2*n, zsnai);
  add_real_block(-1.0, 2*n, 2*n, 2*n, 2*n, kinetic);

}

