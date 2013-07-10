//
// BAGEL - Parallel electron correlation program.
// Filename: reloverlap.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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


#include <src/util/constants.h>
#include <src/rel/reloverlap.h>

using namespace std;
using namespace bagel;

void RelOverlap::compute_() {
  const int n = mol_->nbasis();
  const complex<double> coeff1 (1.0, 0.0);

  auto out = make_shared<ZMatrix>(4*n, 4*n);
  auto ovl = make_shared<Overlap>(*overlap_);
  auto k12 = make_shared<Matrix>(*kinetic_);

  if (half_inverse_) {
    ovl->inverse_half();
    k12->inverse_half();
    k12->scale(c__/sqrt(0.5));
  } else {
    k12->scale(0.5/(c__*c__));
  }


  copy_real_block(coeff1, 0, 0, n, n, ovl);
  copy_real_block(coeff1, n, n, n, n, ovl);
  copy_real_block(coeff1, 2*n, 2*n, n, n, k12);
  copy_real_block(coeff1, 3*n, 3*n, n, n, k12);
}

