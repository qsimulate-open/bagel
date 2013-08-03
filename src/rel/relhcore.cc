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
  // distributed hcore and overlap
  const int n = geom_->nbasis();

  auto znai = make_shared<ZMatrix>(2*n, 2*n);
  auto zkinetic = make_shared<ZMatrix>(2*n, 2*n);

  array<shared_ptr<ZMatrix>,4> zsmallnai;
  for (auto& i : zsmallnai)
    i = znai->clone();

  const complex<double> coeff1 (1.0, 0.0);
  const complex<double> coeffi (0.0, 1.0);

  znai->copy_real_block(coeff1, 0, 0, n, n, nai_);
  znai->copy_real_block(coeff1, n, n, n, n, nai_);
  zkinetic->copy_real_block(coeff1, 0, 0, n, n, kinetic_);
  zkinetic->copy_real_block(coeff1, n, n, n, n, kinetic_);

  zsmallnai[0]->copy_real_block(coeff1, 0, 0, n, n, (*smallnai_)[0]);
  zsmallnai[0]->copy_real_block(coeff1, n, n, n, n, (*smallnai_)[0]);
  zsmallnai[1]->copy_real_block(coeffi, 0, 0, n, n, (*smallnai_)[1]);
  zsmallnai[1]->copy_real_block(-coeffi, n, n, n, n, (*smallnai_)[1]);
  zsmallnai[2]->copy_real_block(coeffi, 0, n, n, n, (*smallnai_)[2]);
  zsmallnai[2]->copy_real_block(coeffi, n, 0, n, n, (*smallnai_)[2]);
  zsmallnai[3]->copy_real_block(coeff1, 0, n, n, n, (*smallnai_)[3]);
  zsmallnai[3]->copy_real_block(-coeff1, n, 0, n, n, (*smallnai_)[3]);

  auto smallnai = make_shared<ZMatrix>(*zsmallnai[0] + *zsmallnai[1] + *zsmallnai[2] + *zsmallnai[3]);

  // RKB hcore: T is off diagonal block matrices, V is first main diagonal, and 1/4m^2c^2W-T is second main diagonal
  const complex<double> w(0.25/(c__*c__), 0.0);
  zero();
  copy_block(0, 0, 2*n, 2*n, znai);
  copy_block(0, 2*n, 2*n, 2*n, zkinetic);
  copy_block(2*n, 0, 2*n, 2*n, zkinetic);
  copy_block(2*n, 2*n, 2*n, 2*n, make_shared<ZMatrix>(*smallnai * w - *zkinetic));

}

