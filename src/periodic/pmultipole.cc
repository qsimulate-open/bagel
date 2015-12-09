//
// BAGEL - Parallel electron correlation program.
// Filename: pmultipole.cc
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


#include <src/periodic/pmultipole.h>
#include <src/periodic/multipolebatch.h>

using namespace std;
using namespace bagel;

PMultipole::PMultipole(const shared_ptr<const Lattice> lattice, const int lmax)
 : PMatrix1eArray<64>(lattice), lmax_(lmax), centre_(lattice->centre()) {

  num_multipoles_ = (lmax + 1) * (lmax + 1);
  assert(num_multipoles_ <= Nblocks());
  init(lattice);
}

void PMultipole::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1,
                              shared_ptr<const Lattice> lattice, const int block) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();

  MultipoleBatch mpole(input, centre_, lmax_);
  mpole.compute();

  for (int i = 0; i != num_multipoles_; ++i) {
    ZMatrix mat(dimb1, dimb0);
    mat.copy_block(0, 0, dimb1, dimb0, mpole.data(i));

    (*pdata_blocks_[i])[block]->add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, mat);
  }
}
