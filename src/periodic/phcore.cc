//
// BAGEL - Parallel electron correlation program.
// Filename: phcore.cc
// Copyright (C) 2014 Toru Shiozaki
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


#include <src/periodic/phcore.h>
#include <src/integral/os/kineticbatch.h>
#include <src/integral/rys/naibatch.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(PHcore)

PHcore::PHcore(const shared_ptr<const Lattice> lattice) : PMatrix1e(lattice) {

  init(lattice);
}

void PHcore::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Lattice> lattice, const int block) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();

  {
    KineticBatch kinetic(input);
    kinetic.compute();
    Matrix k(dimb1, dimb0);
    k.copy_block(0, 0, dimb1, dimb0, kinetic.data());

    pdata_[block]->copy_real_block(1.0, offsetb1, offsetb0, dimb1, dimb0, k);
  }
  {
    auto mol = make_shared<const Geometry>(*(lattice->primitive_cell()), lattice->lattice_vectors(block));
    NAIBatch nai(input, mol);
    nai.compute();
    Matrix n(dimb1, dimb0);
    n.copy_block(0, 0, dimb1, dimb0, nai.data());

    pdata_[block]->add_real_block(1.0, offsetb1, offsetb0, dimb1, dimb0, n);
  }

}


