//
// BAGEL - Parallel electron correlation program.
// Filename: sohcore_base.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/scf/sohcore_base.h>
#include <src/integral/os/kineticbatch.h>
#include <src/integral/os/mmbatch.h>
#include <src/integral/rys/sonaibatch.h>

using namespace std;
using namespace bagel;

SOHcore_base::SOHcore_base(const shared_ptr<const Molecule> mol) : Matrix1e(mol) {

  init();
  fill_upper();

}


void SOHcore_base::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();

  {
    KineticBatch kinetic(input);
    kinetic.compute();

    copy_block(offsetb1, offsetb0, dimb1, dimb0, kinetic.data());
  }
  {
    SONAIBatch nai(input, mol_);
    nai.compute();

    add_block(1.0, offsetb1, offsetb0, dimb1, dimb0, nai.data());
  }

  if (mol_->external()) {
    DipoleBatch dipole(input, mol_);
    dipole.compute();
    const size_t block = dipole.size_block();
    const double* dip = dipole.data();

    int cnt = 0;
    for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
      for (int j = offsetb1; j != dimb1 + offsetb1; ++j, ++cnt) {
        data_[i*ndim_ + j] += dip[cnt        ]*mol_->external(0);
        data_[i*ndim_ + j] += dip[cnt+block  ]*mol_->external(1);
        data_[i*ndim_ + j] += dip[cnt+block*2]*mol_->external(2);
      }
    }
  }
}


