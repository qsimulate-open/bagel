//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: spindipole.cc
// Copyright (C) 2015 Toru Shiozaki
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


#include <src/mat1e/spindipole.h>
#include <src/integral/rys/spindipolebatch.h>
#include <src/util/atommap.h>

using namespace std;
using namespace bagel;

const static AtomMap atommap;

SpinDipole::SpinDipole(shared_ptr<const Molecule> mol, shared_ptr<const Atom> atom, const int s) : Matrix1eArray<6>(mol), atom_(atom) {

  init(mol);
  fill_upper();

  const double sfac = s ? 1.0/(s*0.5) : 0.0;

  scale(atommap.hfcc_pfac(atom->name())/2.0 * sfac);
}


void SpinDipole::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Molecule>) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();
  SpinDipoleBatch sd(input, atom_);
  sd.compute();

  for (int i = 0; i < Nblocks(); ++i) {
    matrices_[i]->copy_block(offsetb1, offsetb0, dimb1, dimb0, sd.data(i));
  }
}
