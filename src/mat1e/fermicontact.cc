//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: fermicontact.cc
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


#include <src/mat1e/fermicontact.h>
#include <src/util/constants.h>
#include <src/util/atommap.h>

using namespace std;
using namespace bagel;

const static AtomMap atommap;

FermiContact::FermiContact(shared_ptr<const Molecule> mol, shared_ptr<const Atom> atom, const int s) : Matrix1e(mol), position_(atom->position()) {

  init(mol);
  fill_upper();

  const double sfac = s ? 1.0/(s*0.5) : 0.0;

  scale(atommap.hfcc_pfac(atom->name())*4.0*pi__/3 * sfac);
}

void FermiContact::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Molecule> dummy) {

  // input = [b1, b0]
  assert(input.size() == 2);
  shared_ptr<const Shell> s1 = input[0];
  shared_ptr<const Shell> s0 = input[1];
  const int dimb1 = s1->nbasis();
  const int dimb0 = s0->nbasis();

  // obtain the value of the basis functions at the nucleus
  VectorB tmp0(dimb0);
  VectorB tmp1(dimb1);

  s0->compute_grid_value(tmp0.data(), nullptr, nullptr, nullptr, position_[0]-s0->position(0), position_[1]-s0->position(1), position_[2]-s0->position(2));
  s1->compute_grid_value(tmp1.data(), nullptr, nullptr, nullptr, position_[0]-s1->position(0), position_[1]-s1->position(1), position_[2]-s1->position(2));

  for (int i = offsetb0; i != dimb0 + offsetb0; ++i)
    for (int j = offsetb1; j != dimb1 + offsetb1; ++j)
      element(j, i) = tmp1[j-offsetb1] * tmp0[i-offsetb0];
}


