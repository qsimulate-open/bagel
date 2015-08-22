//
// BAGEL - Parallel electron correlation program.
// Filename: fermicontact.cc
// Copyright (C) 2015 Toru Shiozaki
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


#include <src/mat1e/fermicontact.h>

using namespace std;
using namespace bagel;

FermiContact::FermiContact(shared_ptr<const Molecule> mol, shared_ptr<const Atom> atom) : Matrix1e(mol), position_(atom->position()) {

  init(mol);
  fill_upper();
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


