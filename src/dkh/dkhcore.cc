//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dkhcore.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Raymond Wang <yiqunwang2021@u.northwestern.edu> 
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


#include <src/dkh/dkhcore.h>
#include <src/util/math/matrix.h>
#include <src/integral/os/kineticbatch.h>
#include <src/integral/rys/naibatch.h>
#include <src/integral/libint/libint.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(DKHcore)

DKHcore::DKHcore(shared_ptr<const Molecule> mol) {

  init(mol);
  fill_upper();
}

void DKHcore::computebatch0(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1) {

  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();

}

void DKHcore::init(shared_ptr<const Molecule> mol) {
  
  Kinetic kinetic_(mol);
  eig_(mol->nbasis());
  kinetic_.diagonalize(eig());
  NAI nai_(mol);

  //loop computebatch0();

}
