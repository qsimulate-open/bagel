//
// BAGEL - Parallel electron correlation program.
// Filename: kinetic.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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


#include <src/scf/kinetic.h>
#include <src/grad/gnaibatch.h>
#include <src/osint/kineticbatch.h>
#include <src/osint/dipolebatch.h>
#include <src/rysint/naibatch.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>

using namespace std;
using namespace bagel;


Kinetic::Kinetic(const shared_ptr<const Geometry> geom) : Matrix1e(geom) {

  init();
  fill_upper();

}


Kinetic::~Kinetic() {

}


void Kinetic::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();

  KineticBatch kinetic(input);
  kinetic.compute();

  copy_block(offsetb1, offsetb0, dimb1, dimb0, kinetic.data());
}


