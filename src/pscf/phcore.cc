//
// BAGEL - Parallel electron correlation program.
// Filename: phcore.cc
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


#include <src/pscf/phcore.h>
#include <src/osint/kineticbatch.h>
#include <src/rysint/naibatch.h>
#include <iostream>

using namespace std;
using namespace bagel;

typedef shared_ptr<PGeometry> RefPGeometry;

PHcore::PHcore(const RefPGeometry g, const bool cabs, const bool kinetic_only)
: PMatrix1e(g), cabs_(cabs), kinetic_only_(kinetic_only) {

  init();

}

PHcore::~PHcore() {

}

void PHcore::computebatch(const array<shared_ptr<const Shell>,2>& input,
    const int offsetb0, const int offsetb1,
    const int nbasis, const int blockoffset) {

  // input = [b1, b0]
  // because of the convention in integral codes.
  assert(input.size() == 2);
  const int dimb0 = input[1]->nbasis(); // b0
  const int dimb1 = input[0]->nbasis(); // b1
  KineticBatch kinetic(input);
  kinetic.compute();
  const double* kdata = kinetic.data();
  NAIBatch nai(input, geom_, L(), A());
  nai.compute();
  const double* ndata = nai.data();

  double scale = cabs_ ? 0.5 : 1.0;
  if (kinetic_only_) scale = 0.0;

  int cnt = 0;
  for (int j = offsetb0; j != dimb0 + offsetb0; ++j) { // atoms in unit cell 0
    for (int i = offsetb1; i != dimb1 + offsetb1; ++i, ++cnt) { // atoms in unit cell m
      (*data_)[blockoffset + j * nbasis + i] = static_cast<complex<double> >(kdata[cnt] + ndata[cnt]*scale);
    }
  }

}

