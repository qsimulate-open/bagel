//
// BAGEL - Parallel electron correlation program.
// Filename: hcore.cc
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


#include <src/scf/hcore.h>
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

typedef shared_ptr<const Geometry> RefGeometry;
typedef shared_ptr<const Shell> RefShell;

Hcore::Hcore(const RefGeometry geom) : Matrix1e(geom) {

  init();
  fill_upper();

}


Hcore::~Hcore() {

}


void Hcore::computebatch(const array<RefShell,2>& input, const int offsetb0, const int offsetb1) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();

  {
    KineticBatch kinetic(input);
    kinetic.compute();
    const double* kdata = kinetic.data();
    int cnt = 0;
    for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
      for (int j = offsetb1; j != dimb1 + offsetb1; ++j, ++cnt) {
        data_[i*nbasis_ + j] = kdata[cnt];
      }
    }
  }
  {
    NAIBatch nai(input, geom_);
    nai.compute();
    const double* ndata = nai.data();
    int cnt = 0;
    for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
      for (int j = offsetb1; j != dimb1 + offsetb1; ++j, ++cnt) {
        data_[i*nbasis_ + j] += ndata[cnt];
      }
    }
  }

  if (geom_->external()) {
    DipoleBatch dipole(input, geom_->charge_center());
    dipole.compute();
    const size_t block = dipole.size_block();
    const double* dip = dipole.data();
    int cnt = 0;
    for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
      for (int j = offsetb1; j != dimb1 + offsetb1; ++j, ++cnt) {
        data_[i*nbasis_ + j] += dip[cnt        ]*geom_->external(0);
        data_[i*nbasis_ + j] += dip[cnt+block  ]*geom_->external(1);
        data_[i*nbasis_ + j] += dip[cnt+block*2]*geom_->external(2);
      }
    }
  }
}


