//
// BAGEL - Parallel electron correlation program.
// Filename: breit.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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


#include <stddef.h>
#include <src/rel/breit.h>
#include <src/rysint/breitbatch.h>

using namespace std;
using namespace bagel;

Breit::Breit(const shared_ptr<const Geometry> geom) : NMatrix1e(geom) {

  for (int i = 0; i != nblocks(); ++i)
    matrix_data_.push_back(shared_ptr<Matrix>(new Matrix(geom->naux(), geom->naux())));

  init();

}


void Breit::print() const {
  int j = 0;
  for (auto& i : matrix_data_) {
    stringstream ss;
    ss << "Breit " << j++;
    i->print(ss.str());
  }
}


void Breit::computebatch(const array<shared_ptr<const Shell>,4>& input, const int offsetb0, const int offsetb1) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[2]->nbasis();
  assert(input[1]->dummy() && input[3]->dummy());
  BreitBatch batch(input, 1.0);
  batch.compute();

  for (int i = 0; i != nblocks(); ++i)
    matrix_data_[i]->copy_block(offsetb1, offsetb0, dimb1, dimb0, batch.data(i));

}


void Breit::init() {

  list<shared_ptr<const Shell>> shells;
  for (auto& i : geom_->aux_atoms())
    shells.insert(shells.end(), i->shells().begin(), i->shells().end());

  const shared_ptr<const Shell> dum(new Shell(shells.front()->spherical()));

  // TODO thread, parallel
  int o0 = 0;
  for (auto& a0 : shells) {
    int o1 = 0;
    for (auto& a1 : shells) {
      array<shared_ptr<const Shell>,4> input = {{a1, dum, a0, dum}};
      computebatch(input, o0, o1);
      o1 += a1->nbasis();
    }
    o0 += a0->nbasis();
  }

}
