//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: breitint.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/mat1e/rel/breitint.h>
#include <src/util/alpha.h>
#include <src/integral/rys/breitbatch.h>

using namespace std;
using namespace bagel;

BreitInt::BreitInt(const shared_ptr<const Molecule> mol) : Matrix1eArray<6>(mol->naux(), mol->naux()) {

  init(mol);

  const vector<int> xyz = { Comp::X, Comp::Y, Comp::Z };
  for (auto& i : xyz) {
    for (auto& j : xyz)
      if (i <= j)
        index_.push_back({i,j});
  }

  localize();
  assert(matrices_.size() == index_.size());
}



void BreitInt::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Molecule>) {

  auto dum0 = make_shared<const Shell>(input[0]->spherical());
  auto dum1 = make_shared<const Shell>(input[1]->spherical());

  // input = [b1, b0]
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();
  const array<shared_ptr<const Shell>,4> in{{input[0], dum0, input[1], dum1}};
  BreitBatch batch(in, 1.0);
  batch.compute();

  for (int i = 0; i != Nblocks(); ++i)
    matrices_[i]->copy_block(offsetb1, offsetb0, dimb1, dimb0, batch.data(i));

}


void BreitInt::init(shared_ptr<const Molecule> mol) {

  list<shared_ptr<const Shell>> shells;
  for (auto& i : mol->aux_atoms())
    shells.insert(shells.end(), i->shells().begin(), i->shells().end());

  // TODO thread, parallel
  int o0 = 0;
  for (auto& a0 : shells) {
    int o1 = 0;
    for (auto& a1 : shells) {
      array<shared_ptr<const Shell>,2> input = {{a1, a0}};
      computebatch(input, o0, o1, nullptr);
      o1 += a1->nbasis();
    }
    o0 += a0->nbasis();
  }

}
