//
// BAGEL - Parallel electron correlation program.
// Filename: smallnai.cc
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
#include <src/rel/smallnai.h>
#include <src/integral/rys/smallnaibatch.h>

using namespace std;
using namespace bagel;

SmallNAI::SmallNAI(const shared_ptr<const Geometry> geom) : Matrix1eArray<4>(geom) {

  init();

}


void SmallNAI::print(const string name) const {
  Matrix1eArray<4>::print(name.empty() ? "SmallNAI" : name);
}


void SmallNAI::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb1 = input[0]->nbasis();
  const int dimb0 = input[1]->nbasis();
  SmallNAIBatch batch(input, geom_);
  batch.compute();

  for (int i = 0; i != nblocks(); ++i)
    matrices_[i]->copy_block(offsetb1, offsetb0, dimb1, dimb0, batch[i]);

}


void SmallNAI::init() {

  list<shared_ptr<const Shell>> shells;
  for (auto& i : geom_->atoms())
    shells.insert(shells.end(), i->shells().begin(), i->shells().end());

  // TODO thread, parallel
  int o0 = 0;
  for (auto& a0 : shells) {
    int o1 = 0;
    for (auto& a1 : shells) {
      array<shared_ptr<const Shell>,2> input = {{a1, a0}};
      computebatch(input, o0, o1);
      o1 += a1->nbasis();
    }
    o0 += a0->nbasis();
  }

}

