//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_gen2.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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


#include <src/smith/CASPT2_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

Task50::Task50(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, a2, c1}}, in, t[0], range));
}

Task51::Task51(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, a2, c1}}, in, t[0], range));
}

