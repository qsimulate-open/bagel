//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_gen13.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/caspt2/MSCASPT2_tasks13.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

Task600::Task600(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[5]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, out, range));
}

Task601::Task601(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task602::Task602(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task603::Task603(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task604::Task604(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task605::Task605(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x1 : *range[1])
          if (t[5]->is_local(x3, x0, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, out, range));
}

Task606::Task606(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x0, x1}}, in, t[0], range));
}

Task607::Task607(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x0, x1}}, in, t[0], range));
}

#endif
