//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks16.h
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

#include <src/smith/relmrci/RelMRCI_tasks16.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

Task750::Task750(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, a3, c2, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c2, a1}}, in, t[0], range));
}

Task751::Task751(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, a3, c2, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c2, a1}}, in, t[0], range));
}

Task752::Task752(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, a2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a2, x0, a1}}, in, t[0], range));
}

Task753::Task753(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, a2, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, a2, x0, x1}}, in, t[0], range));
}

Task754::Task754(vector<shared_ptr<Tensor>> t, const bool reset) : reset_(reset) {
  n_ =  t[0];
}

Task755::Task755(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, x1, c1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, c1, x0}}, in, t[0], range));
}

Task756::Task756(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& c1 : *range[0])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, x1, c1, c2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, c1, c2}}, in, t[0], range));
}

Task757::Task757(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x2, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x2, x0, x1}}, in, t[0], range));
}

Task758::Task758(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x2, x1, x0}}, in, t[0], range));
}

Task759::Task759(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& x0 : *range[1])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, x0, c1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, x0, c1, a2}}, in, t[0], range));
}

Task760::Task760(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x0}}, in, t[0], range));
}

Task761::Task761(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x0}}, in, t[0], range));
}

Task762::Task762(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, x1, c1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, c1, a2}}, in, t[0], range));
}

Task763::Task763(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, a2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a2, x1, x0}}, in, t[0], range));
}

Task764::Task764(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x2, x0, a1}}, in, t[0], range));
}

Task765::Task765(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, x0, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x2, x1}}, in, t[0], range));
}

Task766::Task766(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a4, c1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a4, c1, a2}}, in, t[0], range));
}

Task767::Task767(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, a4, c3, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}

Task768::Task768(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, a3, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a3, x0, a1}}, in, t[0], range));
}

Task769::Task769(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, a3, c2, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c2, a1}}, in, t[0], range));
}

Task770::Task770(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, a3, c2, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c2, a1}}, in, t[0], range));
}

Task771::Task771(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, a2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a2, x0, a1}}, in, t[0], range));
}

Task772::Task772(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, x1, a1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}

#endif
