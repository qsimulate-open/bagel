//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: SPCASPT2_gen6.cc
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

#include <src/smith/caspt2/SPCASPT2_tasks6.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::SPCASPT2;

Task250::Task250(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c4 : *range[0])
    for (auto& c2 : *range[0])
      if (t[0]->is_local(c2, c4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c2, c4}}, in, t[0], range));
}

Task251::Task251(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, a1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task252::Task252(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a1 : *range[2])
      if (t[0]->is_local(a1, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a1, a4}}, in, t[0], range));
}

Task253::Task253(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a1 : *range[2])
      if (t[0]->is_local(a1, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a1, a4}}, in, t[0], range));
}

Task254::Task254(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, a1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task255::Task255(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a1 : *range[2])
      if (t[0]->is_local(a1, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a1, a4}}, in, t[0], range));
}

Task256::Task256(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, a1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task257::Task257(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a3 : *range[2])
      if (t[0]->is_local(a3, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a3, a4}}, in, t[0], range));
}

Task258::Task258(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a3 : *range[2])
      if (t[0]->is_local(a3, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a3, a4}}, in, t[0], range));
}

Task259::Task259(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, a1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task260::Task260(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a3 : *range[2])
      if (t[0]->is_local(a3, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a3, a4}}, in, t[0], range));
}

Task261::Task261(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, a1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task262::Task262(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, c2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c2}}, in, t[0], range));
}

Task263::Task263(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, c2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c2}}, in, t[0], range));
}

Task264::Task264(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a3 : *range[2])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, a3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, a3, x0, x1}}, in, t[0], range));
}

Task265::Task265(vector<shared_ptr<Tensor>> t, const bool reset) : reset_(reset) {
  den1_ =  t[0];
}

Task266::Task266(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, x0}}, in, t[0], range));
}

Task267::Task267(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, x0}}, in, t[0], range));
}

Task268::Task268(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, a2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, a2}}, in, t[0], range));
}

Task269::Task269(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      if (t[0]->is_local(a2, c1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a2, c1}}, in, t[0], range));
}

Task270::Task270(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      if (t[0]->is_local(a2, c1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a2, c1}}, in, t[0], range));
}

Task271::Task271(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, a1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, a1}}, in, t[0], range));
}

Task272::Task272(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      if (t[0]->is_local(a1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a1, x0}}, in, t[0], range));
}

#endif
