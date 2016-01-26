//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASPT2_gen5.cc
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

#include <src/smith/relcaspt2/RelCASPT2_tasks5.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelCASPT2;

Task200::Task200(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, a1, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a1, x1, x0}}, in, t[0], range));
}

Task201::Task201(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task202::Task202(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task203::Task203(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task204::Task204(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, a2, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a2, x3, x2}}, in, t[0], range));
}

Task205::Task205(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task206::Task206(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task207::Task207(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task208::Task208(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task209::Task209(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task210::Task210(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task211::Task211(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c2, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c2, x0, x1}}, in, t[0], range));
}

Task212::Task212(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task213::Task213(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& c1 : *range[0])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, c1, a2, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, c1, a2, x0}}, in, t[0], range));
}

Task214::Task214(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& x1 : *range[1])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, x1, c1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, x1, c1, a2}}, in, t[0], range));
}

Task215::Task215(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task216::Task216(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task217::Task217(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task218::Task218(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, a1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x0}}, in, t[0], range));
}

Task219::Task219(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, a3, c2, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a3, c2, a1}}, in, t[0], range));
}

Task220::Task220(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task221::Task221(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

#endif
