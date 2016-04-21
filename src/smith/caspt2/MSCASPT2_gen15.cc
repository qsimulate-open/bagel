//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_gen15.cc
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

#include <src/smith/caspt2/MSCASPT2_tasks15.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

Task700::Task700(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task701::Task701(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task702::Task702(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task703::Task703(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, c2, a1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, c2, a1, x0}}, in, t[0], range));
}

Task704::Task704(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task705::Task705(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, a2, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, a2, c1}}, in, t[0], range));
}

Task706::Task706(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task707::Task707(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, x3, x2}}, in, t[0], range));
}

Task708::Task708(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, x3, x2}}, in, t[0], range));
}

Task709::Task709(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, x3, x2}}, in, t[0], range));
}

Task710::Task710(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, x3, x2}}, in, t[0], range));
}

Task711::Task711(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task712::Task712(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x3 : *range[1])
              if (t[0]->is_local(x3, x4, x5, x0, x1, x2))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x3, x4, x5, x0, x1, x2}}, in, t[0], range));
}

Task713::Task713(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task714::Task714(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task715::Task715(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task716::Task716(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task717::Task717(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task718::Task718(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task719::Task719(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x3 : *range[1])
              if (t[0]->is_local(x3, x4, x5, x0, x1, x2))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x3, x4, x5, x0, x1, x2}}, in, t[0], range));
}

#endif
