//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_gen18.cc
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

#include <src/smith/mrci/MRCI_tasks18.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

Task850::Task850(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a1, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a1, x3, x0}}, in, t[0], range));
}

Task851::Task851(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a1}}, in, t[0], range));
}

Task852::Task852(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a1, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a1, x3, x0}}, in, t[0], range));
}

Task853::Task853(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a1, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a1, x2, x1}}, in, t[0], range));
}

Task854::Task854(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a1, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a1, x3, x0}}, in, t[0], range));
}

Task855::Task855(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a1, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a1, x3, x0}}, in, t[0], range));
}

Task856::Task856(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a1}}, in, t[0], range));
}

Task857::Task857(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a3 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a3, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a3, x3, x0}}, in, t[0], range));
}

Task858::Task858(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& a3 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a3, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a3, x2, x1}}, in, t[0], range));
}

Task859::Task859(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a3 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a3, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a3, x3, x0}}, in, t[0], range));
}

Task860::Task860(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a3 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a3, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a3, x3, x0}}, in, t[0], range));
}

Task861::Task861(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a1}}, in, t[0], range));
}

Task862::Task862(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a3 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a3, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a3, x3, x0}}, in, t[0], range));
}

Task863::Task863(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& a3 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a3, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a3, x2, x1}}, in, t[0], range));
}

Task864::Task864(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a3 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a3, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a3, x3, x0}}, in, t[0], range));
}

Task865::Task865(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a1}}, in, t[0], range));
}

Task866::Task866(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, x5, x0, x4))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x5, x0, x4}}, in, t[0], range));
}

Task867::Task867(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, x5, x0, x4))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x5, x0, x4}}, in, t[0], range));
}

Task868::Task868(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a1}}, in, t[0], range));
}

Task869::Task869(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a4 : *range[2])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, a4, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, a4, x0, x1}}, in, t[0], range));
}

Task870::Task870(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a1}}, in, t[0], range));
}

Task871::Task871(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a4 : *range[2])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, a4, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, a4, x0, x1}}, in, t[0], range));
}

Task872::Task872(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a1}}, in, t[0], range));
}

Task873::Task873(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& a4 : *range[2])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, a4, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, a4, x1, x0}}, in, t[0], range));
}

Task874::Task874(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a1}}, in, t[0], range));
}

Task875::Task875(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a4 : *range[2])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, a4, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, a4, x0, x1}}, in, t[0], range));
}

Task876::Task876(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a1}}, in, t[0], range));
}

Task877::Task877(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, a3, c2, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, a3, c2, a1}}, in, t[0], range));
}

Task878::Task878(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a1}}, in, t[0], range));
}

Task879::Task879(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x5 : *range[1])
          if (t[0]->is_local(x5, a3, c2, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x5, a3, c2, a1}}, in, t[0], range));
}

Task880::Task880(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task881::Task881(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, x0, x1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, a2}}, in, t[0], range));
}

Task882::Task882(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task883::Task883(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, x0, x1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, a2}}, in, t[0], range));
}

Task884::Task884(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, a1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, a1, a2}}, in, t[0], range));
}

Task885::Task885(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, a1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, a1, a2}}, in, t[0], range));
}

Task886::Task886(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, a1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, a1, a2}}, in, t[0], range));
}

Task887::Task887(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, a1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, a1, a2}}, in, t[0], range));
}

Task888::Task888(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, a1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, a1, a2}}, in, t[0], range));
}

Task889::Task889(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, x0, x1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, a2}}, in, t[0], range));
}

Task890::Task890(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& c3 : *range[0])
            for (auto& a2 : *range[2])
              if (t[0]->is_local(a2, c3, x5, x0, x4, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{a2, c3, x5, x0, x4, x1}}, in, t[0], range));
}

Task891::Task891(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, x0, x1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, a2}}, in, t[0], range));
}

Task892::Task892(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& c3 : *range[0])
            for (auto& a2 : *range[2])
              if (t[0]->is_local(a2, c3, x5, x4, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{a2, c3, x5, x4, x1, x0}}, in, t[0], range));
}

Task893::Task893(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, x0, x1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, a2}}, in, t[0], range));
}

Task894::Task894(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& x6 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x7 : *range[1])
            for (auto& a2 : *range[2])
              if (t[0]->is_local(a2, x7, x0, x6, x5, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{a2, x7, x0, x6, x5, x1}}, in, t[0], range));
}

Task895::Task895(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& x6 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x7 : *range[1])
            for (auto& a2 : *range[2])
              if (t[0]->is_local(a2, x7, x0, x6, x5, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{a2, x7, x0, x6, x5, x1}}, in, t[0], range));
}

Task896::Task896(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, x0, x1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, a2}}, in, t[0], range));
}

Task897::Task897(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, x1, x2, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, x1, x2, x0}}, in, t[0], range));
}

Task898::Task898(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, x0, x1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, a2}}, in, t[0], range));
}

Task899::Task899(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, x5, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, x5, x0, x1}}, in, t[0], range));
}

#endif
