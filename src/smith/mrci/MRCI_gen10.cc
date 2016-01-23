//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_gen10.cc
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

#include <src/smith/mrci/MRCI_tasks10.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

Task450::Task450(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, x1, x0, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, x0, a1}}, in, t[0], range));
}

Task451::Task451(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task452::Task452(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x1, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x1, x0, a2}}, in, t[0], range));
}

Task453::Task453(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task454::Task454(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x1, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x1, x0, a2}}, in, t[0], range));
}

Task455::Task455(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, a2, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, c1, x3}}, in, t[0], range));
}

Task456::Task456(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, a2, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, c1, x3}}, in, t[0], range));
}

Task457::Task457(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, a2, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, c1, x3}}, in, t[0], range));
}

Task458::Task458(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, a2, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, c1, x3}}, in, t[0], range));
}

Task459::Task459(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, a2, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, c1, x3}}, in, t[0], range));
}

Task460::Task460(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, a2, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, c1, x3}}, in, t[0], range));
}

Task461::Task461(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, a2, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, c1, x3}}, in, t[0], range));
}

Task462::Task462(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, a2, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, c1, x3}}, in, t[0], range));
}

Task463::Task463(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, a2, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, c1, x3}}, in, t[0], range));
}

Task464::Task464(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, a2, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, c1, x3}}, in, t[0], range));
}

Task465::Task465(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x1, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x1, x0, a2}}, in, t[0], range));
}

Task466::Task466(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task467::Task467(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task468::Task468(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task469::Task469(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task470::Task470(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task471::Task471(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task472::Task472(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task473::Task473(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task474::Task474(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task475::Task475(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task476::Task476(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task477::Task477(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, c3, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, x3, x2}}, in, t[0], range));
}

Task478::Task478(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task479::Task479(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, c3, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, x3, x2}}, in, t[0], range));
}

Task480::Task480(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task481::Task481(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task482::Task482(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task483::Task483(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& c3 : *range[0])
    for (auto& c1 : *range[0])
      for (auto& x2 : *range[1])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, x2, c1, c3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, x2, c1, c3}}, in, t[0], range));
}

Task484::Task484(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task485::Task485(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& c3 : *range[0])
    for (auto& c1 : *range[0])
      for (auto& x2 : *range[1])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, x2, c1, c3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, x2, c1, c3}}, in, t[0], range));
}

Task486::Task486(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task487::Task487(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x3, a2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, a2, x2}}, in, t[0], range));
}

Task488::Task488(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x1, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x1, x0, a2}}, in, t[0], range));
}

Task489::Task489(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, x2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, x2, x1, x0}}, in, t[0], range));
}

Task490::Task490(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x1, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x1, x0, a2}}, in, t[0], range));
}

Task491::Task491(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, a2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, a2}}, in, t[0], range));
}

Task492::Task492(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, a2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, a2}}, in, t[0], range));
}

Task493::Task493(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, a2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, a2}}, in, t[0], range));
}

Task494::Task494(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, a2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, a2}}, in, t[0], range));
}

Task495::Task495(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, a2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, a2}}, in, t[0], range));
}

Task496::Task496(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, a2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, a2}}, in, t[0], range));
}

Task497::Task497(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x1, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x1, x0, a2}}, in, t[0], range));
}

Task498::Task498(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& c3 : *range[0])
            for (auto& c1 : *range[0])
              if (t[0]->is_local(c1, c3, x3, x2, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{c1, c3, x3, x2, x1, x0}}, in, t[0], range));
}

Task499::Task499(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x1, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x1, x0, a2}}, in, t[0], range));
}

#endif
