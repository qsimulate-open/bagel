//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks6.h
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

#include <src/smith/relmrci/RelMRCI_tasks6.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

Task250::Task250(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a2, x0, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a2, x0, x3}}, in, t[0], range));
}

Task251::Task251(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a2, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a2, x2, x1}}, in, t[0], range));
}

Task252::Task252(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a2, x0, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a2, x0, x3}}, in, t[0], range));
}

Task253::Task253(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a2, x0, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a2, x0, x3}}, in, t[0], range));
}

Task254::Task254(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task255::Task255(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c4 : *range[0])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, c4, x0, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, c4, x0, x3}}, in, t[0], range));
}

Task256::Task256(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& c4 : *range[0])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, c4, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, c4, x2, x1}}, in, t[0], range));
}

Task257::Task257(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c4 : *range[0])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, c4, x0, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, c4, x0, x3}}, in, t[0], range));
}

Task258::Task258(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c4 : *range[0])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, c4, x0, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, c4, x0, x3}}, in, t[0], range));
}

Task259::Task259(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task260::Task260(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a2, x0, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a2, x0, x3}}, in, t[0], range));
}

Task261::Task261(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a2, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a2, x2, x1}}, in, t[0], range));
}

Task262::Task262(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, a2, x0, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, a2, x0, x3}}, in, t[0], range));
}

Task263::Task263(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task264::Task264(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x5, x4, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x5, x4, x0}}, in, t[0], range));
}

Task265::Task265(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x5, x4, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x5, x4, x0}}, in, t[0], range));
}

Task266::Task266(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task267::Task267(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, x5, x4, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, x5, x4, x0}}, in, t[0], range));
}

Task268::Task268(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, x5, x4, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, x5, x4, x0}}, in, t[0], range));
}

Task269::Task269(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task270::Task270(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& c4 : *range[0])
          if (t[0]->is_local(c4, a2, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c4, a2, x0, x1}}, in, t[0], range));
}

Task271::Task271(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task272::Task272(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a4, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a4, x0, x1}}, in, t[0], range));
}

Task273::Task273(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task274::Task274(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& c4 : *range[0])
          if (t[0]->is_local(c4, a2, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c4, a2, x0, x1}}, in, t[0], range));
}

Task275::Task275(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task276::Task276(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a4, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a4, x0, x1}}, in, t[0], range));
}

Task277::Task277(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task278::Task278(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, a4, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, x0, x1}}, in, t[0], range));
}

Task279::Task279(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task280::Task280(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, a4, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, x0, x1}}, in, t[0], range));
}

Task281::Task281(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task282::Task282(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, x2, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, x2, x0, x1}}, in, t[0], range));
}

Task283::Task283(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task284::Task284(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c4 : *range[0])
      for (auto& c1 : *range[0])
        for (auto& a5 : *range[2])
          if (t[0]->is_local(a5, c1, c4, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a5, c1, c4, x0}}, in, t[0], range));
}

Task285::Task285(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task286::Task286(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, a4, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, x3, x0}}, in, t[0], range));
}

Task287::Task287(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task288::Task288(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, a4, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, x3, x0}}, in, t[0], range));
}

Task289::Task289(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task290::Task290(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a4, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a4, x3, x0}}, in, t[0], range));
}

Task291::Task291(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task292::Task292(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a4, x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a4, x3, x0}}, in, t[0], range));
}

Task293::Task293(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task294::Task294(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x3}}, in, t[0], range));
}

Task295::Task295(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c3, x0, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x0, a2}}, in, t[0], range));
}

Task296::Task296(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x5 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x5))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x5}}, in, t[0], range));
}

Task297::Task297(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task298::Task298(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task299::Task299(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

#endif
