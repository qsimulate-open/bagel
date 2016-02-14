//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_gen7.cc
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

#include <src/smith/mrci/MRCI_tasks7.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

Task300::Task300(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task301::Task301(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x5, x0, x4))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x5, x0, x4}}, in, t[0], range));
}

Task302::Task302(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x5, x0, x4))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x5, x0, x4}}, in, t[0], range));
}

Task303::Task303(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task304::Task304(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task305::Task305(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task306::Task306(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task307::Task307(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c4 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c4, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c4, x0, x1}}, in, t[0], range));
}

Task308::Task308(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& c4 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, a2, c4, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, a2, c4, x2}}, in, t[0], range));
}

Task309::Task309(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task310::Task310(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, c3, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, x1, x0}}, in, t[0], range));
}

Task311::Task311(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, c3, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, x1, x0}}, in, t[0], range));
}

Task312::Task312(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task313::Task313(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c4 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c4, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c4, x1, x0}}, in, t[0], range));
}

Task314::Task314(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c4 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c4, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c4, x1, x0}}, in, t[0], range));
}

Task315::Task315(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task316::Task316(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, c3, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, x1, x0}}, in, t[0], range));
}

Task317::Task317(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, c3, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, x1, x0}}, in, t[0], range));
}

Task318::Task318(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task319::Task319(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, c1, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c1, x1, x0}}, in, t[0], range));
}

Task320::Task320(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, c1, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c1, x1, x0}}, in, t[0], range));
}

Task321::Task321(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task322::Task322(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& a4 : *range[2])
          if (t[0]->is_local(a4, c1, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c1, x0, x1}}, in, t[0], range));
}

Task323::Task323(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, a4, c1, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, a4, c1, x2}}, in, t[0], range));
}

Task324::Task324(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task325::Task325(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task326::Task326(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task327::Task327(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task328::Task328(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task329::Task329(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task330::Task330(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task331::Task331(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task332::Task332(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task333::Task333(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task334::Task334(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task335::Task335(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task336::Task336(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task337::Task337(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task338::Task338(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task339::Task339(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task340::Task340(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task341::Task341(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task342::Task342(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task343::Task343(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task344::Task344(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
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

Task345::Task345(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task346::Task346(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, x1, x2, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, x2, x0}}, in, t[0], range));
}

Task347::Task347(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task348::Task348(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c2 : *range[0])
      for (auto& a1 : *range[2])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, a1, c2, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a1, c2, x3}}, in, t[0], range));
}

Task349::Task349(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c2 : *range[0])
      for (auto& a1 : *range[2])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, a1, c2, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a1, c2, x3}}, in, t[0], range));
}

#endif
