//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_tasks3.h
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

#include <src/smith/relmrci/RelMRCI_tasks3.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

Task100::Task100(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, c1, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, c1, x0, x1}}, in, t[0], range));
}

Task101::Task101(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x2 : *range[1])
          for (auto& c2 : *range[0])
            for (auto& x3 : *range[1])
              if (t[0]->is_local(x3, c2, x2, c1, x5, x4))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x3, c2, x2, c1, x5, x4}}, in, t[0], range));
}

Task102::Task102(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, c1, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, c1, x0, x1}}, in, t[0], range));
}

Task103::Task103(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x6 : *range[1])
          for (auto& x7 : *range[1])
            for (auto& c2 : *range[0])
              if (t[0]->is_local(c2, x7, x6, x0, x5, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{c2, x7, x6, x0, x5, x1}}, in, t[0], range));
}

Task104::Task104(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x6 : *range[1])
          for (auto& x7 : *range[1])
            for (auto& c2 : *range[0])
              if (t[0]->is_local(c2, x7, x6, x0, x5, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{c2, x7, x6, x0, x5, x1}}, in, t[0], range));
}

Task105::Task105(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, c1, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, c1, x0, x1}}, in, t[0], range));
}

Task106::Task106(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, x1, x0, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, x1, x0, x2}}, in, t[0], range));
}

Task107::Task107(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, c1, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, c1, x0, x1}}, in, t[0], range));
}

Task108::Task108(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, x1, x5, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, x1, x5, x0}}, in, t[0], range));
}

Task109::Task109(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, x1, x5, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, x1, x5, x0}}, in, t[0], range));
}

Task110::Task110(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, c1, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, c1, x0, x1}}, in, t[0], range));
}

Task111::Task111(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x2 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& c2 : *range[0])
              if (t[0]->is_local(c2, x3, x2, c1, x5, x4))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{c2, x3, x2, c1, x5, x4}}, in, t[0], range));
}

Task112::Task112(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task113::Task113(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task114::Task114(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, c1, x5, x4))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, c1, x5, x4}}, in, t[0], range));
}

Task115::Task115(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, c1, x5, x4))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, c1, x5, x4}}, in, t[0], range));
}

Task116::Task116(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, c1, x5, x4))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, c1, x5, x4}}, in, t[0], range));
}

Task117::Task117(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task118::Task118(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x5, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x5, x4, x3}}, in, t[0], range));
}

Task119::Task119(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x5, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x5, x4, x3}}, in, t[0], range));
}

Task120::Task120(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x5, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x5, x4, x3}}, in, t[0], range));
}

Task121::Task121(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x5, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x5, x4, x3}}, in, t[0], range));
}

Task122::Task122(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task123::Task123(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, x3}}, in, t[0], range));
}

Task124::Task124(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, x3}}, in, t[0], range));
}

Task125::Task125(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, x3}}, in, t[0], range));
}

Task126::Task126(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, x3}}, in, t[0], range));
}

Task127::Task127(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, x3}}, in, t[0], range));
}

Task128::Task128(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, x3}}, in, t[0], range));
}

Task129::Task129(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task130::Task130(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x6 : *range[1])
    for (auto& x7 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              if (t[0]->is_local(x5, x4, x3, c1, x7, x6))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x3, c1, x7, x6}}, in, t[0], range));
}

Task131::Task131(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task132::Task132(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x6 : *range[1])
    for (auto& x7 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              if (t[0]->is_local(x5, x4, x3, c1, x7, x6))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x3, c1, x7, x6}}, in, t[0], range));
}

Task133::Task133(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task134::Task134(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x6 : *range[1])
      for (auto& x7 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& c1 : *range[0])
              if (t[0]->is_local(c1, x4, x3, x7, x6, x5))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{c1, x4, x3, x7, x6, x5}}, in, t[0], range));
}

Task135::Task135(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, c2, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c2, x4, x3}}, in, t[0], range));
}

Task136::Task136(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x6 : *range[1])
      for (auto& x7 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& c1 : *range[0])
              if (t[0]->is_local(c1, x4, x3, x7, x6, x5))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{c1, x4, x3, x7, x6, x5}}, in, t[0], range));
}

Task137::Task137(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task138::Task138(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x6 : *range[1])
      for (auto& x7 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& c1 : *range[0])
              if (t[0]->is_local(c1, x4, x3, x7, x6, x5))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{c1, x4, x3, x7, x6, x5}}, in, t[0], range));
}

Task139::Task139(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task140::Task140(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x6 : *range[1])
      for (auto& x7 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& c1 : *range[0])
            for (auto& x4 : *range[1])
              if (t[0]->is_local(x4, c1, x3, x7, x6, x5))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x4, c1, x3, x7, x6, x5}}, in, t[0], range));
}

Task141::Task141(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task142::Task142(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          if (t[0]->is_local(x4, x3, c1, x5))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x4, x3, c1, x5}}, in, t[0], range));
}

Task143::Task143(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x4, x3}}, in, t[0], range));
}

Task144::Task144(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          if (t[0]->is_local(x4, x3, c1, x5))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x4, x3, c1, x5}}, in, t[0], range));
}

Task145::Task145(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          if (t[0]->is_local(a3, c2, x4, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x4, x3}}, in, t[0], range));
}

Task146::Task146(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          if (t[0]->is_local(x4, x3, c1, x5))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x4, x3, c1, x5}}, in, t[0], range));
}

Task147::Task147(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task148::Task148(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          if (t[0]->is_local(x4, x3, c1, x5))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x4, x3, c1, x5}}, in, t[0], range));
}

Task149::Task149(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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
