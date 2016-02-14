//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_gen20.cc
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

#include <src/smith/mrci/MRCI_tasks20.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

Task950::Task950(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, c1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, c1, a2}}, in, t[0], range));
}

Task951::Task951(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task952::Task952(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, x2, x1, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x2, x1, a1}}, in, t[0], range));
}

Task953::Task953(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, x2, x1, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x2, x1, a1}}, in, t[0], range));
}

Task954::Task954(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, x2, x1, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x2, x1, a1}}, in, t[0], range));
}

Task955::Task955(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task956::Task956(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task957::Task957(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task958::Task958(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& c3 : *range[0])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, c3, c1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c3, c1, a2}}, in, t[0], range));
}

Task959::Task959(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& c3 : *range[0])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, c3, c1, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c3, c1, a2}}, in, t[0], range));
}

Task960::Task960(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task961::Task961(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task962::Task962(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task963::Task963(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task964::Task964(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task965::Task965(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task966::Task966(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task967::Task967(vector<shared_ptr<Tensor>> t, const bool reset) : reset_(reset) {
  n_ =  t[0];
}

Task968::Task968(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task969::Task969(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task970::Task970(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task971::Task971(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x1, x0, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x1, x0, c1}}, in, t[0], range));
}

Task972::Task972(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task973::Task973(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, c3, a2, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c3, a2, c1}}, in, t[0], range));
}

Task974::Task974(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, c3, a2, c1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c3, a2, c1}}, in, t[0], range));
}

Task975::Task975(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task976::Task976(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, x1, a1, c2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, c2}}, in, t[0], range));
}

Task977::Task977(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, x1, a1, c2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, c2}}, in, t[0], range));
}

Task978::Task978(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task979::Task979(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
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

Task980::Task980(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
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

Task981::Task981(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task982::Task982(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, x2, x1, a1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x2, x1, a1}}, in, t[0], range));
}

Task983::Task983(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task984::Task984(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task985::Task985(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task986::Task986(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task987::Task987(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task988::Task988(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task989::Task989(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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
