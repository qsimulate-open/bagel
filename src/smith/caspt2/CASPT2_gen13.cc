//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_gen13.cc
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

#include <src/smith/caspt2/CASPT2_tasks13.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

Task600::Task600(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task601::Task601(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock());
  for (auto& a4 : *range[2])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, a4}}, in, t[0], range));
}

Task602::Task602(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task603::Task603(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x1}}, in, t[0], range, e));
}

Task604::Task604(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task605::Task605(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task606::Task606(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task607::Task607(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task608::Task608(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x2 : *range[1])
          if (t[5]->is_local(x3, x1, x0, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x1, x0, x2}}, in, out, range));
}

Task609::Task609(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, x1, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, x3, x2}}, in, t[0], range));
}

Task610::Task610(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, c3, a2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c3, a2, x1}}, in, t[0], range));
}

Task611::Task611(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x0, x1, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, x3, x2}}, in, t[0], range));
}

Task612::Task612(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x0 : *range[1])
              if (t[5]->is_local(x5, x4, x1, x3, x2, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x1, x3, x2, x0}}, in, out, range));
}

Task613::Task613(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              if (t[0]->is_local(x1, x0, x2, x5, x4, x3))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, x2, x5, x4, x3}}, in, t[0], range));
}

Task614::Task614(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, c2, x0, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, c2, x0, x2}}, in, t[0], range));
}

Task615::Task615(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x0 : *range[1])
          if (t[5]->is_local(x1, x3, x2, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x3, x2, x0}}, in, out, range));
}

Task616::Task616(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, x2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task617::Task617(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a1 : *range[2])
          if (t[0]->is_local(a1, c2, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, c2, x3, x2}}, in, t[0], range));
}

Task618::Task618(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, x2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task619::Task619(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[5], t[6]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x2 : *range[1])
              if (t[5]->is_local(x5, x0, x1, x4))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x1, x4, x3, x2}}, in, out, range));
}

Task620::Task620(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          if (t[0]->is_local(x5, x4, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x5, x4, x1, x0}}, in, t[0], range));
}

Task621::Task621(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[5]->is_local(x3, x0, x1, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x1, x2}}, in, out, range));
}

Task622::Task622(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task623::Task623(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, a1, x0, c3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a1, x0, c3}}, in, t[0], range));
}

Task624::Task624(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task625::Task625(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, c2, x0, a3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, c2, x0, a3}}, in, t[0], range));
}

Task626::Task626(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task627::Task627(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& c2 : *range[0])
      for (auto& a1 : *range[2])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, a1, c2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, a1, c2, x2}}, in, t[0], range, e));
}

Task628::Task628(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& c2 : *range[0])
      for (auto& a1 : *range[2])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, a1, c2, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, a1, c2, x2}}, in, t[0], range));
}

Task629::Task629(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task630::Task630(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, a1, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a1, x0, x1}}, in, t[0], range));
}

Task631::Task631(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task632::Task632(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task633::Task633(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[5], t[6]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x2 : *range[1])
              if (t[5]->is_local(x5, x4, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x1, x0, x3, x2}}, in, out, range));
}

Task634::Task634(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          if (t[0]->is_local(x5, x4, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x5, x4, x1, x0}}, in, t[0], range));
}

Task635::Task635(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          if (t[0]->is_local(x5, x4, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x5, x4, x1, x0}}, in, t[0], range));
}

Task636::Task636(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x5 : *range[1])
          if (t[0]->is_local(x5, a2, c1, x4))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x5, a2, c1, x4}}, in, t[0], range));
}

Task637::Task637(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          if (t[5]->is_local(x3, x2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, out, range));
}

Task638::Task638(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task639::Task639(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, a1, x0, c3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a1, x0, c3}}, in, t[0], range));
}

Task640::Task640(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task641::Task641(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, c2, x0, a3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, c2, x0, a3}}, in, t[0], range));
}

Task642::Task642(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, c2, x0, a3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, c2, x0, a3}}, in, t[0], range));
}

Task643::Task643(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task644::Task644(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, c2, a1, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, c2, a1, x2}}, in, t[0], range, e));
}

Task645::Task645(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, c2, a1, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, c2, a1, x2}}, in, t[0], range));
}

Task646::Task646(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task647::Task647(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, a2, c3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, a2, c3}}, in, t[0], range));
}

Task648::Task648(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task649::Task649(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, c1, a3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, c1, a3}}, in, t[0], range));
}

#endif
