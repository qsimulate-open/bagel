//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_gen8.cc
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

#include <src/smith/caspt2/CASPT2_tasks8.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

Task350::Task350(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c1 : *range[0])
    for (auto& c4 : *range[0])
      if (t[0]->is_local(c4, c1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c4, c1}}, in, t[0], range));
}

Task351::Task351(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c4 : *range[0])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, c4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, c4}}, in, t[0], range));
}

Task352::Task352(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x1}}, in, t[0], range));
}

Task353::Task353(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c4 : *range[0])
    for (auto& c1 : *range[0])
      if (t[0]->is_local(c1, c4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, c4}}, in, t[0], range));
}

Task354::Task354(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x1}}, in, t[0], range));
}

Task355::Task355(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c3 : *range[0])
    for (auto& c4 : *range[0])
      if (t[0]->is_local(c4, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c4, c3}}, in, t[0], range));
}

Task356::Task356(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c4 : *range[0])
    for (auto& c3 : *range[0])
      if (t[0]->is_local(c3, c4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, c4}}, in, t[0], range));
}

Task357::Task357(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x1}}, in, t[0], range));
}

Task358::Task358(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c4 : *range[0])
    for (auto& c3 : *range[0])
      if (t[0]->is_local(c3, c4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, c4}}, in, t[0], range));
}

Task359::Task359(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x1}}, in, t[0], range));
}

Task360::Task360(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a2 : *range[2])
      if (t[0]->is_local(a2, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a2, a4}}, in, t[0], range));
}

Task361::Task361(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a2 : *range[2])
      if (t[0]->is_local(a2, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a2, a4}}, in, t[0], range));
}

Task362::Task362(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x1}}, in, t[0], range));
}

Task363::Task363(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a2 : *range[2])
      if (t[0]->is_local(a2, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a2, a4}}, in, t[0], range));
}

Task364::Task364(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c3 : *range[0])
          if (t[0]->is_local(c3, a2, c1, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a2, c1, x1}}, in, t[0], range));
}

Task365::Task365(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, c1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c1}}, in, t[0], range));
}

Task366::Task366(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, c1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c1}}, in, t[0], range));
}

Task367::Task367(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x1, x0}}, in, t[0], range));
}

Task368::Task368(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c3, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c3, x1, x0}}, in, t[0], range));
}

Task369::Task369(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c3}}, in, t[0], range));
}

Task370::Task370(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c3}}, in, t[0], range));
}

Task371::Task371(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c1, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c1, x0, x1}}, in, t[0], range));
}

Task372::Task372(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x3 : *range[1])
          if (t[0]->is_local(x3, a2, c1, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, a2, c1, x2}}, in, t[0], range));
}

Task373::Task373(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock());
  for (auto& a4 : *range[2])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, a4}}, in, t[0], range));
}

Task374::Task374(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a4 : *range[2])
      if (t[0]->is_local(a4, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, x1}}, in, t[0], range));
}

Task375::Task375(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a4 : *range[2])
      if (t[0]->is_local(a4, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, x0}}, in, t[0], range));
}

Task376::Task376(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task377::Task377(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a1 : *range[2])
      if (t[0]->is_local(a1, x2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a1, x2}}, in, t[0], range));
}

Task378::Task378(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x2 : *range[1])
      if (t[0]->is_local(x2, a1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, a1}}, in, t[0], range));
}

Task379::Task379(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task380::Task380(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& c3 : *range[0])
      if (t[0]->is_local(c3, x2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, x2}}, in, t[0], range));
}

Task381::Task381(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x2 : *range[1])
      if (t[0]->is_local(x2, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, c3}}, in, t[0], range));
}

Task382::Task382(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, a1, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a1, x3, x2}}, in, t[0], range));
}

Task383::Task383(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x2 : *range[1])
      if (t[0]->is_local(x2, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, c3}}, in, t[0], range));
}

Task384::Task384(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, a1, x2, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a1, x2, x3}}, in, t[0], range));
}

Task385::Task385(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x2 : *range[1])
      if (t[0]->is_local(x2, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, c3}}, in, t[0], range));
}

Task386::Task386(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c1, x2, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c1, x2, x3}}, in, t[0], range));
}

Task387::Task387(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x2 : *range[1])
      if (t[0]->is_local(x2, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, c3}}, in, t[0], range));
}

Task388::Task388(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, c1, x2, x3))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, c1, x2, x3}}, in, t[0], range));
}

Task389::Task389(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x2 : *range[1])
      if (t[0]->is_local(x2, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, c3}}, in, t[0], range));
}

Task390::Task390(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& a2 : *range[2])
          if (t[0]->is_local(a2, a1, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a2, a1, x3, x2}}, in, t[0], range));
}

Task391::Task391(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      if (t[0]->is_local(a1, a3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a1, a3}}, in, t[0], range));
}

Task392::Task392(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      if (t[0]->is_local(a1, a3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a1, a3}}, in, t[0], range));
}

Task393::Task393(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, a1, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a1, x3, x2}}, in, t[0], range));
}

Task394::Task394(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      if (t[0]->is_local(a1, a3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a1, a3}}, in, t[0], range));
}

Task395::Task395(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, a1, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a1, x3, x2}}, in, t[0], range));
}

Task396::Task396(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& c3 : *range[0])
      if (t[0]->is_local(c3, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, a4}}, in, t[0], range));
}

Task397::Task397(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a4 : *range[2])
      if (t[0]->is_local(a4, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, c3}}, in, t[0], range));
}

Task398::Task398(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      if (t[0]->is_local(c2, a1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c2, a1}}, in, t[0], range));
}

Task399::Task399(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a4 : *range[2])
      if (t[0]->is_local(a4, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, c3}}, in, t[0], range));
}

#endif
