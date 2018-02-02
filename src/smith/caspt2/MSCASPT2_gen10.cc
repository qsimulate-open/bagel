//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_gen10.cc
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

#include <src/smith/caspt2/MSCASPT2_tasks10.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

Task450::Task450(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[5], t[6]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[5]->is_local())
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, out, range));
}

Task451::Task451(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          if (t[1]->is_local(c1, a4, c3, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}

Task452::Task452(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          if (t[1]->is_local(c1, a2, c3, a4))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a2, c3, a4}}, in, t[0], range));
}

Task453::Task453(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.push_back(make_shared<Task_local>(in, t[0], range));
}

Task454::Task454(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c5 : *range[0])
    for (auto& c3 : *range[0])
      if (t[1]->is_local(c3, c5))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, c5}}, in, t[0], range));
}

Task455::Task455(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c5 : *range[0])
    for (auto& c3 : *range[0])
      if (t[0]->is_local(c3, c5))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, c5}}, in, t[0], range));
}

Task456::Task456(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.push_back(make_shared<Task_local>(in, t[0], range));
}

Task457::Task457(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c5 : *range[0])
    for (auto& c3 : *range[0])
      if (t[1]->is_local(c3, c5))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, c5}}, in, t[0], range));
}

Task458::Task458(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c5 : *range[0])
    for (auto& c3 : *range[0])
      if (t[0]->is_local(c3, c5))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, c5}}, in, t[0], range));
}

Task459::Task459(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.push_back(make_shared<Task_local>(in, t[0], range));
}

Task460::Task460(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a5 : *range[2])
      if (t[1]->is_local(a5, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a5, a4}}, in, t[0], range));
}

Task461::Task461(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a5 : *range[2])
    for (auto& a4 : *range[2])
      if (t[0]->is_local(a4, a5))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, a5}}, in, t[0], range));
}

Task462::Task462(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.push_back(make_shared<Task_local>(in, t[0], range));
}

Task463::Task463(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a5 : *range[2])
      if (t[1]->is_local(a5, a4))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a5, a4}}, in, t[0], range));
}

Task464::Task464(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a5 : *range[2])
    for (auto& a4 : *range[2])
      if (t[0]->is_local(a4, a5))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, a5}}, in, t[0], range));
}

Task465::Task465(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[5], t[6]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x1 : *range[1])
          if (t[5]->is_local(x3, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, out, range));
}

Task466::Task466(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x3}}, in, t[0], range));
}

Task467::Task467(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x3}}, in, t[0], range));
}

Task468::Task468(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[5], t[6]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x2 : *range[1])
              if (t[5]->is_local(x5, x0, x4, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x1, x3, x2}}, in, out, range));
}

Task469::Task469(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, x5, x4))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, x5, x4}}, in, t[0], range));
}

Task470::Task470(vector<shared_ptr<Tensor>> t, const bool reset) : reset_(reset) {
  den0ci_ =  t[0];
  den1ci_ =  t[1];
  den2ci_ =  t[2];
  den3ci_ =  t[3];
  den4ci_ =  t[4];
}

Task472::Task472(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[5]->is_local(x0, x3, x1, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x3, x1, x2}}, in, out, range));
}

Task473::Task473(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, x3, x2}}, in, t[0], range, e));
}

Task474::Task474(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              if (t[5]->is_local(x5, x4, x2, x3, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x2, x3, x1, x0}}, in, out, range));
}

Task475::Task475(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              if (t[0]->is_local(x2, x1, x0, x5, x4, x3))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range, e));
}

Task476::Task476(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[5]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, out, range));
}

Task477::Task477(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range, e));
}

Task478::Task478(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range, e));
}

Task479::Task479(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task480::Task480(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, x3, x2}}, in, t[0], range, e));
}

Task481::Task481(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task482::Task482(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, x3, x2}}, in, t[0], range, e));
}

Task483::Task483(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, x3, x2}}, in, t[0], range, e));
}

Task484::Task484(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, x3, x2}}, in, t[0], range, e));
}

Task485::Task485(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x1 : *range[1])
              if (t[5]->is_local(x5, x0, x4, x3, x2, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x3, x2, x1}}, in, out, range));
}

Task486::Task486(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              if (t[0]->is_local(x2, x1, x0, x5, x4, x3))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range, e));
}

Task487::Task487(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.push_back(make_shared<Task_local>(in, t[0], range));
}

Task488::Task488(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          if (t[1]->is_local(c1, a4, c3, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range, e));
}

Task489::Task489(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.push_back(make_shared<Task_local>(in, t[0], range));
}

Task490::Task490(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          if (t[1]->is_local(c1, a2, c3, a4))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a2, c3, a4}}, in, t[0], range, e));
}

Task491::Task491(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[5]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, out, range));
}

Task492::Task492(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range, e));
}

Task493::Task493(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range, e));
}

Task494::Task494(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x1 : *range[1])
          if (t[5]->is_local(x3, x0, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, out, range));
}

Task495::Task495(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x1, x0, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, x3, x2}}, in, t[0], range, e));
}

Task496::Task496(vector<shared_ptr<Tensor>> t, const bool reset) : reset_(reset) {
  den0ci_ =  t[0];
  den1ci_ =  t[1];
  den2ci_ =  t[2];
  den3ci_ =  t[3];
  den4ci_ =  t[4];
}

Task498::Task498(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[5]->is_local(x0, x3, x1, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x3, x1, x2}}, in, out, range));
}

Task499::Task499(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

#endif
