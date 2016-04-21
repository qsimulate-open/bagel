//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_gen14.cc
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

#include <src/smith/caspt2/MSCASPT2_tasks14.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

Task650::Task650(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task651::Task651(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x3}}, in, t[0], range));
}

Task652::Task652(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x3}}, in, t[0], range));
}

Task653::Task653(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x3}}, in, t[0], range));
}

Task654::Task654(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x3}}, in, t[0], range));
}

Task655::Task655(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task656::Task656(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task657::Task657(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task658::Task658(vector<shared_ptr<Tensor>> t, const bool reset) : reset_(reset) {
  deci_ =  t[0];
}

Task659::Task659(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task660::Task660(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task661::Task661(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task662::Task662(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task663::Task663(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task664::Task664(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range));
}

Task665::Task665(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task666::Task666(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range));
}

Task667::Task667(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range));
}

Task668::Task668(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task669::Task669(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task670::Task670(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task671::Task671(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task672::Task672(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task673::Task673(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task674::Task674(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task675::Task675(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task676::Task676(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task677::Task677(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task678::Task678(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task679::Task679(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task680::Task680(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task681::Task681(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task682::Task682(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task683::Task683(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task684::Task684(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task685::Task685(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task686::Task686(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task687::Task687(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task688::Task688(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task689::Task689(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task690::Task690(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task691::Task691(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task692::Task692(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task693::Task693(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range));
}

Task694::Task694(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range));
}

Task695::Task695(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task696::Task696(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range));
}

Task697::Task697(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task698::Task698(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task699::Task699(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

#endif
