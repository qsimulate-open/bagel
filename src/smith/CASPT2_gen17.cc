//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_gen17.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/CASPT2_tasks17.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

Task800::Task800(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x3 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x3, c2}}, in, t[0], range));
}

Task801::Task801(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task802::Task802(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x6 : *range[1])
      for (auto& x7 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x7, x6, x5}}, in, t[0], range));
}

Task803::Task803(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, c1, x1, x0}}, in, t[0], range));
}

Task804::Task804(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task805::Task805(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range));
}

Task806::Task806(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x1, x0, c2}}, in, t[0], range));
}

Task807::Task807(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, c1, x1, x0}}, in, t[0], range));
}

Task808::Task808(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range));
}

Task809::Task809(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, c1, x5, x4}}, in, t[0], range));
}

Task810::Task810(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock());
  for (auto& x3 : *range[1])
    for (auto& a2 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a2, x3}}, in, t[0], range));
}

Task811::Task811(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range));
}

Task812::Task812(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, c1, x2}}, in, t[0], range, e));
}

Task813::Task813(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, c1, x2}}, in, t[0], range));
}

Task814::Task814(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, a2, c1}}, in, t[0], range));
}

Task815::Task815(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range));
}

Task816::Task816(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, c1, x1, x0}}, in, t[0], range));
}

Task817::Task817(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x1, x0, x5, x4, x3}}, in, t[0], range));
}

Task818::Task818(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, c1, x4, x5}}, in, t[0], range));
}

Task819::Task819(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task820::Task820(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task821::Task821(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, x3}}, in, t[0], range));
}

Task822::Task822(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c2 : *range[0])
    for (auto& a3 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a3, c2}}, in, t[0], range));
}

Task823::Task823(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, x3}}, in, t[0], range));
}

Task824::Task824(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c2 : *range[0])
    for (auto& a3 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a3, c2}}, in, t[0], range));
}

Task825::Task825(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task826::Task826(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, c2, a1, x3}}, in, t[0], range));
}

Task827::Task827(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x2 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, c3}}, in, t[0], range));
}

Task828::Task828(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task829::Task829(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, c1, x3}}, in, t[0], range));
}

Task830::Task830(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x2 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, c3}}, in, t[0], range));
}

Task831::Task831(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, c1, x3}}, in, t[0], range));
}

Task832::Task832(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x2 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, c3}}, in, t[0], range));
}

Task833::Task833(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task834::Task834(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, a2, c1}}, in, t[0], range));
}

Task835::Task835(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task836::Task836(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, c1, x1, x0}}, in, t[0], range));
}

Task837::Task837(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task838::Task838(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x3, x5, x4, x2, x1, x0}}, in, t[0], range));
}

Task839::Task839(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& x5 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x5, c1, x4}}, in, t[0], range));
}

Task840::Task840(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock());
  for (auto& x3 : *range[1])
    for (auto& a2 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a2, x3}}, in, t[0], range));
}

Task841::Task841(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task842::Task842(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x3, x2, x1}}, in, t[0], range));
}

Task843::Task843(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x0 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, c3}}, in, t[0], range));
}

Task844::Task844(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c3, a2, c1}}, in, t[0], range));
}

Task845::Task845(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x3, x2, x1}}, in, t[0], range));
}

Task846::Task846(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x0 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, c1}}, in, t[0], range));
}

Task847::Task847(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c3, a2, c1}}, in, t[0], range));
}

Task848::Task848(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x3, x2, x1}}, in, t[0], range));
}

Task849::Task849(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a2, c1, x1}}, in, t[0], range));
}

#endif
