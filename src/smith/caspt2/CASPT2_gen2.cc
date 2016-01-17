//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_gen2.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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

#include <src/smith/caspt2/CASPT2_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

Task50::Task50(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(array<const Index,7>{{ci0, x5, x4, x1, x3, x2, x0}}, in, t[0], range));
}

Task51::Task51(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(array<const Index,5>{{ci0, x1, x3, x2, x0}}, in, t[0], range));
}

Task52::Task52(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            for (auto& x2 : *range[1])
              for (auto& x3 : *range[1])
                subtasks_.push_back(make_shared<Task_local>(array<const Index,7>{{ci0, x5, x0, x1, x4, x3, x2}}, in, t[0], range));
}

Task53::Task53(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(array<const Index,5>{{ci0, x3, x0, x1, x2}}, in, t[0], range));
}

Task54::Task54(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            for (auto& x2 : *range[1])
              for (auto& x3 : *range[1])
                subtasks_.push_back(make_shared<Task_local>(array<const Index,7>{{ci0, x5, x4, x1, x0, x3, x2}}, in, t[0], range));
}

Task55::Task55(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(array<const Index,5>{{ci0, x3, x2, x1, x0}}, in, t[0], range));
}

Task56::Task56(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(array<const Index,7>{{ci0, x5, x0, x4, x3, x1, x2}}, in, t[0], range));
}

Task57::Task57(array<shared_ptr<Tensor>,2> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}

Task58::Task58(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(array<const Index,7>{{ci0, x5, x2, x4, x3, x1, x0}}, in, t[0], range));
}

Task59::Task59(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(array<const Index,7>{{ci0, x5, x0, x3, x4, x2, x1}}, in, t[0], range));
}

Task60::Task60(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(array<const Index,7>{{ci0, x5, x4, x3, x0, x2, x1}}, in, t[0], range));
}

Task61::Task61(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x6 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x7 : *range[1])
              for (auto& ci0 : *range[3])
                for (auto& x3 : *range[1])
                  for (auto& x4 : *range[1])
                    subtasks_.push_back(make_shared<Task_local>(array<const Index,9>{{ci0, x7, x0, x6, x5, x2, x1, x4, x3}}, in, t[0], range));
}

Task62::Task62(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(array<const Index,7>{{ci0, x5, x0, x4, x3, x2, x1}}, in, t[0], range));
}

Task63::Task63(array<shared_ptr<Tensor>,2> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range));
}

Task64::Task64(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& ci0 : *range[3])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        subtasks_.push_back(make_shared<Task_local>(array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}

Task65::Task65(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& ci0 : *range[3])
        for (auto& x1 : *range[1])
          for (auto& x2 : *range[1])
            subtasks_.push_back(make_shared<Task_local>(array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range));
}

Task66::Task66(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            for (auto& x2 : *range[1])
              for (auto& x3 : *range[1])
                subtasks_.push_back(make_shared<Task_local>(array<const Index,7>{{ci0, x5, x0, x4, x1, x3, x2}}, in, t[0], range));
}

Task67::Task67(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x2 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(array<const Index,7>{{ci0, x2, x5, x4, x3, x1, x0}}, in, t[0], range));
}

Task68::Task68(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(array<const Index,7>{{ci0, x5, x4, x0, x3, x2, x1}}, in, t[0], range));
}

#endif
