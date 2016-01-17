//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_gen2.cc
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

#include <src/smith/relmrci/RelMRCI_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

Task50::Task50(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x3, x0, x2, x1}}, in, t[0], range));
}

Task51::Task51(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x3, x2, x1}}, in, t[0], range));
}

Task52::Task52(array<shared_ptr<Tensor>,2> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}

Task53::Task53(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& x6 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x6, x3, x5, x4, x0, x2, x1}}, in, t[0], range));
}

Task54::Task54(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x4 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x4, x5, x3, x0, x2, x1}}, in, t[0], range));
}

Task55::Task55(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& x6 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x6, x5, x0, x4, x3, x2, x1}}, in, t[0], range));
}

Task56::Task56(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& x6 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x6, x5, x4, x3, x0, x2, x1}}, in, t[0], range));
}

Task57::Task57(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x6 : *range[1])
              for (auto& x0 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x0, x6, x5, x4, x3, x2, x1}}, in, t[0], range));
}

Task58::Task58(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x6 : *range[1])
              for (auto& x4 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x4, x6, x5, x3, x0, x2, x1}}, in, t[0], range));
}

Task59::Task59(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x6 : *range[1])
              for (auto& x3 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x3, x6, x5, x4, x0, x2, x1}}, in, t[0], range));
}

Task60::Task60(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x3, x4, x2, x1}}, in, t[0], range));
}

Task61::Task61(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x3, x4, x0, x2, x1}}, in, t[0], range));
}

Task62::Task62(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x6 : *range[1])
              for (auto& x0 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x0, x6, x3, x5, x4, x2, x1}}, in, t[0], range));
}

Task63::Task63(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x6 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x7 : *range[1])
              for (auto& x3 : *range[1])
                for (auto& x4 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x0, x6, x5, x2, x1, x4, x3}}, in, t[0], range));
}

Task64::Task64(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x7 : *range[1])
        for (auto& x8 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x9 : *range[1])
              for (auto& x3 : *range[1])
                for (auto& x4 : *range[1])
                  for (auto& x5 : *range[1])
                    for (auto& x6 : *range[1])
                      subtasks_.push_back(make_shared<Task_local>(array<const Index,10>{{x9, x0, x8, x7, x2, x1, x6, x5, x4, x3}}, in, t[0], range));
}

Task65::Task65(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x1, x3, x2, x0}}, in, t[0], range));
}

Task66::Task66(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x3, x2, x0}}, in, t[0], range));
}

Task67::Task67(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x3, x1, x2}}, in, t[0], range));
}

Task68::Task68(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x1, x4, x3, x2, x0}}, in, t[0], range));
}

Task69::Task69(array<shared_ptr<Tensor>,2> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x1, x2, x0}}, in, t[0], range));
}

Task70::Task70(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x1, x2}}, in, t[0], range));
}

Task71::Task71(array<shared_ptr<Tensor>,2> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x1, x3, x2}}, in, t[0], range));
}

Task72::Task72(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}

Task73::Task73(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x4 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x3, x2, x1}}, in, t[0], range));
}

Task74::Task74(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x3, x1, x2, x0}}, in, t[0], range));
}

Task75::Task75(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x6 : *range[1])
              for (auto& x0 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x0, x6, x5, x4, x1, x3, x2}}, in, t[0], range));
}

Task76::Task76(array<shared_ptr<Tensor>,2> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x2, x3, x1}}, in, t[0], range));
}

Task77::Task77(array<shared_ptr<Tensor>,6> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,5> in = {{t[1], t[2], t[3], t[4], t[5]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x5, x1, x4, x3, x2}}, in, t[0], range));
}

Task78::Task78(array<shared_ptr<Tensor>,7> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,6> in = {{t[1], t[2], t[3], t[4], t[5], t[6]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x6 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x7 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& x4 : *range[1])
                for (auto& x5 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x0, x7, x1, x6, x5, x4, x3, x2}}, in, t[0], range));
}

Task79::Task79(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task80::Task80(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task81::Task81(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x1, x3, x2}}, in, t[0], range));
}

Task82::Task82(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x6 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x7 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& x4 : *range[1])
                for (auto& x5 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x0, x6, x1, x5, x4, x3, x2}}, in, t[0], range));
}

#endif
