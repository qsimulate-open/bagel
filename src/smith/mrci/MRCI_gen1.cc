//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_gen1.cc
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

#include <src/smith/mrci/MRCI_tasks1.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

Task0::Task0(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x3, x1, x2}}, in, t[0], range));
}

Task1::Task1(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x0, x3, x1, x2}}, in, t[0], range));
}

Task2::Task2(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x3, x0, x2}}, in, t[0], range));
}

Task3::Task3(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x5, x1, x4, x3, x2}}, in, t[0], range));
}

Task4::Task4(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x5, x2, x4, x1, x3}}, in, t[0], range));
}

Task5::Task5(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x5, x3, x4, x1, x2}}, in, t[0], range));
}

Task6::Task6(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x0 : *range[1])
              for (auto& x6 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x6, x0, x5, x1, x4, x3, x2}}, in, t[0], range));
}

Task7::Task7(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x0 : *range[1])
              for (auto& x6 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x6, x0, x5, x4, x3, x1, x2}}, in, t[0], range));
}

Task8::Task8(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x1, x3, x0, x2}}, in, t[0], range));
}

Task9::Task9(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x5, x0, x4, x3, x2}}, in, t[0], range));
}

Task10::Task10(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x5, x4, x3, x0, x2}}, in, t[0], range));
}

Task11::Task11(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x2, x0, x4, x1, x3}}, in, t[0], range));
}

Task12::Task12(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x2 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x5, x3, x4, x1, x0}}, in, t[0], range));
}

Task13::Task13(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x2, x3, x1, x0}}, in, t[0], range));
}

Task14::Task14(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x1, x0}}, in, t[0], range));
}

Task15::Task15(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x3, x2, x4, x1, x0}}, in, t[0], range));
}

Task16::Task16(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x6 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& x7 : *range[1])
                for (auto& x2 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x2, x7, x5, x6, x4, x3, x1, x0}}, in, t[0], range));
}

Task17::Task17(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& x6 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& x7 : *range[1])
                for (auto& x2 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x2, x7, x3, x6, x5, x4, x1, x0}}, in, t[0], range));
}

Task18::Task18(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x2 : *range[1])
              for (auto& x6 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x6, x2, x5, x4, x3, x1, x0}}, in, t[0], range));
}

Task19::Task19(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& x6 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x6, x3, x5, x2, x4, x1, x0}}, in, t[0], range));
}

Task20::Task20(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x4 : *range[1])
              for (auto& x6 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x6, x4, x5, x2, x3, x1, x0}}, in, t[0], range));
}

Task21::Task21(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x2 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x5, x4, x3, x1, x0}}, in, t[0], range));
}

Task22::Task22(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x3, x5, x2, x4, x1, x0}}, in, t[0], range));
}

Task23::Task23(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x4 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x4, x5, x2, x3, x1, x0}}, in, t[0], range));
}

Task24::Task24(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x6 : *range[1])
            for (auto& x2 : *range[1])
              for (auto& x5 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x5, x2, x6, x4, x3, x1, x0}}, in, t[0], range));
}

Task25::Task25(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& x6 : *range[1])
            for (auto& x2 : *range[1])
              for (auto& x3 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x3, x2, x6, x5, x4, x1, x0}}, in, t[0], range));
}

Task26::Task26(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& x6 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x6, x5, x4, x2, x3, x1, x0}}, in, t[0], range));
}

Task27::Task27(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x6 : *range[1])
              for (auto& x3 : *range[1])
                for (auto& x7 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x3, x6, x5, x2, x4, x1, x0}}, in, t[0], range));
}

Task28::Task28(array<shared_ptr<Tensor>,6> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,5> in = {{t[1], t[2], t[3], t[4], t[5]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x6 : *range[1])
            for (auto& x7 : *range[1])
              for (auto& x3 : *range[1])
                for (auto& x4 : *range[1])
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x6, x2, x5, x1, x0, x4, x3}}, in, t[0], range));
}

Task29::Task29(array<shared_ptr<Tensor>,6> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,5> in = {{t[1], t[2], t[3], t[4], t[5]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x7 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x8 : *range[1])
            for (auto& x9 : *range[1])
              for (auto& x3 : *range[1])
                for (auto& x4 : *range[1])
                  for (auto& x5 : *range[1])
                    for (auto& x6 : *range[1])
                      subtasks_.push_back(make_shared<Task_local>(array<const Index,10>{{x9, x8, x2, x7, x1, x0, x6, x5, x4, x3}}, in, t[0], range));
}

Task30::Task30(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x0, x1}}, in, t[0], range));
}

Task31::Task31(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task32::Task32(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x1, x0, x2}}, in, t[0], range));
}

Task33::Task33(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x3, x2, x1}}, in, t[0], range));
}

Task34::Task34(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x5 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x3, x5, x0, x4, x2, x1}}, in, t[0], range));
}

Task35::Task35(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x0, x3, x2, x1}}, in, t[0], range));
}

Task36::Task36(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x2, x3, x0, x1}}, in, t[0], range));
}

Task37::Task37(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x0, x1}}, in, t[0], range));
}

Task38::Task38(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x3, x0, x4, x2, x1}}, in, t[0], range));
}

Task39::Task39(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x1, x0, x4, x3, x2}}, in, t[0], range));
}

Task40::Task40(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x3, x2, x0, x1}}, in, t[0], range));
}

Task41::Task41(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x5 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x2, x4, x3, x0, x1}}, in, t[0], range));
}

Task42::Task42(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x3, x2, x1}}, in, t[0], range));
}

Task43::Task43(array<shared_ptr<Tensor>,6> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,5> in = {{t[1], t[2], t[3], t[4], t[5]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x4 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x5, x4, x3, x2, x1}}, in, t[0], range));
}

Task44::Task44(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task45::Task45(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x3, x2, x0}}, in, t[0], range));
}

Task46::Task46(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x1, x2}}, in, t[0], range));
}

Task47::Task47(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task48::Task48(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task49::Task49(array<shared_ptr<Tensor>,2> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

#endif
