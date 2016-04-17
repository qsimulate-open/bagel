//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_gen1.cc
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

#include <src/smith/caspt2/MSCASPT2_tasks1.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

Task0::Task0(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x2 : *range[1])
              if (t[0]->is_local(x0, x5, x1, x4, x3, x2))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x5, x1, x4, x3, x2}}, in, t[0], range));
}

Task1::Task1(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x2 : *range[1])
              if (t[0]->is_local(x5, x0, x1, x4, x3, x2))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x1, x4, x3, x2}}, in, t[0], range));
}

Task2::Task2(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              if (t[0]->is_local(x5, x4, x3, x2, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x3, x2, x1, x0}}, in, t[0], range));
}

Task3::Task3(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x2 : *range[1])
              if (t[0]->is_local(x5, x0, x4, x1, x3, x2))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x1, x3, x2}}, in, t[0], range));
}

Task4::Task4(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x0, x3, x1, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x3, x1, x2}}, in, t[0], range));
}

Task5::Task5(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x3, x0, x1, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x1, x2}}, in, t[0], range));
}

Task6::Task6(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x3, x2, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task7::Task7(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              if (t[0]->is_local(x5, x4, x0, x3, x1, x2))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x0, x3, x1, x2}}, in, t[0], range));
}

Task8::Task8(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x2 : *range[1])
              if (t[0]->is_local(x5, x0, x4, x3, x1, x2))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x3, x1, x2}}, in, t[0], range));
}

Task9::Task9(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x1, x3, x0, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x3, x0, x2}}, in, t[0], range));
}

Task10::Task10(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              if (t[0]->is_local(x2, x5, x3, x4, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x2, x5, x3, x4, x1, x0}}, in, t[0], range));
}

Task11::Task11(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x1 : *range[1])
              if (t[0]->is_local(x5, x0, x3, x4, x2, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x3, x4, x2, x1}}, in, t[0], range));
}

Task12::Task12(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x1 : *range[1])
              if (t[0]->is_local(x5, x4, x3, x0, x2, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x3, x0, x2, x1}}, in, t[0], range));
}

Task13::Task13(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x7 : *range[1])
    for (auto& x6 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& x1 : *range[1])
                for (auto& x0 : *range[1])
                  if (t[0]->is_local(x7, x6, x2, x5, x4, x3, x1, x0))
                    subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x6, x2, x5, x4, x3, x1, x0}}, in, t[0], range));
}

Task14::Task14(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x7 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x6 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& x2 : *range[1])
                for (auto& x1 : *range[1])
                  if (t[0]->is_local(x7, x0, x6, x5, x4, x3, x2, x1))
                    subtasks_.push_back(make_shared<Task_local>(array<const Index,8>{{x7, x0, x6, x5, x4, x3, x2, x1}}, in, t[0], range));
}

Task15::Task15(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              if (t[0]->is_local(x5, x4, x2, x3, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x2, x3, x1, x0}}, in, t[0], range));
}

Task16::Task16(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x2, x3, x1, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x1, x0}}, in, t[0], range));
}

Task17::Task17(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x3, x0, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}

Task18::Task18(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x4 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              if (t[0]->is_local(x5, x3, x2, x4, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x3, x2, x4, x1, x0}}, in, t[0], range));
}

Task19::Task19(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x1 : *range[1])
              if (t[0]->is_local(x5, x0, x4, x3, x2, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x3, x2, x1}}, in, t[0], range));
}

Task20::Task20(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x3, x2, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x0, x1}}, in, t[0], range));
}

Task21::Task21(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task22::Task22(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x1 : *range[1])
          if (t[0]->is_local(x0, x3, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x3, x2, x1}}, in, t[0], range));
}

Task23::Task23(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x0, x1))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}

Task24::Task24(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x3, x1, x0, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x1, x0, x2}}, in, t[0], range));
}

Task25::Task25(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x0 : *range[1])
              if (t[0]->is_local(x5, x4, x1, x3, x2, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x4, x1, x3, x2, x0}}, in, t[0], range));
}

Task26::Task26(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x0 : *range[1])
          if (t[0]->is_local(x1, x3, x2, x0))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x3, x2, x0}}, in, t[0], range));
}

Task27::Task27(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              if (t[0]->is_local(x5, x2, x4, x3, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x2, x4, x3, x1, x0}}, in, t[0], range));
}

Task28::Task28(vector<shared_ptr<Tensor>> t, const bool reset) : reset_(reset) {
  den2_ =  t[0];
}

Task29::Task29(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x2 : *range[1])
      if (t[0]->is_local(x2, x3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, x3}}, in, t[0], range));
}

Task30::Task30(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      if (t[0]->is_local(x3, x2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x3, x2}}, in, t[0], range));
}

Task31::Task31(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task32::Task32(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      if (t[0]->is_local(x3, x2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x3, x2}}, in, t[0], range));
}

Task33::Task33(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task34::Task34(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      if (t[0]->is_local(x3, x2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x3, x2}}, in, t[0], range));
}

Task35::Task35(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task36::Task36(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task37::Task37(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task38::Task38(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      if (t[0]->is_local(x3, x2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x3, x2}}, in, t[0], range));
}

Task39::Task39(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task40::Task40(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c2 : *range[0])
    for (auto& c3 : *range[0])
      if (t[0]->is_local(c3, c2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, c2}}, in, t[0], range));
}

Task41::Task41(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c3 : *range[0])
    for (auto& c2 : *range[0])
      if (t[0]->is_local(c2, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c2, c3}}, in, t[0], range));
}

Task42::Task42(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& c1 : *range[0])
        for (auto& c2 : *range[0])
          if (t[0]->is_local(c2, c1, x3, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, c1, x3, x2}}, in, t[0], range));
}

Task43::Task43(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c3 : *range[0])
    for (auto& c2 : *range[0])
      if (t[0]->is_local(c2, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c2, c3}}, in, t[0], range));
}

Task44::Task44(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task45::Task45(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c3 : *range[0])
    for (auto& c2 : *range[0])
      if (t[0]->is_local(c2, c3))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c2, c3}}, in, t[0], range));
}

Task46::Task46(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task47::Task47(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x2 : *range[1])
      if (t[0]->is_local(x2, c2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, c2}}, in, t[0], range));
}

Task48::Task48(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x2 : *range[1])
      if (t[0]->is_local(x2, c2))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, c2}}, in, t[0], range));
}

Task49::Task49(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& c1 : *range[0])
          if (t[0]->is_local(c1, x0, x1, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x0, x1, x2}}, in, t[0], range));
}

#endif
