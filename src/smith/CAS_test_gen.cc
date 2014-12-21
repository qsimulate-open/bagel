//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_gen.cc
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


#include <src/smith/CAS_test_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

Task0::Task0(vector<shared_ptr<Tensor>> t) {
  r_ =  t[0];
}


Task1::Task1(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task2::Task2(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task3::Task3(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task4::Task4(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task5::Task5(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x1, x3, x2}}, in, t[0], range));
}


Task6::Task6(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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


Task7::Task7(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[3]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& ci0 : *range[3])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        subtasks_.push_back(make_shared<Task_local>(array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task8::Task8(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task9::Task9(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& ci0 : *range[3])
        for (auto& x1 : *range[1])
          for (auto& x2 : *range[1])
            subtasks_.push_back(make_shared<Task_local>(array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range));
}


Task10::Task10(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range));
}


Task11::Task11(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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


Task12::Task12(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a4, c1, a2}}, in, t[0], range));
}


Task13::Task13(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range, e));
}


Task14::Task14(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}


Task15::Task15(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.push_back(make_shared<Task_local>(in, t[0], range));
}


Task16::Task16(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}


Task17::Task17(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.push_back(make_shared<Task_local>(in, t[0], range));
}


Task18::Task18(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a4, c1, a2}}, in, t[0], range));
}


Task19::Task19(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, a2, c3}}, in, t[0], range));
}


Task20::Task20(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c5 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c5, a2}}, in, t[0], range));
}


Task21::Task21(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, a2, c3}}, in, t[0], range));
}


Task22::Task22(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a5 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a5, c3, a2}}, in, t[0], range));
}


Task23::Task23(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, a2, c3}}, in, t[0], range));
}


Task24::Task24(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c1, a2}}, in, t[0], range));
}


Task25::Task25(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task26::Task26(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c1, a2}}, in, t[0], range));
}


Task27::Task27(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task28::Task28(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a3, x0, a1}}, in, t[0], range));
}


Task29::Task29(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task30::Task30(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, c2, a3, c4, a1}}, in, t[0], range));
}


Task31::Task31(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task32::Task32(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, c2, a3, c4, a1}}, in, t[0], range));
}


Task33::Task33(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task34::Task34(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task35::Task35(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x3, x0}}, in, t[0], range));
}


Task36::Task36(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task37::Task37(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x3, x0}}, in, t[0], range));
}


Task38::Task38(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task39::Task39(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c4, a1}}, in, t[0], range));
}


Task40::Task40(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task41::Task41(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c4, a1}}, in, t[0], range));
}


Task42::Task42(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task43::Task43(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task44::Task44(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c2, a3}}, in, t[0], range));
}


Task45::Task45(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task46::Task46(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c2, a3}}, in, t[0], range));
}


Task47::Task47(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task48::Task48(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task49::Task49(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c2, a1}}, in, t[0], range));
}


Task50::Task50(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task51::Task51(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c2, a1}}, in, t[0], range));
}


Task52::Task52(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task53::Task53(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task54::Task54(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a3}}, in, t[0], range));
}


Task55::Task55(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task56::Task56(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task57::Task57(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range, e));
}


Task58::Task58(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task59::Task59(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range, e));
}


Task60::Task60(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task61::Task61(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task62::Task62(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task63::Task63(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task64::Task64(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a2, x0, a1}}, in, t[0], range));
}


Task65::Task65(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}


Task66::Task66(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a1 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x2, x1, a1, c3, a2}}, in, t[0], range));
}


Task67::Task67(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task68::Task68(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}


Task69::Task69(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a3}}, in, t[0], range));
}


Task70::Task70(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task71::Task71(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a2, x0, a1}}, in, t[0], range));
}


Task72::Task72(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}


Task73::Task73(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x5, x0, x4, x1}}, in, t[0], range));
}


Task74::Task74(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}


Task75::Task75(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range, e));
}


Task76::Task76(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}


Task77::Task77(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task78::Task78(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a2, c3, a4}}, in, t[0], range));
}


Task79::Task79(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range, e));
}


Task80::Task80(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}


Task81::Task81(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.push_back(make_shared<Task_local>(in, t[0], range));
}


Task82::Task82(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}


Task83::Task83(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.push_back(make_shared<Task_local>(in, t[0], range));
}


Task84::Task84(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}


Task85::Task85(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c5 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c5, a2}}, in, t[0], range));
}


Task86::Task86(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}


Task87::Task87(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a5 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a5, c3, a2}}, in, t[0], range));
}


Task88::Task88(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}


Task89::Task89(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c1, a2}}, in, t[0], range));
}


Task90::Task90(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task91::Task91(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c1, a2}}, in, t[0], range));
}


Task92::Task92(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task93::Task93(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a1 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a1, c2, a3}}, in, t[0], range));
}


Task94::Task94(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task95::Task95(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, c2, a3, c4, a1}}, in, t[0], range));
}


Task96::Task96(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task97::Task97(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, c2, a3, c4, a1}}, in, t[0], range));
}


Task98::Task98(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task99::Task99(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task100::Task100(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x3, x0}}, in, t[0], range));
}


Task101::Task101(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task102::Task102(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x3, x0}}, in, t[0], range));
}


Task103::Task103(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task104::Task104(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c4, a1}}, in, t[0], range));
}


Task105::Task105(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task106::Task106(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c4, a1}}, in, t[0], range));
}


Task107::Task107(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task108::Task108(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task109::Task109(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c2, a3}}, in, t[0], range));
}


Task110::Task110(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task111::Task111(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c2, a3}}, in, t[0], range));
}


Task112::Task112(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task113::Task113(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task114::Task114(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c2, a1}}, in, t[0], range));
}


Task115::Task115(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task116::Task116(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c2, a1}}, in, t[0], range));
}


Task117::Task117(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task118::Task118(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task119::Task119(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a3}}, in, t[0], range));
}


Task120::Task120(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task121::Task121(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task122::Task122(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range, e));
}


Task123::Task123(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task124::Task124(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range, e));
}


Task125::Task125(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task126::Task126(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task127::Task127(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, c2, a3, a1}}, in, t[0], range));
}


Task128::Task128(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task129::Task129(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a1, x1, a2}}, in, t[0], range));
}


Task130::Task130(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}


Task131::Task131(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a1 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x2, x1, a1, c3, a2}}, in, t[0], range));
}


Task132::Task132(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task133::Task133(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}


Task134::Task134(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x5, x0, x4, x1}}, in, t[0], range));
}


Task135::Task135(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}


Task136::Task136(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a3}}, in, t[0], range));
}


Task137::Task137(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task138::Task138(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}


Task139::Task139(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range, e));
}


Task140::Task140(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}


Task141::Task141(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task142::Task142(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a2, c3, a4}}, in, t[0], range));
}


Task143::Task143(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}


Task144::Task144(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a1 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a1, c2, a3}}, in, t[0], range));
}


Task145::Task145(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c2, a1}}, in, t[0], range));
}


Task146::Task146(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task147::Task147(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c2, a1}}, in, t[0], range));
}


Task148::Task148(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task149::Task149(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a1, x1, a2}}, in, t[0], range));
}


Task150::Task150(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}


Task151::Task151(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task152::Task152(vector<shared_ptr<Tensor>> t) {
  d_ =  t[0];
}


Task153::Task153(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, x1}}, in, t[0], range));
}


Task154::Task154(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task155::Task155(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, c1, a4, c3, a2}}, in, t[0], range));
}


Task156::Task156(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task157::Task157(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, c1, a4, c3, a2}}, in, t[0], range));
}


Task158::Task158(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task159::Task159(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c3 : *range[0])
    for (auto& c5 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c5, c3}}, in, t[0], range));
}


Task160::Task160(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c3 : *range[0])
    for (auto& c5 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c5, c3}}, in, t[0], range));
}


Task161::Task161(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c5 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c5, a2}}, in, t[0], range));
}


Task162::Task162(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a5 : *range[2])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, a5}}, in, t[0], range));
}


Task163::Task163(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a5 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a5, a4}}, in, t[0], range));
}


Task164::Task164(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a5 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a5, c3, a2}}, in, t[0], range));
}


Task165::Task165(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x0 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, c3}}, in, t[0], range));
}


Task166::Task166(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x0 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x0, c3}}, in, t[0], range));
}


Task167::Task167(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c1, a2}}, in, t[0], range));
}


Task168::Task168(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task169::Task169(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c1, a2}}, in, t[0], range));
}


Task170::Task170(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task171::Task171(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& c4 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c4, x1}}, in, t[0], range));
}


Task172::Task172(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c4}}, in, t[0], range));
}


Task173::Task173(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, c2, a3, c4, a1}}, in, t[0], range));
}


Task174::Task174(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task175::Task175(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, c2, a3, c4, a1}}, in, t[0], range));
}


Task176::Task176(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task177::Task177(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x2}}, in, t[0], range));
}


Task178::Task178(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, x1}}, in, t[0], range));
}


Task179::Task179(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x2, x1, a3, c2, a1}}, in, t[0], range));
}


Task180::Task180(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task181::Task181(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x2, x1, a3, c2, a1}}, in, t[0], range));
}


Task182::Task182(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task183::Task183(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c2 : *range[0])
    for (auto& c4 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c4, c2}}, in, t[0], range));
}


Task184::Task184(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c2 : *range[0])
    for (auto& c4 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c4, c2}}, in, t[0], range));
}


Task185::Task185(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c4, a1}}, in, t[0], range));
}


Task186::Task186(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task187::Task187(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c4, a1}}, in, t[0], range));
}


Task188::Task188(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task189::Task189(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a1 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a1, a4}}, in, t[0], range));
}


Task190::Task190(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, a1}}, in, t[0], range));
}


Task191::Task191(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c2, a3}}, in, t[0], range));
}


Task192::Task192(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task193::Task193(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c2, a3}}, in, t[0], range));
}


Task194::Task194(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task195::Task195(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a3 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a3, a4}}, in, t[0], range));
}


Task196::Task196(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, a3}}, in, t[0], range));
}


Task197::Task197(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c2, a1}}, in, t[0], range));
}


Task198::Task198(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task199::Task199(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a4, c2, a1}}, in, t[0], range));
}


Task200::Task200(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task201::Task201(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c2}}, in, t[0], range));
}


Task202::Task202(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c2}}, in, t[0], range));
}


Task203::Task203(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a3}}, in, t[0], range));
}


Task204::Task204(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task205::Task205(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& c3 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, x2}}, in, t[0], range));
}


Task206::Task206(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x2 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, c3}}, in, t[0], range));
}


Task207::Task207(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a1 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x2, x1, a1, c3, a2}}, in, t[0], range));
}


Task208::Task208(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task209::Task209(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x2 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, x3}}, in, t[0], range));
}


Task210::Task210(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x3, x2}}, in, t[0], range));
}


Task211::Task211(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x1, x3, x2, a1, a2}}, in, t[0], range));
}


Task212::Task212(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
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


Task213::Task213(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a2 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a2, a3}}, in, t[0], range));
}


Task214::Task214(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a3 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a3, a2}}, in, t[0], range));
}


Task215::Task215(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a3}}, in, t[0], range));
}


Task216::Task216(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task217::Task217(vector<shared_ptr<Tensor>> t) {
  d2_ =  t[0];
}


Task218::Task218(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a4, c1, a2}}, in, t[0], range));
}


Task219::Task219(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}


Task220::Task220(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a3, x0, a1}}, in, t[0], range));
}


Task221::Task221(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c2, a1}}, in, t[0], range));
}


Task222::Task222(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task223::Task223(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, a3, c2, a1}}, in, t[0], range));
}


Task224::Task224(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}


Task225::Task225(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a2, x0, a1}}, in, t[0], range));
}


Task226::Task226(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x0, x1, a1, a2}}, in, t[0], range));
}


Task227::Task227(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, t[0], range));
}


Task228::Task228(vector<shared_ptr<Tensor>> t) {
  dec_ =  t[0];
}


Task229::Task229(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}


Task230::Task230(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}


Task231::Task231(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, c1, a4, c3, a2}}, in, t[0], range));
}


Task232::Task232(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}


Task233::Task233(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, c1, a4, c3, a2}}, in, t[0], range));
}


Task234::Task234(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}


Task235::Task235(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, c1, a4, c3, a2}}, in, t[0], range));
}


Task236::Task236(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, a4, c1, a2}}, in, t[0], range));
}


Task237::Task237(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task238::Task238(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, a4, c1, a2}}, in, t[0], range));
}


Task239::Task239(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task240::Task240(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, c1, a4, c3, a2}}, in, t[0], range));
}


Task241::Task241(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a4, c1, a2}}, in, t[0], range));
}


Task242::Task242(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task243::Task243(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a4, c1, a2}}, in, t[0], range));
}


Task244::Task244(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task245::Task245(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}


Task246::Task246(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range));
}


Task247::Task247(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x1, x0, c2, a3, c4, a1}}, in, t[0], range));
}


Task248::Task248(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task249::Task249(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x1, x0, c2, a3, c4, a1}}, in, t[0], range));
}


Task250::Task250(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task251::Task251(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range));
}


Task252::Task252(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x3, x0}}, in, t[0], range));
}


Task253::Task253(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range));
}


Task254::Task254(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x3, x0}}, in, t[0], range));
}


Task255::Task255(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range));
}


Task256::Task256(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, a3, c4, a1}}, in, t[0], range));
}


Task257::Task257(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task258::Task258(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, a3, c4, a1}}, in, t[0], range));
}


Task259::Task259(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task260::Task260(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range));
}


Task261::Task261(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, a4, c2, a3}}, in, t[0], range));
}


Task262::Task262(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task263::Task263(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, a4, c2, a3}}, in, t[0], range));
}


Task264::Task264(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task265::Task265(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range));
}


Task266::Task266(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, a4, c2, a1}}, in, t[0], range));
}


Task267::Task267(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task268::Task268(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, a4, c2, a1}}, in, t[0], range));
}


Task269::Task269(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task270::Task270(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range));
}


Task271::Task271(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, x1, a1, a3}}, in, t[0], range));
}


Task272::Task272(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range));
}


Task273::Task273(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range));
}


Task274::Task274(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range, e));
}


Task275::Task275(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range));
}


Task276::Task276(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range, e));
}


Task277::Task277(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range));
}


Task278::Task278(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task279::Task279(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, c2, a3, a1}}, in, t[0], range));
}


Task280::Task280(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task281::Task281(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}


Task282::Task282(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, x1, a1, a2}}, in, t[0], range));
}


Task283::Task283(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a1 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x0, x2, x1, a1, c3, a2}}, in, t[0], range));
}


Task284::Task284(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range));
}


Task285::Task285(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, x1, a1, a2}}, in, t[0], range));
}


Task286::Task286(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x0, x4, x1}}, in, t[0], range));
}


Task287::Task287(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, x1, a1, a2}}, in, t[0], range));
}


Task288::Task288(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, x1, a1, a3}}, in, t[0], range));
}


Task289::Task289(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range));
}


Task290::Task290(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, x1, a1, a2}}, in, t[0], range));
}


Task291::Task291(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range, e));
}


Task292::Task292(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x0, x1, a1, a2}}, in, t[0], range));
}


Task293::Task293(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range));
}


Task294::Task294(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}


Task295::Task295(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range));
}


Task296::Task296(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x1, x0, c2, a3, c4, a1}}, in, t[0], range));
}


Task297::Task297(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task298::Task298(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& c2 : *range[0])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x1, x0, c2, a3, c4, a1}}, in, t[0], range));
}


Task299::Task299(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task300::Task300(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range));
}


Task301::Task301(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a3, c4, a1}}, in, t[0], range));
}


Task302::Task302(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task303::Task303(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c4 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a3, c4, a1}}, in, t[0], range));
}


Task304::Task304(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task305::Task305(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range));
}


Task306::Task306(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a4, c2, a3}}, in, t[0], range));
}


Task307::Task307(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task308::Task308(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a3 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a4, c2, a3}}, in, t[0], range));
}


Task309::Task309(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task310::Task310(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range));
}


Task311::Task311(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a4, c2, a1}}, in, t[0], range));
}


Task312::Task312(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task313::Task313(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a4, c2, a1}}, in, t[0], range));
}


Task314::Task314(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task315::Task315(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range));
}


Task316::Task316(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range, e));
}


Task317::Task317(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range));
}


Task318::Task318(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range, e));
}


Task319::Task319(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range));
}


Task320::Task320(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task321::Task321(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c2, a3, a1}}, in, t[0], range));
}


Task322::Task322(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x1, x0}}, in, t[0], range));
}


Task323::Task323(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}


Task324::Task324(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, a3, c2, a1}}, in, t[0], range));
}


Task325::Task325(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x3, x0}}, in, t[0], range));
}


Task326::Task326(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, a3, c2, a1}}, in, t[0], range));
}


Task327::Task327(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x3, x0}}, in, t[0], range));
}


Task328::Task328(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& a3 : *range[2])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, a3, c2, a1}}, in, t[0], range));
}


Task329::Task329(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, a1, a3}}, in, t[0], range));
}


Task330::Task330(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range));
}


Task331::Task331(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}


Task332::Task332(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, a1, a2}}, in, t[0], range));
}


Task333::Task333(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a1 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x3, x2, x1, a1, c3, a2}}, in, t[0], range));
}


Task334::Task334(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range));
}


Task335::Task335(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, a1, a2}}, in, t[0], range));
}


Task336::Task336(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, a1, a3}}, in, t[0], range));
}


Task337::Task337(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range));
}


Task338::Task338(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, a1, a2}}, in, t[0], range));
}


Task339::Task339(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range, e));
}


Task340::Task340(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, a1, a2}}, in, t[0], range));
}


Task341::Task341(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x0, x2, x1}}, in, t[0], range));
}


Task342::Task342(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}


Task343::Task343(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};

  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x4, a1, a2}}, in, t[0], range));
}


Task344::Task344(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};

  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x0, x4, x1}}, in, t[0], range));
}


