//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_gen22.cc
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

#include <src/smith/MRCI_tasks22.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

Task1050::Task1050(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c4 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c4, x0}}, in, t[0], range));
}

Task1051::Task1051(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c4}}, in, t[0], range));
}

Task1052::Task1052(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c4 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c4, x0}}, in, t[0], range));
}

Task1053::Task1053(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& c4 : *range[0])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, c4, x2, x1}}, in, t[0], range));
}

Task1054::Task1054(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c4 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c4, x0}}, in, t[0], range));
}

Task1055::Task1055(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, c4}}, in, t[0], range));
}

Task1056::Task1056(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a1}}, in, t[0], range));
}

Task1057::Task1057(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c4 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c4, x0}}, in, t[0], range));
}

Task1058::Task1058(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c4}}, in, t[0], range));
}

Task1059::Task1059(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c4 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c4, x0}}, in, t[0], range));
}

Task1060::Task1060(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& c4 : *range[0])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, c4, x2, x1}}, in, t[0], range));
}

Task1061::Task1061(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c4 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c4, x0}}, in, t[0], range));
}

Task1062::Task1062(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, c4}}, in, t[0], range));
}

Task1063::Task1063(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a1}}, in, t[0], range));
}

Task1064::Task1064(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1065::Task1065(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c4 : *range[0])
    for (auto& c2 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c2, c4}}, in, t[0], range));
}

Task1066::Task1066(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1067::Task1067(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock());
  for (auto& c4 : *range[0])
    for (auto& c2 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c2, c4}}, in, t[0], range));
}

Task1068::Task1068(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1069::Task1069(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, a1}}, in, t[0], range));
}

Task1070::Task1070(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1071::Task1071(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, a1}}, in, t[0], range));
}

Task1072::Task1072(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1073::Task1073(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, a3}}, in, t[0], range));
}

Task1074::Task1074(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1075::Task1075(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, a3}}, in, t[0], range));
}

Task1076::Task1076(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1077::Task1077(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a4 : *range[2])
      for (auto& c5 : *range[0])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, c5, a4, a1}}, in, t[0], range));
}

Task1078::Task1078(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1079::Task1079(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a4 : *range[2])
      for (auto& c5 : *range[0])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, c5, a4, a1}}, in, t[0], range));
}

Task1080::Task1080(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1081::Task1081(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a4 : *range[2])
      for (auto& c5 : *range[0])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, c5, a4, a3}}, in, t[0], range));
}

Task1082::Task1082(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1083::Task1083(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a4 : *range[2])
      for (auto& c5 : *range[0])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, c5, a4, a3}}, in, t[0], range));
}

Task1084::Task1084(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1085::Task1085(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& a5 : *range[2])
      for (auto& a1 : *range[2])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a1, a5, c4}}, in, t[0], range));
}

Task1086::Task1086(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1087::Task1087(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& a5 : *range[2])
      for (auto& a1 : *range[2])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a1, a5, c4}}, in, t[0], range));
}

Task1088::Task1088(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1089::Task1089(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& a5 : *range[2])
      for (auto& a3 : *range[2])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a3, a5, c4}}, in, t[0], range));
}

Task1090::Task1090(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1091::Task1091(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& a5 : *range[2])
      for (auto& a3 : *range[2])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a3, a5, c4}}, in, t[0], range));
}

Task1092::Task1092(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1093::Task1093(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a4 : *range[2])
      for (auto& c5 : *range[0])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, c5, a4, a1}}, in, t[0], range));
}

Task1094::Task1094(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1095::Task1095(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a4 : *range[2])
      for (auto& c5 : *range[0])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, c5, a4, a1}}, in, t[0], range));
}

Task1096::Task1096(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1097::Task1097(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a4 : *range[2])
      for (auto& c5 : *range[0])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, c5, a4, a3}}, in, t[0], range));
}

Task1098::Task1098(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, a3, a1}}, in, t[0], range));
}

Task1099::Task1099(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a3 : *range[2])
    for (auto& a4 : *range[2])
      for (auto& c5 : *range[0])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, c5, a4, a3}}, in, t[0], range));
}

#endif
