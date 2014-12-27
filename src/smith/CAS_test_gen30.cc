//
// BAGEL - Parallel electron correlation program.
// Filename: CAS_test_gen30.cc
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


#include <src/smith/CAS_test_tasks30.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CAS_test;

Task1450::Task1450(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c2 : *range[0])
    for (auto& c1 : *range[0])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x5, x4, x3, x2, c1, c2}}, in, t[0], range));
}

Task1451::Task1451(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x5, x4, x0, x3, x1, x2}}, in, t[0], range));
}

Task1452::Task1452(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x3 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x4, x3, c1}}, in, t[0], range));
}

Task1453::Task1453(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x3 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x4, x3, c2}}, in, t[0], range));
}

Task1454::Task1454(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x5, x4, x2, x3, x1, x0}}, in, t[0], range));
}

Task1455::Task1455(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x3 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x4, x3, c1}}, in, t[0], range));
}

Task1456::Task1456(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x5, x4, x3, x2, a2, c1}}, in, t[0], range));
}

Task1457::Task1457(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x5, x4, x1, x3, x2, x0}}, in, t[0], range));
}

Task1458::Task1458(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x5, x4, x3, x2, a2, c1}}, in, t[0], range));
}

Task1459::Task1459(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x5, x4, x2, x3, x1, x0}}, in, t[0], range));
}

Task1460::Task1460(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x3 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x4, x3, c1}}, in, t[0], range));
}

Task1461::Task1461(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x5, x4, x2, x3, x1, x0}}, in, t[0], range, e));
}

Task1462::Task1462(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x3 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x4, x3, c1}}, in, t[0], range));
}

Task1463::Task1463(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x0 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x5, x4, x0, x3, x2, x1}}, in, t[0], range));
}

Task1464::Task1464(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x3 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x4, x3, c1}}, in, t[0], range));
}

Task1465::Task1465(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x5, x4, x2, x3, x1, x0}}, in, t[0], range));
}

Task1466::Task1466(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}

Task1467::Task1467(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x5 : *range[1])
      for (auto& x6 : *range[1])
        for (auto& x7 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x7, x6, x5, c1}}, in, t[0], range));
}

Task1468::Task1468(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x6 : *range[1])
            for (auto& x7 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x7, x6, x2, x5, x1, x0}}, in, t[0], range));
}

Task1469::Task1469(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}

Task1470::Task1470(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, x1, c1}}, in, t[0], range));
}

Task1471::Task1471(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x3, x2, x1, c2, a3, c1}}, in, t[0], range));
}

Task1472::Task1472(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, x0, x1}}, in, t[0], range));
}

Task1473::Task1473(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a3 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x3, x2, x1, c2, a3, c1}}, in, t[0], range));
}

Task1474::Task1474(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, x0, x1}}, in, t[0], range));
}

Task1475::Task1475(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, x1, c1}}, in, t[0], range));
}

Task1476::Task1476(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, x0, x1}}, in, t[0], range));
}

Task1477::Task1477(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}

Task1478::Task1478(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, c1, c3, a2}}, in, t[0], range));
}

Task1479::Task1479(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c3 : *range[0])
    for (auto& c1 : *range[0])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, c1, c3}}, in, t[0], range));
}

Task1480::Task1480(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, x3, x0, x2}}, in, t[0], range));
}

Task1481::Task1481(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, c1, c3, a2}}, in, t[0], range));
}

Task1482::Task1482(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x3 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x3, c3}}, in, t[0], range));
}

Task1483::Task1483(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x2, x3, x1, x0}}, in, t[0], range));
}

Task1484::Task1484(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, c1, c3, a2}}, in, t[0], range));
}

Task1485::Task1485(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& x3 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x3, c1}}, in, t[0], range));
}

Task1486::Task1486(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x2, x3, x1, x0}}, in, t[0], range));
}

Task1487::Task1487(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, c1, c3, a2}}, in, t[0], range));
}

Task1488::Task1488(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x3}}, in, t[0], range));
}

Task1489::Task1489(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, c1, c3, a2}}, in, t[0], range));
}

Task1490::Task1490(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x3}}, in, t[0], range));
}

Task1491::Task1491(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, c1, c3, a2}}, in, t[0], range));
}

Task1492::Task1492(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, a2, c3}}, in, t[0], range));
}

Task1493::Task1493(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, x3, x2, x0}}, in, t[0], range));
}

Task1494::Task1494(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, a2, c3}}, in, t[0], range));
}

Task1495::Task1495(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x2, x3, x1, x0}}, in, t[0], range));
}

Task1496::Task1496(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& c1 : *range[0])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, c1, c3, a2}}, in, t[0], range));
}

Task1497::Task1497(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x2, x3, a2, c1}}, in, t[0], range));
}

Task1498::Task1498(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x2, x3, x1, x0}}, in, t[0], range));
}

Task1499::Task1499(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x2, x3, a2, c1}}, in, t[0], range));
}

