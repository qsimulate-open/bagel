//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_gen31.cc
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


#include <src/smith/CASPT2_tasks31.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

Task1500::Task1500(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x2, x3, x1, x0}}, in, t[0], range));
}

Task1501::Task1501(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}

Task1502::Task1502(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a2, c3, c1}}, in, t[0], range));
}

Task1503::Task1503(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c4 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c4, a2, c3}}, in, t[0], range));
}

Task1504::Task1504(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x1}}, in, t[0], range));
}

Task1505::Task1505(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c4 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c4, a2, c3}}, in, t[0], range));
}

Task1506::Task1506(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x1}}, in, t[0], range));
}

Task1507::Task1507(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a2, c3, c1}}, in, t[0], range));
}

Task1508::Task1508(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c4 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c4, a2, c1}}, in, t[0], range));
}

Task1509::Task1509(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x1}}, in, t[0], range));
}

Task1510::Task1510(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c4 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c4, a2, c1}}, in, t[0], range));
}

Task1511::Task1511(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x1}}, in, t[0], range));
}

Task1512::Task1512(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a2, c3, c1}}, in, t[0], range));
}

Task1513::Task1513(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a4 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c3, a4, c1}}, in, t[0], range));
}

Task1514::Task1514(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x1}}, in, t[0], range));
}

Task1515::Task1515(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a4 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, c3, a4, c1}}, in, t[0], range));
}

Task1516::Task1516(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x1}}, in, t[0], range));
}

Task1517::Task1517(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a2, c3, c1}}, in, t[0], range));
}

Task1518::Task1518(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x0, x1, c1, a4, c3, a2}}, in, t[0], range));
}

Task1519::Task1519(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x1}}, in, t[0], range));
}

Task1520::Task1520(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x0, x1, c1, a4, c3, a2}}, in, t[0], range));
}

Task1521::Task1521(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x1}}, in, t[0], range));
}

Task1522::Task1522(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a2, c3, c1}}, in, t[0], range));
}

Task1523::Task1523(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x1}}, in, t[0], range, e));
}

Task1524::Task1524(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a2, c3, c1}}, in, t[0], range));
}

Task1525::Task1525(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x1}}, in, t[0], range, e));
}

Task1526::Task1526(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a2, c3, c1}}, in, t[0], range));
}

Task1527::Task1527(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x1}}, in, t[0], range));
}

Task1528::Task1528(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c1 : *range[0])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x1, a2, c3, c1}}, in, t[0], range));
}

Task1529::Task1529(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& ci0 : *range[3])
        subtasks_.push_back(make_shared<Task_local>(std::array<const Index,3>{{ci0, x0, x1}}, in, t[0], range));
}

Task1530::Task1530(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}

Task1531::Task1531(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x4, c2, a1}}, in, t[0], range));
}

Task1532::Task1532(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x4 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x3, x4, c2}}, in, t[0], range));
}

Task1533::Task1533(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x5, x3, x2, x4, x1, x0}}, in, t[0], range));
}

Task1534::Task1534(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x4, c2, a1}}, in, t[0], range));
}

Task1535::Task1535(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x0, x1, x4}}, in, t[0], range));
}

Task1536::Task1536(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x4, c2, a1}}, in, t[0], range));
}

Task1537::Task1537(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x4, x1, x0}}, in, t[0], range));
}

Task1538::Task1538(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      for (auto& x4 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x4, c2, a1}}, in, t[0], range));
}

Task1539::Task1539(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x4 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x5, x3, x4, a1}}, in, t[0], range));
}

Task1540::Task1540(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x5 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x5, x0, x3, x4, x2, x1}}, in, t[0], range));
}

Task1541::Task1541(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(std::array<const Index,1>{{ci0}}, in, t[0], range));
}

Task1542::Task1542(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c2 : *range[0])
    for (auto& a1 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, a1, c2}}, in, t[0], range));
}

Task1543::Task1543(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c2 : *range[0])
    for (auto& a1 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& x2 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x3, x1, x2, c3, a1, c2}}, in, t[0], range));
}

Task1544::Task1544(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x1, x0, x2}}, in, t[0], range));
}

Task1545::Task1545(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c2 : *range[0])
    for (auto& a1 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& x2 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x3 : *range[1])
              for (auto& ci0 : *range[3])
                subtasks_.push_back(make_shared<Task_local>(std::array<const Index,7>{{ci0, x3, x1, x2, c3, a1, c2}}, in, t[0], range));
}

Task1546::Task1546(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, x0, x1}}, in, t[0], range));
}

Task1547::Task1547(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c2 : *range[0])
    for (auto& a1 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, a1, c2}}, in, t[0], range));
}

Task1548::Task1548(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a1 : *range[2])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x2, a1, c3}}, in, t[0], range));
}

Task1549::Task1549(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[3]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& ci0 : *range[3])
            subtasks_.push_back(make_shared<Task_local>(std::array<const Index,5>{{ci0, x3, x0, x1, x2}}, in, t[0], range));
}

