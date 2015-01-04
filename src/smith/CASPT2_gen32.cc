//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_gen32.cc
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


#include <src/smith/CASPT2_tasks32.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

Task1550::Task1550(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a4}}, in, t[0], range));
}

Task1551::Task1551(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x0}}, in, t[0], range));
}

Task1552::Task1552(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1553::Task1553(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x0, a4}}, in, t[0], range));
}

Task1554::Task1554(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x0}}, in, t[0], range));
}

Task1555::Task1555(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1556::Task1556(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a1, x0, a4}}, in, t[0], range));
}

Task1557::Task1557(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x0}}, in, t[0], range));
}

Task1558::Task1558(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1559::Task1559(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a1, x0, a4}}, in, t[0], range));
}

Task1560::Task1560(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x0}}, in, t[0], range));
}

Task1561::Task1561(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1562::Task1562(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c2, a1}}, in, t[0], range));
}

Task1563::Task1563(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, c3}}, in, t[0], range));
}

Task1564::Task1564(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a1 : *range[2])
    for (auto& c2 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c2, a1}}, in, t[0], range));
}

Task1565::Task1565(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, c3}}, in, t[0], range));
}

Task1566::Task1566(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1567::Task1567(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, a2}}, in, t[0], range));
}

Task1568::Task1568(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, c3}}, in, t[0], range));
}

Task1569::Task1569(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c1, a2}}, in, t[0], range));
}

Task1570::Task1570(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c3 : *range[0])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, c3}}, in, t[0], range));
}

Task1571::Task1571(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1572::Task1572(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, c1}}, in, t[0], range));
}

Task1573::Task1573(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, a2, c1}}, in, t[0], range));
}

Task1574::Task1574(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1575::Task1575(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a2, c1}}, in, t[0], range));
}

Task1576::Task1576(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, a2, c1}}, in, t[0], range));
}

Task1577::Task1577(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1578::Task1578(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a4 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a4, c1}}, in, t[0], range));
}

Task1579::Task1579(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, a2, c1}}, in, t[0], range));
}

Task1580::Task1580(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1581::Task1581(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{a2, c1}}, in, t[0], range));
}

Task1582::Task1582(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, a2, c1}}, in, t[0], range));
}

Task1583::Task1583(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1584::Task1584(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c3 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, x0}}, in, t[0], range));
}

Task1585::Task1585(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, a2, c1}}, in, t[0], range));
}

Task1586::Task1586(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c3 : *range[0])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{c3, x0}}, in, t[0], range));
}

Task1587::Task1587(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& a4 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, c3, a2, c1}}, in, t[0], range));
}

Task1588::Task1588(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1589::Task1589(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c4}}, in, t[0], range));
}

Task1590::Task1590(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task1591::Task1591(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c4 : *range[0])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, c4}}, in, t[0], range));
}

Task1592::Task1592(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task1593::Task1593(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1594::Task1594(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& c4 : *range[0])
    for (auto& x1 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, a1, x1, c4}}, in, t[0], range));
}

Task1595::Task1595(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task1596::Task1596(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1597::Task1597(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& c4 : *range[0])
    for (auto& x1 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, a1, x1, c4}}, in, t[0], range));
}

Task1598::Task1598(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task1599::Task1599(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

