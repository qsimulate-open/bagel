//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_gen33.cc
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


#include <src/smith/CASPT2_tasks33.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

Task1600::Task1600(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x1, a4}}, in, t[0], range));
}

Task1601::Task1601(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task1602::Task1602(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1603::Task1603(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, x1, a4}}, in, t[0], range));
}

Task1604::Task1604(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task1605::Task1605(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1606::Task1606(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a1, x1, a4}}, in, t[0], range));
}

Task1607::Task1607(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task1608::Task1608(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1609::Task1609(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, a1, x1, a4}}, in, t[0], range));
}

Task1610::Task1610(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task1611::Task1611(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1612::Task1612(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x0}}, in, t[0], range, e));
}

Task1613::Task1613(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1614::Task1614(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x0}}, in, t[0], range, e));
}

Task1615::Task1615(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1616::Task1616(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range, e));
}

Task1617::Task1617(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1618::Task1618(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, const double e) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range, e));
}

Task1619::Task1619(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1620::Task1620(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x0}}, in, t[0], range));
}

Task1621::Task1621(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1622::Task1622(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x0}}, in, t[0], range));
}

Task1623::Task1623(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1624::Task1624(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task1625::Task1625(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1626::Task1626(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[2]->nblock());
  for (auto& x1 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, c2, a1, x1}}, in, t[0], range));
}

Task1627::Task1627(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1628::Task1628(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, c2, a1, x0}}, in, t[0], range));
}

Task1629::Task1629(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task1630::Task1630(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, a2, c1}}, in, t[0], range));
}

Task1631::Task1631(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task1632::Task1632(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, x2, x5, x4, x3}}, in, t[0], range));
}

Task1633::Task1633(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a2 : *range[2])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, a2, x2}}, in, t[0], range));
}

Task1634::Task1634(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c1 : *range[0])
    for (auto& a2 : *range[2])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, a2, c1}}, in, t[0], range));
}

Task1635::Task1635(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, x2, x5, x4, x3}}, in, t[0], range));
}

Task1636::Task1636(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a1, x0, x1}}, in, t[0], range));
}

Task1637::Task1637(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x2 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x2, c2}}, in, t[0], range));
}

Task1638::Task1638(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x4 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x0 : *range[1])
            for (auto& x1 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x1, x0, x2, x5, x4, x3}}, in, t[0], range));
}

Task1639::Task1639(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x4, a1, x5}}, in, t[0], range));
}

Task1640::Task1640(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task1641::Task1641(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task1642::Task1642(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& x5 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x5, a1, x4}}, in, t[0], range));
}

Task1643::Task1643(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock());
  for (auto& c2 : *range[0])
    for (auto& x3 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x3, c2}}, in, t[0], range));
}

Task1644::Task1644(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task1645::Task1645(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& x4 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x4, a1, x5, x3}}, in, t[0], range));
}

Task1646::Task1646(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& c2 : *range[0])
        for (auto& x4 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x4, c2, a1, x5}}, in, t[0], range));
}

Task1647::Task1647(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task1648::Task1648(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
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

Task1649::Task1649(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& x5 : *range[1])
      for (auto& a1 : *range[2])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, a1, x5, x4}}, in, t[0], range));
}

