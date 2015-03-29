//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_gen27.cc
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

#include <src/smith/MRCI_tasks27.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

Task1300::Task1300(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, a2}}, in, t[0], range));
}

Task1301::Task1301(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c3 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, x5, x0, x1}}, in, t[0], range));
}

Task1302::Task1302(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& x4 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x4, c3, x3, x2}}, in, t[0], range));
}

Task1303::Task1303(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x5 : *range[1])
        for (auto& c3 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, x5, x0, x1}}, in, t[0], range));
}

Task1304::Task1304(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x4 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x4, x3, x2, c3}}, in, t[0], range));
}

Task1305::Task1305(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, a2}}, in, t[0], range));
}

Task1306::Task1306(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x3 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, a2, x3, a1}}, in, t[0], range));
}

Task1307::Task1307(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a3 : *range[2])
      for (auto& c4 : *range[0])
        for (auto& x2 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, c4, a3, a2}}, in, t[0], range));
}

Task1308::Task1308(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, a2}}, in, t[0], range));
}

Task1309::Task1309(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& x4 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& x5 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& a2 : *range[2])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{a2, x3, x2, x5, a1, x4}}, in, t[0], range));
}

Task1310::Task1310(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, a2, x3, x2}}, in, t[0], range));
}

Task1311::Task1311(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, a2}}, in, t[0], range));
}

Task1312::Task1312(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& x5 : *range[1])
        for (auto& a2 : *range[2])
          for (auto& x2 : *range[1])
            for (auto& x3 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x3, x2, a2, x5, a1, x4}}, in, t[0], range));
}

Task1313::Task1313(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x2 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& a3 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a3, x3, x2, a2}}, in, t[0], range));
}

Task1314::Task1314(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& a2 : *range[2])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& a1 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a1, x0, x1, a2}}, in, t[0], range));
}

Task1315::Task1315(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x4 : *range[1])
    for (auto& a1 : *range[2])
      for (auto& x5 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& a2 : *range[2])
            for (auto& x3 : *range[1])
              subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x3, a2, x2, x5, a1, x4}}, in, t[0], range));
}

Task1316::Task1316(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a3 : *range[2])
      for (auto& a2 : *range[2])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, a2, a3, x2}}, in, t[0], range));
}

Task1317::Task1317(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c1 : *range[0])
      for (auto& x1 : *range[1])
        for (auto& c2 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c2, x1, c1, x0}}, in, t[0], range));
}

Task1318::Task1318(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c2, x1, x0}}, in, t[0], range));
}

Task1319::Task1319(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c2, x3, x2}}, in, t[0], range));
}

Task1320::Task1320(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& c3 : *range[0])
    for (auto& c2 : *range[0])
      for (auto& c4 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c4, c2, c3}}, in, t[0], range));
}

Task1321::Task1321(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c2, x3, x2}}, in, t[0], range));
}

Task1322::Task1322(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[2]->nblock()*range[1]->nblock()*range[2]->nblock());
  for (auto& x2 : *range[1])
    for (auto& a3 : *range[2])
      for (auto& x3 : *range[1])
        for (auto& a4 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a4, x3, a3, x2}}, in, t[0], range));
}

Task1323::Task1323(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c2, x1, x0}}, in, t[0], range));
}

Task1324::Task1324(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x4 : *range[1])
    for (auto& c2 : *range[0])
      for (auto& x5 : *range[1])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x5, c2, x4}}, in, t[0], range));
}

Task1325::Task1325(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c2 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c2, x1, x0}}, in, t[0], range));
}

Task1326::Task1326(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x6 : *range[1])
    for (auto& c2 : *range[0])
      for (auto& x7 : *range[1])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x7, c2, x6}}, in, t[0], range));
}

Task1327::Task1327(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c1 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c3 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c3, a4, c1, a2}}, in, t[0], range));
}

Task1328::Task1328(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, a2, a4}}, in, t[0], range));
}

Task1329::Task1329(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x1, x0}}, in, t[0], range));
}

Task1330::Task1330(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& c3 : *range[0])
      for (auto& x3 : *range[1])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x3, c3, x2}}, in, t[0], range));
}

Task1331::Task1331(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, a2, a4}}, in, t[0], range));
}

Task1332::Task1332(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& c5 : *range[0])
    for (auto& c3 : *range[0])
      for (auto& c6 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c6, c3, c5}}, in, t[0], range));
}

Task1333::Task1333(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, a2, a4}}, in, t[0], range));
}

Task1334::Task1334(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[0]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& c5 : *range[0])
    for (auto& c3 : *range[0])
      for (auto& c6 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c6, c3, c5}}, in, t[0], range));
}

Task1335::Task1335(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, a2, a4}}, in, t[0], range));
}

Task1336::Task1336(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a5 : *range[2])
      for (auto& a2 : *range[2])
        for (auto& a6 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a6, a2, a5, a4}}, in, t[0], range));
}

Task1337::Task1337(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, a2, a4}}, in, t[0], range));
}

Task1338::Task1338(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[2]->nblock()*range[2]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a5 : *range[2])
      for (auto& a2 : *range[2])
        for (auto& a6 : *range[2])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{a6, a2, a5, a4}}, in, t[0], range));
}

Task1339::Task1339(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, a2, a4}}, in, t[0], range));
}

Task1340::Task1340(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x3 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, x3, x2}}, in, t[0], range));
}

Task1341::Task1341(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[1]->nblock()*range[0]->nblock()*range[1]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& c3 : *range[0])
      for (auto& x1 : *range[1])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, x1, c3, x0}}, in, t[0], range));
}

Task1342::Task1342(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, a2, a4}}, in, t[0], range));
}

Task1343::Task1343(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}

Task1344::Task1344(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& a2 : *range[2])
      for (auto& c3 : *range[0])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, c3, a2, a4}}, in, t[0], range));
}

Task1345::Task1345(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}

Task1346::Task1346(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[2]->nblock()*range[1]->nblock()*range[2]->nblock()*range[1]->nblock());
  for (auto& a1 : *range[2])
    for (auto& x0 : *range[1])
      for (auto& a2 : *range[2])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, a2, x0, a1}}, in, t[0], range));
}

Task1347::Task1347(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[2]->nblock()*range[2]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& a2 : *range[2])
    for (auto& a1 : *range[2])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x1, x0, a1, a2}}, in, t[0], range));
}

Task1348::Task1348(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[0]->nblock()*range[0]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& c3 : *range[0])
        for (auto& c4 : *range[0])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c4, c3, x1, x0}}, in, t[0], range));
}

Task1349::Task1349(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  subtasks_.reserve(range[0]->nblock()*range[1]->nblock()*range[0]->nblock()*range[1]->nblock());
  for (auto& c3 : *range[0])
    for (auto& x2 : *range[1])
      for (auto& c4 : *range[0])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, c4, x2, c3}}, in, t[0], range));
}

#endif
