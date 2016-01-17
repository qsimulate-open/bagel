//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_gen3.cc
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

#include <src/smith/mrci/MRCI_tasks3.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

Task100::Task100(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task101::Task101(array<shared_ptr<Tensor>,2> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task102::Task102(array<shared_ptr<Tensor>,6> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task103::Task103(array<shared_ptr<Tensor>,7> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task104::Task104(array<shared_ptr<Tensor>,3> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task105::Task105(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x3 : *range[1])
          subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x1, x0}}, in, t[0], range));
}

Task106::Task106(array<shared_ptr<Tensor>,4> t, array<shared_ptr<const IndexRange>,3> range) {
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

Task107::Task107(array<shared_ptr<Tensor>,5> t, array<shared_ptr<const IndexRange>,3> range) {
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
