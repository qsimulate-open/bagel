//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_gen12.cc
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

#include <src/smith/caspt2/MSCASPT2_tasks12.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

Task550::Task550(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x0, x1}}, in, t[0], range));
}

Task551::Task551(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x0, x1}}, in, t[0], range));
}

Task552::Task552(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x0, x1}}, in, t[0], range));
}

Task553::Task553(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x2 : *range[1])
          if (t[5]->is_local(x3, x1, x0, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x1, x0, x2}}, in, out, range));
}

Task554::Task554(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x0, x1}}, in, t[0], range));
}

Task555::Task555(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          if (t[5]->is_local(x3, x0, x1, x2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x1, x2}}, in, out, range));
}

Task556::Task556(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x0, x1}}, in, t[0], range));
}

Task557::Task557(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x1 : *range[1])
          if (t[5]->is_local(x3, x2, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x2, x0, x1}}, in, out, range));
}

Task558::Task558(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x0, x1}}, in, t[0], range));
}

Task559::Task559(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x0, x1}}, in, t[0], range));
}

Task560::Task560(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x2 : *range[1])
            for (auto& x1 : *range[1])
              if (t[5]->is_local(x5, x0, x4, x3, x2, x1))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x0, x4, x3, x2, x1}}, in, out, range));
}

Task561::Task561(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x3 : *range[1])
              if (t[0]->is_local(x3, x4, x5, x0, x1, x2))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x3, x4, x5, x0, x1, x2}}, in, t[0], range));
}

Task562::Task562(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x5 : *range[1])
    for (auto& x2 : *range[1])
      for (auto& x4 : *range[1])
        for (auto& x3 : *range[1])
          for (auto& x1 : *range[1])
            for (auto& x0 : *range[1])
              if (t[5]->is_local(x5, x2, x4, x3, x1, x0))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x5, x2, x4, x3, x1, x0}}, in, out, range));
}

Task563::Task563(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x2 : *range[1])
    for (auto& x1 : *range[1])
      for (auto& x0 : *range[1])
        for (auto& x5 : *range[1])
          for (auto& x4 : *range[1])
            for (auto& x3 : *range[1])
              if (t[0]->is_local(x3, x4, x5, x0, x1, x2))
                subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x3, x4, x5, x0, x1, x2}}, in, t[0], range));
}

Task564::Task564(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.push_back(make_shared<Task_local>(in, t[0], range));
}

Task565::Task565(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a2 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a4 : *range[2])
        for (auto& c1 : *range[0])
          if (t[1]->is_local(c1, a4, c3, a2))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a4, c3, a2}}, in, t[0], range));
}

Task566::Task566(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.push_back(make_shared<Task_local>(in, t[0], range));
}

Task567::Task567(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[2]->nblock()*range[0]->nblock()*range[2]->nblock()*range[0]->nblock());
  for (auto& a4 : *range[2])
    for (auto& c3 : *range[0])
      for (auto& a2 : *range[2])
        for (auto& c1 : *range[0])
          if (t[1]->is_local(c1, a2, c3, a4))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{c1, a2, c3, a4}}, in, t[0], range));
}

Task568::Task568(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      if (t[5]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, out, range));
}

Task569::Task569(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task570::Task570(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task571::Task571(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task572::Task572(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock());
  for (auto& x0 : *range[1])
    for (auto& x1 : *range[1])
      if (t[0]->is_local(x1, x0))
        subtasks_.push_back(make_shared<Task_local>(array<const Index,2>{{x1, x0}}, in, t[0], range));
}

Task573::Task573(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,1> in = {{t[5]}};
  array<shared_ptr<Tensor>,5> out = {{t[0], t[1], t[2], t[3], t[4]}};
  in_ = in;
  out_ = out;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x3 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x2 : *range[1])
        for (auto& x1 : *range[1])
          if (t[5]->is_local(x3, x0, x2, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x3, x0, x2, x1}}, in, out, range));
}

Task574::Task574(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x0, x1}}, in, t[0], range));
}

Task575::Task575(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,3> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;
  subtasks_.reserve(range[1]->nblock()*range[1]->nblock()*range[1]->nblock()*range[1]->nblock());
  for (auto& x1 : *range[1])
    for (auto& x0 : *range[1])
      for (auto& x3 : *range[1])
        for (auto& x2 : *range[1])
          if (t[0]->is_local(x2, x3, x0, x1))
            subtasks_.push_back(make_shared<Task_local>(array<const Index,4>{{x2, x3, x0, x1}}, in, t[0], range));
}

#endif
