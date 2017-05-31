//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_contract_gen.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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


#include <src/smith/caspt2/MSCASPT2.h>
#include <src/smith/caspt2/MSCASPT2_contract_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

Task900::Task900(vector<shared_ptr<Tensor>> t, const bool reset) : reset_(reset) {
  deci_ = t[0];
}

Task901::Task901(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task902::Task902(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task903::Task903(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task914::Task914(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range, shared_ptr<const CIWfn> ciwfn, shared_ptr<VectorB> bdata, const size_t offset, const size_t size) {
  array<shared_ptr<const Tensor>,1> in = {{t[1]}};
  out_ = t[0];
  in_ = in;
  subtasks_.push_back(make_shared<Task_local>(in, t[0], range, ciwfn, bdata, offset, size));
}


Task915::Task915(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,4> in = {{t[1], t[2], t[3], t[4]}};
  out_ = t[0];
  in_ = in;

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[4]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}


Task916::Task916(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}


Task917::Task917(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[3]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}


Task918::Task918(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,3> in = {{t[1], t[2], t[3]}};
  out_ = t[0];
  in_ = in;

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[3]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}

Task921::Task921(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;

  subtasks_.reserve(range[1]->nblock() * range[1]->nblock() * range[1]->nblock() * range[1]->nblock() * range[1]->nblock() * range[1]->nblock());
    for (auto& x0 : *range[1])
      for (auto& x1 : *range[1])
        for (auto& x2 : *range[1])
          for (auto& x3 : *range[1])
            for (auto& x4 : *range[1])
              for (auto& x5 : *range[1])
                if (t[0]->is_local(x0, x1, x2, x3, x4, x5))
                  subtasks_.push_back(make_shared<Task_local>(array<const Index,6>{{x0, x1, x2, x3, x4, x5}}, in, t[0], range));
}

Task923::Task923(vector<shared_ptr<Tensor>> t, array<shared_ptr<const IndexRange>,4> range) {
  array<shared_ptr<const Tensor>,2> in = {{t[1], t[2]}};
  out_ = t[0];
  in_ = in;

  subtasks_.reserve(range[3]->nblock());
  for (auto& ci0 : *range[3])
    if (t[0]->is_local(ci0))
      subtasks_.push_back(make_shared<Task_local>(array<const Index,1>{{ci0}}, in, t[0], range));
}



#endif
