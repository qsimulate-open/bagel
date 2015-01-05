//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_density2qq.cc
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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_density2q() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density2q = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor42 = {Den1};
  auto task42 = make_shared<Task42>(tensor42);
  density2q->add_task(task42);

  vector<IndexRange> I38_index = {closed_, virt_, closed_, virt_};
  auto I38 = make_shared<Tensor>(I38_index);
  vector<shared_ptr<Tensor>> tensor43 = {Den1, I38};
  auto task43 = make_shared<Task43>(tensor43, pindex);
  task43->add_dep(task42);
  density2q->add_task(task43);

  vector<shared_ptr<Tensor>> tensor44 = {I38, t2};
  auto task44 = make_shared<Task44>(tensor44, pindex);
  task43->add_dep(task44);
  task44->add_dep(task42);
  density2q->add_task(task44);

  return density2q;
}


