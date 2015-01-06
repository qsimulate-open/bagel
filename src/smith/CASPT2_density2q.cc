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
  vector<shared_ptr<Tensor>> tensor113 = {Den1};
  auto task113 = make_shared<Task113>(tensor113);
  density2q->add_task(task113);

  vector<IndexRange> I114_index = {closed_, virt_, closed_, virt_};
  auto I114 = make_shared<Tensor>(I114_index);
  vector<shared_ptr<Tensor>> tensor114 = {Den1, I114};
  auto task114 = make_shared<Task114>(tensor114, pindex);
  task114->add_dep(task113);
  density2q->add_task(task114);

  vector<shared_ptr<Tensor>> tensor115 = {I114, t2};
  auto task115 = make_shared<Task115>(tensor115, pindex);
  task114->add_dep(task115);
  task115->add_dep(task113);
  density2q->add_task(task115);

  vector<IndexRange> I116_index = {virt_, closed_, virt_, active_};
  auto I116 = make_shared<Tensor>(I116_index);
  vector<shared_ptr<Tensor>> tensor116 = {Den1, I116};
  auto task116 = make_shared<Task116>(tensor116, pindex);
  task116->add_dep(task113);
  density2q->add_task(task116);

  vector<IndexRange> I117_index = {active_, virt_, closed_, virt_};
  auto I117 = make_shared<Tensor>(I117_index);
  vector<shared_ptr<Tensor>> tensor117 = {I116, Gamma0_(), I117};
  auto task117 = make_shared<Task117>(tensor117, pindex);
  task116->add_dep(task117);
  task117->add_dep(task113);
  density2q->add_task(task117);

  vector<shared_ptr<Tensor>> tensor118 = {I117, t2};
  auto task118 = make_shared<Task118>(tensor118, pindex);
  task117->add_dep(task118);
  task118->add_dep(task113);
  density2q->add_task(task118);

  return density2q;
}


