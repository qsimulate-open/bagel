//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_energyqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_energyq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto energyq = make_shared<Queue>();
  vector<IndexRange> I44_index;
  auto I44 = make_shared<Tensor>(I44_index);
  vector<IndexRange> I45_index = {closed_, virt_, closed_, virt_};
  auto I45 = make_shared<Tensor>(I45_index);
  vector<shared_ptr<Tensor>> tensor39 = {I44, t2, I45};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  energyq->add_task(task39);

  vector<shared_ptr<Tensor>> tensor40 = {I45, v2_};
  auto task40 = make_shared<Task40>(tensor40, pindex);
  task39->add_dep(task40);
  energyq->add_task(task40);

  vector<IndexRange> I47_index = {virt_, closed_, virt_, closed_};
  auto I47 = make_shared<Tensor>(I47_index);
  vector<shared_ptr<Tensor>> tensor41 = {I44, v2_, I47};
  auto task41 = make_shared<Task41>(tensor41, pindex);
  task39->add_dep(task41);
  energyq->add_task(task41);

  vector<shared_ptr<Tensor>> tensor42 = {I47, t2};
  auto task42 = make_shared<Task42>(tensor42, pindex);
  task41->add_dep(task42);
  energyq->add_task(task42);

  vector<IndexRange> I49_index = {active_, active_};
  auto I49 = make_shared<Tensor>(I49_index);
  vector<shared_ptr<Tensor>> tensor43 = {I44, Gamma0_(), I49};
  auto task43 = make_shared<Task43>(tensor43, pindex);
  task39->add_dep(task43);
  energyq->add_task(task43);

  vector<IndexRange> I50_index = {virt_, closed_, virt_, active_};
  auto I50 = make_shared<Tensor>(I50_index);
  vector<shared_ptr<Tensor>> tensor44 = {I49, v2_, I50};
  auto task44 = make_shared<Task44>(tensor44, pindex);
  task43->add_dep(task44);
  energyq->add_task(task44);

  vector<shared_ptr<Tensor>> tensor45 = {I50, t2};
  auto task45 = make_shared<Task45>(tensor45, pindex);
  task44->add_dep(task45);
  energyq->add_task(task45);

  vector<IndexRange> I53_index = {active_, virt_, closed_, virt_};
  auto I53 = make_shared<Tensor>(I53_index);
  vector<shared_ptr<Tensor>> tensor46 = {I49, t2, I53};
  auto task46 = make_shared<Task46>(tensor46, pindex);
  task43->add_dep(task46);
  energyq->add_task(task46);

  vector<shared_ptr<Tensor>> tensor47 = {I53, v2_};
  auto task47 = make_shared<Task47>(tensor47, pindex);
  task46->add_dep(task47);
  energyq->add_task(task47);

  return energyq;
}


