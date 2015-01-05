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
  vector<IndexRange> I16_index;
  auto I16 = make_shared<Tensor>(I16_index);
  vector<IndexRange> I17_index = {virt_, closed_, virt_, closed_};
  auto I17 = make_shared<Tensor>(I17_index);
  vector<shared_ptr<Tensor>> tensor19 = {I16, v2_, I17};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  energyq->add_task(task19);

  vector<shared_ptr<Tensor>> tensor20 = {I17, t2};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  task19->add_dep(task20);
  energyq->add_task(task20);

  vector<IndexRange> I19_index = {closed_, virt_, closed_, virt_};
  auto I19 = make_shared<Tensor>(I19_index);
  vector<shared_ptr<Tensor>> tensor21 = {I16, t2, I19};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  task19->add_dep(task21);
  energyq->add_task(task21);

  vector<shared_ptr<Tensor>> tensor22 = {I19, v2_};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  task21->add_dep(task22);
  energyq->add_task(task22);

  return energyq;
}


