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
  vector<shared_ptr<Tensor>> tensor59 = {Den1};
  auto task59 = make_shared<Task59>(tensor59);
  density2q->add_task(task59);

  vector<IndexRange> I60_index = {closed_, virt_, closed_, virt_};
  auto I60 = make_shared<Tensor>(I60_index);
  vector<shared_ptr<Tensor>> tensor60 = {Den1, I60};
  auto task60 = make_shared<Task60>(tensor60, pindex);
  task60->add_dep(task59);
  density2q->add_task(task60);

  vector<shared_ptr<Tensor>> tensor61 = {I60, t2};
  auto task61 = make_shared<Task61>(tensor61, pindex);
  task60->add_dep(task61);
  task61->add_dep(task59);
  density2q->add_task(task61);

  return density2q;
}


