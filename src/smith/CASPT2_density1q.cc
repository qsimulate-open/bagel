//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_density1qq.cc
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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_density1q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density1q = make_shared<Queue>();
  auto tensor522 = vector<shared_ptr<Tensor>>{den1};
  auto task522 = make_shared<Task522>(tensor522, reset);
  density1q->add_task(task522);

  vector<IndexRange> I734_index = {active_, closed_};
  auto I734 = make_shared<Tensor>(I734_index);
  auto tensor523 = vector<shared_ptr<Tensor>>{den1, I734};
  auto task523 = make_shared<Task523>(tensor523, pindex);
  task523->add_dep(task522);
  density1q->add_task(task523);

  auto tensor524 = vector<shared_ptr<Tensor>>{I734, t2, Gamma12_()};
  auto task524 = make_shared<Task524>(tensor524, pindex);
  task523->add_dep(task524);
  task524->add_dep(task522);
  density1q->add_task(task524);

  vector<IndexRange> I736_index = {virt_, closed_};
  auto I736 = make_shared<Tensor>(I736_index);
  auto tensor525 = vector<shared_ptr<Tensor>>{den1, I736};
  auto task525 = make_shared<Task525>(tensor525, pindex);
  task525->add_dep(task522);
  density1q->add_task(task525);

  auto tensor526 = vector<shared_ptr<Tensor>>{I736, t2, Gamma38_()};
  auto task526 = make_shared<Task526>(tensor526, pindex);
  task525->add_dep(task526);
  task526->add_dep(task522);
  density1q->add_task(task526);

  auto tensor527 = vector<shared_ptr<Tensor>>{I736, Gamma38_(), t2};
  auto task527 = make_shared<Task527>(tensor527, pindex);
  task525->add_dep(task527);
  task527->add_dep(task522);
  density1q->add_task(task527);

  vector<IndexRange> I740_index = {active_, virt_};
  auto I740 = make_shared<Tensor>(I740_index);
  auto tensor528 = vector<shared_ptr<Tensor>>{den1, I740};
  auto task528 = make_shared<Task528>(tensor528, pindex);
  task528->add_dep(task522);
  density1q->add_task(task528);

  auto tensor529 = vector<shared_ptr<Tensor>>{I740, t2, Gamma60_()};
  auto task529 = make_shared<Task529>(tensor529, pindex);
  task528->add_dep(task529);
  task529->add_dep(task522);
  density1q->add_task(task529);

  return density1q;
}


#endif
