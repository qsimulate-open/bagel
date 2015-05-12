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

shared_ptr<Queue> CASPT2::CASPT2::make_density1q() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density1q = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor730 = {den1};
  auto task730 = make_shared<Task730>(tensor730);
  density1q->add_task(task730);

  vector<IndexRange> I734_index = {closed_, active_};
  auto I734 = make_shared<Tensor>(I734_index);
  vector<shared_ptr<Tensor>> tensor731 = {den1, I734};
  auto task731 = make_shared<Task731>(tensor731, pindex);
  task731->add_dep(task730);
  density1q->add_task(task731);

  vector<IndexRange> I735_index = {active_, active_, closed_, active_};
  auto I735 = make_shared<Tensor>(I735_index);
  vector<shared_ptr<Tensor>> tensor732 = {I734, Gamma12_(), I735};
  auto task732 = make_shared<Task732>(tensor732, pindex);
  task731->add_dep(task732);
  task732->add_dep(task730);
  density1q->add_task(task732);

  vector<shared_ptr<Tensor>> tensor733 = {I735, t2};
  auto task733 = make_shared<Task733>(tensor733, pindex);
  task732->add_dep(task733);
  task733->add_dep(task730);
  density1q->add_task(task733);

  vector<IndexRange> I736_index = {virt_, closed_};
  auto I736 = make_shared<Tensor>(I736_index);
  vector<shared_ptr<Tensor>> tensor734 = {den1, I736};
  auto task734 = make_shared<Task734>(tensor734, pindex);
  task734->add_dep(task730);
  density1q->add_task(task734);

  vector<IndexRange> I737_index = {active_, virt_, closed_, active_};
  auto I737 = make_shared<Tensor>(I737_index);
  vector<shared_ptr<Tensor>> tensor735 = {I736, Gamma38_(), I737};
  auto task735 = make_shared<Task735>(tensor735, pindex);
  task734->add_dep(task735);
  task735->add_dep(task730);
  density1q->add_task(task735);

  vector<shared_ptr<Tensor>> tensor736 = {I737, t2};
  auto task736 = make_shared<Task736>(tensor736, pindex);
  task735->add_dep(task736);
  task736->add_dep(task730);
  density1q->add_task(task736);

  vector<IndexRange> I739_index = {active_, active_};
  auto I739 = make_shared<Tensor>(I739_index);
  vector<shared_ptr<Tensor>> tensor737 = {I736, t2, I739};
  auto task737 = make_shared<Task737>(tensor737, pindex);
  task734->add_dep(task737);
  task737->add_dep(task730);
  density1q->add_task(task737);

  vector<shared_ptr<Tensor>> tensor738 = {I739, Gamma38_()};
  auto task738 = make_shared<Task738>(tensor738, pindex);
  task737->add_dep(task738);
  task738->add_dep(task730);
  density1q->add_task(task738);

  vector<IndexRange> I740_index = {active_, virt_};
  auto I740 = make_shared<Tensor>(I740_index);
  vector<shared_ptr<Tensor>> tensor739 = {den1, I740};
  auto task739 = make_shared<Task739>(tensor739, pindex);
  task739->add_dep(task730);
  density1q->add_task(task739);

  vector<IndexRange> I741_index = {active_, active_, active_, active_};
  auto I741 = make_shared<Tensor>(I741_index);
  vector<shared_ptr<Tensor>> tensor740 = {I740, t2, I741};
  auto task740 = make_shared<Task740>(tensor740, pindex);
  task739->add_dep(task740);
  task740->add_dep(task730);
  density1q->add_task(task740);

  vector<shared_ptr<Tensor>> tensor741 = {I741, Gamma60_()};
  auto task741 = make_shared<Task741>(tensor741, pindex);
  task740->add_dep(task741);
  task741->add_dep(task730);
  density1q->add_task(task741);

  return density1q;
}


#endif
