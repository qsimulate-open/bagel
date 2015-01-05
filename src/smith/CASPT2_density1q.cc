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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_density1q() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density1q = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor741 = {den1};
  auto task741 = make_shared<Task741>(tensor741);
  density1q->add_task(task741);

  vector<IndexRange> I748_index = {closed_, active_};
  auto I748 = make_shared<Tensor>(I748_index);
  vector<shared_ptr<Tensor>> tensor742 = {den1, I748};
  auto task742 = make_shared<Task742>(tensor742, pindex);
  task742->add_dep(task741);
  density1q->add_task(task742);

  vector<IndexRange> I749_index = {active_, active_, closed_, active_};
  auto I749 = make_shared<Tensor>(I749_index);
  vector<shared_ptr<Tensor>> tensor743 = {I748, Gamma12_(), I749};
  auto task743 = make_shared<Task743>(tensor743, pindex);
  task742->add_dep(task743);
  task743->add_dep(task741);
  density1q->add_task(task743);

  vector<shared_ptr<Tensor>> tensor744 = {I749, t2};
  auto task744 = make_shared<Task744>(tensor744, pindex);
  task743->add_dep(task744);
  task744->add_dep(task741);
  density1q->add_task(task744);

  vector<IndexRange> I750_index = {virt_, closed_};
  auto I750 = make_shared<Tensor>(I750_index);
  vector<shared_ptr<Tensor>> tensor745 = {den1, I750};
  auto task745 = make_shared<Task745>(tensor745, pindex);
  task745->add_dep(task741);
  density1q->add_task(task745);

  vector<IndexRange> I751_index = {active_, virt_, closed_, active_};
  auto I751 = make_shared<Tensor>(I751_index);
  vector<shared_ptr<Tensor>> tensor746 = {I750, Gamma38_(), I751};
  auto task746 = make_shared<Task746>(tensor746, pindex);
  task745->add_dep(task746);
  task746->add_dep(task741);
  density1q->add_task(task746);

  vector<shared_ptr<Tensor>> tensor747 = {I751, t2};
  auto task747 = make_shared<Task747>(tensor747, pindex);
  task746->add_dep(task747);
  task747->add_dep(task741);
  density1q->add_task(task747);

  vector<IndexRange> I753_index = {active_, active_};
  auto I753 = make_shared<Tensor>(I753_index);
  vector<shared_ptr<Tensor>> tensor748 = {I750, t2, I753};
  auto task748 = make_shared<Task748>(tensor748, pindex);
  task745->add_dep(task748);
  task748->add_dep(task741);
  density1q->add_task(task748);

  vector<shared_ptr<Tensor>> tensor749 = {I753, Gamma38_()};
  auto task749 = make_shared<Task749>(tensor749, pindex);
  task748->add_dep(task749);
  task749->add_dep(task741);
  density1q->add_task(task749);

  vector<IndexRange> I754_index = {active_, virt_};
  auto I754 = make_shared<Tensor>(I754_index);
  vector<shared_ptr<Tensor>> tensor750 = {den1, I754};
  auto task750 = make_shared<Task750>(tensor750, pindex);
  task750->add_dep(task741);
  density1q->add_task(task750);

  vector<IndexRange> I755_index = {active_, active_, active_, active_};
  auto I755 = make_shared<Tensor>(I755_index);
  vector<shared_ptr<Tensor>> tensor751 = {I754, t2, I755};
  auto task751 = make_shared<Task751>(tensor751, pindex);
  task750->add_dep(task751);
  task751->add_dep(task741);
  density1q->add_task(task751);

  vector<shared_ptr<Tensor>> tensor752 = {I755, Gamma60_()};
  auto task752 = make_shared<Task752>(tensor752, pindex);
  task751->add_dep(task752);
  task752->add_dep(task741);
  density1q->add_task(task752);

  return density1q;
}


