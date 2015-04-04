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

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_density2q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density2q = make_shared<Queue>();
  auto tensor530 = vector<shared_ptr<Tensor>>{Den1};
  auto task530 = make_shared<Task530>(tensor530, reset);
  density2q->add_task(task530);

  vector<IndexRange> I742_index = {closed_, closed_, active_, active_};
  auto I742 = make_shared<Tensor>(I742_index);
  auto tensor531 = vector<shared_ptr<Tensor>>{Den1, I742};
  auto task531 = make_shared<Task531>(tensor531, pindex);
  task531->add_dep(task530);
  density2q->add_task(task531);

  auto tensor532 = vector<shared_ptr<Tensor>>{I742, Gamma92_(), t2};
  auto task532 = make_shared<Task532>(tensor532, pindex);
  task531->add_dep(task532);
  task532->add_dep(task530);
  density2q->add_task(task532);

  vector<IndexRange> I744_index = {closed_, active_, active_, active_};
  auto I744 = make_shared<Tensor>(I744_index);
  auto tensor533 = vector<shared_ptr<Tensor>>{Den1, I744};
  auto task533 = make_shared<Task533>(tensor533, pindex);
  task533->add_dep(task530);
  density2q->add_task(task533);

  auto tensor534 = vector<shared_ptr<Tensor>>{I744, Gamma6_(), t2};
  auto task534 = make_shared<Task534>(tensor534, pindex);
  task533->add_dep(task534);
  task534->add_dep(task530);
  density2q->add_task(task534);

  vector<IndexRange> I746_index = {closed_, virt_, closed_, active_};
  auto I746 = make_shared<Tensor>(I746_index);
  auto tensor535 = vector<shared_ptr<Tensor>>{Den1, I746};
  auto task535 = make_shared<Task535>(tensor535, pindex);
  task535->add_dep(task530);
  density2q->add_task(task535);

  vector<IndexRange> I747_index = {closed_, virt_, closed_, active_};
  auto I747 = make_shared<Tensor>(I747_index);
  auto tensor536 = vector<shared_ptr<Tensor>>{I746, Gamma16_(), I747};
  auto task536 = make_shared<Task536>(tensor536, pindex);
  task535->add_dep(task536);
  task536->add_dep(task530);
  density2q->add_task(task536);

  auto tensor537 = vector<shared_ptr<Tensor>>{I747, t2};
  auto task537 = make_shared<Task537>(tensor537, pindex);
  task536->add_dep(task537);
  task537->add_dep(task530);
  density2q->add_task(task537);

  vector<IndexRange> I750_index = {virt_, closed_, active_, active_};
  auto I750 = make_shared<Tensor>(I750_index);
  auto tensor538 = vector<shared_ptr<Tensor>>{Den1, I750};
  auto task538 = make_shared<Task538>(tensor538, pindex);
  task538->add_dep(task530);
  density2q->add_task(task538);

  auto tensor539 = vector<shared_ptr<Tensor>>{I750, Gamma32_(), t2};
  auto task539 = make_shared<Task539>(tensor539, pindex);
  task538->add_dep(task539);
  task539->add_dep(task530);
  density2q->add_task(task539);

  auto tensor540 = vector<shared_ptr<Tensor>>{I750, Gamma35_(), t2};
  auto task540 = make_shared<Task540>(tensor540, pindex);
  task538->add_dep(task540);
  task540->add_dep(task530);
  density2q->add_task(task540);

  vector<IndexRange> I754_index = {active_, active_, virt_, closed_};
  auto I754 = make_shared<Tensor>(I754_index);
  auto tensor541 = vector<shared_ptr<Tensor>>{Den1, I754};
  auto task541 = make_shared<Task541>(tensor541, pindex);
  task541->add_dep(task530);
  density2q->add_task(task541);

  auto tensor542 = vector<shared_ptr<Tensor>>{I754, t2, Gamma35_()};
  auto task542 = make_shared<Task542>(tensor542, pindex);
  task541->add_dep(task542);
  task542->add_dep(task530);
  density2q->add_task(task542);

  auto tensor543 = vector<shared_ptr<Tensor>>{I754, Gamma35_(), t2};
  auto task543 = make_shared<Task543>(tensor543, pindex);
  task541->add_dep(task543);
  task543->add_dep(task530);
  density2q->add_task(task543);

  vector<IndexRange> I758_index = {active_, active_, active_, virt_};
  auto I758 = make_shared<Tensor>(I758_index);
  auto tensor544 = vector<shared_ptr<Tensor>>{Den1, I758};
  auto task544 = make_shared<Task544>(tensor544, pindex);
  task544->add_dep(task530);
  density2q->add_task(task544);

  auto tensor545 = vector<shared_ptr<Tensor>>{I758, t2, Gamma59_()};
  auto task545 = make_shared<Task545>(tensor545, pindex);
  task544->add_dep(task545);
  task545->add_dep(task530);
  density2q->add_task(task545);

  shared_ptr<Tensor> I760;
  if (diagonal) {
    vector<IndexRange> I760_index = {closed_, virt_, closed_, virt_};
    I760 = make_shared<Tensor>(I760_index);
  }
  shared_ptr<Task546> task546;
  if (diagonal) {
    auto tensor546 = vector<shared_ptr<Tensor>>{Den1, I760};
    task546 = make_shared<Task546>(tensor546, pindex);
    task546->add_dep(task530);
    density2q->add_task(task546);
  }

  shared_ptr<Task547> task547;
  if (diagonal) {
    auto tensor547 = vector<shared_ptr<Tensor>>{I760, t2};
    task547 = make_shared<Task547>(tensor547, pindex);
    task546->add_dep(task547);
    task547->add_dep(task530);
    density2q->add_task(task547);
  }

  vector<IndexRange> I762_index = {virt_, closed_, virt_, active_};
  auto I762 = make_shared<Tensor>(I762_index);
  auto tensor548 = vector<shared_ptr<Tensor>>{Den1, I762};
  auto task548 = make_shared<Task548>(tensor548, pindex);
  task548->add_dep(task530);
  density2q->add_task(task548);

  vector<IndexRange> I763_index = {active_, virt_, closed_, virt_};
  auto I763 = make_shared<Tensor>(I763_index);
  auto tensor549 = vector<shared_ptr<Tensor>>{I762, Gamma38_(), I763};
  auto task549 = make_shared<Task549>(tensor549, pindex);
  task548->add_dep(task549);
  task549->add_dep(task530);
  density2q->add_task(task549);

  auto tensor550 = vector<shared_ptr<Tensor>>{I763, t2};
  auto task550 = make_shared<Task550>(tensor550, pindex);
  task549->add_dep(task550);
  task550->add_dep(task530);
  density2q->add_task(task550);

  vector<IndexRange> I766_index = {virt_, virt_, active_, active_};
  auto I766 = make_shared<Tensor>(I766_index);
  auto tensor551 = vector<shared_ptr<Tensor>>{Den1, I766};
  auto task551 = make_shared<Task551>(tensor551, pindex);
  task551->add_dep(task530);
  density2q->add_task(task551);

  auto tensor552 = vector<shared_ptr<Tensor>>{I766, Gamma60_(), t2};
  auto task552 = make_shared<Task552>(tensor552, pindex);
  task551->add_dep(task552);
  task552->add_dep(task530);
  density2q->add_task(task552);

  return density2q;
}


#endif
