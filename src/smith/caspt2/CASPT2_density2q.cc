//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_density2qq.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_density2q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density2q = make_shared<Queue>();
  auto tensor528 = vector<shared_ptr<Tensor>>{Den1};
  auto task528 = make_shared<Task528>(tensor528, reset);
  density2q->add_task(task528);

  vector<IndexRange> I658_index = {closed_, closed_, active_, active_};
  auto I658 = make_shared<Tensor>(I658_index);
  auto tensor529 = vector<shared_ptr<Tensor>>{Den1, I658};
  auto task529 = make_shared<Task529>(tensor529, pindex);
  task529->add_dep(task528);
  density2q->add_task(task529);

  auto tensor530 = vector<shared_ptr<Tensor>>{I658, Gamma92_(), t2};
  auto task530 = make_shared<Task530>(tensor530, pindex);
  task529->add_dep(task530);
  task530->add_dep(task528);
  density2q->add_task(task530);

  vector<IndexRange> I660_index = {closed_, active_, active_, active_};
  auto I660 = make_shared<Tensor>(I660_index);
  auto tensor531 = vector<shared_ptr<Tensor>>{Den1, I660};
  auto task531 = make_shared<Task531>(tensor531, pindex);
  task531->add_dep(task528);
  density2q->add_task(task531);

  auto tensor532 = vector<shared_ptr<Tensor>>{I660, Gamma6_(), t2};
  auto task532 = make_shared<Task532>(tensor532, pindex);
  task531->add_dep(task532);
  task532->add_dep(task528);
  density2q->add_task(task532);

  vector<IndexRange> I662_index = {closed_, virt_, closed_, active_};
  auto I662 = make_shared<Tensor>(I662_index);
  auto tensor533 = vector<shared_ptr<Tensor>>{Den1, I662};
  auto task533 = make_shared<Task533>(tensor533, pindex);
  task533->add_dep(task528);
  density2q->add_task(task533);

  vector<IndexRange> I663_index = {closed_, virt_, closed_, active_};
  auto I663 = make_shared<Tensor>(I663_index);
  auto tensor534 = vector<shared_ptr<Tensor>>{I662, Gamma16_(), I663};
  auto task534 = make_shared<Task534>(tensor534, pindex);
  task533->add_dep(task534);
  task534->add_dep(task528);
  density2q->add_task(task534);

  auto tensor535 = vector<shared_ptr<Tensor>>{I663, t2};
  auto task535 = make_shared<Task535>(tensor535, pindex);
  task534->add_dep(task535);
  task535->add_dep(task528);
  density2q->add_task(task535);

  vector<IndexRange> I666_index = {virt_, closed_, active_, active_};
  auto I666 = make_shared<Tensor>(I666_index);
  auto tensor536 = vector<shared_ptr<Tensor>>{Den1, I666};
  auto task536 = make_shared<Task536>(tensor536, pindex);
  task536->add_dep(task528);
  density2q->add_task(task536);

  auto tensor537 = vector<shared_ptr<Tensor>>{I666, Gamma32_(), t2};
  auto task537 = make_shared<Task537>(tensor537, pindex);
  task536->add_dep(task537);
  task537->add_dep(task528);
  density2q->add_task(task537);

  auto tensor538 = vector<shared_ptr<Tensor>>{I666, Gamma35_(), t2};
  auto task538 = make_shared<Task538>(tensor538, pindex);
  task536->add_dep(task538);
  task538->add_dep(task528);
  density2q->add_task(task538);

  vector<IndexRange> I670_index = {virt_, closed_, active_, active_};
  auto I670 = make_shared<Tensor>(I670_index);
  auto tensor539 = vector<shared_ptr<Tensor>>{Den1, I670};
  auto task539 = make_shared<Task539>(tensor539, pindex);
  task539->add_dep(task528);
  density2q->add_task(task539);

  vector<IndexRange> I671_index = {active_, virt_, closed_, active_};
  auto I671 = make_shared<Tensor>(I671_index);
  auto tensor540 = vector<shared_ptr<Tensor>>{I670, Gamma35_(), I671};
  auto task540 = make_shared<Task540>(tensor540, pindex);
  task539->add_dep(task540);
  task540->add_dep(task528);
  density2q->add_task(task540);

  auto tensor541 = vector<shared_ptr<Tensor>>{I671, t2};
  auto task541 = make_shared<Task541>(tensor541, pindex);
  task540->add_dep(task541);
  task541->add_dep(task528);
  density2q->add_task(task541);

  vector<IndexRange> I674_index = {virt_, active_, active_, active_};
  auto I674 = make_shared<Tensor>(I674_index);
  auto tensor542 = vector<shared_ptr<Tensor>>{Den1, I674};
  auto task542 = make_shared<Task542>(tensor542, pindex);
  task542->add_dep(task528);
  density2q->add_task(task542);

  auto tensor543 = vector<shared_ptr<Tensor>>{I674, Gamma59_(), t2};
  auto task543 = make_shared<Task543>(tensor543, pindex);
  task542->add_dep(task543);
  task543->add_dep(task528);
  density2q->add_task(task543);

  shared_ptr<Tensor> I676;
  if (diagonal) {
    vector<IndexRange> I676_index = {closed_, virt_, closed_, virt_};
    I676 = make_shared<Tensor>(I676_index);
  }
  shared_ptr<Task544> task544;
  if (diagonal) {
    auto tensor544 = vector<shared_ptr<Tensor>>{Den1, I676};
    task544 = make_shared<Task544>(tensor544, pindex);
    task544->add_dep(task528);
    density2q->add_task(task544);
  }

  shared_ptr<Task545> task545;
  if (diagonal) {
    auto tensor545 = vector<shared_ptr<Tensor>>{I676, t2};
    task545 = make_shared<Task545>(tensor545, pindex);
    task544->add_dep(task545);
    task545->add_dep(task528);
    density2q->add_task(task545);
  }

  vector<IndexRange> I678_index = {virt_, closed_, virt_, active_};
  auto I678 = make_shared<Tensor>(I678_index);
  auto tensor546 = vector<shared_ptr<Tensor>>{Den1, I678};
  auto task546 = make_shared<Task546>(tensor546, pindex);
  task546->add_dep(task528);
  density2q->add_task(task546);

  vector<IndexRange> I679_index = {active_, virt_, closed_, virt_};
  auto I679 = make_shared<Tensor>(I679_index);
  auto tensor547 = vector<shared_ptr<Tensor>>{I678, Gamma38_(), I679};
  auto task547 = make_shared<Task547>(tensor547, pindex);
  task546->add_dep(task547);
  task547->add_dep(task528);
  density2q->add_task(task547);

  auto tensor548 = vector<shared_ptr<Tensor>>{I679, t2};
  auto task548 = make_shared<Task548>(tensor548, pindex);
  task547->add_dep(task548);
  task548->add_dep(task528);
  density2q->add_task(task548);

  vector<IndexRange> I682_index = {virt_, virt_, active_, active_};
  auto I682 = make_shared<Tensor>(I682_index);
  auto tensor549 = vector<shared_ptr<Tensor>>{Den1, I682};
  auto task549 = make_shared<Task549>(tensor549, pindex);
  task549->add_dep(task528);
  density2q->add_task(task549);

  auto tensor550 = vector<shared_ptr<Tensor>>{I682, Gamma60_(), t2};
  auto task550 = make_shared<Task550>(tensor550, pindex);
  task549->add_dep(task550);
  task550->add_dep(task528);
  density2q->add_task(task550);

  return density2q;
}


#endif
