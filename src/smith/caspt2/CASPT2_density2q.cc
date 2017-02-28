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
#include <src/smith/caspt2/CASPT2_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_density2q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density2q = make_shared<Queue>();
  auto tensor516 = vector<shared_ptr<Tensor>>{Den1};
  auto task516 = make_shared<Task516>(tensor516, reset);
  density2q->add_task(task516);

  vector<IndexRange> I658_index = {active_, active_, closed_, closed_};
  auto I658 = make_shared<Tensor>(I658_index);
  auto tensor517 = vector<shared_ptr<Tensor>>{Den1, I658};
  auto task517 = make_shared<Task517>(tensor517, pindex);
  task517->add_dep(task516);
  density2q->add_task(task517);

  auto tensor518 = vector<shared_ptr<Tensor>>{I658, t2, Gamma92_()};
  auto task518 = make_shared<Task518>(tensor518, pindex);
  task517->add_dep(task518);
  task518->add_dep(task516);
  density2q->add_task(task518);

  vector<IndexRange> I660_index = {closed_, active_, active_, active_};
  auto I660 = make_shared<Tensor>(I660_index);
  auto tensor519 = vector<shared_ptr<Tensor>>{Den1, I660};
  auto task519 = make_shared<Task519>(tensor519, pindex);
  task519->add_dep(task516);
  density2q->add_task(task519);

  auto tensor520 = vector<shared_ptr<Tensor>>{I660, Gamma6_(), t2};
  auto task520 = make_shared<Task520>(tensor520, pindex);
  task519->add_dep(task520);
  task520->add_dep(task516);
  density2q->add_task(task520);

  vector<IndexRange> I662_index = {closed_, virt_, closed_, active_};
  auto I662 = make_shared<Tensor>(I662_index);
  auto tensor521 = vector<shared_ptr<Tensor>>{Den1, I662};
  auto task521 = make_shared<Task521>(tensor521, pindex);
  task521->add_dep(task516);
  density2q->add_task(task521);

  auto tensor522 = vector<shared_ptr<Tensor>>{I662, Gamma16_(), t2};
  auto task522 = make_shared<Task522>(tensor522, pindex);
  task521->add_dep(task522);
  task522->add_dep(task516);
  density2q->add_task(task522);

  auto tensor523 = vector<shared_ptr<Tensor>>{I662, t2, Gamma16_()};
  auto task523 = make_shared<Task523>(tensor523, pindex);
  task521->add_dep(task523);
  task523->add_dep(task516);
  density2q->add_task(task523);

  vector<IndexRange> I666_index = {virt_, closed_, active_, active_};
  auto I666 = make_shared<Tensor>(I666_index);
  auto tensor524 = vector<shared_ptr<Tensor>>{Den1, I666};
  auto task524 = make_shared<Task524>(tensor524, pindex);
  task524->add_dep(task516);
  density2q->add_task(task524);

  auto tensor525 = vector<shared_ptr<Tensor>>{I666, Gamma32_(), t2};
  auto task525 = make_shared<Task525>(tensor525, pindex);
  task524->add_dep(task525);
  task525->add_dep(task516);
  density2q->add_task(task525);

  auto tensor526 = vector<shared_ptr<Tensor>>{I666, t2, Gamma35_()};
  auto task526 = make_shared<Task526>(tensor526, pindex);
  task524->add_dep(task526);
  task526->add_dep(task516);
  density2q->add_task(task526);

  vector<IndexRange> I670_index = {virt_, closed_, active_, active_};
  auto I670 = make_shared<Tensor>(I670_index);
  auto tensor527 = vector<shared_ptr<Tensor>>{Den1, I670};
  auto task527 = make_shared<Task527>(tensor527, pindex);
  task527->add_dep(task516);
  density2q->add_task(task527);

  vector<IndexRange> I671_index = {active_, virt_, closed_, active_};
  auto I671 = make_shared<Tensor>(I671_index);
  auto tensor528 = vector<shared_ptr<Tensor>>{I670, Gamma35_(), I671};
  auto task528 = make_shared<Task528>(tensor528, pindex);
  task527->add_dep(task528);
  task528->add_dep(task516);
  density2q->add_task(task528);

  auto tensor529 = vector<shared_ptr<Tensor>>{I671, t2};
  auto task529 = make_shared<Task529>(tensor529, pindex);
  task528->add_dep(task529);
  task529->add_dep(task516);
  density2q->add_task(task529);

  vector<IndexRange> I674_index = {virt_, active_, active_, active_};
  auto I674 = make_shared<Tensor>(I674_index);
  auto tensor530 = vector<shared_ptr<Tensor>>{Den1, I674};
  auto task530 = make_shared<Task530>(tensor530, pindex);
  task530->add_dep(task516);
  density2q->add_task(task530);

  auto tensor531 = vector<shared_ptr<Tensor>>{I674, Gamma59_(), t2};
  auto task531 = make_shared<Task531>(tensor531, pindex);
  task530->add_dep(task531);
  task531->add_dep(task516);
  density2q->add_task(task531);

  shared_ptr<Tensor> I676;
  if (diagonal) {
    vector<IndexRange> I676_index = {closed_, virt_, closed_, virt_};
    I676 = make_shared<Tensor>(I676_index);
  }
  shared_ptr<Task532> task532;
  if (diagonal) {
    auto tensor532 = vector<shared_ptr<Tensor>>{Den1, I676};
    task532 = make_shared<Task532>(tensor532, pindex);
    task532->add_dep(task516);
    density2q->add_task(task532);
  }

  shared_ptr<Task533> task533;
  if (diagonal) {
    auto tensor533 = vector<shared_ptr<Tensor>>{I676, t2};
    task533 = make_shared<Task533>(tensor533, pindex);
    task532->add_dep(task533);
    task533->add_dep(task516);
    density2q->add_task(task533);
  }

  vector<IndexRange> I678_index = {virt_, closed_, virt_, active_};
  auto I678 = make_shared<Tensor>(I678_index);
  auto tensor534 = vector<shared_ptr<Tensor>>{Den1, I678};
  auto task534 = make_shared<Task534>(tensor534, pindex);
  task534->add_dep(task516);
  density2q->add_task(task534);

  vector<IndexRange> I679_index = {active_, virt_, closed_, virt_};
  auto I679 = make_shared<Tensor>(I679_index);
  auto tensor535 = vector<shared_ptr<Tensor>>{I678, Gamma38_(), I679};
  auto task535 = make_shared<Task535>(tensor535, pindex);
  task534->add_dep(task535);
  task535->add_dep(task516);
  density2q->add_task(task535);

  auto tensor536 = vector<shared_ptr<Tensor>>{I679, t2};
  auto task536 = make_shared<Task536>(tensor536, pindex);
  task535->add_dep(task536);
  task536->add_dep(task516);
  density2q->add_task(task536);

  vector<IndexRange> I682_index = {virt_, virt_, active_, active_};
  auto I682 = make_shared<Tensor>(I682_index);
  auto tensor537 = vector<shared_ptr<Tensor>>{Den1, I682};
  auto task537 = make_shared<Task537>(tensor537, pindex);
  task537->add_dep(task516);
  density2q->add_task(task537);

  auto tensor538 = vector<shared_ptr<Tensor>>{I682, Gamma60_(), t2};
  auto task538 = make_shared<Task538>(tensor538, pindex);
  task537->add_dep(task538);
  task538->add_dep(task516);
  density2q->add_task(task538);

  return density2q;
}


#endif
