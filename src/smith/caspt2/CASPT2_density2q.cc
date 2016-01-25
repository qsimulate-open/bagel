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
  auto tensor515 = vector<shared_ptr<Tensor>>{Den1};
  auto task515 = make_shared<Task515>(tensor515, reset);
  density2q->add_task(task515);

  vector<IndexRange> I672_index = {active_, active_, closed_, closed_};
  auto I672 = make_shared<Tensor>(I672_index);
  auto tensor516 = vector<shared_ptr<Tensor>>{Den1, I672};
  auto task516 = make_shared<Task516>(tensor516, pindex);
  task516->add_dep(task515);
  density2q->add_task(task516);

  auto tensor517 = vector<shared_ptr<Tensor>>{I672, t2, Gamma92_()};
  auto task517 = make_shared<Task517>(tensor517, pindex);
  task516->add_dep(task517);
  task517->add_dep(task515);
  density2q->add_task(task517);

  vector<IndexRange> I674_index = {closed_, active_, active_, active_};
  auto I674 = make_shared<Tensor>(I674_index);
  auto tensor518 = vector<shared_ptr<Tensor>>{Den1, I674};
  auto task518 = make_shared<Task518>(tensor518, pindex);
  task518->add_dep(task515);
  density2q->add_task(task518);

  auto tensor519 = vector<shared_ptr<Tensor>>{I674, Gamma6_(), t2};
  auto task519 = make_shared<Task519>(tensor519, pindex);
  task518->add_dep(task519);
  task519->add_dep(task515);
  density2q->add_task(task519);

  vector<IndexRange> I676_index = {closed_, virt_, closed_, active_};
  auto I676 = make_shared<Tensor>(I676_index);
  auto tensor520 = vector<shared_ptr<Tensor>>{Den1, I676};
  auto task520 = make_shared<Task520>(tensor520, pindex);
  task520->add_dep(task515);
  density2q->add_task(task520);

  auto tensor521 = vector<shared_ptr<Tensor>>{I676, Gamma16_(), t2};
  auto task521 = make_shared<Task521>(tensor521, pindex);
  task520->add_dep(task521);
  task521->add_dep(task515);
  density2q->add_task(task521);

  auto tensor522 = vector<shared_ptr<Tensor>>{I676, t2, Gamma16_()};
  auto task522 = make_shared<Task522>(tensor522, pindex);
  task520->add_dep(task522);
  task522->add_dep(task515);
  density2q->add_task(task522);

  vector<IndexRange> I680_index = {virt_, closed_, active_, active_};
  auto I680 = make_shared<Tensor>(I680_index);
  auto tensor523 = vector<shared_ptr<Tensor>>{Den1, I680};
  auto task523 = make_shared<Task523>(tensor523, pindex);
  task523->add_dep(task515);
  density2q->add_task(task523);

  auto tensor524 = vector<shared_ptr<Tensor>>{I680, Gamma32_(), t2};
  auto task524 = make_shared<Task524>(tensor524, pindex);
  task523->add_dep(task524);
  task524->add_dep(task515);
  density2q->add_task(task524);

  auto tensor525 = vector<shared_ptr<Tensor>>{I680, t2, Gamma35_()};
  auto task525 = make_shared<Task525>(tensor525, pindex);
  task523->add_dep(task525);
  task525->add_dep(task515);
  density2q->add_task(task525);

  vector<IndexRange> I684_index = {virt_, closed_, active_, active_};
  auto I684 = make_shared<Tensor>(I684_index);
  auto tensor526 = vector<shared_ptr<Tensor>>{Den1, I684};
  auto task526 = make_shared<Task526>(tensor526, pindex);
  task526->add_dep(task515);
  density2q->add_task(task526);

  vector<IndexRange> I685_index = {active_, virt_, closed_, active_};
  auto I685 = make_shared<Tensor>(I685_index);
  auto tensor527 = vector<shared_ptr<Tensor>>{I684, Gamma35_(), I685};
  auto task527 = make_shared<Task527>(tensor527, pindex);
  task526->add_dep(task527);
  task527->add_dep(task515);
  density2q->add_task(task527);

  auto tensor528 = vector<shared_ptr<Tensor>>{I685, t2};
  auto task528 = make_shared<Task528>(tensor528, pindex);
  task527->add_dep(task528);
  task528->add_dep(task515);
  density2q->add_task(task528);

  vector<IndexRange> I688_index = {virt_, active_, active_, active_};
  auto I688 = make_shared<Tensor>(I688_index);
  auto tensor529 = vector<shared_ptr<Tensor>>{Den1, I688};
  auto task529 = make_shared<Task529>(tensor529, pindex);
  task529->add_dep(task515);
  density2q->add_task(task529);

  auto tensor530 = vector<shared_ptr<Tensor>>{I688, Gamma59_(), t2};
  auto task530 = make_shared<Task530>(tensor530, pindex);
  task529->add_dep(task530);
  task530->add_dep(task515);
  density2q->add_task(task530);

  shared_ptr<Tensor> I690;
  if (diagonal) {
    vector<IndexRange> I690_index = {closed_, virt_, closed_, virt_};
    I690 = make_shared<Tensor>(I690_index);
  }
  shared_ptr<Task531> task531;
  if (diagonal) {
    auto tensor531 = vector<shared_ptr<Tensor>>{Den1, I690};
    task531 = make_shared<Task531>(tensor531, pindex);
    task531->add_dep(task515);
    density2q->add_task(task531);
  }

  shared_ptr<Task532> task532;
  if (diagonal) {
    auto tensor532 = vector<shared_ptr<Tensor>>{I690, t2};
    task532 = make_shared<Task532>(tensor532, pindex);
    task531->add_dep(task532);
    task532->add_dep(task515);
    density2q->add_task(task532);
  }

  vector<IndexRange> I692_index = {virt_, closed_, virt_, active_};
  auto I692 = make_shared<Tensor>(I692_index);
  auto tensor533 = vector<shared_ptr<Tensor>>{Den1, I692};
  auto task533 = make_shared<Task533>(tensor533, pindex);
  task533->add_dep(task515);
  density2q->add_task(task533);

  vector<IndexRange> I693_index = {active_, virt_, closed_, virt_};
  auto I693 = make_shared<Tensor>(I693_index);
  auto tensor534 = vector<shared_ptr<Tensor>>{I692, Gamma38_(), I693};
  auto task534 = make_shared<Task534>(tensor534, pindex);
  task533->add_dep(task534);
  task534->add_dep(task515);
  density2q->add_task(task534);

  auto tensor535 = vector<shared_ptr<Tensor>>{I693, t2};
  auto task535 = make_shared<Task535>(tensor535, pindex);
  task534->add_dep(task535);
  task535->add_dep(task515);
  density2q->add_task(task535);

  vector<IndexRange> I696_index = {virt_, virt_, active_, active_};
  auto I696 = make_shared<Tensor>(I696_index);
  auto tensor536 = vector<shared_ptr<Tensor>>{Den1, I696};
  auto task536 = make_shared<Task536>(tensor536, pindex);
  task536->add_dep(task515);
  density2q->add_task(task536);

  auto tensor537 = vector<shared_ptr<Tensor>>{I696, Gamma60_(), t2};
  auto task537 = make_shared<Task537>(tensor537, pindex);
  task536->add_dep(task537);
  task537->add_dep(task515);
  density2q->add_task(task537);

  return density2q;
}


#endif
