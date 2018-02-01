//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_deci4q.cc
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


#include <src/smith/caspt2/MSCASPT2.h>
#include <src/smith/caspt2/MSCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_deci4q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto deci4q = make_shared<Queue>();
  auto tensor516 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci};
  auto task516 = make_shared<Task516>(tensor516, reset);
  deci4q->add_task(task516);

  vector<IndexRange> I789_index = {active_, active_, active_, active_};
  auto I789 = make_shared<Tensor>(I789_index);
  auto tensor518 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I789};
  auto task518 = make_shared<Task518>(tensor518, pindex);
  task518->add_dep(task516);
  deci4q->add_task(task518);

  auto tensor519 = vector<shared_ptr<Tensor>>{I789, v2_, l2};
  auto task519 = make_shared<Task519>(tensor519, pindex);
  task518->add_dep(task519);
  task519->add_dep(task516);
  deci4q->add_task(task519);

  vector<IndexRange> I792_index = {active_, active_, active_, active_, active_, active_};
  auto I792 = make_shared<Tensor>(I792_index);
  auto tensor520 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I792};
  auto task520 = make_shared<Task520>(tensor520, pindex);
  task520->add_dep(task516);
  deci4q->add_task(task520);

  auto tensor521 = vector<shared_ptr<Tensor>>{I792, v2_, l2};
  auto task521 = make_shared<Task521>(tensor521, pindex);
  task520->add_dep(task521);
  task521->add_dep(task516);
  deci4q->add_task(task521);

  vector<IndexRange> I795_index = {active_, active_, active_, active_, active_, active_};
  auto I795 = make_shared<Tensor>(I795_index);
  auto tensor522 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I795};
  auto task522 = make_shared<Task522>(tensor522, pindex);
  task522->add_dep(task516);
  deci4q->add_task(task522);

  auto tensor523 = vector<shared_ptr<Tensor>>{I795, v2_, l2};
  auto task523 = make_shared<Task523>(tensor523, pindex);
  task522->add_dep(task523);
  task523->add_dep(task516);
  deci4q->add_task(task523);

  vector<IndexRange> I798_index = {active_, active_};
  auto I798 = make_shared<Tensor>(I798_index);
  auto tensor524 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I798};
  auto task524 = make_shared<Task524>(tensor524, pindex);
  task524->add_dep(task516);
  deci4q->add_task(task524);

  auto tensor525 = vector<shared_ptr<Tensor>>{I798, v2_, l2};
  auto task525 = make_shared<Task525>(tensor525, pindex);
  task524->add_dep(task525);
  task525->add_dep(task516);
  deci4q->add_task(task525);

  auto tensor526 = vector<shared_ptr<Tensor>>{I798, v2_, l2};
  auto task526 = make_shared<Task526>(tensor526, pindex);
  task524->add_dep(task526);
  task526->add_dep(task516);
  deci4q->add_task(task526);

  vector<IndexRange> I804_index = {active_, active_, active_, active_};
  auto I804 = make_shared<Tensor>(I804_index);
  auto tensor527 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I804};
  auto task527 = make_shared<Task527>(tensor527, pindex);
  task527->add_dep(task516);
  deci4q->add_task(task527);

  auto tensor528 = vector<shared_ptr<Tensor>>{I804, v2_, l2};
  auto task528 = make_shared<Task528>(tensor528, pindex);
  task527->add_dep(task528);
  task528->add_dep(task516);
  deci4q->add_task(task528);

  auto tensor529 = vector<shared_ptr<Tensor>>{I804, v2_, l2};
  auto task529 = make_shared<Task529>(tensor529, pindex);
  task527->add_dep(task529);
  task529->add_dep(task516);
  deci4q->add_task(task529);

  auto tensor530 = vector<shared_ptr<Tensor>>{I804, v2_, l2};
  auto task530 = make_shared<Task530>(tensor530, pindex);
  task527->add_dep(task530);
  task530->add_dep(task516);
  deci4q->add_task(task530);

  auto tensor531 = vector<shared_ptr<Tensor>>{I804, v2_, l2};
  auto task531 = make_shared<Task531>(tensor531, pindex);
  task527->add_dep(task531);
  task531->add_dep(task516);
  deci4q->add_task(task531);

  auto tensor532 = vector<shared_ptr<Tensor>>{I804, v2_, l2};
  auto task532 = make_shared<Task532>(tensor532, pindex);
  task527->add_dep(task532);
  task532->add_dep(task516);
  deci4q->add_task(task532);

  vector<IndexRange> I807_index = {active_, active_, active_, active_};
  auto I807 = make_shared<Tensor>(I807_index);
  auto tensor533 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I807};
  auto task533 = make_shared<Task533>(tensor533, pindex);
  task533->add_dep(task516);
  deci4q->add_task(task533);

  auto tensor534 = vector<shared_ptr<Tensor>>{I807, v2_, l2};
  auto task534 = make_shared<Task534>(tensor534, pindex);
  task533->add_dep(task534);
  task534->add_dep(task516);
  deci4q->add_task(task534);

  vector<IndexRange> I810_index = {active_, active_, active_, active_};
  auto I810 = make_shared<Tensor>(I810_index);
  auto tensor535 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I810};
  auto task535 = make_shared<Task535>(tensor535, pindex);
  task535->add_dep(task516);
  deci4q->add_task(task535);

  auto tensor536 = vector<shared_ptr<Tensor>>{I810, v2_, l2};
  auto task536 = make_shared<Task536>(tensor536, pindex);
  task535->add_dep(task536);
  task536->add_dep(task516);
  deci4q->add_task(task536);

  vector<IndexRange> I819_index = {active_, active_, active_, active_};
  auto I819 = make_shared<Tensor>(I819_index);
  auto tensor537 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I819};
  auto task537 = make_shared<Task537>(tensor537, pindex);
  task537->add_dep(task516);
  deci4q->add_task(task537);

  auto tensor538 = vector<shared_ptr<Tensor>>{I819, v2_, l2};
  auto task538 = make_shared<Task538>(tensor538, pindex);
  task537->add_dep(task538);
  task538->add_dep(task516);
  deci4q->add_task(task538);

  auto tensor539 = vector<shared_ptr<Tensor>>{I819, h1_, l2};
  auto task539 = make_shared<Task539>(tensor539, pindex);
  task537->add_dep(task539);
  task539->add_dep(task516);
  deci4q->add_task(task539);

  vector<IndexRange> I828_index = {active_, active_, active_, active_, active_, active_};
  auto I828 = make_shared<Tensor>(I828_index);
  auto tensor540 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I828};
  auto task540 = make_shared<Task540>(tensor540, pindex);
  task540->add_dep(task516);
  deci4q->add_task(task540);

  auto tensor541 = vector<shared_ptr<Tensor>>{I828, v2_, l2};
  auto task541 = make_shared<Task541>(tensor541, pindex);
  task540->add_dep(task541);
  task541->add_dep(task516);
  deci4q->add_task(task541);

  vector<IndexRange> I831_index = {active_, active_, active_, active_, active_, active_};
  auto I831 = make_shared<Tensor>(I831_index);
  auto tensor542 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I831};
  auto task542 = make_shared<Task542>(tensor542, pindex);
  task542->add_dep(task516);
  deci4q->add_task(task542);

  auto tensor543 = vector<shared_ptr<Tensor>>{I831, v2_, l2};
  auto task543 = make_shared<Task543>(tensor543, pindex);
  task542->add_dep(task543);
  task543->add_dep(task516);
  deci4q->add_task(task543);

  vector<IndexRange> I834_index = {active_, active_};
  auto I834 = make_shared<Tensor>(I834_index);
  auto tensor544 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I834};
  auto task544 = make_shared<Task544>(tensor544, pindex);
  task544->add_dep(task516);
  deci4q->add_task(task544);

  auto tensor545 = vector<shared_ptr<Tensor>>{I834, v2_, l2};
  auto task545 = make_shared<Task545>(tensor545, pindex);
  task544->add_dep(task545);
  task545->add_dep(task516);
  deci4q->add_task(task545);

  auto tensor546 = vector<shared_ptr<Tensor>>{I834, v2_, l2};
  auto task546 = make_shared<Task546>(tensor546, pindex);
  task544->add_dep(task546);
  task546->add_dep(task516);
  deci4q->add_task(task546);

  auto tensor547 = vector<shared_ptr<Tensor>>{I834, h1_, l2};
  auto task547 = make_shared<Task547>(tensor547, pindex);
  task544->add_dep(task547);
  task547->add_dep(task516);
  deci4q->add_task(task547);

  auto tensor548 = vector<shared_ptr<Tensor>>{I834, h1_, l2};
  auto task548 = make_shared<Task548>(tensor548, pindex);
  task544->add_dep(task548);
  task548->add_dep(task516);
  deci4q->add_task(task548);

  vector<IndexRange> I840_index = {active_, active_, active_, active_};
  auto I840 = make_shared<Tensor>(I840_index);
  auto tensor549 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I840};
  auto task549 = make_shared<Task549>(tensor549, pindex);
  task549->add_dep(task516);
  deci4q->add_task(task549);

  auto tensor550 = vector<shared_ptr<Tensor>>{I840, v2_, l2};
  auto task550 = make_shared<Task550>(tensor550, pindex);
  task549->add_dep(task550);
  task550->add_dep(task516);
  deci4q->add_task(task550);

  auto tensor551 = vector<shared_ptr<Tensor>>{I840, h1_, l2};
  auto task551 = make_shared<Task551>(tensor551, pindex);
  task549->add_dep(task551);
  task551->add_dep(task516);
  deci4q->add_task(task551);

  return deci4q;
}


#endif
