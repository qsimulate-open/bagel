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
  auto tensor536 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci};
  auto task536 = make_shared<Task536>(tensor536, reset);
  deci4q->add_task(task536);

  vector<IndexRange> I817_index = {active_, active_, active_, active_};
  auto I817 = make_shared<Tensor>(I817_index);
  auto tensor538 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I817};
  auto task538 = make_shared<Task538>(tensor538, pindex);
  task538->add_dep(task536);
  deci4q->add_task(task538);

  auto tensor539 = vector<shared_ptr<Tensor>>{I817, v2_, l2};
  auto task539 = make_shared<Task539>(tensor539, pindex);
  task538->add_dep(task539);
  task539->add_dep(task536);
  deci4q->add_task(task539);

  vector<IndexRange> I820_index = {active_, active_, active_, active_, active_, active_};
  auto I820 = make_shared<Tensor>(I820_index);
  auto tensor540 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I820};
  auto task540 = make_shared<Task540>(tensor540, pindex);
  task540->add_dep(task536);
  deci4q->add_task(task540);

  auto tensor541 = vector<shared_ptr<Tensor>>{I820, v2_, l2};
  auto task541 = make_shared<Task541>(tensor541, pindex);
  task540->add_dep(task541);
  task541->add_dep(task536);
  deci4q->add_task(task541);

  vector<IndexRange> I823_index = {active_, active_, active_, active_, active_, active_};
  auto I823 = make_shared<Tensor>(I823_index);
  auto tensor542 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I823};
  auto task542 = make_shared<Task542>(tensor542, pindex);
  task542->add_dep(task536);
  deci4q->add_task(task542);

  auto tensor543 = vector<shared_ptr<Tensor>>{I823, v2_, l2};
  auto task543 = make_shared<Task543>(tensor543, pindex);
  task542->add_dep(task543);
  task543->add_dep(task536);
  deci4q->add_task(task543);

  vector<IndexRange> I826_index = {active_, active_};
  auto I826 = make_shared<Tensor>(I826_index);
  auto tensor544 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I826};
  auto task544 = make_shared<Task544>(tensor544, pindex);
  task544->add_dep(task536);
  deci4q->add_task(task544);

  auto tensor545 = vector<shared_ptr<Tensor>>{I826, v2_, l2};
  auto task545 = make_shared<Task545>(tensor545, pindex);
  task544->add_dep(task545);
  task545->add_dep(task536);
  deci4q->add_task(task545);

  auto tensor546 = vector<shared_ptr<Tensor>>{I826, v2_, l2};
  auto task546 = make_shared<Task546>(tensor546, pindex);
  task544->add_dep(task546);
  task546->add_dep(task536);
  deci4q->add_task(task546);

  vector<IndexRange> I832_index = {active_, active_, active_, active_};
  auto I832 = make_shared<Tensor>(I832_index);
  auto tensor547 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I832};
  auto task547 = make_shared<Task547>(tensor547, pindex);
  task547->add_dep(task536);
  deci4q->add_task(task547);

  auto tensor548 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task548 = make_shared<Task548>(tensor548, pindex);
  task547->add_dep(task548);
  task548->add_dep(task536);
  deci4q->add_task(task548);

  auto tensor549 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task549 = make_shared<Task549>(tensor549, pindex);
  task547->add_dep(task549);
  task549->add_dep(task536);
  deci4q->add_task(task549);

  auto tensor550 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task550 = make_shared<Task550>(tensor550, pindex);
  task547->add_dep(task550);
  task550->add_dep(task536);
  deci4q->add_task(task550);

  auto tensor551 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task551 = make_shared<Task551>(tensor551, pindex);
  task547->add_dep(task551);
  task551->add_dep(task536);
  deci4q->add_task(task551);

  auto tensor552 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task552 = make_shared<Task552>(tensor552, pindex);
  task547->add_dep(task552);
  task552->add_dep(task536);
  deci4q->add_task(task552);

  vector<IndexRange> I835_index = {active_, active_, active_, active_};
  auto I835 = make_shared<Tensor>(I835_index);
  auto tensor553 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I835};
  auto task553 = make_shared<Task553>(tensor553, pindex);
  task553->add_dep(task536);
  deci4q->add_task(task553);

  auto tensor554 = vector<shared_ptr<Tensor>>{I835, v2_, l2};
  auto task554 = make_shared<Task554>(tensor554, pindex);
  task553->add_dep(task554);
  task554->add_dep(task536);
  deci4q->add_task(task554);

  vector<IndexRange> I838_index = {active_, active_, active_, active_};
  auto I838 = make_shared<Tensor>(I838_index);
  auto tensor555 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I838};
  auto task555 = make_shared<Task555>(tensor555, pindex);
  task555->add_dep(task536);
  deci4q->add_task(task555);

  auto tensor556 = vector<shared_ptr<Tensor>>{I838, v2_, l2};
  auto task556 = make_shared<Task556>(tensor556, pindex);
  task555->add_dep(task556);
  task556->add_dep(task536);
  deci4q->add_task(task556);

  vector<IndexRange> I847_index = {active_, active_, active_, active_};
  auto I847 = make_shared<Tensor>(I847_index);
  auto tensor557 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I847};
  auto task557 = make_shared<Task557>(tensor557, pindex);
  task557->add_dep(task536);
  deci4q->add_task(task557);

  auto tensor558 = vector<shared_ptr<Tensor>>{I847, v2_, l2};
  auto task558 = make_shared<Task558>(tensor558, pindex);
  task557->add_dep(task558);
  task558->add_dep(task536);
  deci4q->add_task(task558);

  auto tensor559 = vector<shared_ptr<Tensor>>{I847, h1_, l2};
  auto task559 = make_shared<Task559>(tensor559, pindex);
  task557->add_dep(task559);
  task559->add_dep(task536);
  deci4q->add_task(task559);

  vector<IndexRange> I856_index = {active_, active_, active_, active_, active_, active_};
  auto I856 = make_shared<Tensor>(I856_index);
  auto tensor560 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I856};
  auto task560 = make_shared<Task560>(tensor560, pindex);
  task560->add_dep(task536);
  deci4q->add_task(task560);

  auto tensor561 = vector<shared_ptr<Tensor>>{I856, v2_, l2};
  auto task561 = make_shared<Task561>(tensor561, pindex);
  task560->add_dep(task561);
  task561->add_dep(task536);
  deci4q->add_task(task561);

  vector<IndexRange> I859_index = {active_, active_, active_, active_, active_, active_};
  auto I859 = make_shared<Tensor>(I859_index);
  auto tensor562 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I859};
  auto task562 = make_shared<Task562>(tensor562, pindex);
  task562->add_dep(task536);
  deci4q->add_task(task562);

  auto tensor563 = vector<shared_ptr<Tensor>>{I859, v2_, l2};
  auto task563 = make_shared<Task563>(tensor563, pindex);
  task562->add_dep(task563);
  task563->add_dep(task536);
  deci4q->add_task(task563);

  shared_ptr<Tensor> I862;
  if (diagonal) {
    vector<IndexRange> I862_index;
    I862 = make_shared<Tensor>(I862_index);
  }
  shared_ptr<Task564> task564;
  if (diagonal) {
    auto tensor564 = vector<shared_ptr<Tensor>>{den0ci, I862};
    task564 = make_shared<Task564>(tensor564, pindex);
    task564->add_dep(task536);
    deci4q->add_task(task564);
  }

  shared_ptr<Task565> task565;
  if (diagonal) {
    auto tensor565 = vector<shared_ptr<Tensor>>{I862, v2_, l2};
    task565 = make_shared<Task565>(tensor565, pindex);
    task564->add_dep(task565);
    task565->add_dep(task536);
    deci4q->add_task(task565);
  }

  shared_ptr<Tensor> I865;
  if (diagonal) {
    vector<IndexRange> I865_index;
    I865 = make_shared<Tensor>(I865_index);
  }
  shared_ptr<Task566> task566;
  if (diagonal) {
    auto tensor566 = vector<shared_ptr<Tensor>>{den0ci, I865};
    task566 = make_shared<Task566>(tensor566, pindex);
    task566->add_dep(task536);
    deci4q->add_task(task566);
  }

  shared_ptr<Task567> task567;
  if (diagonal) {
    auto tensor567 = vector<shared_ptr<Tensor>>{I865, v2_, l2};
    task567 = make_shared<Task567>(tensor567, pindex);
    task566->add_dep(task567);
    task567->add_dep(task536);
    deci4q->add_task(task567);
  }

  vector<IndexRange> I868_index = {active_, active_};
  auto I868 = make_shared<Tensor>(I868_index);
  auto tensor568 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I868};
  auto task568 = make_shared<Task568>(tensor568, pindex);
  task568->add_dep(task536);
  deci4q->add_task(task568);

  auto tensor569 = vector<shared_ptr<Tensor>>{I868, v2_, l2};
  auto task569 = make_shared<Task569>(tensor569, pindex);
  task568->add_dep(task569);
  task569->add_dep(task536);
  deci4q->add_task(task569);

  auto tensor570 = vector<shared_ptr<Tensor>>{I868, v2_, l2};
  auto task570 = make_shared<Task570>(tensor570, pindex);
  task568->add_dep(task570);
  task570->add_dep(task536);
  deci4q->add_task(task570);

  auto tensor571 = vector<shared_ptr<Tensor>>{I868, h1_, l2};
  auto task571 = make_shared<Task571>(tensor571, pindex);
  task568->add_dep(task571);
  task571->add_dep(task536);
  deci4q->add_task(task571);

  auto tensor572 = vector<shared_ptr<Tensor>>{I868, h1_, l2};
  auto task572 = make_shared<Task572>(tensor572, pindex);
  task568->add_dep(task572);
  task572->add_dep(task536);
  deci4q->add_task(task572);

  vector<IndexRange> I874_index = {active_, active_, active_, active_};
  auto I874 = make_shared<Tensor>(I874_index);
  auto tensor573 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I874};
  auto task573 = make_shared<Task573>(tensor573, pindex);
  task573->add_dep(task536);
  deci4q->add_task(task573);

  auto tensor574 = vector<shared_ptr<Tensor>>{I874, v2_, l2};
  auto task574 = make_shared<Task574>(tensor574, pindex);
  task573->add_dep(task574);
  task574->add_dep(task536);
  deci4q->add_task(task574);

  auto tensor575 = vector<shared_ptr<Tensor>>{I874, h1_, l2};
  auto task575 = make_shared<Task575>(tensor575, pindex);
  task573->add_dep(task575);
  task575->add_dep(task536);
  deci4q->add_task(task575);

  return deci4q;
}


#endif
