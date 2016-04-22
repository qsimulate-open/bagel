//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_deci3qq.cc
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

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_deci3q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deci3q = make_shared<Queue>();
  auto tensor538 = vector<shared_ptr<Tensor>>{deci};
  auto task538 = make_shared<Task538>(tensor538, reset);
  deci3q->add_task(task538);

  vector<IndexRange> I788_index = {ci_};
  auto I788 = make_shared<Tensor>(I788_index);
  auto tensor539 = vector<shared_ptr<Tensor>>{deci, I788};
  auto task539 = make_shared<Task539>(tensor539, cindex);
  task539->add_dep(task538);
  deci3q->add_task(task539);

  vector<IndexRange> I789_index = {active_, active_, active_, active_};
  auto I789 = make_shared<Tensor>(I789_index);
  auto tensor540 = vector<shared_ptr<Tensor>>{I788, Gamma111_(), I789};
  auto task540 = make_shared<Task540>(tensor540, cindex);
  task539->add_dep(task540);
  task540->add_dep(task538);
  deci3q->add_task(task540);

  auto tensor541 = vector<shared_ptr<Tensor>>{I789, v2_, l2};
  auto task541 = make_shared<Task541>(tensor541, cindex);
  task540->add_dep(task541);
  task541->add_dep(task538);
  deci3q->add_task(task541);

  vector<IndexRange> I792_index = {active_, active_, active_, active_, active_, active_};
  auto I792 = make_shared<Tensor>(I792_index);
  auto tensor542 = vector<shared_ptr<Tensor>>{I788, Gamma239_(), I792};
  auto task542 = make_shared<Task542>(tensor542, cindex);
  task539->add_dep(task542);
  task542->add_dep(task538);
  deci3q->add_task(task542);

  auto tensor543 = vector<shared_ptr<Tensor>>{I792, v2_, l2};
  auto task543 = make_shared<Task543>(tensor543, cindex);
  task542->add_dep(task543);
  task543->add_dep(task538);
  deci3q->add_task(task543);

  vector<IndexRange> I795_index = {active_, active_, active_, active_, active_, active_};
  auto I795 = make_shared<Tensor>(I795_index);
  auto tensor544 = vector<shared_ptr<Tensor>>{I788, Gamma116_(), I795};
  auto task544 = make_shared<Task544>(tensor544, cindex);
  task539->add_dep(task544);
  task544->add_dep(task538);
  deci3q->add_task(task544);

  auto tensor545 = vector<shared_ptr<Tensor>>{I795, v2_, l2};
  auto task545 = make_shared<Task545>(tensor545, cindex);
  task544->add_dep(task545);
  task545->add_dep(task538);
  deci3q->add_task(task545);

  vector<IndexRange> I798_index = {active_, active_};
  auto I798 = make_shared<Tensor>(I798_index);
  auto tensor546 = vector<shared_ptr<Tensor>>{I788, Gamma126_(), I798};
  auto task546 = make_shared<Task546>(tensor546, cindex);
  task539->add_dep(task546);
  task546->add_dep(task538);
  deci3q->add_task(task546);

  auto tensor547 = vector<shared_ptr<Tensor>>{I798, v2_, l2};
  auto task547 = make_shared<Task547>(tensor547, cindex);
  task546->add_dep(task547);
  task547->add_dep(task538);
  deci3q->add_task(task547);

  auto tensor548 = vector<shared_ptr<Tensor>>{I798, v2_, l2};
  auto task548 = make_shared<Task548>(tensor548, cindex);
  task546->add_dep(task548);
  task548->add_dep(task538);
  deci3q->add_task(task548);

  vector<IndexRange> I804_index = {active_, active_, active_, active_};
  auto I804 = make_shared<Tensor>(I804_index);
  auto tensor549 = vector<shared_ptr<Tensor>>{I788, Gamma145_(), I804};
  auto task549 = make_shared<Task549>(tensor549, cindex);
  task539->add_dep(task549);
  task549->add_dep(task538);
  deci3q->add_task(task549);

  auto tensor550 = vector<shared_ptr<Tensor>>{I804, v2_, l2};
  auto task550 = make_shared<Task550>(tensor550, cindex);
  task549->add_dep(task550);
  task550->add_dep(task538);
  deci3q->add_task(task550);

  auto tensor551 = vector<shared_ptr<Tensor>>{I804, v2_, l2};
  auto task551 = make_shared<Task551>(tensor551, cindex);
  task549->add_dep(task551);
  task551->add_dep(task538);
  deci3q->add_task(task551);

  auto tensor552 = vector<shared_ptr<Tensor>>{I804, v2_, l2};
  auto task552 = make_shared<Task552>(tensor552, cindex);
  task549->add_dep(task552);
  task552->add_dep(task538);
  deci3q->add_task(task552);

  auto tensor553 = vector<shared_ptr<Tensor>>{I804, v2_, l2};
  auto task553 = make_shared<Task553>(tensor553, cindex);
  task549->add_dep(task553);
  task553->add_dep(task538);
  deci3q->add_task(task553);

  auto tensor554 = vector<shared_ptr<Tensor>>{I804, v2_, l2};
  auto task554 = make_shared<Task554>(tensor554, cindex);
  task549->add_dep(task554);
  task554->add_dep(task538);
  deci3q->add_task(task554);

  vector<IndexRange> I807_index = {active_, active_, active_, active_};
  auto I807 = make_shared<Tensor>(I807_index);
  auto tensor555 = vector<shared_ptr<Tensor>>{I788, Gamma132_(), I807};
  auto task555 = make_shared<Task555>(tensor555, cindex);
  task539->add_dep(task555);
  task555->add_dep(task538);
  deci3q->add_task(task555);

  auto tensor556 = vector<shared_ptr<Tensor>>{I807, v2_, l2};
  auto task556 = make_shared<Task556>(tensor556, cindex);
  task555->add_dep(task556);
  task556->add_dep(task538);
  deci3q->add_task(task556);

  vector<IndexRange> I810_index = {active_, active_, active_, active_};
  auto I810 = make_shared<Tensor>(I810_index);
  auto tensor557 = vector<shared_ptr<Tensor>>{I788, Gamma142_(), I810};
  auto task557 = make_shared<Task557>(tensor557, cindex);
  task539->add_dep(task557);
  task557->add_dep(task538);
  deci3q->add_task(task557);

  auto tensor558 = vector<shared_ptr<Tensor>>{I810, v2_, l2};
  auto task558 = make_shared<Task558>(tensor558, cindex);
  task557->add_dep(task558);
  task558->add_dep(task538);
  deci3q->add_task(task558);

  vector<IndexRange> I819_index = {active_, active_, active_, active_};
  auto I819 = make_shared<Tensor>(I819_index);
  auto tensor559 = vector<shared_ptr<Tensor>>{I788, Gamma122_(), I819};
  auto task559 = make_shared<Task559>(tensor559, cindex);
  task539->add_dep(task559);
  task559->add_dep(task538);
  deci3q->add_task(task559);

  auto tensor560 = vector<shared_ptr<Tensor>>{I819, v2_, l2};
  auto task560 = make_shared<Task560>(tensor560, cindex);
  task559->add_dep(task560);
  task560->add_dep(task538);
  deci3q->add_task(task560);

  auto tensor561 = vector<shared_ptr<Tensor>>{I819, h1_, l2};
  auto task561 = make_shared<Task561>(tensor561, cindex);
  task559->add_dep(task561);
  task561->add_dep(task538);
  deci3q->add_task(task561);

  vector<IndexRange> I828_index = {active_, active_, active_, active_, active_, active_};
  auto I828 = make_shared<Tensor>(I828_index);
  auto tensor562 = vector<shared_ptr<Tensor>>{I788, Gamma169_(), I828};
  auto task562 = make_shared<Task562>(tensor562, cindex);
  task539->add_dep(task562);
  task562->add_dep(task538);
  deci3q->add_task(task562);

  auto tensor563 = vector<shared_ptr<Tensor>>{I828, v2_, l2};
  auto task563 = make_shared<Task563>(tensor563, cindex);
  task562->add_dep(task563);
  task563->add_dep(task538);
  deci3q->add_task(task563);

  vector<IndexRange> I831_index = {active_, active_, active_, active_, active_, active_};
  auto I831 = make_shared<Tensor>(I831_index);
  auto tensor564 = vector<shared_ptr<Tensor>>{I788, Gamma161_(), I831};
  auto task564 = make_shared<Task564>(tensor564, cindex);
  task539->add_dep(task564);
  task564->add_dep(task538);
  deci3q->add_task(task564);

  auto tensor565 = vector<shared_ptr<Tensor>>{I831, v2_, l2};
  auto task565 = make_shared<Task565>(tensor565, cindex);
  task564->add_dep(task565);
  task565->add_dep(task538);
  deci3q->add_task(task565);

  vector<IndexRange> I834_index = {active_, active_};
  auto I834 = make_shared<Tensor>(I834_index);
  auto tensor566 = vector<shared_ptr<Tensor>>{I788, Gamma148_(), I834};
  auto task566 = make_shared<Task566>(tensor566, cindex);
  task539->add_dep(task566);
  task566->add_dep(task538);
  deci3q->add_task(task566);

  auto tensor567 = vector<shared_ptr<Tensor>>{I834, v2_, l2};
  auto task567 = make_shared<Task567>(tensor567, cindex);
  task566->add_dep(task567);
  task567->add_dep(task538);
  deci3q->add_task(task567);

  auto tensor568 = vector<shared_ptr<Tensor>>{I834, v2_, l2};
  auto task568 = make_shared<Task568>(tensor568, cindex);
  task566->add_dep(task568);
  task568->add_dep(task538);
  deci3q->add_task(task568);

  auto tensor569 = vector<shared_ptr<Tensor>>{I834, h1_, l2};
  auto task569 = make_shared<Task569>(tensor569, cindex);
  task566->add_dep(task569);
  task569->add_dep(task538);
  deci3q->add_task(task569);

  auto tensor570 = vector<shared_ptr<Tensor>>{I834, h1_, l2};
  auto task570 = make_shared<Task570>(tensor570, cindex);
  task566->add_dep(task570);
  task570->add_dep(task538);
  deci3q->add_task(task570);

  vector<IndexRange> I840_index = {active_, active_, active_, active_};
  auto I840 = make_shared<Tensor>(I840_index);
  auto tensor571 = vector<shared_ptr<Tensor>>{I788, Gamma170_(), I840};
  auto task571 = make_shared<Task571>(tensor571, cindex);
  task539->add_dep(task571);
  task571->add_dep(task538);
  deci3q->add_task(task571);

  auto tensor572 = vector<shared_ptr<Tensor>>{I840, v2_, l2};
  auto task572 = make_shared<Task572>(tensor572, cindex);
  task571->add_dep(task572);
  task572->add_dep(task538);
  deci3q->add_task(task572);

  auto tensor573 = vector<shared_ptr<Tensor>>{I840, h1_, l2};
  auto task573 = make_shared<Task573>(tensor573, cindex);
  task571->add_dep(task573);
  task573->add_dep(task538);
  deci3q->add_task(task573);

  return deci3q;
}


#endif
