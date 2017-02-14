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

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deci4q = make_shared<Queue>();
  auto tensor568 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci};
  auto task568 = make_shared<Task568>(tensor568, reset);
  deci4q->add_task(task568);

  vector<IndexRange> I817_index = {active_, active_, active_, active_};
  auto I817 = make_shared<Tensor>(I817_index);
  auto tensor570 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I817};
  auto task570 = make_shared<Task570>(tensor570, cindex);
  task570->add_dep(task568);
  deci4q->add_task(task570);

  auto tensor571 = vector<shared_ptr<Tensor>>{I817, v2_, l2};
  auto task571 = make_shared<Task571>(tensor571, cindex);
  task570->add_dep(task571);
  task571->add_dep(task568);
  deci4q->add_task(task571);

  vector<IndexRange> I820_index = {active_, active_, active_, active_, active_, active_};
  auto I820 = make_shared<Tensor>(I820_index);
  auto tensor572 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I820};
  auto task572 = make_shared<Task572>(tensor572, cindex);
  task572->add_dep(task568);
  deci4q->add_task(task572);

  auto tensor573 = vector<shared_ptr<Tensor>>{I820, v2_, l2};
  auto task573 = make_shared<Task573>(tensor573, cindex);
  task572->add_dep(task573);
  task573->add_dep(task568);
  deci4q->add_task(task573);

  vector<IndexRange> I823_index = {active_, active_, active_, active_, active_, active_};
  auto I823 = make_shared<Tensor>(I823_index);
  auto tensor574 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I823};
  auto task574 = make_shared<Task574>(tensor574, cindex);
  task574->add_dep(task568);
  deci4q->add_task(task574);

  auto tensor575 = vector<shared_ptr<Tensor>>{I823, v2_, l2};
  auto task575 = make_shared<Task575>(tensor575, cindex);
  task574->add_dep(task575);
  task575->add_dep(task568);
  deci4q->add_task(task575);

  vector<IndexRange> I826_index = {active_, active_};
  auto I826 = make_shared<Tensor>(I826_index);
  auto tensor576 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I826};
  auto task576 = make_shared<Task576>(tensor576, cindex);
  task576->add_dep(task568);
  deci4q->add_task(task576);

  auto tensor577 = vector<shared_ptr<Tensor>>{I826, v2_, l2};
  auto task577 = make_shared<Task577>(tensor577, cindex);
  task576->add_dep(task577);
  task577->add_dep(task568);
  deci4q->add_task(task577);

  auto tensor578 = vector<shared_ptr<Tensor>>{I826, v2_, l2};
  auto task578 = make_shared<Task578>(tensor578, cindex);
  task576->add_dep(task578);
  task578->add_dep(task568);
  deci4q->add_task(task578);

  vector<IndexRange> I832_index = {active_, active_, active_, active_};
  auto I832 = make_shared<Tensor>(I832_index);
  auto tensor579 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I832};
  auto task579 = make_shared<Task579>(tensor579, cindex);
  task579->add_dep(task568);
  deci4q->add_task(task579);

  auto tensor580 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task580 = make_shared<Task580>(tensor580, cindex);
  task579->add_dep(task580);
  task580->add_dep(task568);
  deci4q->add_task(task580);

  auto tensor581 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task581 = make_shared<Task581>(tensor581, cindex);
  task579->add_dep(task581);
  task581->add_dep(task568);
  deci4q->add_task(task581);

  auto tensor582 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task582 = make_shared<Task582>(tensor582, cindex);
  task579->add_dep(task582);
  task582->add_dep(task568);
  deci4q->add_task(task582);

  auto tensor583 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task583 = make_shared<Task583>(tensor583, cindex);
  task579->add_dep(task583);
  task583->add_dep(task568);
  deci4q->add_task(task583);

  auto tensor584 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task584 = make_shared<Task584>(tensor584, cindex);
  task579->add_dep(task584);
  task584->add_dep(task568);
  deci4q->add_task(task584);

  vector<IndexRange> I835_index = {active_, active_, active_, active_};
  auto I835 = make_shared<Tensor>(I835_index);
  auto tensor585 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I835};
  auto task585 = make_shared<Task585>(tensor585, cindex);
  task585->add_dep(task568);
  deci4q->add_task(task585);

  auto tensor586 = vector<shared_ptr<Tensor>>{I835, v2_, l2};
  auto task586 = make_shared<Task586>(tensor586, cindex);
  task585->add_dep(task586);
  task586->add_dep(task568);
  deci4q->add_task(task586);

  vector<IndexRange> I838_index = {active_, active_, active_, active_};
  auto I838 = make_shared<Tensor>(I838_index);
  auto tensor587 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I838};
  auto task587 = make_shared<Task587>(tensor587, cindex);
  task587->add_dep(task568);
  deci4q->add_task(task587);

  auto tensor588 = vector<shared_ptr<Tensor>>{I838, v2_, l2};
  auto task588 = make_shared<Task588>(tensor588, cindex);
  task587->add_dep(task588);
  task588->add_dep(task568);
  deci4q->add_task(task588);

  vector<IndexRange> I847_index = {active_, active_, active_, active_};
  auto I847 = make_shared<Tensor>(I847_index);
  auto tensor589 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I847};
  auto task589 = make_shared<Task589>(tensor589, cindex);
  task589->add_dep(task568);
  deci4q->add_task(task589);

  auto tensor590 = vector<shared_ptr<Tensor>>{I847, v2_, l2};
  auto task590 = make_shared<Task590>(tensor590, cindex);
  task589->add_dep(task590);
  task590->add_dep(task568);
  deci4q->add_task(task590);

  auto tensor591 = vector<shared_ptr<Tensor>>{I847, h1_, l2};
  auto task591 = make_shared<Task591>(tensor591, cindex);
  task589->add_dep(task591);
  task591->add_dep(task568);
  deci4q->add_task(task591);

  vector<IndexRange> I856_index = {active_, active_, active_, active_, active_, active_};
  auto I856 = make_shared<Tensor>(I856_index);
  auto tensor592 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I856};
  auto task592 = make_shared<Task592>(tensor592, cindex);
  task592->add_dep(task568);
  deci4q->add_task(task592);

  auto tensor593 = vector<shared_ptr<Tensor>>{I856, v2_, l2};
  auto task593 = make_shared<Task593>(tensor593, cindex);
  task592->add_dep(task593);
  task593->add_dep(task568);
  deci4q->add_task(task593);

  vector<IndexRange> I859_index = {active_, active_, active_, active_, active_, active_};
  auto I859 = make_shared<Tensor>(I859_index);
  auto tensor594 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I859};
  auto task594 = make_shared<Task594>(tensor594, cindex);
  task594->add_dep(task568);
  deci4q->add_task(task594);

  auto tensor595 = vector<shared_ptr<Tensor>>{I859, v2_, l2};
  auto task595 = make_shared<Task595>(tensor595, cindex);
  task594->add_dep(task595);
  task595->add_dep(task568);
  deci4q->add_task(task595);

  shared_ptr<Tensor> I862;
  if (diagonal) {
    vector<IndexRange> I862_index;
    I862 = make_shared<Tensor>(I862_index);
  }
  shared_ptr<Task596> task596;
  if (diagonal) {
    auto tensor596 = vector<shared_ptr<Tensor>>{den0ci, I862};
    task596 = make_shared<Task596>(tensor596, cindex);
    task596->add_dep(task568);
    deci4q->add_task(task596);
  }

  shared_ptr<Task597> task597;
  if (diagonal) {
    auto tensor597 = vector<shared_ptr<Tensor>>{I862, v2_, l2};
    task597 = make_shared<Task597>(tensor597, cindex);
    task596->add_dep(task597);
    task597->add_dep(task568);
    deci4q->add_task(task597);
  }

  shared_ptr<Tensor> I865;
  if (diagonal) {
    vector<IndexRange> I865_index;
    I865 = make_shared<Tensor>(I865_index);
  }
  shared_ptr<Task598> task598;
  if (diagonal) {
    auto tensor598 = vector<shared_ptr<Tensor>>{den0ci, I865};
    task598 = make_shared<Task598>(tensor598, cindex);
    task598->add_dep(task568);
    deci4q->add_task(task598);
  }

  shared_ptr<Task599> task599;
  if (diagonal) {
    auto tensor599 = vector<shared_ptr<Tensor>>{I865, v2_, l2};
    task599 = make_shared<Task599>(tensor599, cindex);
    task598->add_dep(task599);
    task599->add_dep(task568);
    deci4q->add_task(task599);
  }

  vector<IndexRange> I868_index = {active_, active_};
  auto I868 = make_shared<Tensor>(I868_index);
  auto tensor600 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I868};
  auto task600 = make_shared<Task600>(tensor600, cindex);
  task600->add_dep(task568);
  deci4q->add_task(task600);

  auto tensor601 = vector<shared_ptr<Tensor>>{I868, v2_, l2};
  auto task601 = make_shared<Task601>(tensor601, cindex);
  task600->add_dep(task601);
  task601->add_dep(task568);
  deci4q->add_task(task601);

  auto tensor602 = vector<shared_ptr<Tensor>>{I868, v2_, l2};
  auto task602 = make_shared<Task602>(tensor602, cindex);
  task600->add_dep(task602);
  task602->add_dep(task568);
  deci4q->add_task(task602);

  auto tensor603 = vector<shared_ptr<Tensor>>{I868, h1_, l2};
  auto task603 = make_shared<Task603>(tensor603, cindex);
  task600->add_dep(task603);
  task603->add_dep(task568);
  deci4q->add_task(task603);

  auto tensor604 = vector<shared_ptr<Tensor>>{I868, h1_, l2};
  auto task604 = make_shared<Task604>(tensor604, cindex);
  task600->add_dep(task604);
  task604->add_dep(task568);
  deci4q->add_task(task604);

  vector<IndexRange> I874_index = {active_, active_, active_, active_};
  auto I874 = make_shared<Tensor>(I874_index);
  auto tensor605 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I874};
  auto task605 = make_shared<Task605>(tensor605, cindex);
  task605->add_dep(task568);
  deci4q->add_task(task605);

  auto tensor606 = vector<shared_ptr<Tensor>>{I874, v2_, l2};
  auto task606 = make_shared<Task606>(tensor606, cindex);
  task605->add_dep(task606);
  task606->add_dep(task568);
  deci4q->add_task(task606);

  auto tensor607 = vector<shared_ptr<Tensor>>{I874, h1_, l2};
  auto task607 = make_shared<Task607>(tensor607, cindex);
  task605->add_dep(task607);
  task607->add_dep(task568);
  deci4q->add_task(task607);

  return deci4q;
}


#endif
