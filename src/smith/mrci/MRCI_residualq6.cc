//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_residualqq.cc
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


#include <src/smith/mrci/MRCI.h>
#include <src/smith/mrci/MRCI_tasks12.h>
#include <src/smith/mrci/MRCI_tasks13.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void MRCI::MRCI::make_residualq6(shared_ptr<Queue> residualq, shared_ptr<Task> task108, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I144_index = {virt_, active_, active_, active_};
  auto I144 = make_shared<Tensor>(I144_index);
  auto tensor552 = vector<shared_ptr<Tensor>>{r, I144};
  auto task552 = make_shared<Task552>(tensor552, pindex);
  task552->add_dep(task108);
  residualq->add_task(task552);

  vector<IndexRange> I145_index = {active_, active_, virt_, active_};
  auto I145 = make_shared<Tensor>(I145_index);
  auto tensor553 = vector<shared_ptr<Tensor>>{I144, Gamma48_(), I145};
  auto task553 = make_shared<Task553>(tensor553, pindex);
  task552->add_dep(task553);
  task553->add_dep(task108);
  residualq->add_task(task553);

  auto tensor554 = vector<shared_ptr<Tensor>>{I145, t2, h1_};
  auto task554 = make_shared<Task554>(tensor554, pindex);
  task553->add_dep(task554);
  task554->add_dep(task108);
  residualq->add_task(task554);

  auto tensor555 = vector<shared_ptr<Tensor>>{I145, t2, v2_};
  auto task555 = make_shared<Task555>(tensor555, pindex);
  task553->add_dep(task555);
  task555->add_dep(task108);
  residualq->add_task(task555);

  auto tensor556 = vector<shared_ptr<Tensor>>{I145, t2, v2_};
  auto task556 = make_shared<Task556>(tensor556, pindex);
  task553->add_dep(task556);
  task556->add_dep(task108);
  residualq->add_task(task556);

  vector<IndexRange> I148_index = {active_, virt_, active_, active_};
  auto I148 = make_shared<Tensor>(I148_index);
  auto tensor557 = vector<shared_ptr<Tensor>>{I144, Gamma49_(), I148};
  auto task557 = make_shared<Task557>(tensor557, pindex);
  task552->add_dep(task557);
  task557->add_dep(task108);
  residualq->add_task(task557);

  auto tensor558 = vector<shared_ptr<Tensor>>{I148, t2, h1_};
  auto task558 = make_shared<Task558>(tensor558, pindex);
  task557->add_dep(task558);
  task558->add_dep(task108);
  residualq->add_task(task558);

  auto tensor559 = vector<shared_ptr<Tensor>>{I148, t2, v2_};
  auto task559 = make_shared<Task559>(tensor559, pindex);
  task557->add_dep(task559);
  task559->add_dep(task108);
  residualq->add_task(task559);

  auto tensor560 = vector<shared_ptr<Tensor>>{I148, t2, v2_};
  auto task560 = make_shared<Task560>(tensor560, pindex);
  task557->add_dep(task560);
  task560->add_dep(task108);
  residualq->add_task(task560);

  auto tensor561 = vector<shared_ptr<Tensor>>{I148, t2, v2_};
  auto task561 = make_shared<Task561>(tensor561, pindex);
  task557->add_dep(task561);
  task561->add_dep(task108);
  residualq->add_task(task561);

  auto tensor562 = vector<shared_ptr<Tensor>>{I148, t2, v2_};
  auto task562 = make_shared<Task562>(tensor562, pindex);
  task557->add_dep(task562);
  task562->add_dep(task108);
  residualq->add_task(task562);

  vector<IndexRange> I151_index = {virt_, active_, active_, active_};
  auto I151 = make_shared<Tensor>(I151_index);
  auto tensor563 = vector<shared_ptr<Tensor>>{I144, Gamma50_(), I151};
  auto task563 = make_shared<Task563>(tensor563, pindex);
  task552->add_dep(task563);
  task563->add_dep(task108);
  residualq->add_task(task563);

  auto tensor564 = vector<shared_ptr<Tensor>>{I151, t2, h1_};
  auto task564 = make_shared<Task564>(tensor564, pindex);
  task563->add_dep(task564);
  task564->add_dep(task108);
  residualq->add_task(task564);

  auto tensor565 = vector<shared_ptr<Tensor>>{I151, t2, h1_};
  auto task565 = make_shared<Task565>(tensor565, pindex);
  task563->add_dep(task565);
  task565->add_dep(task108);
  residualq->add_task(task565);

  vector<IndexRange> I1075_index = {virt_, closed_, active_, active_};
  auto I1075 = make_shared<Tensor>(I1075_index);
  auto tensor566 = vector<shared_ptr<Tensor>>{I151, t2, I1075};
  auto task566 = make_shared<Task566>(tensor566, pindex);
  task563->add_dep(task566);
  task566->add_dep(task108);
  residualq->add_task(task566);

  auto tensor567 = vector<shared_ptr<Tensor>>{I1075, v2_};
  auto task567 = make_shared<Task567>(tensor567, pindex);
  task566->add_dep(task567);
  task567->add_dep(task108);
  residualq->add_task(task567);

  vector<IndexRange> I1078_index = {virt_, closed_, active_, active_};
  auto I1078 = make_shared<Tensor>(I1078_index);
  auto tensor568 = vector<shared_ptr<Tensor>>{I151, t2, I1078};
  auto task568 = make_shared<Task568>(tensor568, pindex);
  task563->add_dep(task568);
  task568->add_dep(task108);
  residualq->add_task(task568);

  auto tensor569 = vector<shared_ptr<Tensor>>{I1078, v2_};
  auto task569 = make_shared<Task569>(tensor569, pindex);
  task568->add_dep(task569);
  task569->add_dep(task108);
  residualq->add_task(task569);

  auto tensor570 = vector<shared_ptr<Tensor>>{I151, t2, v2_};
  auto task570 = make_shared<Task570>(tensor570, pindex);
  task563->add_dep(task570);
  task570->add_dep(task108);
  residualq->add_task(task570);

  auto tensor571 = vector<shared_ptr<Tensor>>{I151, t2, v2_};
  auto task571 = make_shared<Task571>(tensor571, pindex);
  task563->add_dep(task571);
  task571->add_dep(task108);
  residualq->add_task(task571);

  vector<IndexRange> I154_index = {active_, virt_};
  auto I154 = make_shared<Tensor>(I154_index);
  auto tensor572 = vector<shared_ptr<Tensor>>{I144, Gamma51_(), I154};
  auto task572 = make_shared<Task572>(tensor572, pindex);
  task552->add_dep(task572);
  task572->add_dep(task108);
  residualq->add_task(task572);

  auto tensor573 = vector<shared_ptr<Tensor>>{I154, t2, h1_};
  auto task573 = make_shared<Task573>(tensor573, pindex);
  task572->add_dep(task573);
  task573->add_dep(task108);
  residualq->add_task(task573);

  auto tensor574 = vector<shared_ptr<Tensor>>{I154, t2, h1_};
  auto task574 = make_shared<Task574>(tensor574, pindex);
  task572->add_dep(task574);
  task574->add_dep(task108);
  residualq->add_task(task574);

  auto tensor575 = vector<shared_ptr<Tensor>>{I154, t2, v2_};
  auto task575 = make_shared<Task575>(tensor575, pindex);
  task572->add_dep(task575);
  task575->add_dep(task108);
  residualq->add_task(task575);

  auto tensor576 = vector<shared_ptr<Tensor>>{I154, t2, v2_};
  auto task576 = make_shared<Task576>(tensor576, pindex);
  task572->add_dep(task576);
  task576->add_dep(task108);
  residualq->add_task(task576);

  auto tensor577 = vector<shared_ptr<Tensor>>{I154, t2, v2_};
  auto task577 = make_shared<Task577>(tensor577, pindex);
  task572->add_dep(task577);
  task577->add_dep(task108);
  residualq->add_task(task577);

  auto tensor578 = vector<shared_ptr<Tensor>>{I154, t2, v2_};
  auto task578 = make_shared<Task578>(tensor578, pindex);
  task572->add_dep(task578);
  task578->add_dep(task108);
  residualq->add_task(task578);

  vector<IndexRange> I1026_index = {closed_, active_, active_, active_, active_, active_};
  auto I1026 = make_shared<Tensor>(I1026_index);
  auto tensor579 = vector<shared_ptr<Tensor>>{I144, v2_, I1026};
  auto task579 = make_shared<Task579>(tensor579, pindex);
  task552->add_dep(task579);
  task579->add_dep(task108);
  residualq->add_task(task579);

  auto tensor580 = vector<shared_ptr<Tensor>>{I1026, Gamma339_(), t2};
  auto task580 = make_shared<Task580>(tensor580, pindex);
  task579->add_dep(task580);
  task580->add_dep(task108);
  residualq->add_task(task580);

  vector<IndexRange> I1029_index = {active_, active_, virt_, active_};
  auto I1029 = make_shared<Tensor>(I1029_index);
  auto tensor581 = vector<shared_ptr<Tensor>>{I144, Gamma340_(), I1029};
  auto task581 = make_shared<Task581>(tensor581, pindex);
  task552->add_dep(task581);
  task581->add_dep(task108);
  residualq->add_task(task581);

  auto tensor582 = vector<shared_ptr<Tensor>>{I1029, t2, v2_};
  auto task582 = make_shared<Task582>(tensor582, pindex);
  task581->add_dep(task582);
  task582->add_dep(task108);
  residualq->add_task(task582);

  vector<IndexRange> I1032_index = {closed_, active_, active_, active_, active_, active_};
  auto I1032 = make_shared<Tensor>(I1032_index);
  auto tensor583 = vector<shared_ptr<Tensor>>{I144, t2, I1032};
  auto task583 = make_shared<Task583>(tensor583, pindex);
  task552->add_dep(task583);
  task583->add_dep(task108);
  residualq->add_task(task583);

  auto tensor584 = vector<shared_ptr<Tensor>>{I1032, Gamma341_(), v2_};
  auto task584 = make_shared<Task584>(tensor584, pindex);
  task583->add_dep(task584);
  task584->add_dep(task108);
  residualq->add_task(task584);

  auto tensor585 = vector<shared_ptr<Tensor>>{I1032, Gamma342_(), v2_};
  auto task585 = make_shared<Task585>(tensor585, pindex);
  task583->add_dep(task585);
  task585->add_dep(task108);
  residualq->add_task(task585);

  vector<IndexRange> I1044_index = {closed_, active_, active_, active_, active_, active_};
  auto I1044 = make_shared<Tensor>(I1044_index);
  auto tensor586 = vector<shared_ptr<Tensor>>{I144, t2, I1044};
  auto task586 = make_shared<Task586>(tensor586, pindex);
  task552->add_dep(task586);
  task586->add_dep(task108);
  residualq->add_task(task586);

  auto tensor587 = vector<shared_ptr<Tensor>>{I1044, Gamma345_(), v2_};
  auto task587 = make_shared<Task587>(tensor587, pindex);
  task586->add_dep(task587);
  task587->add_dep(task108);
  residualq->add_task(task587);

  auto tensor588 = vector<shared_ptr<Tensor>>{I1044, Gamma346_(), v2_};
  auto task588 = make_shared<Task588>(tensor588, pindex);
  task586->add_dep(task588);
  task588->add_dep(task108);
  residualq->add_task(task588);

  vector<IndexRange> I1056_index = {virt_, active_, active_, active_, active_, active_};
  auto I1056 = make_shared<Tensor>(I1056_index);
  auto tensor589 = vector<shared_ptr<Tensor>>{I144, Gamma349_(), I1056};
  auto task589 = make_shared<Task589>(tensor589, pindex);
  task552->add_dep(task589);
  task589->add_dep(task108);
  residualq->add_task(task589);

  vector<IndexRange> I1057_index = {virt_, virt_, active_, active_};
  auto I1057 = make_shared<Tensor>(I1057_index);
  auto tensor590 = vector<shared_ptr<Tensor>>{I1056, t2, I1057};
  auto task590 = make_shared<Task590>(tensor590, pindex);
  task589->add_dep(task590);
  task590->add_dep(task108);
  residualq->add_task(task590);

  auto tensor591 = vector<shared_ptr<Tensor>>{I1057, v2_};
  auto task591 = make_shared<Task591>(tensor591, pindex);
  task590->add_dep(task591);
  task591->add_dep(task108);
  residualq->add_task(task591);

  auto tensor592 = vector<shared_ptr<Tensor>>{I1056, t2, v2_};
  auto task592 = make_shared<Task592>(tensor592, pindex);
  task589->add_dep(task592);
  task592->add_dep(task108);
  residualq->add_task(task592);

  vector<IndexRange> I1059_index = {active_, active_, virt_, active_, active_, active_};
  auto I1059 = make_shared<Tensor>(I1059_index);
  auto tensor593 = vector<shared_ptr<Tensor>>{I144, Gamma350_(), I1059};
  auto task593 = make_shared<Task593>(tensor593, pindex);
  task552->add_dep(task593);
  task593->add_dep(task108);
  residualq->add_task(task593);

  auto tensor594 = vector<shared_ptr<Tensor>>{I1059, t2, v2_};
  auto task594 = make_shared<Task594>(tensor594, pindex);
  task593->add_dep(task594);
  task594->add_dep(task108);
  residualq->add_task(task594);

  vector<IndexRange> I1062_index = {active_, virt_, active_, active_, active_, active_};
  auto I1062 = make_shared<Tensor>(I1062_index);
  auto tensor595 = vector<shared_ptr<Tensor>>{I144, Gamma351_(), I1062};
  auto task595 = make_shared<Task595>(tensor595, pindex);
  task552->add_dep(task595);
  task595->add_dep(task108);
  residualq->add_task(task595);

  auto tensor596 = vector<shared_ptr<Tensor>>{I1062, t2, v2_};
  auto task596 = make_shared<Task596>(tensor596, pindex);
  task595->add_dep(task596);
  task596->add_dep(task108);
  residualq->add_task(task596);

  vector<IndexRange> I1086_index = {active_, active_, active_, virt_};
  auto I1086 = make_shared<Tensor>(I1086_index);
  auto tensor597 = vector<shared_ptr<Tensor>>{I144, Gamma359_(), I1086};
  auto task597 = make_shared<Task597>(tensor597, pindex);
  task552->add_dep(task597);
  task597->add_dep(task108);
  residualq->add_task(task597);

  auto tensor598 = vector<shared_ptr<Tensor>>{I1086, t2, v2_};
  auto task598 = make_shared<Task598>(tensor598, pindex);
  task597->add_dep(task598);
  task598->add_dep(task108);
  residualq->add_task(task598);

  vector<IndexRange> I1107_index = {active_, active_, active_, active_, virt_, active_};
  auto I1107 = make_shared<Tensor>(I1107_index);
  auto tensor599 = vector<shared_ptr<Tensor>>{I144, Gamma366_(), I1107};
  auto task599 = make_shared<Task599>(tensor599, pindex);
  task552->add_dep(task599);
  task599->add_dep(task108);
  residualq->add_task(task599);

  auto tensor600 = vector<shared_ptr<Tensor>>{I1107, t2, v2_};
  auto task600 = make_shared<Task600>(tensor600, pindex);
  task599->add_dep(task600);
  task600->add_dep(task108);
  residualq->add_task(task600);

  auto tensor601 = vector<shared_ptr<Tensor>>{I144, Gamma556_(), t2};
  auto task601 = make_shared<Task601>(tensor601, pindex);
  task552->add_dep(task601);
  task601->add_dep(task108);
  residualq->add_task(task601);

  auto tensor602 = vector<shared_ptr<Tensor>>{I144, Gamma557_(), t2};
  auto task602 = make_shared<Task602>(tensor602, pindex);
  task552->add_dep(task602);
  task602->add_dep(task108);
  residualq->add_task(task602);
}

#endif
