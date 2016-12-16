//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_residualq5.cc
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


#include <src/smith/relmrci/RelMRCI.h>
#include <src/smith/relmrci/RelMRCI_tasks11.h>
#include <src/smith/relmrci/RelMRCI_tasks12.h>
#include <src/smith/relmrci/RelMRCI_tasks13.h>
#include <src/smith/relmrci/RelMRCI_tasks14.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void RelMRCI::RelMRCI::make_residualq5(shared_ptr<Queue> residualq, shared_ptr<Task> task83, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I134_index = {closed_, virt_, active_, virt_};
  auto I134 = make_shared<Tensor>(I134_index);
  auto tensor535 = vector<shared_ptr<Tensor>>{r, I134};
  auto task535 = make_shared<Task535>(tensor535, pindex);
  task535->add_dep(task83);
  residualq->add_task(task535);

  vector<IndexRange> I135_index = {closed_, virt_, active_, active_};
  auto I135 = make_shared<Tensor>(I135_index);
  auto tensor536 = vector<shared_ptr<Tensor>>{I134, h1_, I135};
  auto task536 = make_shared<Task536>(tensor536, pindex);
  task535->add_dep(task536);
  task536->add_dep(task83);
  residualq->add_task(task536);

  auto tensor537 = vector<shared_ptr<Tensor>>{I135, Gamma24_(), t2};
  auto task537 = make_shared<Task537>(tensor537, pindex);
  task536->add_dep(task537);
  task537->add_dep(task83);
  residualq->add_task(task537);

  vector<IndexRange> I138_index = {closed_, virt_, active_, active_};
  auto I138 = make_shared<Tensor>(I138_index);
  auto tensor538 = vector<shared_ptr<Tensor>>{I134, h1_, I138};
  auto task538 = make_shared<Task538>(tensor538, pindex);
  task535->add_dep(task538);
  task538->add_dep(task83);
  residualq->add_task(task538);

  auto tensor539 = vector<shared_ptr<Tensor>>{I138, Gamma24_(), t2};
  auto task539 = make_shared<Task539>(tensor539, pindex);
  task538->add_dep(task539);
  task539->add_dep(task83);
  residualq->add_task(task539);

  vector<IndexRange> I141_index = {virt_, active_};
  auto I141 = make_shared<Tensor>(I141_index);
  auto tensor540 = vector<shared_ptr<Tensor>>{I134, h1_, I141};
  auto task540 = make_shared<Task540>(tensor540, pindex);
  task535->add_dep(task540);
  task540->add_dep(task83);
  residualq->add_task(task540);

  auto tensor541 = vector<shared_ptr<Tensor>>{I141, Gamma33_(), t2};
  auto task541 = make_shared<Task541>(tensor541, pindex);
  task540->add_dep(task541);
  task541->add_dep(task83);
  residualq->add_task(task541);

  vector<IndexRange> I144_index = {virt_, active_};
  auto I144 = make_shared<Tensor>(I144_index);
  auto tensor542 = vector<shared_ptr<Tensor>>{I134, h1_, I144};
  auto task542 = make_shared<Task542>(tensor542, pindex);
  task535->add_dep(task542);
  task542->add_dep(task83);
  residualq->add_task(task542);

  auto tensor543 = vector<shared_ptr<Tensor>>{I144, Gamma33_(), t2};
  auto task543 = make_shared<Task543>(tensor543, pindex);
  task542->add_dep(task543);
  task543->add_dep(task83);
  residualq->add_task(task543);

  vector<IndexRange> I147_index = {closed_, active_};
  auto I147 = make_shared<Tensor>(I147_index);
  auto tensor544 = vector<shared_ptr<Tensor>>{I134, t2, I147};
  auto task544 = make_shared<Task544>(tensor544, pindex);
  task535->add_dep(task544);
  task544->add_dep(task83);
  residualq->add_task(task544);

  auto tensor545 = vector<shared_ptr<Tensor>>{I147, Gamma27_(), h1_};
  auto task545 = make_shared<Task545>(tensor545, pindex);
  task544->add_dep(task545);
  task545->add_dep(task83);
  residualq->add_task(task545);

  auto tensor546 = vector<shared_ptr<Tensor>>{I147, Gamma33_(), v2_};
  auto task546 = make_shared<Task546>(tensor546, pindex);
  task544->add_dep(task546);
  task546->add_dep(task83);
  residualq->add_task(task546);

  auto tensor547 = vector<shared_ptr<Tensor>>{I147, Gamma24_(), v2_};
  auto task547 = make_shared<Task547>(tensor547, pindex);
  task544->add_dep(task547);
  task547->add_dep(task83);
  residualq->add_task(task547);

  vector<IndexRange> I150_index = {closed_, active_};
  auto I150 = make_shared<Tensor>(I150_index);
  auto tensor548 = vector<shared_ptr<Tensor>>{I134, t2, I150};
  auto task548 = make_shared<Task548>(tensor548, pindex);
  task535->add_dep(task548);
  task548->add_dep(task83);
  residualq->add_task(task548);

  auto tensor549 = vector<shared_ptr<Tensor>>{I150, Gamma27_(), h1_};
  auto task549 = make_shared<Task549>(tensor549, pindex);
  task548->add_dep(task549);
  task549->add_dep(task83);
  residualq->add_task(task549);

  auto tensor550 = vector<shared_ptr<Tensor>>{I150, Gamma33_(), v2_};
  auto task550 = make_shared<Task550>(tensor550, pindex);
  task548->add_dep(task550);
  task550->add_dep(task83);
  residualq->add_task(task550);

  auto tensor551 = vector<shared_ptr<Tensor>>{I150, Gamma24_(), v2_};
  auto task551 = make_shared<Task551>(tensor551, pindex);
  task548->add_dep(task551);
  task551->add_dep(task83);
  residualq->add_task(task551);

  vector<IndexRange> I153_index = {closed_, active_, virt_, virt_};
  auto I153 = make_shared<Tensor>(I153_index);
  auto tensor552 = vector<shared_ptr<Tensor>>{I134, Gamma27_(), I153};
  auto task552 = make_shared<Task552>(tensor552, pindex);
  task535->add_dep(task552);
  task552->add_dep(task83);
  residualq->add_task(task552);

  auto tensor553 = vector<shared_ptr<Tensor>>{I153, t2, h1_};
  auto task553 = make_shared<Task553>(tensor553, pindex);
  task552->add_dep(task553);
  task553->add_dep(task83);
  residualq->add_task(task553);

  auto tensor554 = vector<shared_ptr<Tensor>>{I153, t2, h1_};
  auto task554 = make_shared<Task554>(tensor554, pindex);
  task552->add_dep(task554);
  task554->add_dep(task83);
  residualq->add_task(task554);

  auto tensor555 = vector<shared_ptr<Tensor>>{I153, t2, h1_};
  auto task555 = make_shared<Task555>(tensor555, pindex);
  task552->add_dep(task555);
  task555->add_dep(task83);
  residualq->add_task(task555);

  auto tensor556 = vector<shared_ptr<Tensor>>{I153, t2, h1_};
  auto task556 = make_shared<Task556>(tensor556, pindex);
  task552->add_dep(task556);
  task556->add_dep(task83);
  residualq->add_task(task556);

  auto tensor557 = vector<shared_ptr<Tensor>>{I153, t2, h1_};
  auto task557 = make_shared<Task557>(tensor557, pindex);
  task552->add_dep(task557);
  task557->add_dep(task83);
  residualq->add_task(task557);

  auto tensor558 = vector<shared_ptr<Tensor>>{I153, t2, h1_};
  auto task558 = make_shared<Task558>(tensor558, pindex);
  task552->add_dep(task558);
  task558->add_dep(task83);
  residualq->add_task(task558);

  auto tensor559 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task559 = make_shared<Task559>(tensor559, pindex);
  task552->add_dep(task559);
  task559->add_dep(task83);
  residualq->add_task(task559);

  auto tensor560 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task560 = make_shared<Task560>(tensor560, pindex);
  task552->add_dep(task560);
  task560->add_dep(task83);
  residualq->add_task(task560);

  auto tensor561 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task561 = make_shared<Task561>(tensor561, pindex);
  task552->add_dep(task561);
  task561->add_dep(task83);
  residualq->add_task(task561);

  auto tensor562 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task562 = make_shared<Task562>(tensor562, pindex);
  task552->add_dep(task562);
  task562->add_dep(task83);
  residualq->add_task(task562);

  auto tensor563 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task563 = make_shared<Task563>(tensor563, pindex);
  task552->add_dep(task563);
  task563->add_dep(task83);
  residualq->add_task(task563);

  auto tensor564 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task564 = make_shared<Task564>(tensor564, pindex);
  task552->add_dep(task564);
  task564->add_dep(task83);
  residualq->add_task(task564);

  auto tensor565 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task565 = make_shared<Task565>(tensor565, pindex);
  task552->add_dep(task565);
  task565->add_dep(task83);
  residualq->add_task(task565);

  auto tensor566 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task566 = make_shared<Task566>(tensor566, pindex);
  task552->add_dep(task566);
  task566->add_dep(task83);
  residualq->add_task(task566);

  auto tensor567 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task567 = make_shared<Task567>(tensor567, pindex);
  task552->add_dep(task567);
  task567->add_dep(task83);
  residualq->add_task(task567);

  auto tensor568 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task568 = make_shared<Task568>(tensor568, pindex);
  task552->add_dep(task568);
  task568->add_dep(task83);
  residualq->add_task(task568);

  auto tensor569 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task569 = make_shared<Task569>(tensor569, pindex);
  task552->add_dep(task569);
  task569->add_dep(task83);
  residualq->add_task(task569);

  auto tensor570 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task570 = make_shared<Task570>(tensor570, pindex);
  task552->add_dep(task570);
  task570->add_dep(task83);
  residualq->add_task(task570);

  auto tensor571 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task571 = make_shared<Task571>(tensor571, pindex);
  task552->add_dep(task571);
  task571->add_dep(task83);
  residualq->add_task(task571);

  auto tensor572 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task572 = make_shared<Task572>(tensor572, pindex);
  task552->add_dep(task572);
  task572->add_dep(task83);
  residualq->add_task(task572);

  auto tensor573 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task573 = make_shared<Task573>(tensor573, pindex);
  task552->add_dep(task573);
  task573->add_dep(task83);
  residualq->add_task(task573);

  auto tensor574 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task574 = make_shared<Task574>(tensor574, pindex);
  task552->add_dep(task574);
  task574->add_dep(task83);
  residualq->add_task(task574);

  auto tensor575 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task575 = make_shared<Task575>(tensor575, pindex);
  task552->add_dep(task575);
  task575->add_dep(task83);
  residualq->add_task(task575);

  auto tensor576 = vector<shared_ptr<Tensor>>{I153, t2, v2_};
  auto task576 = make_shared<Task576>(tensor576, pindex);
  task552->add_dep(task576);
  task576->add_dep(task83);
  residualq->add_task(task576);

  vector<IndexRange> I171_index = {virt_, virt_, active_, active_};
  auto I171 = make_shared<Tensor>(I171_index);
  auto tensor577 = vector<shared_ptr<Tensor>>{I134, h1_, I171};
  auto task577 = make_shared<Task577>(tensor577, pindex);
  task535->add_dep(task577);
  task577->add_dep(task83);
  residualq->add_task(task577);

  auto tensor578 = vector<shared_ptr<Tensor>>{I171, Gamma33_(), t2};
  auto task578 = make_shared<Task578>(tensor578, pindex);
  task577->add_dep(task578);
  task578->add_dep(task83);
  residualq->add_task(task578);

  vector<IndexRange> I980_index = {closed_, active_, active_, active_};
  auto I980 = make_shared<Tensor>(I980_index);
  auto tensor579 = vector<shared_ptr<Tensor>>{I134, v2_, I980};
  auto task579 = make_shared<Task579>(tensor579, pindex);
  task535->add_dep(task579);
  task579->add_dep(task83);
  residualq->add_task(task579);

  auto tensor580 = vector<shared_ptr<Tensor>>{I980, Gamma317_(), t2};
  auto task580 = make_shared<Task580>(tensor580, pindex);
  task579->add_dep(task580);
  task580->add_dep(task83);
  residualq->add_task(task580);

  vector<IndexRange> I983_index = {virt_, closed_, active_, active_};
  auto I983 = make_shared<Tensor>(I983_index);
  auto tensor581 = vector<shared_ptr<Tensor>>{I134, t2, I983};
  auto task581 = make_shared<Task581>(tensor581, pindex);
  task535->add_dep(task581);
  task581->add_dep(task83);
  residualq->add_task(task581);

  auto tensor582 = vector<shared_ptr<Tensor>>{I983, Gamma318_(), v2_};
  auto task582 = make_shared<Task582>(tensor582, pindex);
  task581->add_dep(task582);
  task582->add_dep(task83);
  residualq->add_task(task582);

  vector<IndexRange> I986_index = {virt_, closed_, active_, active_};
  auto I986 = make_shared<Tensor>(I986_index);
  auto tensor583 = vector<shared_ptr<Tensor>>{I134, t2, I986};
  auto task583 = make_shared<Task583>(tensor583, pindex);
  task535->add_dep(task583);
  task583->add_dep(task83);
  residualq->add_task(task583);

  auto tensor584 = vector<shared_ptr<Tensor>>{I986, Gamma5_(), v2_};
  auto task584 = make_shared<Task584>(tensor584, pindex);
  task583->add_dep(task584);
  task584->add_dep(task83);
  residualq->add_task(task584);

  vector<IndexRange> I989_index = {virt_, closed_, active_, active_};
  auto I989 = make_shared<Tensor>(I989_index);
  auto tensor585 = vector<shared_ptr<Tensor>>{I134, t2, I989};
  auto task585 = make_shared<Task585>(tensor585, pindex);
  task535->add_dep(task585);
  task585->add_dep(task83);
  residualq->add_task(task585);

  auto tensor586 = vector<shared_ptr<Tensor>>{I989, Gamma318_(), v2_};
  auto task586 = make_shared<Task586>(tensor586, pindex);
  task585->add_dep(task586);
  task586->add_dep(task83);
  residualq->add_task(task586);

  vector<IndexRange> I992_index = {virt_, closed_, active_, active_};
  auto I992 = make_shared<Tensor>(I992_index);
  auto tensor587 = vector<shared_ptr<Tensor>>{I134, t2, I992};
  auto task587 = make_shared<Task587>(tensor587, pindex);
  task535->add_dep(task587);
  task587->add_dep(task83);
  residualq->add_task(task587);

  auto tensor588 = vector<shared_ptr<Tensor>>{I992, Gamma318_(), v2_};
  auto task588 = make_shared<Task588>(tensor588, pindex);
  task587->add_dep(task588);
  task588->add_dep(task83);
  residualq->add_task(task588);

  vector<IndexRange> I995_index = {virt_, active_, active_, active_};
  auto I995 = make_shared<Tensor>(I995_index);
  auto tensor589 = vector<shared_ptr<Tensor>>{I134, t2, I995};
  auto task589 = make_shared<Task589>(tensor589, pindex);
  task535->add_dep(task589);
  task589->add_dep(task83);
  residualq->add_task(task589);

  auto tensor590 = vector<shared_ptr<Tensor>>{I995, Gamma31_(), v2_};
  auto task590 = make_shared<Task590>(tensor590, pindex);
  task589->add_dep(task590);
  task590->add_dep(task83);
  residualq->add_task(task590);

  auto tensor591 = vector<shared_ptr<Tensor>>{I995, Gamma193_(), v2_};
  auto task591 = make_shared<Task591>(tensor591, pindex);
  task589->add_dep(task591);
  task591->add_dep(task83);
  residualq->add_task(task591);

  vector<IndexRange> I998_index = {virt_, active_, active_, active_};
  auto I998 = make_shared<Tensor>(I998_index);
  auto tensor592 = vector<shared_ptr<Tensor>>{I134, t2, I998};
  auto task592 = make_shared<Task592>(tensor592, pindex);
  task535->add_dep(task592);
  task592->add_dep(task83);
  residualq->add_task(task592);

  auto tensor593 = vector<shared_ptr<Tensor>>{I998, Gamma31_(), v2_};
  auto task593 = make_shared<Task593>(tensor593, pindex);
  task592->add_dep(task593);
  task593->add_dep(task83);
  residualq->add_task(task593);

  auto tensor594 = vector<shared_ptr<Tensor>>{I998, Gamma193_(), v2_};
  auto task594 = make_shared<Task594>(tensor594, pindex);
  task592->add_dep(task594);
  task594->add_dep(task83);
  residualq->add_task(task594);

  vector<IndexRange> I1007_index = {closed_, virt_, active_, active_};
  auto I1007 = make_shared<Tensor>(I1007_index);
  auto tensor595 = vector<shared_ptr<Tensor>>{I134, v2_, I1007};
  auto task595 = make_shared<Task595>(tensor595, pindex);
  task535->add_dep(task595);
  task595->add_dep(task83);
  residualq->add_task(task595);

  auto tensor596 = vector<shared_ptr<Tensor>>{I1007, Gamma24_(), t2};
  auto task596 = make_shared<Task596>(tensor596, pindex);
  task595->add_dep(task596);
  task596->add_dep(task83);
  residualq->add_task(task596);

  vector<IndexRange> I1010_index = {closed_, virt_, active_, active_};
  auto I1010 = make_shared<Tensor>(I1010_index);
  auto tensor597 = vector<shared_ptr<Tensor>>{I134, v2_, I1010};
  auto task597 = make_shared<Task597>(tensor597, pindex);
  task535->add_dep(task597);
  task597->add_dep(task83);
  residualq->add_task(task597);

  auto tensor598 = vector<shared_ptr<Tensor>>{I1010, Gamma24_(), t2};
  auto task598 = make_shared<Task598>(tensor598, pindex);
  task597->add_dep(task598);
  task598->add_dep(task83);
  residualq->add_task(task598);

  vector<IndexRange> I1013_index = {closed_, virt_, active_, active_};
  auto I1013 = make_shared<Tensor>(I1013_index);
  auto tensor599 = vector<shared_ptr<Tensor>>{I134, v2_, I1013};
  auto task599 = make_shared<Task599>(tensor599, pindex);
  task535->add_dep(task599);
  task599->add_dep(task83);
  residualq->add_task(task599);

  auto tensor600 = vector<shared_ptr<Tensor>>{I1013, Gamma24_(), t2};
  auto task600 = make_shared<Task600>(tensor600, pindex);
  task599->add_dep(task600);
  task600->add_dep(task83);
  residualq->add_task(task600);

  vector<IndexRange> I1016_index = {closed_, virt_, active_, active_};
  auto I1016 = make_shared<Tensor>(I1016_index);
  auto tensor601 = vector<shared_ptr<Tensor>>{I134, v2_, I1016};
  auto task601 = make_shared<Task601>(tensor601, pindex);
  task535->add_dep(task601);
  task601->add_dep(task83);
  residualq->add_task(task601);

  auto tensor602 = vector<shared_ptr<Tensor>>{I1016, Gamma24_(), t2};
  auto task602 = make_shared<Task602>(tensor602, pindex);
  task601->add_dep(task602);
  task602->add_dep(task83);
  residualq->add_task(task602);

  vector<IndexRange> I1019_index = {closed_, virt_, active_, active_};
  auto I1019 = make_shared<Tensor>(I1019_index);
  auto tensor603 = vector<shared_ptr<Tensor>>{I134, v2_, I1019};
  auto task603 = make_shared<Task603>(tensor603, pindex);
  task535->add_dep(task603);
  task603->add_dep(task83);
  residualq->add_task(task603);

  auto tensor604 = vector<shared_ptr<Tensor>>{I1019, Gamma24_(), t2};
  auto task604 = make_shared<Task604>(tensor604, pindex);
  task603->add_dep(task604);
  task604->add_dep(task83);
  residualq->add_task(task604);

  vector<IndexRange> I1022_index = {closed_, virt_, active_, active_};
  auto I1022 = make_shared<Tensor>(I1022_index);
  auto tensor605 = vector<shared_ptr<Tensor>>{I134, v2_, I1022};
  auto task605 = make_shared<Task605>(tensor605, pindex);
  task535->add_dep(task605);
  task605->add_dep(task83);
  residualq->add_task(task605);

  auto tensor606 = vector<shared_ptr<Tensor>>{I1022, Gamma24_(), t2};
  auto task606 = make_shared<Task606>(tensor606, pindex);
  task605->add_dep(task606);
  task606->add_dep(task83);
  residualq->add_task(task606);

  vector<IndexRange> I1025_index = {virt_, active_, active_, active_};
  auto I1025 = make_shared<Tensor>(I1025_index);
  auto tensor607 = vector<shared_ptr<Tensor>>{I134, v2_, I1025};
  auto task607 = make_shared<Task607>(tensor607, pindex);
  task535->add_dep(task607);
  task607->add_dep(task83);
  residualq->add_task(task607);

  auto tensor608 = vector<shared_ptr<Tensor>>{I1025, Gamma32_(), t2};
  auto task608 = make_shared<Task608>(tensor608, pindex);
  task607->add_dep(task608);
  task608->add_dep(task83);
  residualq->add_task(task608);

  vector<IndexRange> I1028_index = {virt_, active_, active_, active_};
  auto I1028 = make_shared<Tensor>(I1028_index);
  auto tensor609 = vector<shared_ptr<Tensor>>{I134, v2_, I1028};
  auto task609 = make_shared<Task609>(tensor609, pindex);
  task535->add_dep(task609);
  task609->add_dep(task83);
  residualq->add_task(task609);

  auto tensor610 = vector<shared_ptr<Tensor>>{I1028, Gamma32_(), t2};
  auto task610 = make_shared<Task610>(tensor610, pindex);
  task609->add_dep(task610);
  task610->add_dep(task83);
  residualq->add_task(task610);

  vector<IndexRange> I1031_index = {virt_, active_, active_, active_};
  auto I1031 = make_shared<Tensor>(I1031_index);
  auto tensor611 = vector<shared_ptr<Tensor>>{I134, v2_, I1031};
  auto task611 = make_shared<Task611>(tensor611, pindex);
  task535->add_dep(task611);
  task611->add_dep(task83);
  residualq->add_task(task611);

  auto tensor612 = vector<shared_ptr<Tensor>>{I1031, Gamma26_(), t2};
  auto task612 = make_shared<Task612>(tensor612, pindex);
  task611->add_dep(task612);
  task612->add_dep(task83);
  residualq->add_task(task612);

  vector<IndexRange> I1034_index = {virt_, active_, active_, active_};
  auto I1034 = make_shared<Tensor>(I1034_index);
  auto tensor613 = vector<shared_ptr<Tensor>>{I134, v2_, I1034};
  auto task613 = make_shared<Task613>(tensor613, pindex);
  task535->add_dep(task613);
  task613->add_dep(task83);
  residualq->add_task(task613);

  auto tensor614 = vector<shared_ptr<Tensor>>{I1034, Gamma335_(), t2};
  auto task614 = make_shared<Task614>(tensor614, pindex);
  task613->add_dep(task614);
  task614->add_dep(task83);
  residualq->add_task(task614);

  vector<IndexRange> I1037_index = {virt_, active_, active_, active_};
  auto I1037 = make_shared<Tensor>(I1037_index);
  auto tensor615 = vector<shared_ptr<Tensor>>{I134, v2_, I1037};
  auto task615 = make_shared<Task615>(tensor615, pindex);
  task535->add_dep(task615);
  task615->add_dep(task83);
  residualq->add_task(task615);

  auto tensor616 = vector<shared_ptr<Tensor>>{I1037, Gamma336_(), t2};
  auto task616 = make_shared<Task616>(tensor616, pindex);
  task615->add_dep(task616);
  task616->add_dep(task83);
  residualq->add_task(task616);

  vector<IndexRange> I1040_index = {virt_, active_, active_, active_};
  auto I1040 = make_shared<Tensor>(I1040_index);
  auto tensor617 = vector<shared_ptr<Tensor>>{I134, v2_, I1040};
  auto task617 = make_shared<Task617>(tensor617, pindex);
  task535->add_dep(task617);
  task617->add_dep(task83);
  residualq->add_task(task617);

  auto tensor618 = vector<shared_ptr<Tensor>>{I1040, Gamma32_(), t2};
  auto task618 = make_shared<Task618>(tensor618, pindex);
  task617->add_dep(task618);
  task618->add_dep(task83);
  residualq->add_task(task618);

  vector<IndexRange> I1043_index = {virt_, active_, active_, active_};
  auto I1043 = make_shared<Tensor>(I1043_index);
  auto tensor619 = vector<shared_ptr<Tensor>>{I134, v2_, I1043};
  auto task619 = make_shared<Task619>(tensor619, pindex);
  task535->add_dep(task619);
  task619->add_dep(task83);
  residualq->add_task(task619);

  auto tensor620 = vector<shared_ptr<Tensor>>{I1043, Gamma32_(), t2};
  auto task620 = make_shared<Task620>(tensor620, pindex);
  task619->add_dep(task620);
  task620->add_dep(task83);
  residualq->add_task(task620);

  vector<IndexRange> I1046_index = {virt_, active_, active_, active_};
  auto I1046 = make_shared<Tensor>(I1046_index);
  auto tensor621 = vector<shared_ptr<Tensor>>{I134, v2_, I1046};
  auto task621 = make_shared<Task621>(tensor621, pindex);
  task535->add_dep(task621);
  task621->add_dep(task83);
  residualq->add_task(task621);

  auto tensor622 = vector<shared_ptr<Tensor>>{I1046, Gamma32_(), t2};
  auto task622 = make_shared<Task622>(tensor622, pindex);
  task621->add_dep(task622);
  task622->add_dep(task83);
  residualq->add_task(task622);

  vector<IndexRange> I1049_index = {virt_, active_};
  auto I1049 = make_shared<Tensor>(I1049_index);
  auto tensor623 = vector<shared_ptr<Tensor>>{I134, v2_, I1049};
  auto task623 = make_shared<Task623>(tensor623, pindex);
  task535->add_dep(task623);
  task623->add_dep(task83);
  residualq->add_task(task623);

  auto tensor624 = vector<shared_ptr<Tensor>>{I1049, Gamma33_(), t2};
  auto task624 = make_shared<Task624>(tensor624, pindex);
  task623->add_dep(task624);
  task624->add_dep(task83);
  residualq->add_task(task624);

  vector<IndexRange> I1052_index = {virt_, active_};
  auto I1052 = make_shared<Tensor>(I1052_index);
  auto tensor625 = vector<shared_ptr<Tensor>>{I134, v2_, I1052};
  auto task625 = make_shared<Task625>(tensor625, pindex);
  task535->add_dep(task625);
  task625->add_dep(task83);
  residualq->add_task(task625);

  auto tensor626 = vector<shared_ptr<Tensor>>{I1052, Gamma33_(), t2};
  auto task626 = make_shared<Task626>(tensor626, pindex);
  task625->add_dep(task626);
  task626->add_dep(task83);
  residualq->add_task(task626);

  vector<IndexRange> I1067_index = {closed_, closed_, closed_, active_};
  auto I1067 = make_shared<Tensor>(I1067_index);
  auto tensor627 = vector<shared_ptr<Tensor>>{I134, t2, I1067};
  auto task627 = make_shared<Task627>(tensor627, pindex);
  task535->add_dep(task627);
  task627->add_dep(task83);
  residualq->add_task(task627);

  auto tensor628 = vector<shared_ptr<Tensor>>{I1067, Gamma27_(), v2_};
  auto task628 = make_shared<Task628>(tensor628, pindex);
  task627->add_dep(task628);
  task628->add_dep(task83);
  residualq->add_task(task628);

  vector<IndexRange> I1070_index = {closed_, closed_, closed_, active_};
  auto I1070 = make_shared<Tensor>(I1070_index);
  auto tensor629 = vector<shared_ptr<Tensor>>{I134, t2, I1070};
  auto task629 = make_shared<Task629>(tensor629, pindex);
  task535->add_dep(task629);
  task629->add_dep(task83);
  residualq->add_task(task629);

  auto tensor630 = vector<shared_ptr<Tensor>>{I1070, Gamma27_(), v2_};
  auto task630 = make_shared<Task630>(tensor630, pindex);
  task629->add_dep(task630);
  task630->add_dep(task83);
  residualq->add_task(task630);

  vector<IndexRange> I1097_index = {closed_, closed_, active_, active_};
  auto I1097 = make_shared<Tensor>(I1097_index);
  auto tensor631 = vector<shared_ptr<Tensor>>{I134, t2, I1097};
  auto task631 = make_shared<Task631>(tensor631, pindex);
  task535->add_dep(task631);
  task631->add_dep(task83);
  residualq->add_task(task631);

  vector<IndexRange> I1098_index = {closed_, closed_, active_, active_};
  auto I1098 = make_shared<Tensor>(I1098_index);
  auto tensor632 = vector<shared_ptr<Tensor>>{I1097, Gamma33_(), I1098};
  auto task632 = make_shared<Task632>(tensor632, pindex);
  task631->add_dep(task632);
  task632->add_dep(task83);
  residualq->add_task(task632);

  auto tensor633 = vector<shared_ptr<Tensor>>{I1098, v2_};
  auto task633 = make_shared<Task633>(tensor633, pindex);
  task632->add_dep(task633);
  task633->add_dep(task83);
  residualq->add_task(task633);

  auto tensor634 = vector<shared_ptr<Tensor>>{I1097, Gamma24_(), v2_};
  auto task634 = make_shared<Task634>(tensor634, pindex);
  task631->add_dep(task634);
  task634->add_dep(task83);
  residualq->add_task(task634);

  auto tensor635 = vector<shared_ptr<Tensor>>{I1097, Gamma368_(), v2_};
  auto task635 = make_shared<Task635>(tensor635, pindex);
  task631->add_dep(task635);
  task635->add_dep(task83);
  residualq->add_task(task635);

  vector<IndexRange> I1100_index = {closed_, closed_, active_, active_};
  auto I1100 = make_shared<Tensor>(I1100_index);
  auto tensor636 = vector<shared_ptr<Tensor>>{I134, t2, I1100};
  auto task636 = make_shared<Task636>(tensor636, pindex);
  task535->add_dep(task636);
  task636->add_dep(task83);
  residualq->add_task(task636);

  vector<IndexRange> I1101_index = {closed_, closed_, active_, active_};
  auto I1101 = make_shared<Tensor>(I1101_index);
  auto tensor637 = vector<shared_ptr<Tensor>>{I1100, Gamma33_(), I1101};
  auto task637 = make_shared<Task637>(tensor637, pindex);
  task636->add_dep(task637);
  task637->add_dep(task83);
  residualq->add_task(task637);

  auto tensor638 = vector<shared_ptr<Tensor>>{I1101, v2_};
  auto task638 = make_shared<Task638>(tensor638, pindex);
  task637->add_dep(task638);
  task638->add_dep(task83);
  residualq->add_task(task638);

  auto tensor639 = vector<shared_ptr<Tensor>>{I1100, Gamma363_(), v2_};
  auto task639 = make_shared<Task639>(tensor639, pindex);
  task636->add_dep(task639);
  task639->add_dep(task83);
  residualq->add_task(task639);

  vector<IndexRange> I1103_index = {virt_, virt_, active_, active_};
  auto I1103 = make_shared<Tensor>(I1103_index);
  auto tensor640 = vector<shared_ptr<Tensor>>{I134, t2, I1103};
  auto task640 = make_shared<Task640>(tensor640, pindex);
  task535->add_dep(task640);
  task640->add_dep(task83);
  residualq->add_task(task640);

  vector<IndexRange> I1104_index = {virt_, virt_, active_, active_};
  auto I1104 = make_shared<Tensor>(I1104_index);
  auto tensor641 = vector<shared_ptr<Tensor>>{I1103, Gamma33_(), I1104};
  auto task641 = make_shared<Task641>(tensor641, pindex);
  task640->add_dep(task641);
  task641->add_dep(task83);
  residualq->add_task(task641);

  auto tensor642 = vector<shared_ptr<Tensor>>{I1104, v2_};
  auto task642 = make_shared<Task642>(tensor642, pindex);
  task641->add_dep(task642);
  task642->add_dep(task83);
  residualq->add_task(task642);

  auto tensor643 = vector<shared_ptr<Tensor>>{I1103, Gamma24_(), v2_};
  auto task643 = make_shared<Task643>(tensor643, pindex);
  task640->add_dep(task643);
  task643->add_dep(task83);
  residualq->add_task(task643);

  auto tensor644 = vector<shared_ptr<Tensor>>{I1103, Gamma368_(), v2_};
  auto task644 = make_shared<Task644>(tensor644, pindex);
  task640->add_dep(task644);
  task644->add_dep(task83);
  residualq->add_task(task644);

  vector<IndexRange> I1106_index = {virt_, virt_, active_, active_};
  auto I1106 = make_shared<Tensor>(I1106_index);
  auto tensor645 = vector<shared_ptr<Tensor>>{I134, t2, I1106};
  auto task645 = make_shared<Task645>(tensor645, pindex);
  task535->add_dep(task645);
  task645->add_dep(task83);
  residualq->add_task(task645);

  vector<IndexRange> I1107_index = {virt_, virt_, active_, active_};
  auto I1107 = make_shared<Tensor>(I1107_index);
  auto tensor646 = vector<shared_ptr<Tensor>>{I1106, Gamma33_(), I1107};
  auto task646 = make_shared<Task646>(tensor646, pindex);
  task645->add_dep(task646);
  task646->add_dep(task83);
  residualq->add_task(task646);

  auto tensor647 = vector<shared_ptr<Tensor>>{I1107, v2_};
  auto task647 = make_shared<Task647>(tensor647, pindex);
  task646->add_dep(task647);
  task647->add_dep(task83);
  residualq->add_task(task647);

  auto tensor648 = vector<shared_ptr<Tensor>>{I1106, Gamma24_(), v2_};
  auto task648 = make_shared<Task648>(tensor648, pindex);
  task645->add_dep(task648);
  task648->add_dep(task83);
  residualq->add_task(task648);

  auto tensor649 = vector<shared_ptr<Tensor>>{I1106, Gamma368_(), v2_};
  auto task649 = make_shared<Task649>(tensor649, pindex);
  task645->add_dep(task649);
  task649->add_dep(task83);
  residualq->add_task(task649);

  vector<IndexRange> I1109_index = {virt_, virt_, active_, active_};
  auto I1109 = make_shared<Tensor>(I1109_index);
  auto tensor650 = vector<shared_ptr<Tensor>>{I134, t2, I1109};
  auto task650 = make_shared<Task650>(tensor650, pindex);
  task535->add_dep(task650);
  task650->add_dep(task83);
  residualq->add_task(task650);

  vector<IndexRange> I1110_index = {virt_, virt_, active_, active_};
  auto I1110 = make_shared<Tensor>(I1110_index);
  auto tensor651 = vector<shared_ptr<Tensor>>{I1109, Gamma33_(), I1110};
  auto task651 = make_shared<Task651>(tensor651, pindex);
  task650->add_dep(task651);
  task651->add_dep(task83);
  residualq->add_task(task651);

  auto tensor652 = vector<shared_ptr<Tensor>>{I1110, v2_};
  auto task652 = make_shared<Task652>(tensor652, pindex);
  task651->add_dep(task652);
  task652->add_dep(task83);
  residualq->add_task(task652);

  auto tensor653 = vector<shared_ptr<Tensor>>{I1109, Gamma24_(), v2_};
  auto task653 = make_shared<Task653>(tensor653, pindex);
  task650->add_dep(task653);
  task653->add_dep(task83);
  residualq->add_task(task653);

  auto tensor654 = vector<shared_ptr<Tensor>>{I1109, Gamma368_(), v2_};
  auto task654 = make_shared<Task654>(tensor654, pindex);
  task650->add_dep(task654);
  task654->add_dep(task83);
  residualq->add_task(task654);

  vector<IndexRange> I1112_index = {virt_, virt_, active_, active_};
  auto I1112 = make_shared<Tensor>(I1112_index);
  auto tensor655 = vector<shared_ptr<Tensor>>{I134, t2, I1112};
  auto task655 = make_shared<Task655>(tensor655, pindex);
  task535->add_dep(task655);
  task655->add_dep(task83);
  residualq->add_task(task655);

  vector<IndexRange> I1113_index = {virt_, virt_, active_, active_};
  auto I1113 = make_shared<Tensor>(I1113_index);
  auto tensor656 = vector<shared_ptr<Tensor>>{I1112, Gamma33_(), I1113};
  auto task656 = make_shared<Task656>(tensor656, pindex);
  task655->add_dep(task656);
  task656->add_dep(task83);
  residualq->add_task(task656);

  auto tensor657 = vector<shared_ptr<Tensor>>{I1113, v2_};
  auto task657 = make_shared<Task657>(tensor657, pindex);
  task656->add_dep(task657);
  task657->add_dep(task83);
  residualq->add_task(task657);

  auto tensor658 = vector<shared_ptr<Tensor>>{I1112, Gamma363_(), v2_};
  auto task658 = make_shared<Task658>(tensor658, pindex);
  task655->add_dep(task658);
  task658->add_dep(task83);
  residualq->add_task(task658);

  auto tensor659 = vector<shared_ptr<Tensor>>{I1112, v2_, Gamma33_()};
  auto task659 = make_shared<Task659>(tensor659, pindex);
  task655->add_dep(task659);
  task659->add_dep(task83);
  residualq->add_task(task659);

  vector<IndexRange> I1199_index = {closed_, active_, active_, active_};
  auto I1199 = make_shared<Tensor>(I1199_index);
  auto tensor660 = vector<shared_ptr<Tensor>>{I134, t2, I1199};
  auto task660 = make_shared<Task660>(tensor660, pindex);
  task535->add_dep(task660);
  task660->add_dep(task83);
  residualq->add_task(task660);

  auto tensor661 = vector<shared_ptr<Tensor>>{I1199, Gamma32_(), v2_};
  auto task661 = make_shared<Task661>(tensor661, pindex);
  task660->add_dep(task661);
  task661->add_dep(task83);
  residualq->add_task(task661);

  auto tensor662 = vector<shared_ptr<Tensor>>{I1199, Gamma391_(), v2_};
  auto task662 = make_shared<Task662>(tensor662, pindex);
  task660->add_dep(task662);
  task662->add_dep(task83);
  residualq->add_task(task662);

  vector<IndexRange> I1205_index = {virt_, virt_, active_, active_};
  auto I1205 = make_shared<Tensor>(I1205_index);
  auto tensor663 = vector<shared_ptr<Tensor>>{I134, v2_, I1205};
  auto task663 = make_shared<Task663>(tensor663, pindex);
  task535->add_dep(task663);
  task663->add_dep(task83);
  residualq->add_task(task663);

  auto tensor664 = vector<shared_ptr<Tensor>>{I1205, Gamma33_(), t2};
  auto task664 = make_shared<Task664>(tensor664, pindex);
  task663->add_dep(task664);
  task664->add_dep(task83);
  residualq->add_task(task664);

  vector<IndexRange> I1208_index = {virt_, virt_, active_, active_};
  auto I1208 = make_shared<Tensor>(I1208_index);
  auto tensor665 = vector<shared_ptr<Tensor>>{I134, v2_, I1208};
  auto task665 = make_shared<Task665>(tensor665, pindex);
  task535->add_dep(task665);
  task665->add_dep(task83);
  residualq->add_task(task665);

  auto tensor666 = vector<shared_ptr<Tensor>>{I1208, Gamma33_(), t2};
  auto task666 = make_shared<Task666>(tensor666, pindex);
  task665->add_dep(task666);
  task666->add_dep(task83);
  residualq->add_task(task666);

  vector<IndexRange> I1211_index = {virt_, virt_, active_, active_};
  auto I1211 = make_shared<Tensor>(I1211_index);
  auto tensor667 = vector<shared_ptr<Tensor>>{I134, v2_, I1211};
  auto task667 = make_shared<Task667>(tensor667, pindex);
  task535->add_dep(task667);
  task667->add_dep(task83);
  residualq->add_task(task667);

  auto tensor668 = vector<shared_ptr<Tensor>>{I1211, Gamma368_(), t2};
  auto task668 = make_shared<Task668>(tensor668, pindex);
  task667->add_dep(task668);
  task668->add_dep(task83);
  residualq->add_task(task668);

  vector<IndexRange> I1214_index = {virt_, virt_, active_, active_};
  auto I1214 = make_shared<Tensor>(I1214_index);
  auto tensor669 = vector<shared_ptr<Tensor>>{I134, v2_, I1214};
  auto task669 = make_shared<Task669>(tensor669, pindex);
  task535->add_dep(task669);
  task669->add_dep(task83);
  residualq->add_task(task669);

  auto tensor670 = vector<shared_ptr<Tensor>>{I1214, Gamma33_(), t2};
  auto task670 = make_shared<Task670>(tensor670, pindex);
  task669->add_dep(task670);
  task670->add_dep(task83);
  residualq->add_task(task670);

  vector<IndexRange> I1297_index = {active_, virt_, closed_, virt_};
  auto I1297 = make_shared<Tensor>(I1297_index);
  auto tensor671 = vector<shared_ptr<Tensor>>{I134, Gamma428_(), I1297};
  auto task671 = make_shared<Task671>(tensor671, pindex);
  task535->add_dep(task671);
  task671->add_dep(task83);
  residualq->add_task(task671);

  auto tensor672 = vector<shared_ptr<Tensor>>{I1297, t2};
  auto task672 = make_shared<Task672>(tensor672, pindex);
  task671->add_dep(task672);
  task672->add_dep(task83);
  residualq->add_task(task672);

  vector<IndexRange> I1301_index = {active_, virt_, closed_, virt_};
  auto I1301 = make_shared<Tensor>(I1301_index);
  auto tensor673 = vector<shared_ptr<Tensor>>{I134, Gamma430_(), I1301};
  auto task673 = make_shared<Task673>(tensor673, pindex);
  task535->add_dep(task673);
  task673->add_dep(task83);
  residualq->add_task(task673);

  auto tensor674 = vector<shared_ptr<Tensor>>{I1301, t2};
  auto task674 = make_shared<Task674>(tensor674, pindex);
  task673->add_dep(task674);
  task674->add_dep(task83);
  residualq->add_task(task674);
}

#endif
