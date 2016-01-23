//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_deciqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_deciq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deciq = make_shared<Queue>();
  auto tensor537 = vector<shared_ptr<Tensor>>{deci};
  auto task537 = make_shared<Task537>(tensor537, reset);
  deciq->add_task(task537);

  vector<IndexRange> I698_index = {ci_};
  auto I698 = make_shared<Tensor>(I698_index);
  auto tensor538 = vector<shared_ptr<Tensor>>{deci, I698};
  auto task538 = make_shared<Task538>(tensor538, cindex);
  task538->add_dep(task537);
  deciq->add_task(task538);

  vector<IndexRange> I699_index = {active_, active_, active_, active_};
  auto I699 = make_shared<Tensor>(I699_index);
  auto tensor539 = vector<shared_ptr<Tensor>>{I698, Gamma248_(), I699};
  auto task539 = make_shared<Task539>(tensor539, cindex);
  task538->add_dep(task539);
  task539->add_dep(task537);
  deciq->add_task(task539);

  auto tensor540 = vector<shared_ptr<Tensor>>{I699, t2};
  auto task540 = make_shared<Task540>(tensor540, cindex);
  task539->add_dep(task540);
  task540->add_dep(task537);
  deciq->add_task(task540);

  vector<IndexRange> I702_index = {active_, active_, active_, active_};
  auto I702 = make_shared<Tensor>(I702_index);
  auto tensor541 = vector<shared_ptr<Tensor>>{I698, Gamma249_(), I702};
  auto task541 = make_shared<Task541>(tensor541, cindex);
  task538->add_dep(task541);
  task541->add_dep(task537);
  deciq->add_task(task541);

  vector<IndexRange> I703_index = {active_, active_, closed_, closed_};
  auto I703 = make_shared<Tensor>(I703_index);
  auto tensor542 = vector<shared_ptr<Tensor>>{I702, t2, I703};
  auto task542 = make_shared<Task542>(tensor542, cindex);
  task541->add_dep(task542);
  task542->add_dep(task537);
  deciq->add_task(task542);

  auto tensor543 = vector<shared_ptr<Tensor>>{I703, f1_, t2};
  auto task543 = make_shared<Task543>(tensor543, cindex);
  task542->add_dep(task543);
  task543->add_dep(task537);
  deciq->add_task(task543);

  auto tensor544 = vector<shared_ptr<Tensor>>{I702, t2};
  auto task544 = make_shared<Task544>(tensor544, cindex, this->e0_);
  task541->add_dep(task544);
  task544->add_dep(task537);
  deciq->add_task(task544);

  auto tensor545 = vector<shared_ptr<Tensor>>{I702, v2_, t2};
  auto task545 = make_shared<Task545>(tensor545, cindex);
  task541->add_dep(task545);
  task545->add_dep(task537);
  deciq->add_task(task545);

  auto tensor546 = vector<shared_ptr<Tensor>>{I702, v2_, t2};
  auto task546 = make_shared<Task546>(tensor546, cindex);
  task541->add_dep(task546);
  task546->add_dep(task537);
  deciq->add_task(task546);

  vector<IndexRange> I706_index = {active_, active_, active_, active_, active_, active_};
  auto I706 = make_shared<Tensor>(I706_index);
  auto tensor547 = vector<shared_ptr<Tensor>>{I698, Gamma250_(), I706};
  auto task547 = make_shared<Task547>(tensor547, cindex);
  task538->add_dep(task547);
  task547->add_dep(task537);
  deciq->add_task(task547);

  vector<IndexRange> I707_index = {active_, active_, closed_, active_};
  auto I707 = make_shared<Tensor>(I707_index);
  auto tensor548 = vector<shared_ptr<Tensor>>{I706, t2, I707};
  auto task548 = make_shared<Task548>(tensor548, cindex);
  task547->add_dep(task548);
  task548->add_dep(task537);
  deciq->add_task(task548);

  auto tensor549 = vector<shared_ptr<Tensor>>{I707, f1_, t2};
  auto task549 = make_shared<Task549>(tensor549, cindex);
  task548->add_dep(task549);
  task549->add_dep(task537);
  deciq->add_task(task549);

  vector<IndexRange> I710_index = {active_, active_, active_, active_};
  auto I710 = make_shared<Tensor>(I710_index);
  auto tensor550 = vector<shared_ptr<Tensor>>{I698, Gamma251_(), I710};
  auto task550 = make_shared<Task550>(tensor550, cindex);
  task538->add_dep(task550);
  task550->add_dep(task537);
  deciq->add_task(task550);

  vector<IndexRange> I711_index = {active_, closed_, closed_, active_};
  auto I711 = make_shared<Tensor>(I711_index);
  auto tensor551 = vector<shared_ptr<Tensor>>{I710, t2, I711};
  auto task551 = make_shared<Task551>(tensor551, cindex);
  task550->add_dep(task551);
  task551->add_dep(task537);
  deciq->add_task(task551);

  auto tensor552 = vector<shared_ptr<Tensor>>{I711, t2, f1_};
  auto task552 = make_shared<Task552>(tensor552, cindex);
  task551->add_dep(task552);
  task552->add_dep(task537);
  deciq->add_task(task552);

  vector<IndexRange> I742_index = {active_, closed_, closed_, active_};
  auto I742 = make_shared<Tensor>(I742_index);
  auto tensor553 = vector<shared_ptr<Tensor>>{I710, t2, I742};
  auto task553 = make_shared<Task553>(tensor553, cindex);
  task550->add_dep(task553);
  task553->add_dep(task537);
  deciq->add_task(task553);

  auto tensor554 = vector<shared_ptr<Tensor>>{I742, f1_, t2};
  auto task554 = make_shared<Task554>(tensor554, cindex);
  task553->add_dep(task554);
  task554->add_dep(task537);
  deciq->add_task(task554);

  vector<IndexRange> I714_index = {active_, active_, active_, active_, active_, active_};
  auto I714 = make_shared<Tensor>(I714_index);
  auto tensor555 = vector<shared_ptr<Tensor>>{I698, Gamma252_(), I714};
  auto task555 = make_shared<Task555>(tensor555, cindex);
  task538->add_dep(task555);
  task555->add_dep(task537);
  deciq->add_task(task555);

  vector<IndexRange> I715_index = {active_, closed_, active_, active_};
  auto I715 = make_shared<Tensor>(I715_index);
  auto tensor556 = vector<shared_ptr<Tensor>>{I714, t2, I715};
  auto task556 = make_shared<Task556>(tensor556, cindex);
  task555->add_dep(task556);
  task556->add_dep(task537);
  deciq->add_task(task556);

  auto tensor557 = vector<shared_ptr<Tensor>>{I715, t2, f1_};
  auto task557 = make_shared<Task557>(tensor557, cindex);
  task556->add_dep(task557);
  task557->add_dep(task537);
  deciq->add_task(task557);

  vector<IndexRange> I718_index = {active_, active_, active_, active_, active_, active_};
  auto I718 = make_shared<Tensor>(I718_index);
  auto tensor558 = vector<shared_ptr<Tensor>>{I698, Gamma253_(), I718};
  auto task558 = make_shared<Task558>(tensor558, cindex);
  task538->add_dep(task558);
  task558->add_dep(task537);
  deciq->add_task(task558);

  auto tensor559 = vector<shared_ptr<Tensor>>{I718, t2};
  auto task559 = make_shared<Task559>(tensor559, cindex);
  task558->add_dep(task559);
  task559->add_dep(task537);
  deciq->add_task(task559);

  vector<IndexRange> I721_index = {active_, active_, active_, active_, active_, active_};
  auto I721 = make_shared<Tensor>(I721_index);
  auto tensor560 = vector<shared_ptr<Tensor>>{I698, Gamma254_(), I721};
  auto task560 = make_shared<Task560>(tensor560, cindex);
  task538->add_dep(task560);
  task560->add_dep(task537);
  deciq->add_task(task560);

  vector<IndexRange> I722_index = {active_, active_, active_, closed_};
  auto I722 = make_shared<Tensor>(I722_index);
  auto tensor561 = vector<shared_ptr<Tensor>>{I721, t2, I722};
  auto task561 = make_shared<Task561>(tensor561, cindex);
  task560->add_dep(task561);
  task561->add_dep(task537);
  deciq->add_task(task561);

  auto tensor562 = vector<shared_ptr<Tensor>>{I722, f1_, t2};
  auto task562 = make_shared<Task562>(tensor562, cindex);
  task561->add_dep(task562);
  task562->add_dep(task537);
  deciq->add_task(task562);

  auto tensor563 = vector<shared_ptr<Tensor>>{I722, t2, f1_};
  auto task563 = make_shared<Task563>(tensor563, cindex);
  task561->add_dep(task563);
  task563->add_dep(task537);
  deciq->add_task(task563);

  vector<IndexRange> I862_index = {active_, active_, closed_, active_};
  auto I862 = make_shared<Tensor>(I862_index);
  auto tensor564 = vector<shared_ptr<Tensor>>{I721, t2, I862};
  auto task564 = make_shared<Task564>(tensor564, cindex);
  task560->add_dep(task564);
  task564->add_dep(task537);
  deciq->add_task(task564);

  auto tensor565 = vector<shared_ptr<Tensor>>{I862, t2};
  auto task565 = make_shared<Task565>(tensor565, cindex, this->e0_);
  task564->add_dep(task565);
  task565->add_dep(task537);
  deciq->add_task(task565);

  auto tensor566 = vector<shared_ptr<Tensor>>{I862, f1_, t2};
  auto task566 = make_shared<Task566>(tensor566, cindex);
  task564->add_dep(task566);
  task566->add_dep(task537);
  deciq->add_task(task566);

  auto tensor567 = vector<shared_ptr<Tensor>>{I721, v2_, t2};
  auto task567 = make_shared<Task567>(tensor567, cindex);
  task560->add_dep(task567);
  task567->add_dep(task537);
  deciq->add_task(task567);

  auto tensor568 = vector<shared_ptr<Tensor>>{I721, v2_, t2};
  auto task568 = make_shared<Task568>(tensor568, cindex);
  task560->add_dep(task568);
  task568->add_dep(task537);
  deciq->add_task(task568);

  vector<IndexRange> I725_index = {active_, active_, active_, active_};
  auto I725 = make_shared<Tensor>(I725_index);
  auto tensor569 = vector<shared_ptr<Tensor>>{I698, Gamma255_(), I725};
  auto task569 = make_shared<Task569>(tensor569, cindex);
  task538->add_dep(task569);
  task569->add_dep(task537);
  deciq->add_task(task569);

  vector<IndexRange> I726_index = {closed_, active_};
  auto I726 = make_shared<Tensor>(I726_index);
  auto tensor570 = vector<shared_ptr<Tensor>>{I725, t2, I726};
  auto task570 = make_shared<Task570>(tensor570, cindex);
  task569->add_dep(task570);
  task570->add_dep(task537);
  deciq->add_task(task570);

  auto tensor571 = vector<shared_ptr<Tensor>>{I726, t2, f1_};
  auto task571 = make_shared<Task571>(tensor571, cindex);
  task570->add_dep(task571);
  task571->add_dep(task537);
  deciq->add_task(task571);

  auto tensor572 = vector<shared_ptr<Tensor>>{I726, t2, f1_};
  auto task572 = make_shared<Task572>(tensor572, cindex);
  task570->add_dep(task572);
  task572->add_dep(task537);
  deciq->add_task(task572);

  vector<IndexRange> I816_index = {active_, closed_, virt_, active_};
  auto I816 = make_shared<Tensor>(I816_index);
  auto tensor573 = vector<shared_ptr<Tensor>>{I725, t2, I816};
  auto task573 = make_shared<Task573>(tensor573, cindex);
  task569->add_dep(task573);
  task573->add_dep(task537);
  deciq->add_task(task573);

  auto tensor574 = vector<shared_ptr<Tensor>>{I816, t2, f1_};
  auto task574 = make_shared<Task574>(tensor574, cindex);
  task573->add_dep(task574);
  task574->add_dep(task537);
  deciq->add_task(task574);

  vector<IndexRange> I866_index = {active_, virt_, closed_, active_};
  auto I866 = make_shared<Tensor>(I866_index);
  auto tensor575 = vector<shared_ptr<Tensor>>{I725, t2, I866};
  auto task575 = make_shared<Task575>(tensor575, cindex);
  task569->add_dep(task575);
  task575->add_dep(task537);
  deciq->add_task(task575);

  auto tensor576 = vector<shared_ptr<Tensor>>{I866, t2, f1_};
  auto task576 = make_shared<Task576>(tensor576, cindex);
  task575->add_dep(task576);
  task576->add_dep(task537);
  deciq->add_task(task576);

  auto tensor577 = vector<shared_ptr<Tensor>>{I866, t2, f1_};
  auto task577 = make_shared<Task577>(tensor577, cindex);
  task575->add_dep(task577);
  task577->add_dep(task537);
  deciq->add_task(task577);

  auto tensor578 = vector<shared_ptr<Tensor>>{I725, v2_, t2};
  auto task578 = make_shared<Task578>(tensor578, cindex);
  task569->add_dep(task578);
  task578->add_dep(task537);
  deciq->add_task(task578);

  auto tensor579 = vector<shared_ptr<Tensor>>{I725, h1_, t2};
  auto task579 = make_shared<Task579>(tensor579, cindex);
  task569->add_dep(task579);
  task579->add_dep(task537);
  deciq->add_task(task579);

  vector<IndexRange> I733_index = {active_, active_, active_, active_, active_, active_};
  auto I733 = make_shared<Tensor>(I733_index);
  auto tensor580 = vector<shared_ptr<Tensor>>{I698, Gamma257_(), I733};
  auto task580 = make_shared<Task580>(tensor580, cindex);
  task538->add_dep(task580);
  task580->add_dep(task537);
  deciq->add_task(task580);

  vector<IndexRange> I734_index = {active_, active_, closed_, active_};
  auto I734 = make_shared<Tensor>(I734_index);
  auto tensor581 = vector<shared_ptr<Tensor>>{I733, t2, I734};
  auto task581 = make_shared<Task581>(tensor581, cindex);
  task580->add_dep(task581);
  task581->add_dep(task537);
  deciq->add_task(task581);

  auto tensor582 = vector<shared_ptr<Tensor>>{I734, t2, f1_};
  auto task582 = make_shared<Task582>(tensor582, cindex);
  task581->add_dep(task582);
  task582->add_dep(task537);
  deciq->add_task(task582);

  vector<IndexRange> I745_index = {active_, active_, active_, active_};
  auto I745 = make_shared<Tensor>(I745_index);
  auto tensor583 = vector<shared_ptr<Tensor>>{I698, Gamma260_(), I745};
  auto task583 = make_shared<Task583>(tensor583, cindex);
  task538->add_dep(task583);
  task583->add_dep(task537);
  deciq->add_task(task583);

  vector<IndexRange> I746_index = {active_, closed_};
  auto I746 = make_shared<Tensor>(I746_index);
  auto tensor584 = vector<shared_ptr<Tensor>>{I745, t2, I746};
  auto task584 = make_shared<Task584>(tensor584, cindex);
  task583->add_dep(task584);
  task584->add_dep(task537);
  deciq->add_task(task584);

  auto tensor585 = vector<shared_ptr<Tensor>>{I746, f1_, t2};
  auto task585 = make_shared<Task585>(tensor585, cindex);
  task584->add_dep(task585);
  task585->add_dep(task537);
  deciq->add_task(task585);

  vector<IndexRange> I750_index = {active_, closed_};
  auto I750 = make_shared<Tensor>(I750_index);
  auto tensor586 = vector<shared_ptr<Tensor>>{I745, t2, I750};
  auto task586 = make_shared<Task586>(tensor586, cindex);
  task583->add_dep(task586);
  task586->add_dep(task537);
  deciq->add_task(task586);

  auto tensor587 = vector<shared_ptr<Tensor>>{I750, t2, f1_};
  auto task587 = make_shared<Task587>(tensor587, cindex);
  task586->add_dep(task587);
  task587->add_dep(task537);
  deciq->add_task(task587);

  vector<IndexRange> I788_index = {active_, virt_, closed_, active_};
  auto I788 = make_shared<Tensor>(I788_index);
  auto tensor588 = vector<shared_ptr<Tensor>>{I745, t2, I788};
  auto task588 = make_shared<Task588>(tensor588, cindex);
  task583->add_dep(task588);
  task588->add_dep(task537);
  deciq->add_task(task588);

  auto tensor589 = vector<shared_ptr<Tensor>>{I788, f1_, t2};
  auto task589 = make_shared<Task589>(tensor589, cindex);
  task588->add_dep(task589);
  task589->add_dep(task537);
  deciq->add_task(task589);

  vector<IndexRange> I792_index = {active_, closed_, virt_, active_};
  auto I792 = make_shared<Tensor>(I792_index);
  auto tensor590 = vector<shared_ptr<Tensor>>{I745, t2, I792};
  auto task590 = make_shared<Task590>(tensor590, cindex);
  task583->add_dep(task590);
  task590->add_dep(task537);
  deciq->add_task(task590);

  auto tensor591 = vector<shared_ptr<Tensor>>{I792, f1_, t2};
  auto task591 = make_shared<Task591>(tensor591, cindex);
  task590->add_dep(task591);
  task591->add_dep(task537);
  deciq->add_task(task591);

  vector<IndexRange> I796_index = {active_, virt_, closed_, active_};
  auto I796 = make_shared<Tensor>(I796_index);
  auto tensor592 = vector<shared_ptr<Tensor>>{I745, t2, I796};
  auto task592 = make_shared<Task592>(tensor592, cindex);
  task583->add_dep(task592);
  task592->add_dep(task537);
  deciq->add_task(task592);

  auto tensor593 = vector<shared_ptr<Tensor>>{I796, f1_, t2};
  auto task593 = make_shared<Task593>(tensor593, cindex);
  task592->add_dep(task593);
  task593->add_dep(task537);
  deciq->add_task(task593);

  auto tensor594 = vector<shared_ptr<Tensor>>{I745, v2_, t2};
  auto task594 = make_shared<Task594>(tensor594, cindex);
  task583->add_dep(task594);
  task594->add_dep(task537);
  deciq->add_task(task594);

  auto tensor595 = vector<shared_ptr<Tensor>>{I745, h1_, t2};
  auto task595 = make_shared<Task595>(tensor595, cindex);
  task583->add_dep(task595);
  task595->add_dep(task537);
  deciq->add_task(task595);

  vector<IndexRange> I753_index = {active_, active_};
  auto I753 = make_shared<Tensor>(I753_index);
  auto tensor596 = vector<shared_ptr<Tensor>>{I698, Gamma262_(), I753};
  auto task596 = make_shared<Task596>(tensor596, cindex);
  task538->add_dep(task596);
  task596->add_dep(task537);
  deciq->add_task(task596);

  auto tensor597 = vector<shared_ptr<Tensor>>{I753, t2};
  auto task597 = make_shared<Task597>(tensor597, cindex);
  task596->add_dep(task597);
  task597->add_dep(task537);
  deciq->add_task(task597);

  auto tensor598 = vector<shared_ptr<Tensor>>{I753, t2};
  auto task598 = make_shared<Task598>(tensor598, cindex);
  task596->add_dep(task598);
  task598->add_dep(task537);
  deciq->add_task(task598);

  vector<IndexRange> I759_index = {active_, active_};
  auto I759 = make_shared<Tensor>(I759_index);
  auto tensor599 = vector<shared_ptr<Tensor>>{I698, Gamma264_(), I759};
  auto task599 = make_shared<Task599>(tensor599, cindex);
  task538->add_dep(task599);
  task599->add_dep(task537);
  deciq->add_task(task599);

  vector<IndexRange> I760_index = {active_, closed_, virt_, closed_};
  auto I760 = make_shared<Tensor>(I760_index);
  auto tensor600 = vector<shared_ptr<Tensor>>{I759, t2, I760};
  auto task600 = make_shared<Task600>(tensor600, cindex);
  task599->add_dep(task600);
  task600->add_dep(task537);
  deciq->add_task(task600);

  auto tensor601 = vector<shared_ptr<Tensor>>{I760, f1_, t2};
  auto task601 = make_shared<Task601>(tensor601, cindex);
  task600->add_dep(task601);
  task601->add_dep(task537);
  deciq->add_task(task601);

  vector<IndexRange> I764_index = {active_, closed_, virt_, closed_};
  auto I764 = make_shared<Tensor>(I764_index);
  auto tensor602 = vector<shared_ptr<Tensor>>{I759, t2, I764};
  auto task602 = make_shared<Task602>(tensor602, cindex);
  task599->add_dep(task602);
  task602->add_dep(task537);
  deciq->add_task(task602);

  auto tensor603 = vector<shared_ptr<Tensor>>{I764, f1_, t2};
  auto task603 = make_shared<Task603>(tensor603, cindex);
  task602->add_dep(task603);
  task603->add_dep(task537);
  deciq->add_task(task603);

  vector<IndexRange> I768_index = {active_, virt_, closed_, closed_};
  auto I768 = make_shared<Tensor>(I768_index);
  auto tensor604 = vector<shared_ptr<Tensor>>{I759, t2, I768};
  auto task604 = make_shared<Task604>(tensor604, cindex);
  task599->add_dep(task604);
  task604->add_dep(task537);
  deciq->add_task(task604);

  auto tensor605 = vector<shared_ptr<Tensor>>{I768, f1_, t2};
  auto task605 = make_shared<Task605>(tensor605, cindex);
  task604->add_dep(task605);
  task605->add_dep(task537);
  deciq->add_task(task605);

  vector<IndexRange> I772_index = {active_, closed_, closed_, virt_};
  auto I772 = make_shared<Tensor>(I772_index);
  auto tensor606 = vector<shared_ptr<Tensor>>{I759, t2, I772};
  auto task606 = make_shared<Task606>(tensor606, cindex);
  task599->add_dep(task606);
  task606->add_dep(task537);
  deciq->add_task(task606);

  auto tensor607 = vector<shared_ptr<Tensor>>{I772, f1_, t2};
  auto task607 = make_shared<Task607>(tensor607, cindex);
  task606->add_dep(task607);
  task607->add_dep(task537);
  deciq->add_task(task607);

  vector<IndexRange> I776_index = {active_, virt_, closed_, closed_};
  auto I776 = make_shared<Tensor>(I776_index);
  auto tensor608 = vector<shared_ptr<Tensor>>{I759, t2, I776};
  auto task608 = make_shared<Task608>(tensor608, cindex);
  task599->add_dep(task608);
  task608->add_dep(task537);
  deciq->add_task(task608);

  auto tensor609 = vector<shared_ptr<Tensor>>{I776, f1_, t2};
  auto task609 = make_shared<Task609>(tensor609, cindex);
  task608->add_dep(task609);
  task609->add_dep(task537);
  deciq->add_task(task609);

  vector<IndexRange> I780_index = {active_, closed_, closed_, virt_};
  auto I780 = make_shared<Tensor>(I780_index);
  auto tensor610 = vector<shared_ptr<Tensor>>{I759, t2, I780};
  auto task610 = make_shared<Task610>(tensor610, cindex);
  task599->add_dep(task610);
  task610->add_dep(task537);
  deciq->add_task(task610);

  auto tensor611 = vector<shared_ptr<Tensor>>{I780, f1_, t2};
  auto task611 = make_shared<Task611>(tensor611, cindex);
  task610->add_dep(task611);
  task611->add_dep(task537);
  deciq->add_task(task611);

  vector<IndexRange> I800_index = {active_, virt_};
  auto I800 = make_shared<Tensor>(I800_index);
  auto tensor612 = vector<shared_ptr<Tensor>>{I759, f1_, I800};
  auto task612 = make_shared<Task612>(tensor612, cindex);
  task599->add_dep(task612);
  task612->add_dep(task537);
  deciq->add_task(task612);

  auto tensor613 = vector<shared_ptr<Tensor>>{I800, t2};
  auto task613 = make_shared<Task613>(tensor613, cindex);
  task612->add_dep(task613);
  task613->add_dep(task537);
  deciq->add_task(task613);

  auto tensor614 = vector<shared_ptr<Tensor>>{I800, t2};
  auto task614 = make_shared<Task614>(tensor614, cindex);
  task612->add_dep(task614);
  task614->add_dep(task537);
  deciq->add_task(task614);

  vector<IndexRange> I943_index = {virt_, active_};
  auto I943 = make_shared<Tensor>(I943_index);
  auto tensor615 = vector<shared_ptr<Tensor>>{I759, f1_, I943};
  auto task615 = make_shared<Task615>(tensor615, cindex);
  task599->add_dep(task615);
  task615->add_dep(task537);
  deciq->add_task(task615);

  auto tensor616 = vector<shared_ptr<Tensor>>{I943, t2};
  auto task616 = make_shared<Task616>(tensor616, cindex);
  task615->add_dep(task616);
  task616->add_dep(task537);
  deciq->add_task(task616);

  vector<IndexRange> I947_index = {virt_, active_};
  auto I947 = make_shared<Tensor>(I947_index);
  auto tensor617 = vector<shared_ptr<Tensor>>{I759, f1_, I947};
  auto task617 = make_shared<Task617>(tensor617, cindex);
  task599->add_dep(task617);
  task617->add_dep(task537);
  deciq->add_task(task617);

  auto tensor618 = vector<shared_ptr<Tensor>>{I947, t2};
  auto task618 = make_shared<Task618>(tensor618, cindex);
  task617->add_dep(task618);
  task618->add_dep(task537);
  deciq->add_task(task618);

  auto tensor619 = vector<shared_ptr<Tensor>>{I759, t2};
  auto task619 = make_shared<Task619>(tensor619, cindex, this->e0_);
  task599->add_dep(task619);
  task619->add_dep(task537);
  deciq->add_task(task619);

  auto tensor620 = vector<shared_ptr<Tensor>>{I759, t2};
  auto task620 = make_shared<Task620>(tensor620, cindex, this->e0_);
  task599->add_dep(task620);
  task620->add_dep(task537);
  deciq->add_task(task620);

  auto tensor621 = vector<shared_ptr<Tensor>>{I759, v2_, t2};
  auto task621 = make_shared<Task621>(tensor621, cindex);
  task599->add_dep(task621);
  task621->add_dep(task537);
  deciq->add_task(task621);

  auto tensor622 = vector<shared_ptr<Tensor>>{I759, v2_, t2};
  auto task622 = make_shared<Task622>(tensor622, cindex);
  task599->add_dep(task622);
  task622->add_dep(task537);
  deciq->add_task(task622);

  auto tensor623 = vector<shared_ptr<Tensor>>{I759, v2_, t2};
  auto task623 = make_shared<Task623>(tensor623, cindex);
  task599->add_dep(task623);
  task623->add_dep(task537);
  deciq->add_task(task623);

  auto tensor624 = vector<shared_ptr<Tensor>>{I759, v2_, t2};
  auto task624 = make_shared<Task624>(tensor624, cindex);
  task599->add_dep(task624);
  task624->add_dep(task537);
  deciq->add_task(task624);

  vector<IndexRange> I783_index = {active_, active_, active_, active_};
  auto I783 = make_shared<Tensor>(I783_index);
  auto tensor625 = vector<shared_ptr<Tensor>>{I698, Gamma270_(), I783};
  auto task625 = make_shared<Task625>(tensor625, cindex);
  task538->add_dep(task625);
  task625->add_dep(task537);
  deciq->add_task(task625);

  vector<IndexRange> I784_index = {active_, closed_, virt_, active_};
  auto I784 = make_shared<Tensor>(I784_index);
  auto tensor626 = vector<shared_ptr<Tensor>>{I783, t2, I784};
  auto task626 = make_shared<Task626>(tensor626, cindex);
  task625->add_dep(task626);
  task626->add_dep(task537);
  deciq->add_task(task626);

  auto tensor627 = vector<shared_ptr<Tensor>>{I784, f1_, t2};
  auto task627 = make_shared<Task627>(tensor627, cindex);
  task626->add_dep(task627);
  task627->add_dep(task537);
  deciq->add_task(task627);

  auto tensor628 = vector<shared_ptr<Tensor>>{I783, v2_, t2};
  auto task628 = make_shared<Task628>(tensor628, cindex);
  task625->add_dep(task628);
  task628->add_dep(task537);
  deciq->add_task(task628);

  vector<IndexRange> I807_index = {active_, active_, active_, active_, active_, active_};
  auto I807 = make_shared<Tensor>(I807_index);
  auto tensor629 = vector<shared_ptr<Tensor>>{I698, Gamma276_(), I807};
  auto task629 = make_shared<Task629>(tensor629, cindex);
  task538->add_dep(task629);
  task629->add_dep(task537);
  deciq->add_task(task629);

  vector<IndexRange> I808_index = {active_, closed_, active_, active_};
  auto I808 = make_shared<Tensor>(I808_index);
  auto tensor630 = vector<shared_ptr<Tensor>>{I807, t2, I808};
  auto task630 = make_shared<Task630>(tensor630, cindex);
  task629->add_dep(task630);
  task630->add_dep(task537);
  deciq->add_task(task630);

  auto tensor631 = vector<shared_ptr<Tensor>>{I808, f1_, t2};
  auto task631 = make_shared<Task631>(tensor631, cindex);
  task630->add_dep(task631);
  task631->add_dep(task537);
  deciq->add_task(task631);

  vector<IndexRange> I811_index = {active_, active_, active_, active_};
  auto I811 = make_shared<Tensor>(I811_index);
  auto tensor632 = vector<shared_ptr<Tensor>>{I698, Gamma277_(), I811};
  auto task632 = make_shared<Task632>(tensor632, cindex);
  task538->add_dep(task632);
  task632->add_dep(task537);
  deciq->add_task(task632);

  vector<IndexRange> I812_index = {active_, virt_, closed_, active_};
  auto I812 = make_shared<Tensor>(I812_index);
  auto tensor633 = vector<shared_ptr<Tensor>>{I811, t2, I812};
  auto task633 = make_shared<Task633>(tensor633, cindex);
  task632->add_dep(task633);
  task633->add_dep(task537);
  deciq->add_task(task633);

  auto tensor634 = vector<shared_ptr<Tensor>>{I812, t2, f1_};
  auto task634 = make_shared<Task634>(tensor634, cindex);
  task633->add_dep(task634);
  task634->add_dep(task537);
  deciq->add_task(task634);

  auto tensor635 = vector<shared_ptr<Tensor>>{I811, v2_, t2};
  auto task635 = make_shared<Task635>(tensor635, cindex);
  task632->add_dep(task635);
  task635->add_dep(task537);
  deciq->add_task(task635);

  vector<IndexRange> I819_index = {active_, active_, active_, active_};
  auto I819 = make_shared<Tensor>(I819_index);
  auto tensor636 = vector<shared_ptr<Tensor>>{I698, Gamma279_(), I819};
  auto task636 = make_shared<Task636>(tensor636, cindex);
  task538->add_dep(task636);
  task636->add_dep(task537);
  deciq->add_task(task636);

  auto tensor637 = vector<shared_ptr<Tensor>>{I819, t2};
  auto task637 = make_shared<Task637>(tensor637, cindex);
  task636->add_dep(task637);
  task637->add_dep(task537);
  deciq->add_task(task637);

  vector<IndexRange> I822_index = {active_, active_, active_, active_};
  auto I822 = make_shared<Tensor>(I822_index);
  auto tensor638 = vector<shared_ptr<Tensor>>{I698, Gamma280_(), I822};
  auto task638 = make_shared<Task638>(tensor638, cindex);
  task538->add_dep(task638);
  task638->add_dep(task537);
  deciq->add_task(task638);

  vector<IndexRange> I823_index = {active_, virt_, active_, closed_};
  auto I823 = make_shared<Tensor>(I823_index);
  auto tensor639 = vector<shared_ptr<Tensor>>{I822, t2, I823};
  auto task639 = make_shared<Task639>(tensor639, cindex);
  task638->add_dep(task639);
  task639->add_dep(task537);
  deciq->add_task(task639);

  auto tensor640 = vector<shared_ptr<Tensor>>{I823, f1_, t2};
  auto task640 = make_shared<Task640>(tensor640, cindex);
  task639->add_dep(task640);
  task640->add_dep(task537);
  deciq->add_task(task640);

  vector<IndexRange> I827_index = {active_, closed_, active_, virt_};
  auto I827 = make_shared<Tensor>(I827_index);
  auto tensor641 = vector<shared_ptr<Tensor>>{I822, t2, I827};
  auto task641 = make_shared<Task641>(tensor641, cindex);
  task638->add_dep(task641);
  task641->add_dep(task537);
  deciq->add_task(task641);

  auto tensor642 = vector<shared_ptr<Tensor>>{I827, f1_, t2};
  auto task642 = make_shared<Task642>(tensor642, cindex);
  task641->add_dep(task642);
  task642->add_dep(task537);
  deciq->add_task(task642);

  vector<IndexRange> I858_index = {active_, active_, virt_, closed_};
  auto I858 = make_shared<Tensor>(I858_index);
  auto tensor643 = vector<shared_ptr<Tensor>>{I822, t2, I858};
  auto task643 = make_shared<Task643>(tensor643, cindex);
  task638->add_dep(task643);
  task643->add_dep(task537);
  deciq->add_task(task643);

  auto tensor644 = vector<shared_ptr<Tensor>>{I858, t2, f1_};
  auto task644 = make_shared<Task644>(tensor644, cindex);
  task643->add_dep(task644);
  task644->add_dep(task537);
  deciq->add_task(task644);

  vector<IndexRange> I985_index = {closed_, virt_, active_, active_};
  auto I985 = make_shared<Tensor>(I985_index);
  auto tensor645 = vector<shared_ptr<Tensor>>{I822, t2, I985};
  auto task645 = make_shared<Task645>(tensor645, cindex);
  task638->add_dep(task645);
  task645->add_dep(task537);
  deciq->add_task(task645);

  auto tensor646 = vector<shared_ptr<Tensor>>{I985, t2};
  auto task646 = make_shared<Task646>(tensor646, cindex, this->e0_);
  task645->add_dep(task646);
  task646->add_dep(task537);
  deciq->add_task(task646);

  auto tensor647 = vector<shared_ptr<Tensor>>{I985, f1_, t2};
  auto task647 = make_shared<Task647>(tensor647, cindex);
  task645->add_dep(task647);
  task647->add_dep(task537);
  deciq->add_task(task647);

  auto tensor648 = vector<shared_ptr<Tensor>>{I822, v2_, t2};
  auto task648 = make_shared<Task648>(tensor648, cindex);
  task638->add_dep(task648);
  task648->add_dep(task537);
  deciq->add_task(task648);

  auto tensor649 = vector<shared_ptr<Tensor>>{I822, v2_, t2};
  auto task649 = make_shared<Task649>(tensor649, cindex);
  task638->add_dep(task649);
  task649->add_dep(task537);
  deciq->add_task(task649);

  vector<IndexRange> I830_index = {active_, active_, active_, active_};
  auto I830 = make_shared<Tensor>(I830_index);
  auto tensor650 = vector<shared_ptr<Tensor>>{I698, Gamma282_(), I830};
  auto task650 = make_shared<Task650>(tensor650, cindex);
  task538->add_dep(task650);
  task650->add_dep(task537);
  deciq->add_task(task650);

  auto tensor651 = vector<shared_ptr<Tensor>>{I830, t2};
  auto task651 = make_shared<Task651>(tensor651, cindex);
  task650->add_dep(task651);
  task651->add_dep(task537);
  deciq->add_task(task651);

  auto tensor652 = vector<shared_ptr<Tensor>>{I830, t2};
  auto task652 = make_shared<Task652>(tensor652, cindex);
  task650->add_dep(task652);
  task652->add_dep(task537);
  deciq->add_task(task652);

  auto tensor653 = vector<shared_ptr<Tensor>>{I830, t2};
  auto task653 = make_shared<Task653>(tensor653, cindex);
  task650->add_dep(task653);
  task653->add_dep(task537);
  deciq->add_task(task653);

  vector<IndexRange> I833_index = {active_, active_, active_, active_};
  auto I833 = make_shared<Tensor>(I833_index);
  auto tensor654 = vector<shared_ptr<Tensor>>{I698, Gamma283_(), I833};
  auto task654 = make_shared<Task654>(tensor654, cindex);
  task538->add_dep(task654);
  task654->add_dep(task537);
  deciq->add_task(task654);

  vector<IndexRange> I834_index = {active_, virt_, active_, closed_};
  auto I834 = make_shared<Tensor>(I834_index);
  auto tensor655 = vector<shared_ptr<Tensor>>{I833, t2, I834};
  auto task655 = make_shared<Task655>(tensor655, cindex);
  task654->add_dep(task655);
  task655->add_dep(task537);
  deciq->add_task(task655);

  auto tensor656 = vector<shared_ptr<Tensor>>{I834, f1_, t2};
  auto task656 = make_shared<Task656>(tensor656, cindex);
  task655->add_dep(task656);
  task656->add_dep(task537);
  deciq->add_task(task656);

  vector<IndexRange> I838_index = {active_, closed_, active_, virt_};
  auto I838 = make_shared<Tensor>(I838_index);
  auto tensor657 = vector<shared_ptr<Tensor>>{I833, t2, I838};
  auto task657 = make_shared<Task657>(tensor657, cindex);
  task654->add_dep(task657);
  task657->add_dep(task537);
  deciq->add_task(task657);

  auto tensor658 = vector<shared_ptr<Tensor>>{I838, f1_, t2};
  auto task658 = make_shared<Task658>(tensor658, cindex);
  task657->add_dep(task658);
  task658->add_dep(task537);
  deciq->add_task(task658);

  auto tensor659 = vector<shared_ptr<Tensor>>{I838, f1_, t2};
  auto task659 = make_shared<Task659>(tensor659, cindex);
  task657->add_dep(task659);
  task659->add_dep(task537);
  deciq->add_task(task659);

  vector<IndexRange> I854_index = {active_, active_, closed_, virt_};
  auto I854 = make_shared<Tensor>(I854_index);
  auto tensor660 = vector<shared_ptr<Tensor>>{I833, t2, I854};
  auto task660 = make_shared<Task660>(tensor660, cindex);
  task654->add_dep(task660);
  task660->add_dep(task537);
  deciq->add_task(task660);

  auto tensor661 = vector<shared_ptr<Tensor>>{I854, t2, f1_};
  auto task661 = make_shared<Task661>(tensor661, cindex);
  task660->add_dep(task661);
  task661->add_dep(task537);
  deciq->add_task(task661);

  vector<IndexRange> I877_index = {active_, active_, virt_, closed_};
  auto I877 = make_shared<Tensor>(I877_index);
  auto tensor662 = vector<shared_ptr<Tensor>>{I833, t2, I877};
  auto task662 = make_shared<Task662>(tensor662, cindex);
  task654->add_dep(task662);
  task662->add_dep(task537);
  deciq->add_task(task662);

  auto tensor663 = vector<shared_ptr<Tensor>>{I877, f1_, t2};
  auto task663 = make_shared<Task663>(tensor663, cindex);
  task662->add_dep(task663);
  task663->add_dep(task537);
  deciq->add_task(task663);

  vector<IndexRange> I881_index = {active_, active_, closed_, virt_};
  auto I881 = make_shared<Tensor>(I881_index);
  auto tensor664 = vector<shared_ptr<Tensor>>{I833, t2, I881};
  auto task664 = make_shared<Task664>(tensor664, cindex);
  task654->add_dep(task664);
  task664->add_dep(task537);
  deciq->add_task(task664);

  auto tensor665 = vector<shared_ptr<Tensor>>{I881, f1_, t2};
  auto task665 = make_shared<Task665>(tensor665, cindex);
  task664->add_dep(task665);
  task665->add_dep(task537);
  deciq->add_task(task665);

  vector<IndexRange> I888_index = {active_, active_, virt_, closed_};
  auto I888 = make_shared<Tensor>(I888_index);
  auto tensor666 = vector<shared_ptr<Tensor>>{I833, t2, I888};
  auto task666 = make_shared<Task666>(tensor666, cindex);
  task654->add_dep(task666);
  task666->add_dep(task537);
  deciq->add_task(task666);

  auto tensor667 = vector<shared_ptr<Tensor>>{I888, f1_, t2};
  auto task667 = make_shared<Task667>(tensor667, cindex);
  task666->add_dep(task667);
  task667->add_dep(task537);
  deciq->add_task(task667);

  vector<IndexRange> I892_index = {active_, active_, closed_, virt_};
  auto I892 = make_shared<Tensor>(I892_index);
  auto tensor668 = vector<shared_ptr<Tensor>>{I833, t2, I892};
  auto task668 = make_shared<Task668>(tensor668, cindex);
  task654->add_dep(task668);
  task668->add_dep(task537);
  deciq->add_task(task668);

  auto tensor669 = vector<shared_ptr<Tensor>>{I892, f1_, t2};
  auto task669 = make_shared<Task669>(tensor669, cindex);
  task668->add_dep(task669);
  task669->add_dep(task537);
  deciq->add_task(task669);

  vector<IndexRange> I908_index = {active_, active_, closed_, virt_};
  auto I908 = make_shared<Tensor>(I908_index);
  auto tensor670 = vector<shared_ptr<Tensor>>{I833, t2, I908};
  auto task670 = make_shared<Task670>(tensor670, cindex);
  task654->add_dep(task670);
  task670->add_dep(task537);
  deciq->add_task(task670);

  auto tensor671 = vector<shared_ptr<Tensor>>{I908, t2, f1_};
  auto task671 = make_shared<Task671>(tensor671, cindex);
  task670->add_dep(task671);
  task671->add_dep(task537);
  deciq->add_task(task671);

  auto tensor672 = vector<shared_ptr<Tensor>>{I908, t2, f1_};
  auto task672 = make_shared<Task672>(tensor672, cindex);
  task670->add_dep(task672);
  task672->add_dep(task537);
  deciq->add_task(task672);

  vector<IndexRange> I981_index = {virt_, closed_, active_, active_};
  auto I981 = make_shared<Tensor>(I981_index);
  auto tensor673 = vector<shared_ptr<Tensor>>{I833, t2, I981};
  auto task673 = make_shared<Task673>(tensor673, cindex);
  task654->add_dep(task673);
  task673->add_dep(task537);
  deciq->add_task(task673);

  auto tensor674 = vector<shared_ptr<Tensor>>{I981, f1_, t2};
  auto task674 = make_shared<Task674>(tensor674, cindex);
  task673->add_dep(task674);
  task674->add_dep(task537);
  deciq->add_task(task674);

  vector<IndexRange> I993_index = {closed_, virt_, active_, active_};
  auto I993 = make_shared<Tensor>(I993_index);
  auto tensor675 = vector<shared_ptr<Tensor>>{I833, t2, I993};
  auto task675 = make_shared<Task675>(tensor675, cindex);
  task654->add_dep(task675);
  task675->add_dep(task537);
  deciq->add_task(task675);

  auto tensor676 = vector<shared_ptr<Tensor>>{I993, t2};
  auto task676 = make_shared<Task676>(tensor676, cindex, this->e0_);
  task675->add_dep(task676);
  task676->add_dep(task537);
  deciq->add_task(task676);

  auto tensor677 = vector<shared_ptr<Tensor>>{I993, f1_, t2};
  auto task677 = make_shared<Task677>(tensor677, cindex);
  task675->add_dep(task677);
  task677->add_dep(task537);
  deciq->add_task(task677);

  auto tensor678 = vector<shared_ptr<Tensor>>{I833, t2};
  auto task678 = make_shared<Task678>(tensor678, cindex, this->e0_);
  task654->add_dep(task678);
  task678->add_dep(task537);
  deciq->add_task(task678);

  auto tensor679 = vector<shared_ptr<Tensor>>{I833, t2};
  auto task679 = make_shared<Task679>(tensor679, cindex, this->e0_);
  task654->add_dep(task679);
  task679->add_dep(task537);
  deciq->add_task(task679);

  auto tensor680 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task680 = make_shared<Task680>(tensor680, cindex);
  task654->add_dep(task680);
  task680->add_dep(task537);
  deciq->add_task(task680);

  auto tensor681 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task681 = make_shared<Task681>(tensor681, cindex);
  task654->add_dep(task681);
  task681->add_dep(task537);
  deciq->add_task(task681);

  auto tensor682 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task682 = make_shared<Task682>(tensor682, cindex);
  task654->add_dep(task682);
  task682->add_dep(task537);
  deciq->add_task(task682);

  auto tensor683 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task683 = make_shared<Task683>(tensor683, cindex);
  task654->add_dep(task683);
  task683->add_dep(task537);
  deciq->add_task(task683);

  auto tensor684 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task684 = make_shared<Task684>(tensor684, cindex);
  task654->add_dep(task684);
  task684->add_dep(task537);
  deciq->add_task(task684);

  auto tensor685 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task685 = make_shared<Task685>(tensor685, cindex);
  task654->add_dep(task685);
  task685->add_dep(task537);
  deciq->add_task(task685);

  auto tensor686 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task686 = make_shared<Task686>(tensor686, cindex);
  task654->add_dep(task686);
  task686->add_dep(task537);
  deciq->add_task(task686);

  auto tensor687 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task687 = make_shared<Task687>(tensor687, cindex);
  task654->add_dep(task687);
  task687->add_dep(task537);
  deciq->add_task(task687);

  auto tensor688 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task688 = make_shared<Task688>(tensor688, cindex);
  task654->add_dep(task688);
  task688->add_dep(task537);
  deciq->add_task(task688);

  auto tensor689 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task689 = make_shared<Task689>(tensor689, cindex);
  task654->add_dep(task689);
  task689->add_dep(task537);
  deciq->add_task(task689);

  vector<IndexRange> I841_index = {active_, active_, active_, active_, active_, active_};
  auto I841 = make_shared<Tensor>(I841_index);
  auto tensor690 = vector<shared_ptr<Tensor>>{I698, Gamma285_(), I841};
  auto task690 = make_shared<Task690>(tensor690, cindex);
  task538->add_dep(task690);
  task690->add_dep(task537);
  deciq->add_task(task690);

  vector<IndexRange> I842_index = {active_, virt_, active_, active_};
  auto I842 = make_shared<Tensor>(I842_index);
  auto tensor691 = vector<shared_ptr<Tensor>>{I841, t2, I842};
  auto task691 = make_shared<Task691>(tensor691, cindex);
  task690->add_dep(task691);
  task691->add_dep(task537);
  deciq->add_task(task691);

  auto tensor692 = vector<shared_ptr<Tensor>>{I842, f1_, t2};
  auto task692 = make_shared<Task692>(tensor692, cindex);
  task691->add_dep(task692);
  task692->add_dep(task537);
  deciq->add_task(task692);

  vector<IndexRange> I845_index = {active_, active_};
  auto I845 = make_shared<Tensor>(I845_index);
  auto tensor693 = vector<shared_ptr<Tensor>>{I698, Gamma286_(), I845};
  auto task693 = make_shared<Task693>(tensor693, cindex);
  task538->add_dep(task693);
  task693->add_dep(task537);
  deciq->add_task(task693);

  vector<IndexRange> I846_index = {closed_, virt_};
  auto I846 = make_shared<Tensor>(I846_index);
  auto tensor694 = vector<shared_ptr<Tensor>>{I845, t2, I846};
  auto task694 = make_shared<Task694>(tensor694, cindex);
  task693->add_dep(task694);
  task694->add_dep(task537);
  deciq->add_task(task694);

  auto tensor695 = vector<shared_ptr<Tensor>>{I846, t2, f1_};
  auto task695 = make_shared<Task695>(tensor695, cindex);
  task694->add_dep(task695);
  task695->add_dep(task537);
  deciq->add_task(task695);

  auto tensor696 = vector<shared_ptr<Tensor>>{I846, t2, f1_};
  auto task696 = make_shared<Task696>(tensor696, cindex);
  task694->add_dep(task696);
  task696->add_dep(task537);
  deciq->add_task(task696);

  vector<IndexRange> I900_index = {closed_, virt_};
  auto I900 = make_shared<Tensor>(I900_index);
  auto tensor697 = vector<shared_ptr<Tensor>>{I845, t2, I900};
  auto task697 = make_shared<Task697>(tensor697, cindex);
  task693->add_dep(task697);
  task697->add_dep(task537);
  deciq->add_task(task697);

  auto tensor698 = vector<shared_ptr<Tensor>>{I900, t2, f1_};
  auto task698 = make_shared<Task698>(tensor698, cindex);
  task697->add_dep(task698);
  task698->add_dep(task537);
  deciq->add_task(task698);

  auto tensor699 = vector<shared_ptr<Tensor>>{I900, t2, f1_};
  auto task699 = make_shared<Task699>(tensor699, cindex);
  task697->add_dep(task699);
  task699->add_dep(task537);
  deciq->add_task(task699);

  vector<IndexRange> I951_index = {virt_, closed_};
  auto I951 = make_shared<Tensor>(I951_index);
  auto tensor700 = vector<shared_ptr<Tensor>>{I845, t2, I951};
  auto task700 = make_shared<Task700>(tensor700, cindex);
  task693->add_dep(task700);
  task700->add_dep(task537);
  deciq->add_task(task700);

  auto tensor701 = vector<shared_ptr<Tensor>>{I951, f1_, t2};
  auto task701 = make_shared<Task701>(tensor701, cindex);
  task700->add_dep(task701);
  task701->add_dep(task537);
  deciq->add_task(task701);

  vector<IndexRange> I955_index = {virt_, closed_};
  auto I955 = make_shared<Tensor>(I955_index);
  auto tensor702 = vector<shared_ptr<Tensor>>{I845, t2, I955};
  auto task702 = make_shared<Task702>(tensor702, cindex);
  task693->add_dep(task702);
  task702->add_dep(task537);
  deciq->add_task(task702);

  auto tensor703 = vector<shared_ptr<Tensor>>{I955, f1_, t2};
  auto task703 = make_shared<Task703>(tensor703, cindex);
  task702->add_dep(task703);
  task703->add_dep(task537);
  deciq->add_task(task703);

  vector<IndexRange> I959_index = {virt_, closed_};
  auto I959 = make_shared<Tensor>(I959_index);
  auto tensor704 = vector<shared_ptr<Tensor>>{I845, t2, I959};
  auto task704 = make_shared<Task704>(tensor704, cindex);
  task693->add_dep(task704);
  task704->add_dep(task537);
  deciq->add_task(task704);

  auto tensor705 = vector<shared_ptr<Tensor>>{I959, f1_, t2};
  auto task705 = make_shared<Task705>(tensor705, cindex);
  task704->add_dep(task705);
  task705->add_dep(task537);
  deciq->add_task(task705);

  vector<IndexRange> I963_index = {virt_, closed_};
  auto I963 = make_shared<Tensor>(I963_index);
  auto tensor706 = vector<shared_ptr<Tensor>>{I845, t2, I963};
  auto task706 = make_shared<Task706>(tensor706, cindex);
  task693->add_dep(task706);
  task706->add_dep(task537);
  deciq->add_task(task706);

  auto tensor707 = vector<shared_ptr<Tensor>>{I963, f1_, t2};
  auto task707 = make_shared<Task707>(tensor707, cindex);
  task706->add_dep(task707);
  task707->add_dep(task537);
  deciq->add_task(task707);

  vector<IndexRange> I973_index = {closed_, active_};
  auto I973 = make_shared<Tensor>(I973_index);
  auto tensor708 = vector<shared_ptr<Tensor>>{I845, f1_, I973};
  auto task708 = make_shared<Task708>(tensor708, cindex);
  task693->add_dep(task708);
  task708->add_dep(task537);
  deciq->add_task(task708);

  auto tensor709 = vector<shared_ptr<Tensor>>{I973, t2};
  auto task709 = make_shared<Task709>(tensor709, cindex);
  task708->add_dep(task709);
  task709->add_dep(task537);
  deciq->add_task(task709);

  auto tensor710 = vector<shared_ptr<Tensor>>{I973, t2};
  auto task710 = make_shared<Task710>(tensor710, cindex);
  task708->add_dep(task710);
  task710->add_dep(task537);
  deciq->add_task(task710);

  vector<IndexRange> I1005_index = {active_, closed_};
  auto I1005 = make_shared<Tensor>(I1005_index);
  auto tensor711 = vector<shared_ptr<Tensor>>{I845, f1_, I1005};
  auto task711 = make_shared<Task711>(tensor711, cindex);
  task693->add_dep(task711);
  task711->add_dep(task537);
  deciq->add_task(task711);

  auto tensor712 = vector<shared_ptr<Tensor>>{I1005, t2};
  auto task712 = make_shared<Task712>(tensor712, cindex);
  task711->add_dep(task712);
  task712->add_dep(task537);
  deciq->add_task(task712);

  auto tensor713 = vector<shared_ptr<Tensor>>{I1005, t2};
  auto task713 = make_shared<Task713>(tensor713, cindex);
  task711->add_dep(task713);
  task713->add_dep(task537);
  deciq->add_task(task713);

  vector<IndexRange> I1019_index = {virt_, virt_, active_, closed_};
  auto I1019 = make_shared<Tensor>(I1019_index);
  auto tensor714 = vector<shared_ptr<Tensor>>{I845, t2, I1019};
  auto task714 = make_shared<Task714>(tensor714, cindex);
  task693->add_dep(task714);
  task714->add_dep(task537);
  deciq->add_task(task714);

  auto tensor715 = vector<shared_ptr<Tensor>>{I1019, f1_, t2};
  auto task715 = make_shared<Task715>(tensor715, cindex);
  task714->add_dep(task715);
  task715->add_dep(task537);
  deciq->add_task(task715);

  vector<IndexRange> I1023_index = {virt_, virt_, active_, closed_};
  auto I1023 = make_shared<Tensor>(I1023_index);
  auto tensor716 = vector<shared_ptr<Tensor>>{I845, t2, I1023};
  auto task716 = make_shared<Task716>(tensor716, cindex);
  task693->add_dep(task716);
  task716->add_dep(task537);
  deciq->add_task(task716);

  auto tensor717 = vector<shared_ptr<Tensor>>{I1023, f1_, t2};
  auto task717 = make_shared<Task717>(tensor717, cindex);
  task716->add_dep(task717);
  task717->add_dep(task537);
  deciq->add_task(task717);

  vector<IndexRange> I1027_index = {virt_, closed_, active_, virt_};
  auto I1027 = make_shared<Tensor>(I1027_index);
  auto tensor718 = vector<shared_ptr<Tensor>>{I845, t2, I1027};
  auto task718 = make_shared<Task718>(tensor718, cindex);
  task693->add_dep(task718);
  task718->add_dep(task537);
  deciq->add_task(task718);

  auto tensor719 = vector<shared_ptr<Tensor>>{I1027, f1_, t2};
  auto task719 = make_shared<Task719>(tensor719, cindex);
  task718->add_dep(task719);
  task719->add_dep(task537);
  deciq->add_task(task719);

  vector<IndexRange> I1031_index = {virt_, closed_, active_, virt_};
  auto I1031 = make_shared<Tensor>(I1031_index);
  auto tensor720 = vector<shared_ptr<Tensor>>{I845, t2, I1031};
  auto task720 = make_shared<Task720>(tensor720, cindex);
  task693->add_dep(task720);
  task720->add_dep(task537);
  deciq->add_task(task720);

  auto tensor721 = vector<shared_ptr<Tensor>>{I1031, f1_, t2};
  auto task721 = make_shared<Task721>(tensor721, cindex);
  task720->add_dep(task721);
  task721->add_dep(task537);
  deciq->add_task(task721);

  vector<IndexRange> I1035_index = {closed_, virt_, active_, virt_};
  auto I1035 = make_shared<Tensor>(I1035_index);
  auto tensor722 = vector<shared_ptr<Tensor>>{I845, t2, I1035};
  auto task722 = make_shared<Task722>(tensor722, cindex);
  task693->add_dep(task722);
  task722->add_dep(task537);
  deciq->add_task(task722);

  auto tensor723 = vector<shared_ptr<Tensor>>{I1035, f1_, t2};
  auto task723 = make_shared<Task723>(tensor723, cindex);
  task722->add_dep(task723);
  task723->add_dep(task537);
  deciq->add_task(task723);

  vector<IndexRange> I1039_index = {closed_, virt_, active_, virt_};
  auto I1039 = make_shared<Tensor>(I1039_index);
  auto tensor724 = vector<shared_ptr<Tensor>>{I845, t2, I1039};
  auto task724 = make_shared<Task724>(tensor724, cindex);
  task693->add_dep(task724);
  task724->add_dep(task537);
  deciq->add_task(task724);

  auto tensor725 = vector<shared_ptr<Tensor>>{I1039, f1_, t2};
  auto task725 = make_shared<Task725>(tensor725, cindex);
  task724->add_dep(task725);
  task725->add_dep(task537);
  deciq->add_task(task725);

  auto tensor726 = vector<shared_ptr<Tensor>>{I845, t2};
  auto task726 = make_shared<Task726>(tensor726, cindex, this->e0_);
  task693->add_dep(task726);
  task726->add_dep(task537);
  deciq->add_task(task726);

  auto tensor727 = vector<shared_ptr<Tensor>>{I845, t2};
  auto task727 = make_shared<Task727>(tensor727, cindex, this->e0_);
  task693->add_dep(task727);
  task727->add_dep(task537);
  deciq->add_task(task727);

  auto tensor728 = vector<shared_ptr<Tensor>>{I845, v2_, t2};
  auto task728 = make_shared<Task728>(tensor728, cindex);
  task693->add_dep(task728);
  task728->add_dep(task537);
  deciq->add_task(task728);

  auto tensor729 = vector<shared_ptr<Tensor>>{I845, v2_, t2};
  auto task729 = make_shared<Task729>(tensor729, cindex);
  task693->add_dep(task729);
  task729->add_dep(task537);
  deciq->add_task(task729);

  auto tensor730 = vector<shared_ptr<Tensor>>{I845, v2_, t2};
  auto task730 = make_shared<Task730>(tensor730, cindex);
  task693->add_dep(task730);
  task730->add_dep(task537);
  deciq->add_task(task730);

  auto tensor731 = vector<shared_ptr<Tensor>>{I845, v2_, t2};
  auto task731 = make_shared<Task731>(tensor731, cindex);
  task693->add_dep(task731);
  task731->add_dep(task537);
  deciq->add_task(task731);

  vector<IndexRange> I1209_index = {active_, closed_, virt_, active_};
  auto I1209 = make_shared<Tensor>(I1209_index);
  auto tensor732 = vector<shared_ptr<Tensor>>{I845, h1_, I1209};
  auto task732 = make_shared<Task732>(tensor732, cindex);
  task693->add_dep(task732);
  task732->add_dep(task537);
  deciq->add_task(task732);

  auto tensor733 = vector<shared_ptr<Tensor>>{I1209, t2};
  auto task733 = make_shared<Task733>(tensor733, cindex);
  task732->add_dep(task733);
  task733->add_dep(task537);
  deciq->add_task(task733);

  vector<IndexRange> I1212_index = {active_, active_, virt_, closed_};
  auto I1212 = make_shared<Tensor>(I1212_index);
  auto tensor734 = vector<shared_ptr<Tensor>>{I845, h1_, I1212};
  auto task734 = make_shared<Task734>(tensor734, cindex);
  task693->add_dep(task734);
  task734->add_dep(task537);
  deciq->add_task(task734);

  auto tensor735 = vector<shared_ptr<Tensor>>{I1212, t2};
  auto task735 = make_shared<Task735>(tensor735, cindex);
  task734->add_dep(task735);
  task735->add_dep(task537);
  deciq->add_task(task735);

  vector<IndexRange> I895_index = {active_, active_, active_, active_, active_, active_};
  auto I895 = make_shared<Tensor>(I895_index);
  auto tensor736 = vector<shared_ptr<Tensor>>{I698, Gamma299_(), I895};
  auto task736 = make_shared<Task736>(tensor736, cindex);
  task538->add_dep(task736);
  task736->add_dep(task537);
  deciq->add_task(task736);

  vector<IndexRange> I896_index = {active_, active_, virt_, active_};
  auto I896 = make_shared<Tensor>(I896_index);
  auto tensor737 = vector<shared_ptr<Tensor>>{I895, t2, I896};
  auto task737 = make_shared<Task737>(tensor737, cindex);
  task736->add_dep(task737);
  task737->add_dep(task537);
  deciq->add_task(task737);

  auto tensor738 = vector<shared_ptr<Tensor>>{I896, f1_, t2};
  auto task738 = make_shared<Task738>(tensor738, cindex);
  task737->add_dep(task738);
  task738->add_dep(task537);
  deciq->add_task(task738);

  auto tensor739 = vector<shared_ptr<Tensor>>{I895, v2_, t2};
  auto task739 = make_shared<Task739>(tensor739, cindex);
  task736->add_dep(task739);
  task739->add_dep(task537);
  deciq->add_task(task739);

  vector<IndexRange> I915_index = {active_, active_, active_, active_, active_, active_};
  auto I915 = make_shared<Tensor>(I915_index);
  auto tensor740 = vector<shared_ptr<Tensor>>{I698, Gamma304_(), I915};
  auto task740 = make_shared<Task740>(tensor740, cindex);
  task538->add_dep(task740);
  task740->add_dep(task537);
  deciq->add_task(task740);

  vector<IndexRange> I916_index = {active_, active_, virt_, active_};
  auto I916 = make_shared<Tensor>(I916_index);
  auto tensor741 = vector<shared_ptr<Tensor>>{I915, t2, I916};
  auto task741 = make_shared<Task741>(tensor741, cindex);
  task740->add_dep(task741);
  task741->add_dep(task537);
  deciq->add_task(task741);

  auto tensor742 = vector<shared_ptr<Tensor>>{I916, t2, f1_};
  auto task742 = make_shared<Task742>(tensor742, cindex);
  task741->add_dep(task742);
  task742->add_dep(task537);
  deciq->add_task(task742);

  vector<IndexRange> I919_index = {active_, active_, active_, active_, active_, active_};
  auto I919 = make_shared<Tensor>(I919_index);
  auto tensor743 = vector<shared_ptr<Tensor>>{I698, Gamma305_(), I919};
  auto task743 = make_shared<Task743>(tensor743, cindex);
  task538->add_dep(task743);
  task743->add_dep(task537);
  deciq->add_task(task743);

  vector<IndexRange> I920_index = {active_, virt_, active_, active_};
  auto I920 = make_shared<Tensor>(I920_index);
  auto tensor744 = vector<shared_ptr<Tensor>>{I919, t2, I920};
  auto task744 = make_shared<Task744>(tensor744, cindex);
  task743->add_dep(task744);
  task744->add_dep(task537);
  deciq->add_task(task744);

  auto tensor745 = vector<shared_ptr<Tensor>>{I920, t2, f1_};
  auto task745 = make_shared<Task745>(tensor745, cindex);
  task744->add_dep(task745);
  task745->add_dep(task537);
  deciq->add_task(task745);

  auto tensor746 = vector<shared_ptr<Tensor>>{I919, v2_, t2};
  auto task746 = make_shared<Task746>(tensor746, cindex);
  task743->add_dep(task746);
  task746->add_dep(task537);
  deciq->add_task(task746);

  vector<IndexRange> I923_index = {active_, active_, active_, active_, active_, active_};
  auto I923 = make_shared<Tensor>(I923_index);
  auto tensor747 = vector<shared_ptr<Tensor>>{I698, Gamma306_(), I923};
  auto task747 = make_shared<Task747>(tensor747, cindex);
  task538->add_dep(task747);
  task747->add_dep(task537);
  deciq->add_task(task747);

  auto tensor748 = vector<shared_ptr<Tensor>>{I923, t2};
  auto task748 = make_shared<Task748>(tensor748, cindex);
  task747->add_dep(task748);
  task748->add_dep(task537);
  deciq->add_task(task748);

  vector<IndexRange> I926_index = {active_, active_, active_, active_, active_, active_};
  auto I926 = make_shared<Tensor>(I926_index);
  auto tensor749 = vector<shared_ptr<Tensor>>{I698, Gamma307_(), I926};
  auto task749 = make_shared<Task749>(tensor749, cindex);
  task538->add_dep(task749);
  task749->add_dep(task537);
  deciq->add_task(task749);

  vector<IndexRange> I927_index = {active_, active_, active_, virt_};
  auto I927 = make_shared<Tensor>(I927_index);
  auto tensor750 = vector<shared_ptr<Tensor>>{I926, t2, I927};
  auto task750 = make_shared<Task750>(tensor750, cindex);
  task749->add_dep(task750);
  task750->add_dep(task537);
  deciq->add_task(task750);

  auto tensor751 = vector<shared_ptr<Tensor>>{I927, f1_, t2};
  auto task751 = make_shared<Task751>(tensor751, cindex);
  task750->add_dep(task751);
  task751->add_dep(task537);
  deciq->add_task(task751);

  vector<IndexRange> I939_index = {active_, active_, virt_, active_};
  auto I939 = make_shared<Tensor>(I939_index);
  auto tensor752 = vector<shared_ptr<Tensor>>{I926, t2, I939};
  auto task752 = make_shared<Task752>(tensor752, cindex);
  task749->add_dep(task752);
  task752->add_dep(task537);
  deciq->add_task(task752);

  auto tensor753 = vector<shared_ptr<Tensor>>{I939, t2, f1_};
  auto task753 = make_shared<Task753>(tensor753, cindex);
  task752->add_dep(task753);
  task753->add_dep(task537);
  deciq->add_task(task753);

  vector<IndexRange> I1047_index = {active_, virt_, active_, active_};
  auto I1047 = make_shared<Tensor>(I1047_index);
  auto tensor754 = vector<shared_ptr<Tensor>>{I926, t2, I1047};
  auto task754 = make_shared<Task754>(tensor754, cindex);
  task749->add_dep(task754);
  task754->add_dep(task537);
  deciq->add_task(task754);

  auto tensor755 = vector<shared_ptr<Tensor>>{I1047, t2};
  auto task755 = make_shared<Task755>(tensor755, cindex, this->e0_);
  task754->add_dep(task755);
  task755->add_dep(task537);
  deciq->add_task(task755);

  auto tensor756 = vector<shared_ptr<Tensor>>{I1047, f1_, t2};
  auto task756 = make_shared<Task756>(tensor756, cindex);
  task754->add_dep(task756);
  task756->add_dep(task537);
  deciq->add_task(task756);

  auto tensor757 = vector<shared_ptr<Tensor>>{I926, v2_, t2};
  auto task757 = make_shared<Task757>(tensor757, cindex);
  task749->add_dep(task757);
  task757->add_dep(task537);
  deciq->add_task(task757);

  auto tensor758 = vector<shared_ptr<Tensor>>{I926, v2_, t2};
  auto task758 = make_shared<Task758>(tensor758, cindex);
  task749->add_dep(task758);
  task758->add_dep(task537);
  deciq->add_task(task758);

  vector<IndexRange> I930_index = {active_, active_, active_, active_};
  auto I930 = make_shared<Tensor>(I930_index);
  auto tensor759 = vector<shared_ptr<Tensor>>{I698, Gamma308_(), I930};
  auto task759 = make_shared<Task759>(tensor759, cindex);
  task538->add_dep(task759);
  task759->add_dep(task537);
  deciq->add_task(task759);

  vector<IndexRange> I931_index = {active_, virt_};
  auto I931 = make_shared<Tensor>(I931_index);
  auto tensor760 = vector<shared_ptr<Tensor>>{I930, t2, I931};
  auto task760 = make_shared<Task760>(tensor760, cindex);
  task759->add_dep(task760);
  task760->add_dep(task537);
  deciq->add_task(task760);

  auto tensor761 = vector<shared_ptr<Tensor>>{I931, t2, f1_};
  auto task761 = make_shared<Task761>(tensor761, cindex);
  task760->add_dep(task761);
  task761->add_dep(task537);
  deciq->add_task(task761);

  auto tensor762 = vector<shared_ptr<Tensor>>{I931, t2, f1_};
  auto task762 = make_shared<Task762>(tensor762, cindex);
  task760->add_dep(task762);
  task762->add_dep(task537);
  deciq->add_task(task762);

  vector<IndexRange> I997_index = {virt_, active_};
  auto I997 = make_shared<Tensor>(I997_index);
  auto tensor763 = vector<shared_ptr<Tensor>>{I930, t2, I997};
  auto task763 = make_shared<Task763>(tensor763, cindex);
  task759->add_dep(task763);
  task763->add_dep(task537);
  deciq->add_task(task763);

  auto tensor764 = vector<shared_ptr<Tensor>>{I997, f1_, t2};
  auto task764 = make_shared<Task764>(tensor764, cindex);
  task763->add_dep(task764);
  task764->add_dep(task537);
  deciq->add_task(task764);

  vector<IndexRange> I1001_index = {virt_, active_};
  auto I1001 = make_shared<Tensor>(I1001_index);
  auto tensor765 = vector<shared_ptr<Tensor>>{I930, t2, I1001};
  auto task765 = make_shared<Task765>(tensor765, cindex);
  task759->add_dep(task765);
  task765->add_dep(task537);
  deciq->add_task(task765);

  auto tensor766 = vector<shared_ptr<Tensor>>{I1001, f1_, t2};
  auto task766 = make_shared<Task766>(tensor766, cindex);
  task765->add_dep(task766);
  task766->add_dep(task537);
  deciq->add_task(task766);

  vector<IndexRange> I1043_index = {virt_, virt_, active_, active_};
  auto I1043 = make_shared<Tensor>(I1043_index);
  auto tensor767 = vector<shared_ptr<Tensor>>{I930, t2, I1043};
  auto task767 = make_shared<Task767>(tensor767, cindex);
  task759->add_dep(task767);
  task767->add_dep(task537);
  deciq->add_task(task767);

  auto tensor768 = vector<shared_ptr<Tensor>>{I1043, f1_, t2};
  auto task768 = make_shared<Task768>(tensor768, cindex);
  task767->add_dep(task768);
  task768->add_dep(task537);
  deciq->add_task(task768);

  auto tensor769 = vector<shared_ptr<Tensor>>{I1043, f1_, t2};
  auto task769 = make_shared<Task769>(tensor769, cindex);
  task767->add_dep(task769);
  task769->add_dep(task537);
  deciq->add_task(task769);

  vector<IndexRange> I1051_index = {active_, active_, virt_, virt_};
  auto I1051 = make_shared<Tensor>(I1051_index);
  auto tensor770 = vector<shared_ptr<Tensor>>{I930, t2, I1051};
  auto task770 = make_shared<Task770>(tensor770, cindex);
  task759->add_dep(task770);
  task770->add_dep(task537);
  deciq->add_task(task770);

  auto tensor771 = vector<shared_ptr<Tensor>>{I1051, t2, f1_};
  auto task771 = make_shared<Task771>(tensor771, cindex);
  task770->add_dep(task771);
  task771->add_dep(task537);
  deciq->add_task(task771);

  auto tensor772 = vector<shared_ptr<Tensor>>{I930, t2};
  auto task772 = make_shared<Task772>(tensor772, cindex, this->e0_);
  task759->add_dep(task772);
  task772->add_dep(task537);
  deciq->add_task(task772);

  auto tensor773 = vector<shared_ptr<Tensor>>{I930, v2_, t2};
  auto task773 = make_shared<Task773>(tensor773, cindex);
  task759->add_dep(task773);
  task773->add_dep(task537);
  deciq->add_task(task773);

  auto tensor774 = vector<shared_ptr<Tensor>>{I930, v2_, t2};
  auto task774 = make_shared<Task774>(tensor774, cindex);
  task759->add_dep(task774);
  task774->add_dep(task537);
  deciq->add_task(task774);

  auto tensor775 = vector<shared_ptr<Tensor>>{I930, h1_, t2};
  auto task775 = make_shared<Task775>(tensor775, cindex);
  task759->add_dep(task775);
  task775->add_dep(task537);
  deciq->add_task(task775);

  auto tensor776 = vector<shared_ptr<Tensor>>{I930, h1_, t2};
  auto task776 = make_shared<Task776>(tensor776, cindex);
  task759->add_dep(task776);
  task776->add_dep(task537);
  deciq->add_task(task776);

  vector<IndexRange> I966_index;
  auto I966 = make_shared<Tensor>(I966_index);
  auto tensor777 = vector<shared_ptr<Tensor>>{I698, Gamma317_(), I966};
  auto task777 = make_shared<Task777>(tensor777, cindex);
  task538->add_dep(task777);
  task777->add_dep(task537);
  deciq->add_task(task777);

  auto tensor778 = vector<shared_ptr<Tensor>>{I966, t2};
  auto task778 = make_shared<Task778>(tensor778, cindex);
  task777->add_dep(task778);
  task778->add_dep(task537);
  deciq->add_task(task778);

  auto tensor779 = vector<shared_ptr<Tensor>>{I966, t2};
  auto task779 = make_shared<Task779>(tensor779, cindex);
  task777->add_dep(task779);
  task779->add_dep(task537);
  deciq->add_task(task779);

  vector<IndexRange> I1012_index = {active_, active_};
  auto I1012 = make_shared<Tensor>(I1012_index);
  auto tensor780 = vector<shared_ptr<Tensor>>{I698, Gamma329_(), I1012};
  auto task780 = make_shared<Task780>(tensor780, cindex);
  task538->add_dep(task780);
  task780->add_dep(task537);
  deciq->add_task(task780);

  auto tensor781 = vector<shared_ptr<Tensor>>{I1012, t2};
  auto task781 = make_shared<Task781>(tensor781, cindex);
  task780->add_dep(task781);
  task781->add_dep(task537);
  deciq->add_task(task781);

  auto tensor782 = vector<shared_ptr<Tensor>>{I1012, t2};
  auto task782 = make_shared<Task782>(tensor782, cindex);
  task780->add_dep(task782);
  task782->add_dep(task537);
  deciq->add_task(task782);

  vector<IndexRange> I1054_index = {active_, active_, active_, active_};
  auto I1054 = make_shared<Tensor>(I1054_index);
  auto tensor783 = vector<shared_ptr<Tensor>>{I698, Gamma340_(), I1054};
  auto task783 = make_shared<Task783>(tensor783, cindex);
  task538->add_dep(task783);
  task783->add_dep(task537);
  deciq->add_task(task783);

  auto tensor784 = vector<shared_ptr<Tensor>>{I1054, t2};
  auto task784 = make_shared<Task784>(tensor784, cindex);
  task783->add_dep(task784);
  task784->add_dep(task537);
  deciq->add_task(task784);

  vector<IndexRange> I1100_index = {active_, active_, active_, active_, active_, active_};
  auto I1100 = make_shared<Tensor>(I1100_index);
  auto tensor785 = vector<shared_ptr<Tensor>>{I698, Gamma355_(), I1100};
  auto task785 = make_shared<Task785>(tensor785, cindex);
  task538->add_dep(task785);
  task785->add_dep(task537);
  deciq->add_task(task785);

  auto tensor786 = vector<shared_ptr<Tensor>>{I1100, v2_, t2};
  auto task786 = make_shared<Task786>(tensor786, cindex);
  task785->add_dep(task786);
  task786->add_dep(task537);
  deciq->add_task(task786);

  vector<IndexRange> I1154_index = {active_, active_, active_, active_, active_, active_};
  auto I1154 = make_shared<Tensor>(I1154_index);
  auto tensor787 = vector<shared_ptr<Tensor>>{I698, Gamma373_(), I1154};
  auto task787 = make_shared<Task787>(tensor787, cindex);
  task538->add_dep(task787);
  task787->add_dep(task537);
  deciq->add_task(task787);

  auto tensor788 = vector<shared_ptr<Tensor>>{I1154, v2_, t2};
  auto task788 = make_shared<Task788>(tensor788, cindex);
  task787->add_dep(task788);
  task788->add_dep(task537);
  deciq->add_task(task788);

  return deciq;
}


#endif
