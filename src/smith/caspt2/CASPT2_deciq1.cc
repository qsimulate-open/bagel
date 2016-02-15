//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_deciq1.cc
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
#include <src/smith/caspt2/CASPT2_tasks12.h>
#include <src/smith/caspt2/CASPT2_tasks13.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_deciq(const bool reset, const bool diagonal) {
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deciq = make_shared<Queue>();
  auto tensor538 = vector<shared_ptr<Tensor>>{deci};
  auto task538 = make_shared<Task538>(tensor538, reset);
  deciq->add_task(task538);

  vector<IndexRange> I698_index = {ci_};
  auto I698 = make_shared<Tensor>(I698_index);
  auto tensor539 = vector<shared_ptr<Tensor>>{deci, I698};
  auto task539 = make_shared<Task539>(tensor539, cindex);
  task539->add_dep(task538);
  deciq->add_task(task539);

  make_deciq1(deciq, task538, task539, diagonal, I698);
  make_deciq2(deciq, task538, task539, diagonal, I698);
  make_deciq3(deciq, task538, task539, diagonal, I698);

  return deciq;
}


void CASPT2::CASPT2::make_deciq1(shared_ptr<Queue> deciq, shared_ptr<Task> task538, shared_ptr<Task> task539, const bool diagonal, shared_ptr<Tensor> I698) {
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  vector<IndexRange> I699_index = {active_, active_, active_, active_};
  auto I699 = make_shared<Tensor>(I699_index);
  auto tensor540 = vector<shared_ptr<Tensor>>{I698, Gamma248_(), I699};
  auto task540 = make_shared<Task540>(tensor540, cindex);
  task539->add_dep(task540);
  task540->add_dep(task538);
  deciq->add_task(task540);

  auto tensor541 = vector<shared_ptr<Tensor>>{I699, t2};
  auto task541 = make_shared<Task541>(tensor541, cindex);
  task540->add_dep(task541);
  task541->add_dep(task538);
  deciq->add_task(task541);

  vector<IndexRange> I702_index = {active_, active_, active_, active_};
  auto I702 = make_shared<Tensor>(I702_index);
  auto tensor542 = vector<shared_ptr<Tensor>>{I698, Gamma249_(), I702};
  auto task542 = make_shared<Task542>(tensor542, cindex);
  task539->add_dep(task542);
  task542->add_dep(task538);
  deciq->add_task(task542);

  vector<IndexRange> I703_index = {active_, active_, closed_, closed_};
  auto I703 = make_shared<Tensor>(I703_index);
  auto tensor543 = vector<shared_ptr<Tensor>>{I702, t2, I703};
  auto task543 = make_shared<Task543>(tensor543, cindex);
  task542->add_dep(task543);
  task543->add_dep(task538);
  deciq->add_task(task543);

  auto tensor544 = vector<shared_ptr<Tensor>>{I703, f1_, t2};
  auto task544 = make_shared<Task544>(tensor544, cindex);
  task543->add_dep(task544);
  task544->add_dep(task538);
  deciq->add_task(task544);

  auto tensor545 = vector<shared_ptr<Tensor>>{I702, t2};
  auto task545 = make_shared<Task545>(tensor545, cindex, this->e0_);
  task542->add_dep(task545);
  task545->add_dep(task538);
  deciq->add_task(task545);

  auto tensor546 = vector<shared_ptr<Tensor>>{I702, v2_, t2};
  auto task546 = make_shared<Task546>(tensor546, cindex);
  task542->add_dep(task546);
  task546->add_dep(task538);
  deciq->add_task(task546);

  auto tensor547 = vector<shared_ptr<Tensor>>{I702, v2_, t2};
  auto task547 = make_shared<Task547>(tensor547, cindex);
  task542->add_dep(task547);
  task547->add_dep(task538);
  deciq->add_task(task547);

  vector<IndexRange> I706_index = {active_, active_, active_, active_, active_, active_};
  auto I706 = make_shared<Tensor>(I706_index);
  auto tensor548 = vector<shared_ptr<Tensor>>{I698, Gamma250_(), I706};
  auto task548 = make_shared<Task548>(tensor548, cindex);
  task539->add_dep(task548);
  task548->add_dep(task538);
  deciq->add_task(task548);

  vector<IndexRange> I707_index = {active_, active_, closed_, active_};
  auto I707 = make_shared<Tensor>(I707_index);
  auto tensor549 = vector<shared_ptr<Tensor>>{I706, t2, I707};
  auto task549 = make_shared<Task549>(tensor549, cindex);
  task548->add_dep(task549);
  task549->add_dep(task538);
  deciq->add_task(task549);

  auto tensor550 = vector<shared_ptr<Tensor>>{I707, f1_, t2};
  auto task550 = make_shared<Task550>(tensor550, cindex);
  task549->add_dep(task550);
  task550->add_dep(task538);
  deciq->add_task(task550);

  vector<IndexRange> I710_index = {active_, active_, active_, active_};
  auto I710 = make_shared<Tensor>(I710_index);
  auto tensor551 = vector<shared_ptr<Tensor>>{I698, Gamma251_(), I710};
  auto task551 = make_shared<Task551>(tensor551, cindex);
  task539->add_dep(task551);
  task551->add_dep(task538);
  deciq->add_task(task551);

  vector<IndexRange> I711_index = {active_, closed_, closed_, active_};
  auto I711 = make_shared<Tensor>(I711_index);
  auto tensor552 = vector<shared_ptr<Tensor>>{I710, t2, I711};
  auto task552 = make_shared<Task552>(tensor552, cindex);
  task551->add_dep(task552);
  task552->add_dep(task538);
  deciq->add_task(task552);

  auto tensor553 = vector<shared_ptr<Tensor>>{I711, t2, f1_};
  auto task553 = make_shared<Task553>(tensor553, cindex);
  task552->add_dep(task553);
  task553->add_dep(task538);
  deciq->add_task(task553);

  vector<IndexRange> I742_index = {active_, closed_, closed_, active_};
  auto I742 = make_shared<Tensor>(I742_index);
  auto tensor554 = vector<shared_ptr<Tensor>>{I710, t2, I742};
  auto task554 = make_shared<Task554>(tensor554, cindex);
  task551->add_dep(task554);
  task554->add_dep(task538);
  deciq->add_task(task554);

  auto tensor555 = vector<shared_ptr<Tensor>>{I742, f1_, t2};
  auto task555 = make_shared<Task555>(tensor555, cindex);
  task554->add_dep(task555);
  task555->add_dep(task538);
  deciq->add_task(task555);

  vector<IndexRange> I714_index = {active_, active_, active_, active_, active_, active_};
  auto I714 = make_shared<Tensor>(I714_index);
  auto tensor556 = vector<shared_ptr<Tensor>>{I698, Gamma252_(), I714};
  auto task556 = make_shared<Task556>(tensor556, cindex);
  task539->add_dep(task556);
  task556->add_dep(task538);
  deciq->add_task(task556);

  vector<IndexRange> I715_index = {active_, closed_, active_, active_};
  auto I715 = make_shared<Tensor>(I715_index);
  auto tensor557 = vector<shared_ptr<Tensor>>{I714, t2, I715};
  auto task557 = make_shared<Task557>(tensor557, cindex);
  task556->add_dep(task557);
  task557->add_dep(task538);
  deciq->add_task(task557);

  auto tensor558 = vector<shared_ptr<Tensor>>{I715, t2, f1_};
  auto task558 = make_shared<Task558>(tensor558, cindex);
  task557->add_dep(task558);
  task558->add_dep(task538);
  deciq->add_task(task558);

  vector<IndexRange> I718_index = {active_, active_, active_, active_, active_, active_};
  auto I718 = make_shared<Tensor>(I718_index);
  auto tensor559 = vector<shared_ptr<Tensor>>{I698, Gamma253_(), I718};
  auto task559 = make_shared<Task559>(tensor559, cindex);
  task539->add_dep(task559);
  task559->add_dep(task538);
  deciq->add_task(task559);

  auto tensor560 = vector<shared_ptr<Tensor>>{I718, t2};
  auto task560 = make_shared<Task560>(tensor560, cindex);
  task559->add_dep(task560);
  task560->add_dep(task538);
  deciq->add_task(task560);

  vector<IndexRange> I721_index = {active_, active_, active_, active_, active_, active_};
  auto I721 = make_shared<Tensor>(I721_index);
  auto tensor561 = vector<shared_ptr<Tensor>>{I698, Gamma254_(), I721};
  auto task561 = make_shared<Task561>(tensor561, cindex);
  task539->add_dep(task561);
  task561->add_dep(task538);
  deciq->add_task(task561);

  vector<IndexRange> I722_index = {active_, active_, active_, closed_};
  auto I722 = make_shared<Tensor>(I722_index);
  auto tensor562 = vector<shared_ptr<Tensor>>{I721, t2, I722};
  auto task562 = make_shared<Task562>(tensor562, cindex);
  task561->add_dep(task562);
  task562->add_dep(task538);
  deciq->add_task(task562);

  auto tensor563 = vector<shared_ptr<Tensor>>{I722, f1_, t2};
  auto task563 = make_shared<Task563>(tensor563, cindex);
  task562->add_dep(task563);
  task563->add_dep(task538);
  deciq->add_task(task563);

  vector<IndexRange> I738_index = {active_, closed_, active_, active_};
  auto I738 = make_shared<Tensor>(I738_index);
  auto tensor564 = vector<shared_ptr<Tensor>>{I721, t2, I738};
  auto task564 = make_shared<Task564>(tensor564, cindex);
  task561->add_dep(task564);
  task564->add_dep(task538);
  deciq->add_task(task564);

  auto tensor565 = vector<shared_ptr<Tensor>>{I738, t2, f1_};
  auto task565 = make_shared<Task565>(tensor565, cindex);
  task564->add_dep(task565);
  task565->add_dep(task538);
  deciq->add_task(task565);

  vector<IndexRange> I862_index = {active_, active_, closed_, active_};
  auto I862 = make_shared<Tensor>(I862_index);
  auto tensor566 = vector<shared_ptr<Tensor>>{I721, t2, I862};
  auto task566 = make_shared<Task566>(tensor566, cindex);
  task561->add_dep(task566);
  task566->add_dep(task538);
  deciq->add_task(task566);

  auto tensor567 = vector<shared_ptr<Tensor>>{I862, t2};
  auto task567 = make_shared<Task567>(tensor567, cindex, this->e0_);
  task566->add_dep(task567);
  task567->add_dep(task538);
  deciq->add_task(task567);

  auto tensor568 = vector<shared_ptr<Tensor>>{I862, f1_, t2};
  auto task568 = make_shared<Task568>(tensor568, cindex);
  task566->add_dep(task568);
  task568->add_dep(task538);
  deciq->add_task(task568);

  auto tensor569 = vector<shared_ptr<Tensor>>{I721, v2_, t2};
  auto task569 = make_shared<Task569>(tensor569, cindex);
  task561->add_dep(task569);
  task569->add_dep(task538);
  deciq->add_task(task569);

  auto tensor570 = vector<shared_ptr<Tensor>>{I721, v2_, t2};
  auto task570 = make_shared<Task570>(tensor570, cindex);
  task561->add_dep(task570);
  task570->add_dep(task538);
  deciq->add_task(task570);

  vector<IndexRange> I725_index = {active_, active_, active_, active_};
  auto I725 = make_shared<Tensor>(I725_index);
  auto tensor571 = vector<shared_ptr<Tensor>>{I698, Gamma255_(), I725};
  auto task571 = make_shared<Task571>(tensor571, cindex);
  task539->add_dep(task571);
  task571->add_dep(task538);
  deciq->add_task(task571);

  vector<IndexRange> I726_index = {closed_, active_};
  auto I726 = make_shared<Tensor>(I726_index);
  auto tensor572 = vector<shared_ptr<Tensor>>{I725, t2, I726};
  auto task572 = make_shared<Task572>(tensor572, cindex);
  task571->add_dep(task572);
  task572->add_dep(task538);
  deciq->add_task(task572);

  auto tensor573 = vector<shared_ptr<Tensor>>{I726, t2, f1_};
  auto task573 = make_shared<Task573>(tensor573, cindex);
  task572->add_dep(task573);
  task573->add_dep(task538);
  deciq->add_task(task573);

  auto tensor574 = vector<shared_ptr<Tensor>>{I726, t2, f1_};
  auto task574 = make_shared<Task574>(tensor574, cindex);
  task572->add_dep(task574);
  task574->add_dep(task538);
  deciq->add_task(task574);

  vector<IndexRange> I816_index = {active_, closed_, virt_, active_};
  auto I816 = make_shared<Tensor>(I816_index);
  auto tensor575 = vector<shared_ptr<Tensor>>{I725, t2, I816};
  auto task575 = make_shared<Task575>(tensor575, cindex);
  task571->add_dep(task575);
  task575->add_dep(task538);
  deciq->add_task(task575);

  auto tensor576 = vector<shared_ptr<Tensor>>{I816, t2, f1_};
  auto task576 = make_shared<Task576>(tensor576, cindex);
  task575->add_dep(task576);
  task576->add_dep(task538);
  deciq->add_task(task576);

  vector<IndexRange> I866_index = {active_, virt_, closed_, active_};
  auto I866 = make_shared<Tensor>(I866_index);
  auto tensor577 = vector<shared_ptr<Tensor>>{I725, t2, I866};
  auto task577 = make_shared<Task577>(tensor577, cindex);
  task571->add_dep(task577);
  task577->add_dep(task538);
  deciq->add_task(task577);

  auto tensor578 = vector<shared_ptr<Tensor>>{I866, t2, f1_};
  auto task578 = make_shared<Task578>(tensor578, cindex);
  task577->add_dep(task578);
  task578->add_dep(task538);
  deciq->add_task(task578);

  auto tensor579 = vector<shared_ptr<Tensor>>{I866, t2, f1_};
  auto task579 = make_shared<Task579>(tensor579, cindex);
  task577->add_dep(task579);
  task579->add_dep(task538);
  deciq->add_task(task579);

  auto tensor580 = vector<shared_ptr<Tensor>>{I725, v2_, t2};
  auto task580 = make_shared<Task580>(tensor580, cindex);
  task571->add_dep(task580);
  task580->add_dep(task538);
  deciq->add_task(task580);

  auto tensor581 = vector<shared_ptr<Tensor>>{I725, h1_, t2};
  auto task581 = make_shared<Task581>(tensor581, cindex);
  task571->add_dep(task581);
  task581->add_dep(task538);
  deciq->add_task(task581);

  vector<IndexRange> I733_index = {active_, active_, active_, active_, active_, active_};
  auto I733 = make_shared<Tensor>(I733_index);
  auto tensor582 = vector<shared_ptr<Tensor>>{I698, Gamma257_(), I733};
  auto task582 = make_shared<Task582>(tensor582, cindex);
  task539->add_dep(task582);
  task582->add_dep(task538);
  deciq->add_task(task582);

  vector<IndexRange> I734_index = {active_, active_, closed_, active_};
  auto I734 = make_shared<Tensor>(I734_index);
  auto tensor583 = vector<shared_ptr<Tensor>>{I733, t2, I734};
  auto task583 = make_shared<Task583>(tensor583, cindex);
  task582->add_dep(task583);
  task583->add_dep(task538);
  deciq->add_task(task583);

  auto tensor584 = vector<shared_ptr<Tensor>>{I734, t2, f1_};
  auto task584 = make_shared<Task584>(tensor584, cindex);
  task583->add_dep(task584);
  task584->add_dep(task538);
  deciq->add_task(task584);

  vector<IndexRange> I745_index = {active_, active_, active_, active_};
  auto I745 = make_shared<Tensor>(I745_index);
  auto tensor585 = vector<shared_ptr<Tensor>>{I698, Gamma260_(), I745};
  auto task585 = make_shared<Task585>(tensor585, cindex);
  task539->add_dep(task585);
  task585->add_dep(task538);
  deciq->add_task(task585);

  vector<IndexRange> I746_index = {active_, closed_};
  auto I746 = make_shared<Tensor>(I746_index);
  auto tensor586 = vector<shared_ptr<Tensor>>{I745, t2, I746};
  auto task586 = make_shared<Task586>(tensor586, cindex);
  task585->add_dep(task586);
  task586->add_dep(task538);
  deciq->add_task(task586);

  auto tensor587 = vector<shared_ptr<Tensor>>{I746, f1_, t2};
  auto task587 = make_shared<Task587>(tensor587, cindex);
  task586->add_dep(task587);
  task587->add_dep(task538);
  deciq->add_task(task587);

  vector<IndexRange> I750_index = {active_, closed_};
  auto I750 = make_shared<Tensor>(I750_index);
  auto tensor588 = vector<shared_ptr<Tensor>>{I745, t2, I750};
  auto task588 = make_shared<Task588>(tensor588, cindex);
  task585->add_dep(task588);
  task588->add_dep(task538);
  deciq->add_task(task588);

  auto tensor589 = vector<shared_ptr<Tensor>>{I750, t2, f1_};
  auto task589 = make_shared<Task589>(tensor589, cindex);
  task588->add_dep(task589);
  task589->add_dep(task538);
  deciq->add_task(task589);

  vector<IndexRange> I788_index = {active_, virt_, closed_, active_};
  auto I788 = make_shared<Tensor>(I788_index);
  auto tensor590 = vector<shared_ptr<Tensor>>{I745, t2, I788};
  auto task590 = make_shared<Task590>(tensor590, cindex);
  task585->add_dep(task590);
  task590->add_dep(task538);
  deciq->add_task(task590);

  auto tensor591 = vector<shared_ptr<Tensor>>{I788, f1_, t2};
  auto task591 = make_shared<Task591>(tensor591, cindex);
  task590->add_dep(task591);
  task591->add_dep(task538);
  deciq->add_task(task591);

  vector<IndexRange> I792_index = {active_, closed_, virt_, active_};
  auto I792 = make_shared<Tensor>(I792_index);
  auto tensor592 = vector<shared_ptr<Tensor>>{I745, t2, I792};
  auto task592 = make_shared<Task592>(tensor592, cindex);
  task585->add_dep(task592);
  task592->add_dep(task538);
  deciq->add_task(task592);

  auto tensor593 = vector<shared_ptr<Tensor>>{I792, f1_, t2};
  auto task593 = make_shared<Task593>(tensor593, cindex);
  task592->add_dep(task593);
  task593->add_dep(task538);
  deciq->add_task(task593);

  vector<IndexRange> I796_index = {active_, virt_, closed_, active_};
  auto I796 = make_shared<Tensor>(I796_index);
  auto tensor594 = vector<shared_ptr<Tensor>>{I745, t2, I796};
  auto task594 = make_shared<Task594>(tensor594, cindex);
  task585->add_dep(task594);
  task594->add_dep(task538);
  deciq->add_task(task594);

  auto tensor595 = vector<shared_ptr<Tensor>>{I796, f1_, t2};
  auto task595 = make_shared<Task595>(tensor595, cindex);
  task594->add_dep(task595);
  task595->add_dep(task538);
  deciq->add_task(task595);

  auto tensor596 = vector<shared_ptr<Tensor>>{I745, v2_, t2};
  auto task596 = make_shared<Task596>(tensor596, cindex);
  task585->add_dep(task596);
  task596->add_dep(task538);
  deciq->add_task(task596);

  auto tensor597 = vector<shared_ptr<Tensor>>{I745, h1_, t2};
  auto task597 = make_shared<Task597>(tensor597, cindex);
  task585->add_dep(task597);
  task597->add_dep(task538);
  deciq->add_task(task597);

  vector<IndexRange> I753_index = {active_, active_};
  auto I753 = make_shared<Tensor>(I753_index);
  auto tensor598 = vector<shared_ptr<Tensor>>{I698, Gamma262_(), I753};
  auto task598 = make_shared<Task598>(tensor598, cindex);
  task539->add_dep(task598);
  task598->add_dep(task538);
  deciq->add_task(task598);

  auto tensor599 = vector<shared_ptr<Tensor>>{I753, t2};
  auto task599 = make_shared<Task599>(tensor599, cindex);
  task598->add_dep(task599);
  task599->add_dep(task538);
  deciq->add_task(task599);

  auto tensor600 = vector<shared_ptr<Tensor>>{I753, t2};
  auto task600 = make_shared<Task600>(tensor600, cindex);
  task598->add_dep(task600);
  task600->add_dep(task538);
  deciq->add_task(task600);
}

#endif
