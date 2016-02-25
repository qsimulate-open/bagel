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
  auto tensor539 = vector<shared_ptr<Tensor>>{deci};
  auto task539 = make_shared<Task539>(tensor539, reset);
  deciq->add_task(task539);

  vector<IndexRange> I684_index = {ci_};
  auto I684 = make_shared<Tensor>(I684_index);
  auto tensor540 = vector<shared_ptr<Tensor>>{deci, I684};
  auto task540 = make_shared<Task540>(tensor540, cindex);
  task540->add_dep(task539);
  deciq->add_task(task540);

  make_deciq1(deciq, task539, task540, diagonal, I684);
  make_deciq2(deciq, task539, task540, diagonal, I684);
  make_deciq3(deciq, task539, task540, diagonal, I684);

  return deciq;
}


void CASPT2::CASPT2::make_deciq1(shared_ptr<Queue> deciq, shared_ptr<Task> task539, shared_ptr<Task> task540, const bool diagonal, shared_ptr<Tensor> I684) {
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  vector<IndexRange> I685_index = {active_, active_, active_, active_};
  auto I685 = make_shared<Tensor>(I685_index);
  auto tensor541 = vector<shared_ptr<Tensor>>{I684, Gamma248_(), I685};
  auto task541 = make_shared<Task541>(tensor541, cindex);
  task540->add_dep(task541);
  task541->add_dep(task539);
  deciq->add_task(task541);

  auto tensor542 = vector<shared_ptr<Tensor>>{I685, t2};
  auto task542 = make_shared<Task542>(tensor542, cindex);
  task541->add_dep(task542);
  task542->add_dep(task539);
  deciq->add_task(task542);

  vector<IndexRange> I688_index = {active_, active_, active_, active_};
  auto I688 = make_shared<Tensor>(I688_index);
  auto tensor543 = vector<shared_ptr<Tensor>>{I684, Gamma249_(), I688};
  auto task543 = make_shared<Task543>(tensor543, cindex);
  task540->add_dep(task543);
  task543->add_dep(task539);
  deciq->add_task(task543);

  vector<IndexRange> I689_index = {active_, active_, closed_, closed_};
  auto I689 = make_shared<Tensor>(I689_index);
  auto tensor544 = vector<shared_ptr<Tensor>>{I688, t2, I689};
  auto task544 = make_shared<Task544>(tensor544, cindex);
  task543->add_dep(task544);
  task544->add_dep(task539);
  deciq->add_task(task544);

  auto tensor545 = vector<shared_ptr<Tensor>>{I689, f1_, t2};
  auto task545 = make_shared<Task545>(tensor545, cindex);
  task544->add_dep(task545);
  task545->add_dep(task539);
  deciq->add_task(task545);

  auto tensor546 = vector<shared_ptr<Tensor>>{I688, t2};
  auto task546 = make_shared<Task546>(tensor546, cindex, this->e0_);
  task543->add_dep(task546);
  task546->add_dep(task539);
  deciq->add_task(task546);

  auto tensor547 = vector<shared_ptr<Tensor>>{I688, v2_, t2};
  auto task547 = make_shared<Task547>(tensor547, cindex);
  task543->add_dep(task547);
  task547->add_dep(task539);
  deciq->add_task(task547);

  auto tensor548 = vector<shared_ptr<Tensor>>{I688, v2_, t2};
  auto task548 = make_shared<Task548>(tensor548, cindex);
  task543->add_dep(task548);
  task548->add_dep(task539);
  deciq->add_task(task548);

  vector<IndexRange> I692_index = {active_, active_, active_, active_, active_, active_};
  auto I692 = make_shared<Tensor>(I692_index);
  auto tensor549 = vector<shared_ptr<Tensor>>{I684, Gamma250_(), I692};
  auto task549 = make_shared<Task549>(tensor549, cindex);
  task540->add_dep(task549);
  task549->add_dep(task539);
  deciq->add_task(task549);

  vector<IndexRange> I693_index = {active_, active_, closed_, active_};
  auto I693 = make_shared<Tensor>(I693_index);
  auto tensor550 = vector<shared_ptr<Tensor>>{I692, t2, I693};
  auto task550 = make_shared<Task550>(tensor550, cindex);
  task549->add_dep(task550);
  task550->add_dep(task539);
  deciq->add_task(task550);

  auto tensor551 = vector<shared_ptr<Tensor>>{I693, f1_, t2};
  auto task551 = make_shared<Task551>(tensor551, cindex);
  task550->add_dep(task551);
  task551->add_dep(task539);
  deciq->add_task(task551);

  vector<IndexRange> I696_index = {active_, active_, active_, active_};
  auto I696 = make_shared<Tensor>(I696_index);
  auto tensor552 = vector<shared_ptr<Tensor>>{I684, Gamma251_(), I696};
  auto task552 = make_shared<Task552>(tensor552, cindex);
  task540->add_dep(task552);
  task552->add_dep(task539);
  deciq->add_task(task552);

  vector<IndexRange> I697_index = {active_, closed_, closed_, active_};
  auto I697 = make_shared<Tensor>(I697_index);
  auto tensor553 = vector<shared_ptr<Tensor>>{I696, t2, I697};
  auto task553 = make_shared<Task553>(tensor553, cindex);
  task552->add_dep(task553);
  task553->add_dep(task539);
  deciq->add_task(task553);

  auto tensor554 = vector<shared_ptr<Tensor>>{I697, t2, f1_};
  auto task554 = make_shared<Task554>(tensor554, cindex);
  task553->add_dep(task554);
  task554->add_dep(task539);
  deciq->add_task(task554);

  vector<IndexRange> I728_index = {active_, closed_, closed_, active_};
  auto I728 = make_shared<Tensor>(I728_index);
  auto tensor555 = vector<shared_ptr<Tensor>>{I696, t2, I728};
  auto task555 = make_shared<Task555>(tensor555, cindex);
  task552->add_dep(task555);
  task555->add_dep(task539);
  deciq->add_task(task555);

  auto tensor556 = vector<shared_ptr<Tensor>>{I728, f1_, t2};
  auto task556 = make_shared<Task556>(tensor556, cindex);
  task555->add_dep(task556);
  task556->add_dep(task539);
  deciq->add_task(task556);

  vector<IndexRange> I700_index = {active_, active_, active_, active_, active_, active_};
  auto I700 = make_shared<Tensor>(I700_index);
  auto tensor557 = vector<shared_ptr<Tensor>>{I684, Gamma252_(), I700};
  auto task557 = make_shared<Task557>(tensor557, cindex);
  task540->add_dep(task557);
  task557->add_dep(task539);
  deciq->add_task(task557);

  vector<IndexRange> I701_index = {active_, closed_, active_, active_};
  auto I701 = make_shared<Tensor>(I701_index);
  auto tensor558 = vector<shared_ptr<Tensor>>{I700, t2, I701};
  auto task558 = make_shared<Task558>(tensor558, cindex);
  task557->add_dep(task558);
  task558->add_dep(task539);
  deciq->add_task(task558);

  auto tensor559 = vector<shared_ptr<Tensor>>{I701, t2, f1_};
  auto task559 = make_shared<Task559>(tensor559, cindex);
  task558->add_dep(task559);
  task559->add_dep(task539);
  deciq->add_task(task559);

  vector<IndexRange> I704_index = {active_, active_, active_, active_, active_, active_};
  auto I704 = make_shared<Tensor>(I704_index);
  auto tensor560 = vector<shared_ptr<Tensor>>{I684, Gamma253_(), I704};
  auto task560 = make_shared<Task560>(tensor560, cindex);
  task540->add_dep(task560);
  task560->add_dep(task539);
  deciq->add_task(task560);

  auto tensor561 = vector<shared_ptr<Tensor>>{I704, t2};
  auto task561 = make_shared<Task561>(tensor561, cindex);
  task560->add_dep(task561);
  task561->add_dep(task539);
  deciq->add_task(task561);

  vector<IndexRange> I707_index = {active_, active_, active_, active_, active_, active_};
  auto I707 = make_shared<Tensor>(I707_index);
  auto tensor562 = vector<shared_ptr<Tensor>>{I684, Gamma254_(), I707};
  auto task562 = make_shared<Task562>(tensor562, cindex);
  task540->add_dep(task562);
  task562->add_dep(task539);
  deciq->add_task(task562);

  vector<IndexRange> I708_index = {active_, active_, active_, closed_};
  auto I708 = make_shared<Tensor>(I708_index);
  auto tensor563 = vector<shared_ptr<Tensor>>{I707, t2, I708};
  auto task563 = make_shared<Task563>(tensor563, cindex);
  task562->add_dep(task563);
  task563->add_dep(task539);
  deciq->add_task(task563);

  auto tensor564 = vector<shared_ptr<Tensor>>{I708, f1_, t2};
  auto task564 = make_shared<Task564>(tensor564, cindex);
  task563->add_dep(task564);
  task564->add_dep(task539);
  deciq->add_task(task564);

  vector<IndexRange> I724_index = {active_, closed_, active_, active_};
  auto I724 = make_shared<Tensor>(I724_index);
  auto tensor565 = vector<shared_ptr<Tensor>>{I707, t2, I724};
  auto task565 = make_shared<Task565>(tensor565, cindex);
  task562->add_dep(task565);
  task565->add_dep(task539);
  deciq->add_task(task565);

  auto tensor566 = vector<shared_ptr<Tensor>>{I724, t2, f1_};
  auto task566 = make_shared<Task566>(tensor566, cindex);
  task565->add_dep(task566);
  task566->add_dep(task539);
  deciq->add_task(task566);

  vector<IndexRange> I848_index = {active_, active_, closed_, active_};
  auto I848 = make_shared<Tensor>(I848_index);
  auto tensor567 = vector<shared_ptr<Tensor>>{I707, t2, I848};
  auto task567 = make_shared<Task567>(tensor567, cindex);
  task562->add_dep(task567);
  task567->add_dep(task539);
  deciq->add_task(task567);

  auto tensor568 = vector<shared_ptr<Tensor>>{I848, t2};
  auto task568 = make_shared<Task568>(tensor568, cindex, this->e0_);
  task567->add_dep(task568);
  task568->add_dep(task539);
  deciq->add_task(task568);

  auto tensor569 = vector<shared_ptr<Tensor>>{I848, f1_, t2};
  auto task569 = make_shared<Task569>(tensor569, cindex);
  task567->add_dep(task569);
  task569->add_dep(task539);
  deciq->add_task(task569);

  auto tensor570 = vector<shared_ptr<Tensor>>{I707, v2_, t2};
  auto task570 = make_shared<Task570>(tensor570, cindex);
  task562->add_dep(task570);
  task570->add_dep(task539);
  deciq->add_task(task570);

  auto tensor571 = vector<shared_ptr<Tensor>>{I707, v2_, t2};
  auto task571 = make_shared<Task571>(tensor571, cindex);
  task562->add_dep(task571);
  task571->add_dep(task539);
  deciq->add_task(task571);

  vector<IndexRange> I711_index = {active_, active_, active_, active_};
  auto I711 = make_shared<Tensor>(I711_index);
  auto tensor572 = vector<shared_ptr<Tensor>>{I684, Gamma255_(), I711};
  auto task572 = make_shared<Task572>(tensor572, cindex);
  task540->add_dep(task572);
  task572->add_dep(task539);
  deciq->add_task(task572);

  vector<IndexRange> I712_index = {closed_, active_};
  auto I712 = make_shared<Tensor>(I712_index);
  auto tensor573 = vector<shared_ptr<Tensor>>{I711, t2, I712};
  auto task573 = make_shared<Task573>(tensor573, cindex);
  task572->add_dep(task573);
  task573->add_dep(task539);
  deciq->add_task(task573);

  auto tensor574 = vector<shared_ptr<Tensor>>{I712, t2, f1_};
  auto task574 = make_shared<Task574>(tensor574, cindex);
  task573->add_dep(task574);
  task574->add_dep(task539);
  deciq->add_task(task574);

  auto tensor575 = vector<shared_ptr<Tensor>>{I712, t2, f1_};
  auto task575 = make_shared<Task575>(tensor575, cindex);
  task573->add_dep(task575);
  task575->add_dep(task539);
  deciq->add_task(task575);

  vector<IndexRange> I802_index = {active_, closed_, virt_, active_};
  auto I802 = make_shared<Tensor>(I802_index);
  auto tensor576 = vector<shared_ptr<Tensor>>{I711, t2, I802};
  auto task576 = make_shared<Task576>(tensor576, cindex);
  task572->add_dep(task576);
  task576->add_dep(task539);
  deciq->add_task(task576);

  auto tensor577 = vector<shared_ptr<Tensor>>{I802, t2, f1_};
  auto task577 = make_shared<Task577>(tensor577, cindex);
  task576->add_dep(task577);
  task577->add_dep(task539);
  deciq->add_task(task577);

  vector<IndexRange> I852_index = {active_, virt_, closed_, active_};
  auto I852 = make_shared<Tensor>(I852_index);
  auto tensor578 = vector<shared_ptr<Tensor>>{I711, t2, I852};
  auto task578 = make_shared<Task578>(tensor578, cindex);
  task572->add_dep(task578);
  task578->add_dep(task539);
  deciq->add_task(task578);

  auto tensor579 = vector<shared_ptr<Tensor>>{I852, t2, f1_};
  auto task579 = make_shared<Task579>(tensor579, cindex);
  task578->add_dep(task579);
  task579->add_dep(task539);
  deciq->add_task(task579);

  auto tensor580 = vector<shared_ptr<Tensor>>{I852, t2, f1_};
  auto task580 = make_shared<Task580>(tensor580, cindex);
  task578->add_dep(task580);
  task580->add_dep(task539);
  deciq->add_task(task580);

  auto tensor581 = vector<shared_ptr<Tensor>>{I711, v2_, t2};
  auto task581 = make_shared<Task581>(tensor581, cindex);
  task572->add_dep(task581);
  task581->add_dep(task539);
  deciq->add_task(task581);

  auto tensor582 = vector<shared_ptr<Tensor>>{I711, h1_, t2};
  auto task582 = make_shared<Task582>(tensor582, cindex);
  task572->add_dep(task582);
  task582->add_dep(task539);
  deciq->add_task(task582);

  vector<IndexRange> I719_index = {active_, active_, active_, active_, active_, active_};
  auto I719 = make_shared<Tensor>(I719_index);
  auto tensor583 = vector<shared_ptr<Tensor>>{I684, Gamma257_(), I719};
  auto task583 = make_shared<Task583>(tensor583, cindex);
  task540->add_dep(task583);
  task583->add_dep(task539);
  deciq->add_task(task583);

  vector<IndexRange> I720_index = {active_, active_, closed_, active_};
  auto I720 = make_shared<Tensor>(I720_index);
  auto tensor584 = vector<shared_ptr<Tensor>>{I719, t2, I720};
  auto task584 = make_shared<Task584>(tensor584, cindex);
  task583->add_dep(task584);
  task584->add_dep(task539);
  deciq->add_task(task584);

  auto tensor585 = vector<shared_ptr<Tensor>>{I720, t2, f1_};
  auto task585 = make_shared<Task585>(tensor585, cindex);
  task584->add_dep(task585);
  task585->add_dep(task539);
  deciq->add_task(task585);

  vector<IndexRange> I731_index = {active_, active_, active_, active_};
  auto I731 = make_shared<Tensor>(I731_index);
  auto tensor586 = vector<shared_ptr<Tensor>>{I684, Gamma260_(), I731};
  auto task586 = make_shared<Task586>(tensor586, cindex);
  task540->add_dep(task586);
  task586->add_dep(task539);
  deciq->add_task(task586);

  vector<IndexRange> I732_index = {active_, closed_};
  auto I732 = make_shared<Tensor>(I732_index);
  auto tensor587 = vector<shared_ptr<Tensor>>{I731, t2, I732};
  auto task587 = make_shared<Task587>(tensor587, cindex);
  task586->add_dep(task587);
  task587->add_dep(task539);
  deciq->add_task(task587);

  auto tensor588 = vector<shared_ptr<Tensor>>{I732, f1_, t2};
  auto task588 = make_shared<Task588>(tensor588, cindex);
  task587->add_dep(task588);
  task588->add_dep(task539);
  deciq->add_task(task588);

  vector<IndexRange> I736_index = {active_, closed_};
  auto I736 = make_shared<Tensor>(I736_index);
  auto tensor589 = vector<shared_ptr<Tensor>>{I731, t2, I736};
  auto task589 = make_shared<Task589>(tensor589, cindex);
  task586->add_dep(task589);
  task589->add_dep(task539);
  deciq->add_task(task589);

  auto tensor590 = vector<shared_ptr<Tensor>>{I736, t2, f1_};
  auto task590 = make_shared<Task590>(tensor590, cindex);
  task589->add_dep(task590);
  task590->add_dep(task539);
  deciq->add_task(task590);

  vector<IndexRange> I774_index = {active_, virt_, closed_, active_};
  auto I774 = make_shared<Tensor>(I774_index);
  auto tensor591 = vector<shared_ptr<Tensor>>{I731, t2, I774};
  auto task591 = make_shared<Task591>(tensor591, cindex);
  task586->add_dep(task591);
  task591->add_dep(task539);
  deciq->add_task(task591);

  auto tensor592 = vector<shared_ptr<Tensor>>{I774, f1_, t2};
  auto task592 = make_shared<Task592>(tensor592, cindex);
  task591->add_dep(task592);
  task592->add_dep(task539);
  deciq->add_task(task592);

  vector<IndexRange> I778_index = {active_, closed_, virt_, active_};
  auto I778 = make_shared<Tensor>(I778_index);
  auto tensor593 = vector<shared_ptr<Tensor>>{I731, t2, I778};
  auto task593 = make_shared<Task593>(tensor593, cindex);
  task586->add_dep(task593);
  task593->add_dep(task539);
  deciq->add_task(task593);

  auto tensor594 = vector<shared_ptr<Tensor>>{I778, f1_, t2};
  auto task594 = make_shared<Task594>(tensor594, cindex);
  task593->add_dep(task594);
  task594->add_dep(task539);
  deciq->add_task(task594);

  vector<IndexRange> I782_index = {active_, virt_, closed_, active_};
  auto I782 = make_shared<Tensor>(I782_index);
  auto tensor595 = vector<shared_ptr<Tensor>>{I731, t2, I782};
  auto task595 = make_shared<Task595>(tensor595, cindex);
  task586->add_dep(task595);
  task595->add_dep(task539);
  deciq->add_task(task595);

  auto tensor596 = vector<shared_ptr<Tensor>>{I782, f1_, t2};
  auto task596 = make_shared<Task596>(tensor596, cindex);
  task595->add_dep(task596);
  task596->add_dep(task539);
  deciq->add_task(task596);

  auto tensor597 = vector<shared_ptr<Tensor>>{I731, v2_, t2};
  auto task597 = make_shared<Task597>(tensor597, cindex);
  task586->add_dep(task597);
  task597->add_dep(task539);
  deciq->add_task(task597);

  auto tensor598 = vector<shared_ptr<Tensor>>{I731, h1_, t2};
  auto task598 = make_shared<Task598>(tensor598, cindex);
  task586->add_dep(task598);
  task598->add_dep(task539);
  deciq->add_task(task598);

  vector<IndexRange> I739_index = {active_, active_};
  auto I739 = make_shared<Tensor>(I739_index);
  auto tensor599 = vector<shared_ptr<Tensor>>{I684, Gamma262_(), I739};
  auto task599 = make_shared<Task599>(tensor599, cindex);
  task540->add_dep(task599);
  task599->add_dep(task539);
  deciq->add_task(task599);

  auto tensor600 = vector<shared_ptr<Tensor>>{I739, t2};
  auto task600 = make_shared<Task600>(tensor600, cindex);
  task599->add_dep(task600);
  task600->add_dep(task539);
  deciq->add_task(task600);

  auto tensor601 = vector<shared_ptr<Tensor>>{I739, t2};
  auto task601 = make_shared<Task601>(tensor601, cindex);
  task599->add_dep(task601);
  task601->add_dep(task539);
  deciq->add_task(task601);
}

#endif
