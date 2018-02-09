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

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto deciq = make_shared<Queue>();
  auto tensor519 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci};
  auto task519 = make_shared<Task519>(tensor519, reset);
  deciq->add_task(task519);

  vector<IndexRange> I685_index = {active_, active_, active_, active_};
  auto I685 = make_shared<Tensor>(I685_index);
  auto tensor521 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I685, f1_};
  auto task521 = make_shared<Task521>(tensor521, pindex);
  task521->add_dep(task519);
  deciq->add_task(task521);

  auto tensor522 = vector<shared_ptr<Tensor>>{I685, t2};
  auto task522 = make_shared<Task522>(tensor522, pindex);
  task521->add_dep(task522);
  task522->add_dep(task519);
  deciq->add_task(task522);

  vector<IndexRange> I688_index = {active_, active_, active_, active_};
  auto I688 = make_shared<Tensor>(I688_index);
  auto tensor523 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I688};
  auto task523 = make_shared<Task523>(tensor523, pindex);
  task523->add_dep(task519);
  deciq->add_task(task523);

  vector<IndexRange> I689_index = {active_, active_, closed_, closed_};
  auto I689 = make_shared<Tensor>(I689_index);
  auto tensor524 = vector<shared_ptr<Tensor>>{I688, t2, I689};
  auto task524 = make_shared<Task524>(tensor524, pindex);
  task523->add_dep(task524);
  task524->add_dep(task519);
  deciq->add_task(task524);

  auto tensor525 = vector<shared_ptr<Tensor>>{I689, f1_, t2};
  auto task525 = make_shared<Task525>(tensor525, pindex);
  task524->add_dep(task525);
  task525->add_dep(task519);
  deciq->add_task(task525);

  auto tensor526 = vector<shared_ptr<Tensor>>{I688, t2};
  auto task526 = make_shared<Task526>(tensor526, pindex, this->e0_);
  task523->add_dep(task526);
  task526->add_dep(task519);
  deciq->add_task(task526);

  auto tensor527 = vector<shared_ptr<Tensor>>{I688, v2_, t2};
  auto task527 = make_shared<Task527>(tensor527, pindex);
  task523->add_dep(task527);
  task527->add_dep(task519);
  deciq->add_task(task527);

  auto tensor528 = vector<shared_ptr<Tensor>>{I688, v2_, t2};
  auto task528 = make_shared<Task528>(tensor528, pindex);
  task523->add_dep(task528);
  task528->add_dep(task519);
  deciq->add_task(task528);

  vector<IndexRange> I692_index = {active_, active_, active_, active_, active_, active_};
  auto I692 = make_shared<Tensor>(I692_index);
  auto tensor529 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I692};
  auto task529 = make_shared<Task529>(tensor529, pindex);
  task529->add_dep(task519);
  deciq->add_task(task529);

  vector<IndexRange> I693_index = {active_, active_, closed_, active_};
  auto I693 = make_shared<Tensor>(I693_index);
  auto tensor530 = vector<shared_ptr<Tensor>>{I692, t2, I693};
  auto task530 = make_shared<Task530>(tensor530, pindex);
  task529->add_dep(task530);
  task530->add_dep(task519);
  deciq->add_task(task530);

  auto tensor531 = vector<shared_ptr<Tensor>>{I693, f1_, t2};
  auto task531 = make_shared<Task531>(tensor531, pindex);
  task530->add_dep(task531);
  task531->add_dep(task519);
  deciq->add_task(task531);

  vector<IndexRange> I696_index = {active_, active_, active_, active_};
  auto I696 = make_shared<Tensor>(I696_index);
  auto tensor532 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I696};
  auto task532 = make_shared<Task532>(tensor532, pindex);
  task532->add_dep(task519);
  deciq->add_task(task532);

  vector<IndexRange> I697_index = {closed_, closed_, active_, active_};
  auto I697 = make_shared<Tensor>(I697_index);
  auto tensor533 = vector<shared_ptr<Tensor>>{I696, t2, I697};
  auto task533 = make_shared<Task533>(tensor533, pindex);
  task532->add_dep(task533);
  task533->add_dep(task519);
  deciq->add_task(task533);

  auto tensor534 = vector<shared_ptr<Tensor>>{I697, f1_, t2};
  auto task534 = make_shared<Task534>(tensor534, pindex);
  task533->add_dep(task534);
  task534->add_dep(task519);
  deciq->add_task(task534);

  vector<IndexRange> I728_index = {active_, closed_, closed_, active_};
  auto I728 = make_shared<Tensor>(I728_index);
  auto tensor535 = vector<shared_ptr<Tensor>>{I696, t2, I728};
  auto task535 = make_shared<Task535>(tensor535, pindex);
  task532->add_dep(task535);
  task535->add_dep(task519);
  deciq->add_task(task535);

  auto tensor536 = vector<shared_ptr<Tensor>>{I728, f1_, t2};
  auto task536 = make_shared<Task536>(tensor536, pindex);
  task535->add_dep(task536);
  task536->add_dep(task519);
  deciq->add_task(task536);

  vector<IndexRange> I700_index = {active_, active_, active_, active_, active_, active_};
  auto I700 = make_shared<Tensor>(I700_index);
  auto tensor537 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I700};
  auto task537 = make_shared<Task537>(tensor537, pindex);
  task537->add_dep(task519);
  deciq->add_task(task537);

  vector<IndexRange> I701_index = {closed_, active_, active_, active_};
  auto I701 = make_shared<Tensor>(I701_index);
  auto tensor538 = vector<shared_ptr<Tensor>>{I700, t2, I701};
  auto task538 = make_shared<Task538>(tensor538, pindex);
  task537->add_dep(task538);
  task538->add_dep(task519);
  deciq->add_task(task538);

  auto tensor539 = vector<shared_ptr<Tensor>>{I701, f1_, t2};
  auto task539 = make_shared<Task539>(tensor539, pindex);
  task538->add_dep(task539);
  task539->add_dep(task519);
  deciq->add_task(task539);

  vector<IndexRange> I704_index = {active_, active_, active_, active_, active_, active_};
  auto I704 = make_shared<Tensor>(I704_index);
  auto tensor540 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I704, f1_};
  auto task540 = make_shared<Task540>(tensor540, pindex);
  task540->add_dep(task519);
  deciq->add_task(task540);

  auto tensor541 = vector<shared_ptr<Tensor>>{I704, t2};
  auto task541 = make_shared<Task541>(tensor541, pindex);
  task540->add_dep(task541);
  task541->add_dep(task519);
  deciq->add_task(task541);

  vector<IndexRange> I707_index = {active_, active_, active_, active_, active_, active_};
  auto I707 = make_shared<Tensor>(I707_index);
  auto tensor542 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I707};
  auto task542 = make_shared<Task542>(tensor542, pindex);
  task542->add_dep(task519);
  deciq->add_task(task542);

  vector<IndexRange> I708_index = {active_, active_, active_, closed_};
  auto I708 = make_shared<Tensor>(I708_index);
  auto tensor543 = vector<shared_ptr<Tensor>>{I707, t2, I708};
  auto task543 = make_shared<Task543>(tensor543, pindex);
  task542->add_dep(task543);
  task543->add_dep(task519);
  deciq->add_task(task543);

  auto tensor544 = vector<shared_ptr<Tensor>>{I708, f1_, t2};
  auto task544 = make_shared<Task544>(tensor544, pindex);
  task543->add_dep(task544);
  task544->add_dep(task519);
  deciq->add_task(task544);

  vector<IndexRange> I724_index = {closed_, active_, active_, active_};
  auto I724 = make_shared<Tensor>(I724_index);
  auto tensor545 = vector<shared_ptr<Tensor>>{I707, t2, I724};
  auto task545 = make_shared<Task545>(tensor545, pindex);
  task542->add_dep(task545);
  task545->add_dep(task519);
  deciq->add_task(task545);

  auto tensor546 = vector<shared_ptr<Tensor>>{I724, t2};
  auto task546 = make_shared<Task546>(tensor546, pindex, this->e0_);
  task545->add_dep(task546);
  task546->add_dep(task519);
  deciq->add_task(task546);

  auto tensor547 = vector<shared_ptr<Tensor>>{I724, f1_, t2};
  auto task547 = make_shared<Task547>(tensor547, pindex);
  task545->add_dep(task547);
  task547->add_dep(task519);
  deciq->add_task(task547);

  vector<IndexRange> I848_index = {active_, active_, closed_, active_};
  auto I848 = make_shared<Tensor>(I848_index);
  auto tensor548 = vector<shared_ptr<Tensor>>{I707, t2, I848};
  auto task548 = make_shared<Task548>(tensor548, pindex);
  task542->add_dep(task548);
  task548->add_dep(task519);
  deciq->add_task(task548);

  auto tensor549 = vector<shared_ptr<Tensor>>{I848, f1_, t2};
  auto task549 = make_shared<Task549>(tensor549, pindex);
  task548->add_dep(task549);
  task549->add_dep(task519);
  deciq->add_task(task549);

  auto tensor550 = vector<shared_ptr<Tensor>>{I707, v2_, t2};
  auto task550 = make_shared<Task550>(tensor550, pindex);
  task542->add_dep(task550);
  task550->add_dep(task519);
  deciq->add_task(task550);

  auto tensor551 = vector<shared_ptr<Tensor>>{I707, v2_, t2};
  auto task551 = make_shared<Task551>(tensor551, pindex);
  task542->add_dep(task551);
  task551->add_dep(task519);
  deciq->add_task(task551);

  vector<IndexRange> I711_index = {active_, active_, active_, active_};
  auto I711 = make_shared<Tensor>(I711_index);
  auto tensor552 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I711};
  auto task552 = make_shared<Task552>(tensor552, pindex);
  task552->add_dep(task519);
  deciq->add_task(task552);

  vector<IndexRange> I712_index = {closed_, active_};
  auto I712 = make_shared<Tensor>(I712_index);
  auto tensor553 = vector<shared_ptr<Tensor>>{I711, t2, I712};
  auto task553 = make_shared<Task553>(tensor553, pindex);
  task552->add_dep(task553);
  task553->add_dep(task519);
  deciq->add_task(task553);

  vector<IndexRange> I713_index = {closed_, virt_, closed_, active_};
  auto I713 = make_shared<Tensor>(I713_index);
  auto tensor554 = vector<shared_ptr<Tensor>>{I712, f1_, I713};
  auto task554 = make_shared<Task554>(tensor554, pindex);
  task553->add_dep(task554);
  task554->add_dep(task519);
  deciq->add_task(task554);

  auto tensor555 = vector<shared_ptr<Tensor>>{I713, t2};
  auto task555 = make_shared<Task555>(tensor555, pindex);
  task554->add_dep(task555);
  task555->add_dep(task519);
  deciq->add_task(task555);

  vector<IndexRange> I802_index = {closed_, virt_, active_, active_};
  auto I802 = make_shared<Tensor>(I802_index);
  auto tensor556 = vector<shared_ptr<Tensor>>{I711, t2, I802};
  auto task556 = make_shared<Task556>(tensor556, pindex);
  task552->add_dep(task556);
  task556->add_dep(task519);
  deciq->add_task(task556);

  auto tensor557 = vector<shared_ptr<Tensor>>{I802, f1_, t2};
  auto task557 = make_shared<Task557>(tensor557, pindex);
  task556->add_dep(task557);
  task557->add_dep(task519);
  deciq->add_task(task557);

  vector<IndexRange> I852_index = {virt_, closed_, active_, active_};
  auto I852 = make_shared<Tensor>(I852_index);
  auto tensor558 = vector<shared_ptr<Tensor>>{I711, t2, I852};
  auto task558 = make_shared<Task558>(tensor558, pindex);
  task552->add_dep(task558);
  task558->add_dep(task519);
  deciq->add_task(task558);

  vector<IndexRange> I853_index = {closed_, virt_, closed_, active_};
  auto I853 = make_shared<Tensor>(I853_index);
  auto tensor559 = vector<shared_ptr<Tensor>>{I852, f1_, I853};
  auto task559 = make_shared<Task559>(tensor559, pindex);
  task558->add_dep(task559);
  task559->add_dep(task519);
  deciq->add_task(task559);

  auto tensor560 = vector<shared_ptr<Tensor>>{I853, t2};
  auto task560 = make_shared<Task560>(tensor560, pindex);
  task559->add_dep(task560);
  task560->add_dep(task519);
  deciq->add_task(task560);

  auto tensor561 = vector<shared_ptr<Tensor>>{I711, v2_, t2};
  auto task561 = make_shared<Task561>(tensor561, pindex);
  task552->add_dep(task561);
  task561->add_dep(task519);
  deciq->add_task(task561);

  auto tensor562 = vector<shared_ptr<Tensor>>{I711, h1_, t2};
  auto task562 = make_shared<Task562>(tensor562, pindex);
  task552->add_dep(task562);
  task562->add_dep(task519);
  deciq->add_task(task562);

  vector<IndexRange> I719_index = {active_, active_, active_, active_, active_, active_};
  auto I719 = make_shared<Tensor>(I719_index);
  auto tensor563 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I719};
  auto task563 = make_shared<Task563>(tensor563, pindex);
  task563->add_dep(task519);
  deciq->add_task(task563);

  vector<IndexRange> I720_index = {active_, closed_, active_, active_};
  auto I720 = make_shared<Tensor>(I720_index);
  auto tensor564 = vector<shared_ptr<Tensor>>{I719, t2, I720};
  auto task564 = make_shared<Task564>(tensor564, pindex);
  task563->add_dep(task564);
  task564->add_dep(task519);
  deciq->add_task(task564);

  auto tensor565 = vector<shared_ptr<Tensor>>{I720, f1_, t2};
  auto task565 = make_shared<Task565>(tensor565, pindex);
  task564->add_dep(task565);
  task565->add_dep(task519);
  deciq->add_task(task565);

  vector<IndexRange> I731_index = {active_, active_, active_, active_};
  auto I731 = make_shared<Tensor>(I731_index);
  auto tensor566 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I731};
  auto task566 = make_shared<Task566>(tensor566, pindex);
  task566->add_dep(task519);
  deciq->add_task(task566);

  vector<IndexRange> I732_index = {active_, closed_};
  auto I732 = make_shared<Tensor>(I732_index);
  auto tensor567 = vector<shared_ptr<Tensor>>{I731, t2, I732};
  auto task567 = make_shared<Task567>(tensor567, pindex);
  task566->add_dep(task567);
  task567->add_dep(task519);
  deciq->add_task(task567);

  auto tensor568 = vector<shared_ptr<Tensor>>{I732, f1_, t2};
  auto task568 = make_shared<Task568>(tensor568, pindex);
  task567->add_dep(task568);
  task568->add_dep(task519);
  deciq->add_task(task568);

  vector<IndexRange> I736_index = {active_, closed_};
  auto I736 = make_shared<Tensor>(I736_index);
  auto tensor569 = vector<shared_ptr<Tensor>>{I731, t2, I736};
  auto task569 = make_shared<Task569>(tensor569, pindex);
  task566->add_dep(task569);
  task569->add_dep(task519);
  deciq->add_task(task569);

  auto tensor570 = vector<shared_ptr<Tensor>>{I736, f1_, t2};
  auto task570 = make_shared<Task570>(tensor570, pindex);
  task569->add_dep(task570);
  task570->add_dep(task519);
  deciq->add_task(task570);

  vector<IndexRange> I774_index = {active_, virt_, closed_, active_};
  auto I774 = make_shared<Tensor>(I774_index);
  auto tensor571 = vector<shared_ptr<Tensor>>{I731, t2, I774};
  auto task571 = make_shared<Task571>(tensor571, pindex);
  task566->add_dep(task571);
  task571->add_dep(task519);
  deciq->add_task(task571);

  auto tensor572 = vector<shared_ptr<Tensor>>{I774, f1_, t2};
  auto task572 = make_shared<Task572>(tensor572, pindex);
  task571->add_dep(task572);
  task572->add_dep(task519);
  deciq->add_task(task572);

  vector<IndexRange> I778_index = {active_, closed_, virt_, active_};
  auto I778 = make_shared<Tensor>(I778_index);
  auto tensor573 = vector<shared_ptr<Tensor>>{I731, t2, I778};
  auto task573 = make_shared<Task573>(tensor573, pindex);
  task566->add_dep(task573);
  task573->add_dep(task519);
  deciq->add_task(task573);

  auto tensor574 = vector<shared_ptr<Tensor>>{I778, f1_, t2};
  auto task574 = make_shared<Task574>(tensor574, pindex);
  task573->add_dep(task574);
  task574->add_dep(task519);
  deciq->add_task(task574);

  vector<IndexRange> I782_index = {active_, virt_, closed_, active_};
  auto I782 = make_shared<Tensor>(I782_index);
  auto tensor575 = vector<shared_ptr<Tensor>>{I731, t2, I782};
  auto task575 = make_shared<Task575>(tensor575, pindex);
  task566->add_dep(task575);
  task575->add_dep(task519);
  deciq->add_task(task575);

  auto tensor576 = vector<shared_ptr<Tensor>>{I782, f1_, t2};
  auto task576 = make_shared<Task576>(tensor576, pindex);
  task575->add_dep(task576);
  task576->add_dep(task519);
  deciq->add_task(task576);

  auto tensor577 = vector<shared_ptr<Tensor>>{I731, v2_, t2};
  auto task577 = make_shared<Task577>(tensor577, pindex);
  task566->add_dep(task577);
  task577->add_dep(task519);
  deciq->add_task(task577);

  auto tensor578 = vector<shared_ptr<Tensor>>{I731, h1_, t2};
  auto task578 = make_shared<Task578>(tensor578, pindex);
  task566->add_dep(task578);
  task578->add_dep(task519);
  deciq->add_task(task578);

  vector<IndexRange> I739_index = {active_, active_};
  auto I739 = make_shared<Tensor>(I739_index);
  auto tensor579 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I739, f1_};
  auto task579 = make_shared<Task579>(tensor579, pindex);
  task579->add_dep(task519);
  deciq->add_task(task579);

  vector<IndexRange> I740_index = {closed_, virt_, closed_, active_};
  auto I740 = make_shared<Tensor>(I740_index);
  auto tensor580 = vector<shared_ptr<Tensor>>{I739, t2, I740};
  auto task580 = make_shared<Task580>(tensor580, pindex);
  task579->add_dep(task580);
  task580->add_dep(task519);
  deciq->add_task(task580);

  auto tensor581 = vector<shared_ptr<Tensor>>{I740, t2};
  auto task581 = make_shared<Task581>(tensor581, pindex);
  task580->add_dep(task581);
  task581->add_dep(task519);
  deciq->add_task(task581);

  vector<IndexRange> I745_index = {active_, active_};
  auto I745 = make_shared<Tensor>(I745_index);
  auto tensor582 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I745};
  auto task582 = make_shared<Task582>(tensor582, pindex);
  task582->add_dep(task519);
  deciq->add_task(task582);

  vector<IndexRange> I746_index = {active_, closed_, virt_, closed_};
  auto I746 = make_shared<Tensor>(I746_index);
  auto tensor583 = vector<shared_ptr<Tensor>>{I745, t2, I746};
  auto task583 = make_shared<Task583>(tensor583, pindex);
  task582->add_dep(task583);
  task583->add_dep(task519);
  deciq->add_task(task583);

  auto tensor584 = vector<shared_ptr<Tensor>>{I746, f1_, t2};
  auto task584 = make_shared<Task584>(tensor584, pindex);
  task583->add_dep(task584);
  task584->add_dep(task519);
  deciq->add_task(task584);

  vector<IndexRange> I750_index = {active_, closed_, virt_, closed_};
  auto I750 = make_shared<Tensor>(I750_index);
  auto tensor585 = vector<shared_ptr<Tensor>>{I745, t2, I750};
  auto task585 = make_shared<Task585>(tensor585, pindex);
  task582->add_dep(task585);
  task585->add_dep(task519);
  deciq->add_task(task585);

  auto tensor586 = vector<shared_ptr<Tensor>>{I750, f1_, t2};
  auto task586 = make_shared<Task586>(tensor586, pindex);
  task585->add_dep(task586);
  task586->add_dep(task519);
  deciq->add_task(task586);

  vector<IndexRange> I754_index = {active_, virt_, closed_, closed_};
  auto I754 = make_shared<Tensor>(I754_index);
  auto tensor587 = vector<shared_ptr<Tensor>>{I745, t2, I754};
  auto task587 = make_shared<Task587>(tensor587, pindex);
  task582->add_dep(task587);
  task587->add_dep(task519);
  deciq->add_task(task587);

  auto tensor588 = vector<shared_ptr<Tensor>>{I754, f1_, t2};
  auto task588 = make_shared<Task588>(tensor588, pindex);
  task587->add_dep(task588);
  task588->add_dep(task519);
  deciq->add_task(task588);

  vector<IndexRange> I758_index = {active_, closed_, closed_, virt_};
  auto I758 = make_shared<Tensor>(I758_index);
  auto tensor589 = vector<shared_ptr<Tensor>>{I745, t2, I758};
  auto task589 = make_shared<Task589>(tensor589, pindex);
  task582->add_dep(task589);
  task589->add_dep(task519);
  deciq->add_task(task589);

  auto tensor590 = vector<shared_ptr<Tensor>>{I758, f1_, t2};
  auto task590 = make_shared<Task590>(tensor590, pindex);
  task589->add_dep(task590);
  task590->add_dep(task519);
  deciq->add_task(task590);

  vector<IndexRange> I762_index = {active_, virt_, closed_, closed_};
  auto I762 = make_shared<Tensor>(I762_index);
  auto tensor591 = vector<shared_ptr<Tensor>>{I745, t2, I762};
  auto task591 = make_shared<Task591>(tensor591, pindex);
  task582->add_dep(task591);
  task591->add_dep(task519);
  deciq->add_task(task591);

  auto tensor592 = vector<shared_ptr<Tensor>>{I762, f1_, t2};
  auto task592 = make_shared<Task592>(tensor592, pindex);
  task591->add_dep(task592);
  task592->add_dep(task519);
  deciq->add_task(task592);

  vector<IndexRange> I766_index = {active_, closed_, closed_, virt_};
  auto I766 = make_shared<Tensor>(I766_index);
  auto tensor593 = vector<shared_ptr<Tensor>>{I745, t2, I766};
  auto task593 = make_shared<Task593>(tensor593, pindex);
  task582->add_dep(task593);
  task593->add_dep(task519);
  deciq->add_task(task593);

  auto tensor594 = vector<shared_ptr<Tensor>>{I766, f1_, t2};
  auto task594 = make_shared<Task594>(tensor594, pindex);
  task593->add_dep(task594);
  task594->add_dep(task519);
  deciq->add_task(task594);

  vector<IndexRange> I786_index = {virt_, active_};
  auto I786 = make_shared<Tensor>(I786_index);
  auto tensor595 = vector<shared_ptr<Tensor>>{I745, f1_, I786};
  auto task595 = make_shared<Task595>(tensor595, pindex);
  task582->add_dep(task595);
  task595->add_dep(task519);
  deciq->add_task(task595);

  vector<IndexRange> I787_index = {closed_, virt_, closed_, virt_};
  auto I787 = make_shared<Tensor>(I787_index);
  auto tensor596 = vector<shared_ptr<Tensor>>{I786, t2, I787};
  auto task596 = make_shared<Task596>(tensor596, pindex);
  task595->add_dep(task596);
  task596->add_dep(task519);
  deciq->add_task(task596);

  auto tensor597 = vector<shared_ptr<Tensor>>{I787, t2};
  auto task597 = make_shared<Task597>(tensor597, pindex);
  task596->add_dep(task597);
  task597->add_dep(task519);
  deciq->add_task(task597);

  vector<IndexRange> I929_index = {active_, virt_};
  auto I929 = make_shared<Tensor>(I929_index);
  auto tensor598 = vector<shared_ptr<Tensor>>{I745, f1_, I929};
  auto task598 = make_shared<Task598>(tensor598, pindex);
  task582->add_dep(task598);
  task598->add_dep(task519);
  deciq->add_task(task598);

  auto tensor599 = vector<shared_ptr<Tensor>>{I929, t2};
  auto task599 = make_shared<Task599>(tensor599, pindex);
  task598->add_dep(task599);
  task599->add_dep(task519);
  deciq->add_task(task599);

  vector<IndexRange> I933_index = {active_, virt_};
  auto I933 = make_shared<Tensor>(I933_index);
  auto tensor600 = vector<shared_ptr<Tensor>>{I745, f1_, I933};
  auto task600 = make_shared<Task600>(tensor600, pindex);
  task582->add_dep(task600);
  task600->add_dep(task519);
  deciq->add_task(task600);

  auto tensor601 = vector<shared_ptr<Tensor>>{I933, t2};
  auto task601 = make_shared<Task601>(tensor601, pindex);
  task600->add_dep(task601);
  task601->add_dep(task519);
  deciq->add_task(task601);

  vector<IndexRange> I1070_index = {closed_, virt_, closed_, active_};
  auto I1070 = make_shared<Tensor>(I1070_index);
  auto tensor602 = vector<shared_ptr<Tensor>>{I745, t2, I1070};
  auto task602 = make_shared<Task602>(tensor602, pindex);
  task582->add_dep(task602);
  task602->add_dep(task519);
  deciq->add_task(task602);

  auto tensor603 = vector<shared_ptr<Tensor>>{I1070, t2};
  auto task603 = make_shared<Task603>(tensor603, pindex, this->e0_);
  task602->add_dep(task603);
  task603->add_dep(task519);
  deciq->add_task(task603);

  auto tensor604 = vector<shared_ptr<Tensor>>{I745, v2_, t2};
  auto task604 = make_shared<Task604>(tensor604, pindex);
  task582->add_dep(task604);
  task604->add_dep(task519);
  deciq->add_task(task604);

  auto tensor605 = vector<shared_ptr<Tensor>>{I745, v2_, t2};
  auto task605 = make_shared<Task605>(tensor605, pindex);
  task582->add_dep(task605);
  task605->add_dep(task519);
  deciq->add_task(task605);

  auto tensor606 = vector<shared_ptr<Tensor>>{I745, v2_, t2};
  auto task606 = make_shared<Task606>(tensor606, pindex);
  task582->add_dep(task606);
  task606->add_dep(task519);
  deciq->add_task(task606);

  auto tensor607 = vector<shared_ptr<Tensor>>{I745, v2_, t2};
  auto task607 = make_shared<Task607>(tensor607, pindex);
  task582->add_dep(task607);
  task607->add_dep(task519);
  deciq->add_task(task607);

  vector<IndexRange> I769_index = {active_, active_, active_, active_};
  auto I769 = make_shared<Tensor>(I769_index);
  auto tensor608 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I769};
  auto task608 = make_shared<Task608>(tensor608, pindex);
  task608->add_dep(task519);
  deciq->add_task(task608);

  vector<IndexRange> I770_index = {active_, closed_, virt_, active_};
  auto I770 = make_shared<Tensor>(I770_index);
  auto tensor609 = vector<shared_ptr<Tensor>>{I769, t2, I770};
  auto task609 = make_shared<Task609>(tensor609, pindex);
  task608->add_dep(task609);
  task609->add_dep(task519);
  deciq->add_task(task609);

  auto tensor610 = vector<shared_ptr<Tensor>>{I770, f1_, t2};
  auto task610 = make_shared<Task610>(tensor610, pindex);
  task609->add_dep(task610);
  task610->add_dep(task519);
  deciq->add_task(task610);

  auto tensor611 = vector<shared_ptr<Tensor>>{I769, v2_, t2};
  auto task611 = make_shared<Task611>(tensor611, pindex);
  task608->add_dep(task611);
  task611->add_dep(task519);
  deciq->add_task(task611);

  vector<IndexRange> I793_index = {active_, active_, active_, active_, active_, active_};
  auto I793 = make_shared<Tensor>(I793_index);
  auto tensor612 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I793};
  auto task612 = make_shared<Task612>(tensor612, pindex);
  task612->add_dep(task519);
  deciq->add_task(task612);

  vector<IndexRange> I794_index = {active_, closed_, active_, active_};
  auto I794 = make_shared<Tensor>(I794_index);
  auto tensor613 = vector<shared_ptr<Tensor>>{I793, t2, I794};
  auto task613 = make_shared<Task613>(tensor613, pindex);
  task612->add_dep(task613);
  task613->add_dep(task519);
  deciq->add_task(task613);

  auto tensor614 = vector<shared_ptr<Tensor>>{I794, f1_, t2};
  auto task614 = make_shared<Task614>(tensor614, pindex);
  task613->add_dep(task614);
  task614->add_dep(task519);
  deciq->add_task(task614);

  vector<IndexRange> I797_index = {active_, active_, active_, active_};
  auto I797 = make_shared<Tensor>(I797_index);
  auto tensor615 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I797};
  auto task615 = make_shared<Task615>(tensor615, pindex);
  task615->add_dep(task519);
  deciq->add_task(task615);

  vector<IndexRange> I798_index = {virt_, closed_, active_, active_};
  auto I798 = make_shared<Tensor>(I798_index);
  auto tensor616 = vector<shared_ptr<Tensor>>{I797, t2, I798};
  auto task616 = make_shared<Task616>(tensor616, pindex);
  task615->add_dep(task616);
  task616->add_dep(task519);
  deciq->add_task(task616);

  auto tensor617 = vector<shared_ptr<Tensor>>{I798, f1_, t2};
  auto task617 = make_shared<Task617>(tensor617, pindex);
  task616->add_dep(task617);
  task617->add_dep(task519);
  deciq->add_task(task617);

  auto tensor618 = vector<shared_ptr<Tensor>>{I797, v2_, t2};
  auto task618 = make_shared<Task618>(tensor618, pindex);
  task615->add_dep(task618);
  task618->add_dep(task519);
  deciq->add_task(task618);

  vector<IndexRange> I805_index = {active_, active_, active_, active_};
  auto I805 = make_shared<Tensor>(I805_index);
  auto tensor619 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I805, f1_};
  auto task619 = make_shared<Task619>(tensor619, pindex);
  task619->add_dep(task519);
  deciq->add_task(task619);

  auto tensor620 = vector<shared_ptr<Tensor>>{I805, t2};
  auto task620 = make_shared<Task620>(tensor620, pindex);
  task619->add_dep(task620);
  task620->add_dep(task519);
  deciq->add_task(task620);

  vector<IndexRange> I808_index = {active_, active_, active_, active_};
  auto I808 = make_shared<Tensor>(I808_index);
  auto tensor621 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I808};
  auto task621 = make_shared<Task621>(tensor621, pindex);
  task621->add_dep(task519);
  deciq->add_task(task621);

  vector<IndexRange> I809_index = {active_, virt_, active_, closed_};
  auto I809 = make_shared<Tensor>(I809_index);
  auto tensor622 = vector<shared_ptr<Tensor>>{I808, t2, I809};
  auto task622 = make_shared<Task622>(tensor622, pindex);
  task621->add_dep(task622);
  task622->add_dep(task519);
  deciq->add_task(task622);

  auto tensor623 = vector<shared_ptr<Tensor>>{I809, f1_, t2};
  auto task623 = make_shared<Task623>(tensor623, pindex);
  task622->add_dep(task623);
  task623->add_dep(task519);
  deciq->add_task(task623);

  vector<IndexRange> I813_index = {active_, closed_, active_, virt_};
  auto I813 = make_shared<Tensor>(I813_index);
  auto tensor624 = vector<shared_ptr<Tensor>>{I808, t2, I813};
  auto task624 = make_shared<Task624>(tensor624, pindex);
  task621->add_dep(task624);
  task624->add_dep(task519);
  deciq->add_task(task624);

  auto tensor625 = vector<shared_ptr<Tensor>>{I813, f1_, t2};
  auto task625 = make_shared<Task625>(tensor625, pindex);
  task624->add_dep(task625);
  task625->add_dep(task519);
  deciq->add_task(task625);

  vector<IndexRange> I844_index = {active_, virt_, closed_, active_};
  auto I844 = make_shared<Tensor>(I844_index);
  auto tensor626 = vector<shared_ptr<Tensor>>{I808, t2, I844};
  auto task626 = make_shared<Task626>(tensor626, pindex);
  task621->add_dep(task626);
  task626->add_dep(task519);
  deciq->add_task(task626);

  auto tensor627 = vector<shared_ptr<Tensor>>{I844, t2};
  auto task627 = make_shared<Task627>(tensor627, pindex, this->e0_);
  task626->add_dep(task627);
  task627->add_dep(task519);
  deciq->add_task(task627);

  auto tensor628 = vector<shared_ptr<Tensor>>{I844, f1_, t2};
  auto task628 = make_shared<Task628>(tensor628, pindex);
  task626->add_dep(task628);
  task628->add_dep(task519);
  deciq->add_task(task628);

  vector<IndexRange> I987_index = {closed_, virt_, active_, active_};
  auto I987 = make_shared<Tensor>(I987_index);
  auto tensor629 = vector<shared_ptr<Tensor>>{I808, t2, I987};
  auto task629 = make_shared<Task629>(tensor629, pindex);
  task621->add_dep(task629);
  task629->add_dep(task519);
  deciq->add_task(task629);

  auto tensor630 = vector<shared_ptr<Tensor>>{I987, f1_, t2};
  auto task630 = make_shared<Task630>(tensor630, pindex);
  task629->add_dep(task630);
  task630->add_dep(task519);
  deciq->add_task(task630);

  auto tensor631 = vector<shared_ptr<Tensor>>{I808, v2_, t2};
  auto task631 = make_shared<Task631>(tensor631, pindex);
  task621->add_dep(task631);
  task631->add_dep(task519);
  deciq->add_task(task631);

  auto tensor632 = vector<shared_ptr<Tensor>>{I808, v2_, t2};
  auto task632 = make_shared<Task632>(tensor632, pindex);
  task621->add_dep(task632);
  task632->add_dep(task519);
  deciq->add_task(task632);

  vector<IndexRange> I816_index = {active_, active_, active_, active_};
  auto I816 = make_shared<Tensor>(I816_index);
  auto tensor633 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I816, f1_};
  auto task633 = make_shared<Task633>(tensor633, pindex);
  task633->add_dep(task519);
  deciq->add_task(task633);

  auto tensor634 = vector<shared_ptr<Tensor>>{I816, t2};
  auto task634 = make_shared<Task634>(tensor634, pindex);
  task633->add_dep(task634);
  task634->add_dep(task519);
  deciq->add_task(task634);

  vector<IndexRange> I860_index = {active_, virt_, closed_, active_};
  auto I860 = make_shared<Tensor>(I860_index);
  auto tensor635 = vector<shared_ptr<Tensor>>{I816, t2, I860};
  auto task635 = make_shared<Task635>(tensor635, pindex);
  task633->add_dep(task635);
  task635->add_dep(task519);
  deciq->add_task(task635);

  auto tensor636 = vector<shared_ptr<Tensor>>{I860, t2};
  auto task636 = make_shared<Task636>(tensor636, pindex);
  task635->add_dep(task636);
  task636->add_dep(task519);
  deciq->add_task(task636);

  vector<IndexRange> I819_index = {active_, active_, active_, active_};
  auto I819 = make_shared<Tensor>(I819_index);
  auto tensor637 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I819};
  auto task637 = make_shared<Task637>(tensor637, pindex);
  task637->add_dep(task519);
  deciq->add_task(task637);

  vector<IndexRange> I820_index = {active_, virt_, active_, closed_};
  auto I820 = make_shared<Tensor>(I820_index);
  auto tensor638 = vector<shared_ptr<Tensor>>{I819, t2, I820};
  auto task638 = make_shared<Task638>(tensor638, pindex);
  task637->add_dep(task638);
  task638->add_dep(task519);
  deciq->add_task(task638);

  auto tensor639 = vector<shared_ptr<Tensor>>{I820, f1_, t2};
  auto task639 = make_shared<Task639>(tensor639, pindex);
  task638->add_dep(task639);
  task639->add_dep(task519);
  deciq->add_task(task639);

  vector<IndexRange> I824_index = {active_, closed_, active_, virt_};
  auto I824 = make_shared<Tensor>(I824_index);
  auto tensor640 = vector<shared_ptr<Tensor>>{I819, t2, I824};
  auto task640 = make_shared<Task640>(tensor640, pindex);
  task637->add_dep(task640);
  task640->add_dep(task519);
  deciq->add_task(task640);

  auto tensor641 = vector<shared_ptr<Tensor>>{I824, f1_, t2};
  auto task641 = make_shared<Task641>(tensor641, pindex);
  task640->add_dep(task641);
  task641->add_dep(task519);
  deciq->add_task(task641);

  auto tensor642 = vector<shared_ptr<Tensor>>{I824, f1_, t2};
  auto task642 = make_shared<Task642>(tensor642, pindex);
  task640->add_dep(task642);
  task642->add_dep(task519);
  deciq->add_task(task642);

  vector<IndexRange> I840_index = {active_, closed_, virt_, active_};
  auto I840 = make_shared<Tensor>(I840_index);
  auto tensor643 = vector<shared_ptr<Tensor>>{I819, t2, I840};
  auto task643 = make_shared<Task643>(tensor643, pindex);
  task637->add_dep(task643);
  task643->add_dep(task519);
  deciq->add_task(task643);

  auto tensor644 = vector<shared_ptr<Tensor>>{I840, t2};
  auto task644 = make_shared<Task644>(tensor644, pindex, this->e0_);
  task643->add_dep(task644);
  task644->add_dep(task519);
  deciq->add_task(task644);

  auto tensor645 = vector<shared_ptr<Tensor>>{I840, f1_, t2};
  auto task645 = make_shared<Task645>(tensor645, pindex);
  task643->add_dep(task645);
  task645->add_dep(task519);
  deciq->add_task(task645);

  vector<IndexRange> I863_index = {active_, active_, virt_, closed_};
  auto I863 = make_shared<Tensor>(I863_index);
  auto tensor646 = vector<shared_ptr<Tensor>>{I819, t2, I863};
  auto task646 = make_shared<Task646>(tensor646, pindex);
  task637->add_dep(task646);
  task646->add_dep(task519);
  deciq->add_task(task646);

  auto tensor647 = vector<shared_ptr<Tensor>>{I863, f1_, t2};
  auto task647 = make_shared<Task647>(tensor647, pindex);
  task646->add_dep(task647);
  task647->add_dep(task519);
  deciq->add_task(task647);

  vector<IndexRange> I867_index = {active_, active_, closed_, virt_};
  auto I867 = make_shared<Tensor>(I867_index);
  auto tensor648 = vector<shared_ptr<Tensor>>{I819, t2, I867};
  auto task648 = make_shared<Task648>(tensor648, pindex);
  task637->add_dep(task648);
  task648->add_dep(task519);
  deciq->add_task(task648);

  auto tensor649 = vector<shared_ptr<Tensor>>{I867, f1_, t2};
  auto task649 = make_shared<Task649>(tensor649, pindex);
  task648->add_dep(task649);
  task649->add_dep(task519);
  deciq->add_task(task649);

  vector<IndexRange> I874_index = {active_, active_, virt_, closed_};
  auto I874 = make_shared<Tensor>(I874_index);
  auto tensor650 = vector<shared_ptr<Tensor>>{I819, t2, I874};
  auto task650 = make_shared<Task650>(tensor650, pindex);
  task637->add_dep(task650);
  task650->add_dep(task519);
  deciq->add_task(task650);

  auto tensor651 = vector<shared_ptr<Tensor>>{I874, f1_, t2};
  auto task651 = make_shared<Task651>(tensor651, pindex);
  task650->add_dep(task651);
  task651->add_dep(task519);
  deciq->add_task(task651);

  vector<IndexRange> I878_index = {active_, active_, closed_, virt_};
  auto I878 = make_shared<Tensor>(I878_index);
  auto tensor652 = vector<shared_ptr<Tensor>>{I819, t2, I878};
  auto task652 = make_shared<Task652>(tensor652, pindex);
  task637->add_dep(task652);
  task652->add_dep(task519);
  deciq->add_task(task652);

  auto tensor653 = vector<shared_ptr<Tensor>>{I878, f1_, t2};
  auto task653 = make_shared<Task653>(tensor653, pindex);
  task652->add_dep(task653);
  task653->add_dep(task519);
  deciq->add_task(task653);

  vector<IndexRange> I894_index = {active_, closed_, virt_, active_};
  auto I894 = make_shared<Tensor>(I894_index);
  auto tensor654 = vector<shared_ptr<Tensor>>{I819, t2, I894};
  auto task654 = make_shared<Task654>(tensor654, pindex);
  task637->add_dep(task654);
  task654->add_dep(task519);
  deciq->add_task(task654);

  auto tensor655 = vector<shared_ptr<Tensor>>{I894, t2};
  auto task655 = make_shared<Task655>(tensor655, pindex, this->e0_);
  task654->add_dep(task655);
  task655->add_dep(task519);
  deciq->add_task(task655);

  vector<IndexRange> I895_index = {active_, virt_, closed_, virt_};
  auto I895 = make_shared<Tensor>(I895_index);
  auto tensor656 = vector<shared_ptr<Tensor>>{I894, f1_, I895};
  auto task656 = make_shared<Task656>(tensor656, pindex);
  task654->add_dep(task656);
  task656->add_dep(task519);
  deciq->add_task(task656);

  auto tensor657 = vector<shared_ptr<Tensor>>{I895, t2};
  auto task657 = make_shared<Task657>(tensor657, pindex);
  task656->add_dep(task657);
  task657->add_dep(task519);
  deciq->add_task(task657);

  vector<IndexRange> I983_index = {virt_, closed_, active_, active_};
  auto I983 = make_shared<Tensor>(I983_index);
  auto tensor658 = vector<shared_ptr<Tensor>>{I819, t2, I983};
  auto task658 = make_shared<Task658>(tensor658, pindex);
  task637->add_dep(task658);
  task658->add_dep(task519);
  deciq->add_task(task658);

  auto tensor659 = vector<shared_ptr<Tensor>>{I983, f1_, t2};
  auto task659 = make_shared<Task659>(tensor659, pindex);
  task658->add_dep(task659);
  task659->add_dep(task519);
  deciq->add_task(task659);

  vector<IndexRange> I995_index = {closed_, virt_, active_, active_};
  auto I995 = make_shared<Tensor>(I995_index);
  auto tensor660 = vector<shared_ptr<Tensor>>{I819, t2, I995};
  auto task660 = make_shared<Task660>(tensor660, pindex);
  task637->add_dep(task660);
  task660->add_dep(task519);
  deciq->add_task(task660);

  auto tensor661 = vector<shared_ptr<Tensor>>{I995, f1_, t2};
  auto task661 = make_shared<Task661>(tensor661, pindex);
  task660->add_dep(task661);
  task661->add_dep(task519);
  deciq->add_task(task661);

  auto tensor662 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task662 = make_shared<Task662>(tensor662, pindex);
  task637->add_dep(task662);
  task662->add_dep(task519);
  deciq->add_task(task662);

  auto tensor663 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task663 = make_shared<Task663>(tensor663, pindex);
  task637->add_dep(task663);
  task663->add_dep(task519);
  deciq->add_task(task663);

  auto tensor664 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task664 = make_shared<Task664>(tensor664, pindex);
  task637->add_dep(task664);
  task664->add_dep(task519);
  deciq->add_task(task664);

  auto tensor665 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task665 = make_shared<Task665>(tensor665, pindex);
  task637->add_dep(task665);
  task665->add_dep(task519);
  deciq->add_task(task665);

  auto tensor666 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task666 = make_shared<Task666>(tensor666, pindex);
  task637->add_dep(task666);
  task666->add_dep(task519);
  deciq->add_task(task666);

  auto tensor667 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task667 = make_shared<Task667>(tensor667, pindex);
  task637->add_dep(task667);
  task667->add_dep(task519);
  deciq->add_task(task667);

  auto tensor668 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task668 = make_shared<Task668>(tensor668, pindex);
  task637->add_dep(task668);
  task668->add_dep(task519);
  deciq->add_task(task668);

  auto tensor669 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task669 = make_shared<Task669>(tensor669, pindex);
  task637->add_dep(task669);
  task669->add_dep(task519);
  deciq->add_task(task669);

  auto tensor670 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task670 = make_shared<Task670>(tensor670, pindex);
  task637->add_dep(task670);
  task670->add_dep(task519);
  deciq->add_task(task670);

  auto tensor671 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task671 = make_shared<Task671>(tensor671, pindex);
  task637->add_dep(task671);
  task671->add_dep(task519);
  deciq->add_task(task671);

  vector<IndexRange> I827_index = {active_, active_, active_, active_, active_, active_};
  auto I827 = make_shared<Tensor>(I827_index);
  auto tensor672 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I827};
  auto task672 = make_shared<Task672>(tensor672, pindex);
  task672->add_dep(task519);
  deciq->add_task(task672);

  vector<IndexRange> I828_index = {active_, virt_, active_, active_};
  auto I828 = make_shared<Tensor>(I828_index);
  auto tensor673 = vector<shared_ptr<Tensor>>{I827, t2, I828};
  auto task673 = make_shared<Task673>(tensor673, pindex);
  task672->add_dep(task673);
  task673->add_dep(task519);
  deciq->add_task(task673);

  auto tensor674 = vector<shared_ptr<Tensor>>{I828, f1_, t2};
  auto task674 = make_shared<Task674>(tensor674, pindex);
  task673->add_dep(task674);
  task674->add_dep(task519);
  deciq->add_task(task674);

  vector<IndexRange> I831_index = {active_, active_};
  auto I831 = make_shared<Tensor>(I831_index);
  auto tensor675 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I831};
  auto task675 = make_shared<Task675>(tensor675, pindex);
  task675->add_dep(task519);
  deciq->add_task(task675);

  vector<IndexRange> I832_index = {closed_, virt_};
  auto I832 = make_shared<Tensor>(I832_index);
  auto tensor676 = vector<shared_ptr<Tensor>>{I831, t2, I832};
  auto task676 = make_shared<Task676>(tensor676, pindex);
  task675->add_dep(task676);
  task676->add_dep(task519);
  deciq->add_task(task676);

  vector<IndexRange> I833_index = {closed_, virt_, closed_, virt_};
  auto I833 = make_shared<Tensor>(I833_index);
  auto tensor677 = vector<shared_ptr<Tensor>>{I832, f1_, I833};
  auto task677 = make_shared<Task677>(tensor677, pindex);
  task676->add_dep(task677);
  task677->add_dep(task519);
  deciq->add_task(task677);

  auto tensor678 = vector<shared_ptr<Tensor>>{I833, t2};
  auto task678 = make_shared<Task678>(tensor678, pindex);
  task677->add_dep(task678);
  task678->add_dep(task519);
  deciq->add_task(task678);

  vector<IndexRange> I886_index = {closed_, virt_};
  auto I886 = make_shared<Tensor>(I886_index);
  auto tensor679 = vector<shared_ptr<Tensor>>{I831, t2, I886};
  auto task679 = make_shared<Task679>(tensor679, pindex);
  task675->add_dep(task679);
  task679->add_dep(task519);
  deciq->add_task(task679);

  vector<IndexRange> I887_index = {closed_, virt_, closed_, virt_};
  auto I887 = make_shared<Tensor>(I887_index);
  auto tensor680 = vector<shared_ptr<Tensor>>{I886, f1_, I887};
  auto task680 = make_shared<Task680>(tensor680, pindex);
  task679->add_dep(task680);
  task680->add_dep(task519);
  deciq->add_task(task680);

  auto tensor681 = vector<shared_ptr<Tensor>>{I887, t2};
  auto task681 = make_shared<Task681>(tensor681, pindex);
  task680->add_dep(task681);
  task681->add_dep(task519);
  deciq->add_task(task681);

  vector<IndexRange> I937_index = {virt_, closed_};
  auto I937 = make_shared<Tensor>(I937_index);
  auto tensor682 = vector<shared_ptr<Tensor>>{I831, t2, I937};
  auto task682 = make_shared<Task682>(tensor682, pindex);
  task675->add_dep(task682);
  task682->add_dep(task519);
  deciq->add_task(task682);

  auto tensor683 = vector<shared_ptr<Tensor>>{I937, f1_, t2};
  auto task683 = make_shared<Task683>(tensor683, pindex);
  task682->add_dep(task683);
  task683->add_dep(task519);
  deciq->add_task(task683);

  vector<IndexRange> I941_index = {virt_, closed_};
  auto I941 = make_shared<Tensor>(I941_index);
  auto tensor684 = vector<shared_ptr<Tensor>>{I831, t2, I941};
  auto task684 = make_shared<Task684>(tensor684, pindex);
  task675->add_dep(task684);
  task684->add_dep(task519);
  deciq->add_task(task684);

  auto tensor685 = vector<shared_ptr<Tensor>>{I941, f1_, t2};
  auto task685 = make_shared<Task685>(tensor685, pindex);
  task684->add_dep(task685);
  task685->add_dep(task519);
  deciq->add_task(task685);

  vector<IndexRange> I945_index = {virt_, closed_};
  auto I945 = make_shared<Tensor>(I945_index);
  auto tensor686 = vector<shared_ptr<Tensor>>{I831, t2, I945};
  auto task686 = make_shared<Task686>(tensor686, pindex);
  task675->add_dep(task686);
  task686->add_dep(task519);
  deciq->add_task(task686);

  auto tensor687 = vector<shared_ptr<Tensor>>{I945, f1_, t2};
  auto task687 = make_shared<Task687>(tensor687, pindex);
  task686->add_dep(task687);
  task687->add_dep(task519);
  deciq->add_task(task687);

  vector<IndexRange> I949_index = {virt_, closed_};
  auto I949 = make_shared<Tensor>(I949_index);
  auto tensor688 = vector<shared_ptr<Tensor>>{I831, t2, I949};
  auto task688 = make_shared<Task688>(tensor688, pindex);
  task675->add_dep(task688);
  task688->add_dep(task519);
  deciq->add_task(task688);

  auto tensor689 = vector<shared_ptr<Tensor>>{I949, f1_, t2};
  auto task689 = make_shared<Task689>(tensor689, pindex);
  task688->add_dep(task689);
  task689->add_dep(task519);
  deciq->add_task(task689);

  vector<IndexRange> I975_index = {active_, closed_};
  auto I975 = make_shared<Tensor>(I975_index);
  auto tensor690 = vector<shared_ptr<Tensor>>{I831, f1_, I975};
  auto task690 = make_shared<Task690>(tensor690, pindex);
  task675->add_dep(task690);
  task690->add_dep(task519);
  deciq->add_task(task690);

  vector<IndexRange> I976_index = {active_, virt_, closed_, virt_};
  auto I976 = make_shared<Tensor>(I976_index);
  auto tensor691 = vector<shared_ptr<Tensor>>{I975, t2, I976};
  auto task691 = make_shared<Task691>(tensor691, pindex);
  task690->add_dep(task691);
  task691->add_dep(task519);
  deciq->add_task(task691);

  auto tensor692 = vector<shared_ptr<Tensor>>{I976, t2};
  auto task692 = make_shared<Task692>(tensor692, pindex);
  task691->add_dep(task692);
  task692->add_dep(task519);
  deciq->add_task(task692);

  vector<IndexRange> I1007_index = {closed_, active_};
  auto I1007 = make_shared<Tensor>(I1007_index);
  auto tensor693 = vector<shared_ptr<Tensor>>{I831, f1_, I1007};
  auto task693 = make_shared<Task693>(tensor693, pindex);
  task675->add_dep(task693);
  task693->add_dep(task519);
  deciq->add_task(task693);

  vector<IndexRange> I1008_index = {closed_, virt_, closed_, virt_};
  auto I1008 = make_shared<Tensor>(I1008_index);
  auto tensor694 = vector<shared_ptr<Tensor>>{I1007, t2, I1008};
  auto task694 = make_shared<Task694>(tensor694, pindex);
  task693->add_dep(task694);
  task694->add_dep(task519);
  deciq->add_task(task694);

  auto tensor695 = vector<shared_ptr<Tensor>>{I1008, t2};
  auto task695 = make_shared<Task695>(tensor695, pindex);
  task694->add_dep(task695);
  task695->add_dep(task519);
  deciq->add_task(task695);

  vector<IndexRange> I1021_index = {virt_, virt_, active_, closed_};
  auto I1021 = make_shared<Tensor>(I1021_index);
  auto tensor696 = vector<shared_ptr<Tensor>>{I831, t2, I1021};
  auto task696 = make_shared<Task696>(tensor696, pindex);
  task675->add_dep(task696);
  task696->add_dep(task519);
  deciq->add_task(task696);

  auto tensor697 = vector<shared_ptr<Tensor>>{I1021, f1_, t2};
  auto task697 = make_shared<Task697>(tensor697, pindex);
  task696->add_dep(task697);
  task697->add_dep(task519);
  deciq->add_task(task697);

  vector<IndexRange> I1025_index = {virt_, virt_, active_, closed_};
  auto I1025 = make_shared<Tensor>(I1025_index);
  auto tensor698 = vector<shared_ptr<Tensor>>{I831, t2, I1025};
  auto task698 = make_shared<Task698>(tensor698, pindex);
  task675->add_dep(task698);
  task698->add_dep(task519);
  deciq->add_task(task698);

  auto tensor699 = vector<shared_ptr<Tensor>>{I1025, f1_, t2};
  auto task699 = make_shared<Task699>(tensor699, pindex);
  task698->add_dep(task699);
  task699->add_dep(task519);
  deciq->add_task(task699);

  vector<IndexRange> I1029_index = {virt_, closed_, active_, virt_};
  auto I1029 = make_shared<Tensor>(I1029_index);
  auto tensor700 = vector<shared_ptr<Tensor>>{I831, t2, I1029};
  auto task700 = make_shared<Task700>(tensor700, pindex);
  task675->add_dep(task700);
  task700->add_dep(task519);
  deciq->add_task(task700);

  auto tensor701 = vector<shared_ptr<Tensor>>{I1029, f1_, t2};
  auto task701 = make_shared<Task701>(tensor701, pindex);
  task700->add_dep(task701);
  task701->add_dep(task519);
  deciq->add_task(task701);

  vector<IndexRange> I1033_index = {virt_, closed_, active_, virt_};
  auto I1033 = make_shared<Tensor>(I1033_index);
  auto tensor702 = vector<shared_ptr<Tensor>>{I831, t2, I1033};
  auto task702 = make_shared<Task702>(tensor702, pindex);
  task675->add_dep(task702);
  task702->add_dep(task519);
  deciq->add_task(task702);

  auto tensor703 = vector<shared_ptr<Tensor>>{I1033, f1_, t2};
  auto task703 = make_shared<Task703>(tensor703, pindex);
  task702->add_dep(task703);
  task703->add_dep(task519);
  deciq->add_task(task703);

  vector<IndexRange> I1037_index = {closed_, virt_, active_, virt_};
  auto I1037 = make_shared<Tensor>(I1037_index);
  auto tensor704 = vector<shared_ptr<Tensor>>{I831, t2, I1037};
  auto task704 = make_shared<Task704>(tensor704, pindex);
  task675->add_dep(task704);
  task704->add_dep(task519);
  deciq->add_task(task704);

  auto tensor705 = vector<shared_ptr<Tensor>>{I1037, f1_, t2};
  auto task705 = make_shared<Task705>(tensor705, pindex);
  task704->add_dep(task705);
  task705->add_dep(task519);
  deciq->add_task(task705);

  vector<IndexRange> I1041_index = {closed_, virt_, active_, virt_};
  auto I1041 = make_shared<Tensor>(I1041_index);
  auto tensor706 = vector<shared_ptr<Tensor>>{I831, t2, I1041};
  auto task706 = make_shared<Task706>(tensor706, pindex);
  task675->add_dep(task706);
  task706->add_dep(task519);
  deciq->add_task(task706);

  auto tensor707 = vector<shared_ptr<Tensor>>{I1041, f1_, t2};
  auto task707 = make_shared<Task707>(tensor707, pindex);
  task706->add_dep(task707);
  task707->add_dep(task519);
  deciq->add_task(task707);

  vector<IndexRange> I1097_index = {active_, virt_, closed_, virt_};
  auto I1097 = make_shared<Tensor>(I1097_index);
  auto tensor708 = vector<shared_ptr<Tensor>>{I831, t2, I1097};
  auto task708 = make_shared<Task708>(tensor708, pindex);
  task675->add_dep(task708);
  task708->add_dep(task519);
  deciq->add_task(task708);

  auto tensor709 = vector<shared_ptr<Tensor>>{I1097, t2};
  auto task709 = make_shared<Task709>(tensor709, pindex, this->e0_);
  task708->add_dep(task709);
  task709->add_dep(task519);
  deciq->add_task(task709);

  auto tensor710 = vector<shared_ptr<Tensor>>{I831, v2_, t2};
  auto task710 = make_shared<Task710>(tensor710, pindex);
  task675->add_dep(task710);
  task710->add_dep(task519);
  deciq->add_task(task710);

  auto tensor711 = vector<shared_ptr<Tensor>>{I831, v2_, t2};
  auto task711 = make_shared<Task711>(tensor711, pindex);
  task675->add_dep(task711);
  task711->add_dep(task519);
  deciq->add_task(task711);

  auto tensor712 = vector<shared_ptr<Tensor>>{I831, v2_, t2};
  auto task712 = make_shared<Task712>(tensor712, pindex);
  task675->add_dep(task712);
  task712->add_dep(task519);
  deciq->add_task(task712);

  auto tensor713 = vector<shared_ptr<Tensor>>{I831, v2_, t2};
  auto task713 = make_shared<Task713>(tensor713, pindex);
  task675->add_dep(task713);
  task713->add_dep(task519);
  deciq->add_task(task713);

  vector<IndexRange> I1223_index = {active_, closed_, virt_, active_};
  auto I1223 = make_shared<Tensor>(I1223_index);
  auto tensor714 = vector<shared_ptr<Tensor>>{I831, h1_, I1223};
  auto task714 = make_shared<Task714>(tensor714, pindex);
  task675->add_dep(task714);
  task714->add_dep(task519);
  deciq->add_task(task714);

  auto tensor715 = vector<shared_ptr<Tensor>>{I1223, t2};
  auto task715 = make_shared<Task715>(tensor715, pindex);
  task714->add_dep(task715);
  task715->add_dep(task519);
  deciq->add_task(task715);

  vector<IndexRange> I1226_index = {active_, active_, virt_, closed_};
  auto I1226 = make_shared<Tensor>(I1226_index);
  auto tensor716 = vector<shared_ptr<Tensor>>{I831, h1_, I1226};
  auto task716 = make_shared<Task716>(tensor716, pindex);
  task675->add_dep(task716);
  task716->add_dep(task519);
  deciq->add_task(task716);

  auto tensor717 = vector<shared_ptr<Tensor>>{I1226, t2};
  auto task717 = make_shared<Task717>(tensor717, pindex);
  task716->add_dep(task717);
  task717->add_dep(task519);
  deciq->add_task(task717);

  vector<IndexRange> I881_index = {active_, active_, active_, active_, active_, active_};
  auto I881 = make_shared<Tensor>(I881_index);
  auto tensor718 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I881};
  auto task718 = make_shared<Task718>(tensor718, pindex);
  task718->add_dep(task519);
  deciq->add_task(task718);

  vector<IndexRange> I882_index = {active_, active_, virt_, active_};
  auto I882 = make_shared<Tensor>(I882_index);
  auto tensor719 = vector<shared_ptr<Tensor>>{I881, t2, I882};
  auto task719 = make_shared<Task719>(tensor719, pindex);
  task718->add_dep(task719);
  task719->add_dep(task519);
  deciq->add_task(task719);

  auto tensor720 = vector<shared_ptr<Tensor>>{I882, f1_, t2};
  auto task720 = make_shared<Task720>(tensor720, pindex);
  task719->add_dep(task720);
  task720->add_dep(task519);
  deciq->add_task(task720);

  auto tensor721 = vector<shared_ptr<Tensor>>{I881, v2_, t2};
  auto task721 = make_shared<Task721>(tensor721, pindex);
  task718->add_dep(task721);
  task721->add_dep(task519);
  deciq->add_task(task721);

  vector<IndexRange> I901_index = {active_, active_, active_, active_, active_, active_};
  auto I901 = make_shared<Tensor>(I901_index);
  auto tensor722 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I901};
  auto task722 = make_shared<Task722>(tensor722, pindex);
  task722->add_dep(task519);
  deciq->add_task(task722);

  vector<IndexRange> I902_index = {active_, virt_, active_, active_};
  auto I902 = make_shared<Tensor>(I902_index);
  auto tensor723 = vector<shared_ptr<Tensor>>{I901, t2, I902};
  auto task723 = make_shared<Task723>(tensor723, pindex);
  task722->add_dep(task723);
  task723->add_dep(task519);
  deciq->add_task(task723);

  auto tensor724 = vector<shared_ptr<Tensor>>{I902, f1_, t2};
  auto task724 = make_shared<Task724>(tensor724, pindex);
  task723->add_dep(task724);
  task724->add_dep(task519);
  deciq->add_task(task724);

  vector<IndexRange> I905_index = {active_, active_, active_, active_, active_, active_};
  auto I905 = make_shared<Tensor>(I905_index);
  auto tensor725 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I905};
  auto task725 = make_shared<Task725>(tensor725, pindex);
  task725->add_dep(task519);
  deciq->add_task(task725);

  vector<IndexRange> I906_index = {virt_, active_, active_, active_};
  auto I906 = make_shared<Tensor>(I906_index);
  auto tensor726 = vector<shared_ptr<Tensor>>{I905, t2, I906};
  auto task726 = make_shared<Task726>(tensor726, pindex);
  task725->add_dep(task726);
  task726->add_dep(task519);
  deciq->add_task(task726);

  auto tensor727 = vector<shared_ptr<Tensor>>{I906, f1_, t2};
  auto task727 = make_shared<Task727>(tensor727, pindex);
  task726->add_dep(task727);
  task727->add_dep(task519);
  deciq->add_task(task727);

  auto tensor728 = vector<shared_ptr<Tensor>>{I905, v2_, t2};
  auto task728 = make_shared<Task728>(tensor728, pindex);
  task725->add_dep(task728);
  task728->add_dep(task519);
  deciq->add_task(task728);

  vector<IndexRange> I909_index = {active_, active_, active_, active_, active_, active_};
  auto I909 = make_shared<Tensor>(I909_index);
  auto tensor729 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I909, f1_};
  auto task729 = make_shared<Task729>(tensor729, pindex);
  task729->add_dep(task519);
  deciq->add_task(task729);

  auto tensor730 = vector<shared_ptr<Tensor>>{I909, t2};
  auto task730 = make_shared<Task730>(tensor730, pindex);
  task729->add_dep(task730);
  task730->add_dep(task519);
  deciq->add_task(task730);

  vector<IndexRange> I912_index = {active_, active_, active_, active_, active_, active_};
  auto I912 = make_shared<Tensor>(I912_index);
  auto tensor731 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I912};
  auto task731 = make_shared<Task731>(tensor731, pindex);
  task731->add_dep(task519);
  deciq->add_task(task731);

  vector<IndexRange> I913_index = {active_, active_, active_, virt_};
  auto I913 = make_shared<Tensor>(I913_index);
  auto tensor732 = vector<shared_ptr<Tensor>>{I912, t2, I913};
  auto task732 = make_shared<Task732>(tensor732, pindex);
  task731->add_dep(task732);
  task732->add_dep(task519);
  deciq->add_task(task732);

  auto tensor733 = vector<shared_ptr<Tensor>>{I913, f1_, t2};
  auto task733 = make_shared<Task733>(tensor733, pindex);
  task732->add_dep(task733);
  task733->add_dep(task519);
  deciq->add_task(task733);

  vector<IndexRange> I925_index = {active_, virt_, active_, active_};
  auto I925 = make_shared<Tensor>(I925_index);
  auto tensor734 = vector<shared_ptr<Tensor>>{I912, t2, I925};
  auto task734 = make_shared<Task734>(tensor734, pindex);
  task731->add_dep(task734);
  task734->add_dep(task519);
  deciq->add_task(task734);

  auto tensor735 = vector<shared_ptr<Tensor>>{I925, t2};
  auto task735 = make_shared<Task735>(tensor735, pindex, this->e0_);
  task734->add_dep(task735);
  task735->add_dep(task519);
  deciq->add_task(task735);

  auto tensor736 = vector<shared_ptr<Tensor>>{I925, f1_, t2};
  auto task736 = make_shared<Task736>(tensor736, pindex);
  task734->add_dep(task736);
  task736->add_dep(task519);
  deciq->add_task(task736);

  vector<IndexRange> I1049_index = {active_, virt_, active_, active_};
  auto I1049 = make_shared<Tensor>(I1049_index);
  auto tensor737 = vector<shared_ptr<Tensor>>{I912, t2, I1049};
  auto task737 = make_shared<Task737>(tensor737, pindex);
  task731->add_dep(task737);
  task737->add_dep(task519);
  deciq->add_task(task737);

  auto tensor738 = vector<shared_ptr<Tensor>>{I1049, f1_, t2};
  auto task738 = make_shared<Task738>(tensor738, pindex);
  task737->add_dep(task738);
  task738->add_dep(task519);
  deciq->add_task(task738);

  auto tensor739 = vector<shared_ptr<Tensor>>{I912, v2_, t2};
  auto task739 = make_shared<Task739>(tensor739, pindex);
  task731->add_dep(task739);
  task739->add_dep(task519);
  deciq->add_task(task739);

  auto tensor740 = vector<shared_ptr<Tensor>>{I912, v2_, t2};
  auto task740 = make_shared<Task740>(tensor740, pindex);
  task731->add_dep(task740);
  task740->add_dep(task519);
  deciq->add_task(task740);

  vector<IndexRange> I916_index = {active_, active_, active_, active_};
  auto I916 = make_shared<Tensor>(I916_index);
  auto tensor741 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I916};
  auto task741 = make_shared<Task741>(tensor741, pindex);
  task741->add_dep(task519);
  deciq->add_task(task741);

  vector<IndexRange> I917_index = {active_, virt_};
  auto I917 = make_shared<Tensor>(I917_index);
  auto tensor742 = vector<shared_ptr<Tensor>>{I916, t2, I917};
  auto task742 = make_shared<Task742>(tensor742, pindex);
  task741->add_dep(task742);
  task742->add_dep(task519);
  deciq->add_task(task742);

  vector<IndexRange> I918_index = {active_, virt_, closed_, virt_};
  auto I918 = make_shared<Tensor>(I918_index);
  auto tensor743 = vector<shared_ptr<Tensor>>{I917, f1_, I918};
  auto task743 = make_shared<Task743>(tensor743, pindex);
  task742->add_dep(task743);
  task743->add_dep(task519);
  deciq->add_task(task743);

  auto tensor744 = vector<shared_ptr<Tensor>>{I918, t2};
  auto task744 = make_shared<Task744>(tensor744, pindex);
  task743->add_dep(task744);
  task744->add_dep(task519);
  deciq->add_task(task744);

  vector<IndexRange> I999_index = {virt_, active_};
  auto I999 = make_shared<Tensor>(I999_index);
  auto tensor745 = vector<shared_ptr<Tensor>>{I916, t2, I999};
  auto task745 = make_shared<Task745>(tensor745, pindex);
  task741->add_dep(task745);
  task745->add_dep(task519);
  deciq->add_task(task745);

  auto tensor746 = vector<shared_ptr<Tensor>>{I999, f1_, t2};
  auto task746 = make_shared<Task746>(tensor746, pindex);
  task745->add_dep(task746);
  task746->add_dep(task519);
  deciq->add_task(task746);

  vector<IndexRange> I1003_index = {virt_, active_};
  auto I1003 = make_shared<Tensor>(I1003_index);
  auto tensor747 = vector<shared_ptr<Tensor>>{I916, t2, I1003};
  auto task747 = make_shared<Task747>(tensor747, pindex);
  task741->add_dep(task747);
  task747->add_dep(task519);
  deciq->add_task(task747);

  auto tensor748 = vector<shared_ptr<Tensor>>{I1003, f1_, t2};
  auto task748 = make_shared<Task748>(tensor748, pindex);
  task747->add_dep(task748);
  task748->add_dep(task519);
  deciq->add_task(task748);

  vector<IndexRange> I1045_index = {virt_, virt_, active_, active_};
  auto I1045 = make_shared<Tensor>(I1045_index);
  auto tensor749 = vector<shared_ptr<Tensor>>{I916, t2, I1045};
  auto task749 = make_shared<Task749>(tensor749, pindex);
  task741->add_dep(task749);
  task749->add_dep(task519);
  deciq->add_task(task749);

  auto tensor750 = vector<shared_ptr<Tensor>>{I1045, f1_, t2};
  auto task750 = make_shared<Task750>(tensor750, pindex);
  task749->add_dep(task750);
  task750->add_dep(task519);
  deciq->add_task(task750);

  auto tensor751 = vector<shared_ptr<Tensor>>{I1045, f1_, t2};
  auto task751 = make_shared<Task751>(tensor751, pindex);
  task749->add_dep(task751);
  task751->add_dep(task519);
  deciq->add_task(task751);

  vector<IndexRange> I1053_index = {active_, virt_, virt_, active_};
  auto I1053 = make_shared<Tensor>(I1053_index);
  auto tensor752 = vector<shared_ptr<Tensor>>{I916, t2, I1053};
  auto task752 = make_shared<Task752>(tensor752, pindex);
  task741->add_dep(task752);
  task752->add_dep(task519);
  deciq->add_task(task752);

  auto tensor753 = vector<shared_ptr<Tensor>>{I1053, t2};
  auto task753 = make_shared<Task753>(tensor753, pindex, this->e0_);
  task752->add_dep(task753);
  task753->add_dep(task519);
  deciq->add_task(task753);

  auto tensor754 = vector<shared_ptr<Tensor>>{I1053, f1_, t2};
  auto task754 = make_shared<Task754>(tensor754, pindex);
  task752->add_dep(task754);
  task754->add_dep(task519);
  deciq->add_task(task754);

  auto tensor755 = vector<shared_ptr<Tensor>>{I916, v2_, t2};
  auto task755 = make_shared<Task755>(tensor755, pindex);
  task741->add_dep(task755);
  task755->add_dep(task519);
  deciq->add_task(task755);

  auto tensor756 = vector<shared_ptr<Tensor>>{I916, v2_, t2};
  auto task756 = make_shared<Task756>(tensor756, pindex);
  task741->add_dep(task756);
  task756->add_dep(task519);
  deciq->add_task(task756);

  auto tensor757 = vector<shared_ptr<Tensor>>{I916, h1_, t2};
  auto task757 = make_shared<Task757>(tensor757, pindex);
  task741->add_dep(task757);
  task757->add_dep(task519);
  deciq->add_task(task757);

  auto tensor758 = vector<shared_ptr<Tensor>>{I916, h1_, t2};
  auto task758 = make_shared<Task758>(tensor758, pindex);
  task741->add_dep(task758);
  task758->add_dep(task519);
  deciq->add_task(task758);

  vector<IndexRange> I952_index;
  auto I952 = make_shared<Tensor>(I952_index);
  auto tensor759 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I952, f1_};
  auto task759 = make_shared<Task759>(tensor759, pindex);
  task759->add_dep(task519);
  deciq->add_task(task759);

  vector<IndexRange> I953_index = {closed_, virt_, closed_, virt_};
  auto I953 = make_shared<Tensor>(I953_index);
  auto tensor760 = vector<shared_ptr<Tensor>>{I952, t2, I953};
  auto task760 = make_shared<Task760>(tensor760, pindex);
  task759->add_dep(task760);
  task760->add_dep(task519);
  deciq->add_task(task760);

  auto tensor761 = vector<shared_ptr<Tensor>>{I953, t2};
  auto task761 = make_shared<Task761>(tensor761, pindex);
  task760->add_dep(task761);
  task761->add_dep(task519);
  deciq->add_task(task761);

  shared_ptr<Tensor> I958;
  if (diagonal) {
    vector<IndexRange> I958_index;
    I958 = make_shared<Tensor>(I958_index);
  }
  shared_ptr<Task762> task762;
  if (diagonal) {
    auto tensor762 = vector<shared_ptr<Tensor>>{den0ci, I958};
    task762 = make_shared<Task762>(tensor762, pindex);
    task762->add_dep(task519);
    deciq->add_task(task762);
  }

  shared_ptr<Tensor> I959;
  if (diagonal) {
    vector<IndexRange> I959_index = {closed_, closed_};
    I959 = make_shared<Tensor>(I959_index);
  }
  shared_ptr<Task763> task763;
  if (diagonal) {
    auto tensor763 = vector<shared_ptr<Tensor>>{I958, f1_, I959};
    task763 = make_shared<Task763>(tensor763, pindex);
    task762->add_dep(task763);
    task763->add_dep(task519);
    deciq->add_task(task763);
  }

  shared_ptr<Task764> task764;
  if (diagonal) {
    auto tensor764 = vector<shared_ptr<Tensor>>{I959, t2};
    task764 = make_shared<Task764>(tensor764, pindex);
    task763->add_dep(task764);
    task764->add_dep(task519);
    deciq->add_task(task764);
  }

  shared_ptr<Task765> task765;
  if (diagonal) {
    auto tensor765 = vector<shared_ptr<Tensor>>{I958, v2_, t2};
    task765 = make_shared<Task765>(tensor765, pindex);
    task762->add_dep(task765);
    task765->add_dep(task519);
    deciq->add_task(task765);
  }

  shared_ptr<Tensor> I962;
  if (diagonal) {
    vector<IndexRange> I962_index;
    I962 = make_shared<Tensor>(I962_index);
  }
  shared_ptr<Task766> task766;
  if (diagonal) {
    auto tensor766 = vector<shared_ptr<Tensor>>{den0ci, I962};
    task766 = make_shared<Task766>(tensor766, pindex);
    task766->add_dep(task519);
    deciq->add_task(task766);
  }

  shared_ptr<Tensor> I963;
  if (diagonal) {
    vector<IndexRange> I963_index = {closed_, closed_};
    I963 = make_shared<Tensor>(I963_index);
  }
  shared_ptr<Task767> task767;
  if (diagonal) {
    auto tensor767 = vector<shared_ptr<Tensor>>{I962, f1_, I963};
    task767 = make_shared<Task767>(tensor767, pindex);
    task766->add_dep(task767);
    task767->add_dep(task519);
    deciq->add_task(task767);
  }

  shared_ptr<Task768> task768;
  if (diagonal) {
    auto tensor768 = vector<shared_ptr<Tensor>>{I963, t2};
    task768 = make_shared<Task768>(tensor768, pindex);
    task767->add_dep(task768);
    task768->add_dep(task519);
    deciq->add_task(task768);
  }

  shared_ptr<Tensor> I966;
  if (diagonal) {
    vector<IndexRange> I966_index;
    I966 = make_shared<Tensor>(I966_index);
  }
  shared_ptr<Task769> task769;
  if (diagonal) {
    auto tensor769 = vector<shared_ptr<Tensor>>{den0ci, I966};
    task769 = make_shared<Task769>(tensor769, pindex);
    task769->add_dep(task519);
    deciq->add_task(task769);
  }

  shared_ptr<Tensor> I967;
  if (diagonal) {
    vector<IndexRange> I967_index = {virt_, virt_};
    I967 = make_shared<Tensor>(I967_index);
  }
  shared_ptr<Task770> task770;
  if (diagonal) {
    auto tensor770 = vector<shared_ptr<Tensor>>{I966, f1_, I967};
    task770 = make_shared<Task770>(tensor770, pindex);
    task769->add_dep(task770);
    task770->add_dep(task519);
    deciq->add_task(task770);
  }

  shared_ptr<Task771> task771;
  if (diagonal) {
    auto tensor771 = vector<shared_ptr<Tensor>>{I967, t2};
    task771 = make_shared<Task771>(tensor771, pindex);
    task770->add_dep(task771);
    task771->add_dep(task519);
    deciq->add_task(task771);
  }

  shared_ptr<Task772> task772;
  if (diagonal) {
    auto tensor772 = vector<shared_ptr<Tensor>>{I966, t2};
    task772 = make_shared<Task772>(tensor772, pindex, this->e0_);
    task769->add_dep(task772);
    task772->add_dep(task519);
    deciq->add_task(task772);
  }

  shared_ptr<Tensor> I970;
  if (diagonal) {
    vector<IndexRange> I970_index;
    I970 = make_shared<Tensor>(I970_index);
  }
  shared_ptr<Task773> task773;
  if (diagonal) {
    auto tensor773 = vector<shared_ptr<Tensor>>{den0ci, I970};
    task773 = make_shared<Task773>(tensor773, pindex);
    task773->add_dep(task519);
    deciq->add_task(task773);
  }

  shared_ptr<Tensor> I971;
  if (diagonal) {
    vector<IndexRange> I971_index = {virt_, virt_};
    I971 = make_shared<Tensor>(I971_index);
  }
  shared_ptr<Task774> task774;
  if (diagonal) {
    auto tensor774 = vector<shared_ptr<Tensor>>{I970, f1_, I971};
    task774 = make_shared<Task774>(tensor774, pindex);
    task773->add_dep(task774);
    task774->add_dep(task519);
    deciq->add_task(task774);
  }

  shared_ptr<Task775> task775;
  if (diagonal) {
    auto tensor775 = vector<shared_ptr<Tensor>>{I971, t2};
    task775 = make_shared<Task775>(tensor775, pindex);
    task774->add_dep(task775);
    task775->add_dep(task519);
    deciq->add_task(task775);
  }

  vector<IndexRange> I1014_index = {active_, active_};
  auto I1014 = make_shared<Tensor>(I1014_index);
  auto tensor776 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I1014, f1_};
  auto task776 = make_shared<Task776>(tensor776, pindex);
  task776->add_dep(task519);
  deciq->add_task(task776);

  vector<IndexRange> I1015_index = {active_, virt_, closed_, virt_};
  auto I1015 = make_shared<Tensor>(I1015_index);
  auto tensor777 = vector<shared_ptr<Tensor>>{I1014, t2, I1015};
  auto task777 = make_shared<Task777>(tensor777, pindex);
  task776->add_dep(task777);
  task777->add_dep(task519);
  deciq->add_task(task777);

  auto tensor778 = vector<shared_ptr<Tensor>>{I1015, t2};
  auto task778 = make_shared<Task778>(tensor778, pindex);
  task777->add_dep(task778);
  task778->add_dep(task519);
  deciq->add_task(task778);

  vector<IndexRange> I1056_index = {active_, active_, active_, active_};
  auto I1056 = make_shared<Tensor>(I1056_index);
  auto tensor779 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I1056, f1_};
  auto task779 = make_shared<Task779>(tensor779, pindex);
  task779->add_dep(task519);
  deciq->add_task(task779);

  auto tensor780 = vector<shared_ptr<Tensor>>{I1056, t2};
  auto task780 = make_shared<Task780>(tensor780, pindex);
  task779->add_dep(task780);
  task780->add_dep(task519);
  deciq->add_task(task780);

  shared_ptr<Tensor> I1090;
  if (diagonal) {
    vector<IndexRange> I1090_index;
    I1090 = make_shared<Tensor>(I1090_index);
  }
  shared_ptr<Task781> task781;
  if (diagonal) {
    auto tensor781 = vector<shared_ptr<Tensor>>{den0ci, I1090};
    task781 = make_shared<Task781>(tensor781, pindex);
    task781->add_dep(task519);
    deciq->add_task(task781);
  }

  shared_ptr<Task782> task782;
  if (diagonal) {
    auto tensor782 = vector<shared_ptr<Tensor>>{I1090, t2};
    task782 = make_shared<Task782>(tensor782, pindex, this->e0_);
    task781->add_dep(task782);
    task782->add_dep(task519);
    deciq->add_task(task782);
  }

  vector<IndexRange> I1108_index = {active_, active_, active_, active_, active_, active_};
  auto I1108 = make_shared<Tensor>(I1108_index);
  auto tensor783 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I1108};
  auto task783 = make_shared<Task783>(tensor783, pindex);
  task783->add_dep(task519);
  deciq->add_task(task783);

  auto tensor784 = vector<shared_ptr<Tensor>>{I1108, v2_, t2};
  auto task784 = make_shared<Task784>(tensor784, pindex);
  task783->add_dep(task784);
  task784->add_dep(task519);
  deciq->add_task(task784);

  vector<IndexRange> I1162_index = {active_, active_, active_, active_, active_, active_};
  auto I1162 = make_shared<Tensor>(I1162_index);
  auto tensor785 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I1162};
  auto task785 = make_shared<Task785>(tensor785, pindex);
  task785->add_dep(task519);
  deciq->add_task(task785);

  auto tensor786 = vector<shared_ptr<Tensor>>{I1162, v2_, t2};
  auto task786 = make_shared<Task786>(tensor786, pindex);
  task785->add_dep(task786);
  task786->add_dep(task519);
  deciq->add_task(task786);

  shared_ptr<Tensor> I1204;
  if (diagonal) {
    vector<IndexRange> I1204_index;
    I1204 = make_shared<Tensor>(I1204_index);
  }
  shared_ptr<Task787> task787;
  if (diagonal) {
    auto tensor787 = vector<shared_ptr<Tensor>>{den0ci, I1204};
    task787 = make_shared<Task787>(tensor787, pindex);
    task787->add_dep(task519);
    deciq->add_task(task787);
  }

  shared_ptr<Task788> task788;
  if (diagonal) {
    auto tensor788 = vector<shared_ptr<Tensor>>{I1204, v2_, t2};
    task788 = make_shared<Task788>(tensor788, pindex);
    task787->add_dep(task788);
    task788->add_dep(task519);
    deciq->add_task(task788);
  }

  return deciq;
}


#endif
