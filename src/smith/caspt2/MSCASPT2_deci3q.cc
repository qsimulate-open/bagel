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
#include <src/smith/caspt2/MSCASPT2_tasks12.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_deci3q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deci3q = make_shared<Queue>();
  auto tensor557 = vector<shared_ptr<Tensor>>{deci};
  auto task557 = make_shared<Task557>(tensor557, reset);
  deci3q->add_task(task557);

  vector<IndexRange> I816_index = {ci_};
  auto I816 = make_shared<Tensor>(I816_index);
  auto tensor558 = vector<shared_ptr<Tensor>>{deci, I816};
  auto task558 = make_shared<Task558>(tensor558, cindex);
  task558->add_dep(task557);
  deci3q->add_task(task558);

  vector<IndexRange> I817_index = {active_, active_, active_, active_};
  auto I817 = make_shared<Tensor>(I817_index);
  auto tensor559 = vector<shared_ptr<Tensor>>{I816, Gamma111_(), I817};
  auto task559 = make_shared<Task559>(tensor559, cindex);
  task558->add_dep(task559);
  task559->add_dep(task557);
  deci3q->add_task(task559);

  auto tensor560 = vector<shared_ptr<Tensor>>{I817, v2_, l2};
  auto task560 = make_shared<Task560>(tensor560, cindex);
  task559->add_dep(task560);
  task560->add_dep(task557);
  deci3q->add_task(task560);

  vector<IndexRange> I820_index = {active_, active_, active_, active_, active_, active_};
  auto I820 = make_shared<Tensor>(I820_index);
  auto tensor561 = vector<shared_ptr<Tensor>>{I816, Gamma239_(), I820};
  auto task561 = make_shared<Task561>(tensor561, cindex);
  task558->add_dep(task561);
  task561->add_dep(task557);
  deci3q->add_task(task561);

  auto tensor562 = vector<shared_ptr<Tensor>>{I820, v2_, l2};
  auto task562 = make_shared<Task562>(tensor562, cindex);
  task561->add_dep(task562);
  task562->add_dep(task557);
  deci3q->add_task(task562);

  vector<IndexRange> I823_index = {active_, active_, active_, active_, active_, active_};
  auto I823 = make_shared<Tensor>(I823_index);
  auto tensor563 = vector<shared_ptr<Tensor>>{I816, Gamma116_(), I823};
  auto task563 = make_shared<Task563>(tensor563, cindex);
  task558->add_dep(task563);
  task563->add_dep(task557);
  deci3q->add_task(task563);

  auto tensor564 = vector<shared_ptr<Tensor>>{I823, v2_, l2};
  auto task564 = make_shared<Task564>(tensor564, cindex);
  task563->add_dep(task564);
  task564->add_dep(task557);
  deci3q->add_task(task564);

  vector<IndexRange> I826_index = {active_, active_};
  auto I826 = make_shared<Tensor>(I826_index);
  auto tensor565 = vector<shared_ptr<Tensor>>{I816, Gamma126_(), I826};
  auto task565 = make_shared<Task565>(tensor565, cindex);
  task558->add_dep(task565);
  task565->add_dep(task557);
  deci3q->add_task(task565);

  auto tensor566 = vector<shared_ptr<Tensor>>{I826, v2_, l2};
  auto task566 = make_shared<Task566>(tensor566, cindex);
  task565->add_dep(task566);
  task566->add_dep(task557);
  deci3q->add_task(task566);

  auto tensor567 = vector<shared_ptr<Tensor>>{I826, v2_, l2};
  auto task567 = make_shared<Task567>(tensor567, cindex);
  task565->add_dep(task567);
  task567->add_dep(task557);
  deci3q->add_task(task567);

  vector<IndexRange> I832_index = {active_, active_, active_, active_};
  auto I832 = make_shared<Tensor>(I832_index);
  auto tensor568 = vector<shared_ptr<Tensor>>{I816, Gamma145_(), I832};
  auto task568 = make_shared<Task568>(tensor568, cindex);
  task558->add_dep(task568);
  task568->add_dep(task557);
  deci3q->add_task(task568);

  auto tensor569 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task569 = make_shared<Task569>(tensor569, cindex);
  task568->add_dep(task569);
  task569->add_dep(task557);
  deci3q->add_task(task569);

  auto tensor570 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task570 = make_shared<Task570>(tensor570, cindex);
  task568->add_dep(task570);
  task570->add_dep(task557);
  deci3q->add_task(task570);

  auto tensor571 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task571 = make_shared<Task571>(tensor571, cindex);
  task568->add_dep(task571);
  task571->add_dep(task557);
  deci3q->add_task(task571);

  auto tensor572 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task572 = make_shared<Task572>(tensor572, cindex);
  task568->add_dep(task572);
  task572->add_dep(task557);
  deci3q->add_task(task572);

  auto tensor573 = vector<shared_ptr<Tensor>>{I832, v2_, l2};
  auto task573 = make_shared<Task573>(tensor573, cindex);
  task568->add_dep(task573);
  task573->add_dep(task557);
  deci3q->add_task(task573);

  vector<IndexRange> I835_index = {active_, active_, active_, active_};
  auto I835 = make_shared<Tensor>(I835_index);
  auto tensor574 = vector<shared_ptr<Tensor>>{I816, Gamma132_(), I835};
  auto task574 = make_shared<Task574>(tensor574, cindex);
  task558->add_dep(task574);
  task574->add_dep(task557);
  deci3q->add_task(task574);

  auto tensor575 = vector<shared_ptr<Tensor>>{I835, v2_, l2};
  auto task575 = make_shared<Task575>(tensor575, cindex);
  task574->add_dep(task575);
  task575->add_dep(task557);
  deci3q->add_task(task575);

  vector<IndexRange> I838_index = {active_, active_, active_, active_};
  auto I838 = make_shared<Tensor>(I838_index);
  auto tensor576 = vector<shared_ptr<Tensor>>{I816, Gamma142_(), I838};
  auto task576 = make_shared<Task576>(tensor576, cindex);
  task558->add_dep(task576);
  task576->add_dep(task557);
  deci3q->add_task(task576);

  auto tensor577 = vector<shared_ptr<Tensor>>{I838, v2_, l2};
  auto task577 = make_shared<Task577>(tensor577, cindex);
  task576->add_dep(task577);
  task577->add_dep(task557);
  deci3q->add_task(task577);

  vector<IndexRange> I847_index = {active_, active_, active_, active_};
  auto I847 = make_shared<Tensor>(I847_index);
  auto tensor578 = vector<shared_ptr<Tensor>>{I816, Gamma122_(), I847};
  auto task578 = make_shared<Task578>(tensor578, cindex);
  task558->add_dep(task578);
  task578->add_dep(task557);
  deci3q->add_task(task578);

  auto tensor579 = vector<shared_ptr<Tensor>>{I847, v2_, l2};
  auto task579 = make_shared<Task579>(tensor579, cindex);
  task578->add_dep(task579);
  task579->add_dep(task557);
  deci3q->add_task(task579);

  auto tensor580 = vector<shared_ptr<Tensor>>{I847, h1_, l2};
  auto task580 = make_shared<Task580>(tensor580, cindex);
  task578->add_dep(task580);
  task580->add_dep(task557);
  deci3q->add_task(task580);

  vector<IndexRange> I856_index = {active_, active_, active_, active_, active_, active_};
  auto I856 = make_shared<Tensor>(I856_index);
  auto tensor581 = vector<shared_ptr<Tensor>>{I816, Gamma169_(), I856};
  auto task581 = make_shared<Task581>(tensor581, cindex);
  task558->add_dep(task581);
  task581->add_dep(task557);
  deci3q->add_task(task581);

  auto tensor582 = vector<shared_ptr<Tensor>>{I856, v2_, l2};
  auto task582 = make_shared<Task582>(tensor582, cindex);
  task581->add_dep(task582);
  task582->add_dep(task557);
  deci3q->add_task(task582);

  vector<IndexRange> I859_index = {active_, active_, active_, active_, active_, active_};
  auto I859 = make_shared<Tensor>(I859_index);
  auto tensor583 = vector<shared_ptr<Tensor>>{I816, Gamma161_(), I859};
  auto task583 = make_shared<Task583>(tensor583, cindex);
  task558->add_dep(task583);
  task583->add_dep(task557);
  deci3q->add_task(task583);

  auto tensor584 = vector<shared_ptr<Tensor>>{I859, v2_, l2};
  auto task584 = make_shared<Task584>(tensor584, cindex);
  task583->add_dep(task584);
  task584->add_dep(task557);
  deci3q->add_task(task584);

  shared_ptr<Tensor> I862;
  if (diagonal) {
    vector<IndexRange> I862_index;
    I862 = make_shared<Tensor>(I862_index);
  }
  shared_ptr<Task585> task585;
  if (diagonal) {
    auto tensor585 = vector<shared_ptr<Tensor>>{I816, rdm0deriv_, I862};
    task585 = make_shared<Task585>(tensor585, cindex);
    task558->add_dep(task585);
    task585->add_dep(task557);
    deci3q->add_task(task585);
  }

  shared_ptr<Task586> task586;
  if (diagonal) {
    auto tensor586 = vector<shared_ptr<Tensor>>{I862, v2_, l2};
    task586 = make_shared<Task586>(tensor586, cindex);
    task585->add_dep(task586);
    task586->add_dep(task557);
    deci3q->add_task(task586);
  }

  shared_ptr<Tensor> I865;
  if (diagonal) {
    vector<IndexRange> I865_index;
    I865 = make_shared<Tensor>(I865_index);
  }
  shared_ptr<Task587> task587;
  if (diagonal) {
    auto tensor587 = vector<shared_ptr<Tensor>>{I816, rdm0deriv_, I865};
    task587 = make_shared<Task587>(tensor587, cindex);
    task558->add_dep(task587);
    task587->add_dep(task557);
    deci3q->add_task(task587);
  }

  shared_ptr<Task588> task588;
  if (diagonal) {
    auto tensor588 = vector<shared_ptr<Tensor>>{I865, v2_, l2};
    task588 = make_shared<Task588>(tensor588, cindex);
    task587->add_dep(task588);
    task588->add_dep(task557);
    deci3q->add_task(task588);
  }

  vector<IndexRange> I868_index = {active_, active_};
  auto I868 = make_shared<Tensor>(I868_index);
  auto tensor589 = vector<shared_ptr<Tensor>>{I816, Gamma148_(), I868};
  auto task589 = make_shared<Task589>(tensor589, cindex);
  task558->add_dep(task589);
  task589->add_dep(task557);
  deci3q->add_task(task589);

  auto tensor590 = vector<shared_ptr<Tensor>>{I868, v2_, l2};
  auto task590 = make_shared<Task590>(tensor590, cindex);
  task589->add_dep(task590);
  task590->add_dep(task557);
  deci3q->add_task(task590);

  auto tensor591 = vector<shared_ptr<Tensor>>{I868, v2_, l2};
  auto task591 = make_shared<Task591>(tensor591, cindex);
  task589->add_dep(task591);
  task591->add_dep(task557);
  deci3q->add_task(task591);

  auto tensor592 = vector<shared_ptr<Tensor>>{I868, h1_, l2};
  auto task592 = make_shared<Task592>(tensor592, cindex);
  task589->add_dep(task592);
  task592->add_dep(task557);
  deci3q->add_task(task592);

  auto tensor593 = vector<shared_ptr<Tensor>>{I868, h1_, l2};
  auto task593 = make_shared<Task593>(tensor593, cindex);
  task589->add_dep(task593);
  task593->add_dep(task557);
  deci3q->add_task(task593);

  vector<IndexRange> I874_index = {active_, active_, active_, active_};
  auto I874 = make_shared<Tensor>(I874_index);
  auto tensor594 = vector<shared_ptr<Tensor>>{I816, Gamma170_(), I874};
  auto task594 = make_shared<Task594>(tensor594, cindex);
  task558->add_dep(task594);
  task594->add_dep(task557);
  deci3q->add_task(task594);

  auto tensor595 = vector<shared_ptr<Tensor>>{I874, v2_, l2};
  auto task595 = make_shared<Task595>(tensor595, cindex);
  task594->add_dep(task595);
  task595->add_dep(task557);
  deci3q->add_task(task595);

  auto tensor596 = vector<shared_ptr<Tensor>>{I874, h1_, l2};
  auto task596 = make_shared<Task596>(tensor596, cindex);
  task594->add_dep(task596);
  task596->add_dep(task557);
  deci3q->add_task(task596);

  return deci3q;
}


#endif
