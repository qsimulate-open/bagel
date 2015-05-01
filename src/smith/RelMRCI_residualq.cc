//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_residualqq.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/RelMRCI.h>
#include <src/smith/RelMRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelMRCI::RelMRCI::make_residualq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto residualq = make_shared<Queue>();
  auto tensor4 = vector<shared_ptr<Tensor>>{r};
  auto task4 = make_shared<Task4>(tensor4, reset);
  residualq->add_task(task4);

  vector<IndexRange> I0_index = {closed_, closed_, virt_, virt_};
  auto I0 = make_shared<Tensor>(I0_index);
  auto tensor5 = vector<shared_ptr<Tensor>>{r, I0};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  task5->add_dep(task4);
  residualq->add_task(task5);

  vector<IndexRange> I1_index = {closed_, closed_};
  auto I1 = make_shared<Tensor>(I1_index);
  auto tensor6 = vector<shared_ptr<Tensor>>{I0, t2, I1};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  task5->add_dep(task6);
  task6->add_dep(task4);
  residualq->add_task(task6);

  shared_ptr<Task7> task7;
  if (diagonal) {
    auto tensor7 = vector<shared_ptr<Tensor>>{I1, h1_};
    task7 = make_shared<Task7>(tensor7, pindex);
    task6->add_dep(task7);
    task7->add_dep(task4);
    residualq->add_task(task7);
  }

  auto tensor8 = vector<shared_ptr<Tensor>>{I1, v2_, Gamma0_()};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  task6->add_dep(task8);
  task8->add_dep(task4);
  residualq->add_task(task8);

  auto tensor9 = vector<shared_ptr<Tensor>>{I1, Gamma4_(), v2_};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  task6->add_dep(task9);
  task9->add_dep(task4);
  residualq->add_task(task9);

  auto tensor10 = vector<shared_ptr<Tensor>>{I1, v2_, Gamma0_()};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  task6->add_dep(task10);
  task10->add_dep(task4);
  residualq->add_task(task10);

  auto tensor11 = vector<shared_ptr<Tensor>>{I1, v2_, Gamma0_()};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  task6->add_dep(task11);
  task11->add_dep(task4);
  residualq->add_task(task11);

  vector<IndexRange> I3_index = {closed_, closed_};
  auto I3 = make_shared<Tensor>(I3_index);
  auto tensor12 = vector<shared_ptr<Tensor>>{I0, t2, I3};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  task5->add_dep(task12);
  task12->add_dep(task4);
  residualq->add_task(task12);

  shared_ptr<Task13> task13;
  if (diagonal) {
    auto tensor13 = vector<shared_ptr<Tensor>>{I3, h1_};
    task13 = make_shared<Task13>(tensor13, pindex);
    task12->add_dep(task13);
    task13->add_dep(task4);
    residualq->add_task(task13);
  }

  auto tensor14 = vector<shared_ptr<Tensor>>{I3, v2_, Gamma0_()};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  task12->add_dep(task14);
  task14->add_dep(task4);
  residualq->add_task(task14);

  auto tensor15 = vector<shared_ptr<Tensor>>{I3, Gamma4_(), v2_};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  task12->add_dep(task15);
  task15->add_dep(task4);
  residualq->add_task(task15);

  vector<IndexRange> I37_index = {active_, closed_, closed_, active_};
  auto I37 = make_shared<Tensor>(I37_index);
  auto tensor16 = vector<shared_ptr<Tensor>>{I3, Gamma0_(), I37};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  task12->add_dep(task16);
  task16->add_dep(task4);
  residualq->add_task(task16);

  auto tensor17 = vector<shared_ptr<Tensor>>{I37, v2_};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  task16->add_dep(task17);
  task17->add_dep(task4);
  residualq->add_task(task17);

  vector<IndexRange> I5_index = {virt_, virt_};
  auto I5 = make_shared<Tensor>(I5_index);
  auto tensor18 = vector<shared_ptr<Tensor>>{I0, t2, I5};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  task5->add_dep(task18);
  task18->add_dep(task4);
  residualq->add_task(task18);

  shared_ptr<Task19> task19;
  if (diagonal) {
    auto tensor19 = vector<shared_ptr<Tensor>>{I5, h1_};
    task19 = make_shared<Task19>(tensor19, pindex);
    task18->add_dep(task19);
    task19->add_dep(task4);
    residualq->add_task(task19);
  }

  auto tensor20 = vector<shared_ptr<Tensor>>{I5, v2_, Gamma0_()};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  task18->add_dep(task20);
  task20->add_dep(task4);
  residualq->add_task(task20);

  auto tensor21 = vector<shared_ptr<Tensor>>{I5, v2_, Gamma4_()};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  task18->add_dep(task21);
  task21->add_dep(task4);
  residualq->add_task(task21);

  vector<IndexRange> I40_index = {active_, virt_, virt_, active_};
  auto I40 = make_shared<Tensor>(I40_index);
  auto tensor22 = vector<shared_ptr<Tensor>>{I5, Gamma0_(), I40};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  task18->add_dep(task22);
  task22->add_dep(task4);
  residualq->add_task(task22);

  auto tensor23 = vector<shared_ptr<Tensor>>{I40, v2_};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  task22->add_dep(task23);
  task23->add_dep(task4);
  residualq->add_task(task23);

  shared_ptr<Task24> task24;
  if (diagonal) {
    auto tensor24 = vector<shared_ptr<Tensor>>{I0, h1_, t2};
    task24 = make_shared<Task24>(tensor24, pindex);
    task5->add_dep(task24);
    task24->add_dep(task4);
    residualq->add_task(task24);
  }

  vector<IndexRange> I18_index = {virt_, virt_};
  auto I18 = make_shared<Tensor>(I18_index);
  auto tensor25 = vector<shared_ptr<Tensor>>{I0, t2, I18};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  task5->add_dep(task25);
  task25->add_dep(task4);
  residualq->add_task(task25);

  vector<IndexRange> I19_index = {virt_, virt_, active_, active_};
  auto I19 = make_shared<Tensor>(I19_index);
  auto tensor26 = vector<shared_ptr<Tensor>>{I18, Gamma0_(), I19};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  task25->add_dep(task26);
  task26->add_dep(task4);
  residualq->add_task(task26);

  auto tensor27 = vector<shared_ptr<Tensor>>{I19, v2_};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  task26->add_dep(task27);
  task27->add_dep(task4);
  residualq->add_task(task27);

  auto tensor28 = vector<shared_ptr<Tensor>>{I18, Gamma4_(), v2_};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  task25->add_dep(task28);
  task28->add_dep(task4);
  residualq->add_task(task28);

  shared_ptr<Task29> task29;
  if (diagonal) {
    auto tensor29 = vector<shared_ptr<Tensor>>{I0, t2, v2_};
    task29 = make_shared<Task29>(tensor29, pindex);
    task5->add_dep(task29);
    task29->add_dep(task4);
    residualq->add_task(task29);
  }

  shared_ptr<Task30> task30;
  if (diagonal) {
    auto tensor30 = vector<shared_ptr<Tensor>>{I0, t2, v2_};
    task30 = make_shared<Task30>(tensor30, pindex);
    task5->add_dep(task30);
    task30->add_dep(task4);
    residualq->add_task(task30);
  }

  shared_ptr<Task31> task31;
  if (diagonal) {
    auto tensor31 = vector<shared_ptr<Tensor>>{I0, v2_, t2};
    task31 = make_shared<Task31>(tensor31, pindex);
    task5->add_dep(task31);
    task31->add_dep(task4);
    residualq->add_task(task31);
  }

  shared_ptr<Task32> task32;
  if (diagonal) {
    auto tensor32 = vector<shared_ptr<Tensor>>{I0, t2, v2_};
    task32 = make_shared<Task32>(tensor32, pindex);
    task5->add_dep(task32);
    task32->add_dep(task4);
    residualq->add_task(task32);
  }

  shared_ptr<Task33> task33;
  if (diagonal) {
    auto tensor33 = vector<shared_ptr<Tensor>>{I0, t2, v2_};
    task33 = make_shared<Task33>(tensor33, pindex);
    task5->add_dep(task33);
    task33->add_dep(task4);
    residualq->add_task(task33);
  }

  shared_ptr<Task34> task34;
  if (diagonal) {
    auto tensor34 = vector<shared_ptr<Tensor>>{I0, t2, v2_};
    task34 = make_shared<Task34>(tensor34, pindex);
    task5->add_dep(task34);
    task34->add_dep(task4);
    residualq->add_task(task34);
  }

  shared_ptr<Task35> task35;
  if (diagonal) {
    auto tensor35 = vector<shared_ptr<Tensor>>{I0, t2, v2_};
    task35 = make_shared<Task35>(tensor35, pindex);
    task5->add_dep(task35);
    task35->add_dep(task4);
    residualq->add_task(task35);
  }

  shared_ptr<Task36> task36;
  if (diagonal) {
    auto tensor36 = vector<shared_ptr<Tensor>>{I0, t2, v2_};
    task36 = make_shared<Task36>(tensor36, pindex);
    task5->add_dep(task36);
    task36->add_dep(task4);
    residualq->add_task(task36);
  }

  vector<IndexRange> I56_index = {closed_, closed_, virt_, virt_};
  auto I56 = make_shared<Tensor>(I56_index);
  auto tensor37 = vector<shared_ptr<Tensor>>{r, I56};
  auto task37 = make_shared<Task37>(tensor37, pindex);
  task37->add_dep(task4);
  residualq->add_task(task37);

  shared_ptr<Task38> task38;
  if (diagonal) {
    auto tensor38 = vector<shared_ptr<Tensor>>{I56, t2, v2_};
    task38 = make_shared<Task38>(tensor38, pindex);
    task37->add_dep(task38);
    task38->add_dep(task4);
    residualq->add_task(task38);
  }

  shared_ptr<Task39> task39;
  if (diagonal) {
    auto tensor39 = vector<shared_ptr<Tensor>>{I56, t2, v2_};
    task39 = make_shared<Task39>(tensor39, pindex);
    task37->add_dep(task39);
    task39->add_dep(task4);
    residualq->add_task(task39);
  }

  shared_ptr<Task40> task40;
  if (diagonal) {
    auto tensor40 = vector<shared_ptr<Tensor>>{I56, t2, v2_};
    task40 = make_shared<Task40>(tensor40, pindex);
    task37->add_dep(task40);
    task40->add_dep(task4);
    residualq->add_task(task40);
  }

  shared_ptr<Task41> task41;
  if (diagonal) {
    auto tensor41 = vector<shared_ptr<Tensor>>{I56, t2, v2_};
    task41 = make_shared<Task41>(tensor41, pindex);
    task37->add_dep(task41);
    task41->add_dep(task4);
    residualq->add_task(task41);
  }

  vector<IndexRange> I81_index = {closed_, virt_, closed_, virt_};
  auto I81 = make_shared<Tensor>(I81_index);
  auto tensor42 = vector<shared_ptr<Tensor>>{I56, Gamma16_(), I81};
  auto task42 = make_shared<Task42>(tensor42, pindex);
  task37->add_dep(task42);
  task42->add_dep(task4);
  residualq->add_task(task42);

  auto tensor43 = vector<shared_ptr<Tensor>>{I81, t2};
  auto task43 = make_shared<Task43>(tensor43, pindex);
  task42->add_dep(task43);
  task43->add_dep(task4);
  residualq->add_task(task43);

  vector<IndexRange> I85_index = {closed_, virt_, closed_, virt_};
  auto I85 = make_shared<Tensor>(I85_index);
  auto tensor44 = vector<shared_ptr<Tensor>>{I56, Gamma18_(), I85};
  auto task44 = make_shared<Task44>(tensor44, pindex);
  task37->add_dep(task44);
  task44->add_dep(task4);
  residualq->add_task(task44);

  auto tensor45 = vector<shared_ptr<Tensor>>{I85, t2};
  auto task45 = make_shared<Task45>(tensor45, pindex);
  task44->add_dep(task45);
  task45->add_dep(task4);
  residualq->add_task(task45);

  return residualq;
}


#endif
