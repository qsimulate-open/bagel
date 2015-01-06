//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_residualqq.cc
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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_residualq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto residualq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor6 = {r};
  auto task6 = make_shared<Task6>(tensor6);
  residualq->add_task(task6);

  vector<IndexRange> I0_index = {closed_, virt_, closed_, virt_};
  auto I0 = make_shared<Tensor>(I0_index);
  vector<shared_ptr<Tensor>> tensor7 = {r, I0};
  auto task7 = make_shared<Task7>(tensor7, pindex);
  task7->add_dep(task6);
  residualq->add_task(task7);

  vector<IndexRange> I1_index = {closed_, active_};
  auto I1 = make_shared<Tensor>(I1_index);
  vector<shared_ptr<Tensor>> tensor8 = {I0, t2, I1};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  task7->add_dep(task8);
  task8->add_dep(task6);
  residualq->add_task(task8);

  vector<IndexRange> I2_index = {closed_, active_};
  auto I2 = make_shared<Tensor>(I2_index);
  vector<shared_ptr<Tensor>> tensor9 = {I1, Gamma0_(), I2};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  task8->add_dep(task9);
  task9->add_dep(task6);
  residualq->add_task(task9);

  vector<shared_ptr<Tensor>> tensor10 = {I2, f1_};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  task9->add_dep(task10);
  task10->add_dep(task6);
  residualq->add_task(task10);

  vector<IndexRange> I4_index = {closed_, active_};
  auto I4 = make_shared<Tensor>(I4_index);
  vector<shared_ptr<Tensor>> tensor11 = {I0, t2, I4};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  task7->add_dep(task11);
  task11->add_dep(task6);
  residualq->add_task(task11);

  vector<IndexRange> I5_index = {closed_, active_};
  auto I5 = make_shared<Tensor>(I5_index);
  vector<shared_ptr<Tensor>> tensor12 = {I4, Gamma0_(), I5};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  task11->add_dep(task12);
  task12->add_dep(task6);
  residualq->add_task(task12);

  vector<shared_ptr<Tensor>> tensor13 = {I5, f1_};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  task12->add_dep(task13);
  task13->add_dep(task6);
  residualq->add_task(task13);

  vector<IndexRange> I6_index = {active_, closed_, virt_, virt_};
  auto I6 = make_shared<Tensor>(I6_index);
  vector<shared_ptr<Tensor>> tensor14 = {r, I6};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  task14->add_dep(task6);
  residualq->add_task(task14);

  vector<IndexRange> I7_index = {closed_, active_};
  auto I7 = make_shared<Tensor>(I7_index);
  vector<shared_ptr<Tensor>> tensor15 = {I6, t2, I7};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  task14->add_dep(task15);
  task15->add_dep(task6);
  residualq->add_task(task15);

  vector<IndexRange> I8_index = {active_, closed_};
  auto I8 = make_shared<Tensor>(I8_index);
  vector<shared_ptr<Tensor>> tensor16 = {I7, Gamma0_(), I8};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  task15->add_dep(task16);
  task16->add_dep(task6);
  residualq->add_task(task16);

  vector<shared_ptr<Tensor>> tensor17 = {I8, f1_};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  task16->add_dep(task17);
  task17->add_dep(task6);
  residualq->add_task(task17);

  vector<IndexRange> I10_index = {closed_, active_};
  auto I10 = make_shared<Tensor>(I10_index);
  vector<shared_ptr<Tensor>> tensor18 = {I6, t2, I10};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  task14->add_dep(task18);
  task18->add_dep(task6);
  residualq->add_task(task18);

  vector<IndexRange> I11_index = {active_, closed_};
  auto I11 = make_shared<Tensor>(I11_index);
  vector<shared_ptr<Tensor>> tensor19 = {I10, Gamma0_(), I11};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  task18->add_dep(task19);
  task19->add_dep(task6);
  residualq->add_task(task19);

  vector<shared_ptr<Tensor>> tensor20 = {I11, f1_};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  task19->add_dep(task20);
  task20->add_dep(task6);
  residualq->add_task(task20);

  vector<IndexRange> I13_index = {active_, virt_, closed_, virt_};
  auto I13 = make_shared<Tensor>(I13_index);
  vector<shared_ptr<Tensor>> tensor21 = {I6, Gamma4_(), I13};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  task14->add_dep(task21);
  task21->add_dep(task6);
  residualq->add_task(task21);

  vector<shared_ptr<Tensor>> tensor22 = {I13, t2};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  task21->add_dep(task22);
  task22->add_dep(task6);
  residualq->add_task(task22);

  vector<IndexRange> I17_index = {closed_, active_, virt_, virt_};
  auto I17 = make_shared<Tensor>(I17_index);
  vector<shared_ptr<Tensor>> tensor23 = {I6, Gamma0_(), I17};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  task14->add_dep(task23);
  task23->add_dep(task6);
  residualq->add_task(task23);

  vector<shared_ptr<Tensor>> tensor24 = {I17, t2, v2_};
  auto task24 = make_shared<Task24>(tensor24, pindex, this->e0_);
  task23->add_dep(task24);
  task24->add_dep(task6);
  residualq->add_task(task24);

  vector<IndexRange> I18_index = {closed_, closed_};
  auto I18 = make_shared<Tensor>(I18_index);
  vector<shared_ptr<Tensor>> tensor25 = {I17, t2, I18};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  task23->add_dep(task25);
  task25->add_dep(task6);
  residualq->add_task(task25);

  vector<shared_ptr<Tensor>> tensor26 = {I18, f1_};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  task25->add_dep(task26);
  task26->add_dep(task6);
  residualq->add_task(task26);

  vector<IndexRange> I21_index = {closed_, closed_};
  auto I21 = make_shared<Tensor>(I21_index);
  vector<shared_ptr<Tensor>> tensor27 = {I17, t2, I21};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  task23->add_dep(task27);
  task27->add_dep(task6);
  residualq->add_task(task27);

  vector<shared_ptr<Tensor>> tensor28 = {I21, f1_};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  task27->add_dep(task28);
  task28->add_dep(task6);
  residualq->add_task(task28);

  vector<IndexRange> I24_index = {virt_, virt_};
  auto I24 = make_shared<Tensor>(I24_index);
  vector<shared_ptr<Tensor>> tensor29 = {I17, t2, I24};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  task23->add_dep(task29);
  task29->add_dep(task6);
  residualq->add_task(task29);

  vector<shared_ptr<Tensor>> tensor30 = {I24, f1_};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  task29->add_dep(task30);
  task30->add_dep(task6);
  residualq->add_task(task30);

  vector<IndexRange> I27_index = {virt_, virt_};
  auto I27 = make_shared<Tensor>(I27_index);
  vector<shared_ptr<Tensor>> tensor31 = {I17, t2, I27};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  task23->add_dep(task31);
  task31->add_dep(task6);
  residualq->add_task(task31);

  vector<shared_ptr<Tensor>> tensor32 = {I27, f1_};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  task31->add_dep(task32);
  task32->add_dep(task6);
  residualq->add_task(task32);

  vector<IndexRange> I30_index = {virt_, virt_};
  auto I30 = make_shared<Tensor>(I30_index);
  vector<shared_ptr<Tensor>> tensor33 = {I17, t2, I30};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  task23->add_dep(task33);
  task33->add_dep(task6);
  residualq->add_task(task33);

  vector<shared_ptr<Tensor>> tensor34 = {I30, f1_};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  task33->add_dep(task34);
  task34->add_dep(task6);
  residualq->add_task(task34);

  vector<IndexRange> I33_index = {virt_, virt_};
  auto I33 = make_shared<Tensor>(I33_index);
  vector<shared_ptr<Tensor>> tensor35 = {I17, t2, I33};
  auto task35 = make_shared<Task35>(tensor35, pindex);
  task23->add_dep(task35);
  task35->add_dep(task6);
  residualq->add_task(task35);

  vector<shared_ptr<Tensor>> tensor36 = {I33, f1_};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  task35->add_dep(task36);
  task36->add_dep(task6);
  residualq->add_task(task36);

  vector<IndexRange> I38_index = {closed_, virt_, closed_, virt_};
  auto I38 = make_shared<Tensor>(I38_index);
  vector<shared_ptr<Tensor>> tensor37 = {r, I38};
  auto task37 = make_shared<Task37>(tensor37, pindex);
  task37->add_dep(task6);
  residualq->add_task(task37);

  vector<shared_ptr<Tensor>> tensor38 = {I38, v2_};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  task37->add_dep(task38);
  task38->add_dep(task6);
  residualq->add_task(task38);

  return residualq;
}


