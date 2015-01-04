//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_energyqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_energyq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto energyq = make_shared<Queue>();
  vector<IndexRange> I16_index;
  auto I16 = make_shared<Tensor>(I16_index);
  vector<IndexRange> I17_index;
  auto I17 = make_shared<Tensor>(I17_index);
  vector<shared_ptr<Tensor>> tensor19 = {I16, Gamma0_(), I17};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  energyq->add_task(task19);

  vector<IndexRange> I18_index = {virt_, closed_, virt_, closed_};
  auto I18 = make_shared<Tensor>(I18_index);
  vector<shared_ptr<Tensor>> tensor20 = {I17, t2, I18};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  task19->add_dep(task20);
  energyq->add_task(task20);

  vector<shared_ptr<Tensor>> tensor21 = {I18, t2};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  task20->add_dep(task21);
  energyq->add_task(task21);

  vector<IndexRange> I21_index = {closed_, virt_, closed_, virt_};
  auto I21 = make_shared<Tensor>(I21_index);
  vector<shared_ptr<Tensor>> tensor22 = {I17, t2, I21};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  task19->add_dep(task22);
  energyq->add_task(task22);

  vector<shared_ptr<Tensor>> tensor23 = {I21, t2};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  task22->add_dep(task23);
  energyq->add_task(task23);

  vector<IndexRange> I23_index = {closed_, closed_};
  auto I23 = make_shared<Tensor>(I23_index);
  vector<shared_ptr<Tensor>> tensor24 = {I16, f1_, I23};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  task19->add_dep(task24);
  energyq->add_task(task24);

  vector<IndexRange> I24_index = {virt_, closed_, virt_, closed_};
  auto I24 = make_shared<Tensor>(I24_index);
  vector<shared_ptr<Tensor>> tensor25 = {I23, t2, I24};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  task24->add_dep(task25);
  energyq->add_task(task25);

  vector<shared_ptr<Tensor>> tensor26 = {I24, t2};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  task25->add_dep(task26);
  energyq->add_task(task26);

  vector<IndexRange> I27_index = {virt_, closed_, virt_, closed_};
  auto I27 = make_shared<Tensor>(I27_index);
  vector<shared_ptr<Tensor>> tensor27 = {I23, t2, I27};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  task24->add_dep(task27);
  energyq->add_task(task27);

  vector<shared_ptr<Tensor>> tensor28 = {I27, t2};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  task27->add_dep(task28);
  energyq->add_task(task28);

  vector<IndexRange> I29_index = {virt_, virt_};
  auto I29 = make_shared<Tensor>(I29_index);
  vector<shared_ptr<Tensor>> tensor29 = {I16, f1_, I29};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  task19->add_dep(task29);
  energyq->add_task(task29);

  vector<IndexRange> I30_index = {virt_, closed_, virt_, closed_};
  auto I30 = make_shared<Tensor>(I30_index);
  vector<shared_ptr<Tensor>> tensor30 = {I29, t2, I30};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  task29->add_dep(task30);
  energyq->add_task(task30);

  vector<shared_ptr<Tensor>> tensor31 = {I30, t2};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  task30->add_dep(task31);
  energyq->add_task(task31);

  vector<IndexRange> I33_index = {virt_, closed_, virt_, closed_};
  auto I33 = make_shared<Tensor>(I33_index);
  vector<shared_ptr<Tensor>> tensor32 = {I29, t2, I33};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  task29->add_dep(task32);
  energyq->add_task(task32);

  vector<shared_ptr<Tensor>> tensor33 = {I33, t2};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  task32->add_dep(task33);
  energyq->add_task(task33);

  vector<IndexRange> I35_index = {closed_, virt_, closed_, virt_};
  auto I35 = make_shared<Tensor>(I35_index);
  vector<shared_ptr<Tensor>> tensor34 = {I16, t2, I35};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  task19->add_dep(task34);
  energyq->add_task(task34);

  vector<shared_ptr<Tensor>> tensor35 = {I35, t2, v2_};
  auto task35 = make_shared<Task35>(tensor35, pindex, this->e0_);
  task34->add_dep(task35);
  energyq->add_task(task35);

  vector<IndexRange> I37_index = {virt_, closed_, virt_, closed_};
  auto I37 = make_shared<Tensor>(I37_index);
  vector<shared_ptr<Tensor>> tensor36 = {I16, t2, I37};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  task19->add_dep(task36);
  energyq->add_task(task36);

  vector<shared_ptr<Tensor>> tensor37 = {I37, t2};
  auto task37 = make_shared<Task37>(tensor37, pindex, this->e0_);
  task36->add_dep(task37);
  energyq->add_task(task37);

  vector<IndexRange> I39_index = {virt_, closed_, virt_, closed_};
  auto I39 = make_shared<Tensor>(I39_index);
  vector<shared_ptr<Tensor>> tensor38 = {I16, v2_, I39};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  task19->add_dep(task38);
  energyq->add_task(task38);

  vector<shared_ptr<Tensor>> tensor39 = {I39, t2};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  task38->add_dep(task39);
  energyq->add_task(task39);

  return energyq;
}


