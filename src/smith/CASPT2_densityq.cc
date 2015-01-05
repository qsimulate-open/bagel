//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_densityqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_densityq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto densityq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor25 = {den2};
  auto task25 = make_shared<Task25>(tensor25);
  densityq->add_task(task25);

  vector<IndexRange> I24_index = {active_, active_};
  auto I24 = make_shared<Tensor>(I24_index);
  vector<shared_ptr<Tensor>> tensor26 = {den2, I24};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  task26->add_dep(task25);
  densityq->add_task(task26);

  vector<IndexRange> I25_index;
  auto I25 = make_shared<Tensor>(I25_index);
  vector<shared_ptr<Tensor>> tensor27 = {I24, Gamma2_(), I25};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  task26->add_dep(task27);
  task27->add_dep(task25);
  densityq->add_task(task27);

  vector<IndexRange> I26_index = {virt_, closed_, virt_, closed_};
  auto I26 = make_shared<Tensor>(I26_index);
  vector<shared_ptr<Tensor>> tensor28 = {I25, t2, I26};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  task27->add_dep(task28);
  task28->add_dep(task25);
  densityq->add_task(task28);

  vector<shared_ptr<Tensor>> tensor29 = {I26, t2};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  task28->add_dep(task29);
  task29->add_dep(task25);
  densityq->add_task(task29);

  vector<IndexRange> I29_index = {closed_, virt_, closed_, virt_};
  auto I29 = make_shared<Tensor>(I29_index);
  vector<shared_ptr<Tensor>> tensor30 = {I25, t2, I29};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  task27->add_dep(task30);
  task30->add_dep(task25);
  densityq->add_task(task30);

  vector<shared_ptr<Tensor>> tensor31 = {I29, t2};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  task30->add_dep(task31);
  task31->add_dep(task25);
  densityq->add_task(task31);

  vector<IndexRange> I30_index = {closed_, closed_};
  auto I30 = make_shared<Tensor>(I30_index);
  vector<shared_ptr<Tensor>> tensor32 = {den2, I30};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  task32->add_dep(task25);
  densityq->add_task(task32);

  vector<IndexRange> I31_index = {virt_, closed_, virt_, closed_};
  auto I31 = make_shared<Tensor>(I31_index);
  vector<shared_ptr<Tensor>> tensor33 = {I30, t2, I31};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  task32->add_dep(task33);
  task33->add_dep(task25);
  densityq->add_task(task33);

  vector<shared_ptr<Tensor>> tensor34 = {I31, t2};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  task33->add_dep(task34);
  task34->add_dep(task25);
  densityq->add_task(task34);

  vector<IndexRange> I33_index = {virt_, closed_, virt_, closed_};
  auto I33 = make_shared<Tensor>(I33_index);
  vector<shared_ptr<Tensor>> tensor35 = {I30, t2, I33};
  auto task35 = make_shared<Task35>(tensor35, pindex);
  task32->add_dep(task35);
  task35->add_dep(task25);
  densityq->add_task(task35);

  vector<shared_ptr<Tensor>> tensor36 = {I33, t2};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  task35->add_dep(task36);
  task36->add_dep(task25);
  densityq->add_task(task36);

  vector<IndexRange> I34_index = {virt_, virt_};
  auto I34 = make_shared<Tensor>(I34_index);
  vector<shared_ptr<Tensor>> tensor37 = {den2, I34};
  auto task37 = make_shared<Task37>(tensor37, pindex);
  task37->add_dep(task25);
  densityq->add_task(task37);

  vector<IndexRange> I35_index = {virt_, closed_, virt_, closed_};
  auto I35 = make_shared<Tensor>(I35_index);
  vector<shared_ptr<Tensor>> tensor38 = {I34, t2, I35};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  task37->add_dep(task38);
  task38->add_dep(task25);
  densityq->add_task(task38);

  vector<shared_ptr<Tensor>> tensor39 = {I35, t2};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  task38->add_dep(task39);
  task39->add_dep(task25);
  densityq->add_task(task39);

  vector<IndexRange> I37_index = {closed_, virt_, closed_, virt_};
  auto I37 = make_shared<Tensor>(I37_index);
  vector<shared_ptr<Tensor>> tensor40 = {I34, t2, I37};
  auto task40 = make_shared<Task40>(tensor40, pindex);
  task37->add_dep(task40);
  task40->add_dep(task25);
  densityq->add_task(task40);

  vector<shared_ptr<Tensor>> tensor41 = {I37, t2};
  auto task41 = make_shared<Task41>(tensor41, pindex);
  task40->add_dep(task41);
  task41->add_dep(task25);
  densityq->add_task(task41);

  return densityq;
}


