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
  vector<shared_ptr<Tensor>> tensor3 = {r};
  auto task3 = make_shared<Task3>(tensor3);
  residualq->add_task(task3);

  vector<IndexRange> I0_index = {closed_, virt_, closed_, virt_};
  auto I0 = make_shared<Tensor>(I0_index);
  vector<shared_ptr<Tensor>> tensor4 = {r, I0};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  task4->add_dep(task3);
  residualq->add_task(task4);

  vector<shared_ptr<Tensor>> tensor5 = {I0, t2, v2_};
  auto task5 = make_shared<Task5>(tensor5, pindex, this->e0_);
  task4->add_dep(task5);
  task5->add_dep(task3);
  residualq->add_task(task5);

  vector<IndexRange> I1_index = {closed_, virt_, closed_, virt_};
  auto I1 = make_shared<Tensor>(I1_index);
  vector<shared_ptr<Tensor>> tensor6 = {I0, Gamma0_(), I1};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  task4->add_dep(task6);
  task6->add_dep(task3);
  residualq->add_task(task6);

  vector<shared_ptr<Tensor>> tensor7 = {I1, t2};
  auto task7 = make_shared<Task7>(tensor7, pindex);
  task6->add_dep(task7);
  task7->add_dep(task3);
  residualq->add_task(task7);

  vector<IndexRange> I3_index;
  auto I3 = make_shared<Tensor>(I3_index);
  vector<shared_ptr<Tensor>> tensor8 = {I0, t2, I3};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  task4->add_dep(task8);
  task8->add_dep(task3);
  residualq->add_task(task8);

  vector<shared_ptr<Tensor>> tensor9 = {I3, Gamma0_()};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  task8->add_dep(task9);
  task9->add_dep(task3);
  residualq->add_task(task9);

  vector<IndexRange> I4_index = {closed_, virt_, virt_, closed_};
  auto I4 = make_shared<Tensor>(I4_index);
  vector<shared_ptr<Tensor>> tensor10 = {r, I4};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  task10->add_dep(task3);
  residualq->add_task(task10);

  vector<IndexRange> I5_index = {closed_, virt_, closed_, virt_};
  auto I5 = make_shared<Tensor>(I5_index);
  vector<shared_ptr<Tensor>> tensor11 = {I4, f1_, I5};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  task10->add_dep(task11);
  task11->add_dep(task3);
  residualq->add_task(task11);

  vector<shared_ptr<Tensor>> tensor12 = {I5, t2};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  task11->add_dep(task12);
  task12->add_dep(task3);
  residualq->add_task(task12);

  vector<IndexRange> I7_index = {closed_, closed_};
  auto I7 = make_shared<Tensor>(I7_index);
  vector<shared_ptr<Tensor>> tensor13 = {I4, t2, I7};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  task10->add_dep(task13);
  task13->add_dep(task3);
  residualq->add_task(task13);

  vector<shared_ptr<Tensor>> tensor14 = {I7, f1_};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  task13->add_dep(task14);
  task14->add_dep(task3);
  residualq->add_task(task14);

  vector<IndexRange> I9_index = {closed_, virt_, closed_, virt_};
  auto I9 = make_shared<Tensor>(I9_index);
  vector<shared_ptr<Tensor>> tensor15 = {I4, f1_, I9};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  task10->add_dep(task15);
  task15->add_dep(task3);
  residualq->add_task(task15);

  vector<shared_ptr<Tensor>> tensor16 = {I9, t2};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  task15->add_dep(task16);
  task16->add_dep(task3);
  residualq->add_task(task16);

  vector<IndexRange> I11_index = {virt_, virt_};
  auto I11 = make_shared<Tensor>(I11_index);
  vector<shared_ptr<Tensor>> tensor17 = {I4, t2, I11};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  task10->add_dep(task17);
  task17->add_dep(task3);
  residualq->add_task(task17);

  vector<shared_ptr<Tensor>> tensor18 = {I11, f1_};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  task17->add_dep(task18);
  task18->add_dep(task3);
  residualq->add_task(task18);

  return residualq;
}


