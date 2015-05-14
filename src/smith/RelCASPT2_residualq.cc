//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2_residualqq.cc
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


#include <src/smith/RelCASPT2.h>
#include <src/smith/RelCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelCASPT2::RelCASPT2::make_residualq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto residualq = make_shared<Queue>();
  auto tensor2 = vector<shared_ptr<Tensor>>{r};
  auto task2 = make_shared<Task2>(tensor2, reset);
  residualq->add_task(task2);

  vector<IndexRange> I0_index = {virt_, closed_, virt_, active_};
  auto I0 = make_shared<Tensor>(I0_index);
  auto tensor3 = vector<shared_ptr<Tensor>>{r, I0};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  task3->add_dep(task2);
  residualq->add_task(task3);

  vector<IndexRange> I1_index = {active_, virt_, closed_, virt_};
  auto I1 = make_shared<Tensor>(I1_index);
  auto tensor4 = vector<shared_ptr<Tensor>>{I0, Gamma0_(), I1};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  task3->add_dep(task4);
  task4->add_dep(task2);
  residualq->add_task(task4);

  auto tensor5 = vector<shared_ptr<Tensor>>{I1, t2};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  task4->add_dep(task5);
  task5->add_dep(task2);
  residualq->add_task(task5);

  vector<IndexRange> I5_index = {closed_, active_, virt_, virt_};
  auto I5 = make_shared<Tensor>(I5_index);
  auto tensor6 = vector<shared_ptr<Tensor>>{I0, Gamma2_(), I5};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  task3->add_dep(task6);
  task6->add_dep(task2);
  residualq->add_task(task6);

  auto tensor7 = vector<shared_ptr<Tensor>>{I5, t2, v2_};
  auto task7 = make_shared<Task7>(tensor7, pindex, this->e0_);
  task6->add_dep(task7);
  task7->add_dep(task2);
  residualq->add_task(task7);

  auto tensor8 = vector<shared_ptr<Tensor>>{I5, t2, f1_};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  task6->add_dep(task8);
  task8->add_dep(task2);
  residualq->add_task(task8);

  auto tensor9 = vector<shared_ptr<Tensor>>{I5, t2, f1_};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  task6->add_dep(task9);
  task9->add_dep(task2);
  residualq->add_task(task9);

  auto tensor10 = vector<shared_ptr<Tensor>>{I5, t2, f1_};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  task6->add_dep(task10);
  task10->add_dep(task2);
  residualq->add_task(task10);

  auto tensor11 = vector<shared_ptr<Tensor>>{I5, t2, f1_};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  task6->add_dep(task11);
  task11->add_dep(task2);
  residualq->add_task(task11);

  auto tensor12 = vector<shared_ptr<Tensor>>{I5, t2, f1_};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  task6->add_dep(task12);
  task12->add_dep(task2);
  residualq->add_task(task12);

  auto tensor13 = vector<shared_ptr<Tensor>>{I5, t2, f1_};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  task6->add_dep(task13);
  task13->add_dep(task2);
  residualq->add_task(task13);

  return residualq;
}


#endif
