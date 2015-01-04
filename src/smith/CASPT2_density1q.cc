//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_density1qq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_density1q() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density1q = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor1053 = {den1};
  auto task1053 = make_shared<Task1053>(tensor1053);
  density1q->add_task(task1053);

  vector<IndexRange> I1162_index = {closed_, active_};
  auto I1162 = make_shared<Tensor>(I1162_index);
  vector<shared_ptr<Tensor>> tensor1054 = {den1, I1162};
  auto task1054 = make_shared<Task1054>(tensor1054, pindex);
  task1054->add_dep(task1053);
  density1q->add_task(task1054);

  vector<IndexRange> I1163_index = {active_, active_, closed_, active_};
  auto I1163 = make_shared<Tensor>(I1163_index);
  vector<shared_ptr<Tensor>> tensor1055 = {I1162, Gamma12_(), I1163};
  auto task1055 = make_shared<Task1055>(tensor1055, pindex);
  task1054->add_dep(task1055);
  task1055->add_dep(task1053);
  density1q->add_task(task1055);

  vector<shared_ptr<Tensor>> tensor1056 = {I1163, t2};
  auto task1056 = make_shared<Task1056>(tensor1056, pindex);
  task1055->add_dep(task1056);
  task1056->add_dep(task1053);
  density1q->add_task(task1056);

  vector<IndexRange> I1164_index = {virt_, closed_};
  auto I1164 = make_shared<Tensor>(I1164_index);
  vector<shared_ptr<Tensor>> tensor1057 = {den1, I1164};
  auto task1057 = make_shared<Task1057>(tensor1057, pindex);
  task1057->add_dep(task1053);
  density1q->add_task(task1057);

  vector<IndexRange> I1165_index = {active_, virt_, closed_, active_};
  auto I1165 = make_shared<Tensor>(I1165_index);
  vector<shared_ptr<Tensor>> tensor1058 = {I1164, Gamma38_(), I1165};
  auto task1058 = make_shared<Task1058>(tensor1058, pindex);
  task1057->add_dep(task1058);
  task1058->add_dep(task1053);
  density1q->add_task(task1058);

  vector<shared_ptr<Tensor>> tensor1059 = {I1165, t2};
  auto task1059 = make_shared<Task1059>(tensor1059, pindex);
  task1058->add_dep(task1059);
  task1059->add_dep(task1053);
  density1q->add_task(task1059);

  vector<IndexRange> I1167_index = {active_, active_};
  auto I1167 = make_shared<Tensor>(I1167_index);
  vector<shared_ptr<Tensor>> tensor1060 = {I1164, t2, I1167};
  auto task1060 = make_shared<Task1060>(tensor1060, pindex);
  task1057->add_dep(task1060);
  task1060->add_dep(task1053);
  density1q->add_task(task1060);

  vector<shared_ptr<Tensor>> tensor1061 = {I1167, Gamma38_()};
  auto task1061 = make_shared<Task1061>(tensor1061, pindex);
  task1060->add_dep(task1061);
  task1061->add_dep(task1053);
  density1q->add_task(task1061);

  vector<IndexRange> I1168_index = {active_, virt_};
  auto I1168 = make_shared<Tensor>(I1168_index);
  vector<shared_ptr<Tensor>> tensor1062 = {den1, I1168};
  auto task1062 = make_shared<Task1062>(tensor1062, pindex);
  task1062->add_dep(task1053);
  density1q->add_task(task1062);

  vector<IndexRange> I1169_index = {active_, active_, active_, active_};
  auto I1169 = make_shared<Tensor>(I1169_index);
  vector<shared_ptr<Tensor>> tensor1063 = {I1168, t2, I1169};
  auto task1063 = make_shared<Task1063>(tensor1063, pindex);
  task1062->add_dep(task1063);
  task1063->add_dep(task1053);
  density1q->add_task(task1063);

  vector<shared_ptr<Tensor>> tensor1064 = {I1169, Gamma60_()};
  auto task1064 = make_shared<Task1064>(tensor1064, pindex);
  task1063->add_dep(task1064);
  task1064->add_dep(task1053);
  density1q->add_task(task1064);

  return density1q;
}


