//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_density1q.cc
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
#include <src/smith/caspt2/MSCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_density1q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density1q = make_shared<Queue>();
  auto tensor283 = vector<shared_ptr<Tensor>>{den1};
  auto task283 = make_shared<Task283>(tensor283, reset);
  density1q->add_task(task283);

  vector<IndexRange> I290_index = {closed_, active_};
  auto I290 = make_shared<Tensor>(I290_index);
  auto tensor284 = vector<shared_ptr<Tensor>>{den1, I290};
  auto task284 = make_shared<Task284>(tensor284, pindex);
  task284->add_dep(task283);
  density1q->add_task(task284);

  auto tensor285 = vector<shared_ptr<Tensor>>{I290, Gamma12_(), l2};
  auto task285 = make_shared<Task285>(tensor285, pindex);
  task284->add_dep(task285);
  task285->add_dep(task283);
  density1q->add_task(task285);

  vector<IndexRange> I292_index = {virt_, closed_};
  auto I292 = make_shared<Tensor>(I292_index);
  auto tensor286 = vector<shared_ptr<Tensor>>{den1, I292};
  auto task286 = make_shared<Task286>(tensor286, pindex);
  task286->add_dep(task283);
  density1q->add_task(task286);

  vector<IndexRange> I293_index = {active_, virt_, closed_, active_};
  auto I293 = make_shared<Tensor>(I293_index);
  auto tensor287 = vector<shared_ptr<Tensor>>{I292, Gamma65_(), I293};
  auto task287 = make_shared<Task287>(tensor287, pindex);
  task286->add_dep(task287);
  task287->add_dep(task283);
  density1q->add_task(task287);

  auto tensor288 = vector<shared_ptr<Tensor>>{I293, l2};
  auto task288 = make_shared<Task288>(tensor288, pindex);
  task287->add_dep(task288);
  task288->add_dep(task283);
  density1q->add_task(task288);

  vector<IndexRange> I296_index = {virt_, active_};
  auto I296 = make_shared<Tensor>(I296_index);
  auto tensor289 = vector<shared_ptr<Tensor>>{den1, I296};
  auto task289 = make_shared<Task289>(tensor289, pindex);
  task289->add_dep(task283);
  density1q->add_task(task289);

  auto tensor290 = vector<shared_ptr<Tensor>>{I296, Gamma60_(), l2};
  auto task290 = make_shared<Task290>(tensor290, pindex);
  task289->add_dep(task290);
  task290->add_dep(task283);
  density1q->add_task(task290);

  return density1q;
}


#endif
