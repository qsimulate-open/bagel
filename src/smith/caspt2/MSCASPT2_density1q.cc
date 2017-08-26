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
  auto tensor251 = vector<shared_ptr<Tensor>>{den1};
  auto task251 = make_shared<Task251>(tensor251, reset);
  density1q->add_task(task251);

  vector<IndexRange> I290_index = {closed_, active_};
  auto I290 = make_shared<Tensor>(I290_index);
  auto tensor252 = vector<shared_ptr<Tensor>>{den1, I290};
  auto task252 = make_shared<Task252>(tensor252, pindex);
  task252->add_dep(task251);
  density1q->add_task(task252);

  auto tensor253 = vector<shared_ptr<Tensor>>{I290, Gamma12_(), l2};
  auto task253 = make_shared<Task253>(tensor253, pindex);
  task252->add_dep(task253);
  task253->add_dep(task251);
  density1q->add_task(task253);

  vector<IndexRange> I292_index = {virt_, closed_};
  auto I292 = make_shared<Tensor>(I292_index);
  auto tensor254 = vector<shared_ptr<Tensor>>{den1, I292};
  auto task254 = make_shared<Task254>(tensor254, pindex);
  task254->add_dep(task251);
  density1q->add_task(task254);

  vector<IndexRange> I293_index = {active_, virt_, closed_, active_};
  auto I293 = make_shared<Tensor>(I293_index);
  auto tensor255 = vector<shared_ptr<Tensor>>{I292, Gamma65_(), I293};
  auto task255 = make_shared<Task255>(tensor255, pindex);
  task254->add_dep(task255);
  task255->add_dep(task251);
  density1q->add_task(task255);

  auto tensor256 = vector<shared_ptr<Tensor>>{I293, l2};
  auto task256 = make_shared<Task256>(tensor256, pindex);
  task255->add_dep(task256);
  task256->add_dep(task251);
  density1q->add_task(task256);

  vector<IndexRange> I296_index = {virt_, active_};
  auto I296 = make_shared<Tensor>(I296_index);
  auto tensor257 = vector<shared_ptr<Tensor>>{den1, I296};
  auto task257 = make_shared<Task257>(tensor257, pindex);
  task257->add_dep(task251);
  density1q->add_task(task257);

  auto tensor258 = vector<shared_ptr<Tensor>>{I296, Gamma60_(), l2};
  auto task258 = make_shared<Task258>(tensor258, pindex);
  task257->add_dep(task258);
  task258->add_dep(task251);
  density1q->add_task(task258);

  return density1q;
}


#endif
