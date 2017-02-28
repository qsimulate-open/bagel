//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_density1qq.cc
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


#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_density1q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density1q = make_shared<Queue>();
  auto tensor520 = vector<shared_ptr<Tensor>>{den1};
  auto task520 = make_shared<Task520>(tensor520, reset);
  density1q->add_task(task520);

  vector<IndexRange> I650_index = {closed_, active_};
  auto I650 = make_shared<Tensor>(I650_index);
  auto tensor521 = vector<shared_ptr<Tensor>>{den1, I650};
  auto task521 = make_shared<Task521>(tensor521, pindex);
  task521->add_dep(task520);
  density1q->add_task(task521);

  auto tensor522 = vector<shared_ptr<Tensor>>{I650, Gamma12_(), t2};
  auto task522 = make_shared<Task522>(tensor522, pindex);
  task521->add_dep(task522);
  task522->add_dep(task520);
  density1q->add_task(task522);

  vector<IndexRange> I652_index = {virt_, closed_};
  auto I652 = make_shared<Tensor>(I652_index);
  auto tensor523 = vector<shared_ptr<Tensor>>{den1, I652};
  auto task523 = make_shared<Task523>(tensor523, pindex);
  task523->add_dep(task520);
  density1q->add_task(task523);

  vector<IndexRange> I653_index = {active_, virt_, closed_, active_};
  auto I653 = make_shared<Tensor>(I653_index);
  auto tensor524 = vector<shared_ptr<Tensor>>{I652, Gamma38_(), I653};
  auto task524 = make_shared<Task524>(tensor524, pindex);
  task523->add_dep(task524);
  task524->add_dep(task520);
  density1q->add_task(task524);

  auto tensor525 = vector<shared_ptr<Tensor>>{I653, t2};
  auto task525 = make_shared<Task525>(tensor525, pindex);
  task524->add_dep(task525);
  task525->add_dep(task520);
  density1q->add_task(task525);

  vector<IndexRange> I656_index = {virt_, active_};
  auto I656 = make_shared<Tensor>(I656_index);
  auto tensor526 = vector<shared_ptr<Tensor>>{den1, I656};
  auto task526 = make_shared<Task526>(tensor526, pindex);
  task526->add_dep(task520);
  density1q->add_task(task526);

  auto tensor527 = vector<shared_ptr<Tensor>>{I656, Gamma60_(), t2};
  auto task527 = make_shared<Task527>(tensor527, pindex);
  task526->add_dep(task527);
  task527->add_dep(task520);
  density1q->add_task(task527);

  return density1q;
}


#endif
