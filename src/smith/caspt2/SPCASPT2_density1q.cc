//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: SPCASPT2_density1qq.cc
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


#include <src/smith/caspt2/SPCASPT2.h>
#include <src/smith/caspt2/SPCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> SPCASPT2::SPCASPT2::make_density1q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density1q = make_shared<Queue>();
  auto tensor265 = vector<shared_ptr<Tensor>>{den1};
  auto task265 = make_shared<Task265>(tensor265, reset);
  density1q->add_task(task265);

  vector<IndexRange> I290_index = {closed_, active_};
  auto I290 = make_shared<Tensor>(I290_index);
  auto tensor266 = vector<shared_ptr<Tensor>>{den1, I290};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task266->add_dep(task265);
  density1q->add_task(task266);

  auto tensor267 = vector<shared_ptr<Tensor>>{I290, Gamma13_(), t2};
  auto task267 = make_shared<Task267>(tensor267, pindex);
  task266->add_dep(task267);
  task267->add_dep(task265);
  density1q->add_task(task267);

  vector<IndexRange> I292_index = {virt_, closed_};
  auto I292 = make_shared<Tensor>(I292_index);
  auto tensor268 = vector<shared_ptr<Tensor>>{den1, I292};
  auto task268 = make_shared<Task268>(tensor268, pindex);
  task268->add_dep(task265);
  density1q->add_task(task268);

  auto tensor269 = vector<shared_ptr<Tensor>>{I292, Gamma65_(), t2};
  auto task269 = make_shared<Task269>(tensor269, pindex);
  task268->add_dep(task269);
  task269->add_dep(task265);
  density1q->add_task(task269);

  auto tensor270 = vector<shared_ptr<Tensor>>{I292, Gamma67_(), t2};
  auto task270 = make_shared<Task270>(tensor270, pindex);
  task268->add_dep(task270);
  task270->add_dep(task265);
  density1q->add_task(task270);

  vector<IndexRange> I296_index = {virt_, active_};
  auto I296 = make_shared<Tensor>(I296_index);
  auto tensor271 = vector<shared_ptr<Tensor>>{den1, I296};
  auto task271 = make_shared<Task271>(tensor271, pindex);
  task271->add_dep(task265);
  density1q->add_task(task271);

  auto tensor272 = vector<shared_ptr<Tensor>>{I296, Gamma60_(), t2};
  auto task272 = make_shared<Task272>(tensor272, pindex);
  task271->add_dep(task272);
  task272->add_dep(task265);
  density1q->add_task(task272);

  return density1q;
}


#endif
