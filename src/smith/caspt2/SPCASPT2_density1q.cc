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
  auto tensor260 = vector<shared_ptr<Tensor>>{den1};
  auto task260 = make_shared<Task260>(tensor260, reset);
  density1q->add_task(task260);

  vector<IndexRange> I290_index = {closed_, active_};
  auto I290 = make_shared<Tensor>(I290_index);
  auto tensor261 = vector<shared_ptr<Tensor>>{den1, I290};
  auto task261 = make_shared<Task261>(tensor261, pindex);
  task261->add_dep(task260);
  density1q->add_task(task261);

  auto tensor262 = vector<shared_ptr<Tensor>>{I290, Gamma13_(), t2};
  auto task262 = make_shared<Task262>(tensor262, pindex);
  task261->add_dep(task262);
  task262->add_dep(task260);
  density1q->add_task(task262);

  vector<IndexRange> I292_index = {virt_, closed_};
  auto I292 = make_shared<Tensor>(I292_index);
  auto tensor263 = vector<shared_ptr<Tensor>>{den1, I292};
  auto task263 = make_shared<Task263>(tensor263, pindex);
  task263->add_dep(task260);
  density1q->add_task(task263);

  vector<IndexRange> I293_index = {active_, virt_, closed_, active_};
  auto I293 = make_shared<Tensor>(I293_index);
  auto tensor264 = vector<shared_ptr<Tensor>>{I292, Gamma65_(), I293};
  auto task264 = make_shared<Task264>(tensor264, pindex);
  task263->add_dep(task264);
  task264->add_dep(task260);
  density1q->add_task(task264);

  auto tensor265 = vector<shared_ptr<Tensor>>{I293, t2};
  auto task265 = make_shared<Task265>(tensor265, pindex);
  task264->add_dep(task265);
  task265->add_dep(task260);
  density1q->add_task(task265);

  vector<IndexRange> I296_index = {virt_, active_};
  auto I296 = make_shared<Tensor>(I296_index);
  auto tensor266 = vector<shared_ptr<Tensor>>{den1, I296};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task266->add_dep(task260);
  density1q->add_task(task266);

  auto tensor267 = vector<shared_ptr<Tensor>>{I296, Gamma60_(), t2};
  auto task267 = make_shared<Task267>(tensor267, pindex);
  task266->add_dep(task267);
  task267->add_dep(task260);
  density1q->add_task(task267);

  return density1q;
}


#endif
