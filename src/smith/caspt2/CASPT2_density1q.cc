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
  auto tensor488 = vector<shared_ptr<Tensor>>{den1};
  auto task488 = make_shared<Task488>(tensor488, reset);
  density1q->add_task(task488);

  vector<IndexRange> I650_index = {closed_, active_};
  auto I650 = make_shared<Tensor>(I650_index);
  auto tensor489 = vector<shared_ptr<Tensor>>{den1, I650};
  auto task489 = make_shared<Task489>(tensor489, pindex);
  task489->add_dep(task488);
  density1q->add_task(task489);

  auto tensor490 = vector<shared_ptr<Tensor>>{I650, Gamma12_(), t2};
  auto task490 = make_shared<Task490>(tensor490, pindex);
  task489->add_dep(task490);
  task490->add_dep(task488);
  density1q->add_task(task490);

  vector<IndexRange> I652_index = {virt_, closed_};
  auto I652 = make_shared<Tensor>(I652_index);
  auto tensor491 = vector<shared_ptr<Tensor>>{den1, I652};
  auto task491 = make_shared<Task491>(tensor491, pindex);
  task491->add_dep(task488);
  density1q->add_task(task491);

  vector<IndexRange> I653_index = {active_, virt_, closed_, active_};
  auto I653 = make_shared<Tensor>(I653_index);
  auto tensor492 = vector<shared_ptr<Tensor>>{I652, Gamma38_(), I653};
  auto task492 = make_shared<Task492>(tensor492, pindex);
  task491->add_dep(task492);
  task492->add_dep(task488);
  density1q->add_task(task492);

  auto tensor493 = vector<shared_ptr<Tensor>>{I653, t2};
  auto task493 = make_shared<Task493>(tensor493, pindex);
  task492->add_dep(task493);
  task493->add_dep(task488);
  density1q->add_task(task493);

  vector<IndexRange> I656_index = {virt_, active_};
  auto I656 = make_shared<Tensor>(I656_index);
  auto tensor494 = vector<shared_ptr<Tensor>>{den1, I656};
  auto task494 = make_shared<Task494>(tensor494, pindex);
  task494->add_dep(task488);
  density1q->add_task(task494);

  auto tensor495 = vector<shared_ptr<Tensor>>{I656, Gamma60_(), t2};
  auto task495 = make_shared<Task495>(tensor495, pindex);
  task494->add_dep(task495);
  task495->add_dep(task488);
  density1q->add_task(task495);

  return density1q;
}


#endif
