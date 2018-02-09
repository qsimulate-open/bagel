//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_deci2q.cc
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

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_deci2q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto deci2q = make_shared<Queue>();
  auto tensor470 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci};
  auto task470 = make_shared<Task470>(tensor470, reset);
  deci2q->add_task(task470);

  vector<IndexRange> I703_index = {active_, active_, active_, active_};
  auto I703 = make_shared<Tensor>(I703_index);
  auto tensor472 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I703};
  auto task472 = make_shared<Task472>(tensor472, pindex);
  task472->add_dep(task470);
  deci2q->add_task(task472);

  auto tensor473 = vector<shared_ptr<Tensor>>{I703, t2, l2};
  auto task473 = make_shared<Task473>(tensor473, pindex, this->e0_);
  task472->add_dep(task473);
  task473->add_dep(task470);
  deci2q->add_task(task473);

  vector<IndexRange> I706_index = {active_, active_, active_, active_, active_, active_};
  auto I706 = make_shared<Tensor>(I706_index);
  auto tensor474 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I706};
  auto task474 = make_shared<Task474>(tensor474, pindex);
  task474->add_dep(task470);
  deci2q->add_task(task474);

  auto tensor475 = vector<shared_ptr<Tensor>>{I706, t2, l2};
  auto task475 = make_shared<Task475>(tensor475, pindex, this->e0_);
  task474->add_dep(task475);
  task475->add_dep(task470);
  deci2q->add_task(task475);

  vector<IndexRange> I709_index = {active_, active_};
  auto I709 = make_shared<Tensor>(I709_index);
  auto tensor476 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I709};
  auto task476 = make_shared<Task476>(tensor476, pindex);
  task476->add_dep(task470);
  deci2q->add_task(task476);

  auto tensor477 = vector<shared_ptr<Tensor>>{I709, t2, l2};
  auto task477 = make_shared<Task477>(tensor477, pindex, this->e0_);
  task476->add_dep(task477);
  task477->add_dep(task470);
  deci2q->add_task(task477);

  auto tensor478 = vector<shared_ptr<Tensor>>{I709, t2, l2};
  auto task478 = make_shared<Task478>(tensor478, pindex, this->e0_);
  task476->add_dep(task478);
  task478->add_dep(task470);
  deci2q->add_task(task478);

  vector<IndexRange> I715_index = {active_, active_, active_, active_};
  auto I715 = make_shared<Tensor>(I715_index);
  auto tensor479 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I715};
  auto task479 = make_shared<Task479>(tensor479, pindex);
  task479->add_dep(task470);
  deci2q->add_task(task479);

  auto tensor480 = vector<shared_ptr<Tensor>>{I715, t2, l2};
  auto task480 = make_shared<Task480>(tensor480, pindex, this->e0_);
  task479->add_dep(task480);
  task480->add_dep(task470);
  deci2q->add_task(task480);

  vector<IndexRange> I718_index = {active_, active_, active_, active_};
  auto I718 = make_shared<Tensor>(I718_index);
  auto tensor481 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I718};
  auto task481 = make_shared<Task481>(tensor481, pindex);
  task481->add_dep(task470);
  deci2q->add_task(task481);

  auto tensor482 = vector<shared_ptr<Tensor>>{I718, t2, l2};
  auto task482 = make_shared<Task482>(tensor482, pindex, this->e0_);
  task481->add_dep(task482);
  task482->add_dep(task470);
  deci2q->add_task(task482);

  auto tensor483 = vector<shared_ptr<Tensor>>{I718, t2, l2};
  auto task483 = make_shared<Task483>(tensor483, pindex, this->e0_);
  task481->add_dep(task483);
  task483->add_dep(task470);
  deci2q->add_task(task483);

  auto tensor484 = vector<shared_ptr<Tensor>>{I718, t2, l2};
  auto task484 = make_shared<Task484>(tensor484, pindex, this->e0_);
  task481->add_dep(task484);
  task484->add_dep(task470);
  deci2q->add_task(task484);

  vector<IndexRange> I727_index = {active_, active_, active_, active_, active_, active_};
  auto I727 = make_shared<Tensor>(I727_index);
  auto tensor485 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I727};
  auto task485 = make_shared<Task485>(tensor485, pindex);
  task485->add_dep(task470);
  deci2q->add_task(task485);

  auto tensor486 = vector<shared_ptr<Tensor>>{I727, t2, l2};
  auto task486 = make_shared<Task486>(tensor486, pindex, this->e0_);
  task485->add_dep(task486);
  task486->add_dep(task470);
  deci2q->add_task(task486);

  shared_ptr<Tensor> I730;
  if (diagonal) {
    vector<IndexRange> I730_index;
    I730 = make_shared<Tensor>(I730_index);
  }
  shared_ptr<Task487> task487;
  if (diagonal) {
    auto tensor487 = vector<shared_ptr<Tensor>>{den0ci, I730};
    task487 = make_shared<Task487>(tensor487, pindex);
    task487->add_dep(task470);
    deci2q->add_task(task487);
  }

  shared_ptr<Task488> task488;
  if (diagonal) {
    auto tensor488 = vector<shared_ptr<Tensor>>{I730, t2, l2};
    task488 = make_shared<Task488>(tensor488, pindex, this->e0_);
    task487->add_dep(task488);
    task488->add_dep(task470);
    deci2q->add_task(task488);
  }

  shared_ptr<Tensor> I733;
  if (diagonal) {
    vector<IndexRange> I733_index;
    I733 = make_shared<Tensor>(I733_index);
  }
  shared_ptr<Task489> task489;
  if (diagonal) {
    auto tensor489 = vector<shared_ptr<Tensor>>{den0ci, I733};
    task489 = make_shared<Task489>(tensor489, pindex);
    task489->add_dep(task470);
    deci2q->add_task(task489);
  }

  shared_ptr<Task490> task490;
  if (diagonal) {
    auto tensor490 = vector<shared_ptr<Tensor>>{I733, t2, l2};
    task490 = make_shared<Task490>(tensor490, pindex, this->e0_);
    task489->add_dep(task490);
    task490->add_dep(task470);
    deci2q->add_task(task490);
  }

  vector<IndexRange> I736_index = {active_, active_};
  auto I736 = make_shared<Tensor>(I736_index);
  auto tensor491 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I736};
  auto task491 = make_shared<Task491>(tensor491, pindex);
  task491->add_dep(task470);
  deci2q->add_task(task491);

  auto tensor492 = vector<shared_ptr<Tensor>>{I736, t2, l2};
  auto task492 = make_shared<Task492>(tensor492, pindex, this->e0_);
  task491->add_dep(task492);
  task492->add_dep(task470);
  deci2q->add_task(task492);

  auto tensor493 = vector<shared_ptr<Tensor>>{I736, t2, l2};
  auto task493 = make_shared<Task493>(tensor493, pindex, this->e0_);
  task491->add_dep(task493);
  task493->add_dep(task470);
  deci2q->add_task(task493);

  vector<IndexRange> I742_index = {active_, active_, active_, active_};
  auto I742 = make_shared<Tensor>(I742_index);
  auto tensor494 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I742};
  auto task494 = make_shared<Task494>(tensor494, pindex);
  task494->add_dep(task470);
  deci2q->add_task(task494);

  auto tensor495 = vector<shared_ptr<Tensor>>{I742, t2, l2};
  auto task495 = make_shared<Task495>(tensor495, pindex, this->e0_);
  task494->add_dep(task495);
  task495->add_dep(task470);
  deci2q->add_task(task495);

  return deci2q;
}


#endif
