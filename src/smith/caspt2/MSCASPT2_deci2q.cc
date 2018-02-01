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
  auto tensor458 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci};
  auto task458 = make_shared<Task458>(tensor458, reset);
  deci2q->add_task(task458);

  vector<IndexRange> I687_index = {active_, active_, active_, active_};
  auto I687 = make_shared<Tensor>(I687_index);
  auto tensor460 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I687};
  auto task460 = make_shared<Task460>(tensor460, pindex);
  task460->add_dep(task458);
  deci2q->add_task(task460);

  auto tensor461 = vector<shared_ptr<Tensor>>{I687, t2, l2};
  auto task461 = make_shared<Task461>(tensor461, pindex, this->e0_);
  task460->add_dep(task461);
  task461->add_dep(task458);
  deci2q->add_task(task461);

  vector<IndexRange> I690_index = {active_, active_, active_, active_, active_, active_};
  auto I690 = make_shared<Tensor>(I690_index);
  auto tensor462 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I690};
  auto task462 = make_shared<Task462>(tensor462, pindex);
  task462->add_dep(task458);
  deci2q->add_task(task462);

  auto tensor463 = vector<shared_ptr<Tensor>>{I690, t2, l2};
  auto task463 = make_shared<Task463>(tensor463, pindex, this->e0_);
  task462->add_dep(task463);
  task463->add_dep(task458);
  deci2q->add_task(task463);

  vector<IndexRange> I693_index = {active_, active_};
  auto I693 = make_shared<Tensor>(I693_index);
  auto tensor464 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I693};
  auto task464 = make_shared<Task464>(tensor464, pindex);
  task464->add_dep(task458);
  deci2q->add_task(task464);

  auto tensor465 = vector<shared_ptr<Tensor>>{I693, t2, l2};
  auto task465 = make_shared<Task465>(tensor465, pindex, this->e0_);
  task464->add_dep(task465);
  task465->add_dep(task458);
  deci2q->add_task(task465);

  auto tensor466 = vector<shared_ptr<Tensor>>{I693, t2, l2};
  auto task466 = make_shared<Task466>(tensor466, pindex, this->e0_);
  task464->add_dep(task466);
  task466->add_dep(task458);
  deci2q->add_task(task466);

  vector<IndexRange> I699_index = {active_, active_, active_, active_};
  auto I699 = make_shared<Tensor>(I699_index);
  auto tensor467 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I699};
  auto task467 = make_shared<Task467>(tensor467, pindex);
  task467->add_dep(task458);
  deci2q->add_task(task467);

  auto tensor468 = vector<shared_ptr<Tensor>>{I699, t2, l2};
  auto task468 = make_shared<Task468>(tensor468, pindex, this->e0_);
  task467->add_dep(task468);
  task468->add_dep(task458);
  deci2q->add_task(task468);

  vector<IndexRange> I702_index = {active_, active_, active_, active_};
  auto I702 = make_shared<Tensor>(I702_index);
  auto tensor469 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I702};
  auto task469 = make_shared<Task469>(tensor469, pindex);
  task469->add_dep(task458);
  deci2q->add_task(task469);

  auto tensor470 = vector<shared_ptr<Tensor>>{I702, t2, l2};
  auto task470 = make_shared<Task470>(tensor470, pindex, this->e0_);
  task469->add_dep(task470);
  task470->add_dep(task458);
  deci2q->add_task(task470);

  auto tensor471 = vector<shared_ptr<Tensor>>{I702, t2, l2};
  auto task471 = make_shared<Task471>(tensor471, pindex, this->e0_);
  task469->add_dep(task471);
  task471->add_dep(task458);
  deci2q->add_task(task471);

  auto tensor472 = vector<shared_ptr<Tensor>>{I702, t2, l2};
  auto task472 = make_shared<Task472>(tensor472, pindex, this->e0_);
  task469->add_dep(task472);
  task472->add_dep(task458);
  deci2q->add_task(task472);

  vector<IndexRange> I711_index = {active_, active_, active_, active_, active_, active_};
  auto I711 = make_shared<Tensor>(I711_index);
  auto tensor473 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I711};
  auto task473 = make_shared<Task473>(tensor473, pindex);
  task473->add_dep(task458);
  deci2q->add_task(task473);

  auto tensor474 = vector<shared_ptr<Tensor>>{I711, t2, l2};
  auto task474 = make_shared<Task474>(tensor474, pindex, this->e0_);
  task473->add_dep(task474);
  task474->add_dep(task458);
  deci2q->add_task(task474);

  vector<IndexRange> I714_index = {active_, active_};
  auto I714 = make_shared<Tensor>(I714_index);
  auto tensor475 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I714};
  auto task475 = make_shared<Task475>(tensor475, pindex);
  task475->add_dep(task458);
  deci2q->add_task(task475);

  auto tensor476 = vector<shared_ptr<Tensor>>{I714, t2, l2};
  auto task476 = make_shared<Task476>(tensor476, pindex, this->e0_);
  task475->add_dep(task476);
  task476->add_dep(task458);
  deci2q->add_task(task476);

  auto tensor477 = vector<shared_ptr<Tensor>>{I714, t2, l2};
  auto task477 = make_shared<Task477>(tensor477, pindex, this->e0_);
  task475->add_dep(task477);
  task477->add_dep(task458);
  deci2q->add_task(task477);

  vector<IndexRange> I720_index = {active_, active_, active_, active_};
  auto I720 = make_shared<Tensor>(I720_index);
  auto tensor478 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I720};
  auto task478 = make_shared<Task478>(tensor478, pindex);
  task478->add_dep(task458);
  deci2q->add_task(task478);

  auto tensor479 = vector<shared_ptr<Tensor>>{I720, t2, l2};
  auto task479 = make_shared<Task479>(tensor479, pindex, this->e0_);
  task478->add_dep(task479);
  task479->add_dep(task458);
  deci2q->add_task(task479);

  return deci2q;
}


#endif
