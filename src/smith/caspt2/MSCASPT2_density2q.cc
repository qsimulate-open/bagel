//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_density2q.cc
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

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_density2q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density2q = make_shared<Queue>();
  auto tensor291 = vector<shared_ptr<Tensor>>{Den1};
  auto task291 = make_shared<Task291>(tensor291, reset);
  density2q->add_task(task291);

  vector<IndexRange> I298_index = {closed_, closed_, active_, active_};
  auto I298 = make_shared<Tensor>(I298_index);
  auto tensor292 = vector<shared_ptr<Tensor>>{Den1, I298};
  auto task292 = make_shared<Task292>(tensor292, pindex);
  task292->add_dep(task291);
  density2q->add_task(task292);

  auto tensor293 = vector<shared_ptr<Tensor>>{I298, Gamma1_(), l2};
  auto task293 = make_shared<Task293>(tensor293, pindex);
  task292->add_dep(task293);
  task293->add_dep(task291);
  density2q->add_task(task293);

  vector<IndexRange> I300_index = {closed_, active_, active_, active_};
  auto I300 = make_shared<Tensor>(I300_index);
  auto tensor294 = vector<shared_ptr<Tensor>>{Den1, I300};
  auto task294 = make_shared<Task294>(tensor294, pindex);
  task294->add_dep(task291);
  density2q->add_task(task294);

  auto tensor295 = vector<shared_ptr<Tensor>>{I300, Gamma6_(), l2};
  auto task295 = make_shared<Task295>(tensor295, pindex);
  task294->add_dep(task295);
  task295->add_dep(task291);
  density2q->add_task(task295);

  vector<IndexRange> I302_index = {closed_, virt_, closed_, active_};
  auto I302 = make_shared<Tensor>(I302_index);
  auto tensor296 = vector<shared_ptr<Tensor>>{Den1, I302};
  auto task296 = make_shared<Task296>(tensor296, pindex);
  task296->add_dep(task291);
  density2q->add_task(task296);

  vector<IndexRange> I303_index = {closed_, virt_, closed_, active_};
  auto I303 = make_shared<Tensor>(I303_index);
  auto tensor297 = vector<shared_ptr<Tensor>>{I302, Gamma16_(), I303};
  auto task297 = make_shared<Task297>(tensor297, pindex);
  task296->add_dep(task297);
  task297->add_dep(task291);
  density2q->add_task(task297);

  auto tensor298 = vector<shared_ptr<Tensor>>{I303, l2};
  auto task298 = make_shared<Task298>(tensor298, pindex);
  task297->add_dep(task298);
  task298->add_dep(task291);
  density2q->add_task(task298);

  vector<IndexRange> I306_index = {virt_, closed_, active_, active_};
  auto I306 = make_shared<Tensor>(I306_index);
  auto tensor299 = vector<shared_ptr<Tensor>>{Den1, I306};
  auto task299 = make_shared<Task299>(tensor299, pindex);
  task299->add_dep(task291);
  density2q->add_task(task299);

  auto tensor300 = vector<shared_ptr<Tensor>>{I306, Gamma32_(), l2};
  auto task300 = make_shared<Task300>(tensor300, pindex);
  task299->add_dep(task300);
  task300->add_dep(task291);
  density2q->add_task(task300);

  auto tensor301 = vector<shared_ptr<Tensor>>{I306, Gamma35_(), l2};
  auto task301 = make_shared<Task301>(tensor301, pindex);
  task299->add_dep(task301);
  task301->add_dep(task291);
  density2q->add_task(task301);

  vector<IndexRange> I310_index = {virt_, closed_, active_, active_};
  auto I310 = make_shared<Tensor>(I310_index);
  auto tensor302 = vector<shared_ptr<Tensor>>{Den1, I310};
  auto task302 = make_shared<Task302>(tensor302, pindex);
  task302->add_dep(task291);
  density2q->add_task(task302);

  vector<IndexRange> I311_index = {active_, virt_, closed_, active_};
  auto I311 = make_shared<Tensor>(I311_index);
  auto tensor303 = vector<shared_ptr<Tensor>>{I310, Gamma35_(), I311};
  auto task303 = make_shared<Task303>(tensor303, pindex);
  task302->add_dep(task303);
  task303->add_dep(task291);
  density2q->add_task(task303);

  auto tensor304 = vector<shared_ptr<Tensor>>{I311, l2};
  auto task304 = make_shared<Task304>(tensor304, pindex);
  task303->add_dep(task304);
  task304->add_dep(task291);
  density2q->add_task(task304);

  vector<IndexRange> I314_index = {virt_, active_, active_, active_};
  auto I314 = make_shared<Tensor>(I314_index);
  auto tensor305 = vector<shared_ptr<Tensor>>{Den1, I314};
  auto task305 = make_shared<Task305>(tensor305, pindex);
  task305->add_dep(task291);
  density2q->add_task(task305);

  auto tensor306 = vector<shared_ptr<Tensor>>{I314, Gamma62_(), l2};
  auto task306 = make_shared<Task306>(tensor306, pindex);
  task305->add_dep(task306);
  task306->add_dep(task291);
  density2q->add_task(task306);

  shared_ptr<Tensor> I316;
  if (diagonal) {
    vector<IndexRange> I316_index = {closed_, virt_, closed_, virt_};
    I316 = make_shared<Tensor>(I316_index);
  }
  shared_ptr<Task307> task307;
  if (diagonal) {
    auto tensor307 = vector<shared_ptr<Tensor>>{Den1, I316};
    task307 = make_shared<Task307>(tensor307, pindex);
    task307->add_dep(task291);
    density2q->add_task(task307);
  }

  shared_ptr<Task308> task308;
  if (diagonal) {
    auto tensor308 = vector<shared_ptr<Tensor>>{I316, l2};
    task308 = make_shared<Task308>(tensor308, pindex);
    task307->add_dep(task308);
    task308->add_dep(task291);
    density2q->add_task(task308);
  }

  vector<IndexRange> I318_index = {virt_, closed_, virt_, active_};
  auto I318 = make_shared<Tensor>(I318_index);
  auto tensor309 = vector<shared_ptr<Tensor>>{Den1, I318};
  auto task309 = make_shared<Task309>(tensor309, pindex);
  task309->add_dep(task291);
  density2q->add_task(task309);

  vector<IndexRange> I319_index = {active_, virt_, closed_, virt_};
  auto I319 = make_shared<Tensor>(I319_index);
  auto tensor310 = vector<shared_ptr<Tensor>>{I318, Gamma65_(), I319};
  auto task310 = make_shared<Task310>(tensor310, pindex);
  task309->add_dep(task310);
  task310->add_dep(task291);
  density2q->add_task(task310);

  auto tensor311 = vector<shared_ptr<Tensor>>{I319, l2};
  auto task311 = make_shared<Task311>(tensor311, pindex);
  task310->add_dep(task311);
  task311->add_dep(task291);
  density2q->add_task(task311);

  vector<IndexRange> I322_index = {virt_, virt_, active_, active_};
  auto I322 = make_shared<Tensor>(I322_index);
  auto tensor312 = vector<shared_ptr<Tensor>>{Den1, I322};
  auto task312 = make_shared<Task312>(tensor312, pindex);
  task312->add_dep(task291);
  density2q->add_task(task312);

  auto tensor313 = vector<shared_ptr<Tensor>>{I322, Gamma60_(), l2};
  auto task313 = make_shared<Task313>(tensor313, pindex);
  task312->add_dep(task313);
  task313->add_dep(task291);
  density2q->add_task(task313);

  return density2q;
}


#endif
