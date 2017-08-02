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
  auto tensor259 = vector<shared_ptr<Tensor>>{Den1};
  auto task259 = make_shared<Task259>(tensor259, reset);
  density2q->add_task(task259);

  vector<IndexRange> I298_index = {closed_, closed_, active_, active_};
  auto I298 = make_shared<Tensor>(I298_index);
  auto tensor260 = vector<shared_ptr<Tensor>>{Den1, I298};
  auto task260 = make_shared<Task260>(tensor260, pindex);
  task260->add_dep(task259);
  density2q->add_task(task260);

  auto tensor261 = vector<shared_ptr<Tensor>>{I298, Gamma1_(), l2};
  auto task261 = make_shared<Task261>(tensor261, pindex);
  task260->add_dep(task261);
  task261->add_dep(task259);
  density2q->add_task(task261);

  vector<IndexRange> I300_index = {closed_, active_, active_, active_};
  auto I300 = make_shared<Tensor>(I300_index);
  auto tensor262 = vector<shared_ptr<Tensor>>{Den1, I300};
  auto task262 = make_shared<Task262>(tensor262, pindex);
  task262->add_dep(task259);
  density2q->add_task(task262);

  auto tensor263 = vector<shared_ptr<Tensor>>{I300, Gamma6_(), l2};
  auto task263 = make_shared<Task263>(tensor263, pindex);
  task262->add_dep(task263);
  task263->add_dep(task259);
  density2q->add_task(task263);

  vector<IndexRange> I302_index = {closed_, virt_, closed_, active_};
  auto I302 = make_shared<Tensor>(I302_index);
  auto tensor264 = vector<shared_ptr<Tensor>>{Den1, I302};
  auto task264 = make_shared<Task264>(tensor264, pindex);
  task264->add_dep(task259);
  density2q->add_task(task264);

  vector<IndexRange> I303_index = {closed_, virt_, closed_, active_};
  auto I303 = make_shared<Tensor>(I303_index);
  auto tensor265 = vector<shared_ptr<Tensor>>{I302, Gamma16_(), I303};
  auto task265 = make_shared<Task265>(tensor265, pindex);
  task264->add_dep(task265);
  task265->add_dep(task259);
  density2q->add_task(task265);

  auto tensor266 = vector<shared_ptr<Tensor>>{I303, l2};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task265->add_dep(task266);
  task266->add_dep(task259);
  density2q->add_task(task266);

  vector<IndexRange> I306_index = {virt_, closed_, active_, active_};
  auto I306 = make_shared<Tensor>(I306_index);
  auto tensor267 = vector<shared_ptr<Tensor>>{Den1, I306};
  auto task267 = make_shared<Task267>(tensor267, pindex);
  task267->add_dep(task259);
  density2q->add_task(task267);

  auto tensor268 = vector<shared_ptr<Tensor>>{I306, Gamma32_(), l2};
  auto task268 = make_shared<Task268>(tensor268, pindex);
  task267->add_dep(task268);
  task268->add_dep(task259);
  density2q->add_task(task268);

  auto tensor269 = vector<shared_ptr<Tensor>>{I306, Gamma35_(), l2};
  auto task269 = make_shared<Task269>(tensor269, pindex);
  task267->add_dep(task269);
  task269->add_dep(task259);
  density2q->add_task(task269);

  vector<IndexRange> I310_index = {virt_, closed_, active_, active_};
  auto I310 = make_shared<Tensor>(I310_index);
  auto tensor270 = vector<shared_ptr<Tensor>>{Den1, I310};
  auto task270 = make_shared<Task270>(tensor270, pindex);
  task270->add_dep(task259);
  density2q->add_task(task270);

  vector<IndexRange> I311_index = {active_, virt_, closed_, active_};
  auto I311 = make_shared<Tensor>(I311_index);
  auto tensor271 = vector<shared_ptr<Tensor>>{I310, Gamma35_(), I311};
  auto task271 = make_shared<Task271>(tensor271, pindex);
  task270->add_dep(task271);
  task271->add_dep(task259);
  density2q->add_task(task271);

  auto tensor272 = vector<shared_ptr<Tensor>>{I311, l2};
  auto task272 = make_shared<Task272>(tensor272, pindex);
  task271->add_dep(task272);
  task272->add_dep(task259);
  density2q->add_task(task272);

  vector<IndexRange> I314_index = {virt_, active_, active_, active_};
  auto I314 = make_shared<Tensor>(I314_index);
  auto tensor273 = vector<shared_ptr<Tensor>>{Den1, I314};
  auto task273 = make_shared<Task273>(tensor273, pindex);
  task273->add_dep(task259);
  density2q->add_task(task273);

  auto tensor274 = vector<shared_ptr<Tensor>>{I314, Gamma62_(), l2};
  auto task274 = make_shared<Task274>(tensor274, pindex);
  task273->add_dep(task274);
  task274->add_dep(task259);
  density2q->add_task(task274);

  shared_ptr<Tensor> I316;
  if (diagonal) {
    vector<IndexRange> I316_index = {closed_, virt_, closed_, virt_};
    I316 = make_shared<Tensor>(I316_index);
  }
  shared_ptr<Task275> task275;
  if (diagonal) {
    auto tensor275 = vector<shared_ptr<Tensor>>{Den1, I316};
    task275 = make_shared<Task275>(tensor275, pindex);
    task275->add_dep(task259);
    density2q->add_task(task275);
  }

  shared_ptr<Task276> task276;
  if (diagonal) {
    auto tensor276 = vector<shared_ptr<Tensor>>{I316, l2};
    task276 = make_shared<Task276>(tensor276, pindex);
    task275->add_dep(task276);
    task276->add_dep(task259);
    density2q->add_task(task276);
  }

  vector<IndexRange> I318_index = {virt_, closed_, virt_, active_};
  auto I318 = make_shared<Tensor>(I318_index);
  auto tensor277 = vector<shared_ptr<Tensor>>{Den1, I318};
  auto task277 = make_shared<Task277>(tensor277, pindex);
  task277->add_dep(task259);
  density2q->add_task(task277);

  vector<IndexRange> I319_index = {active_, virt_, closed_, virt_};
  auto I319 = make_shared<Tensor>(I319_index);
  auto tensor278 = vector<shared_ptr<Tensor>>{I318, Gamma65_(), I319};
  auto task278 = make_shared<Task278>(tensor278, pindex);
  task277->add_dep(task278);
  task278->add_dep(task259);
  density2q->add_task(task278);

  auto tensor279 = vector<shared_ptr<Tensor>>{I319, l2};
  auto task279 = make_shared<Task279>(tensor279, pindex);
  task278->add_dep(task279);
  task279->add_dep(task259);
  density2q->add_task(task279);

  vector<IndexRange> I322_index = {virt_, virt_, active_, active_};
  auto I322 = make_shared<Tensor>(I322_index);
  auto tensor280 = vector<shared_ptr<Tensor>>{Den1, I322};
  auto task280 = make_shared<Task280>(tensor280, pindex);
  task280->add_dep(task259);
  density2q->add_task(task280);

  auto tensor281 = vector<shared_ptr<Tensor>>{I322, Gamma60_(), l2};
  auto task281 = make_shared<Task281>(tensor281, pindex);
  task280->add_dep(task281);
  task281->add_dep(task259);
  density2q->add_task(task281);

  return density2q;
}


#endif
