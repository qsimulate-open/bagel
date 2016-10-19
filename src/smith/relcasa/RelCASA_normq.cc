//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASA_normqq.cc
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


#include <src/smith/relcasa/RelCASA.h>
#include <src/smith/relcasa/RelCASA_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelCASA::RelCASA::make_normq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto normq = make_shared<Queue>();
  auto tensor267 = vector<shared_ptr<Tensor>>{n};
  auto task267 = make_shared<Task267>(tensor267, reset);
  normq->add_task(task267);

  vector<IndexRange> I376_index = {closed_, closed_, active_, active_};
  auto I376 = make_shared<Tensor>(I376_index);
  auto tensor268 = vector<shared_ptr<Tensor>>{n, I376};
  auto task268 = make_shared<Task268>(tensor268, pindex);
  task268->add_dep(task267);
  normq->add_task(task268);

  auto tensor269 = vector<shared_ptr<Tensor>>{I376, Gamma0_(), t2};
  auto task269 = make_shared<Task269>(tensor269, pindex);
  task268->add_dep(task269);
  task269->add_dep(task267);
  normq->add_task(task269);

  vector<IndexRange> I378_index = {closed_, active_, active_, active_};
  auto I378 = make_shared<Tensor>(I378_index);
  auto tensor270 = vector<shared_ptr<Tensor>>{n, I378};
  auto task270 = make_shared<Task270>(tensor270, pindex);
  task270->add_dep(task267);
  normq->add_task(task270);

  auto tensor271 = vector<shared_ptr<Tensor>>{I378, Gamma3_(), t2};
  auto task271 = make_shared<Task271>(tensor271, pindex);
  task270->add_dep(task271);
  task271->add_dep(task267);
  normq->add_task(task271);

  vector<IndexRange> I380_index = {closed_, virt_, closed_, active_};
  auto I380 = make_shared<Tensor>(I380_index);
  auto tensor272 = vector<shared_ptr<Tensor>>{n, I380};
  auto task272 = make_shared<Task272>(tensor272, pindex);
  task272->add_dep(task267);
  normq->add_task(task272);

  vector<IndexRange> I381_index = {closed_, virt_, closed_, active_};
  auto I381 = make_shared<Tensor>(I381_index);
  auto tensor273 = vector<shared_ptr<Tensor>>{I380, Gamma9_(), I381};
  auto task273 = make_shared<Task273>(tensor273, pindex);
  task272->add_dep(task273);
  task273->add_dep(task267);
  normq->add_task(task273);

  auto tensor274 = vector<shared_ptr<Tensor>>{I381, t2};
  auto task274 = make_shared<Task274>(tensor274, pindex);
  task273->add_dep(task274);
  task274->add_dep(task267);
  normq->add_task(task274);

  vector<IndexRange> I384_index = {virt_, closed_, active_, active_};
  auto I384 = make_shared<Tensor>(I384_index);
  auto tensor275 = vector<shared_ptr<Tensor>>{n, I384};
  auto task275 = make_shared<Task275>(tensor275, pindex);
  task275->add_dep(task267);
  normq->add_task(task275);

  auto tensor276 = vector<shared_ptr<Tensor>>{I384, Gamma24_(), t2};
  auto task276 = make_shared<Task276>(tensor276, pindex);
  task275->add_dep(task276);
  task276->add_dep(task267);
  normq->add_task(task276);

  auto tensor277 = vector<shared_ptr<Tensor>>{I384, Gamma25_(), t2};
  auto task277 = make_shared<Task277>(tensor277, pindex);
  task275->add_dep(task277);
  task277->add_dep(task267);
  normq->add_task(task277);

  vector<IndexRange> I388_index = {virt_, closed_, active_, active_};
  auto I388 = make_shared<Tensor>(I388_index);
  auto tensor278 = vector<shared_ptr<Tensor>>{n, I388};
  auto task278 = make_shared<Task278>(tensor278, pindex);
  task278->add_dep(task267);
  normq->add_task(task278);

  vector<IndexRange> I389_index = {active_, virt_, closed_, active_};
  auto I389 = make_shared<Tensor>(I389_index);
  auto tensor279 = vector<shared_ptr<Tensor>>{I388, Gamma25_(), I389};
  auto task279 = make_shared<Task279>(tensor279, pindex);
  task278->add_dep(task279);
  task279->add_dep(task267);
  normq->add_task(task279);

  auto tensor280 = vector<shared_ptr<Tensor>>{I389, t2};
  auto task280 = make_shared<Task280>(tensor280, pindex);
  task279->add_dep(task280);
  task280->add_dep(task267);
  normq->add_task(task280);

  vector<IndexRange> I392_index = {virt_, active_, active_, active_};
  auto I392 = make_shared<Tensor>(I392_index);
  auto tensor281 = vector<shared_ptr<Tensor>>{n, I392};
  auto task281 = make_shared<Task281>(tensor281, pindex);
  task281->add_dep(task267);
  normq->add_task(task281);

  auto tensor282 = vector<shared_ptr<Tensor>>{I392, Gamma52_(), t2};
  auto task282 = make_shared<Task282>(tensor282, pindex);
  task281->add_dep(task282);
  task282->add_dep(task267);
  normq->add_task(task282);

  shared_ptr<Tensor> I394;
  if (diagonal) {
    vector<IndexRange> I394_index = {closed_, virt_, closed_, virt_};
    I394 = make_shared<Tensor>(I394_index);
  }
  shared_ptr<Task283> task283;
  if (diagonal) {
    auto tensor283 = vector<shared_ptr<Tensor>>{n, I394};
    task283 = make_shared<Task283>(tensor283, pindex);
    task283->add_dep(task267);
    normq->add_task(task283);
  }

  shared_ptr<Task284> task284;
  if (diagonal) {
    auto tensor284 = vector<shared_ptr<Tensor>>{I394, t2};
    task284 = make_shared<Task284>(tensor284, pindex);
    task283->add_dep(task284);
    task284->add_dep(task267);
    normq->add_task(task284);
  }

  vector<IndexRange> I396_index = {virt_, closed_, virt_, active_};
  auto I396 = make_shared<Tensor>(I396_index);
  auto tensor285 = vector<shared_ptr<Tensor>>{n, I396};
  auto task285 = make_shared<Task285>(tensor285, pindex);
  task285->add_dep(task267);
  normq->add_task(task285);

  vector<IndexRange> I397_index = {active_, virt_, closed_, virt_};
  auto I397 = make_shared<Tensor>(I397_index);
  auto tensor286 = vector<shared_ptr<Tensor>>{I396, Gamma29_(), I397};
  auto task286 = make_shared<Task286>(tensor286, pindex);
  task285->add_dep(task286);
  task286->add_dep(task267);
  normq->add_task(task286);

  auto tensor287 = vector<shared_ptr<Tensor>>{I397, t2};
  auto task287 = make_shared<Task287>(tensor287, pindex);
  task286->add_dep(task287);
  task287->add_dep(task267);
  normq->add_task(task287);

  vector<IndexRange> I400_index = {virt_, virt_, active_, active_};
  auto I400 = make_shared<Tensor>(I400_index);
  auto tensor288 = vector<shared_ptr<Tensor>>{n, I400};
  auto task288 = make_shared<Task288>(tensor288, pindex);
  task288->add_dep(task267);
  normq->add_task(task288);

  auto tensor289 = vector<shared_ptr<Tensor>>{I400, Gamma50_(), t2};
  auto task289 = make_shared<Task289>(tensor289, pindex);
  task288->add_dep(task289);
  task289->add_dep(task267);
  normq->add_task(task289);

  return normq;
}


#endif
