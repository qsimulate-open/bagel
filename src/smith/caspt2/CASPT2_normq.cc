//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_normqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_normq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto normq = make_shared<Queue>();
  auto tensor274 = vector<shared_ptr<Tensor>>{n};
  auto task274 = make_shared<Task274>(tensor274, reset);
  normq->add_task(task274);

  vector<IndexRange> I334_index = {closed_, closed_, active_, active_};
  auto I334 = make_shared<Tensor>(I334_index);
  auto tensor275 = vector<shared_ptr<Tensor>>{n, I334};
  auto task275 = make_shared<Task275>(tensor275, pindex);
  task275->add_dep(task274);
  normq->add_task(task275);

  auto tensor276 = vector<shared_ptr<Tensor>>{I334, Gamma92_(), t2};
  auto task276 = make_shared<Task276>(tensor276, pindex);
  task275->add_dep(task276);
  task276->add_dep(task274);
  normq->add_task(task276);

  vector<IndexRange> I336_index = {closed_, active_, active_, active_};
  auto I336 = make_shared<Tensor>(I336_index);
  auto tensor277 = vector<shared_ptr<Tensor>>{n, I336};
  auto task277 = make_shared<Task277>(tensor277, pindex);
  task277->add_dep(task274);
  normq->add_task(task277);

  auto tensor278 = vector<shared_ptr<Tensor>>{I336, Gamma6_(), t2};
  auto task278 = make_shared<Task278>(tensor278, pindex);
  task277->add_dep(task278);
  task278->add_dep(task274);
  normq->add_task(task278);

  vector<IndexRange> I338_index = {closed_, virt_, closed_, active_};
  auto I338 = make_shared<Tensor>(I338_index);
  auto tensor279 = vector<shared_ptr<Tensor>>{n, I338};
  auto task279 = make_shared<Task279>(tensor279, pindex);
  task279->add_dep(task274);
  normq->add_task(task279);

  vector<IndexRange> I339_index = {closed_, virt_, closed_, active_};
  auto I339 = make_shared<Tensor>(I339_index);
  auto tensor280 = vector<shared_ptr<Tensor>>{I338, Gamma16_(), I339};
  auto task280 = make_shared<Task280>(tensor280, pindex);
  task279->add_dep(task280);
  task280->add_dep(task274);
  normq->add_task(task280);

  auto tensor281 = vector<shared_ptr<Tensor>>{I339, t2};
  auto task281 = make_shared<Task281>(tensor281, pindex);
  task280->add_dep(task281);
  task281->add_dep(task274);
  normq->add_task(task281);

  vector<IndexRange> I342_index = {virt_, closed_, active_, active_};
  auto I342 = make_shared<Tensor>(I342_index);
  auto tensor282 = vector<shared_ptr<Tensor>>{n, I342};
  auto task282 = make_shared<Task282>(tensor282, pindex);
  task282->add_dep(task274);
  normq->add_task(task282);

  auto tensor283 = vector<shared_ptr<Tensor>>{I342, Gamma32_(), t2};
  auto task283 = make_shared<Task283>(tensor283, pindex);
  task282->add_dep(task283);
  task283->add_dep(task274);
  normq->add_task(task283);

  auto tensor284 = vector<shared_ptr<Tensor>>{I342, Gamma35_(), t2};
  auto task284 = make_shared<Task284>(tensor284, pindex);
  task282->add_dep(task284);
  task284->add_dep(task274);
  normq->add_task(task284);

  vector<IndexRange> I346_index = {virt_, closed_, active_, active_};
  auto I346 = make_shared<Tensor>(I346_index);
  auto tensor285 = vector<shared_ptr<Tensor>>{n, I346};
  auto task285 = make_shared<Task285>(tensor285, pindex);
  task285->add_dep(task274);
  normq->add_task(task285);

  vector<IndexRange> I347_index = {active_, virt_, closed_, active_};
  auto I347 = make_shared<Tensor>(I347_index);
  auto tensor286 = vector<shared_ptr<Tensor>>{I346, Gamma35_(), I347};
  auto task286 = make_shared<Task286>(tensor286, pindex);
  task285->add_dep(task286);
  task286->add_dep(task274);
  normq->add_task(task286);

  auto tensor287 = vector<shared_ptr<Tensor>>{I347, t2};
  auto task287 = make_shared<Task287>(tensor287, pindex);
  task286->add_dep(task287);
  task287->add_dep(task274);
  normq->add_task(task287);

  vector<IndexRange> I350_index = {virt_, active_, active_, active_};
  auto I350 = make_shared<Tensor>(I350_index);
  auto tensor288 = vector<shared_ptr<Tensor>>{n, I350};
  auto task288 = make_shared<Task288>(tensor288, pindex);
  task288->add_dep(task274);
  normq->add_task(task288);

  auto tensor289 = vector<shared_ptr<Tensor>>{I350, Gamma59_(), t2};
  auto task289 = make_shared<Task289>(tensor289, pindex);
  task288->add_dep(task289);
  task289->add_dep(task274);
  normq->add_task(task289);

  shared_ptr<Tensor> I352;
  if (diagonal) {
    vector<IndexRange> I352_index = {closed_, virt_, closed_, virt_};
    I352 = make_shared<Tensor>(I352_index);
  }
  shared_ptr<Task290> task290;
  if (diagonal) {
    auto tensor290 = vector<shared_ptr<Tensor>>{n, I352};
    task290 = make_shared<Task290>(tensor290, pindex);
    task290->add_dep(task274);
    normq->add_task(task290);
  }

  shared_ptr<Task291> task291;
  if (diagonal) {
    auto tensor291 = vector<shared_ptr<Tensor>>{I352, t2};
    task291 = make_shared<Task291>(tensor291, pindex);
    task290->add_dep(task291);
    task291->add_dep(task274);
    normq->add_task(task291);
  }

  vector<IndexRange> I354_index = {virt_, closed_, virt_, active_};
  auto I354 = make_shared<Tensor>(I354_index);
  auto tensor292 = vector<shared_ptr<Tensor>>{n, I354};
  auto task292 = make_shared<Task292>(tensor292, pindex);
  task292->add_dep(task274);
  normq->add_task(task292);

  vector<IndexRange> I355_index = {active_, virt_, closed_, virt_};
  auto I355 = make_shared<Tensor>(I355_index);
  auto tensor293 = vector<shared_ptr<Tensor>>{I354, Gamma38_(), I355};
  auto task293 = make_shared<Task293>(tensor293, pindex);
  task292->add_dep(task293);
  task293->add_dep(task274);
  normq->add_task(task293);

  auto tensor294 = vector<shared_ptr<Tensor>>{I355, t2};
  auto task294 = make_shared<Task294>(tensor294, pindex);
  task293->add_dep(task294);
  task294->add_dep(task274);
  normq->add_task(task294);

  vector<IndexRange> I358_index = {virt_, virt_, active_, active_};
  auto I358 = make_shared<Tensor>(I358_index);
  auto tensor295 = vector<shared_ptr<Tensor>>{n, I358};
  auto task295 = make_shared<Task295>(tensor295, pindex);
  task295->add_dep(task274);
  normq->add_task(task295);

  auto tensor296 = vector<shared_ptr<Tensor>>{I358, Gamma60_(), t2};
  auto task296 = make_shared<Task296>(tensor296, pindex);
  task295->add_dep(task296);
  task296->add_dep(task274);
  normq->add_task(task296);

  return normq;
}


#endif
