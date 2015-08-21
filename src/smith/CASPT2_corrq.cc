//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_corrqq.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_corrq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto corrq = make_shared<Queue>();
  vector<IndexRange> I404_index;
  auto I404 = make_shared<Tensor>(I404_index);
  vector<IndexRange> I405_index = {active_, active_, active_, active_};
  auto I405 = make_shared<Tensor>(I405_index);
  auto tensor277 = vector<shared_ptr<Tensor>>{I404, Gamma92_(), I405};
  auto task277 = make_shared<Task277>(tensor277, pindex);
  corrq->add_task(task277);

  auto tensor278 = vector<shared_ptr<Tensor>>{I405, t2};
  auto task278 = make_shared<Task278>(tensor278, pindex);
  task277->add_dep(task278);
  corrq->add_task(task278);

  vector<IndexRange> I408_index = {closed_, active_, active_, active_};
  auto I408 = make_shared<Tensor>(I408_index);
  auto tensor279 = vector<shared_ptr<Tensor>>{I404, t2, I408};
  auto task279 = make_shared<Task279>(tensor279, pindex);
  task277->add_dep(task279);
  corrq->add_task(task279);

  auto tensor280 = vector<shared_ptr<Tensor>>{I408, Gamma6_(), t2};
  auto task280 = make_shared<Task280>(tensor280, pindex);
  task279->add_dep(task280);
  corrq->add_task(task280);

  vector<IndexRange> I411_index = {active_, active_};
  auto I411 = make_shared<Tensor>(I411_index);
  auto tensor281 = vector<shared_ptr<Tensor>>{I404, Gamma16_(), I411};
  auto task281 = make_shared<Task281>(tensor281, pindex);
  task277->add_dep(task281);
  corrq->add_task(task281);

  auto tensor282 = vector<shared_ptr<Tensor>>{I411, t2};
  auto task282 = make_shared<Task282>(tensor282, pindex);
  task281->add_dep(task282);
  corrq->add_task(task282);

  auto tensor283 = vector<shared_ptr<Tensor>>{I411, t2};
  auto task283 = make_shared<Task283>(tensor283, pindex);
  task281->add_dep(task283);
  corrq->add_task(task283);

  vector<IndexRange> I417_index = {active_, active_, active_, active_};
  auto I417 = make_shared<Tensor>(I417_index);
  auto tensor284 = vector<shared_ptr<Tensor>>{I404, Gamma32_(), I417};
  auto task284 = make_shared<Task284>(tensor284, pindex);
  task277->add_dep(task284);
  corrq->add_task(task284);

  auto tensor285 = vector<shared_ptr<Tensor>>{I417, t2};
  auto task285 = make_shared<Task285>(tensor285, pindex);
  task284->add_dep(task285);
  corrq->add_task(task285);

  vector<IndexRange> I420_index = {active_, active_, active_, active_};
  auto I420 = make_shared<Tensor>(I420_index);
  auto tensor286 = vector<shared_ptr<Tensor>>{I404, Gamma35_(), I420};
  auto task286 = make_shared<Task286>(tensor286, pindex);
  task277->add_dep(task286);
  corrq->add_task(task286);

  auto tensor287 = vector<shared_ptr<Tensor>>{I420, t2};
  auto task287 = make_shared<Task287>(tensor287, pindex);
  task286->add_dep(task287);
  corrq->add_task(task287);

  auto tensor288 = vector<shared_ptr<Tensor>>{I420, t2};
  auto task288 = make_shared<Task288>(tensor288, pindex);
  task286->add_dep(task288);
  corrq->add_task(task288);

  auto tensor289 = vector<shared_ptr<Tensor>>{I420, t2};
  auto task289 = make_shared<Task289>(tensor289, pindex);
  task286->add_dep(task289);
  corrq->add_task(task289);

  vector<IndexRange> I429_index = {active_, active_, active_, active_, active_, active_};
  auto I429 = make_shared<Tensor>(I429_index);
  auto tensor290 = vector<shared_ptr<Tensor>>{I404, Gamma59_(), I429};
  auto task290 = make_shared<Task290>(tensor290, pindex);
  task277->add_dep(task290);
  corrq->add_task(task290);

  auto tensor291 = vector<shared_ptr<Tensor>>{I429, t2};
  auto task291 = make_shared<Task291>(tensor291, pindex);
  task290->add_dep(task291);
  corrq->add_task(task291);

  auto tensor292 = vector<shared_ptr<Tensor>>{I404, t2};
  auto task292 = make_shared<Task292>(tensor292, pindex);
  task277->add_dep(task292);
  corrq->add_task(task292);

  auto tensor293 = vector<shared_ptr<Tensor>>{I404, t2};
  auto task293 = make_shared<Task293>(tensor293, pindex);
  task277->add_dep(task293);
  corrq->add_task(task293);

  vector<IndexRange> I436_index = {active_, active_};
  auto I436 = make_shared<Tensor>(I436_index);
  auto tensor294 = vector<shared_ptr<Tensor>>{I404, Gamma38_(), I436};
  auto task294 = make_shared<Task294>(tensor294, pindex);
  task277->add_dep(task294);
  corrq->add_task(task294);

  auto tensor295 = vector<shared_ptr<Tensor>>{I436, t2};
  auto task295 = make_shared<Task295>(tensor295, pindex);
  task294->add_dep(task295);
  corrq->add_task(task295);

  auto tensor296 = vector<shared_ptr<Tensor>>{I436, t2};
  auto task296 = make_shared<Task296>(tensor296, pindex);
  task294->add_dep(task296);
  corrq->add_task(task296);

  vector<IndexRange> I442_index = {active_, active_, active_, active_};
  auto I442 = make_shared<Tensor>(I442_index);
  auto tensor297 = vector<shared_ptr<Tensor>>{I404, Gamma60_(), I442};
  auto task297 = make_shared<Task297>(tensor297, pindex);
  task277->add_dep(task297);
  corrq->add_task(task297);

  auto tensor298 = vector<shared_ptr<Tensor>>{I442, t2};
  auto task298 = make_shared<Task298>(tensor298, pindex);
  task297->add_dep(task298);
  corrq->add_task(task298);

  return corrq;
}


#endif
