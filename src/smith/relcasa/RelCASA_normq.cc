//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASA_normq.cc
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
  auto tensor255 = vector<shared_ptr<Tensor>>{n};
  auto task255 = make_shared<Task255>(tensor255, reset);
  normq->add_task(task255);

  vector<IndexRange> I350_index = {closed_, closed_, active_, active_};
  auto I350 = make_shared<Tensor>(I350_index);
  auto tensor256 = vector<shared_ptr<Tensor>>{n, I350};
  auto task256 = make_shared<Task256>(tensor256, pindex);
  task256->add_dep(task255);
  normq->add_task(task256);

  auto tensor257 = vector<shared_ptr<Tensor>>{I350, Gamma0_(), t2};
  auto task257 = make_shared<Task257>(tensor257, pindex);
  task256->add_dep(task257);
  task257->add_dep(task255);
  normq->add_task(task257);

  vector<IndexRange> I352_index = {closed_, active_, active_, active_};
  auto I352 = make_shared<Tensor>(I352_index);
  auto tensor258 = vector<shared_ptr<Tensor>>{n, I352};
  auto task258 = make_shared<Task258>(tensor258, pindex);
  task258->add_dep(task255);
  normq->add_task(task258);

  auto tensor259 = vector<shared_ptr<Tensor>>{I352, Gamma3_(), t2};
  auto task259 = make_shared<Task259>(tensor259, pindex);
  task258->add_dep(task259);
  task259->add_dep(task255);
  normq->add_task(task259);

  vector<IndexRange> I354_index = {closed_, virt_, closed_, active_};
  auto I354 = make_shared<Tensor>(I354_index);
  auto tensor260 = vector<shared_ptr<Tensor>>{n, I354};
  auto task260 = make_shared<Task260>(tensor260, pindex);
  task260->add_dep(task255);
  normq->add_task(task260);

  vector<IndexRange> I355_index = {closed_, virt_, closed_, active_};
  auto I355 = make_shared<Tensor>(I355_index);
  auto tensor261 = vector<shared_ptr<Tensor>>{I354, Gamma9_(), I355};
  auto task261 = make_shared<Task261>(tensor261, pindex);
  task260->add_dep(task261);
  task261->add_dep(task255);
  normq->add_task(task261);

  auto tensor262 = vector<shared_ptr<Tensor>>{I355, t2};
  auto task262 = make_shared<Task262>(tensor262, pindex);
  task261->add_dep(task262);
  task262->add_dep(task255);
  normq->add_task(task262);

  vector<IndexRange> I358_index = {virt_, closed_, active_, active_};
  auto I358 = make_shared<Tensor>(I358_index);
  auto tensor263 = vector<shared_ptr<Tensor>>{n, I358};
  auto task263 = make_shared<Task263>(tensor263, pindex);
  task263->add_dep(task255);
  normq->add_task(task263);

  auto tensor264 = vector<shared_ptr<Tensor>>{I358, Gamma24_(), t2};
  auto task264 = make_shared<Task264>(tensor264, pindex);
  task263->add_dep(task264);
  task264->add_dep(task255);
  normq->add_task(task264);

  auto tensor265 = vector<shared_ptr<Tensor>>{I358, Gamma25_(), t2};
  auto task265 = make_shared<Task265>(tensor265, pindex);
  task263->add_dep(task265);
  task265->add_dep(task255);
  normq->add_task(task265);

  vector<IndexRange> I362_index = {virt_, closed_, active_, active_};
  auto I362 = make_shared<Tensor>(I362_index);
  auto tensor266 = vector<shared_ptr<Tensor>>{n, I362};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task266->add_dep(task255);
  normq->add_task(task266);

  vector<IndexRange> I363_index = {active_, virt_, closed_, active_};
  auto I363 = make_shared<Tensor>(I363_index);
  auto tensor267 = vector<shared_ptr<Tensor>>{I362, Gamma25_(), I363};
  auto task267 = make_shared<Task267>(tensor267, pindex);
  task266->add_dep(task267);
  task267->add_dep(task255);
  normq->add_task(task267);

  auto tensor268 = vector<shared_ptr<Tensor>>{I363, t2};
  auto task268 = make_shared<Task268>(tensor268, pindex);
  task267->add_dep(task268);
  task268->add_dep(task255);
  normq->add_task(task268);

  vector<IndexRange> I366_index = {virt_, active_, active_, active_};
  auto I366 = make_shared<Tensor>(I366_index);
  auto tensor269 = vector<shared_ptr<Tensor>>{n, I366};
  auto task269 = make_shared<Task269>(tensor269, pindex);
  task269->add_dep(task255);
  normq->add_task(task269);

  auto tensor270 = vector<shared_ptr<Tensor>>{I366, Gamma52_(), t2};
  auto task270 = make_shared<Task270>(tensor270, pindex);
  task269->add_dep(task270);
  task270->add_dep(task255);
  normq->add_task(task270);

  shared_ptr<Tensor> I368;
  if (diagonal) {
    vector<IndexRange> I368_index = {closed_, virt_, closed_, virt_};
    I368 = make_shared<Tensor>(I368_index);
  }
  shared_ptr<Task271> task271;
  if (diagonal) {
    auto tensor271 = vector<shared_ptr<Tensor>>{n, I368};
    task271 = make_shared<Task271>(tensor271, pindex);
    task271->add_dep(task255);
    normq->add_task(task271);
  }

  shared_ptr<Task272> task272;
  if (diagonal) {
    auto tensor272 = vector<shared_ptr<Tensor>>{I368, t2};
    task272 = make_shared<Task272>(tensor272, pindex);
    task271->add_dep(task272);
    task272->add_dep(task255);
    normq->add_task(task272);
  }

  vector<IndexRange> I370_index = {virt_, closed_, virt_, active_};
  auto I370 = make_shared<Tensor>(I370_index);
  auto tensor273 = vector<shared_ptr<Tensor>>{n, I370};
  auto task273 = make_shared<Task273>(tensor273, pindex);
  task273->add_dep(task255);
  normq->add_task(task273);

  vector<IndexRange> I371_index = {active_, virt_, closed_, virt_};
  auto I371 = make_shared<Tensor>(I371_index);
  auto tensor274 = vector<shared_ptr<Tensor>>{I370, Gamma29_(), I371};
  auto task274 = make_shared<Task274>(tensor274, pindex);
  task273->add_dep(task274);
  task274->add_dep(task255);
  normq->add_task(task274);

  auto tensor275 = vector<shared_ptr<Tensor>>{I371, t2};
  auto task275 = make_shared<Task275>(tensor275, pindex);
  task274->add_dep(task275);
  task275->add_dep(task255);
  normq->add_task(task275);

  vector<IndexRange> I374_index = {virt_, virt_, active_, active_};
  auto I374 = make_shared<Tensor>(I374_index);
  auto tensor276 = vector<shared_ptr<Tensor>>{n, I374};
  auto task276 = make_shared<Task276>(tensor276, pindex);
  task276->add_dep(task255);
  normq->add_task(task276);

  auto tensor277 = vector<shared_ptr<Tensor>>{I374, Gamma50_(), t2};
  auto task277 = make_shared<Task277>(tensor277, pindex);
  task276->add_dep(task277);
  task277->add_dep(task255);
  normq->add_task(task277);

  return normq;
}


#endif
