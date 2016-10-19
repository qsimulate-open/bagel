//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASA_sourceqq.cc
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


#include <src/smith/casa/CASA.h>
#include <src/smith/casa/CASA_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASA::CASA::make_sourceq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  auto tensor235 = vector<shared_ptr<Tensor>>{s};
  auto task235 = make_shared<Task235>(tensor235, reset);
  sourceq->add_task(task235);

  vector<IndexRange> I330_index = {closed_, closed_, active_, active_};
  auto I330 = make_shared<Tensor>(I330_index);
  auto tensor236 = vector<shared_ptr<Tensor>>{s, I330};
  auto task236 = make_shared<Task236>(tensor236, pindex);
  task236->add_dep(task235);
  sourceq->add_task(task236);

  auto tensor237 = vector<shared_ptr<Tensor>>{I330, Gamma0_(), v2_};
  auto task237 = make_shared<Task237>(tensor237, pindex);
  task236->add_dep(task237);
  task237->add_dep(task235);
  sourceq->add_task(task237);

  vector<IndexRange> I332_index = {closed_, active_, active_, active_};
  auto I332 = make_shared<Tensor>(I332_index);
  auto tensor238 = vector<shared_ptr<Tensor>>{s, I332};
  auto task238 = make_shared<Task238>(tensor238, pindex);
  task238->add_dep(task235);
  sourceq->add_task(task238);

  auto tensor239 = vector<shared_ptr<Tensor>>{I332, Gamma121_(), v2_};
  auto task239 = make_shared<Task239>(tensor239, pindex);
  task238->add_dep(task239);
  task239->add_dep(task235);
  sourceq->add_task(task239);

  auto tensor240 = vector<shared_ptr<Tensor>>{I332, Gamma3_(), v2_};
  auto task240 = make_shared<Task240>(tensor240, pindex);
  task238->add_dep(task240);
  task240->add_dep(task235);
  sourceq->add_task(task240);

  auto tensor241 = vector<shared_ptr<Tensor>>{I332, Gamma5_(), h1_};
  auto task241 = make_shared<Task241>(tensor241, pindex);
  task238->add_dep(task241);
  task241->add_dep(task235);
  sourceq->add_task(task241);

  vector<IndexRange> I336_index = {closed_, closed_, virt_, active_};
  auto I336 = make_shared<Tensor>(I336_index);
  auto tensor242 = vector<shared_ptr<Tensor>>{s, I336};
  auto task242 = make_shared<Task242>(tensor242, pindex);
  task242->add_dep(task235);
  sourceq->add_task(task242);

  vector<IndexRange> I337_index = {closed_, active_, closed_, virt_};
  auto I337 = make_shared<Tensor>(I337_index);
  auto tensor243 = vector<shared_ptr<Tensor>>{I336, Gamma9_(), I337};
  auto task243 = make_shared<Task243>(tensor243, pindex);
  task242->add_dep(task243);
  task243->add_dep(task235);
  sourceq->add_task(task243);

  auto tensor244 = vector<shared_ptr<Tensor>>{I337, v2_};
  auto task244 = make_shared<Task244>(tensor244, pindex);
  task243->add_dep(task244);
  task244->add_dep(task235);
  sourceq->add_task(task244);

  vector<IndexRange> I340_index = {closed_, virt_, active_, active_};
  auto I340 = make_shared<Tensor>(I340_index);
  auto tensor245 = vector<shared_ptr<Tensor>>{s, I340};
  auto task245 = make_shared<Task245>(tensor245, pindex);
  task245->add_dep(task235);
  sourceq->add_task(task245);

  vector<IndexRange> I341_index = {closed_, virt_, active_, active_};
  auto I341 = make_shared<Tensor>(I341_index);
  auto tensor246 = vector<shared_ptr<Tensor>>{I340, Gamma25_(), I341};
  auto task246 = make_shared<Task246>(tensor246, pindex);
  task245->add_dep(task246);
  task246->add_dep(task235);
  sourceq->add_task(task246);

  auto tensor247 = vector<shared_ptr<Tensor>>{I341, v2_};
  auto task247 = make_shared<Task247>(tensor247, pindex);
  task246->add_dep(task247);
  task247->add_dep(task235);
  sourceq->add_task(task247);

  auto tensor248 = vector<shared_ptr<Tensor>>{I340, Gamma26_(), v2_};
  auto task248 = make_shared<Task248>(tensor248, pindex);
  task245->add_dep(task248);
  task248->add_dep(task235);
  sourceq->add_task(task248);

  auto tensor249 = vector<shared_ptr<Tensor>>{I340, Gamma24_(), v2_};
  auto task249 = make_shared<Task249>(tensor249, pindex);
  task245->add_dep(task249);
  task249->add_dep(task235);
  sourceq->add_task(task249);

  auto tensor250 = vector<shared_ptr<Tensor>>{I340, Gamma29_(), h1_};
  auto task250 = make_shared<Task250>(tensor250, pindex);
  task245->add_dep(task250);
  task250->add_dep(task235);
  sourceq->add_task(task250);

  vector<IndexRange> I348_index = {closed_, virt_, active_, active_};
  auto I348 = make_shared<Tensor>(I348_index);
  auto tensor251 = vector<shared_ptr<Tensor>>{s, I348};
  auto task251 = make_shared<Task251>(tensor251, pindex);
  task251->add_dep(task235);
  sourceq->add_task(task251);

  vector<IndexRange> I349_index = {closed_, virt_, active_, active_};
  auto I349 = make_shared<Tensor>(I349_index);
  auto tensor252 = vector<shared_ptr<Tensor>>{I348, Gamma25_(), I349};
  auto task252 = make_shared<Task252>(tensor252, pindex);
  task251->add_dep(task252);
  task252->add_dep(task235);
  sourceq->add_task(task252);

  auto tensor253 = vector<shared_ptr<Tensor>>{I349, v2_};
  auto task253 = make_shared<Task253>(tensor253, pindex);
  task252->add_dep(task253);
  task253->add_dep(task235);
  sourceq->add_task(task253);

  auto tensor254 = vector<shared_ptr<Tensor>>{I348, Gamma5_(), v2_};
  auto task254 = make_shared<Task254>(tensor254, pindex);
  task251->add_dep(task254);
  task254->add_dep(task235);
  sourceq->add_task(task254);

  auto tensor255 = vector<shared_ptr<Tensor>>{I348, Gamma29_(), h1_};
  auto task255 = make_shared<Task255>(tensor255, pindex);
  task251->add_dep(task255);
  task255->add_dep(task235);
  sourceq->add_task(task255);

  vector<IndexRange> I356_index = {virt_, active_, active_, active_};
  auto I356 = make_shared<Tensor>(I356_index);
  auto tensor256 = vector<shared_ptr<Tensor>>{s, I356};
  auto task256 = make_shared<Task256>(tensor256, pindex);
  task256->add_dep(task235);
  sourceq->add_task(task256);

  auto tensor257 = vector<shared_ptr<Tensor>>{I356, Gamma52_(), v2_};
  auto task257 = make_shared<Task257>(tensor257, pindex);
  task256->add_dep(task257);
  task257->add_dep(task235);
  sourceq->add_task(task257);

  auto tensor258 = vector<shared_ptr<Tensor>>{I356, Gamma49_(), v2_};
  auto task258 = make_shared<Task258>(tensor258, pindex);
  task256->add_dep(task258);
  task258->add_dep(task235);
  sourceq->add_task(task258);

  auto tensor259 = vector<shared_ptr<Tensor>>{I356, Gamma50_(), h1_};
  auto task259 = make_shared<Task259>(tensor259, pindex);
  task256->add_dep(task259);
  task259->add_dep(task235);
  sourceq->add_task(task259);

  shared_ptr<Tensor> I360;
  if (diagonal) {
    vector<IndexRange> I360_index = {closed_, virt_, closed_, virt_};
    I360 = make_shared<Tensor>(I360_index);
  }
  shared_ptr<Task260> task260;
  if (diagonal) {
    auto tensor260 = vector<shared_ptr<Tensor>>{s, I360};
    task260 = make_shared<Task260>(tensor260, pindex);
    task260->add_dep(task235);
    sourceq->add_task(task260);
  }

  shared_ptr<Task261> task261;
  if (diagonal) {
    auto tensor261 = vector<shared_ptr<Tensor>>{I360, v2_};
    task261 = make_shared<Task261>(tensor261, pindex);
    task260->add_dep(task261);
    task261->add_dep(task235);
    sourceq->add_task(task261);
  }

  vector<IndexRange> I362_index = {virt_, closed_, virt_, active_};
  auto I362 = make_shared<Tensor>(I362_index);
  auto tensor262 = vector<shared_ptr<Tensor>>{s, I362};
  auto task262 = make_shared<Task262>(tensor262, pindex);
  task262->add_dep(task235);
  sourceq->add_task(task262);

  vector<IndexRange> I363_index = {active_, virt_, closed_, virt_};
  auto I363 = make_shared<Tensor>(I363_index);
  auto tensor263 = vector<shared_ptr<Tensor>>{I362, Gamma29_(), I363};
  auto task263 = make_shared<Task263>(tensor263, pindex);
  task262->add_dep(task263);
  task263->add_dep(task235);
  sourceq->add_task(task263);

  auto tensor264 = vector<shared_ptr<Tensor>>{I363, v2_};
  auto task264 = make_shared<Task264>(tensor264, pindex);
  task263->add_dep(task264);
  task264->add_dep(task235);
  sourceq->add_task(task264);

  vector<IndexRange> I366_index = {virt_, virt_, active_, active_};
  auto I366 = make_shared<Tensor>(I366_index);
  auto tensor265 = vector<shared_ptr<Tensor>>{s, I366};
  auto task265 = make_shared<Task265>(tensor265, pindex);
  task265->add_dep(task235);
  sourceq->add_task(task265);

  auto tensor266 = vector<shared_ptr<Tensor>>{I366, Gamma50_(), v2_};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task265->add_dep(task266);
  task266->add_dep(task235);
  sourceq->add_task(task266);

  return sourceq;
}


#endif
