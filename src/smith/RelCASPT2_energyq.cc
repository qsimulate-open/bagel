//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2_energyqq.cc
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


#include <src/smith/RelCASPT2.h>
#include <src/smith/RelCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelCASPT2::RelCASPT2::make_energyq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto energyq = make_shared<Queue>();
  vector<IndexRange> I348_index;
  auto I348 = make_shared<Tensor>(I348_index);
  vector<IndexRange> I349_index = {active_, active_, active_, active_};
  auto I349 = make_shared<Tensor>(I349_index);
  auto tensor209 = vector<shared_ptr<Tensor>>{I348, Gamma94_(), I349};
  auto task209 = make_shared<Task209>(tensor209, pindex);
  energyq->add_task(task209);

  auto tensor210 = vector<shared_ptr<Tensor>>{I349, t2, v2_};
  auto task210 = make_shared<Task210>(tensor210, pindex);
  task209->add_dep(task210);
  energyq->add_task(task210);

  vector<IndexRange> I352_index = {closed_, active_, active_, active_};
  auto I352 = make_shared<Tensor>(I352_index);
  auto tensor211 = vector<shared_ptr<Tensor>>{I348, t2, I352};
  auto task211 = make_shared<Task211>(tensor211, pindex);
  task209->add_dep(task211);
  energyq->add_task(task211);

  auto tensor212 = vector<shared_ptr<Tensor>>{I352, Gamma107_(), v2_};
  auto task212 = make_shared<Task212>(tensor212, pindex);
  task211->add_dep(task212);
  energyq->add_task(task212);

  vector<IndexRange> I355_index = {closed_, active_, active_, active_};
  auto I355 = make_shared<Tensor>(I355_index);
  auto tensor213 = vector<shared_ptr<Tensor>>{I348, v2_, I355};
  auto task213 = make_shared<Task213>(tensor213, pindex);
  task209->add_dep(task213);
  energyq->add_task(task213);

  auto tensor214 = vector<shared_ptr<Tensor>>{I355, Gamma6_(), t2};
  auto task214 = make_shared<Task214>(tensor214, pindex);
  task213->add_dep(task214);
  energyq->add_task(task214);

  vector<IndexRange> I358_index = {active_, active_};
  auto I358 = make_shared<Tensor>(I358_index);
  auto tensor215 = vector<shared_ptr<Tensor>>{I348, Gamma16_(), I358};
  auto task215 = make_shared<Task215>(tensor215, pindex);
  task209->add_dep(task215);
  energyq->add_task(task215);

  vector<IndexRange> I359_index = {closed_, active_, closed_, virt_};
  auto I359 = make_shared<Tensor>(I359_index);
  auto tensor216 = vector<shared_ptr<Tensor>>{I358, t2, I359};
  auto task216 = make_shared<Task216>(tensor216, pindex);
  task215->add_dep(task216);
  energyq->add_task(task216);

  auto tensor217 = vector<shared_ptr<Tensor>>{I359, v2_};
  auto task217 = make_shared<Task217>(tensor217, pindex);
  task216->add_dep(task217);
  energyq->add_task(task217);

  vector<IndexRange> I364_index = {active_, active_, active_, active_};
  auto I364 = make_shared<Tensor>(I364_index);
  auto tensor218 = vector<shared_ptr<Tensor>>{I348, Gamma35_(), I364};
  auto task218 = make_shared<Task218>(tensor218, pindex);
  task209->add_dep(task218);
  energyq->add_task(task218);

  auto tensor219 = vector<shared_ptr<Tensor>>{I364, v2_, t2};
  auto task219 = make_shared<Task219>(tensor219, pindex);
  task218->add_dep(task219);
  energyq->add_task(task219);

  auto tensor220 = vector<shared_ptr<Tensor>>{I364, v2_, t2};
  auto task220 = make_shared<Task220>(tensor220, pindex);
  task218->add_dep(task220);
  energyq->add_task(task220);

  auto tensor221 = vector<shared_ptr<Tensor>>{I364, v2_, t2};
  auto task221 = make_shared<Task221>(tensor221, pindex);
  task218->add_dep(task221);
  energyq->add_task(task221);

  auto tensor222 = vector<shared_ptr<Tensor>>{I364, v2_, t2};
  auto task222 = make_shared<Task222>(tensor222, pindex);
  task218->add_dep(task222);
  energyq->add_task(task222);

  auto tensor223 = vector<shared_ptr<Tensor>>{I364, v2_, t2};
  auto task223 = make_shared<Task223>(tensor223, pindex);
  task218->add_dep(task223);
  energyq->add_task(task223);

  vector<IndexRange> I367_index = {active_, active_, active_, active_};
  auto I367 = make_shared<Tensor>(I367_index);
  auto tensor224 = vector<shared_ptr<Tensor>>{I348, Gamma29_(), I367};
  auto task224 = make_shared<Task224>(tensor224, pindex);
  task209->add_dep(task224);
  energyq->add_task(task224);

  auto tensor225 = vector<shared_ptr<Tensor>>{I367, v2_, t2};
  auto task225 = make_shared<Task225>(tensor225, pindex);
  task224->add_dep(task225);
  energyq->add_task(task225);

  vector<IndexRange> I370_index = {active_, active_, active_, active_};
  auto I370 = make_shared<Tensor>(I370_index);
  auto tensor226 = vector<shared_ptr<Tensor>>{I348, Gamma32_(), I370};
  auto task226 = make_shared<Task226>(tensor226, pindex);
  task209->add_dep(task226);
  energyq->add_task(task226);

  auto tensor227 = vector<shared_ptr<Tensor>>{I370, t2, v2_};
  auto task227 = make_shared<Task227>(tensor227, pindex);
  task226->add_dep(task227);
  energyq->add_task(task227);

  vector<IndexRange> I379_index = {active_, active_, active_, active_};
  auto I379 = make_shared<Tensor>(I379_index);
  auto tensor228 = vector<shared_ptr<Tensor>>{I348, Gamma7_(), I379};
  auto task228 = make_shared<Task228>(tensor228, pindex);
  task209->add_dep(task228);
  energyq->add_task(task228);

  auto tensor229 = vector<shared_ptr<Tensor>>{I379, t2, v2_};
  auto task229 = make_shared<Task229>(tensor229, pindex);
  task228->add_dep(task229);
  energyq->add_task(task229);

  vector<IndexRange> I388_index = {active_, active_, active_, active_, active_, active_};
  auto I388 = make_shared<Tensor>(I388_index);
  auto tensor230 = vector<shared_ptr<Tensor>>{I348, Gamma59_(), I388};
  auto task230 = make_shared<Task230>(tensor230, pindex);
  task209->add_dep(task230);
  energyq->add_task(task230);

  auto tensor231 = vector<shared_ptr<Tensor>>{I388, t2, v2_};
  auto task231 = make_shared<Task231>(tensor231, pindex);
  task230->add_dep(task231);
  energyq->add_task(task231);

  vector<IndexRange> I391_index = {active_, active_, active_, active_, active_, active_};
  auto I391 = make_shared<Tensor>(I391_index);
  auto tensor232 = vector<shared_ptr<Tensor>>{I348, Gamma57_(), I391};
  auto task232 = make_shared<Task232>(tensor232, pindex);
  task209->add_dep(task232);
  energyq->add_task(task232);

  auto tensor233 = vector<shared_ptr<Tensor>>{I391, v2_, t2};
  auto task233 = make_shared<Task233>(tensor233, pindex);
  task232->add_dep(task233);
  energyq->add_task(task233);

  auto tensor234 = vector<shared_ptr<Tensor>>{I348, v2_, t2};
  auto task234 = make_shared<Task234>(tensor234, pindex);
  task209->add_dep(task234);
  energyq->add_task(task234);

  auto tensor235 = vector<shared_ptr<Tensor>>{I348, v2_, t2};
  auto task235 = make_shared<Task235>(tensor235, pindex);
  task209->add_dep(task235);
  energyq->add_task(task235);

  vector<IndexRange> I398_index = {active_, active_};
  auto I398 = make_shared<Tensor>(I398_index);
  auto tensor236 = vector<shared_ptr<Tensor>>{I348, Gamma38_(), I398};
  auto task236 = make_shared<Task236>(tensor236, pindex);
  task209->add_dep(task236);
  energyq->add_task(task236);

  auto tensor237 = vector<shared_ptr<Tensor>>{I398, t2, v2_};
  auto task237 = make_shared<Task237>(tensor237, pindex);
  task236->add_dep(task237);
  energyq->add_task(task237);

  auto tensor238 = vector<shared_ptr<Tensor>>{I398, v2_, t2};
  auto task238 = make_shared<Task238>(tensor238, pindex);
  task236->add_dep(task238);
  energyq->add_task(task238);

  auto tensor239 = vector<shared_ptr<Tensor>>{I398, t2, h1_};
  auto task239 = make_shared<Task239>(tensor239, pindex);
  task236->add_dep(task239);
  energyq->add_task(task239);

  auto tensor240 = vector<shared_ptr<Tensor>>{I398, t2, h1_};
  auto task240 = make_shared<Task240>(tensor240, pindex);
  task236->add_dep(task240);
  energyq->add_task(task240);

  vector<IndexRange> I404_index = {active_, active_, active_, active_};
  auto I404 = make_shared<Tensor>(I404_index);
  auto tensor241 = vector<shared_ptr<Tensor>>{I348, Gamma60_(), I404};
  auto task241 = make_shared<Task241>(tensor241, pindex);
  task209->add_dep(task241);
  energyq->add_task(task241);

  auto tensor242 = vector<shared_ptr<Tensor>>{I404, t2, v2_};
  auto task242 = make_shared<Task242>(tensor242, pindex);
  task241->add_dep(task242);
  energyq->add_task(task242);

  auto tensor243 = vector<shared_ptr<Tensor>>{I404, h1_, t2};
  auto task243 = make_shared<Task243>(tensor243, pindex);
  task241->add_dep(task243);
  energyq->add_task(task243);

  vector<IndexRange> I407_index = {active_, closed_};
  auto I407 = make_shared<Tensor>(I407_index);
  auto tensor244 = vector<shared_ptr<Tensor>>{I348, h1_, I407};
  auto task244 = make_shared<Task244>(tensor244, pindex);
  task209->add_dep(task244);
  energyq->add_task(task244);

  auto tensor245 = vector<shared_ptr<Tensor>>{I407, t2, Gamma7_()};
  auto task245 = make_shared<Task245>(tensor245, pindex);
  task244->add_dep(task245);
  energyq->add_task(task245);

  return energyq;
}


#endif
