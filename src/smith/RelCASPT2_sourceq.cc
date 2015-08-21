//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2_sourceqq.cc
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

shared_ptr<Queue> RelCASPT2::RelCASPT2::make_sourceq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  auto tensor209 = vector<shared_ptr<Tensor>>{s};
  auto task209 = make_shared<Task209>(tensor209, reset);
  sourceq->add_task(task209);

  vector<IndexRange> I348_index = {active_, active_, active_, closed_};
  auto I348 = make_shared<Tensor>(I348_index);
  auto tensor210 = vector<shared_ptr<Tensor>>{s, I348};
  auto task210 = make_shared<Task210>(tensor210, pindex);
  task210->add_dep(task209);
  sourceq->add_task(task210);

  auto tensor211 = vector<shared_ptr<Tensor>>{I348, h1_, Gamma7_()};
  auto task211 = make_shared<Task211>(tensor211, pindex);
  task210->add_dep(task211);
  task211->add_dep(task209);
  sourceq->add_task(task211);

  auto tensor212 = vector<shared_ptr<Tensor>>{I348, v2_, Gamma107_()};
  auto task212 = make_shared<Task212>(tensor212, pindex);
  task210->add_dep(task212);
  task212->add_dep(task209);
  sourceq->add_task(task212);

  auto tensor213 = vector<shared_ptr<Tensor>>{I348, Gamma6_(), v2_};
  auto task213 = make_shared<Task213>(tensor213, pindex);
  task210->add_dep(task213);
  task213->add_dep(task209);
  sourceq->add_task(task213);

  vector<IndexRange> I350_index = {closed_, virt_, active_, active_};
  auto I350 = make_shared<Tensor>(I350_index);
  auto tensor214 = vector<shared_ptr<Tensor>>{s, I350};
  auto task214 = make_shared<Task214>(tensor214, pindex);
  task214->add_dep(task209);
  sourceq->add_task(task214);

  auto tensor215 = vector<shared_ptr<Tensor>>{I350, Gamma38_(), h1_};
  auto task215 = make_shared<Task215>(tensor215, pindex);
  task214->add_dep(task215);
  task215->add_dep(task209);
  sourceq->add_task(task215);

  auto tensor216 = vector<shared_ptr<Tensor>>{I350, v2_, Gamma35_()};
  auto task216 = make_shared<Task216>(tensor216, pindex);
  task214->add_dep(task216);
  task216->add_dep(task209);
  sourceq->add_task(task216);

  auto tensor217 = vector<shared_ptr<Tensor>>{I350, v2_, Gamma29_()};
  auto task217 = make_shared<Task217>(tensor217, pindex);
  task214->add_dep(task217);
  task217->add_dep(task209);
  sourceq->add_task(task217);

  auto tensor218 = vector<shared_ptr<Tensor>>{I350, Gamma32_(), v2_};
  auto task218 = make_shared<Task218>(tensor218, pindex);
  task214->add_dep(task218);
  task218->add_dep(task209);
  sourceq->add_task(task218);

  auto tensor219 = vector<shared_ptr<Tensor>>{I350, v2_, Gamma35_()};
  auto task219 = make_shared<Task219>(tensor219, pindex);
  task214->add_dep(task219);
  task219->add_dep(task209);
  sourceq->add_task(task219);

  vector<IndexRange> I352_index = {active_, active_, closed_, virt_};
  auto I352 = make_shared<Tensor>(I352_index);
  auto tensor220 = vector<shared_ptr<Tensor>>{s, I352};
  auto task220 = make_shared<Task220>(tensor220, pindex);
  task220->add_dep(task209);
  sourceq->add_task(task220);

  auto tensor221 = vector<shared_ptr<Tensor>>{I352, h1_, Gamma38_()};
  auto task221 = make_shared<Task221>(tensor221, pindex);
  task220->add_dep(task221);
  task221->add_dep(task209);
  sourceq->add_task(task221);

  auto tensor222 = vector<shared_ptr<Tensor>>{I352, v2_, Gamma35_()};
  auto task222 = make_shared<Task222>(tensor222, pindex);
  task220->add_dep(task222);
  task222->add_dep(task209);
  sourceq->add_task(task222);

  auto tensor223 = vector<shared_ptr<Tensor>>{I352, v2_, Gamma7_()};
  auto task223 = make_shared<Task223>(tensor223, pindex);
  task220->add_dep(task223);
  task223->add_dep(task209);
  sourceq->add_task(task223);

  auto tensor224 = vector<shared_ptr<Tensor>>{I352, v2_, Gamma35_()};
  auto task224 = make_shared<Task224>(tensor224, pindex);
  task220->add_dep(task224);
  task224->add_dep(task209);
  sourceq->add_task(task224);

  auto tensor225 = vector<shared_ptr<Tensor>>{I352, v2_, Gamma35_()};
  auto task225 = make_shared<Task225>(tensor225, pindex);
  task220->add_dep(task225);
  task225->add_dep(task209);
  sourceq->add_task(task225);

  vector<IndexRange> I354_index = {virt_, active_, active_, active_};
  auto I354 = make_shared<Tensor>(I354_index);
  auto tensor226 = vector<shared_ptr<Tensor>>{s, I354};
  auto task226 = make_shared<Task226>(tensor226, pindex);
  task226->add_dep(task209);
  sourceq->add_task(task226);

  auto tensor227 = vector<shared_ptr<Tensor>>{I354, Gamma60_(), h1_};
  auto task227 = make_shared<Task227>(tensor227, pindex);
  task226->add_dep(task227);
  task227->add_dep(task209);
  sourceq->add_task(task227);

  auto tensor228 = vector<shared_ptr<Tensor>>{I354, v2_, Gamma59_()};
  auto task228 = make_shared<Task228>(tensor228, pindex);
  task226->add_dep(task228);
  task228->add_dep(task209);
  sourceq->add_task(task228);

  auto tensor229 = vector<shared_ptr<Tensor>>{I354, Gamma57_(), v2_};
  auto task229 = make_shared<Task229>(tensor229, pindex);
  task226->add_dep(task229);
  task229->add_dep(task209);
  sourceq->add_task(task229);

  vector<IndexRange> I356_index = {closed_, closed_, active_, active_};
  auto I356 = make_shared<Tensor>(I356_index);
  auto tensor230 = vector<shared_ptr<Tensor>>{s, I356};
  auto task230 = make_shared<Task230>(tensor230, pindex);
  task230->add_dep(task209);
  sourceq->add_task(task230);

  auto tensor231 = vector<shared_ptr<Tensor>>{I356, Gamma94_(), v2_};
  auto task231 = make_shared<Task231>(tensor231, pindex);
  task230->add_dep(task231);
  task231->add_dep(task209);
  sourceq->add_task(task231);

  vector<IndexRange> I362_index = {active_, closed_, closed_, virt_};
  auto I362 = make_shared<Tensor>(I362_index);
  auto tensor232 = vector<shared_ptr<Tensor>>{s, I362};
  auto task232 = make_shared<Task232>(tensor232, pindex);
  task232->add_dep(task209);
  sourceq->add_task(task232);

  auto tensor233 = vector<shared_ptr<Tensor>>{I362, v2_, Gamma16_()};
  auto task233 = make_shared<Task233>(tensor233, pindex);
  task232->add_dep(task233);
  task233->add_dep(task209);
  sourceq->add_task(task233);

  auto tensor234 = vector<shared_ptr<Tensor>>{I362, Gamma16_(), v2_};
  auto task234 = make_shared<Task234>(tensor234, pindex);
  task232->add_dep(task234);
  task234->add_dep(task209);
  sourceq->add_task(task234);

  shared_ptr<Tensor> I386;
  if (diagonal) {
    vector<IndexRange> I386_index = {closed_, virt_, closed_, virt_};
    I386 = make_shared<Tensor>(I386_index);
  }
  shared_ptr<Task235> task235;
  if (diagonal) {
    auto tensor235 = vector<shared_ptr<Tensor>>{s, I386};
    task235 = make_shared<Task235>(tensor235, pindex);
    task235->add_dep(task209);
    sourceq->add_task(task235);
  }

  shared_ptr<Task236> task236;
  if (diagonal) {
    auto tensor236 = vector<shared_ptr<Tensor>>{I386, v2_};
    task236 = make_shared<Task236>(tensor236, pindex);
    task235->add_dep(task236);
    task236->add_dep(task209);
    sourceq->add_task(task236);
  }

  vector<IndexRange> I388_index = {active_, virt_, closed_, virt_};
  auto I388 = make_shared<Tensor>(I388_index);
  auto tensor237 = vector<shared_ptr<Tensor>>{s, I388};
  auto task237 = make_shared<Task237>(tensor237, pindex);
  task237->add_dep(task209);
  sourceq->add_task(task237);

  auto tensor238 = vector<shared_ptr<Tensor>>{I388, v2_, Gamma38_()};
  auto task238 = make_shared<Task238>(tensor238, pindex);
  task237->add_dep(task238);
  task238->add_dep(task209);
  sourceq->add_task(task238);

  auto tensor239 = vector<shared_ptr<Tensor>>{I388, v2_, Gamma38_()};
  auto task239 = make_shared<Task239>(tensor239, pindex);
  task237->add_dep(task239);
  task239->add_dep(task209);
  sourceq->add_task(task239);

  vector<IndexRange> I392_index = {active_, active_, virt_, virt_};
  auto I392 = make_shared<Tensor>(I392_index);
  auto tensor240 = vector<shared_ptr<Tensor>>{s, I392};
  auto task240 = make_shared<Task240>(tensor240, pindex);
  task240->add_dep(task209);
  sourceq->add_task(task240);

  auto tensor241 = vector<shared_ptr<Tensor>>{I392, v2_, Gamma60_()};
  auto task241 = make_shared<Task241>(tensor241, pindex);
  task240->add_dep(task241);
  task241->add_dep(task209);
  sourceq->add_task(task241);

  return sourceq;
}


#endif
