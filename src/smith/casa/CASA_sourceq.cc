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
  auto tensor223 = vector<shared_ptr<Tensor>>{s};
  auto task223 = make_shared<Task223>(tensor223, reset);
  sourceq->add_task(task223);

  vector<IndexRange> I304_index = {closed_, closed_, active_, active_};
  auto I304 = make_shared<Tensor>(I304_index);
  auto tensor224 = vector<shared_ptr<Tensor>>{s, I304};
  auto task224 = make_shared<Task224>(tensor224, pindex);
  task224->add_dep(task223);
  sourceq->add_task(task224);

  auto tensor225 = vector<shared_ptr<Tensor>>{I304, Gamma0_(), v2_};
  auto task225 = make_shared<Task225>(tensor225, pindex);
  task224->add_dep(task225);
  task225->add_dep(task223);
  sourceq->add_task(task225);

  vector<IndexRange> I306_index = {closed_, active_, active_, active_};
  auto I306 = make_shared<Tensor>(I306_index);
  auto tensor226 = vector<shared_ptr<Tensor>>{s, I306};
  auto task226 = make_shared<Task226>(tensor226, pindex);
  task226->add_dep(task223);
  sourceq->add_task(task226);

  auto tensor227 = vector<shared_ptr<Tensor>>{I306, Gamma109_(), v2_};
  auto task227 = make_shared<Task227>(tensor227, pindex);
  task226->add_dep(task227);
  task227->add_dep(task223);
  sourceq->add_task(task227);

  auto tensor228 = vector<shared_ptr<Tensor>>{I306, Gamma3_(), v2_};
  auto task228 = make_shared<Task228>(tensor228, pindex);
  task226->add_dep(task228);
  task228->add_dep(task223);
  sourceq->add_task(task228);

  auto tensor229 = vector<shared_ptr<Tensor>>{I306, Gamma5_(), h1_};
  auto task229 = make_shared<Task229>(tensor229, pindex);
  task226->add_dep(task229);
  task229->add_dep(task223);
  sourceq->add_task(task229);

  vector<IndexRange> I310_index = {closed_, closed_, virt_, active_};
  auto I310 = make_shared<Tensor>(I310_index);
  auto tensor230 = vector<shared_ptr<Tensor>>{s, I310};
  auto task230 = make_shared<Task230>(tensor230, pindex);
  task230->add_dep(task223);
  sourceq->add_task(task230);

  vector<IndexRange> I311_index = {closed_, active_, closed_, virt_};
  auto I311 = make_shared<Tensor>(I311_index);
  auto tensor231 = vector<shared_ptr<Tensor>>{I310, Gamma9_(), I311};
  auto task231 = make_shared<Task231>(tensor231, pindex);
  task230->add_dep(task231);
  task231->add_dep(task223);
  sourceq->add_task(task231);

  auto tensor232 = vector<shared_ptr<Tensor>>{I311, v2_};
  auto task232 = make_shared<Task232>(tensor232, pindex);
  task231->add_dep(task232);
  task232->add_dep(task223);
  sourceq->add_task(task232);

  vector<IndexRange> I314_index = {closed_, virt_, active_, active_};
  auto I314 = make_shared<Tensor>(I314_index);
  auto tensor233 = vector<shared_ptr<Tensor>>{s, I314};
  auto task233 = make_shared<Task233>(tensor233, pindex);
  task233->add_dep(task223);
  sourceq->add_task(task233);

  vector<IndexRange> I315_index = {closed_, virt_, active_, active_};
  auto I315 = make_shared<Tensor>(I315_index);
  auto tensor234 = vector<shared_ptr<Tensor>>{I314, Gamma25_(), I315};
  auto task234 = make_shared<Task234>(tensor234, pindex);
  task233->add_dep(task234);
  task234->add_dep(task223);
  sourceq->add_task(task234);

  auto tensor235 = vector<shared_ptr<Tensor>>{I315, v2_};
  auto task235 = make_shared<Task235>(tensor235, pindex);
  task234->add_dep(task235);
  task235->add_dep(task223);
  sourceq->add_task(task235);

  auto tensor236 = vector<shared_ptr<Tensor>>{I314, Gamma26_(), v2_};
  auto task236 = make_shared<Task236>(tensor236, pindex);
  task233->add_dep(task236);
  task236->add_dep(task223);
  sourceq->add_task(task236);

  auto tensor237 = vector<shared_ptr<Tensor>>{I314, Gamma24_(), v2_};
  auto task237 = make_shared<Task237>(tensor237, pindex);
  task233->add_dep(task237);
  task237->add_dep(task223);
  sourceq->add_task(task237);

  auto tensor238 = vector<shared_ptr<Tensor>>{I314, Gamma29_(), h1_};
  auto task238 = make_shared<Task238>(tensor238, pindex);
  task233->add_dep(task238);
  task238->add_dep(task223);
  sourceq->add_task(task238);

  vector<IndexRange> I322_index = {closed_, virt_, active_, active_};
  auto I322 = make_shared<Tensor>(I322_index);
  auto tensor239 = vector<shared_ptr<Tensor>>{s, I322};
  auto task239 = make_shared<Task239>(tensor239, pindex);
  task239->add_dep(task223);
  sourceq->add_task(task239);

  vector<IndexRange> I323_index = {closed_, virt_, active_, active_};
  auto I323 = make_shared<Tensor>(I323_index);
  auto tensor240 = vector<shared_ptr<Tensor>>{I322, Gamma25_(), I323};
  auto task240 = make_shared<Task240>(tensor240, pindex);
  task239->add_dep(task240);
  task240->add_dep(task223);
  sourceq->add_task(task240);

  auto tensor241 = vector<shared_ptr<Tensor>>{I323, v2_};
  auto task241 = make_shared<Task241>(tensor241, pindex);
  task240->add_dep(task241);
  task241->add_dep(task223);
  sourceq->add_task(task241);

  auto tensor242 = vector<shared_ptr<Tensor>>{I322, Gamma5_(), v2_};
  auto task242 = make_shared<Task242>(tensor242, pindex);
  task239->add_dep(task242);
  task242->add_dep(task223);
  sourceq->add_task(task242);

  auto tensor243 = vector<shared_ptr<Tensor>>{I322, Gamma29_(), h1_};
  auto task243 = make_shared<Task243>(tensor243, pindex);
  task239->add_dep(task243);
  task243->add_dep(task223);
  sourceq->add_task(task243);

  vector<IndexRange> I330_index = {virt_, active_, active_, active_};
  auto I330 = make_shared<Tensor>(I330_index);
  auto tensor244 = vector<shared_ptr<Tensor>>{s, I330};
  auto task244 = make_shared<Task244>(tensor244, pindex);
  task244->add_dep(task223);
  sourceq->add_task(task244);

  auto tensor245 = vector<shared_ptr<Tensor>>{I330, Gamma52_(), v2_};
  auto task245 = make_shared<Task245>(tensor245, pindex);
  task244->add_dep(task245);
  task245->add_dep(task223);
  sourceq->add_task(task245);

  auto tensor246 = vector<shared_ptr<Tensor>>{I330, Gamma49_(), v2_};
  auto task246 = make_shared<Task246>(tensor246, pindex);
  task244->add_dep(task246);
  task246->add_dep(task223);
  sourceq->add_task(task246);

  auto tensor247 = vector<shared_ptr<Tensor>>{I330, Gamma50_(), h1_};
  auto task247 = make_shared<Task247>(tensor247, pindex);
  task244->add_dep(task247);
  task247->add_dep(task223);
  sourceq->add_task(task247);

  shared_ptr<Tensor> I334;
  if (diagonal) {
    vector<IndexRange> I334_index = {closed_, virt_, closed_, virt_};
    I334 = make_shared<Tensor>(I334_index);
  }
  shared_ptr<Task248> task248;
  if (diagonal) {
    auto tensor248 = vector<shared_ptr<Tensor>>{s, I334};
    task248 = make_shared<Task248>(tensor248, pindex);
    task248->add_dep(task223);
    sourceq->add_task(task248);
  }

  shared_ptr<Task249> task249;
  if (diagonal) {
    auto tensor249 = vector<shared_ptr<Tensor>>{I334, v2_};
    task249 = make_shared<Task249>(tensor249, pindex);
    task248->add_dep(task249);
    task249->add_dep(task223);
    sourceq->add_task(task249);
  }

  vector<IndexRange> I336_index = {virt_, closed_, virt_, active_};
  auto I336 = make_shared<Tensor>(I336_index);
  auto tensor250 = vector<shared_ptr<Tensor>>{s, I336};
  auto task250 = make_shared<Task250>(tensor250, pindex);
  task250->add_dep(task223);
  sourceq->add_task(task250);

  vector<IndexRange> I337_index = {active_, virt_, closed_, virt_};
  auto I337 = make_shared<Tensor>(I337_index);
  auto tensor251 = vector<shared_ptr<Tensor>>{I336, Gamma29_(), I337};
  auto task251 = make_shared<Task251>(tensor251, pindex);
  task250->add_dep(task251);
  task251->add_dep(task223);
  sourceq->add_task(task251);

  auto tensor252 = vector<shared_ptr<Tensor>>{I337, v2_};
  auto task252 = make_shared<Task252>(tensor252, pindex);
  task251->add_dep(task252);
  task252->add_dep(task223);
  sourceq->add_task(task252);

  vector<IndexRange> I340_index = {virt_, virt_, active_, active_};
  auto I340 = make_shared<Tensor>(I340_index);
  auto tensor253 = vector<shared_ptr<Tensor>>{s, I340};
  auto task253 = make_shared<Task253>(tensor253, pindex);
  task253->add_dep(task223);
  sourceq->add_task(task253);

  auto tensor254 = vector<shared_ptr<Tensor>>{I340, Gamma50_(), v2_};
  auto task254 = make_shared<Task254>(tensor254, pindex);
  task253->add_dep(task254);
  task254->add_dep(task223);
  sourceq->add_task(task254);

  return sourceq;
}


#endif
