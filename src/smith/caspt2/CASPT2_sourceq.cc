//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_sourceqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_sourceq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  auto tensor229 = vector<shared_ptr<Tensor>>{s};
  auto task229 = make_shared<Task229>(tensor229, reset);
  sourceq->add_task(task229);

  vector<IndexRange> I288_index = {active_, active_, closed_, closed_};
  auto I288 = make_shared<Tensor>(I288_index);
  auto tensor230 = vector<shared_ptr<Tensor>>{s, I288};
  auto task230 = make_shared<Task230>(tensor230, pindex);
  task230->add_dep(task229);
  sourceq->add_task(task230);

  auto tensor231 = vector<shared_ptr<Tensor>>{I288, v2_, Gamma92_()};
  auto task231 = make_shared<Task231>(tensor231, pindex);
  task230->add_dep(task231);
  task231->add_dep(task229);
  sourceq->add_task(task231);

  vector<IndexRange> I290_index = {closed_, active_, active_, active_};
  auto I290 = make_shared<Tensor>(I290_index);
  auto tensor232 = vector<shared_ptr<Tensor>>{s, I290};
  auto task232 = make_shared<Task232>(tensor232, pindex);
  task232->add_dep(task229);
  sourceq->add_task(task232);

  auto tensor233 = vector<shared_ptr<Tensor>>{I290, Gamma105_(), v2_};
  auto task233 = make_shared<Task233>(tensor233, pindex);
  task232->add_dep(task233);
  task233->add_dep(task229);
  sourceq->add_task(task233);

  auto tensor234 = vector<shared_ptr<Tensor>>{I290, v2_, Gamma6_()};
  auto task234 = make_shared<Task234>(tensor234, pindex);
  task232->add_dep(task234);
  task234->add_dep(task229);
  sourceq->add_task(task234);

  auto tensor235 = vector<shared_ptr<Tensor>>{I290, Gamma7_(), h1_};
  auto task235 = make_shared<Task235>(tensor235, pindex);
  task232->add_dep(task235);
  task235->add_dep(task229);
  sourceq->add_task(task235);

  vector<IndexRange> I294_index = {active_, closed_, closed_, virt_};
  auto I294 = make_shared<Tensor>(I294_index);
  auto tensor236 = vector<shared_ptr<Tensor>>{s, I294};
  auto task236 = make_shared<Task236>(tensor236, pindex);
  task236->add_dep(task229);
  sourceq->add_task(task236);

  auto tensor237 = vector<shared_ptr<Tensor>>{I294, v2_, Gamma16_()};
  auto task237 = make_shared<Task237>(tensor237, pindex);
  task236->add_dep(task237);
  task237->add_dep(task229);
  sourceq->add_task(task237);

  auto tensor238 = vector<shared_ptr<Tensor>>{I294, v2_, Gamma16_()};
  auto task238 = make_shared<Task238>(tensor238, pindex);
  task236->add_dep(task238);
  task238->add_dep(task229);
  sourceq->add_task(task238);

  vector<IndexRange> I298_index = {active_, active_, closed_, virt_};
  auto I298 = make_shared<Tensor>(I298_index);
  auto tensor239 = vector<shared_ptr<Tensor>>{s, I298};
  auto task239 = make_shared<Task239>(tensor239, pindex);
  task239->add_dep(task229);
  sourceq->add_task(task239);

  auto tensor240 = vector<shared_ptr<Tensor>>{I298, v2_, Gamma35_()};
  auto task240 = make_shared<Task240>(tensor240, pindex);
  task239->add_dep(task240);
  task240->add_dep(task229);
  sourceq->add_task(task240);

  auto tensor241 = vector<shared_ptr<Tensor>>{I298, v2_, Gamma29_()};
  auto task241 = make_shared<Task241>(tensor241, pindex);
  task239->add_dep(task241);
  task241->add_dep(task229);
  sourceq->add_task(task241);

  auto tensor242 = vector<shared_ptr<Tensor>>{I298, Gamma32_(), v2_};
  auto task242 = make_shared<Task242>(tensor242, pindex);
  task239->add_dep(task242);
  task242->add_dep(task229);
  sourceq->add_task(task242);

  auto tensor243 = vector<shared_ptr<Tensor>>{I298, v2_, Gamma35_()};
  auto task243 = make_shared<Task243>(tensor243, pindex);
  task239->add_dep(task243);
  task243->add_dep(task229);
  sourceq->add_task(task243);

  auto tensor244 = vector<shared_ptr<Tensor>>{I298, h1_, Gamma38_()};
  auto task244 = make_shared<Task244>(tensor244, pindex);
  task239->add_dep(task244);
  task244->add_dep(task229);
  sourceq->add_task(task244);

  vector<IndexRange> I306_index = {closed_, virt_, active_, active_};
  auto I306 = make_shared<Tensor>(I306_index);
  auto tensor245 = vector<shared_ptr<Tensor>>{s, I306};
  auto task245 = make_shared<Task245>(tensor245, pindex);
  task245->add_dep(task229);
  sourceq->add_task(task245);

  vector<IndexRange> I307_index = {closed_, virt_, active_, active_};
  auto I307 = make_shared<Tensor>(I307_index);
  auto tensor246 = vector<shared_ptr<Tensor>>{I306, Gamma35_(), I307};
  auto task246 = make_shared<Task246>(tensor246, pindex);
  task245->add_dep(task246);
  task246->add_dep(task229);
  sourceq->add_task(task246);

  auto tensor247 = vector<shared_ptr<Tensor>>{I307, v2_};
  auto task247 = make_shared<Task247>(tensor247, pindex);
  task246->add_dep(task247);
  task247->add_dep(task229);
  sourceq->add_task(task247);

  auto tensor248 = vector<shared_ptr<Tensor>>{I306, v2_, Gamma7_()};
  auto task248 = make_shared<Task248>(tensor248, pindex);
  task245->add_dep(task248);
  task248->add_dep(task229);
  sourceq->add_task(task248);

  auto tensor249 = vector<shared_ptr<Tensor>>{I306, v2_, Gamma35_()};
  auto task249 = make_shared<Task249>(tensor249, pindex);
  task245->add_dep(task249);
  task249->add_dep(task229);
  sourceq->add_task(task249);

  auto tensor250 = vector<shared_ptr<Tensor>>{I306, Gamma38_(), h1_};
  auto task250 = make_shared<Task250>(tensor250, pindex);
  task245->add_dep(task250);
  task250->add_dep(task229);
  sourceq->add_task(task250);

  vector<IndexRange> I314_index = {virt_, active_, active_, active_};
  auto I314 = make_shared<Tensor>(I314_index);
  auto tensor251 = vector<shared_ptr<Tensor>>{s, I314};
  auto task251 = make_shared<Task251>(tensor251, pindex);
  task251->add_dep(task229);
  sourceq->add_task(task251);

  auto tensor252 = vector<shared_ptr<Tensor>>{I314, Gamma59_(), v2_};
  auto task252 = make_shared<Task252>(tensor252, pindex);
  task251->add_dep(task252);
  task252->add_dep(task229);
  sourceq->add_task(task252);

  auto tensor253 = vector<shared_ptr<Tensor>>{I314, Gamma57_(), v2_};
  auto task253 = make_shared<Task253>(tensor253, pindex);
  task251->add_dep(task253);
  task253->add_dep(task229);
  sourceq->add_task(task253);

  auto tensor254 = vector<shared_ptr<Tensor>>{I314, h1_, Gamma60_()};
  auto task254 = make_shared<Task254>(tensor254, pindex);
  task251->add_dep(task254);
  task254->add_dep(task229);
  sourceq->add_task(task254);

  shared_ptr<Tensor> I318;
  if (diagonal) {
    vector<IndexRange> I318_index = {closed_, virt_, closed_, virt_};
    I318 = make_shared<Tensor>(I318_index);
  }
  shared_ptr<Task255> task255;
  if (diagonal) {
    auto tensor255 = vector<shared_ptr<Tensor>>{s, I318};
    task255 = make_shared<Task255>(tensor255, pindex);
    task255->add_dep(task229);
    sourceq->add_task(task255);
  }

  shared_ptr<Task256> task256;
  if (diagonal) {
    auto tensor256 = vector<shared_ptr<Tensor>>{I318, v2_};
    task256 = make_shared<Task256>(tensor256, pindex);
    task255->add_dep(task256);
    task256->add_dep(task229);
    sourceq->add_task(task256);
  }

  vector<IndexRange> I320_index = {virt_, closed_, virt_, active_};
  auto I320 = make_shared<Tensor>(I320_index);
  auto tensor257 = vector<shared_ptr<Tensor>>{s, I320};
  auto task257 = make_shared<Task257>(tensor257, pindex);
  task257->add_dep(task229);
  sourceq->add_task(task257);

  vector<IndexRange> I321_index = {active_, virt_, closed_, virt_};
  auto I321 = make_shared<Tensor>(I321_index);
  auto tensor258 = vector<shared_ptr<Tensor>>{I320, Gamma38_(), I321};
  auto task258 = make_shared<Task258>(tensor258, pindex);
  task257->add_dep(task258);
  task258->add_dep(task229);
  sourceq->add_task(task258);

  auto tensor259 = vector<shared_ptr<Tensor>>{I321, v2_};
  auto task259 = make_shared<Task259>(tensor259, pindex);
  task258->add_dep(task259);
  task259->add_dep(task229);
  sourceq->add_task(task259);

  vector<IndexRange> I324_index = {active_, active_, virt_, virt_};
  auto I324 = make_shared<Tensor>(I324_index);
  auto tensor260 = vector<shared_ptr<Tensor>>{s, I324};
  auto task260 = make_shared<Task260>(tensor260, pindex);
  task260->add_dep(task229);
  sourceq->add_task(task260);

  auto tensor261 = vector<shared_ptr<Tensor>>{I324, v2_, Gamma60_()};
  auto task261 = make_shared<Task261>(tensor261, pindex);
  task260->add_dep(task261);
  task261->add_dep(task229);
  sourceq->add_task(task261);

  return sourceq;
}


#endif
