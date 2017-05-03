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
  auto tensor242 = vector<shared_ptr<Tensor>>{s};
  auto task242 = make_shared<Task242>(tensor242, reset);
  sourceq->add_task(task242);

  vector<IndexRange> I288_index = {closed_, closed_, active_, active_};
  auto I288 = make_shared<Tensor>(I288_index);
  auto tensor243 = vector<shared_ptr<Tensor>>{s, I288};
  auto task243 = make_shared<Task243>(tensor243, pindex);
  task243->add_dep(task242);
  sourceq->add_task(task243);

  auto tensor244 = vector<shared_ptr<Tensor>>{I288, Gamma92_(), v2_};
  auto task244 = make_shared<Task244>(tensor244, pindex);
  task243->add_dep(task244);
  task244->add_dep(task242);
  sourceq->add_task(task244);

  vector<IndexRange> I290_index = {closed_, active_, active_, active_};
  auto I290 = make_shared<Tensor>(I290_index);
  auto tensor245 = vector<shared_ptr<Tensor>>{s, I290};
  auto task245 = make_shared<Task245>(tensor245, pindex);
  task245->add_dep(task242);
  sourceq->add_task(task245);

  auto tensor246 = vector<shared_ptr<Tensor>>{I290, Gamma105_(), v2_};
  auto task246 = make_shared<Task246>(tensor246, pindex);
  task245->add_dep(task246);
  task246->add_dep(task242);
  sourceq->add_task(task246);

  auto tensor247 = vector<shared_ptr<Tensor>>{I290, Gamma6_(), v2_};
  auto task247 = make_shared<Task247>(tensor247, pindex);
  task245->add_dep(task247);
  task247->add_dep(task242);
  sourceq->add_task(task247);

  auto tensor248 = vector<shared_ptr<Tensor>>{I290, Gamma7_(), h1_};
  auto task248 = make_shared<Task248>(tensor248, pindex);
  task245->add_dep(task248);
  task248->add_dep(task242);
  sourceq->add_task(task248);

  vector<IndexRange> I294_index = {closed_, closed_, virt_, active_};
  auto I294 = make_shared<Tensor>(I294_index);
  auto tensor249 = vector<shared_ptr<Tensor>>{s, I294};
  auto task249 = make_shared<Task249>(tensor249, pindex);
  task249->add_dep(task242);
  sourceq->add_task(task249);

  vector<IndexRange> I295_index = {closed_, active_, closed_, virt_};
  auto I295 = make_shared<Tensor>(I295_index);
  auto tensor250 = vector<shared_ptr<Tensor>>{I294, Gamma16_(), I295};
  auto task250 = make_shared<Task250>(tensor250, pindex);
  task249->add_dep(task250);
  task250->add_dep(task242);
  sourceq->add_task(task250);

  auto tensor251 = vector<shared_ptr<Tensor>>{I295, v2_};
  auto task251 = make_shared<Task251>(tensor251, pindex);
  task250->add_dep(task251);
  task251->add_dep(task242);
  sourceq->add_task(task251);

  vector<IndexRange> I298_index = {closed_, virt_, active_, active_};
  auto I298 = make_shared<Tensor>(I298_index);
  auto tensor252 = vector<shared_ptr<Tensor>>{s, I298};
  auto task252 = make_shared<Task252>(tensor252, pindex);
  task252->add_dep(task242);
  sourceq->add_task(task252);

  vector<IndexRange> I299_index = {closed_, virt_, active_, active_};
  auto I299 = make_shared<Tensor>(I299_index);
  auto tensor253 = vector<shared_ptr<Tensor>>{I298, Gamma35_(), I299};
  auto task253 = make_shared<Task253>(tensor253, pindex);
  task252->add_dep(task253);
  task253->add_dep(task242);
  sourceq->add_task(task253);

  auto tensor254 = vector<shared_ptr<Tensor>>{I299, v2_};
  auto task254 = make_shared<Task254>(tensor254, pindex);
  task253->add_dep(task254);
  task254->add_dep(task242);
  sourceq->add_task(task254);

  auto tensor255 = vector<shared_ptr<Tensor>>{I298, Gamma29_(), v2_};
  auto task255 = make_shared<Task255>(tensor255, pindex);
  task252->add_dep(task255);
  task255->add_dep(task242);
  sourceq->add_task(task255);

  auto tensor256 = vector<shared_ptr<Tensor>>{I298, Gamma32_(), v2_};
  auto task256 = make_shared<Task256>(tensor256, pindex);
  task252->add_dep(task256);
  task256->add_dep(task242);
  sourceq->add_task(task256);

  auto tensor257 = vector<shared_ptr<Tensor>>{I298, Gamma38_(), h1_};
  auto task257 = make_shared<Task257>(tensor257, pindex);
  task252->add_dep(task257);
  task257->add_dep(task242);
  sourceq->add_task(task257);

  vector<IndexRange> I306_index = {closed_, virt_, active_, active_};
  auto I306 = make_shared<Tensor>(I306_index);
  auto tensor258 = vector<shared_ptr<Tensor>>{s, I306};
  auto task258 = make_shared<Task258>(tensor258, pindex);
  task258->add_dep(task242);
  sourceq->add_task(task258);

  vector<IndexRange> I307_index = {closed_, virt_, active_, active_};
  auto I307 = make_shared<Tensor>(I307_index);
  auto tensor259 = vector<shared_ptr<Tensor>>{I306, Gamma35_(), I307};
  auto task259 = make_shared<Task259>(tensor259, pindex);
  task258->add_dep(task259);
  task259->add_dep(task242);
  sourceq->add_task(task259);

  auto tensor260 = vector<shared_ptr<Tensor>>{I307, v2_};
  auto task260 = make_shared<Task260>(tensor260, pindex);
  task259->add_dep(task260);
  task260->add_dep(task242);
  sourceq->add_task(task260);

  auto tensor261 = vector<shared_ptr<Tensor>>{I306, Gamma7_(), v2_};
  auto task261 = make_shared<Task261>(tensor261, pindex);
  task258->add_dep(task261);
  task261->add_dep(task242);
  sourceq->add_task(task261);

  auto tensor262 = vector<shared_ptr<Tensor>>{I306, Gamma38_(), h1_};
  auto task262 = make_shared<Task262>(tensor262, pindex);
  task258->add_dep(task262);
  task262->add_dep(task242);
  sourceq->add_task(task262);

  vector<IndexRange> I314_index = {virt_, active_, active_, active_};
  auto I314 = make_shared<Tensor>(I314_index);
  auto tensor263 = vector<shared_ptr<Tensor>>{s, I314};
  auto task263 = make_shared<Task263>(tensor263, pindex);
  task263->add_dep(task242);
  sourceq->add_task(task263);

  auto tensor264 = vector<shared_ptr<Tensor>>{I314, Gamma59_(), v2_};
  auto task264 = make_shared<Task264>(tensor264, pindex);
  task263->add_dep(task264);
  task264->add_dep(task242);
  sourceq->add_task(task264);

  auto tensor265 = vector<shared_ptr<Tensor>>{I314, Gamma57_(), v2_};
  auto task265 = make_shared<Task265>(tensor265, pindex);
  task263->add_dep(task265);
  task265->add_dep(task242);
  sourceq->add_task(task265);

  auto tensor266 = vector<shared_ptr<Tensor>>{I314, Gamma60_(), h1_};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task263->add_dep(task266);
  task266->add_dep(task242);
  sourceq->add_task(task266);

  shared_ptr<Tensor> I318;
  if (diagonal) {
    vector<IndexRange> I318_index = {closed_, virt_, closed_, virt_};
    I318 = make_shared<Tensor>(I318_index);
  }
  shared_ptr<Task267> task267;
  if (diagonal) {
    auto tensor267 = vector<shared_ptr<Tensor>>{s, I318};
    task267 = make_shared<Task267>(tensor267, pindex);
    task267->add_dep(task242);
    sourceq->add_task(task267);
  }

  shared_ptr<Task268> task268;
  if (diagonal) {
    auto tensor268 = vector<shared_ptr<Tensor>>{I318, v2_};
    task268 = make_shared<Task268>(tensor268, pindex);
    task267->add_dep(task268);
    task268->add_dep(task242);
    sourceq->add_task(task268);
  }

  vector<IndexRange> I320_index = {virt_, closed_, virt_, active_};
  auto I320 = make_shared<Tensor>(I320_index);
  auto tensor269 = vector<shared_ptr<Tensor>>{s, I320};
  auto task269 = make_shared<Task269>(tensor269, pindex);
  task269->add_dep(task242);
  sourceq->add_task(task269);

  vector<IndexRange> I321_index = {active_, virt_, closed_, virt_};
  auto I321 = make_shared<Tensor>(I321_index);
  auto tensor270 = vector<shared_ptr<Tensor>>{I320, Gamma38_(), I321};
  auto task270 = make_shared<Task270>(tensor270, pindex);
  task269->add_dep(task270);
  task270->add_dep(task242);
  sourceq->add_task(task270);

  auto tensor271 = vector<shared_ptr<Tensor>>{I321, v2_};
  auto task271 = make_shared<Task271>(tensor271, pindex);
  task270->add_dep(task271);
  task271->add_dep(task242);
  sourceq->add_task(task271);

  vector<IndexRange> I324_index = {virt_, virt_, active_, active_};
  auto I324 = make_shared<Tensor>(I324_index);
  auto tensor272 = vector<shared_ptr<Tensor>>{s, I324};
  auto task272 = make_shared<Task272>(tensor272, pindex);
  task272->add_dep(task242);
  sourceq->add_task(task272);

  auto tensor273 = vector<shared_ptr<Tensor>>{I324, Gamma60_(), v2_};
  auto task273 = make_shared<Task273>(tensor273, pindex);
  task272->add_dep(task273);
  task273->add_dep(task242);
  sourceq->add_task(task273);

  return sourceq;
}


#endif
