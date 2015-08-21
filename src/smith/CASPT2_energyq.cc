//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_energyqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_energyq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto energyq = make_shared<Queue>();
  vector<IndexRange> I334_index;
  auto I334 = make_shared<Tensor>(I334_index);
  vector<IndexRange> I335_index = {active_, active_, active_, active_};
  auto I335 = make_shared<Tensor>(I335_index);
  auto tensor241 = vector<shared_ptr<Tensor>>{I334, Gamma92_(), I335};
  auto task241 = make_shared<Task241>(tensor241, pindex);
  energyq->add_task(task241);

  auto tensor242 = vector<shared_ptr<Tensor>>{I335, v2_, t2};
  auto task242 = make_shared<Task242>(tensor242, pindex);
  task241->add_dep(task242);
  energyq->add_task(task242);

  vector<IndexRange> I338_index = {closed_, active_, active_, active_};
  auto I338 = make_shared<Tensor>(I338_index);
  auto tensor243 = vector<shared_ptr<Tensor>>{I334, t2, I338};
  auto task243 = make_shared<Task243>(tensor243, pindex);
  task241->add_dep(task243);
  energyq->add_task(task243);

  auto tensor244 = vector<shared_ptr<Tensor>>{I338, Gamma105_(), v2_};
  auto task244 = make_shared<Task244>(tensor244, pindex);
  task243->add_dep(task244);
  energyq->add_task(task244);

  auto tensor245 = vector<shared_ptr<Tensor>>{I338, Gamma6_(), v2_};
  auto task245 = make_shared<Task245>(tensor245, pindex);
  task243->add_dep(task245);
  energyq->add_task(task245);

  vector<IndexRange> I344_index = {active_, active_};
  auto I344 = make_shared<Tensor>(I344_index);
  auto tensor246 = vector<shared_ptr<Tensor>>{I334, Gamma16_(), I344};
  auto task246 = make_shared<Task246>(tensor246, pindex);
  task241->add_dep(task246);
  energyq->add_task(task246);

  vector<IndexRange> I345_index = {closed_, active_, closed_, virt_};
  auto I345 = make_shared<Tensor>(I345_index);
  auto tensor247 = vector<shared_ptr<Tensor>>{I344, t2, I345};
  auto task247 = make_shared<Task247>(tensor247, pindex);
  task246->add_dep(task247);
  energyq->add_task(task247);

  auto tensor248 = vector<shared_ptr<Tensor>>{I345, v2_};
  auto task248 = make_shared<Task248>(tensor248, pindex);
  task247->add_dep(task248);
  energyq->add_task(task248);

  vector<IndexRange> I350_index = {active_, active_, active_, active_};
  auto I350 = make_shared<Tensor>(I350_index);
  auto tensor249 = vector<shared_ptr<Tensor>>{I334, Gamma35_(), I350};
  auto task249 = make_shared<Task249>(tensor249, pindex);
  task241->add_dep(task249);
  energyq->add_task(task249);

  auto tensor250 = vector<shared_ptr<Tensor>>{I350, v2_, t2};
  auto task250 = make_shared<Task250>(tensor250, pindex);
  task249->add_dep(task250);
  energyq->add_task(task250);

  auto tensor251 = vector<shared_ptr<Tensor>>{I350, v2_, t2};
  auto task251 = make_shared<Task251>(tensor251, pindex);
  task249->add_dep(task251);
  energyq->add_task(task251);

  vector<IndexRange> I363_index = {closed_, virt_, active_, active_};
  auto I363 = make_shared<Tensor>(I363_index);
  auto tensor252 = vector<shared_ptr<Tensor>>{I350, t2, I363};
  auto task252 = make_shared<Task252>(tensor252, pindex);
  task249->add_dep(task252);
  energyq->add_task(task252);

  auto tensor253 = vector<shared_ptr<Tensor>>{I363, v2_};
  auto task253 = make_shared<Task253>(tensor253, pindex);
  task252->add_dep(task253);
  energyq->add_task(task253);

  auto tensor254 = vector<shared_ptr<Tensor>>{I350, v2_, t2};
  auto task254 = make_shared<Task254>(tensor254, pindex);
  task249->add_dep(task254);
  energyq->add_task(task254);

  vector<IndexRange> I353_index = {active_, active_, active_, active_};
  auto I353 = make_shared<Tensor>(I353_index);
  auto tensor255 = vector<shared_ptr<Tensor>>{I334, Gamma29_(), I353};
  auto task255 = make_shared<Task255>(tensor255, pindex);
  task241->add_dep(task255);
  energyq->add_task(task255);

  auto tensor256 = vector<shared_ptr<Tensor>>{I353, v2_, t2};
  auto task256 = make_shared<Task256>(tensor256, pindex);
  task255->add_dep(task256);
  energyq->add_task(task256);

  vector<IndexRange> I356_index = {active_, active_, active_, active_};
  auto I356 = make_shared<Tensor>(I356_index);
  auto tensor257 = vector<shared_ptr<Tensor>>{I334, Gamma32_(), I356};
  auto task257 = make_shared<Task257>(tensor257, pindex);
  task241->add_dep(task257);
  energyq->add_task(task257);

  auto tensor258 = vector<shared_ptr<Tensor>>{I356, v2_, t2};
  auto task258 = make_shared<Task258>(tensor258, pindex);
  task257->add_dep(task258);
  energyq->add_task(task258);

  vector<IndexRange> I365_index = {active_, active_, active_, active_};
  auto I365 = make_shared<Tensor>(I365_index);
  auto tensor259 = vector<shared_ptr<Tensor>>{I334, Gamma7_(), I365};
  auto task259 = make_shared<Task259>(tensor259, pindex);
  task241->add_dep(task259);
  energyq->add_task(task259);

  auto tensor260 = vector<shared_ptr<Tensor>>{I365, v2_, t2};
  auto task260 = make_shared<Task260>(tensor260, pindex);
  task259->add_dep(task260);
  energyq->add_task(task260);

  vector<IndexRange> I374_index = {active_, active_, active_, active_, active_, active_};
  auto I374 = make_shared<Tensor>(I374_index);
  auto tensor261 = vector<shared_ptr<Tensor>>{I334, Gamma59_(), I374};
  auto task261 = make_shared<Task261>(tensor261, pindex);
  task241->add_dep(task261);
  energyq->add_task(task261);

  auto tensor262 = vector<shared_ptr<Tensor>>{I374, v2_, t2};
  auto task262 = make_shared<Task262>(tensor262, pindex);
  task261->add_dep(task262);
  energyq->add_task(task262);

  vector<IndexRange> I377_index = {active_, active_, active_, active_, active_, active_};
  auto I377 = make_shared<Tensor>(I377_index);
  auto tensor263 = vector<shared_ptr<Tensor>>{I334, Gamma57_(), I377};
  auto task263 = make_shared<Task263>(tensor263, pindex);
  task241->add_dep(task263);
  energyq->add_task(task263);

  auto tensor264 = vector<shared_ptr<Tensor>>{I377, v2_, t2};
  auto task264 = make_shared<Task264>(tensor264, pindex);
  task263->add_dep(task264);
  energyq->add_task(task264);

  auto tensor265 = vector<shared_ptr<Tensor>>{I334, t2, v2_};
  auto task265 = make_shared<Task265>(tensor265, pindex);
  task241->add_dep(task265);
  energyq->add_task(task265);

  auto tensor266 = vector<shared_ptr<Tensor>>{I334, v2_, t2};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task241->add_dep(task266);
  energyq->add_task(task266);

  vector<IndexRange> I384_index = {active_, active_};
  auto I384 = make_shared<Tensor>(I384_index);
  auto tensor267 = vector<shared_ptr<Tensor>>{I334, Gamma38_(), I384};
  auto task267 = make_shared<Task267>(tensor267, pindex);
  task241->add_dep(task267);
  energyq->add_task(task267);

  vector<IndexRange> I385_index = {active_, virt_, closed_, virt_};
  auto I385 = make_shared<Tensor>(I385_index);
  auto tensor268 = vector<shared_ptr<Tensor>>{I384, t2, I385};
  auto task268 = make_shared<Task268>(tensor268, pindex);
  task267->add_dep(task268);
  energyq->add_task(task268);

  auto tensor269 = vector<shared_ptr<Tensor>>{I385, v2_};
  auto task269 = make_shared<Task269>(tensor269, pindex);
  task268->add_dep(task269);
  energyq->add_task(task269);

  auto tensor270 = vector<shared_ptr<Tensor>>{I384, h1_, t2};
  auto task270 = make_shared<Task270>(tensor270, pindex);
  task267->add_dep(task270);
  energyq->add_task(task270);

  auto tensor271 = vector<shared_ptr<Tensor>>{I384, t2, h1_};
  auto task271 = make_shared<Task271>(tensor271, pindex);
  task267->add_dep(task271);
  energyq->add_task(task271);

  vector<IndexRange> I390_index = {active_, active_, active_, active_};
  auto I390 = make_shared<Tensor>(I390_index);
  auto tensor272 = vector<shared_ptr<Tensor>>{I334, Gamma60_(), I390};
  auto task272 = make_shared<Task272>(tensor272, pindex);
  task241->add_dep(task272);
  energyq->add_task(task272);

  auto tensor273 = vector<shared_ptr<Tensor>>{I390, v2_, t2};
  auto task273 = make_shared<Task273>(tensor273, pindex);
  task272->add_dep(task273);
  energyq->add_task(task273);

  auto tensor274 = vector<shared_ptr<Tensor>>{I390, t2, h1_};
  auto task274 = make_shared<Task274>(tensor274, pindex);
  task272->add_dep(task274);
  energyq->add_task(task274);

  vector<IndexRange> I393_index = {closed_, active_};
  auto I393 = make_shared<Tensor>(I393_index);
  auto tensor275 = vector<shared_ptr<Tensor>>{I334, h1_, I393};
  auto task275 = make_shared<Task275>(tensor275, pindex);
  task241->add_dep(task275);
  energyq->add_task(task275);

  auto tensor276 = vector<shared_ptr<Tensor>>{I393, Gamma7_(), t2};
  auto task276 = make_shared<Task276>(tensor276, pindex);
  task275->add_dep(task276);
  energyq->add_task(task276);

  return energyq;
}


#endif
