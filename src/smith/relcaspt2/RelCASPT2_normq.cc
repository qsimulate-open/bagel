//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASPT2_normqq.cc
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


#include <src/smith/relcaspt2/RelCASPT2.h>
#include <src/smith/relcaspt2/RelCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelCASPT2::RelCASPT2::make_normq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto normq = make_shared<Queue>();
  auto tensor235 = vector<shared_ptr<Tensor>>{n};
  auto task235 = make_shared<Task235>(tensor235, reset);
  normq->add_task(task235);

  vector<IndexRange> I334_index = {closed_, closed_, active_, active_};
  auto I334 = make_shared<Tensor>(I334_index);
  auto tensor236 = vector<shared_ptr<Tensor>>{n, I334};
  auto task236 = make_shared<Task236>(tensor236, pindex);
  task236->add_dep(task235);
  normq->add_task(task236);

  auto tensor237 = vector<shared_ptr<Tensor>>{I334, Gamma92_(), t2};
  auto task237 = make_shared<Task237>(tensor237, pindex);
  task236->add_dep(task237);
  task237->add_dep(task235);
  normq->add_task(task237);

  vector<IndexRange> I336_index = {closed_, active_, active_, active_};
  auto I336 = make_shared<Tensor>(I336_index);
  auto tensor238 = vector<shared_ptr<Tensor>>{n, I336};
  auto task238 = make_shared<Task238>(tensor238, pindex);
  task238->add_dep(task235);
  normq->add_task(task238);

  auto tensor239 = vector<shared_ptr<Tensor>>{I336, Gamma6_(), t2};
  auto task239 = make_shared<Task239>(tensor239, pindex);
  task238->add_dep(task239);
  task239->add_dep(task235);
  normq->add_task(task239);

  vector<IndexRange> I338_index = {closed_, virt_, closed_, active_};
  auto I338 = make_shared<Tensor>(I338_index);
  auto tensor240 = vector<shared_ptr<Tensor>>{n, I338};
  auto task240 = make_shared<Task240>(tensor240, pindex);
  task240->add_dep(task235);
  normq->add_task(task240);

  vector<IndexRange> I339_index = {closed_, virt_, closed_, active_};
  auto I339 = make_shared<Tensor>(I339_index);
  auto tensor241 = vector<shared_ptr<Tensor>>{I338, Gamma16_(), I339};
  auto task241 = make_shared<Task241>(tensor241, pindex);
  task240->add_dep(task241);
  task241->add_dep(task235);
  normq->add_task(task241);

  auto tensor242 = vector<shared_ptr<Tensor>>{I339, t2};
  auto task242 = make_shared<Task242>(tensor242, pindex);
  task241->add_dep(task242);
  task242->add_dep(task235);
  normq->add_task(task242);

  vector<IndexRange> I342_index = {virt_, closed_, active_, active_};
  auto I342 = make_shared<Tensor>(I342_index);
  auto tensor243 = vector<shared_ptr<Tensor>>{n, I342};
  auto task243 = make_shared<Task243>(tensor243, pindex);
  task243->add_dep(task235);
  normq->add_task(task243);

  auto tensor244 = vector<shared_ptr<Tensor>>{I342, Gamma32_(), t2};
  auto task244 = make_shared<Task244>(tensor244, pindex);
  task243->add_dep(task244);
  task244->add_dep(task235);
  normq->add_task(task244);

  auto tensor245 = vector<shared_ptr<Tensor>>{I342, Gamma35_(), t2};
  auto task245 = make_shared<Task245>(tensor245, pindex);
  task243->add_dep(task245);
  task245->add_dep(task235);
  normq->add_task(task245);

  vector<IndexRange> I346_index = {virt_, closed_, active_, active_};
  auto I346 = make_shared<Tensor>(I346_index);
  auto tensor246 = vector<shared_ptr<Tensor>>{n, I346};
  auto task246 = make_shared<Task246>(tensor246, pindex);
  task246->add_dep(task235);
  normq->add_task(task246);

  vector<IndexRange> I347_index = {active_, virt_, closed_, active_};
  auto I347 = make_shared<Tensor>(I347_index);
  auto tensor247 = vector<shared_ptr<Tensor>>{I346, Gamma35_(), I347};
  auto task247 = make_shared<Task247>(tensor247, pindex);
  task246->add_dep(task247);
  task247->add_dep(task235);
  normq->add_task(task247);

  auto tensor248 = vector<shared_ptr<Tensor>>{I347, t2};
  auto task248 = make_shared<Task248>(tensor248, pindex);
  task247->add_dep(task248);
  task248->add_dep(task235);
  normq->add_task(task248);

  vector<IndexRange> I350_index = {virt_, active_, active_, active_};
  auto I350 = make_shared<Tensor>(I350_index);
  auto tensor249 = vector<shared_ptr<Tensor>>{n, I350};
  auto task249 = make_shared<Task249>(tensor249, pindex);
  task249->add_dep(task235);
  normq->add_task(task249);

  auto tensor250 = vector<shared_ptr<Tensor>>{I350, Gamma59_(), t2};
  auto task250 = make_shared<Task250>(tensor250, pindex);
  task249->add_dep(task250);
  task250->add_dep(task235);
  normq->add_task(task250);

  shared_ptr<Tensor> I352;
  if (diagonal) {
    vector<IndexRange> I352_index = {closed_, virt_, closed_, virt_};
    I352 = make_shared<Tensor>(I352_index);
  }
  shared_ptr<Task251> task251;
  if (diagonal) {
    auto tensor251 = vector<shared_ptr<Tensor>>{n, I352};
    task251 = make_shared<Task251>(tensor251, pindex);
    task251->add_dep(task235);
    normq->add_task(task251);
  }

  shared_ptr<Task252> task252;
  if (diagonal) {
    auto tensor252 = vector<shared_ptr<Tensor>>{I352, t2};
    task252 = make_shared<Task252>(tensor252, pindex);
    task251->add_dep(task252);
    task252->add_dep(task235);
    normq->add_task(task252);
  }

  vector<IndexRange> I354_index = {virt_, closed_, virt_, active_};
  auto I354 = make_shared<Tensor>(I354_index);
  auto tensor253 = vector<shared_ptr<Tensor>>{n, I354};
  auto task253 = make_shared<Task253>(tensor253, pindex);
  task253->add_dep(task235);
  normq->add_task(task253);

  vector<IndexRange> I355_index = {active_, virt_, closed_, virt_};
  auto I355 = make_shared<Tensor>(I355_index);
  auto tensor254 = vector<shared_ptr<Tensor>>{I354, Gamma38_(), I355};
  auto task254 = make_shared<Task254>(tensor254, pindex);
  task253->add_dep(task254);
  task254->add_dep(task235);
  normq->add_task(task254);

  auto tensor255 = vector<shared_ptr<Tensor>>{I355, t2};
  auto task255 = make_shared<Task255>(tensor255, pindex);
  task254->add_dep(task255);
  task255->add_dep(task235);
  normq->add_task(task255);

  vector<IndexRange> I358_index = {virt_, virt_, active_, active_};
  auto I358 = make_shared<Tensor>(I358_index);
  auto tensor256 = vector<shared_ptr<Tensor>>{n, I358};
  auto task256 = make_shared<Task256>(tensor256, pindex);
  task256->add_dep(task235);
  normq->add_task(task256);

  auto tensor257 = vector<shared_ptr<Tensor>>{I358, Gamma60_(), t2};
  auto task257 = make_shared<Task257>(tensor257, pindex);
  task256->add_dep(task257);
  task257->add_dep(task235);
  normq->add_task(task257);

  return normq;
}


#endif
