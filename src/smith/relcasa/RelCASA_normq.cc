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

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelCASA::RelCASA::make_normq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto normq = make_shared<Queue>();
  auto tensor234 = vector<shared_ptr<Tensor>>{n};
  auto task234 = make_shared<Task234>(tensor234, reset);
  normq->add_task(task234);

  vector<IndexRange> I310_index = {closed_, closed_, active_, active_};
  auto I310 = make_shared<Tensor>(I310_index);
  auto tensor235 = vector<shared_ptr<Tensor>>{n, I310};
  auto task235 = make_shared<Task235>(tensor235, pindex);
  task235->add_dep(task234);
  normq->add_task(task235);

  auto tensor236 = vector<shared_ptr<Tensor>>{I310, Gamma96_(), t2};
  auto task236 = make_shared<Task236>(tensor236, pindex);
  task235->add_dep(task236);
  task236->add_dep(task234);
  normq->add_task(task236);

  vector<IndexRange> I312_index = {closed_, active_, active_, active_};
  auto I312 = make_shared<Tensor>(I312_index);
  auto tensor237 = vector<shared_ptr<Tensor>>{n, I312};
  auto task237 = make_shared<Task237>(tensor237, pindex);
  task237->add_dep(task234);
  normq->add_task(task237);

  auto tensor238 = vector<shared_ptr<Tensor>>{I312, Gamma6_(), t2};
  auto task238 = make_shared<Task238>(tensor238, pindex);
  task237->add_dep(task238);
  task238->add_dep(task234);
  normq->add_task(task238);

  vector<IndexRange> I314_index = {closed_, virt_, closed_, active_};
  auto I314 = make_shared<Tensor>(I314_index);
  auto tensor239 = vector<shared_ptr<Tensor>>{n, I314};
  auto task239 = make_shared<Task239>(tensor239, pindex);
  task239->add_dep(task234);
  normq->add_task(task239);

  vector<IndexRange> I315_index = {closed_, virt_, closed_, active_};
  auto I315 = make_shared<Tensor>(I315_index);
  auto tensor240 = vector<shared_ptr<Tensor>>{I314, Gamma13_(), I315};
  auto task240 = make_shared<Task240>(tensor240, pindex);
  task239->add_dep(task240);
  task240->add_dep(task234);
  normq->add_task(task240);

  auto tensor241 = vector<shared_ptr<Tensor>>{I315, t2};
  auto task241 = make_shared<Task241>(tensor241, pindex);
  task240->add_dep(task241);
  task241->add_dep(task234);
  normq->add_task(task241);

  vector<IndexRange> I318_index = {virt_, closed_, active_, active_};
  auto I318 = make_shared<Tensor>(I318_index);
  auto tensor242 = vector<shared_ptr<Tensor>>{n, I318};
  auto task242 = make_shared<Task242>(tensor242, pindex);
  task242->add_dep(task234);
  normq->add_task(task242);

  auto tensor243 = vector<shared_ptr<Tensor>>{I318, Gamma24_(), t2};
  auto task243 = make_shared<Task243>(tensor243, pindex);
  task242->add_dep(task243);
  task243->add_dep(task234);
  normq->add_task(task243);

  auto tensor244 = vector<shared_ptr<Tensor>>{I318, Gamma23_(), t2};
  auto task244 = make_shared<Task244>(tensor244, pindex);
  task242->add_dep(task244);
  task244->add_dep(task234);
  normq->add_task(task244);

  vector<IndexRange> I322_index = {virt_, closed_, active_, active_};
  auto I322 = make_shared<Tensor>(I322_index);
  auto tensor245 = vector<shared_ptr<Tensor>>{n, I322};
  auto task245 = make_shared<Task245>(tensor245, pindex);
  task245->add_dep(task234);
  normq->add_task(task245);

  vector<IndexRange> I323_index = {active_, virt_, closed_, active_};
  auto I323 = make_shared<Tensor>(I323_index);
  auto tensor246 = vector<shared_ptr<Tensor>>{I322, Gamma23_(), I323};
  auto task246 = make_shared<Task246>(tensor246, pindex);
  task245->add_dep(task246);
  task246->add_dep(task234);
  normq->add_task(task246);

  auto tensor247 = vector<shared_ptr<Tensor>>{I323, t2};
  auto task247 = make_shared<Task247>(tensor247, pindex);
  task246->add_dep(task247);
  task247->add_dep(task234);
  normq->add_task(task247);

  vector<IndexRange> I326_index = {virt_, active_, active_, active_};
  auto I326 = make_shared<Tensor>(I326_index);
  auto tensor248 = vector<shared_ptr<Tensor>>{n, I326};
  auto task248 = make_shared<Task248>(tensor248, pindex);
  task248->add_dep(task234);
  normq->add_task(task248);

  auto tensor249 = vector<shared_ptr<Tensor>>{I326, Gamma42_(), t2};
  auto task249 = make_shared<Task249>(tensor249, pindex);
  task248->add_dep(task249);
  task249->add_dep(task234);
  normq->add_task(task249);

  shared_ptr<Tensor> I328;
  if (diagonal) {
    vector<IndexRange> I328_index = {closed_, virt_, closed_, virt_};
    I328 = make_shared<Tensor>(I328_index);
  }
  shared_ptr<Task250> task250;
  if (diagonal) {
    auto tensor250 = vector<shared_ptr<Tensor>>{n, I328};
    task250 = make_shared<Task250>(tensor250, pindex);
    task250->add_dep(task234);
    normq->add_task(task250);
  }

  shared_ptr<Task251> task251;
  if (diagonal) {
    auto tensor251 = vector<shared_ptr<Tensor>>{I328, t2};
    task251 = make_shared<Task251>(tensor251, pindex);
    task250->add_dep(task251);
    task251->add_dep(task234);
    normq->add_task(task251);
  }

  vector<IndexRange> I330_index = {virt_, closed_, virt_, active_};
  auto I330 = make_shared<Tensor>(I330_index);
  auto tensor252 = vector<shared_ptr<Tensor>>{n, I330};
  auto task252 = make_shared<Task252>(tensor252, pindex);
  task252->add_dep(task234);
  normq->add_task(task252);

  vector<IndexRange> I331_index = {active_, virt_, closed_, virt_};
  auto I331 = make_shared<Tensor>(I331_index);
  auto tensor253 = vector<shared_ptr<Tensor>>{I330, Gamma21_(), I331};
  auto task253 = make_shared<Task253>(tensor253, pindex);
  task252->add_dep(task253);
  task253->add_dep(task234);
  normq->add_task(task253);

  auto tensor254 = vector<shared_ptr<Tensor>>{I331, t2};
  auto task254 = make_shared<Task254>(tensor254, pindex);
  task253->add_dep(task254);
  task254->add_dep(task234);
  normq->add_task(task254);

  vector<IndexRange> I334_index = {virt_, virt_, active_, active_};
  auto I334 = make_shared<Tensor>(I334_index);
  auto tensor255 = vector<shared_ptr<Tensor>>{n, I334};
  auto task255 = make_shared<Task255>(tensor255, pindex);
  task255->add_dep(task234);
  normq->add_task(task255);

  auto tensor256 = vector<shared_ptr<Tensor>>{I334, Gamma40_(), t2};
  auto task256 = make_shared<Task256>(tensor256, pindex);
  task255->add_dep(task256);
  task256->add_dep(task234);
  normq->add_task(task256);

  return normq;
}

#endif
