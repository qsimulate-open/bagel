//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASPT2_sourceqq.cc
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

shared_ptr<Queue> RelCASPT2::RelCASPT2::make_sourceq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  auto tensor190 = vector<shared_ptr<Tensor>>{s};
  auto task190 = make_shared<Task190>(tensor190, reset);
  sourceq->add_task(task190);

  vector<IndexRange> I288_index = {closed_, active_, active_, active_};
  auto I288 = make_shared<Tensor>(I288_index);
  auto tensor191 = vector<shared_ptr<Tensor>>{s, I288};
  auto task191 = make_shared<Task191>(tensor191, pindex);
  task191->add_dep(task190);
  sourceq->add_task(task191);

  auto tensor192 = vector<shared_ptr<Tensor>>{I288, Gamma7_(), h1_};
  auto task192 = make_shared<Task192>(tensor192, pindex);
  task191->add_dep(task192);
  task192->add_dep(task190);
  sourceq->add_task(task192);

  auto tensor193 = vector<shared_ptr<Tensor>>{I288, Gamma109_(), v2_};
  auto task193 = make_shared<Task193>(tensor193, pindex);
  task191->add_dep(task193);
  task193->add_dep(task190);
  sourceq->add_task(task193);

  auto tensor194 = vector<shared_ptr<Tensor>>{I288, Gamma6_(), v2_};
  auto task194 = make_shared<Task194>(tensor194, pindex);
  task191->add_dep(task194);
  task194->add_dep(task190);
  sourceq->add_task(task194);

  vector<IndexRange> I290_index = {closed_, virt_, active_, active_};
  auto I290 = make_shared<Tensor>(I290_index);
  auto tensor195 = vector<shared_ptr<Tensor>>{s, I290};
  auto task195 = make_shared<Task195>(tensor195, pindex);
  task195->add_dep(task190);
  sourceq->add_task(task195);

  auto tensor196 = vector<shared_ptr<Tensor>>{I290, Gamma38_(), h1_};
  auto task196 = make_shared<Task196>(tensor196, pindex);
  task195->add_dep(task196);
  task196->add_dep(task190);
  sourceq->add_task(task196);

  vector<IndexRange> I307_index = {closed_, virt_, active_, active_};
  auto I307 = make_shared<Tensor>(I307_index);
  auto tensor197 = vector<shared_ptr<Tensor>>{I290, Gamma35_(), I307};
  auto task197 = make_shared<Task197>(tensor197, pindex);
  task195->add_dep(task197);
  task197->add_dep(task190);
  sourceq->add_task(task197);

  auto tensor198 = vector<shared_ptr<Tensor>>{I307, v2_};
  auto task198 = make_shared<Task198>(tensor198, pindex);
  task197->add_dep(task198);
  task198->add_dep(task190);
  sourceq->add_task(task198);

  auto tensor199 = vector<shared_ptr<Tensor>>{I290, Gamma29_(), v2_};
  auto task199 = make_shared<Task199>(tensor199, pindex);
  task195->add_dep(task199);
  task199->add_dep(task190);
  sourceq->add_task(task199);

  auto tensor200 = vector<shared_ptr<Tensor>>{I290, Gamma32_(), v2_};
  auto task200 = make_shared<Task200>(tensor200, pindex);
  task195->add_dep(task200);
  task200->add_dep(task190);
  sourceq->add_task(task200);

  vector<IndexRange> I292_index = {closed_, virt_, active_, active_};
  auto I292 = make_shared<Tensor>(I292_index);
  auto tensor201 = vector<shared_ptr<Tensor>>{s, I292};
  auto task201 = make_shared<Task201>(tensor201, pindex);
  task201->add_dep(task190);
  sourceq->add_task(task201);

  auto tensor202 = vector<shared_ptr<Tensor>>{I292, Gamma38_(), h1_};
  auto task202 = make_shared<Task202>(tensor202, pindex);
  task201->add_dep(task202);
  task202->add_dep(task190);
  sourceq->add_task(task202);

  vector<IndexRange> I315_index = {closed_, virt_, active_, active_};
  auto I315 = make_shared<Tensor>(I315_index);
  auto tensor203 = vector<shared_ptr<Tensor>>{I292, Gamma35_(), I315};
  auto task203 = make_shared<Task203>(tensor203, pindex);
  task201->add_dep(task203);
  task203->add_dep(task190);
  sourceq->add_task(task203);

  auto tensor204 = vector<shared_ptr<Tensor>>{I315, v2_};
  auto task204 = make_shared<Task204>(tensor204, pindex);
  task203->add_dep(task204);
  task204->add_dep(task190);
  sourceq->add_task(task204);

  auto tensor205 = vector<shared_ptr<Tensor>>{I292, Gamma7_(), v2_};
  auto task205 = make_shared<Task205>(tensor205, pindex);
  task201->add_dep(task205);
  task205->add_dep(task190);
  sourceq->add_task(task205);

  vector<IndexRange> I294_index = {virt_, active_, active_, active_};
  auto I294 = make_shared<Tensor>(I294_index);
  auto tensor206 = vector<shared_ptr<Tensor>>{s, I294};
  auto task206 = make_shared<Task206>(tensor206, pindex);
  task206->add_dep(task190);
  sourceq->add_task(task206);

  auto tensor207 = vector<shared_ptr<Tensor>>{I294, Gamma60_(), h1_};
  auto task207 = make_shared<Task207>(tensor207, pindex);
  task206->add_dep(task207);
  task207->add_dep(task190);
  sourceq->add_task(task207);

  auto tensor208 = vector<shared_ptr<Tensor>>{I294, Gamma59_(), v2_};
  auto task208 = make_shared<Task208>(tensor208, pindex);
  task206->add_dep(task208);
  task208->add_dep(task190);
  sourceq->add_task(task208);

  auto tensor209 = vector<shared_ptr<Tensor>>{I294, Gamma57_(), v2_};
  auto task209 = make_shared<Task209>(tensor209, pindex);
  task206->add_dep(task209);
  task209->add_dep(task190);
  sourceq->add_task(task209);

  vector<IndexRange> I296_index = {closed_, closed_, active_, active_};
  auto I296 = make_shared<Tensor>(I296_index);
  auto tensor210 = vector<shared_ptr<Tensor>>{s, I296};
  auto task210 = make_shared<Task210>(tensor210, pindex);
  task210->add_dep(task190);
  sourceq->add_task(task210);

  auto tensor211 = vector<shared_ptr<Tensor>>{I296, Gamma92_(), v2_};
  auto task211 = make_shared<Task211>(tensor211, pindex);
  task210->add_dep(task211);
  task211->add_dep(task190);
  sourceq->add_task(task211);

  vector<IndexRange> I302_index = {closed_, closed_, virt_, active_};
  auto I302 = make_shared<Tensor>(I302_index);
  auto tensor212 = vector<shared_ptr<Tensor>>{s, I302};
  auto task212 = make_shared<Task212>(tensor212, pindex);
  task212->add_dep(task190);
  sourceq->add_task(task212);

  vector<IndexRange> I303_index = {closed_, active_, closed_, virt_};
  auto I303 = make_shared<Tensor>(I303_index);
  auto tensor213 = vector<shared_ptr<Tensor>>{I302, Gamma16_(), I303};
  auto task213 = make_shared<Task213>(tensor213, pindex);
  task212->add_dep(task213);
  task213->add_dep(task190);
  sourceq->add_task(task213);

  auto tensor214 = vector<shared_ptr<Tensor>>{I303, v2_};
  auto task214 = make_shared<Task214>(tensor214, pindex);
  task213->add_dep(task214);
  task214->add_dep(task190);
  sourceq->add_task(task214);

  shared_ptr<Tensor> I326;
  if (diagonal) {
    vector<IndexRange> I326_index = {closed_, virt_, closed_, virt_};
    I326 = make_shared<Tensor>(I326_index);
  }
  shared_ptr<Task215> task215;
  if (diagonal) {
    auto tensor215 = vector<shared_ptr<Tensor>>{s, I326};
    task215 = make_shared<Task215>(tensor215, pindex);
    task215->add_dep(task190);
    sourceq->add_task(task215);
  }

  shared_ptr<Task216> task216;
  if (diagonal) {
    auto tensor216 = vector<shared_ptr<Tensor>>{I326, v2_};
    task216 = make_shared<Task216>(tensor216, pindex);
    task215->add_dep(task216);
    task216->add_dep(task190);
    sourceq->add_task(task216);
  }

  vector<IndexRange> I328_index = {virt_, closed_, virt_, active_};
  auto I328 = make_shared<Tensor>(I328_index);
  auto tensor217 = vector<shared_ptr<Tensor>>{s, I328};
  auto task217 = make_shared<Task217>(tensor217, pindex);
  task217->add_dep(task190);
  sourceq->add_task(task217);

  vector<IndexRange> I329_index = {active_, virt_, closed_, virt_};
  auto I329 = make_shared<Tensor>(I329_index);
  auto tensor218 = vector<shared_ptr<Tensor>>{I328, Gamma38_(), I329};
  auto task218 = make_shared<Task218>(tensor218, pindex);
  task217->add_dep(task218);
  task218->add_dep(task190);
  sourceq->add_task(task218);

  auto tensor219 = vector<shared_ptr<Tensor>>{I329, v2_};
  auto task219 = make_shared<Task219>(tensor219, pindex);
  task218->add_dep(task219);
  task219->add_dep(task190);
  sourceq->add_task(task219);

  vector<IndexRange> I332_index = {virt_, virt_, active_, active_};
  auto I332 = make_shared<Tensor>(I332_index);
  auto tensor220 = vector<shared_ptr<Tensor>>{s, I332};
  auto task220 = make_shared<Task220>(tensor220, pindex);
  task220->add_dep(task190);
  sourceq->add_task(task220);

  auto tensor221 = vector<shared_ptr<Tensor>>{I332, Gamma60_(), v2_};
  auto task221 = make_shared<Task221>(tensor221, pindex);
  task220->add_dep(task221);
  task221->add_dep(task190);
  sourceq->add_task(task221);

  return sourceq;
}


#endif
