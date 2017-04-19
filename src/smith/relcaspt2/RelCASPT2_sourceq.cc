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
  auto tensor203 = vector<shared_ptr<Tensor>>{s};
  auto task203 = make_shared<Task203>(tensor203, reset);
  sourceq->add_task(task203);

  vector<IndexRange> I288_index = {closed_, closed_, active_, active_};
  auto I288 = make_shared<Tensor>(I288_index);
  auto tensor204 = vector<shared_ptr<Tensor>>{s, I288};
  auto task204 = make_shared<Task204>(tensor204, pindex);
  task204->add_dep(task203);
  sourceq->add_task(task204);

  auto tensor205 = vector<shared_ptr<Tensor>>{I288, Gamma92_(), v2_};
  auto task205 = make_shared<Task205>(tensor205, pindex);
  task204->add_dep(task205);
  task205->add_dep(task203);
  sourceq->add_task(task205);

  vector<IndexRange> I290_index = {closed_, active_, active_, active_};
  auto I290 = make_shared<Tensor>(I290_index);
  auto tensor206 = vector<shared_ptr<Tensor>>{s, I290};
  auto task206 = make_shared<Task206>(tensor206, pindex);
  task206->add_dep(task203);
  sourceq->add_task(task206);

  auto tensor207 = vector<shared_ptr<Tensor>>{I290, Gamma105_(), v2_};
  auto task207 = make_shared<Task207>(tensor207, pindex);
  task206->add_dep(task207);
  task207->add_dep(task203);
  sourceq->add_task(task207);

  auto tensor208 = vector<shared_ptr<Tensor>>{I290, Gamma6_(), v2_};
  auto task208 = make_shared<Task208>(tensor208, pindex);
  task206->add_dep(task208);
  task208->add_dep(task203);
  sourceq->add_task(task208);

  auto tensor209 = vector<shared_ptr<Tensor>>{I290, Gamma7_(), h1_};
  auto task209 = make_shared<Task209>(tensor209, pindex);
  task206->add_dep(task209);
  task209->add_dep(task203);
  sourceq->add_task(task209);

  vector<IndexRange> I294_index = {closed_, closed_, virt_, active_};
  auto I294 = make_shared<Tensor>(I294_index);
  auto tensor210 = vector<shared_ptr<Tensor>>{s, I294};
  auto task210 = make_shared<Task210>(tensor210, pindex);
  task210->add_dep(task203);
  sourceq->add_task(task210);

  vector<IndexRange> I295_index = {closed_, active_, closed_, virt_};
  auto I295 = make_shared<Tensor>(I295_index);
  auto tensor211 = vector<shared_ptr<Tensor>>{I294, Gamma16_(), I295};
  auto task211 = make_shared<Task211>(tensor211, pindex);
  task210->add_dep(task211);
  task211->add_dep(task203);
  sourceq->add_task(task211);

  auto tensor212 = vector<shared_ptr<Tensor>>{I295, v2_};
  auto task212 = make_shared<Task212>(tensor212, pindex);
  task211->add_dep(task212);
  task212->add_dep(task203);
  sourceq->add_task(task212);

  vector<IndexRange> I298_index = {closed_, virt_, active_, active_};
  auto I298 = make_shared<Tensor>(I298_index);
  auto tensor213 = vector<shared_ptr<Tensor>>{s, I298};
  auto task213 = make_shared<Task213>(tensor213, pindex);
  task213->add_dep(task203);
  sourceq->add_task(task213);

  vector<IndexRange> I299_index = {closed_, virt_, active_, active_};
  auto I299 = make_shared<Tensor>(I299_index);
  auto tensor214 = vector<shared_ptr<Tensor>>{I298, Gamma35_(), I299};
  auto task214 = make_shared<Task214>(tensor214, pindex);
  task213->add_dep(task214);
  task214->add_dep(task203);
  sourceq->add_task(task214);

  auto tensor215 = vector<shared_ptr<Tensor>>{I299, v2_};
  auto task215 = make_shared<Task215>(tensor215, pindex);
  task214->add_dep(task215);
  task215->add_dep(task203);
  sourceq->add_task(task215);

  auto tensor216 = vector<shared_ptr<Tensor>>{I298, Gamma29_(), v2_};
  auto task216 = make_shared<Task216>(tensor216, pindex);
  task213->add_dep(task216);
  task216->add_dep(task203);
  sourceq->add_task(task216);

  auto tensor217 = vector<shared_ptr<Tensor>>{I298, Gamma32_(), v2_};
  auto task217 = make_shared<Task217>(tensor217, pindex);
  task213->add_dep(task217);
  task217->add_dep(task203);
  sourceq->add_task(task217);

  auto tensor218 = vector<shared_ptr<Tensor>>{I298, Gamma38_(), h1_};
  auto task218 = make_shared<Task218>(tensor218, pindex);
  task213->add_dep(task218);
  task218->add_dep(task203);
  sourceq->add_task(task218);

  vector<IndexRange> I306_index = {closed_, virt_, active_, active_};
  auto I306 = make_shared<Tensor>(I306_index);
  auto tensor219 = vector<shared_ptr<Tensor>>{s, I306};
  auto task219 = make_shared<Task219>(tensor219, pindex);
  task219->add_dep(task203);
  sourceq->add_task(task219);

  vector<IndexRange> I307_index = {closed_, virt_, active_, active_};
  auto I307 = make_shared<Tensor>(I307_index);
  auto tensor220 = vector<shared_ptr<Tensor>>{I306, Gamma35_(), I307};
  auto task220 = make_shared<Task220>(tensor220, pindex);
  task219->add_dep(task220);
  task220->add_dep(task203);
  sourceq->add_task(task220);

  auto tensor221 = vector<shared_ptr<Tensor>>{I307, v2_};
  auto task221 = make_shared<Task221>(tensor221, pindex);
  task220->add_dep(task221);
  task221->add_dep(task203);
  sourceq->add_task(task221);

  auto tensor222 = vector<shared_ptr<Tensor>>{I306, Gamma7_(), v2_};
  auto task222 = make_shared<Task222>(tensor222, pindex);
  task219->add_dep(task222);
  task222->add_dep(task203);
  sourceq->add_task(task222);

  auto tensor223 = vector<shared_ptr<Tensor>>{I306, Gamma38_(), h1_};
  auto task223 = make_shared<Task223>(tensor223, pindex);
  task219->add_dep(task223);
  task223->add_dep(task203);
  sourceq->add_task(task223);

  vector<IndexRange> I314_index = {virt_, active_, active_, active_};
  auto I314 = make_shared<Tensor>(I314_index);
  auto tensor224 = vector<shared_ptr<Tensor>>{s, I314};
  auto task224 = make_shared<Task224>(tensor224, pindex);
  task224->add_dep(task203);
  sourceq->add_task(task224);

  auto tensor225 = vector<shared_ptr<Tensor>>{I314, Gamma59_(), v2_};
  auto task225 = make_shared<Task225>(tensor225, pindex);
  task224->add_dep(task225);
  task225->add_dep(task203);
  sourceq->add_task(task225);

  auto tensor226 = vector<shared_ptr<Tensor>>{I314, Gamma57_(), v2_};
  auto task226 = make_shared<Task226>(tensor226, pindex);
  task224->add_dep(task226);
  task226->add_dep(task203);
  sourceq->add_task(task226);

  auto tensor227 = vector<shared_ptr<Tensor>>{I314, Gamma60_(), h1_};
  auto task227 = make_shared<Task227>(tensor227, pindex);
  task224->add_dep(task227);
  task227->add_dep(task203);
  sourceq->add_task(task227);

  shared_ptr<Tensor> I318;
  if (diagonal) {
    vector<IndexRange> I318_index = {closed_, virt_, closed_, virt_};
    I318 = make_shared<Tensor>(I318_index);
  }
  shared_ptr<Task228> task228;
  if (diagonal) {
    auto tensor228 = vector<shared_ptr<Tensor>>{s, I318};
    task228 = make_shared<Task228>(tensor228, pindex);
    task228->add_dep(task203);
    sourceq->add_task(task228);
  }

  shared_ptr<Task229> task229;
  if (diagonal) {
    auto tensor229 = vector<shared_ptr<Tensor>>{I318, v2_};
    task229 = make_shared<Task229>(tensor229, pindex);
    task228->add_dep(task229);
    task229->add_dep(task203);
    sourceq->add_task(task229);
  }

  vector<IndexRange> I320_index = {virt_, closed_, virt_, active_};
  auto I320 = make_shared<Tensor>(I320_index);
  auto tensor230 = vector<shared_ptr<Tensor>>{s, I320};
  auto task230 = make_shared<Task230>(tensor230, pindex);
  task230->add_dep(task203);
  sourceq->add_task(task230);

  vector<IndexRange> I321_index = {active_, virt_, closed_, virt_};
  auto I321 = make_shared<Tensor>(I321_index);
  auto tensor231 = vector<shared_ptr<Tensor>>{I320, Gamma38_(), I321};
  auto task231 = make_shared<Task231>(tensor231, pindex);
  task230->add_dep(task231);
  task231->add_dep(task203);
  sourceq->add_task(task231);

  auto tensor232 = vector<shared_ptr<Tensor>>{I321, v2_};
  auto task232 = make_shared<Task232>(tensor232, pindex);
  task231->add_dep(task232);
  task232->add_dep(task203);
  sourceq->add_task(task232);

  vector<IndexRange> I324_index = {virt_, virt_, active_, active_};
  auto I324 = make_shared<Tensor>(I324_index);
  auto tensor233 = vector<shared_ptr<Tensor>>{s, I324};
  auto task233 = make_shared<Task233>(tensor233, pindex);
  task233->add_dep(task203);
  sourceq->add_task(task233);

  auto tensor234 = vector<shared_ptr<Tensor>>{I324, Gamma60_(), v2_};
  auto task234 = make_shared<Task234>(tensor234, pindex);
  task233->add_dep(task234);
  task234->add_dep(task203);
  sourceq->add_task(task234);

  return sourceq;
}


#endif
