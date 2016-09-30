//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASA_sourceqq.cc
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
#include <src/smith/relcasa/RelCASA_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelCASA::RelCASA::make_sourceq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  auto tensor202 = vector<shared_ptr<Tensor>>{s};
  auto task202 = make_shared<Task202>(tensor202, reset);
  sourceq->add_task(task202);

  vector<IndexRange> I264_index = {closed_, closed_, active_, active_};
  auto I264 = make_shared<Tensor>(I264_index);
  auto tensor203 = vector<shared_ptr<Tensor>>{s, I264};
  auto task203 = make_shared<Task203>(tensor203, pindex);
  task203->add_dep(task202);
  sourceq->add_task(task203);

  auto tensor204 = vector<shared_ptr<Tensor>>{I264, Gamma96_(), v2_};
  auto task204 = make_shared<Task204>(tensor204, pindex);
  task203->add_dep(task204);
  task204->add_dep(task202);
  sourceq->add_task(task204);

  vector<IndexRange> I266_index = {closed_, active_, active_, active_};
  auto I266 = make_shared<Tensor>(I266_index);
  auto tensor205 = vector<shared_ptr<Tensor>>{s, I266};
  auto task205 = make_shared<Task205>(tensor205, pindex);
  task205->add_dep(task202);
  sourceq->add_task(task205);

  auto tensor206 = vector<shared_ptr<Tensor>>{I266, Gamma97_(), v2_};
  auto task206 = make_shared<Task206>(tensor206, pindex);
  task205->add_dep(task206);
  task206->add_dep(task202);
  sourceq->add_task(task206);

  auto tensor207 = vector<shared_ptr<Tensor>>{I266, Gamma6_(), v2_};
  auto task207 = make_shared<Task207>(tensor207, pindex);
  task205->add_dep(task207);
  task207->add_dep(task202);
  sourceq->add_task(task207);

  auto tensor208 = vector<shared_ptr<Tensor>>{I266, Gamma3_(), h1_};
  auto task208 = make_shared<Task208>(tensor208, pindex);
  task205->add_dep(task208);
  task208->add_dep(task202);
  sourceq->add_task(task208);

  vector<IndexRange> I270_index = {closed_, closed_, virt_, active_};
  auto I270 = make_shared<Tensor>(I270_index);
  auto tensor209 = vector<shared_ptr<Tensor>>{s, I270};
  auto task209 = make_shared<Task209>(tensor209, pindex);
  task209->add_dep(task202);
  sourceq->add_task(task209);

  vector<IndexRange> I271_index = {closed_, active_, closed_, virt_};
  auto I271 = make_shared<Tensor>(I271_index);
  auto tensor210 = vector<shared_ptr<Tensor>>{I270, Gamma13_(), I271};
  auto task210 = make_shared<Task210>(tensor210, pindex);
  task209->add_dep(task210);
  task210->add_dep(task202);
  sourceq->add_task(task210);

  auto tensor211 = vector<shared_ptr<Tensor>>{I271, v2_};
  auto task211 = make_shared<Task211>(tensor211, pindex);
  task210->add_dep(task211);
  task211->add_dep(task202);
  sourceq->add_task(task211);

  vector<IndexRange> I274_index = {closed_, virt_, active_, active_};
  auto I274 = make_shared<Tensor>(I274_index);
  auto tensor212 = vector<shared_ptr<Tensor>>{s, I274};
  auto task212 = make_shared<Task212>(tensor212, pindex);
  task212->add_dep(task202);
  sourceq->add_task(task212);

  vector<IndexRange> I275_index = {closed_, virt_, active_, active_};
  auto I275 = make_shared<Tensor>(I275_index);
  auto tensor213 = vector<shared_ptr<Tensor>>{I274, Gamma23_(), I275};
  auto task213 = make_shared<Task213>(tensor213, pindex);
  task212->add_dep(task213);
  task213->add_dep(task202);
  sourceq->add_task(task213);

  auto tensor214 = vector<shared_ptr<Tensor>>{I275, v2_};
  auto task214 = make_shared<Task214>(tensor214, pindex);
  task213->add_dep(task214);
  task214->add_dep(task202);
  sourceq->add_task(task214);

  auto tensor215 = vector<shared_ptr<Tensor>>{I274, Gamma18_(), v2_};
  auto task215 = make_shared<Task215>(tensor215, pindex);
  task212->add_dep(task215);
  task215->add_dep(task202);
  sourceq->add_task(task215);

  auto tensor216 = vector<shared_ptr<Tensor>>{I274, Gamma24_(), v2_};
  auto task216 = make_shared<Task216>(tensor216, pindex);
  task212->add_dep(task216);
  task216->add_dep(task202);
  sourceq->add_task(task216);

  auto tensor217 = vector<shared_ptr<Tensor>>{I274, Gamma21_(), h1_};
  auto task217 = make_shared<Task217>(tensor217, pindex);
  task212->add_dep(task217);
  task217->add_dep(task202);
  sourceq->add_task(task217);

  vector<IndexRange> I282_index = {closed_, virt_, active_, active_};
  auto I282 = make_shared<Tensor>(I282_index);
  auto tensor218 = vector<shared_ptr<Tensor>>{s, I282};
  auto task218 = make_shared<Task218>(tensor218, pindex);
  task218->add_dep(task202);
  sourceq->add_task(task218);

  vector<IndexRange> I283_index = {closed_, virt_, active_, active_};
  auto I283 = make_shared<Tensor>(I283_index);
  auto tensor219 = vector<shared_ptr<Tensor>>{I282, Gamma23_(), I283};
  auto task219 = make_shared<Task219>(tensor219, pindex);
  task218->add_dep(task219);
  task219->add_dep(task202);
  sourceq->add_task(task219);

  auto tensor220 = vector<shared_ptr<Tensor>>{I283, v2_};
  auto task220 = make_shared<Task220>(tensor220, pindex);
  task219->add_dep(task220);
  task220->add_dep(task202);
  sourceq->add_task(task220);

  auto tensor221 = vector<shared_ptr<Tensor>>{I282, Gamma3_(), v2_};
  auto task221 = make_shared<Task221>(tensor221, pindex);
  task218->add_dep(task221);
  task221->add_dep(task202);
  sourceq->add_task(task221);

  auto tensor222 = vector<shared_ptr<Tensor>>{I282, Gamma21_(), h1_};
  auto task222 = make_shared<Task222>(tensor222, pindex);
  task218->add_dep(task222);
  task222->add_dep(task202);
  sourceq->add_task(task222);

  vector<IndexRange> I290_index = {virt_, active_, active_, active_};
  auto I290 = make_shared<Tensor>(I290_index);
  auto tensor223 = vector<shared_ptr<Tensor>>{s, I290};
  auto task223 = make_shared<Task223>(tensor223, pindex);
  task223->add_dep(task202);
  sourceq->add_task(task223);

  auto tensor224 = vector<shared_ptr<Tensor>>{I290, Gamma42_(), v2_};
  auto task224 = make_shared<Task224>(tensor224, pindex);
  task223->add_dep(task224);
  task224->add_dep(task202);
  sourceq->add_task(task224);

  auto tensor225 = vector<shared_ptr<Tensor>>{I290, Gamma39_(), v2_};
  auto task225 = make_shared<Task225>(tensor225, pindex);
  task223->add_dep(task225);
  task225->add_dep(task202);
  sourceq->add_task(task225);

  auto tensor226 = vector<shared_ptr<Tensor>>{I290, Gamma40_(), h1_};
  auto task226 = make_shared<Task226>(tensor226, pindex);
  task223->add_dep(task226);
  task226->add_dep(task202);
  sourceq->add_task(task226);

  shared_ptr<Tensor> I294;
  if (diagonal) {
    vector<IndexRange> I294_index = {closed_, virt_, closed_, virt_};
    I294 = make_shared<Tensor>(I294_index);
  }
  shared_ptr<Task227> task227;
  if (diagonal) {
    auto tensor227 = vector<shared_ptr<Tensor>>{s, I294};
    task227 = make_shared<Task227>(tensor227, pindex);
    task227->add_dep(task202);
    sourceq->add_task(task227);
  }

  shared_ptr<Task228> task228;
  if (diagonal) {
    auto tensor228 = vector<shared_ptr<Tensor>>{I294, v2_};
    task228 = make_shared<Task228>(tensor228, pindex);
    task227->add_dep(task228);
    task228->add_dep(task202);
    sourceq->add_task(task228);
  }

  vector<IndexRange> I296_index = {virt_, closed_, virt_, active_};
  auto I296 = make_shared<Tensor>(I296_index);
  auto tensor229 = vector<shared_ptr<Tensor>>{s, I296};
  auto task229 = make_shared<Task229>(tensor229, pindex);
  task229->add_dep(task202);
  sourceq->add_task(task229);

  vector<IndexRange> I297_index = {active_, virt_, closed_, virt_};
  auto I297 = make_shared<Tensor>(I297_index);
  auto tensor230 = vector<shared_ptr<Tensor>>{I296, Gamma21_(), I297};
  auto task230 = make_shared<Task230>(tensor230, pindex);
  task229->add_dep(task230);
  task230->add_dep(task202);
  sourceq->add_task(task230);

  auto tensor231 = vector<shared_ptr<Tensor>>{I297, v2_};
  auto task231 = make_shared<Task231>(tensor231, pindex);
  task230->add_dep(task231);
  task231->add_dep(task202);
  sourceq->add_task(task231);

  vector<IndexRange> I300_index = {virt_, virt_, active_, active_};
  auto I300 = make_shared<Tensor>(I300_index);
  auto tensor232 = vector<shared_ptr<Tensor>>{s, I300};
  auto task232 = make_shared<Task232>(tensor232, pindex);
  task232->add_dep(task202);
  sourceq->add_task(task232);

  auto tensor233 = vector<shared_ptr<Tensor>>{I300, Gamma40_(), v2_};
  auto task233 = make_shared<Task233>(tensor233, pindex);
  task232->add_dep(task233);
  task233->add_dep(task202);
  sourceq->add_task(task233);

  return sourceq;
}


#endif
