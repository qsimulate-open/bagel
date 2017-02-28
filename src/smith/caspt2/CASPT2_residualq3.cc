//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_residualq3.cc
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
#include <src/smith/caspt2/CASPT2_tasks4.h>
#include <src/smith/caspt2/CASPT2_tasks5.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void CASPT2::CASPT2::make_residualq3(shared_ptr<Queue> residualq, shared_ptr<Task> task69, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I204_index = {virt_, closed_, active_, virt_};
  auto I204 = make_shared<Tensor>(I204_index);
  auto tensor193 = vector<shared_ptr<Tensor>>{r, I204};
  auto task193 = make_shared<Task193>(tensor193, pindex);
  task193->add_dep(task69);
  residualq->add_task(task193);

  vector<IndexRange> I205_index = {virt_, closed_, active_, active_};
  auto I205 = make_shared<Tensor>(I205_index);
  auto tensor194 = vector<shared_ptr<Tensor>>{I204, f1_, I205};
  auto task194 = make_shared<Task194>(tensor194, pindex);
  task193->add_dep(task194);
  task194->add_dep(task69);
  residualq->add_task(task194);

  vector<IndexRange> I206_index = {active_, virt_, closed_, active_};
  auto I206 = make_shared<Tensor>(I206_index);
  auto tensor195 = vector<shared_ptr<Tensor>>{I205, Gamma35_(), I206};
  auto task195 = make_shared<Task195>(tensor195, pindex);
  task194->add_dep(task195);
  task195->add_dep(task69);
  residualq->add_task(task195);

  auto tensor196 = vector<shared_ptr<Tensor>>{I206, t2};
  auto task196 = make_shared<Task196>(tensor196, pindex);
  task195->add_dep(task196);
  task196->add_dep(task69);
  residualq->add_task(task196);

  vector<IndexRange> I208_index = {virt_, closed_, active_, active_};
  auto I208 = make_shared<Tensor>(I208_index);
  auto tensor197 = vector<shared_ptr<Tensor>>{I204, f1_, I208};
  auto task197 = make_shared<Task197>(tensor197, pindex);
  task193->add_dep(task197);
  task197->add_dep(task69);
  residualq->add_task(task197);

  auto tensor198 = vector<shared_ptr<Tensor>>{I208, Gamma32_(), t2};
  auto task198 = make_shared<Task198>(tensor198, pindex);
  task197->add_dep(task198);
  task198->add_dep(task69);
  residualq->add_task(task198);

  auto tensor199 = vector<shared_ptr<Tensor>>{I208, Gamma35_(), t2};
  auto task199 = make_shared<Task199>(tensor199, pindex);
  task197->add_dep(task199);
  task199->add_dep(task69);
  residualq->add_task(task199);

  vector<IndexRange> I217_index = {virt_, active_};
  auto I217 = make_shared<Tensor>(I217_index);
  auto tensor200 = vector<shared_ptr<Tensor>>{I204, f1_, I217};
  auto task200 = make_shared<Task200>(tensor200, pindex);
  task193->add_dep(task200);
  task200->add_dep(task69);
  residualq->add_task(task200);

  auto tensor201 = vector<shared_ptr<Tensor>>{I217, Gamma60_(), t2};
  auto task201 = make_shared<Task201>(tensor201, pindex);
  task200->add_dep(task201);
  task201->add_dep(task69);
  residualq->add_task(task201);

  vector<IndexRange> I220_index = {virt_, active_};
  auto I220 = make_shared<Tensor>(I220_index);
  auto tensor202 = vector<shared_ptr<Tensor>>{I204, f1_, I220};
  auto task202 = make_shared<Task202>(tensor202, pindex);
  task193->add_dep(task202);
  task202->add_dep(task69);
  residualq->add_task(task202);

  auto tensor203 = vector<shared_ptr<Tensor>>{I220, Gamma60_(), t2};
  auto task203 = make_shared<Task203>(tensor203, pindex);
  task202->add_dep(task203);
  task203->add_dep(task69);
  residualq->add_task(task203);

  vector<IndexRange> I223_index = {closed_, active_};
  auto I223 = make_shared<Tensor>(I223_index);
  auto tensor204 = vector<shared_ptr<Tensor>>{I204, t2, I223};
  auto task204 = make_shared<Task204>(tensor204, pindex);
  task193->add_dep(task204);
  task204->add_dep(task69);
  residualq->add_task(task204);

  auto tensor205 = vector<shared_ptr<Tensor>>{I223, Gamma38_(), f1_};
  auto task205 = make_shared<Task205>(tensor205, pindex);
  task204->add_dep(task205);
  task205->add_dep(task69);
  residualq->add_task(task205);

  vector<IndexRange> I226_index = {closed_, active_};
  auto I226 = make_shared<Tensor>(I226_index);
  auto tensor206 = vector<shared_ptr<Tensor>>{I204, t2, I226};
  auto task206 = make_shared<Task206>(tensor206, pindex);
  task193->add_dep(task206);
  task206->add_dep(task69);
  residualq->add_task(task206);

  auto tensor207 = vector<shared_ptr<Tensor>>{I226, Gamma38_(), f1_};
  auto task207 = make_shared<Task207>(tensor207, pindex);
  task206->add_dep(task207);
  task207->add_dep(task69);
  residualq->add_task(task207);

  vector<IndexRange> I229_index = {active_, virt_, closed_, virt_};
  auto I229 = make_shared<Tensor>(I229_index);
  auto tensor208 = vector<shared_ptr<Tensor>>{I204, Gamma79_(), I229};
  auto task208 = make_shared<Task208>(tensor208, pindex);
  task193->add_dep(task208);
  task208->add_dep(task69);
  residualq->add_task(task208);

  auto tensor209 = vector<shared_ptr<Tensor>>{I229, t2};
  auto task209 = make_shared<Task209>(tensor209, pindex);
  task208->add_dep(task209);
  task209->add_dep(task69);
  residualq->add_task(task209);

  vector<IndexRange> I233_index = {closed_, active_, virt_, virt_};
  auto I233 = make_shared<Tensor>(I233_index);
  auto tensor210 = vector<shared_ptr<Tensor>>{I204, Gamma38_(), I233};
  auto task210 = make_shared<Task210>(tensor210, pindex);
  task193->add_dep(task210);
  task210->add_dep(task69);
  residualq->add_task(task210);

  auto tensor211 = vector<shared_ptr<Tensor>>{I233, t2};
  auto task211 = make_shared<Task211>(tensor211, pindex, this->e0_);
  task210->add_dep(task211);
  task211->add_dep(task69);
  residualq->add_task(task211);

  auto tensor212 = vector<shared_ptr<Tensor>>{I233, t2, f1_};
  auto task212 = make_shared<Task212>(tensor212, pindex);
  task210->add_dep(task212);
  task212->add_dep(task69);
  residualq->add_task(task212);

  auto tensor213 = vector<shared_ptr<Tensor>>{I233, t2, f1_};
  auto task213 = make_shared<Task213>(tensor213, pindex);
  task210->add_dep(task213);
  task213->add_dep(task69);
  residualq->add_task(task213);

  auto tensor214 = vector<shared_ptr<Tensor>>{I233, t2, f1_};
  auto task214 = make_shared<Task214>(tensor214, pindex);
  task210->add_dep(task214);
  task214->add_dep(task69);
  residualq->add_task(task214);

  auto tensor215 = vector<shared_ptr<Tensor>>{I233, t2, f1_};
  auto task215 = make_shared<Task215>(tensor215, pindex);
  task210->add_dep(task215);
  task215->add_dep(task69);
  residualq->add_task(task215);

  auto tensor216 = vector<shared_ptr<Tensor>>{I233, t2, f1_};
  auto task216 = make_shared<Task216>(tensor216, pindex);
  task210->add_dep(task216);
  task216->add_dep(task69);
  residualq->add_task(task216);

  auto tensor217 = vector<shared_ptr<Tensor>>{I233, t2, f1_};
  auto task217 = make_shared<Task217>(tensor217, pindex);
  task210->add_dep(task217);
  task217->add_dep(task69);
  residualq->add_task(task217);

  vector<IndexRange> I251_index = {virt_, virt_, active_, active_};
  auto I251 = make_shared<Tensor>(I251_index);
  auto tensor218 = vector<shared_ptr<Tensor>>{I204, f1_, I251};
  auto task218 = make_shared<Task218>(tensor218, pindex);
  task193->add_dep(task218);
  task218->add_dep(task69);
  residualq->add_task(task218);

  auto tensor219 = vector<shared_ptr<Tensor>>{I251, Gamma60_(), t2};
  auto task219 = make_shared<Task219>(tensor219, pindex);
  task218->add_dep(task219);
  task219->add_dep(task69);
  residualq->add_task(task219);

  vector<IndexRange> I253_index = {virt_, active_, active_, virt_};
  auto I253 = make_shared<Tensor>(I253_index);
  auto tensor220 = vector<shared_ptr<Tensor>>{r, I253};
  auto task220 = make_shared<Task220>(tensor220, pindex);
  task220->add_dep(task69);
  residualq->add_task(task220);

  vector<IndexRange> I254_index = {virt_, active_, active_, active_};
  auto I254 = make_shared<Tensor>(I254_index);
  auto tensor221 = vector<shared_ptr<Tensor>>{I253, f1_, I254};
  auto task221 = make_shared<Task221>(tensor221, pindex);
  task220->add_dep(task221);
  task221->add_dep(task69);
  residualq->add_task(task221);

  auto tensor222 = vector<shared_ptr<Tensor>>{I254, Gamma59_(), t2};
  auto task222 = make_shared<Task222>(tensor222, pindex);
  task221->add_dep(task222);
  task222->add_dep(task69);
  residualq->add_task(task222);

  vector<IndexRange> I257_index = {active_, active_, virt_, virt_};
  auto I257 = make_shared<Tensor>(I257_index);
  auto tensor223 = vector<shared_ptr<Tensor>>{I253, Gamma60_(), I257};
  auto task223 = make_shared<Task223>(tensor223, pindex);
  task220->add_dep(task223);
  task223->add_dep(task69);
  residualq->add_task(task223);

  auto tensor224 = vector<shared_ptr<Tensor>>{I257, t2, f1_};
  auto task224 = make_shared<Task224>(tensor224, pindex);
  task223->add_dep(task224);
  task224->add_dep(task69);
  residualq->add_task(task224);

  auto tensor225 = vector<shared_ptr<Tensor>>{I257, t2, f1_};
  auto task225 = make_shared<Task225>(tensor225, pindex);
  task223->add_dep(task225);
  task225->add_dep(task69);
  residualq->add_task(task225);

  vector<IndexRange> I259_index = {virt_, virt_, active_, active_};
  auto I259 = make_shared<Tensor>(I259_index);
  auto tensor226 = vector<shared_ptr<Tensor>>{r, I259};
  auto task226 = make_shared<Task226>(tensor226, pindex);
  task226->add_dep(task69);
  residualq->add_task(task226);

  auto tensor227 = vector<shared_ptr<Tensor>>{I259, Gamma90_(), t2};
  auto task227 = make_shared<Task227>(tensor227, pindex);
  task226->add_dep(task227);
  task227->add_dep(task69);
  residualq->add_task(task227);

  auto tensor228 = vector<shared_ptr<Tensor>>{I259, Gamma60_(), t2};
  auto task228 = make_shared<Task228>(tensor228, pindex, this->e0_);
  task226->add_dep(task228);
  task228->add_dep(task69);
  residualq->add_task(task228);
}

#endif
