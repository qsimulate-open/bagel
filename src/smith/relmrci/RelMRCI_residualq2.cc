//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_residualq2.cc
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


#include <src/smith/relmrci/RelMRCI.h>
#include <src/smith/relmrci/RelMRCI_tasks4.h>
#include <src/smith/relmrci/RelMRCI_tasks5.h>
#include <src/smith/relmrci/RelMRCI_tasks6.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void RelMRCI::RelMRCI::make_residualq2(shared_ptr<Queue> residualq, shared_ptr<Task> task83, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I24_index = {closed_, closed_, active_, virt_};
  auto I24 = make_shared<Tensor>(I24_index);
  auto tensor159 = vector<shared_ptr<Tensor>>{r, I24};
  auto task159 = make_shared<Task159>(tensor159, pindex);
  task159->add_dep(task83);
  residualq->add_task(task159);

  vector<IndexRange> I25_index = {closed_, closed_, active_, active_};
  auto I25 = make_shared<Tensor>(I25_index);
  auto tensor160 = vector<shared_ptr<Tensor>>{I24, h1_, I25};
  auto task160 = make_shared<Task160>(tensor160, pindex);
  task159->add_dep(task160);
  task160->add_dep(task83);
  residualq->add_task(task160);

  auto tensor161 = vector<shared_ptr<Tensor>>{I25, Gamma2_(), t2};
  auto task161 = make_shared<Task161>(tensor161, pindex);
  task160->add_dep(task161);
  task161->add_dep(task83);
  residualq->add_task(task161);

  vector<IndexRange> I28_index = {closed_, active_};
  auto I28 = make_shared<Tensor>(I28_index);
  auto tensor162 = vector<shared_ptr<Tensor>>{I24, h1_, I28};
  auto task162 = make_shared<Task162>(tensor162, pindex);
  task159->add_dep(task162);
  task162->add_dep(task83);
  residualq->add_task(task162);

  auto tensor163 = vector<shared_ptr<Tensor>>{I28, Gamma9_(), t2};
  auto task163 = make_shared<Task163>(tensor163, pindex);
  task162->add_dep(task163);
  task163->add_dep(task83);
  residualq->add_task(task163);

  vector<IndexRange> I31_index = {closed_, active_};
  auto I31 = make_shared<Tensor>(I31_index);
  auto tensor164 = vector<shared_ptr<Tensor>>{I24, h1_, I31};
  auto task164 = make_shared<Task164>(tensor164, pindex);
  task159->add_dep(task164);
  task164->add_dep(task83);
  residualq->add_task(task164);

  auto tensor165 = vector<shared_ptr<Tensor>>{I31, Gamma9_(), t2};
  auto task165 = make_shared<Task165>(tensor165, pindex);
  task164->add_dep(task165);
  task165->add_dep(task83);
  residualq->add_task(task165);

  vector<IndexRange> I34_index = {closed_, virt_, closed_, active_};
  auto I34 = make_shared<Tensor>(I34_index);
  auto tensor166 = vector<shared_ptr<Tensor>>{I24, Gamma11_(), I34};
  auto task166 = make_shared<Task166>(tensor166, pindex);
  task159->add_dep(task166);
  task166->add_dep(task83);
  residualq->add_task(task166);

  auto tensor167 = vector<shared_ptr<Tensor>>{I34, t2, h1_};
  auto task167 = make_shared<Task167>(tensor167, pindex);
  task166->add_dep(task167);
  task167->add_dep(task83);
  residualq->add_task(task167);

  auto tensor168 = vector<shared_ptr<Tensor>>{I34, t2, h1_};
  auto task168 = make_shared<Task168>(tensor168, pindex);
  task166->add_dep(task168);
  task168->add_dep(task83);
  residualq->add_task(task168);

  auto tensor169 = vector<shared_ptr<Tensor>>{I34, t2, h1_};
  auto task169 = make_shared<Task169>(tensor169, pindex);
  task166->add_dep(task169);
  task169->add_dep(task83);
  residualq->add_task(task169);

  auto tensor170 = vector<shared_ptr<Tensor>>{I34, t2, h1_};
  auto task170 = make_shared<Task170>(tensor170, pindex);
  task166->add_dep(task170);
  task170->add_dep(task83);
  residualq->add_task(task170);

  auto tensor171 = vector<shared_ptr<Tensor>>{I34, t2, h1_};
  auto task171 = make_shared<Task171>(tensor171, pindex);
  task166->add_dep(task171);
  task171->add_dep(task83);
  residualq->add_task(task171);

  auto tensor172 = vector<shared_ptr<Tensor>>{I34, t2, h1_};
  auto task172 = make_shared<Task172>(tensor172, pindex);
  task166->add_dep(task172);
  task172->add_dep(task83);
  residualq->add_task(task172);

  auto tensor173 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task173 = make_shared<Task173>(tensor173, pindex);
  task166->add_dep(task173);
  task173->add_dep(task83);
  residualq->add_task(task173);

  auto tensor174 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task174 = make_shared<Task174>(tensor174, pindex);
  task166->add_dep(task174);
  task174->add_dep(task83);
  residualq->add_task(task174);

  auto tensor175 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task175 = make_shared<Task175>(tensor175, pindex);
  task166->add_dep(task175);
  task175->add_dep(task83);
  residualq->add_task(task175);

  auto tensor176 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task176 = make_shared<Task176>(tensor176, pindex);
  task166->add_dep(task176);
  task176->add_dep(task83);
  residualq->add_task(task176);

  auto tensor177 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task177 = make_shared<Task177>(tensor177, pindex);
  task166->add_dep(task177);
  task177->add_dep(task83);
  residualq->add_task(task177);

  auto tensor178 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task178 = make_shared<Task178>(tensor178, pindex);
  task166->add_dep(task178);
  task178->add_dep(task83);
  residualq->add_task(task178);

  auto tensor179 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task179 = make_shared<Task179>(tensor179, pindex);
  task166->add_dep(task179);
  task179->add_dep(task83);
  residualq->add_task(task179);

  auto tensor180 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task180 = make_shared<Task180>(tensor180, pindex);
  task166->add_dep(task180);
  task180->add_dep(task83);
  residualq->add_task(task180);

  auto tensor181 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task181 = make_shared<Task181>(tensor181, pindex);
  task166->add_dep(task181);
  task181->add_dep(task83);
  residualq->add_task(task181);

  auto tensor182 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task182 = make_shared<Task182>(tensor182, pindex);
  task166->add_dep(task182);
  task182->add_dep(task83);
  residualq->add_task(task182);

  vector<IndexRange> I505_index = {virt_, active_, closed_, closed_};
  auto I505 = make_shared<Tensor>(I505_index);
  auto tensor183 = vector<shared_ptr<Tensor>>{I34, t2, I505};
  auto task183 = make_shared<Task183>(tensor183, pindex);
  task166->add_dep(task183);
  task183->add_dep(task83);
  residualq->add_task(task183);

  auto tensor184 = vector<shared_ptr<Tensor>>{I505, v2_};
  auto task184 = make_shared<Task184>(tensor184, pindex);
  task183->add_dep(task184);
  task184->add_dep(task83);
  residualq->add_task(task184);

  auto tensor185 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task185 = make_shared<Task185>(tensor185, pindex);
  task166->add_dep(task185);
  task185->add_dep(task83);
  residualq->add_task(task185);

  vector<IndexRange> I514_index = {virt_, active_, closed_, closed_};
  auto I514 = make_shared<Tensor>(I514_index);
  auto tensor186 = vector<shared_ptr<Tensor>>{I34, t2, I514};
  auto task186 = make_shared<Task186>(tensor186, pindex);
  task166->add_dep(task186);
  task186->add_dep(task83);
  residualq->add_task(task186);

  auto tensor187 = vector<shared_ptr<Tensor>>{I514, v2_};
  auto task187 = make_shared<Task187>(tensor187, pindex);
  task186->add_dep(task187);
  task187->add_dep(task83);
  residualq->add_task(task187);

  vector<IndexRange> I517_index = {virt_, active_, closed_, closed_};
  auto I517 = make_shared<Tensor>(I517_index);
  auto tensor188 = vector<shared_ptr<Tensor>>{I34, t2, I517};
  auto task188 = make_shared<Task188>(tensor188, pindex);
  task166->add_dep(task188);
  task188->add_dep(task83);
  residualq->add_task(task188);

  auto tensor189 = vector<shared_ptr<Tensor>>{I517, v2_};
  auto task189 = make_shared<Task189>(tensor189, pindex);
  task188->add_dep(task189);
  task189->add_dep(task83);
  residualq->add_task(task189);

  auto tensor190 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task190 = make_shared<Task190>(tensor190, pindex);
  task166->add_dep(task190);
  task190->add_dep(task83);
  residualq->add_task(task190);

  auto tensor191 = vector<shared_ptr<Tensor>>{I34, t2, v2_};
  auto task191 = make_shared<Task191>(tensor191, pindex);
  task166->add_dep(task191);
  task191->add_dep(task83);
  residualq->add_task(task191);

  vector<IndexRange> I52_index = {closed_, virt_, active_, active_};
  auto I52 = make_shared<Tensor>(I52_index);
  auto tensor192 = vector<shared_ptr<Tensor>>{I24, h1_, I52};
  auto task192 = make_shared<Task192>(tensor192, pindex);
  task159->add_dep(task192);
  task192->add_dep(task83);
  residualq->add_task(task192);

  auto tensor193 = vector<shared_ptr<Tensor>>{I52, Gamma9_(), t2};
  auto task193 = make_shared<Task193>(tensor193, pindex);
  task192->add_dep(task193);
  task193->add_dep(task83);
  residualq->add_task(task193);

  vector<IndexRange> I55_index = {closed_, virt_, active_, active_};
  auto I55 = make_shared<Tensor>(I55_index);
  auto tensor194 = vector<shared_ptr<Tensor>>{I24, h1_, I55};
  auto task194 = make_shared<Task194>(tensor194, pindex);
  task159->add_dep(task194);
  task194->add_dep(task83);
  residualq->add_task(task194);

  auto tensor195 = vector<shared_ptr<Tensor>>{I55, Gamma9_(), t2};
  auto task195 = make_shared<Task195>(tensor195, pindex);
  task194->add_dep(task195);
  task195->add_dep(task83);
  residualq->add_task(task195);

  vector<IndexRange> I58_index = {virt_, active_};
  auto I58 = make_shared<Tensor>(I58_index);
  auto tensor196 = vector<shared_ptr<Tensor>>{I24, t2, I58};
  auto task196 = make_shared<Task196>(tensor196, pindex);
  task159->add_dep(task196);
  task196->add_dep(task83);
  residualq->add_task(task196);

  auto tensor197 = vector<shared_ptr<Tensor>>{I58, Gamma11_(), h1_};
  auto task197 = make_shared<Task197>(tensor197, pindex);
  task196->add_dep(task197);
  task197->add_dep(task83);
  residualq->add_task(task197);

  auto tensor198 = vector<shared_ptr<Tensor>>{I58, Gamma160_(), v2_};
  auto task198 = make_shared<Task198>(tensor198, pindex);
  task196->add_dep(task198);
  task198->add_dep(task83);
  residualq->add_task(task198);

  auto tensor199 = vector<shared_ptr<Tensor>>{I58, Gamma9_(), v2_};
  auto task199 = make_shared<Task199>(tensor199, pindex);
  task196->add_dep(task199);
  task199->add_dep(task83);
  residualq->add_task(task199);

  vector<IndexRange> I61_index = {active_, virt_};
  auto I61 = make_shared<Tensor>(I61_index);
  auto tensor200 = vector<shared_ptr<Tensor>>{I24, t2, I61};
  auto task200 = make_shared<Task200>(tensor200, pindex);
  task159->add_dep(task200);
  task200->add_dep(task83);
  residualq->add_task(task200);

  auto tensor201 = vector<shared_ptr<Tensor>>{I61, h1_, Gamma11_()};
  auto task201 = make_shared<Task201>(tensor201, pindex);
  task200->add_dep(task201);
  task201->add_dep(task83);
  residualq->add_task(task201);

  auto tensor202 = vector<shared_ptr<Tensor>>{I61, Gamma160_(), v2_};
  auto task202 = make_shared<Task202>(tensor202, pindex);
  task200->add_dep(task202);
  task202->add_dep(task83);
  residualq->add_task(task202);

  auto tensor203 = vector<shared_ptr<Tensor>>{I61, Gamma9_(), v2_};
  auto task203 = make_shared<Task203>(tensor203, pindex);
  task200->add_dep(task203);
  task203->add_dep(task83);
  residualq->add_task(task203);

  vector<IndexRange> I306_index = {virt_, active_, active_, active_};
  auto I306 = make_shared<Tensor>(I306_index);
  auto tensor204 = vector<shared_ptr<Tensor>>{I24, t2, I306};
  auto task204 = make_shared<Task204>(tensor204, pindex);
  task159->add_dep(task204);
  task204->add_dep(task83);
  residualq->add_task(task204);

  auto tensor205 = vector<shared_ptr<Tensor>>{I306, Gamma99_(), v2_};
  auto task205 = make_shared<Task205>(tensor205, pindex);
  task204->add_dep(task205);
  task205->add_dep(task83);
  residualq->add_task(task205);

  auto tensor206 = vector<shared_ptr<Tensor>>{I306, Gamma66_(), v2_};
  auto task206 = make_shared<Task206>(tensor206, pindex);
  task204->add_dep(task206);
  task206->add_dep(task83);
  residualq->add_task(task206);

  vector<IndexRange> I312_index = {closed_, closed_, active_, active_};
  auto I312 = make_shared<Tensor>(I312_index);
  auto tensor207 = vector<shared_ptr<Tensor>>{I24, v2_, I312};
  auto task207 = make_shared<Task207>(tensor207, pindex);
  task159->add_dep(task207);
  task207->add_dep(task83);
  residualq->add_task(task207);

  auto tensor208 = vector<shared_ptr<Tensor>>{I312, Gamma0_(), t2};
  auto task208 = make_shared<Task208>(tensor208, pindex);
  task207->add_dep(task208);
  task208->add_dep(task83);
  residualq->add_task(task208);

  vector<IndexRange> I315_index = {closed_, closed_, active_, active_};
  auto I315 = make_shared<Tensor>(I315_index);
  auto tensor209 = vector<shared_ptr<Tensor>>{I24, v2_, I315};
  auto task209 = make_shared<Task209>(tensor209, pindex);
  task159->add_dep(task209);
  task209->add_dep(task83);
  residualq->add_task(task209);

  auto tensor210 = vector<shared_ptr<Tensor>>{I315, Gamma0_(), t2};
  auto task210 = make_shared<Task210>(tensor210, pindex);
  task209->add_dep(task210);
  task210->add_dep(task83);
  residualq->add_task(task210);

  vector<IndexRange> I318_index = {closed_, closed_, active_, active_};
  auto I318 = make_shared<Tensor>(I318_index);
  auto tensor211 = vector<shared_ptr<Tensor>>{I24, v2_, I318};
  auto task211 = make_shared<Task211>(tensor211, pindex);
  task159->add_dep(task211);
  task211->add_dep(task83);
  residualq->add_task(task211);

  auto tensor212 = vector<shared_ptr<Tensor>>{I318, Gamma0_(), t2};
  auto task212 = make_shared<Task212>(tensor212, pindex);
  task211->add_dep(task212);
  task212->add_dep(task83);
  residualq->add_task(task212);

  vector<IndexRange> I321_index = {closed_, closed_, active_, active_};
  auto I321 = make_shared<Tensor>(I321_index);
  auto tensor213 = vector<shared_ptr<Tensor>>{I24, v2_, I321};
  auto task213 = make_shared<Task213>(tensor213, pindex);
  task159->add_dep(task213);
  task213->add_dep(task83);
  residualq->add_task(task213);

  auto tensor214 = vector<shared_ptr<Tensor>>{I321, Gamma2_(), t2};
  auto task214 = make_shared<Task214>(tensor214, pindex);
  task213->add_dep(task214);
  task214->add_dep(task83);
  residualq->add_task(task214);

  vector<IndexRange> I324_index = {closed_, active_, active_, active_};
  auto I324 = make_shared<Tensor>(I324_index);
  auto tensor215 = vector<shared_ptr<Tensor>>{I24, v2_, I324};
  auto task215 = make_shared<Task215>(tensor215, pindex);
  task159->add_dep(task215);
  task215->add_dep(task83);
  residualq->add_task(task215);

  auto tensor216 = vector<shared_ptr<Tensor>>{I324, Gamma105_(), t2};
  auto task216 = make_shared<Task216>(tensor216, pindex);
  task215->add_dep(task216);
  task216->add_dep(task83);
  residualq->add_task(task216);

  vector<IndexRange> I327_index = {closed_, active_, active_, active_};
  auto I327 = make_shared<Tensor>(I327_index);
  auto tensor217 = vector<shared_ptr<Tensor>>{I24, v2_, I327};
  auto task217 = make_shared<Task217>(tensor217, pindex);
  task159->add_dep(task217);
  task217->add_dep(task83);
  residualq->add_task(task217);

  auto tensor218 = vector<shared_ptr<Tensor>>{I327, Gamma105_(), t2};
  auto task218 = make_shared<Task218>(tensor218, pindex);
  task217->add_dep(task218);
  task218->add_dep(task83);
  residualq->add_task(task218);

  vector<IndexRange> I330_index = {closed_, active_, active_, active_};
  auto I330 = make_shared<Tensor>(I330_index);
  auto tensor219 = vector<shared_ptr<Tensor>>{I24, v2_, I330};
  auto task219 = make_shared<Task219>(tensor219, pindex);
  task159->add_dep(task219);
  task219->add_dep(task83);
  residualq->add_task(task219);

  auto tensor220 = vector<shared_ptr<Tensor>>{I330, Gamma1_(), t2};
  auto task220 = make_shared<Task220>(tensor220, pindex);
  task219->add_dep(task220);
  task220->add_dep(task83);
  residualq->add_task(task220);

  vector<IndexRange> I333_index = {closed_, active_, active_, active_};
  auto I333 = make_shared<Tensor>(I333_index);
  auto tensor221 = vector<shared_ptr<Tensor>>{I24, v2_, I333};
  auto task221 = make_shared<Task221>(tensor221, pindex);
  task159->add_dep(task221);
  task221->add_dep(task83);
  residualq->add_task(task221);

  auto tensor222 = vector<shared_ptr<Tensor>>{I333, Gamma65_(), t2};
  auto task222 = make_shared<Task222>(tensor222, pindex);
  task221->add_dep(task222);
  task222->add_dep(task83);
  residualq->add_task(task222);

  vector<IndexRange> I336_index = {closed_, active_, active_, active_};
  auto I336 = make_shared<Tensor>(I336_index);
  auto tensor223 = vector<shared_ptr<Tensor>>{I24, v2_, I336};
  auto task223 = make_shared<Task223>(tensor223, pindex);
  task159->add_dep(task223);
  task223->add_dep(task83);
  residualq->add_task(task223);

  auto tensor224 = vector<shared_ptr<Tensor>>{I336, Gamma105_(), t2};
  auto task224 = make_shared<Task224>(tensor224, pindex);
  task223->add_dep(task224);
  task224->add_dep(task83);
  residualq->add_task(task224);

  vector<IndexRange> I339_index = {closed_, active_, active_, active_};
  auto I339 = make_shared<Tensor>(I339_index);
  auto tensor225 = vector<shared_ptr<Tensor>>{I24, v2_, I339};
  auto task225 = make_shared<Task225>(tensor225, pindex);
  task159->add_dep(task225);
  task225->add_dep(task83);
  residualq->add_task(task225);

  auto tensor226 = vector<shared_ptr<Tensor>>{I339, Gamma110_(), t2};
  auto task226 = make_shared<Task226>(tensor226, pindex);
  task225->add_dep(task226);
  task226->add_dep(task83);
  residualq->add_task(task226);

  vector<IndexRange> I342_index = {closed_, active_, active_, active_};
  auto I342 = make_shared<Tensor>(I342_index);
  auto tensor227 = vector<shared_ptr<Tensor>>{I24, v2_, I342};
  auto task227 = make_shared<Task227>(tensor227, pindex);
  task159->add_dep(task227);
  task227->add_dep(task83);
  residualq->add_task(task227);

  auto tensor228 = vector<shared_ptr<Tensor>>{I342, Gamma105_(), t2};
  auto task228 = make_shared<Task228>(tensor228, pindex);
  task227->add_dep(task228);
  task228->add_dep(task83);
  residualq->add_task(task228);

  vector<IndexRange> I345_index = {closed_, active_, active_, active_};
  auto I345 = make_shared<Tensor>(I345_index);
  auto tensor229 = vector<shared_ptr<Tensor>>{I24, v2_, I345};
  auto task229 = make_shared<Task229>(tensor229, pindex);
  task159->add_dep(task229);
  task229->add_dep(task83);
  residualq->add_task(task229);

  auto tensor230 = vector<shared_ptr<Tensor>>{I345, Gamma105_(), t2};
  auto task230 = make_shared<Task230>(tensor230, pindex);
  task229->add_dep(task230);
  task230->add_dep(task83);
  residualq->add_task(task230);

  vector<IndexRange> I348_index = {closed_, active_};
  auto I348 = make_shared<Tensor>(I348_index);
  auto tensor231 = vector<shared_ptr<Tensor>>{I24, v2_, I348};
  auto task231 = make_shared<Task231>(tensor231, pindex);
  task159->add_dep(task231);
  task231->add_dep(task83);
  residualq->add_task(task231);

  auto tensor232 = vector<shared_ptr<Tensor>>{I348, Gamma9_(), t2};
  auto task232 = make_shared<Task232>(tensor232, pindex);
  task231->add_dep(task232);
  task232->add_dep(task83);
  residualq->add_task(task232);

  vector<IndexRange> I351_index = {closed_, active_};
  auto I351 = make_shared<Tensor>(I351_index);
  auto tensor233 = vector<shared_ptr<Tensor>>{I24, v2_, I351};
  auto task233 = make_shared<Task233>(tensor233, pindex);
  task159->add_dep(task233);
  task233->add_dep(task83);
  residualq->add_task(task233);

  auto tensor234 = vector<shared_ptr<Tensor>>{I351, Gamma9_(), t2};
  auto task234 = make_shared<Task234>(tensor234, pindex);
  task233->add_dep(task234);
  task234->add_dep(task83);
  residualq->add_task(task234);

  vector<IndexRange> I354_index = {closed_, closed_, active_, active_};
  auto I354 = make_shared<Tensor>(I354_index);
  auto tensor235 = vector<shared_ptr<Tensor>>{I24, t2, I354};
  auto task235 = make_shared<Task235>(tensor235, pindex);
  task159->add_dep(task235);
  task235->add_dep(task83);
  residualq->add_task(task235);

  vector<IndexRange> I355_index = {closed_, closed_, active_, active_};
  auto I355 = make_shared<Tensor>(I355_index);
  auto tensor236 = vector<shared_ptr<Tensor>>{I354, Gamma160_(), I355};
  auto task236 = make_shared<Task236>(tensor236, pindex);
  task235->add_dep(task236);
  task236->add_dep(task83);
  residualq->add_task(task236);

  auto tensor237 = vector<shared_ptr<Tensor>>{I355, v2_};
  auto task237 = make_shared<Task237>(tensor237, pindex);
  task236->add_dep(task237);
  task237->add_dep(task83);
  residualq->add_task(task237);

  auto tensor238 = vector<shared_ptr<Tensor>>{I354, Gamma0_(), v2_};
  auto task238 = make_shared<Task238>(tensor238, pindex);
  task235->add_dep(task238);
  task238->add_dep(task83);
  residualq->add_task(task238);

  vector<IndexRange> I357_index = {closed_, closed_, active_, active_};
  auto I357 = make_shared<Tensor>(I357_index);
  auto tensor239 = vector<shared_ptr<Tensor>>{I24, t2, I357};
  auto task239 = make_shared<Task239>(tensor239, pindex);
  task159->add_dep(task239);
  task239->add_dep(task83);
  residualq->add_task(task239);

  vector<IndexRange> I358_index = {closed_, closed_, active_, active_};
  auto I358 = make_shared<Tensor>(I358_index);
  auto tensor240 = vector<shared_ptr<Tensor>>{I357, Gamma160_(), I358};
  auto task240 = make_shared<Task240>(tensor240, pindex);
  task239->add_dep(task240);
  task240->add_dep(task83);
  residualq->add_task(task240);

  auto tensor241 = vector<shared_ptr<Tensor>>{I358, v2_};
  auto task241 = make_shared<Task241>(tensor241, pindex);
  task240->add_dep(task241);
  task241->add_dep(task83);
  residualq->add_task(task241);

  auto tensor242 = vector<shared_ptr<Tensor>>{I357, Gamma2_(), v2_};
  auto task242 = make_shared<Task242>(tensor242, pindex);
  task239->add_dep(task242);
  task242->add_dep(task83);
  residualq->add_task(task242);

  auto tensor243 = vector<shared_ptr<Tensor>>{I357, Gamma128_(), v2_};
  auto task243 = make_shared<Task243>(tensor243, pindex);
  task239->add_dep(task243);
  task243->add_dep(task83);
  residualq->add_task(task243);

  vector<IndexRange> I360_index = {closed_, closed_, active_, active_};
  auto I360 = make_shared<Tensor>(I360_index);
  auto tensor244 = vector<shared_ptr<Tensor>>{I24, t2, I360};
  auto task244 = make_shared<Task244>(tensor244, pindex);
  task159->add_dep(task244);
  task244->add_dep(task83);
  residualq->add_task(task244);

  vector<IndexRange> I361_index = {closed_, closed_, active_, active_};
  auto I361 = make_shared<Tensor>(I361_index);
  auto tensor245 = vector<shared_ptr<Tensor>>{I360, Gamma160_(), I361};
  auto task245 = make_shared<Task245>(tensor245, pindex);
  task244->add_dep(task245);
  task245->add_dep(task83);
  residualq->add_task(task245);

  auto tensor246 = vector<shared_ptr<Tensor>>{I361, v2_};
  auto task246 = make_shared<Task246>(tensor246, pindex);
  task245->add_dep(task246);
  task246->add_dep(task83);
  residualq->add_task(task246);

  auto tensor247 = vector<shared_ptr<Tensor>>{I360, Gamma2_(), v2_};
  auto task247 = make_shared<Task247>(tensor247, pindex);
  task244->add_dep(task247);
  task247->add_dep(task83);
  residualq->add_task(task247);

  auto tensor248 = vector<shared_ptr<Tensor>>{I360, Gamma128_(), v2_};
  auto task248 = make_shared<Task248>(tensor248, pindex);
  task244->add_dep(task248);
  task248->add_dep(task83);
  residualq->add_task(task248);

  vector<IndexRange> I363_index = {virt_, virt_, active_, active_};
  auto I363 = make_shared<Tensor>(I363_index);
  auto tensor249 = vector<shared_ptr<Tensor>>{I24, t2, I363};
  auto task249 = make_shared<Task249>(tensor249, pindex);
  task159->add_dep(task249);
  task249->add_dep(task83);
  residualq->add_task(task249);

  vector<IndexRange> I364_index = {virt_, virt_, active_, active_};
  auto I364 = make_shared<Tensor>(I364_index);
  auto tensor250 = vector<shared_ptr<Tensor>>{I363, Gamma160_(), I364};
  auto task250 = make_shared<Task250>(tensor250, pindex);
  task249->add_dep(task250);
  task250->add_dep(task83);
  residualq->add_task(task250);

  auto tensor251 = vector<shared_ptr<Tensor>>{I364, v2_};
  auto task251 = make_shared<Task251>(tensor251, pindex);
  task250->add_dep(task251);
  task251->add_dep(task83);
  residualq->add_task(task251);

  auto tensor252 = vector<shared_ptr<Tensor>>{I363, Gamma2_(), v2_};
  auto task252 = make_shared<Task252>(tensor252, pindex);
  task249->add_dep(task252);
  task252->add_dep(task83);
  residualq->add_task(task252);

  auto tensor253 = vector<shared_ptr<Tensor>>{I363, Gamma128_(), v2_};
  auto task253 = make_shared<Task253>(tensor253, pindex);
  task249->add_dep(task253);
  task253->add_dep(task83);
  residualq->add_task(task253);

  vector<IndexRange> I366_index = {closed_, closed_, active_, active_};
  auto I366 = make_shared<Tensor>(I366_index);
  auto tensor254 = vector<shared_ptr<Tensor>>{I24, t2, I366};
  auto task254 = make_shared<Task254>(tensor254, pindex);
  task159->add_dep(task254);
  task254->add_dep(task83);
  residualq->add_task(task254);

  vector<IndexRange> I367_index = {closed_, closed_, active_, active_};
  auto I367 = make_shared<Tensor>(I367_index);
  auto tensor255 = vector<shared_ptr<Tensor>>{I366, Gamma160_(), I367};
  auto task255 = make_shared<Task255>(tensor255, pindex);
  task254->add_dep(task255);
  task255->add_dep(task83);
  residualq->add_task(task255);

  auto tensor256 = vector<shared_ptr<Tensor>>{I367, v2_};
  auto task256 = make_shared<Task256>(tensor256, pindex);
  task255->add_dep(task256);
  task256->add_dep(task83);
  residualq->add_task(task256);

  auto tensor257 = vector<shared_ptr<Tensor>>{I366, Gamma2_(), v2_};
  auto task257 = make_shared<Task257>(tensor257, pindex);
  task254->add_dep(task257);
  task257->add_dep(task83);
  residualq->add_task(task257);

  auto tensor258 = vector<shared_ptr<Tensor>>{I366, Gamma128_(), v2_};
  auto task258 = make_shared<Task258>(tensor258, pindex);
  task254->add_dep(task258);
  task258->add_dep(task83);
  residualq->add_task(task258);

  vector<IndexRange> I369_index = {virt_, virt_, active_, active_};
  auto I369 = make_shared<Tensor>(I369_index);
  auto tensor259 = vector<shared_ptr<Tensor>>{I24, t2, I369};
  auto task259 = make_shared<Task259>(tensor259, pindex);
  task159->add_dep(task259);
  task259->add_dep(task83);
  residualq->add_task(task259);

  vector<IndexRange> I370_index = {virt_, virt_, active_, active_};
  auto I370 = make_shared<Tensor>(I370_index);
  auto tensor260 = vector<shared_ptr<Tensor>>{I369, Gamma160_(), I370};
  auto task260 = make_shared<Task260>(tensor260, pindex);
  task259->add_dep(task260);
  task260->add_dep(task83);
  residualq->add_task(task260);

  auto tensor261 = vector<shared_ptr<Tensor>>{I370, v2_};
  auto task261 = make_shared<Task261>(tensor261, pindex);
  task260->add_dep(task261);
  task261->add_dep(task83);
  residualq->add_task(task261);

  auto tensor262 = vector<shared_ptr<Tensor>>{I369, Gamma0_(), v2_};
  auto task262 = make_shared<Task262>(tensor262, pindex);
  task259->add_dep(task262);
  task262->add_dep(task83);
  residualq->add_task(task262);

  vector<IndexRange> I456_index = {closed_, active_, active_, active_};
  auto I456 = make_shared<Tensor>(I456_index);
  auto tensor263 = vector<shared_ptr<Tensor>>{I24, t2, I456};
  auto task263 = make_shared<Task263>(tensor263, pindex);
  task159->add_dep(task263);
  task263->add_dep(task83);
  residualq->add_task(task263);

  auto tensor264 = vector<shared_ptr<Tensor>>{I456, Gamma105_(), v2_};
  auto task264 = make_shared<Task264>(tensor264, pindex);
  task263->add_dep(task264);
  task264->add_dep(task83);
  residualq->add_task(task264);

  auto tensor265 = vector<shared_ptr<Tensor>>{I456, Gamma151_(), v2_};
  auto task265 = make_shared<Task265>(tensor265, pindex);
  task263->add_dep(task265);
  task265->add_dep(task83);
  residualq->add_task(task265);

  vector<IndexRange> I459_index = {closed_, active_, active_, active_};
  auto I459 = make_shared<Tensor>(I459_index);
  auto tensor266 = vector<shared_ptr<Tensor>>{I24, t2, I459};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task159->add_dep(task266);
  task266->add_dep(task83);
  residualq->add_task(task266);

  auto tensor267 = vector<shared_ptr<Tensor>>{I459, Gamma105_(), v2_};
  auto task267 = make_shared<Task267>(tensor267, pindex);
  task266->add_dep(task267);
  task267->add_dep(task83);
  residualq->add_task(task267);

  auto tensor268 = vector<shared_ptr<Tensor>>{I459, Gamma151_(), v2_};
  auto task268 = make_shared<Task268>(tensor268, pindex);
  task266->add_dep(task268);
  task268->add_dep(task83);
  residualq->add_task(task268);

  vector<IndexRange> I468_index = {closed_, virt_, active_, active_};
  auto I468 = make_shared<Tensor>(I468_index);
  auto tensor269 = vector<shared_ptr<Tensor>>{I24, v2_, I468};
  auto task269 = make_shared<Task269>(tensor269, pindex);
  task159->add_dep(task269);
  task269->add_dep(task83);
  residualq->add_task(task269);

  auto tensor270 = vector<shared_ptr<Tensor>>{I468, Gamma9_(), t2};
  auto task270 = make_shared<Task270>(tensor270, pindex);
  task269->add_dep(task270);
  task270->add_dep(task83);
  residualq->add_task(task270);

  vector<IndexRange> I471_index = {closed_, virt_, active_, active_};
  auto I471 = make_shared<Tensor>(I471_index);
  auto tensor271 = vector<shared_ptr<Tensor>>{I24, v2_, I471};
  auto task271 = make_shared<Task271>(tensor271, pindex);
  task159->add_dep(task271);
  task271->add_dep(task83);
  residualq->add_task(task271);

  auto tensor272 = vector<shared_ptr<Tensor>>{I471, Gamma9_(), t2};
  auto task272 = make_shared<Task272>(tensor272, pindex);
  task271->add_dep(task272);
  task272->add_dep(task83);
  residualq->add_task(task272);

  vector<IndexRange> I474_index = {closed_, virt_, active_, active_};
  auto I474 = make_shared<Tensor>(I474_index);
  auto tensor273 = vector<shared_ptr<Tensor>>{I24, v2_, I474};
  auto task273 = make_shared<Task273>(tensor273, pindex);
  task159->add_dep(task273);
  task273->add_dep(task83);
  residualq->add_task(task273);

  auto tensor274 = vector<shared_ptr<Tensor>>{I474, Gamma9_(), t2};
  auto task274 = make_shared<Task274>(tensor274, pindex);
  task273->add_dep(task274);
  task274->add_dep(task83);
  residualq->add_task(task274);

  vector<IndexRange> I477_index = {closed_, virt_, active_, active_};
  auto I477 = make_shared<Tensor>(I477_index);
  auto tensor275 = vector<shared_ptr<Tensor>>{I24, v2_, I477};
  auto task275 = make_shared<Task275>(tensor275, pindex);
  task159->add_dep(task275);
  task275->add_dep(task83);
  residualq->add_task(task275);

  auto tensor276 = vector<shared_ptr<Tensor>>{I477, Gamma9_(), t2};
  auto task276 = make_shared<Task276>(tensor276, pindex);
  task275->add_dep(task276);
  task276->add_dep(task83);
  residualq->add_task(task276);

  vector<IndexRange> I480_index = {closed_, virt_, active_, active_};
  auto I480 = make_shared<Tensor>(I480_index);
  auto tensor277 = vector<shared_ptr<Tensor>>{I24, v2_, I480};
  auto task277 = make_shared<Task277>(tensor277, pindex);
  task159->add_dep(task277);
  task277->add_dep(task83);
  residualq->add_task(task277);

  auto tensor278 = vector<shared_ptr<Tensor>>{I480, Gamma9_(), t2};
  auto task278 = make_shared<Task278>(tensor278, pindex);
  task277->add_dep(task278);
  task278->add_dep(task83);
  residualq->add_task(task278);

  vector<IndexRange> I483_index = {closed_, virt_, active_, active_};
  auto I483 = make_shared<Tensor>(I483_index);
  auto tensor279 = vector<shared_ptr<Tensor>>{I24, v2_, I483};
  auto task279 = make_shared<Task279>(tensor279, pindex);
  task159->add_dep(task279);
  task279->add_dep(task83);
  residualq->add_task(task279);

  auto tensor280 = vector<shared_ptr<Tensor>>{I483, Gamma9_(), t2};
  auto task280 = make_shared<Task280>(tensor280, pindex);
  task279->add_dep(task280);
  task280->add_dep(task83);
  residualq->add_task(task280);

  vector<IndexRange> I486_index = {virt_, active_, active_, active_};
  auto I486 = make_shared<Tensor>(I486_index);
  auto tensor281 = vector<shared_ptr<Tensor>>{I24, v2_, I486};
  auto task281 = make_shared<Task281>(tensor281, pindex);
  task159->add_dep(task281);
  task281->add_dep(task83);
  residualq->add_task(task281);

  auto tensor282 = vector<shared_ptr<Tensor>>{I486, Gamma159_(), t2};
  auto task282 = make_shared<Task282>(tensor282, pindex);
  task281->add_dep(task282);
  task282->add_dep(task83);
  residualq->add_task(task282);

  vector<IndexRange> I501_index = {virt_, closed_, closed_, active_};
  auto I501 = make_shared<Tensor>(I501_index);
  auto tensor283 = vector<shared_ptr<Tensor>>{I24, t2, I501};
  auto task283 = make_shared<Task283>(tensor283, pindex);
  task159->add_dep(task283);
  task283->add_dep(task83);
  residualq->add_task(task283);

  auto tensor284 = vector<shared_ptr<Tensor>>{I501, Gamma11_(), v2_};
  auto task284 = make_shared<Task284>(tensor284, pindex);
  task283->add_dep(task284);
  task284->add_dep(task83);
  residualq->add_task(task284);

  vector<IndexRange> I531_index = {closed_, virt_, active_, active_};
  auto I531 = make_shared<Tensor>(I531_index);
  auto tensor285 = vector<shared_ptr<Tensor>>{I24, t2, I531};
  auto task285 = make_shared<Task285>(tensor285, pindex);
  task159->add_dep(task285);
  task285->add_dep(task83);
  residualq->add_task(task285);

  auto tensor286 = vector<shared_ptr<Tensor>>{I531, Gamma174_(), v2_};
  auto task286 = make_shared<Task286>(tensor286, pindex);
  task285->add_dep(task286);
  task286->add_dep(task83);
  residualq->add_task(task286);

  vector<IndexRange> I534_index = {closed_, virt_, active_, active_};
  auto I534 = make_shared<Tensor>(I534_index);
  auto tensor287 = vector<shared_ptr<Tensor>>{I24, t2, I534};
  auto task287 = make_shared<Task287>(tensor287, pindex);
  task159->add_dep(task287);
  task287->add_dep(task83);
  residualq->add_task(task287);

  auto tensor288 = vector<shared_ptr<Tensor>>{I534, Gamma9_(), v2_};
  auto task288 = make_shared<Task288>(tensor288, pindex);
  task287->add_dep(task288);
  task288->add_dep(task83);
  residualq->add_task(task288);

  vector<IndexRange> I537_index = {closed_, virt_, active_, active_};
  auto I537 = make_shared<Tensor>(I537_index);
  auto tensor289 = vector<shared_ptr<Tensor>>{I24, t2, I537};
  auto task289 = make_shared<Task289>(tensor289, pindex);
  task159->add_dep(task289);
  task289->add_dep(task83);
  residualq->add_task(task289);

  auto tensor290 = vector<shared_ptr<Tensor>>{I537, Gamma174_(), v2_};
  auto task290 = make_shared<Task290>(tensor290, pindex);
  task289->add_dep(task290);
  task290->add_dep(task83);
  residualq->add_task(task290);

  vector<IndexRange> I540_index = {closed_, virt_, active_, active_};
  auto I540 = make_shared<Tensor>(I540_index);
  auto tensor291 = vector<shared_ptr<Tensor>>{I24, t2, I540};
  auto task291 = make_shared<Task291>(tensor291, pindex);
  task159->add_dep(task291);
  task291->add_dep(task83);
  residualq->add_task(task291);

  auto tensor292 = vector<shared_ptr<Tensor>>{I540, Gamma174_(), v2_};
  auto task292 = make_shared<Task292>(tensor292, pindex);
  task291->add_dep(task292);
  task292->add_dep(task83);
  residualq->add_task(task292);

  vector<IndexRange> I1273_index = {closed_, virt_, closed_, active_};
  auto I1273 = make_shared<Tensor>(I1273_index);
  auto tensor293 = vector<shared_ptr<Tensor>>{I24, Gamma416_(), I1273};
  auto task293 = make_shared<Task293>(tensor293, pindex);
  task159->add_dep(task293);
  task293->add_dep(task83);
  residualq->add_task(task293);

  auto tensor294 = vector<shared_ptr<Tensor>>{I1273, t2};
  auto task294 = make_shared<Task294>(tensor294, pindex);
  task293->add_dep(task294);
  task294->add_dep(task83);
  residualq->add_task(task294);

  vector<IndexRange> I1277_index = {closed_, virt_, closed_, active_};
  auto I1277 = make_shared<Tensor>(I1277_index);
  auto tensor295 = vector<shared_ptr<Tensor>>{I24, Gamma418_(), I1277};
  auto task295 = make_shared<Task295>(tensor295, pindex);
  task159->add_dep(task295);
  task295->add_dep(task83);
  residualq->add_task(task295);

  auto tensor296 = vector<shared_ptr<Tensor>>{I1277, t2};
  auto task296 = make_shared<Task296>(tensor296, pindex);
  task295->add_dep(task296);
  task296->add_dep(task83);
  residualq->add_task(task296);
}

#endif
