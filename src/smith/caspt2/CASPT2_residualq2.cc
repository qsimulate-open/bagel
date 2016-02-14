//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_residualq2.cc
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
#include <src/smith/caspt2/CASPT2_tasks3.h>
#include <src/smith/caspt2/CASPT2_tasks4.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void CASPT2::CASPT2::make_residualq2(shared_ptr<Queue> residualq, shared_ptr<Task> task69, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I120_index = {closed_, active_, active_, virt_};
  auto I120 = make_shared<Tensor>(I120_index);
  auto tensor144 = vector<shared_ptr<Tensor>>{r, I120};
  auto task144 = make_shared<Task144>(tensor144, pindex);
  task144->add_dep(task69);
  residualq->add_task(task144);

  vector<IndexRange> I121_index = {closed_, active_, active_, active_};
  auto I121 = make_shared<Tensor>(I121_index);
  auto tensor145 = vector<shared_ptr<Tensor>>{I120, f1_, I121};
  auto task145 = make_shared<Task145>(tensor145, pindex);
  task144->add_dep(task145);
  task145->add_dep(task69);
  residualq->add_task(task145);

  auto tensor146 = vector<shared_ptr<Tensor>>{I121, Gamma6_(), t2};
  auto task146 = make_shared<Task146>(tensor146, pindex);
  task145->add_dep(task146);
  task146->add_dep(task69);
  residualq->add_task(task146);

  vector<IndexRange> I124_index = {active_, virt_, closed_, active_};
  auto I124 = make_shared<Tensor>(I124_index);
  auto tensor147 = vector<shared_ptr<Tensor>>{I120, Gamma7_(), I124};
  auto task147 = make_shared<Task147>(tensor147, pindex);
  task144->add_dep(task147);
  task147->add_dep(task69);
  residualq->add_task(task147);

  auto tensor148 = vector<shared_ptr<Tensor>>{I124, t2, f1_};
  auto task148 = make_shared<Task148>(tensor148, pindex);
  task147->add_dep(task148);
  task148->add_dep(task69);
  residualq->add_task(task148);

  auto tensor149 = vector<shared_ptr<Tensor>>{I124, t2, f1_};
  auto task149 = make_shared<Task149>(tensor149, pindex);
  task147->add_dep(task149);
  task149->add_dep(task69);
  residualq->add_task(task149);

  vector<IndexRange> I130_index = {active_, virt_, closed_, active_};
  auto I130 = make_shared<Tensor>(I130_index);
  auto tensor150 = vector<shared_ptr<Tensor>>{I120, Gamma34_(), I130};
  auto task150 = make_shared<Task150>(tensor150, pindex);
  task144->add_dep(task150);
  task150->add_dep(task69);
  residualq->add_task(task150);

  auto tensor151 = vector<shared_ptr<Tensor>>{I130, t2};
  auto task151 = make_shared<Task151>(tensor151, pindex);
  task150->add_dep(task151);
  task151->add_dep(task69);
  residualq->add_task(task151);

  vector<IndexRange> I132_index = {closed_, active_, virt_, active_};
  auto I132 = make_shared<Tensor>(I132_index);
  auto tensor152 = vector<shared_ptr<Tensor>>{I120, Gamma35_(), I132};
  auto task152 = make_shared<Task152>(tensor152, pindex);
  task144->add_dep(task152);
  task152->add_dep(task69);
  residualq->add_task(task152);

  auto tensor153 = vector<shared_ptr<Tensor>>{I132, t2};
  auto task153 = make_shared<Task153>(tensor153, pindex, this->e0_);
  task152->add_dep(task153);
  task153->add_dep(task69);
  residualq->add_task(task153);

  auto tensor154 = vector<shared_ptr<Tensor>>{I132, t2, f1_};
  auto task154 = make_shared<Task154>(tensor154, pindex);
  task152->add_dep(task154);
  task154->add_dep(task69);
  residualq->add_task(task154);

  auto tensor155 = vector<shared_ptr<Tensor>>{I132, t2, f1_};
  auto task155 = make_shared<Task155>(tensor155, pindex);
  task152->add_dep(task155);
  task155->add_dep(task69);
  residualq->add_task(task155);

  auto tensor156 = vector<shared_ptr<Tensor>>{I132, t2, f1_};
  auto task156 = make_shared<Task156>(tensor156, pindex);
  task152->add_dep(task156);
  task156->add_dep(task69);
  residualq->add_task(task156);

  auto tensor157 = vector<shared_ptr<Tensor>>{I132, t2, f1_};
  auto task157 = make_shared<Task157>(tensor157, pindex);
  task152->add_dep(task157);
  task157->add_dep(task69);
  residualq->add_task(task157);

  auto tensor158 = vector<shared_ptr<Tensor>>{I132, t2, f1_};
  auto task158 = make_shared<Task158>(tensor158, pindex);
  task152->add_dep(task158);
  task158->add_dep(task69);
  residualq->add_task(task158);

  auto tensor159 = vector<shared_ptr<Tensor>>{I132, t2, f1_};
  auto task159 = make_shared<Task159>(tensor159, pindex);
  task152->add_dep(task159);
  task159->add_dep(task69);
  residualq->add_task(task159);

  vector<IndexRange> I146_index = {virt_, active_, active_, active_};
  auto I146 = make_shared<Tensor>(I146_index);
  auto tensor160 = vector<shared_ptr<Tensor>>{I120, f1_, I146};
  auto task160 = make_shared<Task160>(tensor160, pindex);
  task144->add_dep(task160);
  task160->add_dep(task69);
  residualq->add_task(task160);

  auto tensor161 = vector<shared_ptr<Tensor>>{I146, Gamma51_(), t2};
  auto task161 = make_shared<Task161>(tensor161, pindex);
  task160->add_dep(task161);
  task161->add_dep(task69);
  residualq->add_task(task161);

  vector<IndexRange> I149_index = {closed_, virt_};
  auto I149 = make_shared<Tensor>(I149_index);
  auto tensor162 = vector<shared_ptr<Tensor>>{I120, Gamma38_(), I149};
  auto task162 = make_shared<Task162>(tensor162, pindex);
  task144->add_dep(task162);
  task162->add_dep(task69);
  residualq->add_task(task162);

  auto tensor163 = vector<shared_ptr<Tensor>>{I149, t2, f1_};
  auto task163 = make_shared<Task163>(tensor163, pindex);
  task162->add_dep(task163);
  task163->add_dep(task69);
  residualq->add_task(task163);

  auto tensor164 = vector<shared_ptr<Tensor>>{I149, t2, f1_};
  auto task164 = make_shared<Task164>(tensor164, pindex);
  task162->add_dep(task164);
  task164->add_dep(task69);
  residualq->add_task(task164);

  vector<IndexRange> I160_index = {virt_, active_, active_, active_};
  auto I160 = make_shared<Tensor>(I160_index);
  auto tensor165 = vector<shared_ptr<Tensor>>{r, I160};
  auto task165 = make_shared<Task165>(tensor165, pindex);
  task165->add_dep(task69);
  residualq->add_task(task165);

  vector<IndexRange> I161_index = {active_, active_, virt_, active_};
  auto I161 = make_shared<Tensor>(I161_index);
  auto tensor166 = vector<shared_ptr<Tensor>>{I160, Gamma56_(), I161};
  auto task166 = make_shared<Task166>(tensor166, pindex);
  task165->add_dep(task166);
  task166->add_dep(task69);
  residualq->add_task(task166);

  auto tensor167 = vector<shared_ptr<Tensor>>{I161, t2, f1_};
  auto task167 = make_shared<Task167>(tensor167, pindex);
  task166->add_dep(task167);
  task167->add_dep(task69);
  residualq->add_task(task167);

  vector<IndexRange> I164_index = {active_, virt_, active_, active_};
  auto I164 = make_shared<Tensor>(I164_index);
  auto tensor168 = vector<shared_ptr<Tensor>>{I160, Gamma57_(), I164};
  auto task168 = make_shared<Task168>(tensor168, pindex);
  task165->add_dep(task168);
  task168->add_dep(task69);
  residualq->add_task(task168);

  auto tensor169 = vector<shared_ptr<Tensor>>{I164, t2, f1_};
  auto task169 = make_shared<Task169>(tensor169, pindex);
  task168->add_dep(task169);
  task169->add_dep(task69);
  residualq->add_task(task169);

  auto tensor170 = vector<shared_ptr<Tensor>>{I160, Gamma58_(), t2};
  auto task170 = make_shared<Task170>(tensor170, pindex);
  task165->add_dep(task170);
  task170->add_dep(task69);
  residualq->add_task(task170);

  vector<IndexRange> I169_index = {virt_, active_, active_, active_};
  auto I169 = make_shared<Tensor>(I169_index);
  auto tensor171 = vector<shared_ptr<Tensor>>{I160, Gamma59_(), I169};
  auto task171 = make_shared<Task171>(tensor171, pindex);
  task165->add_dep(task171);
  task171->add_dep(task69);
  residualq->add_task(task171);

  auto tensor172 = vector<shared_ptr<Tensor>>{I169, t2};
  auto task172 = make_shared<Task172>(tensor172, pindex, this->e0_);
  task171->add_dep(task172);
  task172->add_dep(task69);
  residualq->add_task(task172);

  auto tensor173 = vector<shared_ptr<Tensor>>{I169, t2, f1_};
  auto task173 = make_shared<Task173>(tensor173, pindex);
  task171->add_dep(task173);
  task173->add_dep(task69);
  residualq->add_task(task173);

  auto tensor174 = vector<shared_ptr<Tensor>>{I169, t2, f1_};
  auto task174 = make_shared<Task174>(tensor174, pindex);
  task171->add_dep(task174);
  task174->add_dep(task69);
  residualq->add_task(task174);

  vector<IndexRange> I172_index = {active_, virt_};
  auto I172 = make_shared<Tensor>(I172_index);
  auto tensor175 = vector<shared_ptr<Tensor>>{I160, Gamma60_(), I172};
  auto task175 = make_shared<Task175>(tensor175, pindex);
  task165->add_dep(task175);
  task175->add_dep(task69);
  residualq->add_task(task175);

  auto tensor176 = vector<shared_ptr<Tensor>>{I172, t2, f1_};
  auto task176 = make_shared<Task176>(tensor176, pindex);
  task175->add_dep(task176);
  task176->add_dep(task69);
  residualq->add_task(task176);

  auto tensor177 = vector<shared_ptr<Tensor>>{I172, t2, f1_};
  auto task177 = make_shared<Task177>(tensor177, pindex);
  task175->add_dep(task177);
  task177->add_dep(task69);
  residualq->add_task(task177);

  vector<IndexRange> I180_index = {virt_, closed_, virt_, closed_};
  auto I180 = make_shared<Tensor>(I180_index);
  auto tensor178 = vector<shared_ptr<Tensor>>{r, I180};
  auto task178 = make_shared<Task178>(tensor178, pindex);
  task178->add_dep(task69);
  residualq->add_task(task178);

  vector<IndexRange> I181_index = {virt_, active_};
  auto I181 = make_shared<Tensor>(I181_index);
  auto tensor179 = vector<shared_ptr<Tensor>>{I180, t2, I181};
  auto task179 = make_shared<Task179>(tensor179, pindex);
  task178->add_dep(task179);
  task179->add_dep(task69);
  residualq->add_task(task179);

  auto tensor180 = vector<shared_ptr<Tensor>>{I181, Gamma16_(), f1_};
  auto task180 = make_shared<Task180>(tensor180, pindex);
  task179->add_dep(task180);
  task180->add_dep(task69);
  residualq->add_task(task180);

  vector<IndexRange> I184_index = {virt_, active_};
  auto I184 = make_shared<Tensor>(I184_index);
  auto tensor181 = vector<shared_ptr<Tensor>>{I180, t2, I184};
  auto task181 = make_shared<Task181>(tensor181, pindex);
  task178->add_dep(task181);
  task181->add_dep(task69);
  residualq->add_task(task181);

  auto tensor182 = vector<shared_ptr<Tensor>>{I184, Gamma16_(), f1_};
  auto task182 = make_shared<Task182>(tensor182, pindex);
  task181->add_dep(task182);
  task182->add_dep(task69);
  residualq->add_task(task182);

  vector<IndexRange> I187_index = {virt_, closed_};
  auto I187 = make_shared<Tensor>(I187_index);
  auto tensor183 = vector<shared_ptr<Tensor>>{I180, f1_, I187};
  auto task183 = make_shared<Task183>(tensor183, pindex);
  task178->add_dep(task183);
  task183->add_dep(task69);
  residualq->add_task(task183);

  vector<IndexRange> I188_index = {active_, virt_, closed_, active_};
  auto I188 = make_shared<Tensor>(I188_index);
  auto tensor184 = vector<shared_ptr<Tensor>>{I187, Gamma38_(), I188};
  auto task184 = make_shared<Task184>(tensor184, pindex);
  task183->add_dep(task184);
  task184->add_dep(task69);
  residualq->add_task(task184);

  auto tensor185 = vector<shared_ptr<Tensor>>{I188, t2};
  auto task185 = make_shared<Task185>(tensor185, pindex);
  task184->add_dep(task185);
  task185->add_dep(task69);
  residualq->add_task(task185);

  vector<IndexRange> I190_index = {virt_, closed_};
  auto I190 = make_shared<Tensor>(I190_index);
  auto tensor186 = vector<shared_ptr<Tensor>>{I180, f1_, I190};
  auto task186 = make_shared<Task186>(tensor186, pindex);
  task178->add_dep(task186);
  task186->add_dep(task69);
  residualq->add_task(task186);

  vector<IndexRange> I191_index = {active_, virt_, closed_, active_};
  auto I191 = make_shared<Tensor>(I191_index);
  auto tensor187 = vector<shared_ptr<Tensor>>{I190, Gamma38_(), I191};
  auto task187 = make_shared<Task187>(tensor187, pindex);
  task186->add_dep(task187);
  task187->add_dep(task69);
  residualq->add_task(task187);

  auto tensor188 = vector<shared_ptr<Tensor>>{I191, t2};
  auto task188 = make_shared<Task188>(tensor188, pindex);
  task187->add_dep(task188);
  task188->add_dep(task69);
  residualq->add_task(task188);

  vector<IndexRange> I199_index = {closed_, active_};
  auto I199 = make_shared<Tensor>(I199_index);
  auto tensor189 = vector<shared_ptr<Tensor>>{I180, t2, I199};
  auto task189 = make_shared<Task189>(tensor189, pindex);
  task178->add_dep(task189);
  task189->add_dep(task69);
  residualq->add_task(task189);

  auto tensor190 = vector<shared_ptr<Tensor>>{I199, Gamma38_(), f1_};
  auto task190 = make_shared<Task190>(tensor190, pindex);
  task189->add_dep(task190);
  task190->add_dep(task69);
  residualq->add_task(task190);

  vector<IndexRange> I202_index = {closed_, active_};
  auto I202 = make_shared<Tensor>(I202_index);
  auto tensor191 = vector<shared_ptr<Tensor>>{I180, t2, I202};
  auto task191 = make_shared<Task191>(tensor191, pindex);
  task178->add_dep(task191);
  task191->add_dep(task69);
  residualq->add_task(task191);

  auto tensor192 = vector<shared_ptr<Tensor>>{I202, Gamma38_(), f1_};
  auto task192 = make_shared<Task192>(tensor192, pindex);
  task191->add_dep(task192);
  task192->add_dep(task69);
  residualq->add_task(task192);
}

#endif
