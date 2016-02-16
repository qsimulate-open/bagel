//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_residualqq.cc
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


#include <src/smith/mrci/MRCI.h>
#include <src/smith/mrci/MRCI_tasks3.h>
#include <src/smith/mrci/MRCI_tasks4.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void MRCI::MRCI::make_residualq2(shared_ptr<Queue> residualq, shared_ptr<Task> task108, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I9_index = {closed_, active_, active_, active_};
  auto I9 = make_shared<Tensor>(I9_index);
  auto tensor139 = vector<shared_ptr<Tensor>>{r, I9};
  auto task139 = make_shared<Task139>(tensor139, pindex);
  task139->add_dep(task108);
  residualq->add_task(task139);

  vector<IndexRange> I10_index = {active_, closed_, active_, active_};
  auto I10 = make_shared<Tensor>(I10_index);
  auto tensor140 = vector<shared_ptr<Tensor>>{I9, Gamma3_(), I10};
  auto task140 = make_shared<Task140>(tensor140, pindex);
  task139->add_dep(task140);
  task140->add_dep(task108);
  residualq->add_task(task140);

  auto tensor141 = vector<shared_ptr<Tensor>>{I10, t2, h1_};
  auto task141 = make_shared<Task141>(tensor141, pindex);
  task140->add_dep(task141);
  task141->add_dep(task108);
  residualq->add_task(task141);

  auto tensor142 = vector<shared_ptr<Tensor>>{I10, t2, v2_};
  auto task142 = make_shared<Task142>(tensor142, pindex);
  task140->add_dep(task142);
  task142->add_dep(task108);
  residualq->add_task(task142);

  auto tensor143 = vector<shared_ptr<Tensor>>{I10, t2, v2_};
  auto task143 = make_shared<Task143>(tensor143, pindex);
  task140->add_dep(task143);
  task143->add_dep(task108);
  residualq->add_task(task143);

  vector<IndexRange> I13_index = {closed_, active_, active_, active_};
  auto I13 = make_shared<Tensor>(I13_index);
  auto tensor144 = vector<shared_ptr<Tensor>>{I9, Gamma4_(), I13};
  auto task144 = make_shared<Task144>(tensor144, pindex);
  task139->add_dep(task144);
  task144->add_dep(task108);
  residualq->add_task(task144);

  auto tensor145 = vector<shared_ptr<Tensor>>{I13, t2, h1_};
  auto task145 = make_shared<Task145>(tensor145, pindex);
  task144->add_dep(task145);
  task145->add_dep(task108);
  residualq->add_task(task145);

  auto tensor146 = vector<shared_ptr<Tensor>>{I13, t2, h1_};
  auto task146 = make_shared<Task146>(tensor146, pindex);
  task144->add_dep(task146);
  task146->add_dep(task108);
  residualq->add_task(task146);

  auto tensor147 = vector<shared_ptr<Tensor>>{I13, t2, v2_};
  auto task147 = make_shared<Task147>(tensor147, pindex);
  task144->add_dep(task147);
  task147->add_dep(task108);
  residualq->add_task(task147);

  vector<IndexRange> I370_index = {virt_, active_, closed_, closed_};
  auto I370 = make_shared<Tensor>(I370_index);
  auto tensor148 = vector<shared_ptr<Tensor>>{I13, t2, I370};
  auto task148 = make_shared<Task148>(tensor148, pindex);
  task144->add_dep(task148);
  task148->add_dep(task108);
  residualq->add_task(task148);

  auto tensor149 = vector<shared_ptr<Tensor>>{I370, v2_};
  auto task149 = make_shared<Task149>(tensor149, pindex);
  task148->add_dep(task149);
  task149->add_dep(task108);
  residualq->add_task(task149);

  vector<IndexRange> I16_index = {closed_, active_};
  auto I16 = make_shared<Tensor>(I16_index);
  auto tensor150 = vector<shared_ptr<Tensor>>{I9, Gamma5_(), I16};
  auto task150 = make_shared<Task150>(tensor150, pindex);
  task139->add_dep(task150);
  task150->add_dep(task108);
  residualq->add_task(task150);

  auto tensor151 = vector<shared_ptr<Tensor>>{I16, t2, h1_};
  auto task151 = make_shared<Task151>(tensor151, pindex);
  task150->add_dep(task151);
  task151->add_dep(task108);
  residualq->add_task(task151);

  auto tensor152 = vector<shared_ptr<Tensor>>{I16, t2, h1_};
  auto task152 = make_shared<Task152>(tensor152, pindex);
  task150->add_dep(task152);
  task152->add_dep(task108);
  residualq->add_task(task152);

  auto tensor153 = vector<shared_ptr<Tensor>>{I16, t2, v2_};
  auto task153 = make_shared<Task153>(tensor153, pindex);
  task150->add_dep(task153);
  task153->add_dep(task108);
  residualq->add_task(task153);

  auto tensor154 = vector<shared_ptr<Tensor>>{I16, t2, v2_};
  auto task154 = make_shared<Task154>(tensor154, pindex);
  task150->add_dep(task154);
  task154->add_dep(task108);
  residualq->add_task(task154);

  auto tensor155 = vector<shared_ptr<Tensor>>{I16, t2, v2_};
  auto task155 = make_shared<Task155>(tensor155, pindex);
  task150->add_dep(task155);
  task155->add_dep(task108);
  residualq->add_task(task155);

  auto tensor156 = vector<shared_ptr<Tensor>>{I16, t2, v2_};
  auto task156 = make_shared<Task156>(tensor156, pindex);
  task150->add_dep(task156);
  task156->add_dep(task108);
  residualq->add_task(task156);

  vector<IndexRange> I22_index = {active_, active_, closed_, active_};
  auto I22 = make_shared<Tensor>(I22_index);
  auto tensor157 = vector<shared_ptr<Tensor>>{I9, Gamma7_(), I22};
  auto task157 = make_shared<Task157>(tensor157, pindex);
  task139->add_dep(task157);
  task157->add_dep(task108);
  residualq->add_task(task157);

  auto tensor158 = vector<shared_ptr<Tensor>>{I22, t2, h1_};
  auto task158 = make_shared<Task158>(tensor158, pindex);
  task157->add_dep(task158);
  task158->add_dep(task108);
  residualq->add_task(task158);

  auto tensor159 = vector<shared_ptr<Tensor>>{I22, t2, v2_};
  auto task159 = make_shared<Task159>(tensor159, pindex);
  task157->add_dep(task159);
  task159->add_dep(task108);
  residualq->add_task(task159);

  auto tensor160 = vector<shared_ptr<Tensor>>{I22, t2, v2_};
  auto task160 = make_shared<Task160>(tensor160, pindex);
  task157->add_dep(task160);
  task160->add_dep(task108);
  residualq->add_task(task160);

  vector<IndexRange> I300_index = {active_, active_, active_, closed_, active_, active_};
  auto I300 = make_shared<Tensor>(I300_index);
  auto tensor161 = vector<shared_ptr<Tensor>>{I9, Gamma97_(), I300};
  auto task161 = make_shared<Task161>(tensor161, pindex);
  task139->add_dep(task161);
  task161->add_dep(task108);
  residualq->add_task(task161);

  auto tensor162 = vector<shared_ptr<Tensor>>{I300, t2, v2_};
  auto task162 = make_shared<Task162>(tensor162, pindex);
  task161->add_dep(task162);
  task162->add_dep(task108);
  residualq->add_task(task162);

  vector<IndexRange> I303_index = {active_, active_, active_, closed_, active_, active_};
  auto I303 = make_shared<Tensor>(I303_index);
  auto tensor163 = vector<shared_ptr<Tensor>>{I9, Gamma98_(), I303};
  auto task163 = make_shared<Task163>(tensor163, pindex);
  task139->add_dep(task163);
  task163->add_dep(task108);
  residualq->add_task(task163);

  auto tensor164 = vector<shared_ptr<Tensor>>{I303, t2, v2_};
  auto task164 = make_shared<Task164>(tensor164, pindex);
  task163->add_dep(task164);
  task164->add_dep(task108);
  residualq->add_task(task164);

  vector<IndexRange> I309_index = {closed_, active_, active_, active_, active_, active_};
  auto I309 = make_shared<Tensor>(I309_index);
  auto tensor165 = vector<shared_ptr<Tensor>>{I9, Gamma100_(), I309};
  auto task165 = make_shared<Task165>(tensor165, pindex);
  task139->add_dep(task165);
  task165->add_dep(task108);
  residualq->add_task(task165);

  vector<IndexRange> I310_index = {closed_, closed_, active_, active_};
  auto I310 = make_shared<Tensor>(I310_index);
  auto tensor166 = vector<shared_ptr<Tensor>>{I309, t2, I310};
  auto task166 = make_shared<Task166>(tensor166, pindex);
  task165->add_dep(task166);
  task166->add_dep(task108);
  residualq->add_task(task166);

  auto tensor167 = vector<shared_ptr<Tensor>>{I310, v2_};
  auto task167 = make_shared<Task167>(tensor167, pindex);
  task166->add_dep(task167);
  task167->add_dep(task108);
  residualq->add_task(task167);

  auto tensor168 = vector<shared_ptr<Tensor>>{I309, t2, v2_};
  auto task168 = make_shared<Task168>(tensor168, pindex);
  task165->add_dep(task168);
  task168->add_dep(task108);
  residualq->add_task(task168);

  vector<IndexRange> I312_index = {closed_, active_, active_, active_, active_, active_};
  auto I312 = make_shared<Tensor>(I312_index);
  auto tensor169 = vector<shared_ptr<Tensor>>{I9, Gamma101_(), I312};
  auto task169 = make_shared<Task169>(tensor169, pindex);
  task139->add_dep(task169);
  task169->add_dep(task108);
  residualq->add_task(task169);

  auto tensor170 = vector<shared_ptr<Tensor>>{I312, t2, v2_};
  auto task170 = make_shared<Task170>(tensor170, pindex);
  task169->add_dep(task170);
  task170->add_dep(task108);
  residualq->add_task(task170);

  vector<IndexRange> I315_index = {active_, closed_, active_, active_, active_, active_};
  auto I315 = make_shared<Tensor>(I315_index);
  auto tensor171 = vector<shared_ptr<Tensor>>{I9, Gamma102_(), I315};
  auto task171 = make_shared<Task171>(tensor171, pindex);
  task139->add_dep(task171);
  task171->add_dep(task108);
  residualq->add_task(task171);

  auto tensor172 = vector<shared_ptr<Tensor>>{I315, t2, v2_};
  auto task172 = make_shared<Task172>(tensor172, pindex);
  task171->add_dep(task172);
  task172->add_dep(task108);
  residualq->add_task(task172);

  vector<IndexRange> I321_index = {active_, active_, closed_, active_};
  auto I321 = make_shared<Tensor>(I321_index);
  auto tensor173 = vector<shared_ptr<Tensor>>{I9, Gamma104_(), I321};
  auto task173 = make_shared<Task173>(tensor173, pindex);
  task139->add_dep(task173);
  task173->add_dep(task108);
  residualq->add_task(task173);

  vector<IndexRange> I322_index = {virt_, closed_, active_, active_};
  auto I322 = make_shared<Tensor>(I322_index);
  auto tensor174 = vector<shared_ptr<Tensor>>{I321, t2, I322};
  auto task174 = make_shared<Task174>(tensor174, pindex);
  task173->add_dep(task174);
  task174->add_dep(task108);
  residualq->add_task(task174);

  auto tensor175 = vector<shared_ptr<Tensor>>{I322, v2_};
  auto task175 = make_shared<Task175>(tensor175, pindex);
  task174->add_dep(task175);
  task175->add_dep(task108);
  residualq->add_task(task175);

  vector<IndexRange> I325_index = {virt_, closed_, active_, active_};
  auto I325 = make_shared<Tensor>(I325_index);
  auto tensor176 = vector<shared_ptr<Tensor>>{I321, t2, I325};
  auto task176 = make_shared<Task176>(tensor176, pindex);
  task173->add_dep(task176);
  task176->add_dep(task108);
  residualq->add_task(task176);

  auto tensor177 = vector<shared_ptr<Tensor>>{I325, v2_};
  auto task177 = make_shared<Task177>(tensor177, pindex);
  task176->add_dep(task177);
  task177->add_dep(task108);
  residualq->add_task(task177);

  auto tensor178 = vector<shared_ptr<Tensor>>{I321, t2, v2_};
  auto task178 = make_shared<Task178>(tensor178, pindex);
  task173->add_dep(task178);
  task178->add_dep(task108);
  residualq->add_task(task178);

  vector<IndexRange> I330_index = {active_, active_, closed_, active_};
  auto I330 = make_shared<Tensor>(I330_index);
  auto tensor179 = vector<shared_ptr<Tensor>>{I9, Gamma107_(), I330};
  auto task179 = make_shared<Task179>(tensor179, pindex);
  task139->add_dep(task179);
  task179->add_dep(task108);
  residualq->add_task(task179);

  auto tensor180 = vector<shared_ptr<Tensor>>{I330, t2, v2_};
  auto task180 = make_shared<Task180>(tensor180, pindex);
  task179->add_dep(task180);
  task180->add_dep(task108);
  residualq->add_task(task180);

  vector<IndexRange> I336_index = {active_, active_, closed_, active_};
  auto I336 = make_shared<Tensor>(I336_index);
  auto tensor181 = vector<shared_ptr<Tensor>>{I9, Gamma109_(), I336};
  auto task181 = make_shared<Task181>(tensor181, pindex);
  task139->add_dep(task181);
  task181->add_dep(task108);
  residualq->add_task(task181);

  auto tensor182 = vector<shared_ptr<Tensor>>{I336, t2, v2_};
  auto task182 = make_shared<Task182>(tensor182, pindex);
  task181->add_dep(task182);
  task182->add_dep(task108);
  residualq->add_task(task182);

  vector<IndexRange> I351_index = {active_, active_, active_, active_, closed_, active_};
  auto I351 = make_shared<Tensor>(I351_index);
  auto tensor183 = vector<shared_ptr<Tensor>>{I9, Gamma114_(), I351};
  auto task183 = make_shared<Task183>(tensor183, pindex);
  task139->add_dep(task183);
  task183->add_dep(task108);
  residualq->add_task(task183);

  auto tensor184 = vector<shared_ptr<Tensor>>{I351, t2, v2_};
  auto task184 = make_shared<Task184>(tensor184, pindex);
  task183->add_dep(task184);
  task184->add_dep(task108);
  residualq->add_task(task184);

  vector<IndexRange> I354_index = {active_, active_, active_, active_, closed_, active_};
  auto I354 = make_shared<Tensor>(I354_index);
  auto tensor185 = vector<shared_ptr<Tensor>>{I9, Gamma115_(), I354};
  auto task185 = make_shared<Task185>(tensor185, pindex);
  task139->add_dep(task185);
  task185->add_dep(task108);
  residualq->add_task(task185);

  auto tensor186 = vector<shared_ptr<Tensor>>{I354, t2, v2_};
  auto task186 = make_shared<Task186>(tensor186, pindex);
  task185->add_dep(task186);
  task186->add_dep(task108);
  residualq->add_task(task186);

  vector<IndexRange> I366_index = {active_, active_, active_, closed_, active_, active_};
  auto I366 = make_shared<Tensor>(I366_index);
  auto tensor187 = vector<shared_ptr<Tensor>>{I9, Gamma119_(), I366};
  auto task187 = make_shared<Task187>(tensor187, pindex);
  task139->add_dep(task187);
  task187->add_dep(task108);
  residualq->add_task(task187);

  auto tensor188 = vector<shared_ptr<Tensor>>{I366, t2, v2_};
  auto task188 = make_shared<Task188>(tensor188, pindex);
  task187->add_dep(task188);
  task188->add_dep(task108);
  residualq->add_task(task188);

  vector<IndexRange> I375_index = {closed_, active_, active_, active_, active_, active_};
  auto I375 = make_shared<Tensor>(I375_index);
  auto tensor189 = vector<shared_ptr<Tensor>>{I9, Gamma122_(), I375};
  auto task189 = make_shared<Task189>(tensor189, pindex);
  task139->add_dep(task189);
  task189->add_dep(task108);
  residualq->add_task(task189);

  auto tensor190 = vector<shared_ptr<Tensor>>{I375, t2, v2_};
  auto task190 = make_shared<Task190>(tensor190, pindex);
  task189->add_dep(task190);
  task190->add_dep(task108);
  residualq->add_task(task190);

  auto tensor191 = vector<shared_ptr<Tensor>>{I9, Gamma550_(), t2};
  auto task191 = make_shared<Task191>(tensor191, pindex);
  task139->add_dep(task191);
  task191->add_dep(task108);
  residualq->add_task(task191);

  auto tensor192 = vector<shared_ptr<Tensor>>{I9, Gamma551_(), t2};
  auto task192 = make_shared<Task192>(tensor192, pindex);
  task139->add_dep(task192);
  task192->add_dep(task108);
  residualq->add_task(task192);
}

#endif
