//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_residualq1.cc
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
#include <src/smith/relmrci/RelMRCI_tasks2.h>
#include <src/smith/relmrci/RelMRCI_tasks3.h>
#include <src/smith/relmrci/RelMRCI_tasks4.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue>RelMRCI::RelMRCI::make_residualq(const bool reset, const bool diagonal) {
  auto residualq = make_shared<Queue>();
  auto tensor83 = vector<shared_ptr<Tensor>>{r};
  auto task83 = make_shared<Task83>(tensor83, reset);
  residualq->add_task(task83);

  make_residualq1(residualq, task83, diagonal);
  make_residualq2(residualq, task83, diagonal);
  make_residualq3(residualq, task83, diagonal);
  make_residualq4(residualq, task83, diagonal);
  make_residualq5(residualq, task83, diagonal);
  make_residualq6(residualq, task83, diagonal);

  return residualq;
}

void RelMRCI::RelMRCI::make_residualq1(shared_ptr<Queue> residualq, shared_ptr<Task> task83, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I0_index = {closed_, closed_, active_, active_};
  auto I0 = make_shared<Tensor>(I0_index);
  auto tensor84 = vector<shared_ptr<Tensor>>{r, I0};
  auto task84 = make_shared<Task84>(tensor84, pindex);
  task84->add_dep(task83);
  residualq->add_task(task84);

  vector<IndexRange> I1_index = {closed_, closed_, active_, active_};
  auto I1 = make_shared<Tensor>(I1_index);
  auto tensor85 = vector<shared_ptr<Tensor>>{I0, Gamma0_(), I1};
  auto task85 = make_shared<Task85>(tensor85, pindex);
  task84->add_dep(task85);
  task85->add_dep(task83);
  residualq->add_task(task85);

  auto tensor86 = vector<shared_ptr<Tensor>>{I1, t2, h1_};
  auto task86 = make_shared<Task86>(tensor86, pindex);
  task85->add_dep(task86);
  task86->add_dep(task83);
  residualq->add_task(task86);

  vector<IndexRange> I214_index = {virt_, active_, closed_, closed_};
  auto I214 = make_shared<Tensor>(I214_index);
  auto tensor87 = vector<shared_ptr<Tensor>>{I1, t2, I214};
  auto task87 = make_shared<Task87>(tensor87, pindex);
  task85->add_dep(task87);
  task87->add_dep(task83);
  residualq->add_task(task87);

  auto tensor88 = vector<shared_ptr<Tensor>>{I214, v2_};
  auto task88 = make_shared<Task88>(tensor88, pindex);
  task87->add_dep(task88);
  task88->add_dep(task83);
  residualq->add_task(task88);

  auto tensor89 = vector<shared_ptr<Tensor>>{I1, t2, v2_};
  auto task89 = make_shared<Task89>(tensor89, pindex);
  task85->add_dep(task89);
  task89->add_dep(task83);
  residualq->add_task(task89);

  vector<IndexRange> I4_index = {active_, active_, active_, closed_};
  auto I4 = make_shared<Tensor>(I4_index);
  auto tensor90 = vector<shared_ptr<Tensor>>{I0, h1_, I4};
  auto task90 = make_shared<Task90>(tensor90, pindex);
  task84->add_dep(task90);
  task90->add_dep(task83);
  residualq->add_task(task90);

  auto tensor91 = vector<shared_ptr<Tensor>>{I4, t2, Gamma1_()};
  auto task91 = make_shared<Task91>(tensor91, pindex);
  task90->add_dep(task91);
  task91->add_dep(task83);
  residualq->add_task(task91);

  vector<IndexRange> I7_index = {closed_, closed_, active_, active_};
  auto I7 = make_shared<Tensor>(I7_index);
  auto tensor92 = vector<shared_ptr<Tensor>>{I0, Gamma2_(), I7};
  auto task92 = make_shared<Task92>(tensor92, pindex);
  task84->add_dep(task92);
  task92->add_dep(task83);
  residualq->add_task(task92);

  auto tensor93 = vector<shared_ptr<Tensor>>{I7, h1_, t2};
  auto task93 = make_shared<Task93>(tensor93, pindex);
  task92->add_dep(task93);
  task93->add_dep(task83);
  residualq->add_task(task93);

  auto tensor94 = vector<shared_ptr<Tensor>>{I7, t2, v2_};
  auto task94 = make_shared<Task94>(tensor94, pindex);
  task92->add_dep(task94);
  task94->add_dep(task83);
  residualq->add_task(task94);

  vector<IndexRange> I183_index = {closed_, active_, active_, closed_, active_, active_};
  auto I183 = make_shared<Tensor>(I183_index);
  auto tensor95 = vector<shared_ptr<Tensor>>{I0, Gamma58_(), I183};
  auto task95 = make_shared<Task95>(tensor95, pindex);
  task84->add_dep(task95);
  task95->add_dep(task83);
  residualq->add_task(task95);

  vector<IndexRange> I184_index = {closed_, closed_, active_, active_};
  auto I184 = make_shared<Tensor>(I184_index);
  auto tensor96 = vector<shared_ptr<Tensor>>{I183, t2, I184};
  auto task96 = make_shared<Task96>(tensor96, pindex);
  task95->add_dep(task96);
  task96->add_dep(task83);
  residualq->add_task(task96);

  auto tensor97 = vector<shared_ptr<Tensor>>{I184, v2_};
  auto task97 = make_shared<Task97>(tensor97, pindex);
  task96->add_dep(task97);
  task97->add_dep(task83);
  residualq->add_task(task97);

  vector<IndexRange> I186_index = {closed_, active_, active_, closed_, active_, active_};
  auto I186 = make_shared<Tensor>(I186_index);
  auto tensor98 = vector<shared_ptr<Tensor>>{I0, Gamma59_(), I186};
  auto task98 = make_shared<Task98>(tensor98, pindex);
  task84->add_dep(task98);
  task98->add_dep(task83);
  residualq->add_task(task98);

  auto tensor99 = vector<shared_ptr<Tensor>>{I186, t2, v2_};
  auto task99 = make_shared<Task99>(tensor99, pindex);
  task98->add_dep(task99);
  task99->add_dep(task83);
  residualq->add_task(task99);

  vector<IndexRange> I189_index = {active_, closed_, active_, closed_, active_, active_};
  auto I189 = make_shared<Tensor>(I189_index);
  auto tensor100 = vector<shared_ptr<Tensor>>{I0, Gamma60_(), I189};
  auto task100 = make_shared<Task100>(tensor100, pindex);
  task84->add_dep(task100);
  task100->add_dep(task83);
  residualq->add_task(task100);

  auto tensor101 = vector<shared_ptr<Tensor>>{I189, t2, v2_};
  auto task101 = make_shared<Task101>(tensor101, pindex);
  task100->add_dep(task101);
  task101->add_dep(task83);
  residualq->add_task(task101);

  vector<IndexRange> I198_index = {closed_, active_, active_, active_, active_, active_};
  auto I198 = make_shared<Tensor>(I198_index);
  auto tensor102 = vector<shared_ptr<Tensor>>{I0, t2, I198};
  auto task102 = make_shared<Task102>(tensor102, pindex);
  task84->add_dep(task102);
  task102->add_dep(task83);
  residualq->add_task(task102);

  auto tensor103 = vector<shared_ptr<Tensor>>{I198, Gamma63_(), v2_};
  auto task103 = make_shared<Task103>(tensor103, pindex);
  task102->add_dep(task103);
  task103->add_dep(task83);
  residualq->add_task(task103);

  auto tensor104 = vector<shared_ptr<Tensor>>{I198, Gamma64_(), v2_};
  auto task104 = make_shared<Task104>(tensor104, pindex);
  task102->add_dep(task104);
  task104->add_dep(task83);
  residualq->add_task(task104);

  vector<IndexRange> I204_index = {closed_, active_, active_, active_};
  auto I204 = make_shared<Tensor>(I204_index);
  auto tensor105 = vector<shared_ptr<Tensor>>{I0, v2_, I204};
  auto task105 = make_shared<Task105>(tensor105, pindex);
  task84->add_dep(task105);
  task105->add_dep(task83);
  residualq->add_task(task105);

  auto tensor106 = vector<shared_ptr<Tensor>>{I204, Gamma65_(), t2};
  auto task106 = make_shared<Task106>(tensor106, pindex);
  task105->add_dep(task106);
  task106->add_dep(task83);
  residualq->add_task(task106);

  vector<IndexRange> I207_index = {virt_, active_, active_, active_};
  auto I207 = make_shared<Tensor>(I207_index);
  auto tensor107 = vector<shared_ptr<Tensor>>{I0, t2, I207};
  auto task107 = make_shared<Task107>(tensor107, pindex);
  task84->add_dep(task107);
  task107->add_dep(task83);
  residualq->add_task(task107);

  auto tensor108 = vector<shared_ptr<Tensor>>{I207, Gamma66_(), v2_};
  auto task108 = make_shared<Task108>(tensor108, pindex);
  task107->add_dep(task108);
  task108->add_dep(task83);
  residualq->add_task(task108);

  auto tensor109 = vector<shared_ptr<Tensor>>{I207, Gamma67_(), v2_};
  auto task109 = make_shared<Task109>(tensor109, pindex);
  task107->add_dep(task109);
  task109->add_dep(task83);
  residualq->add_task(task109);

  vector<IndexRange> I225_index = {closed_, active_, active_, closed_, active_, active_};
  auto I225 = make_shared<Tensor>(I225_index);
  auto tensor110 = vector<shared_ptr<Tensor>>{I0, Gamma65_(), I225};
  auto task110 = make_shared<Task110>(tensor110, pindex);
  task84->add_dep(task110);
  task110->add_dep(task83);
  residualq->add_task(task110);

  auto tensor111 = vector<shared_ptr<Tensor>>{I225, t2, v2_};
  auto task111 = make_shared<Task111>(tensor111, pindex);
  task110->add_dep(task111);
  task111->add_dep(task83);
  residualq->add_task(task111);

  vector<IndexRange> I9_index = {closed_, active_, active_, active_};
  auto I9 = make_shared<Tensor>(I9_index);
  auto tensor112 = vector<shared_ptr<Tensor>>{r, I9};
  auto task112 = make_shared<Task112>(tensor112, pindex);
  task112->add_dep(task83);
  residualq->add_task(task112);

  vector<IndexRange> I10_index = {active_, closed_, active_, active_};
  auto I10 = make_shared<Tensor>(I10_index);
  auto tensor113 = vector<shared_ptr<Tensor>>{I9, Gamma3_(), I10};
  auto task113 = make_shared<Task113>(tensor113, pindex);
  task112->add_dep(task113);
  task113->add_dep(task83);
  residualq->add_task(task113);

  auto tensor114 = vector<shared_ptr<Tensor>>{I10, t2, h1_};
  auto task114 = make_shared<Task114>(tensor114, pindex);
  task113->add_dep(task114);
  task114->add_dep(task83);
  residualq->add_task(task114);

  auto tensor115 = vector<shared_ptr<Tensor>>{I10, t2, v2_};
  auto task115 = make_shared<Task115>(tensor115, pindex);
  task113->add_dep(task115);
  task115->add_dep(task83);
  residualq->add_task(task115);

  auto tensor116 = vector<shared_ptr<Tensor>>{I10, t2, v2_};
  auto task116 = make_shared<Task116>(tensor116, pindex);
  task113->add_dep(task116);
  task116->add_dep(task83);
  residualq->add_task(task116);

  vector<IndexRange> I13_index = {closed_, active_, active_, active_};
  auto I13 = make_shared<Tensor>(I13_index);
  auto tensor117 = vector<shared_ptr<Tensor>>{I9, Gamma4_(), I13};
  auto task117 = make_shared<Task117>(tensor117, pindex);
  task112->add_dep(task117);
  task117->add_dep(task83);
  residualq->add_task(task117);

  auto tensor118 = vector<shared_ptr<Tensor>>{I13, t2, h1_};
  auto task118 = make_shared<Task118>(tensor118, pindex);
  task117->add_dep(task118);
  task118->add_dep(task83);
  residualq->add_task(task118);

  auto tensor119 = vector<shared_ptr<Tensor>>{I13, t2, h1_};
  auto task119 = make_shared<Task119>(tensor119, pindex);
  task117->add_dep(task119);
  task119->add_dep(task83);
  residualq->add_task(task119);

  auto tensor120 = vector<shared_ptr<Tensor>>{I13, v2_, t2};
  auto task120 = make_shared<Task120>(tensor120, pindex);
  task117->add_dep(task120);
  task120->add_dep(task83);
  residualq->add_task(task120);

  auto tensor121 = vector<shared_ptr<Tensor>>{I13, t2, v2_};
  auto task121 = make_shared<Task121>(tensor121, pindex);
  task117->add_dep(task121);
  task121->add_dep(task83);
  residualq->add_task(task121);

  vector<IndexRange> I16_index = {closed_, active_};
  auto I16 = make_shared<Tensor>(I16_index);
  auto tensor122 = vector<shared_ptr<Tensor>>{I9, Gamma5_(), I16};
  auto task122 = make_shared<Task122>(tensor122, pindex);
  task112->add_dep(task122);
  task122->add_dep(task83);
  residualq->add_task(task122);

  auto tensor123 = vector<shared_ptr<Tensor>>{I16, t2, h1_};
  auto task123 = make_shared<Task123>(tensor123, pindex);
  task122->add_dep(task123);
  task123->add_dep(task83);
  residualq->add_task(task123);

  auto tensor124 = vector<shared_ptr<Tensor>>{I16, t2, h1_};
  auto task124 = make_shared<Task124>(tensor124, pindex);
  task122->add_dep(task124);
  task124->add_dep(task83);
  residualq->add_task(task124);

  auto tensor125 = vector<shared_ptr<Tensor>>{I16, t2, v2_};
  auto task125 = make_shared<Task125>(tensor125, pindex);
  task122->add_dep(task125);
  task125->add_dep(task83);
  residualq->add_task(task125);

  auto tensor126 = vector<shared_ptr<Tensor>>{I16, t2, v2_};
  auto task126 = make_shared<Task126>(tensor126, pindex);
  task122->add_dep(task126);
  task126->add_dep(task83);
  residualq->add_task(task126);

  auto tensor127 = vector<shared_ptr<Tensor>>{I16, t2, v2_};
  auto task127 = make_shared<Task127>(tensor127, pindex);
  task122->add_dep(task127);
  task127->add_dep(task83);
  residualq->add_task(task127);

  auto tensor128 = vector<shared_ptr<Tensor>>{I16, t2, v2_};
  auto task128 = make_shared<Task128>(tensor128, pindex);
  task122->add_dep(task128);
  task128->add_dep(task83);
  residualq->add_task(task128);

  vector<IndexRange> I231_index = {active_, active_, active_, closed_, active_, active_};
  auto I231 = make_shared<Tensor>(I231_index);
  auto tensor129 = vector<shared_ptr<Tensor>>{I9, Gamma74_(), I231};
  auto task129 = make_shared<Task129>(tensor129, pindex);
  task112->add_dep(task129);
  task129->add_dep(task83);
  residualq->add_task(task129);

  auto tensor130 = vector<shared_ptr<Tensor>>{I231, t2, v2_};
  auto task130 = make_shared<Task130>(tensor130, pindex);
  task129->add_dep(task130);
  task130->add_dep(task83);
  residualq->add_task(task130);

  vector<IndexRange> I234_index = {active_, active_, active_, closed_, active_, active_};
  auto I234 = make_shared<Tensor>(I234_index);
  auto tensor131 = vector<shared_ptr<Tensor>>{I9, Gamma75_(), I234};
  auto task131 = make_shared<Task131>(tensor131, pindex);
  task112->add_dep(task131);
  task131->add_dep(task83);
  residualq->add_task(task131);

  auto tensor132 = vector<shared_ptr<Tensor>>{I234, t2, v2_};
  auto task132 = make_shared<Task132>(tensor132, pindex);
  task131->add_dep(task132);
  task132->add_dep(task83);
  residualq->add_task(task132);

  vector<IndexRange> I240_index = {closed_, active_, active_, active_, active_, active_};
  auto I240 = make_shared<Tensor>(I240_index);
  auto tensor133 = vector<shared_ptr<Tensor>>{I9, Gamma77_(), I240};
  auto task133 = make_shared<Task133>(tensor133, pindex);
  task112->add_dep(task133);
  task133->add_dep(task83);
  residualq->add_task(task133);

  vector<IndexRange> I241_index = {closed_, closed_, active_, active_};
  auto I241 = make_shared<Tensor>(I241_index);
  auto tensor134 = vector<shared_ptr<Tensor>>{I240, t2, I241};
  auto task134 = make_shared<Task134>(tensor134, pindex);
  task133->add_dep(task134);
  task134->add_dep(task83);
  residualq->add_task(task134);

  auto tensor135 = vector<shared_ptr<Tensor>>{I241, v2_};
  auto task135 = make_shared<Task135>(tensor135, pindex);
  task134->add_dep(task135);
  task135->add_dep(task83);
  residualq->add_task(task135);

  auto tensor136 = vector<shared_ptr<Tensor>>{I240, t2, v2_};
  auto task136 = make_shared<Task136>(tensor136, pindex);
  task133->add_dep(task136);
  task136->add_dep(task83);
  residualq->add_task(task136);

  vector<IndexRange> I243_index = {closed_, active_, active_, active_, active_, active_};
  auto I243 = make_shared<Tensor>(I243_index);
  auto tensor137 = vector<shared_ptr<Tensor>>{I9, Gamma78_(), I243};
  auto task137 = make_shared<Task137>(tensor137, pindex);
  task112->add_dep(task137);
  task137->add_dep(task83);
  residualq->add_task(task137);

  auto tensor138 = vector<shared_ptr<Tensor>>{I243, t2, v2_};
  auto task138 = make_shared<Task138>(tensor138, pindex);
  task137->add_dep(task138);
  task138->add_dep(task83);
  residualq->add_task(task138);

  vector<IndexRange> I246_index = {active_, closed_, active_, active_, active_, active_};
  auto I246 = make_shared<Tensor>(I246_index);
  auto tensor139 = vector<shared_ptr<Tensor>>{I9, Gamma79_(), I246};
  auto task139 = make_shared<Task139>(tensor139, pindex);
  task112->add_dep(task139);
  task139->add_dep(task83);
  residualq->add_task(task139);

  auto tensor140 = vector<shared_ptr<Tensor>>{I246, t2, v2_};
  auto task140 = make_shared<Task140>(tensor140, pindex);
  task139->add_dep(task140);
  task140->add_dep(task83);
  residualq->add_task(task140);

  vector<IndexRange> I252_index = {active_, active_, closed_, active_};
  auto I252 = make_shared<Tensor>(I252_index);
  auto tensor141 = vector<shared_ptr<Tensor>>{I9, Gamma81_(), I252};
  auto task141 = make_shared<Task141>(tensor141, pindex);
  task112->add_dep(task141);
  task141->add_dep(task83);
  residualq->add_task(task141);

  vector<IndexRange> I253_index = {virt_, closed_, active_, active_};
  auto I253 = make_shared<Tensor>(I253_index);
  auto tensor142 = vector<shared_ptr<Tensor>>{I252, t2, I253};
  auto task142 = make_shared<Task142>(tensor142, pindex);
  task141->add_dep(task142);
  task142->add_dep(task83);
  residualq->add_task(task142);

  auto tensor143 = vector<shared_ptr<Tensor>>{I253, v2_};
  auto task143 = make_shared<Task143>(tensor143, pindex);
  task142->add_dep(task143);
  task143->add_dep(task83);
  residualq->add_task(task143);

  vector<IndexRange> I256_index = {virt_, closed_, active_, active_};
  auto I256 = make_shared<Tensor>(I256_index);
  auto tensor144 = vector<shared_ptr<Tensor>>{I252, t2, I256};
  auto task144 = make_shared<Task144>(tensor144, pindex);
  task141->add_dep(task144);
  task144->add_dep(task83);
  residualq->add_task(task144);

  auto tensor145 = vector<shared_ptr<Tensor>>{I256, v2_};
  auto task145 = make_shared<Task145>(tensor145, pindex);
  task144->add_dep(task145);
  task145->add_dep(task83);
  residualq->add_task(task145);

  auto tensor146 = vector<shared_ptr<Tensor>>{I252, t2, v2_};
  auto task146 = make_shared<Task146>(tensor146, pindex);
  task141->add_dep(task146);
  task146->add_dep(task83);
  residualq->add_task(task146);

  vector<IndexRange> I261_index = {active_, active_, closed_, active_};
  auto I261 = make_shared<Tensor>(I261_index);
  auto tensor147 = vector<shared_ptr<Tensor>>{I9, Gamma84_(), I261};
  auto task147 = make_shared<Task147>(tensor147, pindex);
  task112->add_dep(task147);
  task147->add_dep(task83);
  residualq->add_task(task147);

  auto tensor148 = vector<shared_ptr<Tensor>>{I261, t2, v2_};
  auto task148 = make_shared<Task148>(tensor148, pindex);
  task147->add_dep(task148);
  task148->add_dep(task83);
  residualq->add_task(task148);

  vector<IndexRange> I267_index = {active_, active_, closed_, active_};
  auto I267 = make_shared<Tensor>(I267_index);
  auto tensor149 = vector<shared_ptr<Tensor>>{I9, Gamma86_(), I267};
  auto task149 = make_shared<Task149>(tensor149, pindex);
  task112->add_dep(task149);
  task149->add_dep(task83);
  residualq->add_task(task149);

  auto tensor150 = vector<shared_ptr<Tensor>>{I267, t2, v2_};
  auto task150 = make_shared<Task150>(tensor150, pindex);
  task149->add_dep(task150);
  task150->add_dep(task83);
  residualq->add_task(task150);

  vector<IndexRange> I285_index = {active_, active_, active_, closed_, active_, active_};
  auto I285 = make_shared<Tensor>(I285_index);
  auto tensor151 = vector<shared_ptr<Tensor>>{I9, Gamma92_(), I285};
  auto task151 = make_shared<Task151>(tensor151, pindex);
  task112->add_dep(task151);
  task151->add_dep(task83);
  residualq->add_task(task151);

  auto tensor152 = vector<shared_ptr<Tensor>>{I285, t2, v2_};
  auto task152 = make_shared<Task152>(tensor152, pindex);
  task151->add_dep(task152);
  task152->add_dep(task83);
  residualq->add_task(task152);

  vector<IndexRange> I294_index = {closed_, active_, active_, active_, active_, active_};
  auto I294 = make_shared<Tensor>(I294_index);
  auto tensor153 = vector<shared_ptr<Tensor>>{I9, Gamma95_(), I294};
  auto task153 = make_shared<Task153>(tensor153, pindex);
  task112->add_dep(task153);
  task153->add_dep(task83);
  residualq->add_task(task153);

  auto tensor154 = vector<shared_ptr<Tensor>>{I294, t2, v2_};
  auto task154 = make_shared<Task154>(tensor154, pindex);
  task153->add_dep(task154);
  task154->add_dep(task83);
  residualq->add_task(task154);

  vector<IndexRange> I303_index = {active_, active_, active_, closed_};
  auto I303 = make_shared<Tensor>(I303_index);
  auto tensor155 = vector<shared_ptr<Tensor>>{I9, Gamma98_(), I303};
  auto task155 = make_shared<Task155>(tensor155, pindex);
  task112->add_dep(task155);
  task155->add_dep(task83);
  residualq->add_task(task155);

  auto tensor156 = vector<shared_ptr<Tensor>>{I303, t2, v2_};
  auto task156 = make_shared<Task156>(tensor156, pindex);
  task155->add_dep(task156);
  task156->add_dep(task83);
  residualq->add_task(task156);

  auto tensor157 = vector<shared_ptr<Tensor>>{I9, Gamma414_(), t2};
  auto task157 = make_shared<Task157>(tensor157, pindex);
  task112->add_dep(task157);
  task157->add_dep(task83);
  residualq->add_task(task157);

  auto tensor158 = vector<shared_ptr<Tensor>>{I9, Gamma415_(), t2};
  auto task158 = make_shared<Task158>(tensor158, pindex);
  task112->add_dep(task158);
  task158->add_dep(task83);
  residualq->add_task(task158);
}

#endif
