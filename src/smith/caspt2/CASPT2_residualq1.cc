//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_residualq1.cc
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
#include <src/smith/caspt2/CASPT2_tasks2.h>
#include <src/smith/caspt2/CASPT2_tasks3.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_residualq(const bool reset, const bool diagonal) {
  auto residualq = make_shared<Queue>();
  auto tensor69 = vector<shared_ptr<Tensor>>{r};
  auto task69 = make_shared<Task69>(tensor69, reset);
  residualq->add_task(task69);

  make_residualq1(residualq, task69, diagonal);
  make_residualq2(residualq, task69, diagonal);
  make_residualq3(residualq, task69, diagonal);

  return residualq;
}


void CASPT2::CASPT2::make_residualq1(shared_ptr<Queue> residualq, shared_ptr<Task> task69, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I0_index = {closed_, closed_, active_, active_};
  auto I0 = make_shared<Tensor>(I0_index);
  auto tensor70 = vector<shared_ptr<Tensor>>{r, I0};
  auto task70 = make_shared<Task70>(tensor70, pindex);
  task70->add_dep(task69);
  residualq->add_task(task70);

  auto tensor71 = vector<shared_ptr<Tensor>>{I0, Gamma0_(), t2};
  auto task71 = make_shared<Task71>(tensor71, pindex);
  task70->add_dep(task71);
  task71->add_dep(task69);
  residualq->add_task(task71);

  auto tensor72 = vector<shared_ptr<Tensor>>{I0, Gamma92_(), t2};
  auto task72 = make_shared<Task72>(tensor72, pindex, this->e0_);
  task70->add_dep(task72);
  task72->add_dep(task69);
  residualq->add_task(task72);

  vector<IndexRange> I2_index = {closed_, closed_, active_, active_};
  auto I2 = make_shared<Tensor>(I2_index);
  auto tensor73 = vector<shared_ptr<Tensor>>{r, I2};
  auto task73 = make_shared<Task73>(tensor73, pindex);
  task73->add_dep(task69);
  residualq->add_task(task73);

  vector<IndexRange> I3_index = {closed_, closed_, active_, active_};
  auto I3 = make_shared<Tensor>(I3_index);
  auto tensor74 = vector<shared_ptr<Tensor>>{I2, Gamma92_(), I3};
  auto task74 = make_shared<Task74>(tensor74, pindex);
  task73->add_dep(task74);
  task74->add_dep(task69);
  residualq->add_task(task74);

  auto tensor75 = vector<shared_ptr<Tensor>>{I3, t2, f1_};
  auto task75 = make_shared<Task75>(tensor75, pindex);
  task74->add_dep(task75);
  task75->add_dep(task69);
  residualq->add_task(task75);

  vector<IndexRange> I6_index = {closed_, active_, active_, active_};
  auto I6 = make_shared<Tensor>(I6_index);
  auto tensor76 = vector<shared_ptr<Tensor>>{I2, f1_, I6};
  auto task76 = make_shared<Task76>(tensor76, pindex);
  task73->add_dep(task76);
  task76->add_dep(task69);
  residualq->add_task(task76);

  auto tensor77 = vector<shared_ptr<Tensor>>{I6, Gamma2_(), t2};
  auto task77 = make_shared<Task77>(tensor77, pindex);
  task76->add_dep(task77);
  task77->add_dep(task69);
  residualq->add_task(task77);

  vector<IndexRange> I9_index = {active_, closed_, closed_, active_};
  auto I9 = make_shared<Tensor>(I9_index);
  auto tensor78 = vector<shared_ptr<Tensor>>{I2, Gamma3_(), I9};
  auto task78 = make_shared<Task78>(tensor78, pindex);
  task73->add_dep(task78);
  task78->add_dep(task69);
  residualq->add_task(task78);

  auto tensor79 = vector<shared_ptr<Tensor>>{I9, t2, f1_};
  auto task79 = make_shared<Task79>(tensor79, pindex);
  task78->add_dep(task79);
  task79->add_dep(task69);
  residualq->add_task(task79);

  vector<IndexRange> I11_index = {closed_, active_, active_, active_};
  auto I11 = make_shared<Tensor>(I11_index);
  auto tensor80 = vector<shared_ptr<Tensor>>{r, I11};
  auto task80 = make_shared<Task80>(tensor80, pindex);
  task80->add_dep(task69);
  residualq->add_task(task80);

  vector<IndexRange> I12_index = {active_, closed_, active_, active_};
  auto I12 = make_shared<Tensor>(I12_index);
  auto tensor81 = vector<shared_ptr<Tensor>>{I11, Gamma4_(), I12};
  auto task81 = make_shared<Task81>(tensor81, pindex);
  task80->add_dep(task81);
  task81->add_dep(task69);
  residualq->add_task(task81);

  auto tensor82 = vector<shared_ptr<Tensor>>{I12, t2, f1_};
  auto task82 = make_shared<Task82>(tensor82, pindex);
  task81->add_dep(task82);
  task82->add_dep(task69);
  residualq->add_task(task82);

  auto tensor83 = vector<shared_ptr<Tensor>>{I11, Gamma5_(), t2};
  auto task83 = make_shared<Task83>(tensor83, pindex);
  task80->add_dep(task83);
  task83->add_dep(task69);
  residualq->add_task(task83);

  vector<IndexRange> I17_index = {closed_, active_, active_, active_};
  auto I17 = make_shared<Tensor>(I17_index);
  auto tensor84 = vector<shared_ptr<Tensor>>{I11, Gamma6_(), I17};
  auto task84 = make_shared<Task84>(tensor84, pindex);
  task80->add_dep(task84);
  task84->add_dep(task69);
  residualq->add_task(task84);

  auto tensor85 = vector<shared_ptr<Tensor>>{I17, t2};
  auto task85 = make_shared<Task85>(tensor85, pindex, this->e0_);
  task84->add_dep(task85);
  task85->add_dep(task69);
  residualq->add_task(task85);

  auto tensor86 = vector<shared_ptr<Tensor>>{I17, t2, f1_};
  auto task86 = make_shared<Task86>(tensor86, pindex);
  task84->add_dep(task86);
  task86->add_dep(task69);
  residualq->add_task(task86);

  auto tensor87 = vector<shared_ptr<Tensor>>{I17, t2, f1_};
  auto task87 = make_shared<Task87>(tensor87, pindex);
  task84->add_dep(task87);
  task87->add_dep(task69);
  residualq->add_task(task87);

  vector<IndexRange> I20_index = {closed_, active_};
  auto I20 = make_shared<Tensor>(I20_index);
  auto tensor88 = vector<shared_ptr<Tensor>>{I11, Gamma7_(), I20};
  auto task88 = make_shared<Task88>(tensor88, pindex);
  task80->add_dep(task88);
  task88->add_dep(task69);
  residualq->add_task(task88);

  auto tensor89 = vector<shared_ptr<Tensor>>{I20, t2, f1_};
  auto task89 = make_shared<Task89>(tensor89, pindex);
  task88->add_dep(task89);
  task89->add_dep(task69);
  residualq->add_task(task89);

  auto tensor90 = vector<shared_ptr<Tensor>>{I20, t2, f1_};
  auto task90 = make_shared<Task90>(tensor90, pindex);
  task88->add_dep(task90);
  task90->add_dep(task69);
  residualq->add_task(task90);

  vector<IndexRange> I26_index = {active_, active_, closed_, active_};
  auto I26 = make_shared<Tensor>(I26_index);
  auto tensor91 = vector<shared_ptr<Tensor>>{I11, Gamma9_(), I26};
  auto task91 = make_shared<Task91>(tensor91, pindex);
  task80->add_dep(task91);
  task91->add_dep(task69);
  residualq->add_task(task91);

  auto tensor92 = vector<shared_ptr<Tensor>>{I26, t2, f1_};
  auto task92 = make_shared<Task92>(tensor92, pindex);
  task91->add_dep(task92);
  task92->add_dep(task69);
  residualq->add_task(task92);

  vector<IndexRange> I31_index = {closed_, closed_, active_, virt_};
  auto I31 = make_shared<Tensor>(I31_index);
  auto tensor93 = vector<shared_ptr<Tensor>>{r, I31};
  auto task93 = make_shared<Task93>(tensor93, pindex);
  task93->add_dep(task69);
  residualq->add_task(task93);

  vector<IndexRange> I32_index = {closed_, closed_, active_, active_};
  auto I32 = make_shared<Tensor>(I32_index);
  auto tensor94 = vector<shared_ptr<Tensor>>{I31, f1_, I32};
  auto task94 = make_shared<Task94>(tensor94, pindex);
  task93->add_dep(task94);
  task94->add_dep(task69);
  residualq->add_task(task94);

  auto tensor95 = vector<shared_ptr<Tensor>>{I32, Gamma3_(), t2};
  auto task95 = make_shared<Task95>(tensor95, pindex);
  task94->add_dep(task95);
  task95->add_dep(task69);
  residualq->add_task(task95);

  vector<IndexRange> I35_index = {closed_, active_};
  auto I35 = make_shared<Tensor>(I35_index);
  auto tensor96 = vector<shared_ptr<Tensor>>{I31, f1_, I35};
  auto task96 = make_shared<Task96>(tensor96, pindex);
  task93->add_dep(task96);
  task96->add_dep(task69);
  residualq->add_task(task96);

  auto tensor97 = vector<shared_ptr<Tensor>>{I35, Gamma12_(), t2};
  auto task97 = make_shared<Task97>(tensor97, pindex);
  task96->add_dep(task97);
  task97->add_dep(task69);
  residualq->add_task(task97);

  vector<IndexRange> I38_index = {closed_, active_};
  auto I38 = make_shared<Tensor>(I38_index);
  auto tensor98 = vector<shared_ptr<Tensor>>{I31, f1_, I38};
  auto task98 = make_shared<Task98>(tensor98, pindex);
  task93->add_dep(task98);
  task98->add_dep(task69);
  residualq->add_task(task98);

  auto tensor99 = vector<shared_ptr<Tensor>>{I38, Gamma12_(), t2};
  auto task99 = make_shared<Task99>(tensor99, pindex);
  task98->add_dep(task99);
  task99->add_dep(task69);
  residualq->add_task(task99);

  vector<IndexRange> I41_index = {closed_, virt_, closed_, active_};
  auto I41 = make_shared<Tensor>(I41_index);
  auto tensor100 = vector<shared_ptr<Tensor>>{I31, Gamma14_(), I41};
  auto task100 = make_shared<Task100>(tensor100, pindex);
  task93->add_dep(task100);
  task100->add_dep(task69);
  residualq->add_task(task100);

  auto tensor101 = vector<shared_ptr<Tensor>>{I41, t2};
  auto task101 = make_shared<Task101>(tensor101, pindex);
  task100->add_dep(task101);
  task101->add_dep(task69);
  residualq->add_task(task101);

  vector<IndexRange> I45_index = {closed_, virt_, closed_, active_};
  auto I45 = make_shared<Tensor>(I45_index);
  auto tensor102 = vector<shared_ptr<Tensor>>{I31, Gamma16_(), I45};
  auto task102 = make_shared<Task102>(tensor102, pindex);
  task93->add_dep(task102);
  task102->add_dep(task69);
  residualq->add_task(task102);

  auto tensor103 = vector<shared_ptr<Tensor>>{I45, t2};
  auto task103 = make_shared<Task103>(tensor103, pindex, this->e0_);
  task102->add_dep(task103);
  task103->add_dep(task69);
  residualq->add_task(task103);

  auto tensor104 = vector<shared_ptr<Tensor>>{I45, t2, f1_};
  auto task104 = make_shared<Task104>(tensor104, pindex);
  task102->add_dep(task104);
  task104->add_dep(task69);
  residualq->add_task(task104);

  auto tensor105 = vector<shared_ptr<Tensor>>{I45, t2, f1_};
  auto task105 = make_shared<Task105>(tensor105, pindex);
  task102->add_dep(task105);
  task105->add_dep(task69);
  residualq->add_task(task105);

  auto tensor106 = vector<shared_ptr<Tensor>>{I45, t2, f1_};
  auto task106 = make_shared<Task106>(tensor106, pindex);
  task102->add_dep(task106);
  task106->add_dep(task69);
  residualq->add_task(task106);

  auto tensor107 = vector<shared_ptr<Tensor>>{I45, t2, f1_};
  auto task107 = make_shared<Task107>(tensor107, pindex);
  task102->add_dep(task107);
  task107->add_dep(task69);
  residualq->add_task(task107);

  auto tensor108 = vector<shared_ptr<Tensor>>{I45, t2, f1_};
  auto task108 = make_shared<Task108>(tensor108, pindex);
  task102->add_dep(task108);
  task108->add_dep(task69);
  residualq->add_task(task108);

  auto tensor109 = vector<shared_ptr<Tensor>>{I45, t2, f1_};
  auto task109 = make_shared<Task109>(tensor109, pindex);
  task102->add_dep(task109);
  task109->add_dep(task69);
  residualq->add_task(task109);

  vector<IndexRange> I63_index = {virt_, closed_, active_, active_};
  auto I63 = make_shared<Tensor>(I63_index);
  auto tensor110 = vector<shared_ptr<Tensor>>{I31, f1_, I63};
  auto task110 = make_shared<Task110>(tensor110, pindex);
  task93->add_dep(task110);
  task110->add_dep(task69);
  residualq->add_task(task110);

  auto tensor111 = vector<shared_ptr<Tensor>>{I63, Gamma22_(), t2};
  auto task111 = make_shared<Task111>(tensor111, pindex);
  task110->add_dep(task111);
  task111->add_dep(task69);
  residualq->add_task(task111);

  auto tensor112 = vector<shared_ptr<Tensor>>{I63, Gamma12_(), t2};
  auto task112 = make_shared<Task112>(tensor112, pindex);
  task110->add_dep(task112);
  task112->add_dep(task69);
  residualq->add_task(task112);

  vector<IndexRange> I66_index = {virt_, closed_, active_, active_};
  auto I66 = make_shared<Tensor>(I66_index);
  auto tensor113 = vector<shared_ptr<Tensor>>{I31, f1_, I66};
  auto task113 = make_shared<Task113>(tensor113, pindex);
  task93->add_dep(task113);
  task113->add_dep(task69);
  residualq->add_task(task113);

  vector<IndexRange> I67_index = {active_, virt_, closed_, active_};
  auto I67 = make_shared<Tensor>(I67_index);
  auto tensor114 = vector<shared_ptr<Tensor>>{I66, Gamma12_(), I67};
  auto task114 = make_shared<Task114>(tensor114, pindex);
  task113->add_dep(task114);
  task114->add_dep(task69);
  residualq->add_task(task114);

  auto tensor115 = vector<shared_ptr<Tensor>>{I67, t2};
  auto task115 = make_shared<Task115>(tensor115, pindex);
  task114->add_dep(task115);
  task115->add_dep(task69);
  residualq->add_task(task115);

  vector<IndexRange> I75_index = {virt_, active_};
  auto I75 = make_shared<Tensor>(I75_index);
  auto tensor116 = vector<shared_ptr<Tensor>>{I31, t2, I75};
  auto task116 = make_shared<Task116>(tensor116, pindex);
  task93->add_dep(task116);
  task116->add_dep(task69);
  residualq->add_task(task116);

  auto tensor117 = vector<shared_ptr<Tensor>>{I75, Gamma16_(), f1_};
  auto task117 = make_shared<Task117>(tensor117, pindex);
  task116->add_dep(task117);
  task117->add_dep(task69);
  residualq->add_task(task117);

  vector<IndexRange> I78_index = {virt_, active_};
  auto I78 = make_shared<Tensor>(I78_index);
  auto tensor118 = vector<shared_ptr<Tensor>>{I31, t2, I78};
  auto task118 = make_shared<Task118>(tensor118, pindex);
  task93->add_dep(task118);
  task118->add_dep(task69);
  residualq->add_task(task118);

  auto tensor119 = vector<shared_ptr<Tensor>>{I78, Gamma16_(), f1_};
  auto task119 = make_shared<Task119>(tensor119, pindex);
  task118->add_dep(task119);
  task119->add_dep(task69);
  residualq->add_task(task119);

  vector<IndexRange> I80_index = {closed_, active_, active_, virt_};
  auto I80 = make_shared<Tensor>(I80_index);
  auto tensor120 = vector<shared_ptr<Tensor>>{r, I80};
  auto task120 = make_shared<Task120>(tensor120, pindex);
  task120->add_dep(task69);
  residualq->add_task(task120);

  vector<IndexRange> I81_index = {closed_, active_, active_, active_};
  auto I81 = make_shared<Tensor>(I81_index);
  auto tensor121 = vector<shared_ptr<Tensor>>{I80, f1_, I81};
  auto task121 = make_shared<Task121>(tensor121, pindex);
  task120->add_dep(task121);
  task121->add_dep(task69);
  residualq->add_task(task121);

  auto tensor122 = vector<shared_ptr<Tensor>>{I81, Gamma28_(), t2};
  auto task122 = make_shared<Task122>(tensor122, pindex);
  task121->add_dep(task122);
  task122->add_dep(task69);
  residualq->add_task(task122);

  vector<IndexRange> I84_index = {active_, virt_, closed_, active_};
  auto I84 = make_shared<Tensor>(I84_index);
  auto tensor123 = vector<shared_ptr<Tensor>>{I80, Gamma29_(), I84};
  auto task123 = make_shared<Task123>(tensor123, pindex);
  task120->add_dep(task123);
  task123->add_dep(task69);
  residualq->add_task(task123);

  auto tensor124 = vector<shared_ptr<Tensor>>{I84, t2, f1_};
  auto task124 = make_shared<Task124>(tensor124, pindex);
  task123->add_dep(task124);
  task124->add_dep(task69);
  residualq->add_task(task124);

  vector<IndexRange> I87_index = {active_, closed_, virt_, active_};
  auto I87 = make_shared<Tensor>(I87_index);
  auto tensor125 = vector<shared_ptr<Tensor>>{I80, Gamma7_(), I87};
  auto task125 = make_shared<Task125>(tensor125, pindex);
  task120->add_dep(task125);
  task125->add_dep(task69);
  residualq->add_task(task125);

  auto tensor126 = vector<shared_ptr<Tensor>>{I87, t2, f1_};
  auto task126 = make_shared<Task126>(tensor126, pindex);
  task125->add_dep(task126);
  task126->add_dep(task69);
  residualq->add_task(task126);

  auto tensor127 = vector<shared_ptr<Tensor>>{I80, Gamma31_(), t2};
  auto task127 = make_shared<Task127>(tensor127, pindex);
  task120->add_dep(task127);
  task127->add_dep(task69);
  residualq->add_task(task127);

  vector<IndexRange> I92_index = {closed_, active_, virt_, active_};
  auto I92 = make_shared<Tensor>(I92_index);
  auto tensor128 = vector<shared_ptr<Tensor>>{I80, Gamma32_(), I92};
  auto task128 = make_shared<Task128>(tensor128, pindex);
  task120->add_dep(task128);
  task128->add_dep(task69);
  residualq->add_task(task128);

  auto tensor129 = vector<shared_ptr<Tensor>>{I92, t2};
  auto task129 = make_shared<Task129>(tensor129, pindex, this->e0_);
  task128->add_dep(task129);
  task129->add_dep(task69);
  residualq->add_task(task129);

  auto tensor130 = vector<shared_ptr<Tensor>>{I92, t2, f1_};
  auto task130 = make_shared<Task130>(tensor130, pindex);
  task128->add_dep(task130);
  task130->add_dep(task69);
  residualq->add_task(task130);

  auto tensor131 = vector<shared_ptr<Tensor>>{I92, t2, f1_};
  auto task131 = make_shared<Task131>(tensor131, pindex);
  task128->add_dep(task131);
  task131->add_dep(task69);
  residualq->add_task(task131);

  auto tensor132 = vector<shared_ptr<Tensor>>{I92, t2, f1_};
  auto task132 = make_shared<Task132>(tensor132, pindex);
  task128->add_dep(task132);
  task132->add_dep(task69);
  residualq->add_task(task132);

  auto tensor133 = vector<shared_ptr<Tensor>>{I80, Gamma34_(), t2};
  auto task133 = make_shared<Task133>(tensor133, pindex);
  task120->add_dep(task133);
  task133->add_dep(task69);
  residualq->add_task(task133);

  vector<IndexRange> I100_index = {closed_, virt_, active_, active_};
  auto I100 = make_shared<Tensor>(I100_index);
  auto tensor134 = vector<shared_ptr<Tensor>>{I80, Gamma35_(), I100};
  auto task134 = make_shared<Task134>(tensor134, pindex);
  task120->add_dep(task134);
  task134->add_dep(task69);
  residualq->add_task(task134);

  auto tensor135 = vector<shared_ptr<Tensor>>{I100, t2};
  auto task135 = make_shared<Task135>(tensor135, pindex, this->e0_);
  task134->add_dep(task135);
  task135->add_dep(task69);
  residualq->add_task(task135);

  auto tensor136 = vector<shared_ptr<Tensor>>{I100, t2, f1_};
  auto task136 = make_shared<Task136>(tensor136, pindex);
  task134->add_dep(task136);
  task136->add_dep(task69);
  residualq->add_task(task136);

  auto tensor137 = vector<shared_ptr<Tensor>>{I100, t2, f1_};
  auto task137 = make_shared<Task137>(tensor137, pindex);
  task134->add_dep(task137);
  task137->add_dep(task69);
  residualq->add_task(task137);

  auto tensor138 = vector<shared_ptr<Tensor>>{I100, t2, f1_};
  auto task138 = make_shared<Task138>(tensor138, pindex);
  task134->add_dep(task138);
  task138->add_dep(task69);
  residualq->add_task(task138);

  vector<IndexRange> I106_index = {virt_, active_, active_, active_};
  auto I106 = make_shared<Tensor>(I106_index);
  auto tensor139 = vector<shared_ptr<Tensor>>{I80, f1_, I106};
  auto task139 = make_shared<Task139>(tensor139, pindex);
  task120->add_dep(task139);
  task139->add_dep(task69);
  residualq->add_task(task139);

  auto tensor140 = vector<shared_ptr<Tensor>>{I106, Gamma37_(), t2};
  auto task140 = make_shared<Task140>(tensor140, pindex);
  task139->add_dep(task140);
  task140->add_dep(task69);
  residualq->add_task(task140);

  vector<IndexRange> I109_index = {closed_, virt_};
  auto I109 = make_shared<Tensor>(I109_index);
  auto tensor141 = vector<shared_ptr<Tensor>>{I80, Gamma38_(), I109};
  auto task141 = make_shared<Task141>(tensor141, pindex);
  task120->add_dep(task141);
  task141->add_dep(task69);
  residualq->add_task(task141);

  auto tensor142 = vector<shared_ptr<Tensor>>{I109, t2, f1_};
  auto task142 = make_shared<Task142>(tensor142, pindex);
  task141->add_dep(task142);
  task142->add_dep(task69);
  residualq->add_task(task142);

  auto tensor143 = vector<shared_ptr<Tensor>>{I109, t2, f1_};
  auto task143 = make_shared<Task143>(tensor143, pindex);
  task141->add_dep(task143);
  task143->add_dep(task69);
  residualq->add_task(task143);
}

#endif
