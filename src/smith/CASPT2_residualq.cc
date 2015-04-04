//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_residualqq.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_residualq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto residualq = make_shared<Queue>();
  auto tensor69 = vector<shared_ptr<Tensor>>{r};
  auto task69 = make_shared<Task69>(tensor69, reset);
  residualq->add_task(task69);

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

  vector<IndexRange> I265_index = {closed_, active_, closed_, active_};
  auto I265 = make_shared<Tensor>(I265_index);
  auto tensor72 = vector<shared_ptr<Tensor>>{I0, Gamma92_(), I265};
  auto task72 = make_shared<Task72>(tensor72, pindex);
  task70->add_dep(task72);
  task72->add_dep(task69);
  residualq->add_task(task72);

  auto tensor73 = vector<shared_ptr<Tensor>>{I265, t2, v2_};
  auto task73 = make_shared<Task73>(tensor73, pindex, this->e0_);
  task72->add_dep(task73);
  task73->add_dep(task69);
  residualq->add_task(task73);

  vector<IndexRange> I2_index = {closed_, closed_, active_, active_};
  auto I2 = make_shared<Tensor>(I2_index);
  auto tensor74 = vector<shared_ptr<Tensor>>{r, I2};
  auto task74 = make_shared<Task74>(tensor74, pindex);
  task74->add_dep(task69);
  residualq->add_task(task74);

  vector<IndexRange> I3_index = {closed_, closed_, active_, active_};
  auto I3 = make_shared<Tensor>(I3_index);
  auto tensor75 = vector<shared_ptr<Tensor>>{I2, Gamma92_(), I3};
  auto task75 = make_shared<Task75>(tensor75, pindex);
  task74->add_dep(task75);
  task75->add_dep(task69);
  residualq->add_task(task75);

  auto tensor76 = vector<shared_ptr<Tensor>>{I3, t2, f1_};
  auto task76 = make_shared<Task76>(tensor76, pindex);
  task75->add_dep(task76);
  task76->add_dep(task69);
  residualq->add_task(task76);

  vector<IndexRange> I6_index = {closed_, active_, active_, active_};
  auto I6 = make_shared<Tensor>(I6_index);
  auto tensor77 = vector<shared_ptr<Tensor>>{I2, f1_, I6};
  auto task77 = make_shared<Task77>(tensor77, pindex);
  task74->add_dep(task77);
  task77->add_dep(task69);
  residualq->add_task(task77);

  auto tensor78 = vector<shared_ptr<Tensor>>{I6, Gamma2_(), t2};
  auto task78 = make_shared<Task78>(tensor78, pindex);
  task77->add_dep(task78);
  task78->add_dep(task69);
  residualq->add_task(task78);

  vector<IndexRange> I9_index = {active_, closed_, closed_, active_};
  auto I9 = make_shared<Tensor>(I9_index);
  auto tensor79 = vector<shared_ptr<Tensor>>{I2, Gamma3_(), I9};
  auto task79 = make_shared<Task79>(tensor79, pindex);
  task74->add_dep(task79);
  task79->add_dep(task69);
  residualq->add_task(task79);

  auto tensor80 = vector<shared_ptr<Tensor>>{I9, t2, f1_};
  auto task80 = make_shared<Task80>(tensor80, pindex);
  task79->add_dep(task80);
  task80->add_dep(task69);
  residualq->add_task(task80);

  vector<IndexRange> I11_index = {closed_, active_, active_, active_};
  auto I11 = make_shared<Tensor>(I11_index);
  auto tensor81 = vector<shared_ptr<Tensor>>{r, I11};
  auto task81 = make_shared<Task81>(tensor81, pindex);
  task81->add_dep(task69);
  residualq->add_task(task81);

  vector<IndexRange> I12_index = {active_, closed_, active_, active_};
  auto I12 = make_shared<Tensor>(I12_index);
  auto tensor82 = vector<shared_ptr<Tensor>>{I11, Gamma4_(), I12};
  auto task82 = make_shared<Task82>(tensor82, pindex);
  task81->add_dep(task82);
  task82->add_dep(task69);
  residualq->add_task(task82);

  auto tensor83 = vector<shared_ptr<Tensor>>{I12, t2, f1_};
  auto task83 = make_shared<Task83>(tensor83, pindex);
  task82->add_dep(task83);
  task83->add_dep(task69);
  residualq->add_task(task83);

  auto tensor84 = vector<shared_ptr<Tensor>>{I11, Gamma5_(), t2};
  auto task84 = make_shared<Task84>(tensor84, pindex);
  task81->add_dep(task84);
  task84->add_dep(task69);
  residualq->add_task(task84);

  vector<IndexRange> I17_index = {closed_, active_, active_, active_};
  auto I17 = make_shared<Tensor>(I17_index);
  auto tensor85 = vector<shared_ptr<Tensor>>{I11, Gamma6_(), I17};
  auto task85 = make_shared<Task85>(tensor85, pindex);
  task81->add_dep(task85);
  task85->add_dep(task69);
  residualq->add_task(task85);

  auto tensor86 = vector<shared_ptr<Tensor>>{I17, t2, v2_};
  auto task86 = make_shared<Task86>(tensor86, pindex, this->e0_);
  task85->add_dep(task86);
  task86->add_dep(task69);
  residualq->add_task(task86);

  auto tensor87 = vector<shared_ptr<Tensor>>{I17, t2, f1_};
  auto task87 = make_shared<Task87>(tensor87, pindex);
  task85->add_dep(task87);
  task87->add_dep(task69);
  residualq->add_task(task87);

  auto tensor88 = vector<shared_ptr<Tensor>>{I17, t2, f1_};
  auto task88 = make_shared<Task88>(tensor88, pindex);
  task85->add_dep(task88);
  task88->add_dep(task69);
  residualq->add_task(task88);

  vector<IndexRange> I20_index = {closed_, active_};
  auto I20 = make_shared<Tensor>(I20_index);
  auto tensor89 = vector<shared_ptr<Tensor>>{I11, Gamma7_(), I20};
  auto task89 = make_shared<Task89>(tensor89, pindex);
  task81->add_dep(task89);
  task89->add_dep(task69);
  residualq->add_task(task89);

  auto tensor90 = vector<shared_ptr<Tensor>>{I20, h1_};
  auto task90 = make_shared<Task90>(tensor90, pindex);
  task89->add_dep(task90);
  task90->add_dep(task69);
  residualq->add_task(task90);

  auto tensor91 = vector<shared_ptr<Tensor>>{I20, t2, f1_};
  auto task91 = make_shared<Task91>(tensor91, pindex);
  task89->add_dep(task91);
  task91->add_dep(task69);
  residualq->add_task(task91);

  auto tensor92 = vector<shared_ptr<Tensor>>{I20, t2, f1_};
  auto task92 = make_shared<Task92>(tensor92, pindex);
  task89->add_dep(task92);
  task92->add_dep(task69);
  residualq->add_task(task92);

  vector<IndexRange> I26_index = {active_, active_, closed_, active_};
  auto I26 = make_shared<Tensor>(I26_index);
  auto tensor93 = vector<shared_ptr<Tensor>>{I11, Gamma9_(), I26};
  auto task93 = make_shared<Task93>(tensor93, pindex);
  task81->add_dep(task93);
  task93->add_dep(task69);
  residualq->add_task(task93);

  auto tensor94 = vector<shared_ptr<Tensor>>{I26, t2, f1_};
  auto task94 = make_shared<Task94>(tensor94, pindex);
  task93->add_dep(task94);
  task94->add_dep(task69);
  residualq->add_task(task94);

  auto tensor95 = vector<shared_ptr<Tensor>>{I11, Gamma105_(), v2_};
  auto task95 = make_shared<Task95>(tensor95, pindex);
  task81->add_dep(task95);
  task95->add_dep(task69);
  residualq->add_task(task95);

  vector<IndexRange> I31_index = {closed_, closed_, active_, virt_};
  auto I31 = make_shared<Tensor>(I31_index);
  auto tensor96 = vector<shared_ptr<Tensor>>{r, I31};
  auto task96 = make_shared<Task96>(tensor96, pindex);
  task96->add_dep(task69);
  residualq->add_task(task96);

  vector<IndexRange> I32_index = {closed_, closed_, active_, active_};
  auto I32 = make_shared<Tensor>(I32_index);
  auto tensor97 = vector<shared_ptr<Tensor>>{I31, f1_, I32};
  auto task97 = make_shared<Task97>(tensor97, pindex);
  task96->add_dep(task97);
  task97->add_dep(task69);
  residualq->add_task(task97);

  auto tensor98 = vector<shared_ptr<Tensor>>{I32, Gamma3_(), t2};
  auto task98 = make_shared<Task98>(tensor98, pindex);
  task97->add_dep(task98);
  task98->add_dep(task69);
  residualq->add_task(task98);

  vector<IndexRange> I35_index = {closed_, active_};
  auto I35 = make_shared<Tensor>(I35_index);
  auto tensor99 = vector<shared_ptr<Tensor>>{I31, f1_, I35};
  auto task99 = make_shared<Task99>(tensor99, pindex);
  task96->add_dep(task99);
  task99->add_dep(task69);
  residualq->add_task(task99);

  auto tensor100 = vector<shared_ptr<Tensor>>{I35, Gamma12_(), t2};
  auto task100 = make_shared<Task100>(tensor100, pindex);
  task99->add_dep(task100);
  task100->add_dep(task69);
  residualq->add_task(task100);

  vector<IndexRange> I38_index = {closed_, active_};
  auto I38 = make_shared<Tensor>(I38_index);
  auto tensor101 = vector<shared_ptr<Tensor>>{I31, f1_, I38};
  auto task101 = make_shared<Task101>(tensor101, pindex);
  task96->add_dep(task101);
  task101->add_dep(task69);
  residualq->add_task(task101);

  auto tensor102 = vector<shared_ptr<Tensor>>{I38, Gamma12_(), t2};
  auto task102 = make_shared<Task102>(tensor102, pindex);
  task101->add_dep(task102);
  task102->add_dep(task69);
  residualq->add_task(task102);

  vector<IndexRange> I41_index = {closed_, virt_, closed_, active_};
  auto I41 = make_shared<Tensor>(I41_index);
  auto tensor103 = vector<shared_ptr<Tensor>>{I31, Gamma14_(), I41};
  auto task103 = make_shared<Task103>(tensor103, pindex);
  task96->add_dep(task103);
  task103->add_dep(task69);
  residualq->add_task(task103);

  auto tensor104 = vector<shared_ptr<Tensor>>{I41, t2};
  auto task104 = make_shared<Task104>(tensor104, pindex);
  task103->add_dep(task104);
  task104->add_dep(task69);
  residualq->add_task(task104);

  vector<IndexRange> I45_index = {closed_, virt_, closed_, active_};
  auto I45 = make_shared<Tensor>(I45_index);
  auto tensor105 = vector<shared_ptr<Tensor>>{I31, Gamma16_(), I45};
  auto task105 = make_shared<Task105>(tensor105, pindex);
  task96->add_dep(task105);
  task105->add_dep(task69);
  residualq->add_task(task105);

  auto tensor106 = vector<shared_ptr<Tensor>>{I45, t2, v2_};
  auto task106 = make_shared<Task106>(tensor106, pindex, this->e0_);
  task105->add_dep(task106);
  task106->add_dep(task69);
  residualq->add_task(task106);

  auto tensor107 = vector<shared_ptr<Tensor>>{I45, t2, f1_};
  auto task107 = make_shared<Task107>(tensor107, pindex);
  task105->add_dep(task107);
  task107->add_dep(task69);
  residualq->add_task(task107);

  auto tensor108 = vector<shared_ptr<Tensor>>{I45, t2, f1_};
  auto task108 = make_shared<Task108>(tensor108, pindex);
  task105->add_dep(task108);
  task108->add_dep(task69);
  residualq->add_task(task108);

  auto tensor109 = vector<shared_ptr<Tensor>>{I45, t2, f1_};
  auto task109 = make_shared<Task109>(tensor109, pindex);
  task105->add_dep(task109);
  task109->add_dep(task69);
  residualq->add_task(task109);

  auto tensor110 = vector<shared_ptr<Tensor>>{I45, t2, f1_};
  auto task110 = make_shared<Task110>(tensor110, pindex);
  task105->add_dep(task110);
  task110->add_dep(task69);
  residualq->add_task(task110);

  auto tensor111 = vector<shared_ptr<Tensor>>{I45, t2, f1_};
  auto task111 = make_shared<Task111>(tensor111, pindex);
  task105->add_dep(task111);
  task111->add_dep(task69);
  residualq->add_task(task111);

  auto tensor112 = vector<shared_ptr<Tensor>>{I45, t2, f1_};
  auto task112 = make_shared<Task112>(tensor112, pindex);
  task105->add_dep(task112);
  task112->add_dep(task69);
  residualq->add_task(task112);

  vector<IndexRange> I63_index = {virt_, closed_, active_, active_};
  auto I63 = make_shared<Tensor>(I63_index);
  auto tensor113 = vector<shared_ptr<Tensor>>{I31, f1_, I63};
  auto task113 = make_shared<Task113>(tensor113, pindex);
  task96->add_dep(task113);
  task113->add_dep(task69);
  residualq->add_task(task113);

  auto tensor114 = vector<shared_ptr<Tensor>>{I63, Gamma22_(), t2};
  auto task114 = make_shared<Task114>(tensor114, pindex);
  task113->add_dep(task114);
  task114->add_dep(task69);
  residualq->add_task(task114);

  auto tensor115 = vector<shared_ptr<Tensor>>{I63, Gamma12_(), t2};
  auto task115 = make_shared<Task115>(tensor115, pindex);
  task113->add_dep(task115);
  task115->add_dep(task69);
  residualq->add_task(task115);

  vector<IndexRange> I66_index = {virt_, closed_, active_, active_};
  auto I66 = make_shared<Tensor>(I66_index);
  auto tensor116 = vector<shared_ptr<Tensor>>{I31, f1_, I66};
  auto task116 = make_shared<Task116>(tensor116, pindex);
  task96->add_dep(task116);
  task116->add_dep(task69);
  residualq->add_task(task116);

  vector<IndexRange> I67_index = {active_, virt_, closed_, active_};
  auto I67 = make_shared<Tensor>(I67_index);
  auto tensor117 = vector<shared_ptr<Tensor>>{I66, Gamma12_(), I67};
  auto task117 = make_shared<Task117>(tensor117, pindex);
  task116->add_dep(task117);
  task117->add_dep(task69);
  residualq->add_task(task117);

  auto tensor118 = vector<shared_ptr<Tensor>>{I67, t2};
  auto task118 = make_shared<Task118>(tensor118, pindex);
  task117->add_dep(task118);
  task118->add_dep(task69);
  residualq->add_task(task118);

  vector<IndexRange> I75_index = {virt_, active_};
  auto I75 = make_shared<Tensor>(I75_index);
  auto tensor119 = vector<shared_ptr<Tensor>>{I31, t2, I75};
  auto task119 = make_shared<Task119>(tensor119, pindex);
  task96->add_dep(task119);
  task119->add_dep(task69);
  residualq->add_task(task119);

  auto tensor120 = vector<shared_ptr<Tensor>>{I75, Gamma16_(), f1_};
  auto task120 = make_shared<Task120>(tensor120, pindex);
  task119->add_dep(task120);
  task120->add_dep(task69);
  residualq->add_task(task120);

  vector<IndexRange> I78_index = {virt_, active_};
  auto I78 = make_shared<Tensor>(I78_index);
  auto tensor121 = vector<shared_ptr<Tensor>>{I31, t2, I78};
  auto task121 = make_shared<Task121>(tensor121, pindex);
  task96->add_dep(task121);
  task121->add_dep(task69);
  residualq->add_task(task121);

  auto tensor122 = vector<shared_ptr<Tensor>>{I78, Gamma16_(), f1_};
  auto task122 = make_shared<Task122>(tensor122, pindex);
  task121->add_dep(task122);
  task122->add_dep(task69);
  residualq->add_task(task122);

  vector<IndexRange> I80_index = {closed_, active_, active_, virt_};
  auto I80 = make_shared<Tensor>(I80_index);
  auto tensor123 = vector<shared_ptr<Tensor>>{r, I80};
  auto task123 = make_shared<Task123>(tensor123, pindex);
  task123->add_dep(task69);
  residualq->add_task(task123);

  vector<IndexRange> I81_index = {closed_, active_, active_, active_};
  auto I81 = make_shared<Tensor>(I81_index);
  auto tensor124 = vector<shared_ptr<Tensor>>{I80, f1_, I81};
  auto task124 = make_shared<Task124>(tensor124, pindex);
  task123->add_dep(task124);
  task124->add_dep(task69);
  residualq->add_task(task124);

  auto tensor125 = vector<shared_ptr<Tensor>>{I81, Gamma28_(), t2};
  auto task125 = make_shared<Task125>(tensor125, pindex);
  task124->add_dep(task125);
  task125->add_dep(task69);
  residualq->add_task(task125);

  vector<IndexRange> I84_index = {active_, virt_, closed_, active_};
  auto I84 = make_shared<Tensor>(I84_index);
  auto tensor126 = vector<shared_ptr<Tensor>>{I80, Gamma29_(), I84};
  auto task126 = make_shared<Task126>(tensor126, pindex);
  task123->add_dep(task126);
  task126->add_dep(task69);
  residualq->add_task(task126);

  auto tensor127 = vector<shared_ptr<Tensor>>{I84, v2_};
  auto task127 = make_shared<Task127>(tensor127, pindex);
  task126->add_dep(task127);
  task127->add_dep(task69);
  residualq->add_task(task127);

  auto tensor128 = vector<shared_ptr<Tensor>>{I84, t2, f1_};
  auto task128 = make_shared<Task128>(tensor128, pindex);
  task126->add_dep(task128);
  task128->add_dep(task69);
  residualq->add_task(task128);

  vector<IndexRange> I87_index = {active_, closed_, virt_, active_};
  auto I87 = make_shared<Tensor>(I87_index);
  auto tensor129 = vector<shared_ptr<Tensor>>{I80, Gamma7_(), I87};
  auto task129 = make_shared<Task129>(tensor129, pindex);
  task123->add_dep(task129);
  task129->add_dep(task69);
  residualq->add_task(task129);

  auto tensor130 = vector<shared_ptr<Tensor>>{I87, t2, f1_};
  auto task130 = make_shared<Task130>(tensor130, pindex);
  task129->add_dep(task130);
  task130->add_dep(task69);
  residualq->add_task(task130);

  auto tensor131 = vector<shared_ptr<Tensor>>{I80, Gamma31_(), t2};
  auto task131 = make_shared<Task131>(tensor131, pindex);
  task123->add_dep(task131);
  task131->add_dep(task69);
  residualq->add_task(task131);

  vector<IndexRange> I92_index = {closed_, active_, virt_, active_};
  auto I92 = make_shared<Tensor>(I92_index);
  auto tensor132 = vector<shared_ptr<Tensor>>{I80, Gamma32_(), I92};
  auto task132 = make_shared<Task132>(tensor132, pindex);
  task123->add_dep(task132);
  task132->add_dep(task69);
  residualq->add_task(task132);

  auto tensor133 = vector<shared_ptr<Tensor>>{I92, t2, v2_};
  auto task133 = make_shared<Task133>(tensor133, pindex, this->e0_);
  task132->add_dep(task133);
  task133->add_dep(task69);
  residualq->add_task(task133);

  auto tensor134 = vector<shared_ptr<Tensor>>{I92, t2, f1_};
  auto task134 = make_shared<Task134>(tensor134, pindex);
  task132->add_dep(task134);
  task134->add_dep(task69);
  residualq->add_task(task134);

  auto tensor135 = vector<shared_ptr<Tensor>>{I92, t2, f1_};
  auto task135 = make_shared<Task135>(tensor135, pindex);
  task132->add_dep(task135);
  task135->add_dep(task69);
  residualq->add_task(task135);

  auto tensor136 = vector<shared_ptr<Tensor>>{I92, t2, f1_};
  auto task136 = make_shared<Task136>(tensor136, pindex);
  task132->add_dep(task136);
  task136->add_dep(task69);
  residualq->add_task(task136);

  auto tensor137 = vector<shared_ptr<Tensor>>{I80, Gamma34_(), t2};
  auto task137 = make_shared<Task137>(tensor137, pindex);
  task123->add_dep(task137);
  task137->add_dep(task69);
  residualq->add_task(task137);

  vector<IndexRange> I100_index = {closed_, virt_, active_, active_};
  auto I100 = make_shared<Tensor>(I100_index);
  auto tensor138 = vector<shared_ptr<Tensor>>{I80, Gamma35_(), I100};
  auto task138 = make_shared<Task138>(tensor138, pindex);
  task123->add_dep(task138);
  task138->add_dep(task69);
  residualq->add_task(task138);

  auto tensor139 = vector<shared_ptr<Tensor>>{I100, t2, v2_};
  auto task139 = make_shared<Task139>(tensor139, pindex, this->e0_);
  task138->add_dep(task139);
  task139->add_dep(task69);
  residualq->add_task(task139);

  auto tensor140 = vector<shared_ptr<Tensor>>{I100, t2, f1_};
  auto task140 = make_shared<Task140>(tensor140, pindex);
  task138->add_dep(task140);
  task140->add_dep(task69);
  residualq->add_task(task140);

  auto tensor141 = vector<shared_ptr<Tensor>>{I100, t2, f1_};
  auto task141 = make_shared<Task141>(tensor141, pindex);
  task138->add_dep(task141);
  task141->add_dep(task69);
  residualq->add_task(task141);

  auto tensor142 = vector<shared_ptr<Tensor>>{I100, t2, f1_};
  auto task142 = make_shared<Task142>(tensor142, pindex);
  task138->add_dep(task142);
  task142->add_dep(task69);
  residualq->add_task(task142);

  vector<IndexRange> I106_index = {virt_, active_, active_, active_};
  auto I106 = make_shared<Tensor>(I106_index);
  auto tensor143 = vector<shared_ptr<Tensor>>{I80, f1_, I106};
  auto task143 = make_shared<Task143>(tensor143, pindex);
  task123->add_dep(task143);
  task143->add_dep(task69);
  residualq->add_task(task143);

  auto tensor144 = vector<shared_ptr<Tensor>>{I106, Gamma37_(), t2};
  auto task144 = make_shared<Task144>(tensor144, pindex);
  task143->add_dep(task144);
  task144->add_dep(task69);
  residualq->add_task(task144);

  vector<IndexRange> I109_index = {closed_, virt_};
  auto I109 = make_shared<Tensor>(I109_index);
  auto tensor145 = vector<shared_ptr<Tensor>>{I80, Gamma38_(), I109};
  auto task145 = make_shared<Task145>(tensor145, pindex);
  task123->add_dep(task145);
  task145->add_dep(task69);
  residualq->add_task(task145);

  auto tensor146 = vector<shared_ptr<Tensor>>{I109, h1_};
  auto task146 = make_shared<Task146>(tensor146, pindex);
  task145->add_dep(task146);
  task146->add_dep(task69);
  residualq->add_task(task146);

  auto tensor147 = vector<shared_ptr<Tensor>>{I109, t2, f1_};
  auto task147 = make_shared<Task147>(tensor147, pindex);
  task145->add_dep(task147);
  task147->add_dep(task69);
  residualq->add_task(task147);

  auto tensor148 = vector<shared_ptr<Tensor>>{I109, t2, f1_};
  auto task148 = make_shared<Task148>(tensor148, pindex);
  task145->add_dep(task148);
  task148->add_dep(task69);
  residualq->add_task(task148);

  vector<IndexRange> I120_index = {closed_, active_, active_, virt_};
  auto I120 = make_shared<Tensor>(I120_index);
  auto tensor149 = vector<shared_ptr<Tensor>>{r, I120};
  auto task149 = make_shared<Task149>(tensor149, pindex);
  task149->add_dep(task69);
  residualq->add_task(task149);

  vector<IndexRange> I121_index = {closed_, active_, active_, active_};
  auto I121 = make_shared<Tensor>(I121_index);
  auto tensor150 = vector<shared_ptr<Tensor>>{I120, f1_, I121};
  auto task150 = make_shared<Task150>(tensor150, pindex);
  task149->add_dep(task150);
  task150->add_dep(task69);
  residualq->add_task(task150);

  auto tensor151 = vector<shared_ptr<Tensor>>{I121, Gamma6_(), t2};
  auto task151 = make_shared<Task151>(tensor151, pindex);
  task150->add_dep(task151);
  task151->add_dep(task69);
  residualq->add_task(task151);

  vector<IndexRange> I124_index = {active_, virt_, closed_, active_};
  auto I124 = make_shared<Tensor>(I124_index);
  auto tensor152 = vector<shared_ptr<Tensor>>{I120, Gamma7_(), I124};
  auto task152 = make_shared<Task152>(tensor152, pindex);
  task149->add_dep(task152);
  task152->add_dep(task69);
  residualq->add_task(task152);

  auto tensor153 = vector<shared_ptr<Tensor>>{I124, v2_};
  auto task153 = make_shared<Task153>(tensor153, pindex);
  task152->add_dep(task153);
  task153->add_dep(task69);
  residualq->add_task(task153);

  auto tensor154 = vector<shared_ptr<Tensor>>{I124, t2, f1_};
  auto task154 = make_shared<Task154>(tensor154, pindex);
  task152->add_dep(task154);
  task154->add_dep(task69);
  residualq->add_task(task154);

  auto tensor155 = vector<shared_ptr<Tensor>>{I124, t2, f1_};
  auto task155 = make_shared<Task155>(tensor155, pindex);
  task152->add_dep(task155);
  task155->add_dep(task69);
  residualq->add_task(task155);

  vector<IndexRange> I130_index = {active_, virt_, closed_, active_};
  auto I130 = make_shared<Tensor>(I130_index);
  auto tensor156 = vector<shared_ptr<Tensor>>{I120, Gamma34_(), I130};
  auto task156 = make_shared<Task156>(tensor156, pindex);
  task149->add_dep(task156);
  task156->add_dep(task69);
  residualq->add_task(task156);

  auto tensor157 = vector<shared_ptr<Tensor>>{I130, t2};
  auto task157 = make_shared<Task157>(tensor157, pindex);
  task156->add_dep(task157);
  task157->add_dep(task69);
  residualq->add_task(task157);

  vector<IndexRange> I132_index = {closed_, active_, virt_, active_};
  auto I132 = make_shared<Tensor>(I132_index);
  auto tensor158 = vector<shared_ptr<Tensor>>{I120, Gamma35_(), I132};
  auto task158 = make_shared<Task158>(tensor158, pindex);
  task149->add_dep(task158);
  task158->add_dep(task69);
  residualq->add_task(task158);

  auto tensor159 = vector<shared_ptr<Tensor>>{I132, t2, v2_};
  auto task159 = make_shared<Task159>(tensor159, pindex, this->e0_);
  task158->add_dep(task159);
  task159->add_dep(task69);
  residualq->add_task(task159);

  auto tensor160 = vector<shared_ptr<Tensor>>{I132, t2, f1_};
  auto task160 = make_shared<Task160>(tensor160, pindex);
  task158->add_dep(task160);
  task160->add_dep(task69);
  residualq->add_task(task160);

  auto tensor161 = vector<shared_ptr<Tensor>>{I132, t2, f1_};
  auto task161 = make_shared<Task161>(tensor161, pindex);
  task158->add_dep(task161);
  task161->add_dep(task69);
  residualq->add_task(task161);

  auto tensor162 = vector<shared_ptr<Tensor>>{I132, t2, f1_};
  auto task162 = make_shared<Task162>(tensor162, pindex);
  task158->add_dep(task162);
  task162->add_dep(task69);
  residualq->add_task(task162);

  auto tensor163 = vector<shared_ptr<Tensor>>{I132, t2, f1_};
  auto task163 = make_shared<Task163>(tensor163, pindex);
  task158->add_dep(task163);
  task163->add_dep(task69);
  residualq->add_task(task163);

  auto tensor164 = vector<shared_ptr<Tensor>>{I132, t2, f1_};
  auto task164 = make_shared<Task164>(tensor164, pindex);
  task158->add_dep(task164);
  task164->add_dep(task69);
  residualq->add_task(task164);

  auto tensor165 = vector<shared_ptr<Tensor>>{I132, t2, f1_};
  auto task165 = make_shared<Task165>(tensor165, pindex);
  task158->add_dep(task165);
  task165->add_dep(task69);
  residualq->add_task(task165);

  vector<IndexRange> I146_index = {virt_, active_, active_, active_};
  auto I146 = make_shared<Tensor>(I146_index);
  auto tensor166 = vector<shared_ptr<Tensor>>{I120, f1_, I146};
  auto task166 = make_shared<Task166>(tensor166, pindex);
  task149->add_dep(task166);
  task166->add_dep(task69);
  residualq->add_task(task166);

  auto tensor167 = vector<shared_ptr<Tensor>>{I146, Gamma51_(), t2};
  auto task167 = make_shared<Task167>(tensor167, pindex);
  task166->add_dep(task167);
  task167->add_dep(task69);
  residualq->add_task(task167);

  vector<IndexRange> I149_index = {closed_, virt_};
  auto I149 = make_shared<Tensor>(I149_index);
  auto tensor168 = vector<shared_ptr<Tensor>>{I120, Gamma38_(), I149};
  auto task168 = make_shared<Task168>(tensor168, pindex);
  task149->add_dep(task168);
  task168->add_dep(task69);
  residualq->add_task(task168);

  auto tensor169 = vector<shared_ptr<Tensor>>{I149, h1_};
  auto task169 = make_shared<Task169>(tensor169, pindex);
  task168->add_dep(task169);
  task169->add_dep(task69);
  residualq->add_task(task169);

  auto tensor170 = vector<shared_ptr<Tensor>>{I149, t2, f1_};
  auto task170 = make_shared<Task170>(tensor170, pindex);
  task168->add_dep(task170);
  task170->add_dep(task69);
  residualq->add_task(task170);

  auto tensor171 = vector<shared_ptr<Tensor>>{I149, t2, f1_};
  auto task171 = make_shared<Task171>(tensor171, pindex);
  task168->add_dep(task171);
  task171->add_dep(task69);
  residualq->add_task(task171);

  vector<IndexRange> I160_index = {virt_, active_, active_, active_};
  auto I160 = make_shared<Tensor>(I160_index);
  auto tensor172 = vector<shared_ptr<Tensor>>{r, I160};
  auto task172 = make_shared<Task172>(tensor172, pindex);
  task172->add_dep(task69);
  residualq->add_task(task172);

  vector<IndexRange> I161_index = {active_, active_, virt_, active_};
  auto I161 = make_shared<Tensor>(I161_index);
  auto tensor173 = vector<shared_ptr<Tensor>>{I160, Gamma56_(), I161};
  auto task173 = make_shared<Task173>(tensor173, pindex);
  task172->add_dep(task173);
  task173->add_dep(task69);
  residualq->add_task(task173);

  auto tensor174 = vector<shared_ptr<Tensor>>{I161, t2, f1_};
  auto task174 = make_shared<Task174>(tensor174, pindex);
  task173->add_dep(task174);
  task174->add_dep(task69);
  residualq->add_task(task174);

  vector<IndexRange> I164_index = {active_, virt_, active_, active_};
  auto I164 = make_shared<Tensor>(I164_index);
  auto tensor175 = vector<shared_ptr<Tensor>>{I160, Gamma57_(), I164};
  auto task175 = make_shared<Task175>(tensor175, pindex);
  task172->add_dep(task175);
  task175->add_dep(task69);
  residualq->add_task(task175);

  auto tensor176 = vector<shared_ptr<Tensor>>{I164, v2_};
  auto task176 = make_shared<Task176>(tensor176, pindex);
  task175->add_dep(task176);
  task176->add_dep(task69);
  residualq->add_task(task176);

  auto tensor177 = vector<shared_ptr<Tensor>>{I164, t2, f1_};
  auto task177 = make_shared<Task177>(tensor177, pindex);
  task175->add_dep(task177);
  task177->add_dep(task69);
  residualq->add_task(task177);

  auto tensor178 = vector<shared_ptr<Tensor>>{I160, Gamma58_(), t2};
  auto task178 = make_shared<Task178>(tensor178, pindex);
  task172->add_dep(task178);
  task178->add_dep(task69);
  residualq->add_task(task178);

  vector<IndexRange> I169_index = {virt_, active_, active_, active_};
  auto I169 = make_shared<Tensor>(I169_index);
  auto tensor179 = vector<shared_ptr<Tensor>>{I160, Gamma59_(), I169};
  auto task179 = make_shared<Task179>(tensor179, pindex);
  task172->add_dep(task179);
  task179->add_dep(task69);
  residualq->add_task(task179);

  auto tensor180 = vector<shared_ptr<Tensor>>{I169, t2, v2_};
  auto task180 = make_shared<Task180>(tensor180, pindex, this->e0_);
  task179->add_dep(task180);
  task180->add_dep(task69);
  residualq->add_task(task180);

  auto tensor181 = vector<shared_ptr<Tensor>>{I169, t2, f1_};
  auto task181 = make_shared<Task181>(tensor181, pindex);
  task179->add_dep(task181);
  task181->add_dep(task69);
  residualq->add_task(task181);

  auto tensor182 = vector<shared_ptr<Tensor>>{I169, t2, f1_};
  auto task182 = make_shared<Task182>(tensor182, pindex);
  task179->add_dep(task182);
  task182->add_dep(task69);
  residualq->add_task(task182);

  vector<IndexRange> I172_index = {active_, virt_};
  auto I172 = make_shared<Tensor>(I172_index);
  auto tensor183 = vector<shared_ptr<Tensor>>{I160, Gamma60_(), I172};
  auto task183 = make_shared<Task183>(tensor183, pindex);
  task172->add_dep(task183);
  task183->add_dep(task69);
  residualq->add_task(task183);

  auto tensor184 = vector<shared_ptr<Tensor>>{I172, h1_};
  auto task184 = make_shared<Task184>(tensor184, pindex);
  task183->add_dep(task184);
  task184->add_dep(task69);
  residualq->add_task(task184);

  auto tensor185 = vector<shared_ptr<Tensor>>{I172, t2, f1_};
  auto task185 = make_shared<Task185>(tensor185, pindex);
  task183->add_dep(task185);
  task185->add_dep(task69);
  residualq->add_task(task185);

  auto tensor186 = vector<shared_ptr<Tensor>>{I172, t2, f1_};
  auto task186 = make_shared<Task186>(tensor186, pindex);
  task183->add_dep(task186);
  task186->add_dep(task69);
  residualq->add_task(task186);

  vector<IndexRange> I180_index = {virt_, closed_, virt_, closed_};
  auto I180 = make_shared<Tensor>(I180_index);
  auto tensor187 = vector<shared_ptr<Tensor>>{r, I180};
  auto task187 = make_shared<Task187>(tensor187, pindex);
  task187->add_dep(task69);
  residualq->add_task(task187);

  vector<IndexRange> I181_index = {virt_, active_};
  auto I181 = make_shared<Tensor>(I181_index);
  auto tensor188 = vector<shared_ptr<Tensor>>{I180, t2, I181};
  auto task188 = make_shared<Task188>(tensor188, pindex);
  task187->add_dep(task188);
  task188->add_dep(task69);
  residualq->add_task(task188);

  auto tensor189 = vector<shared_ptr<Tensor>>{I181, Gamma16_(), f1_};
  auto task189 = make_shared<Task189>(tensor189, pindex);
  task188->add_dep(task189);
  task189->add_dep(task69);
  residualq->add_task(task189);

  vector<IndexRange> I184_index = {virt_, active_};
  auto I184 = make_shared<Tensor>(I184_index);
  auto tensor190 = vector<shared_ptr<Tensor>>{I180, t2, I184};
  auto task190 = make_shared<Task190>(tensor190, pindex);
  task187->add_dep(task190);
  task190->add_dep(task69);
  residualq->add_task(task190);

  auto tensor191 = vector<shared_ptr<Tensor>>{I184, Gamma16_(), f1_};
  auto task191 = make_shared<Task191>(tensor191, pindex);
  task190->add_dep(task191);
  task191->add_dep(task69);
  residualq->add_task(task191);

  vector<IndexRange> I187_index = {virt_, closed_};
  auto I187 = make_shared<Tensor>(I187_index);
  auto tensor192 = vector<shared_ptr<Tensor>>{I180, f1_, I187};
  auto task192 = make_shared<Task192>(tensor192, pindex);
  task187->add_dep(task192);
  task192->add_dep(task69);
  residualq->add_task(task192);

  vector<IndexRange> I188_index = {active_, virt_, closed_, active_};
  auto I188 = make_shared<Tensor>(I188_index);
  auto tensor193 = vector<shared_ptr<Tensor>>{I187, Gamma38_(), I188};
  auto task193 = make_shared<Task193>(tensor193, pindex);
  task192->add_dep(task193);
  task193->add_dep(task69);
  residualq->add_task(task193);

  auto tensor194 = vector<shared_ptr<Tensor>>{I188, t2};
  auto task194 = make_shared<Task194>(tensor194, pindex);
  task193->add_dep(task194);
  task194->add_dep(task69);
  residualq->add_task(task194);

  vector<IndexRange> I190_index = {virt_, closed_};
  auto I190 = make_shared<Tensor>(I190_index);
  auto tensor195 = vector<shared_ptr<Tensor>>{I180, f1_, I190};
  auto task195 = make_shared<Task195>(tensor195, pindex);
  task187->add_dep(task195);
  task195->add_dep(task69);
  residualq->add_task(task195);

  vector<IndexRange> I191_index = {active_, virt_, closed_, active_};
  auto I191 = make_shared<Tensor>(I191_index);
  auto tensor196 = vector<shared_ptr<Tensor>>{I190, Gamma38_(), I191};
  auto task196 = make_shared<Task196>(tensor196, pindex);
  task195->add_dep(task196);
  task196->add_dep(task69);
  residualq->add_task(task196);

  auto tensor197 = vector<shared_ptr<Tensor>>{I191, t2};
  auto task197 = make_shared<Task197>(tensor197, pindex);
  task196->add_dep(task197);
  task197->add_dep(task69);
  residualq->add_task(task197);

  vector<IndexRange> I199_index = {closed_, active_};
  auto I199 = make_shared<Tensor>(I199_index);
  auto tensor198 = vector<shared_ptr<Tensor>>{I180, t2, I199};
  auto task198 = make_shared<Task198>(tensor198, pindex);
  task187->add_dep(task198);
  task198->add_dep(task69);
  residualq->add_task(task198);

  auto tensor199 = vector<shared_ptr<Tensor>>{I199, Gamma38_(), f1_};
  auto task199 = make_shared<Task199>(tensor199, pindex);
  task198->add_dep(task199);
  task199->add_dep(task69);
  residualq->add_task(task199);

  vector<IndexRange> I202_index = {closed_, active_};
  auto I202 = make_shared<Tensor>(I202_index);
  auto tensor200 = vector<shared_ptr<Tensor>>{I180, t2, I202};
  auto task200 = make_shared<Task200>(tensor200, pindex);
  task187->add_dep(task200);
  task200->add_dep(task69);
  residualq->add_task(task200);

  auto tensor201 = vector<shared_ptr<Tensor>>{I202, Gamma38_(), f1_};
  auto task201 = make_shared<Task201>(tensor201, pindex);
  task200->add_dep(task201);
  task201->add_dep(task69);
  residualq->add_task(task201);

  vector<IndexRange> I204_index = {virt_, closed_, active_, virt_};
  auto I204 = make_shared<Tensor>(I204_index);
  auto tensor202 = vector<shared_ptr<Tensor>>{r, I204};
  auto task202 = make_shared<Task202>(tensor202, pindex);
  task202->add_dep(task69);
  residualq->add_task(task202);

  vector<IndexRange> I205_index = {virt_, closed_, active_, active_};
  auto I205 = make_shared<Tensor>(I205_index);
  auto tensor203 = vector<shared_ptr<Tensor>>{I204, f1_, I205};
  auto task203 = make_shared<Task203>(tensor203, pindex);
  task202->add_dep(task203);
  task203->add_dep(task69);
  residualq->add_task(task203);

  vector<IndexRange> I206_index = {active_, virt_, closed_, active_};
  auto I206 = make_shared<Tensor>(I206_index);
  auto tensor204 = vector<shared_ptr<Tensor>>{I205, Gamma35_(), I206};
  auto task204 = make_shared<Task204>(tensor204, pindex);
  task203->add_dep(task204);
  task204->add_dep(task69);
  residualq->add_task(task204);

  auto tensor205 = vector<shared_ptr<Tensor>>{I206, t2};
  auto task205 = make_shared<Task205>(tensor205, pindex);
  task204->add_dep(task205);
  task205->add_dep(task69);
  residualq->add_task(task205);

  vector<IndexRange> I208_index = {virt_, closed_, active_, active_};
  auto I208 = make_shared<Tensor>(I208_index);
  auto tensor206 = vector<shared_ptr<Tensor>>{I204, f1_, I208};
  auto task206 = make_shared<Task206>(tensor206, pindex);
  task202->add_dep(task206);
  task206->add_dep(task69);
  residualq->add_task(task206);

  auto tensor207 = vector<shared_ptr<Tensor>>{I208, Gamma32_(), t2};
  auto task207 = make_shared<Task207>(tensor207, pindex);
  task206->add_dep(task207);
  task207->add_dep(task69);
  residualq->add_task(task207);

  auto tensor208 = vector<shared_ptr<Tensor>>{I208, Gamma35_(), t2};
  auto task208 = make_shared<Task208>(tensor208, pindex);
  task206->add_dep(task208);
  task208->add_dep(task69);
  residualq->add_task(task208);

  vector<IndexRange> I217_index = {virt_, active_};
  auto I217 = make_shared<Tensor>(I217_index);
  auto tensor209 = vector<shared_ptr<Tensor>>{I204, f1_, I217};
  auto task209 = make_shared<Task209>(tensor209, pindex);
  task202->add_dep(task209);
  task209->add_dep(task69);
  residualq->add_task(task209);

  auto tensor210 = vector<shared_ptr<Tensor>>{I217, Gamma60_(), t2};
  auto task210 = make_shared<Task210>(tensor210, pindex);
  task209->add_dep(task210);
  task210->add_dep(task69);
  residualq->add_task(task210);

  vector<IndexRange> I220_index = {virt_, active_};
  auto I220 = make_shared<Tensor>(I220_index);
  auto tensor211 = vector<shared_ptr<Tensor>>{I204, f1_, I220};
  auto task211 = make_shared<Task211>(tensor211, pindex);
  task202->add_dep(task211);
  task211->add_dep(task69);
  residualq->add_task(task211);

  auto tensor212 = vector<shared_ptr<Tensor>>{I220, Gamma60_(), t2};
  auto task212 = make_shared<Task212>(tensor212, pindex);
  task211->add_dep(task212);
  task212->add_dep(task69);
  residualq->add_task(task212);

  vector<IndexRange> I223_index = {closed_, active_};
  auto I223 = make_shared<Tensor>(I223_index);
  auto tensor213 = vector<shared_ptr<Tensor>>{I204, t2, I223};
  auto task213 = make_shared<Task213>(tensor213, pindex);
  task202->add_dep(task213);
  task213->add_dep(task69);
  residualq->add_task(task213);

  auto tensor214 = vector<shared_ptr<Tensor>>{I223, Gamma38_(), f1_};
  auto task214 = make_shared<Task214>(tensor214, pindex);
  task213->add_dep(task214);
  task214->add_dep(task69);
  residualq->add_task(task214);

  vector<IndexRange> I226_index = {closed_, active_};
  auto I226 = make_shared<Tensor>(I226_index);
  auto tensor215 = vector<shared_ptr<Tensor>>{I204, t2, I226};
  auto task215 = make_shared<Task215>(tensor215, pindex);
  task202->add_dep(task215);
  task215->add_dep(task69);
  residualq->add_task(task215);

  auto tensor216 = vector<shared_ptr<Tensor>>{I226, Gamma38_(), f1_};
  auto task216 = make_shared<Task216>(tensor216, pindex);
  task215->add_dep(task216);
  task216->add_dep(task69);
  residualq->add_task(task216);

  vector<IndexRange> I229_index = {active_, virt_, closed_, virt_};
  auto I229 = make_shared<Tensor>(I229_index);
  auto tensor217 = vector<shared_ptr<Tensor>>{I204, Gamma79_(), I229};
  auto task217 = make_shared<Task217>(tensor217, pindex);
  task202->add_dep(task217);
  task217->add_dep(task69);
  residualq->add_task(task217);

  auto tensor218 = vector<shared_ptr<Tensor>>{I229, t2};
  auto task218 = make_shared<Task218>(tensor218, pindex);
  task217->add_dep(task218);
  task218->add_dep(task69);
  residualq->add_task(task218);

  vector<IndexRange> I233_index = {closed_, active_, virt_, virt_};
  auto I233 = make_shared<Tensor>(I233_index);
  auto tensor219 = vector<shared_ptr<Tensor>>{I204, Gamma38_(), I233};
  auto task219 = make_shared<Task219>(tensor219, pindex);
  task202->add_dep(task219);
  task219->add_dep(task69);
  residualq->add_task(task219);

  auto tensor220 = vector<shared_ptr<Tensor>>{I233, t2, v2_};
  auto task220 = make_shared<Task220>(tensor220, pindex, this->e0_);
  task219->add_dep(task220);
  task220->add_dep(task69);
  residualq->add_task(task220);

  auto tensor221 = vector<shared_ptr<Tensor>>{I233, t2, f1_};
  auto task221 = make_shared<Task221>(tensor221, pindex);
  task219->add_dep(task221);
  task221->add_dep(task69);
  residualq->add_task(task221);

  auto tensor222 = vector<shared_ptr<Tensor>>{I233, t2, f1_};
  auto task222 = make_shared<Task222>(tensor222, pindex);
  task219->add_dep(task222);
  task222->add_dep(task69);
  residualq->add_task(task222);

  auto tensor223 = vector<shared_ptr<Tensor>>{I233, t2, f1_};
  auto task223 = make_shared<Task223>(tensor223, pindex);
  task219->add_dep(task223);
  task223->add_dep(task69);
  residualq->add_task(task223);

  auto tensor224 = vector<shared_ptr<Tensor>>{I233, t2, f1_};
  auto task224 = make_shared<Task224>(tensor224, pindex);
  task219->add_dep(task224);
  task224->add_dep(task69);
  residualq->add_task(task224);

  auto tensor225 = vector<shared_ptr<Tensor>>{I233, t2, f1_};
  auto task225 = make_shared<Task225>(tensor225, pindex);
  task219->add_dep(task225);
  task225->add_dep(task69);
  residualq->add_task(task225);

  auto tensor226 = vector<shared_ptr<Tensor>>{I233, t2, f1_};
  auto task226 = make_shared<Task226>(tensor226, pindex);
  task219->add_dep(task226);
  task226->add_dep(task69);
  residualq->add_task(task226);

  vector<IndexRange> I251_index = {virt_, virt_, active_, active_};
  auto I251 = make_shared<Tensor>(I251_index);
  auto tensor227 = vector<shared_ptr<Tensor>>{I204, f1_, I251};
  auto task227 = make_shared<Task227>(tensor227, pindex);
  task202->add_dep(task227);
  task227->add_dep(task69);
  residualq->add_task(task227);

  auto tensor228 = vector<shared_ptr<Tensor>>{I251, Gamma60_(), t2};
  auto task228 = make_shared<Task228>(tensor228, pindex);
  task227->add_dep(task228);
  task228->add_dep(task69);
  residualq->add_task(task228);

  vector<IndexRange> I253_index = {virt_, active_, active_, virt_};
  auto I253 = make_shared<Tensor>(I253_index);
  auto tensor229 = vector<shared_ptr<Tensor>>{r, I253};
  auto task229 = make_shared<Task229>(tensor229, pindex);
  task229->add_dep(task69);
  residualq->add_task(task229);

  vector<IndexRange> I254_index = {virt_, active_, active_, active_};
  auto I254 = make_shared<Tensor>(I254_index);
  auto tensor230 = vector<shared_ptr<Tensor>>{I253, f1_, I254};
  auto task230 = make_shared<Task230>(tensor230, pindex);
  task229->add_dep(task230);
  task230->add_dep(task69);
  residualq->add_task(task230);

  auto tensor231 = vector<shared_ptr<Tensor>>{I254, Gamma59_(), t2};
  auto task231 = make_shared<Task231>(tensor231, pindex);
  task230->add_dep(task231);
  task231->add_dep(task69);
  residualq->add_task(task231);

  vector<IndexRange> I257_index = {active_, active_, virt_, virt_};
  auto I257 = make_shared<Tensor>(I257_index);
  auto tensor232 = vector<shared_ptr<Tensor>>{I253, Gamma60_(), I257};
  auto task232 = make_shared<Task232>(tensor232, pindex);
  task229->add_dep(task232);
  task232->add_dep(task69);
  residualq->add_task(task232);

  auto tensor233 = vector<shared_ptr<Tensor>>{I257, t2, f1_};
  auto task233 = make_shared<Task233>(tensor233, pindex);
  task232->add_dep(task233);
  task233->add_dep(task69);
  residualq->add_task(task233);

  auto tensor234 = vector<shared_ptr<Tensor>>{I257, t2, f1_};
  auto task234 = make_shared<Task234>(tensor234, pindex);
  task232->add_dep(task234);
  task234->add_dep(task69);
  residualq->add_task(task234);

  vector<IndexRange> I259_index = {virt_, virt_, active_, active_};
  auto I259 = make_shared<Tensor>(I259_index);
  auto tensor235 = vector<shared_ptr<Tensor>>{r, I259};
  auto task235 = make_shared<Task235>(tensor235, pindex);
  task235->add_dep(task69);
  residualq->add_task(task235);

  auto tensor236 = vector<shared_ptr<Tensor>>{I259, Gamma90_(), t2};
  auto task236 = make_shared<Task236>(tensor236, pindex);
  task235->add_dep(task236);
  task236->add_dep(task69);
  residualq->add_task(task236);

  vector<IndexRange> I287_index = {active_, virt_, active_, virt_};
  auto I287 = make_shared<Tensor>(I287_index);
  auto tensor237 = vector<shared_ptr<Tensor>>{I259, Gamma60_(), I287};
  auto task237 = make_shared<Task237>(tensor237, pindex);
  task235->add_dep(task237);
  task237->add_dep(task69);
  residualq->add_task(task237);

  auto tensor238 = vector<shared_ptr<Tensor>>{I287, t2, v2_};
  auto task238 = make_shared<Task238>(tensor238, pindex, this->e0_);
  task237->add_dep(task238);
  task238->add_dep(task69);
  residualq->add_task(task238);

  shared_ptr<Tensor> I318;
  if (diagonal) {
    vector<IndexRange> I318_index = {closed_, virt_, closed_, virt_};
    I318 = make_shared<Tensor>(I318_index);
  }
  shared_ptr<Task239> task239;
  if (diagonal) {
    auto tensor239 = vector<shared_ptr<Tensor>>{r, I318};
    task239 = make_shared<Task239>(tensor239, pindex);
    task239->add_dep(task69);
    residualq->add_task(task239);
  }

  shared_ptr<Task240> task240;
  if (diagonal) {
    auto tensor240 = vector<shared_ptr<Tensor>>{I318, v2_};
    task240 = make_shared<Task240>(tensor240, pindex);
    task239->add_dep(task240);
    task240->add_dep(task69);
    residualq->add_task(task240);
  }

  return residualq;
}


#endif
