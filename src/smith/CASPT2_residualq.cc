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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_residualq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto residualq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor70 = {r};
  auto task70 = make_shared<Task70>(tensor70);
  residualq->add_task(task70);

  vector<IndexRange> I0_index = {closed_, closed_, active_, active_};
  auto I0 = make_shared<Tensor>(I0_index);
  vector<shared_ptr<Tensor>> tensor71 = {r, I0};
  auto task71 = make_shared<Task71>(tensor71, pindex);
  task71->add_dep(task70);
  residualq->add_task(task71);

  vector<IndexRange> I1_index = {closed_, active_, closed_, active_};
  auto I1 = make_shared<Tensor>(I1_index);
  vector<shared_ptr<Tensor>> tensor72 = {I0, Gamma0_(), I1};
  auto task72 = make_shared<Task72>(tensor72, pindex);
  task71->add_dep(task72);
  task72->add_dep(task70);
  residualq->add_task(task72);

  vector<shared_ptr<Tensor>> tensor73 = {I1, t2};
  auto task73 = make_shared<Task73>(tensor73, pindex);
  task72->add_dep(task73);
  task73->add_dep(task70);
  residualq->add_task(task73);

  vector<IndexRange> I277_index = {closed_, active_, closed_, active_};
  auto I277 = make_shared<Tensor>(I277_index);
  vector<shared_ptr<Tensor>> tensor74 = {I0, Gamma94_(), I277};
  auto task74 = make_shared<Task74>(tensor74, pindex);
  task71->add_dep(task74);
  task74->add_dep(task70);
  residualq->add_task(task74);

  vector<shared_ptr<Tensor>> tensor75 = {I277, t2, v2_};
  auto task75 = make_shared<Task75>(tensor75, pindex, this->e0_);
  task74->add_dep(task75);
  task75->add_dep(task70);
  residualq->add_task(task75);

  vector<IndexRange> I2_index = {closed_, closed_, active_, active_};
  auto I2 = make_shared<Tensor>(I2_index);
  vector<shared_ptr<Tensor>> tensor76 = {r, I2};
  auto task76 = make_shared<Task76>(tensor76, pindex);
  task76->add_dep(task70);
  residualq->add_task(task76);

  vector<IndexRange> I3_index = {closed_, closed_, active_, active_};
  auto I3 = make_shared<Tensor>(I3_index);
  vector<shared_ptr<Tensor>> tensor77 = {I2, Gamma94_(), I3};
  auto task77 = make_shared<Task77>(tensor77, pindex);
  task76->add_dep(task77);
  task77->add_dep(task70);
  residualq->add_task(task77);

  vector<IndexRange> I4_index = {closed_, closed_};
  auto I4 = make_shared<Tensor>(I4_index);
  vector<shared_ptr<Tensor>> tensor78 = {I3, t2, I4};
  auto task78 = make_shared<Task78>(tensor78, pindex);
  task77->add_dep(task78);
  task78->add_dep(task70);
  residualq->add_task(task78);

  vector<shared_ptr<Tensor>> tensor79 = {I4, f1_};
  auto task79 = make_shared<Task79>(tensor79, pindex);
  task78->add_dep(task79);
  task79->add_dep(task70);
  residualq->add_task(task79);

  vector<IndexRange> I6_index = {closed_, active_, active_, active_};
  auto I6 = make_shared<Tensor>(I6_index);
  vector<shared_ptr<Tensor>> tensor80 = {I2, f1_, I6};
  auto task80 = make_shared<Task80>(tensor80, pindex);
  task76->add_dep(task80);
  task80->add_dep(task70);
  residualq->add_task(task80);

  vector<IndexRange> I7_index = {active_, active_, closed_, active_};
  auto I7 = make_shared<Tensor>(I7_index);
  vector<shared_ptr<Tensor>> tensor81 = {I6, Gamma2_(), I7};
  auto task81 = make_shared<Task81>(tensor81, pindex);
  task80->add_dep(task81);
  task81->add_dep(task70);
  residualq->add_task(task81);

  vector<shared_ptr<Tensor>> tensor82 = {I7, t2};
  auto task82 = make_shared<Task82>(tensor82, pindex);
  task81->add_dep(task82);
  task82->add_dep(task70);
  residualq->add_task(task82);

  vector<IndexRange> I9_index = {active_, closed_, closed_, active_};
  auto I9 = make_shared<Tensor>(I9_index);
  vector<shared_ptr<Tensor>> tensor83 = {I2, Gamma3_(), I9};
  auto task83 = make_shared<Task83>(tensor83, pindex);
  task76->add_dep(task83);
  task83->add_dep(task70);
  residualq->add_task(task83);

  vector<IndexRange> I10_index = {virt_, active_};
  auto I10 = make_shared<Tensor>(I10_index);
  vector<shared_ptr<Tensor>> tensor84 = {I9, t2, I10};
  auto task84 = make_shared<Task84>(tensor84, pindex);
  task83->add_dep(task84);
  task84->add_dep(task70);
  residualq->add_task(task84);

  vector<shared_ptr<Tensor>> tensor85 = {I10, f1_};
  auto task85 = make_shared<Task85>(tensor85, pindex);
  task84->add_dep(task85);
  task85->add_dep(task70);
  residualq->add_task(task85);

  vector<IndexRange> I11_index = {closed_, active_, active_, active_};
  auto I11 = make_shared<Tensor>(I11_index);
  vector<shared_ptr<Tensor>> tensor86 = {r, I11};
  auto task86 = make_shared<Task86>(tensor86, pindex);
  task86->add_dep(task70);
  residualq->add_task(task86);

  vector<IndexRange> I12_index = {active_, closed_, active_, active_};
  auto I12 = make_shared<Tensor>(I12_index);
  vector<shared_ptr<Tensor>> tensor87 = {I11, Gamma4_(), I12};
  auto task87 = make_shared<Task87>(tensor87, pindex);
  task86->add_dep(task87);
  task87->add_dep(task70);
  residualq->add_task(task87);

  vector<IndexRange> I13_index = {active_, closed_};
  auto I13 = make_shared<Tensor>(I13_index);
  vector<shared_ptr<Tensor>> tensor88 = {I12, t2, I13};
  auto task88 = make_shared<Task88>(tensor88, pindex);
  task87->add_dep(task88);
  task88->add_dep(task70);
  residualq->add_task(task88);

  vector<shared_ptr<Tensor>> tensor89 = {I13, f1_};
  auto task89 = make_shared<Task89>(tensor89, pindex);
  task88->add_dep(task89);
  task89->add_dep(task70);
  residualq->add_task(task89);

  vector<IndexRange> I15_index = {active_, active_, closed_, active_};
  auto I15 = make_shared<Tensor>(I15_index);
  vector<shared_ptr<Tensor>> tensor90 = {I11, Gamma5_(), I15};
  auto task90 = make_shared<Task90>(tensor90, pindex);
  task86->add_dep(task90);
  task90->add_dep(task70);
  residualq->add_task(task90);

  vector<shared_ptr<Tensor>> tensor91 = {I15, t2};
  auto task91 = make_shared<Task91>(tensor91, pindex);
  task90->add_dep(task91);
  task91->add_dep(task70);
  residualq->add_task(task91);

  vector<IndexRange> I17_index = {closed_, active_, active_, active_};
  auto I17 = make_shared<Tensor>(I17_index);
  vector<shared_ptr<Tensor>> tensor92 = {I11, Gamma6_(), I17};
  auto task92 = make_shared<Task92>(tensor92, pindex);
  task86->add_dep(task92);
  task92->add_dep(task70);
  residualq->add_task(task92);

  vector<shared_ptr<Tensor>> tensor93 = {I17, t2, v2_};
  auto task93 = make_shared<Task93>(tensor93, pindex, this->e0_);
  task92->add_dep(task93);
  task93->add_dep(task70);
  residualq->add_task(task93);

  vector<IndexRange> I18_index = {closed_, closed_};
  auto I18 = make_shared<Tensor>(I18_index);
  vector<shared_ptr<Tensor>> tensor94 = {I17, t2, I18};
  auto task94 = make_shared<Task94>(tensor94, pindex);
  task92->add_dep(task94);
  task94->add_dep(task70);
  residualq->add_task(task94);

  vector<shared_ptr<Tensor>> tensor95 = {I18, f1_};
  auto task95 = make_shared<Task95>(tensor95, pindex);
  task94->add_dep(task95);
  task95->add_dep(task70);
  residualq->add_task(task95);

  vector<IndexRange> I30_index = {virt_, active_};
  auto I30 = make_shared<Tensor>(I30_index);
  vector<shared_ptr<Tensor>> tensor96 = {I17, t2, I30};
  auto task96 = make_shared<Task96>(tensor96, pindex);
  task92->add_dep(task96);
  task96->add_dep(task70);
  residualq->add_task(task96);

  vector<shared_ptr<Tensor>> tensor97 = {I30, f1_};
  auto task97 = make_shared<Task97>(tensor97, pindex);
  task96->add_dep(task97);
  task97->add_dep(task70);
  residualq->add_task(task97);

  vector<IndexRange> I20_index = {closed_, active_};
  auto I20 = make_shared<Tensor>(I20_index);
  vector<shared_ptr<Tensor>> tensor98 = {I11, Gamma7_(), I20};
  auto task98 = make_shared<Task98>(tensor98, pindex);
  task86->add_dep(task98);
  task98->add_dep(task70);
  residualq->add_task(task98);

  vector<shared_ptr<Tensor>> tensor99 = {I20, h1_};
  auto task99 = make_shared<Task99>(tensor99, pindex);
  task98->add_dep(task99);
  task99->add_dep(task70);
  residualq->add_task(task99);

  vector<IndexRange> I21_index = {virt_, closed_};
  auto I21 = make_shared<Tensor>(I21_index);
  vector<shared_ptr<Tensor>> tensor100 = {I20, t2, I21};
  auto task100 = make_shared<Task100>(tensor100, pindex);
  task98->add_dep(task100);
  task100->add_dep(task70);
  residualq->add_task(task100);

  vector<shared_ptr<Tensor>> tensor101 = {I21, f1_};
  auto task101 = make_shared<Task101>(tensor101, pindex);
  task100->add_dep(task101);
  task101->add_dep(task70);
  residualq->add_task(task101);

  vector<IndexRange> I24_index = {virt_, closed_};
  auto I24 = make_shared<Tensor>(I24_index);
  vector<shared_ptr<Tensor>> tensor102 = {I20, t2, I24};
  auto task102 = make_shared<Task102>(tensor102, pindex);
  task98->add_dep(task102);
  task102->add_dep(task70);
  residualq->add_task(task102);

  vector<shared_ptr<Tensor>> tensor103 = {I24, f1_};
  auto task103 = make_shared<Task103>(tensor103, pindex);
  task102->add_dep(task103);
  task103->add_dep(task70);
  residualq->add_task(task103);

  vector<IndexRange> I26_index = {active_, active_, closed_, active_};
  auto I26 = make_shared<Tensor>(I26_index);
  vector<shared_ptr<Tensor>> tensor104 = {I11, Gamma9_(), I26};
  auto task104 = make_shared<Task104>(tensor104, pindex);
  task86->add_dep(task104);
  task104->add_dep(task70);
  residualq->add_task(task104);

  vector<IndexRange> I27_index = {virt_, active_};
  auto I27 = make_shared<Tensor>(I27_index);
  vector<shared_ptr<Tensor>> tensor105 = {I26, t2, I27};
  auto task105 = make_shared<Task105>(tensor105, pindex);
  task104->add_dep(task105);
  task105->add_dep(task70);
  residualq->add_task(task105);

  vector<shared_ptr<Tensor>> tensor106 = {I27, f1_};
  auto task106 = make_shared<Task106>(tensor106, pindex);
  task105->add_dep(task106);
  task106->add_dep(task70);
  residualq->add_task(task106);

  vector<IndexRange> I305_index = {closed_, active_, active_, active_};
  auto I305 = make_shared<Tensor>(I305_index);
  vector<shared_ptr<Tensor>> tensor107 = {I11, Gamma107_(), I305};
  auto task107 = make_shared<Task107>(tensor107, pindex);
  task86->add_dep(task107);
  task107->add_dep(task70);
  residualq->add_task(task107);

  vector<shared_ptr<Tensor>> tensor108 = {I305, v2_};
  auto task108 = make_shared<Task108>(tensor108, pindex);
  task107->add_dep(task108);
  task108->add_dep(task70);
  residualq->add_task(task108);

  vector<IndexRange> I31_index = {closed_, closed_, active_, virt_};
  auto I31 = make_shared<Tensor>(I31_index);
  vector<shared_ptr<Tensor>> tensor109 = {r, I31};
  auto task109 = make_shared<Task109>(tensor109, pindex);
  task109->add_dep(task70);
  residualq->add_task(task109);

  vector<IndexRange> I32_index = {closed_, closed_, active_, active_};
  auto I32 = make_shared<Tensor>(I32_index);
  vector<shared_ptr<Tensor>> tensor110 = {I31, f1_, I32};
  auto task110 = make_shared<Task110>(tensor110, pindex);
  task109->add_dep(task110);
  task110->add_dep(task70);
  residualq->add_task(task110);

  vector<IndexRange> I33_index = {closed_, active_, closed_, active_};
  auto I33 = make_shared<Tensor>(I33_index);
  vector<shared_ptr<Tensor>> tensor111 = {I32, Gamma3_(), I33};
  auto task111 = make_shared<Task111>(tensor111, pindex);
  task110->add_dep(task111);
  task111->add_dep(task70);
  residualq->add_task(task111);

  vector<shared_ptr<Tensor>> tensor112 = {I33, t2};
  auto task112 = make_shared<Task112>(tensor112, pindex);
  task111->add_dep(task112);
  task112->add_dep(task70);
  residualq->add_task(task112);

  vector<IndexRange> I35_index = {closed_, active_};
  auto I35 = make_shared<Tensor>(I35_index);
  vector<shared_ptr<Tensor>> tensor113 = {I31, f1_, I35};
  auto task113 = make_shared<Task113>(tensor113, pindex);
  task109->add_dep(task113);
  task113->add_dep(task70);
  residualq->add_task(task113);

  vector<IndexRange> I36_index = {active_, active_, closed_, active_};
  auto I36 = make_shared<Tensor>(I36_index);
  vector<shared_ptr<Tensor>> tensor114 = {I35, Gamma12_(), I36};
  auto task114 = make_shared<Task114>(tensor114, pindex);
  task113->add_dep(task114);
  task114->add_dep(task70);
  residualq->add_task(task114);

  vector<shared_ptr<Tensor>> tensor115 = {I36, t2};
  auto task115 = make_shared<Task115>(tensor115, pindex);
  task114->add_dep(task115);
  task115->add_dep(task70);
  residualq->add_task(task115);

  vector<IndexRange> I38_index = {closed_, active_};
  auto I38 = make_shared<Tensor>(I38_index);
  vector<shared_ptr<Tensor>> tensor116 = {I31, f1_, I38};
  auto task116 = make_shared<Task116>(tensor116, pindex);
  task109->add_dep(task116);
  task116->add_dep(task70);
  residualq->add_task(task116);

  vector<IndexRange> I39_index = {active_, active_, closed_, active_};
  auto I39 = make_shared<Tensor>(I39_index);
  vector<shared_ptr<Tensor>> tensor117 = {I38, Gamma12_(), I39};
  auto task117 = make_shared<Task117>(tensor117, pindex);
  task116->add_dep(task117);
  task117->add_dep(task70);
  residualq->add_task(task117);

  vector<shared_ptr<Tensor>> tensor118 = {I39, t2};
  auto task118 = make_shared<Task118>(tensor118, pindex);
  task117->add_dep(task118);
  task118->add_dep(task70);
  residualq->add_task(task118);

  vector<IndexRange> I41_index = {closed_, virt_, closed_, active_};
  auto I41 = make_shared<Tensor>(I41_index);
  vector<shared_ptr<Tensor>> tensor119 = {I31, Gamma14_(), I41};
  auto task119 = make_shared<Task119>(tensor119, pindex);
  task109->add_dep(task119);
  task119->add_dep(task70);
  residualq->add_task(task119);

  vector<shared_ptr<Tensor>> tensor120 = {I41, t2};
  auto task120 = make_shared<Task120>(tensor120, pindex);
  task119->add_dep(task120);
  task120->add_dep(task70);
  residualq->add_task(task120);

  vector<IndexRange> I45_index = {closed_, virt_, closed_, active_};
  auto I45 = make_shared<Tensor>(I45_index);
  vector<shared_ptr<Tensor>> tensor121 = {I31, Gamma16_(), I45};
  auto task121 = make_shared<Task121>(tensor121, pindex);
  task109->add_dep(task121);
  task121->add_dep(task70);
  residualq->add_task(task121);

  vector<shared_ptr<Tensor>> tensor122 = {I45, t2, v2_};
  auto task122 = make_shared<Task122>(tensor122, pindex, this->e0_);
  task121->add_dep(task122);
  task122->add_dep(task70);
  residualq->add_task(task122);

  vector<IndexRange> I46_index = {closed_, closed_};
  auto I46 = make_shared<Tensor>(I46_index);
  vector<shared_ptr<Tensor>> tensor123 = {I45, t2, I46};
  auto task123 = make_shared<Task123>(tensor123, pindex);
  task121->add_dep(task123);
  task123->add_dep(task70);
  residualq->add_task(task123);

  vector<shared_ptr<Tensor>> tensor124 = {I46, f1_};
  auto task124 = make_shared<Task124>(tensor124, pindex);
  task123->add_dep(task124);
  task124->add_dep(task70);
  residualq->add_task(task124);

  vector<IndexRange> I49_index = {closed_, closed_};
  auto I49 = make_shared<Tensor>(I49_index);
  vector<shared_ptr<Tensor>> tensor125 = {I45, t2, I49};
  auto task125 = make_shared<Task125>(tensor125, pindex);
  task121->add_dep(task125);
  task125->add_dep(task70);
  residualq->add_task(task125);

  vector<shared_ptr<Tensor>> tensor126 = {I49, f1_};
  auto task126 = make_shared<Task126>(tensor126, pindex);
  task125->add_dep(task126);
  task126->add_dep(task70);
  residualq->add_task(task126);

  vector<IndexRange> I52_index = {closed_, closed_};
  auto I52 = make_shared<Tensor>(I52_index);
  vector<shared_ptr<Tensor>> tensor127 = {I45, t2, I52};
  auto task127 = make_shared<Task127>(tensor127, pindex);
  task121->add_dep(task127);
  task127->add_dep(task70);
  residualq->add_task(task127);

  vector<shared_ptr<Tensor>> tensor128 = {I52, f1_};
  auto task128 = make_shared<Task128>(tensor128, pindex);
  task127->add_dep(task128);
  task128->add_dep(task70);
  residualq->add_task(task128);

  vector<IndexRange> I55_index = {virt_, virt_};
  auto I55 = make_shared<Tensor>(I55_index);
  vector<shared_ptr<Tensor>> tensor129 = {I45, t2, I55};
  auto task129 = make_shared<Task129>(tensor129, pindex);
  task121->add_dep(task129);
  task129->add_dep(task70);
  residualq->add_task(task129);

  vector<shared_ptr<Tensor>> tensor130 = {I55, f1_};
  auto task130 = make_shared<Task130>(tensor130, pindex);
  task129->add_dep(task130);
  task130->add_dep(task70);
  residualq->add_task(task130);

  vector<IndexRange> I58_index = {closed_, closed_};
  auto I58 = make_shared<Tensor>(I58_index);
  vector<shared_ptr<Tensor>> tensor131 = {I45, t2, I58};
  auto task131 = make_shared<Task131>(tensor131, pindex);
  task121->add_dep(task131);
  task131->add_dep(task70);
  residualq->add_task(task131);

  vector<shared_ptr<Tensor>> tensor132 = {I58, f1_};
  auto task132 = make_shared<Task132>(tensor132, pindex);
  task131->add_dep(task132);
  task132->add_dep(task70);
  residualq->add_task(task132);

  vector<IndexRange> I61_index = {virt_, virt_};
  auto I61 = make_shared<Tensor>(I61_index);
  vector<shared_ptr<Tensor>> tensor133 = {I45, t2, I61};
  auto task133 = make_shared<Task133>(tensor133, pindex);
  task121->add_dep(task133);
  task133->add_dep(task70);
  residualq->add_task(task133);

  vector<shared_ptr<Tensor>> tensor134 = {I61, f1_};
  auto task134 = make_shared<Task134>(tensor134, pindex);
  task133->add_dep(task134);
  task134->add_dep(task70);
  residualq->add_task(task134);

  vector<IndexRange> I63_index = {virt_, closed_, active_, active_};
  auto I63 = make_shared<Tensor>(I63_index);
  vector<shared_ptr<Tensor>> tensor135 = {I31, f1_, I63};
  auto task135 = make_shared<Task135>(tensor135, pindex);
  task109->add_dep(task135);
  task135->add_dep(task70);
  residualq->add_task(task135);

  vector<IndexRange> I64_index = {active_, virt_, closed_, active_};
  auto I64 = make_shared<Tensor>(I64_index);
  vector<shared_ptr<Tensor>> tensor136 = {I63, Gamma22_(), I64};
  auto task136 = make_shared<Task136>(tensor136, pindex);
  task135->add_dep(task136);
  task136->add_dep(task70);
  residualq->add_task(task136);

  vector<shared_ptr<Tensor>> tensor137 = {I64, t2};
  auto task137 = make_shared<Task137>(tensor137, pindex);
  task136->add_dep(task137);
  task137->add_dep(task70);
  residualq->add_task(task137);

  vector<IndexRange> I70_index = {closed_, virt_, active_, active_};
  auto I70 = make_shared<Tensor>(I70_index);
  vector<shared_ptr<Tensor>> tensor138 = {I63, Gamma12_(), I70};
  auto task138 = make_shared<Task138>(tensor138, pindex);
  task135->add_dep(task138);
  task138->add_dep(task70);
  residualq->add_task(task138);

  vector<shared_ptr<Tensor>> tensor139 = {I70, t2};
  auto task139 = make_shared<Task139>(tensor139, pindex);
  task138->add_dep(task139);
  task139->add_dep(task70);
  residualq->add_task(task139);

  vector<IndexRange> I66_index = {virt_, closed_, active_, active_};
  auto I66 = make_shared<Tensor>(I66_index);
  vector<shared_ptr<Tensor>> tensor140 = {I31, f1_, I66};
  auto task140 = make_shared<Task140>(tensor140, pindex);
  task109->add_dep(task140);
  task140->add_dep(task70);
  residualq->add_task(task140);

  vector<IndexRange> I67_index = {active_, virt_, closed_, active_};
  auto I67 = make_shared<Tensor>(I67_index);
  vector<shared_ptr<Tensor>> tensor141 = {I66, Gamma12_(), I67};
  auto task141 = make_shared<Task141>(tensor141, pindex);
  task140->add_dep(task141);
  task141->add_dep(task70);
  residualq->add_task(task141);

  vector<shared_ptr<Tensor>> tensor142 = {I67, t2};
  auto task142 = make_shared<Task142>(tensor142, pindex);
  task141->add_dep(task142);
  task142->add_dep(task70);
  residualq->add_task(task142);

  vector<IndexRange> I75_index = {virt_, active_};
  auto I75 = make_shared<Tensor>(I75_index);
  vector<shared_ptr<Tensor>> tensor143 = {I31, t2, I75};
  auto task143 = make_shared<Task143>(tensor143, pindex);
  task109->add_dep(task143);
  task143->add_dep(task70);
  residualq->add_task(task143);

  vector<IndexRange> I76_index = {virt_, active_};
  auto I76 = make_shared<Tensor>(I76_index);
  vector<shared_ptr<Tensor>> tensor144 = {I75, Gamma16_(), I76};
  auto task144 = make_shared<Task144>(tensor144, pindex);
  task143->add_dep(task144);
  task144->add_dep(task70);
  residualq->add_task(task144);

  vector<shared_ptr<Tensor>> tensor145 = {I76, f1_};
  auto task145 = make_shared<Task145>(tensor145, pindex);
  task144->add_dep(task145);
  task145->add_dep(task70);
  residualq->add_task(task145);

  vector<IndexRange> I78_index = {virt_, active_};
  auto I78 = make_shared<Tensor>(I78_index);
  vector<shared_ptr<Tensor>> tensor146 = {I31, t2, I78};
  auto task146 = make_shared<Task146>(tensor146, pindex);
  task109->add_dep(task146);
  task146->add_dep(task70);
  residualq->add_task(task146);

  vector<IndexRange> I79_index = {virt_, active_};
  auto I79 = make_shared<Tensor>(I79_index);
  vector<shared_ptr<Tensor>> tensor147 = {I78, Gamma16_(), I79};
  auto task147 = make_shared<Task147>(tensor147, pindex);
  task146->add_dep(task147);
  task147->add_dep(task70);
  residualq->add_task(task147);

  vector<shared_ptr<Tensor>> tensor148 = {I79, f1_};
  auto task148 = make_shared<Task148>(tensor148, pindex);
  task147->add_dep(task148);
  task148->add_dep(task70);
  residualq->add_task(task148);

  vector<IndexRange> I80_index = {closed_, active_, active_, virt_};
  auto I80 = make_shared<Tensor>(I80_index);
  vector<shared_ptr<Tensor>> tensor149 = {r, I80};
  auto task149 = make_shared<Task149>(tensor149, pindex);
  task149->add_dep(task70);
  residualq->add_task(task149);

  vector<IndexRange> I81_index = {closed_, active_, active_, active_};
  auto I81 = make_shared<Tensor>(I81_index);
  vector<shared_ptr<Tensor>> tensor150 = {I80, f1_, I81};
  auto task150 = make_shared<Task150>(tensor150, pindex);
  task149->add_dep(task150);
  task150->add_dep(task70);
  residualq->add_task(task150);

  vector<IndexRange> I82_index = {active_, active_, closed_, active_};
  auto I82 = make_shared<Tensor>(I82_index);
  vector<shared_ptr<Tensor>> tensor151 = {I81, Gamma28_(), I82};
  auto task151 = make_shared<Task151>(tensor151, pindex);
  task150->add_dep(task151);
  task151->add_dep(task70);
  residualq->add_task(task151);

  vector<shared_ptr<Tensor>> tensor152 = {I82, t2};
  auto task152 = make_shared<Task152>(tensor152, pindex);
  task151->add_dep(task152);
  task152->add_dep(task70);
  residualq->add_task(task152);

  vector<IndexRange> I84_index = {active_, virt_, closed_, active_};
  auto I84 = make_shared<Tensor>(I84_index);
  vector<shared_ptr<Tensor>> tensor153 = {I80, Gamma29_(), I84};
  auto task153 = make_shared<Task153>(tensor153, pindex);
  task149->add_dep(task153);
  task153->add_dep(task70);
  residualq->add_task(task153);

  vector<shared_ptr<Tensor>> tensor154 = {I84, v2_};
  auto task154 = make_shared<Task154>(tensor154, pindex);
  task153->add_dep(task154);
  task154->add_dep(task70);
  residualq->add_task(task154);

  vector<IndexRange> I85_index = {active_, closed_};
  auto I85 = make_shared<Tensor>(I85_index);
  vector<shared_ptr<Tensor>> tensor155 = {I84, t2, I85};
  auto task155 = make_shared<Task155>(tensor155, pindex);
  task153->add_dep(task155);
  task155->add_dep(task70);
  residualq->add_task(task155);

  vector<shared_ptr<Tensor>> tensor156 = {I85, f1_};
  auto task156 = make_shared<Task156>(tensor156, pindex);
  task155->add_dep(task156);
  task156->add_dep(task70);
  residualq->add_task(task156);

  vector<IndexRange> I87_index = {active_, closed_, virt_, active_};
  auto I87 = make_shared<Tensor>(I87_index);
  vector<shared_ptr<Tensor>> tensor157 = {I80, Gamma7_(), I87};
  auto task157 = make_shared<Task157>(tensor157, pindex);
  task149->add_dep(task157);
  task157->add_dep(task70);
  residualq->add_task(task157);

  vector<IndexRange> I88_index = {active_, closed_};
  auto I88 = make_shared<Tensor>(I88_index);
  vector<shared_ptr<Tensor>> tensor158 = {I87, t2, I88};
  auto task158 = make_shared<Task158>(tensor158, pindex);
  task157->add_dep(task158);
  task158->add_dep(task70);
  residualq->add_task(task158);

  vector<shared_ptr<Tensor>> tensor159 = {I88, f1_};
  auto task159 = make_shared<Task159>(tensor159, pindex);
  task158->add_dep(task159);
  task159->add_dep(task70);
  residualq->add_task(task159);

  vector<IndexRange> I90_index = {active_, virt_, closed_, active_};
  auto I90 = make_shared<Tensor>(I90_index);
  vector<shared_ptr<Tensor>> tensor160 = {I80, Gamma31_(), I90};
  auto task160 = make_shared<Task160>(tensor160, pindex);
  task149->add_dep(task160);
  task160->add_dep(task70);
  residualq->add_task(task160);

  vector<shared_ptr<Tensor>> tensor161 = {I90, t2};
  auto task161 = make_shared<Task161>(tensor161, pindex);
  task160->add_dep(task161);
  task161->add_dep(task70);
  residualq->add_task(task161);

  vector<IndexRange> I92_index = {closed_, active_, virt_, active_};
  auto I92 = make_shared<Tensor>(I92_index);
  vector<shared_ptr<Tensor>> tensor162 = {I80, Gamma32_(), I92};
  auto task162 = make_shared<Task162>(tensor162, pindex);
  task149->add_dep(task162);
  task162->add_dep(task70);
  residualq->add_task(task162);

  vector<shared_ptr<Tensor>> tensor163 = {I92, t2, v2_};
  auto task163 = make_shared<Task163>(tensor163, pindex, this->e0_);
  task162->add_dep(task163);
  task163->add_dep(task70);
  residualq->add_task(task163);

  vector<IndexRange> I93_index = {closed_, closed_};
  auto I93 = make_shared<Tensor>(I93_index);
  vector<shared_ptr<Tensor>> tensor164 = {I92, t2, I93};
  auto task164 = make_shared<Task164>(tensor164, pindex);
  task162->add_dep(task164);
  task164->add_dep(task70);
  residualq->add_task(task164);

  vector<shared_ptr<Tensor>> tensor165 = {I93, f1_};
  auto task165 = make_shared<Task165>(tensor165, pindex);
  task164->add_dep(task165);
  task165->add_dep(task70);
  residualq->add_task(task165);

  vector<IndexRange> I96_index = {virt_, virt_};
  auto I96 = make_shared<Tensor>(I96_index);
  vector<shared_ptr<Tensor>> tensor166 = {I92, t2, I96};
  auto task166 = make_shared<Task166>(tensor166, pindex);
  task162->add_dep(task166);
  task166->add_dep(task70);
  residualq->add_task(task166);

  vector<shared_ptr<Tensor>> tensor167 = {I96, f1_};
  auto task167 = make_shared<Task167>(tensor167, pindex);
  task166->add_dep(task167);
  task167->add_dep(task70);
  residualq->add_task(task167);

  vector<IndexRange> I119_index = {virt_, active_};
  auto I119 = make_shared<Tensor>(I119_index);
  vector<shared_ptr<Tensor>> tensor168 = {I92, t2, I119};
  auto task168 = make_shared<Task168>(tensor168, pindex);
  task162->add_dep(task168);
  task168->add_dep(task70);
  residualq->add_task(task168);

  vector<shared_ptr<Tensor>> tensor169 = {I119, f1_};
  auto task169 = make_shared<Task169>(tensor169, pindex);
  task168->add_dep(task169);
  task169->add_dep(task70);
  residualq->add_task(task169);

  vector<IndexRange> I98_index = {closed_, virt_, active_, active_};
  auto I98 = make_shared<Tensor>(I98_index);
  vector<shared_ptr<Tensor>> tensor170 = {I80, Gamma34_(), I98};
  auto task170 = make_shared<Task170>(tensor170, pindex);
  task149->add_dep(task170);
  task170->add_dep(task70);
  residualq->add_task(task170);

  vector<shared_ptr<Tensor>> tensor171 = {I98, t2};
  auto task171 = make_shared<Task171>(tensor171, pindex);
  task170->add_dep(task171);
  task171->add_dep(task70);
  residualq->add_task(task171);

  vector<IndexRange> I100_index = {closed_, virt_, active_, active_};
  auto I100 = make_shared<Tensor>(I100_index);
  vector<shared_ptr<Tensor>> tensor172 = {I80, Gamma35_(), I100};
  auto task172 = make_shared<Task172>(tensor172, pindex);
  task149->add_dep(task172);
  task172->add_dep(task70);
  residualq->add_task(task172);

  vector<shared_ptr<Tensor>> tensor173 = {I100, t2, v2_};
  auto task173 = make_shared<Task173>(tensor173, pindex, this->e0_);
  task172->add_dep(task173);
  task173->add_dep(task70);
  residualq->add_task(task173);

  vector<IndexRange> I101_index = {closed_, closed_};
  auto I101 = make_shared<Tensor>(I101_index);
  vector<shared_ptr<Tensor>> tensor174 = {I100, t2, I101};
  auto task174 = make_shared<Task174>(tensor174, pindex);
  task172->add_dep(task174);
  task174->add_dep(task70);
  residualq->add_task(task174);

  vector<shared_ptr<Tensor>> tensor175 = {I101, f1_};
  auto task175 = make_shared<Task175>(tensor175, pindex);
  task174->add_dep(task175);
  task175->add_dep(task70);
  residualq->add_task(task175);

  vector<IndexRange> I104_index = {virt_, virt_};
  auto I104 = make_shared<Tensor>(I104_index);
  vector<shared_ptr<Tensor>> tensor176 = {I100, t2, I104};
  auto task176 = make_shared<Task176>(tensor176, pindex);
  task172->add_dep(task176);
  task176->add_dep(task70);
  residualq->add_task(task176);

  vector<shared_ptr<Tensor>> tensor177 = {I104, f1_};
  auto task177 = make_shared<Task177>(tensor177, pindex);
  task176->add_dep(task177);
  task177->add_dep(task70);
  residualq->add_task(task177);

  vector<IndexRange> I116_index = {virt_, active_};
  auto I116 = make_shared<Tensor>(I116_index);
  vector<shared_ptr<Tensor>> tensor178 = {I100, t2, I116};
  auto task178 = make_shared<Task178>(tensor178, pindex);
  task172->add_dep(task178);
  task178->add_dep(task70);
  residualq->add_task(task178);

  vector<shared_ptr<Tensor>> tensor179 = {I116, f1_};
  auto task179 = make_shared<Task179>(tensor179, pindex);
  task178->add_dep(task179);
  task179->add_dep(task70);
  residualq->add_task(task179);

  vector<IndexRange> I106_index = {virt_, active_, active_, active_};
  auto I106 = make_shared<Tensor>(I106_index);
  vector<shared_ptr<Tensor>> tensor180 = {I80, f1_, I106};
  auto task180 = make_shared<Task180>(tensor180, pindex);
  task149->add_dep(task180);
  task180->add_dep(task70);
  residualq->add_task(task180);

  vector<IndexRange> I107_index = {active_, virt_, active_, active_};
  auto I107 = make_shared<Tensor>(I107_index);
  vector<shared_ptr<Tensor>> tensor181 = {I106, Gamma37_(), I107};
  auto task181 = make_shared<Task181>(tensor181, pindex);
  task180->add_dep(task181);
  task181->add_dep(task70);
  residualq->add_task(task181);

  vector<shared_ptr<Tensor>> tensor182 = {I107, t2};
  auto task182 = make_shared<Task182>(tensor182, pindex);
  task181->add_dep(task182);
  task182->add_dep(task70);
  residualq->add_task(task182);

  vector<IndexRange> I109_index = {closed_, virt_};
  auto I109 = make_shared<Tensor>(I109_index);
  vector<shared_ptr<Tensor>> tensor183 = {I80, Gamma38_(), I109};
  auto task183 = make_shared<Task183>(tensor183, pindex);
  task149->add_dep(task183);
  task183->add_dep(task70);
  residualq->add_task(task183);

  vector<shared_ptr<Tensor>> tensor184 = {I109, h1_};
  auto task184 = make_shared<Task184>(tensor184, pindex);
  task183->add_dep(task184);
  task184->add_dep(task70);
  residualq->add_task(task184);

  vector<IndexRange> I110_index = {virt_, closed_};
  auto I110 = make_shared<Tensor>(I110_index);
  vector<shared_ptr<Tensor>> tensor185 = {I109, t2, I110};
  auto task185 = make_shared<Task185>(tensor185, pindex);
  task183->add_dep(task185);
  task185->add_dep(task70);
  residualq->add_task(task185);

  vector<shared_ptr<Tensor>> tensor186 = {I110, f1_};
  auto task186 = make_shared<Task186>(tensor186, pindex);
  task185->add_dep(task186);
  task186->add_dep(task70);
  residualq->add_task(task186);

  vector<IndexRange> I113_index = {virt_, closed_};
  auto I113 = make_shared<Tensor>(I113_index);
  vector<shared_ptr<Tensor>> tensor187 = {I109, t2, I113};
  auto task187 = make_shared<Task187>(tensor187, pindex);
  task183->add_dep(task187);
  task187->add_dep(task70);
  residualq->add_task(task187);

  vector<shared_ptr<Tensor>> tensor188 = {I113, f1_};
  auto task188 = make_shared<Task188>(tensor188, pindex);
  task187->add_dep(task188);
  task188->add_dep(task70);
  residualq->add_task(task188);

  vector<IndexRange> I120_index = {closed_, active_, active_, virt_};
  auto I120 = make_shared<Tensor>(I120_index);
  vector<shared_ptr<Tensor>> tensor189 = {r, I120};
  auto task189 = make_shared<Task189>(tensor189, pindex);
  task189->add_dep(task70);
  residualq->add_task(task189);

  vector<IndexRange> I121_index = {closed_, active_, active_, active_};
  auto I121 = make_shared<Tensor>(I121_index);
  vector<shared_ptr<Tensor>> tensor190 = {I120, f1_, I121};
  auto task190 = make_shared<Task190>(tensor190, pindex);
  task189->add_dep(task190);
  task190->add_dep(task70);
  residualq->add_task(task190);

  vector<IndexRange> I122_index = {active_, active_, closed_, active_};
  auto I122 = make_shared<Tensor>(I122_index);
  vector<shared_ptr<Tensor>> tensor191 = {I121, Gamma6_(), I122};
  auto task191 = make_shared<Task191>(tensor191, pindex);
  task190->add_dep(task191);
  task191->add_dep(task70);
  residualq->add_task(task191);

  vector<shared_ptr<Tensor>> tensor192 = {I122, t2};
  auto task192 = make_shared<Task192>(tensor192, pindex);
  task191->add_dep(task192);
  task192->add_dep(task70);
  residualq->add_task(task192);

  vector<IndexRange> I124_index = {active_, virt_, closed_, active_};
  auto I124 = make_shared<Tensor>(I124_index);
  vector<shared_ptr<Tensor>> tensor193 = {I120, Gamma7_(), I124};
  auto task193 = make_shared<Task193>(tensor193, pindex);
  task189->add_dep(task193);
  task193->add_dep(task70);
  residualq->add_task(task193);

  vector<shared_ptr<Tensor>> tensor194 = {I124, v2_};
  auto task194 = make_shared<Task194>(tensor194, pindex);
  task193->add_dep(task194);
  task194->add_dep(task70);
  residualq->add_task(task194);

  vector<IndexRange> I125_index = {active_, closed_};
  auto I125 = make_shared<Tensor>(I125_index);
  vector<shared_ptr<Tensor>> tensor195 = {I124, t2, I125};
  auto task195 = make_shared<Task195>(tensor195, pindex);
  task193->add_dep(task195);
  task195->add_dep(task70);
  residualq->add_task(task195);

  vector<shared_ptr<Tensor>> tensor196 = {I125, f1_};
  auto task196 = make_shared<Task196>(tensor196, pindex);
  task195->add_dep(task196);
  task196->add_dep(task70);
  residualq->add_task(task196);

  vector<IndexRange> I128_index = {active_, closed_};
  auto I128 = make_shared<Tensor>(I128_index);
  vector<shared_ptr<Tensor>> tensor197 = {I124, t2, I128};
  auto task197 = make_shared<Task197>(tensor197, pindex);
  task193->add_dep(task197);
  task197->add_dep(task70);
  residualq->add_task(task197);

  vector<shared_ptr<Tensor>> tensor198 = {I128, f1_};
  auto task198 = make_shared<Task198>(tensor198, pindex);
  task197->add_dep(task198);
  task198->add_dep(task70);
  residualq->add_task(task198);

  vector<IndexRange> I130_index = {active_, virt_, closed_, active_};
  auto I130 = make_shared<Tensor>(I130_index);
  vector<shared_ptr<Tensor>> tensor199 = {I120, Gamma34_(), I130};
  auto task199 = make_shared<Task199>(tensor199, pindex);
  task189->add_dep(task199);
  task199->add_dep(task70);
  residualq->add_task(task199);

  vector<shared_ptr<Tensor>> tensor200 = {I130, t2};
  auto task200 = make_shared<Task200>(tensor200, pindex);
  task199->add_dep(task200);
  task200->add_dep(task70);
  residualq->add_task(task200);

  vector<IndexRange> I132_index = {closed_, active_, virt_, active_};
  auto I132 = make_shared<Tensor>(I132_index);
  vector<shared_ptr<Tensor>> tensor201 = {I120, Gamma35_(), I132};
  auto task201 = make_shared<Task201>(tensor201, pindex);
  task189->add_dep(task201);
  task201->add_dep(task70);
  residualq->add_task(task201);

  vector<shared_ptr<Tensor>> tensor202 = {I132, t2, v2_};
  auto task202 = make_shared<Task202>(tensor202, pindex, this->e0_);
  task201->add_dep(task202);
  task202->add_dep(task70);
  residualq->add_task(task202);

  vector<IndexRange> I133_index = {closed_, closed_};
  auto I133 = make_shared<Tensor>(I133_index);
  vector<shared_ptr<Tensor>> tensor203 = {I132, t2, I133};
  auto task203 = make_shared<Task203>(tensor203, pindex);
  task201->add_dep(task203);
  task203->add_dep(task70);
  residualq->add_task(task203);

  vector<shared_ptr<Tensor>> tensor204 = {I133, f1_};
  auto task204 = make_shared<Task204>(tensor204, pindex);
  task203->add_dep(task204);
  task204->add_dep(task70);
  residualq->add_task(task204);

  vector<IndexRange> I136_index = {virt_, virt_};
  auto I136 = make_shared<Tensor>(I136_index);
  vector<shared_ptr<Tensor>> tensor205 = {I132, t2, I136};
  auto task205 = make_shared<Task205>(tensor205, pindex);
  task201->add_dep(task205);
  task205->add_dep(task70);
  residualq->add_task(task205);

  vector<shared_ptr<Tensor>> tensor206 = {I136, f1_};
  auto task206 = make_shared<Task206>(tensor206, pindex);
  task205->add_dep(task206);
  task206->add_dep(task70);
  residualq->add_task(task206);

  vector<IndexRange> I141_index = {closed_, closed_};
  auto I141 = make_shared<Tensor>(I141_index);
  vector<shared_ptr<Tensor>> tensor207 = {I132, t2, I141};
  auto task207 = make_shared<Task207>(tensor207, pindex);
  task201->add_dep(task207);
  task207->add_dep(task70);
  residualq->add_task(task207);

  vector<shared_ptr<Tensor>> tensor208 = {I141, f1_};
  auto task208 = make_shared<Task208>(tensor208, pindex);
  task207->add_dep(task208);
  task208->add_dep(task70);
  residualq->add_task(task208);

  vector<IndexRange> I144_index = {virt_, virt_};
  auto I144 = make_shared<Tensor>(I144_index);
  vector<shared_ptr<Tensor>> tensor209 = {I132, t2, I144};
  auto task209 = make_shared<Task209>(tensor209, pindex);
  task201->add_dep(task209);
  task209->add_dep(task70);
  residualq->add_task(task209);

  vector<shared_ptr<Tensor>> tensor210 = {I144, f1_};
  auto task210 = make_shared<Task210>(tensor210, pindex);
  task209->add_dep(task210);
  task210->add_dep(task70);
  residualq->add_task(task210);

  vector<IndexRange> I156_index = {virt_, active_};
  auto I156 = make_shared<Tensor>(I156_index);
  vector<shared_ptr<Tensor>> tensor211 = {I132, t2, I156};
  auto task211 = make_shared<Task211>(tensor211, pindex);
  task201->add_dep(task211);
  task211->add_dep(task70);
  residualq->add_task(task211);

  vector<shared_ptr<Tensor>> tensor212 = {I156, f1_};
  auto task212 = make_shared<Task212>(tensor212, pindex);
  task211->add_dep(task212);
  task212->add_dep(task70);
  residualq->add_task(task212);

  vector<IndexRange> I159_index = {virt_, active_};
  auto I159 = make_shared<Tensor>(I159_index);
  vector<shared_ptr<Tensor>> tensor213 = {I132, t2, I159};
  auto task213 = make_shared<Task213>(tensor213, pindex);
  task201->add_dep(task213);
  task213->add_dep(task70);
  residualq->add_task(task213);

  vector<shared_ptr<Tensor>> tensor214 = {I159, f1_};
  auto task214 = make_shared<Task214>(tensor214, pindex);
  task213->add_dep(task214);
  task214->add_dep(task70);
  residualq->add_task(task214);

  vector<IndexRange> I146_index = {virt_, active_, active_, active_};
  auto I146 = make_shared<Tensor>(I146_index);
  vector<shared_ptr<Tensor>> tensor215 = {I120, f1_, I146};
  auto task215 = make_shared<Task215>(tensor215, pindex);
  task189->add_dep(task215);
  task215->add_dep(task70);
  residualq->add_task(task215);

  vector<IndexRange> I147_index = {active_, virt_, active_, active_};
  auto I147 = make_shared<Tensor>(I147_index);
  vector<shared_ptr<Tensor>> tensor216 = {I146, Gamma51_(), I147};
  auto task216 = make_shared<Task216>(tensor216, pindex);
  task215->add_dep(task216);
  task216->add_dep(task70);
  residualq->add_task(task216);

  vector<shared_ptr<Tensor>> tensor217 = {I147, t2};
  auto task217 = make_shared<Task217>(tensor217, pindex);
  task216->add_dep(task217);
  task217->add_dep(task70);
  residualq->add_task(task217);

  vector<IndexRange> I149_index = {closed_, virt_};
  auto I149 = make_shared<Tensor>(I149_index);
  vector<shared_ptr<Tensor>> tensor218 = {I120, Gamma38_(), I149};
  auto task218 = make_shared<Task218>(tensor218, pindex);
  task189->add_dep(task218);
  task218->add_dep(task70);
  residualq->add_task(task218);

  vector<shared_ptr<Tensor>> tensor219 = {I149, h1_};
  auto task219 = make_shared<Task219>(tensor219, pindex);
  task218->add_dep(task219);
  task219->add_dep(task70);
  residualq->add_task(task219);

  vector<IndexRange> I150_index = {virt_, closed_};
  auto I150 = make_shared<Tensor>(I150_index);
  vector<shared_ptr<Tensor>> tensor220 = {I149, t2, I150};
  auto task220 = make_shared<Task220>(tensor220, pindex);
  task218->add_dep(task220);
  task220->add_dep(task70);
  residualq->add_task(task220);

  vector<shared_ptr<Tensor>> tensor221 = {I150, f1_};
  auto task221 = make_shared<Task221>(tensor221, pindex);
  task220->add_dep(task221);
  task221->add_dep(task70);
  residualq->add_task(task221);

  vector<IndexRange> I153_index = {virt_, closed_};
  auto I153 = make_shared<Tensor>(I153_index);
  vector<shared_ptr<Tensor>> tensor222 = {I149, t2, I153};
  auto task222 = make_shared<Task222>(tensor222, pindex);
  task218->add_dep(task222);
  task222->add_dep(task70);
  residualq->add_task(task222);

  vector<shared_ptr<Tensor>> tensor223 = {I153, f1_};
  auto task223 = make_shared<Task223>(tensor223, pindex);
  task222->add_dep(task223);
  task223->add_dep(task70);
  residualq->add_task(task223);

  vector<IndexRange> I160_index = {virt_, active_, active_, active_};
  auto I160 = make_shared<Tensor>(I160_index);
  vector<shared_ptr<Tensor>> tensor224 = {r, I160};
  auto task224 = make_shared<Task224>(tensor224, pindex);
  task224->add_dep(task70);
  residualq->add_task(task224);

  vector<IndexRange> I161_index = {active_, active_, virt_, active_};
  auto I161 = make_shared<Tensor>(I161_index);
  vector<shared_ptr<Tensor>> tensor225 = {I160, Gamma56_(), I161};
  auto task225 = make_shared<Task225>(tensor225, pindex);
  task224->add_dep(task225);
  task225->add_dep(task70);
  residualq->add_task(task225);

  vector<IndexRange> I162_index = {active_, closed_};
  auto I162 = make_shared<Tensor>(I162_index);
  vector<shared_ptr<Tensor>> tensor226 = {I161, t2, I162};
  auto task226 = make_shared<Task226>(tensor226, pindex);
  task225->add_dep(task226);
  task226->add_dep(task70);
  residualq->add_task(task226);

  vector<shared_ptr<Tensor>> tensor227 = {I162, f1_};
  auto task227 = make_shared<Task227>(tensor227, pindex);
  task226->add_dep(task227);
  task227->add_dep(task70);
  residualq->add_task(task227);

  vector<IndexRange> I164_index = {active_, virt_, active_, active_};
  auto I164 = make_shared<Tensor>(I164_index);
  vector<shared_ptr<Tensor>> tensor228 = {I160, Gamma57_(), I164};
  auto task228 = make_shared<Task228>(tensor228, pindex);
  task224->add_dep(task228);
  task228->add_dep(task70);
  residualq->add_task(task228);

  vector<shared_ptr<Tensor>> tensor229 = {I164, v2_};
  auto task229 = make_shared<Task229>(tensor229, pindex);
  task228->add_dep(task229);
  task229->add_dep(task70);
  residualq->add_task(task229);

  vector<IndexRange> I165_index = {active_, closed_};
  auto I165 = make_shared<Tensor>(I165_index);
  vector<shared_ptr<Tensor>> tensor230 = {I164, t2, I165};
  auto task230 = make_shared<Task230>(tensor230, pindex);
  task228->add_dep(task230);
  task230->add_dep(task70);
  residualq->add_task(task230);

  vector<shared_ptr<Tensor>> tensor231 = {I165, f1_};
  auto task231 = make_shared<Task231>(tensor231, pindex);
  task230->add_dep(task231);
  task231->add_dep(task70);
  residualq->add_task(task231);

  vector<IndexRange> I167_index = {active_, virt_, active_, active_};
  auto I167 = make_shared<Tensor>(I167_index);
  vector<shared_ptr<Tensor>> tensor232 = {I160, Gamma58_(), I167};
  auto task232 = make_shared<Task232>(tensor232, pindex);
  task224->add_dep(task232);
  task232->add_dep(task70);
  residualq->add_task(task232);

  vector<shared_ptr<Tensor>> tensor233 = {I167, t2};
  auto task233 = make_shared<Task233>(tensor233, pindex);
  task232->add_dep(task233);
  task233->add_dep(task70);
  residualq->add_task(task233);

  vector<IndexRange> I169_index = {virt_, active_, active_, active_};
  auto I169 = make_shared<Tensor>(I169_index);
  vector<shared_ptr<Tensor>> tensor234 = {I160, Gamma59_(), I169};
  auto task234 = make_shared<Task234>(tensor234, pindex);
  task224->add_dep(task234);
  task234->add_dep(task70);
  residualq->add_task(task234);

  vector<shared_ptr<Tensor>> tensor235 = {I169, t2, v2_};
  auto task235 = make_shared<Task235>(tensor235, pindex, this->e0_);
  task234->add_dep(task235);
  task235->add_dep(task70);
  residualq->add_task(task235);

  vector<IndexRange> I170_index = {virt_, virt_};
  auto I170 = make_shared<Tensor>(I170_index);
  vector<shared_ptr<Tensor>> tensor236 = {I169, t2, I170};
  auto task236 = make_shared<Task236>(tensor236, pindex);
  task234->add_dep(task236);
  task236->add_dep(task70);
  residualq->add_task(task236);

  vector<shared_ptr<Tensor>> tensor237 = {I170, f1_};
  auto task237 = make_shared<Task237>(tensor237, pindex);
  task236->add_dep(task237);
  task237->add_dep(task70);
  residualq->add_task(task237);

  vector<IndexRange> I179_index = {virt_, active_};
  auto I179 = make_shared<Tensor>(I179_index);
  vector<shared_ptr<Tensor>> tensor238 = {I169, t2, I179};
  auto task238 = make_shared<Task238>(tensor238, pindex);
  task234->add_dep(task238);
  task238->add_dep(task70);
  residualq->add_task(task238);

  vector<shared_ptr<Tensor>> tensor239 = {I179, f1_};
  auto task239 = make_shared<Task239>(tensor239, pindex);
  task238->add_dep(task239);
  task239->add_dep(task70);
  residualq->add_task(task239);

  vector<IndexRange> I172_index = {active_, virt_};
  auto I172 = make_shared<Tensor>(I172_index);
  vector<shared_ptr<Tensor>> tensor240 = {I160, Gamma60_(), I172};
  auto task240 = make_shared<Task240>(tensor240, pindex);
  task224->add_dep(task240);
  task240->add_dep(task70);
  residualq->add_task(task240);

  vector<shared_ptr<Tensor>> tensor241 = {I172, h1_};
  auto task241 = make_shared<Task241>(tensor241, pindex);
  task240->add_dep(task241);
  task241->add_dep(task70);
  residualq->add_task(task241);

  vector<IndexRange> I173_index = {virt_, closed_};
  auto I173 = make_shared<Tensor>(I173_index);
  vector<shared_ptr<Tensor>> tensor242 = {I172, t2, I173};
  auto task242 = make_shared<Task242>(tensor242, pindex);
  task240->add_dep(task242);
  task242->add_dep(task70);
  residualq->add_task(task242);

  vector<shared_ptr<Tensor>> tensor243 = {I173, f1_};
  auto task243 = make_shared<Task243>(tensor243, pindex);
  task242->add_dep(task243);
  task243->add_dep(task70);
  residualq->add_task(task243);

  vector<IndexRange> I176_index = {virt_, closed_};
  auto I176 = make_shared<Tensor>(I176_index);
  vector<shared_ptr<Tensor>> tensor244 = {I172, t2, I176};
  auto task244 = make_shared<Task244>(tensor244, pindex);
  task240->add_dep(task244);
  task244->add_dep(task70);
  residualq->add_task(task244);

  vector<shared_ptr<Tensor>> tensor245 = {I176, f1_};
  auto task245 = make_shared<Task245>(tensor245, pindex);
  task244->add_dep(task245);
  task245->add_dep(task70);
  residualq->add_task(task245);

  vector<IndexRange> I180_index = {virt_, closed_, virt_, closed_};
  auto I180 = make_shared<Tensor>(I180_index);
  vector<shared_ptr<Tensor>> tensor246 = {r, I180};
  auto task246 = make_shared<Task246>(tensor246, pindex);
  task246->add_dep(task70);
  residualq->add_task(task246);

  vector<IndexRange> I181_index = {virt_, active_};
  auto I181 = make_shared<Tensor>(I181_index);
  vector<shared_ptr<Tensor>> tensor247 = {I180, t2, I181};
  auto task247 = make_shared<Task247>(tensor247, pindex);
  task246->add_dep(task247);
  task247->add_dep(task70);
  residualq->add_task(task247);

  vector<IndexRange> I182_index = {active_, virt_};
  auto I182 = make_shared<Tensor>(I182_index);
  vector<shared_ptr<Tensor>> tensor248 = {I181, Gamma16_(), I182};
  auto task248 = make_shared<Task248>(tensor248, pindex);
  task247->add_dep(task248);
  task248->add_dep(task70);
  residualq->add_task(task248);

  vector<shared_ptr<Tensor>> tensor249 = {I182, f1_};
  auto task249 = make_shared<Task249>(tensor249, pindex);
  task248->add_dep(task249);
  task249->add_dep(task70);
  residualq->add_task(task249);

  vector<IndexRange> I184_index = {virt_, active_};
  auto I184 = make_shared<Tensor>(I184_index);
  vector<shared_ptr<Tensor>> tensor250 = {I180, t2, I184};
  auto task250 = make_shared<Task250>(tensor250, pindex);
  task246->add_dep(task250);
  task250->add_dep(task70);
  residualq->add_task(task250);

  vector<IndexRange> I185_index = {active_, virt_};
  auto I185 = make_shared<Tensor>(I185_index);
  vector<shared_ptr<Tensor>> tensor251 = {I184, Gamma16_(), I185};
  auto task251 = make_shared<Task251>(tensor251, pindex);
  task250->add_dep(task251);
  task251->add_dep(task70);
  residualq->add_task(task251);

  vector<shared_ptr<Tensor>> tensor252 = {I185, f1_};
  auto task252 = make_shared<Task252>(tensor252, pindex);
  task251->add_dep(task252);
  task252->add_dep(task70);
  residualq->add_task(task252);

  vector<IndexRange> I187_index = {virt_, closed_};
  auto I187 = make_shared<Tensor>(I187_index);
  vector<shared_ptr<Tensor>> tensor253 = {I180, f1_, I187};
  auto task253 = make_shared<Task253>(tensor253, pindex);
  task246->add_dep(task253);
  task253->add_dep(task70);
  residualq->add_task(task253);

  vector<IndexRange> I188_index = {active_, virt_, closed_, active_};
  auto I188 = make_shared<Tensor>(I188_index);
  vector<shared_ptr<Tensor>> tensor254 = {I187, Gamma38_(), I188};
  auto task254 = make_shared<Task254>(tensor254, pindex);
  task253->add_dep(task254);
  task254->add_dep(task70);
  residualq->add_task(task254);

  vector<shared_ptr<Tensor>> tensor255 = {I188, t2};
  auto task255 = make_shared<Task255>(tensor255, pindex);
  task254->add_dep(task255);
  task255->add_dep(task70);
  residualq->add_task(task255);

  vector<IndexRange> I190_index = {virt_, closed_};
  auto I190 = make_shared<Tensor>(I190_index);
  vector<shared_ptr<Tensor>> tensor256 = {I180, f1_, I190};
  auto task256 = make_shared<Task256>(tensor256, pindex);
  task246->add_dep(task256);
  task256->add_dep(task70);
  residualq->add_task(task256);

  vector<IndexRange> I191_index = {active_, virt_, closed_, active_};
  auto I191 = make_shared<Tensor>(I191_index);
  vector<shared_ptr<Tensor>> tensor257 = {I190, Gamma38_(), I191};
  auto task257 = make_shared<Task257>(tensor257, pindex);
  task256->add_dep(task257);
  task257->add_dep(task70);
  residualq->add_task(task257);

  vector<shared_ptr<Tensor>> tensor258 = {I191, t2};
  auto task258 = make_shared<Task258>(tensor258, pindex);
  task257->add_dep(task258);
  task258->add_dep(task70);
  residualq->add_task(task258);

  vector<IndexRange> I203_index = {closed_, closed_};
  auto I203 = make_shared<Tensor>(I203_index);
  vector<shared_ptr<Tensor>> tensor259 = {I180, t2, I203};
  auto task259 = make_shared<Task259>(tensor259, pindex);
  task246->add_dep(task259);
  task259->add_dep(task70);
  residualq->add_task(task259);

  vector<shared_ptr<Tensor>> tensor260 = {I203, f1_};
  auto task260 = make_shared<Task260>(tensor260, pindex);
  task259->add_dep(task260);
  task260->add_dep(task70);
  residualq->add_task(task260);

  vector<IndexRange> I205_index = {closed_, closed_};
  auto I205 = make_shared<Tensor>(I205_index);
  vector<shared_ptr<Tensor>> tensor261 = {I180, t2, I205};
  auto task261 = make_shared<Task261>(tensor261, pindex);
  task246->add_dep(task261);
  task261->add_dep(task70);
  residualq->add_task(task261);

  vector<shared_ptr<Tensor>> tensor262 = {I205, f1_};
  auto task262 = make_shared<Task262>(tensor262, pindex);
  task261->add_dep(task262);
  task262->add_dep(task70);
  residualq->add_task(task262);

  vector<IndexRange> I207_index = {virt_, virt_};
  auto I207 = make_shared<Tensor>(I207_index);
  vector<shared_ptr<Tensor>> tensor263 = {I180, t2, I207};
  auto task263 = make_shared<Task263>(tensor263, pindex);
  task246->add_dep(task263);
  task263->add_dep(task70);
  residualq->add_task(task263);

  vector<shared_ptr<Tensor>> tensor264 = {I207, f1_};
  auto task264 = make_shared<Task264>(tensor264, pindex);
  task263->add_dep(task264);
  task264->add_dep(task70);
  residualq->add_task(task264);

  vector<IndexRange> I209_index = {virt_, virt_};
  auto I209 = make_shared<Tensor>(I209_index);
  vector<shared_ptr<Tensor>> tensor265 = {I180, t2, I209};
  auto task265 = make_shared<Task265>(tensor265, pindex);
  task246->add_dep(task265);
  task265->add_dep(task70);
  residualq->add_task(task265);

  vector<shared_ptr<Tensor>> tensor266 = {I209, f1_};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task265->add_dep(task266);
  task266->add_dep(task70);
  residualq->add_task(task266);

  vector<IndexRange> I211_index = {closed_, active_};
  auto I211 = make_shared<Tensor>(I211_index);
  vector<shared_ptr<Tensor>> tensor267 = {I180, t2, I211};
  auto task267 = make_shared<Task267>(tensor267, pindex);
  task246->add_dep(task267);
  task267->add_dep(task70);
  residualq->add_task(task267);

  vector<IndexRange> I212_index = {closed_, active_};
  auto I212 = make_shared<Tensor>(I212_index);
  vector<shared_ptr<Tensor>> tensor268 = {I211, Gamma38_(), I212};
  auto task268 = make_shared<Task268>(tensor268, pindex);
  task267->add_dep(task268);
  task268->add_dep(task70);
  residualq->add_task(task268);

  vector<shared_ptr<Tensor>> tensor269 = {I212, f1_};
  auto task269 = make_shared<Task269>(tensor269, pindex);
  task268->add_dep(task269);
  task269->add_dep(task70);
  residualq->add_task(task269);

  vector<IndexRange> I214_index = {closed_, active_};
  auto I214 = make_shared<Tensor>(I214_index);
  vector<shared_ptr<Tensor>> tensor270 = {I180, t2, I214};
  auto task270 = make_shared<Task270>(tensor270, pindex);
  task246->add_dep(task270);
  task270->add_dep(task70);
  residualq->add_task(task270);

  vector<IndexRange> I215_index = {closed_, active_};
  auto I215 = make_shared<Tensor>(I215_index);
  vector<shared_ptr<Tensor>> tensor271 = {I214, Gamma38_(), I215};
  auto task271 = make_shared<Task271>(tensor271, pindex);
  task270->add_dep(task271);
  task271->add_dep(task70);
  residualq->add_task(task271);

  vector<shared_ptr<Tensor>> tensor272 = {I215, f1_};
  auto task272 = make_shared<Task272>(tensor272, pindex);
  task271->add_dep(task272);
  task272->add_dep(task70);
  residualq->add_task(task272);

  vector<IndexRange> I198_index = {closed_, virt_, closed_, virt_};
  auto I198 = make_shared<Tensor>(I198_index);
  vector<shared_ptr<Tensor>> tensor273 = {r, I198};
  auto task273 = make_shared<Task273>(tensor273, pindex);
  task273->add_dep(task70);
  residualq->add_task(task273);

  vector<shared_ptr<Tensor>> tensor274 = {I198, t2, v2_};
  auto task274 = make_shared<Task274>(tensor274, pindex, this->e0_);
  task273->add_dep(task274);
  task274->add_dep(task70);
  residualq->add_task(task274);

  vector<IndexRange> I199_index = {closed_, virt_, closed_, virt_};
  auto I199 = make_shared<Tensor>(I199_index);
  vector<shared_ptr<Tensor>> tensor275 = {I198, Gamma69_(), I199};
  auto task275 = make_shared<Task275>(tensor275, pindex);
  task273->add_dep(task275);
  task275->add_dep(task70);
  residualq->add_task(task275);

  vector<shared_ptr<Tensor>> tensor276 = {I199, t2};
  auto task276 = make_shared<Task276>(tensor276, pindex);
  task275->add_dep(task276);
  task276->add_dep(task70);
  residualq->add_task(task276);

  vector<IndexRange> I216_index = {virt_, closed_, active_, virt_};
  auto I216 = make_shared<Tensor>(I216_index);
  vector<shared_ptr<Tensor>> tensor277 = {r, I216};
  auto task277 = make_shared<Task277>(tensor277, pindex);
  task277->add_dep(task70);
  residualq->add_task(task277);

  vector<IndexRange> I217_index = {virt_, closed_, active_, active_};
  auto I217 = make_shared<Tensor>(I217_index);
  vector<shared_ptr<Tensor>> tensor278 = {I216, f1_, I217};
  auto task278 = make_shared<Task278>(tensor278, pindex);
  task277->add_dep(task278);
  task278->add_dep(task70);
  residualq->add_task(task278);

  vector<IndexRange> I218_index = {active_, virt_, closed_, active_};
  auto I218 = make_shared<Tensor>(I218_index);
  vector<shared_ptr<Tensor>> tensor279 = {I217, Gamma35_(), I218};
  auto task279 = make_shared<Task279>(tensor279, pindex);
  task278->add_dep(task279);
  task279->add_dep(task70);
  residualq->add_task(task279);

  vector<shared_ptr<Tensor>> tensor280 = {I218, t2};
  auto task280 = make_shared<Task280>(tensor280, pindex);
  task279->add_dep(task280);
  task280->add_dep(task70);
  residualq->add_task(task280);

  vector<IndexRange> I220_index = {virt_, closed_, active_, active_};
  auto I220 = make_shared<Tensor>(I220_index);
  vector<shared_ptr<Tensor>> tensor281 = {I216, f1_, I220};
  auto task281 = make_shared<Task281>(tensor281, pindex);
  task277->add_dep(task281);
  task281->add_dep(task70);
  residualq->add_task(task281);

  vector<IndexRange> I221_index = {active_, virt_, closed_, active_};
  auto I221 = make_shared<Tensor>(I221_index);
  vector<shared_ptr<Tensor>> tensor282 = {I220, Gamma32_(), I221};
  auto task282 = make_shared<Task282>(tensor282, pindex);
  task281->add_dep(task282);
  task282->add_dep(task70);
  residualq->add_task(task282);

  vector<shared_ptr<Tensor>> tensor283 = {I221, t2};
  auto task283 = make_shared<Task283>(tensor283, pindex);
  task282->add_dep(task283);
  task283->add_dep(task70);
  residualq->add_task(task283);

  vector<IndexRange> I227_index = {closed_, virt_, active_, active_};
  auto I227 = make_shared<Tensor>(I227_index);
  vector<shared_ptr<Tensor>> tensor284 = {I220, Gamma35_(), I227};
  auto task284 = make_shared<Task284>(tensor284, pindex);
  task281->add_dep(task284);
  task284->add_dep(task70);
  residualq->add_task(task284);

  vector<shared_ptr<Tensor>> tensor285 = {I227, t2};
  auto task285 = make_shared<Task285>(tensor285, pindex);
  task284->add_dep(task285);
  task285->add_dep(task70);
  residualq->add_task(task285);

  vector<IndexRange> I229_index = {virt_, active_};
  auto I229 = make_shared<Tensor>(I229_index);
  vector<shared_ptr<Tensor>> tensor286 = {I216, f1_, I229};
  auto task286 = make_shared<Task286>(tensor286, pindex);
  task277->add_dep(task286);
  task286->add_dep(task70);
  residualq->add_task(task286);

  vector<IndexRange> I230_index = {active_, virt_, active_, active_};
  auto I230 = make_shared<Tensor>(I230_index);
  vector<shared_ptr<Tensor>> tensor287 = {I229, Gamma60_(), I230};
  auto task287 = make_shared<Task287>(tensor287, pindex);
  task286->add_dep(task287);
  task287->add_dep(task70);
  residualq->add_task(task287);

  vector<shared_ptr<Tensor>> tensor288 = {I230, t2};
  auto task288 = make_shared<Task288>(tensor288, pindex);
  task287->add_dep(task288);
  task288->add_dep(task70);
  residualq->add_task(task288);

  vector<IndexRange> I232_index = {virt_, active_};
  auto I232 = make_shared<Tensor>(I232_index);
  vector<shared_ptr<Tensor>> tensor289 = {I216, f1_, I232};
  auto task289 = make_shared<Task289>(tensor289, pindex);
  task277->add_dep(task289);
  task289->add_dep(task70);
  residualq->add_task(task289);

  vector<IndexRange> I233_index = {active_, virt_, active_, active_};
  auto I233 = make_shared<Tensor>(I233_index);
  vector<shared_ptr<Tensor>> tensor290 = {I232, Gamma60_(), I233};
  auto task290 = make_shared<Task290>(tensor290, pindex);
  task289->add_dep(task290);
  task290->add_dep(task70);
  residualq->add_task(task290);

  vector<shared_ptr<Tensor>> tensor291 = {I233, t2};
  auto task291 = make_shared<Task291>(tensor291, pindex);
  task290->add_dep(task291);
  task291->add_dep(task70);
  residualq->add_task(task291);

  vector<IndexRange> I235_index = {closed_, active_};
  auto I235 = make_shared<Tensor>(I235_index);
  vector<shared_ptr<Tensor>> tensor292 = {I216, t2, I235};
  auto task292 = make_shared<Task292>(tensor292, pindex);
  task277->add_dep(task292);
  task292->add_dep(task70);
  residualq->add_task(task292);

  vector<IndexRange> I236_index = {active_, closed_};
  auto I236 = make_shared<Tensor>(I236_index);
  vector<shared_ptr<Tensor>> tensor293 = {I235, Gamma38_(), I236};
  auto task293 = make_shared<Task293>(tensor293, pindex);
  task292->add_dep(task293);
  task293->add_dep(task70);
  residualq->add_task(task293);

  vector<shared_ptr<Tensor>> tensor294 = {I236, f1_};
  auto task294 = make_shared<Task294>(tensor294, pindex);
  task293->add_dep(task294);
  task294->add_dep(task70);
  residualq->add_task(task294);

  vector<IndexRange> I238_index = {closed_, active_};
  auto I238 = make_shared<Tensor>(I238_index);
  vector<shared_ptr<Tensor>> tensor295 = {I216, t2, I238};
  auto task295 = make_shared<Task295>(tensor295, pindex);
  task277->add_dep(task295);
  task295->add_dep(task70);
  residualq->add_task(task295);

  vector<IndexRange> I239_index = {active_, closed_};
  auto I239 = make_shared<Tensor>(I239_index);
  vector<shared_ptr<Tensor>> tensor296 = {I238, Gamma38_(), I239};
  auto task296 = make_shared<Task296>(tensor296, pindex);
  task295->add_dep(task296);
  task296->add_dep(task70);
  residualq->add_task(task296);

  vector<shared_ptr<Tensor>> tensor297 = {I239, f1_};
  auto task297 = make_shared<Task297>(tensor297, pindex);
  task296->add_dep(task297);
  task297->add_dep(task70);
  residualq->add_task(task297);

  vector<IndexRange> I241_index = {active_, virt_, closed_, virt_};
  auto I241 = make_shared<Tensor>(I241_index);
  vector<shared_ptr<Tensor>> tensor298 = {I216, Gamma81_(), I241};
  auto task298 = make_shared<Task298>(tensor298, pindex);
  task277->add_dep(task298);
  task298->add_dep(task70);
  residualq->add_task(task298);

  vector<shared_ptr<Tensor>> tensor299 = {I241, t2};
  auto task299 = make_shared<Task299>(tensor299, pindex);
  task298->add_dep(task299);
  task299->add_dep(task70);
  residualq->add_task(task299);

  vector<IndexRange> I245_index = {closed_, active_, virt_, virt_};
  auto I245 = make_shared<Tensor>(I245_index);
  vector<shared_ptr<Tensor>> tensor300 = {I216, Gamma38_(), I245};
  auto task300 = make_shared<Task300>(tensor300, pindex);
  task277->add_dep(task300);
  task300->add_dep(task70);
  residualq->add_task(task300);

  vector<shared_ptr<Tensor>> tensor301 = {I245, t2, v2_};
  auto task301 = make_shared<Task301>(tensor301, pindex, this->e0_);
  task300->add_dep(task301);
  task301->add_dep(task70);
  residualq->add_task(task301);

  vector<IndexRange> I246_index = {closed_, closed_};
  auto I246 = make_shared<Tensor>(I246_index);
  vector<shared_ptr<Tensor>> tensor302 = {I245, t2, I246};
  auto task302 = make_shared<Task302>(tensor302, pindex);
  task300->add_dep(task302);
  task302->add_dep(task70);
  residualq->add_task(task302);

  vector<shared_ptr<Tensor>> tensor303 = {I246, f1_};
  auto task303 = make_shared<Task303>(tensor303, pindex);
  task302->add_dep(task303);
  task303->add_dep(task70);
  residualq->add_task(task303);

  vector<IndexRange> I249_index = {closed_, closed_};
  auto I249 = make_shared<Tensor>(I249_index);
  vector<shared_ptr<Tensor>> tensor304 = {I245, t2, I249};
  auto task304 = make_shared<Task304>(tensor304, pindex);
  task300->add_dep(task304);
  task304->add_dep(task70);
  residualq->add_task(task304);

  vector<shared_ptr<Tensor>> tensor305 = {I249, f1_};
  auto task305 = make_shared<Task305>(tensor305, pindex);
  task304->add_dep(task305);
  task305->add_dep(task70);
  residualq->add_task(task305);

  vector<IndexRange> I252_index = {virt_, virt_};
  auto I252 = make_shared<Tensor>(I252_index);
  vector<shared_ptr<Tensor>> tensor306 = {I245, t2, I252};
  auto task306 = make_shared<Task306>(tensor306, pindex);
  task300->add_dep(task306);
  task306->add_dep(task70);
  residualq->add_task(task306);

  vector<shared_ptr<Tensor>> tensor307 = {I252, f1_};
  auto task307 = make_shared<Task307>(tensor307, pindex);
  task306->add_dep(task307);
  task307->add_dep(task70);
  residualq->add_task(task307);

  vector<IndexRange> I255_index = {virt_, virt_};
  auto I255 = make_shared<Tensor>(I255_index);
  vector<shared_ptr<Tensor>> tensor308 = {I245, t2, I255};
  auto task308 = make_shared<Task308>(tensor308, pindex);
  task300->add_dep(task308);
  task308->add_dep(task70);
  residualq->add_task(task308);

  vector<shared_ptr<Tensor>> tensor309 = {I255, f1_};
  auto task309 = make_shared<Task309>(tensor309, pindex);
  task308->add_dep(task309);
  task309->add_dep(task70);
  residualq->add_task(task309);

  vector<IndexRange> I258_index = {virt_, virt_};
  auto I258 = make_shared<Tensor>(I258_index);
  vector<shared_ptr<Tensor>> tensor310 = {I245, t2, I258};
  auto task310 = make_shared<Task310>(tensor310, pindex);
  task300->add_dep(task310);
  task310->add_dep(task70);
  residualq->add_task(task310);

  vector<shared_ptr<Tensor>> tensor311 = {I258, f1_};
  auto task311 = make_shared<Task311>(tensor311, pindex);
  task310->add_dep(task311);
  task311->add_dep(task70);
  residualq->add_task(task311);

  vector<IndexRange> I261_index = {virt_, virt_};
  auto I261 = make_shared<Tensor>(I261_index);
  vector<shared_ptr<Tensor>> tensor312 = {I245, t2, I261};
  auto task312 = make_shared<Task312>(tensor312, pindex);
  task300->add_dep(task312);
  task312->add_dep(task70);
  residualq->add_task(task312);

  vector<shared_ptr<Tensor>> tensor313 = {I261, f1_};
  auto task313 = make_shared<Task313>(tensor313, pindex);
  task312->add_dep(task313);
  task313->add_dep(task70);
  residualq->add_task(task313);

  vector<IndexRange> I263_index = {virt_, virt_, active_, active_};
  auto I263 = make_shared<Tensor>(I263_index);
  vector<shared_ptr<Tensor>> tensor314 = {I216, f1_, I263};
  auto task314 = make_shared<Task314>(tensor314, pindex);
  task277->add_dep(task314);
  task314->add_dep(task70);
  residualq->add_task(task314);

  vector<IndexRange> I264_index = {active_, virt_, active_, virt_};
  auto I264 = make_shared<Tensor>(I264_index);
  vector<shared_ptr<Tensor>> tensor315 = {I263, Gamma60_(), I264};
  auto task315 = make_shared<Task315>(tensor315, pindex);
  task314->add_dep(task315);
  task315->add_dep(task70);
  residualq->add_task(task315);

  vector<shared_ptr<Tensor>> tensor316 = {I264, t2};
  auto task316 = make_shared<Task316>(tensor316, pindex);
  task315->add_dep(task316);
  task316->add_dep(task70);
  residualq->add_task(task316);

  vector<IndexRange> I265_index = {virt_, active_, active_, virt_};
  auto I265 = make_shared<Tensor>(I265_index);
  vector<shared_ptr<Tensor>> tensor317 = {r, I265};
  auto task317 = make_shared<Task317>(tensor317, pindex);
  task317->add_dep(task70);
  residualq->add_task(task317);

  vector<IndexRange> I266_index = {virt_, active_, active_, active_};
  auto I266 = make_shared<Tensor>(I266_index);
  vector<shared_ptr<Tensor>> tensor318 = {I265, f1_, I266};
  auto task318 = make_shared<Task318>(tensor318, pindex);
  task317->add_dep(task318);
  task318->add_dep(task70);
  residualq->add_task(task318);

  vector<IndexRange> I267_index = {active_, virt_, active_, active_};
  auto I267 = make_shared<Tensor>(I267_index);
  vector<shared_ptr<Tensor>> tensor319 = {I266, Gamma59_(), I267};
  auto task319 = make_shared<Task319>(tensor319, pindex);
  task318->add_dep(task319);
  task319->add_dep(task70);
  residualq->add_task(task319);

  vector<shared_ptr<Tensor>> tensor320 = {I267, t2};
  auto task320 = make_shared<Task320>(tensor320, pindex);
  task319->add_dep(task320);
  task320->add_dep(task70);
  residualq->add_task(task320);

  vector<IndexRange> I269_index = {active_, active_, virt_, virt_};
  auto I269 = make_shared<Tensor>(I269_index);
  vector<shared_ptr<Tensor>> tensor321 = {I265, Gamma60_(), I269};
  auto task321 = make_shared<Task321>(tensor321, pindex);
  task317->add_dep(task321);
  task321->add_dep(task70);
  residualq->add_task(task321);

  vector<IndexRange> I270_index = {active_, closed_};
  auto I270 = make_shared<Tensor>(I270_index);
  vector<shared_ptr<Tensor>> tensor322 = {I269, t2, I270};
  auto task322 = make_shared<Task322>(tensor322, pindex);
  task321->add_dep(task322);
  task322->add_dep(task70);
  residualq->add_task(task322);

  vector<shared_ptr<Tensor>> tensor323 = {I270, f1_};
  auto task323 = make_shared<Task323>(tensor323, pindex);
  task322->add_dep(task323);
  task323->add_dep(task70);
  residualq->add_task(task323);

  vector<IndexRange> I275_index = {virt_, virt_};
  auto I275 = make_shared<Tensor>(I275_index);
  vector<shared_ptr<Tensor>> tensor324 = {I269, t2, I275};
  auto task324 = make_shared<Task324>(tensor324, pindex);
  task321->add_dep(task324);
  task324->add_dep(task70);
  residualq->add_task(task324);

  vector<shared_ptr<Tensor>> tensor325 = {I275, f1_};
  auto task325 = make_shared<Task325>(tensor325, pindex);
  task324->add_dep(task325);
  task325->add_dep(task70);
  residualq->add_task(task325);

  vector<IndexRange> I271_index = {virt_, virt_, active_, active_};
  auto I271 = make_shared<Tensor>(I271_index);
  vector<shared_ptr<Tensor>> tensor326 = {r, I271};
  auto task326 = make_shared<Task326>(tensor326, pindex);
  task326->add_dep(task70);
  residualq->add_task(task326);

  vector<IndexRange> I272_index = {active_, virt_, active_, virt_};
  auto I272 = make_shared<Tensor>(I272_index);
  vector<shared_ptr<Tensor>> tensor327 = {I271, Gamma92_(), I272};
  auto task327 = make_shared<Task327>(tensor327, pindex);
  task326->add_dep(task327);
  task327->add_dep(task70);
  residualq->add_task(task327);

  vector<shared_ptr<Tensor>> tensor328 = {I272, t2};
  auto task328 = make_shared<Task328>(tensor328, pindex);
  task327->add_dep(task328);
  task328->add_dep(task70);
  residualq->add_task(task328);

  vector<IndexRange> I301_index = {active_, virt_, active_, virt_};
  auto I301 = make_shared<Tensor>(I301_index);
  vector<shared_ptr<Tensor>> tensor329 = {I271, Gamma60_(), I301};
  auto task329 = make_shared<Task329>(tensor329, pindex);
  task326->add_dep(task329);
  task329->add_dep(task70);
  residualq->add_task(task329);

  vector<shared_ptr<Tensor>> tensor330 = {I301, t2, v2_};
  auto task330 = make_shared<Task330>(tensor330, pindex, this->e0_);
  task329->add_dep(task330);
  task330->add_dep(task70);
  residualq->add_task(task330);

  return residualq;
}


