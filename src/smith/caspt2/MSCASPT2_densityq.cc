//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_densityqq.cc
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


#include <src/smith/caspt2/MSCASPT2.h>
#include <src/smith/caspt2/MSCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_densityq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto densityq = make_shared<Queue>();
  auto tensor60 = vector<shared_ptr<Tensor>>{den2};
  auto task60 = make_shared<Task60>(tensor60, reset);
  densityq->add_task(task60);

  vector<IndexRange> I0_index = {active_, active_};
  auto I0 = make_shared<Tensor>(I0_index);
  auto tensor61 = vector<shared_ptr<Tensor>>{den2, I0};
  auto task61 = make_shared<Task61>(tensor61, pindex);
  task61->add_dep(task60);
  densityq->add_task(task61);

  vector<IndexRange> I1_index = {active_, active_, active_, active_};
  auto I1 = make_shared<Tensor>(I1_index);
  auto tensor62 = vector<shared_ptr<Tensor>>{I0, Gamma0_(), I1};
  auto task62 = make_shared<Task62>(tensor62, pindex);
  task61->add_dep(task62);
  task62->add_dep(task60);
  densityq->add_task(task62);

  auto tensor63 = vector<shared_ptr<Tensor>>{I1, t2, l2};
  auto task63 = make_shared<Task63>(tensor63, pindex);
  task62->add_dep(task63);
  task63->add_dep(task60);
  densityq->add_task(task63);

  vector<IndexRange> I94_index = {active_, active_, active_, active_};
  auto I94 = make_shared<Tensor>(I94_index);
  auto tensor64 = vector<shared_ptr<Tensor>>{I0, Gamma31_(), I94};
  auto task64 = make_shared<Task64>(tensor64, pindex);
  task61->add_dep(task64);
  task64->add_dep(task60);
  densityq->add_task(task64);

  auto tensor65 = vector<shared_ptr<Tensor>>{I94, t2, l2};
  auto task65 = make_shared<Task65>(tensor65, pindex);
  task64->add_dep(task65);
  task65->add_dep(task60);
  densityq->add_task(task65);

  vector<IndexRange> I103_index = {active_, active_, active_, active_};
  auto I103 = make_shared<Tensor>(I103_index);
  auto tensor66 = vector<shared_ptr<Tensor>>{I0, Gamma34_(), I103};
  auto task66 = make_shared<Task66>(tensor66, pindex);
  task61->add_dep(task66);
  task66->add_dep(task60);
  densityq->add_task(task66);

  auto tensor67 = vector<shared_ptr<Tensor>>{I103, t2, l2};
  auto task67 = make_shared<Task67>(tensor67, pindex);
  task66->add_dep(task67);
  task67->add_dep(task60);
  densityq->add_task(task67);

  vector<IndexRange> I137_index = {active_, virt_, closed_, active_};
  auto I137 = make_shared<Tensor>(I137_index);
  auto tensor68 = vector<shared_ptr<Tensor>>{I103, t2, I137};
  auto task68 = make_shared<Task68>(tensor68, pindex);
  task66->add_dep(task68);
  task68->add_dep(task60);
  densityq->add_task(task68);

  auto tensor69 = vector<shared_ptr<Tensor>>{I137, l2};
  auto task69 = make_shared<Task69>(tensor69, pindex);
  task68->add_dep(task69);
  task69->add_dep(task60);
  densityq->add_task(task69);

  vector<IndexRange> I285_index = {active_, active_, active_, active_};
  auto I285 = make_shared<Tensor>(I285_index);
  auto tensor70 = vector<shared_ptr<Tensor>>{I0, Gamma92_(), I285};
  auto task70 = make_shared<Task70>(tensor70, pindex);
  task61->add_dep(task70);
  task70->add_dep(task60);
  densityq->add_task(task70);

  auto tensor71 = vector<shared_ptr<Tensor>>{I285, t2, l2};
  auto task71 = make_shared<Task71>(tensor71, pindex);
  task70->add_dep(task71);
  task71->add_dep(task60);
  densityq->add_task(task71);

  vector<IndexRange> I3_index = {closed_, closed_};
  auto I3 = make_shared<Tensor>(I3_index);
  auto tensor72 = vector<shared_ptr<Tensor>>{den2, I3};
  auto task72 = make_shared<Task72>(tensor72, pindex);
  task72->add_dep(task60);
  densityq->add_task(task72);

  vector<IndexRange> I4_index = {closed_, closed_, active_, active_};
  auto I4 = make_shared<Tensor>(I4_index);
  auto tensor73 = vector<shared_ptr<Tensor>>{I3, l2, I4};
  auto task73 = make_shared<Task73>(tensor73, pindex);
  task72->add_dep(task73);
  task73->add_dep(task60);
  densityq->add_task(task73);

  auto tensor74 = vector<shared_ptr<Tensor>>{I4, Gamma1_(), t2};
  auto task74 = make_shared<Task74>(tensor74, pindex);
  task73->add_dep(task74);
  task74->add_dep(task60);
  densityq->add_task(task74);

  vector<IndexRange> I97_index = {closed_, virt_, active_, active_};
  auto I97 = make_shared<Tensor>(I97_index);
  auto tensor75 = vector<shared_ptr<Tensor>>{I3, l2, I97};
  auto task75 = make_shared<Task75>(tensor75, pindex);
  task72->add_dep(task75);
  task75->add_dep(task60);
  densityq->add_task(task75);

  auto tensor76 = vector<shared_ptr<Tensor>>{I97, Gamma32_(), t2};
  auto task76 = make_shared<Task76>(tensor76, pindex);
  task75->add_dep(task76);
  task76->add_dep(task60);
  densityq->add_task(task76);

  vector<IndexRange> I106_index = {closed_, virt_, active_, active_};
  auto I106 = make_shared<Tensor>(I106_index);
  auto tensor77 = vector<shared_ptr<Tensor>>{I3, l2, I106};
  auto task77 = make_shared<Task77>(tensor77, pindex);
  task72->add_dep(task77);
  task77->add_dep(task60);
  densityq->add_task(task77);

  auto tensor78 = vector<shared_ptr<Tensor>>{I106, Gamma35_(), t2};
  auto task78 = make_shared<Task78>(tensor78, pindex);
  task77->add_dep(task78);
  task78->add_dep(task60);
  densityq->add_task(task78);

  vector<IndexRange> I6_index = {active_, closed_};
  auto I6 = make_shared<Tensor>(I6_index);
  auto tensor79 = vector<shared_ptr<Tensor>>{den2, I6};
  auto task79 = make_shared<Task79>(tensor79, pindex);
  task79->add_dep(task60);
  densityq->add_task(task79);

  vector<IndexRange> I7_index = {closed_, active_, active_, active_};
  auto I7 = make_shared<Tensor>(I7_index);
  auto tensor80 = vector<shared_ptr<Tensor>>{I6, t2, I7};
  auto task80 = make_shared<Task80>(tensor80, pindex);
  task79->add_dep(task80);
  task80->add_dep(task60);
  densityq->add_task(task80);

  auto tensor81 = vector<shared_ptr<Tensor>>{I7, Gamma2_(), l2};
  auto task81 = make_shared<Task81>(tensor81, pindex);
  task80->add_dep(task81);
  task81->add_dep(task60);
  densityq->add_task(task81);

  vector<IndexRange> I112_index = {virt_, active_, active_, active_};
  auto I112 = make_shared<Tensor>(I112_index);
  auto tensor82 = vector<shared_ptr<Tensor>>{I6, t2, I112};
  auto task82 = make_shared<Task82>(tensor82, pindex);
  task79->add_dep(task82);
  task82->add_dep(task60);
  densityq->add_task(task82);

  auto tensor83 = vector<shared_ptr<Tensor>>{I112, Gamma37_(), l2};
  auto task83 = make_shared<Task83>(tensor83, pindex);
  task82->add_dep(task83);
  task83->add_dep(task60);
  densityq->add_task(task83);

  vector<IndexRange> I9_index = {active_, virt_};
  auto I9 = make_shared<Tensor>(I9_index);
  auto tensor84 = vector<shared_ptr<Tensor>>{den2, I9};
  auto task84 = make_shared<Task84>(tensor84, pindex);
  task84->add_dep(task60);
  densityq->add_task(task84);

  vector<IndexRange> I10_index = {closed_, closed_, active_, active_};
  auto I10 = make_shared<Tensor>(I10_index);
  auto tensor85 = vector<shared_ptr<Tensor>>{I9, l2, I10};
  auto task85 = make_shared<Task85>(tensor85, pindex);
  task84->add_dep(task85);
  task85->add_dep(task60);
  densityq->add_task(task85);

  auto tensor86 = vector<shared_ptr<Tensor>>{I10, Gamma3_(), t2};
  auto task86 = make_shared<Task86>(tensor86, pindex);
  task85->add_dep(task86);
  task86->add_dep(task60);
  densityq->add_task(task86);

  vector<IndexRange> I121_index = {closed_, virt_, active_, active_};
  auto I121 = make_shared<Tensor>(I121_index);
  auto tensor87 = vector<shared_ptr<Tensor>>{I9, l2, I121};
  auto task87 = make_shared<Task87>(tensor87, pindex);
  task84->add_dep(task87);
  task87->add_dep(task60);
  densityq->add_task(task87);

  auto tensor88 = vector<shared_ptr<Tensor>>{I121, Gamma35_(), t2};
  auto task88 = make_shared<Task88>(tensor88, pindex);
  task87->add_dep(task88);
  task88->add_dep(task60);
  densityq->add_task(task88);

  vector<IndexRange> I124_index = {closed_, virt_, active_, active_};
  auto I124 = make_shared<Tensor>(I124_index);
  auto tensor89 = vector<shared_ptr<Tensor>>{I9, l2, I124};
  auto task89 = make_shared<Task89>(tensor89, pindex);
  task84->add_dep(task89);
  task89->add_dep(task60);
  densityq->add_task(task89);

  auto tensor90 = vector<shared_ptr<Tensor>>{I124, Gamma32_(), t2};
  auto task90 = make_shared<Task90>(tensor90, pindex);
  task89->add_dep(task90);
  task90->add_dep(task60);
  densityq->add_task(task90);

  vector<IndexRange> I163_index = {virt_, closed_, active_, active_};
  auto I163 = make_shared<Tensor>(I163_index);
  auto tensor91 = vector<shared_ptr<Tensor>>{I9, l2, I163};
  auto task91 = make_shared<Task91>(tensor91, pindex);
  task84->add_dep(task91);
  task91->add_dep(task60);
  densityq->add_task(task91);

  auto tensor92 = vector<shared_ptr<Tensor>>{I163, Gamma35_(), t2};
  auto task92 = make_shared<Task92>(tensor92, pindex);
  task91->add_dep(task92);
  task92->add_dep(task60);
  densityq->add_task(task92);

  vector<IndexRange> I166_index = {virt_, closed_, active_, active_};
  auto I166 = make_shared<Tensor>(I166_index);
  auto tensor93 = vector<shared_ptr<Tensor>>{I9, l2, I166};
  auto task93 = make_shared<Task93>(tensor93, pindex);
  task84->add_dep(task93);
  task93->add_dep(task60);
  densityq->add_task(task93);

  auto tensor94 = vector<shared_ptr<Tensor>>{I166, Gamma35_(), t2};
  auto task94 = make_shared<Task94>(tensor94, pindex);
  task93->add_dep(task94);
  task94->add_dep(task60);
  densityq->add_task(task94);

  vector<IndexRange> I12_index = {active_, closed_};
  auto I12 = make_shared<Tensor>(I12_index);
  auto tensor95 = vector<shared_ptr<Tensor>>{den2, I12};
  auto task95 = make_shared<Task95>(tensor95, pindex);
  task95->add_dep(task60);
  densityq->add_task(task95);

  vector<IndexRange> I13_index = {closed_, active_, active_, active_};
  auto I13 = make_shared<Tensor>(I13_index);
  auto tensor96 = vector<shared_ptr<Tensor>>{I12, l2, I13};
  auto task96 = make_shared<Task96>(tensor96, pindex);
  task95->add_dep(task96);
  task96->add_dep(task60);
  densityq->add_task(task96);

  auto tensor97 = vector<shared_ptr<Tensor>>{I13, Gamma4_(), t2};
  auto task97 = make_shared<Task97>(tensor97, pindex);
  task96->add_dep(task97);
  task97->add_dep(task60);
  densityq->add_task(task97);

  vector<IndexRange> I169_index = {virt_, active_, active_, active_};
  auto I169 = make_shared<Tensor>(I169_index);
  auto tensor98 = vector<shared_ptr<Tensor>>{I12, l2, I169};
  auto task98 = make_shared<Task98>(tensor98, pindex);
  task95->add_dep(task98);
  task98->add_dep(task60);
  densityq->add_task(task98);

  auto tensor99 = vector<shared_ptr<Tensor>>{I169, Gamma56_(), t2};
  auto task99 = make_shared<Task99>(tensor99, pindex);
  task98->add_dep(task99);
  task99->add_dep(task60);
  densityq->add_task(task99);

  vector<IndexRange> I172_index = {virt_, active_, active_, active_};
  auto I172 = make_shared<Tensor>(I172_index);
  auto tensor100 = vector<shared_ptr<Tensor>>{I12, l2, I172};
  auto task100 = make_shared<Task100>(tensor100, pindex);
  task95->add_dep(task100);
  task100->add_dep(task60);
  densityq->add_task(task100);

  auto tensor101 = vector<shared_ptr<Tensor>>{I172, Gamma57_(), t2};
  auto task101 = make_shared<Task101>(tensor101, pindex);
  task100->add_dep(task101);
  task101->add_dep(task60);
  densityq->add_task(task101);

  vector<IndexRange> I15_index = {active_, active_};
  auto I15 = make_shared<Tensor>(I15_index);
  auto tensor102 = vector<shared_ptr<Tensor>>{den2, I15};
  auto task102 = make_shared<Task102>(tensor102, pindex);
  task102->add_dep(task60);
  densityq->add_task(task102);

  vector<IndexRange> I16_index = {active_, active_, active_, active_, active_, active_};
  auto I16 = make_shared<Tensor>(I16_index);
  auto tensor103 = vector<shared_ptr<Tensor>>{I15, Gamma5_(), I16};
  auto task103 = make_shared<Task103>(tensor103, pindex);
  task102->add_dep(task103);
  task103->add_dep(task60);
  densityq->add_task(task103);

  auto tensor104 = vector<shared_ptr<Tensor>>{I16, t2, l2};
  auto task104 = make_shared<Task104>(tensor104, pindex);
  task103->add_dep(task104);
  task104->add_dep(task60);
  densityq->add_task(task104);

  vector<IndexRange> I175_index = {active_, active_, active_, active_, active_, active_};
  auto I175 = make_shared<Tensor>(I175_index);
  auto tensor105 = vector<shared_ptr<Tensor>>{I15, Gamma58_(), I175};
  auto task105 = make_shared<Task105>(tensor105, pindex);
  task102->add_dep(task105);
  task105->add_dep(task60);
  densityq->add_task(task105);

  auto tensor106 = vector<shared_ptr<Tensor>>{I175, t2, l2};
  auto task106 = make_shared<Task106>(tensor106, pindex);
  task105->add_dep(task106);
  task106->add_dep(task60);
  densityq->add_task(task106);

  vector<IndexRange> I18_index = {closed_, closed_};
  auto I18 = make_shared<Tensor>(I18_index);
  auto tensor107 = vector<shared_ptr<Tensor>>{den2, I18};
  auto task107 = make_shared<Task107>(tensor107, pindex);
  task107->add_dep(task60);
  densityq->add_task(task107);

  vector<IndexRange> I19_index = {closed_, active_, active_, active_};
  auto I19 = make_shared<Tensor>(I19_index);
  auto tensor108 = vector<shared_ptr<Tensor>>{I18, l2, I19};
  auto task108 = make_shared<Task108>(tensor108, pindex);
  task107->add_dep(task108);
  task108->add_dep(task60);
  densityq->add_task(task108);

  auto tensor109 = vector<shared_ptr<Tensor>>{I19, Gamma6_(), t2};
  auto task109 = make_shared<Task109>(tensor109, pindex);
  task108->add_dep(task109);
  task109->add_dep(task60);
  densityq->add_task(task109);

  vector<IndexRange> I21_index = {closed_, virt_};
  auto I21 = make_shared<Tensor>(I21_index);
  auto tensor110 = vector<shared_ptr<Tensor>>{den2, I21};
  auto task110 = make_shared<Task110>(tensor110, pindex);
  task110->add_dep(task60);
  densityq->add_task(task110);

  vector<IndexRange> I22_index = {closed_, active_};
  auto I22 = make_shared<Tensor>(I22_index);
  auto tensor111 = vector<shared_ptr<Tensor>>{I21, l2, I22};
  auto task111 = make_shared<Task111>(tensor111, pindex);
  task110->add_dep(task111);
  task111->add_dep(task60);
  densityq->add_task(task111);

  auto tensor112 = vector<shared_ptr<Tensor>>{I22, Gamma7_(), t2};
  auto task112 = make_shared<Task112>(tensor112, pindex);
  task111->add_dep(task112);
  task112->add_dep(task60);
  densityq->add_task(task112);

  vector<IndexRange> I25_index = {closed_, active_};
  auto I25 = make_shared<Tensor>(I25_index);
  auto tensor113 = vector<shared_ptr<Tensor>>{I21, l2, I25};
  auto task113 = make_shared<Task113>(tensor113, pindex);
  task110->add_dep(task113);
  task113->add_dep(task60);
  densityq->add_task(task113);

  auto tensor114 = vector<shared_ptr<Tensor>>{I25, Gamma7_(), t2};
  auto task114 = make_shared<Task114>(tensor114, pindex);
  task113->add_dep(task114);
  task114->add_dep(task60);
  densityq->add_task(task114);

  vector<IndexRange> I181_index = {virt_, active_};
  auto I181 = make_shared<Tensor>(I181_index);
  auto tensor115 = vector<shared_ptr<Tensor>>{I21, l2, I181};
  auto task115 = make_shared<Task115>(tensor115, pindex);
  task110->add_dep(task115);
  task115->add_dep(task60);
  densityq->add_task(task115);

  auto tensor116 = vector<shared_ptr<Tensor>>{I181, Gamma60_(), t2};
  auto task116 = make_shared<Task116>(tensor116, pindex);
  task115->add_dep(task116);
  task116->add_dep(task60);
  densityq->add_task(task116);

  vector<IndexRange> I184_index = {virt_, active_};
  auto I184 = make_shared<Tensor>(I184_index);
  auto tensor117 = vector<shared_ptr<Tensor>>{I21, l2, I184};
  auto task117 = make_shared<Task117>(tensor117, pindex);
  task110->add_dep(task117);
  task117->add_dep(task60);
  densityq->add_task(task117);

  auto tensor118 = vector<shared_ptr<Tensor>>{I184, Gamma60_(), t2};
  auto task118 = make_shared<Task118>(tensor118, pindex);
  task117->add_dep(task118);
  task118->add_dep(task60);
  densityq->add_task(task118);

  vector<IndexRange> I27_index = {active_, virt_};
  auto I27 = make_shared<Tensor>(I27_index);
  auto tensor119 = vector<shared_ptr<Tensor>>{den2, I27};
  auto task119 = make_shared<Task119>(tensor119, pindex);
  task119->add_dep(task60);
  densityq->add_task(task119);

  vector<IndexRange> I28_index = {closed_, active_, active_, active_};
  auto I28 = make_shared<Tensor>(I28_index);
  auto tensor120 = vector<shared_ptr<Tensor>>{I27, l2, I28};
  auto task120 = make_shared<Task120>(tensor120, pindex);
  task119->add_dep(task120);
  task120->add_dep(task60);
  densityq->add_task(task120);

  auto tensor121 = vector<shared_ptr<Tensor>>{I28, Gamma9_(), t2};
  auto task121 = make_shared<Task121>(tensor121, pindex);
  task120->add_dep(task121);
  task121->add_dep(task60);
  densityq->add_task(task121);

  vector<IndexRange> I31_index = {closed_, active_, active_, active_};
  auto I31 = make_shared<Tensor>(I31_index);
  auto tensor122 = vector<shared_ptr<Tensor>>{I27, l2, I31};
  auto task122 = make_shared<Task122>(tensor122, pindex);
  task119->add_dep(task122);
  task122->add_dep(task60);
  densityq->add_task(task122);

  auto tensor123 = vector<shared_ptr<Tensor>>{I31, Gamma6_(), t2};
  auto task123 = make_shared<Task123>(tensor123, pindex);
  task122->add_dep(task123);
  task123->add_dep(task60);
  densityq->add_task(task123);

  vector<IndexRange> I187_index = {virt_, active_, active_, active_};
  auto I187 = make_shared<Tensor>(I187_index);
  auto tensor124 = vector<shared_ptr<Tensor>>{I27, l2, I187};
  auto task124 = make_shared<Task124>(tensor124, pindex);
  task119->add_dep(task124);
  task124->add_dep(task60);
  densityq->add_task(task124);

  auto tensor125 = vector<shared_ptr<Tensor>>{I187, Gamma62_(), t2};
  auto task125 = make_shared<Task125>(tensor125, pindex);
  task124->add_dep(task125);
  task125->add_dep(task60);
  densityq->add_task(task125);

  vector<IndexRange> I33_index = {active_, virt_};
  auto I33 = make_shared<Tensor>(I33_index);
  auto tensor126 = vector<shared_ptr<Tensor>>{den2, I33};
  auto task126 = make_shared<Task126>(tensor126, pindex);
  task126->add_dep(task60);
  densityq->add_task(task126);

  vector<IndexRange> I34_index = {closed_, closed_, active_, active_};
  auto I34 = make_shared<Tensor>(I34_index);
  auto tensor127 = vector<shared_ptr<Tensor>>{I33, t2, I34};
  auto task127 = make_shared<Task127>(tensor127, pindex);
  task126->add_dep(task127);
  task127->add_dep(task60);
  densityq->add_task(task127);

  auto tensor128 = vector<shared_ptr<Tensor>>{I34, Gamma3_(), l2};
  auto task128 = make_shared<Task128>(tensor128, pindex);
  task127->add_dep(task128);
  task128->add_dep(task60);
  densityq->add_task(task128);

  vector<IndexRange> I36_index = {virt_, closed_};
  auto I36 = make_shared<Tensor>(I36_index);
  auto tensor129 = vector<shared_ptr<Tensor>>{den2, I36};
  auto task129 = make_shared<Task129>(tensor129, pindex);
  task129->add_dep(task60);
  densityq->add_task(task129);

  vector<IndexRange> I37_index = {closed_, active_};
  auto I37 = make_shared<Tensor>(I37_index);
  auto tensor130 = vector<shared_ptr<Tensor>>{I36, t2, I37};
  auto task130 = make_shared<Task130>(tensor130, pindex);
  task129->add_dep(task130);
  task130->add_dep(task60);
  densityq->add_task(task130);

  auto tensor131 = vector<shared_ptr<Tensor>>{I37, Gamma12_(), l2};
  auto task131 = make_shared<Task131>(tensor131, pindex);
  task130->add_dep(task131);
  task131->add_dep(task60);
  densityq->add_task(task131);

  vector<IndexRange> I39_index = {closed_, virt_};
  auto I39 = make_shared<Tensor>(I39_index);
  auto tensor132 = vector<shared_ptr<Tensor>>{den2, I39};
  auto task132 = make_shared<Task132>(tensor132, pindex);
  task132->add_dep(task60);
  densityq->add_task(task132);

  vector<IndexRange> I40_index = {closed_, active_};
  auto I40 = make_shared<Tensor>(I40_index);
  auto tensor133 = vector<shared_ptr<Tensor>>{I39, t2, I40};
  auto task133 = make_shared<Task133>(tensor133, pindex);
  task132->add_dep(task133);
  task133->add_dep(task60);
  densityq->add_task(task133);

  auto tensor134 = vector<shared_ptr<Tensor>>{I40, Gamma12_(), l2};
  auto task134 = make_shared<Task134>(tensor134, pindex);
  task133->add_dep(task134);
  task134->add_dep(task60);
  densityq->add_task(task134);

  vector<IndexRange> I196_index = {virt_, closed_};
  auto I196 = make_shared<Tensor>(I196_index);
  auto tensor135 = vector<shared_ptr<Tensor>>{I39, t2, I196};
  auto task135 = make_shared<Task135>(tensor135, pindex);
  task132->add_dep(task135);
  task135->add_dep(task60);
  densityq->add_task(task135);

  vector<IndexRange> I197_index = {active_, virt_, closed_, active_};
  auto I197 = make_shared<Tensor>(I197_index);
  auto tensor136 = vector<shared_ptr<Tensor>>{I196, Gamma65_(), I197};
  auto task136 = make_shared<Task136>(tensor136, pindex);
  task135->add_dep(task136);
  task136->add_dep(task60);
  densityq->add_task(task136);

  auto tensor137 = vector<shared_ptr<Tensor>>{I197, l2};
  auto task137 = make_shared<Task137>(tensor137, pindex);
  task136->add_dep(task137);
  task137->add_dep(task60);
  densityq->add_task(task137);

  vector<IndexRange> I42_index = {active_, active_};
  auto I42 = make_shared<Tensor>(I42_index);
  auto tensor138 = vector<shared_ptr<Tensor>>{den2, I42};
  auto task138 = make_shared<Task138>(tensor138, pindex);
  task138->add_dep(task60);
  densityq->add_task(task138);

  vector<IndexRange> I43_index = {active_, active_};
  auto I43 = make_shared<Tensor>(I43_index);
  auto tensor139 = vector<shared_ptr<Tensor>>{I42, Gamma14_(), I43};
  auto task139 = make_shared<Task139>(tensor139, pindex);
  task138->add_dep(task139);
  task139->add_dep(task60);
  densityq->add_task(task139);

  vector<IndexRange> I44_index = {closed_, virt_, closed_, active_};
  auto I44 = make_shared<Tensor>(I44_index);
  auto tensor140 = vector<shared_ptr<Tensor>>{I43, t2, I44};
  auto task140 = make_shared<Task140>(tensor140, pindex);
  task139->add_dep(task140);
  task140->add_dep(task60);
  densityq->add_task(task140);

  auto tensor141 = vector<shared_ptr<Tensor>>{I44, l2};
  auto task141 = make_shared<Task141>(tensor141, pindex);
  task140->add_dep(task141);
  task141->add_dep(task60);
  densityq->add_task(task141);

  vector<IndexRange> I252_index = {active_, active_};
  auto I252 = make_shared<Tensor>(I252_index);
  auto tensor142 = vector<shared_ptr<Tensor>>{I42, Gamma60_(), I252};
  auto task142 = make_shared<Task142>(tensor142, pindex);
  task138->add_dep(task142);
  task142->add_dep(task60);
  densityq->add_task(task142);

  vector<IndexRange> I253_index = {active_, virt_, closed_, virt_};
  auto I253 = make_shared<Tensor>(I253_index);
  auto tensor143 = vector<shared_ptr<Tensor>>{I252, t2, I253};
  auto task143 = make_shared<Task143>(tensor143, pindex);
  task142->add_dep(task143);
  task143->add_dep(task60);
  densityq->add_task(task143);

  auto tensor144 = vector<shared_ptr<Tensor>>{I253, l2};
  auto task144 = make_shared<Task144>(tensor144, pindex);
  task143->add_dep(task144);
  task144->add_dep(task60);
  densityq->add_task(task144);

  vector<IndexRange> I48_index = {closed_, closed_};
  auto I48 = make_shared<Tensor>(I48_index);
  auto tensor145 = vector<shared_ptr<Tensor>>{den2, I48};
  auto task145 = make_shared<Task145>(tensor145, pindex);
  task145->add_dep(task60);
  densityq->add_task(task145);

  vector<IndexRange> I49_index = {closed_, virt_, closed_, active_};
  auto I49 = make_shared<Tensor>(I49_index);
  auto tensor146 = vector<shared_ptr<Tensor>>{I48, l2, I49};
  auto task146 = make_shared<Task146>(tensor146, pindex);
  task145->add_dep(task146);
  task146->add_dep(task60);
  densityq->add_task(task146);

  auto tensor147 = vector<shared_ptr<Tensor>>{I49, Gamma16_(), t2};
  auto task147 = make_shared<Task147>(tensor147, pindex);
  task146->add_dep(task147);
  task147->add_dep(task60);
  densityq->add_task(task147);

  vector<IndexRange> I52_index = {closed_, virt_, closed_, active_};
  auto I52 = make_shared<Tensor>(I52_index);
  auto tensor148 = vector<shared_ptr<Tensor>>{I48, l2, I52};
  auto task148 = make_shared<Task148>(tensor148, pindex);
  task145->add_dep(task148);
  task148->add_dep(task60);
  densityq->add_task(task148);

  auto tensor149 = vector<shared_ptr<Tensor>>{I52, Gamma16_(), t2};
  auto task149 = make_shared<Task149>(tensor149, pindex);
  task148->add_dep(task149);
  task149->add_dep(task60);
  densityq->add_task(task149);

  vector<IndexRange> I54_index = {closed_, closed_};
  auto I54 = make_shared<Tensor>(I54_index);
  auto tensor150 = vector<shared_ptr<Tensor>>{den2, I54};
  auto task150 = make_shared<Task150>(tensor150, pindex);
  task150->add_dep(task60);
  densityq->add_task(task150);

  vector<IndexRange> I55_index = {closed_, virt_, closed_, active_};
  auto I55 = make_shared<Tensor>(I55_index);
  auto tensor151 = vector<shared_ptr<Tensor>>{I54, l2, I55};
  auto task151 = make_shared<Task151>(tensor151, pindex);
  task150->add_dep(task151);
  task151->add_dep(task60);
  densityq->add_task(task151);

  auto tensor152 = vector<shared_ptr<Tensor>>{I55, Gamma16_(), t2};
  auto task152 = make_shared<Task152>(tensor152, pindex);
  task151->add_dep(task152);
  task152->add_dep(task60);
  densityq->add_task(task152);

  vector<IndexRange> I61_index = {closed_, virt_, closed_, active_};
  auto I61 = make_shared<Tensor>(I61_index);
  auto tensor153 = vector<shared_ptr<Tensor>>{I54, l2, I61};
  auto task153 = make_shared<Task153>(tensor153, pindex);
  task150->add_dep(task153);
  task153->add_dep(task60);
  densityq->add_task(task153);

  auto tensor154 = vector<shared_ptr<Tensor>>{I61, Gamma16_(), t2};
  auto task154 = make_shared<Task154>(tensor154, pindex);
  task153->add_dep(task154);
  task154->add_dep(task60);
  densityq->add_task(task154);

  vector<IndexRange> I57_index = {virt_, virt_};
  auto I57 = make_shared<Tensor>(I57_index);
  auto tensor155 = vector<shared_ptr<Tensor>>{den2, I57};
  auto task155 = make_shared<Task155>(tensor155, pindex);
  task155->add_dep(task60);
  densityq->add_task(task155);

  vector<IndexRange> I58_index = {closed_, virt_, closed_, active_};
  auto I58 = make_shared<Tensor>(I58_index);
  auto tensor156 = vector<shared_ptr<Tensor>>{I57, l2, I58};
  auto task156 = make_shared<Task156>(tensor156, pindex);
  task155->add_dep(task156);
  task156->add_dep(task60);
  densityq->add_task(task156);

  auto tensor157 = vector<shared_ptr<Tensor>>{I58, Gamma16_(), t2};
  auto task157 = make_shared<Task157>(tensor157, pindex);
  task156->add_dep(task157);
  task157->add_dep(task60);
  densityq->add_task(task157);

  vector<IndexRange> I64_index = {closed_, virt_, closed_, active_};
  auto I64 = make_shared<Tensor>(I64_index);
  auto tensor158 = vector<shared_ptr<Tensor>>{I57, l2, I64};
  auto task158 = make_shared<Task158>(tensor158, pindex);
  task155->add_dep(task158);
  task158->add_dep(task60);
  densityq->add_task(task158);

  auto tensor159 = vector<shared_ptr<Tensor>>{I64, Gamma16_(), t2};
  auto task159 = make_shared<Task159>(tensor159, pindex);
  task158->add_dep(task159);
  task159->add_dep(task60);
  densityq->add_task(task159);

  vector<IndexRange> I66_index = {active_, closed_};
  auto I66 = make_shared<Tensor>(I66_index);
  auto tensor160 = vector<shared_ptr<Tensor>>{den2, I66};
  auto task160 = make_shared<Task160>(tensor160, pindex);
  task160->add_dep(task60);
  densityq->add_task(task160);

  vector<IndexRange> I67_index = {virt_, closed_, active_, active_};
  auto I67 = make_shared<Tensor>(I67_index);
  auto tensor161 = vector<shared_ptr<Tensor>>{I66, t2, I67};
  auto task161 = make_shared<Task161>(tensor161, pindex);
  task160->add_dep(task161);
  task161->add_dep(task60);
  densityq->add_task(task161);

  auto tensor162 = vector<shared_ptr<Tensor>>{I67, Gamma22_(), l2};
  auto task162 = make_shared<Task162>(tensor162, pindex);
  task161->add_dep(task162);
  task162->add_dep(task60);
  densityq->add_task(task162);

  auto tensor163 = vector<shared_ptr<Tensor>>{I67, Gamma12_(), l2};
  auto task163 = make_shared<Task163>(tensor163, pindex);
  task161->add_dep(task163);
  task163->add_dep(task60);
  densityq->add_task(task163);

  vector<IndexRange> I69_index = {active_, closed_};
  auto I69 = make_shared<Tensor>(I69_index);
  auto tensor164 = vector<shared_ptr<Tensor>>{den2, I69};
  auto task164 = make_shared<Task164>(tensor164, pindex);
  task164->add_dep(task60);
  densityq->add_task(task164);

  vector<IndexRange> I70_index = {virt_, closed_, active_, active_};
  auto I70 = make_shared<Tensor>(I70_index);
  auto tensor165 = vector<shared_ptr<Tensor>>{I69, t2, I70};
  auto task165 = make_shared<Task165>(tensor165, pindex);
  task164->add_dep(task165);
  task165->add_dep(task60);
  densityq->add_task(task165);

  vector<IndexRange> I71_index = {active_, virt_, closed_, active_};
  auto I71 = make_shared<Tensor>(I71_index);
  auto tensor166 = vector<shared_ptr<Tensor>>{I70, Gamma12_(), I71};
  auto task166 = make_shared<Task166>(tensor166, pindex);
  task165->add_dep(task166);
  task166->add_dep(task60);
  densityq->add_task(task166);

  auto tensor167 = vector<shared_ptr<Tensor>>{I71, l2};
  auto task167 = make_shared<Task167>(tensor167, pindex);
  task166->add_dep(task167);
  task167->add_dep(task60);
  densityq->add_task(task167);

  vector<IndexRange> I78_index = {virt_, active_};
  auto I78 = make_shared<Tensor>(I78_index);
  auto tensor168 = vector<shared_ptr<Tensor>>{den2, I78};
  auto task168 = make_shared<Task168>(tensor168, pindex);
  task168->add_dep(task60);
  densityq->add_task(task168);

  vector<IndexRange> I79_index = {virt_, active_};
  auto I79 = make_shared<Tensor>(I79_index);
  auto tensor169 = vector<shared_ptr<Tensor>>{I78, Gamma16_(), I79};
  auto task169 = make_shared<Task169>(tensor169, pindex);
  task168->add_dep(task169);
  task169->add_dep(task60);
  densityq->add_task(task169);

  vector<IndexRange> I80_index = {closed_, virt_, closed_, virt_};
  auto I80 = make_shared<Tensor>(I80_index);
  auto tensor170 = vector<shared_ptr<Tensor>>{I79, t2, I80};
  auto task170 = make_shared<Task170>(tensor170, pindex);
  task169->add_dep(task170);
  task170->add_dep(task60);
  densityq->add_task(task170);

  auto tensor171 = vector<shared_ptr<Tensor>>{I80, l2};
  auto task171 = make_shared<Task171>(tensor171, pindex);
  task170->add_dep(task171);
  task171->add_dep(task60);
  densityq->add_task(task171);

  vector<IndexRange> I84_index = {active_, virt_};
  auto I84 = make_shared<Tensor>(I84_index);
  auto tensor172 = vector<shared_ptr<Tensor>>{den2, I84};
  auto task172 = make_shared<Task172>(tensor172, pindex);
  task172->add_dep(task60);
  densityq->add_task(task172);

  vector<IndexRange> I85_index = {closed_, active_, active_, active_};
  auto I85 = make_shared<Tensor>(I85_index);
  auto tensor173 = vector<shared_ptr<Tensor>>{I84, t2, I85};
  auto task173 = make_shared<Task173>(tensor173, pindex);
  task172->add_dep(task173);
  task173->add_dep(task60);
  densityq->add_task(task173);

  auto tensor174 = vector<shared_ptr<Tensor>>{I85, Gamma28_(), l2};
  auto task174 = make_shared<Task174>(tensor174, pindex);
  task173->add_dep(task174);
  task174->add_dep(task60);
  densityq->add_task(task174);

  vector<IndexRange> I87_index = {active_, closed_};
  auto I87 = make_shared<Tensor>(I87_index);
  auto tensor175 = vector<shared_ptr<Tensor>>{den2, I87};
  auto task175 = make_shared<Task175>(tensor175, pindex);
  task175->add_dep(task60);
  densityq->add_task(task175);

  vector<IndexRange> I88_index = {closed_, virt_, active_, active_};
  auto I88 = make_shared<Tensor>(I88_index);
  auto tensor176 = vector<shared_ptr<Tensor>>{I87, l2, I88};
  auto task176 = make_shared<Task176>(tensor176, pindex);
  task175->add_dep(task176);
  task176->add_dep(task60);
  densityq->add_task(task176);

  auto tensor177 = vector<shared_ptr<Tensor>>{I88, Gamma29_(), t2};
  auto task177 = make_shared<Task177>(tensor177, pindex);
  task176->add_dep(task177);
  task177->add_dep(task60);
  densityq->add_task(task177);

  vector<IndexRange> I91_index = {closed_, virt_, active_, active_};
  auto I91 = make_shared<Tensor>(I91_index);
  auto tensor178 = vector<shared_ptr<Tensor>>{I87, l2, I91};
  auto task178 = make_shared<Task178>(tensor178, pindex);
  task175->add_dep(task178);
  task178->add_dep(task60);
  densityq->add_task(task178);

  auto tensor179 = vector<shared_ptr<Tensor>>{I91, Gamma7_(), t2};
  auto task179 = make_shared<Task179>(tensor179, pindex);
  task178->add_dep(task179);
  task179->add_dep(task60);
  densityq->add_task(task179);

  vector<IndexRange> I130_index = {virt_, closed_, active_, active_};
  auto I130 = make_shared<Tensor>(I130_index);
  auto tensor180 = vector<shared_ptr<Tensor>>{I87, l2, I130};
  auto task180 = make_shared<Task180>(tensor180, pindex);
  task175->add_dep(task180);
  task180->add_dep(task60);
  densityq->add_task(task180);

  auto tensor181 = vector<shared_ptr<Tensor>>{I130, Gamma7_(), t2};
  auto task181 = make_shared<Task181>(tensor181, pindex);
  task180->add_dep(task181);
  task181->add_dep(task60);
  densityq->add_task(task181);

  vector<IndexRange> I133_index = {virt_, closed_, active_, active_};
  auto I133 = make_shared<Tensor>(I133_index);
  auto tensor182 = vector<shared_ptr<Tensor>>{I87, l2, I133};
  auto task182 = make_shared<Task182>(tensor182, pindex);
  task175->add_dep(task182);
  task182->add_dep(task60);
  densityq->add_task(task182);

  auto tensor183 = vector<shared_ptr<Tensor>>{I133, Gamma7_(), t2};
  auto task183 = make_shared<Task183>(tensor183, pindex);
  task182->add_dep(task183);
  task183->add_dep(task60);
  densityq->add_task(task183);

  vector<IndexRange> I282_index = {virt_, virt_, active_, active_};
  auto I282 = make_shared<Tensor>(I282_index);
  auto tensor184 = vector<shared_ptr<Tensor>>{I87, l2, I282};
  auto task184 = make_shared<Task184>(tensor184, pindex);
  task175->add_dep(task184);
  task184->add_dep(task60);
  densityq->add_task(task184);

  auto tensor185 = vector<shared_ptr<Tensor>>{I282, Gamma60_(), t2};
  auto task185 = make_shared<Task185>(tensor185, pindex);
  task184->add_dep(task185);
  task185->add_dep(task60);
  densityq->add_task(task185);

  vector<IndexRange> I99_index = {virt_, virt_};
  auto I99 = make_shared<Tensor>(I99_index);
  auto tensor186 = vector<shared_ptr<Tensor>>{den2, I99};
  auto task186 = make_shared<Task186>(tensor186, pindex);
  task186->add_dep(task60);
  densityq->add_task(task186);

  vector<IndexRange> I100_index = {closed_, virt_, active_, active_};
  auto I100 = make_shared<Tensor>(I100_index);
  auto tensor187 = vector<shared_ptr<Tensor>>{I99, l2, I100};
  auto task187 = make_shared<Task187>(tensor187, pindex);
  task186->add_dep(task187);
  task187->add_dep(task60);
  densityq->add_task(task187);

  auto tensor188 = vector<shared_ptr<Tensor>>{I100, Gamma32_(), t2};
  auto task188 = make_shared<Task188>(tensor188, pindex);
  task187->add_dep(task188);
  task188->add_dep(task60);
  densityq->add_task(task188);

  vector<IndexRange> I109_index = {closed_, virt_, active_, active_};
  auto I109 = make_shared<Tensor>(I109_index);
  auto tensor189 = vector<shared_ptr<Tensor>>{I99, l2, I109};
  auto task189 = make_shared<Task189>(tensor189, pindex);
  task186->add_dep(task189);
  task189->add_dep(task60);
  densityq->add_task(task189);

  auto tensor190 = vector<shared_ptr<Tensor>>{I109, Gamma35_(), t2};
  auto task190 = make_shared<Task190>(tensor190, pindex);
  task189->add_dep(task190);
  task190->add_dep(task60);
  densityq->add_task(task190);

  vector<IndexRange> I114_index = {virt_, closed_};
  auto I114 = make_shared<Tensor>(I114_index);
  auto tensor191 = vector<shared_ptr<Tensor>>{den2, I114};
  auto task191 = make_shared<Task191>(tensor191, pindex);
  task191->add_dep(task60);
  densityq->add_task(task191);

  vector<IndexRange> I115_index = {closed_, virt_};
  auto I115 = make_shared<Tensor>(I115_index);
  auto tensor192 = vector<shared_ptr<Tensor>>{I114, l2, I115};
  auto task192 = make_shared<Task192>(tensor192, pindex);
  task191->add_dep(task192);
  task192->add_dep(task60);
  densityq->add_task(task192);

  auto tensor193 = vector<shared_ptr<Tensor>>{I115, Gamma65_(), t2};
  auto task193 = make_shared<Task193>(tensor193, pindex);
  task192->add_dep(task193);
  task193->add_dep(task60);
  densityq->add_task(task193);

  vector<IndexRange> I118_index = {closed_, virt_};
  auto I118 = make_shared<Tensor>(I118_index);
  auto tensor194 = vector<shared_ptr<Tensor>>{I114, l2, I118};
  auto task194 = make_shared<Task194>(tensor194, pindex);
  task191->add_dep(task194);
  task194->add_dep(task60);
  densityq->add_task(task194);

  auto tensor195 = vector<shared_ptr<Tensor>>{I118, Gamma65_(), t2};
  auto task195 = make_shared<Task195>(tensor195, pindex);
  task194->add_dep(task195);
  task195->add_dep(task60);
  densityq->add_task(task195);

  vector<IndexRange> I157_index = {virt_, closed_};
  auto I157 = make_shared<Tensor>(I157_index);
  auto tensor196 = vector<shared_ptr<Tensor>>{I114, l2, I157};
  auto task196 = make_shared<Task196>(tensor196, pindex);
  task191->add_dep(task196);
  task196->add_dep(task60);
  densityq->add_task(task196);

  auto tensor197 = vector<shared_ptr<Tensor>>{I157, Gamma65_(), t2};
  auto task197 = make_shared<Task197>(tensor197, pindex);
  task196->add_dep(task197);
  task197->add_dep(task60);
  densityq->add_task(task197);

  vector<IndexRange> I160_index = {virt_, closed_};
  auto I160 = make_shared<Tensor>(I160_index);
  auto tensor198 = vector<shared_ptr<Tensor>>{I114, l2, I160};
  auto task198 = make_shared<Task198>(tensor198, pindex);
  task191->add_dep(task198);
  task198->add_dep(task60);
  densityq->add_task(task198);

  auto tensor199 = vector<shared_ptr<Tensor>>{I160, Gamma65_(), t2};
  auto task199 = make_shared<Task199>(tensor199, pindex);
  task198->add_dep(task199);
  task199->add_dep(task60);
  densityq->add_task(task199);

  vector<IndexRange> I126_index = {active_, virt_};
  auto I126 = make_shared<Tensor>(I126_index);
  auto tensor200 = vector<shared_ptr<Tensor>>{den2, I126};
  auto task200 = make_shared<Task200>(tensor200, pindex);
  task200->add_dep(task60);
  densityq->add_task(task200);

  vector<IndexRange> I127_index = {closed_, active_, active_, active_};
  auto I127 = make_shared<Tensor>(I127_index);
  auto tensor201 = vector<shared_ptr<Tensor>>{I126, t2, I127};
  auto task201 = make_shared<Task201>(tensor201, pindex);
  task200->add_dep(task201);
  task201->add_dep(task60);
  densityq->add_task(task201);

  auto tensor202 = vector<shared_ptr<Tensor>>{I127, Gamma6_(), l2};
  auto task202 = make_shared<Task202>(tensor202, pindex);
  task201->add_dep(task202);
  task202->add_dep(task60);
  densityq->add_task(task202);

  vector<IndexRange> I279_index = {virt_, active_, active_, active_};
  auto I279 = make_shared<Tensor>(I279_index);
  auto tensor203 = vector<shared_ptr<Tensor>>{I126, t2, I279};
  auto task203 = make_shared<Task203>(tensor203, pindex);
  task200->add_dep(task203);
  task203->add_dep(task60);
  densityq->add_task(task203);

  auto tensor204 = vector<shared_ptr<Tensor>>{I279, Gamma62_(), l2};
  auto task204 = make_shared<Task204>(tensor204, pindex);
  task203->add_dep(task204);
  task204->add_dep(task60);
  densityq->add_task(task204);

  vector<IndexRange> I138_index = {closed_, closed_};
  auto I138 = make_shared<Tensor>(I138_index);
  auto tensor205 = vector<shared_ptr<Tensor>>{den2, I138};
  auto task205 = make_shared<Task205>(tensor205, pindex);
  task205->add_dep(task60);
  densityq->add_task(task205);

  vector<IndexRange> I139_index = {virt_, closed_, active_, active_};
  auto I139 = make_shared<Tensor>(I139_index);
  auto tensor206 = vector<shared_ptr<Tensor>>{I138, l2, I139};
  auto task206 = make_shared<Task206>(tensor206, pindex);
  task205->add_dep(task206);
  task206->add_dep(task60);
  densityq->add_task(task206);

  auto tensor207 = vector<shared_ptr<Tensor>>{I139, Gamma35_(), t2};
  auto task207 = make_shared<Task207>(tensor207, pindex);
  task206->add_dep(task207);
  task207->add_dep(task60);
  densityq->add_task(task207);

  vector<IndexRange> I148_index = {virt_, closed_, active_, active_};
  auto I148 = make_shared<Tensor>(I148_index);
  auto tensor208 = vector<shared_ptr<Tensor>>{I138, l2, I148};
  auto task208 = make_shared<Task208>(tensor208, pindex);
  task205->add_dep(task208);
  task208->add_dep(task60);
  densityq->add_task(task208);

  auto tensor209 = vector<shared_ptr<Tensor>>{I148, Gamma35_(), t2};
  auto task209 = make_shared<Task209>(tensor209, pindex);
  task208->add_dep(task209);
  task209->add_dep(task60);
  densityq->add_task(task209);

  vector<IndexRange> I141_index = {virt_, virt_};
  auto I141 = make_shared<Tensor>(I141_index);
  auto tensor210 = vector<shared_ptr<Tensor>>{den2, I141};
  auto task210 = make_shared<Task210>(tensor210, pindex);
  task210->add_dep(task60);
  densityq->add_task(task210);

  vector<IndexRange> I142_index = {virt_, closed_, active_, active_};
  auto I142 = make_shared<Tensor>(I142_index);
  auto tensor211 = vector<shared_ptr<Tensor>>{I141, l2, I142};
  auto task211 = make_shared<Task211>(tensor211, pindex);
  task210->add_dep(task211);
  task211->add_dep(task60);
  densityq->add_task(task211);

  auto tensor212 = vector<shared_ptr<Tensor>>{I142, Gamma35_(), t2};
  auto task212 = make_shared<Task212>(tensor212, pindex);
  task211->add_dep(task212);
  task212->add_dep(task60);
  densityq->add_task(task212);

  vector<IndexRange> I151_index = {virt_, closed_, active_, active_};
  auto I151 = make_shared<Tensor>(I151_index);
  auto tensor213 = vector<shared_ptr<Tensor>>{I141, l2, I151};
  auto task213 = make_shared<Task213>(tensor213, pindex);
  task210->add_dep(task213);
  task213->add_dep(task60);
  densityq->add_task(task213);

  auto tensor214 = vector<shared_ptr<Tensor>>{I151, Gamma35_(), t2};
  auto task214 = make_shared<Task214>(tensor214, pindex);
  task213->add_dep(task214);
  task214->add_dep(task60);
  densityq->add_task(task214);

  vector<IndexRange> I288_index = {virt_, virt_, active_, active_};
  auto I288 = make_shared<Tensor>(I288_index);
  auto tensor215 = vector<shared_ptr<Tensor>>{I141, l2, I288};
  auto task215 = make_shared<Task215>(tensor215, pindex);
  task210->add_dep(task215);
  task215->add_dep(task60);
  densityq->add_task(task215);

  auto tensor216 = vector<shared_ptr<Tensor>>{I288, Gamma60_(), t2};
  auto task216 = make_shared<Task216>(tensor216, pindex);
  task215->add_dep(task216);
  task216->add_dep(task60);
  densityq->add_task(task216);

  vector<IndexRange> I153_index = {active_, closed_};
  auto I153 = make_shared<Tensor>(I153_index);
  auto tensor217 = vector<shared_ptr<Tensor>>{den2, I153};
  auto task217 = make_shared<Task217>(tensor217, pindex);
  task217->add_dep(task60);
  densityq->add_task(task217);

  vector<IndexRange> I154_index = {virt_, active_, active_, active_};
  auto I154 = make_shared<Tensor>(I154_index);
  auto tensor218 = vector<shared_ptr<Tensor>>{I153, t2, I154};
  auto task218 = make_shared<Task218>(tensor218, pindex);
  task217->add_dep(task218);
  task218->add_dep(task60);
  densityq->add_task(task218);

  auto tensor219 = vector<shared_ptr<Tensor>>{I154, Gamma51_(), l2};
  auto task219 = make_shared<Task219>(tensor219, pindex);
  task218->add_dep(task219);
  task219->add_dep(task60);
  densityq->add_task(task219);

  vector<IndexRange> I177_index = {virt_, virt_};
  auto I177 = make_shared<Tensor>(I177_index);
  auto tensor220 = vector<shared_ptr<Tensor>>{den2, I177};
  auto task220 = make_shared<Task220>(tensor220, pindex);
  task220->add_dep(task60);
  densityq->add_task(task220);

  vector<IndexRange> I178_index = {virt_, active_, active_, active_};
  auto I178 = make_shared<Tensor>(I178_index);
  auto tensor221 = vector<shared_ptr<Tensor>>{I177, l2, I178};
  auto task221 = make_shared<Task221>(tensor221, pindex);
  task220->add_dep(task221);
  task221->add_dep(task60);
  densityq->add_task(task221);

  auto tensor222 = vector<shared_ptr<Tensor>>{I178, Gamma62_(), t2};
  auto task222 = make_shared<Task222>(tensor222, pindex);
  task221->add_dep(task222);
  task222->add_dep(task60);
  densityq->add_task(task222);

  vector<IndexRange> I189_index = {virt_, active_};
  auto I189 = make_shared<Tensor>(I189_index);
  auto tensor223 = vector<shared_ptr<Tensor>>{den2, I189};
  auto task223 = make_shared<Task223>(tensor223, pindex);
  task223->add_dep(task60);
  densityq->add_task(task223);

  vector<IndexRange> I190_index = {active_, virt_};
  auto I190 = make_shared<Tensor>(I190_index);
  auto tensor224 = vector<shared_ptr<Tensor>>{I189, Gamma16_(), I190};
  auto task224 = make_shared<Task224>(tensor224, pindex);
  task223->add_dep(task224);
  task224->add_dep(task60);
  densityq->add_task(task224);

  auto tensor225 = vector<shared_ptr<Tensor>>{I190, t2, l2};
  auto task225 = make_shared<Task225>(tensor225, pindex);
  task224->add_dep(task225);
  task225->add_dep(task60);
  densityq->add_task(task225);

  vector<IndexRange> I192_index = {virt_, active_};
  auto I192 = make_shared<Tensor>(I192_index);
  auto tensor226 = vector<shared_ptr<Tensor>>{den2, I192};
  auto task226 = make_shared<Task226>(tensor226, pindex);
  task226->add_dep(task60);
  densityq->add_task(task226);

  vector<IndexRange> I193_index = {active_, virt_};
  auto I193 = make_shared<Tensor>(I193_index);
  auto tensor227 = vector<shared_ptr<Tensor>>{I192, Gamma16_(), I193};
  auto task227 = make_shared<Task227>(tensor227, pindex);
  task226->add_dep(task227);
  task227->add_dep(task60);
  densityq->add_task(task227);

  auto tensor228 = vector<shared_ptr<Tensor>>{I193, t2, l2};
  auto task228 = make_shared<Task228>(tensor228, pindex);
  task227->add_dep(task228);
  task228->add_dep(task60);
  densityq->add_task(task228);

  vector<IndexRange> I198_index = {virt_, closed_};
  auto I198 = make_shared<Tensor>(I198_index);
  auto tensor229 = vector<shared_ptr<Tensor>>{den2, I198};
  auto task229 = make_shared<Task229>(tensor229, pindex);
  task229->add_dep(task60);
  densityq->add_task(task229);

  vector<IndexRange> I199_index = {virt_, closed_};
  auto I199 = make_shared<Tensor>(I199_index);
  auto tensor230 = vector<shared_ptr<Tensor>>{I198, t2, I199};
  auto task230 = make_shared<Task230>(tensor230, pindex);
  task229->add_dep(task230);
  task230->add_dep(task60);
  densityq->add_task(task230);

  vector<IndexRange> I200_index = {active_, virt_, closed_, active_};
  auto I200 = make_shared<Tensor>(I200_index);
  auto tensor231 = vector<shared_ptr<Tensor>>{I199, Gamma65_(), I200};
  auto task231 = make_shared<Task231>(tensor231, pindex);
  task230->add_dep(task231);
  task231->add_dep(task60);
  densityq->add_task(task231);

  auto tensor232 = vector<shared_ptr<Tensor>>{I200, l2};
  auto task232 = make_shared<Task232>(tensor232, pindex);
  task231->add_dep(task232);
  task232->add_dep(task60);
  densityq->add_task(task232);

  vector<IndexRange> I207_index = {active_, active_};
  auto I207 = make_shared<Tensor>(I207_index);
  auto tensor233 = vector<shared_ptr<Tensor>>{den2, I207};
  auto task233 = make_shared<Task233>(tensor233, pindex);
  task233->add_dep(task60);
  densityq->add_task(task233);

  vector<IndexRange> I208_index;
  auto I208 = make_shared<Tensor>(I208_index);
  auto tensor234 = vector<shared_ptr<Tensor>>{I207, Gamma65_(), I208};
  auto task234 = make_shared<Task234>(tensor234, pindex);
  task233->add_dep(task234);
  task234->add_dep(task60);
  densityq->add_task(task234);

  vector<IndexRange> I209_index = {closed_, virt_, closed_, virt_};
  auto I209 = make_shared<Tensor>(I209_index);
  auto tensor235 = vector<shared_ptr<Tensor>>{I208, t2, I209};
  auto task235 = make_shared<Task235>(tensor235, pindex);
  task234->add_dep(task235);
  task235->add_dep(task60);
  densityq->add_task(task235);

  auto tensor236 = vector<shared_ptr<Tensor>>{I209, l2};
  auto task236 = make_shared<Task236>(tensor236, pindex);
  task235->add_dep(task236);
  task236->add_dep(task60);
  densityq->add_task(task236);

  shared_ptr<Tensor> I213;
  if (diagonal) {
    vector<IndexRange> I213_index = {closed_, closed_};
    I213 = make_shared<Tensor>(I213_index);
  }
  shared_ptr<Task237> task237;
  if (diagonal) {
    auto tensor237 = vector<shared_ptr<Tensor>>{den2, I213};
    task237 = make_shared<Task237>(tensor237, pindex);
    task237->add_dep(task60);
    densityq->add_task(task237);
  }

  shared_ptr<Tensor> I214;
  if (diagonal) {
    vector<IndexRange> I214_index = {closed_, virt_, closed_, virt_};
    I214 = make_shared<Tensor>(I214_index);
  }
  shared_ptr<Task238> task238;
  if (diagonal) {
    auto tensor238 = vector<shared_ptr<Tensor>>{I213, t2, I214};
    task238 = make_shared<Task238>(tensor238, pindex);
    task237->add_dep(task238);
    task238->add_dep(task60);
    densityq->add_task(task238);
  }

  shared_ptr<Task239> task239;
  if (diagonal) {
    auto tensor239 = vector<shared_ptr<Tensor>>{I214, l2};
    task239 = make_shared<Task239>(tensor239, pindex);
    task238->add_dep(task239);
    task239->add_dep(task60);
    densityq->add_task(task239);
  }

  shared_ptr<Tensor> I217;
  if (diagonal) {
    vector<IndexRange> I217_index = {virt_, virt_};
    I217 = make_shared<Tensor>(I217_index);
  }
  shared_ptr<Task240> task240;
  if (diagonal) {
    auto tensor240 = vector<shared_ptr<Tensor>>{den2, I217};
    task240 = make_shared<Task240>(tensor240, pindex);
    task240->add_dep(task60);
    densityq->add_task(task240);
  }

  shared_ptr<Tensor> I218;
  if (diagonal) {
    vector<IndexRange> I218_index = {closed_, virt_, closed_, virt_};
    I218 = make_shared<Tensor>(I218_index);
  }
  shared_ptr<Task241> task241;
  if (diagonal) {
    auto tensor241 = vector<shared_ptr<Tensor>>{I217, t2, I218};
    task241 = make_shared<Task241>(tensor241, pindex);
    task240->add_dep(task241);
    task241->add_dep(task60);
    densityq->add_task(task241);
  }

  shared_ptr<Task242> task242;
  if (diagonal) {
    auto tensor242 = vector<shared_ptr<Tensor>>{I218, l2};
    task242 = make_shared<Task242>(tensor242, pindex);
    task241->add_dep(task242);
    task242->add_dep(task60);
    densityq->add_task(task242);
  }

  vector<IndexRange> I221_index = {closed_, active_};
  auto I221 = make_shared<Tensor>(I221_index);
  auto tensor243 = vector<shared_ptr<Tensor>>{den2, I221};
  auto task243 = make_shared<Task243>(tensor243, pindex);
  task243->add_dep(task60);
  densityq->add_task(task243);

  vector<IndexRange> I222_index = {active_, closed_};
  auto I222 = make_shared<Tensor>(I222_index);
  auto tensor244 = vector<shared_ptr<Tensor>>{I221, Gamma65_(), I222};
  auto task244 = make_shared<Task244>(tensor244, pindex);
  task243->add_dep(task244);
  task244->add_dep(task60);
  densityq->add_task(task244);

  vector<IndexRange> I223_index = {active_, virt_, closed_, virt_};
  auto I223 = make_shared<Tensor>(I223_index);
  auto tensor245 = vector<shared_ptr<Tensor>>{I222, t2, I223};
  auto task245 = make_shared<Task245>(tensor245, pindex);
  task244->add_dep(task245);
  task245->add_dep(task60);
  densityq->add_task(task245);

  auto tensor246 = vector<shared_ptr<Tensor>>{I223, l2};
  auto task246 = make_shared<Task246>(tensor246, pindex);
  task245->add_dep(task246);
  task246->add_dep(task60);
  densityq->add_task(task246);

  vector<IndexRange> I227_index = {active_, virt_};
  auto I227 = make_shared<Tensor>(I227_index);
  auto tensor247 = vector<shared_ptr<Tensor>>{den2, I227};
  auto task247 = make_shared<Task247>(tensor247, pindex);
  task247->add_dep(task60);
  densityq->add_task(task247);

  vector<IndexRange> I228_index = {virt_, closed_, active_, active_};
  auto I228 = make_shared<Tensor>(I228_index);
  auto tensor248 = vector<shared_ptr<Tensor>>{I227, t2, I228};
  auto task248 = make_shared<Task248>(tensor248, pindex);
  task247->add_dep(task248);
  task248->add_dep(task60);
  densityq->add_task(task248);

  vector<IndexRange> I229_index = {active_, virt_, closed_, active_};
  auto I229 = make_shared<Tensor>(I229_index);
  auto tensor249 = vector<shared_ptr<Tensor>>{I228, Gamma35_(), I229};
  auto task249 = make_shared<Task249>(tensor249, pindex);
  task248->add_dep(task249);
  task249->add_dep(task60);
  densityq->add_task(task249);

  auto tensor250 = vector<shared_ptr<Tensor>>{I229, l2};
  auto task250 = make_shared<Task250>(tensor250, pindex);
  task249->add_dep(task250);
  task250->add_dep(task60);
  densityq->add_task(task250);

  vector<IndexRange> I230_index = {active_, virt_};
  auto I230 = make_shared<Tensor>(I230_index);
  auto tensor251 = vector<shared_ptr<Tensor>>{den2, I230};
  auto task251 = make_shared<Task251>(tensor251, pindex);
  task251->add_dep(task60);
  densityq->add_task(task251);

  vector<IndexRange> I231_index = {virt_, closed_, active_, active_};
  auto I231 = make_shared<Tensor>(I231_index);
  auto tensor252 = vector<shared_ptr<Tensor>>{I230, t2, I231};
  auto task252 = make_shared<Task252>(tensor252, pindex);
  task251->add_dep(task252);
  task252->add_dep(task60);
  densityq->add_task(task252);

  auto tensor253 = vector<shared_ptr<Tensor>>{I231, Gamma32_(), l2};
  auto task253 = make_shared<Task253>(tensor253, pindex);
  task252->add_dep(task253);
  task253->add_dep(task60);
  densityq->add_task(task253);

  auto tensor254 = vector<shared_ptr<Tensor>>{I231, Gamma35_(), l2};
  auto task254 = make_shared<Task254>(tensor254, pindex);
  task252->add_dep(task254);
  task254->add_dep(task60);
  densityq->add_task(task254);

  vector<IndexRange> I239_index = {closed_, virt_};
  auto I239 = make_shared<Tensor>(I239_index);
  auto tensor255 = vector<shared_ptr<Tensor>>{den2, I239};
  auto task255 = make_shared<Task255>(tensor255, pindex);
  task255->add_dep(task60);
  densityq->add_task(task255);

  vector<IndexRange> I240_index = {virt_, active_};
  auto I240 = make_shared<Tensor>(I240_index);
  auto tensor256 = vector<shared_ptr<Tensor>>{I239, t2, I240};
  auto task256 = make_shared<Task256>(tensor256, pindex);
  task255->add_dep(task256);
  task256->add_dep(task60);
  densityq->add_task(task256);

  auto tensor257 = vector<shared_ptr<Tensor>>{I240, Gamma60_(), l2};
  auto task257 = make_shared<Task257>(tensor257, pindex);
  task256->add_dep(task257);
  task257->add_dep(task60);
  densityq->add_task(task257);

  vector<IndexRange> I242_index = {virt_, closed_};
  auto I242 = make_shared<Tensor>(I242_index);
  auto tensor258 = vector<shared_ptr<Tensor>>{den2, I242};
  auto task258 = make_shared<Task258>(tensor258, pindex);
  task258->add_dep(task60);
  densityq->add_task(task258);

  vector<IndexRange> I243_index = {virt_, active_};
  auto I243 = make_shared<Tensor>(I243_index);
  auto tensor259 = vector<shared_ptr<Tensor>>{I242, t2, I243};
  auto task259 = make_shared<Task259>(tensor259, pindex);
  task258->add_dep(task259);
  task259->add_dep(task60);
  densityq->add_task(task259);

  auto tensor260 = vector<shared_ptr<Tensor>>{I243, Gamma60_(), l2};
  auto task260 = make_shared<Task260>(tensor260, pindex);
  task259->add_dep(task260);
  task260->add_dep(task60);
  densityq->add_task(task260);

  vector<IndexRange> I245_index = {closed_, active_};
  auto I245 = make_shared<Tensor>(I245_index);
  auto tensor261 = vector<shared_ptr<Tensor>>{den2, I245};
  auto task261 = make_shared<Task261>(tensor261, pindex);
  task261->add_dep(task60);
  densityq->add_task(task261);

  vector<IndexRange> I246_index = {closed_, active_};
  auto I246 = make_shared<Tensor>(I246_index);
  auto tensor262 = vector<shared_ptr<Tensor>>{I245, Gamma65_(), I246};
  auto task262 = make_shared<Task262>(tensor262, pindex);
  task261->add_dep(task262);
  task262->add_dep(task60);
  densityq->add_task(task262);

  vector<IndexRange> I247_index = {closed_, virt_, closed_, virt_};
  auto I247 = make_shared<Tensor>(I247_index);
  auto tensor263 = vector<shared_ptr<Tensor>>{I246, t2, I247};
  auto task263 = make_shared<Task263>(tensor263, pindex);
  task262->add_dep(task263);
  task263->add_dep(task60);
  densityq->add_task(task263);

  auto tensor264 = vector<shared_ptr<Tensor>>{I247, l2};
  auto task264 = make_shared<Task264>(tensor264, pindex);
  task263->add_dep(task264);
  task264->add_dep(task60);
  densityq->add_task(task264);

  vector<IndexRange> I257_index = {closed_, closed_};
  auto I257 = make_shared<Tensor>(I257_index);
  auto tensor265 = vector<shared_ptr<Tensor>>{den2, I257};
  auto task265 = make_shared<Task265>(tensor265, pindex);
  task265->add_dep(task60);
  densityq->add_task(task265);

  vector<IndexRange> I258_index = {virt_, closed_, virt_, active_};
  auto I258 = make_shared<Tensor>(I258_index);
  auto tensor266 = vector<shared_ptr<Tensor>>{I257, l2, I258};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task265->add_dep(task266);
  task266->add_dep(task60);
  densityq->add_task(task266);

  auto tensor267 = vector<shared_ptr<Tensor>>{I258, Gamma65_(), t2};
  auto task267 = make_shared<Task267>(tensor267, pindex);
  task266->add_dep(task267);
  task267->add_dep(task60);
  densityq->add_task(task267);

  vector<IndexRange> I261_index = {virt_, closed_, virt_, active_};
  auto I261 = make_shared<Tensor>(I261_index);
  auto tensor268 = vector<shared_ptr<Tensor>>{I257, l2, I261};
  auto task268 = make_shared<Task268>(tensor268, pindex);
  task265->add_dep(task268);
  task268->add_dep(task60);
  densityq->add_task(task268);

  auto tensor269 = vector<shared_ptr<Tensor>>{I261, Gamma65_(), t2};
  auto task269 = make_shared<Task269>(tensor269, pindex);
  task268->add_dep(task269);
  task269->add_dep(task60);
  densityq->add_task(task269);

  vector<IndexRange> I263_index = {virt_, virt_};
  auto I263 = make_shared<Tensor>(I263_index);
  auto tensor270 = vector<shared_ptr<Tensor>>{den2, I263};
  auto task270 = make_shared<Task270>(tensor270, pindex);
  task270->add_dep(task60);
  densityq->add_task(task270);

  vector<IndexRange> I264_index = {virt_, closed_, virt_, active_};
  auto I264 = make_shared<Tensor>(I264_index);
  auto tensor271 = vector<shared_ptr<Tensor>>{I263, l2, I264};
  auto task271 = make_shared<Task271>(tensor271, pindex);
  task270->add_dep(task271);
  task271->add_dep(task60);
  densityq->add_task(task271);

  auto tensor272 = vector<shared_ptr<Tensor>>{I264, Gamma65_(), t2};
  auto task272 = make_shared<Task272>(tensor272, pindex);
  task271->add_dep(task272);
  task272->add_dep(task60);
  densityq->add_task(task272);

  vector<IndexRange> I267_index = {virt_, closed_, virt_, active_};
  auto I267 = make_shared<Tensor>(I267_index);
  auto tensor273 = vector<shared_ptr<Tensor>>{I263, l2, I267};
  auto task273 = make_shared<Task273>(tensor273, pindex);
  task270->add_dep(task273);
  task273->add_dep(task60);
  densityq->add_task(task273);

  auto tensor274 = vector<shared_ptr<Tensor>>{I267, Gamma65_(), t2};
  auto task274 = make_shared<Task274>(tensor274, pindex);
  task273->add_dep(task274);
  task274->add_dep(task60);
  densityq->add_task(task274);

  vector<IndexRange> I269_index = {virt_, virt_};
  auto I269 = make_shared<Tensor>(I269_index);
  auto tensor275 = vector<shared_ptr<Tensor>>{den2, I269};
  auto task275 = make_shared<Task275>(tensor275, pindex);
  task275->add_dep(task60);
  densityq->add_task(task275);

  vector<IndexRange> I270_index = {virt_, closed_, virt_, active_};
  auto I270 = make_shared<Tensor>(I270_index);
  auto tensor276 = vector<shared_ptr<Tensor>>{I269, l2, I270};
  auto task276 = make_shared<Task276>(tensor276, pindex);
  task275->add_dep(task276);
  task276->add_dep(task60);
  densityq->add_task(task276);

  auto tensor277 = vector<shared_ptr<Tensor>>{I270, Gamma65_(), t2};
  auto task277 = make_shared<Task277>(tensor277, pindex);
  task276->add_dep(task277);
  task277->add_dep(task60);
  densityq->add_task(task277);

  vector<IndexRange> I273_index = {virt_, closed_, virt_, active_};
  auto I273 = make_shared<Tensor>(I273_index);
  auto tensor278 = vector<shared_ptr<Tensor>>{I269, l2, I273};
  auto task278 = make_shared<Task278>(tensor278, pindex);
  task275->add_dep(task278);
  task278->add_dep(task60);
  densityq->add_task(task278);

  auto tensor279 = vector<shared_ptr<Tensor>>{I273, Gamma65_(), t2};
  auto task279 = make_shared<Task279>(tensor279, pindex);
  task278->add_dep(task279);
  task279->add_dep(task60);
  densityq->add_task(task279);

  vector<IndexRange> I275_index = {active_, closed_};
  auto I275 = make_shared<Tensor>(I275_index);
  auto tensor280 = vector<shared_ptr<Tensor>>{den2, I275};
  auto task280 = make_shared<Task280>(tensor280, pindex);
  task280->add_dep(task60);
  densityq->add_task(task280);

  vector<IndexRange> I276_index = {virt_, virt_, active_, active_};
  auto I276 = make_shared<Tensor>(I276_index);
  auto tensor281 = vector<shared_ptr<Tensor>>{I275, t2, I276};
  auto task281 = make_shared<Task281>(tensor281, pindex);
  task280->add_dep(task281);
  task281->add_dep(task60);
  densityq->add_task(task281);

  auto tensor282 = vector<shared_ptr<Tensor>>{I276, Gamma60_(), l2};
  auto task282 = make_shared<Task282>(tensor282, pindex);
  task281->add_dep(task282);
  task282->add_dep(task60);
  densityq->add_task(task282);

  return densityq;
}


#endif
