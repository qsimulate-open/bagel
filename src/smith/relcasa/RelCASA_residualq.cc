//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASA_residualqq.cc
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

shared_ptr<Queue> RelCASA::RelCASA::make_residualq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto residualq = make_shared<Queue>();
  auto tensor40 = vector<shared_ptr<Tensor>>{r};
  auto task40 = make_shared<Task40>(tensor40, reset);
  residualq->add_task(task40);

  vector<IndexRange> I0_index = {closed_, active_, active_, closed_};
  auto I0 = make_shared<Tensor>(I0_index);
  auto tensor41 = vector<shared_ptr<Tensor>>{r, I0};
  auto task41 = make_shared<Task41>(tensor41, pindex);
  task41->add_dep(task40);
  residualq->add_task(task41);

  vector<IndexRange> I1_index = {closed_, active_, active_, active_};
  auto I1 = make_shared<Tensor>(I1_index);
  auto tensor42 = vector<shared_ptr<Tensor>>{I0, f1_, I1};
  auto task42 = make_shared<Task42>(tensor42, pindex);
  task41->add_dep(task42);
  task42->add_dep(task40);
  residualq->add_task(task42);

  auto tensor43 = vector<shared_ptr<Tensor>>{I1, Gamma0_(), t2};
  auto task43 = make_shared<Task43>(tensor43, pindex);
  task42->add_dep(task43);
  task43->add_dep(task40);
  residualq->add_task(task43);

  vector<IndexRange> I4_index = {closed_, closed_, active_, active_};
  auto I4 = make_shared<Tensor>(I4_index);
  auto tensor44 = vector<shared_ptr<Tensor>>{I0, Gamma1_(), I4};
  auto task44 = make_shared<Task44>(tensor44, pindex);
  task41->add_dep(task44);
  task44->add_dep(task40);
  residualq->add_task(task44);

  auto tensor45 = vector<shared_ptr<Tensor>>{I4, f1_, t2};
  auto task45 = make_shared<Task45>(tensor45, pindex);
  task44->add_dep(task45);
  task45->add_dep(task40);
  residualq->add_task(task45);

  vector<IndexRange> I6_index = {closed_, active_, active_, active_};
  auto I6 = make_shared<Tensor>(I6_index);
  auto tensor46 = vector<shared_ptr<Tensor>>{r, I6};
  auto task46 = make_shared<Task46>(tensor46, pindex);
  task46->add_dep(task40);
  residualq->add_task(task46);

  vector<IndexRange> I7_index = {closed_, active_, active_, active_};
  auto I7 = make_shared<Tensor>(I7_index);
  auto tensor47 = vector<shared_ptr<Tensor>>{I6, Gamma2_(), I7};
  auto task47 = make_shared<Task47>(tensor47, pindex);
  task46->add_dep(task47);
  task47->add_dep(task40);
  residualq->add_task(task47);

  auto tensor48 = vector<shared_ptr<Tensor>>{I7, f1_, t2};
  auto task48 = make_shared<Task48>(tensor48, pindex);
  task47->add_dep(task48);
  task48->add_dep(task40);
  residualq->add_task(task48);

  vector<IndexRange> I10_index = {closed_, active_};
  auto I10 = make_shared<Tensor>(I10_index);
  auto tensor49 = vector<shared_ptr<Tensor>>{I6, Gamma3_(), I10};
  auto task49 = make_shared<Task49>(tensor49, pindex);
  task46->add_dep(task49);
  task49->add_dep(task40);
  residualq->add_task(task49);

  vector<IndexRange> I11_index = {closed_, virt_, closed_, active_};
  auto I11 = make_shared<Tensor>(I11_index);
  auto tensor50 = vector<shared_ptr<Tensor>>{I10, f1_, I11};
  auto task50 = make_shared<Task50>(tensor50, pindex);
  task49->add_dep(task50);
  task50->add_dep(task40);
  residualq->add_task(task50);

  auto tensor51 = vector<shared_ptr<Tensor>>{I11, t2};
  auto task51 = make_shared<Task51>(tensor51, pindex);
  task50->add_dep(task51);
  task51->add_dep(task40);
  residualq->add_task(task51);

  vector<IndexRange> I16_index = {active_, closed_, active_, active_};
  auto I16 = make_shared<Tensor>(I16_index);
  auto tensor52 = vector<shared_ptr<Tensor>>{I6, Gamma5_(), I16};
  auto task52 = make_shared<Task52>(tensor52, pindex);
  task46->add_dep(task52);
  task52->add_dep(task40);
  residualq->add_task(task52);

  auto tensor53 = vector<shared_ptr<Tensor>>{I16, f1_, t2};
  auto task53 = make_shared<Task53>(tensor53, pindex);
  task52->add_dep(task53);
  task53->add_dep(task40);
  residualq->add_task(task53);

  vector<IndexRange> I19_index = {closed_, active_, active_, active_};
  auto I19 = make_shared<Tensor>(I19_index);
  auto tensor54 = vector<shared_ptr<Tensor>>{I6, Gamma6_(), I19};
  auto task54 = make_shared<Task54>(tensor54, pindex);
  task46->add_dep(task54);
  task54->add_dep(task40);
  residualq->add_task(task54);

  auto tensor55 = vector<shared_ptr<Tensor>>{I19, f1_, t2};
  auto task55 = make_shared<Task55>(tensor55, pindex);
  task54->add_dep(task55);
  task55->add_dep(task40);
  residualq->add_task(task55);

  auto tensor56 = vector<shared_ptr<Tensor>>{I6, Gamma70_(), t2};
  auto task56 = make_shared<Task56>(tensor56, pindex);
  task46->add_dep(task56);
  task56->add_dep(task40);
  residualq->add_task(task56);

  auto tensor57 = vector<shared_ptr<Tensor>>{I6, Gamma71_(), t2};
  auto task57 = make_shared<Task57>(tensor57, pindex);
  task46->add_dep(task57);
  task57->add_dep(task40);
  residualq->add_task(task57);

  vector<IndexRange> I21_index = {virt_, closed_, active_, closed_};
  auto I21 = make_shared<Tensor>(I21_index);
  auto tensor58 = vector<shared_ptr<Tensor>>{r, I21};
  auto task58 = make_shared<Task58>(tensor58, pindex);
  task58->add_dep(task40);
  residualq->add_task(task58);

  vector<IndexRange> I22_index = {virt_, closed_, active_, active_};
  auto I22 = make_shared<Tensor>(I22_index);
  auto tensor59 = vector<shared_ptr<Tensor>>{I21, f1_, I22};
  auto task59 = make_shared<Task59>(tensor59, pindex);
  task58->add_dep(task59);
  task59->add_dep(task40);
  residualq->add_task(task59);

  auto tensor60 = vector<shared_ptr<Tensor>>{I22, Gamma7_(), t2};
  auto task60 = make_shared<Task60>(tensor60, pindex);
  task59->add_dep(task60);
  task60->add_dep(task40);
  residualq->add_task(task60);

  auto tensor61 = vector<shared_ptr<Tensor>>{I22, Gamma9_(), t2};
  auto task61 = make_shared<Task61>(tensor61, pindex);
  task59->add_dep(task61);
  task61->add_dep(task40);
  residualq->add_task(task61);

  vector<IndexRange> I25_index = {virt_, closed_, active_, active_};
  auto I25 = make_shared<Tensor>(I25_index);
  auto tensor62 = vector<shared_ptr<Tensor>>{I21, f1_, I25};
  auto task62 = make_shared<Task62>(tensor62, pindex);
  task58->add_dep(task62);
  task62->add_dep(task40);
  residualq->add_task(task62);

  vector<IndexRange> I26_index = {active_, virt_, closed_, active_};
  auto I26 = make_shared<Tensor>(I26_index);
  auto tensor63 = vector<shared_ptr<Tensor>>{I25, Gamma9_(), I26};
  auto task63 = make_shared<Task63>(tensor63, pindex);
  task62->add_dep(task63);
  task63->add_dep(task40);
  residualq->add_task(task63);

  auto tensor64 = vector<shared_ptr<Tensor>>{I26, t2};
  auto task64 = make_shared<Task64>(tensor64, pindex);
  task63->add_dep(task64);
  task64->add_dep(task40);
  residualq->add_task(task64);

  vector<IndexRange> I34_index = {closed_, active_};
  auto I34 = make_shared<Tensor>(I34_index);
  auto tensor65 = vector<shared_ptr<Tensor>>{I21, f1_, I34};
  auto task65 = make_shared<Task65>(tensor65, pindex);
  task58->add_dep(task65);
  task65->add_dep(task40);
  residualq->add_task(task65);

  auto tensor66 = vector<shared_ptr<Tensor>>{I34, Gamma9_(), t2};
  auto task66 = make_shared<Task66>(tensor66, pindex);
  task65->add_dep(task66);
  task66->add_dep(task40);
  residualq->add_task(task66);

  vector<IndexRange> I37_index = {closed_, active_};
  auto I37 = make_shared<Tensor>(I37_index);
  auto tensor67 = vector<shared_ptr<Tensor>>{I21, f1_, I37};
  auto task67 = make_shared<Task67>(tensor67, pindex);
  task58->add_dep(task67);
  task67->add_dep(task40);
  residualq->add_task(task67);

  auto tensor68 = vector<shared_ptr<Tensor>>{I37, Gamma9_(), t2};
  auto task68 = make_shared<Task68>(tensor68, pindex);
  task67->add_dep(task68);
  task68->add_dep(task40);
  residualq->add_task(task68);

  vector<IndexRange> I40_index = {virt_, active_};
  auto I40 = make_shared<Tensor>(I40_index);
  auto tensor69 = vector<shared_ptr<Tensor>>{I21, t2, I40};
  auto task69 = make_shared<Task69>(tensor69, pindex);
  task58->add_dep(task69);
  task69->add_dep(task40);
  residualq->add_task(task69);

  auto tensor70 = vector<shared_ptr<Tensor>>{I40, Gamma13_(), f1_};
  auto task70 = make_shared<Task70>(tensor70, pindex);
  task69->add_dep(task70);
  task70->add_dep(task40);
  residualq->add_task(task70);

  vector<IndexRange> I43_index = {virt_, active_};
  auto I43 = make_shared<Tensor>(I43_index);
  auto tensor71 = vector<shared_ptr<Tensor>>{I21, t2, I43};
  auto task71 = make_shared<Task71>(tensor71, pindex);
  task58->add_dep(task71);
  task71->add_dep(task40);
  residualq->add_task(task71);

  auto tensor72 = vector<shared_ptr<Tensor>>{I43, Gamma13_(), f1_};
  auto task72 = make_shared<Task72>(tensor72, pindex);
  task71->add_dep(task72);
  task72->add_dep(task40);
  residualq->add_task(task72);

  vector<IndexRange> I46_index = {closed_, closed_, active_, active_};
  auto I46 = make_shared<Tensor>(I46_index);
  auto tensor73 = vector<shared_ptr<Tensor>>{I21, f1_, I46};
  auto task73 = make_shared<Task73>(tensor73, pindex);
  task58->add_dep(task73);
  task73->add_dep(task40);
  residualq->add_task(task73);

  auto tensor74 = vector<shared_ptr<Tensor>>{I46, Gamma1_(), t2};
  auto task74 = make_shared<Task74>(tensor74, pindex);
  task73->add_dep(task74);
  task74->add_dep(task40);
  residualq->add_task(task74);

  vector<IndexRange> I49_index = {closed_, virt_, closed_, active_};
  auto I49 = make_shared<Tensor>(I49_index);
  auto tensor75 = vector<shared_ptr<Tensor>>{I21, f1_, I49};
  auto task75 = make_shared<Task75>(tensor75, pindex);
  task58->add_dep(task75);
  task75->add_dep(task40);
  residualq->add_task(task75);

  vector<IndexRange> I50_index = {closed_, virt_, closed_, active_};
  auto I50 = make_shared<Tensor>(I50_index);
  auto tensor76 = vector<shared_ptr<Tensor>>{I49, Gamma13_(), I50};
  auto task76 = make_shared<Task76>(tensor76, pindex);
  task75->add_dep(task76);
  task76->add_dep(task40);
  residualq->add_task(task76);

  auto tensor77 = vector<shared_ptr<Tensor>>{I50, t2};
  auto task77 = make_shared<Task77>(tensor77, pindex);
  task76->add_dep(task77);
  task77->add_dep(task40);
  residualq->add_task(task77);

  vector<IndexRange> I217_index = {closed_, virt_, closed_, active_};
  auto I217 = make_shared<Tensor>(I217_index);
  auto tensor78 = vector<shared_ptr<Tensor>>{I21, Gamma72_(), I217};
  auto task78 = make_shared<Task78>(tensor78, pindex);
  task58->add_dep(task78);
  task78->add_dep(task40);
  residualq->add_task(task78);

  auto tensor79 = vector<shared_ptr<Tensor>>{I217, t2};
  auto task79 = make_shared<Task79>(tensor79, pindex);
  task78->add_dep(task79);
  task79->add_dep(task40);
  residualq->add_task(task79);

  vector<IndexRange> I221_index = {closed_, virt_, closed_, active_};
  auto I221 = make_shared<Tensor>(I221_index);
  auto tensor80 = vector<shared_ptr<Tensor>>{I21, Gamma74_(), I221};
  auto task80 = make_shared<Task80>(tensor80, pindex);
  task58->add_dep(task80);
  task80->add_dep(task40);
  residualq->add_task(task80);

  auto tensor81 = vector<shared_ptr<Tensor>>{I221, t2};
  auto task81 = make_shared<Task81>(tensor81, pindex);
  task80->add_dep(task81);
  task81->add_dep(task40);
  residualq->add_task(task81);

  vector<IndexRange> I54_index = {virt_, closed_, active_, active_};
  auto I54 = make_shared<Tensor>(I54_index);
  auto tensor82 = vector<shared_ptr<Tensor>>{r, I54};
  auto task82 = make_shared<Task82>(tensor82, pindex);
  task82->add_dep(task40);
  residualq->add_task(task82);

  vector<IndexRange> I55_index = {virt_, closed_, active_, active_};
  auto I55 = make_shared<Tensor>(I55_index);
  auto tensor83 = vector<shared_ptr<Tensor>>{I54, Gamma18_(), I55};
  auto task83 = make_shared<Task83>(tensor83, pindex);
  task82->add_dep(task83);
  task83->add_dep(task40);
  residualq->add_task(task83);

  auto tensor84 = vector<shared_ptr<Tensor>>{I55, f1_, t2};
  auto task84 = make_shared<Task84>(tensor84, pindex);
  task83->add_dep(task84);
  task84->add_dep(task40);
  residualq->add_task(task84);

  vector<IndexRange> I58_index = {closed_, virt_, active_, active_};
  auto I58 = make_shared<Tensor>(I58_index);
  auto tensor85 = vector<shared_ptr<Tensor>>{I54, Gamma3_(), I58};
  auto task85 = make_shared<Task85>(tensor85, pindex);
  task82->add_dep(task85);
  task85->add_dep(task40);
  residualq->add_task(task85);

  auto tensor86 = vector<shared_ptr<Tensor>>{I58, f1_, t2};
  auto task86 = make_shared<Task86>(tensor86, pindex);
  task85->add_dep(task86);
  task86->add_dep(task40);
  residualq->add_task(task86);

  vector<IndexRange> I61_index = {virt_, active_, active_, active_};
  auto I61 = make_shared<Tensor>(I61_index);
  auto tensor87 = vector<shared_ptr<Tensor>>{I54, f1_, I61};
  auto task87 = make_shared<Task87>(tensor87, pindex);
  task82->add_dep(task87);
  task87->add_dep(task40);
  residualq->add_task(task87);

  auto tensor88 = vector<shared_ptr<Tensor>>{I61, Gamma20_(), t2};
  auto task88 = make_shared<Task88>(tensor88, pindex);
  task87->add_dep(task88);
  task88->add_dep(task40);
  residualq->add_task(task88);

  vector<IndexRange> I64_index = {closed_, virt_};
  auto I64 = make_shared<Tensor>(I64_index);
  auto tensor89 = vector<shared_ptr<Tensor>>{I54, Gamma21_(), I64};
  auto task89 = make_shared<Task89>(tensor89, pindex);
  task82->add_dep(task89);
  task89->add_dep(task40);
  residualq->add_task(task89);

  vector<IndexRange> I65_index = {closed_, virt_, closed_, virt_};
  auto I65 = make_shared<Tensor>(I65_index);
  auto tensor90 = vector<shared_ptr<Tensor>>{I64, f1_, I65};
  auto task90 = make_shared<Task90>(tensor90, pindex);
  task89->add_dep(task90);
  task90->add_dep(task40);
  residualq->add_task(task90);

  auto tensor91 = vector<shared_ptr<Tensor>>{I65, t2};
  auto task91 = make_shared<Task91>(tensor91, pindex);
  task90->add_dep(task91);
  task91->add_dep(task40);
  residualq->add_task(task91);

  vector<IndexRange> I70_index = {active_, closed_, virt_, active_};
  auto I70 = make_shared<Tensor>(I70_index);
  auto tensor92 = vector<shared_ptr<Tensor>>{I54, Gamma23_(), I70};
  auto task92 = make_shared<Task92>(tensor92, pindex);
  task82->add_dep(task92);
  task92->add_dep(task40);
  residualq->add_task(task92);

  auto tensor93 = vector<shared_ptr<Tensor>>{I70, f1_, t2};
  auto task93 = make_shared<Task93>(tensor93, pindex);
  task92->add_dep(task93);
  task93->add_dep(task40);
  residualq->add_task(task93);

  vector<IndexRange> I73_index = {active_, virt_, closed_, active_};
  auto I73 = make_shared<Tensor>(I73_index);
  auto tensor94 = vector<shared_ptr<Tensor>>{I54, Gamma24_(), I73};
  auto task94 = make_shared<Task94>(tensor94, pindex);
  task82->add_dep(task94);
  task94->add_dep(task40);
  residualq->add_task(task94);

  auto tensor95 = vector<shared_ptr<Tensor>>{I73, f1_, t2};
  auto task95 = make_shared<Task95>(tensor95, pindex);
  task94->add_dep(task95);
  task95->add_dep(task40);
  residualq->add_task(task95);

  vector<IndexRange> I76_index = {closed_, active_, active_, active_};
  auto I76 = make_shared<Tensor>(I76_index);
  auto tensor96 = vector<shared_ptr<Tensor>>{I54, f1_, I76};
  auto task96 = make_shared<Task96>(tensor96, pindex);
  task82->add_dep(task96);
  task96->add_dep(task40);
  residualq->add_task(task96);

  auto tensor97 = vector<shared_ptr<Tensor>>{I76, Gamma25_(), t2};
  auto task97 = make_shared<Task97>(tensor97, pindex);
  task96->add_dep(task97);
  task97->add_dep(task40);
  residualq->add_task(task97);

  vector<IndexRange> I79_index = {virt_, closed_, active_, active_};
  auto I79 = make_shared<Tensor>(I79_index);
  auto tensor98 = vector<shared_ptr<Tensor>>{I54, f1_, I79};
  auto task98 = make_shared<Task98>(tensor98, pindex);
  task82->add_dep(task98);
  task98->add_dep(task40);
  residualq->add_task(task98);

  auto tensor99 = vector<shared_ptr<Tensor>>{I79, Gamma24_(), t2};
  auto task99 = make_shared<Task99>(tensor99, pindex);
  task98->add_dep(task99);
  task99->add_dep(task40);
  residualq->add_task(task99);

  auto tensor100 = vector<shared_ptr<Tensor>>{I79, Gamma23_(), t2};
  auto task100 = make_shared<Task100>(tensor100, pindex);
  task98->add_dep(task100);
  task100->add_dep(task40);
  residualq->add_task(task100);

  auto tensor101 = vector<shared_ptr<Tensor>>{I54, Gamma88_(), t2};
  auto task101 = make_shared<Task101>(tensor101, pindex);
  task82->add_dep(task101);
  task101->add_dep(task40);
  residualq->add_task(task101);

  auto tensor102 = vector<shared_ptr<Tensor>>{I54, Gamma89_(), t2};
  auto task102 = make_shared<Task102>(tensor102, pindex);
  task82->add_dep(task102);
  task102->add_dep(task40);
  residualq->add_task(task102);

  auto tensor103 = vector<shared_ptr<Tensor>>{I54, Gamma92_(), t2};
  auto task103 = make_shared<Task103>(tensor103, pindex);
  task82->add_dep(task103);
  task103->add_dep(task40);
  residualq->add_task(task103);

  auto tensor104 = vector<shared_ptr<Tensor>>{I54, Gamma93_(), t2};
  auto task104 = make_shared<Task104>(tensor104, pindex);
  task82->add_dep(task104);
  task104->add_dep(task40);
  residualq->add_task(task104);

  vector<IndexRange> I84_index = {virt_, closed_, active_, active_};
  auto I84 = make_shared<Tensor>(I84_index);
  auto tensor105 = vector<shared_ptr<Tensor>>{r, I84};
  auto task105 = make_shared<Task105>(tensor105, pindex);
  task105->add_dep(task40);
  residualq->add_task(task105);

  vector<IndexRange> I85_index = {virt_, closed_, active_, active_};
  auto I85 = make_shared<Tensor>(I85_index);
  auto tensor106 = vector<shared_ptr<Tensor>>{I84, Gamma3_(), I85};
  auto task106 = make_shared<Task106>(tensor106, pindex);
  task105->add_dep(task106);
  task106->add_dep(task40);
  residualq->add_task(task106);

  vector<IndexRange> I86_index = {closed_, virt_, closed_, active_};
  auto I86 = make_shared<Tensor>(I86_index);
  auto tensor107 = vector<shared_ptr<Tensor>>{I85, f1_, I86};
  auto task107 = make_shared<Task107>(tensor107, pindex);
  task106->add_dep(task107);
  task107->add_dep(task40);
  residualq->add_task(task107);

  auto tensor108 = vector<shared_ptr<Tensor>>{I86, t2};
  auto task108 = make_shared<Task108>(tensor108, pindex);
  task107->add_dep(task108);
  task108->add_dep(task40);
  residualq->add_task(task108);

  vector<IndexRange> I91_index = {virt_, active_, active_, active_};
  auto I91 = make_shared<Tensor>(I91_index);
  auto tensor109 = vector<shared_ptr<Tensor>>{I84, f1_, I91};
  auto task109 = make_shared<Task109>(tensor109, pindex);
  task105->add_dep(task109);
  task109->add_dep(task40);
  residualq->add_task(task109);

  auto tensor110 = vector<shared_ptr<Tensor>>{I91, Gamma30_(), t2};
  auto task110 = make_shared<Task110>(tensor110, pindex);
  task109->add_dep(task110);
  task110->add_dep(task40);
  residualq->add_task(task110);

  vector<IndexRange> I94_index = {closed_, virt_};
  auto I94 = make_shared<Tensor>(I94_index);
  auto tensor111 = vector<shared_ptr<Tensor>>{I84, Gamma21_(), I94};
  auto task111 = make_shared<Task111>(tensor111, pindex);
  task105->add_dep(task111);
  task111->add_dep(task40);
  residualq->add_task(task111);

  vector<IndexRange> I95_index = {closed_, virt_, closed_, virt_};
  auto I95 = make_shared<Tensor>(I95_index);
  auto tensor112 = vector<shared_ptr<Tensor>>{I94, f1_, I95};
  auto task112 = make_shared<Task112>(tensor112, pindex);
  task111->add_dep(task112);
  task112->add_dep(task40);
  residualq->add_task(task112);

  auto tensor113 = vector<shared_ptr<Tensor>>{I95, t2};
  auto task113 = make_shared<Task113>(tensor113, pindex);
  task112->add_dep(task113);
  task113->add_dep(task40);
  residualq->add_task(task113);

  vector<IndexRange> I100_index = {active_, closed_, virt_, active_};
  auto I100 = make_shared<Tensor>(I100_index);
  auto tensor114 = vector<shared_ptr<Tensor>>{I84, Gamma23_(), I100};
  auto task114 = make_shared<Task114>(tensor114, pindex);
  task105->add_dep(task114);
  task114->add_dep(task40);
  residualq->add_task(task114);

  vector<IndexRange> I101_index = {active_, virt_, closed_, virt_};
  auto I101 = make_shared<Tensor>(I101_index);
  auto tensor115 = vector<shared_ptr<Tensor>>{I100, f1_, I101};
  auto task115 = make_shared<Task115>(tensor115, pindex);
  task114->add_dep(task115);
  task115->add_dep(task40);
  residualq->add_task(task115);

  auto tensor116 = vector<shared_ptr<Tensor>>{I101, t2};
  auto task116 = make_shared<Task116>(tensor116, pindex);
  task115->add_dep(task116);
  task116->add_dep(task40);
  residualq->add_task(task116);

  vector<IndexRange> I106_index = {closed_, active_, active_, active_};
  auto I106 = make_shared<Tensor>(I106_index);
  auto tensor117 = vector<shared_ptr<Tensor>>{I84, f1_, I106};
  auto task117 = make_shared<Task117>(tensor117, pindex);
  task105->add_dep(task117);
  task117->add_dep(task40);
  residualq->add_task(task117);

  auto tensor118 = vector<shared_ptr<Tensor>>{I106, Gamma6_(), t2};
  auto task118 = make_shared<Task118>(tensor118, pindex);
  task117->add_dep(task118);
  task118->add_dep(task40);
  residualq->add_task(task118);

  vector<IndexRange> I109_index = {virt_, closed_, active_, active_};
  auto I109 = make_shared<Tensor>(I109_index);
  auto tensor119 = vector<shared_ptr<Tensor>>{I84, f1_, I109};
  auto task119 = make_shared<Task119>(tensor119, pindex);
  task105->add_dep(task119);
  task119->add_dep(task40);
  residualq->add_task(task119);

  vector<IndexRange> I110_index = {active_, virt_, closed_, active_};
  auto I110 = make_shared<Tensor>(I110_index);
  auto tensor120 = vector<shared_ptr<Tensor>>{I109, Gamma23_(), I110};
  auto task120 = make_shared<Task120>(tensor120, pindex);
  task119->add_dep(task120);
  task120->add_dep(task40);
  residualq->add_task(task120);

  auto tensor121 = vector<shared_ptr<Tensor>>{I110, t2};
  auto task121 = make_shared<Task121>(tensor121, pindex);
  task120->add_dep(task121);
  task121->add_dep(task40);
  residualq->add_task(task121);

  vector<IndexRange> I253_index = {active_, virt_, closed_, active_};
  auto I253 = make_shared<Tensor>(I253_index);
  auto tensor122 = vector<shared_ptr<Tensor>>{I84, Gamma89_(), I253};
  auto task122 = make_shared<Task122>(tensor122, pindex);
  task105->add_dep(task122);
  task122->add_dep(task40);
  residualq->add_task(task122);

  auto tensor123 = vector<shared_ptr<Tensor>>{I253, t2};
  auto task123 = make_shared<Task123>(tensor123, pindex);
  task122->add_dep(task123);
  task123->add_dep(task40);
  residualq->add_task(task123);

  vector<IndexRange> I261_index = {active_, virt_, closed_, active_};
  auto I261 = make_shared<Tensor>(I261_index);
  auto tensor124 = vector<shared_ptr<Tensor>>{I84, Gamma93_(), I261};
  auto task124 = make_shared<Task124>(tensor124, pindex);
  task105->add_dep(task124);
  task124->add_dep(task40);
  residualq->add_task(task124);

  auto tensor125 = vector<shared_ptr<Tensor>>{I261, t2};
  auto task125 = make_shared<Task125>(tensor125, pindex);
  task124->add_dep(task125);
  task125->add_dep(task40);
  residualq->add_task(task125);

  vector<IndexRange> I114_index = {virt_, active_, active_, active_};
  auto I114 = make_shared<Tensor>(I114_index);
  auto tensor126 = vector<shared_ptr<Tensor>>{r, I114};
  auto task126 = make_shared<Task126>(tensor126, pindex);
  task126->add_dep(task40);
  residualq->add_task(task126);

  vector<IndexRange> I115_index = {active_, virt_, active_, active_};
  auto I115 = make_shared<Tensor>(I115_index);
  auto tensor127 = vector<shared_ptr<Tensor>>{I114, Gamma38_(), I115};
  auto task127 = make_shared<Task127>(tensor127, pindex);
  task126->add_dep(task127);
  task127->add_dep(task40);
  residualq->add_task(task127);

  auto tensor128 = vector<shared_ptr<Tensor>>{I115, f1_, t2};
  auto task128 = make_shared<Task128>(tensor128, pindex);
  task127->add_dep(task128);
  task128->add_dep(task40);
  residualq->add_task(task128);

  vector<IndexRange> I118_index = {virt_, active_, active_, active_};
  auto I118 = make_shared<Tensor>(I118_index);
  auto tensor129 = vector<shared_ptr<Tensor>>{I114, Gamma39_(), I118};
  auto task129 = make_shared<Task129>(tensor129, pindex);
  task126->add_dep(task129);
  task129->add_dep(task40);
  residualq->add_task(task129);

  auto tensor130 = vector<shared_ptr<Tensor>>{I118, f1_, t2};
  auto task130 = make_shared<Task130>(tensor130, pindex);
  task129->add_dep(task130);
  task130->add_dep(task40);
  residualq->add_task(task130);

  vector<IndexRange> I121_index = {active_, virt_};
  auto I121 = make_shared<Tensor>(I121_index);
  auto tensor131 = vector<shared_ptr<Tensor>>{I114, Gamma40_(), I121};
  auto task131 = make_shared<Task131>(tensor131, pindex);
  task126->add_dep(task131);
  task131->add_dep(task40);
  residualq->add_task(task131);

  vector<IndexRange> I122_index = {active_, virt_, closed_, virt_};
  auto I122 = make_shared<Tensor>(I122_index);
  auto tensor132 = vector<shared_ptr<Tensor>>{I121, f1_, I122};
  auto task132 = make_shared<Task132>(tensor132, pindex);
  task131->add_dep(task132);
  task132->add_dep(task40);
  residualq->add_task(task132);

  auto tensor133 = vector<shared_ptr<Tensor>>{I122, t2};
  auto task133 = make_shared<Task133>(tensor133, pindex);
  task132->add_dep(task133);
  task133->add_dep(task40);
  residualq->add_task(task133);

  vector<IndexRange> I127_index = {active_, virt_, active_, active_};
  auto I127 = make_shared<Tensor>(I127_index);
  auto tensor134 = vector<shared_ptr<Tensor>>{I114, Gamma42_(), I127};
  auto task134 = make_shared<Task134>(tensor134, pindex);
  task126->add_dep(task134);
  task134->add_dep(task40);
  residualq->add_task(task134);

  auto tensor135 = vector<shared_ptr<Tensor>>{I127, f1_, t2};
  auto task135 = make_shared<Task135>(tensor135, pindex);
  task134->add_dep(task135);
  task135->add_dep(task40);
  residualq->add_task(task135);

  vector<IndexRange> I130_index = {virt_, active_, active_, active_};
  auto I130 = make_shared<Tensor>(I130_index);
  auto tensor136 = vector<shared_ptr<Tensor>>{I114, f1_, I130};
  auto task136 = make_shared<Task136>(tensor136, pindex);
  task126->add_dep(task136);
  task136->add_dep(task40);
  residualq->add_task(task136);

  auto tensor137 = vector<shared_ptr<Tensor>>{I130, Gamma42_(), t2};
  auto task137 = make_shared<Task137>(tensor137, pindex);
  task136->add_dep(task137);
  task137->add_dep(task40);
  residualq->add_task(task137);

  auto tensor138 = vector<shared_ptr<Tensor>>{I114, Gamma76_(), t2};
  auto task138 = make_shared<Task138>(tensor138, pindex);
  task126->add_dep(task138);
  task138->add_dep(task40);
  residualq->add_task(task138);

  auto tensor139 = vector<shared_ptr<Tensor>>{I114, Gamma77_(), t2};
  auto task139 = make_shared<Task139>(tensor139, pindex);
  task126->add_dep(task139);
  task139->add_dep(task40);
  residualq->add_task(task139);

  vector<IndexRange> I132_index = {closed_, virt_, closed_, virt_};
  auto I132 = make_shared<Tensor>(I132_index);
  auto tensor140 = vector<shared_ptr<Tensor>>{r, I132};
  auto task140 = make_shared<Task140>(tensor140, pindex);
  task140->add_dep(task40);
  residualq->add_task(task140);

  vector<IndexRange> I133_index = {closed_, active_};
  auto I133 = make_shared<Tensor>(I133_index);
  auto tensor141 = vector<shared_ptr<Tensor>>{I132, t2, I133};
  auto task141 = make_shared<Task141>(tensor141, pindex);
  task140->add_dep(task141);
  task141->add_dep(task40);
  residualq->add_task(task141);

  auto tensor142 = vector<shared_ptr<Tensor>>{I133, Gamma21_(), f1_};
  auto task142 = make_shared<Task142>(tensor142, pindex);
  task141->add_dep(task142);
  task142->add_dep(task40);
  residualq->add_task(task142);

  vector<IndexRange> I136_index = {closed_, active_};
  auto I136 = make_shared<Tensor>(I136_index);
  auto tensor143 = vector<shared_ptr<Tensor>>{I132, t2, I136};
  auto task143 = make_shared<Task143>(tensor143, pindex);
  task140->add_dep(task143);
  task143->add_dep(task40);
  residualq->add_task(task143);

  auto tensor144 = vector<shared_ptr<Tensor>>{I136, Gamma21_(), f1_};
  auto task144 = make_shared<Task144>(tensor144, pindex);
  task143->add_dep(task144);
  task144->add_dep(task40);
  residualq->add_task(task144);

  vector<IndexRange> I139_index = {virt_, closed_};
  auto I139 = make_shared<Tensor>(I139_index);
  auto tensor145 = vector<shared_ptr<Tensor>>{I132, f1_, I139};
  auto task145 = make_shared<Task145>(tensor145, pindex);
  task140->add_dep(task145);
  task145->add_dep(task40);
  residualq->add_task(task145);

  vector<IndexRange> I140_index = {active_, virt_, closed_, active_};
  auto I140 = make_shared<Tensor>(I140_index);
  auto tensor146 = vector<shared_ptr<Tensor>>{I139, Gamma21_(), I140};
  auto task146 = make_shared<Task146>(tensor146, pindex);
  task145->add_dep(task146);
  task146->add_dep(task40);
  residualq->add_task(task146);

  auto tensor147 = vector<shared_ptr<Tensor>>{I140, t2};
  auto task147 = make_shared<Task147>(tensor147, pindex);
  task146->add_dep(task147);
  task147->add_dep(task40);
  residualq->add_task(task147);

  vector<IndexRange> I142_index = {virt_, closed_};
  auto I142 = make_shared<Tensor>(I142_index);
  auto tensor148 = vector<shared_ptr<Tensor>>{I132, f1_, I142};
  auto task148 = make_shared<Task148>(tensor148, pindex);
  task140->add_dep(task148);
  task148->add_dep(task40);
  residualq->add_task(task148);

  vector<IndexRange> I143_index = {active_, virt_, closed_, active_};
  auto I143 = make_shared<Tensor>(I143_index);
  auto tensor149 = vector<shared_ptr<Tensor>>{I142, Gamma21_(), I143};
  auto task149 = make_shared<Task149>(tensor149, pindex);
  task148->add_dep(task149);
  task149->add_dep(task40);
  residualq->add_task(task149);

  auto tensor150 = vector<shared_ptr<Tensor>>{I143, t2};
  auto task150 = make_shared<Task150>(tensor150, pindex);
  task149->add_dep(task150);
  task150->add_dep(task40);
  residualq->add_task(task150);

  vector<IndexRange> I151_index = {virt_, active_};
  auto I151 = make_shared<Tensor>(I151_index);
  auto tensor151 = vector<shared_ptr<Tensor>>{I132, t2, I151};
  auto task151 = make_shared<Task151>(tensor151, pindex);
  task140->add_dep(task151);
  task151->add_dep(task40);
  residualq->add_task(task151);

  auto tensor152 = vector<shared_ptr<Tensor>>{I151, Gamma13_(), f1_};
  auto task152 = make_shared<Task152>(tensor152, pindex);
  task151->add_dep(task152);
  task152->add_dep(task40);
  residualq->add_task(task152);

  vector<IndexRange> I154_index = {virt_, active_};
  auto I154 = make_shared<Tensor>(I154_index);
  auto tensor153 = vector<shared_ptr<Tensor>>{I132, t2, I154};
  auto task153 = make_shared<Task153>(tensor153, pindex);
  task140->add_dep(task153);
  task153->add_dep(task40);
  residualq->add_task(task153);

  auto tensor154 = vector<shared_ptr<Tensor>>{I154, Gamma13_(), f1_};
  auto task154 = make_shared<Task154>(tensor154, pindex);
  task153->add_dep(task154);
  task154->add_dep(task40);
  residualq->add_task(task154);

  shared_ptr<Tensor> I157;
  if (diagonal) {
    vector<IndexRange> I157_index = {closed_, virt_, closed_, virt_};
    I157 = make_shared<Tensor>(I157_index);
  }
  shared_ptr<Task155> task155;
  if (diagonal) {
    auto tensor155 = vector<shared_ptr<Tensor>>{I132, f1_, I157};
    task155 = make_shared<Task155>(tensor155, pindex);
    task140->add_dep(task155);
    task155->add_dep(task40);
    residualq->add_task(task155);
  }

  shared_ptr<Task156> task156;
  if (diagonal) {
    auto tensor156 = vector<shared_ptr<Tensor>>{I157, t2};
    task156 = make_shared<Task156>(tensor156, pindex);
    task155->add_dep(task156);
    task156->add_dep(task40);
    residualq->add_task(task156);
  }

  vector<IndexRange> I160_index = {active_, closed_, virt_, virt_};
  auto I160 = make_shared<Tensor>(I160_index);
  auto tensor157 = vector<shared_ptr<Tensor>>{r, I160};
  auto task157 = make_shared<Task157>(tensor157, pindex);
  task157->add_dep(task40);
  residualq->add_task(task157);

  vector<IndexRange> I161_index = {closed_, active_};
  auto I161 = make_shared<Tensor>(I161_index);
  auto tensor158 = vector<shared_ptr<Tensor>>{I160, t2, I161};
  auto task158 = make_shared<Task158>(tensor158, pindex);
  task157->add_dep(task158);
  task158->add_dep(task40);
  residualq->add_task(task158);

  auto tensor159 = vector<shared_ptr<Tensor>>{I161, Gamma21_(), f1_};
  auto task159 = make_shared<Task159>(tensor159, pindex);
  task158->add_dep(task159);
  task159->add_dep(task40);
  residualq->add_task(task159);

  vector<IndexRange> I164_index = {closed_, active_};
  auto I164 = make_shared<Tensor>(I164_index);
  auto tensor160 = vector<shared_ptr<Tensor>>{I160, t2, I164};
  auto task160 = make_shared<Task160>(tensor160, pindex);
  task157->add_dep(task160);
  task160->add_dep(task40);
  residualq->add_task(task160);

  auto tensor161 = vector<shared_ptr<Tensor>>{I164, Gamma21_(), f1_};
  auto task161 = make_shared<Task161>(tensor161, pindex);
  task160->add_dep(task161);
  task161->add_dep(task40);
  residualq->add_task(task161);

  vector<IndexRange> I167_index = {virt_, virt_, active_, active_};
  auto I167 = make_shared<Tensor>(I167_index);
  auto tensor162 = vector<shared_ptr<Tensor>>{I160, f1_, I167};
  auto task162 = make_shared<Task162>(tensor162, pindex);
  task157->add_dep(task162);
  task162->add_dep(task40);
  residualq->add_task(task162);

  auto tensor163 = vector<shared_ptr<Tensor>>{I167, Gamma40_(), t2};
  auto task163 = make_shared<Task163>(tensor163, pindex);
  task162->add_dep(task163);
  task163->add_dep(task40);
  residualq->add_task(task163);

  vector<IndexRange> I170_index = {virt_, active_};
  auto I170 = make_shared<Tensor>(I170_index);
  auto tensor164 = vector<shared_ptr<Tensor>>{I160, f1_, I170};
  auto task164 = make_shared<Task164>(tensor164, pindex);
  task157->add_dep(task164);
  task164->add_dep(task40);
  residualq->add_task(task164);

  auto tensor165 = vector<shared_ptr<Tensor>>{I170, Gamma40_(), t2};
  auto task165 = make_shared<Task165>(tensor165, pindex);
  task164->add_dep(task165);
  task165->add_dep(task40);
  residualq->add_task(task165);

  vector<IndexRange> I173_index = {virt_, active_};
  auto I173 = make_shared<Tensor>(I173_index);
  auto tensor166 = vector<shared_ptr<Tensor>>{I160, f1_, I173};
  auto task166 = make_shared<Task166>(tensor166, pindex);
  task157->add_dep(task166);
  task166->add_dep(task40);
  residualq->add_task(task166);

  auto tensor167 = vector<shared_ptr<Tensor>>{I173, Gamma40_(), t2};
  auto task167 = make_shared<Task167>(tensor167, pindex);
  task166->add_dep(task167);
  task167->add_dep(task40);
  residualq->add_task(task167);

  vector<IndexRange> I176_index = {virt_, closed_, active_, active_};
  auto I176 = make_shared<Tensor>(I176_index);
  auto tensor168 = vector<shared_ptr<Tensor>>{I160, f1_, I176};
  auto task168 = make_shared<Task168>(tensor168, pindex);
  task157->add_dep(task168);
  task168->add_dep(task40);
  residualq->add_task(task168);

  vector<IndexRange> I177_index = {active_, virt_, closed_, active_};
  auto I177 = make_shared<Tensor>(I177_index);
  auto tensor169 = vector<shared_ptr<Tensor>>{I176, Gamma23_(), I177};
  auto task169 = make_shared<Task169>(tensor169, pindex);
  task168->add_dep(task169);
  task169->add_dep(task40);
  residualq->add_task(task169);

  auto tensor170 = vector<shared_ptr<Tensor>>{I177, t2};
  auto task170 = make_shared<Task170>(tensor170, pindex);
  task169->add_dep(task170);
  task170->add_dep(task40);
  residualq->add_task(task170);

  vector<IndexRange> I179_index = {virt_, closed_, active_, active_};
  auto I179 = make_shared<Tensor>(I179_index);
  auto tensor171 = vector<shared_ptr<Tensor>>{I160, f1_, I179};
  auto task171 = make_shared<Task171>(tensor171, pindex);
  task157->add_dep(task171);
  task171->add_dep(task40);
  residualq->add_task(task171);

  auto tensor172 = vector<shared_ptr<Tensor>>{I179, Gamma24_(), t2};
  auto task172 = make_shared<Task172>(tensor172, pindex);
  task171->add_dep(task172);
  task172->add_dep(task40);
  residualq->add_task(task172);

  auto tensor173 = vector<shared_ptr<Tensor>>{I179, Gamma23_(), t2};
  auto task173 = make_shared<Task173>(tensor173, pindex);
  task171->add_dep(task173);
  task173->add_dep(task40);
  residualq->add_task(task173);

  vector<IndexRange> I188_index = {virt_, closed_, virt_, active_};
  auto I188 = make_shared<Tensor>(I188_index);
  auto tensor174 = vector<shared_ptr<Tensor>>{I160, f1_, I188};
  auto task174 = make_shared<Task174>(tensor174, pindex);
  task157->add_dep(task174);
  task174->add_dep(task40);
  residualq->add_task(task174);

  vector<IndexRange> I189_index = {active_, virt_, closed_, virt_};
  auto I189 = make_shared<Tensor>(I189_index);
  auto tensor175 = vector<shared_ptr<Tensor>>{I188, Gamma21_(), I189};
  auto task175 = make_shared<Task175>(tensor175, pindex);
  task174->add_dep(task175);
  task175->add_dep(task40);
  residualq->add_task(task175);

  auto tensor176 = vector<shared_ptr<Tensor>>{I189, t2};
  auto task176 = make_shared<Task176>(tensor176, pindex);
  task175->add_dep(task176);
  task176->add_dep(task40);
  residualq->add_task(task176);

  vector<IndexRange> I194_index = {virt_, closed_, virt_, active_};
  auto I194 = make_shared<Tensor>(I194_index);
  auto tensor177 = vector<shared_ptr<Tensor>>{I160, f1_, I194};
  auto task177 = make_shared<Task177>(tensor177, pindex);
  task157->add_dep(task177);
  task177->add_dep(task40);
  residualq->add_task(task177);

  vector<IndexRange> I195_index = {active_, virt_, closed_, virt_};
  auto I195 = make_shared<Tensor>(I195_index);
  auto tensor178 = vector<shared_ptr<Tensor>>{I194, Gamma21_(), I195};
  auto task178 = make_shared<Task178>(tensor178, pindex);
  task177->add_dep(task178);
  task178->add_dep(task40);
  residualq->add_task(task178);

  auto tensor179 = vector<shared_ptr<Tensor>>{I195, t2};
  auto task179 = make_shared<Task179>(tensor179, pindex);
  task178->add_dep(task179);
  task179->add_dep(task40);
  residualq->add_task(task179);

  vector<IndexRange> I237_index = {active_, virt_, closed_, virt_};
  auto I237 = make_shared<Tensor>(I237_index);
  auto tensor180 = vector<shared_ptr<Tensor>>{I160, Gamma82_(), I237};
  auto task180 = make_shared<Task180>(tensor180, pindex);
  task157->add_dep(task180);
  task180->add_dep(task40);
  residualq->add_task(task180);

  auto tensor181 = vector<shared_ptr<Tensor>>{I237, t2};
  auto task181 = make_shared<Task181>(tensor181, pindex);
  task180->add_dep(task181);
  task181->add_dep(task40);
  residualq->add_task(task181);

  vector<IndexRange> I241_index = {active_, virt_, closed_, virt_};
  auto I241 = make_shared<Tensor>(I241_index);
  auto tensor182 = vector<shared_ptr<Tensor>>{I160, Gamma84_(), I241};
  auto task182 = make_shared<Task182>(tensor182, pindex);
  task157->add_dep(task182);
  task182->add_dep(task40);
  residualq->add_task(task182);

  auto tensor183 = vector<shared_ptr<Tensor>>{I241, t2};
  auto task183 = make_shared<Task183>(tensor183, pindex);
  task182->add_dep(task183);
  task183->add_dep(task40);
  residualq->add_task(task183);

  vector<IndexRange> I199_index = {virt_, virt_, active_, active_};
  auto I199 = make_shared<Tensor>(I199_index);
  auto tensor184 = vector<shared_ptr<Tensor>>{r, I199};
  auto task184 = make_shared<Task184>(tensor184, pindex);
  task184->add_dep(task40);
  residualq->add_task(task184);

  vector<IndexRange> I200_index = {active_, virt_, virt_, active_};
  auto I200 = make_shared<Tensor>(I200_index);
  auto tensor185 = vector<shared_ptr<Tensor>>{I199, Gamma40_(), I200};
  auto task185 = make_shared<Task185>(tensor185, pindex);
  task184->add_dep(task185);
  task185->add_dep(task40);
  residualq->add_task(task185);

  auto tensor186 = vector<shared_ptr<Tensor>>{I200, f1_, t2};
  auto task186 = make_shared<Task186>(tensor186, pindex);
  task185->add_dep(task186);
  task186->add_dep(task40);
  residualq->add_task(task186);

  vector<IndexRange> I203_index = {virt_, active_, active_, active_};
  auto I203 = make_shared<Tensor>(I203_index);
  auto tensor187 = vector<shared_ptr<Tensor>>{I199, f1_, I203};
  auto task187 = make_shared<Task187>(tensor187, pindex);
  task184->add_dep(task187);
  task187->add_dep(task40);
  residualq->add_task(task187);

  auto tensor188 = vector<shared_ptr<Tensor>>{I203, Gamma42_(), t2};
  auto task188 = make_shared<Task188>(tensor188, pindex);
  task187->add_dep(task188);
  task188->add_dep(task40);
  residualq->add_task(task188);

  vector<IndexRange> I206_index = {virt_, virt_, active_, active_};
  auto I206 = make_shared<Tensor>(I206_index);
  auto tensor189 = vector<shared_ptr<Tensor>>{I199, f1_, I206};
  auto task189 = make_shared<Task189>(tensor189, pindex);
  task184->add_dep(task189);
  task189->add_dep(task40);
  residualq->add_task(task189);

  auto tensor190 = vector<shared_ptr<Tensor>>{I206, Gamma40_(), t2};
  auto task190 = make_shared<Task190>(tensor190, pindex);
  task189->add_dep(task190);
  task190->add_dep(task40);
  residualq->add_task(task190);

  vector<IndexRange> I208_index = {closed_, closed_, active_, active_};
  auto I208 = make_shared<Tensor>(I208_index);
  auto tensor191 = vector<shared_ptr<Tensor>>{r, I208};
  auto task191 = make_shared<Task191>(tensor191, pindex);
  task191->add_dep(task40);
  residualq->add_task(task191);

  auto tensor192 = vector<shared_ptr<Tensor>>{I208, Gamma68_(), t2};
  auto task192 = make_shared<Task192>(tensor192, pindex);
  task191->add_dep(task192);
  task192->add_dep(task40);
  residualq->add_task(task192);

  auto tensor193 = vector<shared_ptr<Tensor>>{I208, Gamma69_(), t2};
  auto task193 = make_shared<Task193>(tensor193, pindex);
  task191->add_dep(task193);
  task193->add_dep(task40);
  residualq->add_task(task193);

  vector<IndexRange> I228_index = {closed_, virt_, closed_, virt_};
  auto I228 = make_shared<Tensor>(I228_index);
  auto tensor194 = vector<shared_ptr<Tensor>>{r, I228};
  auto task194 = make_shared<Task194>(tensor194, pindex);
  task194->add_dep(task40);
  residualq->add_task(task194);

  vector<IndexRange> I229_index = {closed_, virt_, closed_, virt_};
  auto I229 = make_shared<Tensor>(I229_index);
  auto tensor195 = vector<shared_ptr<Tensor>>{I228, Gamma78_(), I229};
  auto task195 = make_shared<Task195>(tensor195, pindex);
  task194->add_dep(task195);
  task195->add_dep(task40);
  residualq->add_task(task195);

  auto tensor196 = vector<shared_ptr<Tensor>>{I229, t2};
  auto task196 = make_shared<Task196>(tensor196, pindex);
  task195->add_dep(task196);
  task196->add_dep(task40);
  residualq->add_task(task196);

  vector<IndexRange> I233_index = {closed_, virt_, closed_, virt_};
  auto I233 = make_shared<Tensor>(I233_index);
  auto tensor197 = vector<shared_ptr<Tensor>>{I228, Gamma80_(), I233};
  auto task197 = make_shared<Task197>(tensor197, pindex);
  task194->add_dep(task197);
  task197->add_dep(task40);
  residualq->add_task(task197);

  auto tensor198 = vector<shared_ptr<Tensor>>{I233, t2};
  auto task198 = make_shared<Task198>(tensor198, pindex);
  task197->add_dep(task198);
  task198->add_dep(task40);
  residualq->add_task(task198);

  vector<IndexRange> I244_index = {virt_, virt_, active_, active_};
  auto I244 = make_shared<Tensor>(I244_index);
  auto tensor199 = vector<shared_ptr<Tensor>>{r, I244};
  auto task199 = make_shared<Task199>(tensor199, pindex);
  task199->add_dep(task40);
  residualq->add_task(task199);

  auto tensor200 = vector<shared_ptr<Tensor>>{I244, Gamma86_(), t2};
  auto task200 = make_shared<Task200>(tensor200, pindex);
  task199->add_dep(task200);
  task200->add_dep(task40);
  residualq->add_task(task200);

  auto tensor201 = vector<shared_ptr<Tensor>>{I244, Gamma87_(), t2};
  auto task201 = make_shared<Task201>(tensor201, pindex);
  task199->add_dep(task201);
  task201->add_dep(task40);
  residualq->add_task(task201);

  return residualq;
}


#endif
