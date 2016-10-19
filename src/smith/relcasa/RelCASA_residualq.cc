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

  vector<IndexRange> I1_index = {closed_, closed_, active_, active_};
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

  vector<IndexRange> I4_index = {closed_, active_, active_, active_};
  auto I4 = make_shared<Tensor>(I4_index);
  auto tensor44 = vector<shared_ptr<Tensor>>{I0, f1_, I4};
  auto task44 = make_shared<Task44>(tensor44, pindex);
  task41->add_dep(task44);
  task44->add_dep(task40);
  residualq->add_task(task44);

  auto tensor45 = vector<shared_ptr<Tensor>>{I4, Gamma1_(), t2};
  auto task45 = make_shared<Task45>(tensor45, pindex);
  task44->add_dep(task45);
  task45->add_dep(task40);
  residualq->add_task(task45);

  vector<IndexRange> I7_index = {closed_, closed_, active_, active_};
  auto I7 = make_shared<Tensor>(I7_index);
  auto tensor46 = vector<shared_ptr<Tensor>>{I0, Gamma2_(), I7};
  auto task46 = make_shared<Task46>(tensor46, pindex);
  task41->add_dep(task46);
  task46->add_dep(task40);
  residualq->add_task(task46);

  auto tensor47 = vector<shared_ptr<Tensor>>{I7, f1_, t2};
  auto task47 = make_shared<Task47>(tensor47, pindex);
  task46->add_dep(task47);
  task47->add_dep(task40);
  residualq->add_task(task47);

  vector<IndexRange> I9_index = {active_, active_, active_, closed_};
  auto I9 = make_shared<Tensor>(I9_index);
  auto tensor48 = vector<shared_ptr<Tensor>>{r, I9};
  auto task48 = make_shared<Task48>(tensor48, pindex);
  task48->add_dep(task40);
  residualq->add_task(task48);

  vector<IndexRange> I10_index = {closed_, active_, active_, active_};
  auto I10 = make_shared<Tensor>(I10_index);
  auto tensor49 = vector<shared_ptr<Tensor>>{I9, f1_, I10};
  auto task49 = make_shared<Task49>(tensor49, pindex);
  task48->add_dep(task49);
  task49->add_dep(task40);
  residualq->add_task(task49);

  auto tensor50 = vector<shared_ptr<Tensor>>{I10, Gamma3_(), t2};
  auto task50 = make_shared<Task50>(tensor50, pindex);
  task49->add_dep(task50);
  task50->add_dep(task40);
  residualq->add_task(task50);

  vector<IndexRange> I13_index = {closed_, active_, active_, active_};
  auto I13 = make_shared<Tensor>(I13_index);
  auto tensor51 = vector<shared_ptr<Tensor>>{I9, Gamma4_(), I13};
  auto task51 = make_shared<Task51>(tensor51, pindex);
  task48->add_dep(task51);
  task51->add_dep(task40);
  residualq->add_task(task51);

  auto tensor52 = vector<shared_ptr<Tensor>>{I13, f1_, t2};
  auto task52 = make_shared<Task52>(tensor52, pindex);
  task51->add_dep(task52);
  task52->add_dep(task40);
  residualq->add_task(task52);

  vector<IndexRange> I16_index = {closed_, active_};
  auto I16 = make_shared<Tensor>(I16_index);
  auto tensor53 = vector<shared_ptr<Tensor>>{I9, Gamma5_(), I16};
  auto task53 = make_shared<Task53>(tensor53, pindex);
  task48->add_dep(task53);
  task53->add_dep(task40);
  residualq->add_task(task53);

  vector<IndexRange> I17_index = {closed_, virt_, closed_, active_};
  auto I17 = make_shared<Tensor>(I17_index);
  auto tensor54 = vector<shared_ptr<Tensor>>{I16, f1_, I17};
  auto task54 = make_shared<Task54>(tensor54, pindex);
  task53->add_dep(task54);
  task54->add_dep(task40);
  residualq->add_task(task54);

  auto tensor55 = vector<shared_ptr<Tensor>>{I17, t2};
  auto task55 = make_shared<Task55>(tensor55, pindex);
  task54->add_dep(task55);
  task55->add_dep(task40);
  residualq->add_task(task55);

  vector<IndexRange> I22_index = {active_, closed_, active_, active_};
  auto I22 = make_shared<Tensor>(I22_index);
  auto tensor56 = vector<shared_ptr<Tensor>>{I9, Gamma7_(), I22};
  auto task56 = make_shared<Task56>(tensor56, pindex);
  task48->add_dep(task56);
  task56->add_dep(task40);
  residualq->add_task(task56);

  auto tensor57 = vector<shared_ptr<Tensor>>{I22, f1_, t2};
  auto task57 = make_shared<Task57>(tensor57, pindex);
  task56->add_dep(task57);
  task57->add_dep(task40);
  residualq->add_task(task57);

  vector<IndexRange> I25_index = {closed_, active_, active_, active_};
  auto I25 = make_shared<Tensor>(I25_index);
  auto tensor58 = vector<shared_ptr<Tensor>>{I9, Gamma3_(), I25};
  auto task58 = make_shared<Task58>(tensor58, pindex);
  task48->add_dep(task58);
  task58->add_dep(task40);
  residualq->add_task(task58);

  auto tensor59 = vector<shared_ptr<Tensor>>{I25, t2};
  auto task59 = make_shared<Task59>(tensor59, pindex, this->e0_);
  task58->add_dep(task59);
  task59->add_dep(task40);
  residualq->add_task(task59);

  auto tensor60 = vector<shared_ptr<Tensor>>{I25, f1_, t2};
  auto task60 = make_shared<Task60>(tensor60, pindex);
  task58->add_dep(task60);
  task60->add_dep(task40);
  residualq->add_task(task60);

  auto tensor61 = vector<shared_ptr<Tensor>>{I9, Gamma94_(), t2};
  auto task61 = make_shared<Task61>(tensor61, pindex);
  task48->add_dep(task61);
  task61->add_dep(task40);
  residualq->add_task(task61);

  auto tensor62 = vector<shared_ptr<Tensor>>{I9, Gamma95_(), t2};
  auto task62 = make_shared<Task62>(tensor62, pindex);
  task48->add_dep(task62);
  task62->add_dep(task40);
  residualq->add_task(task62);

  vector<IndexRange> I27_index = {virt_, closed_, active_, closed_};
  auto I27 = make_shared<Tensor>(I27_index);
  auto tensor63 = vector<shared_ptr<Tensor>>{r, I27};
  auto task63 = make_shared<Task63>(tensor63, pindex);
  task63->add_dep(task40);
  residualq->add_task(task63);

  vector<IndexRange> I28_index = {closed_, virt_, closed_, active_};
  auto I28 = make_shared<Tensor>(I28_index);
  auto tensor64 = vector<shared_ptr<Tensor>>{I27, f1_, I28};
  auto task64 = make_shared<Task64>(tensor64, pindex);
  task63->add_dep(task64);
  task64->add_dep(task40);
  residualq->add_task(task64);

  vector<IndexRange> I29_index = {closed_, virt_, closed_, active_};
  auto I29 = make_shared<Tensor>(I29_index);
  auto tensor65 = vector<shared_ptr<Tensor>>{I28, Gamma9_(), I29};
  auto task65 = make_shared<Task65>(tensor65, pindex);
  task64->add_dep(task65);
  task65->add_dep(task40);
  residualq->add_task(task65);

  auto tensor66 = vector<shared_ptr<Tensor>>{I29, t2};
  auto task66 = make_shared<Task66>(tensor66, pindex);
  task65->add_dep(task66);
  task66->add_dep(task40);
  residualq->add_task(task66);

  vector<IndexRange> I34_index = {closed_, virt_, closed_, active_};
  auto I34 = make_shared<Tensor>(I34_index);
  auto tensor67 = vector<shared_ptr<Tensor>>{I27, f1_, I34};
  auto task67 = make_shared<Task67>(tensor67, pindex);
  task63->add_dep(task67);
  task67->add_dep(task40);
  residualq->add_task(task67);

  vector<IndexRange> I35_index = {closed_, virt_, closed_, active_};
  auto I35 = make_shared<Tensor>(I35_index);
  auto tensor68 = vector<shared_ptr<Tensor>>{I34, Gamma9_(), I35};
  auto task68 = make_shared<Task68>(tensor68, pindex);
  task67->add_dep(task68);
  task68->add_dep(task40);
  residualq->add_task(task68);

  auto tensor69 = vector<shared_ptr<Tensor>>{I35, t2};
  auto task69 = make_shared<Task69>(tensor69, pindex);
  task68->add_dep(task69);
  task69->add_dep(task40);
  residualq->add_task(task69);

  vector<IndexRange> I40_index = {virt_, closed_, active_, active_};
  auto I40 = make_shared<Tensor>(I40_index);
  auto tensor70 = vector<shared_ptr<Tensor>>{I27, f1_, I40};
  auto task70 = make_shared<Task70>(tensor70, pindex);
  task63->add_dep(task70);
  task70->add_dep(task40);
  residualq->add_task(task70);

  auto tensor71 = vector<shared_ptr<Tensor>>{I40, Gamma13_(), t2};
  auto task71 = make_shared<Task71>(tensor71, pindex);
  task70->add_dep(task71);
  task71->add_dep(task40);
  residualq->add_task(task71);

  auto tensor72 = vector<shared_ptr<Tensor>>{I40, Gamma15_(), t2};
  auto task72 = make_shared<Task72>(tensor72, pindex);
  task70->add_dep(task72);
  task72->add_dep(task40);
  residualq->add_task(task72);

  vector<IndexRange> I43_index = {virt_, closed_, active_, active_};
  auto I43 = make_shared<Tensor>(I43_index);
  auto tensor73 = vector<shared_ptr<Tensor>>{I27, f1_, I43};
  auto task73 = make_shared<Task73>(tensor73, pindex);
  task63->add_dep(task73);
  task73->add_dep(task40);
  residualq->add_task(task73);

  vector<IndexRange> I44_index = {active_, virt_, closed_, active_};
  auto I44 = make_shared<Tensor>(I44_index);
  auto tensor74 = vector<shared_ptr<Tensor>>{I43, Gamma15_(), I44};
  auto task74 = make_shared<Task74>(tensor74, pindex);
  task73->add_dep(task74);
  task74->add_dep(task40);
  residualq->add_task(task74);

  auto tensor75 = vector<shared_ptr<Tensor>>{I44, t2};
  auto task75 = make_shared<Task75>(tensor75, pindex);
  task74->add_dep(task75);
  task75->add_dep(task40);
  residualq->add_task(task75);

  vector<IndexRange> I52_index = {closed_, active_};
  auto I52 = make_shared<Tensor>(I52_index);
  auto tensor76 = vector<shared_ptr<Tensor>>{I27, f1_, I52};
  auto task76 = make_shared<Task76>(tensor76, pindex);
  task63->add_dep(task76);
  task76->add_dep(task40);
  residualq->add_task(task76);

  auto tensor77 = vector<shared_ptr<Tensor>>{I52, Gamma15_(), t2};
  auto task77 = make_shared<Task77>(tensor77, pindex);
  task76->add_dep(task77);
  task77->add_dep(task40);
  residualq->add_task(task77);

  vector<IndexRange> I55_index = {closed_, active_};
  auto I55 = make_shared<Tensor>(I55_index);
  auto tensor78 = vector<shared_ptr<Tensor>>{I27, f1_, I55};
  auto task78 = make_shared<Task78>(tensor78, pindex);
  task63->add_dep(task78);
  task78->add_dep(task40);
  residualq->add_task(task78);

  auto tensor79 = vector<shared_ptr<Tensor>>{I55, Gamma15_(), t2};
  auto task79 = make_shared<Task79>(tensor79, pindex);
  task78->add_dep(task79);
  task79->add_dep(task40);
  residualq->add_task(task79);

  vector<IndexRange> I58_index = {virt_, active_};
  auto I58 = make_shared<Tensor>(I58_index);
  auto tensor80 = vector<shared_ptr<Tensor>>{I27, t2, I58};
  auto task80 = make_shared<Task80>(tensor80, pindex);
  task63->add_dep(task80);
  task80->add_dep(task40);
  residualq->add_task(task80);

  auto tensor81 = vector<shared_ptr<Tensor>>{I58, Gamma9_(), f1_};
  auto task81 = make_shared<Task81>(tensor81, pindex);
  task80->add_dep(task81);
  task81->add_dep(task40);
  residualq->add_task(task81);

  vector<IndexRange> I61_index = {virt_, active_};
  auto I61 = make_shared<Tensor>(I61_index);
  auto tensor82 = vector<shared_ptr<Tensor>>{I27, t2, I61};
  auto task82 = make_shared<Task82>(tensor82, pindex);
  task63->add_dep(task82);
  task82->add_dep(task40);
  residualq->add_task(task82);

  auto tensor83 = vector<shared_ptr<Tensor>>{I61, Gamma9_(), f1_};
  auto task83 = make_shared<Task83>(tensor83, pindex);
  task82->add_dep(task83);
  task83->add_dep(task40);
  residualq->add_task(task83);

  vector<IndexRange> I64_index = {closed_, closed_, active_, active_};
  auto I64 = make_shared<Tensor>(I64_index);
  auto tensor84 = vector<shared_ptr<Tensor>>{I27, f1_, I64};
  auto task84 = make_shared<Task84>(tensor84, pindex);
  task63->add_dep(task84);
  task84->add_dep(task40);
  residualq->add_task(task84);

  auto tensor85 = vector<shared_ptr<Tensor>>{I64, Gamma2_(), t2};
  auto task85 = make_shared<Task85>(tensor85, pindex);
  task84->add_dep(task85);
  task85->add_dep(task40);
  residualq->add_task(task85);

  vector<IndexRange> I67_index = {closed_, virt_, closed_, active_};
  auto I67 = make_shared<Tensor>(I67_index);
  auto tensor86 = vector<shared_ptr<Tensor>>{I27, f1_, I67};
  auto task86 = make_shared<Task86>(tensor86, pindex);
  task63->add_dep(task86);
  task86->add_dep(task40);
  residualq->add_task(task86);

  vector<IndexRange> I68_index = {closed_, virt_, closed_, active_};
  auto I68 = make_shared<Tensor>(I68_index);
  auto tensor87 = vector<shared_ptr<Tensor>>{I67, Gamma9_(), I68};
  auto task87 = make_shared<Task87>(tensor87, pindex);
  task86->add_dep(task87);
  task87->add_dep(task40);
  residualq->add_task(task87);

  auto tensor88 = vector<shared_ptr<Tensor>>{I68, t2};
  auto task88 = make_shared<Task88>(tensor88, pindex);
  task87->add_dep(task88);
  task88->add_dep(task40);
  residualq->add_task(task88);

  vector<IndexRange> I253_index = {closed_, virt_, closed_, active_};
  auto I253 = make_shared<Tensor>(I253_index);
  auto tensor89 = vector<shared_ptr<Tensor>>{I27, Gamma9_(), I253};
  auto task89 = make_shared<Task89>(tensor89, pindex);
  task63->add_dep(task89);
  task89->add_dep(task40);
  residualq->add_task(task89);

  auto tensor90 = vector<shared_ptr<Tensor>>{I253, t2};
  auto task90 = make_shared<Task90>(tensor90, pindex, this->e0_);
  task89->add_dep(task90);
  task90->add_dep(task40);
  residualq->add_task(task90);

  vector<IndexRange> I283_index = {closed_, virt_, closed_, active_};
  auto I283 = make_shared<Tensor>(I283_index);
  auto tensor91 = vector<shared_ptr<Tensor>>{I27, Gamma96_(), I283};
  auto task91 = make_shared<Task91>(tensor91, pindex);
  task63->add_dep(task91);
  task91->add_dep(task40);
  residualq->add_task(task91);

  auto tensor92 = vector<shared_ptr<Tensor>>{I283, t2};
  auto task92 = make_shared<Task92>(tensor92, pindex);
  task91->add_dep(task92);
  task92->add_dep(task40);
  residualq->add_task(task92);

  vector<IndexRange> I287_index = {closed_, virt_, closed_, active_};
  auto I287 = make_shared<Tensor>(I287_index);
  auto tensor93 = vector<shared_ptr<Tensor>>{I27, Gamma98_(), I287};
  auto task93 = make_shared<Task93>(tensor93, pindex);
  task63->add_dep(task93);
  task93->add_dep(task40);
  residualq->add_task(task93);

  auto tensor94 = vector<shared_ptr<Tensor>>{I287, t2};
  auto task94 = make_shared<Task94>(tensor94, pindex);
  task93->add_dep(task94);
  task94->add_dep(task40);
  residualq->add_task(task94);

  vector<IndexRange> I72_index = {virt_, active_, active_, closed_};
  auto I72 = make_shared<Tensor>(I72_index);
  auto tensor95 = vector<shared_ptr<Tensor>>{r, I72};
  auto task95 = make_shared<Task95>(tensor95, pindex);
  task95->add_dep(task40);
  residualq->add_task(task95);

  vector<IndexRange> I73_index = {virt_, closed_, active_, active_};
  auto I73 = make_shared<Tensor>(I73_index);
  auto tensor96 = vector<shared_ptr<Tensor>>{I72, f1_, I73};
  auto task96 = make_shared<Task96>(tensor96, pindex);
  task95->add_dep(task96);
  task96->add_dep(task40);
  residualq->add_task(task96);

  auto tensor97 = vector<shared_ptr<Tensor>>{I73, Gamma24_(), t2};
  auto task97 = make_shared<Task97>(tensor97, pindex);
  task96->add_dep(task97);
  task97->add_dep(task40);
  residualq->add_task(task97);

  auto tensor98 = vector<shared_ptr<Tensor>>{I73, Gamma25_(), t2};
  auto task98 = make_shared<Task98>(tensor98, pindex);
  task96->add_dep(task98);
  task98->add_dep(task40);
  residualq->add_task(task98);

  vector<IndexRange> I79_index = {virt_, closed_, active_, active_};
  auto I79 = make_shared<Tensor>(I79_index);
  auto tensor99 = vector<shared_ptr<Tensor>>{I72, Gamma26_(), I79};
  auto task99 = make_shared<Task99>(tensor99, pindex);
  task95->add_dep(task99);
  task99->add_dep(task40);
  residualq->add_task(task99);

  auto tensor100 = vector<shared_ptr<Tensor>>{I79, f1_, t2};
  auto task100 = make_shared<Task100>(tensor100, pindex);
  task99->add_dep(task100);
  task100->add_dep(task40);
  residualq->add_task(task100);

  vector<IndexRange> I82_index = {closed_, virt_, active_, active_};
  auto I82 = make_shared<Tensor>(I82_index);
  auto tensor101 = vector<shared_ptr<Tensor>>{I72, Gamma5_(), I82};
  auto task101 = make_shared<Task101>(tensor101, pindex);
  task95->add_dep(task101);
  task101->add_dep(task40);
  residualq->add_task(task101);

  auto tensor102 = vector<shared_ptr<Tensor>>{I82, f1_, t2};
  auto task102 = make_shared<Task102>(tensor102, pindex);
  task101->add_dep(task102);
  task102->add_dep(task40);
  residualq->add_task(task102);

  vector<IndexRange> I85_index = {virt_, active_, active_, active_};
  auto I85 = make_shared<Tensor>(I85_index);
  auto tensor103 = vector<shared_ptr<Tensor>>{I72, f1_, I85};
  auto task103 = make_shared<Task103>(tensor103, pindex);
  task95->add_dep(task103);
  task103->add_dep(task40);
  residualq->add_task(task103);

  auto tensor104 = vector<shared_ptr<Tensor>>{I85, Gamma28_(), t2};
  auto task104 = make_shared<Task104>(tensor104, pindex);
  task103->add_dep(task104);
  task104->add_dep(task40);
  residualq->add_task(task104);

  vector<IndexRange> I88_index = {closed_, virt_};
  auto I88 = make_shared<Tensor>(I88_index);
  auto tensor105 = vector<shared_ptr<Tensor>>{I72, Gamma29_(), I88};
  auto task105 = make_shared<Task105>(tensor105, pindex);
  task95->add_dep(task105);
  task105->add_dep(task40);
  residualq->add_task(task105);

  vector<IndexRange> I89_index = {closed_, virt_, closed_, virt_};
  auto I89 = make_shared<Tensor>(I89_index);
  auto tensor106 = vector<shared_ptr<Tensor>>{I88, f1_, I89};
  auto task106 = make_shared<Task106>(tensor106, pindex);
  task105->add_dep(task106);
  task106->add_dep(task40);
  residualq->add_task(task106);

  auto tensor107 = vector<shared_ptr<Tensor>>{I89, t2};
  auto task107 = make_shared<Task107>(tensor107, pindex);
  task106->add_dep(task107);
  task107->add_dep(task40);
  residualq->add_task(task107);

  vector<IndexRange> I94_index = {active_, closed_, virt_, active_};
  auto I94 = make_shared<Tensor>(I94_index);
  auto tensor108 = vector<shared_ptr<Tensor>>{I72, Gamma25_(), I94};
  auto task108 = make_shared<Task108>(tensor108, pindex);
  task95->add_dep(task108);
  task108->add_dep(task40);
  residualq->add_task(task108);

  auto tensor109 = vector<shared_ptr<Tensor>>{I94, t2};
  auto task109 = make_shared<Task109>(tensor109, pindex, this->e0_);
  task108->add_dep(task109);
  task109->add_dep(task40);
  residualq->add_task(task109);

  auto tensor110 = vector<shared_ptr<Tensor>>{I94, f1_, t2};
  auto task110 = make_shared<Task110>(tensor110, pindex);
  task108->add_dep(task110);
  task110->add_dep(task40);
  residualq->add_task(task110);

  vector<IndexRange> I97_index = {active_, virt_, closed_, active_};
  auto I97 = make_shared<Tensor>(I97_index);
  auto tensor111 = vector<shared_ptr<Tensor>>{I72, Gamma24_(), I97};
  auto task111 = make_shared<Task111>(tensor111, pindex);
  task95->add_dep(task111);
  task111->add_dep(task40);
  residualq->add_task(task111);

  auto tensor112 = vector<shared_ptr<Tensor>>{I97, t2};
  auto task112 = make_shared<Task112>(tensor112, pindex, this->e0_);
  task111->add_dep(task112);
  task112->add_dep(task40);
  residualq->add_task(task112);

  auto tensor113 = vector<shared_ptr<Tensor>>{I97, f1_, t2};
  auto task113 = make_shared<Task113>(tensor113, pindex);
  task111->add_dep(task113);
  task113->add_dep(task40);
  residualq->add_task(task113);

  vector<IndexRange> I100_index = {closed_, active_, active_, active_};
  auto I100 = make_shared<Tensor>(I100_index);
  auto tensor114 = vector<shared_ptr<Tensor>>{I72, f1_, I100};
  auto task114 = make_shared<Task114>(tensor114, pindex);
  task95->add_dep(task114);
  task114->add_dep(task40);
  residualq->add_task(task114);

  auto tensor115 = vector<shared_ptr<Tensor>>{I100, Gamma33_(), t2};
  auto task115 = make_shared<Task115>(tensor115, pindex);
  task114->add_dep(task115);
  task115->add_dep(task40);
  residualq->add_task(task115);

  vector<IndexRange> I103_index = {virt_, closed_, active_, active_};
  auto I103 = make_shared<Tensor>(I103_index);
  auto tensor116 = vector<shared_ptr<Tensor>>{I72, f1_, I103};
  auto task116 = make_shared<Task116>(tensor116, pindex);
  task95->add_dep(task116);
  task116->add_dep(task40);
  residualq->add_task(task116);

  auto tensor117 = vector<shared_ptr<Tensor>>{I103, Gamma24_(), t2};
  auto task117 = make_shared<Task117>(tensor117, pindex);
  task116->add_dep(task117);
  task117->add_dep(task40);
  residualq->add_task(task117);

  auto tensor118 = vector<shared_ptr<Tensor>>{I103, Gamma25_(), t2};
  auto task118 = make_shared<Task118>(tensor118, pindex);
  task116->add_dep(task118);
  task118->add_dep(task40);
  residualq->add_task(task118);

  auto tensor119 = vector<shared_ptr<Tensor>>{I72, Gamma112_(), t2};
  auto task119 = make_shared<Task119>(tensor119, pindex);
  task95->add_dep(task119);
  task119->add_dep(task40);
  residualq->add_task(task119);

  auto tensor120 = vector<shared_ptr<Tensor>>{I72, Gamma113_(), t2};
  auto task120 = make_shared<Task120>(tensor120, pindex);
  task95->add_dep(task120);
  task120->add_dep(task40);
  residualq->add_task(task120);

  auto tensor121 = vector<shared_ptr<Tensor>>{I72, Gamma116_(), t2};
  auto task121 = make_shared<Task121>(tensor121, pindex);
  task95->add_dep(task121);
  task121->add_dep(task40);
  residualq->add_task(task121);

  auto tensor122 = vector<shared_ptr<Tensor>>{I72, Gamma117_(), t2};
  auto task122 = make_shared<Task122>(tensor122, pindex);
  task95->add_dep(task122);
  task122->add_dep(task40);
  residualq->add_task(task122);

  vector<IndexRange> I108_index = {virt_, active_, active_, closed_};
  auto I108 = make_shared<Tensor>(I108_index);
  auto tensor123 = vector<shared_ptr<Tensor>>{r, I108};
  auto task123 = make_shared<Task123>(tensor123, pindex);
  task123->add_dep(task40);
  residualq->add_task(task123);

  vector<IndexRange> I109_index = {virt_, closed_, active_, active_};
  auto I109 = make_shared<Tensor>(I109_index);
  auto tensor124 = vector<shared_ptr<Tensor>>{I108, f1_, I109};
  auto task124 = make_shared<Task124>(tensor124, pindex);
  task123->add_dep(task124);
  task124->add_dep(task40);
  residualq->add_task(task124);

  vector<IndexRange> I110_index = {active_, virt_, closed_, active_};
  auto I110 = make_shared<Tensor>(I110_index);
  auto tensor125 = vector<shared_ptr<Tensor>>{I109, Gamma25_(), I110};
  auto task125 = make_shared<Task125>(tensor125, pindex);
  task124->add_dep(task125);
  task125->add_dep(task40);
  residualq->add_task(task125);

  auto tensor126 = vector<shared_ptr<Tensor>>{I110, t2};
  auto task126 = make_shared<Task126>(tensor126, pindex);
  task125->add_dep(task126);
  task126->add_dep(task40);
  residualq->add_task(task126);

  vector<IndexRange> I115_index = {virt_, closed_, active_, active_};
  auto I115 = make_shared<Tensor>(I115_index);
  auto tensor127 = vector<shared_ptr<Tensor>>{I108, Gamma5_(), I115};
  auto task127 = make_shared<Task127>(tensor127, pindex);
  task123->add_dep(task127);
  task127->add_dep(task40);
  residualq->add_task(task127);

  vector<IndexRange> I116_index = {closed_, virt_, closed_, active_};
  auto I116 = make_shared<Tensor>(I116_index);
  auto tensor128 = vector<shared_ptr<Tensor>>{I115, f1_, I116};
  auto task128 = make_shared<Task128>(tensor128, pindex);
  task127->add_dep(task128);
  task128->add_dep(task40);
  residualq->add_task(task128);

  auto tensor129 = vector<shared_ptr<Tensor>>{I116, t2};
  auto task129 = make_shared<Task129>(tensor129, pindex);
  task128->add_dep(task129);
  task129->add_dep(task40);
  residualq->add_task(task129);

  vector<IndexRange> I121_index = {virt_, active_, active_, active_};
  auto I121 = make_shared<Tensor>(I121_index);
  auto tensor130 = vector<shared_ptr<Tensor>>{I108, f1_, I121};
  auto task130 = make_shared<Task130>(tensor130, pindex);
  task123->add_dep(task130);
  task130->add_dep(task40);
  residualq->add_task(task130);

  auto tensor131 = vector<shared_ptr<Tensor>>{I121, Gamma40_(), t2};
  auto task131 = make_shared<Task131>(tensor131, pindex);
  task130->add_dep(task131);
  task131->add_dep(task40);
  residualq->add_task(task131);

  vector<IndexRange> I124_index = {closed_, virt_};
  auto I124 = make_shared<Tensor>(I124_index);
  auto tensor132 = vector<shared_ptr<Tensor>>{I108, Gamma29_(), I124};
  auto task132 = make_shared<Task132>(tensor132, pindex);
  task123->add_dep(task132);
  task132->add_dep(task40);
  residualq->add_task(task132);

  vector<IndexRange> I125_index = {closed_, virt_, closed_, virt_};
  auto I125 = make_shared<Tensor>(I125_index);
  auto tensor133 = vector<shared_ptr<Tensor>>{I124, f1_, I125};
  auto task133 = make_shared<Task133>(tensor133, pindex);
  task132->add_dep(task133);
  task133->add_dep(task40);
  residualq->add_task(task133);

  auto tensor134 = vector<shared_ptr<Tensor>>{I125, t2};
  auto task134 = make_shared<Task134>(tensor134, pindex);
  task133->add_dep(task134);
  task134->add_dep(task40);
  residualq->add_task(task134);

  vector<IndexRange> I130_index = {active_, closed_, virt_, active_};
  auto I130 = make_shared<Tensor>(I130_index);
  auto tensor135 = vector<shared_ptr<Tensor>>{I108, Gamma25_(), I130};
  auto task135 = make_shared<Task135>(tensor135, pindex);
  task123->add_dep(task135);
  task135->add_dep(task40);
  residualq->add_task(task135);

  auto tensor136 = vector<shared_ptr<Tensor>>{I130, t2};
  auto task136 = make_shared<Task136>(tensor136, pindex, this->e0_);
  task135->add_dep(task136);
  task136->add_dep(task40);
  residualq->add_task(task136);

  vector<IndexRange> I131_index = {active_, virt_, closed_, virt_};
  auto I131 = make_shared<Tensor>(I131_index);
  auto tensor137 = vector<shared_ptr<Tensor>>{I130, f1_, I131};
  auto task137 = make_shared<Task137>(tensor137, pindex);
  task135->add_dep(task137);
  task137->add_dep(task40);
  residualq->add_task(task137);

  auto tensor138 = vector<shared_ptr<Tensor>>{I131, t2};
  auto task138 = make_shared<Task138>(tensor138, pindex);
  task137->add_dep(task138);
  task138->add_dep(task40);
  residualq->add_task(task138);

  vector<IndexRange> I136_index = {closed_, active_, active_, active_};
  auto I136 = make_shared<Tensor>(I136_index);
  auto tensor139 = vector<shared_ptr<Tensor>>{I108, f1_, I136};
  auto task139 = make_shared<Task139>(tensor139, pindex);
  task123->add_dep(task139);
  task139->add_dep(task40);
  residualq->add_task(task139);

  auto tensor140 = vector<shared_ptr<Tensor>>{I136, Gamma3_(), t2};
  auto task140 = make_shared<Task140>(tensor140, pindex);
  task139->add_dep(task140);
  task140->add_dep(task40);
  residualq->add_task(task140);

  vector<IndexRange> I139_index = {virt_, closed_, active_, active_};
  auto I139 = make_shared<Tensor>(I139_index);
  auto tensor141 = vector<shared_ptr<Tensor>>{I108, f1_, I139};
  auto task141 = make_shared<Task141>(tensor141, pindex);
  task123->add_dep(task141);
  task141->add_dep(task40);
  residualq->add_task(task141);

  vector<IndexRange> I140_index = {active_, virt_, closed_, active_};
  auto I140 = make_shared<Tensor>(I140_index);
  auto tensor142 = vector<shared_ptr<Tensor>>{I139, Gamma25_(), I140};
  auto task142 = make_shared<Task142>(tensor142, pindex);
  task141->add_dep(task142);
  task142->add_dep(task40);
  residualq->add_task(task142);

  auto tensor143 = vector<shared_ptr<Tensor>>{I140, t2};
  auto task143 = make_shared<Task143>(tensor143, pindex);
  task142->add_dep(task143);
  task143->add_dep(task40);
  residualq->add_task(task143);

  vector<IndexRange> I319_index = {active_, virt_, closed_, active_};
  auto I319 = make_shared<Tensor>(I319_index);
  auto tensor144 = vector<shared_ptr<Tensor>>{I108, Gamma113_(), I319};
  auto task144 = make_shared<Task144>(tensor144, pindex);
  task123->add_dep(task144);
  task144->add_dep(task40);
  residualq->add_task(task144);

  auto tensor145 = vector<shared_ptr<Tensor>>{I319, t2};
  auto task145 = make_shared<Task145>(tensor145, pindex);
  task144->add_dep(task145);
  task145->add_dep(task40);
  residualq->add_task(task145);

  vector<IndexRange> I327_index = {active_, virt_, closed_, active_};
  auto I327 = make_shared<Tensor>(I327_index);
  auto tensor146 = vector<shared_ptr<Tensor>>{I108, Gamma117_(), I327};
  auto task146 = make_shared<Task146>(tensor146, pindex);
  task123->add_dep(task146);
  task146->add_dep(task40);
  residualq->add_task(task146);

  auto tensor147 = vector<shared_ptr<Tensor>>{I327, t2};
  auto task147 = make_shared<Task147>(tensor147, pindex);
  task146->add_dep(task147);
  task147->add_dep(task40);
  residualq->add_task(task147);

  vector<IndexRange> I144_index = {virt_, active_, active_, active_};
  auto I144 = make_shared<Tensor>(I144_index);
  auto tensor148 = vector<shared_ptr<Tensor>>{r, I144};
  auto task148 = make_shared<Task148>(tensor148, pindex);
  task148->add_dep(task40);
  residualq->add_task(task148);

  vector<IndexRange> I145_index = {active_, virt_, active_, active_};
  auto I145 = make_shared<Tensor>(I145_index);
  auto tensor149 = vector<shared_ptr<Tensor>>{I144, Gamma48_(), I145};
  auto task149 = make_shared<Task149>(tensor149, pindex);
  task148->add_dep(task149);
  task149->add_dep(task40);
  residualq->add_task(task149);

  auto tensor150 = vector<shared_ptr<Tensor>>{I145, f1_, t2};
  auto task150 = make_shared<Task150>(tensor150, pindex);
  task149->add_dep(task150);
  task150->add_dep(task40);
  residualq->add_task(task150);

  vector<IndexRange> I148_index = {virt_, active_, active_, active_};
  auto I148 = make_shared<Tensor>(I148_index);
  auto tensor151 = vector<shared_ptr<Tensor>>{I144, Gamma49_(), I148};
  auto task151 = make_shared<Task151>(tensor151, pindex);
  task148->add_dep(task151);
  task151->add_dep(task40);
  residualq->add_task(task151);

  auto tensor152 = vector<shared_ptr<Tensor>>{I148, f1_, t2};
  auto task152 = make_shared<Task152>(tensor152, pindex);
  task151->add_dep(task152);
  task152->add_dep(task40);
  residualq->add_task(task152);

  vector<IndexRange> I151_index = {active_, virt_};
  auto I151 = make_shared<Tensor>(I151_index);
  auto tensor153 = vector<shared_ptr<Tensor>>{I144, Gamma50_(), I151};
  auto task153 = make_shared<Task153>(tensor153, pindex);
  task148->add_dep(task153);
  task153->add_dep(task40);
  residualq->add_task(task153);

  vector<IndexRange> I152_index = {active_, virt_, closed_, virt_};
  auto I152 = make_shared<Tensor>(I152_index);
  auto tensor154 = vector<shared_ptr<Tensor>>{I151, f1_, I152};
  auto task154 = make_shared<Task154>(tensor154, pindex);
  task153->add_dep(task154);
  task154->add_dep(task40);
  residualq->add_task(task154);

  auto tensor155 = vector<shared_ptr<Tensor>>{I152, t2};
  auto task155 = make_shared<Task155>(tensor155, pindex);
  task154->add_dep(task155);
  task155->add_dep(task40);
  residualq->add_task(task155);

  vector<IndexRange> I157_index = {active_, virt_, active_, active_};
  auto I157 = make_shared<Tensor>(I157_index);
  auto tensor156 = vector<shared_ptr<Tensor>>{I144, Gamma52_(), I157};
  auto task156 = make_shared<Task156>(tensor156, pindex);
  task148->add_dep(task156);
  task156->add_dep(task40);
  residualq->add_task(task156);

  auto tensor157 = vector<shared_ptr<Tensor>>{I157, t2};
  auto task157 = make_shared<Task157>(tensor157, pindex, this->e0_);
  task156->add_dep(task157);
  task157->add_dep(task40);
  residualq->add_task(task157);

  auto tensor158 = vector<shared_ptr<Tensor>>{I157, f1_, t2};
  auto task158 = make_shared<Task158>(tensor158, pindex);
  task156->add_dep(task158);
  task158->add_dep(task40);
  residualq->add_task(task158);

  vector<IndexRange> I160_index = {virt_, active_, active_, active_};
  auto I160 = make_shared<Tensor>(I160_index);
  auto tensor159 = vector<shared_ptr<Tensor>>{I144, f1_, I160};
  auto task159 = make_shared<Task159>(tensor159, pindex);
  task148->add_dep(task159);
  task159->add_dep(task40);
  residualq->add_task(task159);

  auto tensor160 = vector<shared_ptr<Tensor>>{I160, Gamma52_(), t2};
  auto task160 = make_shared<Task160>(tensor160, pindex);
  task159->add_dep(task160);
  task160->add_dep(task40);
  residualq->add_task(task160);

  auto tensor161 = vector<shared_ptr<Tensor>>{I144, Gamma100_(), t2};
  auto task161 = make_shared<Task161>(tensor161, pindex);
  task148->add_dep(task161);
  task161->add_dep(task40);
  residualq->add_task(task161);

  auto tensor162 = vector<shared_ptr<Tensor>>{I144, Gamma101_(), t2};
  auto task162 = make_shared<Task162>(tensor162, pindex);
  task148->add_dep(task162);
  task162->add_dep(task40);
  residualq->add_task(task162);

  vector<IndexRange> I162_index = {closed_, virt_, virt_, closed_};
  auto I162 = make_shared<Tensor>(I162_index);
  auto tensor163 = vector<shared_ptr<Tensor>>{r, I162};
  auto task163 = make_shared<Task163>(tensor163, pindex);
  task163->add_dep(task40);
  residualq->add_task(task163);

  shared_ptr<Tensor> I163;
  if (diagonal) {
    vector<IndexRange> I163_index = {closed_, virt_, closed_, virt_};
    I163 = make_shared<Tensor>(I163_index);
  }
  shared_ptr<Task164> task164;
  if (diagonal) {
    auto tensor164 = vector<shared_ptr<Tensor>>{I162, f1_, I163};
    task164 = make_shared<Task164>(tensor164, pindex);
    task163->add_dep(task164);
    task164->add_dep(task40);
    residualq->add_task(task164);
  }

  shared_ptr<Task165> task165;
  if (diagonal) {
    auto tensor165 = vector<shared_ptr<Tensor>>{I163, t2};
    task165 = make_shared<Task165>(tensor165, pindex);
    task164->add_dep(task165);
    task165->add_dep(task40);
    residualq->add_task(task165);
  }

  vector<IndexRange> I167_index = {closed_, active_};
  auto I167 = make_shared<Tensor>(I167_index);
  auto tensor166 = vector<shared_ptr<Tensor>>{I162, t2, I167};
  auto task166 = make_shared<Task166>(tensor166, pindex);
  task163->add_dep(task166);
  task166->add_dep(task40);
  residualq->add_task(task166);

  auto tensor167 = vector<shared_ptr<Tensor>>{I167, Gamma29_(), f1_};
  auto task167 = make_shared<Task167>(tensor167, pindex);
  task166->add_dep(task167);
  task167->add_dep(task40);
  residualq->add_task(task167);

  vector<IndexRange> I170_index = {closed_, active_};
  auto I170 = make_shared<Tensor>(I170_index);
  auto tensor168 = vector<shared_ptr<Tensor>>{I162, t2, I170};
  auto task168 = make_shared<Task168>(tensor168, pindex);
  task163->add_dep(task168);
  task168->add_dep(task40);
  residualq->add_task(task168);

  auto tensor169 = vector<shared_ptr<Tensor>>{I170, Gamma29_(), f1_};
  auto task169 = make_shared<Task169>(tensor169, pindex);
  task168->add_dep(task169);
  task169->add_dep(task40);
  residualq->add_task(task169);

  vector<IndexRange> I173_index = {virt_, closed_};
  auto I173 = make_shared<Tensor>(I173_index);
  auto tensor170 = vector<shared_ptr<Tensor>>{I162, f1_, I173};
  auto task170 = make_shared<Task170>(tensor170, pindex);
  task163->add_dep(task170);
  task170->add_dep(task40);
  residualq->add_task(task170);

  vector<IndexRange> I174_index = {active_, virt_, closed_, active_};
  auto I174 = make_shared<Tensor>(I174_index);
  auto tensor171 = vector<shared_ptr<Tensor>>{I173, Gamma29_(), I174};
  auto task171 = make_shared<Task171>(tensor171, pindex);
  task170->add_dep(task171);
  task171->add_dep(task40);
  residualq->add_task(task171);

  auto tensor172 = vector<shared_ptr<Tensor>>{I174, t2};
  auto task172 = make_shared<Task172>(tensor172, pindex);
  task171->add_dep(task172);
  task172->add_dep(task40);
  residualq->add_task(task172);

  vector<IndexRange> I176_index = {virt_, closed_};
  auto I176 = make_shared<Tensor>(I176_index);
  auto tensor173 = vector<shared_ptr<Tensor>>{I162, f1_, I176};
  auto task173 = make_shared<Task173>(tensor173, pindex);
  task163->add_dep(task173);
  task173->add_dep(task40);
  residualq->add_task(task173);

  vector<IndexRange> I177_index = {active_, virt_, closed_, active_};
  auto I177 = make_shared<Tensor>(I177_index);
  auto tensor174 = vector<shared_ptr<Tensor>>{I176, Gamma29_(), I177};
  auto task174 = make_shared<Task174>(tensor174, pindex);
  task173->add_dep(task174);
  task174->add_dep(task40);
  residualq->add_task(task174);

  auto tensor175 = vector<shared_ptr<Tensor>>{I177, t2};
  auto task175 = make_shared<Task175>(tensor175, pindex);
  task174->add_dep(task175);
  task175->add_dep(task40);
  residualq->add_task(task175);

  vector<IndexRange> I185_index = {virt_, active_};
  auto I185 = make_shared<Tensor>(I185_index);
  auto tensor176 = vector<shared_ptr<Tensor>>{I162, t2, I185};
  auto task176 = make_shared<Task176>(tensor176, pindex);
  task163->add_dep(task176);
  task176->add_dep(task40);
  residualq->add_task(task176);

  auto tensor177 = vector<shared_ptr<Tensor>>{I185, Gamma9_(), f1_};
  auto task177 = make_shared<Task177>(tensor177, pindex);
  task176->add_dep(task177);
  task177->add_dep(task40);
  residualq->add_task(task177);

  vector<IndexRange> I188_index = {virt_, active_};
  auto I188 = make_shared<Tensor>(I188_index);
  auto tensor178 = vector<shared_ptr<Tensor>>{I162, t2, I188};
  auto task178 = make_shared<Task178>(tensor178, pindex);
  task163->add_dep(task178);
  task178->add_dep(task40);
  residualq->add_task(task178);

  auto tensor179 = vector<shared_ptr<Tensor>>{I188, Gamma9_(), f1_};
  auto task179 = make_shared<Task179>(tensor179, pindex);
  task178->add_dep(task179);
  task179->add_dep(task40);
  residualq->add_task(task179);

  shared_ptr<Tensor> I191;
  if (diagonal) {
    vector<IndexRange> I191_index = {closed_, virt_, closed_, virt_};
    I191 = make_shared<Tensor>(I191_index);
  }
  shared_ptr<Task180> task180;
  if (diagonal) {
    auto tensor180 = vector<shared_ptr<Tensor>>{I162, f1_, I191};
    task180 = make_shared<Task180>(tensor180, pindex);
    task163->add_dep(task180);
    task180->add_dep(task40);
    residualq->add_task(task180);
  }

  shared_ptr<Task181> task181;
  if (diagonal) {
    auto tensor181 = vector<shared_ptr<Tensor>>{I191, t2};
    task181 = make_shared<Task181>(tensor181, pindex);
    task180->add_dep(task181);
    task181->add_dep(task40);
    residualq->add_task(task181);
  }

  vector<IndexRange> I194_index = {virt_, virt_, active_, closed_};
  auto I194 = make_shared<Tensor>(I194_index);
  auto tensor182 = vector<shared_ptr<Tensor>>{r, I194};
  auto task182 = make_shared<Task182>(tensor182, pindex);
  task182->add_dep(task40);
  residualq->add_task(task182);

  vector<IndexRange> I195_index = {virt_, closed_, virt_, active_};
  auto I195 = make_shared<Tensor>(I195_index);
  auto tensor183 = vector<shared_ptr<Tensor>>{I194, f1_, I195};
  auto task183 = make_shared<Task183>(tensor183, pindex);
  task182->add_dep(task183);
  task183->add_dep(task40);
  residualq->add_task(task183);

  vector<IndexRange> I196_index = {active_, virt_, closed_, virt_};
  auto I196 = make_shared<Tensor>(I196_index);
  auto tensor184 = vector<shared_ptr<Tensor>>{I195, Gamma29_(), I196};
  auto task184 = make_shared<Task184>(tensor184, pindex);
  task183->add_dep(task184);
  task184->add_dep(task40);
  residualq->add_task(task184);

  auto tensor185 = vector<shared_ptr<Tensor>>{I196, t2};
  auto task185 = make_shared<Task185>(tensor185, pindex);
  task184->add_dep(task185);
  task185->add_dep(task40);
  residualq->add_task(task185);

  vector<IndexRange> I201_index = {closed_, active_};
  auto I201 = make_shared<Tensor>(I201_index);
  auto tensor186 = vector<shared_ptr<Tensor>>{I194, t2, I201};
  auto task186 = make_shared<Task186>(tensor186, pindex);
  task182->add_dep(task186);
  task186->add_dep(task40);
  residualq->add_task(task186);

  auto tensor187 = vector<shared_ptr<Tensor>>{I201, Gamma29_(), f1_};
  auto task187 = make_shared<Task187>(tensor187, pindex);
  task186->add_dep(task187);
  task187->add_dep(task40);
  residualq->add_task(task187);

  vector<IndexRange> I204_index = {closed_, active_};
  auto I204 = make_shared<Tensor>(I204_index);
  auto tensor188 = vector<shared_ptr<Tensor>>{I194, t2, I204};
  auto task188 = make_shared<Task188>(tensor188, pindex);
  task182->add_dep(task188);
  task188->add_dep(task40);
  residualq->add_task(task188);

  auto tensor189 = vector<shared_ptr<Tensor>>{I204, Gamma29_(), f1_};
  auto task189 = make_shared<Task189>(tensor189, pindex);
  task188->add_dep(task189);
  task189->add_dep(task40);
  residualq->add_task(task189);

  vector<IndexRange> I207_index = {virt_, virt_, active_, active_};
  auto I207 = make_shared<Tensor>(I207_index);
  auto tensor190 = vector<shared_ptr<Tensor>>{I194, f1_, I207};
  auto task190 = make_shared<Task190>(tensor190, pindex);
  task182->add_dep(task190);
  task190->add_dep(task40);
  residualq->add_task(task190);

  auto tensor191 = vector<shared_ptr<Tensor>>{I207, Gamma50_(), t2};
  auto task191 = make_shared<Task191>(tensor191, pindex);
  task190->add_dep(task191);
  task191->add_dep(task40);
  residualq->add_task(task191);

  vector<IndexRange> I210_index = {virt_, active_};
  auto I210 = make_shared<Tensor>(I210_index);
  auto tensor192 = vector<shared_ptr<Tensor>>{I194, f1_, I210};
  auto task192 = make_shared<Task192>(tensor192, pindex);
  task182->add_dep(task192);
  task192->add_dep(task40);
  residualq->add_task(task192);

  auto tensor193 = vector<shared_ptr<Tensor>>{I210, Gamma50_(), t2};
  auto task193 = make_shared<Task193>(tensor193, pindex);
  task192->add_dep(task193);
  task193->add_dep(task40);
  residualq->add_task(task193);

  vector<IndexRange> I213_index = {virt_, active_};
  auto I213 = make_shared<Tensor>(I213_index);
  auto tensor194 = vector<shared_ptr<Tensor>>{I194, f1_, I213};
  auto task194 = make_shared<Task194>(tensor194, pindex);
  task182->add_dep(task194);
  task194->add_dep(task40);
  residualq->add_task(task194);

  auto tensor195 = vector<shared_ptr<Tensor>>{I213, Gamma50_(), t2};
  auto task195 = make_shared<Task195>(tensor195, pindex);
  task194->add_dep(task195);
  task195->add_dep(task40);
  residualq->add_task(task195);

  vector<IndexRange> I216_index = {virt_, closed_, active_, active_};
  auto I216 = make_shared<Tensor>(I216_index);
  auto tensor196 = vector<shared_ptr<Tensor>>{I194, f1_, I216};
  auto task196 = make_shared<Task196>(tensor196, pindex);
  task182->add_dep(task196);
  task196->add_dep(task40);
  residualq->add_task(task196);

  vector<IndexRange> I217_index = {active_, virt_, closed_, active_};
  auto I217 = make_shared<Tensor>(I217_index);
  auto tensor197 = vector<shared_ptr<Tensor>>{I216, Gamma25_(), I217};
  auto task197 = make_shared<Task197>(tensor197, pindex);
  task196->add_dep(task197);
  task197->add_dep(task40);
  residualq->add_task(task197);

  auto tensor198 = vector<shared_ptr<Tensor>>{I217, t2};
  auto task198 = make_shared<Task198>(tensor198, pindex);
  task197->add_dep(task198);
  task198->add_dep(task40);
  residualq->add_task(task198);

  vector<IndexRange> I219_index = {virt_, closed_, active_, active_};
  auto I219 = make_shared<Tensor>(I219_index);
  auto tensor199 = vector<shared_ptr<Tensor>>{I194, f1_, I219};
  auto task199 = make_shared<Task199>(tensor199, pindex);
  task182->add_dep(task199);
  task199->add_dep(task40);
  residualq->add_task(task199);

  auto tensor200 = vector<shared_ptr<Tensor>>{I219, Gamma24_(), t2};
  auto task200 = make_shared<Task200>(tensor200, pindex);
  task199->add_dep(task200);
  task200->add_dep(task40);
  residualq->add_task(task200);

  auto tensor201 = vector<shared_ptr<Tensor>>{I219, Gamma25_(), t2};
  auto task201 = make_shared<Task201>(tensor201, pindex);
  task199->add_dep(task201);
  task201->add_dep(task40);
  residualq->add_task(task201);

  vector<IndexRange> I228_index = {virt_, closed_, virt_, active_};
  auto I228 = make_shared<Tensor>(I228_index);
  auto tensor202 = vector<shared_ptr<Tensor>>{I194, f1_, I228};
  auto task202 = make_shared<Task202>(tensor202, pindex);
  task182->add_dep(task202);
  task202->add_dep(task40);
  residualq->add_task(task202);

  vector<IndexRange> I229_index = {active_, virt_, closed_, virt_};
  auto I229 = make_shared<Tensor>(I229_index);
  auto tensor203 = vector<shared_ptr<Tensor>>{I228, Gamma29_(), I229};
  auto task203 = make_shared<Task203>(tensor203, pindex);
  task202->add_dep(task203);
  task203->add_dep(task40);
  residualq->add_task(task203);

  auto tensor204 = vector<shared_ptr<Tensor>>{I229, t2};
  auto task204 = make_shared<Task204>(tensor204, pindex);
  task203->add_dep(task204);
  task204->add_dep(task40);
  residualq->add_task(task204);

  vector<IndexRange> I234_index = {virt_, closed_, virt_, active_};
  auto I234 = make_shared<Tensor>(I234_index);
  auto tensor205 = vector<shared_ptr<Tensor>>{I194, f1_, I234};
  auto task205 = make_shared<Task205>(tensor205, pindex);
  task182->add_dep(task205);
  task205->add_dep(task40);
  residualq->add_task(task205);

  vector<IndexRange> I235_index = {active_, virt_, closed_, virt_};
  auto I235 = make_shared<Tensor>(I235_index);
  auto tensor206 = vector<shared_ptr<Tensor>>{I234, Gamma29_(), I235};
  auto task206 = make_shared<Task206>(tensor206, pindex);
  task205->add_dep(task206);
  task206->add_dep(task40);
  residualq->add_task(task206);

  auto tensor207 = vector<shared_ptr<Tensor>>{I235, t2};
  auto task207 = make_shared<Task207>(tensor207, pindex);
  task206->add_dep(task207);
  task207->add_dep(task40);
  residualq->add_task(task207);

  vector<IndexRange> I269_index = {active_, virt_, closed_, virt_};
  auto I269 = make_shared<Tensor>(I269_index);
  auto tensor208 = vector<shared_ptr<Tensor>>{I194, Gamma29_(), I269};
  auto task208 = make_shared<Task208>(tensor208, pindex);
  task182->add_dep(task208);
  task208->add_dep(task40);
  residualq->add_task(task208);

  auto tensor209 = vector<shared_ptr<Tensor>>{I269, t2};
  auto task209 = make_shared<Task209>(tensor209, pindex, this->e0_);
  task208->add_dep(task209);
  task209->add_dep(task40);
  residualq->add_task(task209);

  vector<IndexRange> I303_index = {active_, virt_, closed_, virt_};
  auto I303 = make_shared<Tensor>(I303_index);
  auto tensor210 = vector<shared_ptr<Tensor>>{I194, Gamma106_(), I303};
  auto task210 = make_shared<Task210>(tensor210, pindex);
  task182->add_dep(task210);
  task210->add_dep(task40);
  residualq->add_task(task210);

  auto tensor211 = vector<shared_ptr<Tensor>>{I303, t2};
  auto task211 = make_shared<Task211>(tensor211, pindex);
  task210->add_dep(task211);
  task211->add_dep(task40);
  residualq->add_task(task211);

  vector<IndexRange> I307_index = {active_, virt_, closed_, virt_};
  auto I307 = make_shared<Tensor>(I307_index);
  auto tensor212 = vector<shared_ptr<Tensor>>{I194, Gamma108_(), I307};
  auto task212 = make_shared<Task212>(tensor212, pindex);
  task182->add_dep(task212);
  task212->add_dep(task40);
  residualq->add_task(task212);

  auto tensor213 = vector<shared_ptr<Tensor>>{I307, t2};
  auto task213 = make_shared<Task213>(tensor213, pindex);
  task212->add_dep(task213);
  task213->add_dep(task40);
  residualq->add_task(task213);

  vector<IndexRange> I239_index = {virt_, virt_, active_, active_};
  auto I239 = make_shared<Tensor>(I239_index);
  auto tensor214 = vector<shared_ptr<Tensor>>{r, I239};
  auto task214 = make_shared<Task214>(tensor214, pindex);
  task214->add_dep(task40);
  residualq->add_task(task214);

  vector<IndexRange> I240_index = {active_, virt_, virt_, active_};
  auto I240 = make_shared<Tensor>(I240_index);
  auto tensor215 = vector<shared_ptr<Tensor>>{I239, Gamma50_(), I240};
  auto task215 = make_shared<Task215>(tensor215, pindex);
  task214->add_dep(task215);
  task215->add_dep(task40);
  residualq->add_task(task215);

  auto tensor216 = vector<shared_ptr<Tensor>>{I240, f1_, t2};
  auto task216 = make_shared<Task216>(tensor216, pindex);
  task215->add_dep(task216);
  task216->add_dep(task40);
  residualq->add_task(task216);

  vector<IndexRange> I243_index = {virt_, active_, active_, active_};
  auto I243 = make_shared<Tensor>(I243_index);
  auto tensor217 = vector<shared_ptr<Tensor>>{I239, f1_, I243};
  auto task217 = make_shared<Task217>(tensor217, pindex);
  task214->add_dep(task217);
  task217->add_dep(task40);
  residualq->add_task(task217);

  auto tensor218 = vector<shared_ptr<Tensor>>{I243, Gamma52_(), t2};
  auto task218 = make_shared<Task218>(tensor218, pindex);
  task217->add_dep(task218);
  task218->add_dep(task40);
  residualq->add_task(task218);

  vector<IndexRange> I246_index = {virt_, virt_, active_, active_};
  auto I246 = make_shared<Tensor>(I246_index);
  auto tensor219 = vector<shared_ptr<Tensor>>{I239, f1_, I246};
  auto task219 = make_shared<Task219>(tensor219, pindex);
  task214->add_dep(task219);
  task219->add_dep(task40);
  residualq->add_task(task219);

  auto tensor220 = vector<shared_ptr<Tensor>>{I246, Gamma50_(), t2};
  auto task220 = make_shared<Task220>(tensor220, pindex);
  task219->add_dep(task220);
  task220->add_dep(task40);
  residualq->add_task(task220);

  vector<IndexRange> I248_index = {closed_, closed_, active_, active_};
  auto I248 = make_shared<Tensor>(I248_index);
  auto tensor221 = vector<shared_ptr<Tensor>>{r, I248};
  auto task221 = make_shared<Task221>(tensor221, pindex);
  task221->add_dep(task40);
  residualq->add_task(task221);

  auto tensor222 = vector<shared_ptr<Tensor>>{I248, Gamma0_(), t2};
  auto task222 = make_shared<Task222>(tensor222, pindex, this->e0_);
  task221->add_dep(task222);
  task222->add_dep(task40);
  residualq->add_task(task222);

  auto tensor223 = vector<shared_ptr<Tensor>>{I248, Gamma92_(), t2};
  auto task223 = make_shared<Task223>(tensor223, pindex);
  task221->add_dep(task223);
  task223->add_dep(task40);
  residualq->add_task(task223);

  auto tensor224 = vector<shared_ptr<Tensor>>{I248, Gamma93_(), t2};
  auto task224 = make_shared<Task224>(tensor224, pindex);
  task221->add_dep(task224);
  task224->add_dep(task40);
  residualq->add_task(task224);

  vector<IndexRange> I266_index = {closed_, virt_, closed_, virt_};
  auto I266 = make_shared<Tensor>(I266_index);
  auto tensor225 = vector<shared_ptr<Tensor>>{r, I266};
  auto task225 = make_shared<Task225>(tensor225, pindex);
  task225->add_dep(task40);
  residualq->add_task(task225);

  shared_ptr<Task226> task226;
  if (diagonal) {
    auto tensor226 = vector<shared_ptr<Tensor>>{I266, t2};
    task226 = make_shared<Task226>(tensor226, pindex, this->e0_);
    task225->add_dep(task226);
    task226->add_dep(task40);
    residualq->add_task(task226);
  }

  vector<IndexRange> I295_index = {closed_, virt_, closed_, virt_};
  auto I295 = make_shared<Tensor>(I295_index);
  auto tensor227 = vector<shared_ptr<Tensor>>{I266, Gamma102_(), I295};
  auto task227 = make_shared<Task227>(tensor227, pindex);
  task225->add_dep(task227);
  task227->add_dep(task40);
  residualq->add_task(task227);

  auto tensor228 = vector<shared_ptr<Tensor>>{I295, t2};
  auto task228 = make_shared<Task228>(tensor228, pindex);
  task227->add_dep(task228);
  task228->add_dep(task40);
  residualq->add_task(task228);

  vector<IndexRange> I299_index = {closed_, virt_, closed_, virt_};
  auto I299 = make_shared<Tensor>(I299_index);
  auto tensor229 = vector<shared_ptr<Tensor>>{I266, Gamma104_(), I299};
  auto task229 = make_shared<Task229>(tensor229, pindex);
  task225->add_dep(task229);
  task229->add_dep(task40);
  residualq->add_task(task229);

  auto tensor230 = vector<shared_ptr<Tensor>>{I299, t2};
  auto task230 = make_shared<Task230>(tensor230, pindex);
  task229->add_dep(task230);
  task230->add_dep(task40);
  residualq->add_task(task230);

  vector<IndexRange> I272_index = {virt_, virt_, active_, active_};
  auto I272 = make_shared<Tensor>(I272_index);
  auto tensor231 = vector<shared_ptr<Tensor>>{r, I272};
  auto task231 = make_shared<Task231>(tensor231, pindex);
  task231->add_dep(task40);
  residualq->add_task(task231);

  auto tensor232 = vector<shared_ptr<Tensor>>{I272, Gamma50_(), t2};
  auto task232 = make_shared<Task232>(tensor232, pindex, this->e0_);
  task231->add_dep(task232);
  task232->add_dep(task40);
  residualq->add_task(task232);

  auto tensor233 = vector<shared_ptr<Tensor>>{I272, Gamma110_(), t2};
  auto task233 = make_shared<Task233>(tensor233, pindex);
  task231->add_dep(task233);
  task233->add_dep(task40);
  residualq->add_task(task233);

  auto tensor234 = vector<shared_ptr<Tensor>>{I272, Gamma111_(), t2};
  auto task234 = make_shared<Task234>(tensor234, pindex);
  task231->add_dep(task234);
  task234->add_dep(task40);
  residualq->add_task(task234);

  return residualq;
}


#endif
