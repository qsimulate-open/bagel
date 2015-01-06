//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_densityqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_densityq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto densityq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor57 = {den2};
  auto task57 = make_shared<Task57>(tensor57);
  densityq->add_task(task57);

  vector<IndexRange> I64_index = {active_, active_};
  auto I64 = make_shared<Tensor>(I64_index);
  vector<shared_ptr<Tensor>> tensor58 = {den2, I64};
  auto task58 = make_shared<Task58>(tensor58, pindex);
  task58->add_dep(task57);
  densityq->add_task(task58);

  vector<IndexRange> I65_index;
  auto I65 = make_shared<Tensor>(I65_index);
  vector<shared_ptr<Tensor>> tensor59 = {I64, Gamma0_(), I65};
  auto task59 = make_shared<Task59>(tensor59, pindex);
  task58->add_dep(task59);
  task59->add_dep(task57);
  densityq->add_task(task59);

  vector<IndexRange> I66_index = {virt_, closed_, virt_, closed_};
  auto I66 = make_shared<Tensor>(I66_index);
  vector<shared_ptr<Tensor>> tensor60 = {I65, t2, I66};
  auto task60 = make_shared<Task60>(tensor60, pindex);
  task59->add_dep(task60);
  task60->add_dep(task57);
  densityq->add_task(task60);

  vector<shared_ptr<Tensor>> tensor61 = {I66, t2};
  auto task61 = make_shared<Task61>(tensor61, pindex);
  task60->add_dep(task61);
  task61->add_dep(task57);
  densityq->add_task(task61);

  vector<IndexRange> I69_index = {closed_, virt_, closed_, virt_};
  auto I69 = make_shared<Tensor>(I69_index);
  vector<shared_ptr<Tensor>> tensor62 = {I65, t2, I69};
  auto task62 = make_shared<Task62>(tensor62, pindex);
  task59->add_dep(task62);
  task62->add_dep(task57);
  densityq->add_task(task62);

  vector<shared_ptr<Tensor>> tensor63 = {I69, t2};
  auto task63 = make_shared<Task63>(tensor63, pindex);
  task62->add_dep(task63);
  task63->add_dep(task57);
  densityq->add_task(task63);

  vector<IndexRange> I70_index = {closed_, closed_};
  auto I70 = make_shared<Tensor>(I70_index);
  vector<shared_ptr<Tensor>> tensor64 = {den2, I70};
  auto task64 = make_shared<Task64>(tensor64, pindex);
  task64->add_dep(task57);
  densityq->add_task(task64);

  vector<IndexRange> I71_index = {closed_, virt_, closed_, virt_};
  auto I71 = make_shared<Tensor>(I71_index);
  vector<shared_ptr<Tensor>> tensor65 = {I70, t2, I71};
  auto task65 = make_shared<Task65>(tensor65, pindex);
  task64->add_dep(task65);
  task65->add_dep(task57);
  densityq->add_task(task65);

  vector<shared_ptr<Tensor>> tensor66 = {I71, t2};
  auto task66 = make_shared<Task66>(tensor66, pindex);
  task65->add_dep(task66);
  task66->add_dep(task57);
  densityq->add_task(task66);

  vector<IndexRange> I73_index = {virt_, closed_, virt_, closed_};
  auto I73 = make_shared<Tensor>(I73_index);
  vector<shared_ptr<Tensor>> tensor67 = {I70, t2, I73};
  auto task67 = make_shared<Task67>(tensor67, pindex);
  task64->add_dep(task67);
  task67->add_dep(task57);
  densityq->add_task(task67);

  vector<shared_ptr<Tensor>> tensor68 = {I73, t2};
  auto task68 = make_shared<Task68>(tensor68, pindex);
  task67->add_dep(task68);
  task68->add_dep(task57);
  densityq->add_task(task68);

  vector<IndexRange> I74_index = {virt_, virt_};
  auto I74 = make_shared<Tensor>(I74_index);
  vector<shared_ptr<Tensor>> tensor69 = {den2, I74};
  auto task69 = make_shared<Task69>(tensor69, pindex);
  task69->add_dep(task57);
  densityq->add_task(task69);

  vector<IndexRange> I75_index = {virt_, closed_, virt_, closed_};
  auto I75 = make_shared<Tensor>(I75_index);
  vector<shared_ptr<Tensor>> tensor70 = {I74, t2, I75};
  auto task70 = make_shared<Task70>(tensor70, pindex);
  task69->add_dep(task70);
  task70->add_dep(task57);
  densityq->add_task(task70);

  vector<shared_ptr<Tensor>> tensor71 = {I75, t2};
  auto task71 = make_shared<Task71>(tensor71, pindex);
  task70->add_dep(task71);
  task71->add_dep(task57);
  densityq->add_task(task71);

  vector<IndexRange> I77_index = {virt_, closed_, virt_, closed_};
  auto I77 = make_shared<Tensor>(I77_index);
  vector<shared_ptr<Tensor>> tensor72 = {I74, t2, I77};
  auto task72 = make_shared<Task72>(tensor72, pindex);
  task69->add_dep(task72);
  task72->add_dep(task57);
  densityq->add_task(task72);

  vector<shared_ptr<Tensor>> tensor73 = {I77, t2};
  auto task73 = make_shared<Task73>(tensor73, pindex);
  task72->add_dep(task73);
  task73->add_dep(task57);
  densityq->add_task(task73);

  vector<IndexRange> I78_index = {closed_, active_};
  auto I78 = make_shared<Tensor>(I78_index);
  vector<shared_ptr<Tensor>> tensor74 = {den2, I78};
  auto task74 = make_shared<Task74>(tensor74, pindex);
  task74->add_dep(task57);
  densityq->add_task(task74);

  vector<IndexRange> I79_index = {closed_, active_};
  auto I79 = make_shared<Tensor>(I79_index);
  vector<shared_ptr<Tensor>> tensor75 = {I78, Gamma0_(), I79};
  auto task75 = make_shared<Task75>(tensor75, pindex);
  task74->add_dep(task75);
  task75->add_dep(task57);
  densityq->add_task(task75);

  vector<IndexRange> I80_index = {virt_, closed_, virt_, closed_};
  auto I80 = make_shared<Tensor>(I80_index);
  vector<shared_ptr<Tensor>> tensor76 = {I79, t2, I80};
  auto task76 = make_shared<Task76>(tensor76, pindex);
  task75->add_dep(task76);
  task76->add_dep(task57);
  densityq->add_task(task76);

  vector<shared_ptr<Tensor>> tensor77 = {I80, t2};
  auto task77 = make_shared<Task77>(tensor77, pindex);
  task76->add_dep(task77);
  task77->add_dep(task57);
  densityq->add_task(task77);

  vector<IndexRange> I83_index = {virt_, closed_, virt_, closed_};
  auto I83 = make_shared<Tensor>(I83_index);
  vector<shared_ptr<Tensor>> tensor78 = {I79, t2, I83};
  auto task78 = make_shared<Task78>(tensor78, pindex);
  task75->add_dep(task78);
  task78->add_dep(task57);
  densityq->add_task(task78);

  vector<shared_ptr<Tensor>> tensor79 = {I83, t2};
  auto task79 = make_shared<Task79>(tensor79, pindex);
  task78->add_dep(task79);
  task79->add_dep(task57);
  densityq->add_task(task79);

  vector<IndexRange> I84_index = {closed_, active_};
  auto I84 = make_shared<Tensor>(I84_index);
  vector<shared_ptr<Tensor>> tensor80 = {den2, I84};
  auto task80 = make_shared<Task80>(tensor80, pindex);
  task80->add_dep(task57);
  densityq->add_task(task80);

  vector<IndexRange> I85_index = {active_, closed_};
  auto I85 = make_shared<Tensor>(I85_index);
  vector<shared_ptr<Tensor>> tensor81 = {I84, Gamma0_(), I85};
  auto task81 = make_shared<Task81>(tensor81, pindex);
  task80->add_dep(task81);
  task81->add_dep(task57);
  densityq->add_task(task81);

  vector<IndexRange> I86_index = {virt_, closed_, virt_, active_};
  auto I86 = make_shared<Tensor>(I86_index);
  vector<shared_ptr<Tensor>> tensor82 = {I85, t2, I86};
  auto task82 = make_shared<Task82>(tensor82, pindex);
  task81->add_dep(task82);
  task82->add_dep(task57);
  densityq->add_task(task82);

  vector<shared_ptr<Tensor>> tensor83 = {I86, t2};
  auto task83 = make_shared<Task83>(tensor83, pindex);
  task82->add_dep(task83);
  task83->add_dep(task57);
  densityq->add_task(task83);

  vector<IndexRange> I89_index = {virt_, closed_, virt_, active_};
  auto I89 = make_shared<Tensor>(I89_index);
  vector<shared_ptr<Tensor>> tensor84 = {I85, t2, I89};
  auto task84 = make_shared<Task84>(tensor84, pindex);
  task81->add_dep(task84);
  task84->add_dep(task57);
  densityq->add_task(task84);

  vector<shared_ptr<Tensor>> tensor85 = {I89, t2};
  auto task85 = make_shared<Task85>(tensor85, pindex);
  task84->add_dep(task85);
  task85->add_dep(task57);
  densityq->add_task(task85);

  vector<IndexRange> I90_index = {active_, active_};
  auto I90 = make_shared<Tensor>(I90_index);
  vector<shared_ptr<Tensor>> tensor86 = {den2, I90};
  auto task86 = make_shared<Task86>(tensor86, pindex);
  task86->add_dep(task57);
  densityq->add_task(task86);

  vector<IndexRange> I91_index = {active_, active_};
  auto I91 = make_shared<Tensor>(I91_index);
  vector<shared_ptr<Tensor>> tensor87 = {I90, Gamma26_(), I91};
  auto task87 = make_shared<Task87>(tensor87, pindex);
  task86->add_dep(task87);
  task87->add_dep(task57);
  densityq->add_task(task87);

  vector<IndexRange> I92_index = {virt_, closed_, virt_, active_};
  auto I92 = make_shared<Tensor>(I92_index);
  vector<shared_ptr<Tensor>> tensor88 = {I91, t2, I92};
  auto task88 = make_shared<Task88>(tensor88, pindex);
  task87->add_dep(task88);
  task88->add_dep(task57);
  densityq->add_task(task88);

  vector<shared_ptr<Tensor>> tensor89 = {I92, t2};
  auto task89 = make_shared<Task89>(tensor89, pindex);
  task88->add_dep(task89);
  task89->add_dep(task57);
  densityq->add_task(task89);

  vector<IndexRange> I95_index = {virt_, closed_, virt_, active_};
  auto I95 = make_shared<Tensor>(I95_index);
  vector<shared_ptr<Tensor>> tensor90 = {I91, t2, I95};
  auto task90 = make_shared<Task90>(tensor90, pindex);
  task87->add_dep(task90);
  task90->add_dep(task57);
  densityq->add_task(task90);

  vector<shared_ptr<Tensor>> tensor91 = {I95, t2};
  auto task91 = make_shared<Task91>(tensor91, pindex);
  task90->add_dep(task91);
  task91->add_dep(task57);
  densityq->add_task(task91);

  vector<IndexRange> I96_index = {closed_, closed_};
  auto I96 = make_shared<Tensor>(I96_index);
  vector<shared_ptr<Tensor>> tensor92 = {den2, I96};
  auto task92 = make_shared<Task92>(tensor92, pindex);
  task92->add_dep(task57);
  densityq->add_task(task92);

  vector<IndexRange> I97_index = {virt_, closed_, virt_, active_};
  auto I97 = make_shared<Tensor>(I97_index);
  vector<shared_ptr<Tensor>> tensor93 = {I96, t2, I97};
  auto task93 = make_shared<Task93>(tensor93, pindex);
  task92->add_dep(task93);
  task93->add_dep(task57);
  densityq->add_task(task93);

  vector<IndexRange> I98_index = {virt_, closed_, virt_, active_};
  auto I98 = make_shared<Tensor>(I98_index);
  vector<shared_ptr<Tensor>> tensor94 = {I97, Gamma0_(), I98};
  auto task94 = make_shared<Task94>(tensor94, pindex);
  task93->add_dep(task94);
  task94->add_dep(task57);
  densityq->add_task(task94);

  vector<shared_ptr<Tensor>> tensor95 = {I98, t2};
  auto task95 = make_shared<Task95>(tensor95, pindex);
  task94->add_dep(task95);
  task95->add_dep(task57);
  densityq->add_task(task95);

  vector<IndexRange> I100_index = {virt_, closed_, virt_, active_};
  auto I100 = make_shared<Tensor>(I100_index);
  vector<shared_ptr<Tensor>> tensor96 = {I96, t2, I100};
  auto task96 = make_shared<Task96>(tensor96, pindex);
  task92->add_dep(task96);
  task96->add_dep(task57);
  densityq->add_task(task96);

  vector<IndexRange> I101_index = {virt_, closed_, virt_, active_};
  auto I101 = make_shared<Tensor>(I101_index);
  vector<shared_ptr<Tensor>> tensor97 = {I100, Gamma0_(), I101};
  auto task97 = make_shared<Task97>(tensor97, pindex);
  task96->add_dep(task97);
  task97->add_dep(task57);
  densityq->add_task(task97);

  vector<shared_ptr<Tensor>> tensor98 = {I101, t2};
  auto task98 = make_shared<Task98>(tensor98, pindex);
  task97->add_dep(task98);
  task98->add_dep(task57);
  densityq->add_task(task98);

  vector<IndexRange> I102_index = {virt_, virt_};
  auto I102 = make_shared<Tensor>(I102_index);
  vector<shared_ptr<Tensor>> tensor99 = {den2, I102};
  auto task99 = make_shared<Task99>(tensor99, pindex);
  task99->add_dep(task57);
  densityq->add_task(task99);

  vector<IndexRange> I103_index = {virt_, closed_, virt_, active_};
  auto I103 = make_shared<Tensor>(I103_index);
  vector<shared_ptr<Tensor>> tensor100 = {I102, t2, I103};
  auto task100 = make_shared<Task100>(tensor100, pindex);
  task99->add_dep(task100);
  task100->add_dep(task57);
  densityq->add_task(task100);

  vector<IndexRange> I104_index = {virt_, closed_, virt_, active_};
  auto I104 = make_shared<Tensor>(I104_index);
  vector<shared_ptr<Tensor>> tensor101 = {I103, Gamma0_(), I104};
  auto task101 = make_shared<Task101>(tensor101, pindex);
  task100->add_dep(task101);
  task101->add_dep(task57);
  densityq->add_task(task101);

  vector<shared_ptr<Tensor>> tensor102 = {I104, t2};
  auto task102 = make_shared<Task102>(tensor102, pindex);
  task101->add_dep(task102);
  task102->add_dep(task57);
  densityq->add_task(task102);

  vector<IndexRange> I106_index = {virt_, closed_, virt_, active_};
  auto I106 = make_shared<Tensor>(I106_index);
  vector<shared_ptr<Tensor>> tensor103 = {I102, t2, I106};
  auto task103 = make_shared<Task103>(tensor103, pindex);
  task99->add_dep(task103);
  task103->add_dep(task57);
  densityq->add_task(task103);

  vector<IndexRange> I107_index = {virt_, closed_, virt_, active_};
  auto I107 = make_shared<Tensor>(I107_index);
  vector<shared_ptr<Tensor>> tensor104 = {I106, Gamma0_(), I107};
  auto task104 = make_shared<Task104>(tensor104, pindex);
  task103->add_dep(task104);
  task104->add_dep(task57);
  densityq->add_task(task104);

  vector<shared_ptr<Tensor>> tensor105 = {I107, t2};
  auto task105 = make_shared<Task105>(tensor105, pindex);
  task104->add_dep(task105);
  task105->add_dep(task57);
  densityq->add_task(task105);

  vector<IndexRange> I108_index = {virt_, virt_};
  auto I108 = make_shared<Tensor>(I108_index);
  vector<shared_ptr<Tensor>> tensor106 = {den2, I108};
  auto task106 = make_shared<Task106>(tensor106, pindex);
  task106->add_dep(task57);
  densityq->add_task(task106);

  vector<IndexRange> I109_index = {virt_, closed_, virt_, active_};
  auto I109 = make_shared<Tensor>(I109_index);
  vector<shared_ptr<Tensor>> tensor107 = {I108, t2, I109};
  auto task107 = make_shared<Task107>(tensor107, pindex);
  task106->add_dep(task107);
  task107->add_dep(task57);
  densityq->add_task(task107);

  vector<IndexRange> I110_index = {virt_, closed_, virt_, active_};
  auto I110 = make_shared<Tensor>(I110_index);
  vector<shared_ptr<Tensor>> tensor108 = {I109, Gamma0_(), I110};
  auto task108 = make_shared<Task108>(tensor108, pindex);
  task107->add_dep(task108);
  task108->add_dep(task57);
  densityq->add_task(task108);

  vector<shared_ptr<Tensor>> tensor109 = {I110, t2};
  auto task109 = make_shared<Task109>(tensor109, pindex);
  task108->add_dep(task109);
  task109->add_dep(task57);
  densityq->add_task(task109);

  vector<IndexRange> I112_index = {virt_, closed_, virt_, active_};
  auto I112 = make_shared<Tensor>(I112_index);
  vector<shared_ptr<Tensor>> tensor110 = {I108, t2, I112};
  auto task110 = make_shared<Task110>(tensor110, pindex);
  task106->add_dep(task110);
  task110->add_dep(task57);
  densityq->add_task(task110);

  vector<IndexRange> I113_index = {virt_, closed_, virt_, active_};
  auto I113 = make_shared<Tensor>(I113_index);
  vector<shared_ptr<Tensor>> tensor111 = {I112, Gamma0_(), I113};
  auto task111 = make_shared<Task111>(tensor111, pindex);
  task110->add_dep(task111);
  task111->add_dep(task57);
  densityq->add_task(task111);

  vector<shared_ptr<Tensor>> tensor112 = {I113, t2};
  auto task112 = make_shared<Task112>(tensor112, pindex);
  task111->add_dep(task112);
  task112->add_dep(task57);
  densityq->add_task(task112);

  return densityq;
}


