//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2_residualqq.cc
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


#include <src/smith/RelCASPT2.h>
#include <src/smith/RelCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelCASPT2::RelCASPT2::make_residualq(const bool reset, const bool diagonal) {

  auto residualq = make_shared<Queue>();
  auto task31 = make_shared<Task31>(r, reset);
  residualq->add_task(task31);

  auto I0 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task32 = make_shared<Task32>(r, I0);
  task32->add_dep(task31);
  residualq->add_task(task32);

  auto task33 = make_shared<Task33>(I0, Gamma0_(), t2);
  task32->add_dep(task33);
  task33->add_dep(task31);
  residualq->add_task(task33);

  auto I277 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, closed_, active_});
  auto task34 = make_shared<Task34>(I0, Gamma94_(), I277);
  task32->add_dep(task34);
  task34->add_dep(task31);
  residualq->add_task(task34);

  auto task35 = make_shared<Task35>(I277, t2, v2_, this->e0_);
  task34->add_dep(task35);
  task35->add_dep(task31);
  residualq->add_task(task35);

  auto I2 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task36 = make_shared<Task36>(r, I2);
  task36->add_dep(task31);
  residualq->add_task(task36);

  auto I3 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task37 = make_shared<Task37>(I2, Gamma94_(), I3);
  task36->add_dep(task37);
  task37->add_dep(task31);
  residualq->add_task(task37);

  auto task38 = make_shared<Task38>(I3, t2, f1_);
  task37->add_dep(task38);
  task38->add_dep(task31);
  residualq->add_task(task38);

  auto I6 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task39 = make_shared<Task39>(I2, f1_, I6);
  task36->add_dep(task39);
  task39->add_dep(task31);
  residualq->add_task(task39);

  auto task40 = make_shared<Task40>(I6, Gamma2_(), t2);
  task39->add_dep(task40);
  task40->add_dep(task31);
  residualq->add_task(task40);

  auto I9 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, closed_, closed_, active_});
  auto task41 = make_shared<Task41>(I2, Gamma3_(), I9);
  task36->add_dep(task41);
  task41->add_dep(task31);
  residualq->add_task(task41);

  auto task42 = make_shared<Task42>(I9, t2, f1_);
  task41->add_dep(task42);
  task42->add_dep(task31);
  residualq->add_task(task42);

  auto I11 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task43 = make_shared<Task43>(r, I11);
  task43->add_dep(task31);
  residualq->add_task(task43);

  auto I12 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, closed_, active_, active_});
  auto task44 = make_shared<Task44>(I11, Gamma4_(), I12);
  task43->add_dep(task44);
  task44->add_dep(task31);
  residualq->add_task(task44);

  auto task45 = make_shared<Task45>(I12, t2, f1_);
  task44->add_dep(task45);
  task45->add_dep(task31);
  residualq->add_task(task45);

  auto task46 = make_shared<Task46>(I11, Gamma5_(), t2);
  task43->add_dep(task46);
  task46->add_dep(task31);
  residualq->add_task(task46);

  auto I17 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task47 = make_shared<Task47>(I11, Gamma6_(), I17);
  task43->add_dep(task47);
  task47->add_dep(task31);
  residualq->add_task(task47);

  auto task48 = make_shared<Task48>(I17, t2, v2_, this->e0_);
  task47->add_dep(task48);
  task48->add_dep(task31);
  residualq->add_task(task48);

  auto task49 = make_shared<Task49>(I17, t2, f1_);
  task47->add_dep(task49);
  task49->add_dep(task31);
  residualq->add_task(task49);

  auto task50 = make_shared<Task50>(I17, t2, f1_);
  task47->add_dep(task50);
  task50->add_dep(task31);
  residualq->add_task(task50);

  auto I20 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{closed_, active_});
  auto task51 = make_shared<Task51>(I11, Gamma7_(), I20);
  task43->add_dep(task51);
  task51->add_dep(task31);
  residualq->add_task(task51);

  auto task52 = make_shared<Task52>(I20, h1_);
  task51->add_dep(task52);
  task52->add_dep(task31);
  residualq->add_task(task52);

  auto task53 = make_shared<Task53>(I20, t2, f1_);
  task51->add_dep(task53);
  task53->add_dep(task31);
  residualq->add_task(task53);

  auto task54 = make_shared<Task54>(I20, t2, f1_);
  task51->add_dep(task54);
  task54->add_dep(task31);
  residualq->add_task(task54);

  auto I26 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, active_, closed_, active_});
  auto task55 = make_shared<Task55>(I11, Gamma9_(), I26);
  task43->add_dep(task55);
  task55->add_dep(task31);
  residualq->add_task(task55);

  auto task56 = make_shared<Task56>(I26, t2, f1_);
  task55->add_dep(task56);
  task56->add_dep(task31);
  residualq->add_task(task56);

  auto task57 = make_shared<Task57>(I11, Gamma107_(), v2_);
  task43->add_dep(task57);
  task57->add_dep(task31);
  residualq->add_task(task57);

  auto I31 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, active_, virt_});
  auto task58 = make_shared<Task58>(r, I31);
  task58->add_dep(task31);
  residualq->add_task(task58);

  auto I32 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task59 = make_shared<Task59>(I31, f1_, I32);
  task58->add_dep(task59);
  task59->add_dep(task31);
  residualq->add_task(task59);

  auto task60 = make_shared<Task60>(I32, Gamma3_(), t2);
  task59->add_dep(task60);
  task60->add_dep(task31);
  residualq->add_task(task60);

  auto I35 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{closed_, active_});
  auto task61 = make_shared<Task61>(I31, f1_, I35);
  task58->add_dep(task61);
  task61->add_dep(task31);
  residualq->add_task(task61);

  auto task62 = make_shared<Task62>(I35, Gamma12_(), t2);
  task61->add_dep(task62);
  task62->add_dep(task31);
  residualq->add_task(task62);

  auto I38 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{closed_, active_});
  auto task63 = make_shared<Task63>(I31, f1_, I38);
  task58->add_dep(task63);
  task63->add_dep(task31);
  residualq->add_task(task63);

  auto task64 = make_shared<Task64>(I38, Gamma12_(), t2);
  task63->add_dep(task64);
  task64->add_dep(task31);
  residualq->add_task(task64);

  auto I41 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task65 = make_shared<Task65>(I31, Gamma14_(), I41);
  task58->add_dep(task65);
  task65->add_dep(task31);
  residualq->add_task(task65);

  auto task66 = make_shared<Task66>(I41, t2);
  task65->add_dep(task66);
  task66->add_dep(task31);
  residualq->add_task(task66);

  auto I45 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task67 = make_shared<Task67>(I31, Gamma16_(), I45);
  task58->add_dep(task67);
  task67->add_dep(task31);
  residualq->add_task(task67);

  auto task68 = make_shared<Task68>(I45, t2, v2_, this->e0_);
  task67->add_dep(task68);
  task68->add_dep(task31);
  residualq->add_task(task68);

  auto task69 = make_shared<Task69>(I45, t2, f1_);
  task67->add_dep(task69);
  task69->add_dep(task31);
  residualq->add_task(task69);

  auto task70 = make_shared<Task70>(I45, t2, f1_);
  task67->add_dep(task70);
  task70->add_dep(task31);
  residualq->add_task(task70);

  auto task71 = make_shared<Task71>(I45, t2, f1_);
  task67->add_dep(task71);
  task71->add_dep(task31);
  residualq->add_task(task71);

  auto task72 = make_shared<Task72>(I45, t2, f1_);
  task67->add_dep(task72);
  task72->add_dep(task31);
  residualq->add_task(task72);

  auto task73 = make_shared<Task73>(I45, t2, f1_);
  task67->add_dep(task73);
  task73->add_dep(task31);
  residualq->add_task(task73);

  auto task74 = make_shared<Task74>(I45, t2, f1_);
  task67->add_dep(task74);
  task74->add_dep(task31);
  residualq->add_task(task74);

  auto I63 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task75 = make_shared<Task75>(I31, f1_, I63);
  task58->add_dep(task75);
  task75->add_dep(task31);
  residualq->add_task(task75);

  auto task76 = make_shared<Task76>(I63, Gamma22_(), t2);
  task75->add_dep(task76);
  task76->add_dep(task31);
  residualq->add_task(task76);

  auto task77 = make_shared<Task77>(I63, Gamma12_(), t2);
  task75->add_dep(task77);
  task77->add_dep(task31);
  residualq->add_task(task77);

  auto I66 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task78 = make_shared<Task78>(I31, f1_, I66);
  task58->add_dep(task78);
  task78->add_dep(task31);
  residualq->add_task(task78);

  auto I67 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task79 = make_shared<Task79>(I66, Gamma12_(), I67);
  task78->add_dep(task79);
  task79->add_dep(task31);
  residualq->add_task(task79);

  auto task80 = make_shared<Task80>(I67, t2);
  task79->add_dep(task80);
  task80->add_dep(task31);
  residualq->add_task(task80);

  auto I75 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{virt_, active_});
  auto task81 = make_shared<Task81>(I31, t2, I75);
  task58->add_dep(task81);
  task81->add_dep(task31);
  residualq->add_task(task81);

  auto task82 = make_shared<Task82>(I75, Gamma16_(), f1_);
  task81->add_dep(task82);
  task82->add_dep(task31);
  residualq->add_task(task82);

  auto I78 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{virt_, active_});
  auto task83 = make_shared<Task83>(I31, t2, I78);
  task58->add_dep(task83);
  task83->add_dep(task31);
  residualq->add_task(task83);

  auto task84 = make_shared<Task84>(I78, Gamma16_(), f1_);
  task83->add_dep(task84);
  task84->add_dep(task31);
  residualq->add_task(task84);

  auto I80 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, virt_});
  auto task85 = make_shared<Task85>(r, I80);
  task85->add_dep(task31);
  residualq->add_task(task85);

  auto I81 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task86 = make_shared<Task86>(I80, f1_, I81);
  task85->add_dep(task86);
  task86->add_dep(task31);
  residualq->add_task(task86);

  auto task87 = make_shared<Task87>(I81, Gamma28_(), t2);
  task86->add_dep(task87);
  task87->add_dep(task31);
  residualq->add_task(task87);

  auto I84 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task88 = make_shared<Task88>(I80, Gamma29_(), I84);
  task85->add_dep(task88);
  task88->add_dep(task31);
  residualq->add_task(task88);

  auto task89 = make_shared<Task89>(I84, v2_);
  task88->add_dep(task89);
  task89->add_dep(task31);
  residualq->add_task(task89);

  auto task90 = make_shared<Task90>(I84, t2, f1_);
  task88->add_dep(task90);
  task90->add_dep(task31);
  residualq->add_task(task90);

  auto I87 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, closed_, virt_, active_});
  auto task91 = make_shared<Task91>(I80, Gamma7_(), I87);
  task85->add_dep(task91);
  task91->add_dep(task31);
  residualq->add_task(task91);

  auto task92 = make_shared<Task92>(I87, t2, f1_);
  task91->add_dep(task92);
  task92->add_dep(task31);
  residualq->add_task(task92);

  auto task93 = make_shared<Task93>(I80, Gamma31_(), t2);
  task85->add_dep(task93);
  task93->add_dep(task31);
  residualq->add_task(task93);

  auto I92 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, virt_, active_});
  auto task94 = make_shared<Task94>(I80, Gamma32_(), I92);
  task85->add_dep(task94);
  task94->add_dep(task31);
  residualq->add_task(task94);

  auto task95 = make_shared<Task95>(I92, t2, v2_, this->e0_);
  task94->add_dep(task95);
  task95->add_dep(task31);
  residualq->add_task(task95);

  auto task96 = make_shared<Task96>(I92, t2, f1_);
  task94->add_dep(task96);
  task96->add_dep(task31);
  residualq->add_task(task96);

  auto task97 = make_shared<Task97>(I92, t2, f1_);
  task94->add_dep(task97);
  task97->add_dep(task31);
  residualq->add_task(task97);

  auto task98 = make_shared<Task98>(I92, t2, f1_);
  task94->add_dep(task98);
  task98->add_dep(task31);
  residualq->add_task(task98);

  auto task99 = make_shared<Task99>(I80, Gamma34_(), t2);
  task85->add_dep(task99);
  task99->add_dep(task31);
  residualq->add_task(task99);

  auto I100 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task100 = make_shared<Task100>(I80, Gamma35_(), I100);
  task85->add_dep(task100);
  task100->add_dep(task31);
  residualq->add_task(task100);

  auto task101 = make_shared<Task101>(I100, t2, v2_, this->e0_);
  task100->add_dep(task101);
  task101->add_dep(task31);
  residualq->add_task(task101);

  auto task102 = make_shared<Task102>(I100, t2, f1_);
  task100->add_dep(task102);
  task102->add_dep(task31);
  residualq->add_task(task102);

  auto task103 = make_shared<Task103>(I100, t2, f1_);
  task100->add_dep(task103);
  task103->add_dep(task31);
  residualq->add_task(task103);

  auto I103 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task104 = make_shared<Task104>(I80, f1_, I103);
  task85->add_dep(task104);
  task104->add_dep(task31);
  residualq->add_task(task104);

  auto task105 = make_shared<Task105>(I103, Gamma35_(), t2);
  task104->add_dep(task105);
  task105->add_dep(task31);
  residualq->add_task(task105);

  auto I106 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task106 = make_shared<Task106>(I80, f1_, I106);
  task85->add_dep(task106);
  task106->add_dep(task31);
  residualq->add_task(task106);

  auto task107 = make_shared<Task107>(I106, Gamma37_(), t2);
  task106->add_dep(task107);
  task107->add_dep(task31);
  residualq->add_task(task107);

  auto I109 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{closed_, virt_});
  auto task108 = make_shared<Task108>(I80, Gamma38_(), I109);
  task85->add_dep(task108);
  task108->add_dep(task31);
  residualq->add_task(task108);

  auto task109 = make_shared<Task109>(I109, h1_);
  task108->add_dep(task109);
  task109->add_dep(task31);
  residualq->add_task(task109);

  auto task110 = make_shared<Task110>(I109, t2, f1_);
  task108->add_dep(task110);
  task110->add_dep(task31);
  residualq->add_task(task110);

  auto task111 = make_shared<Task111>(I109, t2, f1_);
  task108->add_dep(task111);
  task111->add_dep(task31);
  residualq->add_task(task111);

  auto I120 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, virt_});
  auto task112 = make_shared<Task112>(r, I120);
  task112->add_dep(task31);
  residualq->add_task(task112);

  auto I121 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task113 = make_shared<Task113>(I120, f1_, I121);
  task112->add_dep(task113);
  task113->add_dep(task31);
  residualq->add_task(task113);

  auto task114 = make_shared<Task114>(I121, Gamma6_(), t2);
  task113->add_dep(task114);
  task114->add_dep(task31);
  residualq->add_task(task114);

  auto I124 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task115 = make_shared<Task115>(I120, Gamma7_(), I124);
  task112->add_dep(task115);
  task115->add_dep(task31);
  residualq->add_task(task115);

  auto task116 = make_shared<Task116>(I124, v2_);
  task115->add_dep(task116);
  task116->add_dep(task31);
  residualq->add_task(task116);

  auto task117 = make_shared<Task117>(I124, t2, f1_);
  task115->add_dep(task117);
  task117->add_dep(task31);
  residualq->add_task(task117);

  auto task118 = make_shared<Task118>(I124, t2, f1_);
  task115->add_dep(task118);
  task118->add_dep(task31);
  residualq->add_task(task118);

  auto I130 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task119 = make_shared<Task119>(I120, Gamma34_(), I130);
  task112->add_dep(task119);
  task119->add_dep(task31);
  residualq->add_task(task119);

  auto task120 = make_shared<Task120>(I130, t2);
  task119->add_dep(task120);
  task120->add_dep(task31);
  residualq->add_task(task120);

  auto I132 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, virt_, active_});
  auto task121 = make_shared<Task121>(I120, Gamma35_(), I132);
  task112->add_dep(task121);
  task121->add_dep(task31);
  residualq->add_task(task121);

  auto task122 = make_shared<Task122>(I132, t2, v2_, this->e0_);
  task121->add_dep(task122);
  task122->add_dep(task31);
  residualq->add_task(task122);

  auto task123 = make_shared<Task123>(I132, t2, f1_);
  task121->add_dep(task123);
  task123->add_dep(task31);
  residualq->add_task(task123);

  auto task124 = make_shared<Task124>(I132, t2, f1_);
  task121->add_dep(task124);
  task124->add_dep(task31);
  residualq->add_task(task124);

  auto task125 = make_shared<Task125>(I132, t2, f1_);
  task121->add_dep(task125);
  task125->add_dep(task31);
  residualq->add_task(task125);

  auto task126 = make_shared<Task126>(I132, t2, f1_);
  task121->add_dep(task126);
  task126->add_dep(task31);
  residualq->add_task(task126);

  auto task127 = make_shared<Task127>(I132, t2, f1_);
  task121->add_dep(task127);
  task127->add_dep(task31);
  residualq->add_task(task127);

  auto task128 = make_shared<Task128>(I132, t2, f1_);
  task121->add_dep(task128);
  task128->add_dep(task31);
  residualq->add_task(task128);

  auto I146 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task129 = make_shared<Task129>(I120, f1_, I146);
  task112->add_dep(task129);
  task129->add_dep(task31);
  residualq->add_task(task129);

  auto task130 = make_shared<Task130>(I146, Gamma51_(), t2);
  task129->add_dep(task130);
  task130->add_dep(task31);
  residualq->add_task(task130);

  auto I149 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{closed_, virt_});
  auto task131 = make_shared<Task131>(I120, Gamma38_(), I149);
  task112->add_dep(task131);
  task131->add_dep(task31);
  residualq->add_task(task131);

  auto task132 = make_shared<Task132>(I149, h1_);
  task131->add_dep(task132);
  task132->add_dep(task31);
  residualq->add_task(task132);

  auto task133 = make_shared<Task133>(I149, t2, f1_);
  task131->add_dep(task133);
  task133->add_dep(task31);
  residualq->add_task(task133);

  auto task134 = make_shared<Task134>(I149, t2, f1_);
  task131->add_dep(task134);
  task134->add_dep(task31);
  residualq->add_task(task134);

  auto I160 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task135 = make_shared<Task135>(r, I160);
  task135->add_dep(task31);
  residualq->add_task(task135);

  auto I161 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, active_, virt_, active_});
  auto task136 = make_shared<Task136>(I160, Gamma56_(), I161);
  task135->add_dep(task136);
  task136->add_dep(task31);
  residualq->add_task(task136);

  auto task137 = make_shared<Task137>(I161, t2, f1_);
  task136->add_dep(task137);
  task137->add_dep(task31);
  residualq->add_task(task137);

  auto I164 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, active_, active_});
  auto task138 = make_shared<Task138>(I160, Gamma57_(), I164);
  task135->add_dep(task138);
  task138->add_dep(task31);
  residualq->add_task(task138);

  auto task139 = make_shared<Task139>(I164, v2_);
  task138->add_dep(task139);
  task139->add_dep(task31);
  residualq->add_task(task139);

  auto task140 = make_shared<Task140>(I164, t2, f1_);
  task138->add_dep(task140);
  task140->add_dep(task31);
  residualq->add_task(task140);

  auto task141 = make_shared<Task141>(I160, Gamma58_(), t2);
  task135->add_dep(task141);
  task141->add_dep(task31);
  residualq->add_task(task141);

  auto I169 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task142 = make_shared<Task142>(I160, Gamma59_(), I169);
  task135->add_dep(task142);
  task142->add_dep(task31);
  residualq->add_task(task142);

  auto task143 = make_shared<Task143>(I169, t2, v2_, this->e0_);
  task142->add_dep(task143);
  task143->add_dep(task31);
  residualq->add_task(task143);

  auto task144 = make_shared<Task144>(I169, t2, f1_);
  task142->add_dep(task144);
  task144->add_dep(task31);
  residualq->add_task(task144);

  auto task145 = make_shared<Task145>(I169, t2, f1_);
  task142->add_dep(task145);
  task145->add_dep(task31);
  residualq->add_task(task145);

  auto I172 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{active_, virt_});
  auto task146 = make_shared<Task146>(I160, Gamma60_(), I172);
  task135->add_dep(task146);
  task146->add_dep(task31);
  residualq->add_task(task146);

  auto task147 = make_shared<Task147>(I172, h1_);
  task146->add_dep(task147);
  task147->add_dep(task31);
  residualq->add_task(task147);

  auto task148 = make_shared<Task148>(I172, t2, f1_);
  task146->add_dep(task148);
  task148->add_dep(task31);
  residualq->add_task(task148);

  auto task149 = make_shared<Task149>(I172, t2, f1_);
  task146->add_dep(task149);
  task149->add_dep(task31);
  residualq->add_task(task149);

  auto I180 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, closed_, virt_, closed_});
  auto task150 = make_shared<Task150>(r, I180);
  task150->add_dep(task31);
  residualq->add_task(task150);

  auto I181 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{virt_, active_});
  auto task151 = make_shared<Task151>(I180, t2, I181);
  task150->add_dep(task151);
  task151->add_dep(task31);
  residualq->add_task(task151);

  auto task152 = make_shared<Task152>(I181, Gamma16_(), f1_);
  task151->add_dep(task152);
  task152->add_dep(task31);
  residualq->add_task(task152);

  auto I184 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{virt_, active_});
  auto task153 = make_shared<Task153>(I180, t2, I184);
  task150->add_dep(task153);
  task153->add_dep(task31);
  residualq->add_task(task153);

  auto task154 = make_shared<Task154>(I184, Gamma16_(), f1_);
  task153->add_dep(task154);
  task154->add_dep(task31);
  residualq->add_task(task154);

  auto I187 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task155 = make_shared<Task155>(I180, f1_, I187);
  task150->add_dep(task155);
  task155->add_dep(task31);
  residualq->add_task(task155);

  auto I188 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task156 = make_shared<Task156>(I187, Gamma38_(), I188);
  task155->add_dep(task156);
  task156->add_dep(task31);
  residualq->add_task(task156);

  auto task157 = make_shared<Task157>(I188, t2);
  task156->add_dep(task157);
  task157->add_dep(task31);
  residualq->add_task(task157);

  auto I190 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task158 = make_shared<Task158>(I180, f1_, I190);
  task150->add_dep(task158);
  task158->add_dep(task31);
  residualq->add_task(task158);

  auto I191 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task159 = make_shared<Task159>(I190, Gamma38_(), I191);
  task158->add_dep(task159);
  task159->add_dep(task31);
  residualq->add_task(task159);

  auto task160 = make_shared<Task160>(I191, t2);
  task159->add_dep(task160);
  task160->add_dep(task31);
  residualq->add_task(task160);

  shared_ptr<Task161> task161;
  if (diagonal) {
    task161 = make_shared<Task161>(I180, t2, f1_);
    task150->add_dep(task161);
    task161->add_dep(task31);
    residualq->add_task(task161);
  }

  shared_ptr<Task162> task162;
  if (diagonal) {
    task162 = make_shared<Task162>(I180, t2, f1_);
    task150->add_dep(task162);
    task162->add_dep(task31);
    residualq->add_task(task162);
  }

  shared_ptr<Task163> task163;
  if (diagonal) {
    task163 = make_shared<Task163>(I180, t2, f1_);
    task150->add_dep(task163);
    task163->add_dep(task31);
    residualq->add_task(task163);
  }

  shared_ptr<Task164> task164;
  if (diagonal) {
    task164 = make_shared<Task164>(I180, t2, f1_);
    task150->add_dep(task164);
    task164->add_dep(task31);
    residualq->add_task(task164);
  }

  auto I211 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{closed_, active_});
  auto task165 = make_shared<Task165>(I180, t2, I211);
  task150->add_dep(task165);
  task165->add_dep(task31);
  residualq->add_task(task165);

  auto task166 = make_shared<Task166>(I211, Gamma38_(), f1_);
  task165->add_dep(task166);
  task166->add_dep(task31);
  residualq->add_task(task166);

  auto I214 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{closed_, active_});
  auto task167 = make_shared<Task167>(I180, t2, I214);
  task150->add_dep(task167);
  task167->add_dep(task31);
  residualq->add_task(task167);

  auto task168 = make_shared<Task168>(I214, Gamma38_(), f1_);
  task167->add_dep(task168);
  task168->add_dep(task31);
  residualq->add_task(task168);

  auto I198 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, virt_});
  auto task169 = make_shared<Task169>(r, I198);
  task169->add_dep(task31);
  residualq->add_task(task169);

  shared_ptr<Task170> task170;
  if (diagonal) {
    task170 = make_shared<Task170>(I198, t2, v2_, this->e0_);
    task169->add_dep(task170);
    task170->add_dep(task31);
    residualq->add_task(task170);
  }

  auto I199 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, virt_, closed_, virt_});
  auto task171 = make_shared<Task171>(I198, Gamma69_(), I199);
  task169->add_dep(task171);
  task171->add_dep(task31);
  residualq->add_task(task171);

  auto task172 = make_shared<Task172>(I199, t2);
  task171->add_dep(task172);
  task172->add_dep(task31);
  residualq->add_task(task172);

  auto I216 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, closed_, active_, virt_});
  auto task173 = make_shared<Task173>(r, I216);
  task173->add_dep(task31);
  residualq->add_task(task173);

  auto I217 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task174 = make_shared<Task174>(I216, f1_, I217);
  task173->add_dep(task174);
  task174->add_dep(task31);
  residualq->add_task(task174);

  auto I218 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task175 = make_shared<Task175>(I217, Gamma35_(), I218);
  task174->add_dep(task175);
  task175->add_dep(task31);
  residualq->add_task(task175);

  auto task176 = make_shared<Task176>(I218, t2);
  task175->add_dep(task176);
  task176->add_dep(task31);
  residualq->add_task(task176);

  auto I220 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task177 = make_shared<Task177>(I216, f1_, I220);
  task173->add_dep(task177);
  task177->add_dep(task31);
  residualq->add_task(task177);

  auto task178 = make_shared<Task178>(I220, Gamma32_(), t2);
  task177->add_dep(task178);
  task178->add_dep(task31);
  residualq->add_task(task178);

  auto task179 = make_shared<Task179>(I220, Gamma35_(), t2);
  task177->add_dep(task179);
  task179->add_dep(task31);
  residualq->add_task(task179);

  auto I229 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{virt_, active_});
  auto task180 = make_shared<Task180>(I216, f1_, I229);
  task173->add_dep(task180);
  task180->add_dep(task31);
  residualq->add_task(task180);

  auto task181 = make_shared<Task181>(I229, Gamma60_(), t2);
  task180->add_dep(task181);
  task181->add_dep(task31);
  residualq->add_task(task181);

  auto I232 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{virt_, active_});
  auto task182 = make_shared<Task182>(I216, f1_, I232);
  task173->add_dep(task182);
  task182->add_dep(task31);
  residualq->add_task(task182);

  auto task183 = make_shared<Task183>(I232, Gamma60_(), t2);
  task182->add_dep(task183);
  task183->add_dep(task31);
  residualq->add_task(task183);

  auto I235 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{closed_, active_});
  auto task184 = make_shared<Task184>(I216, t2, I235);
  task173->add_dep(task184);
  task184->add_dep(task31);
  residualq->add_task(task184);

  auto task185 = make_shared<Task185>(I235, Gamma38_(), f1_);
  task184->add_dep(task185);
  task185->add_dep(task31);
  residualq->add_task(task185);

  auto I238 = make_shared<TATensor<std::complex<double>,2>>(std::vector<IndexRange>{closed_, active_});
  auto task186 = make_shared<Task186>(I216, t2, I238);
  task173->add_dep(task186);
  task186->add_dep(task31);
  residualq->add_task(task186);

  auto task187 = make_shared<Task187>(I238, Gamma38_(), f1_);
  task186->add_dep(task187);
  task187->add_dep(task31);
  residualq->add_task(task187);

  auto I241 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, closed_, virt_});
  auto task188 = make_shared<Task188>(I216, Gamma81_(), I241);
  task173->add_dep(task188);
  task188->add_dep(task31);
  residualq->add_task(task188);

  auto task189 = make_shared<Task189>(I241, t2);
  task188->add_dep(task189);
  task189->add_dep(task31);
  residualq->add_task(task189);

  auto I245 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{closed_, active_, virt_, virt_});
  auto task190 = make_shared<Task190>(I216, Gamma38_(), I245);
  task173->add_dep(task190);
  task190->add_dep(task31);
  residualq->add_task(task190);

  auto task191 = make_shared<Task191>(I245, t2, v2_, this->e0_);
  task190->add_dep(task191);
  task191->add_dep(task31);
  residualq->add_task(task191);

  auto task192 = make_shared<Task192>(I245, t2, f1_);
  task190->add_dep(task192);
  task192->add_dep(task31);
  residualq->add_task(task192);

  auto task193 = make_shared<Task193>(I245, t2, f1_);
  task190->add_dep(task193);
  task193->add_dep(task31);
  residualq->add_task(task193);

  auto task194 = make_shared<Task194>(I245, t2, f1_);
  task190->add_dep(task194);
  task194->add_dep(task31);
  residualq->add_task(task194);

  auto task195 = make_shared<Task195>(I245, t2, f1_);
  task190->add_dep(task195);
  task195->add_dep(task31);
  residualq->add_task(task195);

  auto task196 = make_shared<Task196>(I245, t2, f1_);
  task190->add_dep(task196);
  task196->add_dep(task31);
  residualq->add_task(task196);

  auto task197 = make_shared<Task197>(I245, t2, f1_);
  task190->add_dep(task197);
  task197->add_dep(task31);
  residualq->add_task(task197);

  auto I263 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task198 = make_shared<Task198>(I216, f1_, I263);
  task173->add_dep(task198);
  task198->add_dep(task31);
  residualq->add_task(task198);

  auto task199 = make_shared<Task199>(I263, Gamma60_(), t2);
  task198->add_dep(task199);
  task199->add_dep(task31);
  residualq->add_task(task199);

  auto I265 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, active_, active_, virt_});
  auto task200 = make_shared<Task200>(r, I265);
  task200->add_dep(task31);
  residualq->add_task(task200);

  auto I266 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task201 = make_shared<Task201>(I265, f1_, I266);
  task200->add_dep(task201);
  task201->add_dep(task31);
  residualq->add_task(task201);

  auto task202 = make_shared<Task202>(I266, Gamma59_(), t2);
  task201->add_dep(task202);
  task202->add_dep(task31);
  residualq->add_task(task202);

  auto I269 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, active_, virt_, virt_});
  auto task203 = make_shared<Task203>(I265, Gamma60_(), I269);
  task200->add_dep(task203);
  task203->add_dep(task31);
  residualq->add_task(task203);

  auto task204 = make_shared<Task204>(I269, t2, f1_);
  task203->add_dep(task204);
  task204->add_dep(task31);
  residualq->add_task(task204);

  auto task205 = make_shared<Task205>(I269, t2, f1_);
  task203->add_dep(task205);
  task205->add_dep(task31);
  residualq->add_task(task205);

  auto I271 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task206 = make_shared<Task206>(r, I271);
  task206->add_dep(task31);
  residualq->add_task(task206);

  auto task207 = make_shared<Task207>(I271, Gamma92_(), t2);
  task206->add_dep(task207);
  task207->add_dep(task31);
  residualq->add_task(task207);

  auto I301 = make_shared<TATensor<std::complex<double>,4>>(std::vector<IndexRange>{active_, virt_, active_, virt_});
  auto task208 = make_shared<Task208>(I271, Gamma60_(), I301);
  task206->add_dep(task208);
  task208->add_dep(task31);
  residualq->add_task(task208);

  auto task209 = make_shared<Task209>(I301, t2, v2_, this->e0_);
  task208->add_dep(task209);
  task209->add_dep(task31);
  residualq->add_task(task209);

  return residualq;
}


#endif
