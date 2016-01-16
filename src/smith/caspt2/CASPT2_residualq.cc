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


#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_residualq(const bool reset, const bool diagonal) {

  auto residualq = make_shared<Queue>();
  auto task69 = make_shared<Task69>(r, reset);
  residualq->add_task(task69);

  auto I0 = make_shared<TATensor<double,4>>({closed_, closed_, active_, active_});
  auto task70 = make_shared<Task70>(r, I0);
  task70->add_dep(task69);
  residualq->add_task(task70);

  auto task71 = make_shared<Task71>(I0, Gamma0_(), t2);
  task70->add_dep(task71);
  task71->add_dep(task69);
  residualq->add_task(task71);

  auto task72 = make_shared<Task72>(I0, Gamma92_(), t2, this->e0_);
  task70->add_dep(task72);
  task72->add_dep(task69);
  residualq->add_task(task72);

  auto I2 = make_shared<TATensor<double,4>>({closed_, closed_, active_, active_});
  auto task73 = make_shared<Task73>(r, I2);
  task73->add_dep(task69);
  residualq->add_task(task73);

  auto I3 = make_shared<TATensor<double,4>>({closed_, closed_, active_, active_});
  auto task74 = make_shared<Task74>(I2, Gamma92_(), I3);
  task73->add_dep(task74);
  task74->add_dep(task69);
  residualq->add_task(task74);

  auto task75 = make_shared<Task75>(I3, t2, f1_);
  task74->add_dep(task75);
  task75->add_dep(task69);
  residualq->add_task(task75);

  auto I6 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task76 = make_shared<Task76>(I2, f1_, I6);
  task73->add_dep(task76);
  task76->add_dep(task69);
  residualq->add_task(task76);

  auto task77 = make_shared<Task77>(I6, Gamma2_(), t2);
  task76->add_dep(task77);
  task77->add_dep(task69);
  residualq->add_task(task77);

  auto I9 = make_shared<TATensor<double,4>>({active_, closed_, closed_, active_});
  auto task78 = make_shared<Task78>(I2, Gamma3_(), I9);
  task73->add_dep(task78);
  task78->add_dep(task69);
  residualq->add_task(task78);

  auto task79 = make_shared<Task79>(I9, t2, f1_);
  task78->add_dep(task79);
  task79->add_dep(task69);
  residualq->add_task(task79);

  auto I11 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task80 = make_shared<Task80>(r, I11);
  task80->add_dep(task69);
  residualq->add_task(task80);

  auto I12 = make_shared<TATensor<double,4>>({active_, closed_, active_, active_});
  auto task81 = make_shared<Task81>(I11, Gamma4_(), I12);
  task80->add_dep(task81);
  task81->add_dep(task69);
  residualq->add_task(task81);

  auto task82 = make_shared<Task82>(I12, t2, f1_);
  task81->add_dep(task82);
  task82->add_dep(task69);
  residualq->add_task(task82);

  auto task83 = make_shared<Task83>(I11, Gamma5_(), t2);
  task80->add_dep(task83);
  task83->add_dep(task69);
  residualq->add_task(task83);

  auto I17 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task84 = make_shared<Task84>(I11, Gamma6_(), I17);
  task80->add_dep(task84);
  task84->add_dep(task69);
  residualq->add_task(task84);

  auto task85 = make_shared<Task85>(I17, t2, this->e0_);
  task84->add_dep(task85);
  task85->add_dep(task69);
  residualq->add_task(task85);

  auto task86 = make_shared<Task86>(I17, t2, f1_);
  task84->add_dep(task86);
  task86->add_dep(task69);
  residualq->add_task(task86);

  auto task87 = make_shared<Task87>(I17, t2, f1_);
  task84->add_dep(task87);
  task87->add_dep(task69);
  residualq->add_task(task87);

  auto I20 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task88 = make_shared<Task88>(I11, Gamma7_(), I20);
  task80->add_dep(task88);
  task88->add_dep(task69);
  residualq->add_task(task88);

  auto task89 = make_shared<Task89>(I20, t2, f1_);
  task88->add_dep(task89);
  task89->add_dep(task69);
  residualq->add_task(task89);

  auto task90 = make_shared<Task90>(I20, t2, f1_);
  task88->add_dep(task90);
  task90->add_dep(task69);
  residualq->add_task(task90);

  auto I26 = make_shared<TATensor<double,4>>({active_, active_, closed_, active_});
  auto task91 = make_shared<Task91>(I11, Gamma9_(), I26);
  task80->add_dep(task91);
  task91->add_dep(task69);
  residualq->add_task(task91);

  auto task92 = make_shared<Task92>(I26, t2, f1_);
  task91->add_dep(task92);
  task92->add_dep(task69);
  residualq->add_task(task92);

  auto I31 = make_shared<TATensor<double,4>>({closed_, closed_, active_, virt_});
  auto task93 = make_shared<Task93>(r, I31);
  task93->add_dep(task69);
  residualq->add_task(task93);

  auto I32 = make_shared<TATensor<double,4>>({closed_, closed_, active_, active_});
  auto task94 = make_shared<Task94>(I31, f1_, I32);
  task93->add_dep(task94);
  task94->add_dep(task69);
  residualq->add_task(task94);

  auto task95 = make_shared<Task95>(I32, Gamma3_(), t2);
  task94->add_dep(task95);
  task95->add_dep(task69);
  residualq->add_task(task95);

  auto I35 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task96 = make_shared<Task96>(I31, f1_, I35);
  task93->add_dep(task96);
  task96->add_dep(task69);
  residualq->add_task(task96);

  auto task97 = make_shared<Task97>(I35, Gamma12_(), t2);
  task96->add_dep(task97);
  task97->add_dep(task69);
  residualq->add_task(task97);

  auto I38 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task98 = make_shared<Task98>(I31, f1_, I38);
  task93->add_dep(task98);
  task98->add_dep(task69);
  residualq->add_task(task98);

  auto task99 = make_shared<Task99>(I38, Gamma12_(), t2);
  task98->add_dep(task99);
  task99->add_dep(task69);
  residualq->add_task(task99);

  auto I41 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, active_});
  auto task100 = make_shared<Task100>(I31, Gamma14_(), I41);
  task93->add_dep(task100);
  task100->add_dep(task69);
  residualq->add_task(task100);

  auto task101 = make_shared<Task101>(I41, t2);
  task100->add_dep(task101);
  task101->add_dep(task69);
  residualq->add_task(task101);

  auto I45 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, active_});
  auto task102 = make_shared<Task102>(I31, Gamma16_(), I45);
  task93->add_dep(task102);
  task102->add_dep(task69);
  residualq->add_task(task102);

  auto task103 = make_shared<Task103>(I45, t2, this->e0_);
  task102->add_dep(task103);
  task103->add_dep(task69);
  residualq->add_task(task103);

  auto task104 = make_shared<Task104>(I45, t2, f1_);
  task102->add_dep(task104);
  task104->add_dep(task69);
  residualq->add_task(task104);

  auto task105 = make_shared<Task105>(I45, t2, f1_);
  task102->add_dep(task105);
  task105->add_dep(task69);
  residualq->add_task(task105);

  auto task106 = make_shared<Task106>(I45, t2, f1_);
  task102->add_dep(task106);
  task106->add_dep(task69);
  residualq->add_task(task106);

  auto task107 = make_shared<Task107>(I45, t2, f1_);
  task102->add_dep(task107);
  task107->add_dep(task69);
  residualq->add_task(task107);

  auto task108 = make_shared<Task108>(I45, t2, f1_);
  task102->add_dep(task108);
  task108->add_dep(task69);
  residualq->add_task(task108);

  auto task109 = make_shared<Task109>(I45, t2, f1_);
  task102->add_dep(task109);
  task109->add_dep(task69);
  residualq->add_task(task109);

  auto I63 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task110 = make_shared<Task110>(I31, f1_, I63);
  task93->add_dep(task110);
  task110->add_dep(task69);
  residualq->add_task(task110);

  auto task111 = make_shared<Task111>(I63, Gamma22_(), t2);
  task110->add_dep(task111);
  task111->add_dep(task69);
  residualq->add_task(task111);

  auto task112 = make_shared<Task112>(I63, Gamma12_(), t2);
  task110->add_dep(task112);
  task112->add_dep(task69);
  residualq->add_task(task112);

  auto I66 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task113 = make_shared<Task113>(I31, f1_, I66);
  task93->add_dep(task113);
  task113->add_dep(task69);
  residualq->add_task(task113);

  auto I67 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task114 = make_shared<Task114>(I66, Gamma12_(), I67);
  task113->add_dep(task114);
  task114->add_dep(task69);
  residualq->add_task(task114);

  auto task115 = make_shared<Task115>(I67, t2);
  task114->add_dep(task115);
  task115->add_dep(task69);
  residualq->add_task(task115);

  auto I75 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task116 = make_shared<Task116>(I31, t2, I75);
  task93->add_dep(task116);
  task116->add_dep(task69);
  residualq->add_task(task116);

  auto task117 = make_shared<Task117>(I75, Gamma16_(), f1_);
  task116->add_dep(task117);
  task117->add_dep(task69);
  residualq->add_task(task117);

  auto I78 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task118 = make_shared<Task118>(I31, t2, I78);
  task93->add_dep(task118);
  task118->add_dep(task69);
  residualq->add_task(task118);

  auto task119 = make_shared<Task119>(I78, Gamma16_(), f1_);
  task118->add_dep(task119);
  task119->add_dep(task69);
  residualq->add_task(task119);

  auto I80 = make_shared<TATensor<double,4>>({closed_, active_, active_, virt_});
  auto task120 = make_shared<Task120>(r, I80);
  task120->add_dep(task69);
  residualq->add_task(task120);

  auto I81 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task121 = make_shared<Task121>(I80, f1_, I81);
  task120->add_dep(task121);
  task121->add_dep(task69);
  residualq->add_task(task121);

  auto task122 = make_shared<Task122>(I81, Gamma28_(), t2);
  task121->add_dep(task122);
  task122->add_dep(task69);
  residualq->add_task(task122);

  auto I84 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task123 = make_shared<Task123>(I80, Gamma29_(), I84);
  task120->add_dep(task123);
  task123->add_dep(task69);
  residualq->add_task(task123);

  auto task124 = make_shared<Task124>(I84, t2, f1_);
  task123->add_dep(task124);
  task124->add_dep(task69);
  residualq->add_task(task124);

  auto I87 = make_shared<TATensor<double,4>>({active_, closed_, virt_, active_});
  auto task125 = make_shared<Task125>(I80, Gamma7_(), I87);
  task120->add_dep(task125);
  task125->add_dep(task69);
  residualq->add_task(task125);

  auto task126 = make_shared<Task126>(I87, t2, f1_);
  task125->add_dep(task126);
  task126->add_dep(task69);
  residualq->add_task(task126);

  auto task127 = make_shared<Task127>(I80, Gamma31_(), t2);
  task120->add_dep(task127);
  task127->add_dep(task69);
  residualq->add_task(task127);

  auto I92 = make_shared<TATensor<double,4>>({closed_, active_, virt_, active_});
  auto task128 = make_shared<Task128>(I80, Gamma32_(), I92);
  task120->add_dep(task128);
  task128->add_dep(task69);
  residualq->add_task(task128);

  auto task129 = make_shared<Task129>(I92, t2, this->e0_);
  task128->add_dep(task129);
  task129->add_dep(task69);
  residualq->add_task(task129);

  auto task130 = make_shared<Task130>(I92, t2, f1_);
  task128->add_dep(task130);
  task130->add_dep(task69);
  residualq->add_task(task130);

  auto task131 = make_shared<Task131>(I92, t2, f1_);
  task128->add_dep(task131);
  task131->add_dep(task69);
  residualq->add_task(task131);

  auto task132 = make_shared<Task132>(I92, t2, f1_);
  task128->add_dep(task132);
  task132->add_dep(task69);
  residualq->add_task(task132);

  auto task133 = make_shared<Task133>(I80, Gamma34_(), t2);
  task120->add_dep(task133);
  task133->add_dep(task69);
  residualq->add_task(task133);

  auto I100 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task134 = make_shared<Task134>(I80, Gamma35_(), I100);
  task120->add_dep(task134);
  task134->add_dep(task69);
  residualq->add_task(task134);

  auto task135 = make_shared<Task135>(I100, t2, this->e0_);
  task134->add_dep(task135);
  task135->add_dep(task69);
  residualq->add_task(task135);

  auto task136 = make_shared<Task136>(I100, t2, f1_);
  task134->add_dep(task136);
  task136->add_dep(task69);
  residualq->add_task(task136);

  auto task137 = make_shared<Task137>(I100, t2, f1_);
  task134->add_dep(task137);
  task137->add_dep(task69);
  residualq->add_task(task137);

  auto task138 = make_shared<Task138>(I100, t2, f1_);
  task134->add_dep(task138);
  task138->add_dep(task69);
  residualq->add_task(task138);

  auto I106 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task139 = make_shared<Task139>(I80, f1_, I106);
  task120->add_dep(task139);
  task139->add_dep(task69);
  residualq->add_task(task139);

  auto task140 = make_shared<Task140>(I106, Gamma37_(), t2);
  task139->add_dep(task140);
  task140->add_dep(task69);
  residualq->add_task(task140);

  auto I109 = make_shared<TATensor<double,2>>({closed_, virt_});
  auto task141 = make_shared<Task141>(I80, Gamma38_(), I109);
  task120->add_dep(task141);
  task141->add_dep(task69);
  residualq->add_task(task141);

  auto task142 = make_shared<Task142>(I109, t2, f1_);
  task141->add_dep(task142);
  task142->add_dep(task69);
  residualq->add_task(task142);

  auto task143 = make_shared<Task143>(I109, t2, f1_);
  task141->add_dep(task143);
  task143->add_dep(task69);
  residualq->add_task(task143);

  auto I120 = make_shared<TATensor<double,4>>({closed_, active_, active_, virt_});
  auto task144 = make_shared<Task144>(r, I120);
  task144->add_dep(task69);
  residualq->add_task(task144);

  auto I121 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task145 = make_shared<Task145>(I120, f1_, I121);
  task144->add_dep(task145);
  task145->add_dep(task69);
  residualq->add_task(task145);

  auto task146 = make_shared<Task146>(I121, Gamma6_(), t2);
  task145->add_dep(task146);
  task146->add_dep(task69);
  residualq->add_task(task146);

  auto I124 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task147 = make_shared<Task147>(I120, Gamma7_(), I124);
  task144->add_dep(task147);
  task147->add_dep(task69);
  residualq->add_task(task147);

  auto task148 = make_shared<Task148>(I124, t2, f1_);
  task147->add_dep(task148);
  task148->add_dep(task69);
  residualq->add_task(task148);

  auto task149 = make_shared<Task149>(I124, t2, f1_);
  task147->add_dep(task149);
  task149->add_dep(task69);
  residualq->add_task(task149);

  auto I130 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task150 = make_shared<Task150>(I120, Gamma34_(), I130);
  task144->add_dep(task150);
  task150->add_dep(task69);
  residualq->add_task(task150);

  auto task151 = make_shared<Task151>(I130, t2);
  task150->add_dep(task151);
  task151->add_dep(task69);
  residualq->add_task(task151);

  auto I132 = make_shared<TATensor<double,4>>({closed_, active_, virt_, active_});
  auto task152 = make_shared<Task152>(I120, Gamma35_(), I132);
  task144->add_dep(task152);
  task152->add_dep(task69);
  residualq->add_task(task152);

  auto task153 = make_shared<Task153>(I132, t2, this->e0_);
  task152->add_dep(task153);
  task153->add_dep(task69);
  residualq->add_task(task153);

  auto task154 = make_shared<Task154>(I132, t2, f1_);
  task152->add_dep(task154);
  task154->add_dep(task69);
  residualq->add_task(task154);

  auto task155 = make_shared<Task155>(I132, t2, f1_);
  task152->add_dep(task155);
  task155->add_dep(task69);
  residualq->add_task(task155);

  auto task156 = make_shared<Task156>(I132, t2, f1_);
  task152->add_dep(task156);
  task156->add_dep(task69);
  residualq->add_task(task156);

  auto task157 = make_shared<Task157>(I132, t2, f1_);
  task152->add_dep(task157);
  task157->add_dep(task69);
  residualq->add_task(task157);

  auto task158 = make_shared<Task158>(I132, t2, f1_);
  task152->add_dep(task158);
  task158->add_dep(task69);
  residualq->add_task(task158);

  auto task159 = make_shared<Task159>(I132, t2, f1_);
  task152->add_dep(task159);
  task159->add_dep(task69);
  residualq->add_task(task159);

  auto I146 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task160 = make_shared<Task160>(I120, f1_, I146);
  task144->add_dep(task160);
  task160->add_dep(task69);
  residualq->add_task(task160);

  auto task161 = make_shared<Task161>(I146, Gamma51_(), t2);
  task160->add_dep(task161);
  task161->add_dep(task69);
  residualq->add_task(task161);

  auto I149 = make_shared<TATensor<double,2>>({closed_, virt_});
  auto task162 = make_shared<Task162>(I120, Gamma38_(), I149);
  task144->add_dep(task162);
  task162->add_dep(task69);
  residualq->add_task(task162);

  auto task163 = make_shared<Task163>(I149, t2, f1_);
  task162->add_dep(task163);
  task163->add_dep(task69);
  residualq->add_task(task163);

  auto task164 = make_shared<Task164>(I149, t2, f1_);
  task162->add_dep(task164);
  task164->add_dep(task69);
  residualq->add_task(task164);

  auto I160 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task165 = make_shared<Task165>(r, I160);
  task165->add_dep(task69);
  residualq->add_task(task165);

  auto I161 = make_shared<TATensor<double,4>>({active_, active_, virt_, active_});
  auto task166 = make_shared<Task166>(I160, Gamma56_(), I161);
  task165->add_dep(task166);
  task166->add_dep(task69);
  residualq->add_task(task166);

  auto task167 = make_shared<Task167>(I161, t2, f1_);
  task166->add_dep(task167);
  task167->add_dep(task69);
  residualq->add_task(task167);

  auto I164 = make_shared<TATensor<double,4>>({active_, virt_, active_, active_});
  auto task168 = make_shared<Task168>(I160, Gamma57_(), I164);
  task165->add_dep(task168);
  task168->add_dep(task69);
  residualq->add_task(task168);

  auto task169 = make_shared<Task169>(I164, t2, f1_);
  task168->add_dep(task169);
  task169->add_dep(task69);
  residualq->add_task(task169);

  auto task170 = make_shared<Task170>(I160, Gamma58_(), t2);
  task165->add_dep(task170);
  task170->add_dep(task69);
  residualq->add_task(task170);

  auto I169 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task171 = make_shared<Task171>(I160, Gamma59_(), I169);
  task165->add_dep(task171);
  task171->add_dep(task69);
  residualq->add_task(task171);

  auto task172 = make_shared<Task172>(I169, t2, this->e0_);
  task171->add_dep(task172);
  task172->add_dep(task69);
  residualq->add_task(task172);

  auto task173 = make_shared<Task173>(I169, t2, f1_);
  task171->add_dep(task173);
  task173->add_dep(task69);
  residualq->add_task(task173);

  auto task174 = make_shared<Task174>(I169, t2, f1_);
  task171->add_dep(task174);
  task174->add_dep(task69);
  residualq->add_task(task174);

  auto I172 = make_shared<TATensor<double,2>>({active_, virt_});
  auto task175 = make_shared<Task175>(I160, Gamma60_(), I172);
  task165->add_dep(task175);
  task175->add_dep(task69);
  residualq->add_task(task175);

  auto task176 = make_shared<Task176>(I172, t2, f1_);
  task175->add_dep(task176);
  task176->add_dep(task69);
  residualq->add_task(task176);

  auto task177 = make_shared<Task177>(I172, t2, f1_);
  task175->add_dep(task177);
  task177->add_dep(task69);
  residualq->add_task(task177);

  auto I180 = make_shared<TATensor<double,4>>({virt_, closed_, virt_, closed_});
  auto task178 = make_shared<Task178>(r, I180);
  task178->add_dep(task69);
  residualq->add_task(task178);

  auto I181 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task179 = make_shared<Task179>(I180, t2, I181);
  task178->add_dep(task179);
  task179->add_dep(task69);
  residualq->add_task(task179);

  auto task180 = make_shared<Task180>(I181, Gamma16_(), f1_);
  task179->add_dep(task180);
  task180->add_dep(task69);
  residualq->add_task(task180);

  auto I184 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task181 = make_shared<Task181>(I180, t2, I184);
  task178->add_dep(task181);
  task181->add_dep(task69);
  residualq->add_task(task181);

  auto task182 = make_shared<Task182>(I184, Gamma16_(), f1_);
  task181->add_dep(task182);
  task182->add_dep(task69);
  residualq->add_task(task182);

  auto I187 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task183 = make_shared<Task183>(I180, f1_, I187);
  task178->add_dep(task183);
  task183->add_dep(task69);
  residualq->add_task(task183);

  auto I188 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task184 = make_shared<Task184>(I187, Gamma38_(), I188);
  task183->add_dep(task184);
  task184->add_dep(task69);
  residualq->add_task(task184);

  auto task185 = make_shared<Task185>(I188, t2);
  task184->add_dep(task185);
  task185->add_dep(task69);
  residualq->add_task(task185);

  auto I190 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task186 = make_shared<Task186>(I180, f1_, I190);
  task178->add_dep(task186);
  task186->add_dep(task69);
  residualq->add_task(task186);

  auto I191 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task187 = make_shared<Task187>(I190, Gamma38_(), I191);
  task186->add_dep(task187);
  task187->add_dep(task69);
  residualq->add_task(task187);

  auto task188 = make_shared<Task188>(I191, t2);
  task187->add_dep(task188);
  task188->add_dep(task69);
  residualq->add_task(task188);

  auto I199 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task189 = make_shared<Task189>(I180, t2, I199);
  task178->add_dep(task189);
  task189->add_dep(task69);
  residualq->add_task(task189);

  auto task190 = make_shared<Task190>(I199, Gamma38_(), f1_);
  task189->add_dep(task190);
  task190->add_dep(task69);
  residualq->add_task(task190);

  auto I202 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task191 = make_shared<Task191>(I180, t2, I202);
  task178->add_dep(task191);
  task191->add_dep(task69);
  residualq->add_task(task191);

  auto task192 = make_shared<Task192>(I202, Gamma38_(), f1_);
  task191->add_dep(task192);
  task192->add_dep(task69);
  residualq->add_task(task192);

  auto I204 = make_shared<TATensor<double,4>>({virt_, closed_, active_, virt_});
  auto task193 = make_shared<Task193>(r, I204);
  task193->add_dep(task69);
  residualq->add_task(task193);

  auto I205 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task194 = make_shared<Task194>(I204, f1_, I205);
  task193->add_dep(task194);
  task194->add_dep(task69);
  residualq->add_task(task194);

  auto I206 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task195 = make_shared<Task195>(I205, Gamma35_(), I206);
  task194->add_dep(task195);
  task195->add_dep(task69);
  residualq->add_task(task195);

  auto task196 = make_shared<Task196>(I206, t2);
  task195->add_dep(task196);
  task196->add_dep(task69);
  residualq->add_task(task196);

  auto I208 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task197 = make_shared<Task197>(I204, f1_, I208);
  task193->add_dep(task197);
  task197->add_dep(task69);
  residualq->add_task(task197);

  auto task198 = make_shared<Task198>(I208, Gamma32_(), t2);
  task197->add_dep(task198);
  task198->add_dep(task69);
  residualq->add_task(task198);

  auto task199 = make_shared<Task199>(I208, Gamma35_(), t2);
  task197->add_dep(task199);
  task199->add_dep(task69);
  residualq->add_task(task199);

  auto I217 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task200 = make_shared<Task200>(I204, f1_, I217);
  task193->add_dep(task200);
  task200->add_dep(task69);
  residualq->add_task(task200);

  auto task201 = make_shared<Task201>(I217, Gamma60_(), t2);
  task200->add_dep(task201);
  task201->add_dep(task69);
  residualq->add_task(task201);

  auto I220 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task202 = make_shared<Task202>(I204, f1_, I220);
  task193->add_dep(task202);
  task202->add_dep(task69);
  residualq->add_task(task202);

  auto task203 = make_shared<Task203>(I220, Gamma60_(), t2);
  task202->add_dep(task203);
  task203->add_dep(task69);
  residualq->add_task(task203);

  auto I223 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task204 = make_shared<Task204>(I204, t2, I223);
  task193->add_dep(task204);
  task204->add_dep(task69);
  residualq->add_task(task204);

  auto task205 = make_shared<Task205>(I223, Gamma38_(), f1_);
  task204->add_dep(task205);
  task205->add_dep(task69);
  residualq->add_task(task205);

  auto I226 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task206 = make_shared<Task206>(I204, t2, I226);
  task193->add_dep(task206);
  task206->add_dep(task69);
  residualq->add_task(task206);

  auto task207 = make_shared<Task207>(I226, Gamma38_(), f1_);
  task206->add_dep(task207);
  task207->add_dep(task69);
  residualq->add_task(task207);

  auto I229 = make_shared<TATensor<double,4>>({active_, virt_, closed_, virt_});
  auto task208 = make_shared<Task208>(I204, Gamma79_(), I229);
  task193->add_dep(task208);
  task208->add_dep(task69);
  residualq->add_task(task208);

  auto task209 = make_shared<Task209>(I229, t2);
  task208->add_dep(task209);
  task209->add_dep(task69);
  residualq->add_task(task209);

  auto I233 = make_shared<TATensor<double,4>>({closed_, active_, virt_, virt_});
  auto task210 = make_shared<Task210>(I204, Gamma38_(), I233);
  task193->add_dep(task210);
  task210->add_dep(task69);
  residualq->add_task(task210);

  auto task211 = make_shared<Task211>(I233, t2, this->e0_);
  task210->add_dep(task211);
  task211->add_dep(task69);
  residualq->add_task(task211);

  auto task212 = make_shared<Task212>(I233, t2, f1_);
  task210->add_dep(task212);
  task212->add_dep(task69);
  residualq->add_task(task212);

  auto task213 = make_shared<Task213>(I233, t2, f1_);
  task210->add_dep(task213);
  task213->add_dep(task69);
  residualq->add_task(task213);

  auto task214 = make_shared<Task214>(I233, t2, f1_);
  task210->add_dep(task214);
  task214->add_dep(task69);
  residualq->add_task(task214);

  auto task215 = make_shared<Task215>(I233, t2, f1_);
  task210->add_dep(task215);
  task215->add_dep(task69);
  residualq->add_task(task215);

  auto task216 = make_shared<Task216>(I233, t2, f1_);
  task210->add_dep(task216);
  task216->add_dep(task69);
  residualq->add_task(task216);

  auto task217 = make_shared<Task217>(I233, t2, f1_);
  task210->add_dep(task217);
  task217->add_dep(task69);
  residualq->add_task(task217);

  auto I251 = make_shared<TATensor<double,4>>({virt_, virt_, active_, active_});
  auto task218 = make_shared<Task218>(I204, f1_, I251);
  task193->add_dep(task218);
  task218->add_dep(task69);
  residualq->add_task(task218);

  auto task219 = make_shared<Task219>(I251, Gamma60_(), t2);
  task218->add_dep(task219);
  task219->add_dep(task69);
  residualq->add_task(task219);

  auto I253 = make_shared<TATensor<double,4>>({virt_, active_, active_, virt_});
  auto task220 = make_shared<Task220>(r, I253);
  task220->add_dep(task69);
  residualq->add_task(task220);

  auto I254 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task221 = make_shared<Task221>(I253, f1_, I254);
  task220->add_dep(task221);
  task221->add_dep(task69);
  residualq->add_task(task221);

  auto task222 = make_shared<Task222>(I254, Gamma59_(), t2);
  task221->add_dep(task222);
  task222->add_dep(task69);
  residualq->add_task(task222);

  auto I257 = make_shared<TATensor<double,4>>({active_, active_, virt_, virt_});
  auto task223 = make_shared<Task223>(I253, Gamma60_(), I257);
  task220->add_dep(task223);
  task223->add_dep(task69);
  residualq->add_task(task223);

  auto task224 = make_shared<Task224>(I257, t2, f1_);
  task223->add_dep(task224);
  task224->add_dep(task69);
  residualq->add_task(task224);

  auto task225 = make_shared<Task225>(I257, t2, f1_);
  task223->add_dep(task225);
  task225->add_dep(task69);
  residualq->add_task(task225);

  auto I259 = make_shared<TATensor<double,4>>({virt_, virt_, active_, active_});
  auto task226 = make_shared<Task226>(r, I259);
  task226->add_dep(task69);
  residualq->add_task(task226);

  auto task227 = make_shared<Task227>(I259, Gamma90_(), t2);
  task226->add_dep(task227);
  task227->add_dep(task69);
  residualq->add_task(task227);

  auto task228 = make_shared<Task228>(I259, Gamma60_(), t2, this->e0_);
  task226->add_dep(task228);
  task228->add_dep(task69);
  residualq->add_task(task228);

  return residualq;
}


#endif
