//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_deciqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_deciq() {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deciq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor119 = {deci};
  auto task119 = make_shared<Task119>(tensor119);
  deciq->add_task(task119);

  vector<IndexRange> I120_index = {ci_};
  auto I120 = make_shared<Tensor>(I120_index);
  vector<shared_ptr<Tensor>> tensor120 = {deci, I120};
  auto task120 = make_shared<Task120>(tensor120, cindex);
  task120->add_dep(task119);
  deciq->add_task(task120);

  vector<IndexRange> I121_index;
  auto I121 = make_shared<Tensor>(I121_index);
  vector<shared_ptr<Tensor>> tensor121 = {I120, Gamma36_(), I121};
  auto task121 = make_shared<Task121>(tensor121, cindex);
  task120->add_dep(task121);
  task121->add_dep(task119);
  deciq->add_task(task121);

  vector<IndexRange> I122_index = {closed_, virt_, closed_, virt_};
  auto I122 = make_shared<Tensor>(I122_index);
  vector<shared_ptr<Tensor>> tensor122 = {I121, t2, I122};
  auto task122 = make_shared<Task122>(tensor122, cindex);
  task121->add_dep(task122);
  task122->add_dep(task119);
  deciq->add_task(task122);

  vector<shared_ptr<Tensor>> tensor123 = {I122, t2};
  auto task123 = make_shared<Task123>(tensor123, cindex);
  task122->add_dep(task123);
  task123->add_dep(task119);
  deciq->add_task(task123);

  vector<IndexRange> I125_index = {virt_, closed_, virt_, closed_};
  auto I125 = make_shared<Tensor>(I125_index);
  vector<shared_ptr<Tensor>> tensor124 = {I121, t2, I125};
  auto task124 = make_shared<Task124>(tensor124, cindex);
  task121->add_dep(task124);
  task124->add_dep(task119);
  deciq->add_task(task124);

  vector<shared_ptr<Tensor>> tensor125 = {I125, t2};
  auto task125 = make_shared<Task125>(tensor125, cindex);
  task124->add_dep(task125);
  task125->add_dep(task119);
  deciq->add_task(task125);

  vector<IndexRange> I127_index = {active_, active_};
  auto I127 = make_shared<Tensor>(I127_index);
  vector<shared_ptr<Tensor>> tensor126 = {I120, Gamma38_(), I127};
  auto task126 = make_shared<Task126>(tensor126, cindex);
  task120->add_dep(task126);
  task126->add_dep(task119);
  deciq->add_task(task126);

  vector<IndexRange> I128_index = {closed_, active_};
  auto I128 = make_shared<Tensor>(I128_index);
  vector<shared_ptr<Tensor>> tensor127 = {I127, f1_, I128};
  auto task127 = make_shared<Task127>(tensor127, cindex);
  task126->add_dep(task127);
  task127->add_dep(task119);
  deciq->add_task(task127);

  vector<IndexRange> I129_index = {virt_, closed_, virt_, closed_};
  auto I129 = make_shared<Tensor>(I129_index);
  vector<shared_ptr<Tensor>> tensor128 = {I128, t2, I129};
  auto task128 = make_shared<Task128>(tensor128, cindex);
  task127->add_dep(task128);
  task128->add_dep(task119);
  deciq->add_task(task128);

  vector<shared_ptr<Tensor>> tensor129 = {I129, t2};
  auto task129 = make_shared<Task129>(tensor129, cindex);
  task128->add_dep(task129);
  task129->add_dep(task119);
  deciq->add_task(task129);

  vector<IndexRange> I133_index = {virt_, closed_, virt_, closed_};
  auto I133 = make_shared<Tensor>(I133_index);
  vector<shared_ptr<Tensor>> tensor130 = {I128, t2, I133};
  auto task130 = make_shared<Task130>(tensor130, cindex);
  task127->add_dep(task130);
  task130->add_dep(task119);
  deciq->add_task(task130);

  vector<shared_ptr<Tensor>> tensor131 = {I133, t2};
  auto task131 = make_shared<Task131>(tensor131, cindex);
  task130->add_dep(task131);
  task131->add_dep(task119);
  deciq->add_task(task131);

  vector<IndexRange> I136_index = {active_, closed_};
  auto I136 = make_shared<Tensor>(I136_index);
  vector<shared_ptr<Tensor>> tensor132 = {I127, f1_, I136};
  auto task132 = make_shared<Task132>(tensor132, cindex);
  task126->add_dep(task132);
  task132->add_dep(task119);
  deciq->add_task(task132);

  vector<IndexRange> I137_index = {virt_, closed_, virt_, active_};
  auto I137 = make_shared<Tensor>(I137_index);
  vector<shared_ptr<Tensor>> tensor133 = {I136, t2, I137};
  auto task133 = make_shared<Task133>(tensor133, cindex);
  task132->add_dep(task133);
  task133->add_dep(task119);
  deciq->add_task(task133);

  vector<shared_ptr<Tensor>> tensor134 = {I137, t2};
  auto task134 = make_shared<Task134>(tensor134, cindex);
  task133->add_dep(task134);
  task134->add_dep(task119);
  deciq->add_task(task134);

  vector<IndexRange> I141_index = {virt_, closed_, virt_, active_};
  auto I141 = make_shared<Tensor>(I141_index);
  vector<shared_ptr<Tensor>> tensor135 = {I136, t2, I141};
  auto task135 = make_shared<Task135>(tensor135, cindex);
  task132->add_dep(task135);
  task135->add_dep(task119);
  deciq->add_task(task135);

  vector<shared_ptr<Tensor>> tensor136 = {I141, t2};
  auto task136 = make_shared<Task136>(tensor136, cindex);
  task135->add_dep(task136);
  task136->add_dep(task119);
  deciq->add_task(task136);

  vector<IndexRange> I150_index = {virt_, virt_, active_, closed_};
  auto I150 = make_shared<Tensor>(I150_index);
  vector<shared_ptr<Tensor>> tensor137 = {I127, t2, I150};
  auto task137 = make_shared<Task137>(tensor137, cindex);
  task126->add_dep(task137);
  task137->add_dep(task119);
  deciq->add_task(task137);

  vector<IndexRange> I151_index = {virt_, closed_, virt_, active_};
  auto I151 = make_shared<Tensor>(I151_index);
  vector<shared_ptr<Tensor>> tensor138 = {I150, f1_, I151};
  auto task138 = make_shared<Task138>(tensor138, cindex);
  task137->add_dep(task138);
  task138->add_dep(task119);
  deciq->add_task(task138);

  vector<shared_ptr<Tensor>> tensor139 = {I151, t2};
  auto task139 = make_shared<Task139>(tensor139, cindex);
  task138->add_dep(task139);
  task139->add_dep(task119);
  deciq->add_task(task139);

  vector<IndexRange> I154_index = {virt_, virt_, active_, closed_};
  auto I154 = make_shared<Tensor>(I154_index);
  vector<shared_ptr<Tensor>> tensor140 = {I127, t2, I154};
  auto task140 = make_shared<Task140>(tensor140, cindex);
  task126->add_dep(task140);
  task140->add_dep(task119);
  deciq->add_task(task140);

  vector<IndexRange> I155_index = {virt_, closed_, virt_, active_};
  auto I155 = make_shared<Tensor>(I155_index);
  vector<shared_ptr<Tensor>> tensor141 = {I154, f1_, I155};
  auto task141 = make_shared<Task141>(tensor141, cindex);
  task140->add_dep(task141);
  task141->add_dep(task119);
  deciq->add_task(task141);

  vector<shared_ptr<Tensor>> tensor142 = {I155, t2};
  auto task142 = make_shared<Task142>(tensor142, cindex);
  task141->add_dep(task142);
  task142->add_dep(task119);
  deciq->add_task(task142);

  vector<IndexRange> I158_index = {virt_, closed_, active_, virt_};
  auto I158 = make_shared<Tensor>(I158_index);
  vector<shared_ptr<Tensor>> tensor143 = {I127, t2, I158};
  auto task143 = make_shared<Task143>(tensor143, cindex);
  task126->add_dep(task143);
  task143->add_dep(task119);
  deciq->add_task(task143);

  vector<IndexRange> I159_index = {virt_, closed_, virt_, active_};
  auto I159 = make_shared<Tensor>(I159_index);
  vector<shared_ptr<Tensor>> tensor144 = {I158, f1_, I159};
  auto task144 = make_shared<Task144>(tensor144, cindex);
  task143->add_dep(task144);
  task144->add_dep(task119);
  deciq->add_task(task144);

  vector<shared_ptr<Tensor>> tensor145 = {I159, t2};
  auto task145 = make_shared<Task145>(tensor145, cindex);
  task144->add_dep(task145);
  task145->add_dep(task119);
  deciq->add_task(task145);

  vector<IndexRange> I162_index = {virt_, closed_, active_, virt_};
  auto I162 = make_shared<Tensor>(I162_index);
  vector<shared_ptr<Tensor>> tensor146 = {I127, t2, I162};
  auto task146 = make_shared<Task146>(tensor146, cindex);
  task126->add_dep(task146);
  task146->add_dep(task119);
  deciq->add_task(task146);

  vector<IndexRange> I163_index = {virt_, closed_, virt_, active_};
  auto I163 = make_shared<Tensor>(I163_index);
  vector<shared_ptr<Tensor>> tensor147 = {I162, f1_, I163};
  auto task147 = make_shared<Task147>(tensor147, cindex);
  task146->add_dep(task147);
  task147->add_dep(task119);
  deciq->add_task(task147);

  vector<shared_ptr<Tensor>> tensor148 = {I163, t2};
  auto task148 = make_shared<Task148>(tensor148, cindex);
  task147->add_dep(task148);
  task148->add_dep(task119);
  deciq->add_task(task148);

  vector<IndexRange> I166_index = {closed_, virt_, active_, virt_};
  auto I166 = make_shared<Tensor>(I166_index);
  vector<shared_ptr<Tensor>> tensor149 = {I127, t2, I166};
  auto task149 = make_shared<Task149>(tensor149, cindex);
  task126->add_dep(task149);
  task149->add_dep(task119);
  deciq->add_task(task149);

  vector<IndexRange> I167_index = {virt_, closed_, virt_, active_};
  auto I167 = make_shared<Tensor>(I167_index);
  vector<shared_ptr<Tensor>> tensor150 = {I166, f1_, I167};
  auto task150 = make_shared<Task150>(tensor150, cindex);
  task149->add_dep(task150);
  task150->add_dep(task119);
  deciq->add_task(task150);

  vector<shared_ptr<Tensor>> tensor151 = {I167, t2};
  auto task151 = make_shared<Task151>(tensor151, cindex);
  task150->add_dep(task151);
  task151->add_dep(task119);
  deciq->add_task(task151);

  vector<IndexRange> I170_index = {closed_, virt_, active_, virt_};
  auto I170 = make_shared<Tensor>(I170_index);
  vector<shared_ptr<Tensor>> tensor152 = {I127, t2, I170};
  auto task152 = make_shared<Task152>(tensor152, cindex);
  task126->add_dep(task152);
  task152->add_dep(task119);
  deciq->add_task(task152);

  vector<IndexRange> I171_index = {virt_, closed_, virt_, active_};
  auto I171 = make_shared<Tensor>(I171_index);
  vector<shared_ptr<Tensor>> tensor153 = {I170, f1_, I171};
  auto task153 = make_shared<Task153>(tensor153, cindex);
  task152->add_dep(task153);
  task153->add_dep(task119);
  deciq->add_task(task153);

  vector<shared_ptr<Tensor>> tensor154 = {I171, t2};
  auto task154 = make_shared<Task154>(tensor154, cindex);
  task153->add_dep(task154);
  task154->add_dep(task119);
  deciq->add_task(task154);

  vector<IndexRange> I174_index = {virt_, closed_, virt_, active_};
  auto I174 = make_shared<Tensor>(I174_index);
  vector<shared_ptr<Tensor>> tensor155 = {I127, t2, I174};
  auto task155 = make_shared<Task155>(tensor155, cindex);
  task126->add_dep(task155);
  task155->add_dep(task119);
  deciq->add_task(task155);

  vector<shared_ptr<Tensor>> tensor156 = {I174, t2};
  auto task156 = make_shared<Task156>(tensor156, cindex, this->e0_);
  task155->add_dep(task156);
  task156->add_dep(task119);
  deciq->add_task(task156);

  vector<IndexRange> I177_index = {virt_, closed_, virt_, active_};
  auto I177 = make_shared<Tensor>(I177_index);
  vector<shared_ptr<Tensor>> tensor157 = {I127, t2, I177};
  auto task157 = make_shared<Task157>(tensor157, cindex);
  task126->add_dep(task157);
  task157->add_dep(task119);
  deciq->add_task(task157);

  vector<shared_ptr<Tensor>> tensor158 = {I177, t2};
  auto task158 = make_shared<Task158>(tensor158, cindex, this->e0_);
  task157->add_dep(task158);
  task158->add_dep(task119);
  deciq->add_task(task158);

  vector<IndexRange> I180_index = {virt_, closed_, virt_, active_};
  auto I180 = make_shared<Tensor>(I180_index);
  vector<shared_ptr<Tensor>> tensor159 = {I127, v2_, I180};
  auto task159 = make_shared<Task159>(tensor159, cindex);
  task126->add_dep(task159);
  task159->add_dep(task119);
  deciq->add_task(task159);

  vector<shared_ptr<Tensor>> tensor160 = {I180, t2};
  auto task160 = make_shared<Task160>(tensor160, cindex);
  task159->add_dep(task160);
  task160->add_dep(task119);
  deciq->add_task(task160);

  vector<IndexRange> I183_index = {virt_, closed_, virt_, active_};
  auto I183 = make_shared<Tensor>(I183_index);
  vector<shared_ptr<Tensor>> tensor161 = {I127, v2_, I183};
  auto task161 = make_shared<Task161>(tensor161, cindex);
  task126->add_dep(task161);
  task161->add_dep(task119);
  deciq->add_task(task161);

  vector<shared_ptr<Tensor>> tensor162 = {I183, t2};
  auto task162 = make_shared<Task162>(tensor162, cindex);
  task161->add_dep(task162);
  task162->add_dep(task119);
  deciq->add_task(task162);

  vector<IndexRange> I186_index = {virt_, closed_, virt_, active_};
  auto I186 = make_shared<Tensor>(I186_index);
  vector<shared_ptr<Tensor>> tensor163 = {I127, v2_, I186};
  auto task163 = make_shared<Task163>(tensor163, cindex);
  task126->add_dep(task163);
  task163->add_dep(task119);
  deciq->add_task(task163);

  vector<shared_ptr<Tensor>> tensor164 = {I186, t2};
  auto task164 = make_shared<Task164>(tensor164, cindex);
  task163->add_dep(task164);
  task164->add_dep(task119);
  deciq->add_task(task164);

  vector<IndexRange> I189_index = {virt_, closed_, virt_, active_};
  auto I189 = make_shared<Tensor>(I189_index);
  vector<shared_ptr<Tensor>> tensor165 = {I127, v2_, I189};
  auto task165 = make_shared<Task165>(tensor165, cindex);
  task126->add_dep(task165);
  task165->add_dep(task119);
  deciq->add_task(task165);

  vector<shared_ptr<Tensor>> tensor166 = {I189, t2};
  auto task166 = make_shared<Task166>(tensor166, cindex);
  task165->add_dep(task166);
  task166->add_dep(task119);
  deciq->add_task(task166);

  vector<IndexRange> I143_index = {active_, active_};
  auto I143 = make_shared<Tensor>(I143_index);
  vector<shared_ptr<Tensor>> tensor167 = {I120, Gamma42_(), I143};
  auto task167 = make_shared<Task167>(tensor167, cindex);
  task120->add_dep(task167);
  task167->add_dep(task119);
  deciq->add_task(task167);

  vector<IndexRange> I144_index = {virt_, closed_, virt_, active_};
  auto I144 = make_shared<Tensor>(I144_index);
  vector<shared_ptr<Tensor>> tensor168 = {I143, t2, I144};
  auto task168 = make_shared<Task168>(tensor168, cindex);
  task167->add_dep(task168);
  task168->add_dep(task119);
  deciq->add_task(task168);

  vector<shared_ptr<Tensor>> tensor169 = {I144, t2};
  auto task169 = make_shared<Task169>(tensor169, cindex);
  task168->add_dep(task169);
  task169->add_dep(task119);
  deciq->add_task(task169);

  vector<IndexRange> I147_index = {virt_, closed_, virt_, active_};
  auto I147 = make_shared<Tensor>(I147_index);
  vector<shared_ptr<Tensor>> tensor170 = {I143, t2, I147};
  auto task170 = make_shared<Task170>(tensor170, cindex);
  task167->add_dep(task170);
  task170->add_dep(task119);
  deciq->add_task(task170);

  vector<shared_ptr<Tensor>> tensor171 = {I147, t2};
  auto task171 = make_shared<Task171>(tensor171, cindex);
  task170->add_dep(task171);
  task171->add_dep(task119);
  deciq->add_task(task171);

  return deciq;
}


