//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_residualqq.cc
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


#include <src/smith/mrci/MRCI.h>
#include <src/smith/mrci/MRCI_tasks4.h>
#include <src/smith/mrci/MRCI_tasks5.h>
#include <src/smith/mrci/MRCI_tasks6.h>
#include <src/smith/mrci/MRCI_tasks7.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void MRCI::MRCI::make_residualq3(shared_ptr<Queue> residualq, shared_ptr<Task> task108, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I27_index = {closed_, closed_, active_, virt_};
  auto I27 = make_shared<Tensor>(I27_index);
  auto tensor193 = vector<shared_ptr<Tensor>>{r, I27};
  auto task193 = make_shared<Task193>(tensor193, pindex);
  task193->add_dep(task108);
  residualq->add_task(task193);

  vector<IndexRange> I28_index = {closed_, closed_, active_, active_};
  auto I28 = make_shared<Tensor>(I28_index);
  auto tensor194 = vector<shared_ptr<Tensor>>{I27, h1_, I28};
  auto task194 = make_shared<Task194>(tensor194, pindex);
  task193->add_dep(task194);
  task194->add_dep(task108);
  residualq->add_task(task194);

  auto tensor195 = vector<shared_ptr<Tensor>>{I28, Gamma2_(), t2};
  auto task195 = make_shared<Task195>(tensor195, pindex);
  task194->add_dep(task195);
  task195->add_dep(task108);
  residualq->add_task(task195);

  vector<IndexRange> I31_index = {closed_, active_};
  auto I31 = make_shared<Tensor>(I31_index);
  auto tensor196 = vector<shared_ptr<Tensor>>{I27, h1_, I31};
  auto task196 = make_shared<Task196>(tensor196, pindex);
  task193->add_dep(task196);
  task196->add_dep(task108);
  residualq->add_task(task196);

  auto tensor197 = vector<shared_ptr<Tensor>>{I31, Gamma10_(), t2};
  auto task197 = make_shared<Task197>(tensor197, pindex);
  task196->add_dep(task197);
  task197->add_dep(task108);
  residualq->add_task(task197);

  vector<IndexRange> I34_index = {closed_, active_};
  auto I34 = make_shared<Tensor>(I34_index);
  auto tensor198 = vector<shared_ptr<Tensor>>{I27, h1_, I34};
  auto task198 = make_shared<Task198>(tensor198, pindex);
  task193->add_dep(task198);
  task198->add_dep(task108);
  residualq->add_task(task198);

  auto tensor199 = vector<shared_ptr<Tensor>>{I34, Gamma10_(), t2};
  auto task199 = make_shared<Task199>(tensor199, pindex);
  task198->add_dep(task199);
  task199->add_dep(task108);
  residualq->add_task(task199);

  vector<IndexRange> I37_index = {closed_, virt_, closed_, active_};
  auto I37 = make_shared<Tensor>(I37_index);
  auto tensor200 = vector<shared_ptr<Tensor>>{I27, Gamma12_(), I37};
  auto task200 = make_shared<Task200>(tensor200, pindex);
  task193->add_dep(task200);
  task200->add_dep(task108);
  residualq->add_task(task200);

  auto tensor201 = vector<shared_ptr<Tensor>>{I37, t2, h1_};
  auto task201 = make_shared<Task201>(tensor201, pindex);
  task200->add_dep(task201);
  task201->add_dep(task108);
  residualq->add_task(task201);

  auto tensor202 = vector<shared_ptr<Tensor>>{I37, t2, h1_};
  auto task202 = make_shared<Task202>(tensor202, pindex);
  task200->add_dep(task202);
  task202->add_dep(task108);
  residualq->add_task(task202);

  auto tensor203 = vector<shared_ptr<Tensor>>{I37, t2, h1_};
  auto task203 = make_shared<Task203>(tensor203, pindex);
  task200->add_dep(task203);
  task203->add_dep(task108);
  residualq->add_task(task203);

  auto tensor204 = vector<shared_ptr<Tensor>>{I37, t2, h1_};
  auto task204 = make_shared<Task204>(tensor204, pindex);
  task200->add_dep(task204);
  task204->add_dep(task108);
  residualq->add_task(task204);

  auto tensor205 = vector<shared_ptr<Tensor>>{I37, t2, h1_};
  auto task205 = make_shared<Task205>(tensor205, pindex);
  task200->add_dep(task205);
  task205->add_dep(task108);
  residualq->add_task(task205);

  auto tensor206 = vector<shared_ptr<Tensor>>{I37, t2, h1_};
  auto task206 = make_shared<Task206>(tensor206, pindex);
  task200->add_dep(task206);
  task206->add_dep(task108);
  residualq->add_task(task206);

  auto tensor207 = vector<shared_ptr<Tensor>>{I37, t2, v2_};
  auto task207 = make_shared<Task207>(tensor207, pindex);
  task200->add_dep(task207);
  task207->add_dep(task108);
  residualq->add_task(task207);

  auto tensor208 = vector<shared_ptr<Tensor>>{I37, t2, v2_};
  auto task208 = make_shared<Task208>(tensor208, pindex);
  task200->add_dep(task208);
  task208->add_dep(task108);
  residualq->add_task(task208);

  auto tensor209 = vector<shared_ptr<Tensor>>{I37, t2, v2_};
  auto task209 = make_shared<Task209>(tensor209, pindex);
  task200->add_dep(task209);
  task209->add_dep(task108);
  residualq->add_task(task209);

  auto tensor210 = vector<shared_ptr<Tensor>>{I37, t2, v2_};
  auto task210 = make_shared<Task210>(tensor210, pindex);
  task200->add_dep(task210);
  task210->add_dep(task108);
  residualq->add_task(task210);

  auto tensor211 = vector<shared_ptr<Tensor>>{I37, t2, v2_};
  auto task211 = make_shared<Task211>(tensor211, pindex);
  task200->add_dep(task211);
  task211->add_dep(task108);
  residualq->add_task(task211);

  auto tensor212 = vector<shared_ptr<Tensor>>{I37, t2, v2_};
  auto task212 = make_shared<Task212>(tensor212, pindex);
  task200->add_dep(task212);
  task212->add_dep(task108);
  residualq->add_task(task212);

  auto tensor213 = vector<shared_ptr<Tensor>>{I37, t2, v2_};
  auto task213 = make_shared<Task213>(tensor213, pindex);
  task200->add_dep(task213);
  task213->add_dep(task108);
  residualq->add_task(task213);

  auto tensor214 = vector<shared_ptr<Tensor>>{I37, t2, v2_};
  auto task214 = make_shared<Task214>(tensor214, pindex);
  task200->add_dep(task214);
  task214->add_dep(task108);
  residualq->add_task(task214);

  auto tensor215 = vector<shared_ptr<Tensor>>{I37, t2, v2_};
  auto task215 = make_shared<Task215>(tensor215, pindex);
  task200->add_dep(task215);
  task215->add_dep(task108);
  residualq->add_task(task215);

  auto tensor216 = vector<shared_ptr<Tensor>>{I37, t2, v2_};
  auto task216 = make_shared<Task216>(tensor216, pindex);
  task200->add_dep(task216);
  task216->add_dep(task108);
  residualq->add_task(task216);

  vector<IndexRange> I613_index = {virt_, active_, closed_, closed_};
  auto I613 = make_shared<Tensor>(I613_index);
  auto tensor217 = vector<shared_ptr<Tensor>>{I37, t2, I613};
  auto task217 = make_shared<Task217>(tensor217, pindex);
  task200->add_dep(task217);
  task217->add_dep(task108);
  residualq->add_task(task217);

  auto tensor218 = vector<shared_ptr<Tensor>>{I613, v2_};
  auto task218 = make_shared<Task218>(tensor218, pindex);
  task217->add_dep(task218);
  task218->add_dep(task108);
  residualq->add_task(task218);

  vector<IndexRange> I616_index = {virt_, active_, closed_, closed_};
  auto I616 = make_shared<Tensor>(I616_index);
  auto tensor219 = vector<shared_ptr<Tensor>>{I37, t2, I616};
  auto task219 = make_shared<Task219>(tensor219, pindex);
  task200->add_dep(task219);
  task219->add_dep(task108);
  residualq->add_task(task219);

  auto tensor220 = vector<shared_ptr<Tensor>>{I616, v2_};
  auto task220 = make_shared<Task220>(tensor220, pindex);
  task219->add_dep(task220);
  task220->add_dep(task108);
  residualq->add_task(task220);

  vector<IndexRange> I625_index = {virt_, active_, closed_, closed_};
  auto I625 = make_shared<Tensor>(I625_index);
  auto tensor221 = vector<shared_ptr<Tensor>>{I37, t2, I625};
  auto task221 = make_shared<Task221>(tensor221, pindex);
  task200->add_dep(task221);
  task221->add_dep(task108);
  residualq->add_task(task221);

  auto tensor222 = vector<shared_ptr<Tensor>>{I625, v2_};
  auto task222 = make_shared<Task222>(tensor222, pindex);
  task221->add_dep(task222);
  task222->add_dep(task108);
  residualq->add_task(task222);

  vector<IndexRange> I628_index = {virt_, active_, closed_, closed_};
  auto I628 = make_shared<Tensor>(I628_index);
  auto tensor223 = vector<shared_ptr<Tensor>>{I37, t2, I628};
  auto task223 = make_shared<Task223>(tensor223, pindex);
  task200->add_dep(task223);
  task223->add_dep(task108);
  residualq->add_task(task223);

  auto tensor224 = vector<shared_ptr<Tensor>>{I628, v2_};
  auto task224 = make_shared<Task224>(tensor224, pindex);
  task223->add_dep(task224);
  task224->add_dep(task108);
  residualq->add_task(task224);

  auto tensor225 = vector<shared_ptr<Tensor>>{I37, t2, v2_};
  auto task225 = make_shared<Task225>(tensor225, pindex);
  task200->add_dep(task225);
  task225->add_dep(task108);
  residualq->add_task(task225);

  auto tensor226 = vector<shared_ptr<Tensor>>{I37, t2, v2_};
  auto task226 = make_shared<Task226>(tensor226, pindex);
  task200->add_dep(task226);
  task226->add_dep(task108);
  residualq->add_task(task226);

  vector<IndexRange> I55_index = {virt_, closed_, active_, active_};
  auto I55 = make_shared<Tensor>(I55_index);
  auto tensor227 = vector<shared_ptr<Tensor>>{I27, h1_, I55};
  auto task227 = make_shared<Task227>(tensor227, pindex);
  task193->add_dep(task227);
  task227->add_dep(task108);
  residualq->add_task(task227);

  auto tensor228 = vector<shared_ptr<Tensor>>{I55, Gamma18_(), t2};
  auto task228 = make_shared<Task228>(tensor228, pindex);
  task227->add_dep(task228);
  task228->add_dep(task108);
  residualq->add_task(task228);

  auto tensor229 = vector<shared_ptr<Tensor>>{I55, Gamma10_(), t2};
  auto task229 = make_shared<Task229>(tensor229, pindex);
  task227->add_dep(task229);
  task229->add_dep(task108);
  residualq->add_task(task229);

  vector<IndexRange> I58_index = {virt_, closed_, active_, active_};
  auto I58 = make_shared<Tensor>(I58_index);
  auto tensor230 = vector<shared_ptr<Tensor>>{I27, h1_, I58};
  auto task230 = make_shared<Task230>(tensor230, pindex);
  task193->add_dep(task230);
  task230->add_dep(task108);
  residualq->add_task(task230);

  vector<IndexRange> I59_index = {active_, virt_, closed_, active_};
  auto I59 = make_shared<Tensor>(I59_index);
  auto tensor231 = vector<shared_ptr<Tensor>>{I58, Gamma10_(), I59};
  auto task231 = make_shared<Task231>(tensor231, pindex);
  task230->add_dep(task231);
  task231->add_dep(task108);
  residualq->add_task(task231);

  auto tensor232 = vector<shared_ptr<Tensor>>{I59, t2};
  auto task232 = make_shared<Task232>(tensor232, pindex);
  task231->add_dep(task232);
  task232->add_dep(task108);
  residualq->add_task(task232);

  vector<IndexRange> I67_index = {virt_, active_};
  auto I67 = make_shared<Tensor>(I67_index);
  auto tensor233 = vector<shared_ptr<Tensor>>{I27, t2, I67};
  auto task233 = make_shared<Task233>(tensor233, pindex);
  task193->add_dep(task233);
  task233->add_dep(task108);
  residualq->add_task(task233);

  auto tensor234 = vector<shared_ptr<Tensor>>{I67, Gamma12_(), h1_};
  auto task234 = make_shared<Task234>(tensor234, pindex);
  task233->add_dep(task234);
  task234->add_dep(task108);
  residualq->add_task(task234);

  auto tensor235 = vector<shared_ptr<Tensor>>{I67, Gamma197_(), v2_};
  auto task235 = make_shared<Task235>(tensor235, pindex);
  task233->add_dep(task235);
  task235->add_dep(task108);
  residualq->add_task(task235);

  auto tensor236 = vector<shared_ptr<Tensor>>{I67, Gamma10_(), v2_};
  auto task236 = make_shared<Task236>(tensor236, pindex);
  task233->add_dep(task236);
  task236->add_dep(task108);
  residualq->add_task(task236);

  vector<IndexRange> I70_index = {virt_, active_};
  auto I70 = make_shared<Tensor>(I70_index);
  auto tensor237 = vector<shared_ptr<Tensor>>{I27, t2, I70};
  auto task237 = make_shared<Task237>(tensor237, pindex);
  task193->add_dep(task237);
  task237->add_dep(task108);
  residualq->add_task(task237);

  auto tensor238 = vector<shared_ptr<Tensor>>{I70, Gamma12_(), h1_};
  auto task238 = make_shared<Task238>(tensor238, pindex);
  task237->add_dep(task238);
  task238->add_dep(task108);
  residualq->add_task(task238);

  auto tensor239 = vector<shared_ptr<Tensor>>{I70, Gamma197_(), v2_};
  auto task239 = make_shared<Task239>(tensor239, pindex);
  task237->add_dep(task239);
  task239->add_dep(task108);
  residualq->add_task(task239);

  auto tensor240 = vector<shared_ptr<Tensor>>{I70, Gamma10_(), v2_};
  auto task240 = make_shared<Task240>(tensor240, pindex);
  task237->add_dep(task240);
  task240->add_dep(task108);
  residualq->add_task(task240);

  vector<IndexRange> I387_index = {virt_, active_, active_, active_};
  auto I387 = make_shared<Tensor>(I387_index);
  auto tensor241 = vector<shared_ptr<Tensor>>{I27, t2, I387};
  auto task241 = make_shared<Task241>(tensor241, pindex);
  task193->add_dep(task241);
  task241->add_dep(task108);
  residualq->add_task(task241);

  auto tensor242 = vector<shared_ptr<Tensor>>{I387, Gamma126_(), v2_};
  auto task242 = make_shared<Task242>(tensor242, pindex);
  task241->add_dep(task242);
  task242->add_dep(task108);
  residualq->add_task(task242);

  auto tensor243 = vector<shared_ptr<Tensor>>{I387, Gamma88_(), v2_};
  auto task243 = make_shared<Task243>(tensor243, pindex);
  task241->add_dep(task243);
  task243->add_dep(task108);
  residualq->add_task(task243);

  vector<IndexRange> I393_index = {closed_, closed_, active_, active_};
  auto I393 = make_shared<Tensor>(I393_index);
  auto tensor244 = vector<shared_ptr<Tensor>>{I27, v2_, I393};
  auto task244 = make_shared<Task244>(tensor244, pindex);
  task193->add_dep(task244);
  task244->add_dep(task108);
  residualq->add_task(task244);

  auto tensor245 = vector<shared_ptr<Tensor>>{I393, Gamma0_(), t2};
  auto task245 = make_shared<Task245>(tensor245, pindex);
  task244->add_dep(task245);
  task245->add_dep(task108);
  residualq->add_task(task245);

  vector<IndexRange> I396_index = {closed_, closed_, active_, active_};
  auto I396 = make_shared<Tensor>(I396_index);
  auto tensor246 = vector<shared_ptr<Tensor>>{I27, v2_, I396};
  auto task246 = make_shared<Task246>(tensor246, pindex);
  task193->add_dep(task246);
  task246->add_dep(task108);
  residualq->add_task(task246);

  auto tensor247 = vector<shared_ptr<Tensor>>{I396, Gamma0_(), t2};
  auto task247 = make_shared<Task247>(tensor247, pindex);
  task246->add_dep(task247);
  task247->add_dep(task108);
  residualq->add_task(task247);

  vector<IndexRange> I399_index = {closed_, closed_, active_, active_};
  auto I399 = make_shared<Tensor>(I399_index);
  auto tensor248 = vector<shared_ptr<Tensor>>{I27, v2_, I399};
  auto task248 = make_shared<Task248>(tensor248, pindex);
  task193->add_dep(task248);
  task248->add_dep(task108);
  residualq->add_task(task248);

  auto tensor249 = vector<shared_ptr<Tensor>>{I399, Gamma0_(), t2};
  auto task249 = make_shared<Task249>(tensor249, pindex);
  task248->add_dep(task249);
  task249->add_dep(task108);
  residualq->add_task(task249);

  vector<IndexRange> I402_index = {closed_, closed_, active_, active_};
  auto I402 = make_shared<Tensor>(I402_index);
  auto tensor250 = vector<shared_ptr<Tensor>>{I27, v2_, I402};
  auto task250 = make_shared<Task250>(tensor250, pindex);
  task193->add_dep(task250);
  task250->add_dep(task108);
  residualq->add_task(task250);

  auto tensor251 = vector<shared_ptr<Tensor>>{I402, Gamma2_(), t2};
  auto task251 = make_shared<Task251>(tensor251, pindex);
  task250->add_dep(task251);
  task251->add_dep(task108);
  residualq->add_task(task251);

  vector<IndexRange> I405_index = {closed_, active_, active_, active_};
  auto I405 = make_shared<Tensor>(I405_index);
  auto tensor252 = vector<shared_ptr<Tensor>>{I27, v2_, I405};
  auto task252 = make_shared<Task252>(tensor252, pindex);
  task193->add_dep(task252);
  task252->add_dep(task108);
  residualq->add_task(task252);

  auto tensor253 = vector<shared_ptr<Tensor>>{I405, Gamma132_(), t2};
  auto task253 = make_shared<Task253>(tensor253, pindex);
  task252->add_dep(task253);
  task253->add_dep(task108);
  residualq->add_task(task253);

  vector<IndexRange> I408_index = {closed_, active_, active_, active_};
  auto I408 = make_shared<Tensor>(I408_index);
  auto tensor254 = vector<shared_ptr<Tensor>>{I27, v2_, I408};
  auto task254 = make_shared<Task254>(tensor254, pindex);
  task193->add_dep(task254);
  task254->add_dep(task108);
  residualq->add_task(task254);

  auto tensor255 = vector<shared_ptr<Tensor>>{I408, Gamma132_(), t2};
  auto task255 = make_shared<Task255>(tensor255, pindex);
  task254->add_dep(task255);
  task255->add_dep(task108);
  residualq->add_task(task255);

  vector<IndexRange> I411_index = {closed_, active_, active_, active_};
  auto I411 = make_shared<Tensor>(I411_index);
  auto tensor256 = vector<shared_ptr<Tensor>>{I27, v2_, I411};
  auto task256 = make_shared<Task256>(tensor256, pindex);
  task193->add_dep(task256);
  task256->add_dep(task108);
  residualq->add_task(task256);

  auto tensor257 = vector<shared_ptr<Tensor>>{I411, Gamma1_(), t2};
  auto task257 = make_shared<Task257>(tensor257, pindex);
  task256->add_dep(task257);
  task257->add_dep(task108);
  residualq->add_task(task257);

  vector<IndexRange> I414_index = {closed_, active_, active_, active_};
  auto I414 = make_shared<Tensor>(I414_index);
  auto tensor258 = vector<shared_ptr<Tensor>>{I27, v2_, I414};
  auto task258 = make_shared<Task258>(tensor258, pindex);
  task193->add_dep(task258);
  task258->add_dep(task108);
  residualq->add_task(task258);

  auto tensor259 = vector<shared_ptr<Tensor>>{I414, Gamma87_(), t2};
  auto task259 = make_shared<Task259>(tensor259, pindex);
  task258->add_dep(task259);
  task259->add_dep(task108);
  residualq->add_task(task259);

  vector<IndexRange> I417_index = {closed_, active_, active_, active_};
  auto I417 = make_shared<Tensor>(I417_index);
  auto tensor260 = vector<shared_ptr<Tensor>>{I27, v2_, I417};
  auto task260 = make_shared<Task260>(tensor260, pindex);
  task193->add_dep(task260);
  task260->add_dep(task108);
  residualq->add_task(task260);

  auto tensor261 = vector<shared_ptr<Tensor>>{I417, Gamma132_(), t2};
  auto task261 = make_shared<Task261>(tensor261, pindex);
  task260->add_dep(task261);
  task261->add_dep(task108);
  residualq->add_task(task261);

  vector<IndexRange> I420_index = {closed_, active_, active_, active_};
  auto I420 = make_shared<Tensor>(I420_index);
  auto tensor262 = vector<shared_ptr<Tensor>>{I27, v2_, I420};
  auto task262 = make_shared<Task262>(tensor262, pindex);
  task193->add_dep(task262);
  task262->add_dep(task108);
  residualq->add_task(task262);

  auto tensor263 = vector<shared_ptr<Tensor>>{I420, Gamma137_(), t2};
  auto task263 = make_shared<Task263>(tensor263, pindex);
  task262->add_dep(task263);
  task263->add_dep(task108);
  residualq->add_task(task263);

  vector<IndexRange> I423_index = {closed_, active_, active_, active_};
  auto I423 = make_shared<Tensor>(I423_index);
  auto tensor264 = vector<shared_ptr<Tensor>>{I27, v2_, I423};
  auto task264 = make_shared<Task264>(tensor264, pindex);
  task193->add_dep(task264);
  task264->add_dep(task108);
  residualq->add_task(task264);

  auto tensor265 = vector<shared_ptr<Tensor>>{I423, Gamma132_(), t2};
  auto task265 = make_shared<Task265>(tensor265, pindex);
  task264->add_dep(task265);
  task265->add_dep(task108);
  residualq->add_task(task265);

  vector<IndexRange> I426_index = {closed_, active_, active_, active_};
  auto I426 = make_shared<Tensor>(I426_index);
  auto tensor266 = vector<shared_ptr<Tensor>>{I27, v2_, I426};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task193->add_dep(task266);
  task266->add_dep(task108);
  residualq->add_task(task266);

  auto tensor267 = vector<shared_ptr<Tensor>>{I426, Gamma132_(), t2};
  auto task267 = make_shared<Task267>(tensor267, pindex);
  task266->add_dep(task267);
  task267->add_dep(task108);
  residualq->add_task(task267);

  vector<IndexRange> I429_index = {closed_, active_};
  auto I429 = make_shared<Tensor>(I429_index);
  auto tensor268 = vector<shared_ptr<Tensor>>{I27, v2_, I429};
  auto task268 = make_shared<Task268>(tensor268, pindex);
  task193->add_dep(task268);
  task268->add_dep(task108);
  residualq->add_task(task268);

  auto tensor269 = vector<shared_ptr<Tensor>>{I429, Gamma10_(), t2};
  auto task269 = make_shared<Task269>(tensor269, pindex);
  task268->add_dep(task269);
  task269->add_dep(task108);
  residualq->add_task(task269);

  vector<IndexRange> I432_index = {closed_, active_};
  auto I432 = make_shared<Tensor>(I432_index);
  auto tensor270 = vector<shared_ptr<Tensor>>{I27, v2_, I432};
  auto task270 = make_shared<Task270>(tensor270, pindex);
  task193->add_dep(task270);
  task270->add_dep(task108);
  residualq->add_task(task270);

  auto tensor271 = vector<shared_ptr<Tensor>>{I432, Gamma10_(), t2};
  auto task271 = make_shared<Task271>(tensor271, pindex);
  task270->add_dep(task271);
  task271->add_dep(task108);
  residualq->add_task(task271);

  vector<IndexRange> I435_index = {closed_, closed_, active_, active_};
  auto I435 = make_shared<Tensor>(I435_index);
  auto tensor272 = vector<shared_ptr<Tensor>>{I27, t2, I435};
  auto task272 = make_shared<Task272>(tensor272, pindex);
  task193->add_dep(task272);
  task272->add_dep(task108);
  residualq->add_task(task272);

  vector<IndexRange> I436_index = {closed_, closed_, active_, active_};
  auto I436 = make_shared<Tensor>(I436_index);
  auto tensor273 = vector<shared_ptr<Tensor>>{I435, Gamma197_(), I436};
  auto task273 = make_shared<Task273>(tensor273, pindex);
  task272->add_dep(task273);
  task273->add_dep(task108);
  residualq->add_task(task273);

  auto tensor274 = vector<shared_ptr<Tensor>>{I436, v2_};
  auto task274 = make_shared<Task274>(tensor274, pindex);
  task273->add_dep(task274);
  task274->add_dep(task108);
  residualq->add_task(task274);

  auto tensor275 = vector<shared_ptr<Tensor>>{I435, Gamma0_(), v2_};
  auto task275 = make_shared<Task275>(tensor275, pindex);
  task272->add_dep(task275);
  task275->add_dep(task108);
  residualq->add_task(task275);

  vector<IndexRange> I438_index = {closed_, closed_, active_, active_};
  auto I438 = make_shared<Tensor>(I438_index);
  auto tensor276 = vector<shared_ptr<Tensor>>{I27, t2, I438};
  auto task276 = make_shared<Task276>(tensor276, pindex);
  task193->add_dep(task276);
  task276->add_dep(task108);
  residualq->add_task(task276);

  vector<IndexRange> I439_index = {closed_, closed_, active_, active_};
  auto I439 = make_shared<Tensor>(I439_index);
  auto tensor277 = vector<shared_ptr<Tensor>>{I438, Gamma197_(), I439};
  auto task277 = make_shared<Task277>(tensor277, pindex);
  task276->add_dep(task277);
  task277->add_dep(task108);
  residualq->add_task(task277);

  auto tensor278 = vector<shared_ptr<Tensor>>{I439, v2_};
  auto task278 = make_shared<Task278>(tensor278, pindex);
  task277->add_dep(task278);
  task278->add_dep(task108);
  residualq->add_task(task278);

  auto tensor279 = vector<shared_ptr<Tensor>>{I438, Gamma2_(), v2_};
  auto task279 = make_shared<Task279>(tensor279, pindex);
  task276->add_dep(task279);
  task279->add_dep(task108);
  residualq->add_task(task279);

  auto tensor280 = vector<shared_ptr<Tensor>>{I438, Gamma155_(), v2_};
  auto task280 = make_shared<Task280>(tensor280, pindex);
  task276->add_dep(task280);
  task280->add_dep(task108);
  residualq->add_task(task280);

  vector<IndexRange> I441_index = {closed_, closed_, active_, active_};
  auto I441 = make_shared<Tensor>(I441_index);
  auto tensor281 = vector<shared_ptr<Tensor>>{I27, t2, I441};
  auto task281 = make_shared<Task281>(tensor281, pindex);
  task193->add_dep(task281);
  task281->add_dep(task108);
  residualq->add_task(task281);

  vector<IndexRange> I442_index = {closed_, closed_, active_, active_};
  auto I442 = make_shared<Tensor>(I442_index);
  auto tensor282 = vector<shared_ptr<Tensor>>{I441, Gamma197_(), I442};
  auto task282 = make_shared<Task282>(tensor282, pindex);
  task281->add_dep(task282);
  task282->add_dep(task108);
  residualq->add_task(task282);

  auto tensor283 = vector<shared_ptr<Tensor>>{I442, v2_};
  auto task283 = make_shared<Task283>(tensor283, pindex);
  task282->add_dep(task283);
  task283->add_dep(task108);
  residualq->add_task(task283);

  auto tensor284 = vector<shared_ptr<Tensor>>{I441, Gamma2_(), v2_};
  auto task284 = make_shared<Task284>(tensor284, pindex);
  task281->add_dep(task284);
  task284->add_dep(task108);
  residualq->add_task(task284);

  auto tensor285 = vector<shared_ptr<Tensor>>{I441, Gamma155_(), v2_};
  auto task285 = make_shared<Task285>(tensor285, pindex);
  task281->add_dep(task285);
  task285->add_dep(task108);
  residualq->add_task(task285);

  vector<IndexRange> I444_index = {virt_, virt_, active_, active_};
  auto I444 = make_shared<Tensor>(I444_index);
  auto tensor286 = vector<shared_ptr<Tensor>>{I27, t2, I444};
  auto task286 = make_shared<Task286>(tensor286, pindex);
  task193->add_dep(task286);
  task286->add_dep(task108);
  residualq->add_task(task286);

  vector<IndexRange> I445_index = {virt_, virt_, active_, active_};
  auto I445 = make_shared<Tensor>(I445_index);
  auto tensor287 = vector<shared_ptr<Tensor>>{I444, Gamma197_(), I445};
  auto task287 = make_shared<Task287>(tensor287, pindex);
  task286->add_dep(task287);
  task287->add_dep(task108);
  residualq->add_task(task287);

  auto tensor288 = vector<shared_ptr<Tensor>>{I445, v2_};
  auto task288 = make_shared<Task288>(tensor288, pindex);
  task287->add_dep(task288);
  task288->add_dep(task108);
  residualq->add_task(task288);

  auto tensor289 = vector<shared_ptr<Tensor>>{I444, Gamma2_(), v2_};
  auto task289 = make_shared<Task289>(tensor289, pindex);
  task286->add_dep(task289);
  task289->add_dep(task108);
  residualq->add_task(task289);

  auto tensor290 = vector<shared_ptr<Tensor>>{I444, Gamma155_(), v2_};
  auto task290 = make_shared<Task290>(tensor290, pindex);
  task286->add_dep(task290);
  task290->add_dep(task108);
  residualq->add_task(task290);

  vector<IndexRange> I447_index = {closed_, closed_, active_, active_};
  auto I447 = make_shared<Tensor>(I447_index);
  auto tensor291 = vector<shared_ptr<Tensor>>{I27, t2, I447};
  auto task291 = make_shared<Task291>(tensor291, pindex);
  task193->add_dep(task291);
  task291->add_dep(task108);
  residualq->add_task(task291);

  vector<IndexRange> I448_index = {closed_, closed_, active_, active_};
  auto I448 = make_shared<Tensor>(I448_index);
  auto tensor292 = vector<shared_ptr<Tensor>>{I447, Gamma197_(), I448};
  auto task292 = make_shared<Task292>(tensor292, pindex);
  task291->add_dep(task292);
  task292->add_dep(task108);
  residualq->add_task(task292);

  auto tensor293 = vector<shared_ptr<Tensor>>{I448, v2_};
  auto task293 = make_shared<Task293>(tensor293, pindex);
  task292->add_dep(task293);
  task293->add_dep(task108);
  residualq->add_task(task293);

  auto tensor294 = vector<shared_ptr<Tensor>>{I447, Gamma2_(), v2_};
  auto task294 = make_shared<Task294>(tensor294, pindex);
  task291->add_dep(task294);
  task294->add_dep(task108);
  residualq->add_task(task294);

  auto tensor295 = vector<shared_ptr<Tensor>>{I447, Gamma155_(), v2_};
  auto task295 = make_shared<Task295>(tensor295, pindex);
  task291->add_dep(task295);
  task295->add_dep(task108);
  residualq->add_task(task295);

  vector<IndexRange> I450_index = {virt_, virt_, active_, active_};
  auto I450 = make_shared<Tensor>(I450_index);
  auto tensor296 = vector<shared_ptr<Tensor>>{I27, t2, I450};
  auto task296 = make_shared<Task296>(tensor296, pindex);
  task193->add_dep(task296);
  task296->add_dep(task108);
  residualq->add_task(task296);

  vector<IndexRange> I451_index = {virt_, virt_, active_, active_};
  auto I451 = make_shared<Tensor>(I451_index);
  auto tensor297 = vector<shared_ptr<Tensor>>{I450, Gamma197_(), I451};
  auto task297 = make_shared<Task297>(tensor297, pindex);
  task296->add_dep(task297);
  task297->add_dep(task108);
  residualq->add_task(task297);

  auto tensor298 = vector<shared_ptr<Tensor>>{I451, v2_};
  auto task298 = make_shared<Task298>(tensor298, pindex);
  task297->add_dep(task298);
  task298->add_dep(task108);
  residualq->add_task(task298);

  auto tensor299 = vector<shared_ptr<Tensor>>{I450, Gamma0_(), v2_};
  auto task299 = make_shared<Task299>(tensor299, pindex);
  task296->add_dep(task299);
  task299->add_dep(task108);
  residualq->add_task(task299);

  vector<IndexRange> I537_index = {closed_, active_, active_, active_};
  auto I537 = make_shared<Tensor>(I537_index);
  auto tensor300 = vector<shared_ptr<Tensor>>{I27, t2, I537};
  auto task300 = make_shared<Task300>(tensor300, pindex);
  task193->add_dep(task300);
  task300->add_dep(task108);
  residualq->add_task(task300);

  auto tensor301 = vector<shared_ptr<Tensor>>{I537, Gamma176_(), v2_};
  auto task301 = make_shared<Task301>(tensor301, pindex);
  task300->add_dep(task301);
  task301->add_dep(task108);
  residualq->add_task(task301);

  auto tensor302 = vector<shared_ptr<Tensor>>{I537, Gamma178_(), v2_};
  auto task302 = make_shared<Task302>(tensor302, pindex);
  task300->add_dep(task302);
  task302->add_dep(task108);
  residualq->add_task(task302);

  vector<IndexRange> I540_index = {closed_, active_, active_, active_};
  auto I540 = make_shared<Tensor>(I540_index);
  auto tensor303 = vector<shared_ptr<Tensor>>{I27, t2, I540};
  auto task303 = make_shared<Task303>(tensor303, pindex);
  task193->add_dep(task303);
  task303->add_dep(task108);
  residualq->add_task(task303);

  auto tensor304 = vector<shared_ptr<Tensor>>{I540, Gamma132_(), v2_};
  auto task304 = make_shared<Task304>(tensor304, pindex);
  task303->add_dep(task304);
  task304->add_dep(task108);
  residualq->add_task(task304);

  auto tensor305 = vector<shared_ptr<Tensor>>{I540, Gamma179_(), v2_};
  auto task305 = make_shared<Task305>(tensor305, pindex);
  task303->add_dep(task305);
  task305->add_dep(task108);
  residualq->add_task(task305);

  vector<IndexRange> I549_index = {virt_, closed_, active_, active_};
  auto I549 = make_shared<Tensor>(I549_index);
  auto tensor306 = vector<shared_ptr<Tensor>>{I27, v2_, I549};
  auto task306 = make_shared<Task306>(tensor306, pindex);
  task193->add_dep(task306);
  task306->add_dep(task108);
  residualq->add_task(task306);

  vector<IndexRange> I550_index = {active_, virt_, closed_, active_};
  auto I550 = make_shared<Tensor>(I550_index);
  auto tensor307 = vector<shared_ptr<Tensor>>{I549, Gamma10_(), I550};
  auto task307 = make_shared<Task307>(tensor307, pindex);
  task306->add_dep(task307);
  task307->add_dep(task108);
  residualq->add_task(task307);

  auto tensor308 = vector<shared_ptr<Tensor>>{I550, t2};
  auto task308 = make_shared<Task308>(tensor308, pindex);
  task307->add_dep(task308);
  task308->add_dep(task108);
  residualq->add_task(task308);

  vector<IndexRange> I552_index = {virt_, closed_, active_, active_};
  auto I552 = make_shared<Tensor>(I552_index);
  auto tensor309 = vector<shared_ptr<Tensor>>{I27, v2_, I552};
  auto task309 = make_shared<Task309>(tensor309, pindex);
  task193->add_dep(task309);
  task309->add_dep(task108);
  residualq->add_task(task309);

  auto tensor310 = vector<shared_ptr<Tensor>>{I552, Gamma18_(), t2};
  auto task310 = make_shared<Task310>(tensor310, pindex);
  task309->add_dep(task310);
  task310->add_dep(task108);
  residualq->add_task(task310);

  auto tensor311 = vector<shared_ptr<Tensor>>{I552, Gamma10_(), t2};
  auto task311 = make_shared<Task311>(tensor311, pindex);
  task309->add_dep(task311);
  task311->add_dep(task108);
  residualq->add_task(task311);

  vector<IndexRange> I555_index = {virt_, closed_, active_, active_};
  auto I555 = make_shared<Tensor>(I555_index);
  auto tensor312 = vector<shared_ptr<Tensor>>{I27, v2_, I555};
  auto task312 = make_shared<Task312>(tensor312, pindex);
  task193->add_dep(task312);
  task312->add_dep(task108);
  residualq->add_task(task312);

  auto tensor313 = vector<shared_ptr<Tensor>>{I555, Gamma18_(), t2};
  auto task313 = make_shared<Task313>(tensor313, pindex);
  task312->add_dep(task313);
  task313->add_dep(task108);
  residualq->add_task(task313);

  auto tensor314 = vector<shared_ptr<Tensor>>{I555, Gamma10_(), t2};
  auto task314 = make_shared<Task314>(tensor314, pindex);
  task312->add_dep(task314);
  task314->add_dep(task108);
  residualq->add_task(task314);

  vector<IndexRange> I558_index = {virt_, closed_, active_, active_};
  auto I558 = make_shared<Tensor>(I558_index);
  auto tensor315 = vector<shared_ptr<Tensor>>{I27, v2_, I558};
  auto task315 = make_shared<Task315>(tensor315, pindex);
  task193->add_dep(task315);
  task315->add_dep(task108);
  residualq->add_task(task315);

  auto tensor316 = vector<shared_ptr<Tensor>>{I558, Gamma18_(), t2};
  auto task316 = make_shared<Task316>(tensor316, pindex);
  task315->add_dep(task316);
  task316->add_dep(task108);
  residualq->add_task(task316);

  auto tensor317 = vector<shared_ptr<Tensor>>{I558, Gamma10_(), t2};
  auto task317 = make_shared<Task317>(tensor317, pindex);
  task315->add_dep(task317);
  task317->add_dep(task108);
  residualq->add_task(task317);

  vector<IndexRange> I561_index = {virt_, closed_, active_, active_};
  auto I561 = make_shared<Tensor>(I561_index);
  auto tensor318 = vector<shared_ptr<Tensor>>{I27, v2_, I561};
  auto task318 = make_shared<Task318>(tensor318, pindex);
  task193->add_dep(task318);
  task318->add_dep(task108);
  residualq->add_task(task318);

  auto tensor319 = vector<shared_ptr<Tensor>>{I561, Gamma18_(), t2};
  auto task319 = make_shared<Task319>(tensor319, pindex);
  task318->add_dep(task319);
  task319->add_dep(task108);
  residualq->add_task(task319);

  auto tensor320 = vector<shared_ptr<Tensor>>{I561, Gamma10_(), t2};
  auto task320 = make_shared<Task320>(tensor320, pindex);
  task318->add_dep(task320);
  task320->add_dep(task108);
  residualq->add_task(task320);

  vector<IndexRange> I564_index = {virt_, closed_, active_, active_};
  auto I564 = make_shared<Tensor>(I564_index);
  auto tensor321 = vector<shared_ptr<Tensor>>{I27, v2_, I564};
  auto task321 = make_shared<Task321>(tensor321, pindex);
  task193->add_dep(task321);
  task321->add_dep(task108);
  residualq->add_task(task321);

  vector<IndexRange> I565_index = {active_, virt_, closed_, active_};
  auto I565 = make_shared<Tensor>(I565_index);
  auto tensor322 = vector<shared_ptr<Tensor>>{I564, Gamma10_(), I565};
  auto task322 = make_shared<Task322>(tensor322, pindex);
  task321->add_dep(task322);
  task322->add_dep(task108);
  residualq->add_task(task322);

  auto tensor323 = vector<shared_ptr<Tensor>>{I565, t2};
  auto task323 = make_shared<Task323>(tensor323, pindex);
  task322->add_dep(task323);
  task323->add_dep(task108);
  residualq->add_task(task323);

  vector<IndexRange> I567_index = {closed_, active_, active_, active_};
  auto I567 = make_shared<Tensor>(I567_index);
  auto tensor324 = vector<shared_ptr<Tensor>>{I27, t2, I567};
  auto task324 = make_shared<Task324>(tensor324, pindex);
  task193->add_dep(task324);
  task324->add_dep(task108);
  residualq->add_task(task324);

  auto tensor325 = vector<shared_ptr<Tensor>>{I567, Gamma132_(), v2_};
  auto task325 = make_shared<Task325>(tensor325, pindex);
  task324->add_dep(task325);
  task325->add_dep(task108);
  residualq->add_task(task325);

  auto tensor326 = vector<shared_ptr<Tensor>>{I567, Gamma179_(), v2_};
  auto task326 = make_shared<Task326>(tensor326, pindex);
  task324->add_dep(task326);
  task326->add_dep(task108);
  residualq->add_task(task326);

  vector<IndexRange> I570_index = {closed_, active_, active_, active_};
  auto I570 = make_shared<Tensor>(I570_index);
  auto tensor327 = vector<shared_ptr<Tensor>>{I27, t2, I570};
  auto task327 = make_shared<Task327>(tensor327, pindex);
  task193->add_dep(task327);
  task327->add_dep(task108);
  residualq->add_task(task327);

  auto tensor328 = vector<shared_ptr<Tensor>>{I570, Gamma132_(), v2_};
  auto task328 = make_shared<Task328>(tensor328, pindex);
  task327->add_dep(task328);
  task328->add_dep(task108);
  residualq->add_task(task328);

  auto tensor329 = vector<shared_ptr<Tensor>>{I570, Gamma179_(), v2_};
  auto task329 = make_shared<Task329>(tensor329, pindex);
  task327->add_dep(task329);
  task329->add_dep(task108);
  residualq->add_task(task329);

  vector<IndexRange> I597_index = {virt_, active_, active_, active_};
  auto I597 = make_shared<Tensor>(I597_index);
  auto tensor330 = vector<shared_ptr<Tensor>>{I27, v2_, I597};
  auto task330 = make_shared<Task330>(tensor330, pindex);
  task193->add_dep(task330);
  task330->add_dep(task108);
  residualq->add_task(task330);

  auto tensor331 = vector<shared_ptr<Tensor>>{I597, Gamma196_(), t2};
  auto task331 = make_shared<Task331>(tensor331, pindex);
  task330->add_dep(task331);
  task331->add_dep(task108);
  residualq->add_task(task331);

  vector<IndexRange> I642_index = {closed_, virt_, active_, active_};
  auto I642 = make_shared<Tensor>(I642_index);
  auto tensor332 = vector<shared_ptr<Tensor>>{I27, t2, I642};
  auto task332 = make_shared<Task332>(tensor332, pindex);
  task193->add_dep(task332);
  task332->add_dep(task108);
  residualq->add_task(task332);

  auto tensor333 = vector<shared_ptr<Tensor>>{I642, Gamma18_(), v2_};
  auto task333 = make_shared<Task333>(tensor333, pindex);
  task332->add_dep(task333);
  task333->add_dep(task108);
  residualq->add_task(task333);

  vector<IndexRange> I645_index = {closed_, virt_, active_, active_};
  auto I645 = make_shared<Tensor>(I645_index);
  auto tensor334 = vector<shared_ptr<Tensor>>{I27, t2, I645};
  auto task334 = make_shared<Task334>(tensor334, pindex);
  task193->add_dep(task334);
  task334->add_dep(task108);
  residualq->add_task(task334);

  auto tensor335 = vector<shared_ptr<Tensor>>{I645, Gamma10_(), v2_};
  auto task335 = make_shared<Task335>(tensor335, pindex);
  task334->add_dep(task335);
  task335->add_dep(task108);
  residualq->add_task(task335);

  vector<IndexRange> I648_index = {closed_, virt_, active_, active_};
  auto I648 = make_shared<Tensor>(I648_index);
  auto tensor336 = vector<shared_ptr<Tensor>>{I27, t2, I648};
  auto task336 = make_shared<Task336>(tensor336, pindex);
  task193->add_dep(task336);
  task336->add_dep(task108);
  residualq->add_task(task336);

  auto tensor337 = vector<shared_ptr<Tensor>>{I648, Gamma18_(), v2_};
  auto task337 = make_shared<Task337>(tensor337, pindex);
  task336->add_dep(task337);
  task337->add_dep(task108);
  residualq->add_task(task337);

  vector<IndexRange> I651_index = {closed_, virt_, active_, active_};
  auto I651 = make_shared<Tensor>(I651_index);
  auto tensor338 = vector<shared_ptr<Tensor>>{I27, t2, I651};
  auto task338 = make_shared<Task338>(tensor338, pindex);
  task193->add_dep(task338);
  task338->add_dep(task108);
  residualq->add_task(task338);

  auto tensor339 = vector<shared_ptr<Tensor>>{I651, Gamma18_(), v2_};
  auto task339 = make_shared<Task339>(tensor339, pindex);
  task338->add_dep(task339);
  task339->add_dep(task108);
  residualq->add_task(task339);

  vector<IndexRange> I1681_index = {closed_, virt_, closed_, active_};
  auto I1681 = make_shared<Tensor>(I1681_index);
  auto tensor340 = vector<shared_ptr<Tensor>>{I27, Gamma552_(), I1681};
  auto task340 = make_shared<Task340>(tensor340, pindex);
  task193->add_dep(task340);
  task340->add_dep(task108);
  residualq->add_task(task340);

  auto tensor341 = vector<shared_ptr<Tensor>>{I1681, t2};
  auto task341 = make_shared<Task341>(tensor341, pindex);
  task340->add_dep(task341);
  task341->add_dep(task108);
  residualq->add_task(task341);

  vector<IndexRange> I1685_index = {closed_, virt_, closed_, active_};
  auto I1685 = make_shared<Tensor>(I1685_index);
  auto tensor342 = vector<shared_ptr<Tensor>>{I27, Gamma554_(), I1685};
  auto task342 = make_shared<Task342>(tensor342, pindex);
  task193->add_dep(task342);
  task342->add_dep(task108);
  residualq->add_task(task342);

  auto tensor343 = vector<shared_ptr<Tensor>>{I1685, t2};
  auto task343 = make_shared<Task343>(tensor343, pindex);
  task342->add_dep(task343);
  task343->add_dep(task108);
  residualq->add_task(task343);
}

#endif
