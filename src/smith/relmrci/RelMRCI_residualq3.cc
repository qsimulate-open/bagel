//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_residualq3.cc
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


#include <src/smith/relmrci/RelMRCI.h>
#include <src/smith/relmrci/RelMRCI_tasks6.h>
#include <src/smith/relmrci/RelMRCI_tasks7.h>
#include <src/smith/relmrci/RelMRCI_tasks8.h>
#include <src/smith/relmrci/RelMRCI_tasks9.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void RelMRCI::RelMRCI::make_residualq3(shared_ptr<Queue> residualq, shared_ptr<Task> task83, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I63_index = {closed_, active_, active_, virt_};
  auto I63 = make_shared<Tensor>(I63_index);
  auto tensor297 = vector<shared_ptr<Tensor>>{r, I63};
  auto task297 = make_shared<Task297>(tensor297, pindex);
  task297->add_dep(task83);
  residualq->add_task(task297);

  vector<IndexRange> I64_index = {closed_, active_, active_, active_};
  auto I64 = make_shared<Tensor>(I64_index);
  auto tensor298 = vector<shared_ptr<Tensor>>{I63, h1_, I64};
  auto task298 = make_shared<Task298>(tensor298, pindex);
  task297->add_dep(task298);
  task298->add_dep(task83);
  residualq->add_task(task298);

  auto tensor299 = vector<shared_ptr<Tensor>>{I64, Gamma4_(), t2};
  auto task299 = make_shared<Task299>(tensor299, pindex);
  task298->add_dep(task299);
  task299->add_dep(task83);
  residualq->add_task(task299);

  vector<IndexRange> I67_index = {active_, virt_, closed_, active_};
  auto I67 = make_shared<Tensor>(I67_index);
  auto tensor300 = vector<shared_ptr<Tensor>>{I63, Gamma5_(), I67};
  auto task300 = make_shared<Task300>(tensor300, pindex);
  task297->add_dep(task300);
  task300->add_dep(task83);
  residualq->add_task(task300);

  auto tensor301 = vector<shared_ptr<Tensor>>{I67, t2, h1_};
  auto task301 = make_shared<Task301>(tensor301, pindex);
  task300->add_dep(task301);
  task301->add_dep(task83);
  residualq->add_task(task301);

  auto tensor302 = vector<shared_ptr<Tensor>>{I67, t2, h1_};
  auto task302 = make_shared<Task302>(tensor302, pindex);
  task300->add_dep(task302);
  task302->add_dep(task83);
  residualq->add_task(task302);

  auto tensor303 = vector<shared_ptr<Tensor>>{I67, t2, v2_};
  auto task303 = make_shared<Task303>(tensor303, pindex);
  task300->add_dep(task303);
  task303->add_dep(task83);
  residualq->add_task(task303);

  auto tensor304 = vector<shared_ptr<Tensor>>{I67, t2, v2_};
  auto task304 = make_shared<Task304>(tensor304, pindex);
  task300->add_dep(task304);
  task304->add_dep(task83);
  residualq->add_task(task304);

  auto tensor305 = vector<shared_ptr<Tensor>>{I67, t2, v2_};
  auto task305 = make_shared<Task305>(tensor305, pindex);
  task300->add_dep(task305);
  task305->add_dep(task83);
  residualq->add_task(task305);

  auto tensor306 = vector<shared_ptr<Tensor>>{I67, t2, v2_};
  auto task306 = make_shared<Task306>(tensor306, pindex);
  task300->add_dep(task306);
  task306->add_dep(task83);
  residualq->add_task(task306);

  auto tensor307 = vector<shared_ptr<Tensor>>{I67, t2, v2_};
  auto task307 = make_shared<Task307>(tensor307, pindex);
  task300->add_dep(task307);
  task307->add_dep(task83);
  residualq->add_task(task307);

  auto tensor308 = vector<shared_ptr<Tensor>>{I67, t2, v2_};
  auto task308 = make_shared<Task308>(tensor308, pindex);
  task300->add_dep(task308);
  task308->add_dep(task83);
  residualq->add_task(task308);

  auto tensor309 = vector<shared_ptr<Tensor>>{I67, t2, v2_};
  auto task309 = make_shared<Task309>(tensor309, pindex);
  task300->add_dep(task309);
  task309->add_dep(task83);
  residualq->add_task(task309);

  auto tensor310 = vector<shared_ptr<Tensor>>{I67, t2, v2_};
  auto task310 = make_shared<Task310>(tensor310, pindex);
  task300->add_dep(task310);
  task310->add_dep(task83);
  residualq->add_task(task310);

  vector<IndexRange> I73_index = {closed_, virt_, active_, active_};
  auto I73 = make_shared<Tensor>(I73_index);
  auto tensor311 = vector<shared_ptr<Tensor>>{I63, Gamma24_(), I73};
  auto task311 = make_shared<Task311>(tensor311, pindex);
  task297->add_dep(task311);
  task311->add_dep(task83);
  residualq->add_task(task311);

  auto tensor312 = vector<shared_ptr<Tensor>>{I73, t2, h1_};
  auto task312 = make_shared<Task312>(tensor312, pindex);
  task311->add_dep(task312);
  task312->add_dep(task83);
  residualq->add_task(task312);

  auto tensor313 = vector<shared_ptr<Tensor>>{I73, t2, h1_};
  auto task313 = make_shared<Task313>(tensor313, pindex);
  task311->add_dep(task313);
  task313->add_dep(task83);
  residualq->add_task(task313);

  vector<IndexRange> I89_index = {active_, virt_, closed_, virt_};
  auto I89 = make_shared<Tensor>(I89_index);
  auto tensor314 = vector<shared_ptr<Tensor>>{I73, h1_, I89};
  auto task314 = make_shared<Task314>(tensor314, pindex);
  task311->add_dep(task314);
  task314->add_dep(task83);
  residualq->add_task(task314);

  auto tensor315 = vector<shared_ptr<Tensor>>{I89, t2};
  auto task315 = make_shared<Task315>(tensor315, pindex);
  task314->add_dep(task315);
  task315->add_dep(task83);
  residualq->add_task(task315);

  auto tensor316 = vector<shared_ptr<Tensor>>{I73, t2, v2_};
  auto task316 = make_shared<Task316>(tensor316, pindex);
  task311->add_dep(task316);
  task316->add_dep(task83);
  residualq->add_task(task316);

  auto tensor317 = vector<shared_ptr<Tensor>>{I73, t2, v2_};
  auto task317 = make_shared<Task317>(tensor317, pindex);
  task311->add_dep(task317);
  task317->add_dep(task83);
  residualq->add_task(task317);

  vector<IndexRange> I631_index = {virt_, closed_, active_, active_};
  auto I631 = make_shared<Tensor>(I631_index);
  auto tensor318 = vector<shared_ptr<Tensor>>{I73, t2, I631};
  auto task318 = make_shared<Task318>(tensor318, pindex);
  task311->add_dep(task318);
  task318->add_dep(task83);
  residualq->add_task(task318);

  auto tensor319 = vector<shared_ptr<Tensor>>{I631, v2_};
  auto task319 = make_shared<Task319>(tensor319, pindex);
  task318->add_dep(task319);
  task319->add_dep(task83);
  residualq->add_task(task319);

  vector<IndexRange> I634_index = {virt_, closed_, active_, active_};
  auto I634 = make_shared<Tensor>(I634_index);
  auto tensor320 = vector<shared_ptr<Tensor>>{I73, t2, I634};
  auto task320 = make_shared<Task320>(tensor320, pindex);
  task311->add_dep(task320);
  task320->add_dep(task83);
  residualq->add_task(task320);

  auto tensor321 = vector<shared_ptr<Tensor>>{I634, v2_};
  auto task321 = make_shared<Task321>(tensor321, pindex);
  task320->add_dep(task321);
  task321->add_dep(task83);
  residualq->add_task(task321);

  auto tensor322 = vector<shared_ptr<Tensor>>{I73, t2, v2_};
  auto task322 = make_shared<Task322>(tensor322, pindex);
  task311->add_dep(task322);
  task322->add_dep(task83);
  residualq->add_task(task322);

  auto tensor323 = vector<shared_ptr<Tensor>>{I73, t2, v2_};
  auto task323 = make_shared<Task323>(tensor323, pindex);
  task311->add_dep(task323);
  task323->add_dep(task83);
  residualq->add_task(task323);

  vector<IndexRange> I679_index = {virt_, active_, closed_, closed_};
  auto I679 = make_shared<Tensor>(I679_index);
  auto tensor324 = vector<shared_ptr<Tensor>>{I73, t2, I679};
  auto task324 = make_shared<Task324>(tensor324, pindex);
  task311->add_dep(task324);
  task324->add_dep(task83);
  residualq->add_task(task324);

  auto tensor325 = vector<shared_ptr<Tensor>>{I679, v2_};
  auto task325 = make_shared<Task325>(tensor325, pindex);
  task324->add_dep(task325);
  task325->add_dep(task83);
  residualq->add_task(task325);

  vector<IndexRange> I682_index = {virt_, active_, closed_, closed_};
  auto I682 = make_shared<Tensor>(I682_index);
  auto tensor326 = vector<shared_ptr<Tensor>>{I73, t2, I682};
  auto task326 = make_shared<Task326>(tensor326, pindex);
  task311->add_dep(task326);
  task326->add_dep(task83);
  residualq->add_task(task326);

  auto tensor327 = vector<shared_ptr<Tensor>>{I682, v2_};
  auto task327 = make_shared<Task327>(tensor327, pindex);
  task326->add_dep(task327);
  task327->add_dep(task83);
  residualq->add_task(task327);

  auto tensor328 = vector<shared_ptr<Tensor>>{I73, t2, v2_};
  auto task328 = make_shared<Task328>(tensor328, pindex);
  task311->add_dep(task328);
  task328->add_dep(task83);
  residualq->add_task(task328);

  auto tensor329 = vector<shared_ptr<Tensor>>{I73, t2, v2_};
  auto task329 = make_shared<Task329>(tensor329, pindex);
  task311->add_dep(task329);
  task329->add_dep(task83);
  residualq->add_task(task329);

  vector<IndexRange> I79_index = {virt_, active_, active_, active_};
  auto I79 = make_shared<Tensor>(I79_index);
  auto tensor330 = vector<shared_ptr<Tensor>>{I63, h1_, I79};
  auto task330 = make_shared<Task330>(tensor330, pindex);
  task297->add_dep(task330);
  task330->add_dep(task83);
  residualq->add_task(task330);

  auto tensor331 = vector<shared_ptr<Tensor>>{I79, Gamma26_(), t2};
  auto task331 = make_shared<Task331>(tensor331, pindex);
  task330->add_dep(task331);
  task331->add_dep(task83);
  residualq->add_task(task331);

  vector<IndexRange> I82_index = {closed_, virt_};
  auto I82 = make_shared<Tensor>(I82_index);
  auto tensor332 = vector<shared_ptr<Tensor>>{I63, Gamma27_(), I82};
  auto task332 = make_shared<Task332>(tensor332, pindex);
  task297->add_dep(task332);
  task332->add_dep(task83);
  residualq->add_task(task332);

  auto tensor333 = vector<shared_ptr<Tensor>>{I82, t2, h1_};
  auto task333 = make_shared<Task333>(tensor333, pindex);
  task332->add_dep(task333);
  task333->add_dep(task83);
  residualq->add_task(task333);

  auto tensor334 = vector<shared_ptr<Tensor>>{I82, t2, h1_};
  auto task334 = make_shared<Task334>(tensor334, pindex);
  task332->add_dep(task334);
  task334->add_dep(task83);
  residualq->add_task(task334);

  auto tensor335 = vector<shared_ptr<Tensor>>{I82, t2, v2_};
  auto task335 = make_shared<Task335>(tensor335, pindex);
  task332->add_dep(task335);
  task335->add_dep(task83);
  residualq->add_task(task335);

  auto tensor336 = vector<shared_ptr<Tensor>>{I82, t2, v2_};
  auto task336 = make_shared<Task336>(tensor336, pindex);
  task332->add_dep(task336);
  task336->add_dep(task83);
  residualq->add_task(task336);

  auto tensor337 = vector<shared_ptr<Tensor>>{I82, t2, v2_};
  auto task337 = make_shared<Task337>(tensor337, pindex);
  task332->add_dep(task337);
  task337->add_dep(task83);
  residualq->add_task(task337);

  auto tensor338 = vector<shared_ptr<Tensor>>{I82, t2, v2_};
  auto task338 = make_shared<Task338>(tensor338, pindex);
  task332->add_dep(task338);
  task338->add_dep(task83);
  residualq->add_task(task338);

  vector<IndexRange> I543_index = {closed_, closed_, active_, active_, active_, active_};
  auto I543 = make_shared<Tensor>(I543_index);
  auto tensor339 = vector<shared_ptr<Tensor>>{I63, v2_, I543};
  auto task339 = make_shared<Task339>(tensor339, pindex);
  task297->add_dep(task339);
  task339->add_dep(task83);
  residualq->add_task(task339);

  auto tensor340 = vector<shared_ptr<Tensor>>{I543, Gamma84_(), t2};
  auto task340 = make_shared<Task340>(tensor340, pindex);
  task339->add_dep(task340);
  task340->add_dep(task83);
  residualq->add_task(task340);

  vector<IndexRange> I546_index = {closed_, active_, active_, active_, active_, active_};
  auto I546 = make_shared<Tensor>(I546_index);
  auto tensor341 = vector<shared_ptr<Tensor>>{I63, v2_, I546};
  auto task341 = make_shared<Task341>(tensor341, pindex);
  task297->add_dep(task341);
  task341->add_dep(task83);
  residualq->add_task(task341);

  auto tensor342 = vector<shared_ptr<Tensor>>{I546, Gamma179_(), t2};
  auto task342 = make_shared<Task342>(tensor342, pindex);
  task341->add_dep(task342);
  task342->add_dep(task83);
  residualq->add_task(task342);

  vector<IndexRange> I549_index = {closed_, active_, active_, active_, active_, active_};
  auto I549 = make_shared<Tensor>(I549_index);
  auto tensor343 = vector<shared_ptr<Tensor>>{I63, v2_, I549};
  auto task343 = make_shared<Task343>(tensor343, pindex);
  task297->add_dep(task343);
  task343->add_dep(task83);
  residualq->add_task(task343);

  auto tensor344 = vector<shared_ptr<Tensor>>{I549, Gamma77_(), t2};
  auto task344 = make_shared<Task344>(tensor344, pindex);
  task343->add_dep(task344);
  task344->add_dep(task83);
  residualq->add_task(task344);

  vector<IndexRange> I552_index = {closed_, active_, active_, active_};
  auto I552 = make_shared<Tensor>(I552_index);
  auto tensor345 = vector<shared_ptr<Tensor>>{I63, v2_, I552};
  auto task345 = make_shared<Task345>(tensor345, pindex);
  task297->add_dep(task345);
  task345->add_dep(task83);
  residualq->add_task(task345);

  auto tensor346 = vector<shared_ptr<Tensor>>{I552, Gamma4_(), t2};
  auto task346 = make_shared<Task346>(tensor346, pindex);
  task345->add_dep(task346);
  task346->add_dep(task83);
  residualq->add_task(task346);

  vector<IndexRange> I555_index = {closed_, active_, active_, active_};
  auto I555 = make_shared<Tensor>(I555_index);
  auto tensor347 = vector<shared_ptr<Tensor>>{I63, v2_, I555};
  auto task347 = make_shared<Task347>(tensor347, pindex);
  task297->add_dep(task347);
  task347->add_dep(task83);
  residualq->add_task(task347);

  auto tensor348 = vector<shared_ptr<Tensor>>{I555, Gamma4_(), t2};
  auto task348 = make_shared<Task348>(tensor348, pindex);
  task347->add_dep(task348);
  task348->add_dep(task83);
  residualq->add_task(task348);

  vector<IndexRange> I558_index = {closed_, active_, active_, active_};
  auto I558 = make_shared<Tensor>(I558_index);
  auto tensor349 = vector<shared_ptr<Tensor>>{I63, t2, I558};
  auto task349 = make_shared<Task349>(tensor349, pindex);
  task297->add_dep(task349);
  task349->add_dep(task83);
  residualq->add_task(task349);

  auto tensor350 = vector<shared_ptr<Tensor>>{I558, Gamma183_(), v2_};
  auto task350 = make_shared<Task350>(tensor350, pindex);
  task349->add_dep(task350);
  task350->add_dep(task83);
  residualq->add_task(task350);

  auto tensor351 = vector<shared_ptr<Tensor>>{I558, Gamma81_(), v2_};
  auto task351 = make_shared<Task351>(tensor351, pindex);
  task349->add_dep(task351);
  task351->add_dep(task83);
  residualq->add_task(task351);

  vector<IndexRange> I561_index = {closed_, active_, active_, active_};
  auto I561 = make_shared<Tensor>(I561_index);
  auto tensor352 = vector<shared_ptr<Tensor>>{I63, t2, I561};
  auto task352 = make_shared<Task352>(tensor352, pindex);
  task297->add_dep(task352);
  task352->add_dep(task83);
  residualq->add_task(task352);

  auto tensor353 = vector<shared_ptr<Tensor>>{I561, Gamma183_(), v2_};
  auto task353 = make_shared<Task353>(tensor353, pindex);
  task352->add_dep(task353);
  task353->add_dep(task83);
  residualq->add_task(task353);

  auto tensor354 = vector<shared_ptr<Tensor>>{I561, Gamma81_(), v2_};
  auto task354 = make_shared<Task354>(tensor354, pindex);
  task352->add_dep(task354);
  task354->add_dep(task83);
  residualq->add_task(task354);

  vector<IndexRange> I588_index = {closed_, closed_, active_, active_, active_, active_};
  auto I588 = make_shared<Tensor>(I588_index);
  auto tensor355 = vector<shared_ptr<Tensor>>{I63, t2, I588};
  auto task355 = make_shared<Task355>(tensor355, pindex);
  task297->add_dep(task355);
  task355->add_dep(task83);
  residualq->add_task(task355);

  vector<IndexRange> I589_index = {closed_, closed_, active_, active_};
  auto I589 = make_shared<Tensor>(I589_index);
  auto tensor356 = vector<shared_ptr<Tensor>>{I588, Gamma193_(), I589};
  auto task356 = make_shared<Task356>(tensor356, pindex);
  task355->add_dep(task356);
  task356->add_dep(task83);
  residualq->add_task(task356);

  auto tensor357 = vector<shared_ptr<Tensor>>{I589, v2_};
  auto task357 = make_shared<Task357>(tensor357, pindex);
  task356->add_dep(task357);
  task357->add_dep(task83);
  residualq->add_task(task357);

  auto tensor358 = vector<shared_ptr<Tensor>>{I588, Gamma4_(), v2_};
  auto task358 = make_shared<Task358>(tensor358, pindex);
  task355->add_dep(task358);
  task358->add_dep(task83);
  residualq->add_task(task358);

  vector<IndexRange> I591_index = {virt_, active_, active_, closed_, active_, active_};
  auto I591 = make_shared<Tensor>(I591_index);
  auto tensor359 = vector<shared_ptr<Tensor>>{I63, Gamma193_(), I591};
  auto task359 = make_shared<Task359>(tensor359, pindex);
  task297->add_dep(task359);
  task359->add_dep(task83);
  residualq->add_task(task359);

  vector<IndexRange> I592_index = {virt_, virt_, active_, active_};
  auto I592 = make_shared<Tensor>(I592_index);
  auto tensor360 = vector<shared_ptr<Tensor>>{I591, t2, I592};
  auto task360 = make_shared<Task360>(tensor360, pindex);
  task359->add_dep(task360);
  task360->add_dep(task83);
  residualq->add_task(task360);

  auto tensor361 = vector<shared_ptr<Tensor>>{I592, v2_};
  auto task361 = make_shared<Task361>(tensor361, pindex);
  task360->add_dep(task361);
  task361->add_dep(task83);
  residualq->add_task(task361);

  vector<IndexRange> I597_index = {active_, active_, virt_, closed_, active_, active_};
  auto I597 = make_shared<Tensor>(I597_index);
  auto tensor362 = vector<shared_ptr<Tensor>>{I63, Gamma4_(), I597};
  auto task362 = make_shared<Task362>(tensor362, pindex);
  task297->add_dep(task362);
  task362->add_dep(task83);
  residualq->add_task(task362);

  auto tensor363 = vector<shared_ptr<Tensor>>{I597, t2, v2_};
  auto task363 = make_shared<Task363>(tensor363, pindex);
  task362->add_dep(task363);
  task363->add_dep(task83);
  residualq->add_task(task363);

  vector<IndexRange> I618_index = {closed_, active_, active_, active_, active_, active_};
  auto I618 = make_shared<Tensor>(I618_index);
  auto tensor364 = vector<shared_ptr<Tensor>>{I63, t2, I618};
  auto task364 = make_shared<Task364>(tensor364, pindex);
  task297->add_dep(task364);
  task364->add_dep(task83);
  residualq->add_task(task364);

  auto tensor365 = vector<shared_ptr<Tensor>>{I618, Gamma203_(), v2_};
  auto task365 = make_shared<Task365>(tensor365, pindex);
  task364->add_dep(task365);
  task365->add_dep(task83);
  residualq->add_task(task365);

  auto tensor366 = vector<shared_ptr<Tensor>>{I618, Gamma204_(), v2_};
  auto task366 = make_shared<Task366>(tensor366, pindex);
  task364->add_dep(task366);
  task366->add_dep(task83);
  residualq->add_task(task366);

  vector<IndexRange> I624_index = {virt_, active_, active_, active_};
  auto I624 = make_shared<Tensor>(I624_index);
  auto tensor367 = vector<shared_ptr<Tensor>>{I63, v2_, I624};
  auto task367 = make_shared<Task367>(tensor367, pindex);
  task297->add_dep(task367);
  task367->add_dep(task83);
  residualq->add_task(task367);

  auto tensor368 = vector<shared_ptr<Tensor>>{I624, Gamma26_(), t2};
  auto task368 = make_shared<Task368>(tensor368, pindex);
  task367->add_dep(task368);
  task368->add_dep(task83);
  residualq->add_task(task368);

  vector<IndexRange> I627_index = {virt_, active_, active_, active_};
  auto I627 = make_shared<Tensor>(I627_index);
  auto tensor369 = vector<shared_ptr<Tensor>>{I63, v2_, I627};
  auto task369 = make_shared<Task369>(tensor369, pindex);
  task297->add_dep(task369);
  task369->add_dep(task83);
  residualq->add_task(task369);

  auto tensor370 = vector<shared_ptr<Tensor>>{I627, Gamma26_(), t2};
  auto task370 = make_shared<Task370>(tensor370, pindex);
  task369->add_dep(task370);
  task370->add_dep(task83);
  residualq->add_task(task370);

  vector<IndexRange> I666_index = {virt_, active_, active_, active_};
  auto I666 = make_shared<Tensor>(I666_index);
  auto tensor371 = vector<shared_ptr<Tensor>>{I63, t2, I666};
  auto task371 = make_shared<Task371>(tensor371, pindex);
  task297->add_dep(task371);
  task371->add_dep(task83);
  residualq->add_task(task371);

  auto tensor372 = vector<shared_ptr<Tensor>>{I666, Gamma193_(), v2_};
  auto task372 = make_shared<Task372>(tensor372, pindex);
  task371->add_dep(task372);
  task372->add_dep(task83);
  residualq->add_task(task372);

  auto tensor373 = vector<shared_ptr<Tensor>>{I666, Gamma26_(), v2_};
  auto task373 = make_shared<Task373>(tensor373, pindex);
  task371->add_dep(task373);
  task373->add_dep(task83);
  residualq->add_task(task373);

  vector<IndexRange> I669_index = {virt_, active_, active_, active_};
  auto I669 = make_shared<Tensor>(I669_index);
  auto tensor374 = vector<shared_ptr<Tensor>>{I63, t2, I669};
  auto task374 = make_shared<Task374>(tensor374, pindex);
  task297->add_dep(task374);
  task374->add_dep(task83);
  residualq->add_task(task374);

  auto tensor375 = vector<shared_ptr<Tensor>>{I669, Gamma193_(), v2_};
  auto task375 = make_shared<Task375>(tensor375, pindex);
  task374->add_dep(task375);
  task375->add_dep(task83);
  residualq->add_task(task375);

  auto tensor376 = vector<shared_ptr<Tensor>>{I669, Gamma26_(), v2_};
  auto task376 = make_shared<Task376>(tensor376, pindex);
  task374->add_dep(task376);
  task376->add_dep(task83);
  residualq->add_task(task376);

  vector<IndexRange> I696_index = {closed_, active_, active_, active_, virt_, active_};
  auto I696 = make_shared<Tensor>(I696_index);
  auto tensor377 = vector<shared_ptr<Tensor>>{I63, Gamma229_(), I696};
  auto task377 = make_shared<Task377>(tensor377, pindex);
  task297->add_dep(task377);
  task377->add_dep(task83);
  residualq->add_task(task377);

  auto tensor378 = vector<shared_ptr<Tensor>>{I696, t2, v2_};
  auto task378 = make_shared<Task378>(tensor378, pindex);
  task377->add_dep(task378);
  task378->add_dep(task83);
  residualq->add_task(task378);

  auto tensor379 = vector<shared_ptr<Tensor>>{I63, Gamma420_(), t2};
  auto task379 = make_shared<Task379>(tensor379, pindex);
  task297->add_dep(task379);
  task379->add_dep(task83);
  residualq->add_task(task379);

  auto tensor380 = vector<shared_ptr<Tensor>>{I63, Gamma421_(), t2};
  auto task380 = make_shared<Task380>(tensor380, pindex);
  task297->add_dep(task380);
  task380->add_dep(task83);
  residualq->add_task(task380);

  vector<IndexRange> I93_index = {virt_, active_, active_, active_};
  auto I93 = make_shared<Tensor>(I93_index);
  auto tensor381 = vector<shared_ptr<Tensor>>{r, I93};
  auto task381 = make_shared<Task381>(tensor381, pindex);
  task381->add_dep(task83);
  residualq->add_task(task381);

  vector<IndexRange> I94_index = {virt_, active_, active_, active_};
  auto I94 = make_shared<Tensor>(I94_index);
  auto tensor382 = vector<shared_ptr<Tensor>>{I93, Gamma31_(), I94};
  auto task382 = make_shared<Task382>(tensor382, pindex);
  task381->add_dep(task382);
  task382->add_dep(task83);
  residualq->add_task(task382);

  auto tensor383 = vector<shared_ptr<Tensor>>{I94, h1_, t2};
  auto task383 = make_shared<Task383>(tensor383, pindex);
  task382->add_dep(task383);
  task383->add_dep(task83);
  residualq->add_task(task383);

  auto tensor384 = vector<shared_ptr<Tensor>>{I94, t2, v2_};
  auto task384 = make_shared<Task384>(tensor384, pindex);
  task382->add_dep(task384);
  task384->add_dep(task83);
  residualq->add_task(task384);

  auto tensor385 = vector<shared_ptr<Tensor>>{I94, t2, v2_};
  auto task385 = make_shared<Task385>(tensor385, pindex);
  task382->add_dep(task385);
  task385->add_dep(task83);
  residualq->add_task(task385);

  auto tensor386 = vector<shared_ptr<Tensor>>{I94, t2, v2_};
  auto task386 = make_shared<Task386>(tensor386, pindex);
  task382->add_dep(task386);
  task386->add_dep(task83);
  residualq->add_task(task386);

  vector<IndexRange> I97_index = {virt_, active_, active_, active_};
  auto I97 = make_shared<Tensor>(I97_index);
  auto tensor387 = vector<shared_ptr<Tensor>>{I93, Gamma32_(), I97};
  auto task387 = make_shared<Task387>(tensor387, pindex);
  task381->add_dep(task387);
  task387->add_dep(task83);
  residualq->add_task(task387);

  auto tensor388 = vector<shared_ptr<Tensor>>{I97, t2, h1_};
  auto task388 = make_shared<Task388>(tensor388, pindex);
  task387->add_dep(task388);
  task388->add_dep(task83);
  residualq->add_task(task388);

  auto tensor389 = vector<shared_ptr<Tensor>>{I97, t2, h1_};
  auto task389 = make_shared<Task389>(tensor389, pindex);
  task387->add_dep(task389);
  task389->add_dep(task83);
  residualq->add_task(task389);

  vector<IndexRange> I736_index = {virt_, closed_, active_, active_};
  auto I736 = make_shared<Tensor>(I736_index);
  auto tensor390 = vector<shared_ptr<Tensor>>{I97, t2, I736};
  auto task390 = make_shared<Task390>(tensor390, pindex);
  task387->add_dep(task390);
  task390->add_dep(task83);
  residualq->add_task(task390);

  auto tensor391 = vector<shared_ptr<Tensor>>{I736, v2_};
  auto task391 = make_shared<Task391>(tensor391, pindex);
  task390->add_dep(task391);
  task391->add_dep(task83);
  residualq->add_task(task391);

  vector<IndexRange> I739_index = {virt_, closed_, active_, active_};
  auto I739 = make_shared<Tensor>(I739_index);
  auto tensor392 = vector<shared_ptr<Tensor>>{I97, t2, I739};
  auto task392 = make_shared<Task392>(tensor392, pindex);
  task387->add_dep(task392);
  task392->add_dep(task83);
  residualq->add_task(task392);

  auto tensor393 = vector<shared_ptr<Tensor>>{I739, v2_};
  auto task393 = make_shared<Task393>(tensor393, pindex);
  task392->add_dep(task393);
  task393->add_dep(task83);
  residualq->add_task(task393);

  auto tensor394 = vector<shared_ptr<Tensor>>{I97, t2, v2_};
  auto task394 = make_shared<Task394>(tensor394, pindex);
  task387->add_dep(task394);
  task394->add_dep(task83);
  residualq->add_task(task394);

  auto tensor395 = vector<shared_ptr<Tensor>>{I97, t2, v2_};
  auto task395 = make_shared<Task395>(tensor395, pindex);
  task387->add_dep(task395);
  task395->add_dep(task83);
  residualq->add_task(task395);

  vector<IndexRange> I100_index = {active_, virt_};
  auto I100 = make_shared<Tensor>(I100_index);
  auto tensor396 = vector<shared_ptr<Tensor>>{I93, Gamma33_(), I100};
  auto task396 = make_shared<Task396>(tensor396, pindex);
  task381->add_dep(task396);
  task396->add_dep(task83);
  residualq->add_task(task396);

  auto tensor397 = vector<shared_ptr<Tensor>>{I100, t2, h1_};
  auto task397 = make_shared<Task397>(tensor397, pindex);
  task396->add_dep(task397);
  task397->add_dep(task83);
  residualq->add_task(task397);

  auto tensor398 = vector<shared_ptr<Tensor>>{I100, h1_, t2};
  auto task398 = make_shared<Task398>(tensor398, pindex);
  task396->add_dep(task398);
  task398->add_dep(task83);
  residualq->add_task(task398);

  auto tensor399 = vector<shared_ptr<Tensor>>{I100, t2, v2_};
  auto task399 = make_shared<Task399>(tensor399, pindex);
  task396->add_dep(task399);
  task399->add_dep(task83);
  residualq->add_task(task399);

  auto tensor400 = vector<shared_ptr<Tensor>>{I100, t2, v2_};
  auto task400 = make_shared<Task400>(tensor400, pindex);
  task396->add_dep(task400);
  task400->add_dep(task83);
  residualq->add_task(task400);

  auto tensor401 = vector<shared_ptr<Tensor>>{I100, t2, v2_};
  auto task401 = make_shared<Task401>(tensor401, pindex);
  task396->add_dep(task401);
  task401->add_dep(task83);
  residualq->add_task(task401);

  auto tensor402 = vector<shared_ptr<Tensor>>{I100, t2, v2_};
  auto task402 = make_shared<Task402>(tensor402, pindex);
  task396->add_dep(task402);
  task402->add_dep(task83);
  residualq->add_task(task402);

  vector<IndexRange> I699_index = {active_, active_, active_, active_, active_, closed_};
  auto I699 = make_shared<Tensor>(I699_index);
  auto tensor403 = vector<shared_ptr<Tensor>>{I93, v2_, I699};
  auto task403 = make_shared<Task403>(tensor403, pindex);
  task381->add_dep(task403);
  task403->add_dep(task83);
  residualq->add_task(task403);

  auto tensor404 = vector<shared_ptr<Tensor>>{I699, t2, Gamma230_()};
  auto task404 = make_shared<Task404>(tensor404, pindex);
  task403->add_dep(task404);
  task404->add_dep(task83);
  residualq->add_task(task404);

  vector<IndexRange> I702_index = {active_, active_, virt_, active_};
  auto I702 = make_shared<Tensor>(I702_index);
  auto tensor405 = vector<shared_ptr<Tensor>>{I93, Gamma231_(), I702};
  auto task405 = make_shared<Task405>(tensor405, pindex);
  task381->add_dep(task405);
  task405->add_dep(task83);
  residualq->add_task(task405);

  auto tensor406 = vector<shared_ptr<Tensor>>{I702, t2, v2_};
  auto task406 = make_shared<Task406>(tensor406, pindex);
  task405->add_dep(task406);
  task406->add_dep(task83);
  residualq->add_task(task406);

  vector<IndexRange> I705_index = {closed_, active_, active_, active_, active_, active_};
  auto I705 = make_shared<Tensor>(I705_index);
  auto tensor407 = vector<shared_ptr<Tensor>>{I93, t2, I705};
  auto task407 = make_shared<Task407>(tensor407, pindex);
  task381->add_dep(task407);
  task407->add_dep(task83);
  residualq->add_task(task407);

  auto tensor408 = vector<shared_ptr<Tensor>>{I705, Gamma232_(), v2_};
  auto task408 = make_shared<Task408>(tensor408, pindex);
  task407->add_dep(task408);
  task408->add_dep(task83);
  residualq->add_task(task408);

  auto tensor409 = vector<shared_ptr<Tensor>>{I705, Gamma233_(), v2_};
  auto task409 = make_shared<Task409>(tensor409, pindex);
  task407->add_dep(task409);
  task409->add_dep(task83);
  residualq->add_task(task409);

  vector<IndexRange> I717_index = {virt_, active_, active_, active_, active_, active_};
  auto I717 = make_shared<Tensor>(I717_index);
  auto tensor410 = vector<shared_ptr<Tensor>>{I93, Gamma236_(), I717};
  auto task410 = make_shared<Task410>(tensor410, pindex);
  task381->add_dep(task410);
  task410->add_dep(task83);
  residualq->add_task(task410);

  vector<IndexRange> I718_index = {virt_, virt_, active_, active_};
  auto I718 = make_shared<Tensor>(I718_index);
  auto tensor411 = vector<shared_ptr<Tensor>>{I717, t2, I718};
  auto task411 = make_shared<Task411>(tensor411, pindex);
  task410->add_dep(task411);
  task411->add_dep(task83);
  residualq->add_task(task411);

  auto tensor412 = vector<shared_ptr<Tensor>>{I718, v2_};
  auto task412 = make_shared<Task412>(tensor412, pindex);
  task411->add_dep(task412);
  task412->add_dep(task83);
  residualq->add_task(task412);

  auto tensor413 = vector<shared_ptr<Tensor>>{I717, t2, v2_};
  auto task413 = make_shared<Task413>(tensor413, pindex);
  task410->add_dep(task413);
  task413->add_dep(task83);
  residualq->add_task(task413);

  vector<IndexRange> I720_index = {active_, active_, virt_, active_, active_, active_};
  auto I720 = make_shared<Tensor>(I720_index);
  auto tensor414 = vector<shared_ptr<Tensor>>{I93, Gamma237_(), I720};
  auto task414 = make_shared<Task414>(tensor414, pindex);
  task381->add_dep(task414);
  task414->add_dep(task83);
  residualq->add_task(task414);

  auto tensor415 = vector<shared_ptr<Tensor>>{I720, t2, v2_};
  auto task415 = make_shared<Task415>(tensor415, pindex);
  task414->add_dep(task415);
  task415->add_dep(task83);
  residualq->add_task(task415);

  vector<IndexRange> I723_index = {active_, virt_, active_, active_, active_, active_};
  auto I723 = make_shared<Tensor>(I723_index);
  auto tensor416 = vector<shared_ptr<Tensor>>{I93, Gamma238_(), I723};
  auto task416 = make_shared<Task416>(tensor416, pindex);
  task381->add_dep(task416);
  task416->add_dep(task83);
  residualq->add_task(task416);

  auto tensor417 = vector<shared_ptr<Tensor>>{I723, t2, v2_};
  auto task417 = make_shared<Task417>(tensor417, pindex);
  task416->add_dep(task417);
  task417->add_dep(task83);
  residualq->add_task(task417);

  vector<IndexRange> I744_index = {active_, active_, active_, virt_};
  auto I744 = make_shared<Tensor>(I744_index);
  auto tensor418 = vector<shared_ptr<Tensor>>{I93, Gamma245_(), I744};
  auto task418 = make_shared<Task418>(tensor418, pindex);
  task381->add_dep(task418);
  task418->add_dep(task83);
  residualq->add_task(task418);

  auto tensor419 = vector<shared_ptr<Tensor>>{I744, t2, v2_};
  auto task419 = make_shared<Task419>(tensor419, pindex);
  task418->add_dep(task419);
  task419->add_dep(task83);
  residualq->add_task(task419);

  vector<IndexRange> I747_index = {active_, active_, active_, virt_};
  auto I747 = make_shared<Tensor>(I747_index);
  auto tensor420 = vector<shared_ptr<Tensor>>{I93, Gamma246_(), I747};
  auto task420 = make_shared<Task420>(tensor420, pindex);
  task381->add_dep(task420);
  task420->add_dep(task83);
  residualq->add_task(task420);

  auto tensor421 = vector<shared_ptr<Tensor>>{I747, t2, v2_};
  auto task421 = make_shared<Task421>(tensor421, pindex);
  task420->add_dep(task421);
  task421->add_dep(task83);
  residualq->add_task(task421);

  vector<IndexRange> I768_index = {active_, active_, active_, active_, virt_, active_};
  auto I768 = make_shared<Tensor>(I768_index);
  auto tensor422 = vector<shared_ptr<Tensor>>{I93, Gamma253_(), I768};
  auto task422 = make_shared<Task422>(tensor422, pindex);
  task381->add_dep(task422);
  task422->add_dep(task83);
  residualq->add_task(task422);

  auto tensor423 = vector<shared_ptr<Tensor>>{I768, t2, v2_};
  auto task423 = make_shared<Task423>(tensor423, pindex);
  task422->add_dep(task423);
  task423->add_dep(task83);
  residualq->add_task(task423);

  auto tensor424 = vector<shared_ptr<Tensor>>{I93, Gamma422_(), t2};
  auto task424 = make_shared<Task424>(tensor424, pindex);
  task381->add_dep(task424);
  task424->add_dep(task83);
  residualq->add_task(task424);

  auto tensor425 = vector<shared_ptr<Tensor>>{I93, Gamma423_(), t2};
  auto task425 = make_shared<Task425>(tensor425, pindex);
  task381->add_dep(task425);
  task425->add_dep(task83);
  residualq->add_task(task425);
}

#endif
