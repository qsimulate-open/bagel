//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_densityqq.cc
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


#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_densityq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto densityq = make_shared<Queue>();
  auto tensor283 = vector<shared_ptr<Tensor>>{den2};
  auto task283 = make_shared<Task283>(tensor283, reset);
  densityq->add_task(task283);

  vector<IndexRange> I374_index = {active_, active_};
  auto I374 = make_shared<Tensor>(I374_index);
  auto tensor284 = vector<shared_ptr<Tensor>>{den2, I374};
  auto task284 = make_shared<Task284>(tensor284, pindex);
  task284->add_dep(task283);
  densityq->add_task(task284);

  vector<IndexRange> I375_index = {active_, active_, active_, active_};
  auto I375 = make_shared<Tensor>(I375_index);
  auto tensor285 = vector<shared_ptr<Tensor>>{I374, Gamma138_(), I375};
  auto task285 = make_shared<Task285>(tensor285, pindex);
  task284->add_dep(task285);
  task285->add_dep(task283);
  densityq->add_task(task285);

  auto tensor286 = vector<shared_ptr<Tensor>>{I375, t2};
  auto task286 = make_shared<Task286>(tensor286, pindex);
  task285->add_dep(task286);
  task286->add_dep(task283);
  densityq->add_task(task286);

  vector<IndexRange> I468_index = {active_, active_, active_, active_};
  auto I468 = make_shared<Tensor>(I468_index);
  auto tensor287 = vector<shared_ptr<Tensor>>{I374, Gamma169_(), I468};
  auto task287 = make_shared<Task287>(tensor287, pindex);
  task284->add_dep(task287);
  task287->add_dep(task283);
  densityq->add_task(task287);

  auto tensor288 = vector<shared_ptr<Tensor>>{I468, t2};
  auto task288 = make_shared<Task288>(tensor288, pindex);
  task287->add_dep(task288);
  task288->add_dep(task283);
  densityq->add_task(task288);

  vector<IndexRange> I477_index = {active_, active_, active_, active_};
  auto I477 = make_shared<Tensor>(I477_index);
  auto tensor289 = vector<shared_ptr<Tensor>>{I374, Gamma172_(), I477};
  auto task289 = make_shared<Task289>(tensor289, pindex);
  task284->add_dep(task289);
  task289->add_dep(task283);
  densityq->add_task(task289);

  auto tensor290 = vector<shared_ptr<Tensor>>{I477, t2};
  auto task290 = make_shared<Task290>(tensor290, pindex);
  task289->add_dep(task290);
  task290->add_dep(task283);
  densityq->add_task(task290);

  auto tensor291 = vector<shared_ptr<Tensor>>{I477, t2};
  auto task291 = make_shared<Task291>(tensor291, pindex);
  task289->add_dep(task291);
  task291->add_dep(task283);
  densityq->add_task(task291);

  auto tensor292 = vector<shared_ptr<Tensor>>{I477, t2};
  auto task292 = make_shared<Task292>(tensor292, pindex);
  task289->add_dep(task292);
  task292->add_dep(task283);
  densityq->add_task(task292);

  vector<IndexRange> I659_index = {active_, active_, active_, active_};
  auto I659 = make_shared<Tensor>(I659_index);
  auto tensor293 = vector<shared_ptr<Tensor>>{I374, Gamma230_(), I659};
  auto task293 = make_shared<Task293>(tensor293, pindex);
  task284->add_dep(task293);
  task293->add_dep(task283);
  densityq->add_task(task293);

  auto tensor294 = vector<shared_ptr<Tensor>>{I659, t2};
  auto task294 = make_shared<Task294>(tensor294, pindex);
  task293->add_dep(task294);
  task294->add_dep(task283);
  densityq->add_task(task294);

  vector<IndexRange> I377_index = {closed_, closed_};
  auto I377 = make_shared<Tensor>(I377_index);
  auto tensor295 = vector<shared_ptr<Tensor>>{den2, I377};
  auto task295 = make_shared<Task295>(tensor295, pindex);
  task295->add_dep(task283);
  densityq->add_task(task295);

  vector<IndexRange> I378_index = {closed_, closed_, active_, active_};
  auto I378 = make_shared<Tensor>(I378_index);
  auto tensor296 = vector<shared_ptr<Tensor>>{I377, t2, I378};
  auto task296 = make_shared<Task296>(tensor296, pindex);
  task295->add_dep(task296);
  task296->add_dep(task283);
  densityq->add_task(task296);

  auto tensor297 = vector<shared_ptr<Tensor>>{I378, Gamma92_(), t2};
  auto task297 = make_shared<Task297>(tensor297, pindex);
  task296->add_dep(task297);
  task297->add_dep(task283);
  densityq->add_task(task297);

  vector<IndexRange> I471_index = {closed_, virt_, active_, active_};
  auto I471 = make_shared<Tensor>(I471_index);
  auto tensor298 = vector<shared_ptr<Tensor>>{I377, t2, I471};
  auto task298 = make_shared<Task298>(tensor298, pindex);
  task295->add_dep(task298);
  task298->add_dep(task283);
  densityq->add_task(task298);

  auto tensor299 = vector<shared_ptr<Tensor>>{I471, Gamma32_(), t2};
  auto task299 = make_shared<Task299>(tensor299, pindex);
  task298->add_dep(task299);
  task299->add_dep(task283);
  densityq->add_task(task299);

  vector<IndexRange> I480_index = {closed_, virt_, active_, active_};
  auto I480 = make_shared<Tensor>(I480_index);
  auto tensor300 = vector<shared_ptr<Tensor>>{I377, t2, I480};
  auto task300 = make_shared<Task300>(tensor300, pindex);
  task295->add_dep(task300);
  task300->add_dep(task283);
  densityq->add_task(task300);

  auto tensor301 = vector<shared_ptr<Tensor>>{I480, Gamma35_(), t2};
  auto task301 = make_shared<Task301>(tensor301, pindex);
  task300->add_dep(task301);
  task301->add_dep(task283);
  densityq->add_task(task301);

  vector<IndexRange> I380_index = {active_, closed_};
  auto I380 = make_shared<Tensor>(I380_index);
  auto tensor302 = vector<shared_ptr<Tensor>>{den2, I380};
  auto task302 = make_shared<Task302>(tensor302, pindex);
  task302->add_dep(task283);
  densityq->add_task(task302);

  vector<IndexRange> I381_index = {closed_, active_, active_, active_};
  auto I381 = make_shared<Tensor>(I381_index);
  auto tensor303 = vector<shared_ptr<Tensor>>{I380, t2, I381};
  auto task303 = make_shared<Task303>(tensor303, pindex);
  task302->add_dep(task303);
  task303->add_dep(task283);
  densityq->add_task(task303);

  auto tensor304 = vector<shared_ptr<Tensor>>{I381, Gamma2_(), t2};
  auto task304 = make_shared<Task304>(tensor304, pindex);
  task303->add_dep(task304);
  task304->add_dep(task283);
  densityq->add_task(task304);

  vector<IndexRange> I486_index = {virt_, active_, active_, active_};
  auto I486 = make_shared<Tensor>(I486_index);
  auto tensor305 = vector<shared_ptr<Tensor>>{I380, t2, I486};
  auto task305 = make_shared<Task305>(tensor305, pindex);
  task302->add_dep(task305);
  task305->add_dep(task283);
  densityq->add_task(task305);

  auto tensor306 = vector<shared_ptr<Tensor>>{I486, Gamma37_(), t2};
  auto task306 = make_shared<Task306>(tensor306, pindex);
  task305->add_dep(task306);
  task306->add_dep(task283);
  densityq->add_task(task306);

  vector<IndexRange> I383_index = {active_, virt_};
  auto I383 = make_shared<Tensor>(I383_index);
  auto tensor307 = vector<shared_ptr<Tensor>>{den2, I383};
  auto task307 = make_shared<Task307>(tensor307, pindex);
  task307->add_dep(task283);
  densityq->add_task(task307);

  vector<IndexRange> I384_index = {closed_, closed_, active_, active_};
  auto I384 = make_shared<Tensor>(I384_index);
  auto tensor308 = vector<shared_ptr<Tensor>>{I383, t2, I384};
  auto task308 = make_shared<Task308>(tensor308, pindex);
  task307->add_dep(task308);
  task308->add_dep(task283);
  densityq->add_task(task308);

  auto tensor309 = vector<shared_ptr<Tensor>>{I384, Gamma3_(), t2};
  auto task309 = make_shared<Task309>(tensor309, pindex);
  task308->add_dep(task309);
  task309->add_dep(task283);
  densityq->add_task(task309);

  vector<IndexRange> I495_index = {closed_, virt_, active_, active_};
  auto I495 = make_shared<Tensor>(I495_index);
  auto tensor310 = vector<shared_ptr<Tensor>>{I383, t2, I495};
  auto task310 = make_shared<Task310>(tensor310, pindex);
  task307->add_dep(task310);
  task310->add_dep(task283);
  densityq->add_task(task310);

  auto tensor311 = vector<shared_ptr<Tensor>>{I495, Gamma35_(), t2};
  auto task311 = make_shared<Task311>(tensor311, pindex);
  task310->add_dep(task311);
  task311->add_dep(task283);
  densityq->add_task(task311);

  vector<IndexRange> I498_index = {closed_, virt_, active_, active_};
  auto I498 = make_shared<Tensor>(I498_index);
  auto tensor312 = vector<shared_ptr<Tensor>>{I383, t2, I498};
  auto task312 = make_shared<Task312>(tensor312, pindex);
  task307->add_dep(task312);
  task312->add_dep(task283);
  densityq->add_task(task312);

  auto tensor313 = vector<shared_ptr<Tensor>>{I498, Gamma32_(), t2};
  auto task313 = make_shared<Task313>(tensor313, pindex);
  task312->add_dep(task313);
  task313->add_dep(task283);
  densityq->add_task(task313);

  vector<IndexRange> I537_index = {virt_, closed_, active_, active_};
  auto I537 = make_shared<Tensor>(I537_index);
  auto tensor314 = vector<shared_ptr<Tensor>>{I383, t2, I537};
  auto task314 = make_shared<Task314>(tensor314, pindex);
  task307->add_dep(task314);
  task314->add_dep(task283);
  densityq->add_task(task314);

  auto tensor315 = vector<shared_ptr<Tensor>>{I537, Gamma35_(), t2};
  auto task315 = make_shared<Task315>(tensor315, pindex);
  task314->add_dep(task315);
  task315->add_dep(task283);
  densityq->add_task(task315);

  vector<IndexRange> I540_index = {virt_, closed_, active_, active_};
  auto I540 = make_shared<Tensor>(I540_index);
  auto tensor316 = vector<shared_ptr<Tensor>>{I383, t2, I540};
  auto task316 = make_shared<Task316>(tensor316, pindex);
  task307->add_dep(task316);
  task316->add_dep(task283);
  densityq->add_task(task316);

  auto tensor317 = vector<shared_ptr<Tensor>>{I540, Gamma35_(), t2};
  auto task317 = make_shared<Task317>(tensor317, pindex);
  task316->add_dep(task317);
  task317->add_dep(task283);
  densityq->add_task(task317);

  vector<IndexRange> I386_index = {active_, closed_};
  auto I386 = make_shared<Tensor>(I386_index);
  auto tensor318 = vector<shared_ptr<Tensor>>{den2, I386};
  auto task318 = make_shared<Task318>(tensor318, pindex);
  task318->add_dep(task283);
  densityq->add_task(task318);

  vector<IndexRange> I387_index = {closed_, active_, active_, active_};
  auto I387 = make_shared<Tensor>(I387_index);
  auto tensor319 = vector<shared_ptr<Tensor>>{I386, t2, I387};
  auto task319 = make_shared<Task319>(tensor319, pindex);
  task318->add_dep(task319);
  task319->add_dep(task283);
  densityq->add_task(task319);

  auto tensor320 = vector<shared_ptr<Tensor>>{I387, Gamma4_(), t2};
  auto task320 = make_shared<Task320>(tensor320, pindex);
  task319->add_dep(task320);
  task320->add_dep(task283);
  densityq->add_task(task320);

  vector<IndexRange> I543_index = {virt_, active_, active_, active_};
  auto I543 = make_shared<Tensor>(I543_index);
  auto tensor321 = vector<shared_ptr<Tensor>>{I386, t2, I543};
  auto task321 = make_shared<Task321>(tensor321, pindex);
  task318->add_dep(task321);
  task321->add_dep(task283);
  densityq->add_task(task321);

  auto tensor322 = vector<shared_ptr<Tensor>>{I543, Gamma56_(), t2};
  auto task322 = make_shared<Task322>(tensor322, pindex);
  task321->add_dep(task322);
  task322->add_dep(task283);
  densityq->add_task(task322);

  vector<IndexRange> I546_index = {virt_, active_, active_, active_};
  auto I546 = make_shared<Tensor>(I546_index);
  auto tensor323 = vector<shared_ptr<Tensor>>{I386, t2, I546};
  auto task323 = make_shared<Task323>(tensor323, pindex);
  task318->add_dep(task323);
  task323->add_dep(task283);
  densityq->add_task(task323);

  auto tensor324 = vector<shared_ptr<Tensor>>{I546, Gamma57_(), t2};
  auto task324 = make_shared<Task324>(tensor324, pindex);
  task323->add_dep(task324);
  task324->add_dep(task283);
  densityq->add_task(task324);

  vector<IndexRange> I389_index = {active_, active_};
  auto I389 = make_shared<Tensor>(I389_index);
  auto tensor325 = vector<shared_ptr<Tensor>>{den2, I389};
  auto task325 = make_shared<Task325>(tensor325, pindex);
  task325->add_dep(task283);
  densityq->add_task(task325);

  vector<IndexRange> I390_index = {active_, active_, active_, active_, active_, active_};
  auto I390 = make_shared<Tensor>(I390_index);
  auto tensor326 = vector<shared_ptr<Tensor>>{I389, Gamma143_(), I390};
  auto task326 = make_shared<Task326>(tensor326, pindex);
  task325->add_dep(task326);
  task326->add_dep(task283);
  densityq->add_task(task326);

  auto tensor327 = vector<shared_ptr<Tensor>>{I390, t2};
  auto task327 = make_shared<Task327>(tensor327, pindex);
  task326->add_dep(task327);
  task327->add_dep(task283);
  densityq->add_task(task327);

  vector<IndexRange> I549_index = {active_, active_, active_, active_, active_, active_};
  auto I549 = make_shared<Tensor>(I549_index);
  auto tensor328 = vector<shared_ptr<Tensor>>{I389, Gamma196_(), I549};
  auto task328 = make_shared<Task328>(tensor328, pindex);
  task325->add_dep(task328);
  task328->add_dep(task283);
  densityq->add_task(task328);

  auto tensor329 = vector<shared_ptr<Tensor>>{I549, t2};
  auto task329 = make_shared<Task329>(tensor329, pindex);
  task328->add_dep(task329);
  task329->add_dep(task283);
  densityq->add_task(task329);

  vector<IndexRange> I392_index = {closed_, closed_};
  auto I392 = make_shared<Tensor>(I392_index);
  auto tensor330 = vector<shared_ptr<Tensor>>{den2, I392};
  auto task330 = make_shared<Task330>(tensor330, pindex);
  task330->add_dep(task283);
  densityq->add_task(task330);

  vector<IndexRange> I393_index = {closed_, active_, active_, active_};
  auto I393 = make_shared<Tensor>(I393_index);
  auto tensor331 = vector<shared_ptr<Tensor>>{I392, t2, I393};
  auto task331 = make_shared<Task331>(tensor331, pindex);
  task330->add_dep(task331);
  task331->add_dep(task283);
  densityq->add_task(task331);

  auto tensor332 = vector<shared_ptr<Tensor>>{I393, Gamma6_(), t2};
  auto task332 = make_shared<Task332>(tensor332, pindex);
  task331->add_dep(task332);
  task332->add_dep(task283);
  densityq->add_task(task332);

  vector<IndexRange> I395_index = {closed_, virt_};
  auto I395 = make_shared<Tensor>(I395_index);
  auto tensor333 = vector<shared_ptr<Tensor>>{den2, I395};
  auto task333 = make_shared<Task333>(tensor333, pindex);
  task333->add_dep(task283);
  densityq->add_task(task333);

  vector<IndexRange> I396_index = {closed_, active_};
  auto I396 = make_shared<Tensor>(I396_index);
  auto tensor334 = vector<shared_ptr<Tensor>>{I395, t2, I396};
  auto task334 = make_shared<Task334>(tensor334, pindex);
  task333->add_dep(task334);
  task334->add_dep(task283);
  densityq->add_task(task334);

  auto tensor335 = vector<shared_ptr<Tensor>>{I396, Gamma7_(), t2};
  auto task335 = make_shared<Task335>(tensor335, pindex);
  task334->add_dep(task335);
  task335->add_dep(task283);
  densityq->add_task(task335);

  vector<IndexRange> I399_index = {closed_, active_};
  auto I399 = make_shared<Tensor>(I399_index);
  auto tensor336 = vector<shared_ptr<Tensor>>{I395, t2, I399};
  auto task336 = make_shared<Task336>(tensor336, pindex);
  task333->add_dep(task336);
  task336->add_dep(task283);
  densityq->add_task(task336);

  auto tensor337 = vector<shared_ptr<Tensor>>{I399, Gamma7_(), t2};
  auto task337 = make_shared<Task337>(tensor337, pindex);
  task336->add_dep(task337);
  task337->add_dep(task283);
  densityq->add_task(task337);

  vector<IndexRange> I555_index = {virt_, active_};
  auto I555 = make_shared<Tensor>(I555_index);
  auto tensor338 = vector<shared_ptr<Tensor>>{I395, t2, I555};
  auto task338 = make_shared<Task338>(tensor338, pindex);
  task333->add_dep(task338);
  task338->add_dep(task283);
  densityq->add_task(task338);

  auto tensor339 = vector<shared_ptr<Tensor>>{I555, Gamma60_(), t2};
  auto task339 = make_shared<Task339>(tensor339, pindex);
  task338->add_dep(task339);
  task339->add_dep(task283);
  densityq->add_task(task339);

  vector<IndexRange> I558_index = {virt_, active_};
  auto I558 = make_shared<Tensor>(I558_index);
  auto tensor340 = vector<shared_ptr<Tensor>>{I395, t2, I558};
  auto task340 = make_shared<Task340>(tensor340, pindex);
  task333->add_dep(task340);
  task340->add_dep(task283);
  densityq->add_task(task340);

  auto tensor341 = vector<shared_ptr<Tensor>>{I558, Gamma60_(), t2};
  auto task341 = make_shared<Task341>(tensor341, pindex);
  task340->add_dep(task341);
  task341->add_dep(task283);
  densityq->add_task(task341);

  vector<IndexRange> I401_index = {active_, virt_};
  auto I401 = make_shared<Tensor>(I401_index);
  auto tensor342 = vector<shared_ptr<Tensor>>{den2, I401};
  auto task342 = make_shared<Task342>(tensor342, pindex);
  task342->add_dep(task283);
  densityq->add_task(task342);

  vector<IndexRange> I402_index = {closed_, active_, active_, active_};
  auto I402 = make_shared<Tensor>(I402_index);
  auto tensor343 = vector<shared_ptr<Tensor>>{I401, t2, I402};
  auto task343 = make_shared<Task343>(tensor343, pindex);
  task342->add_dep(task343);
  task343->add_dep(task283);
  densityq->add_task(task343);

  auto tensor344 = vector<shared_ptr<Tensor>>{I402, Gamma9_(), t2};
  auto task344 = make_shared<Task344>(tensor344, pindex);
  task343->add_dep(task344);
  task344->add_dep(task283);
  densityq->add_task(task344);

  vector<IndexRange> I405_index = {closed_, active_, active_, active_};
  auto I405 = make_shared<Tensor>(I405_index);
  auto tensor345 = vector<shared_ptr<Tensor>>{I401, t2, I405};
  auto task345 = make_shared<Task345>(tensor345, pindex);
  task342->add_dep(task345);
  task345->add_dep(task283);
  densityq->add_task(task345);

  auto tensor346 = vector<shared_ptr<Tensor>>{I405, Gamma6_(), t2};
  auto task346 = make_shared<Task346>(tensor346, pindex);
  task345->add_dep(task346);
  task346->add_dep(task283);
  densityq->add_task(task346);

  vector<IndexRange> I561_index = {virt_, active_, active_, active_};
  auto I561 = make_shared<Tensor>(I561_index);
  auto tensor347 = vector<shared_ptr<Tensor>>{I401, t2, I561};
  auto task347 = make_shared<Task347>(tensor347, pindex);
  task342->add_dep(task347);
  task347->add_dep(task283);
  densityq->add_task(task347);

  auto tensor348 = vector<shared_ptr<Tensor>>{I561, Gamma59_(), t2};
  auto task348 = make_shared<Task348>(tensor348, pindex);
  task347->add_dep(task348);
  task348->add_dep(task283);
  densityq->add_task(task348);

  vector<IndexRange> I407_index = {active_, virt_};
  auto I407 = make_shared<Tensor>(I407_index);
  auto tensor349 = vector<shared_ptr<Tensor>>{den2, I407};
  auto task349 = make_shared<Task349>(tensor349, pindex);
  task349->add_dep(task283);
  densityq->add_task(task349);

  vector<IndexRange> I408_index = {closed_, closed_, active_, active_};
  auto I408 = make_shared<Tensor>(I408_index);
  auto tensor350 = vector<shared_ptr<Tensor>>{I407, t2, I408};
  auto task350 = make_shared<Task350>(tensor350, pindex);
  task349->add_dep(task350);
  task350->add_dep(task283);
  densityq->add_task(task350);

  auto tensor351 = vector<shared_ptr<Tensor>>{I408, Gamma3_(), t2};
  auto task351 = make_shared<Task351>(tensor351, pindex);
  task350->add_dep(task351);
  task351->add_dep(task283);
  densityq->add_task(task351);

  vector<IndexRange> I410_index = {virt_, closed_};
  auto I410 = make_shared<Tensor>(I410_index);
  auto tensor352 = vector<shared_ptr<Tensor>>{den2, I410};
  auto task352 = make_shared<Task352>(tensor352, pindex);
  task352->add_dep(task283);
  densityq->add_task(task352);

  vector<IndexRange> I411_index = {closed_, active_};
  auto I411 = make_shared<Tensor>(I411_index);
  auto tensor353 = vector<shared_ptr<Tensor>>{I410, t2, I411};
  auto task353 = make_shared<Task353>(tensor353, pindex);
  task352->add_dep(task353);
  task353->add_dep(task283);
  densityq->add_task(task353);

  auto tensor354 = vector<shared_ptr<Tensor>>{I411, Gamma12_(), t2};
  auto task354 = make_shared<Task354>(tensor354, pindex);
  task353->add_dep(task354);
  task354->add_dep(task283);
  densityq->add_task(task354);

  vector<IndexRange> I413_index = {closed_, virt_};
  auto I413 = make_shared<Tensor>(I413_index);
  auto tensor355 = vector<shared_ptr<Tensor>>{den2, I413};
  auto task355 = make_shared<Task355>(tensor355, pindex);
  task355->add_dep(task283);
  densityq->add_task(task355);

  vector<IndexRange> I414_index = {closed_, active_};
  auto I414 = make_shared<Tensor>(I414_index);
  auto tensor356 = vector<shared_ptr<Tensor>>{I413, t2, I414};
  auto task356 = make_shared<Task356>(tensor356, pindex);
  task355->add_dep(task356);
  task356->add_dep(task283);
  densityq->add_task(task356);

  auto tensor357 = vector<shared_ptr<Tensor>>{I414, Gamma12_(), t2};
  auto task357 = make_shared<Task357>(tensor357, pindex);
  task356->add_dep(task357);
  task357->add_dep(task283);
  densityq->add_task(task357);

  vector<IndexRange> I570_index = {virt_, closed_};
  auto I570 = make_shared<Tensor>(I570_index);
  auto tensor358 = vector<shared_ptr<Tensor>>{I413, t2, I570};
  auto task358 = make_shared<Task358>(tensor358, pindex);
  task355->add_dep(task358);
  task358->add_dep(task283);
  densityq->add_task(task358);

  vector<IndexRange> I571_index = {active_, virt_, closed_, active_};
  auto I571 = make_shared<Tensor>(I571_index);
  auto tensor359 = vector<shared_ptr<Tensor>>{I570, Gamma38_(), I571};
  auto task359 = make_shared<Task359>(tensor359, pindex);
  task358->add_dep(task359);
  task359->add_dep(task283);
  densityq->add_task(task359);

  auto tensor360 = vector<shared_ptr<Tensor>>{I571, t2};
  auto task360 = make_shared<Task360>(tensor360, pindex);
  task359->add_dep(task360);
  task360->add_dep(task283);
  densityq->add_task(task360);

  vector<IndexRange> I416_index = {active_, active_};
  auto I416 = make_shared<Tensor>(I416_index);
  auto tensor361 = vector<shared_ptr<Tensor>>{den2, I416};
  auto task361 = make_shared<Task361>(tensor361, pindex);
  task361->add_dep(task283);
  densityq->add_task(task361);

  vector<IndexRange> I417_index = {active_, active_};
  auto I417 = make_shared<Tensor>(I417_index);
  auto tensor362 = vector<shared_ptr<Tensor>>{I416, Gamma152_(), I417};
  auto task362 = make_shared<Task362>(tensor362, pindex);
  task361->add_dep(task362);
  task362->add_dep(task283);
  densityq->add_task(task362);

  auto tensor363 = vector<shared_ptr<Tensor>>{I417, t2};
  auto task363 = make_shared<Task363>(tensor363, pindex);
  task362->add_dep(task363);
  task363->add_dep(task283);
  densityq->add_task(task363);

  auto tensor364 = vector<shared_ptr<Tensor>>{I417, t2};
  auto task364 = make_shared<Task364>(tensor364, pindex);
  task362->add_dep(task364);
  task364->add_dep(task283);
  densityq->add_task(task364);

  vector<IndexRange> I626_index = {active_, active_};
  auto I626 = make_shared<Tensor>(I626_index);
  auto tensor365 = vector<shared_ptr<Tensor>>{I416, Gamma60_(), I626};
  auto task365 = make_shared<Task365>(tensor365, pindex);
  task361->add_dep(task365);
  task365->add_dep(task283);
  densityq->add_task(task365);

  auto tensor366 = vector<shared_ptr<Tensor>>{I626, t2};
  auto task366 = make_shared<Task366>(tensor366, pindex);
  task365->add_dep(task366);
  task366->add_dep(task283);
  densityq->add_task(task366);

  auto tensor367 = vector<shared_ptr<Tensor>>{I626, t2};
  auto task367 = make_shared<Task367>(tensor367, pindex);
  task365->add_dep(task367);
  task367->add_dep(task283);
  densityq->add_task(task367);

  vector<IndexRange> I422_index = {closed_, closed_};
  auto I422 = make_shared<Tensor>(I422_index);
  auto tensor368 = vector<shared_ptr<Tensor>>{den2, I422};
  auto task368 = make_shared<Task368>(tensor368, pindex);
  task368->add_dep(task283);
  densityq->add_task(task368);

  vector<IndexRange> I423_index = {closed_, virt_, closed_, active_};
  auto I423 = make_shared<Tensor>(I423_index);
  auto tensor369 = vector<shared_ptr<Tensor>>{I422, t2, I423};
  auto task369 = make_shared<Task369>(tensor369, pindex);
  task368->add_dep(task369);
  task369->add_dep(task283);
  densityq->add_task(task369);

  auto tensor370 = vector<shared_ptr<Tensor>>{I423, Gamma16_(), t2};
  auto task370 = make_shared<Task370>(tensor370, pindex);
  task369->add_dep(task370);
  task370->add_dep(task283);
  densityq->add_task(task370);

  vector<IndexRange> I426_index = {closed_, virt_, closed_, active_};
  auto I426 = make_shared<Tensor>(I426_index);
  auto tensor371 = vector<shared_ptr<Tensor>>{I422, t2, I426};
  auto task371 = make_shared<Task371>(tensor371, pindex);
  task368->add_dep(task371);
  task371->add_dep(task283);
  densityq->add_task(task371);

  auto tensor372 = vector<shared_ptr<Tensor>>{I426, Gamma16_(), t2};
  auto task372 = make_shared<Task372>(tensor372, pindex);
  task371->add_dep(task372);
  task372->add_dep(task283);
  densityq->add_task(task372);

  vector<IndexRange> I428_index = {closed_, closed_};
  auto I428 = make_shared<Tensor>(I428_index);
  auto tensor373 = vector<shared_ptr<Tensor>>{den2, I428};
  auto task373 = make_shared<Task373>(tensor373, pindex);
  task373->add_dep(task283);
  densityq->add_task(task373);

  vector<IndexRange> I429_index = {closed_, virt_, closed_, active_};
  auto I429 = make_shared<Tensor>(I429_index);
  auto tensor374 = vector<shared_ptr<Tensor>>{I428, t2, I429};
  auto task374 = make_shared<Task374>(tensor374, pindex);
  task373->add_dep(task374);
  task374->add_dep(task283);
  densityq->add_task(task374);

  auto tensor375 = vector<shared_ptr<Tensor>>{I429, Gamma16_(), t2};
  auto task375 = make_shared<Task375>(tensor375, pindex);
  task374->add_dep(task375);
  task375->add_dep(task283);
  densityq->add_task(task375);

  vector<IndexRange> I435_index = {closed_, virt_, closed_, active_};
  auto I435 = make_shared<Tensor>(I435_index);
  auto tensor376 = vector<shared_ptr<Tensor>>{I428, t2, I435};
  auto task376 = make_shared<Task376>(tensor376, pindex);
  task373->add_dep(task376);
  task376->add_dep(task283);
  densityq->add_task(task376);

  auto tensor377 = vector<shared_ptr<Tensor>>{I435, Gamma16_(), t2};
  auto task377 = make_shared<Task377>(tensor377, pindex);
  task376->add_dep(task377);
  task377->add_dep(task283);
  densityq->add_task(task377);

  vector<IndexRange> I431_index = {virt_, virt_};
  auto I431 = make_shared<Tensor>(I431_index);
  auto tensor378 = vector<shared_ptr<Tensor>>{den2, I431};
  auto task378 = make_shared<Task378>(tensor378, pindex);
  task378->add_dep(task283);
  densityq->add_task(task378);

  vector<IndexRange> I432_index = {closed_, virt_, closed_, active_};
  auto I432 = make_shared<Tensor>(I432_index);
  auto tensor379 = vector<shared_ptr<Tensor>>{I431, t2, I432};
  auto task379 = make_shared<Task379>(tensor379, pindex);
  task378->add_dep(task379);
  task379->add_dep(task283);
  densityq->add_task(task379);

  auto tensor380 = vector<shared_ptr<Tensor>>{I432, Gamma16_(), t2};
  auto task380 = make_shared<Task380>(tensor380, pindex);
  task379->add_dep(task380);
  task380->add_dep(task283);
  densityq->add_task(task380);

  vector<IndexRange> I438_index = {closed_, virt_, closed_, active_};
  auto I438 = make_shared<Tensor>(I438_index);
  auto tensor381 = vector<shared_ptr<Tensor>>{I431, t2, I438};
  auto task381 = make_shared<Task381>(tensor381, pindex);
  task378->add_dep(task381);
  task381->add_dep(task283);
  densityq->add_task(task381);

  auto tensor382 = vector<shared_ptr<Tensor>>{I438, Gamma16_(), t2};
  auto task382 = make_shared<Task382>(tensor382, pindex);
  task381->add_dep(task382);
  task382->add_dep(task283);
  densityq->add_task(task382);

  vector<IndexRange> I440_index = {active_, closed_};
  auto I440 = make_shared<Tensor>(I440_index);
  auto tensor383 = vector<shared_ptr<Tensor>>{den2, I440};
  auto task383 = make_shared<Task383>(tensor383, pindex);
  task383->add_dep(task283);
  densityq->add_task(task383);

  vector<IndexRange> I441_index = {virt_, closed_, active_, active_};
  auto I441 = make_shared<Tensor>(I441_index);
  auto tensor384 = vector<shared_ptr<Tensor>>{I440, t2, I441};
  auto task384 = make_shared<Task384>(tensor384, pindex);
  task383->add_dep(task384);
  task384->add_dep(task283);
  densityq->add_task(task384);

  auto tensor385 = vector<shared_ptr<Tensor>>{I441, Gamma22_(), t2};
  auto task385 = make_shared<Task385>(tensor385, pindex);
  task384->add_dep(task385);
  task385->add_dep(task283);
  densityq->add_task(task385);

  auto tensor386 = vector<shared_ptr<Tensor>>{I441, Gamma12_(), t2};
  auto task386 = make_shared<Task386>(tensor386, pindex);
  task384->add_dep(task386);
  task386->add_dep(task283);
  densityq->add_task(task386);

  vector<IndexRange> I443_index = {active_, closed_};
  auto I443 = make_shared<Tensor>(I443_index);
  auto tensor387 = vector<shared_ptr<Tensor>>{den2, I443};
  auto task387 = make_shared<Task387>(tensor387, pindex);
  task387->add_dep(task283);
  densityq->add_task(task387);

  vector<IndexRange> I444_index = {virt_, closed_, active_, active_};
  auto I444 = make_shared<Tensor>(I444_index);
  auto tensor388 = vector<shared_ptr<Tensor>>{I443, t2, I444};
  auto task388 = make_shared<Task388>(tensor388, pindex);
  task387->add_dep(task388);
  task388->add_dep(task283);
  densityq->add_task(task388);

  vector<IndexRange> I445_index = {active_, virt_, closed_, active_};
  auto I445 = make_shared<Tensor>(I445_index);
  auto tensor389 = vector<shared_ptr<Tensor>>{I444, Gamma12_(), I445};
  auto task389 = make_shared<Task389>(tensor389, pindex);
  task388->add_dep(task389);
  task389->add_dep(task283);
  densityq->add_task(task389);

  auto tensor390 = vector<shared_ptr<Tensor>>{I445, t2};
  auto task390 = make_shared<Task390>(tensor390, pindex);
  task389->add_dep(task390);
  task390->add_dep(task283);
  densityq->add_task(task390);

  vector<IndexRange> I452_index = {virt_, active_};
  auto I452 = make_shared<Tensor>(I452_index);
  auto tensor391 = vector<shared_ptr<Tensor>>{den2, I452};
  auto task391 = make_shared<Task391>(tensor391, pindex);
  task391->add_dep(task283);
  densityq->add_task(task391);

  vector<IndexRange> I453_index = {active_, virt_};
  auto I453 = make_shared<Tensor>(I453_index);
  auto tensor392 = vector<shared_ptr<Tensor>>{I452, Gamma16_(), I453};
  auto task392 = make_shared<Task392>(tensor392, pindex);
  task391->add_dep(task392);
  task392->add_dep(task283);
  densityq->add_task(task392);

  auto tensor393 = vector<shared_ptr<Tensor>>{I453, t2};
  auto task393 = make_shared<Task393>(tensor393, pindex);
  task392->add_dep(task393);
  task393->add_dep(task283);
  densityq->add_task(task393);

  auto tensor394 = vector<shared_ptr<Tensor>>{I453, t2};
  auto task394 = make_shared<Task394>(tensor394, pindex);
  task392->add_dep(task394);
  task394->add_dep(task283);
  densityq->add_task(task394);

  vector<IndexRange> I458_index = {active_, virt_};
  auto I458 = make_shared<Tensor>(I458_index);
  auto tensor395 = vector<shared_ptr<Tensor>>{den2, I458};
  auto task395 = make_shared<Task395>(tensor395, pindex);
  task395->add_dep(task283);
  densityq->add_task(task395);

  vector<IndexRange> I459_index = {closed_, active_, active_, active_};
  auto I459 = make_shared<Tensor>(I459_index);
  auto tensor396 = vector<shared_ptr<Tensor>>{I458, t2, I459};
  auto task396 = make_shared<Task396>(tensor396, pindex);
  task395->add_dep(task396);
  task396->add_dep(task283);
  densityq->add_task(task396);

  auto tensor397 = vector<shared_ptr<Tensor>>{I459, Gamma28_(), t2};
  auto task397 = make_shared<Task397>(tensor397, pindex);
  task396->add_dep(task397);
  task397->add_dep(task283);
  densityq->add_task(task397);

  vector<IndexRange> I461_index = {active_, closed_};
  auto I461 = make_shared<Tensor>(I461_index);
  auto tensor398 = vector<shared_ptr<Tensor>>{den2, I461};
  auto task398 = make_shared<Task398>(tensor398, pindex);
  task398->add_dep(task283);
  densityq->add_task(task398);

  vector<IndexRange> I462_index = {closed_, virt_, active_, active_};
  auto I462 = make_shared<Tensor>(I462_index);
  auto tensor399 = vector<shared_ptr<Tensor>>{I461, t2, I462};
  auto task399 = make_shared<Task399>(tensor399, pindex);
  task398->add_dep(task399);
  task399->add_dep(task283);
  densityq->add_task(task399);

  auto tensor400 = vector<shared_ptr<Tensor>>{I462, Gamma29_(), t2};
  auto task400 = make_shared<Task400>(tensor400, pindex);
  task399->add_dep(task400);
  task400->add_dep(task283);
  densityq->add_task(task400);

  vector<IndexRange> I465_index = {closed_, virt_, active_, active_};
  auto I465 = make_shared<Tensor>(I465_index);
  auto tensor401 = vector<shared_ptr<Tensor>>{I461, t2, I465};
  auto task401 = make_shared<Task401>(tensor401, pindex);
  task398->add_dep(task401);
  task401->add_dep(task283);
  densityq->add_task(task401);

  auto tensor402 = vector<shared_ptr<Tensor>>{I465, Gamma7_(), t2};
  auto task402 = make_shared<Task402>(tensor402, pindex);
  task401->add_dep(task402);
  task402->add_dep(task283);
  densityq->add_task(task402);

  vector<IndexRange> I504_index = {virt_, closed_, active_, active_};
  auto I504 = make_shared<Tensor>(I504_index);
  auto tensor403 = vector<shared_ptr<Tensor>>{I461, t2, I504};
  auto task403 = make_shared<Task403>(tensor403, pindex);
  task398->add_dep(task403);
  task403->add_dep(task283);
  densityq->add_task(task403);

  auto tensor404 = vector<shared_ptr<Tensor>>{I504, Gamma7_(), t2};
  auto task404 = make_shared<Task404>(tensor404, pindex);
  task403->add_dep(task404);
  task404->add_dep(task283);
  densityq->add_task(task404);

  vector<IndexRange> I507_index = {virt_, closed_, active_, active_};
  auto I507 = make_shared<Tensor>(I507_index);
  auto tensor405 = vector<shared_ptr<Tensor>>{I461, t2, I507};
  auto task405 = make_shared<Task405>(tensor405, pindex);
  task398->add_dep(task405);
  task405->add_dep(task283);
  densityq->add_task(task405);

  auto tensor406 = vector<shared_ptr<Tensor>>{I507, Gamma7_(), t2};
  auto task406 = make_shared<Task406>(tensor406, pindex);
  task405->add_dep(task406);
  task406->add_dep(task283);
  densityq->add_task(task406);

  vector<IndexRange> I656_index = {virt_, virt_, active_, active_};
  auto I656 = make_shared<Tensor>(I656_index);
  auto tensor407 = vector<shared_ptr<Tensor>>{I461, t2, I656};
  auto task407 = make_shared<Task407>(tensor407, pindex);
  task398->add_dep(task407);
  task407->add_dep(task283);
  densityq->add_task(task407);

  auto tensor408 = vector<shared_ptr<Tensor>>{I656, Gamma60_(), t2};
  auto task408 = make_shared<Task408>(tensor408, pindex);
  task407->add_dep(task408);
  task408->add_dep(task283);
  densityq->add_task(task408);

  vector<IndexRange> I473_index = {virt_, virt_};
  auto I473 = make_shared<Tensor>(I473_index);
  auto tensor409 = vector<shared_ptr<Tensor>>{den2, I473};
  auto task409 = make_shared<Task409>(tensor409, pindex);
  task409->add_dep(task283);
  densityq->add_task(task409);

  vector<IndexRange> I474_index = {closed_, virt_, active_, active_};
  auto I474 = make_shared<Tensor>(I474_index);
  auto tensor410 = vector<shared_ptr<Tensor>>{I473, t2, I474};
  auto task410 = make_shared<Task410>(tensor410, pindex);
  task409->add_dep(task410);
  task410->add_dep(task283);
  densityq->add_task(task410);

  auto tensor411 = vector<shared_ptr<Tensor>>{I474, Gamma32_(), t2};
  auto task411 = make_shared<Task411>(tensor411, pindex);
  task410->add_dep(task411);
  task411->add_dep(task283);
  densityq->add_task(task411);

  vector<IndexRange> I483_index = {closed_, virt_, active_, active_};
  auto I483 = make_shared<Tensor>(I483_index);
  auto tensor412 = vector<shared_ptr<Tensor>>{I473, t2, I483};
  auto task412 = make_shared<Task412>(tensor412, pindex);
  task409->add_dep(task412);
  task412->add_dep(task283);
  densityq->add_task(task412);

  auto tensor413 = vector<shared_ptr<Tensor>>{I483, Gamma35_(), t2};
  auto task413 = make_shared<Task413>(tensor413, pindex);
  task412->add_dep(task413);
  task413->add_dep(task283);
  densityq->add_task(task413);

  vector<IndexRange> I488_index = {virt_, closed_};
  auto I488 = make_shared<Tensor>(I488_index);
  auto tensor414 = vector<shared_ptr<Tensor>>{den2, I488};
  auto task414 = make_shared<Task414>(tensor414, pindex);
  task414->add_dep(task283);
  densityq->add_task(task414);

  vector<IndexRange> I489_index = {closed_, virt_};
  auto I489 = make_shared<Tensor>(I489_index);
  auto tensor415 = vector<shared_ptr<Tensor>>{I488, t2, I489};
  auto task415 = make_shared<Task415>(tensor415, pindex);
  task414->add_dep(task415);
  task415->add_dep(task283);
  densityq->add_task(task415);

  auto tensor416 = vector<shared_ptr<Tensor>>{I489, Gamma38_(), t2};
  auto task416 = make_shared<Task416>(tensor416, pindex);
  task415->add_dep(task416);
  task416->add_dep(task283);
  densityq->add_task(task416);

  vector<IndexRange> I492_index = {closed_, virt_};
  auto I492 = make_shared<Tensor>(I492_index);
  auto tensor417 = vector<shared_ptr<Tensor>>{I488, t2, I492};
  auto task417 = make_shared<Task417>(tensor417, pindex);
  task414->add_dep(task417);
  task417->add_dep(task283);
  densityq->add_task(task417);

  auto tensor418 = vector<shared_ptr<Tensor>>{I492, Gamma38_(), t2};
  auto task418 = make_shared<Task418>(tensor418, pindex);
  task417->add_dep(task418);
  task418->add_dep(task283);
  densityq->add_task(task418);

  vector<IndexRange> I531_index = {virt_, closed_};
  auto I531 = make_shared<Tensor>(I531_index);
  auto tensor419 = vector<shared_ptr<Tensor>>{I488, t2, I531};
  auto task419 = make_shared<Task419>(tensor419, pindex);
  task414->add_dep(task419);
  task419->add_dep(task283);
  densityq->add_task(task419);

  auto tensor420 = vector<shared_ptr<Tensor>>{I531, Gamma38_(), t2};
  auto task420 = make_shared<Task420>(tensor420, pindex);
  task419->add_dep(task420);
  task420->add_dep(task283);
  densityq->add_task(task420);

  vector<IndexRange> I534_index = {virt_, closed_};
  auto I534 = make_shared<Tensor>(I534_index);
  auto tensor421 = vector<shared_ptr<Tensor>>{I488, t2, I534};
  auto task421 = make_shared<Task421>(tensor421, pindex);
  task414->add_dep(task421);
  task421->add_dep(task283);
  densityq->add_task(task421);

  auto tensor422 = vector<shared_ptr<Tensor>>{I534, Gamma38_(), t2};
  auto task422 = make_shared<Task422>(tensor422, pindex);
  task421->add_dep(task422);
  task422->add_dep(task283);
  densityq->add_task(task422);

  vector<IndexRange> I500_index = {active_, virt_};
  auto I500 = make_shared<Tensor>(I500_index);
  auto tensor423 = vector<shared_ptr<Tensor>>{den2, I500};
  auto task423 = make_shared<Task423>(tensor423, pindex);
  task423->add_dep(task283);
  densityq->add_task(task423);

  vector<IndexRange> I501_index = {closed_, active_, active_, active_};
  auto I501 = make_shared<Tensor>(I501_index);
  auto tensor424 = vector<shared_ptr<Tensor>>{I500, t2, I501};
  auto task424 = make_shared<Task424>(tensor424, pindex);
  task423->add_dep(task424);
  task424->add_dep(task283);
  densityq->add_task(task424);

  auto tensor425 = vector<shared_ptr<Tensor>>{I501, Gamma6_(), t2};
  auto task425 = make_shared<Task425>(tensor425, pindex);
  task424->add_dep(task425);
  task425->add_dep(task283);
  densityq->add_task(task425);

  vector<IndexRange> I653_index = {virt_, active_, active_, active_};
  auto I653 = make_shared<Tensor>(I653_index);
  auto tensor426 = vector<shared_ptr<Tensor>>{I500, t2, I653};
  auto task426 = make_shared<Task426>(tensor426, pindex);
  task423->add_dep(task426);
  task426->add_dep(task283);
  densityq->add_task(task426);

  auto tensor427 = vector<shared_ptr<Tensor>>{I653, Gamma59_(), t2};
  auto task427 = make_shared<Task427>(tensor427, pindex);
  task426->add_dep(task427);
  task427->add_dep(task283);
  densityq->add_task(task427);

  vector<IndexRange> I512_index = {closed_, closed_};
  auto I512 = make_shared<Tensor>(I512_index);
  auto tensor428 = vector<shared_ptr<Tensor>>{den2, I512};
  auto task428 = make_shared<Task428>(tensor428, pindex);
  task428->add_dep(task283);
  densityq->add_task(task428);

  vector<IndexRange> I513_index = {virt_, closed_, active_, active_};
  auto I513 = make_shared<Tensor>(I513_index);
  auto tensor429 = vector<shared_ptr<Tensor>>{I512, t2, I513};
  auto task429 = make_shared<Task429>(tensor429, pindex);
  task428->add_dep(task429);
  task429->add_dep(task283);
  densityq->add_task(task429);

  auto tensor430 = vector<shared_ptr<Tensor>>{I513, Gamma35_(), t2};
  auto task430 = make_shared<Task430>(tensor430, pindex);
  task429->add_dep(task430);
  task430->add_dep(task283);
  densityq->add_task(task430);

  vector<IndexRange> I522_index = {virt_, closed_, active_, active_};
  auto I522 = make_shared<Tensor>(I522_index);
  auto tensor431 = vector<shared_ptr<Tensor>>{I512, t2, I522};
  auto task431 = make_shared<Task431>(tensor431, pindex);
  task428->add_dep(task431);
  task431->add_dep(task283);
  densityq->add_task(task431);

  auto tensor432 = vector<shared_ptr<Tensor>>{I522, Gamma35_(), t2};
  auto task432 = make_shared<Task432>(tensor432, pindex);
  task431->add_dep(task432);
  task432->add_dep(task283);
  densityq->add_task(task432);

  vector<IndexRange> I515_index = {virt_, virt_};
  auto I515 = make_shared<Tensor>(I515_index);
  auto tensor433 = vector<shared_ptr<Tensor>>{den2, I515};
  auto task433 = make_shared<Task433>(tensor433, pindex);
  task433->add_dep(task283);
  densityq->add_task(task433);

  vector<IndexRange> I516_index = {virt_, closed_, active_, active_};
  auto I516 = make_shared<Tensor>(I516_index);
  auto tensor434 = vector<shared_ptr<Tensor>>{I515, t2, I516};
  auto task434 = make_shared<Task434>(tensor434, pindex);
  task433->add_dep(task434);
  task434->add_dep(task283);
  densityq->add_task(task434);

  auto tensor435 = vector<shared_ptr<Tensor>>{I516, Gamma35_(), t2};
  auto task435 = make_shared<Task435>(tensor435, pindex);
  task434->add_dep(task435);
  task435->add_dep(task283);
  densityq->add_task(task435);

  vector<IndexRange> I525_index = {virt_, closed_, active_, active_};
  auto I525 = make_shared<Tensor>(I525_index);
  auto tensor436 = vector<shared_ptr<Tensor>>{I515, t2, I525};
  auto task436 = make_shared<Task436>(tensor436, pindex);
  task433->add_dep(task436);
  task436->add_dep(task283);
  densityq->add_task(task436);

  auto tensor437 = vector<shared_ptr<Tensor>>{I525, Gamma35_(), t2};
  auto task437 = make_shared<Task437>(tensor437, pindex);
  task436->add_dep(task437);
  task437->add_dep(task283);
  densityq->add_task(task437);

  vector<IndexRange> I662_index = {virt_, virt_, active_, active_};
  auto I662 = make_shared<Tensor>(I662_index);
  auto tensor438 = vector<shared_ptr<Tensor>>{I515, t2, I662};
  auto task438 = make_shared<Task438>(tensor438, pindex);
  task433->add_dep(task438);
  task438->add_dep(task283);
  densityq->add_task(task438);

  auto tensor439 = vector<shared_ptr<Tensor>>{I662, Gamma60_(), t2};
  auto task439 = make_shared<Task439>(tensor439, pindex);
  task438->add_dep(task439);
  task439->add_dep(task283);
  densityq->add_task(task439);

  vector<IndexRange> I527_index = {active_, closed_};
  auto I527 = make_shared<Tensor>(I527_index);
  auto tensor440 = vector<shared_ptr<Tensor>>{den2, I527};
  auto task440 = make_shared<Task440>(tensor440, pindex);
  task440->add_dep(task283);
  densityq->add_task(task440);

  vector<IndexRange> I528_index = {virt_, active_, active_, active_};
  auto I528 = make_shared<Tensor>(I528_index);
  auto tensor441 = vector<shared_ptr<Tensor>>{I527, t2, I528};
  auto task441 = make_shared<Task441>(tensor441, pindex);
  task440->add_dep(task441);
  task441->add_dep(task283);
  densityq->add_task(task441);

  auto tensor442 = vector<shared_ptr<Tensor>>{I528, Gamma51_(), t2};
  auto task442 = make_shared<Task442>(tensor442, pindex);
  task441->add_dep(task442);
  task442->add_dep(task283);
  densityq->add_task(task442);

  vector<IndexRange> I551_index = {virt_, virt_};
  auto I551 = make_shared<Tensor>(I551_index);
  auto tensor443 = vector<shared_ptr<Tensor>>{den2, I551};
  auto task443 = make_shared<Task443>(tensor443, pindex);
  task443->add_dep(task283);
  densityq->add_task(task443);

  vector<IndexRange> I552_index = {virt_, active_, active_, active_};
  auto I552 = make_shared<Tensor>(I552_index);
  auto tensor444 = vector<shared_ptr<Tensor>>{I551, t2, I552};
  auto task444 = make_shared<Task444>(tensor444, pindex);
  task443->add_dep(task444);
  task444->add_dep(task283);
  densityq->add_task(task444);

  auto tensor445 = vector<shared_ptr<Tensor>>{I552, Gamma59_(), t2};
  auto task445 = make_shared<Task445>(tensor445, pindex);
  task444->add_dep(task445);
  task445->add_dep(task283);
  densityq->add_task(task445);

  vector<IndexRange> I563_index = {virt_, active_};
  auto I563 = make_shared<Tensor>(I563_index);
  auto tensor446 = vector<shared_ptr<Tensor>>{den2, I563};
  auto task446 = make_shared<Task446>(tensor446, pindex);
  task446->add_dep(task283);
  densityq->add_task(task446);

  vector<IndexRange> I564_index = {virt_, active_};
  auto I564 = make_shared<Tensor>(I564_index);
  auto tensor447 = vector<shared_ptr<Tensor>>{I563, Gamma16_(), I564};
  auto task447 = make_shared<Task447>(tensor447, pindex);
  task446->add_dep(task447);
  task447->add_dep(task283);
  densityq->add_task(task447);

  auto tensor448 = vector<shared_ptr<Tensor>>{I564, t2};
  auto task448 = make_shared<Task448>(tensor448, pindex);
  task447->add_dep(task448);
  task448->add_dep(task283);
  densityq->add_task(task448);

  vector<IndexRange> I566_index = {virt_, active_};
  auto I566 = make_shared<Tensor>(I566_index);
  auto tensor449 = vector<shared_ptr<Tensor>>{den2, I566};
  auto task449 = make_shared<Task449>(tensor449, pindex);
  task449->add_dep(task283);
  densityq->add_task(task449);

  vector<IndexRange> I567_index = {virt_, active_};
  auto I567 = make_shared<Tensor>(I567_index);
  auto tensor450 = vector<shared_ptr<Tensor>>{I566, Gamma16_(), I567};
  auto task450 = make_shared<Task450>(tensor450, pindex);
  task449->add_dep(task450);
  task450->add_dep(task283);
  densityq->add_task(task450);

  auto tensor451 = vector<shared_ptr<Tensor>>{I567, t2};
  auto task451 = make_shared<Task451>(tensor451, pindex);
  task450->add_dep(task451);
  task451->add_dep(task283);
  densityq->add_task(task451);

  vector<IndexRange> I572_index = {virt_, closed_};
  auto I572 = make_shared<Tensor>(I572_index);
  auto tensor452 = vector<shared_ptr<Tensor>>{den2, I572};
  auto task452 = make_shared<Task452>(tensor452, pindex);
  task452->add_dep(task283);
  densityq->add_task(task452);

  vector<IndexRange> I573_index = {virt_, closed_};
  auto I573 = make_shared<Tensor>(I573_index);
  auto tensor453 = vector<shared_ptr<Tensor>>{I572, t2, I573};
  auto task453 = make_shared<Task453>(tensor453, pindex);
  task452->add_dep(task453);
  task453->add_dep(task283);
  densityq->add_task(task453);

  vector<IndexRange> I574_index = {active_, virt_, closed_, active_};
  auto I574 = make_shared<Tensor>(I574_index);
  auto tensor454 = vector<shared_ptr<Tensor>>{I573, Gamma38_(), I574};
  auto task454 = make_shared<Task454>(tensor454, pindex);
  task453->add_dep(task454);
  task454->add_dep(task283);
  densityq->add_task(task454);

  auto tensor455 = vector<shared_ptr<Tensor>>{I574, t2};
  auto task455 = make_shared<Task455>(tensor455, pindex);
  task454->add_dep(task455);
  task455->add_dep(task283);
  densityq->add_task(task455);

  vector<IndexRange> I581_index = {active_, active_};
  auto I581 = make_shared<Tensor>(I581_index);
  auto tensor456 = vector<shared_ptr<Tensor>>{den2, I581};
  auto task456 = make_shared<Task456>(tensor456, pindex);
  task456->add_dep(task283);
  densityq->add_task(task456);

  vector<IndexRange> I582_index;
  auto I582 = make_shared<Tensor>(I582_index);
  auto tensor457 = vector<shared_ptr<Tensor>>{I581, Gamma38_(), I582};
  auto task457 = make_shared<Task457>(tensor457, pindex);
  task456->add_dep(task457);
  task457->add_dep(task283);
  densityq->add_task(task457);

  auto tensor458 = vector<shared_ptr<Tensor>>{I582, t2};
  auto task458 = make_shared<Task458>(tensor458, pindex);
  task457->add_dep(task458);
  task458->add_dep(task283);
  densityq->add_task(task458);

  auto tensor459 = vector<shared_ptr<Tensor>>{I582, t2};
  auto task459 = make_shared<Task459>(tensor459, pindex);
  task457->add_dep(task459);
  task459->add_dep(task283);
  densityq->add_task(task459);

  shared_ptr<Tensor> I587;
  if (diagonal) {
    vector<IndexRange> I587_index = {closed_, closed_};
    I587 = make_shared<Tensor>(I587_index);
  }
  shared_ptr<Task460> task460;
  if (diagonal) {
    auto tensor460 = vector<shared_ptr<Tensor>>{den2, I587};
    task460 = make_shared<Task460>(tensor460, pindex);
    task460->add_dep(task283);
    densityq->add_task(task460);
  }

  shared_ptr<Task461> task461;
  if (diagonal) {
    auto tensor461 = vector<shared_ptr<Tensor>>{I587, t2};
    task461 = make_shared<Task461>(tensor461, pindex);
    task460->add_dep(task461);
    task461->add_dep(task283);
    densityq->add_task(task461);
  }

  shared_ptr<Task462> task462;
  if (diagonal) {
    auto tensor462 = vector<shared_ptr<Tensor>>{I587, t2};
    task462 = make_shared<Task462>(tensor462, pindex);
    task460->add_dep(task462);
    task462->add_dep(task283);
    densityq->add_task(task462);
  }

  shared_ptr<Tensor> I591;
  if (diagonal) {
    vector<IndexRange> I591_index = {virt_, virt_};
    I591 = make_shared<Tensor>(I591_index);
  }
  shared_ptr<Task463> task463;
  if (diagonal) {
    auto tensor463 = vector<shared_ptr<Tensor>>{den2, I591};
    task463 = make_shared<Task463>(tensor463, pindex);
    task463->add_dep(task283);
    densityq->add_task(task463);
  }

  shared_ptr<Task464> task464;
  if (diagonal) {
    auto tensor464 = vector<shared_ptr<Tensor>>{I591, t2};
    task464 = make_shared<Task464>(tensor464, pindex);
    task463->add_dep(task464);
    task464->add_dep(task283);
    densityq->add_task(task464);
  }

  shared_ptr<Task465> task465;
  if (diagonal) {
    auto tensor465 = vector<shared_ptr<Tensor>>{I591, t2};
    task465 = make_shared<Task465>(tensor465, pindex);
    task463->add_dep(task465);
    task465->add_dep(task283);
    densityq->add_task(task465);
  }

  vector<IndexRange> I595_index = {closed_, active_};
  auto I595 = make_shared<Tensor>(I595_index);
  auto tensor466 = vector<shared_ptr<Tensor>>{den2, I595};
  auto task466 = make_shared<Task466>(tensor466, pindex);
  task466->add_dep(task283);
  densityq->add_task(task466);

  vector<IndexRange> I596_index = {closed_, active_};
  auto I596 = make_shared<Tensor>(I596_index);
  auto tensor467 = vector<shared_ptr<Tensor>>{I595, Gamma38_(), I596};
  auto task467 = make_shared<Task467>(tensor467, pindex);
  task466->add_dep(task467);
  task467->add_dep(task283);
  densityq->add_task(task467);

  auto tensor468 = vector<shared_ptr<Tensor>>{I596, t2};
  auto task468 = make_shared<Task468>(tensor468, pindex);
  task467->add_dep(task468);
  task468->add_dep(task283);
  densityq->add_task(task468);

  auto tensor469 = vector<shared_ptr<Tensor>>{I596, t2};
  auto task469 = make_shared<Task469>(tensor469, pindex);
  task467->add_dep(task469);
  task469->add_dep(task283);
  densityq->add_task(task469);

  vector<IndexRange> I601_index = {active_, virt_};
  auto I601 = make_shared<Tensor>(I601_index);
  auto tensor470 = vector<shared_ptr<Tensor>>{den2, I601};
  auto task470 = make_shared<Task470>(tensor470, pindex);
  task470->add_dep(task283);
  densityq->add_task(task470);

  vector<IndexRange> I602_index = {virt_, closed_, active_, active_};
  auto I602 = make_shared<Tensor>(I602_index);
  auto tensor471 = vector<shared_ptr<Tensor>>{I601, t2, I602};
  auto task471 = make_shared<Task471>(tensor471, pindex);
  task470->add_dep(task471);
  task471->add_dep(task283);
  densityq->add_task(task471);

  vector<IndexRange> I603_index = {active_, virt_, closed_, active_};
  auto I603 = make_shared<Tensor>(I603_index);
  auto tensor472 = vector<shared_ptr<Tensor>>{I602, Gamma35_(), I603};
  auto task472 = make_shared<Task472>(tensor472, pindex);
  task471->add_dep(task472);
  task472->add_dep(task283);
  densityq->add_task(task472);

  auto tensor473 = vector<shared_ptr<Tensor>>{I603, t2};
  auto task473 = make_shared<Task473>(tensor473, pindex);
  task472->add_dep(task473);
  task473->add_dep(task283);
  densityq->add_task(task473);

  vector<IndexRange> I604_index = {active_, virt_};
  auto I604 = make_shared<Tensor>(I604_index);
  auto tensor474 = vector<shared_ptr<Tensor>>{den2, I604};
  auto task474 = make_shared<Task474>(tensor474, pindex);
  task474->add_dep(task283);
  densityq->add_task(task474);

  vector<IndexRange> I605_index = {virt_, closed_, active_, active_};
  auto I605 = make_shared<Tensor>(I605_index);
  auto tensor475 = vector<shared_ptr<Tensor>>{I604, t2, I605};
  auto task475 = make_shared<Task475>(tensor475, pindex);
  task474->add_dep(task475);
  task475->add_dep(task283);
  densityq->add_task(task475);

  auto tensor476 = vector<shared_ptr<Tensor>>{I605, Gamma32_(), t2};
  auto task476 = make_shared<Task476>(tensor476, pindex);
  task475->add_dep(task476);
  task476->add_dep(task283);
  densityq->add_task(task476);

  auto tensor477 = vector<shared_ptr<Tensor>>{I605, Gamma35_(), t2};
  auto task477 = make_shared<Task477>(tensor477, pindex);
  task475->add_dep(task477);
  task477->add_dep(task283);
  densityq->add_task(task477);

  vector<IndexRange> I613_index = {closed_, virt_};
  auto I613 = make_shared<Tensor>(I613_index);
  auto tensor478 = vector<shared_ptr<Tensor>>{den2, I613};
  auto task478 = make_shared<Task478>(tensor478, pindex);
  task478->add_dep(task283);
  densityq->add_task(task478);

  vector<IndexRange> I614_index = {virt_, active_};
  auto I614 = make_shared<Tensor>(I614_index);
  auto tensor479 = vector<shared_ptr<Tensor>>{I613, t2, I614};
  auto task479 = make_shared<Task479>(tensor479, pindex);
  task478->add_dep(task479);
  task479->add_dep(task283);
  densityq->add_task(task479);

  auto tensor480 = vector<shared_ptr<Tensor>>{I614, Gamma60_(), t2};
  auto task480 = make_shared<Task480>(tensor480, pindex);
  task479->add_dep(task480);
  task480->add_dep(task283);
  densityq->add_task(task480);

  vector<IndexRange> I616_index = {virt_, closed_};
  auto I616 = make_shared<Tensor>(I616_index);
  auto tensor481 = vector<shared_ptr<Tensor>>{den2, I616};
  auto task481 = make_shared<Task481>(tensor481, pindex);
  task481->add_dep(task283);
  densityq->add_task(task481);

  vector<IndexRange> I617_index = {virt_, active_};
  auto I617 = make_shared<Tensor>(I617_index);
  auto tensor482 = vector<shared_ptr<Tensor>>{I616, t2, I617};
  auto task482 = make_shared<Task482>(tensor482, pindex);
  task481->add_dep(task482);
  task482->add_dep(task283);
  densityq->add_task(task482);

  auto tensor483 = vector<shared_ptr<Tensor>>{I617, Gamma60_(), t2};
  auto task483 = make_shared<Task483>(tensor483, pindex);
  task482->add_dep(task483);
  task483->add_dep(task283);
  densityq->add_task(task483);

  vector<IndexRange> I619_index = {closed_, active_};
  auto I619 = make_shared<Tensor>(I619_index);
  auto tensor484 = vector<shared_ptr<Tensor>>{den2, I619};
  auto task484 = make_shared<Task484>(tensor484, pindex);
  task484->add_dep(task283);
  densityq->add_task(task484);

  vector<IndexRange> I620_index = {active_, closed_};
  auto I620 = make_shared<Tensor>(I620_index);
  auto tensor485 = vector<shared_ptr<Tensor>>{I619, Gamma38_(), I620};
  auto task485 = make_shared<Task485>(tensor485, pindex);
  task484->add_dep(task485);
  task485->add_dep(task283);
  densityq->add_task(task485);

  auto tensor486 = vector<shared_ptr<Tensor>>{I620, t2};
  auto task486 = make_shared<Task486>(tensor486, pindex);
  task485->add_dep(task486);
  task486->add_dep(task283);
  densityq->add_task(task486);

  auto tensor487 = vector<shared_ptr<Tensor>>{I620, t2};
  auto task487 = make_shared<Task487>(tensor487, pindex);
  task485->add_dep(task487);
  task487->add_dep(task283);
  densityq->add_task(task487);

  vector<IndexRange> I631_index = {closed_, closed_};
  auto I631 = make_shared<Tensor>(I631_index);
  auto tensor488 = vector<shared_ptr<Tensor>>{den2, I631};
  auto task488 = make_shared<Task488>(tensor488, pindex);
  task488->add_dep(task283);
  densityq->add_task(task488);

  vector<IndexRange> I632_index = {virt_, closed_, virt_, active_};
  auto I632 = make_shared<Tensor>(I632_index);
  auto tensor489 = vector<shared_ptr<Tensor>>{I631, t2, I632};
  auto task489 = make_shared<Task489>(tensor489, pindex);
  task488->add_dep(task489);
  task489->add_dep(task283);
  densityq->add_task(task489);

  auto tensor490 = vector<shared_ptr<Tensor>>{I632, Gamma38_(), t2};
  auto task490 = make_shared<Task490>(tensor490, pindex);
  task489->add_dep(task490);
  task490->add_dep(task283);
  densityq->add_task(task490);

  vector<IndexRange> I635_index = {virt_, closed_, virt_, active_};
  auto I635 = make_shared<Tensor>(I635_index);
  auto tensor491 = vector<shared_ptr<Tensor>>{I631, t2, I635};
  auto task491 = make_shared<Task491>(tensor491, pindex);
  task488->add_dep(task491);
  task491->add_dep(task283);
  densityq->add_task(task491);

  auto tensor492 = vector<shared_ptr<Tensor>>{I635, Gamma38_(), t2};
  auto task492 = make_shared<Task492>(tensor492, pindex);
  task491->add_dep(task492);
  task492->add_dep(task283);
  densityq->add_task(task492);

  vector<IndexRange> I637_index = {virt_, virt_};
  auto I637 = make_shared<Tensor>(I637_index);
  auto tensor493 = vector<shared_ptr<Tensor>>{den2, I637};
  auto task493 = make_shared<Task493>(tensor493, pindex);
  task493->add_dep(task283);
  densityq->add_task(task493);

  vector<IndexRange> I638_index = {virt_, closed_, virt_, active_};
  auto I638 = make_shared<Tensor>(I638_index);
  auto tensor494 = vector<shared_ptr<Tensor>>{I637, t2, I638};
  auto task494 = make_shared<Task494>(tensor494, pindex);
  task493->add_dep(task494);
  task494->add_dep(task283);
  densityq->add_task(task494);

  auto tensor495 = vector<shared_ptr<Tensor>>{I638, Gamma38_(), t2};
  auto task495 = make_shared<Task495>(tensor495, pindex);
  task494->add_dep(task495);
  task495->add_dep(task283);
  densityq->add_task(task495);

  vector<IndexRange> I641_index = {virt_, closed_, virt_, active_};
  auto I641 = make_shared<Tensor>(I641_index);
  auto tensor496 = vector<shared_ptr<Tensor>>{I637, t2, I641};
  auto task496 = make_shared<Task496>(tensor496, pindex);
  task493->add_dep(task496);
  task496->add_dep(task283);
  densityq->add_task(task496);

  auto tensor497 = vector<shared_ptr<Tensor>>{I641, Gamma38_(), t2};
  auto task497 = make_shared<Task497>(tensor497, pindex);
  task496->add_dep(task497);
  task497->add_dep(task283);
  densityq->add_task(task497);

  vector<IndexRange> I643_index = {virt_, virt_};
  auto I643 = make_shared<Tensor>(I643_index);
  auto tensor498 = vector<shared_ptr<Tensor>>{den2, I643};
  auto task498 = make_shared<Task498>(tensor498, pindex);
  task498->add_dep(task283);
  densityq->add_task(task498);

  vector<IndexRange> I644_index = {virt_, closed_, virt_, active_};
  auto I644 = make_shared<Tensor>(I644_index);
  auto tensor499 = vector<shared_ptr<Tensor>>{I643, t2, I644};
  auto task499 = make_shared<Task499>(tensor499, pindex);
  task498->add_dep(task499);
  task499->add_dep(task283);
  densityq->add_task(task499);

  auto tensor500 = vector<shared_ptr<Tensor>>{I644, Gamma38_(), t2};
  auto task500 = make_shared<Task500>(tensor500, pindex);
  task499->add_dep(task500);
  task500->add_dep(task283);
  densityq->add_task(task500);

  vector<IndexRange> I647_index = {virt_, closed_, virt_, active_};
  auto I647 = make_shared<Tensor>(I647_index);
  auto tensor501 = vector<shared_ptr<Tensor>>{I643, t2, I647};
  auto task501 = make_shared<Task501>(tensor501, pindex);
  task498->add_dep(task501);
  task501->add_dep(task283);
  densityq->add_task(task501);

  auto tensor502 = vector<shared_ptr<Tensor>>{I647, Gamma38_(), t2};
  auto task502 = make_shared<Task502>(tensor502, pindex);
  task501->add_dep(task502);
  task502->add_dep(task283);
  densityq->add_task(task502);

  vector<IndexRange> I649_index = {active_, closed_};
  auto I649 = make_shared<Tensor>(I649_index);
  auto tensor503 = vector<shared_ptr<Tensor>>{den2, I649};
  auto task503 = make_shared<Task503>(tensor503, pindex);
  task503->add_dep(task283);
  densityq->add_task(task503);

  vector<IndexRange> I650_index = {virt_, virt_, active_, active_};
  auto I650 = make_shared<Tensor>(I650_index);
  auto tensor504 = vector<shared_ptr<Tensor>>{I649, t2, I650};
  auto task504 = make_shared<Task504>(tensor504, pindex);
  task503->add_dep(task504);
  task504->add_dep(task283);
  densityq->add_task(task504);

  auto tensor505 = vector<shared_ptr<Tensor>>{I650, Gamma60_(), t2};
  auto task505 = make_shared<Task505>(tensor505, pindex);
  task504->add_dep(task505);
  task505->add_dep(task283);
  densityq->add_task(task505);

  return densityq;
}


#endif
