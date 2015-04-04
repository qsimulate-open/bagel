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

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_densityq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto densityq = make_shared<Queue>();
  auto tensor299 = vector<shared_ptr<Tensor>>{den2};
  auto task299 = make_shared<Task299>(tensor299, reset);
  densityq->add_task(task299);

  vector<IndexRange> I444_index = {active_, active_};
  auto I444 = make_shared<Tensor>(I444_index);
  auto tensor300 = vector<shared_ptr<Tensor>>{den2, I444};
  auto task300 = make_shared<Task300>(tensor300, pindex);
  task300->add_dep(task299);
  densityq->add_task(task300);

  vector<IndexRange> I445_index = {active_, active_, active_, active_};
  auto I445 = make_shared<Tensor>(I445_index);
  auto tensor301 = vector<shared_ptr<Tensor>>{I444, Gamma160_(), I445};
  auto task301 = make_shared<Task301>(tensor301, pindex);
  task300->add_dep(task301);
  task301->add_dep(task299);
  densityq->add_task(task301);

  auto tensor302 = vector<shared_ptr<Tensor>>{I445, t2};
  auto task302 = make_shared<Task302>(tensor302, pindex);
  task301->add_dep(task302);
  task302->add_dep(task299);
  densityq->add_task(task302);

  vector<IndexRange> I538_index = {active_, active_, active_, active_};
  auto I538 = make_shared<Tensor>(I538_index);
  auto tensor303 = vector<shared_ptr<Tensor>>{I444, Gamma191_(), I538};
  auto task303 = make_shared<Task303>(tensor303, pindex);
  task300->add_dep(task303);
  task303->add_dep(task299);
  densityq->add_task(task303);

  auto tensor304 = vector<shared_ptr<Tensor>>{I538, t2};
  auto task304 = make_shared<Task304>(tensor304, pindex);
  task303->add_dep(task304);
  task304->add_dep(task299);
  densityq->add_task(task304);

  vector<IndexRange> I547_index = {active_, active_, active_, active_};
  auto I547 = make_shared<Tensor>(I547_index);
  auto tensor305 = vector<shared_ptr<Tensor>>{I444, Gamma194_(), I547};
  auto task305 = make_shared<Task305>(tensor305, pindex);
  task300->add_dep(task305);
  task305->add_dep(task299);
  densityq->add_task(task305);

  auto tensor306 = vector<shared_ptr<Tensor>>{I547, t2};
  auto task306 = make_shared<Task306>(tensor306, pindex);
  task305->add_dep(task306);
  task306->add_dep(task299);
  densityq->add_task(task306);

  auto tensor307 = vector<shared_ptr<Tensor>>{I547, t2};
  auto task307 = make_shared<Task307>(tensor307, pindex);
  task305->add_dep(task307);
  task307->add_dep(task299);
  densityq->add_task(task307);

  auto tensor308 = vector<shared_ptr<Tensor>>{I547, t2};
  auto task308 = make_shared<Task308>(tensor308, pindex);
  task305->add_dep(task308);
  task308->add_dep(task299);
  densityq->add_task(task308);

  vector<IndexRange> I729_index = {active_, active_, active_, active_};
  auto I729 = make_shared<Tensor>(I729_index);
  auto tensor309 = vector<shared_ptr<Tensor>>{I444, Gamma252_(), I729};
  auto task309 = make_shared<Task309>(tensor309, pindex);
  task300->add_dep(task309);
  task309->add_dep(task299);
  densityq->add_task(task309);

  auto tensor310 = vector<shared_ptr<Tensor>>{I729, t2};
  auto task310 = make_shared<Task310>(tensor310, pindex);
  task309->add_dep(task310);
  task310->add_dep(task299);
  densityq->add_task(task310);

  vector<IndexRange> I447_index = {closed_, closed_};
  auto I447 = make_shared<Tensor>(I447_index);
  auto tensor311 = vector<shared_ptr<Tensor>>{den2, I447};
  auto task311 = make_shared<Task311>(tensor311, pindex);
  task311->add_dep(task299);
  densityq->add_task(task311);

  vector<IndexRange> I448_index = {closed_, closed_, active_, active_};
  auto I448 = make_shared<Tensor>(I448_index);
  auto tensor312 = vector<shared_ptr<Tensor>>{I447, t2, I448};
  auto task312 = make_shared<Task312>(tensor312, pindex);
  task311->add_dep(task312);
  task312->add_dep(task299);
  densityq->add_task(task312);

  auto tensor313 = vector<shared_ptr<Tensor>>{I448, Gamma92_(), t2};
  auto task313 = make_shared<Task313>(tensor313, pindex);
  task312->add_dep(task313);
  task313->add_dep(task299);
  densityq->add_task(task313);

  vector<IndexRange> I541_index = {closed_, virt_, active_, active_};
  auto I541 = make_shared<Tensor>(I541_index);
  auto tensor314 = vector<shared_ptr<Tensor>>{I447, t2, I541};
  auto task314 = make_shared<Task314>(tensor314, pindex);
  task311->add_dep(task314);
  task314->add_dep(task299);
  densityq->add_task(task314);

  auto tensor315 = vector<shared_ptr<Tensor>>{I541, Gamma32_(), t2};
  auto task315 = make_shared<Task315>(tensor315, pindex);
  task314->add_dep(task315);
  task315->add_dep(task299);
  densityq->add_task(task315);

  vector<IndexRange> I550_index = {closed_, virt_, active_, active_};
  auto I550 = make_shared<Tensor>(I550_index);
  auto tensor316 = vector<shared_ptr<Tensor>>{I447, t2, I550};
  auto task316 = make_shared<Task316>(tensor316, pindex);
  task311->add_dep(task316);
  task316->add_dep(task299);
  densityq->add_task(task316);

  auto tensor317 = vector<shared_ptr<Tensor>>{I550, Gamma35_(), t2};
  auto task317 = make_shared<Task317>(tensor317, pindex);
  task316->add_dep(task317);
  task317->add_dep(task299);
  densityq->add_task(task317);

  vector<IndexRange> I450_index = {active_, closed_};
  auto I450 = make_shared<Tensor>(I450_index);
  auto tensor318 = vector<shared_ptr<Tensor>>{den2, I450};
  auto task318 = make_shared<Task318>(tensor318, pindex);
  task318->add_dep(task299);
  densityq->add_task(task318);

  vector<IndexRange> I451_index = {closed_, active_, active_, active_};
  auto I451 = make_shared<Tensor>(I451_index);
  auto tensor319 = vector<shared_ptr<Tensor>>{I450, t2, I451};
  auto task319 = make_shared<Task319>(tensor319, pindex);
  task318->add_dep(task319);
  task319->add_dep(task299);
  densityq->add_task(task319);

  auto tensor320 = vector<shared_ptr<Tensor>>{I451, Gamma2_(), t2};
  auto task320 = make_shared<Task320>(tensor320, pindex);
  task319->add_dep(task320);
  task320->add_dep(task299);
  densityq->add_task(task320);

  vector<IndexRange> I556_index = {virt_, active_, active_, active_};
  auto I556 = make_shared<Tensor>(I556_index);
  auto tensor321 = vector<shared_ptr<Tensor>>{I450, t2, I556};
  auto task321 = make_shared<Task321>(tensor321, pindex);
  task318->add_dep(task321);
  task321->add_dep(task299);
  densityq->add_task(task321);

  auto tensor322 = vector<shared_ptr<Tensor>>{I556, Gamma37_(), t2};
  auto task322 = make_shared<Task322>(tensor322, pindex);
  task321->add_dep(task322);
  task322->add_dep(task299);
  densityq->add_task(task322);

  vector<IndexRange> I453_index = {active_, virt_};
  auto I453 = make_shared<Tensor>(I453_index);
  auto tensor323 = vector<shared_ptr<Tensor>>{den2, I453};
  auto task323 = make_shared<Task323>(tensor323, pindex);
  task323->add_dep(task299);
  densityq->add_task(task323);

  vector<IndexRange> I454_index = {closed_, closed_, active_, active_};
  auto I454 = make_shared<Tensor>(I454_index);
  auto tensor324 = vector<shared_ptr<Tensor>>{I453, t2, I454};
  auto task324 = make_shared<Task324>(tensor324, pindex);
  task323->add_dep(task324);
  task324->add_dep(task299);
  densityq->add_task(task324);

  auto tensor325 = vector<shared_ptr<Tensor>>{I454, Gamma3_(), t2};
  auto task325 = make_shared<Task325>(tensor325, pindex);
  task324->add_dep(task325);
  task325->add_dep(task299);
  densityq->add_task(task325);

  vector<IndexRange> I565_index = {closed_, virt_, active_, active_};
  auto I565 = make_shared<Tensor>(I565_index);
  auto tensor326 = vector<shared_ptr<Tensor>>{I453, t2, I565};
  auto task326 = make_shared<Task326>(tensor326, pindex);
  task323->add_dep(task326);
  task326->add_dep(task299);
  densityq->add_task(task326);

  auto tensor327 = vector<shared_ptr<Tensor>>{I565, Gamma35_(), t2};
  auto task327 = make_shared<Task327>(tensor327, pindex);
  task326->add_dep(task327);
  task327->add_dep(task299);
  densityq->add_task(task327);

  vector<IndexRange> I568_index = {closed_, virt_, active_, active_};
  auto I568 = make_shared<Tensor>(I568_index);
  auto tensor328 = vector<shared_ptr<Tensor>>{I453, t2, I568};
  auto task328 = make_shared<Task328>(tensor328, pindex);
  task323->add_dep(task328);
  task328->add_dep(task299);
  densityq->add_task(task328);

  auto tensor329 = vector<shared_ptr<Tensor>>{I568, Gamma32_(), t2};
  auto task329 = make_shared<Task329>(tensor329, pindex);
  task328->add_dep(task329);
  task329->add_dep(task299);
  densityq->add_task(task329);

  vector<IndexRange> I607_index = {virt_, closed_, active_, active_};
  auto I607 = make_shared<Tensor>(I607_index);
  auto tensor330 = vector<shared_ptr<Tensor>>{I453, t2, I607};
  auto task330 = make_shared<Task330>(tensor330, pindex);
  task323->add_dep(task330);
  task330->add_dep(task299);
  densityq->add_task(task330);

  auto tensor331 = vector<shared_ptr<Tensor>>{I607, Gamma35_(), t2};
  auto task331 = make_shared<Task331>(tensor331, pindex);
  task330->add_dep(task331);
  task331->add_dep(task299);
  densityq->add_task(task331);

  vector<IndexRange> I610_index = {virt_, closed_, active_, active_};
  auto I610 = make_shared<Tensor>(I610_index);
  auto tensor332 = vector<shared_ptr<Tensor>>{I453, t2, I610};
  auto task332 = make_shared<Task332>(tensor332, pindex);
  task323->add_dep(task332);
  task332->add_dep(task299);
  densityq->add_task(task332);

  auto tensor333 = vector<shared_ptr<Tensor>>{I610, Gamma35_(), t2};
  auto task333 = make_shared<Task333>(tensor333, pindex);
  task332->add_dep(task333);
  task333->add_dep(task299);
  densityq->add_task(task333);

  vector<IndexRange> I456_index = {active_, closed_};
  auto I456 = make_shared<Tensor>(I456_index);
  auto tensor334 = vector<shared_ptr<Tensor>>{den2, I456};
  auto task334 = make_shared<Task334>(tensor334, pindex);
  task334->add_dep(task299);
  densityq->add_task(task334);

  vector<IndexRange> I457_index = {closed_, active_, active_, active_};
  auto I457 = make_shared<Tensor>(I457_index);
  auto tensor335 = vector<shared_ptr<Tensor>>{I456, t2, I457};
  auto task335 = make_shared<Task335>(tensor335, pindex);
  task334->add_dep(task335);
  task335->add_dep(task299);
  densityq->add_task(task335);

  auto tensor336 = vector<shared_ptr<Tensor>>{I457, Gamma4_(), t2};
  auto task336 = make_shared<Task336>(tensor336, pindex);
  task335->add_dep(task336);
  task336->add_dep(task299);
  densityq->add_task(task336);

  vector<IndexRange> I613_index = {virt_, active_, active_, active_};
  auto I613 = make_shared<Tensor>(I613_index);
  auto tensor337 = vector<shared_ptr<Tensor>>{I456, t2, I613};
  auto task337 = make_shared<Task337>(tensor337, pindex);
  task334->add_dep(task337);
  task337->add_dep(task299);
  densityq->add_task(task337);

  auto tensor338 = vector<shared_ptr<Tensor>>{I613, Gamma56_(), t2};
  auto task338 = make_shared<Task338>(tensor338, pindex);
  task337->add_dep(task338);
  task338->add_dep(task299);
  densityq->add_task(task338);

  vector<IndexRange> I616_index = {virt_, active_, active_, active_};
  auto I616 = make_shared<Tensor>(I616_index);
  auto tensor339 = vector<shared_ptr<Tensor>>{I456, t2, I616};
  auto task339 = make_shared<Task339>(tensor339, pindex);
  task334->add_dep(task339);
  task339->add_dep(task299);
  densityq->add_task(task339);

  auto tensor340 = vector<shared_ptr<Tensor>>{I616, Gamma57_(), t2};
  auto task340 = make_shared<Task340>(tensor340, pindex);
  task339->add_dep(task340);
  task340->add_dep(task299);
  densityq->add_task(task340);

  vector<IndexRange> I459_index = {active_, active_};
  auto I459 = make_shared<Tensor>(I459_index);
  auto tensor341 = vector<shared_ptr<Tensor>>{den2, I459};
  auto task341 = make_shared<Task341>(tensor341, pindex);
  task341->add_dep(task299);
  densityq->add_task(task341);

  vector<IndexRange> I460_index = {active_, active_, active_, active_, active_, active_};
  auto I460 = make_shared<Tensor>(I460_index);
  auto tensor342 = vector<shared_ptr<Tensor>>{I459, Gamma165_(), I460};
  auto task342 = make_shared<Task342>(tensor342, pindex);
  task341->add_dep(task342);
  task342->add_dep(task299);
  densityq->add_task(task342);

  auto tensor343 = vector<shared_ptr<Tensor>>{I460, t2};
  auto task343 = make_shared<Task343>(tensor343, pindex);
  task342->add_dep(task343);
  task343->add_dep(task299);
  densityq->add_task(task343);

  vector<IndexRange> I619_index = {active_, active_, active_, active_, active_, active_};
  auto I619 = make_shared<Tensor>(I619_index);
  auto tensor344 = vector<shared_ptr<Tensor>>{I459, Gamma218_(), I619};
  auto task344 = make_shared<Task344>(tensor344, pindex);
  task341->add_dep(task344);
  task344->add_dep(task299);
  densityq->add_task(task344);

  auto tensor345 = vector<shared_ptr<Tensor>>{I619, t2};
  auto task345 = make_shared<Task345>(tensor345, pindex);
  task344->add_dep(task345);
  task345->add_dep(task299);
  densityq->add_task(task345);

  vector<IndexRange> I462_index = {closed_, closed_};
  auto I462 = make_shared<Tensor>(I462_index);
  auto tensor346 = vector<shared_ptr<Tensor>>{den2, I462};
  auto task346 = make_shared<Task346>(tensor346, pindex);
  task346->add_dep(task299);
  densityq->add_task(task346);

  vector<IndexRange> I463_index = {closed_, active_, active_, active_};
  auto I463 = make_shared<Tensor>(I463_index);
  auto tensor347 = vector<shared_ptr<Tensor>>{I462, t2, I463};
  auto task347 = make_shared<Task347>(tensor347, pindex);
  task346->add_dep(task347);
  task347->add_dep(task299);
  densityq->add_task(task347);

  auto tensor348 = vector<shared_ptr<Tensor>>{I463, Gamma6_(), t2};
  auto task348 = make_shared<Task348>(tensor348, pindex);
  task347->add_dep(task348);
  task348->add_dep(task299);
  densityq->add_task(task348);

  vector<IndexRange> I465_index = {closed_, virt_};
  auto I465 = make_shared<Tensor>(I465_index);
  auto tensor349 = vector<shared_ptr<Tensor>>{den2, I465};
  auto task349 = make_shared<Task349>(tensor349, pindex);
  task349->add_dep(task299);
  densityq->add_task(task349);

  vector<IndexRange> I466_index = {closed_, active_};
  auto I466 = make_shared<Tensor>(I466_index);
  auto tensor350 = vector<shared_ptr<Tensor>>{I465, t2, I466};
  auto task350 = make_shared<Task350>(tensor350, pindex);
  task349->add_dep(task350);
  task350->add_dep(task299);
  densityq->add_task(task350);

  auto tensor351 = vector<shared_ptr<Tensor>>{I466, Gamma7_(), t2};
  auto task351 = make_shared<Task351>(tensor351, pindex);
  task350->add_dep(task351);
  task351->add_dep(task299);
  densityq->add_task(task351);

  vector<IndexRange> I469_index = {closed_, active_};
  auto I469 = make_shared<Tensor>(I469_index);
  auto tensor352 = vector<shared_ptr<Tensor>>{I465, t2, I469};
  auto task352 = make_shared<Task352>(tensor352, pindex);
  task349->add_dep(task352);
  task352->add_dep(task299);
  densityq->add_task(task352);

  auto tensor353 = vector<shared_ptr<Tensor>>{I469, Gamma7_(), t2};
  auto task353 = make_shared<Task353>(tensor353, pindex);
  task352->add_dep(task353);
  task353->add_dep(task299);
  densityq->add_task(task353);

  vector<IndexRange> I625_index = {virt_, active_};
  auto I625 = make_shared<Tensor>(I625_index);
  auto tensor354 = vector<shared_ptr<Tensor>>{I465, t2, I625};
  auto task354 = make_shared<Task354>(tensor354, pindex);
  task349->add_dep(task354);
  task354->add_dep(task299);
  densityq->add_task(task354);

  auto tensor355 = vector<shared_ptr<Tensor>>{I625, Gamma60_(), t2};
  auto task355 = make_shared<Task355>(tensor355, pindex);
  task354->add_dep(task355);
  task355->add_dep(task299);
  densityq->add_task(task355);

  vector<IndexRange> I628_index = {virt_, active_};
  auto I628 = make_shared<Tensor>(I628_index);
  auto tensor356 = vector<shared_ptr<Tensor>>{I465, t2, I628};
  auto task356 = make_shared<Task356>(tensor356, pindex);
  task349->add_dep(task356);
  task356->add_dep(task299);
  densityq->add_task(task356);

  auto tensor357 = vector<shared_ptr<Tensor>>{I628, Gamma60_(), t2};
  auto task357 = make_shared<Task357>(tensor357, pindex);
  task356->add_dep(task357);
  task357->add_dep(task299);
  densityq->add_task(task357);

  vector<IndexRange> I471_index = {active_, virt_};
  auto I471 = make_shared<Tensor>(I471_index);
  auto tensor358 = vector<shared_ptr<Tensor>>{den2, I471};
  auto task358 = make_shared<Task358>(tensor358, pindex);
  task358->add_dep(task299);
  densityq->add_task(task358);

  vector<IndexRange> I472_index = {closed_, active_, active_, active_};
  auto I472 = make_shared<Tensor>(I472_index);
  auto tensor359 = vector<shared_ptr<Tensor>>{I471, t2, I472};
  auto task359 = make_shared<Task359>(tensor359, pindex);
  task358->add_dep(task359);
  task359->add_dep(task299);
  densityq->add_task(task359);

  auto tensor360 = vector<shared_ptr<Tensor>>{I472, Gamma9_(), t2};
  auto task360 = make_shared<Task360>(tensor360, pindex);
  task359->add_dep(task360);
  task360->add_dep(task299);
  densityq->add_task(task360);

  vector<IndexRange> I475_index = {closed_, active_, active_, active_};
  auto I475 = make_shared<Tensor>(I475_index);
  auto tensor361 = vector<shared_ptr<Tensor>>{I471, t2, I475};
  auto task361 = make_shared<Task361>(tensor361, pindex);
  task358->add_dep(task361);
  task361->add_dep(task299);
  densityq->add_task(task361);

  auto tensor362 = vector<shared_ptr<Tensor>>{I475, Gamma6_(), t2};
  auto task362 = make_shared<Task362>(tensor362, pindex);
  task361->add_dep(task362);
  task362->add_dep(task299);
  densityq->add_task(task362);

  vector<IndexRange> I631_index = {virt_, active_, active_, active_};
  auto I631 = make_shared<Tensor>(I631_index);
  auto tensor363 = vector<shared_ptr<Tensor>>{I471, t2, I631};
  auto task363 = make_shared<Task363>(tensor363, pindex);
  task358->add_dep(task363);
  task363->add_dep(task299);
  densityq->add_task(task363);

  auto tensor364 = vector<shared_ptr<Tensor>>{I631, Gamma59_(), t2};
  auto task364 = make_shared<Task364>(tensor364, pindex);
  task363->add_dep(task364);
  task364->add_dep(task299);
  densityq->add_task(task364);

  vector<IndexRange> I477_index = {active_, virt_};
  auto I477 = make_shared<Tensor>(I477_index);
  auto tensor365 = vector<shared_ptr<Tensor>>{den2, I477};
  auto task365 = make_shared<Task365>(tensor365, pindex);
  task365->add_dep(task299);
  densityq->add_task(task365);

  vector<IndexRange> I478_index = {closed_, closed_, active_, active_};
  auto I478 = make_shared<Tensor>(I478_index);
  auto tensor366 = vector<shared_ptr<Tensor>>{I477, t2, I478};
  auto task366 = make_shared<Task366>(tensor366, pindex);
  task365->add_dep(task366);
  task366->add_dep(task299);
  densityq->add_task(task366);

  auto tensor367 = vector<shared_ptr<Tensor>>{I478, Gamma3_(), t2};
  auto task367 = make_shared<Task367>(tensor367, pindex);
  task366->add_dep(task367);
  task367->add_dep(task299);
  densityq->add_task(task367);

  vector<IndexRange> I480_index = {virt_, closed_};
  auto I480 = make_shared<Tensor>(I480_index);
  auto tensor368 = vector<shared_ptr<Tensor>>{den2, I480};
  auto task368 = make_shared<Task368>(tensor368, pindex);
  task368->add_dep(task299);
  densityq->add_task(task368);

  vector<IndexRange> I481_index = {closed_, active_};
  auto I481 = make_shared<Tensor>(I481_index);
  auto tensor369 = vector<shared_ptr<Tensor>>{I480, t2, I481};
  auto task369 = make_shared<Task369>(tensor369, pindex);
  task368->add_dep(task369);
  task369->add_dep(task299);
  densityq->add_task(task369);

  auto tensor370 = vector<shared_ptr<Tensor>>{I481, Gamma12_(), t2};
  auto task370 = make_shared<Task370>(tensor370, pindex);
  task369->add_dep(task370);
  task370->add_dep(task299);
  densityq->add_task(task370);

  vector<IndexRange> I483_index = {closed_, virt_};
  auto I483 = make_shared<Tensor>(I483_index);
  auto tensor371 = vector<shared_ptr<Tensor>>{den2, I483};
  auto task371 = make_shared<Task371>(tensor371, pindex);
  task371->add_dep(task299);
  densityq->add_task(task371);

  vector<IndexRange> I484_index = {closed_, active_};
  auto I484 = make_shared<Tensor>(I484_index);
  auto tensor372 = vector<shared_ptr<Tensor>>{I483, t2, I484};
  auto task372 = make_shared<Task372>(tensor372, pindex);
  task371->add_dep(task372);
  task372->add_dep(task299);
  densityq->add_task(task372);

  auto tensor373 = vector<shared_ptr<Tensor>>{I484, Gamma12_(), t2};
  auto task373 = make_shared<Task373>(tensor373, pindex);
  task372->add_dep(task373);
  task373->add_dep(task299);
  densityq->add_task(task373);

  vector<IndexRange> I640_index = {virt_, closed_};
  auto I640 = make_shared<Tensor>(I640_index);
  auto tensor374 = vector<shared_ptr<Tensor>>{I483, t2, I640};
  auto task374 = make_shared<Task374>(tensor374, pindex);
  task371->add_dep(task374);
  task374->add_dep(task299);
  densityq->add_task(task374);

  vector<IndexRange> I641_index = {active_, virt_, closed_, active_};
  auto I641 = make_shared<Tensor>(I641_index);
  auto tensor375 = vector<shared_ptr<Tensor>>{I640, Gamma38_(), I641};
  auto task375 = make_shared<Task375>(tensor375, pindex);
  task374->add_dep(task375);
  task375->add_dep(task299);
  densityq->add_task(task375);

  auto tensor376 = vector<shared_ptr<Tensor>>{I641, t2};
  auto task376 = make_shared<Task376>(tensor376, pindex);
  task375->add_dep(task376);
  task376->add_dep(task299);
  densityq->add_task(task376);

  vector<IndexRange> I486_index = {active_, active_};
  auto I486 = make_shared<Tensor>(I486_index);
  auto tensor377 = vector<shared_ptr<Tensor>>{den2, I486};
  auto task377 = make_shared<Task377>(tensor377, pindex);
  task377->add_dep(task299);
  densityq->add_task(task377);

  vector<IndexRange> I487_index = {active_, active_};
  auto I487 = make_shared<Tensor>(I487_index);
  auto tensor378 = vector<shared_ptr<Tensor>>{I486, Gamma174_(), I487};
  auto task378 = make_shared<Task378>(tensor378, pindex);
  task377->add_dep(task378);
  task378->add_dep(task299);
  densityq->add_task(task378);

  auto tensor379 = vector<shared_ptr<Tensor>>{I487, t2};
  auto task379 = make_shared<Task379>(tensor379, pindex);
  task378->add_dep(task379);
  task379->add_dep(task299);
  densityq->add_task(task379);

  auto tensor380 = vector<shared_ptr<Tensor>>{I487, t2};
  auto task380 = make_shared<Task380>(tensor380, pindex);
  task378->add_dep(task380);
  task380->add_dep(task299);
  densityq->add_task(task380);

  vector<IndexRange> I696_index = {active_, active_};
  auto I696 = make_shared<Tensor>(I696_index);
  auto tensor381 = vector<shared_ptr<Tensor>>{I486, Gamma60_(), I696};
  auto task381 = make_shared<Task381>(tensor381, pindex);
  task377->add_dep(task381);
  task381->add_dep(task299);
  densityq->add_task(task381);

  auto tensor382 = vector<shared_ptr<Tensor>>{I696, t2};
  auto task382 = make_shared<Task382>(tensor382, pindex);
  task381->add_dep(task382);
  task382->add_dep(task299);
  densityq->add_task(task382);

  auto tensor383 = vector<shared_ptr<Tensor>>{I696, t2};
  auto task383 = make_shared<Task383>(tensor383, pindex);
  task381->add_dep(task383);
  task383->add_dep(task299);
  densityq->add_task(task383);

  vector<IndexRange> I492_index = {closed_, closed_};
  auto I492 = make_shared<Tensor>(I492_index);
  auto tensor384 = vector<shared_ptr<Tensor>>{den2, I492};
  auto task384 = make_shared<Task384>(tensor384, pindex);
  task384->add_dep(task299);
  densityq->add_task(task384);

  vector<IndexRange> I493_index = {closed_, virt_, closed_, active_};
  auto I493 = make_shared<Tensor>(I493_index);
  auto tensor385 = vector<shared_ptr<Tensor>>{I492, t2, I493};
  auto task385 = make_shared<Task385>(tensor385, pindex);
  task384->add_dep(task385);
  task385->add_dep(task299);
  densityq->add_task(task385);

  auto tensor386 = vector<shared_ptr<Tensor>>{I493, Gamma16_(), t2};
  auto task386 = make_shared<Task386>(tensor386, pindex);
  task385->add_dep(task386);
  task386->add_dep(task299);
  densityq->add_task(task386);

  vector<IndexRange> I496_index = {closed_, virt_, closed_, active_};
  auto I496 = make_shared<Tensor>(I496_index);
  auto tensor387 = vector<shared_ptr<Tensor>>{I492, t2, I496};
  auto task387 = make_shared<Task387>(tensor387, pindex);
  task384->add_dep(task387);
  task387->add_dep(task299);
  densityq->add_task(task387);

  auto tensor388 = vector<shared_ptr<Tensor>>{I496, Gamma16_(), t2};
  auto task388 = make_shared<Task388>(tensor388, pindex);
  task387->add_dep(task388);
  task388->add_dep(task299);
  densityq->add_task(task388);

  vector<IndexRange> I498_index = {closed_, closed_};
  auto I498 = make_shared<Tensor>(I498_index);
  auto tensor389 = vector<shared_ptr<Tensor>>{den2, I498};
  auto task389 = make_shared<Task389>(tensor389, pindex);
  task389->add_dep(task299);
  densityq->add_task(task389);

  vector<IndexRange> I499_index = {closed_, virt_, closed_, active_};
  auto I499 = make_shared<Tensor>(I499_index);
  auto tensor390 = vector<shared_ptr<Tensor>>{I498, t2, I499};
  auto task390 = make_shared<Task390>(tensor390, pindex);
  task389->add_dep(task390);
  task390->add_dep(task299);
  densityq->add_task(task390);

  auto tensor391 = vector<shared_ptr<Tensor>>{I499, Gamma16_(), t2};
  auto task391 = make_shared<Task391>(tensor391, pindex);
  task390->add_dep(task391);
  task391->add_dep(task299);
  densityq->add_task(task391);

  vector<IndexRange> I505_index = {closed_, virt_, closed_, active_};
  auto I505 = make_shared<Tensor>(I505_index);
  auto tensor392 = vector<shared_ptr<Tensor>>{I498, t2, I505};
  auto task392 = make_shared<Task392>(tensor392, pindex);
  task389->add_dep(task392);
  task392->add_dep(task299);
  densityq->add_task(task392);

  auto tensor393 = vector<shared_ptr<Tensor>>{I505, Gamma16_(), t2};
  auto task393 = make_shared<Task393>(tensor393, pindex);
  task392->add_dep(task393);
  task393->add_dep(task299);
  densityq->add_task(task393);

  vector<IndexRange> I501_index = {virt_, virt_};
  auto I501 = make_shared<Tensor>(I501_index);
  auto tensor394 = vector<shared_ptr<Tensor>>{den2, I501};
  auto task394 = make_shared<Task394>(tensor394, pindex);
  task394->add_dep(task299);
  densityq->add_task(task394);

  vector<IndexRange> I502_index = {closed_, virt_, closed_, active_};
  auto I502 = make_shared<Tensor>(I502_index);
  auto tensor395 = vector<shared_ptr<Tensor>>{I501, t2, I502};
  auto task395 = make_shared<Task395>(tensor395, pindex);
  task394->add_dep(task395);
  task395->add_dep(task299);
  densityq->add_task(task395);

  auto tensor396 = vector<shared_ptr<Tensor>>{I502, Gamma16_(), t2};
  auto task396 = make_shared<Task396>(tensor396, pindex);
  task395->add_dep(task396);
  task396->add_dep(task299);
  densityq->add_task(task396);

  vector<IndexRange> I508_index = {closed_, virt_, closed_, active_};
  auto I508 = make_shared<Tensor>(I508_index);
  auto tensor397 = vector<shared_ptr<Tensor>>{I501, t2, I508};
  auto task397 = make_shared<Task397>(tensor397, pindex);
  task394->add_dep(task397);
  task397->add_dep(task299);
  densityq->add_task(task397);

  auto tensor398 = vector<shared_ptr<Tensor>>{I508, Gamma16_(), t2};
  auto task398 = make_shared<Task398>(tensor398, pindex);
  task397->add_dep(task398);
  task398->add_dep(task299);
  densityq->add_task(task398);

  vector<IndexRange> I510_index = {active_, closed_};
  auto I510 = make_shared<Tensor>(I510_index);
  auto tensor399 = vector<shared_ptr<Tensor>>{den2, I510};
  auto task399 = make_shared<Task399>(tensor399, pindex);
  task399->add_dep(task299);
  densityq->add_task(task399);

  vector<IndexRange> I511_index = {virt_, closed_, active_, active_};
  auto I511 = make_shared<Tensor>(I511_index);
  auto tensor400 = vector<shared_ptr<Tensor>>{I510, t2, I511};
  auto task400 = make_shared<Task400>(tensor400, pindex);
  task399->add_dep(task400);
  task400->add_dep(task299);
  densityq->add_task(task400);

  auto tensor401 = vector<shared_ptr<Tensor>>{I511, Gamma22_(), t2};
  auto task401 = make_shared<Task401>(tensor401, pindex);
  task400->add_dep(task401);
  task401->add_dep(task299);
  densityq->add_task(task401);

  auto tensor402 = vector<shared_ptr<Tensor>>{I511, Gamma12_(), t2};
  auto task402 = make_shared<Task402>(tensor402, pindex);
  task400->add_dep(task402);
  task402->add_dep(task299);
  densityq->add_task(task402);

  vector<IndexRange> I513_index = {active_, closed_};
  auto I513 = make_shared<Tensor>(I513_index);
  auto tensor403 = vector<shared_ptr<Tensor>>{den2, I513};
  auto task403 = make_shared<Task403>(tensor403, pindex);
  task403->add_dep(task299);
  densityq->add_task(task403);

  vector<IndexRange> I514_index = {virt_, closed_, active_, active_};
  auto I514 = make_shared<Tensor>(I514_index);
  auto tensor404 = vector<shared_ptr<Tensor>>{I513, t2, I514};
  auto task404 = make_shared<Task404>(tensor404, pindex);
  task403->add_dep(task404);
  task404->add_dep(task299);
  densityq->add_task(task404);

  vector<IndexRange> I515_index = {active_, virt_, closed_, active_};
  auto I515 = make_shared<Tensor>(I515_index);
  auto tensor405 = vector<shared_ptr<Tensor>>{I514, Gamma12_(), I515};
  auto task405 = make_shared<Task405>(tensor405, pindex);
  task404->add_dep(task405);
  task405->add_dep(task299);
  densityq->add_task(task405);

  auto tensor406 = vector<shared_ptr<Tensor>>{I515, t2};
  auto task406 = make_shared<Task406>(tensor406, pindex);
  task405->add_dep(task406);
  task406->add_dep(task299);
  densityq->add_task(task406);

  vector<IndexRange> I522_index = {virt_, active_};
  auto I522 = make_shared<Tensor>(I522_index);
  auto tensor407 = vector<shared_ptr<Tensor>>{den2, I522};
  auto task407 = make_shared<Task407>(tensor407, pindex);
  task407->add_dep(task299);
  densityq->add_task(task407);

  vector<IndexRange> I523_index = {active_, virt_};
  auto I523 = make_shared<Tensor>(I523_index);
  auto tensor408 = vector<shared_ptr<Tensor>>{I522, Gamma16_(), I523};
  auto task408 = make_shared<Task408>(tensor408, pindex);
  task407->add_dep(task408);
  task408->add_dep(task299);
  densityq->add_task(task408);

  auto tensor409 = vector<shared_ptr<Tensor>>{I523, t2};
  auto task409 = make_shared<Task409>(tensor409, pindex);
  task408->add_dep(task409);
  task409->add_dep(task299);
  densityq->add_task(task409);

  auto tensor410 = vector<shared_ptr<Tensor>>{I523, t2};
  auto task410 = make_shared<Task410>(tensor410, pindex);
  task408->add_dep(task410);
  task410->add_dep(task299);
  densityq->add_task(task410);

  vector<IndexRange> I528_index = {active_, virt_};
  auto I528 = make_shared<Tensor>(I528_index);
  auto tensor411 = vector<shared_ptr<Tensor>>{den2, I528};
  auto task411 = make_shared<Task411>(tensor411, pindex);
  task411->add_dep(task299);
  densityq->add_task(task411);

  vector<IndexRange> I529_index = {closed_, active_, active_, active_};
  auto I529 = make_shared<Tensor>(I529_index);
  auto tensor412 = vector<shared_ptr<Tensor>>{I528, t2, I529};
  auto task412 = make_shared<Task412>(tensor412, pindex);
  task411->add_dep(task412);
  task412->add_dep(task299);
  densityq->add_task(task412);

  auto tensor413 = vector<shared_ptr<Tensor>>{I529, Gamma28_(), t2};
  auto task413 = make_shared<Task413>(tensor413, pindex);
  task412->add_dep(task413);
  task413->add_dep(task299);
  densityq->add_task(task413);

  vector<IndexRange> I531_index = {active_, closed_};
  auto I531 = make_shared<Tensor>(I531_index);
  auto tensor414 = vector<shared_ptr<Tensor>>{den2, I531};
  auto task414 = make_shared<Task414>(tensor414, pindex);
  task414->add_dep(task299);
  densityq->add_task(task414);

  vector<IndexRange> I532_index = {closed_, virt_, active_, active_};
  auto I532 = make_shared<Tensor>(I532_index);
  auto tensor415 = vector<shared_ptr<Tensor>>{I531, t2, I532};
  auto task415 = make_shared<Task415>(tensor415, pindex);
  task414->add_dep(task415);
  task415->add_dep(task299);
  densityq->add_task(task415);

  auto tensor416 = vector<shared_ptr<Tensor>>{I532, Gamma29_(), t2};
  auto task416 = make_shared<Task416>(tensor416, pindex);
  task415->add_dep(task416);
  task416->add_dep(task299);
  densityq->add_task(task416);

  vector<IndexRange> I535_index = {closed_, virt_, active_, active_};
  auto I535 = make_shared<Tensor>(I535_index);
  auto tensor417 = vector<shared_ptr<Tensor>>{I531, t2, I535};
  auto task417 = make_shared<Task417>(tensor417, pindex);
  task414->add_dep(task417);
  task417->add_dep(task299);
  densityq->add_task(task417);

  auto tensor418 = vector<shared_ptr<Tensor>>{I535, Gamma7_(), t2};
  auto task418 = make_shared<Task418>(tensor418, pindex);
  task417->add_dep(task418);
  task418->add_dep(task299);
  densityq->add_task(task418);

  vector<IndexRange> I574_index = {virt_, closed_, active_, active_};
  auto I574 = make_shared<Tensor>(I574_index);
  auto tensor419 = vector<shared_ptr<Tensor>>{I531, t2, I574};
  auto task419 = make_shared<Task419>(tensor419, pindex);
  task414->add_dep(task419);
  task419->add_dep(task299);
  densityq->add_task(task419);

  auto tensor420 = vector<shared_ptr<Tensor>>{I574, Gamma7_(), t2};
  auto task420 = make_shared<Task420>(tensor420, pindex);
  task419->add_dep(task420);
  task420->add_dep(task299);
  densityq->add_task(task420);

  vector<IndexRange> I577_index = {virt_, closed_, active_, active_};
  auto I577 = make_shared<Tensor>(I577_index);
  auto tensor421 = vector<shared_ptr<Tensor>>{I531, t2, I577};
  auto task421 = make_shared<Task421>(tensor421, pindex);
  task414->add_dep(task421);
  task421->add_dep(task299);
  densityq->add_task(task421);

  auto tensor422 = vector<shared_ptr<Tensor>>{I577, Gamma7_(), t2};
  auto task422 = make_shared<Task422>(tensor422, pindex);
  task421->add_dep(task422);
  task422->add_dep(task299);
  densityq->add_task(task422);

  vector<IndexRange> I726_index = {virt_, virt_, active_, active_};
  auto I726 = make_shared<Tensor>(I726_index);
  auto tensor423 = vector<shared_ptr<Tensor>>{I531, t2, I726};
  auto task423 = make_shared<Task423>(tensor423, pindex);
  task414->add_dep(task423);
  task423->add_dep(task299);
  densityq->add_task(task423);

  auto tensor424 = vector<shared_ptr<Tensor>>{I726, Gamma60_(), t2};
  auto task424 = make_shared<Task424>(tensor424, pindex);
  task423->add_dep(task424);
  task424->add_dep(task299);
  densityq->add_task(task424);

  vector<IndexRange> I543_index = {virt_, virt_};
  auto I543 = make_shared<Tensor>(I543_index);
  auto tensor425 = vector<shared_ptr<Tensor>>{den2, I543};
  auto task425 = make_shared<Task425>(tensor425, pindex);
  task425->add_dep(task299);
  densityq->add_task(task425);

  vector<IndexRange> I544_index = {closed_, virt_, active_, active_};
  auto I544 = make_shared<Tensor>(I544_index);
  auto tensor426 = vector<shared_ptr<Tensor>>{I543, t2, I544};
  auto task426 = make_shared<Task426>(tensor426, pindex);
  task425->add_dep(task426);
  task426->add_dep(task299);
  densityq->add_task(task426);

  auto tensor427 = vector<shared_ptr<Tensor>>{I544, Gamma32_(), t2};
  auto task427 = make_shared<Task427>(tensor427, pindex);
  task426->add_dep(task427);
  task427->add_dep(task299);
  densityq->add_task(task427);

  vector<IndexRange> I553_index = {closed_, virt_, active_, active_};
  auto I553 = make_shared<Tensor>(I553_index);
  auto tensor428 = vector<shared_ptr<Tensor>>{I543, t2, I553};
  auto task428 = make_shared<Task428>(tensor428, pindex);
  task425->add_dep(task428);
  task428->add_dep(task299);
  densityq->add_task(task428);

  auto tensor429 = vector<shared_ptr<Tensor>>{I553, Gamma35_(), t2};
  auto task429 = make_shared<Task429>(tensor429, pindex);
  task428->add_dep(task429);
  task429->add_dep(task299);
  densityq->add_task(task429);

  vector<IndexRange> I558_index = {virt_, closed_};
  auto I558 = make_shared<Tensor>(I558_index);
  auto tensor430 = vector<shared_ptr<Tensor>>{den2, I558};
  auto task430 = make_shared<Task430>(tensor430, pindex);
  task430->add_dep(task299);
  densityq->add_task(task430);

  vector<IndexRange> I559_index = {closed_, virt_};
  auto I559 = make_shared<Tensor>(I559_index);
  auto tensor431 = vector<shared_ptr<Tensor>>{I558, t2, I559};
  auto task431 = make_shared<Task431>(tensor431, pindex);
  task430->add_dep(task431);
  task431->add_dep(task299);
  densityq->add_task(task431);

  auto tensor432 = vector<shared_ptr<Tensor>>{I559, Gamma38_(), t2};
  auto task432 = make_shared<Task432>(tensor432, pindex);
  task431->add_dep(task432);
  task432->add_dep(task299);
  densityq->add_task(task432);

  vector<IndexRange> I562_index = {closed_, virt_};
  auto I562 = make_shared<Tensor>(I562_index);
  auto tensor433 = vector<shared_ptr<Tensor>>{I558, t2, I562};
  auto task433 = make_shared<Task433>(tensor433, pindex);
  task430->add_dep(task433);
  task433->add_dep(task299);
  densityq->add_task(task433);

  auto tensor434 = vector<shared_ptr<Tensor>>{I562, Gamma38_(), t2};
  auto task434 = make_shared<Task434>(tensor434, pindex);
  task433->add_dep(task434);
  task434->add_dep(task299);
  densityq->add_task(task434);

  vector<IndexRange> I601_index = {virt_, closed_};
  auto I601 = make_shared<Tensor>(I601_index);
  auto tensor435 = vector<shared_ptr<Tensor>>{I558, t2, I601};
  auto task435 = make_shared<Task435>(tensor435, pindex);
  task430->add_dep(task435);
  task435->add_dep(task299);
  densityq->add_task(task435);

  auto tensor436 = vector<shared_ptr<Tensor>>{I601, Gamma38_(), t2};
  auto task436 = make_shared<Task436>(tensor436, pindex);
  task435->add_dep(task436);
  task436->add_dep(task299);
  densityq->add_task(task436);

  vector<IndexRange> I604_index = {virt_, closed_};
  auto I604 = make_shared<Tensor>(I604_index);
  auto tensor437 = vector<shared_ptr<Tensor>>{I558, t2, I604};
  auto task437 = make_shared<Task437>(tensor437, pindex);
  task430->add_dep(task437);
  task437->add_dep(task299);
  densityq->add_task(task437);

  auto tensor438 = vector<shared_ptr<Tensor>>{I604, Gamma38_(), t2};
  auto task438 = make_shared<Task438>(tensor438, pindex);
  task437->add_dep(task438);
  task438->add_dep(task299);
  densityq->add_task(task438);

  vector<IndexRange> I570_index = {active_, virt_};
  auto I570 = make_shared<Tensor>(I570_index);
  auto tensor439 = vector<shared_ptr<Tensor>>{den2, I570};
  auto task439 = make_shared<Task439>(tensor439, pindex);
  task439->add_dep(task299);
  densityq->add_task(task439);

  vector<IndexRange> I571_index = {closed_, active_, active_, active_};
  auto I571 = make_shared<Tensor>(I571_index);
  auto tensor440 = vector<shared_ptr<Tensor>>{I570, t2, I571};
  auto task440 = make_shared<Task440>(tensor440, pindex);
  task439->add_dep(task440);
  task440->add_dep(task299);
  densityq->add_task(task440);

  auto tensor441 = vector<shared_ptr<Tensor>>{I571, Gamma6_(), t2};
  auto task441 = make_shared<Task441>(tensor441, pindex);
  task440->add_dep(task441);
  task441->add_dep(task299);
  densityq->add_task(task441);

  vector<IndexRange> I723_index = {virt_, active_, active_, active_};
  auto I723 = make_shared<Tensor>(I723_index);
  auto tensor442 = vector<shared_ptr<Tensor>>{I570, t2, I723};
  auto task442 = make_shared<Task442>(tensor442, pindex);
  task439->add_dep(task442);
  task442->add_dep(task299);
  densityq->add_task(task442);

  auto tensor443 = vector<shared_ptr<Tensor>>{I723, Gamma59_(), t2};
  auto task443 = make_shared<Task443>(tensor443, pindex);
  task442->add_dep(task443);
  task443->add_dep(task299);
  densityq->add_task(task443);

  vector<IndexRange> I582_index = {closed_, closed_};
  auto I582 = make_shared<Tensor>(I582_index);
  auto tensor444 = vector<shared_ptr<Tensor>>{den2, I582};
  auto task444 = make_shared<Task444>(tensor444, pindex);
  task444->add_dep(task299);
  densityq->add_task(task444);

  vector<IndexRange> I583_index = {virt_, closed_, active_, active_};
  auto I583 = make_shared<Tensor>(I583_index);
  auto tensor445 = vector<shared_ptr<Tensor>>{I582, t2, I583};
  auto task445 = make_shared<Task445>(tensor445, pindex);
  task444->add_dep(task445);
  task445->add_dep(task299);
  densityq->add_task(task445);

  auto tensor446 = vector<shared_ptr<Tensor>>{I583, Gamma35_(), t2};
  auto task446 = make_shared<Task446>(tensor446, pindex);
  task445->add_dep(task446);
  task446->add_dep(task299);
  densityq->add_task(task446);

  vector<IndexRange> I592_index = {virt_, closed_, active_, active_};
  auto I592 = make_shared<Tensor>(I592_index);
  auto tensor447 = vector<shared_ptr<Tensor>>{I582, t2, I592};
  auto task447 = make_shared<Task447>(tensor447, pindex);
  task444->add_dep(task447);
  task447->add_dep(task299);
  densityq->add_task(task447);

  auto tensor448 = vector<shared_ptr<Tensor>>{I592, Gamma35_(), t2};
  auto task448 = make_shared<Task448>(tensor448, pindex);
  task447->add_dep(task448);
  task448->add_dep(task299);
  densityq->add_task(task448);

  vector<IndexRange> I585_index = {virt_, virt_};
  auto I585 = make_shared<Tensor>(I585_index);
  auto tensor449 = vector<shared_ptr<Tensor>>{den2, I585};
  auto task449 = make_shared<Task449>(tensor449, pindex);
  task449->add_dep(task299);
  densityq->add_task(task449);

  vector<IndexRange> I586_index = {virt_, closed_, active_, active_};
  auto I586 = make_shared<Tensor>(I586_index);
  auto tensor450 = vector<shared_ptr<Tensor>>{I585, t2, I586};
  auto task450 = make_shared<Task450>(tensor450, pindex);
  task449->add_dep(task450);
  task450->add_dep(task299);
  densityq->add_task(task450);

  auto tensor451 = vector<shared_ptr<Tensor>>{I586, Gamma35_(), t2};
  auto task451 = make_shared<Task451>(tensor451, pindex);
  task450->add_dep(task451);
  task451->add_dep(task299);
  densityq->add_task(task451);

  vector<IndexRange> I595_index = {virt_, closed_, active_, active_};
  auto I595 = make_shared<Tensor>(I595_index);
  auto tensor452 = vector<shared_ptr<Tensor>>{I585, t2, I595};
  auto task452 = make_shared<Task452>(tensor452, pindex);
  task449->add_dep(task452);
  task452->add_dep(task299);
  densityq->add_task(task452);

  auto tensor453 = vector<shared_ptr<Tensor>>{I595, Gamma35_(), t2};
  auto task453 = make_shared<Task453>(tensor453, pindex);
  task452->add_dep(task453);
  task453->add_dep(task299);
  densityq->add_task(task453);

  vector<IndexRange> I732_index = {virt_, virt_, active_, active_};
  auto I732 = make_shared<Tensor>(I732_index);
  auto tensor454 = vector<shared_ptr<Tensor>>{I585, t2, I732};
  auto task454 = make_shared<Task454>(tensor454, pindex);
  task449->add_dep(task454);
  task454->add_dep(task299);
  densityq->add_task(task454);

  auto tensor455 = vector<shared_ptr<Tensor>>{I732, Gamma60_(), t2};
  auto task455 = make_shared<Task455>(tensor455, pindex);
  task454->add_dep(task455);
  task455->add_dep(task299);
  densityq->add_task(task455);

  vector<IndexRange> I597_index = {active_, closed_};
  auto I597 = make_shared<Tensor>(I597_index);
  auto tensor456 = vector<shared_ptr<Tensor>>{den2, I597};
  auto task456 = make_shared<Task456>(tensor456, pindex);
  task456->add_dep(task299);
  densityq->add_task(task456);

  vector<IndexRange> I598_index = {virt_, active_, active_, active_};
  auto I598 = make_shared<Tensor>(I598_index);
  auto tensor457 = vector<shared_ptr<Tensor>>{I597, t2, I598};
  auto task457 = make_shared<Task457>(tensor457, pindex);
  task456->add_dep(task457);
  task457->add_dep(task299);
  densityq->add_task(task457);

  auto tensor458 = vector<shared_ptr<Tensor>>{I598, Gamma51_(), t2};
  auto task458 = make_shared<Task458>(tensor458, pindex);
  task457->add_dep(task458);
  task458->add_dep(task299);
  densityq->add_task(task458);

  vector<IndexRange> I621_index = {virt_, virt_};
  auto I621 = make_shared<Tensor>(I621_index);
  auto tensor459 = vector<shared_ptr<Tensor>>{den2, I621};
  auto task459 = make_shared<Task459>(tensor459, pindex);
  task459->add_dep(task299);
  densityq->add_task(task459);

  vector<IndexRange> I622_index = {virt_, active_, active_, active_};
  auto I622 = make_shared<Tensor>(I622_index);
  auto tensor460 = vector<shared_ptr<Tensor>>{I621, t2, I622};
  auto task460 = make_shared<Task460>(tensor460, pindex);
  task459->add_dep(task460);
  task460->add_dep(task299);
  densityq->add_task(task460);

  auto tensor461 = vector<shared_ptr<Tensor>>{I622, Gamma59_(), t2};
  auto task461 = make_shared<Task461>(tensor461, pindex);
  task460->add_dep(task461);
  task461->add_dep(task299);
  densityq->add_task(task461);

  vector<IndexRange> I633_index = {virt_, active_};
  auto I633 = make_shared<Tensor>(I633_index);
  auto tensor462 = vector<shared_ptr<Tensor>>{den2, I633};
  auto task462 = make_shared<Task462>(tensor462, pindex);
  task462->add_dep(task299);
  densityq->add_task(task462);

  vector<IndexRange> I634_index = {virt_, active_};
  auto I634 = make_shared<Tensor>(I634_index);
  auto tensor463 = vector<shared_ptr<Tensor>>{I633, Gamma16_(), I634};
  auto task463 = make_shared<Task463>(tensor463, pindex);
  task462->add_dep(task463);
  task463->add_dep(task299);
  densityq->add_task(task463);

  auto tensor464 = vector<shared_ptr<Tensor>>{I634, t2};
  auto task464 = make_shared<Task464>(tensor464, pindex);
  task463->add_dep(task464);
  task464->add_dep(task299);
  densityq->add_task(task464);

  vector<IndexRange> I636_index = {virt_, active_};
  auto I636 = make_shared<Tensor>(I636_index);
  auto tensor465 = vector<shared_ptr<Tensor>>{den2, I636};
  auto task465 = make_shared<Task465>(tensor465, pindex);
  task465->add_dep(task299);
  densityq->add_task(task465);

  vector<IndexRange> I637_index = {virt_, active_};
  auto I637 = make_shared<Tensor>(I637_index);
  auto tensor466 = vector<shared_ptr<Tensor>>{I636, Gamma16_(), I637};
  auto task466 = make_shared<Task466>(tensor466, pindex);
  task465->add_dep(task466);
  task466->add_dep(task299);
  densityq->add_task(task466);

  auto tensor467 = vector<shared_ptr<Tensor>>{I637, t2};
  auto task467 = make_shared<Task467>(tensor467, pindex);
  task466->add_dep(task467);
  task467->add_dep(task299);
  densityq->add_task(task467);

  vector<IndexRange> I642_index = {virt_, closed_};
  auto I642 = make_shared<Tensor>(I642_index);
  auto tensor468 = vector<shared_ptr<Tensor>>{den2, I642};
  auto task468 = make_shared<Task468>(tensor468, pindex);
  task468->add_dep(task299);
  densityq->add_task(task468);

  vector<IndexRange> I643_index = {virt_, closed_};
  auto I643 = make_shared<Tensor>(I643_index);
  auto tensor469 = vector<shared_ptr<Tensor>>{I642, t2, I643};
  auto task469 = make_shared<Task469>(tensor469, pindex);
  task468->add_dep(task469);
  task469->add_dep(task299);
  densityq->add_task(task469);

  vector<IndexRange> I644_index = {active_, virt_, closed_, active_};
  auto I644 = make_shared<Tensor>(I644_index);
  auto tensor470 = vector<shared_ptr<Tensor>>{I643, Gamma38_(), I644};
  auto task470 = make_shared<Task470>(tensor470, pindex);
  task469->add_dep(task470);
  task470->add_dep(task299);
  densityq->add_task(task470);

  auto tensor471 = vector<shared_ptr<Tensor>>{I644, t2};
  auto task471 = make_shared<Task471>(tensor471, pindex);
  task470->add_dep(task471);
  task471->add_dep(task299);
  densityq->add_task(task471);

  vector<IndexRange> I651_index = {active_, active_};
  auto I651 = make_shared<Tensor>(I651_index);
  auto tensor472 = vector<shared_ptr<Tensor>>{den2, I651};
  auto task472 = make_shared<Task472>(tensor472, pindex);
  task472->add_dep(task299);
  densityq->add_task(task472);

  vector<IndexRange> I652_index;
  auto I652 = make_shared<Tensor>(I652_index);
  auto tensor473 = vector<shared_ptr<Tensor>>{I651, Gamma38_(), I652};
  auto task473 = make_shared<Task473>(tensor473, pindex);
  task472->add_dep(task473);
  task473->add_dep(task299);
  densityq->add_task(task473);

  auto tensor474 = vector<shared_ptr<Tensor>>{I652, t2};
  auto task474 = make_shared<Task474>(tensor474, pindex);
  task473->add_dep(task474);
  task474->add_dep(task299);
  densityq->add_task(task474);

  auto tensor475 = vector<shared_ptr<Tensor>>{I652, t2};
  auto task475 = make_shared<Task475>(tensor475, pindex);
  task473->add_dep(task475);
  task475->add_dep(task299);
  densityq->add_task(task475);

  shared_ptr<Tensor> I657;
  if (diagonal) {
    vector<IndexRange> I657_index = {closed_, closed_};
    I657 = make_shared<Tensor>(I657_index);
  }
  shared_ptr<Task476> task476;
  if (diagonal) {
    auto tensor476 = vector<shared_ptr<Tensor>>{den2, I657};
    task476 = make_shared<Task476>(tensor476, pindex);
    task476->add_dep(task299);
    densityq->add_task(task476);
  }

  shared_ptr<Task477> task477;
  if (diagonal) {
    auto tensor477 = vector<shared_ptr<Tensor>>{I657, t2};
    task477 = make_shared<Task477>(tensor477, pindex);
    task476->add_dep(task477);
    task477->add_dep(task299);
    densityq->add_task(task477);
  }

  shared_ptr<Task478> task478;
  if (diagonal) {
    auto tensor478 = vector<shared_ptr<Tensor>>{I657, t2};
    task478 = make_shared<Task478>(tensor478, pindex);
    task476->add_dep(task478);
    task478->add_dep(task299);
    densityq->add_task(task478);
  }

  shared_ptr<Tensor> I661;
  if (diagonal) {
    vector<IndexRange> I661_index = {virt_, virt_};
    I661 = make_shared<Tensor>(I661_index);
  }
  shared_ptr<Task479> task479;
  if (diagonal) {
    auto tensor479 = vector<shared_ptr<Tensor>>{den2, I661};
    task479 = make_shared<Task479>(tensor479, pindex);
    task479->add_dep(task299);
    densityq->add_task(task479);
  }

  shared_ptr<Task480> task480;
  if (diagonal) {
    auto tensor480 = vector<shared_ptr<Tensor>>{I661, t2};
    task480 = make_shared<Task480>(tensor480, pindex);
    task479->add_dep(task480);
    task480->add_dep(task299);
    densityq->add_task(task480);
  }

  shared_ptr<Task481> task481;
  if (diagonal) {
    auto tensor481 = vector<shared_ptr<Tensor>>{I661, t2};
    task481 = make_shared<Task481>(tensor481, pindex);
    task479->add_dep(task481);
    task481->add_dep(task299);
    densityq->add_task(task481);
  }

  vector<IndexRange> I665_index = {closed_, active_};
  auto I665 = make_shared<Tensor>(I665_index);
  auto tensor482 = vector<shared_ptr<Tensor>>{den2, I665};
  auto task482 = make_shared<Task482>(tensor482, pindex);
  task482->add_dep(task299);
  densityq->add_task(task482);

  vector<IndexRange> I666_index = {closed_, active_};
  auto I666 = make_shared<Tensor>(I666_index);
  auto tensor483 = vector<shared_ptr<Tensor>>{I665, Gamma38_(), I666};
  auto task483 = make_shared<Task483>(tensor483, pindex);
  task482->add_dep(task483);
  task483->add_dep(task299);
  densityq->add_task(task483);

  auto tensor484 = vector<shared_ptr<Tensor>>{I666, t2};
  auto task484 = make_shared<Task484>(tensor484, pindex);
  task483->add_dep(task484);
  task484->add_dep(task299);
  densityq->add_task(task484);

  auto tensor485 = vector<shared_ptr<Tensor>>{I666, t2};
  auto task485 = make_shared<Task485>(tensor485, pindex);
  task483->add_dep(task485);
  task485->add_dep(task299);
  densityq->add_task(task485);

  vector<IndexRange> I671_index = {active_, virt_};
  auto I671 = make_shared<Tensor>(I671_index);
  auto tensor486 = vector<shared_ptr<Tensor>>{den2, I671};
  auto task486 = make_shared<Task486>(tensor486, pindex);
  task486->add_dep(task299);
  densityq->add_task(task486);

  vector<IndexRange> I672_index = {virt_, closed_, active_, active_};
  auto I672 = make_shared<Tensor>(I672_index);
  auto tensor487 = vector<shared_ptr<Tensor>>{I671, t2, I672};
  auto task487 = make_shared<Task487>(tensor487, pindex);
  task486->add_dep(task487);
  task487->add_dep(task299);
  densityq->add_task(task487);

  vector<IndexRange> I673_index = {active_, virt_, closed_, active_};
  auto I673 = make_shared<Tensor>(I673_index);
  auto tensor488 = vector<shared_ptr<Tensor>>{I672, Gamma35_(), I673};
  auto task488 = make_shared<Task488>(tensor488, pindex);
  task487->add_dep(task488);
  task488->add_dep(task299);
  densityq->add_task(task488);

  auto tensor489 = vector<shared_ptr<Tensor>>{I673, t2};
  auto task489 = make_shared<Task489>(tensor489, pindex);
  task488->add_dep(task489);
  task489->add_dep(task299);
  densityq->add_task(task489);

  vector<IndexRange> I674_index = {active_, virt_};
  auto I674 = make_shared<Tensor>(I674_index);
  auto tensor490 = vector<shared_ptr<Tensor>>{den2, I674};
  auto task490 = make_shared<Task490>(tensor490, pindex);
  task490->add_dep(task299);
  densityq->add_task(task490);

  vector<IndexRange> I675_index = {virt_, closed_, active_, active_};
  auto I675 = make_shared<Tensor>(I675_index);
  auto tensor491 = vector<shared_ptr<Tensor>>{I674, t2, I675};
  auto task491 = make_shared<Task491>(tensor491, pindex);
  task490->add_dep(task491);
  task491->add_dep(task299);
  densityq->add_task(task491);

  auto tensor492 = vector<shared_ptr<Tensor>>{I675, Gamma32_(), t2};
  auto task492 = make_shared<Task492>(tensor492, pindex);
  task491->add_dep(task492);
  task492->add_dep(task299);
  densityq->add_task(task492);

  auto tensor493 = vector<shared_ptr<Tensor>>{I675, Gamma35_(), t2};
  auto task493 = make_shared<Task493>(tensor493, pindex);
  task491->add_dep(task493);
  task493->add_dep(task299);
  densityq->add_task(task493);

  vector<IndexRange> I683_index = {closed_, virt_};
  auto I683 = make_shared<Tensor>(I683_index);
  auto tensor494 = vector<shared_ptr<Tensor>>{den2, I683};
  auto task494 = make_shared<Task494>(tensor494, pindex);
  task494->add_dep(task299);
  densityq->add_task(task494);

  vector<IndexRange> I684_index = {virt_, active_};
  auto I684 = make_shared<Tensor>(I684_index);
  auto tensor495 = vector<shared_ptr<Tensor>>{I683, t2, I684};
  auto task495 = make_shared<Task495>(tensor495, pindex);
  task494->add_dep(task495);
  task495->add_dep(task299);
  densityq->add_task(task495);

  auto tensor496 = vector<shared_ptr<Tensor>>{I684, Gamma60_(), t2};
  auto task496 = make_shared<Task496>(tensor496, pindex);
  task495->add_dep(task496);
  task496->add_dep(task299);
  densityq->add_task(task496);

  vector<IndexRange> I686_index = {virt_, closed_};
  auto I686 = make_shared<Tensor>(I686_index);
  auto tensor497 = vector<shared_ptr<Tensor>>{den2, I686};
  auto task497 = make_shared<Task497>(tensor497, pindex);
  task497->add_dep(task299);
  densityq->add_task(task497);

  vector<IndexRange> I687_index = {virt_, active_};
  auto I687 = make_shared<Tensor>(I687_index);
  auto tensor498 = vector<shared_ptr<Tensor>>{I686, t2, I687};
  auto task498 = make_shared<Task498>(tensor498, pindex);
  task497->add_dep(task498);
  task498->add_dep(task299);
  densityq->add_task(task498);

  auto tensor499 = vector<shared_ptr<Tensor>>{I687, Gamma60_(), t2};
  auto task499 = make_shared<Task499>(tensor499, pindex);
  task498->add_dep(task499);
  task499->add_dep(task299);
  densityq->add_task(task499);

  vector<IndexRange> I689_index = {closed_, active_};
  auto I689 = make_shared<Tensor>(I689_index);
  auto tensor500 = vector<shared_ptr<Tensor>>{den2, I689};
  auto task500 = make_shared<Task500>(tensor500, pindex);
  task500->add_dep(task299);
  densityq->add_task(task500);

  vector<IndexRange> I690_index = {active_, closed_};
  auto I690 = make_shared<Tensor>(I690_index);
  auto tensor501 = vector<shared_ptr<Tensor>>{I689, Gamma38_(), I690};
  auto task501 = make_shared<Task501>(tensor501, pindex);
  task500->add_dep(task501);
  task501->add_dep(task299);
  densityq->add_task(task501);

  auto tensor502 = vector<shared_ptr<Tensor>>{I690, t2};
  auto task502 = make_shared<Task502>(tensor502, pindex);
  task501->add_dep(task502);
  task502->add_dep(task299);
  densityq->add_task(task502);

  auto tensor503 = vector<shared_ptr<Tensor>>{I690, t2};
  auto task503 = make_shared<Task503>(tensor503, pindex);
  task501->add_dep(task503);
  task503->add_dep(task299);
  densityq->add_task(task503);

  vector<IndexRange> I701_index = {closed_, closed_};
  auto I701 = make_shared<Tensor>(I701_index);
  auto tensor504 = vector<shared_ptr<Tensor>>{den2, I701};
  auto task504 = make_shared<Task504>(tensor504, pindex);
  task504->add_dep(task299);
  densityq->add_task(task504);

  vector<IndexRange> I702_index = {virt_, closed_, virt_, active_};
  auto I702 = make_shared<Tensor>(I702_index);
  auto tensor505 = vector<shared_ptr<Tensor>>{I701, t2, I702};
  auto task505 = make_shared<Task505>(tensor505, pindex);
  task504->add_dep(task505);
  task505->add_dep(task299);
  densityq->add_task(task505);

  auto tensor506 = vector<shared_ptr<Tensor>>{I702, Gamma38_(), t2};
  auto task506 = make_shared<Task506>(tensor506, pindex);
  task505->add_dep(task506);
  task506->add_dep(task299);
  densityq->add_task(task506);

  vector<IndexRange> I705_index = {virt_, closed_, virt_, active_};
  auto I705 = make_shared<Tensor>(I705_index);
  auto tensor507 = vector<shared_ptr<Tensor>>{I701, t2, I705};
  auto task507 = make_shared<Task507>(tensor507, pindex);
  task504->add_dep(task507);
  task507->add_dep(task299);
  densityq->add_task(task507);

  auto tensor508 = vector<shared_ptr<Tensor>>{I705, Gamma38_(), t2};
  auto task508 = make_shared<Task508>(tensor508, pindex);
  task507->add_dep(task508);
  task508->add_dep(task299);
  densityq->add_task(task508);

  vector<IndexRange> I707_index = {virt_, virt_};
  auto I707 = make_shared<Tensor>(I707_index);
  auto tensor509 = vector<shared_ptr<Tensor>>{den2, I707};
  auto task509 = make_shared<Task509>(tensor509, pindex);
  task509->add_dep(task299);
  densityq->add_task(task509);

  vector<IndexRange> I708_index = {virt_, closed_, virt_, active_};
  auto I708 = make_shared<Tensor>(I708_index);
  auto tensor510 = vector<shared_ptr<Tensor>>{I707, t2, I708};
  auto task510 = make_shared<Task510>(tensor510, pindex);
  task509->add_dep(task510);
  task510->add_dep(task299);
  densityq->add_task(task510);

  auto tensor511 = vector<shared_ptr<Tensor>>{I708, Gamma38_(), t2};
  auto task511 = make_shared<Task511>(tensor511, pindex);
  task510->add_dep(task511);
  task511->add_dep(task299);
  densityq->add_task(task511);

  vector<IndexRange> I711_index = {virt_, closed_, virt_, active_};
  auto I711 = make_shared<Tensor>(I711_index);
  auto tensor512 = vector<shared_ptr<Tensor>>{I707, t2, I711};
  auto task512 = make_shared<Task512>(tensor512, pindex);
  task509->add_dep(task512);
  task512->add_dep(task299);
  densityq->add_task(task512);

  auto tensor513 = vector<shared_ptr<Tensor>>{I711, Gamma38_(), t2};
  auto task513 = make_shared<Task513>(tensor513, pindex);
  task512->add_dep(task513);
  task513->add_dep(task299);
  densityq->add_task(task513);

  vector<IndexRange> I713_index = {virt_, virt_};
  auto I713 = make_shared<Tensor>(I713_index);
  auto tensor514 = vector<shared_ptr<Tensor>>{den2, I713};
  auto task514 = make_shared<Task514>(tensor514, pindex);
  task514->add_dep(task299);
  densityq->add_task(task514);

  vector<IndexRange> I714_index = {virt_, closed_, virt_, active_};
  auto I714 = make_shared<Tensor>(I714_index);
  auto tensor515 = vector<shared_ptr<Tensor>>{I713, t2, I714};
  auto task515 = make_shared<Task515>(tensor515, pindex);
  task514->add_dep(task515);
  task515->add_dep(task299);
  densityq->add_task(task515);

  auto tensor516 = vector<shared_ptr<Tensor>>{I714, Gamma38_(), t2};
  auto task516 = make_shared<Task516>(tensor516, pindex);
  task515->add_dep(task516);
  task516->add_dep(task299);
  densityq->add_task(task516);

  vector<IndexRange> I717_index = {virt_, closed_, virt_, active_};
  auto I717 = make_shared<Tensor>(I717_index);
  auto tensor517 = vector<shared_ptr<Tensor>>{I713, t2, I717};
  auto task517 = make_shared<Task517>(tensor517, pindex);
  task514->add_dep(task517);
  task517->add_dep(task299);
  densityq->add_task(task517);

  auto tensor518 = vector<shared_ptr<Tensor>>{I717, Gamma38_(), t2};
  auto task518 = make_shared<Task518>(tensor518, pindex);
  task517->add_dep(task518);
  task518->add_dep(task299);
  densityq->add_task(task518);

  vector<IndexRange> I719_index = {active_, closed_};
  auto I719 = make_shared<Tensor>(I719_index);
  auto tensor519 = vector<shared_ptr<Tensor>>{den2, I719};
  auto task519 = make_shared<Task519>(tensor519, pindex);
  task519->add_dep(task299);
  densityq->add_task(task519);

  vector<IndexRange> I720_index = {virt_, virt_, active_, active_};
  auto I720 = make_shared<Tensor>(I720_index);
  auto tensor520 = vector<shared_ptr<Tensor>>{I719, t2, I720};
  auto task520 = make_shared<Task520>(tensor520, pindex);
  task519->add_dep(task520);
  task520->add_dep(task299);
  densityq->add_task(task520);

  auto tensor521 = vector<shared_ptr<Tensor>>{I720, Gamma60_(), t2};
  auto task521 = make_shared<Task521>(tensor521, pindex);
  task520->add_dep(task521);
  task521->add_dep(task299);
  densityq->add_task(task521);

  return densityq;
}


#endif
