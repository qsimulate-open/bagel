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

  auto densityq = make_shared<Queue>();
  auto task300 = make_shared<Task300>(den2, reset);
  densityq->add_task(task300);

  auto I444 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task301 = make_shared<Task301>(den2, I444);
  task301->add_dep(task300);
  densityq->add_task(task301);

  auto I445 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task302 = make_shared<Task302>(I444, Gamma160_(), I445);
  task301->add_dep(task302);
  task302->add_dep(task300);
  densityq->add_task(task302);

  auto task303 = make_shared<Task303>(I445, t2);
  task302->add_dep(task303);
  task303->add_dep(task300);
  densityq->add_task(task303);

  auto I538 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task304 = make_shared<Task304>(I444, Gamma191_(), I538);
  task301->add_dep(task304);
  task304->add_dep(task300);
  densityq->add_task(task304);

  auto task305 = make_shared<Task305>(I538, t2);
  task304->add_dep(task305);
  task305->add_dep(task300);
  densityq->add_task(task305);

  auto I547 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task306 = make_shared<Task306>(I444, Gamma194_(), I547);
  task301->add_dep(task306);
  task306->add_dep(task300);
  densityq->add_task(task306);

  auto task307 = make_shared<Task307>(I547, t2);
  task306->add_dep(task307);
  task307->add_dep(task300);
  densityq->add_task(task307);

  auto task308 = make_shared<Task308>(I547, t2);
  task306->add_dep(task308);
  task308->add_dep(task300);
  densityq->add_task(task308);

  auto task309 = make_shared<Task309>(I547, t2);
  task306->add_dep(task309);
  task309->add_dep(task300);
  densityq->add_task(task309);

  auto I729 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task310 = make_shared<Task310>(I444, Gamma252_(), I729);
  task301->add_dep(task310);
  task310->add_dep(task300);
  densityq->add_task(task310);

  auto task311 = make_shared<Task311>(I729, t2);
  task310->add_dep(task311);
  task311->add_dep(task300);
  densityq->add_task(task311);

  auto I447 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, closed_});
  auto task312 = make_shared<Task312>(den2, I447);
  task312->add_dep(task300);
  densityq->add_task(task312);

  auto I448 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task313 = make_shared<Task313>(I447, t2, I448);
  task312->add_dep(task313);
  task313->add_dep(task300);
  densityq->add_task(task313);

  auto task314 = make_shared<Task314>(I448, Gamma92_(), t2);
  task313->add_dep(task314);
  task314->add_dep(task300);
  densityq->add_task(task314);

  auto I541 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task315 = make_shared<Task315>(I447, t2, I541);
  task312->add_dep(task315);
  task315->add_dep(task300);
  densityq->add_task(task315);

  auto task316 = make_shared<Task316>(I541, Gamma32_(), t2);
  task315->add_dep(task316);
  task316->add_dep(task300);
  densityq->add_task(task316);

  auto I550 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task317 = make_shared<Task317>(I447, t2, I550);
  task312->add_dep(task317);
  task317->add_dep(task300);
  densityq->add_task(task317);

  auto task318 = make_shared<Task318>(I550, Gamma35_(), t2);
  task317->add_dep(task318);
  task318->add_dep(task300);
  densityq->add_task(task318);

  auto I450 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, closed_});
  auto task319 = make_shared<Task319>(den2, I450);
  task319->add_dep(task300);
  densityq->add_task(task319);

  auto I451 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task320 = make_shared<Task320>(I450, t2, I451);
  task319->add_dep(task320);
  task320->add_dep(task300);
  densityq->add_task(task320);

  auto task321 = make_shared<Task321>(I451, Gamma2_(), t2);
  task320->add_dep(task321);
  task321->add_dep(task300);
  densityq->add_task(task321);

  auto I556 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task322 = make_shared<Task322>(I450, t2, I556);
  task319->add_dep(task322);
  task322->add_dep(task300);
  densityq->add_task(task322);

  auto task323 = make_shared<Task323>(I556, Gamma37_(), t2);
  task322->add_dep(task323);
  task323->add_dep(task300);
  densityq->add_task(task323);

  auto I453 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, virt_});
  auto task324 = make_shared<Task324>(den2, I453);
  task324->add_dep(task300);
  densityq->add_task(task324);

  auto I454 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task325 = make_shared<Task325>(I453, t2, I454);
  task324->add_dep(task325);
  task325->add_dep(task300);
  densityq->add_task(task325);

  auto task326 = make_shared<Task326>(I454, Gamma3_(), t2);
  task325->add_dep(task326);
  task326->add_dep(task300);
  densityq->add_task(task326);

  auto I565 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task327 = make_shared<Task327>(I453, t2, I565);
  task324->add_dep(task327);
  task327->add_dep(task300);
  densityq->add_task(task327);

  auto task328 = make_shared<Task328>(I565, Gamma35_(), t2);
  task327->add_dep(task328);
  task328->add_dep(task300);
  densityq->add_task(task328);

  auto I568 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task329 = make_shared<Task329>(I453, t2, I568);
  task324->add_dep(task329);
  task329->add_dep(task300);
  densityq->add_task(task329);

  auto task330 = make_shared<Task330>(I568, Gamma32_(), t2);
  task329->add_dep(task330);
  task330->add_dep(task300);
  densityq->add_task(task330);

  auto I607 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task331 = make_shared<Task331>(I453, t2, I607);
  task324->add_dep(task331);
  task331->add_dep(task300);
  densityq->add_task(task331);

  auto task332 = make_shared<Task332>(I607, Gamma35_(), t2);
  task331->add_dep(task332);
  task332->add_dep(task300);
  densityq->add_task(task332);

  auto I610 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task333 = make_shared<Task333>(I453, t2, I610);
  task324->add_dep(task333);
  task333->add_dep(task300);
  densityq->add_task(task333);

  auto task334 = make_shared<Task334>(I610, Gamma35_(), t2);
  task333->add_dep(task334);
  task334->add_dep(task300);
  densityq->add_task(task334);

  auto I456 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, closed_});
  auto task335 = make_shared<Task335>(den2, I456);
  task335->add_dep(task300);
  densityq->add_task(task335);

  auto I457 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task336 = make_shared<Task336>(I456, t2, I457);
  task335->add_dep(task336);
  task336->add_dep(task300);
  densityq->add_task(task336);

  auto task337 = make_shared<Task337>(I457, Gamma4_(), t2);
  task336->add_dep(task337);
  task337->add_dep(task300);
  densityq->add_task(task337);

  auto I613 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task338 = make_shared<Task338>(I456, t2, I613);
  task335->add_dep(task338);
  task338->add_dep(task300);
  densityq->add_task(task338);

  auto task339 = make_shared<Task339>(I613, Gamma56_(), t2);
  task338->add_dep(task339);
  task339->add_dep(task300);
  densityq->add_task(task339);

  auto I616 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task340 = make_shared<Task340>(I456, t2, I616);
  task335->add_dep(task340);
  task340->add_dep(task300);
  densityq->add_task(task340);

  auto task341 = make_shared<Task341>(I616, Gamma57_(), t2);
  task340->add_dep(task341);
  task341->add_dep(task300);
  densityq->add_task(task341);

  auto I459 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task342 = make_shared<Task342>(den2, I459);
  task342->add_dep(task300);
  densityq->add_task(task342);

  auto I460 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task343 = make_shared<Task343>(I459, Gamma165_(), I460);
  task342->add_dep(task343);
  task343->add_dep(task300);
  densityq->add_task(task343);

  auto task344 = make_shared<Task344>(I460, t2);
  task343->add_dep(task344);
  task344->add_dep(task300);
  densityq->add_task(task344);

  auto I619 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task345 = make_shared<Task345>(I459, Gamma218_(), I619);
  task342->add_dep(task345);
  task345->add_dep(task300);
  densityq->add_task(task345);

  auto task346 = make_shared<Task346>(I619, t2);
  task345->add_dep(task346);
  task346->add_dep(task300);
  densityq->add_task(task346);

  auto I462 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, closed_});
  auto task347 = make_shared<Task347>(den2, I462);
  task347->add_dep(task300);
  densityq->add_task(task347);

  auto I463 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task348 = make_shared<Task348>(I462, t2, I463);
  task347->add_dep(task348);
  task348->add_dep(task300);
  densityq->add_task(task348);

  auto task349 = make_shared<Task349>(I463, Gamma6_(), t2);
  task348->add_dep(task349);
  task349->add_dep(task300);
  densityq->add_task(task349);

  auto I465 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, virt_});
  auto task350 = make_shared<Task350>(den2, I465);
  task350->add_dep(task300);
  densityq->add_task(task350);

  auto I466 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, active_});
  auto task351 = make_shared<Task351>(I465, t2, I466);
  task350->add_dep(task351);
  task351->add_dep(task300);
  densityq->add_task(task351);

  auto task352 = make_shared<Task352>(I466, Gamma7_(), t2);
  task351->add_dep(task352);
  task352->add_dep(task300);
  densityq->add_task(task352);

  auto I469 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, active_});
  auto task353 = make_shared<Task353>(I465, t2, I469);
  task350->add_dep(task353);
  task353->add_dep(task300);
  densityq->add_task(task353);

  auto task354 = make_shared<Task354>(I469, Gamma7_(), t2);
  task353->add_dep(task354);
  task354->add_dep(task300);
  densityq->add_task(task354);

  auto I625 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task355 = make_shared<Task355>(I465, t2, I625);
  task350->add_dep(task355);
  task355->add_dep(task300);
  densityq->add_task(task355);

  auto task356 = make_shared<Task356>(I625, Gamma60_(), t2);
  task355->add_dep(task356);
  task356->add_dep(task300);
  densityq->add_task(task356);

  auto I628 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task357 = make_shared<Task357>(I465, t2, I628);
  task350->add_dep(task357);
  task357->add_dep(task300);
  densityq->add_task(task357);

  auto task358 = make_shared<Task358>(I628, Gamma60_(), t2);
  task357->add_dep(task358);
  task358->add_dep(task300);
  densityq->add_task(task358);

  auto I471 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, virt_});
  auto task359 = make_shared<Task359>(den2, I471);
  task359->add_dep(task300);
  densityq->add_task(task359);

  auto I472 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task360 = make_shared<Task360>(I471, t2, I472);
  task359->add_dep(task360);
  task360->add_dep(task300);
  densityq->add_task(task360);

  auto task361 = make_shared<Task361>(I472, Gamma9_(), t2);
  task360->add_dep(task361);
  task361->add_dep(task300);
  densityq->add_task(task361);

  auto I475 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task362 = make_shared<Task362>(I471, t2, I475);
  task359->add_dep(task362);
  task362->add_dep(task300);
  densityq->add_task(task362);

  auto task363 = make_shared<Task363>(I475, Gamma6_(), t2);
  task362->add_dep(task363);
  task363->add_dep(task300);
  densityq->add_task(task363);

  auto I631 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task364 = make_shared<Task364>(I471, t2, I631);
  task359->add_dep(task364);
  task364->add_dep(task300);
  densityq->add_task(task364);

  auto task365 = make_shared<Task365>(I631, Gamma59_(), t2);
  task364->add_dep(task365);
  task365->add_dep(task300);
  densityq->add_task(task365);

  auto I477 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, virt_});
  auto task366 = make_shared<Task366>(den2, I477);
  task366->add_dep(task300);
  densityq->add_task(task366);

  auto I478 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task367 = make_shared<Task367>(I477, t2, I478);
  task366->add_dep(task367);
  task367->add_dep(task300);
  densityq->add_task(task367);

  auto task368 = make_shared<Task368>(I478, Gamma3_(), t2);
  task367->add_dep(task368);
  task368->add_dep(task300);
  densityq->add_task(task368);

  auto I480 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task369 = make_shared<Task369>(den2, I480);
  task369->add_dep(task300);
  densityq->add_task(task369);

  auto I481 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, active_});
  auto task370 = make_shared<Task370>(I480, t2, I481);
  task369->add_dep(task370);
  task370->add_dep(task300);
  densityq->add_task(task370);

  auto task371 = make_shared<Task371>(I481, Gamma12_(), t2);
  task370->add_dep(task371);
  task371->add_dep(task300);
  densityq->add_task(task371);

  auto I483 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, virt_});
  auto task372 = make_shared<Task372>(den2, I483);
  task372->add_dep(task300);
  densityq->add_task(task372);

  auto I484 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, active_});
  auto task373 = make_shared<Task373>(I483, t2, I484);
  task372->add_dep(task373);
  task373->add_dep(task300);
  densityq->add_task(task373);

  auto task374 = make_shared<Task374>(I484, Gamma12_(), t2);
  task373->add_dep(task374);
  task374->add_dep(task300);
  densityq->add_task(task374);

  auto I640 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task375 = make_shared<Task375>(I483, t2, I640);
  task372->add_dep(task375);
  task375->add_dep(task300);
  densityq->add_task(task375);

  auto I641 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task376 = make_shared<Task376>(I640, Gamma38_(), I641);
  task375->add_dep(task376);
  task376->add_dep(task300);
  densityq->add_task(task376);

  auto task377 = make_shared<Task377>(I641, t2);
  task376->add_dep(task377);
  task377->add_dep(task300);
  densityq->add_task(task377);

  auto I486 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task378 = make_shared<Task378>(den2, I486);
  task378->add_dep(task300);
  densityq->add_task(task378);

  auto I487 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task379 = make_shared<Task379>(I486, Gamma174_(), I487);
  task378->add_dep(task379);
  task379->add_dep(task300);
  densityq->add_task(task379);

  auto task380 = make_shared<Task380>(I487, t2);
  task379->add_dep(task380);
  task380->add_dep(task300);
  densityq->add_task(task380);

  auto task381 = make_shared<Task381>(I487, t2);
  task379->add_dep(task381);
  task381->add_dep(task300);
  densityq->add_task(task381);

  auto I696 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task382 = make_shared<Task382>(I486, Gamma60_(), I696);
  task378->add_dep(task382);
  task382->add_dep(task300);
  densityq->add_task(task382);

  auto task383 = make_shared<Task383>(I696, t2);
  task382->add_dep(task383);
  task383->add_dep(task300);
  densityq->add_task(task383);

  auto task384 = make_shared<Task384>(I696, t2);
  task382->add_dep(task384);
  task384->add_dep(task300);
  densityq->add_task(task384);

  auto I492 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, closed_});
  auto task385 = make_shared<Task385>(den2, I492);
  task385->add_dep(task300);
  densityq->add_task(task385);

  auto I493 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task386 = make_shared<Task386>(I492, t2, I493);
  task385->add_dep(task386);
  task386->add_dep(task300);
  densityq->add_task(task386);

  auto task387 = make_shared<Task387>(I493, Gamma16_(), t2);
  task386->add_dep(task387);
  task387->add_dep(task300);
  densityq->add_task(task387);

  auto I496 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task388 = make_shared<Task388>(I492, t2, I496);
  task385->add_dep(task388);
  task388->add_dep(task300);
  densityq->add_task(task388);

  auto task389 = make_shared<Task389>(I496, Gamma16_(), t2);
  task388->add_dep(task389);
  task389->add_dep(task300);
  densityq->add_task(task389);

  auto I498 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, closed_});
  auto task390 = make_shared<Task390>(den2, I498);
  task390->add_dep(task300);
  densityq->add_task(task390);

  auto I499 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task391 = make_shared<Task391>(I498, t2, I499);
  task390->add_dep(task391);
  task391->add_dep(task300);
  densityq->add_task(task391);

  auto task392 = make_shared<Task392>(I499, Gamma16_(), t2);
  task391->add_dep(task392);
  task392->add_dep(task300);
  densityq->add_task(task392);

  auto I505 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task393 = make_shared<Task393>(I498, t2, I505);
  task390->add_dep(task393);
  task393->add_dep(task300);
  densityq->add_task(task393);

  auto task394 = make_shared<Task394>(I505, Gamma16_(), t2);
  task393->add_dep(task394);
  task394->add_dep(task300);
  densityq->add_task(task394);

  auto I501 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, virt_});
  auto task395 = make_shared<Task395>(den2, I501);
  task395->add_dep(task300);
  densityq->add_task(task395);

  auto I502 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task396 = make_shared<Task396>(I501, t2, I502);
  task395->add_dep(task396);
  task396->add_dep(task300);
  densityq->add_task(task396);

  auto task397 = make_shared<Task397>(I502, Gamma16_(), t2);
  task396->add_dep(task397);
  task397->add_dep(task300);
  densityq->add_task(task397);

  auto I508 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task398 = make_shared<Task398>(I501, t2, I508);
  task395->add_dep(task398);
  task398->add_dep(task300);
  densityq->add_task(task398);

  auto task399 = make_shared<Task399>(I508, Gamma16_(), t2);
  task398->add_dep(task399);
  task399->add_dep(task300);
  densityq->add_task(task399);

  auto I510 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, closed_});
  auto task400 = make_shared<Task400>(den2, I510);
  task400->add_dep(task300);
  densityq->add_task(task400);

  auto I511 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task401 = make_shared<Task401>(I510, t2, I511);
  task400->add_dep(task401);
  task401->add_dep(task300);
  densityq->add_task(task401);

  auto task402 = make_shared<Task402>(I511, Gamma22_(), t2);
  task401->add_dep(task402);
  task402->add_dep(task300);
  densityq->add_task(task402);

  auto task403 = make_shared<Task403>(I511, Gamma12_(), t2);
  task401->add_dep(task403);
  task403->add_dep(task300);
  densityq->add_task(task403);

  auto I513 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, closed_});
  auto task404 = make_shared<Task404>(den2, I513);
  task404->add_dep(task300);
  densityq->add_task(task404);

  auto I514 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task405 = make_shared<Task405>(I513, t2, I514);
  task404->add_dep(task405);
  task405->add_dep(task300);
  densityq->add_task(task405);

  auto I515 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task406 = make_shared<Task406>(I514, Gamma12_(), I515);
  task405->add_dep(task406);
  task406->add_dep(task300);
  densityq->add_task(task406);

  auto task407 = make_shared<Task407>(I515, t2);
  task406->add_dep(task407);
  task407->add_dep(task300);
  densityq->add_task(task407);

  auto I522 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task408 = make_shared<Task408>(den2, I522);
  task408->add_dep(task300);
  densityq->add_task(task408);

  auto I523 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, virt_});
  auto task409 = make_shared<Task409>(I522, Gamma16_(), I523);
  task408->add_dep(task409);
  task409->add_dep(task300);
  densityq->add_task(task409);

  auto task410 = make_shared<Task410>(I523, t2);
  task409->add_dep(task410);
  task410->add_dep(task300);
  densityq->add_task(task410);

  auto task411 = make_shared<Task411>(I523, t2);
  task409->add_dep(task411);
  task411->add_dep(task300);
  densityq->add_task(task411);

  auto I528 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, virt_});
  auto task412 = make_shared<Task412>(den2, I528);
  task412->add_dep(task300);
  densityq->add_task(task412);

  auto I529 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task413 = make_shared<Task413>(I528, t2, I529);
  task412->add_dep(task413);
  task413->add_dep(task300);
  densityq->add_task(task413);

  auto task414 = make_shared<Task414>(I529, Gamma28_(), t2);
  task413->add_dep(task414);
  task414->add_dep(task300);
  densityq->add_task(task414);

  auto I531 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, closed_});
  auto task415 = make_shared<Task415>(den2, I531);
  task415->add_dep(task300);
  densityq->add_task(task415);

  auto I532 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task416 = make_shared<Task416>(I531, t2, I532);
  task415->add_dep(task416);
  task416->add_dep(task300);
  densityq->add_task(task416);

  auto task417 = make_shared<Task417>(I532, Gamma29_(), t2);
  task416->add_dep(task417);
  task417->add_dep(task300);
  densityq->add_task(task417);

  auto I535 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task418 = make_shared<Task418>(I531, t2, I535);
  task415->add_dep(task418);
  task418->add_dep(task300);
  densityq->add_task(task418);

  auto task419 = make_shared<Task419>(I535, Gamma7_(), t2);
  task418->add_dep(task419);
  task419->add_dep(task300);
  densityq->add_task(task419);

  auto I574 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task420 = make_shared<Task420>(I531, t2, I574);
  task415->add_dep(task420);
  task420->add_dep(task300);
  densityq->add_task(task420);

  auto task421 = make_shared<Task421>(I574, Gamma7_(), t2);
  task420->add_dep(task421);
  task421->add_dep(task300);
  densityq->add_task(task421);

  auto I577 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task422 = make_shared<Task422>(I531, t2, I577);
  task415->add_dep(task422);
  task422->add_dep(task300);
  densityq->add_task(task422);

  auto task423 = make_shared<Task423>(I577, Gamma7_(), t2);
  task422->add_dep(task423);
  task423->add_dep(task300);
  densityq->add_task(task423);

  auto I726 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task424 = make_shared<Task424>(I531, t2, I726);
  task415->add_dep(task424);
  task424->add_dep(task300);
  densityq->add_task(task424);

  auto task425 = make_shared<Task425>(I726, Gamma60_(), t2);
  task424->add_dep(task425);
  task425->add_dep(task300);
  densityq->add_task(task425);

  auto I543 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, virt_});
  auto task426 = make_shared<Task426>(den2, I543);
  task426->add_dep(task300);
  densityq->add_task(task426);

  auto I544 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task427 = make_shared<Task427>(I543, t2, I544);
  task426->add_dep(task427);
  task427->add_dep(task300);
  densityq->add_task(task427);

  auto task428 = make_shared<Task428>(I544, Gamma32_(), t2);
  task427->add_dep(task428);
  task428->add_dep(task300);
  densityq->add_task(task428);

  auto I553 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, active_, active_});
  auto task429 = make_shared<Task429>(I543, t2, I553);
  task426->add_dep(task429);
  task429->add_dep(task300);
  densityq->add_task(task429);

  auto task430 = make_shared<Task430>(I553, Gamma35_(), t2);
  task429->add_dep(task430);
  task430->add_dep(task300);
  densityq->add_task(task430);

  auto I558 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task431 = make_shared<Task431>(den2, I558);
  task431->add_dep(task300);
  densityq->add_task(task431);

  auto I559 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, virt_});
  auto task432 = make_shared<Task432>(I558, t2, I559);
  task431->add_dep(task432);
  task432->add_dep(task300);
  densityq->add_task(task432);

  auto task433 = make_shared<Task433>(I559, Gamma38_(), t2);
  task432->add_dep(task433);
  task433->add_dep(task300);
  densityq->add_task(task433);

  auto I562 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, virt_});
  auto task434 = make_shared<Task434>(I558, t2, I562);
  task431->add_dep(task434);
  task434->add_dep(task300);
  densityq->add_task(task434);

  auto task435 = make_shared<Task435>(I562, Gamma38_(), t2);
  task434->add_dep(task435);
  task435->add_dep(task300);
  densityq->add_task(task435);

  auto I601 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task436 = make_shared<Task436>(I558, t2, I601);
  task431->add_dep(task436);
  task436->add_dep(task300);
  densityq->add_task(task436);

  auto task437 = make_shared<Task437>(I601, Gamma38_(), t2);
  task436->add_dep(task437);
  task437->add_dep(task300);
  densityq->add_task(task437);

  auto I604 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task438 = make_shared<Task438>(I558, t2, I604);
  task431->add_dep(task438);
  task438->add_dep(task300);
  densityq->add_task(task438);

  auto task439 = make_shared<Task439>(I604, Gamma38_(), t2);
  task438->add_dep(task439);
  task439->add_dep(task300);
  densityq->add_task(task439);

  auto I570 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, virt_});
  auto task440 = make_shared<Task440>(den2, I570);
  task440->add_dep(task300);
  densityq->add_task(task440);

  auto I571 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task441 = make_shared<Task441>(I570, t2, I571);
  task440->add_dep(task441);
  task441->add_dep(task300);
  densityq->add_task(task441);

  auto task442 = make_shared<Task442>(I571, Gamma6_(), t2);
  task441->add_dep(task442);
  task442->add_dep(task300);
  densityq->add_task(task442);

  auto I723 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task443 = make_shared<Task443>(I570, t2, I723);
  task440->add_dep(task443);
  task443->add_dep(task300);
  densityq->add_task(task443);

  auto task444 = make_shared<Task444>(I723, Gamma59_(), t2);
  task443->add_dep(task444);
  task444->add_dep(task300);
  densityq->add_task(task444);

  auto I582 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, closed_});
  auto task445 = make_shared<Task445>(den2, I582);
  task445->add_dep(task300);
  densityq->add_task(task445);

  auto I583 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task446 = make_shared<Task446>(I582, t2, I583);
  task445->add_dep(task446);
  task446->add_dep(task300);
  densityq->add_task(task446);

  auto task447 = make_shared<Task447>(I583, Gamma35_(), t2);
  task446->add_dep(task447);
  task447->add_dep(task300);
  densityq->add_task(task447);

  auto I592 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task448 = make_shared<Task448>(I582, t2, I592);
  task445->add_dep(task448);
  task448->add_dep(task300);
  densityq->add_task(task448);

  auto task449 = make_shared<Task449>(I592, Gamma35_(), t2);
  task448->add_dep(task449);
  task449->add_dep(task300);
  densityq->add_task(task449);

  auto I585 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, virt_});
  auto task450 = make_shared<Task450>(den2, I585);
  task450->add_dep(task300);
  densityq->add_task(task450);

  auto I586 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task451 = make_shared<Task451>(I585, t2, I586);
  task450->add_dep(task451);
  task451->add_dep(task300);
  densityq->add_task(task451);

  auto task452 = make_shared<Task452>(I586, Gamma35_(), t2);
  task451->add_dep(task452);
  task452->add_dep(task300);
  densityq->add_task(task452);

  auto I595 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task453 = make_shared<Task453>(I585, t2, I595);
  task450->add_dep(task453);
  task453->add_dep(task300);
  densityq->add_task(task453);

  auto task454 = make_shared<Task454>(I595, Gamma35_(), t2);
  task453->add_dep(task454);
  task454->add_dep(task300);
  densityq->add_task(task454);

  auto I732 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task455 = make_shared<Task455>(I585, t2, I732);
  task450->add_dep(task455);
  task455->add_dep(task300);
  densityq->add_task(task455);

  auto task456 = make_shared<Task456>(I732, Gamma60_(), t2);
  task455->add_dep(task456);
  task456->add_dep(task300);
  densityq->add_task(task456);

  auto I597 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, closed_});
  auto task457 = make_shared<Task457>(den2, I597);
  task457->add_dep(task300);
  densityq->add_task(task457);

  auto I598 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task458 = make_shared<Task458>(I597, t2, I598);
  task457->add_dep(task458);
  task458->add_dep(task300);
  densityq->add_task(task458);

  auto task459 = make_shared<Task459>(I598, Gamma51_(), t2);
  task458->add_dep(task459);
  task459->add_dep(task300);
  densityq->add_task(task459);

  auto I621 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, virt_});
  auto task460 = make_shared<Task460>(den2, I621);
  task460->add_dep(task300);
  densityq->add_task(task460);

  auto I622 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task461 = make_shared<Task461>(I621, t2, I622);
  task460->add_dep(task461);
  task461->add_dep(task300);
  densityq->add_task(task461);

  auto task462 = make_shared<Task462>(I622, Gamma59_(), t2);
  task461->add_dep(task462);
  task462->add_dep(task300);
  densityq->add_task(task462);

  auto I633 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task463 = make_shared<Task463>(den2, I633);
  task463->add_dep(task300);
  densityq->add_task(task463);

  auto I634 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task464 = make_shared<Task464>(I633, Gamma16_(), I634);
  task463->add_dep(task464);
  task464->add_dep(task300);
  densityq->add_task(task464);

  auto task465 = make_shared<Task465>(I634, t2);
  task464->add_dep(task465);
  task465->add_dep(task300);
  densityq->add_task(task465);

  auto I636 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task466 = make_shared<Task466>(den2, I636);
  task466->add_dep(task300);
  densityq->add_task(task466);

  auto I637 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task467 = make_shared<Task467>(I636, Gamma16_(), I637);
  task466->add_dep(task467);
  task467->add_dep(task300);
  densityq->add_task(task467);

  auto task468 = make_shared<Task468>(I637, t2);
  task467->add_dep(task468);
  task468->add_dep(task300);
  densityq->add_task(task468);

  auto I642 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task469 = make_shared<Task469>(den2, I642);
  task469->add_dep(task300);
  densityq->add_task(task469);

  auto I643 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task470 = make_shared<Task470>(I642, t2, I643);
  task469->add_dep(task470);
  task470->add_dep(task300);
  densityq->add_task(task470);

  auto I644 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task471 = make_shared<Task471>(I643, Gamma38_(), I644);
  task470->add_dep(task471);
  task471->add_dep(task300);
  densityq->add_task(task471);

  auto task472 = make_shared<Task472>(I644, t2);
  task471->add_dep(task472);
  task472->add_dep(task300);
  densityq->add_task(task472);

  auto I651 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task473 = make_shared<Task473>(den2, I651);
  task473->add_dep(task300);
  densityq->add_task(task473);

  auto I652 = make_shared<TATensor<double,0>>(std::vector<IndexRange>{});
  auto task474 = make_shared<Task474>(I651, Gamma38_(), I652);
  task473->add_dep(task474);
  task474->add_dep(task300);
  densityq->add_task(task474);

  auto task475 = make_shared<Task475>(I652, t2);
  task474->add_dep(task475);
  task475->add_dep(task300);
  densityq->add_task(task475);

  auto task476 = make_shared<Task476>(I652, t2);
  task474->add_dep(task476);
  task476->add_dep(task300);
  densityq->add_task(task476);

  shared_ptr<TATensor<double,2>> I657;
  if (diagonal) {
    I657 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, closed_});
  }
  shared_ptr<Task477> task477;
  if (diagonal) {
    task477 = make_shared<Task477>(den2, I657);
    task477->add_dep(task300);
    densityq->add_task(task477);
  }

  shared_ptr<Task478> task478;
  if (diagonal) {
    task478 = make_shared<Task478>(I657, t2);
    task477->add_dep(task478);
    task478->add_dep(task300);
    densityq->add_task(task478);
  }

  shared_ptr<Task479> task479;
  if (diagonal) {
    task479 = make_shared<Task479>(I657, t2);
    task477->add_dep(task479);
    task479->add_dep(task300);
    densityq->add_task(task479);
  }

  shared_ptr<TATensor<double,2>> I661;
  if (diagonal) {
    I661 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, virt_});
  }
  shared_ptr<Task480> task480;
  if (diagonal) {
    task480 = make_shared<Task480>(den2, I661);
    task480->add_dep(task300);
    densityq->add_task(task480);
  }

  shared_ptr<Task481> task481;
  if (diagonal) {
    task481 = make_shared<Task481>(I661, t2);
    task480->add_dep(task481);
    task481->add_dep(task300);
    densityq->add_task(task481);
  }

  shared_ptr<Task482> task482;
  if (diagonal) {
    task482 = make_shared<Task482>(I661, t2);
    task480->add_dep(task482);
    task482->add_dep(task300);
    densityq->add_task(task482);
  }

  auto I665 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, active_});
  auto task483 = make_shared<Task483>(den2, I665);
  task483->add_dep(task300);
  densityq->add_task(task483);

  auto I666 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, active_});
  auto task484 = make_shared<Task484>(I665, Gamma38_(), I666);
  task483->add_dep(task484);
  task484->add_dep(task300);
  densityq->add_task(task484);

  auto task485 = make_shared<Task485>(I666, t2);
  task484->add_dep(task485);
  task485->add_dep(task300);
  densityq->add_task(task485);

  auto task486 = make_shared<Task486>(I666, t2);
  task484->add_dep(task486);
  task486->add_dep(task300);
  densityq->add_task(task486);

  auto I671 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, virt_});
  auto task487 = make_shared<Task487>(den2, I671);
  task487->add_dep(task300);
  densityq->add_task(task487);

  auto I672 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task488 = make_shared<Task488>(I671, t2, I672);
  task487->add_dep(task488);
  task488->add_dep(task300);
  densityq->add_task(task488);

  auto I673 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task489 = make_shared<Task489>(I672, Gamma35_(), I673);
  task488->add_dep(task489);
  task489->add_dep(task300);
  densityq->add_task(task489);

  auto task490 = make_shared<Task490>(I673, t2);
  task489->add_dep(task490);
  task490->add_dep(task300);
  densityq->add_task(task490);

  auto I674 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, virt_});
  auto task491 = make_shared<Task491>(den2, I674);
  task491->add_dep(task300);
  densityq->add_task(task491);

  auto I675 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task492 = make_shared<Task492>(I674, t2, I675);
  task491->add_dep(task492);
  task492->add_dep(task300);
  densityq->add_task(task492);

  auto task493 = make_shared<Task493>(I675, Gamma32_(), t2);
  task492->add_dep(task493);
  task493->add_dep(task300);
  densityq->add_task(task493);

  auto task494 = make_shared<Task494>(I675, Gamma35_(), t2);
  task492->add_dep(task494);
  task494->add_dep(task300);
  densityq->add_task(task494);

  auto I683 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, virt_});
  auto task495 = make_shared<Task495>(den2, I683);
  task495->add_dep(task300);
  densityq->add_task(task495);

  auto I684 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task496 = make_shared<Task496>(I683, t2, I684);
  task495->add_dep(task496);
  task496->add_dep(task300);
  densityq->add_task(task496);

  auto task497 = make_shared<Task497>(I684, Gamma60_(), t2);
  task496->add_dep(task497);
  task497->add_dep(task300);
  densityq->add_task(task497);

  auto I686 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task498 = make_shared<Task498>(den2, I686);
  task498->add_dep(task300);
  densityq->add_task(task498);

  auto I687 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task499 = make_shared<Task499>(I686, t2, I687);
  task498->add_dep(task499);
  task499->add_dep(task300);
  densityq->add_task(task499);

  auto task500 = make_shared<Task500>(I687, Gamma60_(), t2);
  task499->add_dep(task500);
  task500->add_dep(task300);
  densityq->add_task(task500);

  auto I689 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, active_});
  auto task501 = make_shared<Task501>(den2, I689);
  task501->add_dep(task300);
  densityq->add_task(task501);

  auto I690 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, closed_});
  auto task502 = make_shared<Task502>(I689, Gamma38_(), I690);
  task501->add_dep(task502);
  task502->add_dep(task300);
  densityq->add_task(task502);

  auto task503 = make_shared<Task503>(I690, t2);
  task502->add_dep(task503);
  task503->add_dep(task300);
  densityq->add_task(task503);

  auto task504 = make_shared<Task504>(I690, t2);
  task502->add_dep(task504);
  task504->add_dep(task300);
  densityq->add_task(task504);

  auto I701 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, closed_});
  auto task505 = make_shared<Task505>(den2, I701);
  task505->add_dep(task300);
  densityq->add_task(task505);

  auto I702 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task506 = make_shared<Task506>(I701, t2, I702);
  task505->add_dep(task506);
  task506->add_dep(task300);
  densityq->add_task(task506);

  auto task507 = make_shared<Task507>(I702, Gamma38_(), t2);
  task506->add_dep(task507);
  task507->add_dep(task300);
  densityq->add_task(task507);

  auto I705 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task508 = make_shared<Task508>(I701, t2, I705);
  task505->add_dep(task508);
  task508->add_dep(task300);
  densityq->add_task(task508);

  auto task509 = make_shared<Task509>(I705, Gamma38_(), t2);
  task508->add_dep(task509);
  task509->add_dep(task300);
  densityq->add_task(task509);

  auto I707 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, virt_});
  auto task510 = make_shared<Task510>(den2, I707);
  task510->add_dep(task300);
  densityq->add_task(task510);

  auto I708 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task511 = make_shared<Task511>(I707, t2, I708);
  task510->add_dep(task511);
  task511->add_dep(task300);
  densityq->add_task(task511);

  auto task512 = make_shared<Task512>(I708, Gamma38_(), t2);
  task511->add_dep(task512);
  task512->add_dep(task300);
  densityq->add_task(task512);

  auto I711 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task513 = make_shared<Task513>(I707, t2, I711);
  task510->add_dep(task513);
  task513->add_dep(task300);
  densityq->add_task(task513);

  auto task514 = make_shared<Task514>(I711, Gamma38_(), t2);
  task513->add_dep(task514);
  task514->add_dep(task300);
  densityq->add_task(task514);

  auto I713 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, virt_});
  auto task515 = make_shared<Task515>(den2, I713);
  task515->add_dep(task300);
  densityq->add_task(task515);

  auto I714 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task516 = make_shared<Task516>(I713, t2, I714);
  task515->add_dep(task516);
  task516->add_dep(task300);
  densityq->add_task(task516);

  auto task517 = make_shared<Task517>(I714, Gamma38_(), t2);
  task516->add_dep(task517);
  task517->add_dep(task300);
  densityq->add_task(task517);

  auto I717 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task518 = make_shared<Task518>(I713, t2, I717);
  task515->add_dep(task518);
  task518->add_dep(task300);
  densityq->add_task(task518);

  auto task519 = make_shared<Task519>(I717, Gamma38_(), t2);
  task518->add_dep(task519);
  task519->add_dep(task300);
  densityq->add_task(task519);

  auto I719 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, closed_});
  auto task520 = make_shared<Task520>(den2, I719);
  task520->add_dep(task300);
  densityq->add_task(task520);

  auto I720 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task521 = make_shared<Task521>(I719, t2, I720);
  task520->add_dep(task521);
  task521->add_dep(task300);
  densityq->add_task(task521);

  auto task522 = make_shared<Task522>(I720, Gamma60_(), t2);
  task521->add_dep(task522);
  task522->add_dep(task300);
  densityq->add_task(task522);

  return densityq;
}


#endif
