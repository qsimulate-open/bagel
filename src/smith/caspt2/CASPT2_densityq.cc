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


#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_densityq(const bool reset, const bool diagonal) {

  auto densityq = make_shared<Queue>();
  auto task283 = make_shared<Task283>(den2, reset);
  densityq->add_task(task283);

  auto I374 = make_shared<TATensor<double,2>>({active_, active_});
  auto task284 = make_shared<Task284>(den2, I374);
  task284->add_dep(task283);
  densityq->add_task(task284);

  auto I375 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task285 = make_shared<Task285>(I374, Gamma138_(), I375);
  task284->add_dep(task285);
  task285->add_dep(task283);
  densityq->add_task(task285);

  auto task286 = make_shared<Task286>(I375, t2);
  task285->add_dep(task286);
  task286->add_dep(task283);
  densityq->add_task(task286);

  auto I468 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task287 = make_shared<Task287>(I374, Gamma169_(), I468);
  task284->add_dep(task287);
  task287->add_dep(task283);
  densityq->add_task(task287);

  auto task288 = make_shared<Task288>(I468, t2);
  task287->add_dep(task288);
  task288->add_dep(task283);
  densityq->add_task(task288);

  auto I477 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task289 = make_shared<Task289>(I374, Gamma172_(), I477);
  task284->add_dep(task289);
  task289->add_dep(task283);
  densityq->add_task(task289);

  auto task290 = make_shared<Task290>(I477, t2);
  task289->add_dep(task290);
  task290->add_dep(task283);
  densityq->add_task(task290);

  auto task291 = make_shared<Task291>(I477, t2);
  task289->add_dep(task291);
  task291->add_dep(task283);
  densityq->add_task(task291);

  auto task292 = make_shared<Task292>(I477, t2);
  task289->add_dep(task292);
  task292->add_dep(task283);
  densityq->add_task(task292);

  auto I659 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task293 = make_shared<Task293>(I374, Gamma230_(), I659);
  task284->add_dep(task293);
  task293->add_dep(task283);
  densityq->add_task(task293);

  auto task294 = make_shared<Task294>(I659, t2);
  task293->add_dep(task294);
  task294->add_dep(task283);
  densityq->add_task(task294);

  auto I377 = make_shared<TATensor<double,2>>({closed_, closed_});
  auto task295 = make_shared<Task295>(den2, I377);
  task295->add_dep(task283);
  densityq->add_task(task295);

  auto I378 = make_shared<TATensor<double,4>>({closed_, closed_, active_, active_});
  auto task296 = make_shared<Task296>(I377, t2, I378);
  task295->add_dep(task296);
  task296->add_dep(task283);
  densityq->add_task(task296);

  auto task297 = make_shared<Task297>(I378, Gamma92_(), t2);
  task296->add_dep(task297);
  task297->add_dep(task283);
  densityq->add_task(task297);

  auto I471 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task298 = make_shared<Task298>(I377, t2, I471);
  task295->add_dep(task298);
  task298->add_dep(task283);
  densityq->add_task(task298);

  auto task299 = make_shared<Task299>(I471, Gamma32_(), t2);
  task298->add_dep(task299);
  task299->add_dep(task283);
  densityq->add_task(task299);

  auto I480 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task300 = make_shared<Task300>(I377, t2, I480);
  task295->add_dep(task300);
  task300->add_dep(task283);
  densityq->add_task(task300);

  auto task301 = make_shared<Task301>(I480, Gamma35_(), t2);
  task300->add_dep(task301);
  task301->add_dep(task283);
  densityq->add_task(task301);

  auto I380 = make_shared<TATensor<double,2>>({active_, closed_});
  auto task302 = make_shared<Task302>(den2, I380);
  task302->add_dep(task283);
  densityq->add_task(task302);

  auto I381 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task303 = make_shared<Task303>(I380, t2, I381);
  task302->add_dep(task303);
  task303->add_dep(task283);
  densityq->add_task(task303);

  auto task304 = make_shared<Task304>(I381, Gamma2_(), t2);
  task303->add_dep(task304);
  task304->add_dep(task283);
  densityq->add_task(task304);

  auto I486 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task305 = make_shared<Task305>(I380, t2, I486);
  task302->add_dep(task305);
  task305->add_dep(task283);
  densityq->add_task(task305);

  auto task306 = make_shared<Task306>(I486, Gamma37_(), t2);
  task305->add_dep(task306);
  task306->add_dep(task283);
  densityq->add_task(task306);

  auto I383 = make_shared<TATensor<double,2>>({active_, virt_});
  auto task307 = make_shared<Task307>(den2, I383);
  task307->add_dep(task283);
  densityq->add_task(task307);

  auto I384 = make_shared<TATensor<double,4>>({closed_, closed_, active_, active_});
  auto task308 = make_shared<Task308>(I383, t2, I384);
  task307->add_dep(task308);
  task308->add_dep(task283);
  densityq->add_task(task308);

  auto task309 = make_shared<Task309>(I384, Gamma3_(), t2);
  task308->add_dep(task309);
  task309->add_dep(task283);
  densityq->add_task(task309);

  auto I495 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task310 = make_shared<Task310>(I383, t2, I495);
  task307->add_dep(task310);
  task310->add_dep(task283);
  densityq->add_task(task310);

  auto task311 = make_shared<Task311>(I495, Gamma35_(), t2);
  task310->add_dep(task311);
  task311->add_dep(task283);
  densityq->add_task(task311);

  auto I498 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task312 = make_shared<Task312>(I383, t2, I498);
  task307->add_dep(task312);
  task312->add_dep(task283);
  densityq->add_task(task312);

  auto task313 = make_shared<Task313>(I498, Gamma32_(), t2);
  task312->add_dep(task313);
  task313->add_dep(task283);
  densityq->add_task(task313);

  auto I537 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task314 = make_shared<Task314>(I383, t2, I537);
  task307->add_dep(task314);
  task314->add_dep(task283);
  densityq->add_task(task314);

  auto task315 = make_shared<Task315>(I537, Gamma35_(), t2);
  task314->add_dep(task315);
  task315->add_dep(task283);
  densityq->add_task(task315);

  auto I540 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task316 = make_shared<Task316>(I383, t2, I540);
  task307->add_dep(task316);
  task316->add_dep(task283);
  densityq->add_task(task316);

  auto task317 = make_shared<Task317>(I540, Gamma35_(), t2);
  task316->add_dep(task317);
  task317->add_dep(task283);
  densityq->add_task(task317);

  auto I386 = make_shared<TATensor<double,2>>({active_, closed_});
  auto task318 = make_shared<Task318>(den2, I386);
  task318->add_dep(task283);
  densityq->add_task(task318);

  auto I387 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task319 = make_shared<Task319>(I386, t2, I387);
  task318->add_dep(task319);
  task319->add_dep(task283);
  densityq->add_task(task319);

  auto task320 = make_shared<Task320>(I387, Gamma4_(), t2);
  task319->add_dep(task320);
  task320->add_dep(task283);
  densityq->add_task(task320);

  auto I543 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task321 = make_shared<Task321>(I386, t2, I543);
  task318->add_dep(task321);
  task321->add_dep(task283);
  densityq->add_task(task321);

  auto task322 = make_shared<Task322>(I543, Gamma56_(), t2);
  task321->add_dep(task322);
  task322->add_dep(task283);
  densityq->add_task(task322);

  auto I546 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task323 = make_shared<Task323>(I386, t2, I546);
  task318->add_dep(task323);
  task323->add_dep(task283);
  densityq->add_task(task323);

  auto task324 = make_shared<Task324>(I546, Gamma57_(), t2);
  task323->add_dep(task324);
  task324->add_dep(task283);
  densityq->add_task(task324);

  auto I389 = make_shared<TATensor<double,2>>({active_, active_});
  auto task325 = make_shared<Task325>(den2, I389);
  task325->add_dep(task283);
  densityq->add_task(task325);

  auto I390 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task326 = make_shared<Task326>(I389, Gamma143_(), I390);
  task325->add_dep(task326);
  task326->add_dep(task283);
  densityq->add_task(task326);

  auto task327 = make_shared<Task327>(I390, t2);
  task326->add_dep(task327);
  task327->add_dep(task283);
  densityq->add_task(task327);

  auto I549 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task328 = make_shared<Task328>(I389, Gamma196_(), I549);
  task325->add_dep(task328);
  task328->add_dep(task283);
  densityq->add_task(task328);

  auto task329 = make_shared<Task329>(I549, t2);
  task328->add_dep(task329);
  task329->add_dep(task283);
  densityq->add_task(task329);

  auto I392 = make_shared<TATensor<double,2>>({closed_, closed_});
  auto task330 = make_shared<Task330>(den2, I392);
  task330->add_dep(task283);
  densityq->add_task(task330);

  auto I393 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task331 = make_shared<Task331>(I392, t2, I393);
  task330->add_dep(task331);
  task331->add_dep(task283);
  densityq->add_task(task331);

  auto task332 = make_shared<Task332>(I393, Gamma6_(), t2);
  task331->add_dep(task332);
  task332->add_dep(task283);
  densityq->add_task(task332);

  auto I395 = make_shared<TATensor<double,2>>({closed_, virt_});
  auto task333 = make_shared<Task333>(den2, I395);
  task333->add_dep(task283);
  densityq->add_task(task333);

  auto I396 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task334 = make_shared<Task334>(I395, t2, I396);
  task333->add_dep(task334);
  task334->add_dep(task283);
  densityq->add_task(task334);

  auto task335 = make_shared<Task335>(I396, Gamma7_(), t2);
  task334->add_dep(task335);
  task335->add_dep(task283);
  densityq->add_task(task335);

  auto I399 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task336 = make_shared<Task336>(I395, t2, I399);
  task333->add_dep(task336);
  task336->add_dep(task283);
  densityq->add_task(task336);

  auto task337 = make_shared<Task337>(I399, Gamma7_(), t2);
  task336->add_dep(task337);
  task337->add_dep(task283);
  densityq->add_task(task337);

  auto I555 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task338 = make_shared<Task338>(I395, t2, I555);
  task333->add_dep(task338);
  task338->add_dep(task283);
  densityq->add_task(task338);

  auto task339 = make_shared<Task339>(I555, Gamma60_(), t2);
  task338->add_dep(task339);
  task339->add_dep(task283);
  densityq->add_task(task339);

  auto I558 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task340 = make_shared<Task340>(I395, t2, I558);
  task333->add_dep(task340);
  task340->add_dep(task283);
  densityq->add_task(task340);

  auto task341 = make_shared<Task341>(I558, Gamma60_(), t2);
  task340->add_dep(task341);
  task341->add_dep(task283);
  densityq->add_task(task341);

  auto I401 = make_shared<TATensor<double,2>>({active_, virt_});
  auto task342 = make_shared<Task342>(den2, I401);
  task342->add_dep(task283);
  densityq->add_task(task342);

  auto I402 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task343 = make_shared<Task343>(I401, t2, I402);
  task342->add_dep(task343);
  task343->add_dep(task283);
  densityq->add_task(task343);

  auto task344 = make_shared<Task344>(I402, Gamma9_(), t2);
  task343->add_dep(task344);
  task344->add_dep(task283);
  densityq->add_task(task344);

  auto I405 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task345 = make_shared<Task345>(I401, t2, I405);
  task342->add_dep(task345);
  task345->add_dep(task283);
  densityq->add_task(task345);

  auto task346 = make_shared<Task346>(I405, Gamma6_(), t2);
  task345->add_dep(task346);
  task346->add_dep(task283);
  densityq->add_task(task346);

  auto I561 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task347 = make_shared<Task347>(I401, t2, I561);
  task342->add_dep(task347);
  task347->add_dep(task283);
  densityq->add_task(task347);

  auto task348 = make_shared<Task348>(I561, Gamma59_(), t2);
  task347->add_dep(task348);
  task348->add_dep(task283);
  densityq->add_task(task348);

  auto I407 = make_shared<TATensor<double,2>>({active_, virt_});
  auto task349 = make_shared<Task349>(den2, I407);
  task349->add_dep(task283);
  densityq->add_task(task349);

  auto I408 = make_shared<TATensor<double,4>>({closed_, closed_, active_, active_});
  auto task350 = make_shared<Task350>(I407, t2, I408);
  task349->add_dep(task350);
  task350->add_dep(task283);
  densityq->add_task(task350);

  auto task351 = make_shared<Task351>(I408, Gamma3_(), t2);
  task350->add_dep(task351);
  task351->add_dep(task283);
  densityq->add_task(task351);

  auto I410 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task352 = make_shared<Task352>(den2, I410);
  task352->add_dep(task283);
  densityq->add_task(task352);

  auto I411 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task353 = make_shared<Task353>(I410, t2, I411);
  task352->add_dep(task353);
  task353->add_dep(task283);
  densityq->add_task(task353);

  auto task354 = make_shared<Task354>(I411, Gamma12_(), t2);
  task353->add_dep(task354);
  task354->add_dep(task283);
  densityq->add_task(task354);

  auto I413 = make_shared<TATensor<double,2>>({closed_, virt_});
  auto task355 = make_shared<Task355>(den2, I413);
  task355->add_dep(task283);
  densityq->add_task(task355);

  auto I414 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task356 = make_shared<Task356>(I413, t2, I414);
  task355->add_dep(task356);
  task356->add_dep(task283);
  densityq->add_task(task356);

  auto task357 = make_shared<Task357>(I414, Gamma12_(), t2);
  task356->add_dep(task357);
  task357->add_dep(task283);
  densityq->add_task(task357);

  auto I570 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task358 = make_shared<Task358>(I413, t2, I570);
  task355->add_dep(task358);
  task358->add_dep(task283);
  densityq->add_task(task358);

  auto I571 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task359 = make_shared<Task359>(I570, Gamma38_(), I571);
  task358->add_dep(task359);
  task359->add_dep(task283);
  densityq->add_task(task359);

  auto task360 = make_shared<Task360>(I571, t2);
  task359->add_dep(task360);
  task360->add_dep(task283);
  densityq->add_task(task360);

  auto I416 = make_shared<TATensor<double,2>>({active_, active_});
  auto task361 = make_shared<Task361>(den2, I416);
  task361->add_dep(task283);
  densityq->add_task(task361);

  auto I417 = make_shared<TATensor<double,2>>({active_, active_});
  auto task362 = make_shared<Task362>(I416, Gamma152_(), I417);
  task361->add_dep(task362);
  task362->add_dep(task283);
  densityq->add_task(task362);

  auto task363 = make_shared<Task363>(I417, t2);
  task362->add_dep(task363);
  task363->add_dep(task283);
  densityq->add_task(task363);

  auto task364 = make_shared<Task364>(I417, t2);
  task362->add_dep(task364);
  task364->add_dep(task283);
  densityq->add_task(task364);

  auto I626 = make_shared<TATensor<double,2>>({active_, active_});
  auto task365 = make_shared<Task365>(I416, Gamma60_(), I626);
  task361->add_dep(task365);
  task365->add_dep(task283);
  densityq->add_task(task365);

  auto task366 = make_shared<Task366>(I626, t2);
  task365->add_dep(task366);
  task366->add_dep(task283);
  densityq->add_task(task366);

  auto task367 = make_shared<Task367>(I626, t2);
  task365->add_dep(task367);
  task367->add_dep(task283);
  densityq->add_task(task367);

  auto I422 = make_shared<TATensor<double,2>>({closed_, closed_});
  auto task368 = make_shared<Task368>(den2, I422);
  task368->add_dep(task283);
  densityq->add_task(task368);

  auto I423 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, active_});
  auto task369 = make_shared<Task369>(I422, t2, I423);
  task368->add_dep(task369);
  task369->add_dep(task283);
  densityq->add_task(task369);

  auto task370 = make_shared<Task370>(I423, Gamma16_(), t2);
  task369->add_dep(task370);
  task370->add_dep(task283);
  densityq->add_task(task370);

  auto I426 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, active_});
  auto task371 = make_shared<Task371>(I422, t2, I426);
  task368->add_dep(task371);
  task371->add_dep(task283);
  densityq->add_task(task371);

  auto task372 = make_shared<Task372>(I426, Gamma16_(), t2);
  task371->add_dep(task372);
  task372->add_dep(task283);
  densityq->add_task(task372);

  auto I428 = make_shared<TATensor<double,2>>({closed_, closed_});
  auto task373 = make_shared<Task373>(den2, I428);
  task373->add_dep(task283);
  densityq->add_task(task373);

  auto I429 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, active_});
  auto task374 = make_shared<Task374>(I428, t2, I429);
  task373->add_dep(task374);
  task374->add_dep(task283);
  densityq->add_task(task374);

  auto task375 = make_shared<Task375>(I429, Gamma16_(), t2);
  task374->add_dep(task375);
  task375->add_dep(task283);
  densityq->add_task(task375);

  auto I435 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, active_});
  auto task376 = make_shared<Task376>(I428, t2, I435);
  task373->add_dep(task376);
  task376->add_dep(task283);
  densityq->add_task(task376);

  auto task377 = make_shared<Task377>(I435, Gamma16_(), t2);
  task376->add_dep(task377);
  task377->add_dep(task283);
  densityq->add_task(task377);

  auto I431 = make_shared<TATensor<double,2>>({virt_, virt_});
  auto task378 = make_shared<Task378>(den2, I431);
  task378->add_dep(task283);
  densityq->add_task(task378);

  auto I432 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, active_});
  auto task379 = make_shared<Task379>(I431, t2, I432);
  task378->add_dep(task379);
  task379->add_dep(task283);
  densityq->add_task(task379);

  auto task380 = make_shared<Task380>(I432, Gamma16_(), t2);
  task379->add_dep(task380);
  task380->add_dep(task283);
  densityq->add_task(task380);

  auto I438 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, active_});
  auto task381 = make_shared<Task381>(I431, t2, I438);
  task378->add_dep(task381);
  task381->add_dep(task283);
  densityq->add_task(task381);

  auto task382 = make_shared<Task382>(I438, Gamma16_(), t2);
  task381->add_dep(task382);
  task382->add_dep(task283);
  densityq->add_task(task382);

  auto I440 = make_shared<TATensor<double,2>>({active_, closed_});
  auto task383 = make_shared<Task383>(den2, I440);
  task383->add_dep(task283);
  densityq->add_task(task383);

  auto I441 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task384 = make_shared<Task384>(I440, t2, I441);
  task383->add_dep(task384);
  task384->add_dep(task283);
  densityq->add_task(task384);

  auto task385 = make_shared<Task385>(I441, Gamma22_(), t2);
  task384->add_dep(task385);
  task385->add_dep(task283);
  densityq->add_task(task385);

  auto task386 = make_shared<Task386>(I441, Gamma12_(), t2);
  task384->add_dep(task386);
  task386->add_dep(task283);
  densityq->add_task(task386);

  auto I443 = make_shared<TATensor<double,2>>({active_, closed_});
  auto task387 = make_shared<Task387>(den2, I443);
  task387->add_dep(task283);
  densityq->add_task(task387);

  auto I444 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task388 = make_shared<Task388>(I443, t2, I444);
  task387->add_dep(task388);
  task388->add_dep(task283);
  densityq->add_task(task388);

  auto I445 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task389 = make_shared<Task389>(I444, Gamma12_(), I445);
  task388->add_dep(task389);
  task389->add_dep(task283);
  densityq->add_task(task389);

  auto task390 = make_shared<Task390>(I445, t2);
  task389->add_dep(task390);
  task390->add_dep(task283);
  densityq->add_task(task390);

  auto I452 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task391 = make_shared<Task391>(den2, I452);
  task391->add_dep(task283);
  densityq->add_task(task391);

  auto I453 = make_shared<TATensor<double,2>>({active_, virt_});
  auto task392 = make_shared<Task392>(I452, Gamma16_(), I453);
  task391->add_dep(task392);
  task392->add_dep(task283);
  densityq->add_task(task392);

  auto task393 = make_shared<Task393>(I453, t2);
  task392->add_dep(task393);
  task393->add_dep(task283);
  densityq->add_task(task393);

  auto task394 = make_shared<Task394>(I453, t2);
  task392->add_dep(task394);
  task394->add_dep(task283);
  densityq->add_task(task394);

  auto I458 = make_shared<TATensor<double,2>>({active_, virt_});
  auto task395 = make_shared<Task395>(den2, I458);
  task395->add_dep(task283);
  densityq->add_task(task395);

  auto I459 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task396 = make_shared<Task396>(I458, t2, I459);
  task395->add_dep(task396);
  task396->add_dep(task283);
  densityq->add_task(task396);

  auto task397 = make_shared<Task397>(I459, Gamma28_(), t2);
  task396->add_dep(task397);
  task397->add_dep(task283);
  densityq->add_task(task397);

  auto I461 = make_shared<TATensor<double,2>>({active_, closed_});
  auto task398 = make_shared<Task398>(den2, I461);
  task398->add_dep(task283);
  densityq->add_task(task398);

  auto I462 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task399 = make_shared<Task399>(I461, t2, I462);
  task398->add_dep(task399);
  task399->add_dep(task283);
  densityq->add_task(task399);

  auto task400 = make_shared<Task400>(I462, Gamma29_(), t2);
  task399->add_dep(task400);
  task400->add_dep(task283);
  densityq->add_task(task400);

  auto I465 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task401 = make_shared<Task401>(I461, t2, I465);
  task398->add_dep(task401);
  task401->add_dep(task283);
  densityq->add_task(task401);

  auto task402 = make_shared<Task402>(I465, Gamma7_(), t2);
  task401->add_dep(task402);
  task402->add_dep(task283);
  densityq->add_task(task402);

  auto I504 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task403 = make_shared<Task403>(I461, t2, I504);
  task398->add_dep(task403);
  task403->add_dep(task283);
  densityq->add_task(task403);

  auto task404 = make_shared<Task404>(I504, Gamma7_(), t2);
  task403->add_dep(task404);
  task404->add_dep(task283);
  densityq->add_task(task404);

  auto I507 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task405 = make_shared<Task405>(I461, t2, I507);
  task398->add_dep(task405);
  task405->add_dep(task283);
  densityq->add_task(task405);

  auto task406 = make_shared<Task406>(I507, Gamma7_(), t2);
  task405->add_dep(task406);
  task406->add_dep(task283);
  densityq->add_task(task406);

  auto I656 = make_shared<TATensor<double,4>>({virt_, virt_, active_, active_});
  auto task407 = make_shared<Task407>(I461, t2, I656);
  task398->add_dep(task407);
  task407->add_dep(task283);
  densityq->add_task(task407);

  auto task408 = make_shared<Task408>(I656, Gamma60_(), t2);
  task407->add_dep(task408);
  task408->add_dep(task283);
  densityq->add_task(task408);

  auto I473 = make_shared<TATensor<double,2>>({virt_, virt_});
  auto task409 = make_shared<Task409>(den2, I473);
  task409->add_dep(task283);
  densityq->add_task(task409);

  auto I474 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task410 = make_shared<Task410>(I473, t2, I474);
  task409->add_dep(task410);
  task410->add_dep(task283);
  densityq->add_task(task410);

  auto task411 = make_shared<Task411>(I474, Gamma32_(), t2);
  task410->add_dep(task411);
  task411->add_dep(task283);
  densityq->add_task(task411);

  auto I483 = make_shared<TATensor<double,4>>({closed_, virt_, active_, active_});
  auto task412 = make_shared<Task412>(I473, t2, I483);
  task409->add_dep(task412);
  task412->add_dep(task283);
  densityq->add_task(task412);

  auto task413 = make_shared<Task413>(I483, Gamma35_(), t2);
  task412->add_dep(task413);
  task413->add_dep(task283);
  densityq->add_task(task413);

  auto I488 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task414 = make_shared<Task414>(den2, I488);
  task414->add_dep(task283);
  densityq->add_task(task414);

  auto I489 = make_shared<TATensor<double,2>>({closed_, virt_});
  auto task415 = make_shared<Task415>(I488, t2, I489);
  task414->add_dep(task415);
  task415->add_dep(task283);
  densityq->add_task(task415);

  auto task416 = make_shared<Task416>(I489, Gamma38_(), t2);
  task415->add_dep(task416);
  task416->add_dep(task283);
  densityq->add_task(task416);

  auto I492 = make_shared<TATensor<double,2>>({closed_, virt_});
  auto task417 = make_shared<Task417>(I488, t2, I492);
  task414->add_dep(task417);
  task417->add_dep(task283);
  densityq->add_task(task417);

  auto task418 = make_shared<Task418>(I492, Gamma38_(), t2);
  task417->add_dep(task418);
  task418->add_dep(task283);
  densityq->add_task(task418);

  auto I531 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task419 = make_shared<Task419>(I488, t2, I531);
  task414->add_dep(task419);
  task419->add_dep(task283);
  densityq->add_task(task419);

  auto task420 = make_shared<Task420>(I531, Gamma38_(), t2);
  task419->add_dep(task420);
  task420->add_dep(task283);
  densityq->add_task(task420);

  auto I534 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task421 = make_shared<Task421>(I488, t2, I534);
  task414->add_dep(task421);
  task421->add_dep(task283);
  densityq->add_task(task421);

  auto task422 = make_shared<Task422>(I534, Gamma38_(), t2);
  task421->add_dep(task422);
  task422->add_dep(task283);
  densityq->add_task(task422);

  auto I500 = make_shared<TATensor<double,2>>({active_, virt_});
  auto task423 = make_shared<Task423>(den2, I500);
  task423->add_dep(task283);
  densityq->add_task(task423);

  auto I501 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task424 = make_shared<Task424>(I500, t2, I501);
  task423->add_dep(task424);
  task424->add_dep(task283);
  densityq->add_task(task424);

  auto task425 = make_shared<Task425>(I501, Gamma6_(), t2);
  task424->add_dep(task425);
  task425->add_dep(task283);
  densityq->add_task(task425);

  auto I653 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task426 = make_shared<Task426>(I500, t2, I653);
  task423->add_dep(task426);
  task426->add_dep(task283);
  densityq->add_task(task426);

  auto task427 = make_shared<Task427>(I653, Gamma59_(), t2);
  task426->add_dep(task427);
  task427->add_dep(task283);
  densityq->add_task(task427);

  auto I512 = make_shared<TATensor<double,2>>({closed_, closed_});
  auto task428 = make_shared<Task428>(den2, I512);
  task428->add_dep(task283);
  densityq->add_task(task428);

  auto I513 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task429 = make_shared<Task429>(I512, t2, I513);
  task428->add_dep(task429);
  task429->add_dep(task283);
  densityq->add_task(task429);

  auto task430 = make_shared<Task430>(I513, Gamma35_(), t2);
  task429->add_dep(task430);
  task430->add_dep(task283);
  densityq->add_task(task430);

  auto I522 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task431 = make_shared<Task431>(I512, t2, I522);
  task428->add_dep(task431);
  task431->add_dep(task283);
  densityq->add_task(task431);

  auto task432 = make_shared<Task432>(I522, Gamma35_(), t2);
  task431->add_dep(task432);
  task432->add_dep(task283);
  densityq->add_task(task432);

  auto I515 = make_shared<TATensor<double,2>>({virt_, virt_});
  auto task433 = make_shared<Task433>(den2, I515);
  task433->add_dep(task283);
  densityq->add_task(task433);

  auto I516 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task434 = make_shared<Task434>(I515, t2, I516);
  task433->add_dep(task434);
  task434->add_dep(task283);
  densityq->add_task(task434);

  auto task435 = make_shared<Task435>(I516, Gamma35_(), t2);
  task434->add_dep(task435);
  task435->add_dep(task283);
  densityq->add_task(task435);

  auto I525 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task436 = make_shared<Task436>(I515, t2, I525);
  task433->add_dep(task436);
  task436->add_dep(task283);
  densityq->add_task(task436);

  auto task437 = make_shared<Task437>(I525, Gamma35_(), t2);
  task436->add_dep(task437);
  task437->add_dep(task283);
  densityq->add_task(task437);

  auto I662 = make_shared<TATensor<double,4>>({virt_, virt_, active_, active_});
  auto task438 = make_shared<Task438>(I515, t2, I662);
  task433->add_dep(task438);
  task438->add_dep(task283);
  densityq->add_task(task438);

  auto task439 = make_shared<Task439>(I662, Gamma60_(), t2);
  task438->add_dep(task439);
  task439->add_dep(task283);
  densityq->add_task(task439);

  auto I527 = make_shared<TATensor<double,2>>({active_, closed_});
  auto task440 = make_shared<Task440>(den2, I527);
  task440->add_dep(task283);
  densityq->add_task(task440);

  auto I528 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task441 = make_shared<Task441>(I527, t2, I528);
  task440->add_dep(task441);
  task441->add_dep(task283);
  densityq->add_task(task441);

  auto task442 = make_shared<Task442>(I528, Gamma51_(), t2);
  task441->add_dep(task442);
  task442->add_dep(task283);
  densityq->add_task(task442);

  auto I551 = make_shared<TATensor<double,2>>({virt_, virt_});
  auto task443 = make_shared<Task443>(den2, I551);
  task443->add_dep(task283);
  densityq->add_task(task443);

  auto I552 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task444 = make_shared<Task444>(I551, t2, I552);
  task443->add_dep(task444);
  task444->add_dep(task283);
  densityq->add_task(task444);

  auto task445 = make_shared<Task445>(I552, Gamma59_(), t2);
  task444->add_dep(task445);
  task445->add_dep(task283);
  densityq->add_task(task445);

  auto I563 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task446 = make_shared<Task446>(den2, I563);
  task446->add_dep(task283);
  densityq->add_task(task446);

  auto I564 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task447 = make_shared<Task447>(I563, Gamma16_(), I564);
  task446->add_dep(task447);
  task447->add_dep(task283);
  densityq->add_task(task447);

  auto task448 = make_shared<Task448>(I564, t2);
  task447->add_dep(task448);
  task448->add_dep(task283);
  densityq->add_task(task448);

  auto I566 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task449 = make_shared<Task449>(den2, I566);
  task449->add_dep(task283);
  densityq->add_task(task449);

  auto I567 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task450 = make_shared<Task450>(I566, Gamma16_(), I567);
  task449->add_dep(task450);
  task450->add_dep(task283);
  densityq->add_task(task450);

  auto task451 = make_shared<Task451>(I567, t2);
  task450->add_dep(task451);
  task451->add_dep(task283);
  densityq->add_task(task451);

  auto I572 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task452 = make_shared<Task452>(den2, I572);
  task452->add_dep(task283);
  densityq->add_task(task452);

  auto I573 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task453 = make_shared<Task453>(I572, t2, I573);
  task452->add_dep(task453);
  task453->add_dep(task283);
  densityq->add_task(task453);

  auto I574 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task454 = make_shared<Task454>(I573, Gamma38_(), I574);
  task453->add_dep(task454);
  task454->add_dep(task283);
  densityq->add_task(task454);

  auto task455 = make_shared<Task455>(I574, t2);
  task454->add_dep(task455);
  task455->add_dep(task283);
  densityq->add_task(task455);

  auto I581 = make_shared<TATensor<double,2>>({active_, active_});
  auto task456 = make_shared<Task456>(den2, I581);
  task456->add_dep(task283);
  densityq->add_task(task456);

  auto I582 = make_shared<TATensor<double,0>>({});
  auto task457 = make_shared<Task457>(I581, Gamma38_(), I582);
  task456->add_dep(task457);
  task457->add_dep(task283);
  densityq->add_task(task457);

  auto task458 = make_shared<Task458>(I582, t2);
  task457->add_dep(task458);
  task458->add_dep(task283);
  densityq->add_task(task458);

  auto task459 = make_shared<Task459>(I582, t2);
  task457->add_dep(task459);
  task459->add_dep(task283);
  densityq->add_task(task459);

  shared_ptr<TATensor<double,2>> I587;
  if (diagonal) {
    I587 = make_shared<TATensor<double,2>>({closed_, closed_});
  }
  shared_ptr<Task460> task460;
  if (diagonal) {
    task460 = make_shared<Task460>(den2, I587);
    task460->add_dep(task283);
    densityq->add_task(task460);
  }

  shared_ptr<Task461> task461;
  if (diagonal) {
    task461 = make_shared<Task461>(I587, t2);
    task460->add_dep(task461);
    task461->add_dep(task283);
    densityq->add_task(task461);
  }

  shared_ptr<Task462> task462;
  if (diagonal) {
    task462 = make_shared<Task462>(I587, t2);
    task460->add_dep(task462);
    task462->add_dep(task283);
    densityq->add_task(task462);
  }

  shared_ptr<TATensor<double,2>> I591;
  if (diagonal) {
    I591 = make_shared<TATensor<double,2>>({virt_, virt_});
  }
  shared_ptr<Task463> task463;
  if (diagonal) {
    task463 = make_shared<Task463>(den2, I591);
    task463->add_dep(task283);
    densityq->add_task(task463);
  }

  shared_ptr<Task464> task464;
  if (diagonal) {
    task464 = make_shared<Task464>(I591, t2);
    task463->add_dep(task464);
    task464->add_dep(task283);
    densityq->add_task(task464);
  }

  shared_ptr<Task465> task465;
  if (diagonal) {
    task465 = make_shared<Task465>(I591, t2);
    task463->add_dep(task465);
    task465->add_dep(task283);
    densityq->add_task(task465);
  }

  auto I595 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task466 = make_shared<Task466>(den2, I595);
  task466->add_dep(task283);
  densityq->add_task(task466);

  auto I596 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task467 = make_shared<Task467>(I595, Gamma38_(), I596);
  task466->add_dep(task467);
  task467->add_dep(task283);
  densityq->add_task(task467);

  auto task468 = make_shared<Task468>(I596, t2);
  task467->add_dep(task468);
  task468->add_dep(task283);
  densityq->add_task(task468);

  auto task469 = make_shared<Task469>(I596, t2);
  task467->add_dep(task469);
  task469->add_dep(task283);
  densityq->add_task(task469);

  auto I601 = make_shared<TATensor<double,2>>({active_, virt_});
  auto task470 = make_shared<Task470>(den2, I601);
  task470->add_dep(task283);
  densityq->add_task(task470);

  auto I602 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task471 = make_shared<Task471>(I601, t2, I602);
  task470->add_dep(task471);
  task471->add_dep(task283);
  densityq->add_task(task471);

  auto I603 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task472 = make_shared<Task472>(I602, Gamma35_(), I603);
  task471->add_dep(task472);
  task472->add_dep(task283);
  densityq->add_task(task472);

  auto task473 = make_shared<Task473>(I603, t2);
  task472->add_dep(task473);
  task473->add_dep(task283);
  densityq->add_task(task473);

  auto I604 = make_shared<TATensor<double,2>>({active_, virt_});
  auto task474 = make_shared<Task474>(den2, I604);
  task474->add_dep(task283);
  densityq->add_task(task474);

  auto I605 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task475 = make_shared<Task475>(I604, t2, I605);
  task474->add_dep(task475);
  task475->add_dep(task283);
  densityq->add_task(task475);

  auto task476 = make_shared<Task476>(I605, Gamma32_(), t2);
  task475->add_dep(task476);
  task476->add_dep(task283);
  densityq->add_task(task476);

  auto task477 = make_shared<Task477>(I605, Gamma35_(), t2);
  task475->add_dep(task477);
  task477->add_dep(task283);
  densityq->add_task(task477);

  auto I613 = make_shared<TATensor<double,2>>({closed_, virt_});
  auto task478 = make_shared<Task478>(den2, I613);
  task478->add_dep(task283);
  densityq->add_task(task478);

  auto I614 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task479 = make_shared<Task479>(I613, t2, I614);
  task478->add_dep(task479);
  task479->add_dep(task283);
  densityq->add_task(task479);

  auto task480 = make_shared<Task480>(I614, Gamma60_(), t2);
  task479->add_dep(task480);
  task480->add_dep(task283);
  densityq->add_task(task480);

  auto I616 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task481 = make_shared<Task481>(den2, I616);
  task481->add_dep(task283);
  densityq->add_task(task481);

  auto I617 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task482 = make_shared<Task482>(I616, t2, I617);
  task481->add_dep(task482);
  task482->add_dep(task283);
  densityq->add_task(task482);

  auto task483 = make_shared<Task483>(I617, Gamma60_(), t2);
  task482->add_dep(task483);
  task483->add_dep(task283);
  densityq->add_task(task483);

  auto I619 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task484 = make_shared<Task484>(den2, I619);
  task484->add_dep(task283);
  densityq->add_task(task484);

  auto I620 = make_shared<TATensor<double,2>>({active_, closed_});
  auto task485 = make_shared<Task485>(I619, Gamma38_(), I620);
  task484->add_dep(task485);
  task485->add_dep(task283);
  densityq->add_task(task485);

  auto task486 = make_shared<Task486>(I620, t2);
  task485->add_dep(task486);
  task486->add_dep(task283);
  densityq->add_task(task486);

  auto task487 = make_shared<Task487>(I620, t2);
  task485->add_dep(task487);
  task487->add_dep(task283);
  densityq->add_task(task487);

  auto I631 = make_shared<TATensor<double,2>>({closed_, closed_});
  auto task488 = make_shared<Task488>(den2, I631);
  task488->add_dep(task283);
  densityq->add_task(task488);

  auto I632 = make_shared<TATensor<double,4>>({virt_, closed_, virt_, active_});
  auto task489 = make_shared<Task489>(I631, t2, I632);
  task488->add_dep(task489);
  task489->add_dep(task283);
  densityq->add_task(task489);

  auto task490 = make_shared<Task490>(I632, Gamma38_(), t2);
  task489->add_dep(task490);
  task490->add_dep(task283);
  densityq->add_task(task490);

  auto I635 = make_shared<TATensor<double,4>>({virt_, closed_, virt_, active_});
  auto task491 = make_shared<Task491>(I631, t2, I635);
  task488->add_dep(task491);
  task491->add_dep(task283);
  densityq->add_task(task491);

  auto task492 = make_shared<Task492>(I635, Gamma38_(), t2);
  task491->add_dep(task492);
  task492->add_dep(task283);
  densityq->add_task(task492);

  auto I637 = make_shared<TATensor<double,2>>({virt_, virt_});
  auto task493 = make_shared<Task493>(den2, I637);
  task493->add_dep(task283);
  densityq->add_task(task493);

  auto I638 = make_shared<TATensor<double,4>>({virt_, closed_, virt_, active_});
  auto task494 = make_shared<Task494>(I637, t2, I638);
  task493->add_dep(task494);
  task494->add_dep(task283);
  densityq->add_task(task494);

  auto task495 = make_shared<Task495>(I638, Gamma38_(), t2);
  task494->add_dep(task495);
  task495->add_dep(task283);
  densityq->add_task(task495);

  auto I641 = make_shared<TATensor<double,4>>({virt_, closed_, virt_, active_});
  auto task496 = make_shared<Task496>(I637, t2, I641);
  task493->add_dep(task496);
  task496->add_dep(task283);
  densityq->add_task(task496);

  auto task497 = make_shared<Task497>(I641, Gamma38_(), t2);
  task496->add_dep(task497);
  task497->add_dep(task283);
  densityq->add_task(task497);

  auto I643 = make_shared<TATensor<double,2>>({virt_, virt_});
  auto task498 = make_shared<Task498>(den2, I643);
  task498->add_dep(task283);
  densityq->add_task(task498);

  auto I644 = make_shared<TATensor<double,4>>({virt_, closed_, virt_, active_});
  auto task499 = make_shared<Task499>(I643, t2, I644);
  task498->add_dep(task499);
  task499->add_dep(task283);
  densityq->add_task(task499);

  auto task500 = make_shared<Task500>(I644, Gamma38_(), t2);
  task499->add_dep(task500);
  task500->add_dep(task283);
  densityq->add_task(task500);

  auto I647 = make_shared<TATensor<double,4>>({virt_, closed_, virt_, active_});
  auto task501 = make_shared<Task501>(I643, t2, I647);
  task498->add_dep(task501);
  task501->add_dep(task283);
  densityq->add_task(task501);

  auto task502 = make_shared<Task502>(I647, Gamma38_(), t2);
  task501->add_dep(task502);
  task502->add_dep(task283);
  densityq->add_task(task502);

  auto I649 = make_shared<TATensor<double,2>>({active_, closed_});
  auto task503 = make_shared<Task503>(den2, I649);
  task503->add_dep(task283);
  densityq->add_task(task503);

  auto I650 = make_shared<TATensor<double,4>>({virt_, virt_, active_, active_});
  auto task504 = make_shared<Task504>(I649, t2, I650);
  task503->add_dep(task504);
  task504->add_dep(task283);
  densityq->add_task(task504);

  auto task505 = make_shared<Task505>(I650, Gamma60_(), t2);
  task504->add_dep(task505);
  task505->add_dep(task283);
  densityq->add_task(task505);

  return densityq;
}


#endif
