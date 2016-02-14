//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_densityq2.cc
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
#include <src/smith/caspt2/CASPT2_tasks7.h>
#include <src/smith/caspt2/CASPT2_tasks8.h>
#include <src/smith/caspt2/CASPT2_tasks9.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void  CASPT2::CASPT2::make_densityq2(shared_ptr<Queue> densityq, shared_ptr<Task> task284, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I401_index = {active_, virt_};
  auto I401 = make_shared<Tensor>(I401_index);
  auto tensor343 = vector<shared_ptr<Tensor>>{den2, I401};
  auto task343 = make_shared<Task343>(tensor343, pindex);
  task343->add_dep(task284);
  densityq->add_task(task343);

  vector<IndexRange> I402_index = {closed_, active_, active_, active_};
  auto I402 = make_shared<Tensor>(I402_index);
  auto tensor344 = vector<shared_ptr<Tensor>>{I401, t2, I402};
  auto task344 = make_shared<Task344>(tensor344, pindex);
  task343->add_dep(task344);
  task344->add_dep(task284);
  densityq->add_task(task344);

  auto tensor345 = vector<shared_ptr<Tensor>>{I402, Gamma9_(), t2};
  auto task345 = make_shared<Task345>(tensor345, pindex);
  task344->add_dep(task345);
  task345->add_dep(task284);
  densityq->add_task(task345);

  vector<IndexRange> I405_index = {closed_, active_, active_, active_};
  auto I405 = make_shared<Tensor>(I405_index);
  auto tensor346 = vector<shared_ptr<Tensor>>{I401, t2, I405};
  auto task346 = make_shared<Task346>(tensor346, pindex);
  task343->add_dep(task346);
  task346->add_dep(task284);
  densityq->add_task(task346);

  auto tensor347 = vector<shared_ptr<Tensor>>{I405, Gamma6_(), t2};
  auto task347 = make_shared<Task347>(tensor347, pindex);
  task346->add_dep(task347);
  task347->add_dep(task284);
  densityq->add_task(task347);

  vector<IndexRange> I561_index = {virt_, active_, active_, active_};
  auto I561 = make_shared<Tensor>(I561_index);
  auto tensor348 = vector<shared_ptr<Tensor>>{I401, t2, I561};
  auto task348 = make_shared<Task348>(tensor348, pindex);
  task343->add_dep(task348);
  task348->add_dep(task284);
  densityq->add_task(task348);

  auto tensor349 = vector<shared_ptr<Tensor>>{I561, Gamma59_(), t2};
  auto task349 = make_shared<Task349>(tensor349, pindex);
  task348->add_dep(task349);
  task349->add_dep(task284);
  densityq->add_task(task349);

  vector<IndexRange> I407_index = {active_, virt_};
  auto I407 = make_shared<Tensor>(I407_index);
  auto tensor350 = vector<shared_ptr<Tensor>>{den2, I407};
  auto task350 = make_shared<Task350>(tensor350, pindex);
  task350->add_dep(task284);
  densityq->add_task(task350);

  vector<IndexRange> I408_index = {closed_, closed_, active_, active_};
  auto I408 = make_shared<Tensor>(I408_index);
  auto tensor351 = vector<shared_ptr<Tensor>>{I407, t2, I408};
  auto task351 = make_shared<Task351>(tensor351, pindex);
  task350->add_dep(task351);
  task351->add_dep(task284);
  densityq->add_task(task351);

  auto tensor352 = vector<shared_ptr<Tensor>>{I408, Gamma3_(), t2};
  auto task352 = make_shared<Task352>(tensor352, pindex);
  task351->add_dep(task352);
  task352->add_dep(task284);
  densityq->add_task(task352);

  vector<IndexRange> I410_index = {virt_, closed_};
  auto I410 = make_shared<Tensor>(I410_index);
  auto tensor353 = vector<shared_ptr<Tensor>>{den2, I410};
  auto task353 = make_shared<Task353>(tensor353, pindex);
  task353->add_dep(task284);
  densityq->add_task(task353);

  vector<IndexRange> I411_index = {closed_, active_};
  auto I411 = make_shared<Tensor>(I411_index);
  auto tensor354 = vector<shared_ptr<Tensor>>{I410, t2, I411};
  auto task354 = make_shared<Task354>(tensor354, pindex);
  task353->add_dep(task354);
  task354->add_dep(task284);
  densityq->add_task(task354);

  auto tensor355 = vector<shared_ptr<Tensor>>{I411, Gamma12_(), t2};
  auto task355 = make_shared<Task355>(tensor355, pindex);
  task354->add_dep(task355);
  task355->add_dep(task284);
  densityq->add_task(task355);

  vector<IndexRange> I413_index = {closed_, virt_};
  auto I413 = make_shared<Tensor>(I413_index);
  auto tensor356 = vector<shared_ptr<Tensor>>{den2, I413};
  auto task356 = make_shared<Task356>(tensor356, pindex);
  task356->add_dep(task284);
  densityq->add_task(task356);

  vector<IndexRange> I414_index = {closed_, active_};
  auto I414 = make_shared<Tensor>(I414_index);
  auto tensor357 = vector<shared_ptr<Tensor>>{I413, t2, I414};
  auto task357 = make_shared<Task357>(tensor357, pindex);
  task356->add_dep(task357);
  task357->add_dep(task284);
  densityq->add_task(task357);

  auto tensor358 = vector<shared_ptr<Tensor>>{I414, Gamma12_(), t2};
  auto task358 = make_shared<Task358>(tensor358, pindex);
  task357->add_dep(task358);
  task358->add_dep(task284);
  densityq->add_task(task358);

  vector<IndexRange> I570_index = {virt_, closed_};
  auto I570 = make_shared<Tensor>(I570_index);
  auto tensor359 = vector<shared_ptr<Tensor>>{I413, t2, I570};
  auto task359 = make_shared<Task359>(tensor359, pindex);
  task356->add_dep(task359);
  task359->add_dep(task284);
  densityq->add_task(task359);

  vector<IndexRange> I571_index = {active_, virt_, closed_, active_};
  auto I571 = make_shared<Tensor>(I571_index);
  auto tensor360 = vector<shared_ptr<Tensor>>{I570, Gamma38_(), I571};
  auto task360 = make_shared<Task360>(tensor360, pindex);
  task359->add_dep(task360);
  task360->add_dep(task284);
  densityq->add_task(task360);

  auto tensor361 = vector<shared_ptr<Tensor>>{I571, t2};
  auto task361 = make_shared<Task361>(tensor361, pindex);
  task360->add_dep(task361);
  task361->add_dep(task284);
  densityq->add_task(task361);

  vector<IndexRange> I416_index = {active_, active_};
  auto I416 = make_shared<Tensor>(I416_index);
  auto tensor362 = vector<shared_ptr<Tensor>>{den2, I416};
  auto task362 = make_shared<Task362>(tensor362, pindex);
  task362->add_dep(task284);
  densityq->add_task(task362);

  vector<IndexRange> I417_index = {active_, active_};
  auto I417 = make_shared<Tensor>(I417_index);
  auto tensor363 = vector<shared_ptr<Tensor>>{I416, Gamma152_(), I417};
  auto task363 = make_shared<Task363>(tensor363, pindex);
  task362->add_dep(task363);
  task363->add_dep(task284);
  densityq->add_task(task363);

  auto tensor364 = vector<shared_ptr<Tensor>>{I417, t2};
  auto task364 = make_shared<Task364>(tensor364, pindex);
  task363->add_dep(task364);
  task364->add_dep(task284);
  densityq->add_task(task364);

  auto tensor365 = vector<shared_ptr<Tensor>>{I417, t2};
  auto task365 = make_shared<Task365>(tensor365, pindex);
  task363->add_dep(task365);
  task365->add_dep(task284);
  densityq->add_task(task365);

  vector<IndexRange> I626_index = {active_, active_};
  auto I626 = make_shared<Tensor>(I626_index);
  auto tensor366 = vector<shared_ptr<Tensor>>{I416, Gamma60_(), I626};
  auto task366 = make_shared<Task366>(tensor366, pindex);
  task362->add_dep(task366);
  task366->add_dep(task284);
  densityq->add_task(task366);

  auto tensor367 = vector<shared_ptr<Tensor>>{I626, t2};
  auto task367 = make_shared<Task367>(tensor367, pindex);
  task366->add_dep(task367);
  task367->add_dep(task284);
  densityq->add_task(task367);

  auto tensor368 = vector<shared_ptr<Tensor>>{I626, t2};
  auto task368 = make_shared<Task368>(tensor368, pindex);
  task366->add_dep(task368);
  task368->add_dep(task284);
  densityq->add_task(task368);

  vector<IndexRange> I422_index = {closed_, closed_};
  auto I422 = make_shared<Tensor>(I422_index);
  auto tensor369 = vector<shared_ptr<Tensor>>{den2, I422};
  auto task369 = make_shared<Task369>(tensor369, pindex);
  task369->add_dep(task284);
  densityq->add_task(task369);

  vector<IndexRange> I423_index = {closed_, virt_, closed_, active_};
  auto I423 = make_shared<Tensor>(I423_index);
  auto tensor370 = vector<shared_ptr<Tensor>>{I422, t2, I423};
  auto task370 = make_shared<Task370>(tensor370, pindex);
  task369->add_dep(task370);
  task370->add_dep(task284);
  densityq->add_task(task370);

  auto tensor371 = vector<shared_ptr<Tensor>>{I423, Gamma16_(), t2};
  auto task371 = make_shared<Task371>(tensor371, pindex);
  task370->add_dep(task371);
  task371->add_dep(task284);
  densityq->add_task(task371);

  vector<IndexRange> I426_index = {closed_, virt_, closed_, active_};
  auto I426 = make_shared<Tensor>(I426_index);
  auto tensor372 = vector<shared_ptr<Tensor>>{I422, t2, I426};
  auto task372 = make_shared<Task372>(tensor372, pindex);
  task369->add_dep(task372);
  task372->add_dep(task284);
  densityq->add_task(task372);

  auto tensor373 = vector<shared_ptr<Tensor>>{I426, Gamma16_(), t2};
  auto task373 = make_shared<Task373>(tensor373, pindex);
  task372->add_dep(task373);
  task373->add_dep(task284);
  densityq->add_task(task373);

  vector<IndexRange> I428_index = {closed_, closed_};
  auto I428 = make_shared<Tensor>(I428_index);
  auto tensor374 = vector<shared_ptr<Tensor>>{den2, I428};
  auto task374 = make_shared<Task374>(tensor374, pindex);
  task374->add_dep(task284);
  densityq->add_task(task374);

  vector<IndexRange> I429_index = {closed_, virt_, closed_, active_};
  auto I429 = make_shared<Tensor>(I429_index);
  auto tensor375 = vector<shared_ptr<Tensor>>{I428, t2, I429};
  auto task375 = make_shared<Task375>(tensor375, pindex);
  task374->add_dep(task375);
  task375->add_dep(task284);
  densityq->add_task(task375);

  auto tensor376 = vector<shared_ptr<Tensor>>{I429, Gamma16_(), t2};
  auto task376 = make_shared<Task376>(tensor376, pindex);
  task375->add_dep(task376);
  task376->add_dep(task284);
  densityq->add_task(task376);

  vector<IndexRange> I435_index = {closed_, virt_, closed_, active_};
  auto I435 = make_shared<Tensor>(I435_index);
  auto tensor377 = vector<shared_ptr<Tensor>>{I428, t2, I435};
  auto task377 = make_shared<Task377>(tensor377, pindex);
  task374->add_dep(task377);
  task377->add_dep(task284);
  densityq->add_task(task377);

  auto tensor378 = vector<shared_ptr<Tensor>>{I435, Gamma16_(), t2};
  auto task378 = make_shared<Task378>(tensor378, pindex);
  task377->add_dep(task378);
  task378->add_dep(task284);
  densityq->add_task(task378);

  vector<IndexRange> I431_index = {virt_, virt_};
  auto I431 = make_shared<Tensor>(I431_index);
  auto tensor379 = vector<shared_ptr<Tensor>>{den2, I431};
  auto task379 = make_shared<Task379>(tensor379, pindex);
  task379->add_dep(task284);
  densityq->add_task(task379);

  vector<IndexRange> I432_index = {closed_, virt_, closed_, active_};
  auto I432 = make_shared<Tensor>(I432_index);
  auto tensor380 = vector<shared_ptr<Tensor>>{I431, t2, I432};
  auto task380 = make_shared<Task380>(tensor380, pindex);
  task379->add_dep(task380);
  task380->add_dep(task284);
  densityq->add_task(task380);

  auto tensor381 = vector<shared_ptr<Tensor>>{I432, Gamma16_(), t2};
  auto task381 = make_shared<Task381>(tensor381, pindex);
  task380->add_dep(task381);
  task381->add_dep(task284);
  densityq->add_task(task381);

  vector<IndexRange> I438_index = {closed_, virt_, closed_, active_};
  auto I438 = make_shared<Tensor>(I438_index);
  auto tensor382 = vector<shared_ptr<Tensor>>{I431, t2, I438};
  auto task382 = make_shared<Task382>(tensor382, pindex);
  task379->add_dep(task382);
  task382->add_dep(task284);
  densityq->add_task(task382);

  auto tensor383 = vector<shared_ptr<Tensor>>{I438, Gamma16_(), t2};
  auto task383 = make_shared<Task383>(tensor383, pindex);
  task382->add_dep(task383);
  task383->add_dep(task284);
  densityq->add_task(task383);

  vector<IndexRange> I440_index = {active_, closed_};
  auto I440 = make_shared<Tensor>(I440_index);
  auto tensor384 = vector<shared_ptr<Tensor>>{den2, I440};
  auto task384 = make_shared<Task384>(tensor384, pindex);
  task384->add_dep(task284);
  densityq->add_task(task384);

  vector<IndexRange> I441_index = {virt_, closed_, active_, active_};
  auto I441 = make_shared<Tensor>(I441_index);
  auto tensor385 = vector<shared_ptr<Tensor>>{I440, t2, I441};
  auto task385 = make_shared<Task385>(tensor385, pindex);
  task384->add_dep(task385);
  task385->add_dep(task284);
  densityq->add_task(task385);

  auto tensor386 = vector<shared_ptr<Tensor>>{I441, Gamma22_(), t2};
  auto task386 = make_shared<Task386>(tensor386, pindex);
  task385->add_dep(task386);
  task386->add_dep(task284);
  densityq->add_task(task386);

  auto tensor387 = vector<shared_ptr<Tensor>>{I441, Gamma12_(), t2};
  auto task387 = make_shared<Task387>(tensor387, pindex);
  task385->add_dep(task387);
  task387->add_dep(task284);
  densityq->add_task(task387);

  vector<IndexRange> I443_index = {active_, closed_};
  auto I443 = make_shared<Tensor>(I443_index);
  auto tensor388 = vector<shared_ptr<Tensor>>{den2, I443};
  auto task388 = make_shared<Task388>(tensor388, pindex);
  task388->add_dep(task284);
  densityq->add_task(task388);

  vector<IndexRange> I444_index = {virt_, closed_, active_, active_};
  auto I444 = make_shared<Tensor>(I444_index);
  auto tensor389 = vector<shared_ptr<Tensor>>{I443, t2, I444};
  auto task389 = make_shared<Task389>(tensor389, pindex);
  task388->add_dep(task389);
  task389->add_dep(task284);
  densityq->add_task(task389);

  vector<IndexRange> I445_index = {active_, virt_, closed_, active_};
  auto I445 = make_shared<Tensor>(I445_index);
  auto tensor390 = vector<shared_ptr<Tensor>>{I444, Gamma12_(), I445};
  auto task390 = make_shared<Task390>(tensor390, pindex);
  task389->add_dep(task390);
  task390->add_dep(task284);
  densityq->add_task(task390);

  auto tensor391 = vector<shared_ptr<Tensor>>{I445, t2};
  auto task391 = make_shared<Task391>(tensor391, pindex);
  task390->add_dep(task391);
  task391->add_dep(task284);
  densityq->add_task(task391);

  vector<IndexRange> I452_index = {virt_, active_};
  auto I452 = make_shared<Tensor>(I452_index);
  auto tensor392 = vector<shared_ptr<Tensor>>{den2, I452};
  auto task392 = make_shared<Task392>(tensor392, pindex);
  task392->add_dep(task284);
  densityq->add_task(task392);

  vector<IndexRange> I453_index = {active_, virt_};
  auto I453 = make_shared<Tensor>(I453_index);
  auto tensor393 = vector<shared_ptr<Tensor>>{I452, Gamma16_(), I453};
  auto task393 = make_shared<Task393>(tensor393, pindex);
  task392->add_dep(task393);
  task393->add_dep(task284);
  densityq->add_task(task393);

  auto tensor394 = vector<shared_ptr<Tensor>>{I453, t2};
  auto task394 = make_shared<Task394>(tensor394, pindex);
  task393->add_dep(task394);
  task394->add_dep(task284);
  densityq->add_task(task394);

  auto tensor395 = vector<shared_ptr<Tensor>>{I453, t2};
  auto task395 = make_shared<Task395>(tensor395, pindex);
  task393->add_dep(task395);
  task395->add_dep(task284);
  densityq->add_task(task395);

  vector<IndexRange> I458_index = {active_, virt_};
  auto I458 = make_shared<Tensor>(I458_index);
  auto tensor396 = vector<shared_ptr<Tensor>>{den2, I458};
  auto task396 = make_shared<Task396>(tensor396, pindex);
  task396->add_dep(task284);
  densityq->add_task(task396);

  vector<IndexRange> I459_index = {closed_, active_, active_, active_};
  auto I459 = make_shared<Tensor>(I459_index);
  auto tensor397 = vector<shared_ptr<Tensor>>{I458, t2, I459};
  auto task397 = make_shared<Task397>(tensor397, pindex);
  task396->add_dep(task397);
  task397->add_dep(task284);
  densityq->add_task(task397);

  auto tensor398 = vector<shared_ptr<Tensor>>{I459, Gamma28_(), t2};
  auto task398 = make_shared<Task398>(tensor398, pindex);
  task397->add_dep(task398);
  task398->add_dep(task284);
  densityq->add_task(task398);

  vector<IndexRange> I461_index = {active_, closed_};
  auto I461 = make_shared<Tensor>(I461_index);
  auto tensor399 = vector<shared_ptr<Tensor>>{den2, I461};
  auto task399 = make_shared<Task399>(tensor399, pindex);
  task399->add_dep(task284);
  densityq->add_task(task399);

  vector<IndexRange> I462_index = {closed_, virt_, active_, active_};
  auto I462 = make_shared<Tensor>(I462_index);
  auto tensor400 = vector<shared_ptr<Tensor>>{I461, t2, I462};
  auto task400 = make_shared<Task400>(tensor400, pindex);
  task399->add_dep(task400);
  task400->add_dep(task284);
  densityq->add_task(task400);

  auto tensor401 = vector<shared_ptr<Tensor>>{I462, Gamma29_(), t2};
  auto task401 = make_shared<Task401>(tensor401, pindex);
  task400->add_dep(task401);
  task401->add_dep(task284);
  densityq->add_task(task401);

  vector<IndexRange> I465_index = {closed_, virt_, active_, active_};
  auto I465 = make_shared<Tensor>(I465_index);
  auto tensor402 = vector<shared_ptr<Tensor>>{I461, t2, I465};
  auto task402 = make_shared<Task402>(tensor402, pindex);
  task399->add_dep(task402);
  task402->add_dep(task284);
  densityq->add_task(task402);

  auto tensor403 = vector<shared_ptr<Tensor>>{I465, Gamma7_(), t2};
  auto task403 = make_shared<Task403>(tensor403, pindex);
  task402->add_dep(task403);
  task403->add_dep(task284);
  densityq->add_task(task403);

  vector<IndexRange> I504_index = {virt_, closed_, active_, active_};
  auto I504 = make_shared<Tensor>(I504_index);
  auto tensor404 = vector<shared_ptr<Tensor>>{I461, t2, I504};
  auto task404 = make_shared<Task404>(tensor404, pindex);
  task399->add_dep(task404);
  task404->add_dep(task284);
  densityq->add_task(task404);

  auto tensor405 = vector<shared_ptr<Tensor>>{I504, Gamma7_(), t2};
  auto task405 = make_shared<Task405>(tensor405, pindex);
  task404->add_dep(task405);
  task405->add_dep(task284);
  densityq->add_task(task405);

  vector<IndexRange> I507_index = {virt_, closed_, active_, active_};
  auto I507 = make_shared<Tensor>(I507_index);
  auto tensor406 = vector<shared_ptr<Tensor>>{I461, t2, I507};
  auto task406 = make_shared<Task406>(tensor406, pindex);
  task399->add_dep(task406);
  task406->add_dep(task284);
  densityq->add_task(task406);

  auto tensor407 = vector<shared_ptr<Tensor>>{I507, Gamma7_(), t2};
  auto task407 = make_shared<Task407>(tensor407, pindex);
  task406->add_dep(task407);
  task407->add_dep(task284);
  densityq->add_task(task407);

  vector<IndexRange> I656_index = {virt_, virt_, active_, active_};
  auto I656 = make_shared<Tensor>(I656_index);
  auto tensor408 = vector<shared_ptr<Tensor>>{I461, t2, I656};
  auto task408 = make_shared<Task408>(tensor408, pindex);
  task399->add_dep(task408);
  task408->add_dep(task284);
  densityq->add_task(task408);

  auto tensor409 = vector<shared_ptr<Tensor>>{I656, Gamma60_(), t2};
  auto task409 = make_shared<Task409>(tensor409, pindex);
  task408->add_dep(task409);
  task409->add_dep(task284);
  densityq->add_task(task409);

  vector<IndexRange> I473_index = {virt_, virt_};
  auto I473 = make_shared<Tensor>(I473_index);
  auto tensor410 = vector<shared_ptr<Tensor>>{den2, I473};
  auto task410 = make_shared<Task410>(tensor410, pindex);
  task410->add_dep(task284);
  densityq->add_task(task410);

  vector<IndexRange> I474_index = {closed_, virt_, active_, active_};
  auto I474 = make_shared<Tensor>(I474_index);
  auto tensor411 = vector<shared_ptr<Tensor>>{I473, t2, I474};
  auto task411 = make_shared<Task411>(tensor411, pindex);
  task410->add_dep(task411);
  task411->add_dep(task284);
  densityq->add_task(task411);

  auto tensor412 = vector<shared_ptr<Tensor>>{I474, Gamma32_(), t2};
  auto task412 = make_shared<Task412>(tensor412, pindex);
  task411->add_dep(task412);
  task412->add_dep(task284);
  densityq->add_task(task412);

  vector<IndexRange> I483_index = {closed_, virt_, active_, active_};
  auto I483 = make_shared<Tensor>(I483_index);
  auto tensor413 = vector<shared_ptr<Tensor>>{I473, t2, I483};
  auto task413 = make_shared<Task413>(tensor413, pindex);
  task410->add_dep(task413);
  task413->add_dep(task284);
  densityq->add_task(task413);

  auto tensor414 = vector<shared_ptr<Tensor>>{I483, Gamma35_(), t2};
  auto task414 = make_shared<Task414>(tensor414, pindex);
  task413->add_dep(task414);
  task414->add_dep(task284);
  densityq->add_task(task414);

  vector<IndexRange> I488_index = {virt_, closed_};
  auto I488 = make_shared<Tensor>(I488_index);
  auto tensor415 = vector<shared_ptr<Tensor>>{den2, I488};
  auto task415 = make_shared<Task415>(tensor415, pindex);
  task415->add_dep(task284);
  densityq->add_task(task415);

  vector<IndexRange> I489_index = {closed_, virt_};
  auto I489 = make_shared<Tensor>(I489_index);
  auto tensor416 = vector<shared_ptr<Tensor>>{I488, t2, I489};
  auto task416 = make_shared<Task416>(tensor416, pindex);
  task415->add_dep(task416);
  task416->add_dep(task284);
  densityq->add_task(task416);

  auto tensor417 = vector<shared_ptr<Tensor>>{I489, Gamma38_(), t2};
  auto task417 = make_shared<Task417>(tensor417, pindex);
  task416->add_dep(task417);
  task417->add_dep(task284);
  densityq->add_task(task417);

  vector<IndexRange> I492_index = {closed_, virt_};
  auto I492 = make_shared<Tensor>(I492_index);
  auto tensor418 = vector<shared_ptr<Tensor>>{I488, t2, I492};
  auto task418 = make_shared<Task418>(tensor418, pindex);
  task415->add_dep(task418);
  task418->add_dep(task284);
  densityq->add_task(task418);

  auto tensor419 = vector<shared_ptr<Tensor>>{I492, Gamma38_(), t2};
  auto task419 = make_shared<Task419>(tensor419, pindex);
  task418->add_dep(task419);
  task419->add_dep(task284);
  densityq->add_task(task419);

  vector<IndexRange> I531_index = {virt_, closed_};
  auto I531 = make_shared<Tensor>(I531_index);
  auto tensor420 = vector<shared_ptr<Tensor>>{I488, t2, I531};
  auto task420 = make_shared<Task420>(tensor420, pindex);
  task415->add_dep(task420);
  task420->add_dep(task284);
  densityq->add_task(task420);

  auto tensor421 = vector<shared_ptr<Tensor>>{I531, Gamma38_(), t2};
  auto task421 = make_shared<Task421>(tensor421, pindex);
  task420->add_dep(task421);
  task421->add_dep(task284);
  densityq->add_task(task421);

  vector<IndexRange> I534_index = {virt_, closed_};
  auto I534 = make_shared<Tensor>(I534_index);
  auto tensor422 = vector<shared_ptr<Tensor>>{I488, t2, I534};
  auto task422 = make_shared<Task422>(tensor422, pindex);
  task415->add_dep(task422);
  task422->add_dep(task284);
  densityq->add_task(task422);

  auto tensor423 = vector<shared_ptr<Tensor>>{I534, Gamma38_(), t2};
  auto task423 = make_shared<Task423>(tensor423, pindex);
  task422->add_dep(task423);
  task423->add_dep(task284);
  densityq->add_task(task423);
}

#endif
