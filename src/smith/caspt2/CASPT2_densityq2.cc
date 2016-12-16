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

void CASPT2::CASPT2::make_densityq2(shared_ptr<Queue> densityq, shared_ptr<Task> task285, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I387_index = {active_, virt_};
  auto I387 = make_shared<Tensor>(I387_index);
  auto tensor344 = vector<shared_ptr<Tensor>>{den2, I387};
  auto task344 = make_shared<Task344>(tensor344, pindex);
  task344->add_dep(task285);
  densityq->add_task(task344);

  vector<IndexRange> I388_index = {closed_, active_, active_, active_};
  auto I388 = make_shared<Tensor>(I388_index);
  auto tensor345 = vector<shared_ptr<Tensor>>{I387, t2, I388};
  auto task345 = make_shared<Task345>(tensor345, pindex);
  task344->add_dep(task345);
  task345->add_dep(task285);
  densityq->add_task(task345);

  auto tensor346 = vector<shared_ptr<Tensor>>{I388, Gamma9_(), t2};
  auto task346 = make_shared<Task346>(tensor346, pindex);
  task345->add_dep(task346);
  task346->add_dep(task285);
  densityq->add_task(task346);

  vector<IndexRange> I391_index = {closed_, active_, active_, active_};
  auto I391 = make_shared<Tensor>(I391_index);
  auto tensor347 = vector<shared_ptr<Tensor>>{I387, t2, I391};
  auto task347 = make_shared<Task347>(tensor347, pindex);
  task344->add_dep(task347);
  task347->add_dep(task285);
  densityq->add_task(task347);

  auto tensor348 = vector<shared_ptr<Tensor>>{I391, Gamma6_(), t2};
  auto task348 = make_shared<Task348>(tensor348, pindex);
  task347->add_dep(task348);
  task348->add_dep(task285);
  densityq->add_task(task348);

  vector<IndexRange> I547_index = {virt_, active_, active_, active_};
  auto I547 = make_shared<Tensor>(I547_index);
  auto tensor349 = vector<shared_ptr<Tensor>>{I387, t2, I547};
  auto task349 = make_shared<Task349>(tensor349, pindex);
  task344->add_dep(task349);
  task349->add_dep(task285);
  densityq->add_task(task349);

  auto tensor350 = vector<shared_ptr<Tensor>>{I547, Gamma59_(), t2};
  auto task350 = make_shared<Task350>(tensor350, pindex);
  task349->add_dep(task350);
  task350->add_dep(task285);
  densityq->add_task(task350);

  vector<IndexRange> I393_index = {active_, virt_};
  auto I393 = make_shared<Tensor>(I393_index);
  auto tensor351 = vector<shared_ptr<Tensor>>{den2, I393};
  auto task351 = make_shared<Task351>(tensor351, pindex);
  task351->add_dep(task285);
  densityq->add_task(task351);

  vector<IndexRange> I394_index = {closed_, closed_, active_, active_};
  auto I394 = make_shared<Tensor>(I394_index);
  auto tensor352 = vector<shared_ptr<Tensor>>{I393, t2, I394};
  auto task352 = make_shared<Task352>(tensor352, pindex);
  task351->add_dep(task352);
  task352->add_dep(task285);
  densityq->add_task(task352);

  auto tensor353 = vector<shared_ptr<Tensor>>{I394, Gamma3_(), t2};
  auto task353 = make_shared<Task353>(tensor353, pindex);
  task352->add_dep(task353);
  task353->add_dep(task285);
  densityq->add_task(task353);

  vector<IndexRange> I396_index = {virt_, closed_};
  auto I396 = make_shared<Tensor>(I396_index);
  auto tensor354 = vector<shared_ptr<Tensor>>{den2, I396};
  auto task354 = make_shared<Task354>(tensor354, pindex);
  task354->add_dep(task285);
  densityq->add_task(task354);

  vector<IndexRange> I397_index = {closed_, active_};
  auto I397 = make_shared<Tensor>(I397_index);
  auto tensor355 = vector<shared_ptr<Tensor>>{I396, t2, I397};
  auto task355 = make_shared<Task355>(tensor355, pindex);
  task354->add_dep(task355);
  task355->add_dep(task285);
  densityq->add_task(task355);

  auto tensor356 = vector<shared_ptr<Tensor>>{I397, Gamma12_(), t2};
  auto task356 = make_shared<Task356>(tensor356, pindex);
  task355->add_dep(task356);
  task356->add_dep(task285);
  densityq->add_task(task356);

  vector<IndexRange> I399_index = {closed_, virt_};
  auto I399 = make_shared<Tensor>(I399_index);
  auto tensor357 = vector<shared_ptr<Tensor>>{den2, I399};
  auto task357 = make_shared<Task357>(tensor357, pindex);
  task357->add_dep(task285);
  densityq->add_task(task357);

  vector<IndexRange> I400_index = {closed_, active_};
  auto I400 = make_shared<Tensor>(I400_index);
  auto tensor358 = vector<shared_ptr<Tensor>>{I399, t2, I400};
  auto task358 = make_shared<Task358>(tensor358, pindex);
  task357->add_dep(task358);
  task358->add_dep(task285);
  densityq->add_task(task358);

  auto tensor359 = vector<shared_ptr<Tensor>>{I400, Gamma12_(), t2};
  auto task359 = make_shared<Task359>(tensor359, pindex);
  task358->add_dep(task359);
  task359->add_dep(task285);
  densityq->add_task(task359);

  vector<IndexRange> I556_index = {virt_, closed_};
  auto I556 = make_shared<Tensor>(I556_index);
  auto tensor360 = vector<shared_ptr<Tensor>>{I399, t2, I556};
  auto task360 = make_shared<Task360>(tensor360, pindex);
  task357->add_dep(task360);
  task360->add_dep(task285);
  densityq->add_task(task360);

  vector<IndexRange> I557_index = {active_, virt_, closed_, active_};
  auto I557 = make_shared<Tensor>(I557_index);
  auto tensor361 = vector<shared_ptr<Tensor>>{I556, Gamma38_(), I557};
  auto task361 = make_shared<Task361>(tensor361, pindex);
  task360->add_dep(task361);
  task361->add_dep(task285);
  densityq->add_task(task361);

  auto tensor362 = vector<shared_ptr<Tensor>>{I557, t2};
  auto task362 = make_shared<Task362>(tensor362, pindex);
  task361->add_dep(task362);
  task362->add_dep(task285);
  densityq->add_task(task362);

  vector<IndexRange> I402_index = {active_, active_};
  auto I402 = make_shared<Tensor>(I402_index);
  auto tensor363 = vector<shared_ptr<Tensor>>{den2, I402};
  auto task363 = make_shared<Task363>(tensor363, pindex);
  task363->add_dep(task285);
  densityq->add_task(task363);

  vector<IndexRange> I403_index = {active_, active_};
  auto I403 = make_shared<Tensor>(I403_index);
  auto tensor364 = vector<shared_ptr<Tensor>>{I402, Gamma152_(), I403};
  auto task364 = make_shared<Task364>(tensor364, pindex);
  task363->add_dep(task364);
  task364->add_dep(task285);
  densityq->add_task(task364);

  auto tensor365 = vector<shared_ptr<Tensor>>{I403, t2};
  auto task365 = make_shared<Task365>(tensor365, pindex);
  task364->add_dep(task365);
  task365->add_dep(task285);
  densityq->add_task(task365);

  auto tensor366 = vector<shared_ptr<Tensor>>{I403, t2};
  auto task366 = make_shared<Task366>(tensor366, pindex);
  task364->add_dep(task366);
  task366->add_dep(task285);
  densityq->add_task(task366);

  vector<IndexRange> I612_index = {active_, active_};
  auto I612 = make_shared<Tensor>(I612_index);
  auto tensor367 = vector<shared_ptr<Tensor>>{I402, Gamma60_(), I612};
  auto task367 = make_shared<Task367>(tensor367, pindex);
  task363->add_dep(task367);
  task367->add_dep(task285);
  densityq->add_task(task367);

  auto tensor368 = vector<shared_ptr<Tensor>>{I612, t2};
  auto task368 = make_shared<Task368>(tensor368, pindex);
  task367->add_dep(task368);
  task368->add_dep(task285);
  densityq->add_task(task368);

  auto tensor369 = vector<shared_ptr<Tensor>>{I612, t2};
  auto task369 = make_shared<Task369>(tensor369, pindex);
  task367->add_dep(task369);
  task369->add_dep(task285);
  densityq->add_task(task369);

  vector<IndexRange> I408_index = {closed_, closed_};
  auto I408 = make_shared<Tensor>(I408_index);
  auto tensor370 = vector<shared_ptr<Tensor>>{den2, I408};
  auto task370 = make_shared<Task370>(tensor370, pindex);
  task370->add_dep(task285);
  densityq->add_task(task370);

  vector<IndexRange> I409_index = {closed_, virt_, closed_, active_};
  auto I409 = make_shared<Tensor>(I409_index);
  auto tensor371 = vector<shared_ptr<Tensor>>{I408, t2, I409};
  auto task371 = make_shared<Task371>(tensor371, pindex);
  task370->add_dep(task371);
  task371->add_dep(task285);
  densityq->add_task(task371);

  auto tensor372 = vector<shared_ptr<Tensor>>{I409, Gamma16_(), t2};
  auto task372 = make_shared<Task372>(tensor372, pindex);
  task371->add_dep(task372);
  task372->add_dep(task285);
  densityq->add_task(task372);

  vector<IndexRange> I412_index = {closed_, virt_, closed_, active_};
  auto I412 = make_shared<Tensor>(I412_index);
  auto tensor373 = vector<shared_ptr<Tensor>>{I408, t2, I412};
  auto task373 = make_shared<Task373>(tensor373, pindex);
  task370->add_dep(task373);
  task373->add_dep(task285);
  densityq->add_task(task373);

  auto tensor374 = vector<shared_ptr<Tensor>>{I412, Gamma16_(), t2};
  auto task374 = make_shared<Task374>(tensor374, pindex);
  task373->add_dep(task374);
  task374->add_dep(task285);
  densityq->add_task(task374);

  vector<IndexRange> I414_index = {closed_, closed_};
  auto I414 = make_shared<Tensor>(I414_index);
  auto tensor375 = vector<shared_ptr<Tensor>>{den2, I414};
  auto task375 = make_shared<Task375>(tensor375, pindex);
  task375->add_dep(task285);
  densityq->add_task(task375);

  vector<IndexRange> I415_index = {closed_, virt_, closed_, active_};
  auto I415 = make_shared<Tensor>(I415_index);
  auto tensor376 = vector<shared_ptr<Tensor>>{I414, t2, I415};
  auto task376 = make_shared<Task376>(tensor376, pindex);
  task375->add_dep(task376);
  task376->add_dep(task285);
  densityq->add_task(task376);

  auto tensor377 = vector<shared_ptr<Tensor>>{I415, Gamma16_(), t2};
  auto task377 = make_shared<Task377>(tensor377, pindex);
  task376->add_dep(task377);
  task377->add_dep(task285);
  densityq->add_task(task377);

  vector<IndexRange> I421_index = {closed_, virt_, closed_, active_};
  auto I421 = make_shared<Tensor>(I421_index);
  auto tensor378 = vector<shared_ptr<Tensor>>{I414, t2, I421};
  auto task378 = make_shared<Task378>(tensor378, pindex);
  task375->add_dep(task378);
  task378->add_dep(task285);
  densityq->add_task(task378);

  auto tensor379 = vector<shared_ptr<Tensor>>{I421, Gamma16_(), t2};
  auto task379 = make_shared<Task379>(tensor379, pindex);
  task378->add_dep(task379);
  task379->add_dep(task285);
  densityq->add_task(task379);

  vector<IndexRange> I417_index = {virt_, virt_};
  auto I417 = make_shared<Tensor>(I417_index);
  auto tensor380 = vector<shared_ptr<Tensor>>{den2, I417};
  auto task380 = make_shared<Task380>(tensor380, pindex);
  task380->add_dep(task285);
  densityq->add_task(task380);

  vector<IndexRange> I418_index = {closed_, virt_, closed_, active_};
  auto I418 = make_shared<Tensor>(I418_index);
  auto tensor381 = vector<shared_ptr<Tensor>>{I417, t2, I418};
  auto task381 = make_shared<Task381>(tensor381, pindex);
  task380->add_dep(task381);
  task381->add_dep(task285);
  densityq->add_task(task381);

  auto tensor382 = vector<shared_ptr<Tensor>>{I418, Gamma16_(), t2};
  auto task382 = make_shared<Task382>(tensor382, pindex);
  task381->add_dep(task382);
  task382->add_dep(task285);
  densityq->add_task(task382);

  vector<IndexRange> I424_index = {closed_, virt_, closed_, active_};
  auto I424 = make_shared<Tensor>(I424_index);
  auto tensor383 = vector<shared_ptr<Tensor>>{I417, t2, I424};
  auto task383 = make_shared<Task383>(tensor383, pindex);
  task380->add_dep(task383);
  task383->add_dep(task285);
  densityq->add_task(task383);

  auto tensor384 = vector<shared_ptr<Tensor>>{I424, Gamma16_(), t2};
  auto task384 = make_shared<Task384>(tensor384, pindex);
  task383->add_dep(task384);
  task384->add_dep(task285);
  densityq->add_task(task384);

  vector<IndexRange> I426_index = {active_, closed_};
  auto I426 = make_shared<Tensor>(I426_index);
  auto tensor385 = vector<shared_ptr<Tensor>>{den2, I426};
  auto task385 = make_shared<Task385>(tensor385, pindex);
  task385->add_dep(task285);
  densityq->add_task(task385);

  vector<IndexRange> I427_index = {virt_, closed_, active_, active_};
  auto I427 = make_shared<Tensor>(I427_index);
  auto tensor386 = vector<shared_ptr<Tensor>>{I426, t2, I427};
  auto task386 = make_shared<Task386>(tensor386, pindex);
  task385->add_dep(task386);
  task386->add_dep(task285);
  densityq->add_task(task386);

  auto tensor387 = vector<shared_ptr<Tensor>>{I427, Gamma22_(), t2};
  auto task387 = make_shared<Task387>(tensor387, pindex);
  task386->add_dep(task387);
  task387->add_dep(task285);
  densityq->add_task(task387);

  auto tensor388 = vector<shared_ptr<Tensor>>{I427, Gamma12_(), t2};
  auto task388 = make_shared<Task388>(tensor388, pindex);
  task386->add_dep(task388);
  task388->add_dep(task285);
  densityq->add_task(task388);

  vector<IndexRange> I429_index = {active_, closed_};
  auto I429 = make_shared<Tensor>(I429_index);
  auto tensor389 = vector<shared_ptr<Tensor>>{den2, I429};
  auto task389 = make_shared<Task389>(tensor389, pindex);
  task389->add_dep(task285);
  densityq->add_task(task389);

  vector<IndexRange> I430_index = {virt_, closed_, active_, active_};
  auto I430 = make_shared<Tensor>(I430_index);
  auto tensor390 = vector<shared_ptr<Tensor>>{I429, t2, I430};
  auto task390 = make_shared<Task390>(tensor390, pindex);
  task389->add_dep(task390);
  task390->add_dep(task285);
  densityq->add_task(task390);

  vector<IndexRange> I431_index = {active_, virt_, closed_, active_};
  auto I431 = make_shared<Tensor>(I431_index);
  auto tensor391 = vector<shared_ptr<Tensor>>{I430, Gamma12_(), I431};
  auto task391 = make_shared<Task391>(tensor391, pindex);
  task390->add_dep(task391);
  task391->add_dep(task285);
  densityq->add_task(task391);

  auto tensor392 = vector<shared_ptr<Tensor>>{I431, t2};
  auto task392 = make_shared<Task392>(tensor392, pindex);
  task391->add_dep(task392);
  task392->add_dep(task285);
  densityq->add_task(task392);

  vector<IndexRange> I438_index = {virt_, active_};
  auto I438 = make_shared<Tensor>(I438_index);
  auto tensor393 = vector<shared_ptr<Tensor>>{den2, I438};
  auto task393 = make_shared<Task393>(tensor393, pindex);
  task393->add_dep(task285);
  densityq->add_task(task393);

  vector<IndexRange> I439_index = {active_, virt_};
  auto I439 = make_shared<Tensor>(I439_index);
  auto tensor394 = vector<shared_ptr<Tensor>>{I438, Gamma16_(), I439};
  auto task394 = make_shared<Task394>(tensor394, pindex);
  task393->add_dep(task394);
  task394->add_dep(task285);
  densityq->add_task(task394);

  auto tensor395 = vector<shared_ptr<Tensor>>{I439, t2};
  auto task395 = make_shared<Task395>(tensor395, pindex);
  task394->add_dep(task395);
  task395->add_dep(task285);
  densityq->add_task(task395);

  auto tensor396 = vector<shared_ptr<Tensor>>{I439, t2};
  auto task396 = make_shared<Task396>(tensor396, pindex);
  task394->add_dep(task396);
  task396->add_dep(task285);
  densityq->add_task(task396);

  vector<IndexRange> I444_index = {active_, virt_};
  auto I444 = make_shared<Tensor>(I444_index);
  auto tensor397 = vector<shared_ptr<Tensor>>{den2, I444};
  auto task397 = make_shared<Task397>(tensor397, pindex);
  task397->add_dep(task285);
  densityq->add_task(task397);

  vector<IndexRange> I445_index = {closed_, active_, active_, active_};
  auto I445 = make_shared<Tensor>(I445_index);
  auto tensor398 = vector<shared_ptr<Tensor>>{I444, t2, I445};
  auto task398 = make_shared<Task398>(tensor398, pindex);
  task397->add_dep(task398);
  task398->add_dep(task285);
  densityq->add_task(task398);

  auto tensor399 = vector<shared_ptr<Tensor>>{I445, Gamma28_(), t2};
  auto task399 = make_shared<Task399>(tensor399, pindex);
  task398->add_dep(task399);
  task399->add_dep(task285);
  densityq->add_task(task399);

  vector<IndexRange> I447_index = {active_, closed_};
  auto I447 = make_shared<Tensor>(I447_index);
  auto tensor400 = vector<shared_ptr<Tensor>>{den2, I447};
  auto task400 = make_shared<Task400>(tensor400, pindex);
  task400->add_dep(task285);
  densityq->add_task(task400);

  vector<IndexRange> I448_index = {closed_, virt_, active_, active_};
  auto I448 = make_shared<Tensor>(I448_index);
  auto tensor401 = vector<shared_ptr<Tensor>>{I447, t2, I448};
  auto task401 = make_shared<Task401>(tensor401, pindex);
  task400->add_dep(task401);
  task401->add_dep(task285);
  densityq->add_task(task401);

  auto tensor402 = vector<shared_ptr<Tensor>>{I448, Gamma29_(), t2};
  auto task402 = make_shared<Task402>(tensor402, pindex);
  task401->add_dep(task402);
  task402->add_dep(task285);
  densityq->add_task(task402);

  vector<IndexRange> I451_index = {closed_, virt_, active_, active_};
  auto I451 = make_shared<Tensor>(I451_index);
  auto tensor403 = vector<shared_ptr<Tensor>>{I447, t2, I451};
  auto task403 = make_shared<Task403>(tensor403, pindex);
  task400->add_dep(task403);
  task403->add_dep(task285);
  densityq->add_task(task403);

  auto tensor404 = vector<shared_ptr<Tensor>>{I451, Gamma7_(), t2};
  auto task404 = make_shared<Task404>(tensor404, pindex);
  task403->add_dep(task404);
  task404->add_dep(task285);
  densityq->add_task(task404);

  vector<IndexRange> I490_index = {virt_, closed_, active_, active_};
  auto I490 = make_shared<Tensor>(I490_index);
  auto tensor405 = vector<shared_ptr<Tensor>>{I447, t2, I490};
  auto task405 = make_shared<Task405>(tensor405, pindex);
  task400->add_dep(task405);
  task405->add_dep(task285);
  densityq->add_task(task405);

  auto tensor406 = vector<shared_ptr<Tensor>>{I490, Gamma7_(), t2};
  auto task406 = make_shared<Task406>(tensor406, pindex);
  task405->add_dep(task406);
  task406->add_dep(task285);
  densityq->add_task(task406);

  vector<IndexRange> I493_index = {virt_, closed_, active_, active_};
  auto I493 = make_shared<Tensor>(I493_index);
  auto tensor407 = vector<shared_ptr<Tensor>>{I447, t2, I493};
  auto task407 = make_shared<Task407>(tensor407, pindex);
  task400->add_dep(task407);
  task407->add_dep(task285);
  densityq->add_task(task407);

  auto tensor408 = vector<shared_ptr<Tensor>>{I493, Gamma7_(), t2};
  auto task408 = make_shared<Task408>(tensor408, pindex);
  task407->add_dep(task408);
  task408->add_dep(task285);
  densityq->add_task(task408);

  vector<IndexRange> I642_index = {virt_, virt_, active_, active_};
  auto I642 = make_shared<Tensor>(I642_index);
  auto tensor409 = vector<shared_ptr<Tensor>>{I447, t2, I642};
  auto task409 = make_shared<Task409>(tensor409, pindex);
  task400->add_dep(task409);
  task409->add_dep(task285);
  densityq->add_task(task409);

  auto tensor410 = vector<shared_ptr<Tensor>>{I642, Gamma60_(), t2};
  auto task410 = make_shared<Task410>(tensor410, pindex);
  task409->add_dep(task410);
  task410->add_dep(task285);
  densityq->add_task(task410);

  vector<IndexRange> I459_index = {virt_, virt_};
  auto I459 = make_shared<Tensor>(I459_index);
  auto tensor411 = vector<shared_ptr<Tensor>>{den2, I459};
  auto task411 = make_shared<Task411>(tensor411, pindex);
  task411->add_dep(task285);
  densityq->add_task(task411);

  vector<IndexRange> I460_index = {closed_, virt_, active_, active_};
  auto I460 = make_shared<Tensor>(I460_index);
  auto tensor412 = vector<shared_ptr<Tensor>>{I459, t2, I460};
  auto task412 = make_shared<Task412>(tensor412, pindex);
  task411->add_dep(task412);
  task412->add_dep(task285);
  densityq->add_task(task412);

  auto tensor413 = vector<shared_ptr<Tensor>>{I460, Gamma32_(), t2};
  auto task413 = make_shared<Task413>(tensor413, pindex);
  task412->add_dep(task413);
  task413->add_dep(task285);
  densityq->add_task(task413);

  vector<IndexRange> I469_index = {closed_, virt_, active_, active_};
  auto I469 = make_shared<Tensor>(I469_index);
  auto tensor414 = vector<shared_ptr<Tensor>>{I459, t2, I469};
  auto task414 = make_shared<Task414>(tensor414, pindex);
  task411->add_dep(task414);
  task414->add_dep(task285);
  densityq->add_task(task414);

  auto tensor415 = vector<shared_ptr<Tensor>>{I469, Gamma35_(), t2};
  auto task415 = make_shared<Task415>(tensor415, pindex);
  task414->add_dep(task415);
  task415->add_dep(task285);
  densityq->add_task(task415);

  vector<IndexRange> I474_index = {virt_, closed_};
  auto I474 = make_shared<Tensor>(I474_index);
  auto tensor416 = vector<shared_ptr<Tensor>>{den2, I474};
  auto task416 = make_shared<Task416>(tensor416, pindex);
  task416->add_dep(task285);
  densityq->add_task(task416);

  vector<IndexRange> I475_index = {closed_, virt_};
  auto I475 = make_shared<Tensor>(I475_index);
  auto tensor417 = vector<shared_ptr<Tensor>>{I474, t2, I475};
  auto task417 = make_shared<Task417>(tensor417, pindex);
  task416->add_dep(task417);
  task417->add_dep(task285);
  densityq->add_task(task417);

  auto tensor418 = vector<shared_ptr<Tensor>>{I475, Gamma38_(), t2};
  auto task418 = make_shared<Task418>(tensor418, pindex);
  task417->add_dep(task418);
  task418->add_dep(task285);
  densityq->add_task(task418);

  vector<IndexRange> I478_index = {closed_, virt_};
  auto I478 = make_shared<Tensor>(I478_index);
  auto tensor419 = vector<shared_ptr<Tensor>>{I474, t2, I478};
  auto task419 = make_shared<Task419>(tensor419, pindex);
  task416->add_dep(task419);
  task419->add_dep(task285);
  densityq->add_task(task419);

  auto tensor420 = vector<shared_ptr<Tensor>>{I478, Gamma38_(), t2};
  auto task420 = make_shared<Task420>(tensor420, pindex);
  task419->add_dep(task420);
  task420->add_dep(task285);
  densityq->add_task(task420);

  vector<IndexRange> I517_index = {virt_, closed_};
  auto I517 = make_shared<Tensor>(I517_index);
  auto tensor421 = vector<shared_ptr<Tensor>>{I474, t2, I517};
  auto task421 = make_shared<Task421>(tensor421, pindex);
  task416->add_dep(task421);
  task421->add_dep(task285);
  densityq->add_task(task421);

  auto tensor422 = vector<shared_ptr<Tensor>>{I517, Gamma38_(), t2};
  auto task422 = make_shared<Task422>(tensor422, pindex);
  task421->add_dep(task422);
  task422->add_dep(task285);
  densityq->add_task(task422);

  vector<IndexRange> I520_index = {virt_, closed_};
  auto I520 = make_shared<Tensor>(I520_index);
  auto tensor423 = vector<shared_ptr<Tensor>>{I474, t2, I520};
  auto task423 = make_shared<Task423>(tensor423, pindex);
  task416->add_dep(task423);
  task423->add_dep(task285);
  densityq->add_task(task423);

  auto tensor424 = vector<shared_ptr<Tensor>>{I520, Gamma38_(), t2};
  auto task424 = make_shared<Task424>(tensor424, pindex);
  task423->add_dep(task424);
  task424->add_dep(task285);
  densityq->add_task(task424);
}

#endif
