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
#include <src/smith/mrci/MRCI_tasks7.h>
#include <src/smith/mrci/MRCI_tasks8.h>
#include <src/smith/mrci/MRCI_tasks9.h>
#include <src/smith/mrci/MRCI_tasks10.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void MRCI::MRCI::make_residualq4(shared_ptr<Queue> residualq, shared_ptr<Task> task108, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I72_index = {closed_, active_, active_, virt_};
  auto I72 = make_shared<Tensor>(I72_index);
  auto tensor344 = vector<shared_ptr<Tensor>>{r, I72};
  auto task344 = make_shared<Task344>(tensor344, pindex);
  task344->add_dep(task108);
  residualq->add_task(task344);

  vector<IndexRange> I73_index = {closed_, active_, active_, active_};
  auto I73 = make_shared<Tensor>(I73_index);
  auto tensor345 = vector<shared_ptr<Tensor>>{I72, h1_, I73};
  auto task345 = make_shared<Task345>(tensor345, pindex);
  task344->add_dep(task345);
  task345->add_dep(task108);
  residualq->add_task(task345);

  auto tensor346 = vector<shared_ptr<Tensor>>{I73, Gamma24_(), t2};
  auto task346 = make_shared<Task346>(tensor346, pindex);
  task345->add_dep(task346);
  task346->add_dep(task108);
  residualq->add_task(task346);

  vector<IndexRange> I76_index = {active_, virt_, closed_, active_};
  auto I76 = make_shared<Tensor>(I76_index);
  auto tensor347 = vector<shared_ptr<Tensor>>{I72, Gamma25_(), I76};
  auto task347 = make_shared<Task347>(tensor347, pindex);
  task344->add_dep(task347);
  task347->add_dep(task108);
  residualq->add_task(task347);

  auto tensor348 = vector<shared_ptr<Tensor>>{I76, t2, h1_};
  auto task348 = make_shared<Task348>(tensor348, pindex);
  task347->add_dep(task348);
  task348->add_dep(task108);
  residualq->add_task(task348);

  auto tensor349 = vector<shared_ptr<Tensor>>{I76, t2, v2_};
  auto task349 = make_shared<Task349>(tensor349, pindex);
  task347->add_dep(task349);
  task349->add_dep(task108);
  residualq->add_task(task349);

  auto tensor350 = vector<shared_ptr<Tensor>>{I76, t2, v2_};
  auto task350 = make_shared<Task350>(tensor350, pindex);
  task347->add_dep(task350);
  task350->add_dep(task108);
  residualq->add_task(task350);

  auto tensor351 = vector<shared_ptr<Tensor>>{I76, t2, v2_};
  auto task351 = make_shared<Task351>(tensor351, pindex);
  task347->add_dep(task351);
  task351->add_dep(task108);
  residualq->add_task(task351);

  auto tensor352 = vector<shared_ptr<Tensor>>{I76, t2, v2_};
  auto task352 = make_shared<Task352>(tensor352, pindex);
  task347->add_dep(task352);
  task352->add_dep(task108);
  residualq->add_task(task352);

  auto tensor353 = vector<shared_ptr<Tensor>>{I76, t2, v2_};
  auto task353 = make_shared<Task353>(tensor353, pindex);
  task347->add_dep(task353);
  task353->add_dep(task108);
  residualq->add_task(task353);

  vector<IndexRange> I79_index = {active_, closed_, virt_, active_};
  auto I79 = make_shared<Tensor>(I79_index);
  auto tensor354 = vector<shared_ptr<Tensor>>{I72, Gamma5_(), I79};
  auto task354 = make_shared<Task354>(tensor354, pindex);
  task344->add_dep(task354);
  task354->add_dep(task108);
  residualq->add_task(task354);

  auto tensor355 = vector<shared_ptr<Tensor>>{I79, t2, h1_};
  auto task355 = make_shared<Task355>(tensor355, pindex);
  task354->add_dep(task355);
  task355->add_dep(task108);
  residualq->add_task(task355);

  auto tensor356 = vector<shared_ptr<Tensor>>{I79, t2, v2_};
  auto task356 = make_shared<Task356>(tensor356, pindex);
  task354->add_dep(task356);
  task356->add_dep(task108);
  residualq->add_task(task356);

  auto tensor357 = vector<shared_ptr<Tensor>>{I79, t2, v2_};
  auto task357 = make_shared<Task357>(tensor357, pindex);
  task354->add_dep(task357);
  task357->add_dep(task108);
  residualq->add_task(task357);

  auto tensor358 = vector<shared_ptr<Tensor>>{I79, t2, v2_};
  auto task358 = make_shared<Task358>(tensor358, pindex);
  task354->add_dep(task358);
  task358->add_dep(task108);
  residualq->add_task(task358);

  vector<IndexRange> I82_index = {closed_, active_, virt_, active_};
  auto I82 = make_shared<Tensor>(I82_index);
  auto tensor359 = vector<shared_ptr<Tensor>>{I72, Gamma27_(), I82};
  auto task359 = make_shared<Task359>(tensor359, pindex);
  task344->add_dep(task359);
  task359->add_dep(task108);
  residualq->add_task(task359);

  auto tensor360 = vector<shared_ptr<Tensor>>{I82, t2, h1_};
  auto task360 = make_shared<Task360>(tensor360, pindex);
  task359->add_dep(task360);
  task360->add_dep(task108);
  residualq->add_task(task360);

  auto tensor361 = vector<shared_ptr<Tensor>>{I82, t2, h1_};
  auto task361 = make_shared<Task361>(tensor361, pindex);
  task359->add_dep(task361);
  task361->add_dep(task108);
  residualq->add_task(task361);

  auto tensor362 = vector<shared_ptr<Tensor>>{I82, t2, h1_};
  auto task362 = make_shared<Task362>(tensor362, pindex);
  task359->add_dep(task362);
  task362->add_dep(task108);
  residualq->add_task(task362);

  auto tensor363 = vector<shared_ptr<Tensor>>{I82, t2, v2_};
  auto task363 = make_shared<Task363>(tensor363, pindex);
  task359->add_dep(task363);
  task363->add_dep(task108);
  residualq->add_task(task363);

  auto tensor364 = vector<shared_ptr<Tensor>>{I82, t2, v2_};
  auto task364 = make_shared<Task364>(tensor364, pindex);
  task359->add_dep(task364);
  task364->add_dep(task108);
  residualq->add_task(task364);

  vector<IndexRange> I823_index = {virt_, active_, closed_, closed_};
  auto I823 = make_shared<Tensor>(I823_index);
  auto tensor365 = vector<shared_ptr<Tensor>>{I82, t2, I823};
  auto task365 = make_shared<Task365>(tensor365, pindex);
  task359->add_dep(task365);
  task365->add_dep(task108);
  residualq->add_task(task365);

  auto tensor366 = vector<shared_ptr<Tensor>>{I823, v2_};
  auto task366 = make_shared<Task366>(tensor366, pindex);
  task365->add_dep(task366);
  task366->add_dep(task108);
  residualq->add_task(task366);

  auto tensor367 = vector<shared_ptr<Tensor>>{I82, t2, v2_};
  auto task367 = make_shared<Task367>(tensor367, pindex);
  task359->add_dep(task367);
  task367->add_dep(task108);
  residualq->add_task(task367);

  auto tensor368 = vector<shared_ptr<Tensor>>{I82, t2, v2_};
  auto task368 = make_shared<Task368>(tensor368, pindex);
  task359->add_dep(task368);
  task368->add_dep(task108);
  residualq->add_task(task368);

  vector<IndexRange> I88_index = {closed_, virt_, active_, active_};
  auto I88 = make_shared<Tensor>(I88_index);
  auto tensor369 = vector<shared_ptr<Tensor>>{I72, Gamma29_(), I88};
  auto task369 = make_shared<Task369>(tensor369, pindex);
  task344->add_dep(task369);
  task369->add_dep(task108);
  residualq->add_task(task369);

  auto tensor370 = vector<shared_ptr<Tensor>>{I88, t2, h1_};
  auto task370 = make_shared<Task370>(tensor370, pindex);
  task369->add_dep(task370);
  task370->add_dep(task108);
  residualq->add_task(task370);

  auto tensor371 = vector<shared_ptr<Tensor>>{I88, t2, h1_};
  auto task371 = make_shared<Task371>(tensor371, pindex);
  task369->add_dep(task371);
  task371->add_dep(task108);
  residualq->add_task(task371);

  auto tensor372 = vector<shared_ptr<Tensor>>{I88, t2, h1_};
  auto task372 = make_shared<Task372>(tensor372, pindex);
  task369->add_dep(task372);
  task372->add_dep(task108);
  residualq->add_task(task372);

  auto tensor373 = vector<shared_ptr<Tensor>>{I88, t2, v2_};
  auto task373 = make_shared<Task373>(tensor373, pindex);
  task369->add_dep(task373);
  task373->add_dep(task108);
  residualq->add_task(task373);

  auto tensor374 = vector<shared_ptr<Tensor>>{I88, t2, v2_};
  auto task374 = make_shared<Task374>(tensor374, pindex);
  task369->add_dep(task374);
  task374->add_dep(task108);
  residualq->add_task(task374);

  auto tensor375 = vector<shared_ptr<Tensor>>{I88, t2, v2_};
  auto task375 = make_shared<Task375>(tensor375, pindex);
  task369->add_dep(task375);
  task375->add_dep(task108);
  residualq->add_task(task375);

  vector<IndexRange> I772_index = {virt_, closed_, active_, active_};
  auto I772 = make_shared<Tensor>(I772_index);
  auto tensor376 = vector<shared_ptr<Tensor>>{I88, t2, I772};
  auto task376 = make_shared<Task376>(tensor376, pindex);
  task369->add_dep(task376);
  task376->add_dep(task108);
  residualq->add_task(task376);

  auto tensor377 = vector<shared_ptr<Tensor>>{I772, v2_};
  auto task377 = make_shared<Task377>(tensor377, pindex);
  task376->add_dep(task377);
  task377->add_dep(task108);
  residualq->add_task(task377);

  vector<IndexRange> I775_index = {virt_, closed_, active_, active_};
  auto I775 = make_shared<Tensor>(I775_index);
  auto tensor378 = vector<shared_ptr<Tensor>>{I88, t2, I775};
  auto task378 = make_shared<Task378>(tensor378, pindex);
  task369->add_dep(task378);
  task378->add_dep(task108);
  residualq->add_task(task378);

  auto tensor379 = vector<shared_ptr<Tensor>>{I775, v2_};
  auto task379 = make_shared<Task379>(tensor379, pindex);
  task378->add_dep(task379);
  task379->add_dep(task108);
  residualq->add_task(task379);

  auto tensor380 = vector<shared_ptr<Tensor>>{I88, t2, v2_};
  auto task380 = make_shared<Task380>(tensor380, pindex);
  task369->add_dep(task380);
  task380->add_dep(task108);
  residualq->add_task(task380);

  auto tensor381 = vector<shared_ptr<Tensor>>{I88, t2, v2_};
  auto task381 = make_shared<Task381>(tensor381, pindex);
  task369->add_dep(task381);
  task381->add_dep(task108);
  residualq->add_task(task381);

  auto tensor382 = vector<shared_ptr<Tensor>>{I88, t2, v2_};
  auto task382 = make_shared<Task382>(tensor382, pindex);
  task369->add_dep(task382);
  task382->add_dep(task108);
  residualq->add_task(task382);

  vector<IndexRange> I94_index = {virt_, active_, active_, active_};
  auto I94 = make_shared<Tensor>(I94_index);
  auto tensor383 = vector<shared_ptr<Tensor>>{I72, h1_, I94};
  auto task383 = make_shared<Task383>(tensor383, pindex);
  task344->add_dep(task383);
  task383->add_dep(task108);
  residualq->add_task(task383);

  auto tensor384 = vector<shared_ptr<Tensor>>{I94, Gamma31_(), t2};
  auto task384 = make_shared<Task384>(tensor384, pindex);
  task383->add_dep(task384);
  task384->add_dep(task108);
  residualq->add_task(task384);

  vector<IndexRange> I97_index = {closed_, virt_};
  auto I97 = make_shared<Tensor>(I97_index);
  auto tensor385 = vector<shared_ptr<Tensor>>{I72, Gamma32_(), I97};
  auto task385 = make_shared<Task385>(tensor385, pindex);
  task344->add_dep(task385);
  task385->add_dep(task108);
  residualq->add_task(task385);

  auto tensor386 = vector<shared_ptr<Tensor>>{I97, t2, h1_};
  auto task386 = make_shared<Task386>(tensor386, pindex);
  task385->add_dep(task386);
  task386->add_dep(task108);
  residualq->add_task(task386);

  auto tensor387 = vector<shared_ptr<Tensor>>{I97, t2, h1_};
  auto task387 = make_shared<Task387>(tensor387, pindex);
  task385->add_dep(task387);
  task387->add_dep(task108);
  residualq->add_task(task387);

  auto tensor388 = vector<shared_ptr<Tensor>>{I97, t2, v2_};
  auto task388 = make_shared<Task388>(tensor388, pindex);
  task385->add_dep(task388);
  task388->add_dep(task108);
  residualq->add_task(task388);

  auto tensor389 = vector<shared_ptr<Tensor>>{I97, t2, v2_};
  auto task389 = make_shared<Task389>(tensor389, pindex);
  task385->add_dep(task389);
  task389->add_dep(task108);
  residualq->add_task(task389);

  auto tensor390 = vector<shared_ptr<Tensor>>{I97, t2, v2_};
  auto task390 = make_shared<Task390>(tensor390, pindex);
  task385->add_dep(task390);
  task390->add_dep(task108);
  residualq->add_task(task390);

  auto tensor391 = vector<shared_ptr<Tensor>>{I97, t2, v2_};
  auto task391 = make_shared<Task391>(tensor391, pindex);
  task385->add_dep(task391);
  task391->add_dep(task108);
  residualq->add_task(task391);

  vector<IndexRange> I654_index = {closed_, closed_, active_, active_, active_, active_};
  auto I654 = make_shared<Tensor>(I654_index);
  auto tensor392 = vector<shared_ptr<Tensor>>{I72, v2_, I654};
  auto task392 = make_shared<Task392>(tensor392, pindex);
  task344->add_dep(task392);
  task392->add_dep(task108);
  residualq->add_task(task392);

  auto tensor393 = vector<shared_ptr<Tensor>>{I654, Gamma215_(), t2};
  auto task393 = make_shared<Task393>(tensor393, pindex);
  task392->add_dep(task393);
  task393->add_dep(task108);
  residualq->add_task(task393);

  vector<IndexRange> I657_index = {closed_, active_, active_, active_, active_, active_};
  auto I657 = make_shared<Tensor>(I657_index);
  auto tensor394 = vector<shared_ptr<Tensor>>{I72, v2_, I657};
  auto task394 = make_shared<Task394>(tensor394, pindex);
  task344->add_dep(task394);
  task394->add_dep(task108);
  residualq->add_task(task394);

  auto tensor395 = vector<shared_ptr<Tensor>>{I657, Gamma216_(), t2};
  auto task395 = make_shared<Task395>(tensor395, pindex);
  task394->add_dep(task395);
  task395->add_dep(task108);
  residualq->add_task(task395);

  vector<IndexRange> I660_index = {closed_, active_, active_, active_, active_, active_};
  auto I660 = make_shared<Tensor>(I660_index);
  auto tensor396 = vector<shared_ptr<Tensor>>{I72, v2_, I660};
  auto task396 = make_shared<Task396>(tensor396, pindex);
  task344->add_dep(task396);
  task396->add_dep(task108);
  residualq->add_task(task396);

  auto tensor397 = vector<shared_ptr<Tensor>>{I660, Gamma217_(), t2};
  auto task397 = make_shared<Task397>(tensor397, pindex);
  task396->add_dep(task397);
  task397->add_dep(task108);
  residualq->add_task(task397);

  vector<IndexRange> I663_index = {closed_, active_, active_, active_};
  auto I663 = make_shared<Tensor>(I663_index);
  auto tensor398 = vector<shared_ptr<Tensor>>{I72, v2_, I663};
  auto task398 = make_shared<Task398>(tensor398, pindex);
  task344->add_dep(task398);
  task398->add_dep(task108);
  residualq->add_task(task398);

  auto tensor399 = vector<shared_ptr<Tensor>>{I663, Gamma4_(), t2};
  auto task399 = make_shared<Task399>(tensor399, pindex);
  task398->add_dep(task399);
  task399->add_dep(task108);
  residualq->add_task(task399);

  vector<IndexRange> I666_index = {closed_, active_, active_, active_};
  auto I666 = make_shared<Tensor>(I666_index);
  auto tensor400 = vector<shared_ptr<Tensor>>{I72, v2_, I666};
  auto task400 = make_shared<Task400>(tensor400, pindex);
  task344->add_dep(task400);
  task400->add_dep(task108);
  residualq->add_task(task400);

  auto tensor401 = vector<shared_ptr<Tensor>>{I666, Gamma24_(), t2};
  auto task401 = make_shared<Task401>(tensor401, pindex);
  task400->add_dep(task401);
  task401->add_dep(task108);
  residualq->add_task(task401);

  vector<IndexRange> I669_index = {closed_, active_, active_, active_};
  auto I669 = make_shared<Tensor>(I669_index);
  auto tensor402 = vector<shared_ptr<Tensor>>{I72, t2, I669};
  auto task402 = make_shared<Task402>(tensor402, pindex);
  task344->add_dep(task402);
  task402->add_dep(task108);
  residualq->add_task(task402);

  auto tensor403 = vector<shared_ptr<Tensor>>{I669, Gamma220_(), v2_};
  auto task403 = make_shared<Task403>(tensor403, pindex);
  task402->add_dep(task403);
  task403->add_dep(task108);
  residualq->add_task(task403);

  auto tensor404 = vector<shared_ptr<Tensor>>{I669, Gamma222_(), v2_};
  auto task404 = make_shared<Task404>(tensor404, pindex);
  task402->add_dep(task404);
  task404->add_dep(task108);
  residualq->add_task(task404);

  vector<IndexRange> I672_index = {closed_, active_, active_, active_};
  auto I672 = make_shared<Tensor>(I672_index);
  auto tensor405 = vector<shared_ptr<Tensor>>{I72, t2, I672};
  auto task405 = make_shared<Task405>(tensor405, pindex);
  task344->add_dep(task405);
  task405->add_dep(task108);
  residualq->add_task(task405);

  auto tensor406 = vector<shared_ptr<Tensor>>{I672, Gamma221_(), v2_};
  auto task406 = make_shared<Task406>(tensor406, pindex);
  task405->add_dep(task406);
  task406->add_dep(task108);
  residualq->add_task(task406);

  auto tensor407 = vector<shared_ptr<Tensor>>{I672, Gamma104_(), v2_};
  auto task407 = make_shared<Task407>(tensor407, pindex);
  task405->add_dep(task407);
  task407->add_dep(task108);
  residualq->add_task(task407);

  vector<IndexRange> I699_index = {closed_, closed_, active_, active_, active_, active_};
  auto I699 = make_shared<Tensor>(I699_index);
  auto tensor408 = vector<shared_ptr<Tensor>>{I72, t2, I699};
  auto task408 = make_shared<Task408>(tensor408, pindex);
  task344->add_dep(task408);
  task408->add_dep(task108);
  residualq->add_task(task408);

  vector<IndexRange> I700_index = {closed_, closed_, active_, active_};
  auto I700 = make_shared<Tensor>(I700_index);
  auto tensor409 = vector<shared_ptr<Tensor>>{I699, Gamma230_(), I700};
  auto task409 = make_shared<Task409>(tensor409, pindex);
  task408->add_dep(task409);
  task409->add_dep(task108);
  residualq->add_task(task409);

  auto tensor410 = vector<shared_ptr<Tensor>>{I700, v2_};
  auto task410 = make_shared<Task410>(tensor410, pindex);
  task409->add_dep(task410);
  task410->add_dep(task108);
  residualq->add_task(task410);

  auto tensor411 = vector<shared_ptr<Tensor>>{I699, Gamma232_(), v2_};
  auto task411 = make_shared<Task411>(tensor411, pindex);
  task408->add_dep(task411);
  task411->add_dep(task108);
  residualq->add_task(task411);

  auto tensor412 = vector<shared_ptr<Tensor>>{I699, Gamma234_(), v2_};
  auto task412 = make_shared<Task412>(tensor412, pindex);
  task408->add_dep(task412);
  task412->add_dep(task108);
  residualq->add_task(task412);

  vector<IndexRange> I702_index = {virt_, active_, active_, active_, closed_, active_};
  auto I702 = make_shared<Tensor>(I702_index);
  auto tensor413 = vector<shared_ptr<Tensor>>{I72, Gamma230_(), I702};
  auto task413 = make_shared<Task413>(tensor413, pindex);
  task344->add_dep(task413);
  task413->add_dep(task108);
  residualq->add_task(task413);

  vector<IndexRange> I703_index = {virt_, virt_, active_, active_};
  auto I703 = make_shared<Tensor>(I703_index);
  auto tensor414 = vector<shared_ptr<Tensor>>{I702, t2, I703};
  auto task414 = make_shared<Task414>(tensor414, pindex);
  task413->add_dep(task414);
  task414->add_dep(task108);
  residualq->add_task(task414);

  auto tensor415 = vector<shared_ptr<Tensor>>{I703, v2_};
  auto task415 = make_shared<Task415>(tensor415, pindex);
  task414->add_dep(task415);
  task415->add_dep(task108);
  residualq->add_task(task415);

  vector<IndexRange> I708_index = {active_, active_, virt_, active_, closed_, active_};
  auto I708 = make_shared<Tensor>(I708_index);
  auto tensor416 = vector<shared_ptr<Tensor>>{I72, Gamma233_(), I708};
  auto task416 = make_shared<Task416>(tensor416, pindex);
  task344->add_dep(task416);
  task416->add_dep(task108);
  residualq->add_task(task416);

  auto tensor417 = vector<shared_ptr<Tensor>>{I708, t2, v2_};
  auto task417 = make_shared<Task417>(tensor417, pindex);
  task416->add_dep(task417);
  task417->add_dep(task108);
  residualq->add_task(task417);

  vector<IndexRange> I714_index = {active_, virt_, active_, active_, closed_, active_};
  auto I714 = make_shared<Tensor>(I714_index);
  auto tensor418 = vector<shared_ptr<Tensor>>{I72, Gamma235_(), I714};
  auto task418 = make_shared<Task418>(tensor418, pindex);
  task344->add_dep(task418);
  task418->add_dep(task108);
  residualq->add_task(task418);

  auto tensor419 = vector<shared_ptr<Tensor>>{I714, t2, v2_};
  auto task419 = make_shared<Task419>(tensor419, pindex);
  task418->add_dep(task419);
  task419->add_dep(task108);
  residualq->add_task(task419);

  vector<IndexRange> I729_index = {closed_, closed_, active_, active_, active_, active_};
  auto I729 = make_shared<Tensor>(I729_index);
  auto tensor420 = vector<shared_ptr<Tensor>>{I72, t2, I729};
  auto task420 = make_shared<Task420>(tensor420, pindex);
  task344->add_dep(task420);
  task420->add_dep(task108);
  residualq->add_task(task420);

  vector<IndexRange> I730_index = {closed_, closed_, active_, active_};
  auto I730 = make_shared<Tensor>(I730_index);
  auto tensor421 = vector<shared_ptr<Tensor>>{I729, Gamma240_(), I730};
  auto task421 = make_shared<Task421>(tensor421, pindex);
  task420->add_dep(task421);
  task421->add_dep(task108);
  residualq->add_task(task421);

  auto tensor422 = vector<shared_ptr<Tensor>>{I730, v2_};
  auto task422 = make_shared<Task422>(tensor422, pindex);
  task421->add_dep(task422);
  task422->add_dep(task108);
  residualq->add_task(task422);

  auto tensor423 = vector<shared_ptr<Tensor>>{I729, Gamma24_(), v2_};
  auto task423 = make_shared<Task423>(tensor423, pindex);
  task420->add_dep(task423);
  task423->add_dep(task108);
  residualq->add_task(task423);

  auto tensor424 = vector<shared_ptr<Tensor>>{I729, Gamma244_(), v2_};
  auto task424 = make_shared<Task424>(tensor424, pindex);
  task420->add_dep(task424);
  task424->add_dep(task108);
  residualq->add_task(task424);

  vector<IndexRange> I732_index = {virt_, active_, active_, closed_, active_, active_};
  auto I732 = make_shared<Tensor>(I732_index);
  auto tensor425 = vector<shared_ptr<Tensor>>{I72, Gamma240_(), I732};
  auto task425 = make_shared<Task425>(tensor425, pindex);
  task344->add_dep(task425);
  task425->add_dep(task108);
  residualq->add_task(task425);

  vector<IndexRange> I733_index = {virt_, virt_, active_, active_};
  auto I733 = make_shared<Tensor>(I733_index);
  auto tensor426 = vector<shared_ptr<Tensor>>{I732, t2, I733};
  auto task426 = make_shared<Task426>(tensor426, pindex);
  task425->add_dep(task426);
  task426->add_dep(task108);
  residualq->add_task(task426);

  auto tensor427 = vector<shared_ptr<Tensor>>{I733, v2_};
  auto task427 = make_shared<Task427>(tensor427, pindex);
  task426->add_dep(task427);
  task427->add_dep(task108);
  residualq->add_task(task427);

  vector<IndexRange> I738_index = {active_, active_, virt_, closed_, active_, active_};
  auto I738 = make_shared<Tensor>(I738_index);
  auto tensor428 = vector<shared_ptr<Tensor>>{I72, Gamma24_(), I738};
  auto task428 = make_shared<Task428>(tensor428, pindex);
  task344->add_dep(task428);
  task428->add_dep(task108);
  residualq->add_task(task428);

  auto tensor429 = vector<shared_ptr<Tensor>>{I738, t2, v2_};
  auto task429 = make_shared<Task429>(tensor429, pindex);
  task428->add_dep(task429);
  task429->add_dep(task108);
  residualq->add_task(task429);

  vector<IndexRange> I744_index = {active_, virt_, active_, closed_, active_, active_};
  auto I744 = make_shared<Tensor>(I744_index);
  auto tensor430 = vector<shared_ptr<Tensor>>{I72, Gamma244_(), I744};
  auto task430 = make_shared<Task430>(tensor430, pindex);
  task344->add_dep(task430);
  task430->add_dep(task108);
  residualq->add_task(task430);

  auto tensor431 = vector<shared_ptr<Tensor>>{I744, t2, v2_};
  auto task431 = make_shared<Task431>(tensor431, pindex);
  task430->add_dep(task431);
  task431->add_dep(task108);
  residualq->add_task(task431);

  vector<IndexRange> I759_index = {closed_, active_, active_, active_, active_, active_};
  auto I759 = make_shared<Tensor>(I759_index);
  auto tensor432 = vector<shared_ptr<Tensor>>{I72, t2, I759};
  auto task432 = make_shared<Task432>(tensor432, pindex);
  task344->add_dep(task432);
  task432->add_dep(task108);
  residualq->add_task(task432);

  auto tensor433 = vector<shared_ptr<Tensor>>{I759, Gamma250_(), v2_};
  auto task433 = make_shared<Task433>(tensor433, pindex);
  task432->add_dep(task433);
  task433->add_dep(task108);
  residualq->add_task(task433);

  auto tensor434 = vector<shared_ptr<Tensor>>{I759, Gamma251_(), v2_};
  auto task434 = make_shared<Task434>(tensor434, pindex);
  task432->add_dep(task434);
  task434->add_dep(task108);
  residualq->add_task(task434);

  vector<IndexRange> I765_index = {virt_, active_, active_, active_};
  auto I765 = make_shared<Tensor>(I765_index);
  auto tensor435 = vector<shared_ptr<Tensor>>{I72, v2_, I765};
  auto task435 = make_shared<Task435>(tensor435, pindex);
  task344->add_dep(task435);
  task435->add_dep(task108);
  residualq->add_task(task435);

  auto tensor436 = vector<shared_ptr<Tensor>>{I765, Gamma252_(), t2};
  auto task436 = make_shared<Task436>(tensor436, pindex);
  task435->add_dep(task436);
  task436->add_dep(task108);
  residualq->add_task(task436);

  vector<IndexRange> I768_index = {virt_, active_, active_, active_};
  auto I768 = make_shared<Tensor>(I768_index);
  auto tensor437 = vector<shared_ptr<Tensor>>{I72, v2_, I768};
  auto task437 = make_shared<Task437>(tensor437, pindex);
  task344->add_dep(task437);
  task437->add_dep(task108);
  residualq->add_task(task437);

  auto tensor438 = vector<shared_ptr<Tensor>>{I768, Gamma31_(), t2};
  auto task438 = make_shared<Task438>(tensor438, pindex);
  task437->add_dep(task438);
  task438->add_dep(task108);
  residualq->add_task(task438);

  vector<IndexRange> I807_index = {virt_, active_, active_, active_};
  auto I807 = make_shared<Tensor>(I807_index);
  auto tensor439 = vector<shared_ptr<Tensor>>{I72, t2, I807};
  auto task439 = make_shared<Task439>(tensor439, pindex);
  task344->add_dep(task439);
  task439->add_dep(task108);
  residualq->add_task(task439);

  auto tensor440 = vector<shared_ptr<Tensor>>{I807, Gamma240_(), v2_};
  auto task440 = make_shared<Task440>(tensor440, pindex);
  task439->add_dep(task440);
  task440->add_dep(task108);
  residualq->add_task(task440);

  auto tensor441 = vector<shared_ptr<Tensor>>{I807, Gamma252_(), v2_};
  auto task441 = make_shared<Task441>(tensor441, pindex);
  task439->add_dep(task441);
  task441->add_dep(task108);
  residualq->add_task(task441);

  vector<IndexRange> I810_index = {virt_, active_, active_, active_};
  auto I810 = make_shared<Tensor>(I810_index);
  auto tensor442 = vector<shared_ptr<Tensor>>{I72, t2, I810};
  auto task442 = make_shared<Task442>(tensor442, pindex);
  task344->add_dep(task442);
  task442->add_dep(task108);
  residualq->add_task(task442);

  auto tensor443 = vector<shared_ptr<Tensor>>{I810, Gamma230_(), v2_};
  auto task443 = make_shared<Task443>(tensor443, pindex);
  task442->add_dep(task443);
  task443->add_dep(task108);
  residualq->add_task(task443);

  auto tensor444 = vector<shared_ptr<Tensor>>{I810, Gamma31_(), v2_};
  auto task444 = make_shared<Task444>(tensor444, pindex);
  task442->add_dep(task444);
  task444->add_dep(task108);
  residualq->add_task(task444);

  vector<IndexRange> I837_index = {closed_, active_, active_, active_, virt_, active_};
  auto I837 = make_shared<Tensor>(I837_index);
  auto tensor445 = vector<shared_ptr<Tensor>>{I72, Gamma276_(), I837};
  auto task445 = make_shared<Task445>(tensor445, pindex);
  task344->add_dep(task445);
  task445->add_dep(task108);
  residualq->add_task(task445);

  auto tensor446 = vector<shared_ptr<Tensor>>{I837, t2, v2_};
  auto task446 = make_shared<Task446>(tensor446, pindex);
  task445->add_dep(task446);
  task446->add_dep(task108);
  residualq->add_task(task446);

  auto tensor447 = vector<shared_ptr<Tensor>>{I72, Gamma568_(), t2};
  auto task447 = make_shared<Task447>(tensor447, pindex);
  task344->add_dep(task447);
  task447->add_dep(task108);
  residualq->add_task(task447);

  auto tensor448 = vector<shared_ptr<Tensor>>{I72, Gamma569_(), t2};
  auto task448 = make_shared<Task448>(tensor448, pindex);
  task344->add_dep(task448);
  task448->add_dep(task108);
  residualq->add_task(task448);

  auto tensor449 = vector<shared_ptr<Tensor>>{I72, Gamma572_(), t2};
  auto task449 = make_shared<Task449>(tensor449, pindex);
  task344->add_dep(task449);
  task449->add_dep(task108);
  residualq->add_task(task449);

  auto tensor450 = vector<shared_ptr<Tensor>>{I72, Gamma573_(), t2};
  auto task450 = make_shared<Task450>(tensor450, pindex);
  task344->add_dep(task450);
  task450->add_dep(task108);
  residualq->add_task(task450);
}

#endif
