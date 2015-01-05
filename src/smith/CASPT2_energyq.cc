//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_energyqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_energyq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto energyq = make_shared<Queue>();
  vector<IndexRange> I348_index;
  auto I348 = make_shared<Tensor>(I348_index);
  vector<IndexRange> I349_index = {active_, active_, active_, active_};
  auto I349 = make_shared<Tensor>(I349_index);
  vector<shared_ptr<Tensor>> tensor331 = {I348, Gamma94_(), I349};
  auto task331 = make_shared<Task331>(tensor331, pindex);
  energyq->add_task(task331);

  vector<IndexRange> I350_index = {closed_, active_, closed_, active_};
  auto I350 = make_shared<Tensor>(I350_index);
  vector<shared_ptr<Tensor>> tensor332 = {I349, t2, I350};
  auto task332 = make_shared<Task332>(tensor332, pindex);
  task331->add_dep(task332);
  energyq->add_task(task332);

  vector<shared_ptr<Tensor>> tensor333 = {I350, v2_};
  auto task333 = make_shared<Task333>(tensor333, pindex);
  task332->add_dep(task333);
  energyq->add_task(task333);

  vector<IndexRange> I352_index = {active_, active_, active_, closed_};
  auto I352 = make_shared<Tensor>(I352_index);
  vector<shared_ptr<Tensor>> tensor334 = {I348, t2, I352};
  auto task334 = make_shared<Task334>(tensor334, pindex);
  task331->add_dep(task334);
  energyq->add_task(task334);

  vector<IndexRange> I353_index = {active_, active_, active_, active_, active_, active_};
  auto I353 = make_shared<Tensor>(I353_index);
  vector<shared_ptr<Tensor>> tensor335 = {I352, v2_, I353};
  auto task335 = make_shared<Task335>(tensor335, pindex);
  task334->add_dep(task335);
  energyq->add_task(task335);

  vector<shared_ptr<Tensor>> tensor336 = {I353, Gamma107_()};
  auto task336 = make_shared<Task336>(tensor336, pindex);
  task335->add_dep(task336);
  energyq->add_task(task336);

  vector<IndexRange> I355_index = {closed_, active_, active_, active_};
  auto I355 = make_shared<Tensor>(I355_index);
  vector<shared_ptr<Tensor>> tensor337 = {I348, v2_, I355};
  auto task337 = make_shared<Task337>(tensor337, pindex);
  task331->add_dep(task337);
  energyq->add_task(task337);

  vector<IndexRange> I356_index = {active_, closed_, active_, active_};
  auto I356 = make_shared<Tensor>(I356_index);
  vector<shared_ptr<Tensor>> tensor338 = {I355, Gamma6_(), I356};
  auto task338 = make_shared<Task338>(tensor338, pindex);
  task337->add_dep(task338);
  energyq->add_task(task338);

  vector<shared_ptr<Tensor>> tensor339 = {I356, t2};
  auto task339 = make_shared<Task339>(tensor339, pindex);
  task338->add_dep(task339);
  energyq->add_task(task339);

  vector<IndexRange> I358_index = {active_, active_};
  auto I358 = make_shared<Tensor>(I358_index);
  vector<shared_ptr<Tensor>> tensor340 = {I348, Gamma16_(), I358};
  auto task340 = make_shared<Task340>(tensor340, pindex);
  task331->add_dep(task340);
  energyq->add_task(task340);

  vector<IndexRange> I359_index = {active_, closed_, virt_, closed_};
  auto I359 = make_shared<Tensor>(I359_index);
  vector<shared_ptr<Tensor>> tensor341 = {I358, v2_, I359};
  auto task341 = make_shared<Task341>(tensor341, pindex);
  task340->add_dep(task341);
  energyq->add_task(task341);

  vector<shared_ptr<Tensor>> tensor342 = {I359, t2};
  auto task342 = make_shared<Task342>(tensor342, pindex);
  task341->add_dep(task342);
  energyq->add_task(task342);

  vector<IndexRange> I362_index = {active_, closed_, virt_, closed_};
  auto I362 = make_shared<Tensor>(I362_index);
  vector<shared_ptr<Tensor>> tensor343 = {I358, v2_, I362};
  auto task343 = make_shared<Task343>(tensor343, pindex);
  task340->add_dep(task343);
  energyq->add_task(task343);

  vector<shared_ptr<Tensor>> tensor344 = {I362, t2};
  auto task344 = make_shared<Task344>(tensor344, pindex);
  task343->add_dep(task344);
  energyq->add_task(task344);

  vector<IndexRange> I364_index = {active_, active_, active_, active_};
  auto I364 = make_shared<Tensor>(I364_index);
  vector<shared_ptr<Tensor>> tensor345 = {I348, Gamma35_(), I364};
  auto task345 = make_shared<Task345>(tensor345, pindex);
  task331->add_dep(task345);
  energyq->add_task(task345);

  vector<IndexRange> I365_index = {active_, closed_, virt_, active_};
  auto I365 = make_shared<Tensor>(I365_index);
  vector<shared_ptr<Tensor>> tensor346 = {I364, v2_, I365};
  auto task346 = make_shared<Task346>(tensor346, pindex);
  task345->add_dep(task346);
  energyq->add_task(task346);

  vector<shared_ptr<Tensor>> tensor347 = {I365, t2};
  auto task347 = make_shared<Task347>(tensor347, pindex);
  task346->add_dep(task347);
  energyq->add_task(task347);

  vector<IndexRange> I374_index = {active_, closed_, virt_, active_};
  auto I374 = make_shared<Tensor>(I374_index);
  vector<shared_ptr<Tensor>> tensor348 = {I364, v2_, I374};
  auto task348 = make_shared<Task348>(tensor348, pindex);
  task345->add_dep(task348);
  energyq->add_task(task348);

  vector<shared_ptr<Tensor>> tensor349 = {I374, t2};
  auto task349 = make_shared<Task349>(tensor349, pindex);
  task348->add_dep(task349);
  energyq->add_task(task349);

  vector<IndexRange> I377_index = {active_, active_, virt_, closed_};
  auto I377 = make_shared<Tensor>(I377_index);
  vector<shared_ptr<Tensor>> tensor350 = {I364, v2_, I377};
  auto task350 = make_shared<Task350>(tensor350, pindex);
  task345->add_dep(task350);
  energyq->add_task(task350);

  vector<shared_ptr<Tensor>> tensor351 = {I377, t2};
  auto task351 = make_shared<Task351>(tensor351, pindex);
  task350->add_dep(task351);
  energyq->add_task(task351);

  vector<IndexRange> I383_index = {active_, active_, virt_, closed_};
  auto I383 = make_shared<Tensor>(I383_index);
  vector<shared_ptr<Tensor>> tensor352 = {I364, v2_, I383};
  auto task352 = make_shared<Task352>(tensor352, pindex);
  task345->add_dep(task352);
  energyq->add_task(task352);

  vector<shared_ptr<Tensor>> tensor353 = {I383, t2};
  auto task353 = make_shared<Task353>(tensor353, pindex);
  task352->add_dep(task353);
  energyq->add_task(task353);

  vector<IndexRange> I386_index = {active_, active_, virt_, closed_};
  auto I386 = make_shared<Tensor>(I386_index);
  vector<shared_ptr<Tensor>> tensor354 = {I364, v2_, I386};
  auto task354 = make_shared<Task354>(tensor354, pindex);
  task345->add_dep(task354);
  energyq->add_task(task354);

  vector<shared_ptr<Tensor>> tensor355 = {I386, t2};
  auto task355 = make_shared<Task355>(tensor355, pindex);
  task354->add_dep(task355);
  energyq->add_task(task355);

  vector<IndexRange> I367_index = {active_, active_, active_, active_};
  auto I367 = make_shared<Tensor>(I367_index);
  vector<shared_ptr<Tensor>> tensor356 = {I348, Gamma29_(), I367};
  auto task356 = make_shared<Task356>(tensor356, pindex);
  task331->add_dep(task356);
  energyq->add_task(task356);

  vector<IndexRange> I368_index = {active_, closed_, virt_, active_};
  auto I368 = make_shared<Tensor>(I368_index);
  vector<shared_ptr<Tensor>> tensor357 = {I367, v2_, I368};
  auto task357 = make_shared<Task357>(tensor357, pindex);
  task356->add_dep(task357);
  energyq->add_task(task357);

  vector<shared_ptr<Tensor>> tensor358 = {I368, t2};
  auto task358 = make_shared<Task358>(tensor358, pindex);
  task357->add_dep(task358);
  energyq->add_task(task358);

  vector<IndexRange> I370_index = {active_, active_, active_, active_};
  auto I370 = make_shared<Tensor>(I370_index);
  vector<shared_ptr<Tensor>> tensor359 = {I348, Gamma32_(), I370};
  auto task359 = make_shared<Task359>(tensor359, pindex);
  task331->add_dep(task359);
  energyq->add_task(task359);

  vector<IndexRange> I371_index = {active_, closed_, virt_, active_};
  auto I371 = make_shared<Tensor>(I371_index);
  vector<shared_ptr<Tensor>> tensor360 = {I370, v2_, I371};
  auto task360 = make_shared<Task360>(tensor360, pindex);
  task359->add_dep(task360);
  energyq->add_task(task360);

  vector<shared_ptr<Tensor>> tensor361 = {I371, t2};
  auto task361 = make_shared<Task361>(tensor361, pindex);
  task360->add_dep(task361);
  energyq->add_task(task361);

  vector<IndexRange> I379_index = {active_, active_, active_, active_};
  auto I379 = make_shared<Tensor>(I379_index);
  vector<shared_ptr<Tensor>> tensor362 = {I348, Gamma7_(), I379};
  auto task362 = make_shared<Task362>(tensor362, pindex);
  task331->add_dep(task362);
  energyq->add_task(task362);

  vector<IndexRange> I380_index = {closed_, active_, active_, virt_};
  auto I380 = make_shared<Tensor>(I380_index);
  vector<shared_ptr<Tensor>> tensor363 = {I379, t2, I380};
  auto task363 = make_shared<Task363>(tensor363, pindex);
  task362->add_dep(task363);
  energyq->add_task(task363);

  vector<shared_ptr<Tensor>> tensor364 = {I380, v2_};
  auto task364 = make_shared<Task364>(tensor364, pindex);
  task363->add_dep(task364);
  energyq->add_task(task364);

  vector<IndexRange> I388_index = {active_, active_, active_, active_, active_, active_};
  auto I388 = make_shared<Tensor>(I388_index);
  vector<shared_ptr<Tensor>> tensor365 = {I348, Gamma59_(), I388};
  auto task365 = make_shared<Task365>(tensor365, pindex);
  task331->add_dep(task365);
  energyq->add_task(task365);

  vector<IndexRange> I389_index = {active_, active_, virt_, active_};
  auto I389 = make_shared<Tensor>(I389_index);
  vector<shared_ptr<Tensor>> tensor366 = {I388, v2_, I389};
  auto task366 = make_shared<Task366>(tensor366, pindex);
  task365->add_dep(task366);
  energyq->add_task(task366);

  vector<shared_ptr<Tensor>> tensor367 = {I389, t2};
  auto task367 = make_shared<Task367>(tensor367, pindex);
  task366->add_dep(task367);
  energyq->add_task(task367);

  vector<IndexRange> I391_index = {active_, active_, active_, active_, active_, active_};
  auto I391 = make_shared<Tensor>(I391_index);
  vector<shared_ptr<Tensor>> tensor368 = {I348, Gamma57_(), I391};
  auto task368 = make_shared<Task368>(tensor368, pindex);
  task331->add_dep(task368);
  energyq->add_task(task368);

  vector<IndexRange> I392_index = {active_, active_, virt_, active_};
  auto I392 = make_shared<Tensor>(I392_index);
  vector<shared_ptr<Tensor>> tensor369 = {I391, v2_, I392};
  auto task369 = make_shared<Task369>(tensor369, pindex);
  task368->add_dep(task369);
  energyq->add_task(task369);

  vector<shared_ptr<Tensor>> tensor370 = {I392, t2};
  auto task370 = make_shared<Task370>(tensor370, pindex);
  task369->add_dep(task370);
  energyq->add_task(task370);

  vector<IndexRange> I394_index = {closed_, virt_, closed_, virt_};
  auto I394 = make_shared<Tensor>(I394_index);
  vector<shared_ptr<Tensor>> tensor371 = {I348, t2, I394};
  auto task371 = make_shared<Task371>(tensor371, pindex);
  task331->add_dep(task371);
  energyq->add_task(task371);

  vector<shared_ptr<Tensor>> tensor372 = {I394, v2_};
  auto task372 = make_shared<Task372>(tensor372, pindex);
  task371->add_dep(task372);
  energyq->add_task(task372);

  vector<IndexRange> I396_index = {virt_, closed_, virt_, closed_};
  auto I396 = make_shared<Tensor>(I396_index);
  vector<shared_ptr<Tensor>> tensor373 = {I348, v2_, I396};
  auto task373 = make_shared<Task373>(tensor373, pindex);
  task331->add_dep(task373);
  energyq->add_task(task373);

  vector<shared_ptr<Tensor>> tensor374 = {I396, t2};
  auto task374 = make_shared<Task374>(tensor374, pindex);
  task373->add_dep(task374);
  energyq->add_task(task374);

  vector<IndexRange> I398_index = {active_, active_};
  auto I398 = make_shared<Tensor>(I398_index);
  vector<shared_ptr<Tensor>> tensor375 = {I348, Gamma38_(), I398};
  auto task375 = make_shared<Task375>(tensor375, pindex);
  task331->add_dep(task375);
  energyq->add_task(task375);

  vector<IndexRange> I399_index = {virt_, closed_, virt_, active_};
  auto I399 = make_shared<Tensor>(I399_index);
  vector<shared_ptr<Tensor>> tensor376 = {I398, v2_, I399};
  auto task376 = make_shared<Task376>(tensor376, pindex);
  task375->add_dep(task376);
  energyq->add_task(task376);

  vector<shared_ptr<Tensor>> tensor377 = {I399, t2};
  auto task377 = make_shared<Task377>(tensor377, pindex);
  task376->add_dep(task377);
  energyq->add_task(task377);

  vector<IndexRange> I402_index = {virt_, closed_, virt_, active_};
  auto I402 = make_shared<Tensor>(I402_index);
  vector<shared_ptr<Tensor>> tensor378 = {I398, v2_, I402};
  auto task378 = make_shared<Task378>(tensor378, pindex);
  task375->add_dep(task378);
  energyq->add_task(task378);

  vector<shared_ptr<Tensor>> tensor379 = {I402, t2};
  auto task379 = make_shared<Task379>(tensor379, pindex);
  task378->add_dep(task379);
  energyq->add_task(task379);

  vector<IndexRange> I411_index = {closed_, virt_};
  auto I411 = make_shared<Tensor>(I411_index);
  vector<shared_ptr<Tensor>> tensor380 = {I398, t2, I411};
  auto task380 = make_shared<Task380>(tensor380, pindex);
  task375->add_dep(task380);
  energyq->add_task(task380);

  vector<shared_ptr<Tensor>> tensor381 = {I411, h1_};
  auto task381 = make_shared<Task381>(tensor381, pindex);
  task380->add_dep(task381);
  energyq->add_task(task381);

  vector<IndexRange> I414_index = {active_, active_, virt_, closed_};
  auto I414 = make_shared<Tensor>(I414_index);
  vector<shared_ptr<Tensor>> tensor382 = {I398, h1_, I414};
  auto task382 = make_shared<Task382>(tensor382, pindex);
  task375->add_dep(task382);
  energyq->add_task(task382);

  vector<shared_ptr<Tensor>> tensor383 = {I414, t2};
  auto task383 = make_shared<Task383>(tensor383, pindex);
  task382->add_dep(task383);
  energyq->add_task(task383);

  vector<IndexRange> I404_index = {active_, active_, active_, active_};
  auto I404 = make_shared<Tensor>(I404_index);
  vector<shared_ptr<Tensor>> tensor384 = {I348, Gamma60_(), I404};
  auto task384 = make_shared<Task384>(tensor384, pindex);
  task331->add_dep(task384);
  energyq->add_task(task384);

  vector<IndexRange> I405_index = {virt_, active_, virt_, active_};
  auto I405 = make_shared<Tensor>(I405_index);
  vector<shared_ptr<Tensor>> tensor385 = {I404, v2_, I405};
  auto task385 = make_shared<Task385>(tensor385, pindex);
  task384->add_dep(task385);
  energyq->add_task(task385);

  vector<shared_ptr<Tensor>> tensor386 = {I405, t2};
  auto task386 = make_shared<Task386>(tensor386, pindex);
  task385->add_dep(task386);
  energyq->add_task(task386);

  vector<IndexRange> I417_index = {active_, active_, virt_, active_};
  auto I417 = make_shared<Tensor>(I417_index);
  vector<shared_ptr<Tensor>> tensor387 = {I404, h1_, I417};
  auto task387 = make_shared<Task387>(tensor387, pindex);
  task384->add_dep(task387);
  energyq->add_task(task387);

  vector<shared_ptr<Tensor>> tensor388 = {I417, t2};
  auto task388 = make_shared<Task388>(tensor388, pindex);
  task387->add_dep(task388);
  energyq->add_task(task388);

  vector<IndexRange> I407_index = {closed_, active_};
  auto I407 = make_shared<Tensor>(I407_index);
  vector<shared_ptr<Tensor>> tensor389 = {I348, h1_, I407};
  auto task389 = make_shared<Task389>(tensor389, pindex);
  task331->add_dep(task389);
  energyq->add_task(task389);

  vector<IndexRange> I408_index = {active_, closed_, active_, active_};
  auto I408 = make_shared<Tensor>(I408_index);
  vector<shared_ptr<Tensor>> tensor390 = {I407, Gamma7_(), I408};
  auto task390 = make_shared<Task390>(tensor390, pindex);
  task389->add_dep(task390);
  energyq->add_task(task390);

  vector<shared_ptr<Tensor>> tensor391 = {I408, t2};
  auto task391 = make_shared<Task391>(tensor391, pindex);
  task390->add_dep(task391);
  energyq->add_task(task391);

  return energyq;
}


