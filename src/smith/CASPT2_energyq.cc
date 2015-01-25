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

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_energyq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto energyq = make_shared<Queue>();
  vector<IndexRange> I334_index;
  auto I334 = make_shared<Tensor>(I334_index);
  vector<IndexRange> I335_index = {active_, active_, active_, active_};
  auto I335 = make_shared<Tensor>(I335_index);
  vector<shared_ptr<Tensor>> tensor320 = {I334, Gamma92_(), I335};
  auto task320 = make_shared<Task320>(tensor320, pindex);
  energyq->add_task(task320);

  vector<IndexRange> I336_index = {closed_, active_, closed_, active_};
  auto I336 = make_shared<Tensor>(I336_index);
  vector<shared_ptr<Tensor>> tensor321 = {I335, t2, I336};
  auto task321 = make_shared<Task321>(tensor321, pindex);
  task320->add_dep(task321);
  energyq->add_task(task321);

  vector<shared_ptr<Tensor>> tensor322 = {I336, v2_};
  auto task322 = make_shared<Task322>(tensor322, pindex);
  task321->add_dep(task322);
  energyq->add_task(task322);

  vector<IndexRange> I338_index = {closed_, active_, active_, active_};
  auto I338 = make_shared<Tensor>(I338_index);
  vector<shared_ptr<Tensor>> tensor323 = {I334, v2_, I338};
  auto task323 = make_shared<Task323>(tensor323, pindex);
  task320->add_dep(task323);
  energyq->add_task(task323);

  vector<IndexRange> I339_index = {active_, closed_, active_, active_};
  auto I339 = make_shared<Tensor>(I339_index);
  vector<shared_ptr<Tensor>> tensor324 = {I338, Gamma105_(), I339};
  auto task324 = make_shared<Task324>(tensor324, pindex);
  task323->add_dep(task324);
  energyq->add_task(task324);

  vector<shared_ptr<Tensor>> tensor325 = {I339, t2};
  auto task325 = make_shared<Task325>(tensor325, pindex);
  task324->add_dep(task325);
  energyq->add_task(task325);

  vector<IndexRange> I341_index = {closed_, active_, active_, active_};
  auto I341 = make_shared<Tensor>(I341_index);
  vector<shared_ptr<Tensor>> tensor326 = {I334, t2, I341};
  auto task326 = make_shared<Task326>(tensor326, pindex);
  task320->add_dep(task326);
  energyq->add_task(task326);

  vector<IndexRange> I342_index = {active_, active_, closed_, active_};
  auto I342 = make_shared<Tensor>(I342_index);
  vector<shared_ptr<Tensor>> tensor327 = {I341, Gamma6_(), I342};
  auto task327 = make_shared<Task327>(tensor327, pindex);
  task326->add_dep(task327);
  energyq->add_task(task327);

  vector<shared_ptr<Tensor>> tensor328 = {I342, v2_};
  auto task328 = make_shared<Task328>(tensor328, pindex);
  task327->add_dep(task328);
  energyq->add_task(task328);

  vector<IndexRange> I344_index = {active_, active_};
  auto I344 = make_shared<Tensor>(I344_index);
  vector<shared_ptr<Tensor>> tensor329 = {I334, Gamma16_(), I344};
  auto task329 = make_shared<Task329>(tensor329, pindex);
  task320->add_dep(task329);
  energyq->add_task(task329);

  vector<IndexRange> I345_index = {active_, closed_, virt_, closed_};
  auto I345 = make_shared<Tensor>(I345_index);
  vector<shared_ptr<Tensor>> tensor330 = {I344, v2_, I345};
  auto task330 = make_shared<Task330>(tensor330, pindex);
  task329->add_dep(task330);
  energyq->add_task(task330);

  vector<shared_ptr<Tensor>> tensor331 = {I345, t2};
  auto task331 = make_shared<Task331>(tensor331, pindex);
  task330->add_dep(task331);
  energyq->add_task(task331);

  vector<IndexRange> I348_index = {closed_, active_, closed_, virt_};
  auto I348 = make_shared<Tensor>(I348_index);
  vector<shared_ptr<Tensor>> tensor332 = {I344, t2, I348};
  auto task332 = make_shared<Task332>(tensor332, pindex);
  task329->add_dep(task332);
  energyq->add_task(task332);

  vector<shared_ptr<Tensor>> tensor333 = {I348, v2_};
  auto task333 = make_shared<Task333>(tensor333, pindex);
  task332->add_dep(task333);
  energyq->add_task(task333);

  vector<IndexRange> I350_index = {active_, active_, active_, active_};
  auto I350 = make_shared<Tensor>(I350_index);
  vector<shared_ptr<Tensor>> tensor334 = {I334, Gamma35_(), I350};
  auto task334 = make_shared<Task334>(tensor334, pindex);
  task320->add_dep(task334);
  energyq->add_task(task334);

  vector<IndexRange> I351_index = {closed_, virt_, active_, active_};
  auto I351 = make_shared<Tensor>(I351_index);
  vector<shared_ptr<Tensor>> tensor335 = {I350, t2, I351};
  auto task335 = make_shared<Task335>(tensor335, pindex);
  task334->add_dep(task335);
  energyq->add_task(task335);

  vector<shared_ptr<Tensor>> tensor336 = {I351, v2_};
  auto task336 = make_shared<Task336>(tensor336, pindex);
  task335->add_dep(task336);
  energyq->add_task(task336);

  vector<IndexRange> I360_index = {active_, closed_, virt_, active_};
  auto I360 = make_shared<Tensor>(I360_index);
  vector<shared_ptr<Tensor>> tensor337 = {I350, v2_, I360};
  auto task337 = make_shared<Task337>(tensor337, pindex);
  task334->add_dep(task337);
  energyq->add_task(task337);

  vector<shared_ptr<Tensor>> tensor338 = {I360, t2};
  auto task338 = make_shared<Task338>(tensor338, pindex);
  task337->add_dep(task338);
  energyq->add_task(task338);

  vector<IndexRange> I363_index = {active_, active_, virt_, closed_};
  auto I363 = make_shared<Tensor>(I363_index);
  vector<shared_ptr<Tensor>> tensor339 = {I350, v2_, I363};
  auto task339 = make_shared<Task339>(tensor339, pindex);
  task334->add_dep(task339);
  energyq->add_task(task339);

  vector<shared_ptr<Tensor>> tensor340 = {I363, t2};
  auto task340 = make_shared<Task340>(tensor340, pindex);
  task339->add_dep(task340);
  energyq->add_task(task340);

  vector<IndexRange> I369_index = {active_, active_, virt_, closed_};
  auto I369 = make_shared<Tensor>(I369_index);
  vector<shared_ptr<Tensor>> tensor341 = {I350, v2_, I369};
  auto task341 = make_shared<Task341>(tensor341, pindex);
  task334->add_dep(task341);
  energyq->add_task(task341);

  vector<shared_ptr<Tensor>> tensor342 = {I369, t2};
  auto task342 = make_shared<Task342>(tensor342, pindex);
  task341->add_dep(task342);
  energyq->add_task(task342);

  vector<IndexRange> I372_index = {active_, active_, virt_, closed_};
  auto I372 = make_shared<Tensor>(I372_index);
  vector<shared_ptr<Tensor>> tensor343 = {I350, v2_, I372};
  auto task343 = make_shared<Task343>(tensor343, pindex);
  task334->add_dep(task343);
  energyq->add_task(task343);

  vector<shared_ptr<Tensor>> tensor344 = {I372, t2};
  auto task344 = make_shared<Task344>(tensor344, pindex);
  task343->add_dep(task344);
  energyq->add_task(task344);

  vector<IndexRange> I353_index = {active_, active_, active_, active_};
  auto I353 = make_shared<Tensor>(I353_index);
  vector<shared_ptr<Tensor>> tensor345 = {I334, Gamma29_(), I353};
  auto task345 = make_shared<Task345>(tensor345, pindex);
  task320->add_dep(task345);
  energyq->add_task(task345);

  vector<IndexRange> I354_index = {active_, closed_, virt_, active_};
  auto I354 = make_shared<Tensor>(I354_index);
  vector<shared_ptr<Tensor>> tensor346 = {I353, v2_, I354};
  auto task346 = make_shared<Task346>(tensor346, pindex);
  task345->add_dep(task346);
  energyq->add_task(task346);

  vector<shared_ptr<Tensor>> tensor347 = {I354, t2};
  auto task347 = make_shared<Task347>(tensor347, pindex);
  task346->add_dep(task347);
  energyq->add_task(task347);

  vector<IndexRange> I356_index = {active_, active_, active_, active_};
  auto I356 = make_shared<Tensor>(I356_index);
  vector<shared_ptr<Tensor>> tensor348 = {I334, Gamma32_(), I356};
  auto task348 = make_shared<Task348>(tensor348, pindex);
  task320->add_dep(task348);
  energyq->add_task(task348);

  vector<IndexRange> I357_index = {active_, closed_, virt_, active_};
  auto I357 = make_shared<Tensor>(I357_index);
  vector<shared_ptr<Tensor>> tensor349 = {I356, v2_, I357};
  auto task349 = make_shared<Task349>(tensor349, pindex);
  task348->add_dep(task349);
  energyq->add_task(task349);

  vector<shared_ptr<Tensor>> tensor350 = {I357, t2};
  auto task350 = make_shared<Task350>(tensor350, pindex);
  task349->add_dep(task350);
  energyq->add_task(task350);

  vector<IndexRange> I365_index = {active_, active_, active_, active_};
  auto I365 = make_shared<Tensor>(I365_index);
  vector<shared_ptr<Tensor>> tensor351 = {I334, Gamma7_(), I365};
  auto task351 = make_shared<Task351>(tensor351, pindex);
  task320->add_dep(task351);
  energyq->add_task(task351);

  vector<IndexRange> I366_index = {active_, active_, virt_, closed_};
  auto I366 = make_shared<Tensor>(I366_index);
  vector<shared_ptr<Tensor>> tensor352 = {I365, v2_, I366};
  auto task352 = make_shared<Task352>(tensor352, pindex);
  task351->add_dep(task352);
  energyq->add_task(task352);

  vector<shared_ptr<Tensor>> tensor353 = {I366, t2};
  auto task353 = make_shared<Task353>(tensor353, pindex);
  task352->add_dep(task353);
  energyq->add_task(task353);

  vector<IndexRange> I374_index = {active_, active_, active_, active_, active_, active_};
  auto I374 = make_shared<Tensor>(I374_index);
  vector<shared_ptr<Tensor>> tensor354 = {I334, Gamma59_(), I374};
  auto task354 = make_shared<Task354>(tensor354, pindex);
  task320->add_dep(task354);
  energyq->add_task(task354);

  vector<IndexRange> I375_index = {active_, virt_, active_, active_};
  auto I375 = make_shared<Tensor>(I375_index);
  vector<shared_ptr<Tensor>> tensor355 = {I374, t2, I375};
  auto task355 = make_shared<Task355>(tensor355, pindex);
  task354->add_dep(task355);
  energyq->add_task(task355);

  vector<shared_ptr<Tensor>> tensor356 = {I375, v2_};
  auto task356 = make_shared<Task356>(tensor356, pindex);
  task355->add_dep(task356);
  energyq->add_task(task356);

  vector<IndexRange> I377_index = {active_, active_, active_, active_, active_, active_};
  auto I377 = make_shared<Tensor>(I377_index);
  vector<shared_ptr<Tensor>> tensor357 = {I334, Gamma57_(), I377};
  auto task357 = make_shared<Task357>(tensor357, pindex);
  task320->add_dep(task357);
  energyq->add_task(task357);

  vector<IndexRange> I378_index = {active_, active_, virt_, active_};
  auto I378 = make_shared<Tensor>(I378_index);
  vector<shared_ptr<Tensor>> tensor358 = {I377, v2_, I378};
  auto task358 = make_shared<Task358>(tensor358, pindex);
  task357->add_dep(task358);
  energyq->add_task(task358);

  vector<shared_ptr<Tensor>> tensor359 = {I378, t2};
  auto task359 = make_shared<Task359>(tensor359, pindex);
  task358->add_dep(task359);
  energyq->add_task(task359);

  vector<IndexRange> I380_index = {virt_, closed_, virt_, closed_};
  auto I380 = make_shared<Tensor>(I380_index);
  vector<shared_ptr<Tensor>> tensor360 = {I334, v2_, I380};
  auto task360 = make_shared<Task360>(tensor360, pindex);
  task320->add_dep(task360);
  energyq->add_task(task360);

  vector<shared_ptr<Tensor>> tensor361 = {I380, t2};
  auto task361 = make_shared<Task361>(tensor361, pindex);
  task360->add_dep(task361);
  energyq->add_task(task361);

  vector<IndexRange> I382_index = {virt_, closed_, virt_, closed_};
  auto I382 = make_shared<Tensor>(I382_index);
  vector<shared_ptr<Tensor>> tensor362 = {I334, v2_, I382};
  auto task362 = make_shared<Task362>(tensor362, pindex);
  task320->add_dep(task362);
  energyq->add_task(task362);

  vector<shared_ptr<Tensor>> tensor363 = {I382, t2};
  auto task363 = make_shared<Task363>(tensor363, pindex);
  task362->add_dep(task363);
  energyq->add_task(task363);

  vector<IndexRange> I384_index = {active_, active_};
  auto I384 = make_shared<Tensor>(I384_index);
  vector<shared_ptr<Tensor>> tensor364 = {I334, Gamma38_(), I384};
  auto task364 = make_shared<Task364>(tensor364, pindex);
  task320->add_dep(task364);
  energyq->add_task(task364);

  vector<IndexRange> I385_index = {virt_, closed_, virt_, active_};
  auto I385 = make_shared<Tensor>(I385_index);
  vector<shared_ptr<Tensor>> tensor365 = {I384, v2_, I385};
  auto task365 = make_shared<Task365>(tensor365, pindex);
  task364->add_dep(task365);
  energyq->add_task(task365);

  vector<shared_ptr<Tensor>> tensor366 = {I385, t2};
  auto task366 = make_shared<Task366>(tensor366, pindex);
  task365->add_dep(task366);
  energyq->add_task(task366);

  vector<IndexRange> I388_index = {active_, virt_, closed_, virt_};
  auto I388 = make_shared<Tensor>(I388_index);
  vector<shared_ptr<Tensor>> tensor367 = {I384, t2, I388};
  auto task367 = make_shared<Task367>(tensor367, pindex);
  task364->add_dep(task367);
  energyq->add_task(task367);

  vector<shared_ptr<Tensor>> tensor368 = {I388, v2_};
  auto task368 = make_shared<Task368>(tensor368, pindex);
  task367->add_dep(task368);
  energyq->add_task(task368);

  vector<IndexRange> I397_index = {active_, closed_, virt_, active_};
  auto I397 = make_shared<Tensor>(I397_index);
  vector<shared_ptr<Tensor>> tensor369 = {I384, h1_, I397};
  auto task369 = make_shared<Task369>(tensor369, pindex);
  task364->add_dep(task369);
  energyq->add_task(task369);

  vector<shared_ptr<Tensor>> tensor370 = {I397, t2};
  auto task370 = make_shared<Task370>(tensor370, pindex);
  task369->add_dep(task370);
  energyq->add_task(task370);

  vector<IndexRange> I400_index = {active_, active_, virt_, closed_};
  auto I400 = make_shared<Tensor>(I400_index);
  vector<shared_ptr<Tensor>> tensor371 = {I384, h1_, I400};
  auto task371 = make_shared<Task371>(tensor371, pindex);
  task364->add_dep(task371);
  energyq->add_task(task371);

  vector<shared_ptr<Tensor>> tensor372 = {I400, t2};
  auto task372 = make_shared<Task372>(tensor372, pindex);
  task371->add_dep(task372);
  energyq->add_task(task372);

  vector<IndexRange> I390_index = {active_, active_, active_, active_};
  auto I390 = make_shared<Tensor>(I390_index);
  vector<shared_ptr<Tensor>> tensor373 = {I334, Gamma60_(), I390};
  auto task373 = make_shared<Task373>(tensor373, pindex);
  task320->add_dep(task373);
  energyq->add_task(task373);

  vector<IndexRange> I391_index = {active_, virt_, active_, virt_};
  auto I391 = make_shared<Tensor>(I391_index);
  vector<shared_ptr<Tensor>> tensor374 = {I390, t2, I391};
  auto task374 = make_shared<Task374>(tensor374, pindex);
  task373->add_dep(task374);
  energyq->add_task(task374);

  vector<shared_ptr<Tensor>> tensor375 = {I391, v2_};
  auto task375 = make_shared<Task375>(tensor375, pindex);
  task374->add_dep(task375);
  energyq->add_task(task375);

  vector<IndexRange> I403_index = {active_, active_, virt_, active_};
  auto I403 = make_shared<Tensor>(I403_index);
  vector<shared_ptr<Tensor>> tensor376 = {I390, h1_, I403};
  auto task376 = make_shared<Task376>(tensor376, pindex);
  task373->add_dep(task376);
  energyq->add_task(task376);

  vector<shared_ptr<Tensor>> tensor377 = {I403, t2};
  auto task377 = make_shared<Task377>(tensor377, pindex);
  task376->add_dep(task377);
  energyq->add_task(task377);

  vector<IndexRange> I393_index = {closed_, active_};
  auto I393 = make_shared<Tensor>(I393_index);
  vector<shared_ptr<Tensor>> tensor378 = {I334, h1_, I393};
  auto task378 = make_shared<Task378>(tensor378, pindex);
  task320->add_dep(task378);
  energyq->add_task(task378);

  vector<IndexRange> I394_index = {active_, closed_, active_, active_};
  auto I394 = make_shared<Tensor>(I394_index);
  vector<shared_ptr<Tensor>> tensor379 = {I393, Gamma7_(), I394};
  auto task379 = make_shared<Task379>(tensor379, pindex);
  task378->add_dep(task379);
  energyq->add_task(task379);

  vector<shared_ptr<Tensor>> tensor380 = {I394, t2};
  auto task380 = make_shared<Task380>(tensor380, pindex);
  task379->add_dep(task380);
  energyq->add_task(task380);

  return energyq;
}


#endif
