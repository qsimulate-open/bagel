//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_deciqq.cc
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


#include <src/smith/caspt2/MSCASPT2.h>
#include <src/smith/caspt2/MSCASPT2_tasks7.h>
#include <src/smith/caspt2/MSCASPT2_tasks8.h>
#include <src/smith/caspt2/MSCASPT2_tasks9.h>
#include <src/smith/caspt2/MSCASPT2_tasks10.h>
#include <src/smith/caspt2/MSCASPT2_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_deciq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deciq = make_shared<Queue>();
  auto tensor314 = vector<shared_ptr<Tensor>>{deci};
  auto task314 = make_shared<Task314>(tensor314, reset);
  deciq->add_task(task314);

  vector<IndexRange> I324_index = {ci_};
  auto I324 = make_shared<Tensor>(I324_index);
  auto tensor315 = vector<shared_ptr<Tensor>>{deci, I324};
  auto task315 = make_shared<Task315>(tensor315, cindex);
  task315->add_dep(task314);
  deciq->add_task(task315);

  vector<IndexRange> I325_index = {active_, active_, active_, active_};
  auto I325 = make_shared<Tensor>(I325_index);
  auto tensor316 = vector<shared_ptr<Tensor>>{I324, Gamma110_(), I325};
  auto task316 = make_shared<Task316>(tensor316, cindex);
  task315->add_dep(task316);
  task316->add_dep(task314);
  deciq->add_task(task316);

  auto tensor317 = vector<shared_ptr<Tensor>>{I325, t2, l2};
  auto task317 = make_shared<Task317>(tensor317, cindex);
  task316->add_dep(task317);
  task317->add_dep(task314);
  deciq->add_task(task317);

  vector<IndexRange> I328_index = {active_, active_, active_, active_};
  auto I328 = make_shared<Tensor>(I328_index);
  auto tensor318 = vector<shared_ptr<Tensor>>{I324, Gamma111_(), I328};
  auto task318 = make_shared<Task318>(tensor318, cindex);
  task315->add_dep(task318);
  task318->add_dep(task314);
  deciq->add_task(task318);

  vector<IndexRange> I329_index = {closed_, active_, active_, closed_};
  auto I329 = make_shared<Tensor>(I329_index);
  auto tensor319 = vector<shared_ptr<Tensor>>{I328, l2, I329};
  auto task319 = make_shared<Task319>(tensor319, cindex);
  task318->add_dep(task319);
  task319->add_dep(task314);
  deciq->add_task(task319);

  auto tensor320 = vector<shared_ptr<Tensor>>{I329, f1_, t2};
  auto task320 = make_shared<Task320>(tensor320, cindex);
  task319->add_dep(task320);
  task320->add_dep(task314);
  deciq->add_task(task320);

  auto tensor321 = vector<shared_ptr<Tensor>>{I328, t2, l2};
  auto task321 = make_shared<Task321>(tensor321, cindex, this->e0_);
  task318->add_dep(task321);
  task321->add_dep(task314);
  deciq->add_task(task321);

  vector<IndexRange> I332_index = {active_, active_, active_, active_, active_, active_};
  auto I332 = make_shared<Tensor>(I332_index);
  auto tensor322 = vector<shared_ptr<Tensor>>{I324, Gamma112_(), I332};
  auto task322 = make_shared<Task322>(tensor322, cindex);
  task315->add_dep(task322);
  task322->add_dep(task314);
  deciq->add_task(task322);

  vector<IndexRange> I333_index = {active_, active_, closed_, active_};
  auto I333 = make_shared<Tensor>(I333_index);
  auto tensor323 = vector<shared_ptr<Tensor>>{I332, t2, I333};
  auto task323 = make_shared<Task323>(tensor323, cindex);
  task322->add_dep(task323);
  task323->add_dep(task314);
  deciq->add_task(task323);

  auto tensor324 = vector<shared_ptr<Tensor>>{I333, f1_, l2};
  auto task324 = make_shared<Task324>(tensor324, cindex);
  task323->add_dep(task324);
  task324->add_dep(task314);
  deciq->add_task(task324);

  vector<IndexRange> I336_index = {active_, active_, active_, active_};
  auto I336 = make_shared<Tensor>(I336_index);
  auto tensor325 = vector<shared_ptr<Tensor>>{I324, Gamma113_(), I336};
  auto task325 = make_shared<Task325>(tensor325, cindex);
  task315->add_dep(task325);
  task325->add_dep(task314);
  deciq->add_task(task325);

  vector<IndexRange> I337_index = {closed_, closed_, active_, active_};
  auto I337 = make_shared<Tensor>(I337_index);
  auto tensor326 = vector<shared_ptr<Tensor>>{I336, l2, I337};
  auto task326 = make_shared<Task326>(tensor326, cindex);
  task325->add_dep(task326);
  task326->add_dep(task314);
  deciq->add_task(task326);

  auto tensor327 = vector<shared_ptr<Tensor>>{I337, f1_, t2};
  auto task327 = make_shared<Task327>(tensor327, cindex);
  task326->add_dep(task327);
  task327->add_dep(task314);
  deciq->add_task(task327);

  vector<IndexRange> I368_index = {active_, closed_, closed_, active_};
  auto I368 = make_shared<Tensor>(I368_index);
  auto tensor328 = vector<shared_ptr<Tensor>>{I336, t2, I368};
  auto task328 = make_shared<Task328>(tensor328, cindex);
  task325->add_dep(task328);
  task328->add_dep(task314);
  deciq->add_task(task328);

  auto tensor329 = vector<shared_ptr<Tensor>>{I368, f1_, l2};
  auto task329 = make_shared<Task329>(tensor329, cindex);
  task328->add_dep(task329);
  task329->add_dep(task314);
  deciq->add_task(task329);

  vector<IndexRange> I340_index = {active_, active_, active_, active_, active_, active_};
  auto I340 = make_shared<Tensor>(I340_index);
  auto tensor330 = vector<shared_ptr<Tensor>>{I324, Gamma114_(), I340};
  auto task330 = make_shared<Task330>(tensor330, cindex);
  task315->add_dep(task330);
  task330->add_dep(task314);
  deciq->add_task(task330);

  vector<IndexRange> I341_index = {closed_, active_, active_, active_};
  auto I341 = make_shared<Tensor>(I341_index);
  auto tensor331 = vector<shared_ptr<Tensor>>{I340, l2, I341};
  auto task331 = make_shared<Task331>(tensor331, cindex);
  task330->add_dep(task331);
  task331->add_dep(task314);
  deciq->add_task(task331);

  auto tensor332 = vector<shared_ptr<Tensor>>{I341, f1_, t2};
  auto task332 = make_shared<Task332>(tensor332, cindex);
  task331->add_dep(task332);
  task332->add_dep(task314);
  deciq->add_task(task332);

  vector<IndexRange> I344_index = {active_, active_, active_, active_, active_, active_};
  auto I344 = make_shared<Tensor>(I344_index);
  auto tensor333 = vector<shared_ptr<Tensor>>{I324, Gamma115_(), I344};
  auto task333 = make_shared<Task333>(tensor333, cindex);
  task315->add_dep(task333);
  task333->add_dep(task314);
  deciq->add_task(task333);

  auto tensor334 = vector<shared_ptr<Tensor>>{I344, t2, l2};
  auto task334 = make_shared<Task334>(tensor334, cindex);
  task333->add_dep(task334);
  task334->add_dep(task314);
  deciq->add_task(task334);

  vector<IndexRange> I347_index = {active_, active_, active_, active_, active_, active_};
  auto I347 = make_shared<Tensor>(I347_index);
  auto tensor335 = vector<shared_ptr<Tensor>>{I324, Gamma116_(), I347};
  auto task335 = make_shared<Task335>(tensor335, cindex);
  task315->add_dep(task335);
  task335->add_dep(task314);
  deciq->add_task(task335);

  vector<IndexRange> I348_index = {active_, active_, active_, closed_};
  auto I348 = make_shared<Tensor>(I348_index);
  auto tensor336 = vector<shared_ptr<Tensor>>{I347, l2, I348};
  auto task336 = make_shared<Task336>(tensor336, cindex);
  task335->add_dep(task336);
  task336->add_dep(task314);
  deciq->add_task(task336);

  auto tensor337 = vector<shared_ptr<Tensor>>{I348, f1_, t2};
  auto task337 = make_shared<Task337>(tensor337, cindex);
  task336->add_dep(task337);
  task337->add_dep(task314);
  deciq->add_task(task337);

  auto tensor338 = vector<shared_ptr<Tensor>>{I348, f1_, t2};
  auto task338 = make_shared<Task338>(tensor338, cindex);
  task336->add_dep(task338);
  task338->add_dep(task314);
  deciq->add_task(task338);

  vector<IndexRange> I488_index = {active_, active_, closed_, active_};
  auto I488 = make_shared<Tensor>(I488_index);
  auto tensor339 = vector<shared_ptr<Tensor>>{I347, t2, I488};
  auto task339 = make_shared<Task339>(tensor339, cindex);
  task335->add_dep(task339);
  task339->add_dep(task314);
  deciq->add_task(task339);

  auto tensor340 = vector<shared_ptr<Tensor>>{I488, l2};
  auto task340 = make_shared<Task340>(tensor340, cindex, this->e0_);
  task339->add_dep(task340);
  task340->add_dep(task314);
  deciq->add_task(task340);

  auto tensor341 = vector<shared_ptr<Tensor>>{I488, f1_, l2};
  auto task341 = make_shared<Task341>(tensor341, cindex);
  task339->add_dep(task341);
  task341->add_dep(task314);
  deciq->add_task(task341);

  vector<IndexRange> I351_index = {active_, active_, active_, active_};
  auto I351 = make_shared<Tensor>(I351_index);
  auto tensor342 = vector<shared_ptr<Tensor>>{I324, Gamma117_(), I351};
  auto task342 = make_shared<Task342>(tensor342, cindex);
  task315->add_dep(task342);
  task342->add_dep(task314);
  deciq->add_task(task342);

  vector<IndexRange> I352_index = {closed_, active_};
  auto I352 = make_shared<Tensor>(I352_index);
  auto tensor343 = vector<shared_ptr<Tensor>>{I351, l2, I352};
  auto task343 = make_shared<Task343>(tensor343, cindex);
  task342->add_dep(task343);
  task343->add_dep(task314);
  deciq->add_task(task343);

  vector<IndexRange> I353_index = {closed_, virt_, closed_, active_};
  auto I353 = make_shared<Tensor>(I353_index);
  auto tensor344 = vector<shared_ptr<Tensor>>{I352, f1_, I353};
  auto task344 = make_shared<Task344>(tensor344, cindex);
  task343->add_dep(task344);
  task344->add_dep(task314);
  deciq->add_task(task344);

  auto tensor345 = vector<shared_ptr<Tensor>>{I353, t2};
  auto task345 = make_shared<Task345>(tensor345, cindex);
  task344->add_dep(task345);
  task345->add_dep(task314);
  deciq->add_task(task345);

  vector<IndexRange> I442_index = {closed_, virt_, active_, active_};
  auto I442 = make_shared<Tensor>(I442_index);
  auto tensor346 = vector<shared_ptr<Tensor>>{I351, l2, I442};
  auto task346 = make_shared<Task346>(tensor346, cindex);
  task342->add_dep(task346);
  task346->add_dep(task314);
  deciq->add_task(task346);

  auto tensor347 = vector<shared_ptr<Tensor>>{I442, f1_, t2};
  auto task347 = make_shared<Task347>(tensor347, cindex);
  task346->add_dep(task347);
  task347->add_dep(task314);
  deciq->add_task(task347);

  vector<IndexRange> I492_index = {virt_, closed_, active_, active_};
  auto I492 = make_shared<Tensor>(I492_index);
  auto tensor348 = vector<shared_ptr<Tensor>>{I351, l2, I492};
  auto task348 = make_shared<Task348>(tensor348, cindex);
  task342->add_dep(task348);
  task348->add_dep(task314);
  deciq->add_task(task348);

  vector<IndexRange> I493_index = {closed_, virt_, closed_, active_};
  auto I493 = make_shared<Tensor>(I493_index);
  auto tensor349 = vector<shared_ptr<Tensor>>{I492, f1_, I493};
  auto task349 = make_shared<Task349>(tensor349, cindex);
  task348->add_dep(task349);
  task349->add_dep(task314);
  deciq->add_task(task349);

  auto tensor350 = vector<shared_ptr<Tensor>>{I493, t2};
  auto task350 = make_shared<Task350>(tensor350, cindex);
  task349->add_dep(task350);
  task350->add_dep(task314);
  deciq->add_task(task350);

  vector<IndexRange> I359_index = {active_, active_, active_, active_, active_, active_};
  auto I359 = make_shared<Tensor>(I359_index);
  auto tensor351 = vector<shared_ptr<Tensor>>{I324, Gamma119_(), I359};
  auto task351 = make_shared<Task351>(tensor351, cindex);
  task315->add_dep(task351);
  task351->add_dep(task314);
  deciq->add_task(task351);

  vector<IndexRange> I360_index = {active_, closed_, active_, active_};
  auto I360 = make_shared<Tensor>(I360_index);
  auto tensor352 = vector<shared_ptr<Tensor>>{I359, l2, I360};
  auto task352 = make_shared<Task352>(tensor352, cindex);
  task351->add_dep(task352);
  task352->add_dep(task314);
  deciq->add_task(task352);

  auto tensor353 = vector<shared_ptr<Tensor>>{I360, f1_, t2};
  auto task353 = make_shared<Task353>(tensor353, cindex);
  task352->add_dep(task353);
  task353->add_dep(task314);
  deciq->add_task(task353);

  vector<IndexRange> I371_index = {active_, active_, active_, active_};
  auto I371 = make_shared<Tensor>(I371_index);
  auto tensor354 = vector<shared_ptr<Tensor>>{I324, Gamma122_(), I371};
  auto task354 = make_shared<Task354>(tensor354, cindex);
  task315->add_dep(task354);
  task354->add_dep(task314);
  deciq->add_task(task354);

  vector<IndexRange> I372_index = {active_, closed_};
  auto I372 = make_shared<Tensor>(I372_index);
  auto tensor355 = vector<shared_ptr<Tensor>>{I371, t2, I372};
  auto task355 = make_shared<Task355>(tensor355, cindex);
  task354->add_dep(task355);
  task355->add_dep(task314);
  deciq->add_task(task355);

  auto tensor356 = vector<shared_ptr<Tensor>>{I372, f1_, l2};
  auto task356 = make_shared<Task356>(tensor356, cindex);
  task355->add_dep(task356);
  task356->add_dep(task314);
  deciq->add_task(task356);

  vector<IndexRange> I376_index = {active_, closed_};
  auto I376 = make_shared<Tensor>(I376_index);
  auto tensor357 = vector<shared_ptr<Tensor>>{I371, t2, I376};
  auto task357 = make_shared<Task357>(tensor357, cindex);
  task354->add_dep(task357);
  task357->add_dep(task314);
  deciq->add_task(task357);

  auto tensor358 = vector<shared_ptr<Tensor>>{I376, f1_, l2};
  auto task358 = make_shared<Task358>(tensor358, cindex);
  task357->add_dep(task358);
  task358->add_dep(task314);
  deciq->add_task(task358);

  vector<IndexRange> I414_index = {active_, virt_, closed_, active_};
  auto I414 = make_shared<Tensor>(I414_index);
  auto tensor359 = vector<shared_ptr<Tensor>>{I371, t2, I414};
  auto task359 = make_shared<Task359>(tensor359, cindex);
  task354->add_dep(task359);
  task359->add_dep(task314);
  deciq->add_task(task359);

  auto tensor360 = vector<shared_ptr<Tensor>>{I414, f1_, l2};
  auto task360 = make_shared<Task360>(tensor360, cindex);
  task359->add_dep(task360);
  task360->add_dep(task314);
  deciq->add_task(task360);

  vector<IndexRange> I418_index = {active_, closed_, virt_, active_};
  auto I418 = make_shared<Tensor>(I418_index);
  auto tensor361 = vector<shared_ptr<Tensor>>{I371, t2, I418};
  auto task361 = make_shared<Task361>(tensor361, cindex);
  task354->add_dep(task361);
  task361->add_dep(task314);
  deciq->add_task(task361);

  auto tensor362 = vector<shared_ptr<Tensor>>{I418, f1_, l2};
  auto task362 = make_shared<Task362>(tensor362, cindex);
  task361->add_dep(task362);
  task362->add_dep(task314);
  deciq->add_task(task362);

  vector<IndexRange> I422_index = {active_, virt_, closed_, active_};
  auto I422 = make_shared<Tensor>(I422_index);
  auto tensor363 = vector<shared_ptr<Tensor>>{I371, t2, I422};
  auto task363 = make_shared<Task363>(tensor363, cindex);
  task354->add_dep(task363);
  task363->add_dep(task314);
  deciq->add_task(task363);

  auto tensor364 = vector<shared_ptr<Tensor>>{I422, f1_, l2};
  auto task364 = make_shared<Task364>(tensor364, cindex);
  task363->add_dep(task364);
  task364->add_dep(task314);
  deciq->add_task(task364);

  vector<IndexRange> I379_index = {active_, active_};
  auto I379 = make_shared<Tensor>(I379_index);
  auto tensor365 = vector<shared_ptr<Tensor>>{I324, Gamma124_(), I379};
  auto task365 = make_shared<Task365>(tensor365, cindex);
  task315->add_dep(task365);
  task365->add_dep(task314);
  deciq->add_task(task365);

  auto tensor366 = vector<shared_ptr<Tensor>>{I379, t2, l2};
  auto task366 = make_shared<Task366>(tensor366, cindex);
  task365->add_dep(task366);
  task366->add_dep(task314);
  deciq->add_task(task366);

  auto tensor367 = vector<shared_ptr<Tensor>>{I379, t2, l2};
  auto task367 = make_shared<Task367>(tensor367, cindex);
  task365->add_dep(task367);
  task367->add_dep(task314);
  deciq->add_task(task367);

  vector<IndexRange> I385_index = {active_, active_};
  auto I385 = make_shared<Tensor>(I385_index);
  auto tensor368 = vector<shared_ptr<Tensor>>{I324, Gamma126_(), I385};
  auto task368 = make_shared<Task368>(tensor368, cindex);
  task315->add_dep(task368);
  task368->add_dep(task314);
  deciq->add_task(task368);

  vector<IndexRange> I386_index = {virt_, closed_, active_, closed_};
  auto I386 = make_shared<Tensor>(I386_index);
  auto tensor369 = vector<shared_ptr<Tensor>>{I385, l2, I386};
  auto task369 = make_shared<Task369>(tensor369, cindex);
  task368->add_dep(task369);
  task369->add_dep(task314);
  deciq->add_task(task369);

  vector<IndexRange> I387_index = {closed_, virt_, closed_, active_};
  auto I387 = make_shared<Tensor>(I387_index);
  auto tensor370 = vector<shared_ptr<Tensor>>{I386, f1_, I387};
  auto task370 = make_shared<Task370>(tensor370, cindex);
  task369->add_dep(task370);
  task370->add_dep(task314);
  deciq->add_task(task370);

  auto tensor371 = vector<shared_ptr<Tensor>>{I387, t2};
  auto task371 = make_shared<Task371>(tensor371, cindex);
  task370->add_dep(task371);
  task371->add_dep(task314);
  deciq->add_task(task371);

  vector<IndexRange> I395_index = {closed_, virt_, closed_, active_};
  auto I395 = make_shared<Tensor>(I395_index);
  auto tensor372 = vector<shared_ptr<Tensor>>{I386, f1_, I395};
  auto task372 = make_shared<Task372>(tensor372, cindex);
  task369->add_dep(task372);
  task372->add_dep(task314);
  deciq->add_task(task372);

  auto tensor373 = vector<shared_ptr<Tensor>>{I395, t2};
  auto task373 = make_shared<Task373>(tensor373, cindex);
  task372->add_dep(task373);
  task373->add_dep(task314);
  deciq->add_task(task373);

  vector<IndexRange> I399_index = {closed_, virt_, closed_, active_};
  auto I399 = make_shared<Tensor>(I399_index);
  auto tensor374 = vector<shared_ptr<Tensor>>{I386, f1_, I399};
  auto task374 = make_shared<Task374>(tensor374, cindex);
  task369->add_dep(task374);
  task374->add_dep(task314);
  deciq->add_task(task374);

  auto tensor375 = vector<shared_ptr<Tensor>>{I399, t2};
  auto task375 = make_shared<Task375>(tensor375, cindex);
  task374->add_dep(task375);
  task375->add_dep(task314);
  deciq->add_task(task375);

  vector<IndexRange> I426_index = {active_, virt_};
  auto I426 = make_shared<Tensor>(I426_index);
  auto tensor376 = vector<shared_ptr<Tensor>>{I385, f1_, I426};
  auto task376 = make_shared<Task376>(tensor376, cindex);
  task368->add_dep(task376);
  task376->add_dep(task314);
  deciq->add_task(task376);

  auto tensor377 = vector<shared_ptr<Tensor>>{I426, t2, l2};
  auto task377 = make_shared<Task377>(tensor377, cindex);
  task376->add_dep(task377);
  task377->add_dep(task314);
  deciq->add_task(task377);

  auto tensor378 = vector<shared_ptr<Tensor>>{I426, t2, l2};
  auto task378 = make_shared<Task378>(tensor378, cindex);
  task376->add_dep(task378);
  task378->add_dep(task314);
  deciq->add_task(task378);

  vector<IndexRange> I569_index = {virt_, active_};
  auto I569 = make_shared<Tensor>(I569_index);
  auto tensor379 = vector<shared_ptr<Tensor>>{I385, f1_, I569};
  auto task379 = make_shared<Task379>(tensor379, cindex);
  task368->add_dep(task379);
  task379->add_dep(task314);
  deciq->add_task(task379);

  auto tensor380 = vector<shared_ptr<Tensor>>{I569, t2, l2};
  auto task380 = make_shared<Task380>(tensor380, cindex);
  task379->add_dep(task380);
  task380->add_dep(task314);
  deciq->add_task(task380);

  vector<IndexRange> I573_index = {virt_, active_};
  auto I573 = make_shared<Tensor>(I573_index);
  auto tensor381 = vector<shared_ptr<Tensor>>{I385, f1_, I573};
  auto task381 = make_shared<Task381>(tensor381, cindex);
  task368->add_dep(task381);
  task381->add_dep(task314);
  deciq->add_task(task381);

  auto tensor382 = vector<shared_ptr<Tensor>>{I573, t2, l2};
  auto task382 = make_shared<Task382>(tensor382, cindex);
  task381->add_dep(task382);
  task382->add_dep(task314);
  deciq->add_task(task382);

  auto tensor383 = vector<shared_ptr<Tensor>>{I385, t2, l2};
  auto task383 = make_shared<Task383>(tensor383, cindex, this->e0_);
  task368->add_dep(task383);
  task383->add_dep(task314);
  deciq->add_task(task383);

  auto tensor384 = vector<shared_ptr<Tensor>>{I385, t2, l2};
  auto task384 = make_shared<Task384>(tensor384, cindex, this->e0_);
  task368->add_dep(task384);
  task384->add_dep(task314);
  deciq->add_task(task384);

  vector<IndexRange> I409_index = {active_, active_, active_, active_};
  auto I409 = make_shared<Tensor>(I409_index);
  auto tensor385 = vector<shared_ptr<Tensor>>{I324, Gamma132_(), I409};
  auto task385 = make_shared<Task385>(tensor385, cindex);
  task315->add_dep(task385);
  task385->add_dep(task314);
  deciq->add_task(task385);

  vector<IndexRange> I410_index = {active_, closed_, virt_, active_};
  auto I410 = make_shared<Tensor>(I410_index);
  auto tensor386 = vector<shared_ptr<Tensor>>{I409, t2, I410};
  auto task386 = make_shared<Task386>(tensor386, cindex);
  task385->add_dep(task386);
  task386->add_dep(task314);
  deciq->add_task(task386);

  auto tensor387 = vector<shared_ptr<Tensor>>{I410, f1_, l2};
  auto task387 = make_shared<Task387>(tensor387, cindex);
  task386->add_dep(task387);
  task387->add_dep(task314);
  deciq->add_task(task387);

  vector<IndexRange> I433_index = {active_, active_, active_, active_, active_, active_};
  auto I433 = make_shared<Tensor>(I433_index);
  auto tensor388 = vector<shared_ptr<Tensor>>{I324, Gamma138_(), I433};
  auto task388 = make_shared<Task388>(tensor388, cindex);
  task315->add_dep(task388);
  task388->add_dep(task314);
  deciq->add_task(task388);

  vector<IndexRange> I434_index = {active_, closed_, active_, active_};
  auto I434 = make_shared<Tensor>(I434_index);
  auto tensor389 = vector<shared_ptr<Tensor>>{I433, t2, I434};
  auto task389 = make_shared<Task389>(tensor389, cindex);
  task388->add_dep(task389);
  task389->add_dep(task314);
  deciq->add_task(task389);

  auto tensor390 = vector<shared_ptr<Tensor>>{I434, f1_, l2};
  auto task390 = make_shared<Task390>(tensor390, cindex);
  task389->add_dep(task390);
  task390->add_dep(task314);
  deciq->add_task(task390);

  vector<IndexRange> I437_index = {active_, active_, active_, active_};
  auto I437 = make_shared<Tensor>(I437_index);
  auto tensor391 = vector<shared_ptr<Tensor>>{I324, Gamma139_(), I437};
  auto task391 = make_shared<Task391>(tensor391, cindex);
  task315->add_dep(task391);
  task391->add_dep(task314);
  deciq->add_task(task391);

  vector<IndexRange> I438_index = {virt_, closed_, active_, active_};
  auto I438 = make_shared<Tensor>(I438_index);
  auto tensor392 = vector<shared_ptr<Tensor>>{I437, l2, I438};
  auto task392 = make_shared<Task392>(tensor392, cindex);
  task391->add_dep(task392);
  task392->add_dep(task314);
  deciq->add_task(task392);

  auto tensor393 = vector<shared_ptr<Tensor>>{I438, f1_, t2};
  auto task393 = make_shared<Task393>(tensor393, cindex);
  task392->add_dep(task393);
  task393->add_dep(task314);
  deciq->add_task(task393);

  vector<IndexRange> I445_index = {active_, active_, active_, active_};
  auto I445 = make_shared<Tensor>(I445_index);
  auto tensor394 = vector<shared_ptr<Tensor>>{I324, Gamma141_(), I445};
  auto task394 = make_shared<Task394>(tensor394, cindex);
  task315->add_dep(task394);
  task394->add_dep(task314);
  deciq->add_task(task394);

  auto tensor395 = vector<shared_ptr<Tensor>>{I445, t2, l2};
  auto task395 = make_shared<Task395>(tensor395, cindex);
  task394->add_dep(task395);
  task395->add_dep(task314);
  deciq->add_task(task395);

  vector<IndexRange> I448_index = {active_, active_, active_, active_};
  auto I448 = make_shared<Tensor>(I448_index);
  auto tensor396 = vector<shared_ptr<Tensor>>{I324, Gamma142_(), I448};
  auto task396 = make_shared<Task396>(tensor396, cindex);
  task315->add_dep(task396);
  task396->add_dep(task314);
  deciq->add_task(task396);

  vector<IndexRange> I449_index = {active_, virt_, active_, closed_};
  auto I449 = make_shared<Tensor>(I449_index);
  auto tensor397 = vector<shared_ptr<Tensor>>{I448, l2, I449};
  auto task397 = make_shared<Task397>(tensor397, cindex);
  task396->add_dep(task397);
  task397->add_dep(task314);
  deciq->add_task(task397);

  auto tensor398 = vector<shared_ptr<Tensor>>{I449, f1_, t2};
  auto task398 = make_shared<Task398>(tensor398, cindex);
  task397->add_dep(task398);
  task398->add_dep(task314);
  deciq->add_task(task398);

  auto tensor399 = vector<shared_ptr<Tensor>>{I449, f1_, t2};
  auto task399 = make_shared<Task399>(tensor399, cindex);
  task397->add_dep(task399);
  task399->add_dep(task314);
  deciq->add_task(task399);

  auto tensor400 = vector<shared_ptr<Tensor>>{I449, f1_, t2};
  auto task400 = make_shared<Task400>(tensor400, cindex);
  task397->add_dep(task400);
  task400->add_dep(task314);
  deciq->add_task(task400);

  vector<IndexRange> I627_index = {closed_, virt_, active_, active_};
  auto I627 = make_shared<Tensor>(I627_index);
  auto tensor401 = vector<shared_ptr<Tensor>>{I448, t2, I627};
  auto task401 = make_shared<Task401>(tensor401, cindex);
  task396->add_dep(task401);
  task401->add_dep(task314);
  deciq->add_task(task401);

  auto tensor402 = vector<shared_ptr<Tensor>>{I627, l2};
  auto task402 = make_shared<Task402>(tensor402, cindex, this->e0_);
  task401->add_dep(task402);
  task402->add_dep(task314);
  deciq->add_task(task402);

  auto tensor403 = vector<shared_ptr<Tensor>>{I627, f1_, l2};
  auto task403 = make_shared<Task403>(tensor403, cindex);
  task401->add_dep(task403);
  task403->add_dep(task314);
  deciq->add_task(task403);

  vector<IndexRange> I456_index = {active_, active_, active_, active_};
  auto I456 = make_shared<Tensor>(I456_index);
  auto tensor404 = vector<shared_ptr<Tensor>>{I324, Gamma144_(), I456};
  auto task404 = make_shared<Task404>(tensor404, cindex);
  task315->add_dep(task404);
  task404->add_dep(task314);
  deciq->add_task(task404);

  auto tensor405 = vector<shared_ptr<Tensor>>{I456, t2, l2};
  auto task405 = make_shared<Task405>(tensor405, cindex);
  task404->add_dep(task405);
  task405->add_dep(task314);
  deciq->add_task(task405);

  auto tensor406 = vector<shared_ptr<Tensor>>{I456, t2, l2};
  auto task406 = make_shared<Task406>(tensor406, cindex);
  task404->add_dep(task406);
  task406->add_dep(task314);
  deciq->add_task(task406);

  auto tensor407 = vector<shared_ptr<Tensor>>{I456, t2, l2};
  auto task407 = make_shared<Task407>(tensor407, cindex);
  task404->add_dep(task407);
  task407->add_dep(task314);
  deciq->add_task(task407);

  vector<IndexRange> I459_index = {active_, active_, active_, active_};
  auto I459 = make_shared<Tensor>(I459_index);
  auto tensor408 = vector<shared_ptr<Tensor>>{I324, Gamma145_(), I459};
  auto task408 = make_shared<Task408>(tensor408, cindex);
  task315->add_dep(task408);
  task408->add_dep(task314);
  deciq->add_task(task408);

  vector<IndexRange> I460_index = {virt_, active_, active_, closed_};
  auto I460 = make_shared<Tensor>(I460_index);
  auto tensor409 = vector<shared_ptr<Tensor>>{I459, l2, I460};
  auto task409 = make_shared<Task409>(tensor409, cindex);
  task408->add_dep(task409);
  task409->add_dep(task314);
  deciq->add_task(task409);

  auto tensor410 = vector<shared_ptr<Tensor>>{I460, f1_, t2};
  auto task410 = make_shared<Task410>(tensor410, cindex);
  task409->add_dep(task410);
  task410->add_dep(task314);
  deciq->add_task(task410);

  auto tensor411 = vector<shared_ptr<Tensor>>{I460, f1_, t2};
  auto task411 = make_shared<Task411>(tensor411, cindex);
  task409->add_dep(task411);
  task411->add_dep(task314);
  deciq->add_task(task411);

  auto tensor412 = vector<shared_ptr<Tensor>>{I460, f1_, t2};
  auto task412 = make_shared<Task412>(tensor412, cindex);
  task409->add_dep(task412);
  task412->add_dep(task314);
  deciq->add_task(task412);

  vector<IndexRange> I503_index = {active_, virt_, active_, closed_};
  auto I503 = make_shared<Tensor>(I503_index);
  auto tensor413 = vector<shared_ptr<Tensor>>{I459, l2, I503};
  auto task413 = make_shared<Task413>(tensor413, cindex);
  task408->add_dep(task413);
  task413->add_dep(task314);
  deciq->add_task(task413);

  vector<IndexRange> I504_index = {active_, virt_, closed_, active_};
  auto I504 = make_shared<Tensor>(I504_index);
  auto tensor414 = vector<shared_ptr<Tensor>>{I503, f1_, I504};
  auto task414 = make_shared<Task414>(tensor414, cindex);
  task413->add_dep(task414);
  task414->add_dep(task314);
  deciq->add_task(task414);

  auto tensor415 = vector<shared_ptr<Tensor>>{I504, t2};
  auto task415 = make_shared<Task415>(tensor415, cindex);
  task414->add_dep(task415);
  task415->add_dep(task314);
  deciq->add_task(task415);

  vector<IndexRange> I508_index = {active_, virt_, closed_, active_};
  auto I508 = make_shared<Tensor>(I508_index);
  auto tensor416 = vector<shared_ptr<Tensor>>{I503, f1_, I508};
  auto task416 = make_shared<Task416>(tensor416, cindex);
  task413->add_dep(task416);
  task416->add_dep(task314);
  deciq->add_task(task416);

  auto tensor417 = vector<shared_ptr<Tensor>>{I508, t2};
  auto task417 = make_shared<Task417>(tensor417, cindex);
  task416->add_dep(task417);
  task417->add_dep(task314);
  deciq->add_task(task417);

  vector<IndexRange> I535_index = {active_, virt_, closed_, virt_};
  auto I535 = make_shared<Tensor>(I535_index);
  auto tensor418 = vector<shared_ptr<Tensor>>{I503, f1_, I535};
  auto task418 = make_shared<Task418>(tensor418, cindex);
  task413->add_dep(task418);
  task418->add_dep(task314);
  deciq->add_task(task418);

  auto tensor419 = vector<shared_ptr<Tensor>>{I535, t2};
  auto task419 = make_shared<Task419>(tensor419, cindex);
  task418->add_dep(task419);
  task419->add_dep(task314);
  deciq->add_task(task419);

  vector<IndexRange> I623_index = {virt_, closed_, active_, active_};
  auto I623 = make_shared<Tensor>(I623_index);
  auto tensor420 = vector<shared_ptr<Tensor>>{I459, t2, I623};
  auto task420 = make_shared<Task420>(tensor420, cindex);
  task408->add_dep(task420);
  task420->add_dep(task314);
  deciq->add_task(task420);

  auto tensor421 = vector<shared_ptr<Tensor>>{I623, f1_, l2};
  auto task421 = make_shared<Task421>(tensor421, cindex);
  task420->add_dep(task421);
  task421->add_dep(task314);
  deciq->add_task(task421);

  vector<IndexRange> I631_index = {virt_, closed_, active_, active_};
  auto I631 = make_shared<Tensor>(I631_index);
  auto tensor422 = vector<shared_ptr<Tensor>>{I459, t2, I631};
  auto task422 = make_shared<Task422>(tensor422, cindex);
  task408->add_dep(task422);
  task422->add_dep(task314);
  deciq->add_task(task422);

  auto tensor423 = vector<shared_ptr<Tensor>>{I631, f1_, l2};
  auto task423 = make_shared<Task423>(tensor423, cindex);
  task422->add_dep(task423);
  task423->add_dep(task314);
  deciq->add_task(task423);

  vector<IndexRange> I635_index = {closed_, virt_, active_, active_};
  auto I635 = make_shared<Tensor>(I635_index);
  auto tensor424 = vector<shared_ptr<Tensor>>{I459, t2, I635};
  auto task424 = make_shared<Task424>(tensor424, cindex);
  task408->add_dep(task424);
  task424->add_dep(task314);
  deciq->add_task(task424);

  auto tensor425 = vector<shared_ptr<Tensor>>{I635, l2};
  auto task425 = make_shared<Task425>(tensor425, cindex, this->e0_);
  task424->add_dep(task425);
  task425->add_dep(task314);
  deciq->add_task(task425);

  auto tensor426 = vector<shared_ptr<Tensor>>{I635, f1_, l2};
  auto task426 = make_shared<Task426>(tensor426, cindex);
  task424->add_dep(task426);
  task426->add_dep(task314);
  deciq->add_task(task426);

  auto tensor427 = vector<shared_ptr<Tensor>>{I459, t2, l2};
  auto task427 = make_shared<Task427>(tensor427, cindex, this->e0_);
  task408->add_dep(task427);
  task427->add_dep(task314);
  deciq->add_task(task427);

  auto tensor428 = vector<shared_ptr<Tensor>>{I459, t2, l2};
  auto task428 = make_shared<Task428>(tensor428, cindex, this->e0_);
  task408->add_dep(task428);
  task428->add_dep(task314);
  deciq->add_task(task428);

  vector<IndexRange> I467_index = {active_, active_, active_, active_, active_, active_};
  auto I467 = make_shared<Tensor>(I467_index);
  auto tensor429 = vector<shared_ptr<Tensor>>{I324, Gamma147_(), I467};
  auto task429 = make_shared<Task429>(tensor429, cindex);
  task315->add_dep(task429);
  task429->add_dep(task314);
  deciq->add_task(task429);

  vector<IndexRange> I468_index = {active_, virt_, active_, active_};
  auto I468 = make_shared<Tensor>(I468_index);
  auto tensor430 = vector<shared_ptr<Tensor>>{I467, t2, I468};
  auto task430 = make_shared<Task430>(tensor430, cindex);
  task429->add_dep(task430);
  task430->add_dep(task314);
  deciq->add_task(task430);

  auto tensor431 = vector<shared_ptr<Tensor>>{I468, f1_, l2};
  auto task431 = make_shared<Task431>(tensor431, cindex);
  task430->add_dep(task431);
  task431->add_dep(task314);
  deciq->add_task(task431);

  vector<IndexRange> I471_index = {active_, active_};
  auto I471 = make_shared<Tensor>(I471_index);
  auto tensor432 = vector<shared_ptr<Tensor>>{I324, Gamma148_(), I471};
  auto task432 = make_shared<Task432>(tensor432, cindex);
  task315->add_dep(task432);
  task432->add_dep(task314);
  deciq->add_task(task432);

  vector<IndexRange> I472_index = {closed_, virt_};
  auto I472 = make_shared<Tensor>(I472_index);
  auto tensor433 = vector<shared_ptr<Tensor>>{I471, l2, I472};
  auto task433 = make_shared<Task433>(tensor433, cindex);
  task432->add_dep(task433);
  task433->add_dep(task314);
  deciq->add_task(task433);

  vector<IndexRange> I473_index = {closed_, virt_, closed_, virt_};
  auto I473 = make_shared<Tensor>(I473_index);
  auto tensor434 = vector<shared_ptr<Tensor>>{I472, f1_, I473};
  auto task434 = make_shared<Task434>(tensor434, cindex);
  task433->add_dep(task434);
  task434->add_dep(task314);
  deciq->add_task(task434);

  auto tensor435 = vector<shared_ptr<Tensor>>{I473, t2};
  auto task435 = make_shared<Task435>(tensor435, cindex);
  task434->add_dep(task435);
  task435->add_dep(task314);
  deciq->add_task(task435);

  vector<IndexRange> I526_index = {closed_, virt_};
  auto I526 = make_shared<Tensor>(I526_index);
  auto tensor436 = vector<shared_ptr<Tensor>>{I471, l2, I526};
  auto task436 = make_shared<Task436>(tensor436, cindex);
  task432->add_dep(task436);
  task436->add_dep(task314);
  deciq->add_task(task436);

  vector<IndexRange> I527_index = {closed_, virt_, closed_, virt_};
  auto I527 = make_shared<Tensor>(I527_index);
  auto tensor437 = vector<shared_ptr<Tensor>>{I526, f1_, I527};
  auto task437 = make_shared<Task437>(tensor437, cindex);
  task436->add_dep(task437);
  task437->add_dep(task314);
  deciq->add_task(task437);

  auto tensor438 = vector<shared_ptr<Tensor>>{I527, t2};
  auto task438 = make_shared<Task438>(tensor438, cindex);
  task437->add_dep(task438);
  task438->add_dep(task314);
  deciq->add_task(task438);

  vector<IndexRange> I577_index = {virt_, closed_};
  auto I577 = make_shared<Tensor>(I577_index);
  auto tensor439 = vector<shared_ptr<Tensor>>{I471, t2, I577};
  auto task439 = make_shared<Task439>(tensor439, cindex);
  task432->add_dep(task439);
  task439->add_dep(task314);
  deciq->add_task(task439);

  auto tensor440 = vector<shared_ptr<Tensor>>{I577, f1_, l2};
  auto task440 = make_shared<Task440>(tensor440, cindex);
  task439->add_dep(task440);
  task440->add_dep(task314);
  deciq->add_task(task440);

  vector<IndexRange> I581_index = {virt_, closed_};
  auto I581 = make_shared<Tensor>(I581_index);
  auto tensor441 = vector<shared_ptr<Tensor>>{I471, t2, I581};
  auto task441 = make_shared<Task441>(tensor441, cindex);
  task432->add_dep(task441);
  task441->add_dep(task314);
  deciq->add_task(task441);

  auto tensor442 = vector<shared_ptr<Tensor>>{I581, f1_, l2};
  auto task442 = make_shared<Task442>(tensor442, cindex);
  task441->add_dep(task442);
  task442->add_dep(task314);
  deciq->add_task(task442);

  vector<IndexRange> I585_index = {virt_, closed_};
  auto I585 = make_shared<Tensor>(I585_index);
  auto tensor443 = vector<shared_ptr<Tensor>>{I471, t2, I585};
  auto task443 = make_shared<Task443>(tensor443, cindex);
  task432->add_dep(task443);
  task443->add_dep(task314);
  deciq->add_task(task443);

  auto tensor444 = vector<shared_ptr<Tensor>>{I585, f1_, l2};
  auto task444 = make_shared<Task444>(tensor444, cindex);
  task443->add_dep(task444);
  task444->add_dep(task314);
  deciq->add_task(task444);

  vector<IndexRange> I589_index = {virt_, closed_};
  auto I589 = make_shared<Tensor>(I589_index);
  auto tensor445 = vector<shared_ptr<Tensor>>{I471, t2, I589};
  auto task445 = make_shared<Task445>(tensor445, cindex);
  task432->add_dep(task445);
  task445->add_dep(task314);
  deciq->add_task(task445);

  auto tensor446 = vector<shared_ptr<Tensor>>{I589, f1_, l2};
  auto task446 = make_shared<Task446>(tensor446, cindex);
  task445->add_dep(task446);
  task446->add_dep(task314);
  deciq->add_task(task446);

  vector<IndexRange> I615_index = {closed_, active_};
  auto I615 = make_shared<Tensor>(I615_index);
  auto tensor447 = vector<shared_ptr<Tensor>>{I471, f1_, I615};
  auto task447 = make_shared<Task447>(tensor447, cindex);
  task432->add_dep(task447);
  task447->add_dep(task314);
  deciq->add_task(task447);

  auto tensor448 = vector<shared_ptr<Tensor>>{I615, t2, l2};
  auto task448 = make_shared<Task448>(tensor448, cindex);
  task447->add_dep(task448);
  task448->add_dep(task314);
  deciq->add_task(task448);

  auto tensor449 = vector<shared_ptr<Tensor>>{I615, t2, l2};
  auto task449 = make_shared<Task449>(tensor449, cindex);
  task447->add_dep(task449);
  task449->add_dep(task314);
  deciq->add_task(task449);

  vector<IndexRange> I647_index = {active_, closed_};
  auto I647 = make_shared<Tensor>(I647_index);
  auto tensor450 = vector<shared_ptr<Tensor>>{I471, f1_, I647};
  auto task450 = make_shared<Task450>(tensor450, cindex);
  task432->add_dep(task450);
  task450->add_dep(task314);
  deciq->add_task(task450);

  auto tensor451 = vector<shared_ptr<Tensor>>{I647, t2, l2};
  auto task451 = make_shared<Task451>(tensor451, cindex);
  task450->add_dep(task451);
  task451->add_dep(task314);
  deciq->add_task(task451);

  auto tensor452 = vector<shared_ptr<Tensor>>{I647, t2, l2};
  auto task452 = make_shared<Task452>(tensor452, cindex);
  task450->add_dep(task452);
  task452->add_dep(task314);
  deciq->add_task(task452);

  vector<IndexRange> I661_index = {active_, virt_, virt_, closed_};
  auto I661 = make_shared<Tensor>(I661_index);
  auto tensor453 = vector<shared_ptr<Tensor>>{I471, l2, I661};
  auto task453 = make_shared<Task453>(tensor453, cindex);
  task432->add_dep(task453);
  task453->add_dep(task314);
  deciq->add_task(task453);

  vector<IndexRange> I662_index = {active_, virt_, closed_, virt_};
  auto I662 = make_shared<Tensor>(I662_index);
  auto tensor454 = vector<shared_ptr<Tensor>>{I661, f1_, I662};
  auto task454 = make_shared<Task454>(tensor454, cindex);
  task453->add_dep(task454);
  task454->add_dep(task314);
  deciq->add_task(task454);

  auto tensor455 = vector<shared_ptr<Tensor>>{I662, t2};
  auto task455 = make_shared<Task455>(tensor455, cindex);
  task454->add_dep(task455);
  task455->add_dep(task314);
  deciq->add_task(task455);

  vector<IndexRange> I670_index = {active_, virt_, closed_, virt_};
  auto I670 = make_shared<Tensor>(I670_index);
  auto tensor456 = vector<shared_ptr<Tensor>>{I661, f1_, I670};
  auto task456 = make_shared<Task456>(tensor456, cindex);
  task453->add_dep(task456);
  task456->add_dep(task314);
  deciq->add_task(task456);

  auto tensor457 = vector<shared_ptr<Tensor>>{I670, t2};
  auto task457 = make_shared<Task457>(tensor457, cindex);
  task456->add_dep(task457);
  task457->add_dep(task314);
  deciq->add_task(task457);

  vector<IndexRange> I678_index = {active_, virt_, closed_, virt_};
  auto I678 = make_shared<Tensor>(I678_index);
  auto tensor458 = vector<shared_ptr<Tensor>>{I661, f1_, I678};
  auto task458 = make_shared<Task458>(tensor458, cindex);
  task453->add_dep(task458);
  task458->add_dep(task314);
  deciq->add_task(task458);

  auto tensor459 = vector<shared_ptr<Tensor>>{I678, t2};
  auto task459 = make_shared<Task459>(tensor459, cindex);
  task458->add_dep(task459);
  task459->add_dep(task314);
  deciq->add_task(task459);

  auto tensor460 = vector<shared_ptr<Tensor>>{I471, t2, l2};
  auto task460 = make_shared<Task460>(tensor460, cindex, this->e0_);
  task432->add_dep(task460);
  task460->add_dep(task314);
  deciq->add_task(task460);

  auto tensor461 = vector<shared_ptr<Tensor>>{I471, t2, l2};
  auto task461 = make_shared<Task461>(tensor461, cindex, this->e0_);
  task432->add_dep(task461);
  task461->add_dep(task314);
  deciq->add_task(task461);

  vector<IndexRange> I521_index = {active_, active_, active_, active_, active_, active_};
  auto I521 = make_shared<Tensor>(I521_index);
  auto tensor462 = vector<shared_ptr<Tensor>>{I324, Gamma161_(), I521};
  auto task462 = make_shared<Task462>(tensor462, cindex);
  task315->add_dep(task462);
  task462->add_dep(task314);
  deciq->add_task(task462);

  vector<IndexRange> I522_index = {active_, active_, virt_, active_};
  auto I522 = make_shared<Tensor>(I522_index);
  auto tensor463 = vector<shared_ptr<Tensor>>{I521, t2, I522};
  auto task463 = make_shared<Task463>(tensor463, cindex);
  task462->add_dep(task463);
  task463->add_dep(task314);
  deciq->add_task(task463);

  auto tensor464 = vector<shared_ptr<Tensor>>{I522, f1_, l2};
  auto task464 = make_shared<Task464>(tensor464, cindex);
  task463->add_dep(task464);
  task464->add_dep(task314);
  deciq->add_task(task464);

  vector<IndexRange> I541_index = {active_, active_, active_, active_, active_, active_};
  auto I541 = make_shared<Tensor>(I541_index);
  auto tensor465 = vector<shared_ptr<Tensor>>{I324, Gamma166_(), I541};
  auto task465 = make_shared<Task465>(tensor465, cindex);
  task315->add_dep(task465);
  task465->add_dep(task314);
  deciq->add_task(task465);

  vector<IndexRange> I542_index = {active_, virt_, active_, active_};
  auto I542 = make_shared<Tensor>(I542_index);
  auto tensor466 = vector<shared_ptr<Tensor>>{I541, l2, I542};
  auto task466 = make_shared<Task466>(tensor466, cindex);
  task465->add_dep(task466);
  task466->add_dep(task314);
  deciq->add_task(task466);

  auto tensor467 = vector<shared_ptr<Tensor>>{I542, f1_, t2};
  auto task467 = make_shared<Task467>(tensor467, cindex);
  task466->add_dep(task467);
  task467->add_dep(task314);
  deciq->add_task(task467);

  vector<IndexRange> I545_index = {active_, active_, active_, active_, active_, active_};
  auto I545 = make_shared<Tensor>(I545_index);
  auto tensor468 = vector<shared_ptr<Tensor>>{I324, Gamma167_(), I545};
  auto task468 = make_shared<Task468>(tensor468, cindex);
  task315->add_dep(task468);
  task468->add_dep(task314);
  deciq->add_task(task468);

  vector<IndexRange> I546_index = {virt_, active_, active_, active_};
  auto I546 = make_shared<Tensor>(I546_index);
  auto tensor469 = vector<shared_ptr<Tensor>>{I545, l2, I546};
  auto task469 = make_shared<Task469>(tensor469, cindex);
  task468->add_dep(task469);
  task469->add_dep(task314);
  deciq->add_task(task469);

  auto tensor470 = vector<shared_ptr<Tensor>>{I546, f1_, t2};
  auto task470 = make_shared<Task470>(tensor470, cindex);
  task469->add_dep(task470);
  task470->add_dep(task314);
  deciq->add_task(task470);

  vector<IndexRange> I549_index = {active_, active_, active_, active_, active_, active_};
  auto I549 = make_shared<Tensor>(I549_index);
  auto tensor471 = vector<shared_ptr<Tensor>>{I324, Gamma168_(), I549};
  auto task471 = make_shared<Task471>(tensor471, cindex);
  task315->add_dep(task471);
  task471->add_dep(task314);
  deciq->add_task(task471);

  auto tensor472 = vector<shared_ptr<Tensor>>{I549, t2, l2};
  auto task472 = make_shared<Task472>(tensor472, cindex);
  task471->add_dep(task472);
  task472->add_dep(task314);
  deciq->add_task(task472);

  vector<IndexRange> I552_index = {active_, active_, active_, active_, active_, active_};
  auto I552 = make_shared<Tensor>(I552_index);
  auto tensor473 = vector<shared_ptr<Tensor>>{I324, Gamma169_(), I552};
  auto task473 = make_shared<Task473>(tensor473, cindex);
  task315->add_dep(task473);
  task473->add_dep(task314);
  deciq->add_task(task473);

  vector<IndexRange> I553_index = {active_, active_, active_, virt_};
  auto I553 = make_shared<Tensor>(I553_index);
  auto tensor474 = vector<shared_ptr<Tensor>>{I552, l2, I553};
  auto task474 = make_shared<Task474>(tensor474, cindex);
  task473->add_dep(task474);
  task474->add_dep(task314);
  deciq->add_task(task474);

  auto tensor475 = vector<shared_ptr<Tensor>>{I553, f1_, t2};
  auto task475 = make_shared<Task475>(tensor475, cindex);
  task474->add_dep(task475);
  task475->add_dep(task314);
  deciq->add_task(task475);

  auto tensor476 = vector<shared_ptr<Tensor>>{I553, f1_, t2};
  auto task476 = make_shared<Task476>(tensor476, cindex);
  task474->add_dep(task476);
  task476->add_dep(task314);
  deciq->add_task(task476);

  vector<IndexRange> I689_index = {active_, virt_, active_, active_};
  auto I689 = make_shared<Tensor>(I689_index);
  auto tensor477 = vector<shared_ptr<Tensor>>{I552, t2, I689};
  auto task477 = make_shared<Task477>(tensor477, cindex);
  task473->add_dep(task477);
  task477->add_dep(task314);
  deciq->add_task(task477);

  auto tensor478 = vector<shared_ptr<Tensor>>{I689, l2};
  auto task478 = make_shared<Task478>(tensor478, cindex, this->e0_);
  task477->add_dep(task478);
  task478->add_dep(task314);
  deciq->add_task(task478);

  auto tensor479 = vector<shared_ptr<Tensor>>{I689, f1_, l2};
  auto task479 = make_shared<Task479>(tensor479, cindex);
  task477->add_dep(task479);
  task479->add_dep(task314);
  deciq->add_task(task479);

  vector<IndexRange> I556_index = {active_, active_, active_, active_};
  auto I556 = make_shared<Tensor>(I556_index);
  auto tensor480 = vector<shared_ptr<Tensor>>{I324, Gamma170_(), I556};
  auto task480 = make_shared<Task480>(tensor480, cindex);
  task315->add_dep(task480);
  task480->add_dep(task314);
  deciq->add_task(task480);

  vector<IndexRange> I557_index = {active_, virt_};
  auto I557 = make_shared<Tensor>(I557_index);
  auto tensor481 = vector<shared_ptr<Tensor>>{I556, l2, I557};
  auto task481 = make_shared<Task481>(tensor481, cindex);
  task480->add_dep(task481);
  task481->add_dep(task314);
  deciq->add_task(task481);

  vector<IndexRange> I558_index = {active_, virt_, closed_, virt_};
  auto I558 = make_shared<Tensor>(I558_index);
  auto tensor482 = vector<shared_ptr<Tensor>>{I557, f1_, I558};
  auto task482 = make_shared<Task482>(tensor482, cindex);
  task481->add_dep(task482);
  task482->add_dep(task314);
  deciq->add_task(task482);

  auto tensor483 = vector<shared_ptr<Tensor>>{I558, t2};
  auto task483 = make_shared<Task483>(tensor483, cindex);
  task482->add_dep(task483);
  task483->add_dep(task314);
  deciq->add_task(task483);

  vector<IndexRange> I639_index = {virt_, active_};
  auto I639 = make_shared<Tensor>(I639_index);
  auto tensor484 = vector<shared_ptr<Tensor>>{I556, t2, I639};
  auto task484 = make_shared<Task484>(tensor484, cindex);
  task480->add_dep(task484);
  task484->add_dep(task314);
  deciq->add_task(task484);

  auto tensor485 = vector<shared_ptr<Tensor>>{I639, f1_, l2};
  auto task485 = make_shared<Task485>(tensor485, cindex);
  task484->add_dep(task485);
  task485->add_dep(task314);
  deciq->add_task(task485);

  vector<IndexRange> I643_index = {virt_, active_};
  auto I643 = make_shared<Tensor>(I643_index);
  auto tensor486 = vector<shared_ptr<Tensor>>{I556, t2, I643};
  auto task486 = make_shared<Task486>(tensor486, cindex);
  task480->add_dep(task486);
  task486->add_dep(task314);
  deciq->add_task(task486);

  auto tensor487 = vector<shared_ptr<Tensor>>{I643, f1_, l2};
  auto task487 = make_shared<Task487>(tensor487, cindex);
  task486->add_dep(task487);
  task487->add_dep(task314);
  deciq->add_task(task487);

  vector<IndexRange> I685_index = {virt_, virt_, active_, active_};
  auto I685 = make_shared<Tensor>(I685_index);
  auto tensor488 = vector<shared_ptr<Tensor>>{I556, t2, I685};
  auto task488 = make_shared<Task488>(tensor488, cindex);
  task480->add_dep(task488);
  task488->add_dep(task314);
  deciq->add_task(task488);

  auto tensor489 = vector<shared_ptr<Tensor>>{I685, f1_, l2};
  auto task489 = make_shared<Task489>(tensor489, cindex);
  task488->add_dep(task489);
  task489->add_dep(task314);
  deciq->add_task(task489);

  vector<IndexRange> I693_index = {active_, virt_, virt_, active_};
  auto I693 = make_shared<Tensor>(I693_index);
  auto tensor490 = vector<shared_ptr<Tensor>>{I556, l2, I693};
  auto task490 = make_shared<Task490>(tensor490, cindex);
  task480->add_dep(task490);
  task490->add_dep(task314);
  deciq->add_task(task490);

  auto tensor491 = vector<shared_ptr<Tensor>>{I693, f1_, t2};
  auto task491 = make_shared<Task491>(tensor491, cindex);
  task490->add_dep(task491);
  task491->add_dep(task314);
  deciq->add_task(task491);

  auto tensor492 = vector<shared_ptr<Tensor>>{I693, f1_, t2};
  auto task492 = make_shared<Task492>(tensor492, cindex);
  task490->add_dep(task492);
  task492->add_dep(task314);
  deciq->add_task(task492);

  auto tensor493 = vector<shared_ptr<Tensor>>{I556, t2, l2};
  auto task493 = make_shared<Task493>(tensor493, cindex, this->e0_);
  task480->add_dep(task493);
  task493->add_dep(task314);
  deciq->add_task(task493);

  vector<IndexRange> I592_index;
  auto I592 = make_shared<Tensor>(I592_index);
  auto tensor494 = vector<shared_ptr<Tensor>>{I324, Gamma179_(), I592};
  auto task494 = make_shared<Task494>(tensor494, cindex);
  task315->add_dep(task494);
  task494->add_dep(task314);
  deciq->add_task(task494);

  auto tensor495 = vector<shared_ptr<Tensor>>{I592, t2, l2};
  auto task495 = make_shared<Task495>(tensor495, cindex);
  task494->add_dep(task495);
  task495->add_dep(task314);
  deciq->add_task(task495);

  auto tensor496 = vector<shared_ptr<Tensor>>{I592, t2, l2};
  auto task496 = make_shared<Task496>(tensor496, cindex);
  task494->add_dep(task496);
  task496->add_dep(task314);
  deciq->add_task(task496);

  shared_ptr<Tensor> I598;
  if (diagonal) {
    vector<IndexRange> I598_index;
    I598 = make_shared<Tensor>(I598_index);
  }
  shared_ptr<Task497> task497;
  if (diagonal) {
    auto tensor497 = vector<shared_ptr<Tensor>>{I324, rdm0deriv_, I598};
    task497 = make_shared<Task497>(tensor497, cindex);
    task315->add_dep(task497);
    task497->add_dep(task314);
    deciq->add_task(task497);
  }

  shared_ptr<Tensor> I599;
  if (diagonal) {
    vector<IndexRange> I599_index = {closed_, closed_};
    I599 = make_shared<Tensor>(I599_index);
  }
  shared_ptr<Task498> task498;
  if (diagonal) {
    auto tensor498 = vector<shared_ptr<Tensor>>{I598, f1_, I599};
    task498 = make_shared<Task498>(tensor498, cindex);
    task497->add_dep(task498);
    task498->add_dep(task314);
    deciq->add_task(task498);
  }

  shared_ptr<Task499> task499;
  if (diagonal) {
    auto tensor499 = vector<shared_ptr<Tensor>>{I599, t2, l2};
    task499 = make_shared<Task499>(tensor499, cindex);
    task498->add_dep(task499);
    task499->add_dep(task314);
    deciq->add_task(task499);
  }

  shared_ptr<Tensor> I602;
  if (diagonal) {
    vector<IndexRange> I602_index;
    I602 = make_shared<Tensor>(I602_index);
  }
  shared_ptr<Task500> task500;
  if (diagonal) {
    auto tensor500 = vector<shared_ptr<Tensor>>{I324, rdm0deriv_, I602};
    task500 = make_shared<Task500>(tensor500, cindex);
    task315->add_dep(task500);
    task500->add_dep(task314);
    deciq->add_task(task500);
  }

  shared_ptr<Tensor> I603;
  if (diagonal) {
    vector<IndexRange> I603_index = {closed_, closed_};
    I603 = make_shared<Tensor>(I603_index);
  }
  shared_ptr<Task501> task501;
  if (diagonal) {
    auto tensor501 = vector<shared_ptr<Tensor>>{I602, f1_, I603};
    task501 = make_shared<Task501>(tensor501, cindex);
    task500->add_dep(task501);
    task501->add_dep(task314);
    deciq->add_task(task501);
  }

  shared_ptr<Task502> task502;
  if (diagonal) {
    auto tensor502 = vector<shared_ptr<Tensor>>{I603, t2, l2};
    task502 = make_shared<Task502>(tensor502, cindex);
    task501->add_dep(task502);
    task502->add_dep(task314);
    deciq->add_task(task502);
  }

  shared_ptr<Tensor> I606;
  if (diagonal) {
    vector<IndexRange> I606_index;
    I606 = make_shared<Tensor>(I606_index);
  }
  shared_ptr<Task503> task503;
  if (diagonal) {
    auto tensor503 = vector<shared_ptr<Tensor>>{I324, rdm0deriv_, I606};
    task503 = make_shared<Task503>(tensor503, cindex);
    task315->add_dep(task503);
    task503->add_dep(task314);
    deciq->add_task(task503);
  }

  shared_ptr<Tensor> I607;
  if (diagonal) {
    vector<IndexRange> I607_index = {virt_, virt_};
    I607 = make_shared<Tensor>(I607_index);
  }
  shared_ptr<Task504> task504;
  if (diagonal) {
    auto tensor504 = vector<shared_ptr<Tensor>>{I606, f1_, I607};
    task504 = make_shared<Task504>(tensor504, cindex);
    task503->add_dep(task504);
    task504->add_dep(task314);
    deciq->add_task(task504);
  }

  shared_ptr<Task505> task505;
  if (diagonal) {
    auto tensor505 = vector<shared_ptr<Tensor>>{I607, t2, l2};
    task505 = make_shared<Task505>(tensor505, cindex);
    task504->add_dep(task505);
    task505->add_dep(task314);
    deciq->add_task(task505);
  }

  shared_ptr<Task506> task506;
  if (diagonal) {
    auto tensor506 = vector<shared_ptr<Tensor>>{I606, t2, l2};
    task506 = make_shared<Task506>(tensor506, cindex, this->e0_);
    task503->add_dep(task506);
    task506->add_dep(task314);
    deciq->add_task(task506);
  }

  shared_ptr<Tensor> I610;
  if (diagonal) {
    vector<IndexRange> I610_index;
    I610 = make_shared<Tensor>(I610_index);
  }
  shared_ptr<Task507> task507;
  if (diagonal) {
    auto tensor507 = vector<shared_ptr<Tensor>>{I324, rdm0deriv_, I610};
    task507 = make_shared<Task507>(tensor507, cindex);
    task315->add_dep(task507);
    task507->add_dep(task314);
    deciq->add_task(task507);
  }

  shared_ptr<Tensor> I611;
  if (diagonal) {
    vector<IndexRange> I611_index = {virt_, virt_};
    I611 = make_shared<Tensor>(I611_index);
  }
  shared_ptr<Task508> task508;
  if (diagonal) {
    auto tensor508 = vector<shared_ptr<Tensor>>{I610, f1_, I611};
    task508 = make_shared<Task508>(tensor508, cindex);
    task507->add_dep(task508);
    task508->add_dep(task314);
    deciq->add_task(task508);
  }

  shared_ptr<Task509> task509;
  if (diagonal) {
    auto tensor509 = vector<shared_ptr<Tensor>>{I611, t2, l2};
    task509 = make_shared<Task509>(tensor509, cindex);
    task508->add_dep(task509);
    task509->add_dep(task314);
    deciq->add_task(task509);
  }

  vector<IndexRange> I654_index = {active_, active_};
  auto I654 = make_shared<Tensor>(I654_index);
  auto tensor510 = vector<shared_ptr<Tensor>>{I324, Gamma191_(), I654};
  auto task510 = make_shared<Task510>(tensor510, cindex);
  task315->add_dep(task510);
  task510->add_dep(task314);
  deciq->add_task(task510);

  auto tensor511 = vector<shared_ptr<Tensor>>{I654, t2, l2};
  auto task511 = make_shared<Task511>(tensor511, cindex);
  task510->add_dep(task511);
  task511->add_dep(task314);
  deciq->add_task(task511);

  auto tensor512 = vector<shared_ptr<Tensor>>{I654, t2, l2};
  auto task512 = make_shared<Task512>(tensor512, cindex);
  task510->add_dep(task512);
  task512->add_dep(task314);
  deciq->add_task(task512);

  vector<IndexRange> I696_index = {active_, active_, active_, active_};
  auto I696 = make_shared<Tensor>(I696_index);
  auto tensor513 = vector<shared_ptr<Tensor>>{I324, Gamma202_(), I696};
  auto task513 = make_shared<Task513>(tensor513, cindex);
  task315->add_dep(task513);
  task513->add_dep(task314);
  deciq->add_task(task513);

  auto tensor514 = vector<shared_ptr<Tensor>>{I696, t2, l2};
  auto task514 = make_shared<Task514>(tensor514, cindex);
  task513->add_dep(task514);
  task514->add_dep(task314);
  deciq->add_task(task514);

  shared_ptr<Tensor> I730;
  if (diagonal) {
    vector<IndexRange> I730_index;
    I730 = make_shared<Tensor>(I730_index);
  }
  shared_ptr<Task515> task515;
  if (diagonal) {
    auto tensor515 = vector<shared_ptr<Tensor>>{I324, rdm0deriv_, I730};
    task515 = make_shared<Task515>(tensor515, cindex);
    task315->add_dep(task515);
    task515->add_dep(task314);
    deciq->add_task(task515);
  }

  shared_ptr<Task516> task516;
  if (diagonal) {
    auto tensor516 = vector<shared_ptr<Tensor>>{I730, t2, l2};
    task516 = make_shared<Task516>(tensor516, cindex, this->e0_);
    task515->add_dep(task516);
    task516->add_dep(task314);
    deciq->add_task(task516);
  }

  return deciq;
}


#endif
