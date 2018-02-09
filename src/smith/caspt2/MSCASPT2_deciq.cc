//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_deciq.cc
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
#include <src/smith/caspt2/MSCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_deciq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto deciq = make_shared<Queue>();
  auto tensor282 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci};
  auto task282 = make_shared<Task282>(tensor282, reset);
  deciq->add_task(task282);

  vector<IndexRange> I325_index = {active_, active_, active_, active_};
  auto I325 = make_shared<Tensor>(I325_index);
  auto tensor284 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I325, f1_};
  auto task284 = make_shared<Task284>(tensor284, pindex);
  task284->add_dep(task282);
  deciq->add_task(task284);

  auto tensor285 = vector<shared_ptr<Tensor>>{I325, t2, l2};
  auto task285 = make_shared<Task285>(tensor285, pindex);
  task284->add_dep(task285);
  task285->add_dep(task282);
  deciq->add_task(task285);

  vector<IndexRange> I328_index = {active_, active_, active_, active_};
  auto I328 = make_shared<Tensor>(I328_index);
  auto tensor286 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I328};
  auto task286 = make_shared<Task286>(tensor286, pindex);
  task286->add_dep(task282);
  deciq->add_task(task286);

  vector<IndexRange> I329_index = {closed_, active_, active_, closed_};
  auto I329 = make_shared<Tensor>(I329_index);
  auto tensor287 = vector<shared_ptr<Tensor>>{I328, l2, I329};
  auto task287 = make_shared<Task287>(tensor287, pindex);
  task286->add_dep(task287);
  task287->add_dep(task282);
  deciq->add_task(task287);

  auto tensor288 = vector<shared_ptr<Tensor>>{I329, f1_, t2};
  auto task288 = make_shared<Task288>(tensor288, pindex);
  task287->add_dep(task288);
  task288->add_dep(task282);
  deciq->add_task(task288);

  vector<IndexRange> I332_index = {active_, active_, active_, active_, active_, active_};
  auto I332 = make_shared<Tensor>(I332_index);
  auto tensor289 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I332};
  auto task289 = make_shared<Task289>(tensor289, pindex);
  task289->add_dep(task282);
  deciq->add_task(task289);

  vector<IndexRange> I333_index = {active_, active_, closed_, active_};
  auto I333 = make_shared<Tensor>(I333_index);
  auto tensor290 = vector<shared_ptr<Tensor>>{I332, t2, I333};
  auto task290 = make_shared<Task290>(tensor290, pindex);
  task289->add_dep(task290);
  task290->add_dep(task282);
  deciq->add_task(task290);

  auto tensor291 = vector<shared_ptr<Tensor>>{I333, f1_, l2};
  auto task291 = make_shared<Task291>(tensor291, pindex);
  task290->add_dep(task291);
  task291->add_dep(task282);
  deciq->add_task(task291);

  vector<IndexRange> I336_index = {active_, active_, active_, active_};
  auto I336 = make_shared<Tensor>(I336_index);
  auto tensor292 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I336};
  auto task292 = make_shared<Task292>(tensor292, pindex);
  task292->add_dep(task282);
  deciq->add_task(task292);

  vector<IndexRange> I337_index = {closed_, closed_, active_, active_};
  auto I337 = make_shared<Tensor>(I337_index);
  auto tensor293 = vector<shared_ptr<Tensor>>{I336, l2, I337};
  auto task293 = make_shared<Task293>(tensor293, pindex);
  task292->add_dep(task293);
  task293->add_dep(task282);
  deciq->add_task(task293);

  auto tensor294 = vector<shared_ptr<Tensor>>{I337, f1_, t2};
  auto task294 = make_shared<Task294>(tensor294, pindex);
  task293->add_dep(task294);
  task294->add_dep(task282);
  deciq->add_task(task294);

  vector<IndexRange> I368_index = {active_, closed_, closed_, active_};
  auto I368 = make_shared<Tensor>(I368_index);
  auto tensor295 = vector<shared_ptr<Tensor>>{I336, t2, I368};
  auto task295 = make_shared<Task295>(tensor295, pindex);
  task292->add_dep(task295);
  task295->add_dep(task282);
  deciq->add_task(task295);

  auto tensor296 = vector<shared_ptr<Tensor>>{I368, f1_, l2};
  auto task296 = make_shared<Task296>(tensor296, pindex);
  task295->add_dep(task296);
  task296->add_dep(task282);
  deciq->add_task(task296);

  vector<IndexRange> I340_index = {active_, active_, active_, active_, active_, active_};
  auto I340 = make_shared<Tensor>(I340_index);
  auto tensor297 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I340};
  auto task297 = make_shared<Task297>(tensor297, pindex);
  task297->add_dep(task282);
  deciq->add_task(task297);

  vector<IndexRange> I341_index = {closed_, active_, active_, active_};
  auto I341 = make_shared<Tensor>(I341_index);
  auto tensor298 = vector<shared_ptr<Tensor>>{I340, l2, I341};
  auto task298 = make_shared<Task298>(tensor298, pindex);
  task297->add_dep(task298);
  task298->add_dep(task282);
  deciq->add_task(task298);

  auto tensor299 = vector<shared_ptr<Tensor>>{I341, f1_, t2};
  auto task299 = make_shared<Task299>(tensor299, pindex);
  task298->add_dep(task299);
  task299->add_dep(task282);
  deciq->add_task(task299);

  vector<IndexRange> I344_index = {active_, active_, active_, active_, active_, active_};
  auto I344 = make_shared<Tensor>(I344_index);
  auto tensor300 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I344, f1_};
  auto task300 = make_shared<Task300>(tensor300, pindex);
  task300->add_dep(task282);
  deciq->add_task(task300);

  auto tensor301 = vector<shared_ptr<Tensor>>{I344, t2, l2};
  auto task301 = make_shared<Task301>(tensor301, pindex);
  task300->add_dep(task301);
  task301->add_dep(task282);
  deciq->add_task(task301);

  vector<IndexRange> I347_index = {active_, active_, active_, active_, active_, active_};
  auto I347 = make_shared<Tensor>(I347_index);
  auto tensor302 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I347};
  auto task302 = make_shared<Task302>(tensor302, pindex);
  task302->add_dep(task282);
  deciq->add_task(task302);

  vector<IndexRange> I348_index = {active_, active_, active_, closed_};
  auto I348 = make_shared<Tensor>(I348_index);
  auto tensor303 = vector<shared_ptr<Tensor>>{I347, l2, I348};
  auto task303 = make_shared<Task303>(tensor303, pindex);
  task302->add_dep(task303);
  task303->add_dep(task282);
  deciq->add_task(task303);

  auto tensor304 = vector<shared_ptr<Tensor>>{I348, f1_, t2};
  auto task304 = make_shared<Task304>(tensor304, pindex);
  task303->add_dep(task304);
  task304->add_dep(task282);
  deciq->add_task(task304);

  auto tensor305 = vector<shared_ptr<Tensor>>{I348, f1_, t2};
  auto task305 = make_shared<Task305>(tensor305, pindex);
  task303->add_dep(task305);
  task305->add_dep(task282);
  deciq->add_task(task305);

  vector<IndexRange> I488_index = {active_, active_, closed_, active_};
  auto I488 = make_shared<Tensor>(I488_index);
  auto tensor306 = vector<shared_ptr<Tensor>>{I347, t2, I488};
  auto task306 = make_shared<Task306>(tensor306, pindex);
  task302->add_dep(task306);
  task306->add_dep(task282);
  deciq->add_task(task306);

  auto tensor307 = vector<shared_ptr<Tensor>>{I488, f1_, l2};
  auto task307 = make_shared<Task307>(tensor307, pindex);
  task306->add_dep(task307);
  task307->add_dep(task282);
  deciq->add_task(task307);

  vector<IndexRange> I351_index = {active_, active_, active_, active_};
  auto I351 = make_shared<Tensor>(I351_index);
  auto tensor308 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I351};
  auto task308 = make_shared<Task308>(tensor308, pindex);
  task308->add_dep(task282);
  deciq->add_task(task308);

  vector<IndexRange> I352_index = {closed_, active_};
  auto I352 = make_shared<Tensor>(I352_index);
  auto tensor309 = vector<shared_ptr<Tensor>>{I351, l2, I352};
  auto task309 = make_shared<Task309>(tensor309, pindex);
  task308->add_dep(task309);
  task309->add_dep(task282);
  deciq->add_task(task309);

  vector<IndexRange> I353_index = {closed_, virt_, closed_, active_};
  auto I353 = make_shared<Tensor>(I353_index);
  auto tensor310 = vector<shared_ptr<Tensor>>{I352, f1_, I353};
  auto task310 = make_shared<Task310>(tensor310, pindex);
  task309->add_dep(task310);
  task310->add_dep(task282);
  deciq->add_task(task310);

  auto tensor311 = vector<shared_ptr<Tensor>>{I353, t2};
  auto task311 = make_shared<Task311>(tensor311, pindex);
  task310->add_dep(task311);
  task311->add_dep(task282);
  deciq->add_task(task311);

  vector<IndexRange> I442_index = {closed_, virt_, active_, active_};
  auto I442 = make_shared<Tensor>(I442_index);
  auto tensor312 = vector<shared_ptr<Tensor>>{I351, l2, I442};
  auto task312 = make_shared<Task312>(tensor312, pindex);
  task308->add_dep(task312);
  task312->add_dep(task282);
  deciq->add_task(task312);

  auto tensor313 = vector<shared_ptr<Tensor>>{I442, f1_, t2};
  auto task313 = make_shared<Task313>(tensor313, pindex);
  task312->add_dep(task313);
  task313->add_dep(task282);
  deciq->add_task(task313);

  vector<IndexRange> I492_index = {virt_, closed_, active_, active_};
  auto I492 = make_shared<Tensor>(I492_index);
  auto tensor314 = vector<shared_ptr<Tensor>>{I351, l2, I492};
  auto task314 = make_shared<Task314>(tensor314, pindex);
  task308->add_dep(task314);
  task314->add_dep(task282);
  deciq->add_task(task314);

  vector<IndexRange> I493_index = {closed_, virt_, closed_, active_};
  auto I493 = make_shared<Tensor>(I493_index);
  auto tensor315 = vector<shared_ptr<Tensor>>{I492, f1_, I493};
  auto task315 = make_shared<Task315>(tensor315, pindex);
  task314->add_dep(task315);
  task315->add_dep(task282);
  deciq->add_task(task315);

  auto tensor316 = vector<shared_ptr<Tensor>>{I493, t2};
  auto task316 = make_shared<Task316>(tensor316, pindex);
  task315->add_dep(task316);
  task316->add_dep(task282);
  deciq->add_task(task316);

  vector<IndexRange> I359_index = {active_, active_, active_, active_, active_, active_};
  auto I359 = make_shared<Tensor>(I359_index);
  auto tensor317 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I359};
  auto task317 = make_shared<Task317>(tensor317, pindex);
  task317->add_dep(task282);
  deciq->add_task(task317);

  vector<IndexRange> I360_index = {active_, closed_, active_, active_};
  auto I360 = make_shared<Tensor>(I360_index);
  auto tensor318 = vector<shared_ptr<Tensor>>{I359, l2, I360};
  auto task318 = make_shared<Task318>(tensor318, pindex);
  task317->add_dep(task318);
  task318->add_dep(task282);
  deciq->add_task(task318);

  auto tensor319 = vector<shared_ptr<Tensor>>{I360, f1_, t2};
  auto task319 = make_shared<Task319>(tensor319, pindex);
  task318->add_dep(task319);
  task319->add_dep(task282);
  deciq->add_task(task319);

  vector<IndexRange> I371_index = {active_, active_, active_, active_};
  auto I371 = make_shared<Tensor>(I371_index);
  auto tensor320 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I371};
  auto task320 = make_shared<Task320>(tensor320, pindex);
  task320->add_dep(task282);
  deciq->add_task(task320);

  vector<IndexRange> I372_index = {active_, closed_};
  auto I372 = make_shared<Tensor>(I372_index);
  auto tensor321 = vector<shared_ptr<Tensor>>{I371, t2, I372};
  auto task321 = make_shared<Task321>(tensor321, pindex);
  task320->add_dep(task321);
  task321->add_dep(task282);
  deciq->add_task(task321);

  auto tensor322 = vector<shared_ptr<Tensor>>{I372, f1_, l2};
  auto task322 = make_shared<Task322>(tensor322, pindex);
  task321->add_dep(task322);
  task322->add_dep(task282);
  deciq->add_task(task322);

  vector<IndexRange> I376_index = {active_, closed_};
  auto I376 = make_shared<Tensor>(I376_index);
  auto tensor323 = vector<shared_ptr<Tensor>>{I371, t2, I376};
  auto task323 = make_shared<Task323>(tensor323, pindex);
  task320->add_dep(task323);
  task323->add_dep(task282);
  deciq->add_task(task323);

  auto tensor324 = vector<shared_ptr<Tensor>>{I376, f1_, l2};
  auto task324 = make_shared<Task324>(tensor324, pindex);
  task323->add_dep(task324);
  task324->add_dep(task282);
  deciq->add_task(task324);

  vector<IndexRange> I414_index = {active_, virt_, closed_, active_};
  auto I414 = make_shared<Tensor>(I414_index);
  auto tensor325 = vector<shared_ptr<Tensor>>{I371, t2, I414};
  auto task325 = make_shared<Task325>(tensor325, pindex);
  task320->add_dep(task325);
  task325->add_dep(task282);
  deciq->add_task(task325);

  auto tensor326 = vector<shared_ptr<Tensor>>{I414, f1_, l2};
  auto task326 = make_shared<Task326>(tensor326, pindex);
  task325->add_dep(task326);
  task326->add_dep(task282);
  deciq->add_task(task326);

  vector<IndexRange> I418_index = {active_, closed_, virt_, active_};
  auto I418 = make_shared<Tensor>(I418_index);
  auto tensor327 = vector<shared_ptr<Tensor>>{I371, t2, I418};
  auto task327 = make_shared<Task327>(tensor327, pindex);
  task320->add_dep(task327);
  task327->add_dep(task282);
  deciq->add_task(task327);

  auto tensor328 = vector<shared_ptr<Tensor>>{I418, f1_, l2};
  auto task328 = make_shared<Task328>(tensor328, pindex);
  task327->add_dep(task328);
  task328->add_dep(task282);
  deciq->add_task(task328);

  vector<IndexRange> I422_index = {active_, virt_, closed_, active_};
  auto I422 = make_shared<Tensor>(I422_index);
  auto tensor329 = vector<shared_ptr<Tensor>>{I371, t2, I422};
  auto task329 = make_shared<Task329>(tensor329, pindex);
  task320->add_dep(task329);
  task329->add_dep(task282);
  deciq->add_task(task329);

  auto tensor330 = vector<shared_ptr<Tensor>>{I422, f1_, l2};
  auto task330 = make_shared<Task330>(tensor330, pindex);
  task329->add_dep(task330);
  task330->add_dep(task282);
  deciq->add_task(task330);

  vector<IndexRange> I379_index = {active_, active_};
  auto I379 = make_shared<Tensor>(I379_index);
  auto tensor331 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I379, f1_};
  auto task331 = make_shared<Task331>(tensor331, pindex);
  task331->add_dep(task282);
  deciq->add_task(task331);

  auto tensor332 = vector<shared_ptr<Tensor>>{I379, t2, l2};
  auto task332 = make_shared<Task332>(tensor332, pindex);
  task331->add_dep(task332);
  task332->add_dep(task282);
  deciq->add_task(task332);

  auto tensor333 = vector<shared_ptr<Tensor>>{I379, t2, l2};
  auto task333 = make_shared<Task333>(tensor333, pindex);
  task331->add_dep(task333);
  task333->add_dep(task282);
  deciq->add_task(task333);

  vector<IndexRange> I385_index = {active_, active_};
  auto I385 = make_shared<Tensor>(I385_index);
  auto tensor334 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I385};
  auto task334 = make_shared<Task334>(tensor334, pindex);
  task334->add_dep(task282);
  deciq->add_task(task334);

  vector<IndexRange> I386_index = {virt_, closed_, active_, closed_};
  auto I386 = make_shared<Tensor>(I386_index);
  auto tensor335 = vector<shared_ptr<Tensor>>{I385, l2, I386};
  auto task335 = make_shared<Task335>(tensor335, pindex);
  task334->add_dep(task335);
  task335->add_dep(task282);
  deciq->add_task(task335);

  vector<IndexRange> I387_index = {closed_, virt_, closed_, active_};
  auto I387 = make_shared<Tensor>(I387_index);
  auto tensor336 = vector<shared_ptr<Tensor>>{I386, f1_, I387};
  auto task336 = make_shared<Task336>(tensor336, pindex);
  task335->add_dep(task336);
  task336->add_dep(task282);
  deciq->add_task(task336);

  auto tensor337 = vector<shared_ptr<Tensor>>{I387, t2};
  auto task337 = make_shared<Task337>(tensor337, pindex);
  task336->add_dep(task337);
  task337->add_dep(task282);
  deciq->add_task(task337);

  vector<IndexRange> I395_index = {closed_, virt_, closed_, active_};
  auto I395 = make_shared<Tensor>(I395_index);
  auto tensor338 = vector<shared_ptr<Tensor>>{I386, f1_, I395};
  auto task338 = make_shared<Task338>(tensor338, pindex);
  task335->add_dep(task338);
  task338->add_dep(task282);
  deciq->add_task(task338);

  auto tensor339 = vector<shared_ptr<Tensor>>{I395, t2};
  auto task339 = make_shared<Task339>(tensor339, pindex);
  task338->add_dep(task339);
  task339->add_dep(task282);
  deciq->add_task(task339);

  vector<IndexRange> I399_index = {closed_, virt_, closed_, active_};
  auto I399 = make_shared<Tensor>(I399_index);
  auto tensor340 = vector<shared_ptr<Tensor>>{I386, f1_, I399};
  auto task340 = make_shared<Task340>(tensor340, pindex);
  task335->add_dep(task340);
  task340->add_dep(task282);
  deciq->add_task(task340);

  auto tensor341 = vector<shared_ptr<Tensor>>{I399, t2};
  auto task341 = make_shared<Task341>(tensor341, pindex);
  task340->add_dep(task341);
  task341->add_dep(task282);
  deciq->add_task(task341);

  vector<IndexRange> I426_index = {active_, virt_};
  auto I426 = make_shared<Tensor>(I426_index);
  auto tensor342 = vector<shared_ptr<Tensor>>{I385, f1_, I426};
  auto task342 = make_shared<Task342>(tensor342, pindex);
  task334->add_dep(task342);
  task342->add_dep(task282);
  deciq->add_task(task342);

  auto tensor343 = vector<shared_ptr<Tensor>>{I426, t2, l2};
  auto task343 = make_shared<Task343>(tensor343, pindex);
  task342->add_dep(task343);
  task343->add_dep(task282);
  deciq->add_task(task343);

  auto tensor344 = vector<shared_ptr<Tensor>>{I426, t2, l2};
  auto task344 = make_shared<Task344>(tensor344, pindex);
  task342->add_dep(task344);
  task344->add_dep(task282);
  deciq->add_task(task344);

  vector<IndexRange> I569_index = {virt_, active_};
  auto I569 = make_shared<Tensor>(I569_index);
  auto tensor345 = vector<shared_ptr<Tensor>>{I385, f1_, I569};
  auto task345 = make_shared<Task345>(tensor345, pindex);
  task334->add_dep(task345);
  task345->add_dep(task282);
  deciq->add_task(task345);

  auto tensor346 = vector<shared_ptr<Tensor>>{I569, t2, l2};
  auto task346 = make_shared<Task346>(tensor346, pindex);
  task345->add_dep(task346);
  task346->add_dep(task282);
  deciq->add_task(task346);

  vector<IndexRange> I573_index = {virt_, active_};
  auto I573 = make_shared<Tensor>(I573_index);
  auto tensor347 = vector<shared_ptr<Tensor>>{I385, f1_, I573};
  auto task347 = make_shared<Task347>(tensor347, pindex);
  task334->add_dep(task347);
  task347->add_dep(task282);
  deciq->add_task(task347);

  auto tensor348 = vector<shared_ptr<Tensor>>{I573, t2, l2};
  auto task348 = make_shared<Task348>(tensor348, pindex);
  task347->add_dep(task348);
  task348->add_dep(task282);
  deciq->add_task(task348);

  vector<IndexRange> I409_index = {active_, active_, active_, active_};
  auto I409 = make_shared<Tensor>(I409_index);
  auto tensor349 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I409};
  auto task349 = make_shared<Task349>(tensor349, pindex);
  task349->add_dep(task282);
  deciq->add_task(task349);

  vector<IndexRange> I410_index = {active_, closed_, virt_, active_};
  auto I410 = make_shared<Tensor>(I410_index);
  auto tensor350 = vector<shared_ptr<Tensor>>{I409, t2, I410};
  auto task350 = make_shared<Task350>(tensor350, pindex);
  task349->add_dep(task350);
  task350->add_dep(task282);
  deciq->add_task(task350);

  auto tensor351 = vector<shared_ptr<Tensor>>{I410, f1_, l2};
  auto task351 = make_shared<Task351>(tensor351, pindex);
  task350->add_dep(task351);
  task351->add_dep(task282);
  deciq->add_task(task351);

  vector<IndexRange> I433_index = {active_, active_, active_, active_, active_, active_};
  auto I433 = make_shared<Tensor>(I433_index);
  auto tensor352 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I433};
  auto task352 = make_shared<Task352>(tensor352, pindex);
  task352->add_dep(task282);
  deciq->add_task(task352);

  vector<IndexRange> I434_index = {active_, closed_, active_, active_};
  auto I434 = make_shared<Tensor>(I434_index);
  auto tensor353 = vector<shared_ptr<Tensor>>{I433, t2, I434};
  auto task353 = make_shared<Task353>(tensor353, pindex);
  task352->add_dep(task353);
  task353->add_dep(task282);
  deciq->add_task(task353);

  auto tensor354 = vector<shared_ptr<Tensor>>{I434, f1_, l2};
  auto task354 = make_shared<Task354>(tensor354, pindex);
  task353->add_dep(task354);
  task354->add_dep(task282);
  deciq->add_task(task354);

  vector<IndexRange> I437_index = {active_, active_, active_, active_};
  auto I437 = make_shared<Tensor>(I437_index);
  auto tensor355 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I437};
  auto task355 = make_shared<Task355>(tensor355, pindex);
  task355->add_dep(task282);
  deciq->add_task(task355);

  vector<IndexRange> I438_index = {virt_, closed_, active_, active_};
  auto I438 = make_shared<Tensor>(I438_index);
  auto tensor356 = vector<shared_ptr<Tensor>>{I437, l2, I438};
  auto task356 = make_shared<Task356>(tensor356, pindex);
  task355->add_dep(task356);
  task356->add_dep(task282);
  deciq->add_task(task356);

  auto tensor357 = vector<shared_ptr<Tensor>>{I438, f1_, t2};
  auto task357 = make_shared<Task357>(tensor357, pindex);
  task356->add_dep(task357);
  task357->add_dep(task282);
  deciq->add_task(task357);

  vector<IndexRange> I445_index = {active_, active_, active_, active_};
  auto I445 = make_shared<Tensor>(I445_index);
  auto tensor358 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I445, f1_};
  auto task358 = make_shared<Task358>(tensor358, pindex);
  task358->add_dep(task282);
  deciq->add_task(task358);

  auto tensor359 = vector<shared_ptr<Tensor>>{I445, t2, l2};
  auto task359 = make_shared<Task359>(tensor359, pindex);
  task358->add_dep(task359);
  task359->add_dep(task282);
  deciq->add_task(task359);

  vector<IndexRange> I448_index = {active_, active_, active_, active_};
  auto I448 = make_shared<Tensor>(I448_index);
  auto tensor360 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I448};
  auto task360 = make_shared<Task360>(tensor360, pindex);
  task360->add_dep(task282);
  deciq->add_task(task360);

  vector<IndexRange> I449_index = {active_, virt_, active_, closed_};
  auto I449 = make_shared<Tensor>(I449_index);
  auto tensor361 = vector<shared_ptr<Tensor>>{I448, l2, I449};
  auto task361 = make_shared<Task361>(tensor361, pindex);
  task360->add_dep(task361);
  task361->add_dep(task282);
  deciq->add_task(task361);

  auto tensor362 = vector<shared_ptr<Tensor>>{I449, f1_, t2};
  auto task362 = make_shared<Task362>(tensor362, pindex);
  task361->add_dep(task362);
  task362->add_dep(task282);
  deciq->add_task(task362);

  auto tensor363 = vector<shared_ptr<Tensor>>{I449, f1_, t2};
  auto task363 = make_shared<Task363>(tensor363, pindex);
  task361->add_dep(task363);
  task363->add_dep(task282);
  deciq->add_task(task363);

  auto tensor364 = vector<shared_ptr<Tensor>>{I449, f1_, t2};
  auto task364 = make_shared<Task364>(tensor364, pindex);
  task361->add_dep(task364);
  task364->add_dep(task282);
  deciq->add_task(task364);

  vector<IndexRange> I627_index = {closed_, virt_, active_, active_};
  auto I627 = make_shared<Tensor>(I627_index);
  auto tensor365 = vector<shared_ptr<Tensor>>{I448, t2, I627};
  auto task365 = make_shared<Task365>(tensor365, pindex);
  task360->add_dep(task365);
  task365->add_dep(task282);
  deciq->add_task(task365);

  auto tensor366 = vector<shared_ptr<Tensor>>{I627, f1_, l2};
  auto task366 = make_shared<Task366>(tensor366, pindex);
  task365->add_dep(task366);
  task366->add_dep(task282);
  deciq->add_task(task366);

  vector<IndexRange> I456_index = {active_, active_, active_, active_};
  auto I456 = make_shared<Tensor>(I456_index);
  auto tensor367 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I456, f1_};
  auto task367 = make_shared<Task367>(tensor367, pindex);
  task367->add_dep(task282);
  deciq->add_task(task367);

  auto tensor368 = vector<shared_ptr<Tensor>>{I456, t2, l2};
  auto task368 = make_shared<Task368>(tensor368, pindex);
  task367->add_dep(task368);
  task368->add_dep(task282);
  deciq->add_task(task368);

  auto tensor369 = vector<shared_ptr<Tensor>>{I456, t2, l2};
  auto task369 = make_shared<Task369>(tensor369, pindex);
  task367->add_dep(task369);
  task369->add_dep(task282);
  deciq->add_task(task369);

  auto tensor370 = vector<shared_ptr<Tensor>>{I456, t2, l2};
  auto task370 = make_shared<Task370>(tensor370, pindex);
  task367->add_dep(task370);
  task370->add_dep(task282);
  deciq->add_task(task370);

  vector<IndexRange> I459_index = {active_, active_, active_, active_};
  auto I459 = make_shared<Tensor>(I459_index);
  auto tensor371 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I459};
  auto task371 = make_shared<Task371>(tensor371, pindex);
  task371->add_dep(task282);
  deciq->add_task(task371);

  vector<IndexRange> I460_index = {virt_, active_, active_, closed_};
  auto I460 = make_shared<Tensor>(I460_index);
  auto tensor372 = vector<shared_ptr<Tensor>>{I459, l2, I460};
  auto task372 = make_shared<Task372>(tensor372, pindex);
  task371->add_dep(task372);
  task372->add_dep(task282);
  deciq->add_task(task372);

  auto tensor373 = vector<shared_ptr<Tensor>>{I460, f1_, t2};
  auto task373 = make_shared<Task373>(tensor373, pindex);
  task372->add_dep(task373);
  task373->add_dep(task282);
  deciq->add_task(task373);

  auto tensor374 = vector<shared_ptr<Tensor>>{I460, f1_, t2};
  auto task374 = make_shared<Task374>(tensor374, pindex);
  task372->add_dep(task374);
  task374->add_dep(task282);
  deciq->add_task(task374);

  auto tensor375 = vector<shared_ptr<Tensor>>{I460, f1_, t2};
  auto task375 = make_shared<Task375>(tensor375, pindex);
  task372->add_dep(task375);
  task375->add_dep(task282);
  deciq->add_task(task375);

  vector<IndexRange> I503_index = {active_, virt_, active_, closed_};
  auto I503 = make_shared<Tensor>(I503_index);
  auto tensor376 = vector<shared_ptr<Tensor>>{I459, l2, I503};
  auto task376 = make_shared<Task376>(tensor376, pindex);
  task371->add_dep(task376);
  task376->add_dep(task282);
  deciq->add_task(task376);

  vector<IndexRange> I504_index = {active_, virt_, closed_, active_};
  auto I504 = make_shared<Tensor>(I504_index);
  auto tensor377 = vector<shared_ptr<Tensor>>{I503, f1_, I504};
  auto task377 = make_shared<Task377>(tensor377, pindex);
  task376->add_dep(task377);
  task377->add_dep(task282);
  deciq->add_task(task377);

  auto tensor378 = vector<shared_ptr<Tensor>>{I504, t2};
  auto task378 = make_shared<Task378>(tensor378, pindex);
  task377->add_dep(task378);
  task378->add_dep(task282);
  deciq->add_task(task378);

  vector<IndexRange> I508_index = {active_, virt_, closed_, active_};
  auto I508 = make_shared<Tensor>(I508_index);
  auto tensor379 = vector<shared_ptr<Tensor>>{I503, f1_, I508};
  auto task379 = make_shared<Task379>(tensor379, pindex);
  task376->add_dep(task379);
  task379->add_dep(task282);
  deciq->add_task(task379);

  auto tensor380 = vector<shared_ptr<Tensor>>{I508, t2};
  auto task380 = make_shared<Task380>(tensor380, pindex);
  task379->add_dep(task380);
  task380->add_dep(task282);
  deciq->add_task(task380);

  vector<IndexRange> I535_index = {active_, virt_, closed_, virt_};
  auto I535 = make_shared<Tensor>(I535_index);
  auto tensor381 = vector<shared_ptr<Tensor>>{I503, f1_, I535};
  auto task381 = make_shared<Task381>(tensor381, pindex);
  task376->add_dep(task381);
  task381->add_dep(task282);
  deciq->add_task(task381);

  auto tensor382 = vector<shared_ptr<Tensor>>{I535, t2};
  auto task382 = make_shared<Task382>(tensor382, pindex);
  task381->add_dep(task382);
  task382->add_dep(task282);
  deciq->add_task(task382);

  vector<IndexRange> I623_index = {virt_, closed_, active_, active_};
  auto I623 = make_shared<Tensor>(I623_index);
  auto tensor383 = vector<shared_ptr<Tensor>>{I459, t2, I623};
  auto task383 = make_shared<Task383>(tensor383, pindex);
  task371->add_dep(task383);
  task383->add_dep(task282);
  deciq->add_task(task383);

  auto tensor384 = vector<shared_ptr<Tensor>>{I623, f1_, l2};
  auto task384 = make_shared<Task384>(tensor384, pindex);
  task383->add_dep(task384);
  task384->add_dep(task282);
  deciq->add_task(task384);

  vector<IndexRange> I631_index = {virt_, closed_, active_, active_};
  auto I631 = make_shared<Tensor>(I631_index);
  auto tensor385 = vector<shared_ptr<Tensor>>{I459, t2, I631};
  auto task385 = make_shared<Task385>(tensor385, pindex);
  task371->add_dep(task385);
  task385->add_dep(task282);
  deciq->add_task(task385);

  auto tensor386 = vector<shared_ptr<Tensor>>{I631, f1_, l2};
  auto task386 = make_shared<Task386>(tensor386, pindex);
  task385->add_dep(task386);
  task386->add_dep(task282);
  deciq->add_task(task386);

  vector<IndexRange> I635_index = {closed_, virt_, active_, active_};
  auto I635 = make_shared<Tensor>(I635_index);
  auto tensor387 = vector<shared_ptr<Tensor>>{I459, t2, I635};
  auto task387 = make_shared<Task387>(tensor387, pindex);
  task371->add_dep(task387);
  task387->add_dep(task282);
  deciq->add_task(task387);

  auto tensor388 = vector<shared_ptr<Tensor>>{I635, f1_, l2};
  auto task388 = make_shared<Task388>(tensor388, pindex);
  task387->add_dep(task388);
  task388->add_dep(task282);
  deciq->add_task(task388);

  vector<IndexRange> I467_index = {active_, active_, active_, active_, active_, active_};
  auto I467 = make_shared<Tensor>(I467_index);
  auto tensor389 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I467};
  auto task389 = make_shared<Task389>(tensor389, pindex);
  task389->add_dep(task282);
  deciq->add_task(task389);

  vector<IndexRange> I468_index = {active_, virt_, active_, active_};
  auto I468 = make_shared<Tensor>(I468_index);
  auto tensor390 = vector<shared_ptr<Tensor>>{I467, t2, I468};
  auto task390 = make_shared<Task390>(tensor390, pindex);
  task389->add_dep(task390);
  task390->add_dep(task282);
  deciq->add_task(task390);

  auto tensor391 = vector<shared_ptr<Tensor>>{I468, f1_, l2};
  auto task391 = make_shared<Task391>(tensor391, pindex);
  task390->add_dep(task391);
  task391->add_dep(task282);
  deciq->add_task(task391);

  vector<IndexRange> I471_index = {active_, active_};
  auto I471 = make_shared<Tensor>(I471_index);
  auto tensor392 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I471};
  auto task392 = make_shared<Task392>(tensor392, pindex);
  task392->add_dep(task282);
  deciq->add_task(task392);

  vector<IndexRange> I472_index = {closed_, virt_};
  auto I472 = make_shared<Tensor>(I472_index);
  auto tensor393 = vector<shared_ptr<Tensor>>{I471, l2, I472};
  auto task393 = make_shared<Task393>(tensor393, pindex);
  task392->add_dep(task393);
  task393->add_dep(task282);
  deciq->add_task(task393);

  vector<IndexRange> I473_index = {closed_, virt_, closed_, virt_};
  auto I473 = make_shared<Tensor>(I473_index);
  auto tensor394 = vector<shared_ptr<Tensor>>{I472, f1_, I473};
  auto task394 = make_shared<Task394>(tensor394, pindex);
  task393->add_dep(task394);
  task394->add_dep(task282);
  deciq->add_task(task394);

  auto tensor395 = vector<shared_ptr<Tensor>>{I473, t2};
  auto task395 = make_shared<Task395>(tensor395, pindex);
  task394->add_dep(task395);
  task395->add_dep(task282);
  deciq->add_task(task395);

  vector<IndexRange> I526_index = {closed_, virt_};
  auto I526 = make_shared<Tensor>(I526_index);
  auto tensor396 = vector<shared_ptr<Tensor>>{I471, l2, I526};
  auto task396 = make_shared<Task396>(tensor396, pindex);
  task392->add_dep(task396);
  task396->add_dep(task282);
  deciq->add_task(task396);

  vector<IndexRange> I527_index = {closed_, virt_, closed_, virt_};
  auto I527 = make_shared<Tensor>(I527_index);
  auto tensor397 = vector<shared_ptr<Tensor>>{I526, f1_, I527};
  auto task397 = make_shared<Task397>(tensor397, pindex);
  task396->add_dep(task397);
  task397->add_dep(task282);
  deciq->add_task(task397);

  auto tensor398 = vector<shared_ptr<Tensor>>{I527, t2};
  auto task398 = make_shared<Task398>(tensor398, pindex);
  task397->add_dep(task398);
  task398->add_dep(task282);
  deciq->add_task(task398);

  vector<IndexRange> I577_index = {virt_, closed_};
  auto I577 = make_shared<Tensor>(I577_index);
  auto tensor399 = vector<shared_ptr<Tensor>>{I471, t2, I577};
  auto task399 = make_shared<Task399>(tensor399, pindex);
  task392->add_dep(task399);
  task399->add_dep(task282);
  deciq->add_task(task399);

  auto tensor400 = vector<shared_ptr<Tensor>>{I577, f1_, l2};
  auto task400 = make_shared<Task400>(tensor400, pindex);
  task399->add_dep(task400);
  task400->add_dep(task282);
  deciq->add_task(task400);

  vector<IndexRange> I581_index = {virt_, closed_};
  auto I581 = make_shared<Tensor>(I581_index);
  auto tensor401 = vector<shared_ptr<Tensor>>{I471, t2, I581};
  auto task401 = make_shared<Task401>(tensor401, pindex);
  task392->add_dep(task401);
  task401->add_dep(task282);
  deciq->add_task(task401);

  auto tensor402 = vector<shared_ptr<Tensor>>{I581, f1_, l2};
  auto task402 = make_shared<Task402>(tensor402, pindex);
  task401->add_dep(task402);
  task402->add_dep(task282);
  deciq->add_task(task402);

  vector<IndexRange> I585_index = {virt_, closed_};
  auto I585 = make_shared<Tensor>(I585_index);
  auto tensor403 = vector<shared_ptr<Tensor>>{I471, t2, I585};
  auto task403 = make_shared<Task403>(tensor403, pindex);
  task392->add_dep(task403);
  task403->add_dep(task282);
  deciq->add_task(task403);

  auto tensor404 = vector<shared_ptr<Tensor>>{I585, f1_, l2};
  auto task404 = make_shared<Task404>(tensor404, pindex);
  task403->add_dep(task404);
  task404->add_dep(task282);
  deciq->add_task(task404);

  vector<IndexRange> I589_index = {virt_, closed_};
  auto I589 = make_shared<Tensor>(I589_index);
  auto tensor405 = vector<shared_ptr<Tensor>>{I471, t2, I589};
  auto task405 = make_shared<Task405>(tensor405, pindex);
  task392->add_dep(task405);
  task405->add_dep(task282);
  deciq->add_task(task405);

  auto tensor406 = vector<shared_ptr<Tensor>>{I589, f1_, l2};
  auto task406 = make_shared<Task406>(tensor406, pindex);
  task405->add_dep(task406);
  task406->add_dep(task282);
  deciq->add_task(task406);

  vector<IndexRange> I615_index = {closed_, active_};
  auto I615 = make_shared<Tensor>(I615_index);
  auto tensor407 = vector<shared_ptr<Tensor>>{I471, f1_, I615};
  auto task407 = make_shared<Task407>(tensor407, pindex);
  task392->add_dep(task407);
  task407->add_dep(task282);
  deciq->add_task(task407);

  auto tensor408 = vector<shared_ptr<Tensor>>{I615, t2, l2};
  auto task408 = make_shared<Task408>(tensor408, pindex);
  task407->add_dep(task408);
  task408->add_dep(task282);
  deciq->add_task(task408);

  auto tensor409 = vector<shared_ptr<Tensor>>{I615, t2, l2};
  auto task409 = make_shared<Task409>(tensor409, pindex);
  task407->add_dep(task409);
  task409->add_dep(task282);
  deciq->add_task(task409);

  vector<IndexRange> I647_index = {active_, closed_};
  auto I647 = make_shared<Tensor>(I647_index);
  auto tensor410 = vector<shared_ptr<Tensor>>{I471, f1_, I647};
  auto task410 = make_shared<Task410>(tensor410, pindex);
  task392->add_dep(task410);
  task410->add_dep(task282);
  deciq->add_task(task410);

  auto tensor411 = vector<shared_ptr<Tensor>>{I647, t2, l2};
  auto task411 = make_shared<Task411>(tensor411, pindex);
  task410->add_dep(task411);
  task411->add_dep(task282);
  deciq->add_task(task411);

  auto tensor412 = vector<shared_ptr<Tensor>>{I647, t2, l2};
  auto task412 = make_shared<Task412>(tensor412, pindex);
  task410->add_dep(task412);
  task412->add_dep(task282);
  deciq->add_task(task412);

  vector<IndexRange> I661_index = {active_, virt_, virt_, closed_};
  auto I661 = make_shared<Tensor>(I661_index);
  auto tensor413 = vector<shared_ptr<Tensor>>{I471, l2, I661};
  auto task413 = make_shared<Task413>(tensor413, pindex);
  task392->add_dep(task413);
  task413->add_dep(task282);
  deciq->add_task(task413);

  vector<IndexRange> I662_index = {active_, virt_, closed_, virt_};
  auto I662 = make_shared<Tensor>(I662_index);
  auto tensor414 = vector<shared_ptr<Tensor>>{I661, f1_, I662};
  auto task414 = make_shared<Task414>(tensor414, pindex);
  task413->add_dep(task414);
  task414->add_dep(task282);
  deciq->add_task(task414);

  auto tensor415 = vector<shared_ptr<Tensor>>{I662, t2};
  auto task415 = make_shared<Task415>(tensor415, pindex);
  task414->add_dep(task415);
  task415->add_dep(task282);
  deciq->add_task(task415);

  vector<IndexRange> I670_index = {active_, virt_, closed_, virt_};
  auto I670 = make_shared<Tensor>(I670_index);
  auto tensor416 = vector<shared_ptr<Tensor>>{I661, f1_, I670};
  auto task416 = make_shared<Task416>(tensor416, pindex);
  task413->add_dep(task416);
  task416->add_dep(task282);
  deciq->add_task(task416);

  auto tensor417 = vector<shared_ptr<Tensor>>{I670, t2};
  auto task417 = make_shared<Task417>(tensor417, pindex);
  task416->add_dep(task417);
  task417->add_dep(task282);
  deciq->add_task(task417);

  vector<IndexRange> I678_index = {active_, virt_, closed_, virt_};
  auto I678 = make_shared<Tensor>(I678_index);
  auto tensor418 = vector<shared_ptr<Tensor>>{I661, f1_, I678};
  auto task418 = make_shared<Task418>(tensor418, pindex);
  task413->add_dep(task418);
  task418->add_dep(task282);
  deciq->add_task(task418);

  auto tensor419 = vector<shared_ptr<Tensor>>{I678, t2};
  auto task419 = make_shared<Task419>(tensor419, pindex);
  task418->add_dep(task419);
  task419->add_dep(task282);
  deciq->add_task(task419);

  vector<IndexRange> I521_index = {active_, active_, active_, active_, active_, active_};
  auto I521 = make_shared<Tensor>(I521_index);
  auto tensor420 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I521};
  auto task420 = make_shared<Task420>(tensor420, pindex);
  task420->add_dep(task282);
  deciq->add_task(task420);

  vector<IndexRange> I522_index = {active_, active_, virt_, active_};
  auto I522 = make_shared<Tensor>(I522_index);
  auto tensor421 = vector<shared_ptr<Tensor>>{I521, t2, I522};
  auto task421 = make_shared<Task421>(tensor421, pindex);
  task420->add_dep(task421);
  task421->add_dep(task282);
  deciq->add_task(task421);

  auto tensor422 = vector<shared_ptr<Tensor>>{I522, f1_, l2};
  auto task422 = make_shared<Task422>(tensor422, pindex);
  task421->add_dep(task422);
  task422->add_dep(task282);
  deciq->add_task(task422);

  vector<IndexRange> I541_index = {active_, active_, active_, active_, active_, active_};
  auto I541 = make_shared<Tensor>(I541_index);
  auto tensor423 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I541};
  auto task423 = make_shared<Task423>(tensor423, pindex);
  task423->add_dep(task282);
  deciq->add_task(task423);

  vector<IndexRange> I542_index = {active_, virt_, active_, active_};
  auto I542 = make_shared<Tensor>(I542_index);
  auto tensor424 = vector<shared_ptr<Tensor>>{I541, l2, I542};
  auto task424 = make_shared<Task424>(tensor424, pindex);
  task423->add_dep(task424);
  task424->add_dep(task282);
  deciq->add_task(task424);

  auto tensor425 = vector<shared_ptr<Tensor>>{I542, f1_, t2};
  auto task425 = make_shared<Task425>(tensor425, pindex);
  task424->add_dep(task425);
  task425->add_dep(task282);
  deciq->add_task(task425);

  vector<IndexRange> I545_index = {active_, active_, active_, active_, active_, active_};
  auto I545 = make_shared<Tensor>(I545_index);
  auto tensor426 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I545};
  auto task426 = make_shared<Task426>(tensor426, pindex);
  task426->add_dep(task282);
  deciq->add_task(task426);

  vector<IndexRange> I546_index = {virt_, active_, active_, active_};
  auto I546 = make_shared<Tensor>(I546_index);
  auto tensor427 = vector<shared_ptr<Tensor>>{I545, l2, I546};
  auto task427 = make_shared<Task427>(tensor427, pindex);
  task426->add_dep(task427);
  task427->add_dep(task282);
  deciq->add_task(task427);

  auto tensor428 = vector<shared_ptr<Tensor>>{I546, f1_, t2};
  auto task428 = make_shared<Task428>(tensor428, pindex);
  task427->add_dep(task428);
  task428->add_dep(task282);
  deciq->add_task(task428);

  vector<IndexRange> I549_index = {active_, active_, active_, active_, active_, active_};
  auto I549 = make_shared<Tensor>(I549_index);
  auto tensor429 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I549, f1_};
  auto task429 = make_shared<Task429>(tensor429, pindex);
  task429->add_dep(task282);
  deciq->add_task(task429);

  auto tensor430 = vector<shared_ptr<Tensor>>{I549, t2, l2};
  auto task430 = make_shared<Task430>(tensor430, pindex);
  task429->add_dep(task430);
  task430->add_dep(task282);
  deciq->add_task(task430);

  vector<IndexRange> I552_index = {active_, active_, active_, active_, active_, active_};
  auto I552 = make_shared<Tensor>(I552_index);
  auto tensor431 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I552};
  auto task431 = make_shared<Task431>(tensor431, pindex);
  task431->add_dep(task282);
  deciq->add_task(task431);

  vector<IndexRange> I553_index = {active_, active_, active_, virt_};
  auto I553 = make_shared<Tensor>(I553_index);
  auto tensor432 = vector<shared_ptr<Tensor>>{I552, l2, I553};
  auto task432 = make_shared<Task432>(tensor432, pindex);
  task431->add_dep(task432);
  task432->add_dep(task282);
  deciq->add_task(task432);

  auto tensor433 = vector<shared_ptr<Tensor>>{I553, f1_, t2};
  auto task433 = make_shared<Task433>(tensor433, pindex);
  task432->add_dep(task433);
  task433->add_dep(task282);
  deciq->add_task(task433);

  auto tensor434 = vector<shared_ptr<Tensor>>{I553, f1_, t2};
  auto task434 = make_shared<Task434>(tensor434, pindex);
  task432->add_dep(task434);
  task434->add_dep(task282);
  deciq->add_task(task434);

  vector<IndexRange> I689_index = {active_, virt_, active_, active_};
  auto I689 = make_shared<Tensor>(I689_index);
  auto tensor435 = vector<shared_ptr<Tensor>>{I552, t2, I689};
  auto task435 = make_shared<Task435>(tensor435, pindex);
  task431->add_dep(task435);
  task435->add_dep(task282);
  deciq->add_task(task435);

  auto tensor436 = vector<shared_ptr<Tensor>>{I689, f1_, l2};
  auto task436 = make_shared<Task436>(tensor436, pindex);
  task435->add_dep(task436);
  task436->add_dep(task282);
  deciq->add_task(task436);

  vector<IndexRange> I556_index = {active_, active_, active_, active_};
  auto I556 = make_shared<Tensor>(I556_index);
  auto tensor437 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I556};
  auto task437 = make_shared<Task437>(tensor437, pindex);
  task437->add_dep(task282);
  deciq->add_task(task437);

  vector<IndexRange> I557_index = {active_, virt_};
  auto I557 = make_shared<Tensor>(I557_index);
  auto tensor438 = vector<shared_ptr<Tensor>>{I556, l2, I557};
  auto task438 = make_shared<Task438>(tensor438, pindex);
  task437->add_dep(task438);
  task438->add_dep(task282);
  deciq->add_task(task438);

  vector<IndexRange> I558_index = {active_, virt_, closed_, virt_};
  auto I558 = make_shared<Tensor>(I558_index);
  auto tensor439 = vector<shared_ptr<Tensor>>{I557, f1_, I558};
  auto task439 = make_shared<Task439>(tensor439, pindex);
  task438->add_dep(task439);
  task439->add_dep(task282);
  deciq->add_task(task439);

  auto tensor440 = vector<shared_ptr<Tensor>>{I558, t2};
  auto task440 = make_shared<Task440>(tensor440, pindex);
  task439->add_dep(task440);
  task440->add_dep(task282);
  deciq->add_task(task440);

  vector<IndexRange> I639_index = {virt_, active_};
  auto I639 = make_shared<Tensor>(I639_index);
  auto tensor441 = vector<shared_ptr<Tensor>>{I556, t2, I639};
  auto task441 = make_shared<Task441>(tensor441, pindex);
  task437->add_dep(task441);
  task441->add_dep(task282);
  deciq->add_task(task441);

  auto tensor442 = vector<shared_ptr<Tensor>>{I639, f1_, l2};
  auto task442 = make_shared<Task442>(tensor442, pindex);
  task441->add_dep(task442);
  task442->add_dep(task282);
  deciq->add_task(task442);

  vector<IndexRange> I643_index = {virt_, active_};
  auto I643 = make_shared<Tensor>(I643_index);
  auto tensor443 = vector<shared_ptr<Tensor>>{I556, t2, I643};
  auto task443 = make_shared<Task443>(tensor443, pindex);
  task437->add_dep(task443);
  task443->add_dep(task282);
  deciq->add_task(task443);

  auto tensor444 = vector<shared_ptr<Tensor>>{I643, f1_, l2};
  auto task444 = make_shared<Task444>(tensor444, pindex);
  task443->add_dep(task444);
  task444->add_dep(task282);
  deciq->add_task(task444);

  vector<IndexRange> I685_index = {virt_, virt_, active_, active_};
  auto I685 = make_shared<Tensor>(I685_index);
  auto tensor445 = vector<shared_ptr<Tensor>>{I556, t2, I685};
  auto task445 = make_shared<Task445>(tensor445, pindex);
  task437->add_dep(task445);
  task445->add_dep(task282);
  deciq->add_task(task445);

  auto tensor446 = vector<shared_ptr<Tensor>>{I685, f1_, l2};
  auto task446 = make_shared<Task446>(tensor446, pindex);
  task445->add_dep(task446);
  task446->add_dep(task282);
  deciq->add_task(task446);

  vector<IndexRange> I693_index = {active_, virt_, virt_, active_};
  auto I693 = make_shared<Tensor>(I693_index);
  auto tensor447 = vector<shared_ptr<Tensor>>{I556, l2, I693};
  auto task447 = make_shared<Task447>(tensor447, pindex);
  task437->add_dep(task447);
  task447->add_dep(task282);
  deciq->add_task(task447);

  auto tensor448 = vector<shared_ptr<Tensor>>{I693, f1_, t2};
  auto task448 = make_shared<Task448>(tensor448, pindex);
  task447->add_dep(task448);
  task448->add_dep(task282);
  deciq->add_task(task448);

  auto tensor449 = vector<shared_ptr<Tensor>>{I693, f1_, t2};
  auto task449 = make_shared<Task449>(tensor449, pindex);
  task447->add_dep(task449);
  task449->add_dep(task282);
  deciq->add_task(task449);

  vector<IndexRange> I592_index;
  auto I592 = make_shared<Tensor>(I592_index);
  auto tensor450 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I592, f1_};
  auto task450 = make_shared<Task450>(tensor450, pindex);
  task450->add_dep(task282);
  deciq->add_task(task450);

  auto tensor451 = vector<shared_ptr<Tensor>>{I592, t2, l2};
  auto task451 = make_shared<Task451>(tensor451, pindex);
  task450->add_dep(task451);
  task451->add_dep(task282);
  deciq->add_task(task451);

  auto tensor452 = vector<shared_ptr<Tensor>>{I592, t2, l2};
  auto task452 = make_shared<Task452>(tensor452, pindex);
  task450->add_dep(task452);
  task452->add_dep(task282);
  deciq->add_task(task452);

  shared_ptr<Tensor> I598;
  if (diagonal) {
    vector<IndexRange> I598_index;
    I598 = make_shared<Tensor>(I598_index);
  }
  shared_ptr<Task453> task453;
  if (diagonal) {
    auto tensor453 = vector<shared_ptr<Tensor>>{den0ci, I598};
    task453 = make_shared<Task453>(tensor453, pindex);
    task453->add_dep(task282);
    deciq->add_task(task453);
  }

  shared_ptr<Tensor> I599;
  if (diagonal) {
    vector<IndexRange> I599_index = {closed_, closed_};
    I599 = make_shared<Tensor>(I599_index);
  }
  shared_ptr<Task454> task454;
  if (diagonal) {
    auto tensor454 = vector<shared_ptr<Tensor>>{I598, f1_, I599};
    task454 = make_shared<Task454>(tensor454, pindex);
    task453->add_dep(task454);
    task454->add_dep(task282);
    deciq->add_task(task454);
  }

  shared_ptr<Task455> task455;
  if (diagonal) {
    auto tensor455 = vector<shared_ptr<Tensor>>{I599, t2, l2};
    task455 = make_shared<Task455>(tensor455, pindex);
    task454->add_dep(task455);
    task455->add_dep(task282);
    deciq->add_task(task455);
  }

  shared_ptr<Tensor> I602;
  if (diagonal) {
    vector<IndexRange> I602_index;
    I602 = make_shared<Tensor>(I602_index);
  }
  shared_ptr<Task456> task456;
  if (diagonal) {
    auto tensor456 = vector<shared_ptr<Tensor>>{den0ci, I602};
    task456 = make_shared<Task456>(tensor456, pindex);
    task456->add_dep(task282);
    deciq->add_task(task456);
  }

  shared_ptr<Tensor> I603;
  if (diagonal) {
    vector<IndexRange> I603_index = {closed_, closed_};
    I603 = make_shared<Tensor>(I603_index);
  }
  shared_ptr<Task457> task457;
  if (diagonal) {
    auto tensor457 = vector<shared_ptr<Tensor>>{I602, f1_, I603};
    task457 = make_shared<Task457>(tensor457, pindex);
    task456->add_dep(task457);
    task457->add_dep(task282);
    deciq->add_task(task457);
  }

  shared_ptr<Task458> task458;
  if (diagonal) {
    auto tensor458 = vector<shared_ptr<Tensor>>{I603, t2, l2};
    task458 = make_shared<Task458>(tensor458, pindex);
    task457->add_dep(task458);
    task458->add_dep(task282);
    deciq->add_task(task458);
  }

  shared_ptr<Tensor> I606;
  if (diagonal) {
    vector<IndexRange> I606_index;
    I606 = make_shared<Tensor>(I606_index);
  }
  shared_ptr<Task459> task459;
  if (diagonal) {
    auto tensor459 = vector<shared_ptr<Tensor>>{den0ci, I606};
    task459 = make_shared<Task459>(tensor459, pindex);
    task459->add_dep(task282);
    deciq->add_task(task459);
  }

  shared_ptr<Tensor> I607;
  if (diagonal) {
    vector<IndexRange> I607_index = {virt_, virt_};
    I607 = make_shared<Tensor>(I607_index);
  }
  shared_ptr<Task460> task460;
  if (diagonal) {
    auto tensor460 = vector<shared_ptr<Tensor>>{I606, f1_, I607};
    task460 = make_shared<Task460>(tensor460, pindex);
    task459->add_dep(task460);
    task460->add_dep(task282);
    deciq->add_task(task460);
  }

  shared_ptr<Task461> task461;
  if (diagonal) {
    auto tensor461 = vector<shared_ptr<Tensor>>{I607, t2, l2};
    task461 = make_shared<Task461>(tensor461, pindex);
    task460->add_dep(task461);
    task461->add_dep(task282);
    deciq->add_task(task461);
  }

  shared_ptr<Tensor> I610;
  if (diagonal) {
    vector<IndexRange> I610_index;
    I610 = make_shared<Tensor>(I610_index);
  }
  shared_ptr<Task462> task462;
  if (diagonal) {
    auto tensor462 = vector<shared_ptr<Tensor>>{den0ci, I610};
    task462 = make_shared<Task462>(tensor462, pindex);
    task462->add_dep(task282);
    deciq->add_task(task462);
  }

  shared_ptr<Tensor> I611;
  if (diagonal) {
    vector<IndexRange> I611_index = {virt_, virt_};
    I611 = make_shared<Tensor>(I611_index);
  }
  shared_ptr<Task463> task463;
  if (diagonal) {
    auto tensor463 = vector<shared_ptr<Tensor>>{I610, f1_, I611};
    task463 = make_shared<Task463>(tensor463, pindex);
    task462->add_dep(task463);
    task463->add_dep(task282);
    deciq->add_task(task463);
  }

  shared_ptr<Task464> task464;
  if (diagonal) {
    auto tensor464 = vector<shared_ptr<Tensor>>{I611, t2, l2};
    task464 = make_shared<Task464>(tensor464, pindex);
    task463->add_dep(task464);
    task464->add_dep(task282);
    deciq->add_task(task464);
  }

  vector<IndexRange> I654_index = {active_, active_};
  auto I654 = make_shared<Tensor>(I654_index);
  auto tensor465 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I654, f1_};
  auto task465 = make_shared<Task465>(tensor465, pindex);
  task465->add_dep(task282);
  deciq->add_task(task465);

  auto tensor466 = vector<shared_ptr<Tensor>>{I654, t2, l2};
  auto task466 = make_shared<Task466>(tensor466, pindex);
  task465->add_dep(task466);
  task466->add_dep(task282);
  deciq->add_task(task466);

  auto tensor467 = vector<shared_ptr<Tensor>>{I654, t2, l2};
  auto task467 = make_shared<Task467>(tensor467, pindex);
  task465->add_dep(task467);
  task467->add_dep(task282);
  deciq->add_task(task467);

  vector<IndexRange> I696_index = {active_, active_, active_, active_};
  auto I696 = make_shared<Tensor>(I696_index);
  auto tensor468 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I696, f1_};
  auto task468 = make_shared<Task468>(tensor468, pindex);
  task468->add_dep(task282);
  deciq->add_task(task468);

  auto tensor469 = vector<shared_ptr<Tensor>>{I696, t2, l2};
  auto task469 = make_shared<Task469>(tensor469, pindex);
  task468->add_dep(task469);
  task469->add_dep(task282);
  deciq->add_task(task469);

  return deciq;
}


#endif
