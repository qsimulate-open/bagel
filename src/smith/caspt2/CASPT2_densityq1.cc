//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_densityq1.cc
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
#include <src/smith/caspt2/CASPT2_tasks6.h>
#include <src/smith/caspt2/CASPT2_tasks7.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_densityq(const bool reset, const bool diagonal) {
  auto densityq = make_shared<Queue>();
  auto tensor284 = vector<shared_ptr<Tensor>>{den2};
  auto task284 = make_shared<Task284>(tensor284, reset);
  densityq->add_task(task284);

  make_densityq1(densityq, task284, diagonal);
  make_densityq2(densityq, task284, diagonal);
  make_densityq3(densityq, task284, diagonal);

  return densityq;
}


void  CASPT2::CASPT2::make_densityq1(shared_ptr<Queue> densityq, shared_ptr<Task> task284, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I374_index = {active_, active_};
  auto I374 = make_shared<Tensor>(I374_index);
  auto tensor285 = vector<shared_ptr<Tensor>>{den2, I374};
  auto task285 = make_shared<Task285>(tensor285, pindex);
  task285->add_dep(task284);
  densityq->add_task(task285);

  vector<IndexRange> I375_index = {active_, active_, active_, active_};
  auto I375 = make_shared<Tensor>(I375_index);
  auto tensor286 = vector<shared_ptr<Tensor>>{I374, Gamma138_(), I375};
  auto task286 = make_shared<Task286>(tensor286, pindex);
  task285->add_dep(task286);
  task286->add_dep(task284);
  densityq->add_task(task286);

  auto tensor287 = vector<shared_ptr<Tensor>>{I375, t2};
  auto task287 = make_shared<Task287>(tensor287, pindex);
  task286->add_dep(task287);
  task287->add_dep(task284);
  densityq->add_task(task287);

  vector<IndexRange> I468_index = {active_, active_, active_, active_};
  auto I468 = make_shared<Tensor>(I468_index);
  auto tensor288 = vector<shared_ptr<Tensor>>{I374, Gamma169_(), I468};
  auto task288 = make_shared<Task288>(tensor288, pindex);
  task285->add_dep(task288);
  task288->add_dep(task284);
  densityq->add_task(task288);

  auto tensor289 = vector<shared_ptr<Tensor>>{I468, t2};
  auto task289 = make_shared<Task289>(tensor289, pindex);
  task288->add_dep(task289);
  task289->add_dep(task284);
  densityq->add_task(task289);

  vector<IndexRange> I477_index = {active_, active_, active_, active_};
  auto I477 = make_shared<Tensor>(I477_index);
  auto tensor290 = vector<shared_ptr<Tensor>>{I374, Gamma172_(), I477};
  auto task290 = make_shared<Task290>(tensor290, pindex);
  task285->add_dep(task290);
  task290->add_dep(task284);
  densityq->add_task(task290);

  auto tensor291 = vector<shared_ptr<Tensor>>{I477, t2};
  auto task291 = make_shared<Task291>(tensor291, pindex);
  task290->add_dep(task291);
  task291->add_dep(task284);
  densityq->add_task(task291);

  auto tensor292 = vector<shared_ptr<Tensor>>{I477, t2};
  auto task292 = make_shared<Task292>(tensor292, pindex);
  task290->add_dep(task292);
  task292->add_dep(task284);
  densityq->add_task(task292);

  auto tensor293 = vector<shared_ptr<Tensor>>{I477, t2};
  auto task293 = make_shared<Task293>(tensor293, pindex);
  task290->add_dep(task293);
  task293->add_dep(task284);
  densityq->add_task(task293);

  vector<IndexRange> I659_index = {active_, active_, active_, active_};
  auto I659 = make_shared<Tensor>(I659_index);
  auto tensor294 = vector<shared_ptr<Tensor>>{I374, Gamma230_(), I659};
  auto task294 = make_shared<Task294>(tensor294, pindex);
  task285->add_dep(task294);
  task294->add_dep(task284);
  densityq->add_task(task294);

  auto tensor295 = vector<shared_ptr<Tensor>>{I659, t2};
  auto task295 = make_shared<Task295>(tensor295, pindex);
  task294->add_dep(task295);
  task295->add_dep(task284);
  densityq->add_task(task295);

  vector<IndexRange> I377_index = {closed_, closed_};
  auto I377 = make_shared<Tensor>(I377_index);
  auto tensor296 = vector<shared_ptr<Tensor>>{den2, I377};
  auto task296 = make_shared<Task296>(tensor296, pindex);
  task296->add_dep(task284);
  densityq->add_task(task296);

  vector<IndexRange> I378_index = {closed_, closed_, active_, active_};
  auto I378 = make_shared<Tensor>(I378_index);
  auto tensor297 = vector<shared_ptr<Tensor>>{I377, t2, I378};
  auto task297 = make_shared<Task297>(tensor297, pindex);
  task296->add_dep(task297);
  task297->add_dep(task284);
  densityq->add_task(task297);

  auto tensor298 = vector<shared_ptr<Tensor>>{I378, Gamma92_(), t2};
  auto task298 = make_shared<Task298>(tensor298, pindex);
  task297->add_dep(task298);
  task298->add_dep(task284);
  densityq->add_task(task298);

  vector<IndexRange> I471_index = {closed_, virt_, active_, active_};
  auto I471 = make_shared<Tensor>(I471_index);
  auto tensor299 = vector<shared_ptr<Tensor>>{I377, t2, I471};
  auto task299 = make_shared<Task299>(tensor299, pindex);
  task296->add_dep(task299);
  task299->add_dep(task284);
  densityq->add_task(task299);

  auto tensor300 = vector<shared_ptr<Tensor>>{I471, Gamma32_(), t2};
  auto task300 = make_shared<Task300>(tensor300, pindex);
  task299->add_dep(task300);
  task300->add_dep(task284);
  densityq->add_task(task300);

  vector<IndexRange> I480_index = {closed_, virt_, active_, active_};
  auto I480 = make_shared<Tensor>(I480_index);
  auto tensor301 = vector<shared_ptr<Tensor>>{I377, t2, I480};
  auto task301 = make_shared<Task301>(tensor301, pindex);
  task296->add_dep(task301);
  task301->add_dep(task284);
  densityq->add_task(task301);

  auto tensor302 = vector<shared_ptr<Tensor>>{I480, Gamma35_(), t2};
  auto task302 = make_shared<Task302>(tensor302, pindex);
  task301->add_dep(task302);
  task302->add_dep(task284);
  densityq->add_task(task302);

  vector<IndexRange> I380_index = {active_, closed_};
  auto I380 = make_shared<Tensor>(I380_index);
  auto tensor303 = vector<shared_ptr<Tensor>>{den2, I380};
  auto task303 = make_shared<Task303>(tensor303, pindex);
  task303->add_dep(task284);
  densityq->add_task(task303);

  vector<IndexRange> I381_index = {closed_, active_, active_, active_};
  auto I381 = make_shared<Tensor>(I381_index);
  auto tensor304 = vector<shared_ptr<Tensor>>{I380, t2, I381};
  auto task304 = make_shared<Task304>(tensor304, pindex);
  task303->add_dep(task304);
  task304->add_dep(task284);
  densityq->add_task(task304);

  auto tensor305 = vector<shared_ptr<Tensor>>{I381, Gamma2_(), t2};
  auto task305 = make_shared<Task305>(tensor305, pindex);
  task304->add_dep(task305);
  task305->add_dep(task284);
  densityq->add_task(task305);

  vector<IndexRange> I486_index = {virt_, active_, active_, active_};
  auto I486 = make_shared<Tensor>(I486_index);
  auto tensor306 = vector<shared_ptr<Tensor>>{I380, t2, I486};
  auto task306 = make_shared<Task306>(tensor306, pindex);
  task303->add_dep(task306);
  task306->add_dep(task284);
  densityq->add_task(task306);

  auto tensor307 = vector<shared_ptr<Tensor>>{I486, Gamma37_(), t2};
  auto task307 = make_shared<Task307>(tensor307, pindex);
  task306->add_dep(task307);
  task307->add_dep(task284);
  densityq->add_task(task307);

  vector<IndexRange> I383_index = {active_, virt_};
  auto I383 = make_shared<Tensor>(I383_index);
  auto tensor308 = vector<shared_ptr<Tensor>>{den2, I383};
  auto task308 = make_shared<Task308>(tensor308, pindex);
  task308->add_dep(task284);
  densityq->add_task(task308);

  vector<IndexRange> I384_index = {closed_, closed_, active_, active_};
  auto I384 = make_shared<Tensor>(I384_index);
  auto tensor309 = vector<shared_ptr<Tensor>>{I383, t2, I384};
  auto task309 = make_shared<Task309>(tensor309, pindex);
  task308->add_dep(task309);
  task309->add_dep(task284);
  densityq->add_task(task309);

  auto tensor310 = vector<shared_ptr<Tensor>>{I384, Gamma3_(), t2};
  auto task310 = make_shared<Task310>(tensor310, pindex);
  task309->add_dep(task310);
  task310->add_dep(task284);
  densityq->add_task(task310);

  vector<IndexRange> I495_index = {closed_, virt_, active_, active_};
  auto I495 = make_shared<Tensor>(I495_index);
  auto tensor311 = vector<shared_ptr<Tensor>>{I383, t2, I495};
  auto task311 = make_shared<Task311>(tensor311, pindex);
  task308->add_dep(task311);
  task311->add_dep(task284);
  densityq->add_task(task311);

  auto tensor312 = vector<shared_ptr<Tensor>>{I495, Gamma35_(), t2};
  auto task312 = make_shared<Task312>(tensor312, pindex);
  task311->add_dep(task312);
  task312->add_dep(task284);
  densityq->add_task(task312);

  vector<IndexRange> I498_index = {closed_, virt_, active_, active_};
  auto I498 = make_shared<Tensor>(I498_index);
  auto tensor313 = vector<shared_ptr<Tensor>>{I383, t2, I498};
  auto task313 = make_shared<Task313>(tensor313, pindex);
  task308->add_dep(task313);
  task313->add_dep(task284);
  densityq->add_task(task313);

  auto tensor314 = vector<shared_ptr<Tensor>>{I498, Gamma32_(), t2};
  auto task314 = make_shared<Task314>(tensor314, pindex);
  task313->add_dep(task314);
  task314->add_dep(task284);
  densityq->add_task(task314);

  vector<IndexRange> I537_index = {virt_, closed_, active_, active_};
  auto I537 = make_shared<Tensor>(I537_index);
  auto tensor315 = vector<shared_ptr<Tensor>>{I383, t2, I537};
  auto task315 = make_shared<Task315>(tensor315, pindex);
  task308->add_dep(task315);
  task315->add_dep(task284);
  densityq->add_task(task315);

  auto tensor316 = vector<shared_ptr<Tensor>>{I537, Gamma35_(), t2};
  auto task316 = make_shared<Task316>(tensor316, pindex);
  task315->add_dep(task316);
  task316->add_dep(task284);
  densityq->add_task(task316);

  vector<IndexRange> I540_index = {virt_, closed_, active_, active_};
  auto I540 = make_shared<Tensor>(I540_index);
  auto tensor317 = vector<shared_ptr<Tensor>>{I383, t2, I540};
  auto task317 = make_shared<Task317>(tensor317, pindex);
  task308->add_dep(task317);
  task317->add_dep(task284);
  densityq->add_task(task317);

  auto tensor318 = vector<shared_ptr<Tensor>>{I540, Gamma35_(), t2};
  auto task318 = make_shared<Task318>(tensor318, pindex);
  task317->add_dep(task318);
  task318->add_dep(task284);
  densityq->add_task(task318);

  vector<IndexRange> I386_index = {active_, closed_};
  auto I386 = make_shared<Tensor>(I386_index);
  auto tensor319 = vector<shared_ptr<Tensor>>{den2, I386};
  auto task319 = make_shared<Task319>(tensor319, pindex);
  task319->add_dep(task284);
  densityq->add_task(task319);

  vector<IndexRange> I387_index = {closed_, active_, active_, active_};
  auto I387 = make_shared<Tensor>(I387_index);
  auto tensor320 = vector<shared_ptr<Tensor>>{I386, t2, I387};
  auto task320 = make_shared<Task320>(tensor320, pindex);
  task319->add_dep(task320);
  task320->add_dep(task284);
  densityq->add_task(task320);

  auto tensor321 = vector<shared_ptr<Tensor>>{I387, Gamma4_(), t2};
  auto task321 = make_shared<Task321>(tensor321, pindex);
  task320->add_dep(task321);
  task321->add_dep(task284);
  densityq->add_task(task321);

  vector<IndexRange> I543_index = {virt_, active_, active_, active_};
  auto I543 = make_shared<Tensor>(I543_index);
  auto tensor322 = vector<shared_ptr<Tensor>>{I386, t2, I543};
  auto task322 = make_shared<Task322>(tensor322, pindex);
  task319->add_dep(task322);
  task322->add_dep(task284);
  densityq->add_task(task322);

  auto tensor323 = vector<shared_ptr<Tensor>>{I543, Gamma56_(), t2};
  auto task323 = make_shared<Task323>(tensor323, pindex);
  task322->add_dep(task323);
  task323->add_dep(task284);
  densityq->add_task(task323);

  vector<IndexRange> I546_index = {virt_, active_, active_, active_};
  auto I546 = make_shared<Tensor>(I546_index);
  auto tensor324 = vector<shared_ptr<Tensor>>{I386, t2, I546};
  auto task324 = make_shared<Task324>(tensor324, pindex);
  task319->add_dep(task324);
  task324->add_dep(task284);
  densityq->add_task(task324);

  auto tensor325 = vector<shared_ptr<Tensor>>{I546, Gamma57_(), t2};
  auto task325 = make_shared<Task325>(tensor325, pindex);
  task324->add_dep(task325);
  task325->add_dep(task284);
  densityq->add_task(task325);

  vector<IndexRange> I389_index = {active_, active_};
  auto I389 = make_shared<Tensor>(I389_index);
  auto tensor326 = vector<shared_ptr<Tensor>>{den2, I389};
  auto task326 = make_shared<Task326>(tensor326, pindex);
  task326->add_dep(task284);
  densityq->add_task(task326);

  vector<IndexRange> I390_index = {active_, active_, active_, active_, active_, active_};
  auto I390 = make_shared<Tensor>(I390_index);
  auto tensor327 = vector<shared_ptr<Tensor>>{I389, Gamma143_(), I390};
  auto task327 = make_shared<Task327>(tensor327, pindex);
  task326->add_dep(task327);
  task327->add_dep(task284);
  densityq->add_task(task327);

  auto tensor328 = vector<shared_ptr<Tensor>>{I390, t2};
  auto task328 = make_shared<Task328>(tensor328, pindex);
  task327->add_dep(task328);
  task328->add_dep(task284);
  densityq->add_task(task328);

  vector<IndexRange> I549_index = {active_, active_, active_, active_, active_, active_};
  auto I549 = make_shared<Tensor>(I549_index);
  auto tensor329 = vector<shared_ptr<Tensor>>{I389, Gamma196_(), I549};
  auto task329 = make_shared<Task329>(tensor329, pindex);
  task326->add_dep(task329);
  task329->add_dep(task284);
  densityq->add_task(task329);

  auto tensor330 = vector<shared_ptr<Tensor>>{I549, t2};
  auto task330 = make_shared<Task330>(tensor330, pindex);
  task329->add_dep(task330);
  task330->add_dep(task284);
  densityq->add_task(task330);

  vector<IndexRange> I392_index = {closed_, closed_};
  auto I392 = make_shared<Tensor>(I392_index);
  auto tensor331 = vector<shared_ptr<Tensor>>{den2, I392};
  auto task331 = make_shared<Task331>(tensor331, pindex);
  task331->add_dep(task284);
  densityq->add_task(task331);

  vector<IndexRange> I393_index = {closed_, active_, active_, active_};
  auto I393 = make_shared<Tensor>(I393_index);
  auto tensor332 = vector<shared_ptr<Tensor>>{I392, t2, I393};
  auto task332 = make_shared<Task332>(tensor332, pindex);
  task331->add_dep(task332);
  task332->add_dep(task284);
  densityq->add_task(task332);

  auto tensor333 = vector<shared_ptr<Tensor>>{I393, Gamma6_(), t2};
  auto task333 = make_shared<Task333>(tensor333, pindex);
  task332->add_dep(task333);
  task333->add_dep(task284);
  densityq->add_task(task333);

  vector<IndexRange> I395_index = {closed_, virt_};
  auto I395 = make_shared<Tensor>(I395_index);
  auto tensor334 = vector<shared_ptr<Tensor>>{den2, I395};
  auto task334 = make_shared<Task334>(tensor334, pindex);
  task334->add_dep(task284);
  densityq->add_task(task334);

  vector<IndexRange> I396_index = {closed_, active_};
  auto I396 = make_shared<Tensor>(I396_index);
  auto tensor335 = vector<shared_ptr<Tensor>>{I395, t2, I396};
  auto task335 = make_shared<Task335>(tensor335, pindex);
  task334->add_dep(task335);
  task335->add_dep(task284);
  densityq->add_task(task335);

  auto tensor336 = vector<shared_ptr<Tensor>>{I396, Gamma7_(), t2};
  auto task336 = make_shared<Task336>(tensor336, pindex);
  task335->add_dep(task336);
  task336->add_dep(task284);
  densityq->add_task(task336);

  vector<IndexRange> I399_index = {closed_, active_};
  auto I399 = make_shared<Tensor>(I399_index);
  auto tensor337 = vector<shared_ptr<Tensor>>{I395, t2, I399};
  auto task337 = make_shared<Task337>(tensor337, pindex);
  task334->add_dep(task337);
  task337->add_dep(task284);
  densityq->add_task(task337);

  auto tensor338 = vector<shared_ptr<Tensor>>{I399, Gamma7_(), t2};
  auto task338 = make_shared<Task338>(tensor338, pindex);
  task337->add_dep(task338);
  task338->add_dep(task284);
  densityq->add_task(task338);

  vector<IndexRange> I555_index = {virt_, active_};
  auto I555 = make_shared<Tensor>(I555_index);
  auto tensor339 = vector<shared_ptr<Tensor>>{I395, t2, I555};
  auto task339 = make_shared<Task339>(tensor339, pindex);
  task334->add_dep(task339);
  task339->add_dep(task284);
  densityq->add_task(task339);

  auto tensor340 = vector<shared_ptr<Tensor>>{I555, Gamma60_(), t2};
  auto task340 = make_shared<Task340>(tensor340, pindex);
  task339->add_dep(task340);
  task340->add_dep(task284);
  densityq->add_task(task340);

  vector<IndexRange> I558_index = {virt_, active_};
  auto I558 = make_shared<Tensor>(I558_index);
  auto tensor341 = vector<shared_ptr<Tensor>>{I395, t2, I558};
  auto task341 = make_shared<Task341>(tensor341, pindex);
  task334->add_dep(task341);
  task341->add_dep(task284);
  densityq->add_task(task341);

  auto tensor342 = vector<shared_ptr<Tensor>>{I558, Gamma60_(), t2};
  auto task342 = make_shared<Task342>(tensor342, pindex);
  task341->add_dep(task342);
  task342->add_dep(task284);
  densityq->add_task(task342);
}

#endif
