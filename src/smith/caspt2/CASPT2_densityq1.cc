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
  auto tensor285 = vector<shared_ptr<Tensor>>{den2};
  auto task285 = make_shared<Task285>(tensor285, reset);
  densityq->add_task(task285);

  make_densityq1(densityq, task285, diagonal);
  make_densityq2(densityq, task285, diagonal);
  make_densityq3(densityq, task285, diagonal);

  return densityq;
}

void CASPT2::CASPT2::make_densityq1(shared_ptr<Queue> densityq, shared_ptr<Task> task285, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I360_index = {active_, active_};
  auto I360 = make_shared<Tensor>(I360_index);
  auto tensor286 = vector<shared_ptr<Tensor>>{den2, I360};
  auto task286 = make_shared<Task286>(tensor286, pindex);
  task286->add_dep(task285);
  densityq->add_task(task286);

  vector<IndexRange> I361_index = {active_, active_, active_, active_};
  auto I361 = make_shared<Tensor>(I361_index);
  auto tensor287 = vector<shared_ptr<Tensor>>{I360, Gamma138_(), I361};
  auto task287 = make_shared<Task287>(tensor287, pindex);
  task286->add_dep(task287);
  task287->add_dep(task285);
  densityq->add_task(task287);

  auto tensor288 = vector<shared_ptr<Tensor>>{I361, t2};
  auto task288 = make_shared<Task288>(tensor288, pindex);
  task287->add_dep(task288);
  task288->add_dep(task285);
  densityq->add_task(task288);

  vector<IndexRange> I454_index = {active_, active_, active_, active_};
  auto I454 = make_shared<Tensor>(I454_index);
  auto tensor289 = vector<shared_ptr<Tensor>>{I360, Gamma169_(), I454};
  auto task289 = make_shared<Task289>(tensor289, pindex);
  task286->add_dep(task289);
  task289->add_dep(task285);
  densityq->add_task(task289);

  auto tensor290 = vector<shared_ptr<Tensor>>{I454, t2};
  auto task290 = make_shared<Task290>(tensor290, pindex);
  task289->add_dep(task290);
  task290->add_dep(task285);
  densityq->add_task(task290);

  vector<IndexRange> I463_index = {active_, active_, active_, active_};
  auto I463 = make_shared<Tensor>(I463_index);
  auto tensor291 = vector<shared_ptr<Tensor>>{I360, Gamma172_(), I463};
  auto task291 = make_shared<Task291>(tensor291, pindex);
  task286->add_dep(task291);
  task291->add_dep(task285);
  densityq->add_task(task291);

  auto tensor292 = vector<shared_ptr<Tensor>>{I463, t2};
  auto task292 = make_shared<Task292>(tensor292, pindex);
  task291->add_dep(task292);
  task292->add_dep(task285);
  densityq->add_task(task292);

  auto tensor293 = vector<shared_ptr<Tensor>>{I463, t2};
  auto task293 = make_shared<Task293>(tensor293, pindex);
  task291->add_dep(task293);
  task293->add_dep(task285);
  densityq->add_task(task293);

  auto tensor294 = vector<shared_ptr<Tensor>>{I463, t2};
  auto task294 = make_shared<Task294>(tensor294, pindex);
  task291->add_dep(task294);
  task294->add_dep(task285);
  densityq->add_task(task294);

  vector<IndexRange> I645_index = {active_, active_, active_, active_};
  auto I645 = make_shared<Tensor>(I645_index);
  auto tensor295 = vector<shared_ptr<Tensor>>{I360, Gamma230_(), I645};
  auto task295 = make_shared<Task295>(tensor295, pindex);
  task286->add_dep(task295);
  task295->add_dep(task285);
  densityq->add_task(task295);

  auto tensor296 = vector<shared_ptr<Tensor>>{I645, t2};
  auto task296 = make_shared<Task296>(tensor296, pindex);
  task295->add_dep(task296);
  task296->add_dep(task285);
  densityq->add_task(task296);

  vector<IndexRange> I363_index = {closed_, closed_};
  auto I363 = make_shared<Tensor>(I363_index);
  auto tensor297 = vector<shared_ptr<Tensor>>{den2, I363};
  auto task297 = make_shared<Task297>(tensor297, pindex);
  task297->add_dep(task285);
  densityq->add_task(task297);

  vector<IndexRange> I364_index = {closed_, closed_, active_, active_};
  auto I364 = make_shared<Tensor>(I364_index);
  auto tensor298 = vector<shared_ptr<Tensor>>{I363, t2, I364};
  auto task298 = make_shared<Task298>(tensor298, pindex);
  task297->add_dep(task298);
  task298->add_dep(task285);
  densityq->add_task(task298);

  auto tensor299 = vector<shared_ptr<Tensor>>{I364, Gamma92_(), t2};
  auto task299 = make_shared<Task299>(tensor299, pindex);
  task298->add_dep(task299);
  task299->add_dep(task285);
  densityq->add_task(task299);

  vector<IndexRange> I457_index = {closed_, virt_, active_, active_};
  auto I457 = make_shared<Tensor>(I457_index);
  auto tensor300 = vector<shared_ptr<Tensor>>{I363, t2, I457};
  auto task300 = make_shared<Task300>(tensor300, pindex);
  task297->add_dep(task300);
  task300->add_dep(task285);
  densityq->add_task(task300);

  auto tensor301 = vector<shared_ptr<Tensor>>{I457, Gamma32_(), t2};
  auto task301 = make_shared<Task301>(tensor301, pindex);
  task300->add_dep(task301);
  task301->add_dep(task285);
  densityq->add_task(task301);

  vector<IndexRange> I466_index = {closed_, virt_, active_, active_};
  auto I466 = make_shared<Tensor>(I466_index);
  auto tensor302 = vector<shared_ptr<Tensor>>{I363, t2, I466};
  auto task302 = make_shared<Task302>(tensor302, pindex);
  task297->add_dep(task302);
  task302->add_dep(task285);
  densityq->add_task(task302);

  auto tensor303 = vector<shared_ptr<Tensor>>{I466, Gamma35_(), t2};
  auto task303 = make_shared<Task303>(tensor303, pindex);
  task302->add_dep(task303);
  task303->add_dep(task285);
  densityq->add_task(task303);

  vector<IndexRange> I366_index = {active_, closed_};
  auto I366 = make_shared<Tensor>(I366_index);
  auto tensor304 = vector<shared_ptr<Tensor>>{den2, I366};
  auto task304 = make_shared<Task304>(tensor304, pindex);
  task304->add_dep(task285);
  densityq->add_task(task304);

  vector<IndexRange> I367_index = {closed_, active_, active_, active_};
  auto I367 = make_shared<Tensor>(I367_index);
  auto tensor305 = vector<shared_ptr<Tensor>>{I366, t2, I367};
  auto task305 = make_shared<Task305>(tensor305, pindex);
  task304->add_dep(task305);
  task305->add_dep(task285);
  densityq->add_task(task305);

  auto tensor306 = vector<shared_ptr<Tensor>>{I367, Gamma2_(), t2};
  auto task306 = make_shared<Task306>(tensor306, pindex);
  task305->add_dep(task306);
  task306->add_dep(task285);
  densityq->add_task(task306);

  vector<IndexRange> I472_index = {virt_, active_, active_, active_};
  auto I472 = make_shared<Tensor>(I472_index);
  auto tensor307 = vector<shared_ptr<Tensor>>{I366, t2, I472};
  auto task307 = make_shared<Task307>(tensor307, pindex);
  task304->add_dep(task307);
  task307->add_dep(task285);
  densityq->add_task(task307);

  auto tensor308 = vector<shared_ptr<Tensor>>{I472, Gamma37_(), t2};
  auto task308 = make_shared<Task308>(tensor308, pindex);
  task307->add_dep(task308);
  task308->add_dep(task285);
  densityq->add_task(task308);

  vector<IndexRange> I369_index = {active_, virt_};
  auto I369 = make_shared<Tensor>(I369_index);
  auto tensor309 = vector<shared_ptr<Tensor>>{den2, I369};
  auto task309 = make_shared<Task309>(tensor309, pindex);
  task309->add_dep(task285);
  densityq->add_task(task309);

  vector<IndexRange> I370_index = {closed_, closed_, active_, active_};
  auto I370 = make_shared<Tensor>(I370_index);
  auto tensor310 = vector<shared_ptr<Tensor>>{I369, t2, I370};
  auto task310 = make_shared<Task310>(tensor310, pindex);
  task309->add_dep(task310);
  task310->add_dep(task285);
  densityq->add_task(task310);

  auto tensor311 = vector<shared_ptr<Tensor>>{I370, Gamma3_(), t2};
  auto task311 = make_shared<Task311>(tensor311, pindex);
  task310->add_dep(task311);
  task311->add_dep(task285);
  densityq->add_task(task311);

  vector<IndexRange> I481_index = {closed_, virt_, active_, active_};
  auto I481 = make_shared<Tensor>(I481_index);
  auto tensor312 = vector<shared_ptr<Tensor>>{I369, t2, I481};
  auto task312 = make_shared<Task312>(tensor312, pindex);
  task309->add_dep(task312);
  task312->add_dep(task285);
  densityq->add_task(task312);

  auto tensor313 = vector<shared_ptr<Tensor>>{I481, Gamma35_(), t2};
  auto task313 = make_shared<Task313>(tensor313, pindex);
  task312->add_dep(task313);
  task313->add_dep(task285);
  densityq->add_task(task313);

  vector<IndexRange> I484_index = {closed_, virt_, active_, active_};
  auto I484 = make_shared<Tensor>(I484_index);
  auto tensor314 = vector<shared_ptr<Tensor>>{I369, t2, I484};
  auto task314 = make_shared<Task314>(tensor314, pindex);
  task309->add_dep(task314);
  task314->add_dep(task285);
  densityq->add_task(task314);

  auto tensor315 = vector<shared_ptr<Tensor>>{I484, Gamma32_(), t2};
  auto task315 = make_shared<Task315>(tensor315, pindex);
  task314->add_dep(task315);
  task315->add_dep(task285);
  densityq->add_task(task315);

  vector<IndexRange> I523_index = {virt_, closed_, active_, active_};
  auto I523 = make_shared<Tensor>(I523_index);
  auto tensor316 = vector<shared_ptr<Tensor>>{I369, t2, I523};
  auto task316 = make_shared<Task316>(tensor316, pindex);
  task309->add_dep(task316);
  task316->add_dep(task285);
  densityq->add_task(task316);

  auto tensor317 = vector<shared_ptr<Tensor>>{I523, Gamma35_(), t2};
  auto task317 = make_shared<Task317>(tensor317, pindex);
  task316->add_dep(task317);
  task317->add_dep(task285);
  densityq->add_task(task317);

  vector<IndexRange> I526_index = {virt_, closed_, active_, active_};
  auto I526 = make_shared<Tensor>(I526_index);
  auto tensor318 = vector<shared_ptr<Tensor>>{I369, t2, I526};
  auto task318 = make_shared<Task318>(tensor318, pindex);
  task309->add_dep(task318);
  task318->add_dep(task285);
  densityq->add_task(task318);

  auto tensor319 = vector<shared_ptr<Tensor>>{I526, Gamma35_(), t2};
  auto task319 = make_shared<Task319>(tensor319, pindex);
  task318->add_dep(task319);
  task319->add_dep(task285);
  densityq->add_task(task319);

  vector<IndexRange> I372_index = {active_, closed_};
  auto I372 = make_shared<Tensor>(I372_index);
  auto tensor320 = vector<shared_ptr<Tensor>>{den2, I372};
  auto task320 = make_shared<Task320>(tensor320, pindex);
  task320->add_dep(task285);
  densityq->add_task(task320);

  vector<IndexRange> I373_index = {closed_, active_, active_, active_};
  auto I373 = make_shared<Tensor>(I373_index);
  auto tensor321 = vector<shared_ptr<Tensor>>{I372, t2, I373};
  auto task321 = make_shared<Task321>(tensor321, pindex);
  task320->add_dep(task321);
  task321->add_dep(task285);
  densityq->add_task(task321);

  auto tensor322 = vector<shared_ptr<Tensor>>{I373, Gamma4_(), t2};
  auto task322 = make_shared<Task322>(tensor322, pindex);
  task321->add_dep(task322);
  task322->add_dep(task285);
  densityq->add_task(task322);

  vector<IndexRange> I529_index = {virt_, active_, active_, active_};
  auto I529 = make_shared<Tensor>(I529_index);
  auto tensor323 = vector<shared_ptr<Tensor>>{I372, t2, I529};
  auto task323 = make_shared<Task323>(tensor323, pindex);
  task320->add_dep(task323);
  task323->add_dep(task285);
  densityq->add_task(task323);

  auto tensor324 = vector<shared_ptr<Tensor>>{I529, Gamma56_(), t2};
  auto task324 = make_shared<Task324>(tensor324, pindex);
  task323->add_dep(task324);
  task324->add_dep(task285);
  densityq->add_task(task324);

  vector<IndexRange> I532_index = {virt_, active_, active_, active_};
  auto I532 = make_shared<Tensor>(I532_index);
  auto tensor325 = vector<shared_ptr<Tensor>>{I372, t2, I532};
  auto task325 = make_shared<Task325>(tensor325, pindex);
  task320->add_dep(task325);
  task325->add_dep(task285);
  densityq->add_task(task325);

  auto tensor326 = vector<shared_ptr<Tensor>>{I532, Gamma57_(), t2};
  auto task326 = make_shared<Task326>(tensor326, pindex);
  task325->add_dep(task326);
  task326->add_dep(task285);
  densityq->add_task(task326);

  vector<IndexRange> I375_index = {active_, active_};
  auto I375 = make_shared<Tensor>(I375_index);
  auto tensor327 = vector<shared_ptr<Tensor>>{den2, I375};
  auto task327 = make_shared<Task327>(tensor327, pindex);
  task327->add_dep(task285);
  densityq->add_task(task327);

  vector<IndexRange> I376_index = {active_, active_, active_, active_, active_, active_};
  auto I376 = make_shared<Tensor>(I376_index);
  auto tensor328 = vector<shared_ptr<Tensor>>{I375, Gamma143_(), I376};
  auto task328 = make_shared<Task328>(tensor328, pindex);
  task327->add_dep(task328);
  task328->add_dep(task285);
  densityq->add_task(task328);

  auto tensor329 = vector<shared_ptr<Tensor>>{I376, t2};
  auto task329 = make_shared<Task329>(tensor329, pindex);
  task328->add_dep(task329);
  task329->add_dep(task285);
  densityq->add_task(task329);

  vector<IndexRange> I535_index = {active_, active_, active_, active_, active_, active_};
  auto I535 = make_shared<Tensor>(I535_index);
  auto tensor330 = vector<shared_ptr<Tensor>>{I375, Gamma196_(), I535};
  auto task330 = make_shared<Task330>(tensor330, pindex);
  task327->add_dep(task330);
  task330->add_dep(task285);
  densityq->add_task(task330);

  auto tensor331 = vector<shared_ptr<Tensor>>{I535, t2};
  auto task331 = make_shared<Task331>(tensor331, pindex);
  task330->add_dep(task331);
  task331->add_dep(task285);
  densityq->add_task(task331);

  vector<IndexRange> I378_index = {closed_, closed_};
  auto I378 = make_shared<Tensor>(I378_index);
  auto tensor332 = vector<shared_ptr<Tensor>>{den2, I378};
  auto task332 = make_shared<Task332>(tensor332, pindex);
  task332->add_dep(task285);
  densityq->add_task(task332);

  vector<IndexRange> I379_index = {closed_, active_, active_, active_};
  auto I379 = make_shared<Tensor>(I379_index);
  auto tensor333 = vector<shared_ptr<Tensor>>{I378, t2, I379};
  auto task333 = make_shared<Task333>(tensor333, pindex);
  task332->add_dep(task333);
  task333->add_dep(task285);
  densityq->add_task(task333);

  auto tensor334 = vector<shared_ptr<Tensor>>{I379, Gamma6_(), t2};
  auto task334 = make_shared<Task334>(tensor334, pindex);
  task333->add_dep(task334);
  task334->add_dep(task285);
  densityq->add_task(task334);

  vector<IndexRange> I381_index = {closed_, virt_};
  auto I381 = make_shared<Tensor>(I381_index);
  auto tensor335 = vector<shared_ptr<Tensor>>{den2, I381};
  auto task335 = make_shared<Task335>(tensor335, pindex);
  task335->add_dep(task285);
  densityq->add_task(task335);

  vector<IndexRange> I382_index = {closed_, active_};
  auto I382 = make_shared<Tensor>(I382_index);
  auto tensor336 = vector<shared_ptr<Tensor>>{I381, t2, I382};
  auto task336 = make_shared<Task336>(tensor336, pindex);
  task335->add_dep(task336);
  task336->add_dep(task285);
  densityq->add_task(task336);

  auto tensor337 = vector<shared_ptr<Tensor>>{I382, Gamma7_(), t2};
  auto task337 = make_shared<Task337>(tensor337, pindex);
  task336->add_dep(task337);
  task337->add_dep(task285);
  densityq->add_task(task337);

  vector<IndexRange> I385_index = {closed_, active_};
  auto I385 = make_shared<Tensor>(I385_index);
  auto tensor338 = vector<shared_ptr<Tensor>>{I381, t2, I385};
  auto task338 = make_shared<Task338>(tensor338, pindex);
  task335->add_dep(task338);
  task338->add_dep(task285);
  densityq->add_task(task338);

  auto tensor339 = vector<shared_ptr<Tensor>>{I385, Gamma7_(), t2};
  auto task339 = make_shared<Task339>(tensor339, pindex);
  task338->add_dep(task339);
  task339->add_dep(task285);
  densityq->add_task(task339);

  vector<IndexRange> I541_index = {virt_, active_};
  auto I541 = make_shared<Tensor>(I541_index);
  auto tensor340 = vector<shared_ptr<Tensor>>{I381, t2, I541};
  auto task340 = make_shared<Task340>(tensor340, pindex);
  task335->add_dep(task340);
  task340->add_dep(task285);
  densityq->add_task(task340);

  auto tensor341 = vector<shared_ptr<Tensor>>{I541, Gamma60_(), t2};
  auto task341 = make_shared<Task341>(tensor341, pindex);
  task340->add_dep(task341);
  task341->add_dep(task285);
  densityq->add_task(task341);

  vector<IndexRange> I544_index = {virt_, active_};
  auto I544 = make_shared<Tensor>(I544_index);
  auto tensor342 = vector<shared_ptr<Tensor>>{I381, t2, I544};
  auto task342 = make_shared<Task342>(tensor342, pindex);
  task335->add_dep(task342);
  task342->add_dep(task285);
  densityq->add_task(task342);

  auto tensor343 = vector<shared_ptr<Tensor>>{I544, Gamma60_(), t2};
  auto task343 = make_shared<Task343>(tensor343, pindex);
  task342->add_dep(task343);
  task343->add_dep(task285);
  densityq->add_task(task343);
}

#endif
