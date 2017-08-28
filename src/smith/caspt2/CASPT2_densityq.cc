//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_densityqq.cc
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
#include <src/smith/caspt2/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_densityq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto densityq = make_shared<Queue>();
  auto tensor265 = vector<shared_ptr<Tensor>>{den2};
  auto task265 = make_shared<Task265>(tensor265, reset);
  densityq->add_task(task265);

  vector<IndexRange> I360_index = {active_, active_};
  auto I360 = make_shared<Tensor>(I360_index);
  auto tensor266 = vector<shared_ptr<Tensor>>{den2, I360};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task266->add_dep(task265);
  densityq->add_task(task266);

  vector<IndexRange> I361_index = {active_, active_, active_, active_};
  auto I361 = make_shared<Tensor>(I361_index);
  auto tensor267 = vector<shared_ptr<Tensor>>{I360, Gamma138_(), I361};
  auto task267 = make_shared<Task267>(tensor267, pindex);
  task266->add_dep(task267);
  task267->add_dep(task265);
  densityq->add_task(task267);

  auto tensor268 = vector<shared_ptr<Tensor>>{I361, t2};
  auto task268 = make_shared<Task268>(tensor268, pindex);
  task267->add_dep(task268);
  task268->add_dep(task265);
  densityq->add_task(task268);

  vector<IndexRange> I454_index = {active_, active_, active_, active_};
  auto I454 = make_shared<Tensor>(I454_index);
  auto tensor269 = vector<shared_ptr<Tensor>>{I360, Gamma169_(), I454};
  auto task269 = make_shared<Task269>(tensor269, pindex);
  task266->add_dep(task269);
  task269->add_dep(task265);
  densityq->add_task(task269);

  auto tensor270 = vector<shared_ptr<Tensor>>{I454, t2};
  auto task270 = make_shared<Task270>(tensor270, pindex);
  task269->add_dep(task270);
  task270->add_dep(task265);
  densityq->add_task(task270);

  vector<IndexRange> I463_index = {active_, active_, active_, active_};
  auto I463 = make_shared<Tensor>(I463_index);
  auto tensor271 = vector<shared_ptr<Tensor>>{I360, Gamma172_(), I463};
  auto task271 = make_shared<Task271>(tensor271, pindex);
  task266->add_dep(task271);
  task271->add_dep(task265);
  densityq->add_task(task271);

  auto tensor272 = vector<shared_ptr<Tensor>>{I463, t2};
  auto task272 = make_shared<Task272>(tensor272, pindex);
  task271->add_dep(task272);
  task272->add_dep(task265);
  densityq->add_task(task272);

  vector<IndexRange> I497_index = {active_, virt_, closed_, active_};
  auto I497 = make_shared<Tensor>(I497_index);
  auto tensor273 = vector<shared_ptr<Tensor>>{I463, t2, I497};
  auto task273 = make_shared<Task273>(tensor273, pindex);
  task271->add_dep(task273);
  task273->add_dep(task265);
  densityq->add_task(task273);

  auto tensor274 = vector<shared_ptr<Tensor>>{I497, t2};
  auto task274 = make_shared<Task274>(tensor274, pindex);
  task273->add_dep(task274);
  task274->add_dep(task265);
  densityq->add_task(task274);

  vector<IndexRange> I645_index = {active_, active_, active_, active_};
  auto I645 = make_shared<Tensor>(I645_index);
  auto tensor275 = vector<shared_ptr<Tensor>>{I360, Gamma230_(), I645};
  auto task275 = make_shared<Task275>(tensor275, pindex);
  task266->add_dep(task275);
  task275->add_dep(task265);
  densityq->add_task(task275);

  auto tensor276 = vector<shared_ptr<Tensor>>{I645, t2};
  auto task276 = make_shared<Task276>(tensor276, pindex);
  task275->add_dep(task276);
  task276->add_dep(task265);
  densityq->add_task(task276);

  vector<IndexRange> I363_index = {closed_, closed_};
  auto I363 = make_shared<Tensor>(I363_index);
  auto tensor277 = vector<shared_ptr<Tensor>>{den2, I363};
  auto task277 = make_shared<Task277>(tensor277, pindex);
  task277->add_dep(task265);
  densityq->add_task(task277);

  vector<IndexRange> I364_index = {closed_, closed_, active_, active_};
  auto I364 = make_shared<Tensor>(I364_index);
  auto tensor278 = vector<shared_ptr<Tensor>>{I363, t2, I364};
  auto task278 = make_shared<Task278>(tensor278, pindex);
  task277->add_dep(task278);
  task278->add_dep(task265);
  densityq->add_task(task278);

  auto tensor279 = vector<shared_ptr<Tensor>>{I364, Gamma92_(), t2};
  auto task279 = make_shared<Task279>(tensor279, pindex);
  task278->add_dep(task279);
  task279->add_dep(task265);
  densityq->add_task(task279);

  vector<IndexRange> I457_index = {closed_, virt_, active_, active_};
  auto I457 = make_shared<Tensor>(I457_index);
  auto tensor280 = vector<shared_ptr<Tensor>>{I363, t2, I457};
  auto task280 = make_shared<Task280>(tensor280, pindex);
  task277->add_dep(task280);
  task280->add_dep(task265);
  densityq->add_task(task280);

  auto tensor281 = vector<shared_ptr<Tensor>>{I457, Gamma32_(), t2};
  auto task281 = make_shared<Task281>(tensor281, pindex);
  task280->add_dep(task281);
  task281->add_dep(task265);
  densityq->add_task(task281);

  vector<IndexRange> I466_index = {closed_, virt_, active_, active_};
  auto I466 = make_shared<Tensor>(I466_index);
  auto tensor282 = vector<shared_ptr<Tensor>>{I363, t2, I466};
  auto task282 = make_shared<Task282>(tensor282, pindex);
  task277->add_dep(task282);
  task282->add_dep(task265);
  densityq->add_task(task282);

  auto tensor283 = vector<shared_ptr<Tensor>>{I466, Gamma35_(), t2};
  auto task283 = make_shared<Task283>(tensor283, pindex);
  task282->add_dep(task283);
  task283->add_dep(task265);
  densityq->add_task(task283);

  vector<IndexRange> I366_index = {active_, closed_};
  auto I366 = make_shared<Tensor>(I366_index);
  auto tensor284 = vector<shared_ptr<Tensor>>{den2, I366};
  auto task284 = make_shared<Task284>(tensor284, pindex);
  task284->add_dep(task265);
  densityq->add_task(task284);

  vector<IndexRange> I367_index = {closed_, active_, active_, active_};
  auto I367 = make_shared<Tensor>(I367_index);
  auto tensor285 = vector<shared_ptr<Tensor>>{I366, t2, I367};
  auto task285 = make_shared<Task285>(tensor285, pindex);
  task284->add_dep(task285);
  task285->add_dep(task265);
  densityq->add_task(task285);

  auto tensor286 = vector<shared_ptr<Tensor>>{I367, Gamma2_(), t2};
  auto task286 = make_shared<Task286>(tensor286, pindex);
  task285->add_dep(task286);
  task286->add_dep(task265);
  densityq->add_task(task286);

  vector<IndexRange> I472_index = {virt_, active_, active_, active_};
  auto I472 = make_shared<Tensor>(I472_index);
  auto tensor287 = vector<shared_ptr<Tensor>>{I366, t2, I472};
  auto task287 = make_shared<Task287>(tensor287, pindex);
  task284->add_dep(task287);
  task287->add_dep(task265);
  densityq->add_task(task287);

  auto tensor288 = vector<shared_ptr<Tensor>>{I472, Gamma37_(), t2};
  auto task288 = make_shared<Task288>(tensor288, pindex);
  task287->add_dep(task288);
  task288->add_dep(task265);
  densityq->add_task(task288);

  vector<IndexRange> I369_index = {active_, virt_};
  auto I369 = make_shared<Tensor>(I369_index);
  auto tensor289 = vector<shared_ptr<Tensor>>{den2, I369};
  auto task289 = make_shared<Task289>(tensor289, pindex);
  task289->add_dep(task265);
  densityq->add_task(task289);

  vector<IndexRange> I370_index = {closed_, closed_, active_, active_};
  auto I370 = make_shared<Tensor>(I370_index);
  auto tensor290 = vector<shared_ptr<Tensor>>{I369, t2, I370};
  auto task290 = make_shared<Task290>(tensor290, pindex);
  task289->add_dep(task290);
  task290->add_dep(task265);
  densityq->add_task(task290);

  auto tensor291 = vector<shared_ptr<Tensor>>{I370, Gamma3_(), t2};
  auto task291 = make_shared<Task291>(tensor291, pindex);
  task290->add_dep(task291);
  task291->add_dep(task265);
  densityq->add_task(task291);

  vector<IndexRange> I481_index = {closed_, virt_, active_, active_};
  auto I481 = make_shared<Tensor>(I481_index);
  auto tensor292 = vector<shared_ptr<Tensor>>{I369, t2, I481};
  auto task292 = make_shared<Task292>(tensor292, pindex);
  task289->add_dep(task292);
  task292->add_dep(task265);
  densityq->add_task(task292);

  auto tensor293 = vector<shared_ptr<Tensor>>{I481, Gamma35_(), t2};
  auto task293 = make_shared<Task293>(tensor293, pindex);
  task292->add_dep(task293);
  task293->add_dep(task265);
  densityq->add_task(task293);

  vector<IndexRange> I484_index = {closed_, virt_, active_, active_};
  auto I484 = make_shared<Tensor>(I484_index);
  auto tensor294 = vector<shared_ptr<Tensor>>{I369, t2, I484};
  auto task294 = make_shared<Task294>(tensor294, pindex);
  task289->add_dep(task294);
  task294->add_dep(task265);
  densityq->add_task(task294);

  auto tensor295 = vector<shared_ptr<Tensor>>{I484, Gamma32_(), t2};
  auto task295 = make_shared<Task295>(tensor295, pindex);
  task294->add_dep(task295);
  task295->add_dep(task265);
  densityq->add_task(task295);

  vector<IndexRange> I523_index = {virt_, closed_, active_, active_};
  auto I523 = make_shared<Tensor>(I523_index);
  auto tensor296 = vector<shared_ptr<Tensor>>{I369, t2, I523};
  auto task296 = make_shared<Task296>(tensor296, pindex);
  task289->add_dep(task296);
  task296->add_dep(task265);
  densityq->add_task(task296);

  auto tensor297 = vector<shared_ptr<Tensor>>{I523, Gamma35_(), t2};
  auto task297 = make_shared<Task297>(tensor297, pindex);
  task296->add_dep(task297);
  task297->add_dep(task265);
  densityq->add_task(task297);

  vector<IndexRange> I526_index = {virt_, closed_, active_, active_};
  auto I526 = make_shared<Tensor>(I526_index);
  auto tensor298 = vector<shared_ptr<Tensor>>{I369, t2, I526};
  auto task298 = make_shared<Task298>(tensor298, pindex);
  task289->add_dep(task298);
  task298->add_dep(task265);
  densityq->add_task(task298);

  auto tensor299 = vector<shared_ptr<Tensor>>{I526, Gamma35_(), t2};
  auto task299 = make_shared<Task299>(tensor299, pindex);
  task298->add_dep(task299);
  task299->add_dep(task265);
  densityq->add_task(task299);

  vector<IndexRange> I372_index = {active_, closed_};
  auto I372 = make_shared<Tensor>(I372_index);
  auto tensor300 = vector<shared_ptr<Tensor>>{den2, I372};
  auto task300 = make_shared<Task300>(tensor300, pindex);
  task300->add_dep(task265);
  densityq->add_task(task300);

  vector<IndexRange> I373_index = {closed_, active_, active_, active_};
  auto I373 = make_shared<Tensor>(I373_index);
  auto tensor301 = vector<shared_ptr<Tensor>>{I372, t2, I373};
  auto task301 = make_shared<Task301>(tensor301, pindex);
  task300->add_dep(task301);
  task301->add_dep(task265);
  densityq->add_task(task301);

  auto tensor302 = vector<shared_ptr<Tensor>>{I373, Gamma4_(), t2};
  auto task302 = make_shared<Task302>(tensor302, pindex);
  task301->add_dep(task302);
  task302->add_dep(task265);
  densityq->add_task(task302);

  vector<IndexRange> I529_index = {virt_, active_, active_, active_};
  auto I529 = make_shared<Tensor>(I529_index);
  auto tensor303 = vector<shared_ptr<Tensor>>{I372, t2, I529};
  auto task303 = make_shared<Task303>(tensor303, pindex);
  task300->add_dep(task303);
  task303->add_dep(task265);
  densityq->add_task(task303);

  auto tensor304 = vector<shared_ptr<Tensor>>{I529, Gamma56_(), t2};
  auto task304 = make_shared<Task304>(tensor304, pindex);
  task303->add_dep(task304);
  task304->add_dep(task265);
  densityq->add_task(task304);

  vector<IndexRange> I532_index = {virt_, active_, active_, active_};
  auto I532 = make_shared<Tensor>(I532_index);
  auto tensor305 = vector<shared_ptr<Tensor>>{I372, t2, I532};
  auto task305 = make_shared<Task305>(tensor305, pindex);
  task300->add_dep(task305);
  task305->add_dep(task265);
  densityq->add_task(task305);

  auto tensor306 = vector<shared_ptr<Tensor>>{I532, Gamma57_(), t2};
  auto task306 = make_shared<Task306>(tensor306, pindex);
  task305->add_dep(task306);
  task306->add_dep(task265);
  densityq->add_task(task306);

  vector<IndexRange> I375_index = {active_, active_};
  auto I375 = make_shared<Tensor>(I375_index);
  auto tensor307 = vector<shared_ptr<Tensor>>{den2, I375};
  auto task307 = make_shared<Task307>(tensor307, pindex);
  task307->add_dep(task265);
  densityq->add_task(task307);

  vector<IndexRange> I376_index = {active_, active_, active_, active_, active_, active_};
  auto I376 = make_shared<Tensor>(I376_index);
  auto tensor308 = vector<shared_ptr<Tensor>>{I375, Gamma143_(), I376};
  auto task308 = make_shared<Task308>(tensor308, pindex);
  task307->add_dep(task308);
  task308->add_dep(task265);
  densityq->add_task(task308);

  auto tensor309 = vector<shared_ptr<Tensor>>{I376, t2};
  auto task309 = make_shared<Task309>(tensor309, pindex);
  task308->add_dep(task309);
  task309->add_dep(task265);
  densityq->add_task(task309);

  vector<IndexRange> I535_index = {active_, active_, active_, active_, active_, active_};
  auto I535 = make_shared<Tensor>(I535_index);
  auto tensor310 = vector<shared_ptr<Tensor>>{I375, Gamma196_(), I535};
  auto task310 = make_shared<Task310>(tensor310, pindex);
  task307->add_dep(task310);
  task310->add_dep(task265);
  densityq->add_task(task310);

  auto tensor311 = vector<shared_ptr<Tensor>>{I535, t2};
  auto task311 = make_shared<Task311>(tensor311, pindex);
  task310->add_dep(task311);
  task311->add_dep(task265);
  densityq->add_task(task311);

  vector<IndexRange> I378_index = {closed_, closed_};
  auto I378 = make_shared<Tensor>(I378_index);
  auto tensor312 = vector<shared_ptr<Tensor>>{den2, I378};
  auto task312 = make_shared<Task312>(tensor312, pindex);
  task312->add_dep(task265);
  densityq->add_task(task312);

  vector<IndexRange> I379_index = {closed_, active_, active_, active_};
  auto I379 = make_shared<Tensor>(I379_index);
  auto tensor313 = vector<shared_ptr<Tensor>>{I378, t2, I379};
  auto task313 = make_shared<Task313>(tensor313, pindex);
  task312->add_dep(task313);
  task313->add_dep(task265);
  densityq->add_task(task313);

  auto tensor314 = vector<shared_ptr<Tensor>>{I379, Gamma6_(), t2};
  auto task314 = make_shared<Task314>(tensor314, pindex);
  task313->add_dep(task314);
  task314->add_dep(task265);
  densityq->add_task(task314);

  vector<IndexRange> I381_index = {closed_, virt_};
  auto I381 = make_shared<Tensor>(I381_index);
  auto tensor315 = vector<shared_ptr<Tensor>>{den2, I381};
  auto task315 = make_shared<Task315>(tensor315, pindex);
  task315->add_dep(task265);
  densityq->add_task(task315);

  vector<IndexRange> I382_index = {closed_, active_};
  auto I382 = make_shared<Tensor>(I382_index);
  auto tensor316 = vector<shared_ptr<Tensor>>{I381, t2, I382};
  auto task316 = make_shared<Task316>(tensor316, pindex);
  task315->add_dep(task316);
  task316->add_dep(task265);
  densityq->add_task(task316);

  auto tensor317 = vector<shared_ptr<Tensor>>{I382, Gamma7_(), t2};
  auto task317 = make_shared<Task317>(tensor317, pindex);
  task316->add_dep(task317);
  task317->add_dep(task265);
  densityq->add_task(task317);

  vector<IndexRange> I385_index = {closed_, active_};
  auto I385 = make_shared<Tensor>(I385_index);
  auto tensor318 = vector<shared_ptr<Tensor>>{I381, t2, I385};
  auto task318 = make_shared<Task318>(tensor318, pindex);
  task315->add_dep(task318);
  task318->add_dep(task265);
  densityq->add_task(task318);

  auto tensor319 = vector<shared_ptr<Tensor>>{I385, Gamma7_(), t2};
  auto task319 = make_shared<Task319>(tensor319, pindex);
  task318->add_dep(task319);
  task319->add_dep(task265);
  densityq->add_task(task319);

  vector<IndexRange> I541_index = {virt_, active_};
  auto I541 = make_shared<Tensor>(I541_index);
  auto tensor320 = vector<shared_ptr<Tensor>>{I381, t2, I541};
  auto task320 = make_shared<Task320>(tensor320, pindex);
  task315->add_dep(task320);
  task320->add_dep(task265);
  densityq->add_task(task320);

  auto tensor321 = vector<shared_ptr<Tensor>>{I541, Gamma60_(), t2};
  auto task321 = make_shared<Task321>(tensor321, pindex);
  task320->add_dep(task321);
  task321->add_dep(task265);
  densityq->add_task(task321);

  vector<IndexRange> I544_index = {virt_, active_};
  auto I544 = make_shared<Tensor>(I544_index);
  auto tensor322 = vector<shared_ptr<Tensor>>{I381, t2, I544};
  auto task322 = make_shared<Task322>(tensor322, pindex);
  task315->add_dep(task322);
  task322->add_dep(task265);
  densityq->add_task(task322);

  auto tensor323 = vector<shared_ptr<Tensor>>{I544, Gamma60_(), t2};
  auto task323 = make_shared<Task323>(tensor323, pindex);
  task322->add_dep(task323);
  task323->add_dep(task265);
  densityq->add_task(task323);

  vector<IndexRange> I387_index = {active_, virt_};
  auto I387 = make_shared<Tensor>(I387_index);
  auto tensor324 = vector<shared_ptr<Tensor>>{den2, I387};
  auto task324 = make_shared<Task324>(tensor324, pindex);
  task324->add_dep(task265);
  densityq->add_task(task324);

  vector<IndexRange> I388_index = {closed_, active_, active_, active_};
  auto I388 = make_shared<Tensor>(I388_index);
  auto tensor325 = vector<shared_ptr<Tensor>>{I387, t2, I388};
  auto task325 = make_shared<Task325>(tensor325, pindex);
  task324->add_dep(task325);
  task325->add_dep(task265);
  densityq->add_task(task325);

  auto tensor326 = vector<shared_ptr<Tensor>>{I388, Gamma9_(), t2};
  auto task326 = make_shared<Task326>(tensor326, pindex);
  task325->add_dep(task326);
  task326->add_dep(task265);
  densityq->add_task(task326);

  vector<IndexRange> I391_index = {closed_, active_, active_, active_};
  auto I391 = make_shared<Tensor>(I391_index);
  auto tensor327 = vector<shared_ptr<Tensor>>{I387, t2, I391};
  auto task327 = make_shared<Task327>(tensor327, pindex);
  task324->add_dep(task327);
  task327->add_dep(task265);
  densityq->add_task(task327);

  auto tensor328 = vector<shared_ptr<Tensor>>{I391, Gamma6_(), t2};
  auto task328 = make_shared<Task328>(tensor328, pindex);
  task327->add_dep(task328);
  task328->add_dep(task265);
  densityq->add_task(task328);

  vector<IndexRange> I547_index = {virt_, active_, active_, active_};
  auto I547 = make_shared<Tensor>(I547_index);
  auto tensor329 = vector<shared_ptr<Tensor>>{I387, t2, I547};
  auto task329 = make_shared<Task329>(tensor329, pindex);
  task324->add_dep(task329);
  task329->add_dep(task265);
  densityq->add_task(task329);

  auto tensor330 = vector<shared_ptr<Tensor>>{I547, Gamma59_(), t2};
  auto task330 = make_shared<Task330>(tensor330, pindex);
  task329->add_dep(task330);
  task330->add_dep(task265);
  densityq->add_task(task330);

  vector<IndexRange> I393_index = {active_, virt_};
  auto I393 = make_shared<Tensor>(I393_index);
  auto tensor331 = vector<shared_ptr<Tensor>>{den2, I393};
  auto task331 = make_shared<Task331>(tensor331, pindex);
  task331->add_dep(task265);
  densityq->add_task(task331);

  vector<IndexRange> I394_index = {closed_, closed_, active_, active_};
  auto I394 = make_shared<Tensor>(I394_index);
  auto tensor332 = vector<shared_ptr<Tensor>>{I393, t2, I394};
  auto task332 = make_shared<Task332>(tensor332, pindex);
  task331->add_dep(task332);
  task332->add_dep(task265);
  densityq->add_task(task332);

  auto tensor333 = vector<shared_ptr<Tensor>>{I394, Gamma3_(), t2};
  auto task333 = make_shared<Task333>(tensor333, pindex);
  task332->add_dep(task333);
  task333->add_dep(task265);
  densityq->add_task(task333);

  vector<IndexRange> I396_index = {virt_, closed_};
  auto I396 = make_shared<Tensor>(I396_index);
  auto tensor334 = vector<shared_ptr<Tensor>>{den2, I396};
  auto task334 = make_shared<Task334>(tensor334, pindex);
  task334->add_dep(task265);
  densityq->add_task(task334);

  vector<IndexRange> I397_index = {closed_, active_};
  auto I397 = make_shared<Tensor>(I397_index);
  auto tensor335 = vector<shared_ptr<Tensor>>{I396, t2, I397};
  auto task335 = make_shared<Task335>(tensor335, pindex);
  task334->add_dep(task335);
  task335->add_dep(task265);
  densityq->add_task(task335);

  auto tensor336 = vector<shared_ptr<Tensor>>{I397, Gamma12_(), t2};
  auto task336 = make_shared<Task336>(tensor336, pindex);
  task335->add_dep(task336);
  task336->add_dep(task265);
  densityq->add_task(task336);

  vector<IndexRange> I399_index = {closed_, virt_};
  auto I399 = make_shared<Tensor>(I399_index);
  auto tensor337 = vector<shared_ptr<Tensor>>{den2, I399};
  auto task337 = make_shared<Task337>(tensor337, pindex);
  task337->add_dep(task265);
  densityq->add_task(task337);

  vector<IndexRange> I400_index = {closed_, active_};
  auto I400 = make_shared<Tensor>(I400_index);
  auto tensor338 = vector<shared_ptr<Tensor>>{I399, t2, I400};
  auto task338 = make_shared<Task338>(tensor338, pindex);
  task337->add_dep(task338);
  task338->add_dep(task265);
  densityq->add_task(task338);

  auto tensor339 = vector<shared_ptr<Tensor>>{I400, Gamma12_(), t2};
  auto task339 = make_shared<Task339>(tensor339, pindex);
  task338->add_dep(task339);
  task339->add_dep(task265);
  densityq->add_task(task339);

  vector<IndexRange> I556_index = {virt_, closed_};
  auto I556 = make_shared<Tensor>(I556_index);
  auto tensor340 = vector<shared_ptr<Tensor>>{I399, t2, I556};
  auto task340 = make_shared<Task340>(tensor340, pindex);
  task337->add_dep(task340);
  task340->add_dep(task265);
  densityq->add_task(task340);

  vector<IndexRange> I557_index = {active_, virt_, closed_, active_};
  auto I557 = make_shared<Tensor>(I557_index);
  auto tensor341 = vector<shared_ptr<Tensor>>{I556, Gamma38_(), I557};
  auto task341 = make_shared<Task341>(tensor341, pindex);
  task340->add_dep(task341);
  task341->add_dep(task265);
  densityq->add_task(task341);

  auto tensor342 = vector<shared_ptr<Tensor>>{I557, t2};
  auto task342 = make_shared<Task342>(tensor342, pindex);
  task341->add_dep(task342);
  task342->add_dep(task265);
  densityq->add_task(task342);

  vector<IndexRange> I402_index = {active_, active_};
  auto I402 = make_shared<Tensor>(I402_index);
  auto tensor343 = vector<shared_ptr<Tensor>>{den2, I402};
  auto task343 = make_shared<Task343>(tensor343, pindex);
  task343->add_dep(task265);
  densityq->add_task(task343);

  vector<IndexRange> I403_index = {active_, active_};
  auto I403 = make_shared<Tensor>(I403_index);
  auto tensor344 = vector<shared_ptr<Tensor>>{I402, Gamma152_(), I403};
  auto task344 = make_shared<Task344>(tensor344, pindex);
  task343->add_dep(task344);
  task344->add_dep(task265);
  densityq->add_task(task344);

  vector<IndexRange> I404_index = {closed_, virt_, closed_, active_};
  auto I404 = make_shared<Tensor>(I404_index);
  auto tensor345 = vector<shared_ptr<Tensor>>{I403, t2, I404};
  auto task345 = make_shared<Task345>(tensor345, pindex);
  task344->add_dep(task345);
  task345->add_dep(task265);
  densityq->add_task(task345);

  auto tensor346 = vector<shared_ptr<Tensor>>{I404, t2};
  auto task346 = make_shared<Task346>(tensor346, pindex);
  task345->add_dep(task346);
  task346->add_dep(task265);
  densityq->add_task(task346);

  vector<IndexRange> I612_index = {active_, active_};
  auto I612 = make_shared<Tensor>(I612_index);
  auto tensor347 = vector<shared_ptr<Tensor>>{I402, Gamma60_(), I612};
  auto task347 = make_shared<Task347>(tensor347, pindex);
  task343->add_dep(task347);
  task347->add_dep(task265);
  densityq->add_task(task347);

  vector<IndexRange> I613_index = {active_, virt_, closed_, virt_};
  auto I613 = make_shared<Tensor>(I613_index);
  auto tensor348 = vector<shared_ptr<Tensor>>{I612, t2, I613};
  auto task348 = make_shared<Task348>(tensor348, pindex);
  task347->add_dep(task348);
  task348->add_dep(task265);
  densityq->add_task(task348);

  auto tensor349 = vector<shared_ptr<Tensor>>{I613, t2};
  auto task349 = make_shared<Task349>(tensor349, pindex);
  task348->add_dep(task349);
  task349->add_dep(task265);
  densityq->add_task(task349);

  vector<IndexRange> I408_index = {closed_, closed_};
  auto I408 = make_shared<Tensor>(I408_index);
  auto tensor350 = vector<shared_ptr<Tensor>>{den2, I408};
  auto task350 = make_shared<Task350>(tensor350, pindex);
  task350->add_dep(task265);
  densityq->add_task(task350);

  vector<IndexRange> I409_index = {closed_, virt_, closed_, active_};
  auto I409 = make_shared<Tensor>(I409_index);
  auto tensor351 = vector<shared_ptr<Tensor>>{I408, t2, I409};
  auto task351 = make_shared<Task351>(tensor351, pindex);
  task350->add_dep(task351);
  task351->add_dep(task265);
  densityq->add_task(task351);

  auto tensor352 = vector<shared_ptr<Tensor>>{I409, Gamma16_(), t2};
  auto task352 = make_shared<Task352>(tensor352, pindex);
  task351->add_dep(task352);
  task352->add_dep(task265);
  densityq->add_task(task352);

  vector<IndexRange> I412_index = {closed_, virt_, closed_, active_};
  auto I412 = make_shared<Tensor>(I412_index);
  auto tensor353 = vector<shared_ptr<Tensor>>{I408, t2, I412};
  auto task353 = make_shared<Task353>(tensor353, pindex);
  task350->add_dep(task353);
  task353->add_dep(task265);
  densityq->add_task(task353);

  auto tensor354 = vector<shared_ptr<Tensor>>{I412, Gamma16_(), t2};
  auto task354 = make_shared<Task354>(tensor354, pindex);
  task353->add_dep(task354);
  task354->add_dep(task265);
  densityq->add_task(task354);

  vector<IndexRange> I414_index = {closed_, closed_};
  auto I414 = make_shared<Tensor>(I414_index);
  auto tensor355 = vector<shared_ptr<Tensor>>{den2, I414};
  auto task355 = make_shared<Task355>(tensor355, pindex);
  task355->add_dep(task265);
  densityq->add_task(task355);

  vector<IndexRange> I415_index = {closed_, virt_, closed_, active_};
  auto I415 = make_shared<Tensor>(I415_index);
  auto tensor356 = vector<shared_ptr<Tensor>>{I414, t2, I415};
  auto task356 = make_shared<Task356>(tensor356, pindex);
  task355->add_dep(task356);
  task356->add_dep(task265);
  densityq->add_task(task356);

  auto tensor357 = vector<shared_ptr<Tensor>>{I415, Gamma16_(), t2};
  auto task357 = make_shared<Task357>(tensor357, pindex);
  task356->add_dep(task357);
  task357->add_dep(task265);
  densityq->add_task(task357);

  vector<IndexRange> I421_index = {closed_, virt_, closed_, active_};
  auto I421 = make_shared<Tensor>(I421_index);
  auto tensor358 = vector<shared_ptr<Tensor>>{I414, t2, I421};
  auto task358 = make_shared<Task358>(tensor358, pindex);
  task355->add_dep(task358);
  task358->add_dep(task265);
  densityq->add_task(task358);

  auto tensor359 = vector<shared_ptr<Tensor>>{I421, Gamma16_(), t2};
  auto task359 = make_shared<Task359>(tensor359, pindex);
  task358->add_dep(task359);
  task359->add_dep(task265);
  densityq->add_task(task359);

  vector<IndexRange> I417_index = {virt_, virt_};
  auto I417 = make_shared<Tensor>(I417_index);
  auto tensor360 = vector<shared_ptr<Tensor>>{den2, I417};
  auto task360 = make_shared<Task360>(tensor360, pindex);
  task360->add_dep(task265);
  densityq->add_task(task360);

  vector<IndexRange> I418_index = {closed_, virt_, closed_, active_};
  auto I418 = make_shared<Tensor>(I418_index);
  auto tensor361 = vector<shared_ptr<Tensor>>{I417, t2, I418};
  auto task361 = make_shared<Task361>(tensor361, pindex);
  task360->add_dep(task361);
  task361->add_dep(task265);
  densityq->add_task(task361);

  auto tensor362 = vector<shared_ptr<Tensor>>{I418, Gamma16_(), t2};
  auto task362 = make_shared<Task362>(tensor362, pindex);
  task361->add_dep(task362);
  task362->add_dep(task265);
  densityq->add_task(task362);

  vector<IndexRange> I424_index = {closed_, virt_, closed_, active_};
  auto I424 = make_shared<Tensor>(I424_index);
  auto tensor363 = vector<shared_ptr<Tensor>>{I417, t2, I424};
  auto task363 = make_shared<Task363>(tensor363, pindex);
  task360->add_dep(task363);
  task363->add_dep(task265);
  densityq->add_task(task363);

  auto tensor364 = vector<shared_ptr<Tensor>>{I424, Gamma16_(), t2};
  auto task364 = make_shared<Task364>(tensor364, pindex);
  task363->add_dep(task364);
  task364->add_dep(task265);
  densityq->add_task(task364);

  vector<IndexRange> I426_index = {active_, closed_};
  auto I426 = make_shared<Tensor>(I426_index);
  auto tensor365 = vector<shared_ptr<Tensor>>{den2, I426};
  auto task365 = make_shared<Task365>(tensor365, pindex);
  task365->add_dep(task265);
  densityq->add_task(task365);

  vector<IndexRange> I427_index = {virt_, closed_, active_, active_};
  auto I427 = make_shared<Tensor>(I427_index);
  auto tensor366 = vector<shared_ptr<Tensor>>{I426, t2, I427};
  auto task366 = make_shared<Task366>(tensor366, pindex);
  task365->add_dep(task366);
  task366->add_dep(task265);
  densityq->add_task(task366);

  auto tensor367 = vector<shared_ptr<Tensor>>{I427, Gamma22_(), t2};
  auto task367 = make_shared<Task367>(tensor367, pindex);
  task366->add_dep(task367);
  task367->add_dep(task265);
  densityq->add_task(task367);

  auto tensor368 = vector<shared_ptr<Tensor>>{I427, Gamma12_(), t2};
  auto task368 = make_shared<Task368>(tensor368, pindex);
  task366->add_dep(task368);
  task368->add_dep(task265);
  densityq->add_task(task368);

  vector<IndexRange> I429_index = {active_, closed_};
  auto I429 = make_shared<Tensor>(I429_index);
  auto tensor369 = vector<shared_ptr<Tensor>>{den2, I429};
  auto task369 = make_shared<Task369>(tensor369, pindex);
  task369->add_dep(task265);
  densityq->add_task(task369);

  vector<IndexRange> I430_index = {virt_, closed_, active_, active_};
  auto I430 = make_shared<Tensor>(I430_index);
  auto tensor370 = vector<shared_ptr<Tensor>>{I429, t2, I430};
  auto task370 = make_shared<Task370>(tensor370, pindex);
  task369->add_dep(task370);
  task370->add_dep(task265);
  densityq->add_task(task370);

  vector<IndexRange> I431_index = {active_, virt_, closed_, active_};
  auto I431 = make_shared<Tensor>(I431_index);
  auto tensor371 = vector<shared_ptr<Tensor>>{I430, Gamma12_(), I431};
  auto task371 = make_shared<Task371>(tensor371, pindex);
  task370->add_dep(task371);
  task371->add_dep(task265);
  densityq->add_task(task371);

  auto tensor372 = vector<shared_ptr<Tensor>>{I431, t2};
  auto task372 = make_shared<Task372>(tensor372, pindex);
  task371->add_dep(task372);
  task372->add_dep(task265);
  densityq->add_task(task372);

  vector<IndexRange> I438_index = {virt_, active_};
  auto I438 = make_shared<Tensor>(I438_index);
  auto tensor373 = vector<shared_ptr<Tensor>>{den2, I438};
  auto task373 = make_shared<Task373>(tensor373, pindex);
  task373->add_dep(task265);
  densityq->add_task(task373);

  vector<IndexRange> I439_index = {virt_, active_};
  auto I439 = make_shared<Tensor>(I439_index);
  auto tensor374 = vector<shared_ptr<Tensor>>{I438, Gamma16_(), I439};
  auto task374 = make_shared<Task374>(tensor374, pindex);
  task373->add_dep(task374);
  task374->add_dep(task265);
  densityq->add_task(task374);

  vector<IndexRange> I440_index = {closed_, virt_, closed_, virt_};
  auto I440 = make_shared<Tensor>(I440_index);
  auto tensor375 = vector<shared_ptr<Tensor>>{I439, t2, I440};
  auto task375 = make_shared<Task375>(tensor375, pindex);
  task374->add_dep(task375);
  task375->add_dep(task265);
  densityq->add_task(task375);

  auto tensor376 = vector<shared_ptr<Tensor>>{I440, t2};
  auto task376 = make_shared<Task376>(tensor376, pindex);
  task375->add_dep(task376);
  task376->add_dep(task265);
  densityq->add_task(task376);

  vector<IndexRange> I444_index = {active_, virt_};
  auto I444 = make_shared<Tensor>(I444_index);
  auto tensor377 = vector<shared_ptr<Tensor>>{den2, I444};
  auto task377 = make_shared<Task377>(tensor377, pindex);
  task377->add_dep(task265);
  densityq->add_task(task377);

  vector<IndexRange> I445_index = {closed_, active_, active_, active_};
  auto I445 = make_shared<Tensor>(I445_index);
  auto tensor378 = vector<shared_ptr<Tensor>>{I444, t2, I445};
  auto task378 = make_shared<Task378>(tensor378, pindex);
  task377->add_dep(task378);
  task378->add_dep(task265);
  densityq->add_task(task378);

  auto tensor379 = vector<shared_ptr<Tensor>>{I445, Gamma28_(), t2};
  auto task379 = make_shared<Task379>(tensor379, pindex);
  task378->add_dep(task379);
  task379->add_dep(task265);
  densityq->add_task(task379);

  vector<IndexRange> I447_index = {active_, closed_};
  auto I447 = make_shared<Tensor>(I447_index);
  auto tensor380 = vector<shared_ptr<Tensor>>{den2, I447};
  auto task380 = make_shared<Task380>(tensor380, pindex);
  task380->add_dep(task265);
  densityq->add_task(task380);

  vector<IndexRange> I448_index = {closed_, virt_, active_, active_};
  auto I448 = make_shared<Tensor>(I448_index);
  auto tensor381 = vector<shared_ptr<Tensor>>{I447, t2, I448};
  auto task381 = make_shared<Task381>(tensor381, pindex);
  task380->add_dep(task381);
  task381->add_dep(task265);
  densityq->add_task(task381);

  auto tensor382 = vector<shared_ptr<Tensor>>{I448, Gamma29_(), t2};
  auto task382 = make_shared<Task382>(tensor382, pindex);
  task381->add_dep(task382);
  task382->add_dep(task265);
  densityq->add_task(task382);

  vector<IndexRange> I451_index = {closed_, virt_, active_, active_};
  auto I451 = make_shared<Tensor>(I451_index);
  auto tensor383 = vector<shared_ptr<Tensor>>{I447, t2, I451};
  auto task383 = make_shared<Task383>(tensor383, pindex);
  task380->add_dep(task383);
  task383->add_dep(task265);
  densityq->add_task(task383);

  auto tensor384 = vector<shared_ptr<Tensor>>{I451, Gamma7_(), t2};
  auto task384 = make_shared<Task384>(tensor384, pindex);
  task383->add_dep(task384);
  task384->add_dep(task265);
  densityq->add_task(task384);

  vector<IndexRange> I490_index = {virt_, closed_, active_, active_};
  auto I490 = make_shared<Tensor>(I490_index);
  auto tensor385 = vector<shared_ptr<Tensor>>{I447, t2, I490};
  auto task385 = make_shared<Task385>(tensor385, pindex);
  task380->add_dep(task385);
  task385->add_dep(task265);
  densityq->add_task(task385);

  auto tensor386 = vector<shared_ptr<Tensor>>{I490, Gamma7_(), t2};
  auto task386 = make_shared<Task386>(tensor386, pindex);
  task385->add_dep(task386);
  task386->add_dep(task265);
  densityq->add_task(task386);

  vector<IndexRange> I493_index = {virt_, closed_, active_, active_};
  auto I493 = make_shared<Tensor>(I493_index);
  auto tensor387 = vector<shared_ptr<Tensor>>{I447, t2, I493};
  auto task387 = make_shared<Task387>(tensor387, pindex);
  task380->add_dep(task387);
  task387->add_dep(task265);
  densityq->add_task(task387);

  auto tensor388 = vector<shared_ptr<Tensor>>{I493, Gamma7_(), t2};
  auto task388 = make_shared<Task388>(tensor388, pindex);
  task387->add_dep(task388);
  task388->add_dep(task265);
  densityq->add_task(task388);

  vector<IndexRange> I642_index = {virt_, virt_, active_, active_};
  auto I642 = make_shared<Tensor>(I642_index);
  auto tensor389 = vector<shared_ptr<Tensor>>{I447, t2, I642};
  auto task389 = make_shared<Task389>(tensor389, pindex);
  task380->add_dep(task389);
  task389->add_dep(task265);
  densityq->add_task(task389);

  auto tensor390 = vector<shared_ptr<Tensor>>{I642, Gamma60_(), t2};
  auto task390 = make_shared<Task390>(tensor390, pindex);
  task389->add_dep(task390);
  task390->add_dep(task265);
  densityq->add_task(task390);

  vector<IndexRange> I459_index = {virt_, virt_};
  auto I459 = make_shared<Tensor>(I459_index);
  auto tensor391 = vector<shared_ptr<Tensor>>{den2, I459};
  auto task391 = make_shared<Task391>(tensor391, pindex);
  task391->add_dep(task265);
  densityq->add_task(task391);

  vector<IndexRange> I460_index = {closed_, virt_, active_, active_};
  auto I460 = make_shared<Tensor>(I460_index);
  auto tensor392 = vector<shared_ptr<Tensor>>{I459, t2, I460};
  auto task392 = make_shared<Task392>(tensor392, pindex);
  task391->add_dep(task392);
  task392->add_dep(task265);
  densityq->add_task(task392);

  auto tensor393 = vector<shared_ptr<Tensor>>{I460, Gamma32_(), t2};
  auto task393 = make_shared<Task393>(tensor393, pindex);
  task392->add_dep(task393);
  task393->add_dep(task265);
  densityq->add_task(task393);

  vector<IndexRange> I469_index = {closed_, virt_, active_, active_};
  auto I469 = make_shared<Tensor>(I469_index);
  auto tensor394 = vector<shared_ptr<Tensor>>{I459, t2, I469};
  auto task394 = make_shared<Task394>(tensor394, pindex);
  task391->add_dep(task394);
  task394->add_dep(task265);
  densityq->add_task(task394);

  auto tensor395 = vector<shared_ptr<Tensor>>{I469, Gamma35_(), t2};
  auto task395 = make_shared<Task395>(tensor395, pindex);
  task394->add_dep(task395);
  task395->add_dep(task265);
  densityq->add_task(task395);

  vector<IndexRange> I474_index = {virt_, closed_};
  auto I474 = make_shared<Tensor>(I474_index);
  auto tensor396 = vector<shared_ptr<Tensor>>{den2, I474};
  auto task396 = make_shared<Task396>(tensor396, pindex);
  task396->add_dep(task265);
  densityq->add_task(task396);

  vector<IndexRange> I475_index = {closed_, virt_};
  auto I475 = make_shared<Tensor>(I475_index);
  auto tensor397 = vector<shared_ptr<Tensor>>{I474, t2, I475};
  auto task397 = make_shared<Task397>(tensor397, pindex);
  task396->add_dep(task397);
  task397->add_dep(task265);
  densityq->add_task(task397);

  auto tensor398 = vector<shared_ptr<Tensor>>{I475, Gamma38_(), t2};
  auto task398 = make_shared<Task398>(tensor398, pindex);
  task397->add_dep(task398);
  task398->add_dep(task265);
  densityq->add_task(task398);

  vector<IndexRange> I478_index = {closed_, virt_};
  auto I478 = make_shared<Tensor>(I478_index);
  auto tensor399 = vector<shared_ptr<Tensor>>{I474, t2, I478};
  auto task399 = make_shared<Task399>(tensor399, pindex);
  task396->add_dep(task399);
  task399->add_dep(task265);
  densityq->add_task(task399);

  auto tensor400 = vector<shared_ptr<Tensor>>{I478, Gamma38_(), t2};
  auto task400 = make_shared<Task400>(tensor400, pindex);
  task399->add_dep(task400);
  task400->add_dep(task265);
  densityq->add_task(task400);

  vector<IndexRange> I517_index = {virt_, closed_};
  auto I517 = make_shared<Tensor>(I517_index);
  auto tensor401 = vector<shared_ptr<Tensor>>{I474, t2, I517};
  auto task401 = make_shared<Task401>(tensor401, pindex);
  task396->add_dep(task401);
  task401->add_dep(task265);
  densityq->add_task(task401);

  auto tensor402 = vector<shared_ptr<Tensor>>{I517, Gamma38_(), t2};
  auto task402 = make_shared<Task402>(tensor402, pindex);
  task401->add_dep(task402);
  task402->add_dep(task265);
  densityq->add_task(task402);

  vector<IndexRange> I520_index = {virt_, closed_};
  auto I520 = make_shared<Tensor>(I520_index);
  auto tensor403 = vector<shared_ptr<Tensor>>{I474, t2, I520};
  auto task403 = make_shared<Task403>(tensor403, pindex);
  task396->add_dep(task403);
  task403->add_dep(task265);
  densityq->add_task(task403);

  auto tensor404 = vector<shared_ptr<Tensor>>{I520, Gamma38_(), t2};
  auto task404 = make_shared<Task404>(tensor404, pindex);
  task403->add_dep(task404);
  task404->add_dep(task265);
  densityq->add_task(task404);

  vector<IndexRange> I486_index = {active_, virt_};
  auto I486 = make_shared<Tensor>(I486_index);
  auto tensor405 = vector<shared_ptr<Tensor>>{den2, I486};
  auto task405 = make_shared<Task405>(tensor405, pindex);
  task405->add_dep(task265);
  densityq->add_task(task405);

  vector<IndexRange> I487_index = {closed_, active_, active_, active_};
  auto I487 = make_shared<Tensor>(I487_index);
  auto tensor406 = vector<shared_ptr<Tensor>>{I486, t2, I487};
  auto task406 = make_shared<Task406>(tensor406, pindex);
  task405->add_dep(task406);
  task406->add_dep(task265);
  densityq->add_task(task406);

  auto tensor407 = vector<shared_ptr<Tensor>>{I487, Gamma6_(), t2};
  auto task407 = make_shared<Task407>(tensor407, pindex);
  task406->add_dep(task407);
  task407->add_dep(task265);
  densityq->add_task(task407);

  vector<IndexRange> I639_index = {virt_, active_, active_, active_};
  auto I639 = make_shared<Tensor>(I639_index);
  auto tensor408 = vector<shared_ptr<Tensor>>{I486, t2, I639};
  auto task408 = make_shared<Task408>(tensor408, pindex);
  task405->add_dep(task408);
  task408->add_dep(task265);
  densityq->add_task(task408);

  auto tensor409 = vector<shared_ptr<Tensor>>{I639, Gamma59_(), t2};
  auto task409 = make_shared<Task409>(tensor409, pindex);
  task408->add_dep(task409);
  task409->add_dep(task265);
  densityq->add_task(task409);

  vector<IndexRange> I498_index = {closed_, closed_};
  auto I498 = make_shared<Tensor>(I498_index);
  auto tensor410 = vector<shared_ptr<Tensor>>{den2, I498};
  auto task410 = make_shared<Task410>(tensor410, pindex);
  task410->add_dep(task265);
  densityq->add_task(task410);

  vector<IndexRange> I499_index = {virt_, closed_, active_, active_};
  auto I499 = make_shared<Tensor>(I499_index);
  auto tensor411 = vector<shared_ptr<Tensor>>{I498, t2, I499};
  auto task411 = make_shared<Task411>(tensor411, pindex);
  task410->add_dep(task411);
  task411->add_dep(task265);
  densityq->add_task(task411);

  auto tensor412 = vector<shared_ptr<Tensor>>{I499, Gamma35_(), t2};
  auto task412 = make_shared<Task412>(tensor412, pindex);
  task411->add_dep(task412);
  task412->add_dep(task265);
  densityq->add_task(task412);

  vector<IndexRange> I508_index = {virt_, closed_, active_, active_};
  auto I508 = make_shared<Tensor>(I508_index);
  auto tensor413 = vector<shared_ptr<Tensor>>{I498, t2, I508};
  auto task413 = make_shared<Task413>(tensor413, pindex);
  task410->add_dep(task413);
  task413->add_dep(task265);
  densityq->add_task(task413);

  auto tensor414 = vector<shared_ptr<Tensor>>{I508, Gamma35_(), t2};
  auto task414 = make_shared<Task414>(tensor414, pindex);
  task413->add_dep(task414);
  task414->add_dep(task265);
  densityq->add_task(task414);

  vector<IndexRange> I501_index = {virt_, virt_};
  auto I501 = make_shared<Tensor>(I501_index);
  auto tensor415 = vector<shared_ptr<Tensor>>{den2, I501};
  auto task415 = make_shared<Task415>(tensor415, pindex);
  task415->add_dep(task265);
  densityq->add_task(task415);

  vector<IndexRange> I502_index = {virt_, closed_, active_, active_};
  auto I502 = make_shared<Tensor>(I502_index);
  auto tensor416 = vector<shared_ptr<Tensor>>{I501, t2, I502};
  auto task416 = make_shared<Task416>(tensor416, pindex);
  task415->add_dep(task416);
  task416->add_dep(task265);
  densityq->add_task(task416);

  auto tensor417 = vector<shared_ptr<Tensor>>{I502, Gamma35_(), t2};
  auto task417 = make_shared<Task417>(tensor417, pindex);
  task416->add_dep(task417);
  task417->add_dep(task265);
  densityq->add_task(task417);

  vector<IndexRange> I511_index = {virt_, closed_, active_, active_};
  auto I511 = make_shared<Tensor>(I511_index);
  auto tensor418 = vector<shared_ptr<Tensor>>{I501, t2, I511};
  auto task418 = make_shared<Task418>(tensor418, pindex);
  task415->add_dep(task418);
  task418->add_dep(task265);
  densityq->add_task(task418);

  auto tensor419 = vector<shared_ptr<Tensor>>{I511, Gamma35_(), t2};
  auto task419 = make_shared<Task419>(tensor419, pindex);
  task418->add_dep(task419);
  task419->add_dep(task265);
  densityq->add_task(task419);

  vector<IndexRange> I648_index = {virt_, virt_, active_, active_};
  auto I648 = make_shared<Tensor>(I648_index);
  auto tensor420 = vector<shared_ptr<Tensor>>{I501, t2, I648};
  auto task420 = make_shared<Task420>(tensor420, pindex);
  task415->add_dep(task420);
  task420->add_dep(task265);
  densityq->add_task(task420);

  auto tensor421 = vector<shared_ptr<Tensor>>{I648, Gamma60_(), t2};
  auto task421 = make_shared<Task421>(tensor421, pindex);
  task420->add_dep(task421);
  task421->add_dep(task265);
  densityq->add_task(task421);

  vector<IndexRange> I513_index = {active_, closed_};
  auto I513 = make_shared<Tensor>(I513_index);
  auto tensor422 = vector<shared_ptr<Tensor>>{den2, I513};
  auto task422 = make_shared<Task422>(tensor422, pindex);
  task422->add_dep(task265);
  densityq->add_task(task422);

  vector<IndexRange> I514_index = {virt_, active_, active_, active_};
  auto I514 = make_shared<Tensor>(I514_index);
  auto tensor423 = vector<shared_ptr<Tensor>>{I513, t2, I514};
  auto task423 = make_shared<Task423>(tensor423, pindex);
  task422->add_dep(task423);
  task423->add_dep(task265);
  densityq->add_task(task423);

  auto tensor424 = vector<shared_ptr<Tensor>>{I514, Gamma51_(), t2};
  auto task424 = make_shared<Task424>(tensor424, pindex);
  task423->add_dep(task424);
  task424->add_dep(task265);
  densityq->add_task(task424);

  vector<IndexRange> I537_index = {virt_, virt_};
  auto I537 = make_shared<Tensor>(I537_index);
  auto tensor425 = vector<shared_ptr<Tensor>>{den2, I537};
  auto task425 = make_shared<Task425>(tensor425, pindex);
  task425->add_dep(task265);
  densityq->add_task(task425);

  vector<IndexRange> I538_index = {virt_, active_, active_, active_};
  auto I538 = make_shared<Tensor>(I538_index);
  auto tensor426 = vector<shared_ptr<Tensor>>{I537, t2, I538};
  auto task426 = make_shared<Task426>(tensor426, pindex);
  task425->add_dep(task426);
  task426->add_dep(task265);
  densityq->add_task(task426);

  auto tensor427 = vector<shared_ptr<Tensor>>{I538, Gamma59_(), t2};
  auto task427 = make_shared<Task427>(tensor427, pindex);
  task426->add_dep(task427);
  task427->add_dep(task265);
  densityq->add_task(task427);

  vector<IndexRange> I549_index = {virt_, active_};
  auto I549 = make_shared<Tensor>(I549_index);
  auto tensor428 = vector<shared_ptr<Tensor>>{den2, I549};
  auto task428 = make_shared<Task428>(tensor428, pindex);
  task428->add_dep(task265);
  densityq->add_task(task428);

  vector<IndexRange> I550_index = {active_, virt_};
  auto I550 = make_shared<Tensor>(I550_index);
  auto tensor429 = vector<shared_ptr<Tensor>>{I549, Gamma16_(), I550};
  auto task429 = make_shared<Task429>(tensor429, pindex);
  task428->add_dep(task429);
  task429->add_dep(task265);
  densityq->add_task(task429);

  auto tensor430 = vector<shared_ptr<Tensor>>{I550, t2};
  auto task430 = make_shared<Task430>(tensor430, pindex);
  task429->add_dep(task430);
  task430->add_dep(task265);
  densityq->add_task(task430);

  vector<IndexRange> I552_index = {virt_, active_};
  auto I552 = make_shared<Tensor>(I552_index);
  auto tensor431 = vector<shared_ptr<Tensor>>{den2, I552};
  auto task431 = make_shared<Task431>(tensor431, pindex);
  task431->add_dep(task265);
  densityq->add_task(task431);

  vector<IndexRange> I553_index = {active_, virt_};
  auto I553 = make_shared<Tensor>(I553_index);
  auto tensor432 = vector<shared_ptr<Tensor>>{I552, Gamma16_(), I553};
  auto task432 = make_shared<Task432>(tensor432, pindex);
  task431->add_dep(task432);
  task432->add_dep(task265);
  densityq->add_task(task432);

  auto tensor433 = vector<shared_ptr<Tensor>>{I553, t2};
  auto task433 = make_shared<Task433>(tensor433, pindex);
  task432->add_dep(task433);
  task433->add_dep(task265);
  densityq->add_task(task433);

  vector<IndexRange> I558_index = {virt_, closed_};
  auto I558 = make_shared<Tensor>(I558_index);
  auto tensor434 = vector<shared_ptr<Tensor>>{den2, I558};
  auto task434 = make_shared<Task434>(tensor434, pindex);
  task434->add_dep(task265);
  densityq->add_task(task434);

  vector<IndexRange> I559_index = {virt_, closed_};
  auto I559 = make_shared<Tensor>(I559_index);
  auto tensor435 = vector<shared_ptr<Tensor>>{I558, t2, I559};
  auto task435 = make_shared<Task435>(tensor435, pindex);
  task434->add_dep(task435);
  task435->add_dep(task265);
  densityq->add_task(task435);

  vector<IndexRange> I560_index = {active_, virt_, closed_, active_};
  auto I560 = make_shared<Tensor>(I560_index);
  auto tensor436 = vector<shared_ptr<Tensor>>{I559, Gamma38_(), I560};
  auto task436 = make_shared<Task436>(tensor436, pindex);
  task435->add_dep(task436);
  task436->add_dep(task265);
  densityq->add_task(task436);

  auto tensor437 = vector<shared_ptr<Tensor>>{I560, t2};
  auto task437 = make_shared<Task437>(tensor437, pindex);
  task436->add_dep(task437);
  task437->add_dep(task265);
  densityq->add_task(task437);

  vector<IndexRange> I567_index = {active_, active_};
  auto I567 = make_shared<Tensor>(I567_index);
  auto tensor438 = vector<shared_ptr<Tensor>>{den2, I567};
  auto task438 = make_shared<Task438>(tensor438, pindex);
  task438->add_dep(task265);
  densityq->add_task(task438);

  vector<IndexRange> I568_index;
  auto I568 = make_shared<Tensor>(I568_index);
  auto tensor439 = vector<shared_ptr<Tensor>>{I567, Gamma38_(), I568};
  auto task439 = make_shared<Task439>(tensor439, pindex);
  task438->add_dep(task439);
  task439->add_dep(task265);
  densityq->add_task(task439);

  vector<IndexRange> I569_index = {closed_, virt_, closed_, virt_};
  auto I569 = make_shared<Tensor>(I569_index);
  auto tensor440 = vector<shared_ptr<Tensor>>{I568, t2, I569};
  auto task440 = make_shared<Task440>(tensor440, pindex);
  task439->add_dep(task440);
  task440->add_dep(task265);
  densityq->add_task(task440);

  auto tensor441 = vector<shared_ptr<Tensor>>{I569, t2};
  auto task441 = make_shared<Task441>(tensor441, pindex);
  task440->add_dep(task441);
  task441->add_dep(task265);
  densityq->add_task(task441);

  shared_ptr<Tensor> I573;
  if (diagonal) {
    vector<IndexRange> I573_index = {closed_, closed_};
    I573 = make_shared<Tensor>(I573_index);
  }
  shared_ptr<Task442> task442;
  if (diagonal) {
    auto tensor442 = vector<shared_ptr<Tensor>>{den2, I573};
    task442 = make_shared<Task442>(tensor442, pindex);
    task442->add_dep(task265);
    densityq->add_task(task442);
  }

  shared_ptr<Tensor> I574;
  if (diagonal) {
    vector<IndexRange> I574_index = {closed_, virt_, closed_, virt_};
    I574 = make_shared<Tensor>(I574_index);
  }
  shared_ptr<Task443> task443;
  if (diagonal) {
    auto tensor443 = vector<shared_ptr<Tensor>>{I573, t2, I574};
    task443 = make_shared<Task443>(tensor443, pindex);
    task442->add_dep(task443);
    task443->add_dep(task265);
    densityq->add_task(task443);
  }

  shared_ptr<Task444> task444;
  if (diagonal) {
    auto tensor444 = vector<shared_ptr<Tensor>>{I574, t2};
    task444 = make_shared<Task444>(tensor444, pindex);
    task443->add_dep(task444);
    task444->add_dep(task265);
    densityq->add_task(task444);
  }

  shared_ptr<Tensor> I577;
  if (diagonal) {
    vector<IndexRange> I577_index = {virt_, virt_};
    I577 = make_shared<Tensor>(I577_index);
  }
  shared_ptr<Task445> task445;
  if (diagonal) {
    auto tensor445 = vector<shared_ptr<Tensor>>{den2, I577};
    task445 = make_shared<Task445>(tensor445, pindex);
    task445->add_dep(task265);
    densityq->add_task(task445);
  }

  shared_ptr<Tensor> I578;
  if (diagonal) {
    vector<IndexRange> I578_index = {closed_, virt_, closed_, virt_};
    I578 = make_shared<Tensor>(I578_index);
  }
  shared_ptr<Task446> task446;
  if (diagonal) {
    auto tensor446 = vector<shared_ptr<Tensor>>{I577, t2, I578};
    task446 = make_shared<Task446>(tensor446, pindex);
    task445->add_dep(task446);
    task446->add_dep(task265);
    densityq->add_task(task446);
  }

  shared_ptr<Task447> task447;
  if (diagonal) {
    auto tensor447 = vector<shared_ptr<Tensor>>{I578, t2};
    task447 = make_shared<Task447>(tensor447, pindex);
    task446->add_dep(task447);
    task447->add_dep(task265);
    densityq->add_task(task447);
  }

  vector<IndexRange> I581_index = {closed_, active_};
  auto I581 = make_shared<Tensor>(I581_index);
  auto tensor448 = vector<shared_ptr<Tensor>>{den2, I581};
  auto task448 = make_shared<Task448>(tensor448, pindex);
  task448->add_dep(task265);
  densityq->add_task(task448);

  vector<IndexRange> I582_index = {active_, closed_};
  auto I582 = make_shared<Tensor>(I582_index);
  auto tensor449 = vector<shared_ptr<Tensor>>{I581, Gamma38_(), I582};
  auto task449 = make_shared<Task449>(tensor449, pindex);
  task448->add_dep(task449);
  task449->add_dep(task265);
  densityq->add_task(task449);

  vector<IndexRange> I583_index = {active_, virt_, closed_, virt_};
  auto I583 = make_shared<Tensor>(I583_index);
  auto tensor450 = vector<shared_ptr<Tensor>>{I582, t2, I583};
  auto task450 = make_shared<Task450>(tensor450, pindex);
  task449->add_dep(task450);
  task450->add_dep(task265);
  densityq->add_task(task450);

  auto tensor451 = vector<shared_ptr<Tensor>>{I583, t2};
  auto task451 = make_shared<Task451>(tensor451, pindex);
  task450->add_dep(task451);
  task451->add_dep(task265);
  densityq->add_task(task451);

  vector<IndexRange> I587_index = {active_, virt_};
  auto I587 = make_shared<Tensor>(I587_index);
  auto tensor452 = vector<shared_ptr<Tensor>>{den2, I587};
  auto task452 = make_shared<Task452>(tensor452, pindex);
  task452->add_dep(task265);
  densityq->add_task(task452);

  vector<IndexRange> I588_index = {virt_, closed_, active_, active_};
  auto I588 = make_shared<Tensor>(I588_index);
  auto tensor453 = vector<shared_ptr<Tensor>>{I587, t2, I588};
  auto task453 = make_shared<Task453>(tensor453, pindex);
  task452->add_dep(task453);
  task453->add_dep(task265);
  densityq->add_task(task453);

  vector<IndexRange> I589_index = {active_, virt_, closed_, active_};
  auto I589 = make_shared<Tensor>(I589_index);
  auto tensor454 = vector<shared_ptr<Tensor>>{I588, Gamma35_(), I589};
  auto task454 = make_shared<Task454>(tensor454, pindex);
  task453->add_dep(task454);
  task454->add_dep(task265);
  densityq->add_task(task454);

  auto tensor455 = vector<shared_ptr<Tensor>>{I589, t2};
  auto task455 = make_shared<Task455>(tensor455, pindex);
  task454->add_dep(task455);
  task455->add_dep(task265);
  densityq->add_task(task455);

  vector<IndexRange> I590_index = {active_, virt_};
  auto I590 = make_shared<Tensor>(I590_index);
  auto tensor456 = vector<shared_ptr<Tensor>>{den2, I590};
  auto task456 = make_shared<Task456>(tensor456, pindex);
  task456->add_dep(task265);
  densityq->add_task(task456);

  vector<IndexRange> I591_index = {virt_, closed_, active_, active_};
  auto I591 = make_shared<Tensor>(I591_index);
  auto tensor457 = vector<shared_ptr<Tensor>>{I590, t2, I591};
  auto task457 = make_shared<Task457>(tensor457, pindex);
  task456->add_dep(task457);
  task457->add_dep(task265);
  densityq->add_task(task457);

  auto tensor458 = vector<shared_ptr<Tensor>>{I591, Gamma32_(), t2};
  auto task458 = make_shared<Task458>(tensor458, pindex);
  task457->add_dep(task458);
  task458->add_dep(task265);
  densityq->add_task(task458);

  auto tensor459 = vector<shared_ptr<Tensor>>{I591, Gamma35_(), t2};
  auto task459 = make_shared<Task459>(tensor459, pindex);
  task457->add_dep(task459);
  task459->add_dep(task265);
  densityq->add_task(task459);

  vector<IndexRange> I599_index = {closed_, virt_};
  auto I599 = make_shared<Tensor>(I599_index);
  auto tensor460 = vector<shared_ptr<Tensor>>{den2, I599};
  auto task460 = make_shared<Task460>(tensor460, pindex);
  task460->add_dep(task265);
  densityq->add_task(task460);

  vector<IndexRange> I600_index = {virt_, active_};
  auto I600 = make_shared<Tensor>(I600_index);
  auto tensor461 = vector<shared_ptr<Tensor>>{I599, t2, I600};
  auto task461 = make_shared<Task461>(tensor461, pindex);
  task460->add_dep(task461);
  task461->add_dep(task265);
  densityq->add_task(task461);

  auto tensor462 = vector<shared_ptr<Tensor>>{I600, Gamma60_(), t2};
  auto task462 = make_shared<Task462>(tensor462, pindex);
  task461->add_dep(task462);
  task462->add_dep(task265);
  densityq->add_task(task462);

  vector<IndexRange> I602_index = {virt_, closed_};
  auto I602 = make_shared<Tensor>(I602_index);
  auto tensor463 = vector<shared_ptr<Tensor>>{den2, I602};
  auto task463 = make_shared<Task463>(tensor463, pindex);
  task463->add_dep(task265);
  densityq->add_task(task463);

  vector<IndexRange> I603_index = {virt_, active_};
  auto I603 = make_shared<Tensor>(I603_index);
  auto tensor464 = vector<shared_ptr<Tensor>>{I602, t2, I603};
  auto task464 = make_shared<Task464>(tensor464, pindex);
  task463->add_dep(task464);
  task464->add_dep(task265);
  densityq->add_task(task464);

  auto tensor465 = vector<shared_ptr<Tensor>>{I603, Gamma60_(), t2};
  auto task465 = make_shared<Task465>(tensor465, pindex);
  task464->add_dep(task465);
  task465->add_dep(task265);
  densityq->add_task(task465);

  vector<IndexRange> I605_index = {closed_, active_};
  auto I605 = make_shared<Tensor>(I605_index);
  auto tensor466 = vector<shared_ptr<Tensor>>{den2, I605};
  auto task466 = make_shared<Task466>(tensor466, pindex);
  task466->add_dep(task265);
  densityq->add_task(task466);

  vector<IndexRange> I606_index = {closed_, active_};
  auto I606 = make_shared<Tensor>(I606_index);
  auto tensor467 = vector<shared_ptr<Tensor>>{I605, Gamma38_(), I606};
  auto task467 = make_shared<Task467>(tensor467, pindex);
  task466->add_dep(task467);
  task467->add_dep(task265);
  densityq->add_task(task467);

  vector<IndexRange> I607_index = {closed_, virt_, closed_, virt_};
  auto I607 = make_shared<Tensor>(I607_index);
  auto tensor468 = vector<shared_ptr<Tensor>>{I606, t2, I607};
  auto task468 = make_shared<Task468>(tensor468, pindex);
  task467->add_dep(task468);
  task468->add_dep(task265);
  densityq->add_task(task468);

  auto tensor469 = vector<shared_ptr<Tensor>>{I607, t2};
  auto task469 = make_shared<Task469>(tensor469, pindex);
  task468->add_dep(task469);
  task469->add_dep(task265);
  densityq->add_task(task469);

  vector<IndexRange> I617_index = {closed_, closed_};
  auto I617 = make_shared<Tensor>(I617_index);
  auto tensor470 = vector<shared_ptr<Tensor>>{den2, I617};
  auto task470 = make_shared<Task470>(tensor470, pindex);
  task470->add_dep(task265);
  densityq->add_task(task470);

  vector<IndexRange> I618_index = {virt_, closed_, virt_, active_};
  auto I618 = make_shared<Tensor>(I618_index);
  auto tensor471 = vector<shared_ptr<Tensor>>{I617, t2, I618};
  auto task471 = make_shared<Task471>(tensor471, pindex);
  task470->add_dep(task471);
  task471->add_dep(task265);
  densityq->add_task(task471);

  auto tensor472 = vector<shared_ptr<Tensor>>{I618, Gamma38_(), t2};
  auto task472 = make_shared<Task472>(tensor472, pindex);
  task471->add_dep(task472);
  task472->add_dep(task265);
  densityq->add_task(task472);

  vector<IndexRange> I621_index = {virt_, closed_, virt_, active_};
  auto I621 = make_shared<Tensor>(I621_index);
  auto tensor473 = vector<shared_ptr<Tensor>>{I617, t2, I621};
  auto task473 = make_shared<Task473>(tensor473, pindex);
  task470->add_dep(task473);
  task473->add_dep(task265);
  densityq->add_task(task473);

  auto tensor474 = vector<shared_ptr<Tensor>>{I621, Gamma38_(), t2};
  auto task474 = make_shared<Task474>(tensor474, pindex);
  task473->add_dep(task474);
  task474->add_dep(task265);
  densityq->add_task(task474);

  vector<IndexRange> I623_index = {virt_, virt_};
  auto I623 = make_shared<Tensor>(I623_index);
  auto tensor475 = vector<shared_ptr<Tensor>>{den2, I623};
  auto task475 = make_shared<Task475>(tensor475, pindex);
  task475->add_dep(task265);
  densityq->add_task(task475);

  vector<IndexRange> I624_index = {virt_, closed_, virt_, active_};
  auto I624 = make_shared<Tensor>(I624_index);
  auto tensor476 = vector<shared_ptr<Tensor>>{I623, t2, I624};
  auto task476 = make_shared<Task476>(tensor476, pindex);
  task475->add_dep(task476);
  task476->add_dep(task265);
  densityq->add_task(task476);

  auto tensor477 = vector<shared_ptr<Tensor>>{I624, Gamma38_(), t2};
  auto task477 = make_shared<Task477>(tensor477, pindex);
  task476->add_dep(task477);
  task477->add_dep(task265);
  densityq->add_task(task477);

  vector<IndexRange> I627_index = {virt_, closed_, virt_, active_};
  auto I627 = make_shared<Tensor>(I627_index);
  auto tensor478 = vector<shared_ptr<Tensor>>{I623, t2, I627};
  auto task478 = make_shared<Task478>(tensor478, pindex);
  task475->add_dep(task478);
  task478->add_dep(task265);
  densityq->add_task(task478);

  auto tensor479 = vector<shared_ptr<Tensor>>{I627, Gamma38_(), t2};
  auto task479 = make_shared<Task479>(tensor479, pindex);
  task478->add_dep(task479);
  task479->add_dep(task265);
  densityq->add_task(task479);

  vector<IndexRange> I629_index = {virt_, virt_};
  auto I629 = make_shared<Tensor>(I629_index);
  auto tensor480 = vector<shared_ptr<Tensor>>{den2, I629};
  auto task480 = make_shared<Task480>(tensor480, pindex);
  task480->add_dep(task265);
  densityq->add_task(task480);

  vector<IndexRange> I630_index = {virt_, closed_, virt_, active_};
  auto I630 = make_shared<Tensor>(I630_index);
  auto tensor481 = vector<shared_ptr<Tensor>>{I629, t2, I630};
  auto task481 = make_shared<Task481>(tensor481, pindex);
  task480->add_dep(task481);
  task481->add_dep(task265);
  densityq->add_task(task481);

  auto tensor482 = vector<shared_ptr<Tensor>>{I630, Gamma38_(), t2};
  auto task482 = make_shared<Task482>(tensor482, pindex);
  task481->add_dep(task482);
  task482->add_dep(task265);
  densityq->add_task(task482);

  vector<IndexRange> I633_index = {virt_, closed_, virt_, active_};
  auto I633 = make_shared<Tensor>(I633_index);
  auto tensor483 = vector<shared_ptr<Tensor>>{I629, t2, I633};
  auto task483 = make_shared<Task483>(tensor483, pindex);
  task480->add_dep(task483);
  task483->add_dep(task265);
  densityq->add_task(task483);

  auto tensor484 = vector<shared_ptr<Tensor>>{I633, Gamma38_(), t2};
  auto task484 = make_shared<Task484>(tensor484, pindex);
  task483->add_dep(task484);
  task484->add_dep(task265);
  densityq->add_task(task484);

  vector<IndexRange> I635_index = {active_, closed_};
  auto I635 = make_shared<Tensor>(I635_index);
  auto tensor485 = vector<shared_ptr<Tensor>>{den2, I635};
  auto task485 = make_shared<Task485>(tensor485, pindex);
  task485->add_dep(task265);
  densityq->add_task(task485);

  vector<IndexRange> I636_index = {virt_, virt_, active_, active_};
  auto I636 = make_shared<Tensor>(I636_index);
  auto tensor486 = vector<shared_ptr<Tensor>>{I635, t2, I636};
  auto task486 = make_shared<Task486>(tensor486, pindex);
  task485->add_dep(task486);
  task486->add_dep(task265);
  densityq->add_task(task486);

  auto tensor487 = vector<shared_ptr<Tensor>>{I636, Gamma60_(), t2};
  auto task487 = make_shared<Task487>(tensor487, pindex);
  task486->add_dep(task487);
  task487->add_dep(task265);
  densityq->add_task(task487);

  return densityq;
}


#endif
