//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_corrqq.cc
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

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_corrq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto corrq = make_shared<Queue>();
  vector<IndexRange> I334_index;
  auto I334 = make_shared<Tensor>(I334_index);
  vector<IndexRange> I335_index = {active_, active_, active_, active_};
  auto I335 = make_shared<Tensor>(I335_index);
  auto tensor262 = vector<shared_ptr<Tensor>>{I334, Gamma92_(), I335};
  auto task262 = make_shared<Task262>(tensor262, pindex);
  corrq->add_task(task262);

  auto tensor263 = vector<shared_ptr<Tensor>>{I335, t2};
  auto task263 = make_shared<Task263>(tensor263, pindex);
  task262->add_dep(task263);
  corrq->add_task(task263);

  vector<IndexRange> I338_index = {closed_, active_, active_, active_};
  auto I338 = make_shared<Tensor>(I338_index);
  auto tensor264 = vector<shared_ptr<Tensor>>{I334, t2, I338};
  auto task264 = make_shared<Task264>(tensor264, pindex);
  task262->add_dep(task264);
  corrq->add_task(task264);

  auto tensor265 = vector<shared_ptr<Tensor>>{I338, Gamma6_(), t2};
  auto task265 = make_shared<Task265>(tensor265, pindex);
  task264->add_dep(task265);
  corrq->add_task(task265);

  vector<IndexRange> I341_index = {active_, active_};
  auto I341 = make_shared<Tensor>(I341_index);
  auto tensor266 = vector<shared_ptr<Tensor>>{I334, Gamma16_(), I341};
  auto task266 = make_shared<Task266>(tensor266, pindex);
  task262->add_dep(task266);
  corrq->add_task(task266);

  auto tensor267 = vector<shared_ptr<Tensor>>{I341, t2};
  auto task267 = make_shared<Task267>(tensor267, pindex);
  task266->add_dep(task267);
  corrq->add_task(task267);

  auto tensor268 = vector<shared_ptr<Tensor>>{I341, t2};
  auto task268 = make_shared<Task268>(tensor268, pindex);
  task266->add_dep(task268);
  corrq->add_task(task268);

  vector<IndexRange> I347_index = {active_, active_, active_, active_};
  auto I347 = make_shared<Tensor>(I347_index);
  auto tensor269 = vector<shared_ptr<Tensor>>{I334, Gamma32_(), I347};
  auto task269 = make_shared<Task269>(tensor269, pindex);
  task262->add_dep(task269);
  corrq->add_task(task269);

  auto tensor270 = vector<shared_ptr<Tensor>>{I347, t2};
  auto task270 = make_shared<Task270>(tensor270, pindex);
  task269->add_dep(task270);
  corrq->add_task(task270);

  vector<IndexRange> I350_index = {active_, active_, active_, active_};
  auto I350 = make_shared<Tensor>(I350_index);
  auto tensor271 = vector<shared_ptr<Tensor>>{I334, Gamma35_(), I350};
  auto task271 = make_shared<Task271>(tensor271, pindex);
  task262->add_dep(task271);
  corrq->add_task(task271);

  auto tensor272 = vector<shared_ptr<Tensor>>{I350, t2};
  auto task272 = make_shared<Task272>(tensor272, pindex);
  task271->add_dep(task272);
  corrq->add_task(task272);

  auto tensor273 = vector<shared_ptr<Tensor>>{I350, t2};
  auto task273 = make_shared<Task273>(tensor273, pindex);
  task271->add_dep(task273);
  corrq->add_task(task273);

  auto tensor274 = vector<shared_ptr<Tensor>>{I350, t2};
  auto task274 = make_shared<Task274>(tensor274, pindex);
  task271->add_dep(task274);
  corrq->add_task(task274);

  vector<IndexRange> I359_index = {active_, active_, active_, active_, active_, active_};
  auto I359 = make_shared<Tensor>(I359_index);
  auto tensor275 = vector<shared_ptr<Tensor>>{I334, Gamma59_(), I359};
  auto task275 = make_shared<Task275>(tensor275, pindex);
  task262->add_dep(task275);
  corrq->add_task(task275);

  auto tensor276 = vector<shared_ptr<Tensor>>{I359, t2};
  auto task276 = make_shared<Task276>(tensor276, pindex);
  task275->add_dep(task276);
  corrq->add_task(task276);

  auto tensor277 = vector<shared_ptr<Tensor>>{I334, t2};
  auto task277 = make_shared<Task277>(tensor277, pindex);
  task262->add_dep(task277);
  corrq->add_task(task277);

  auto tensor278 = vector<shared_ptr<Tensor>>{I334, t2};
  auto task278 = make_shared<Task278>(tensor278, pindex);
  task262->add_dep(task278);
  corrq->add_task(task278);

  vector<IndexRange> I366_index = {active_, active_};
  auto I366 = make_shared<Tensor>(I366_index);
  auto tensor279 = vector<shared_ptr<Tensor>>{I334, Gamma38_(), I366};
  auto task279 = make_shared<Task279>(tensor279, pindex);
  task262->add_dep(task279);
  corrq->add_task(task279);

  auto tensor280 = vector<shared_ptr<Tensor>>{I366, t2};
  auto task280 = make_shared<Task280>(tensor280, pindex);
  task279->add_dep(task280);
  corrq->add_task(task280);

  auto tensor281 = vector<shared_ptr<Tensor>>{I366, t2};
  auto task281 = make_shared<Task281>(tensor281, pindex);
  task279->add_dep(task281);
  corrq->add_task(task281);

  vector<IndexRange> I372_index = {active_, active_, active_, active_};
  auto I372 = make_shared<Tensor>(I372_index);
  auto tensor282 = vector<shared_ptr<Tensor>>{I334, Gamma60_(), I372};
  auto task282 = make_shared<Task282>(tensor282, pindex);
  task262->add_dep(task282);
  corrq->add_task(task282);

  auto tensor283 = vector<shared_ptr<Tensor>>{I372, t2};
  auto task283 = make_shared<Task283>(tensor283, pindex);
  task282->add_dep(task283);
  corrq->add_task(task283);

  return corrq;
}


#endif
