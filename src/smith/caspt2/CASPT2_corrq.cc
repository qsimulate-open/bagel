//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_corrqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_corrq(const bool reset, const bool diagonal) {

  auto corrq = make_shared<Queue>();
  auto I405 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task277 = make_shared<Task277>(Gamma92_(), I405);
  corrq->add_task(task277);

  auto task278 = make_shared<Task278>(I405, t2);
  task277->add_dep(task278);
  corrq->add_task(task278);

  auto I408 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task279 = make_shared<Task279>(t2, I408);
  task277->add_dep(task279);
  corrq->add_task(task279);

  auto task280 = make_shared<Task280>(I408, Gamma6_(), t2);
  task279->add_dep(task280);
  corrq->add_task(task280);

  auto I411 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task281 = make_shared<Task281>(Gamma16_(), I411);
  task277->add_dep(task281);
  corrq->add_task(task281);

  auto task282 = make_shared<Task282>(I411, t2);
  task281->add_dep(task282);
  corrq->add_task(task282);

  auto task283 = make_shared<Task283>(I411, t2);
  task281->add_dep(task283);
  corrq->add_task(task283);

  auto I417 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task284 = make_shared<Task284>(Gamma32_(), I417);
  task277->add_dep(task284);
  corrq->add_task(task284);

  auto task285 = make_shared<Task285>(I417, t2);
  task284->add_dep(task285);
  corrq->add_task(task285);

  auto I420 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task286 = make_shared<Task286>(Gamma35_(), I420);
  task277->add_dep(task286);
  corrq->add_task(task286);

  auto task287 = make_shared<Task287>(I420, t2);
  task286->add_dep(task287);
  corrq->add_task(task287);

  auto task288 = make_shared<Task288>(I420, t2);
  task286->add_dep(task288);
  corrq->add_task(task288);

  auto task289 = make_shared<Task289>(I420, t2);
  task286->add_dep(task289);
  corrq->add_task(task289);

  auto I429 = make_shared<TATensor<double,6>>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto task290 = make_shared<Task290>(Gamma59_(), I429);
  task277->add_dep(task290);
  corrq->add_task(task290);

  auto task291 = make_shared<Task291>(I429, t2);
  task290->add_dep(task291);
  corrq->add_task(task291);

  shared_ptr<Task292> task292;
  if (diagonal) {
    task292 = make_shared<Task292>(t2);
    task277->add_dep(task292);
    corrq->add_task(task292);
  }

  shared_ptr<Task293> task293;
  if (diagonal) {
    task293 = make_shared<Task293>(t2);
    task277->add_dep(task293);
    corrq->add_task(task293);
  }

  auto I436 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{active_, active_});
  auto task294 = make_shared<Task294>(Gamma38_(), I436);
  task277->add_dep(task294);
  corrq->add_task(task294);

  auto task295 = make_shared<Task295>(I436, t2);
  task294->add_dep(task295);
  corrq->add_task(task295);

  auto task296 = make_shared<Task296>(I436, t2);
  task294->add_dep(task296);
  corrq->add_task(task296);

  auto I442 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto task297 = make_shared<Task297>(Gamma60_(), I442);
  task277->add_dep(task297);
  corrq->add_task(task297);

  auto task298 = make_shared<Task298>(I442, t2);
  task297->add_dep(task298);
  corrq->add_task(task298);

  return corrq;
}


#endif
