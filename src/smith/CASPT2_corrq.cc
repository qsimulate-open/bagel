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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_corrq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto corrq = make_shared<Queue>();
  vector<IndexRange> I418_index;
  auto I418 = make_shared<Tensor>(I418_index);
  vector<IndexRange> I419_index = {active_, active_, active_, active_};
  auto I419 = make_shared<Tensor>(I419_index);
  vector<shared_ptr<Tensor>> tensor392 = {I418, Gamma94_(), I419};
  auto task392 = make_shared<Task392>(tensor392, pindex);
  corrq->add_task(task392);

  vector<IndexRange> I420_index = {active_, closed_, active_, closed_};
  auto I420 = make_shared<Tensor>(I420_index);
  vector<shared_ptr<Tensor>> tensor393 = {I419, t2, I420};
  auto task393 = make_shared<Task393>(tensor393, pindex);
  task392->add_dep(task393);
  corrq->add_task(task393);

  vector<shared_ptr<Tensor>> tensor394 = {I420, t2};
  auto task394 = make_shared<Task394>(tensor394, pindex);
  task393->add_dep(task394);
  corrq->add_task(task394);

  vector<IndexRange> I422_index = {closed_, active_, active_, active_};
  auto I422 = make_shared<Tensor>(I422_index);
  vector<shared_ptr<Tensor>> tensor395 = {I418, t2, I422};
  auto task395 = make_shared<Task395>(tensor395, pindex);
  task392->add_dep(task395);
  corrq->add_task(task395);

  vector<IndexRange> I423_index = {active_, closed_, active_, active_};
  auto I423 = make_shared<Tensor>(I423_index);
  vector<shared_ptr<Tensor>> tensor396 = {I422, Gamma6_(), I423};
  auto task396 = make_shared<Task396>(tensor396, pindex);
  task395->add_dep(task396);
  corrq->add_task(task396);

  vector<shared_ptr<Tensor>> tensor397 = {I423, t2};
  auto task397 = make_shared<Task397>(tensor397, pindex);
  task396->add_dep(task397);
  corrq->add_task(task397);

  vector<IndexRange> I425_index = {active_, active_};
  auto I425 = make_shared<Tensor>(I425_index);
  vector<shared_ptr<Tensor>> tensor398 = {I418, Gamma16_(), I425};
  auto task398 = make_shared<Task398>(tensor398, pindex);
  task392->add_dep(task398);
  corrq->add_task(task398);

  vector<IndexRange> I426_index = {active_, closed_, virt_, closed_};
  auto I426 = make_shared<Tensor>(I426_index);
  vector<shared_ptr<Tensor>> tensor399 = {I425, t2, I426};
  auto task399 = make_shared<Task399>(tensor399, pindex);
  task398->add_dep(task399);
  corrq->add_task(task399);

  vector<shared_ptr<Tensor>> tensor400 = {I426, t2};
  auto task400 = make_shared<Task400>(tensor400, pindex);
  task399->add_dep(task400);
  corrq->add_task(task400);

  vector<IndexRange> I429_index = {active_, closed_, virt_, closed_};
  auto I429 = make_shared<Tensor>(I429_index);
  vector<shared_ptr<Tensor>> tensor401 = {I425, t2, I429};
  auto task401 = make_shared<Task401>(tensor401, pindex);
  task398->add_dep(task401);
  corrq->add_task(task401);

  vector<shared_ptr<Tensor>> tensor402 = {I429, t2};
  auto task402 = make_shared<Task402>(tensor402, pindex);
  task401->add_dep(task402);
  corrq->add_task(task402);

  vector<IndexRange> I431_index = {active_, active_, active_, active_};
  auto I431 = make_shared<Tensor>(I431_index);
  vector<shared_ptr<Tensor>> tensor403 = {I418, Gamma32_(), I431};
  auto task403 = make_shared<Task403>(tensor403, pindex);
  task392->add_dep(task403);
  corrq->add_task(task403);

  vector<IndexRange> I432_index = {active_, closed_, virt_, active_};
  auto I432 = make_shared<Tensor>(I432_index);
  vector<shared_ptr<Tensor>> tensor404 = {I431, t2, I432};
  auto task404 = make_shared<Task404>(tensor404, pindex);
  task403->add_dep(task404);
  corrq->add_task(task404);

  vector<shared_ptr<Tensor>> tensor405 = {I432, t2};
  auto task405 = make_shared<Task405>(tensor405, pindex);
  task404->add_dep(task405);
  corrq->add_task(task405);

  vector<IndexRange> I434_index = {active_, active_, active_, active_};
  auto I434 = make_shared<Tensor>(I434_index);
  vector<shared_ptr<Tensor>> tensor406 = {I418, Gamma35_(), I434};
  auto task406 = make_shared<Task406>(tensor406, pindex);
  task392->add_dep(task406);
  corrq->add_task(task406);

  vector<IndexRange> I435_index = {active_, closed_, virt_, active_};
  auto I435 = make_shared<Tensor>(I435_index);
  vector<shared_ptr<Tensor>> tensor407 = {I434, t2, I435};
  auto task407 = make_shared<Task407>(tensor407, pindex);
  task406->add_dep(task407);
  corrq->add_task(task407);

  vector<shared_ptr<Tensor>> tensor408 = {I435, t2};
  auto task408 = make_shared<Task408>(tensor408, pindex);
  task407->add_dep(task408);
  corrq->add_task(task408);

  vector<IndexRange> I438_index = {active_, active_, virt_, closed_};
  auto I438 = make_shared<Tensor>(I438_index);
  vector<shared_ptr<Tensor>> tensor409 = {I434, t2, I438};
  auto task409 = make_shared<Task409>(tensor409, pindex);
  task406->add_dep(task409);
  corrq->add_task(task409);

  vector<shared_ptr<Tensor>> tensor410 = {I438, t2};
  auto task410 = make_shared<Task410>(tensor410, pindex);
  task409->add_dep(task410);
  corrq->add_task(task410);

  vector<IndexRange> I441_index = {active_, active_, virt_, closed_};
  auto I441 = make_shared<Tensor>(I441_index);
  vector<shared_ptr<Tensor>> tensor411 = {I434, t2, I441};
  auto task411 = make_shared<Task411>(tensor411, pindex);
  task406->add_dep(task411);
  corrq->add_task(task411);

  vector<shared_ptr<Tensor>> tensor412 = {I441, t2};
  auto task412 = make_shared<Task412>(tensor412, pindex);
  task411->add_dep(task412);
  corrq->add_task(task412);

  vector<IndexRange> I443_index = {active_, active_, active_, active_, active_, active_};
  auto I443 = make_shared<Tensor>(I443_index);
  vector<shared_ptr<Tensor>> tensor413 = {I418, Gamma59_(), I443};
  auto task413 = make_shared<Task413>(tensor413, pindex);
  task392->add_dep(task413);
  corrq->add_task(task413);

  vector<IndexRange> I444_index = {active_, active_, virt_, active_};
  auto I444 = make_shared<Tensor>(I444_index);
  vector<shared_ptr<Tensor>> tensor414 = {I443, t2, I444};
  auto task414 = make_shared<Task414>(tensor414, pindex);
  task413->add_dep(task414);
  corrq->add_task(task414);

  vector<shared_ptr<Tensor>> tensor415 = {I444, t2};
  auto task415 = make_shared<Task415>(tensor415, pindex);
  task414->add_dep(task415);
  corrq->add_task(task415);

  vector<IndexRange> I446_index = {virt_, closed_, virt_, closed_};
  auto I446 = make_shared<Tensor>(I446_index);
  vector<shared_ptr<Tensor>> tensor416 = {I418, t2, I446};
  auto task416 = make_shared<Task416>(tensor416, pindex);
  task392->add_dep(task416);
  corrq->add_task(task416);

  vector<shared_ptr<Tensor>> tensor417 = {I446, t2};
  auto task417 = make_shared<Task417>(tensor417, pindex);
  task416->add_dep(task417);
  corrq->add_task(task417);

  vector<IndexRange> I448_index = {virt_, closed_, virt_, closed_};
  auto I448 = make_shared<Tensor>(I448_index);
  vector<shared_ptr<Tensor>> tensor418 = {I418, t2, I448};
  auto task418 = make_shared<Task418>(tensor418, pindex);
  task392->add_dep(task418);
  corrq->add_task(task418);

  vector<shared_ptr<Tensor>> tensor419 = {I448, t2};
  auto task419 = make_shared<Task419>(tensor419, pindex);
  task418->add_dep(task419);
  corrq->add_task(task419);

  vector<IndexRange> I450_index = {active_, active_};
  auto I450 = make_shared<Tensor>(I450_index);
  vector<shared_ptr<Tensor>> tensor420 = {I418, Gamma38_(), I450};
  auto task420 = make_shared<Task420>(tensor420, pindex);
  task392->add_dep(task420);
  corrq->add_task(task420);

  vector<IndexRange> I451_index = {virt_, closed_, virt_, active_};
  auto I451 = make_shared<Tensor>(I451_index);
  vector<shared_ptr<Tensor>> tensor421 = {I450, t2, I451};
  auto task421 = make_shared<Task421>(tensor421, pindex);
  task420->add_dep(task421);
  corrq->add_task(task421);

  vector<shared_ptr<Tensor>> tensor422 = {I451, t2};
  auto task422 = make_shared<Task422>(tensor422, pindex);
  task421->add_dep(task422);
  corrq->add_task(task422);

  vector<IndexRange> I454_index = {virt_, closed_, virt_, active_};
  auto I454 = make_shared<Tensor>(I454_index);
  vector<shared_ptr<Tensor>> tensor423 = {I450, t2, I454};
  auto task423 = make_shared<Task423>(tensor423, pindex);
  task420->add_dep(task423);
  corrq->add_task(task423);

  vector<shared_ptr<Tensor>> tensor424 = {I454, t2};
  auto task424 = make_shared<Task424>(tensor424, pindex);
  task423->add_dep(task424);
  corrq->add_task(task424);

  vector<IndexRange> I456_index = {active_, active_, active_, active_};
  auto I456 = make_shared<Tensor>(I456_index);
  vector<shared_ptr<Tensor>> tensor425 = {I418, Gamma60_(), I456};
  auto task425 = make_shared<Task425>(tensor425, pindex);
  task392->add_dep(task425);
  corrq->add_task(task425);

  vector<IndexRange> I457_index = {virt_, active_, virt_, active_};
  auto I457 = make_shared<Tensor>(I457_index);
  vector<shared_ptr<Tensor>> tensor426 = {I456, t2, I457};
  auto task426 = make_shared<Task426>(tensor426, pindex);
  task425->add_dep(task426);
  corrq->add_task(task426);

  vector<shared_ptr<Tensor>> tensor427 = {I457, t2};
  auto task427 = make_shared<Task427>(tensor427, pindex);
  task426->add_dep(task427);
  corrq->add_task(task427);

  return corrq;
}


