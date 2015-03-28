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

shared_ptr<Queue> CASPT2::CASPT2::make_corrq(const bool reset) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto corrq = make_shared<Queue>();
  vector<IndexRange> I404_index;
  auto I404 = make_shared<Tensor>(I404_index);
  vector<IndexRange> I405_index = {active_, active_, active_, active_};
  auto I405 = make_shared<Tensor>(I405_index);
  vector<shared_ptr<Tensor>> tensor379 = {I404, Gamma92_(), I405};
  auto task379 = make_shared<Task379>(tensor379, pindex);
  corrq->add_task(task379);

  vector<IndexRange> I406_index = {active_, closed_, active_, closed_};
  auto I406 = make_shared<Tensor>(I406_index);
  vector<shared_ptr<Tensor>> tensor380 = {I405, t2, I406};
  auto task380 = make_shared<Task380>(tensor380, pindex);
  task379->add_dep(task380);
  corrq->add_task(task380);

  vector<shared_ptr<Tensor>> tensor381 = {I406, t2};
  auto task381 = make_shared<Task381>(tensor381, pindex);
  task380->add_dep(task381);
  corrq->add_task(task381);

  vector<IndexRange> I408_index = {closed_, active_, active_, active_};
  auto I408 = make_shared<Tensor>(I408_index);
  vector<shared_ptr<Tensor>> tensor382 = {I404, t2, I408};
  auto task382 = make_shared<Task382>(tensor382, pindex);
  task379->add_dep(task382);
  corrq->add_task(task382);

  vector<IndexRange> I409_index = {active_, closed_, active_, active_};
  auto I409 = make_shared<Tensor>(I409_index);
  vector<shared_ptr<Tensor>> tensor383 = {I408, Gamma6_(), I409};
  auto task383 = make_shared<Task383>(tensor383, pindex);
  task382->add_dep(task383);
  corrq->add_task(task383);

  vector<shared_ptr<Tensor>> tensor384 = {I409, t2};
  auto task384 = make_shared<Task384>(tensor384, pindex);
  task383->add_dep(task384);
  corrq->add_task(task384);

  vector<IndexRange> I411_index = {active_, active_};
  auto I411 = make_shared<Tensor>(I411_index);
  vector<shared_ptr<Tensor>> tensor385 = {I404, Gamma16_(), I411};
  auto task385 = make_shared<Task385>(tensor385, pindex);
  task379->add_dep(task385);
  corrq->add_task(task385);

  vector<IndexRange> I412_index = {active_, closed_, virt_, closed_};
  auto I412 = make_shared<Tensor>(I412_index);
  vector<shared_ptr<Tensor>> tensor386 = {I411, t2, I412};
  auto task386 = make_shared<Task386>(tensor386, pindex);
  task385->add_dep(task386);
  corrq->add_task(task386);

  vector<shared_ptr<Tensor>> tensor387 = {I412, t2};
  auto task387 = make_shared<Task387>(tensor387, pindex);
  task386->add_dep(task387);
  corrq->add_task(task387);

  vector<IndexRange> I415_index = {active_, closed_, virt_, closed_};
  auto I415 = make_shared<Tensor>(I415_index);
  vector<shared_ptr<Tensor>> tensor388 = {I411, t2, I415};
  auto task388 = make_shared<Task388>(tensor388, pindex);
  task385->add_dep(task388);
  corrq->add_task(task388);

  vector<shared_ptr<Tensor>> tensor389 = {I415, t2};
  auto task389 = make_shared<Task389>(tensor389, pindex);
  task388->add_dep(task389);
  corrq->add_task(task389);

  vector<IndexRange> I417_index = {active_, active_, active_, active_};
  auto I417 = make_shared<Tensor>(I417_index);
  vector<shared_ptr<Tensor>> tensor390 = {I404, Gamma32_(), I417};
  auto task390 = make_shared<Task390>(tensor390, pindex);
  task379->add_dep(task390);
  corrq->add_task(task390);

  vector<IndexRange> I418_index = {active_, closed_, virt_, active_};
  auto I418 = make_shared<Tensor>(I418_index);
  vector<shared_ptr<Tensor>> tensor391 = {I417, t2, I418};
  auto task391 = make_shared<Task391>(tensor391, pindex);
  task390->add_dep(task391);
  corrq->add_task(task391);

  vector<shared_ptr<Tensor>> tensor392 = {I418, t2};
  auto task392 = make_shared<Task392>(tensor392, pindex);
  task391->add_dep(task392);
  corrq->add_task(task392);

  vector<IndexRange> I420_index = {active_, active_, active_, active_};
  auto I420 = make_shared<Tensor>(I420_index);
  vector<shared_ptr<Tensor>> tensor393 = {I404, Gamma35_(), I420};
  auto task393 = make_shared<Task393>(tensor393, pindex);
  task379->add_dep(task393);
  corrq->add_task(task393);

  vector<IndexRange> I421_index = {active_, closed_, virt_, active_};
  auto I421 = make_shared<Tensor>(I421_index);
  vector<shared_ptr<Tensor>> tensor394 = {I420, t2, I421};
  auto task394 = make_shared<Task394>(tensor394, pindex);
  task393->add_dep(task394);
  corrq->add_task(task394);

  vector<shared_ptr<Tensor>> tensor395 = {I421, t2};
  auto task395 = make_shared<Task395>(tensor395, pindex);
  task394->add_dep(task395);
  corrq->add_task(task395);

  vector<IndexRange> I424_index = {active_, active_, virt_, closed_};
  auto I424 = make_shared<Tensor>(I424_index);
  vector<shared_ptr<Tensor>> tensor396 = {I420, t2, I424};
  auto task396 = make_shared<Task396>(tensor396, pindex);
  task393->add_dep(task396);
  corrq->add_task(task396);

  vector<shared_ptr<Tensor>> tensor397 = {I424, t2};
  auto task397 = make_shared<Task397>(tensor397, pindex);
  task396->add_dep(task397);
  corrq->add_task(task397);

  vector<IndexRange> I427_index = {active_, active_, virt_, closed_};
  auto I427 = make_shared<Tensor>(I427_index);
  vector<shared_ptr<Tensor>> tensor398 = {I420, t2, I427};
  auto task398 = make_shared<Task398>(tensor398, pindex);
  task393->add_dep(task398);
  corrq->add_task(task398);

  vector<shared_ptr<Tensor>> tensor399 = {I427, t2};
  auto task399 = make_shared<Task399>(tensor399, pindex);
  task398->add_dep(task399);
  corrq->add_task(task399);

  vector<IndexRange> I429_index = {active_, active_, active_, active_, active_, active_};
  auto I429 = make_shared<Tensor>(I429_index);
  vector<shared_ptr<Tensor>> tensor400 = {I404, Gamma59_(), I429};
  auto task400 = make_shared<Task400>(tensor400, pindex);
  task379->add_dep(task400);
  corrq->add_task(task400);

  vector<IndexRange> I430_index = {active_, active_, virt_, active_};
  auto I430 = make_shared<Tensor>(I430_index);
  vector<shared_ptr<Tensor>> tensor401 = {I429, t2, I430};
  auto task401 = make_shared<Task401>(tensor401, pindex);
  task400->add_dep(task401);
  corrq->add_task(task401);

  vector<shared_ptr<Tensor>> tensor402 = {I430, t2};
  auto task402 = make_shared<Task402>(tensor402, pindex);
  task401->add_dep(task402);
  corrq->add_task(task402);

  vector<IndexRange> I432_index = {virt_, closed_, virt_, closed_};
  auto I432 = make_shared<Tensor>(I432_index);
  vector<shared_ptr<Tensor>> tensor403 = {I404, t2, I432};
  auto task403 = make_shared<Task403>(tensor403, pindex);
  task379->add_dep(task403);
  corrq->add_task(task403);

  vector<shared_ptr<Tensor>> tensor404 = {I432, t2};
  auto task404 = make_shared<Task404>(tensor404, pindex);
  task403->add_dep(task404);
  corrq->add_task(task404);

  vector<IndexRange> I434_index = {virt_, closed_, virt_, closed_};
  auto I434 = make_shared<Tensor>(I434_index);
  vector<shared_ptr<Tensor>> tensor405 = {I404, t2, I434};
  auto task405 = make_shared<Task405>(tensor405, pindex);
  task379->add_dep(task405);
  corrq->add_task(task405);

  vector<shared_ptr<Tensor>> tensor406 = {I434, t2};
  auto task406 = make_shared<Task406>(tensor406, pindex);
  task405->add_dep(task406);
  corrq->add_task(task406);

  vector<IndexRange> I436_index = {active_, active_};
  auto I436 = make_shared<Tensor>(I436_index);
  vector<shared_ptr<Tensor>> tensor407 = {I404, Gamma38_(), I436};
  auto task407 = make_shared<Task407>(tensor407, pindex);
  task379->add_dep(task407);
  corrq->add_task(task407);

  vector<IndexRange> I437_index = {virt_, closed_, virt_, active_};
  auto I437 = make_shared<Tensor>(I437_index);
  vector<shared_ptr<Tensor>> tensor408 = {I436, t2, I437};
  auto task408 = make_shared<Task408>(tensor408, pindex);
  task407->add_dep(task408);
  corrq->add_task(task408);

  vector<shared_ptr<Tensor>> tensor409 = {I437, t2};
  auto task409 = make_shared<Task409>(tensor409, pindex);
  task408->add_dep(task409);
  corrq->add_task(task409);

  vector<IndexRange> I440_index = {virt_, closed_, virt_, active_};
  auto I440 = make_shared<Tensor>(I440_index);
  vector<shared_ptr<Tensor>> tensor410 = {I436, t2, I440};
  auto task410 = make_shared<Task410>(tensor410, pindex);
  task407->add_dep(task410);
  corrq->add_task(task410);

  vector<shared_ptr<Tensor>> tensor411 = {I440, t2};
  auto task411 = make_shared<Task411>(tensor411, pindex);
  task410->add_dep(task411);
  corrq->add_task(task411);

  vector<IndexRange> I442_index = {active_, active_, active_, active_};
  auto I442 = make_shared<Tensor>(I442_index);
  vector<shared_ptr<Tensor>> tensor412 = {I404, Gamma60_(), I442};
  auto task412 = make_shared<Task412>(tensor412, pindex);
  task379->add_dep(task412);
  corrq->add_task(task412);

  vector<IndexRange> I443_index = {virt_, active_, virt_, active_};
  auto I443 = make_shared<Tensor>(I443_index);
  vector<shared_ptr<Tensor>> tensor413 = {I442, t2, I443};
  auto task413 = make_shared<Task413>(tensor413, pindex);
  task412->add_dep(task413);
  corrq->add_task(task413);

  vector<shared_ptr<Tensor>> tensor414 = {I443, t2};
  auto task414 = make_shared<Task414>(tensor414, pindex);
  task413->add_dep(task414);
  corrq->add_task(task414);

  return corrq;
}


#endif
