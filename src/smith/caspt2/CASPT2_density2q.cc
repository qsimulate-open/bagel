//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_density2qq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_density2q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density2q = make_shared<Queue>();
  auto tensor496 = vector<shared_ptr<Tensor>>{Den1};
  auto task496 = make_shared<Task496>(tensor496, reset);
  density2q->add_task(task496);

  vector<IndexRange> I658_index = {closed_, closed_, active_, active_};
  auto I658 = make_shared<Tensor>(I658_index);
  auto tensor497 = vector<shared_ptr<Tensor>>{Den1, I658};
  auto task497 = make_shared<Task497>(tensor497, pindex);
  task497->add_dep(task496);
  density2q->add_task(task497);

  auto tensor498 = vector<shared_ptr<Tensor>>{I658, Gamma92_(), t2};
  auto task498 = make_shared<Task498>(tensor498, pindex);
  task497->add_dep(task498);
  task498->add_dep(task496);
  density2q->add_task(task498);

  vector<IndexRange> I660_index = {closed_, active_, active_, active_};
  auto I660 = make_shared<Tensor>(I660_index);
  auto tensor499 = vector<shared_ptr<Tensor>>{Den1, I660};
  auto task499 = make_shared<Task499>(tensor499, pindex);
  task499->add_dep(task496);
  density2q->add_task(task499);

  auto tensor500 = vector<shared_ptr<Tensor>>{I660, Gamma6_(), t2};
  auto task500 = make_shared<Task500>(tensor500, pindex);
  task499->add_dep(task500);
  task500->add_dep(task496);
  density2q->add_task(task500);

  vector<IndexRange> I662_index = {closed_, virt_, closed_, active_};
  auto I662 = make_shared<Tensor>(I662_index);
  auto tensor501 = vector<shared_ptr<Tensor>>{Den1, I662};
  auto task501 = make_shared<Task501>(tensor501, pindex);
  task501->add_dep(task496);
  density2q->add_task(task501);

  vector<IndexRange> I663_index = {closed_, virt_, closed_, active_};
  auto I663 = make_shared<Tensor>(I663_index);
  auto tensor502 = vector<shared_ptr<Tensor>>{I662, Gamma16_(), I663};
  auto task502 = make_shared<Task502>(tensor502, pindex);
  task501->add_dep(task502);
  task502->add_dep(task496);
  density2q->add_task(task502);

  auto tensor503 = vector<shared_ptr<Tensor>>{I663, t2};
  auto task503 = make_shared<Task503>(tensor503, pindex);
  task502->add_dep(task503);
  task503->add_dep(task496);
  density2q->add_task(task503);

  vector<IndexRange> I666_index = {virt_, closed_, active_, active_};
  auto I666 = make_shared<Tensor>(I666_index);
  auto tensor504 = vector<shared_ptr<Tensor>>{Den1, I666};
  auto task504 = make_shared<Task504>(tensor504, pindex);
  task504->add_dep(task496);
  density2q->add_task(task504);

  auto tensor505 = vector<shared_ptr<Tensor>>{I666, Gamma32_(), t2};
  auto task505 = make_shared<Task505>(tensor505, pindex);
  task504->add_dep(task505);
  task505->add_dep(task496);
  density2q->add_task(task505);

  auto tensor506 = vector<shared_ptr<Tensor>>{I666, Gamma35_(), t2};
  auto task506 = make_shared<Task506>(tensor506, pindex);
  task504->add_dep(task506);
  task506->add_dep(task496);
  density2q->add_task(task506);

  vector<IndexRange> I670_index = {virt_, closed_, active_, active_};
  auto I670 = make_shared<Tensor>(I670_index);
  auto tensor507 = vector<shared_ptr<Tensor>>{Den1, I670};
  auto task507 = make_shared<Task507>(tensor507, pindex);
  task507->add_dep(task496);
  density2q->add_task(task507);

  vector<IndexRange> I671_index = {active_, virt_, closed_, active_};
  auto I671 = make_shared<Tensor>(I671_index);
  auto tensor508 = vector<shared_ptr<Tensor>>{I670, Gamma35_(), I671};
  auto task508 = make_shared<Task508>(tensor508, pindex);
  task507->add_dep(task508);
  task508->add_dep(task496);
  density2q->add_task(task508);

  auto tensor509 = vector<shared_ptr<Tensor>>{I671, t2};
  auto task509 = make_shared<Task509>(tensor509, pindex);
  task508->add_dep(task509);
  task509->add_dep(task496);
  density2q->add_task(task509);

  vector<IndexRange> I674_index = {virt_, active_, active_, active_};
  auto I674 = make_shared<Tensor>(I674_index);
  auto tensor510 = vector<shared_ptr<Tensor>>{Den1, I674};
  auto task510 = make_shared<Task510>(tensor510, pindex);
  task510->add_dep(task496);
  density2q->add_task(task510);

  auto tensor511 = vector<shared_ptr<Tensor>>{I674, Gamma59_(), t2};
  auto task511 = make_shared<Task511>(tensor511, pindex);
  task510->add_dep(task511);
  task511->add_dep(task496);
  density2q->add_task(task511);

  shared_ptr<Tensor> I676;
  if (diagonal) {
    vector<IndexRange> I676_index = {closed_, virt_, closed_, virt_};
    I676 = make_shared<Tensor>(I676_index);
  }
  shared_ptr<Task512> task512;
  if (diagonal) {
    auto tensor512 = vector<shared_ptr<Tensor>>{Den1, I676};
    task512 = make_shared<Task512>(tensor512, pindex);
    task512->add_dep(task496);
    density2q->add_task(task512);
  }

  shared_ptr<Task513> task513;
  if (diagonal) {
    auto tensor513 = vector<shared_ptr<Tensor>>{I676, t2};
    task513 = make_shared<Task513>(tensor513, pindex);
    task512->add_dep(task513);
    task513->add_dep(task496);
    density2q->add_task(task513);
  }

  vector<IndexRange> I678_index = {virt_, closed_, virt_, active_};
  auto I678 = make_shared<Tensor>(I678_index);
  auto tensor514 = vector<shared_ptr<Tensor>>{Den1, I678};
  auto task514 = make_shared<Task514>(tensor514, pindex);
  task514->add_dep(task496);
  density2q->add_task(task514);

  vector<IndexRange> I679_index = {active_, virt_, closed_, virt_};
  auto I679 = make_shared<Tensor>(I679_index);
  auto tensor515 = vector<shared_ptr<Tensor>>{I678, Gamma38_(), I679};
  auto task515 = make_shared<Task515>(tensor515, pindex);
  task514->add_dep(task515);
  task515->add_dep(task496);
  density2q->add_task(task515);

  auto tensor516 = vector<shared_ptr<Tensor>>{I679, t2};
  auto task516 = make_shared<Task516>(tensor516, pindex);
  task515->add_dep(task516);
  task516->add_dep(task496);
  density2q->add_task(task516);

  vector<IndexRange> I682_index = {virt_, virt_, active_, active_};
  auto I682 = make_shared<Tensor>(I682_index);
  auto tensor517 = vector<shared_ptr<Tensor>>{Den1, I682};
  auto task517 = make_shared<Task517>(tensor517, pindex);
  task517->add_dep(task496);
  density2q->add_task(task517);

  auto tensor518 = vector<shared_ptr<Tensor>>{I682, Gamma60_(), t2};
  auto task518 = make_shared<Task518>(tensor518, pindex);
  task517->add_dep(task518);
  task518->add_dep(task496);
  density2q->add_task(task518);

  return density2q;
}


#endif
