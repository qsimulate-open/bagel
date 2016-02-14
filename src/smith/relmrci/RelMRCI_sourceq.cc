//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_sourceqq.cc
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


#include <src/smith/relmrci/RelMRCI.h>
#include <src/smith/relmrci/RelMRCI_tasks15.h>
#include <src/smith/relmrci/RelMRCI_tasks16.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelMRCI::RelMRCI::make_sourceq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  auto tensor727 = vector<shared_ptr<Tensor>>{s};
  auto task727 = make_shared<Task727>(tensor727, reset);
  sourceq->add_task(task727);

  vector<IndexRange> I1308_index = {closed_, active_, active_, active_};
  auto I1308 = make_shared<Tensor>(I1308_index);
  auto tensor728 = vector<shared_ptr<Tensor>>{s, I1308};
  auto task728 = make_shared<Task728>(tensor728, pindex);
  task728->add_dep(task727);
  sourceq->add_task(task728);

  auto tensor729 = vector<shared_ptr<Tensor>>{I1308, Gamma5_(), h1_};
  auto task729 = make_shared<Task729>(tensor729, pindex);
  task728->add_dep(task729);
  task729->add_dep(task727);
  sourceq->add_task(task729);

  auto tensor730 = vector<shared_ptr<Tensor>>{I1308, v2_, Gamma81_()};
  auto task730 = make_shared<Task730>(tensor730, pindex);
  task728->add_dep(task730);
  task730->add_dep(task727);
  sourceq->add_task(task730);

  auto tensor731 = vector<shared_ptr<Tensor>>{I1308, v2_, Gamma4_()};
  auto task731 = make_shared<Task731>(tensor731, pindex);
  task728->add_dep(task731);
  task731->add_dep(task727);
  sourceq->add_task(task731);

  vector<IndexRange> I1310_index = {active_, active_, closed_, virt_};
  auto I1310 = make_shared<Tensor>(I1310_index);
  auto tensor732 = vector<shared_ptr<Tensor>>{s, I1310};
  auto task732 = make_shared<Task732>(tensor732, pindex);
  task732->add_dep(task727);
  sourceq->add_task(task732);

  auto tensor733 = vector<shared_ptr<Tensor>>{I1310, h1_, Gamma27_()};
  auto task733 = make_shared<Task733>(tensor733, pindex);
  task732->add_dep(task733);
  task733->add_dep(task727);
  sourceq->add_task(task733);

  auto tensor734 = vector<shared_ptr<Tensor>>{I1310, Gamma24_(), v2_};
  auto task734 = make_shared<Task734>(tensor734, pindex);
  task732->add_dep(task734);
  task734->add_dep(task727);
  sourceq->add_task(task734);

  auto tensor735 = vector<shared_ptr<Tensor>>{I1310, v2_, Gamma5_()};
  auto task735 = make_shared<Task735>(tensor735, pindex);
  task732->add_dep(task735);
  task735->add_dep(task727);
  sourceq->add_task(task735);

  auto tensor736 = vector<shared_ptr<Tensor>>{I1310, v2_, Gamma24_()};
  auto task736 = make_shared<Task736>(tensor736, pindex);
  task732->add_dep(task736);
  task736->add_dep(task727);
  sourceq->add_task(task736);

  auto tensor737 = vector<shared_ptr<Tensor>>{I1310, v2_, Gamma24_()};
  auto task737 = make_shared<Task737>(tensor737, pindex);
  task732->add_dep(task737);
  task737->add_dep(task727);
  sourceq->add_task(task737);

  vector<IndexRange> I1312_index = {active_, active_, active_, virt_};
  auto I1312 = make_shared<Tensor>(I1312_index);
  auto tensor738 = vector<shared_ptr<Tensor>>{s, I1312};
  auto task738 = make_shared<Task738>(tensor738, pindex);
  task738->add_dep(task727);
  sourceq->add_task(task738);

  auto tensor739 = vector<shared_ptr<Tensor>>{I1312, h1_, Gamma33_()};
  auto task739 = make_shared<Task739>(tensor739, pindex);
  task738->add_dep(task739);
  task739->add_dep(task727);
  sourceq->add_task(task739);

  auto tensor740 = vector<shared_ptr<Tensor>>{I1312, v2_, Gamma32_()};
  auto task740 = make_shared<Task740>(tensor740, pindex);
  task738->add_dep(task740);
  task740->add_dep(task727);
  sourceq->add_task(task740);

  auto tensor741 = vector<shared_ptr<Tensor>>{I1312, Gamma31_(), v2_};
  auto task741 = make_shared<Task741>(tensor741, pindex);
  task738->add_dep(task741);
  task741->add_dep(task727);
  sourceq->add_task(task741);

  vector<IndexRange> I1314_index = {closed_, closed_, active_, active_};
  auto I1314 = make_shared<Tensor>(I1314_index);
  auto tensor742 = vector<shared_ptr<Tensor>>{s, I1314};
  auto task742 = make_shared<Task742>(tensor742, pindex);
  task742->add_dep(task727);
  sourceq->add_task(task742);

  auto tensor743 = vector<shared_ptr<Tensor>>{I1314, Gamma0_(), v2_};
  auto task743 = make_shared<Task743>(tensor743, pindex);
  task742->add_dep(task743);
  task743->add_dep(task727);
  sourceq->add_task(task743);

  vector<IndexRange> I1320_index = {active_, closed_, closed_, virt_};
  auto I1320 = make_shared<Tensor>(I1320_index);
  auto tensor744 = vector<shared_ptr<Tensor>>{s, I1320};
  auto task744 = make_shared<Task744>(tensor744, pindex);
  task744->add_dep(task727);
  sourceq->add_task(task744);

  auto tensor745 = vector<shared_ptr<Tensor>>{I1320, v2_, Gamma11_()};
  auto task745 = make_shared<Task745>(tensor745, pindex);
  task744->add_dep(task745);
  task745->add_dep(task727);
  sourceq->add_task(task745);

  auto tensor746 = vector<shared_ptr<Tensor>>{I1320, v2_, Gamma11_()};
  auto task746 = make_shared<Task746>(tensor746, pindex);
  task744->add_dep(task746);
  task746->add_dep(task727);
  sourceq->add_task(task746);

  shared_ptr<Tensor> I1336;
  if (diagonal) {
    vector<IndexRange> I1336_index = {closed_, virt_, closed_, virt_};
    I1336 = make_shared<Tensor>(I1336_index);
  }
  shared_ptr<Task747> task747;
  if (diagonal) {
    auto tensor747 = vector<shared_ptr<Tensor>>{s, I1336};
    task747 = make_shared<Task747>(tensor747, pindex);
    task747->add_dep(task727);
    sourceq->add_task(task747);
  }

  shared_ptr<Task748> task748;
  if (diagonal) {
    auto tensor748 = vector<shared_ptr<Tensor>>{I1336, v2_};
    task748 = make_shared<Task748>(tensor748, pindex);
    task747->add_dep(task748);
    task748->add_dep(task727);
    sourceq->add_task(task748);
  }

  vector<IndexRange> I1338_index = {active_, virt_, closed_, virt_};
  auto I1338 = make_shared<Tensor>(I1338_index);
  auto tensor749 = vector<shared_ptr<Tensor>>{s, I1338};
  auto task749 = make_shared<Task749>(tensor749, pindex);
  task749->add_dep(task727);
  sourceq->add_task(task749);

  auto tensor750 = vector<shared_ptr<Tensor>>{I1338, v2_, Gamma27_()};
  auto task750 = make_shared<Task750>(tensor750, pindex);
  task749->add_dep(task750);
  task750->add_dep(task727);
  sourceq->add_task(task750);

  auto tensor751 = vector<shared_ptr<Tensor>>{I1338, v2_, Gamma27_()};
  auto task751 = make_shared<Task751>(tensor751, pindex);
  task749->add_dep(task751);
  task751->add_dep(task727);
  sourceq->add_task(task751);

  vector<IndexRange> I1342_index = {virt_, virt_, active_, active_};
  auto I1342 = make_shared<Tensor>(I1342_index);
  auto tensor752 = vector<shared_ptr<Tensor>>{s, I1342};
  auto task752 = make_shared<Task752>(tensor752, pindex);
  task752->add_dep(task727);
  sourceq->add_task(task752);

  auto tensor753 = vector<shared_ptr<Tensor>>{I1342, Gamma33_(), v2_};
  auto task753 = make_shared<Task753>(tensor753, pindex);
  task752->add_dep(task753);
  task753->add_dep(task727);
  sourceq->add_task(task753);

  return sourceq;
}


#endif
