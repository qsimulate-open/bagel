//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_sourceqq.cc
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


#include <src/smith/RelMRCI.h>
#include <src/smith/RelMRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelMRCI::RelMRCI::make_sourceq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  auto tensor728 = vector<shared_ptr<Tensor>>{s};
  auto task728 = make_shared<Task728>(tensor728, reset);
  sourceq->add_task(task728);

  vector<IndexRange> I1312_index = {active_, active_, active_, closed_};
  auto I1312 = make_shared<Tensor>(I1312_index);
  auto tensor729 = vector<shared_ptr<Tensor>>{s, I1312};
  auto task729 = make_shared<Task729>(tensor729, pindex);
  task729->add_dep(task728);
  sourceq->add_task(task729);

  auto tensor730 = vector<shared_ptr<Tensor>>{I1312, h1_, Gamma5_()};
  auto task730 = make_shared<Task730>(tensor730, pindex);
  task729->add_dep(task730);
  task730->add_dep(task728);
  sourceq->add_task(task730);

  auto tensor731 = vector<shared_ptr<Tensor>>{I1312, v2_, Gamma81_()};
  auto task731 = make_shared<Task731>(tensor731, pindex);
  task729->add_dep(task731);
  task731->add_dep(task728);
  sourceq->add_task(task731);

  auto tensor732 = vector<shared_ptr<Tensor>>{I1312, Gamma4_(), v2_};
  auto task732 = make_shared<Task732>(tensor732, pindex);
  task729->add_dep(task732);
  task732->add_dep(task728);
  sourceq->add_task(task732);

  vector<IndexRange> I1314_index = {active_, active_, closed_, virt_};
  auto I1314 = make_shared<Tensor>(I1314_index);
  auto tensor733 = vector<shared_ptr<Tensor>>{s, I1314};
  auto task733 = make_shared<Task733>(tensor733, pindex);
  task733->add_dep(task728);
  sourceq->add_task(task733);

  auto tensor734 = vector<shared_ptr<Tensor>>{I1314, h1_, Gamma27_()};
  auto task734 = make_shared<Task734>(tensor734, pindex);
  task733->add_dep(task734);
  task734->add_dep(task728);
  sourceq->add_task(task734);

  auto tensor735 = vector<shared_ptr<Tensor>>{I1314, Gamma24_(), v2_};
  auto task735 = make_shared<Task735>(tensor735, pindex);
  task733->add_dep(task735);
  task735->add_dep(task728);
  sourceq->add_task(task735);

  auto tensor736 = vector<shared_ptr<Tensor>>{I1314, v2_, Gamma5_()};
  auto task736 = make_shared<Task736>(tensor736, pindex);
  task733->add_dep(task736);
  task736->add_dep(task728);
  sourceq->add_task(task736);

  auto tensor737 = vector<shared_ptr<Tensor>>{I1314, v2_, Gamma24_()};
  auto task737 = make_shared<Task737>(tensor737, pindex);
  task733->add_dep(task737);
  task737->add_dep(task728);
  sourceq->add_task(task737);

  auto tensor738 = vector<shared_ptr<Tensor>>{I1314, v2_, Gamma24_()};
  auto task738 = make_shared<Task738>(tensor738, pindex);
  task733->add_dep(task738);
  task738->add_dep(task728);
  sourceq->add_task(task738);

  vector<IndexRange> I1316_index = {active_, active_, active_, virt_};
  auto I1316 = make_shared<Tensor>(I1316_index);
  auto tensor739 = vector<shared_ptr<Tensor>>{s, I1316};
  auto task739 = make_shared<Task739>(tensor739, pindex);
  task739->add_dep(task728);
  sourceq->add_task(task739);

  auto tensor740 = vector<shared_ptr<Tensor>>{I1316, h1_, Gamma33_()};
  auto task740 = make_shared<Task740>(tensor740, pindex);
  task739->add_dep(task740);
  task740->add_dep(task728);
  sourceq->add_task(task740);

  auto tensor741 = vector<shared_ptr<Tensor>>{I1316, v2_, Gamma32_()};
  auto task741 = make_shared<Task741>(tensor741, pindex);
  task739->add_dep(task741);
  task741->add_dep(task728);
  sourceq->add_task(task741);

  auto tensor742 = vector<shared_ptr<Tensor>>{I1316, v2_, Gamma31_()};
  auto task742 = make_shared<Task742>(tensor742, pindex);
  task739->add_dep(task742);
  task742->add_dep(task728);
  sourceq->add_task(task742);

  vector<IndexRange> I1318_index = {active_, active_, closed_, closed_};
  auto I1318 = make_shared<Tensor>(I1318_index);
  auto tensor743 = vector<shared_ptr<Tensor>>{s, I1318};
  auto task743 = make_shared<Task743>(tensor743, pindex);
  task743->add_dep(task728);
  sourceq->add_task(task743);

  auto tensor744 = vector<shared_ptr<Tensor>>{I1318, v2_, Gamma0_()};
  auto task744 = make_shared<Task744>(tensor744, pindex);
  task743->add_dep(task744);
  task744->add_dep(task728);
  sourceq->add_task(task744);

  vector<IndexRange> I1324_index = {closed_, closed_, virt_, active_};
  auto I1324 = make_shared<Tensor>(I1324_index);
  auto tensor745 = vector<shared_ptr<Tensor>>{s, I1324};
  auto task745 = make_shared<Task745>(tensor745, pindex);
  task745->add_dep(task728);
  sourceq->add_task(task745);

  auto tensor746 = vector<shared_ptr<Tensor>>{I1324, Gamma11_(), v2_};
  auto task746 = make_shared<Task746>(tensor746, pindex);
  task745->add_dep(task746);
  task746->add_dep(task728);
  sourceq->add_task(task746);

  auto tensor747 = vector<shared_ptr<Tensor>>{I1324, v2_, Gamma11_()};
  auto task747 = make_shared<Task747>(tensor747, pindex);
  task745->add_dep(task747);
  task747->add_dep(task728);
  sourceq->add_task(task747);

  shared_ptr<Tensor> I1340;
  if (diagonal) {
    vector<IndexRange> I1340_index = {closed_, virt_, closed_, virt_};
    I1340 = make_shared<Tensor>(I1340_index);
  }
  shared_ptr<Task748> task748;
  if (diagonal) {
    auto tensor748 = vector<shared_ptr<Tensor>>{s, I1340};
    task748 = make_shared<Task748>(tensor748, pindex);
    task748->add_dep(task728);
    sourceq->add_task(task748);
  }

  shared_ptr<Task749> task749;
  if (diagonal) {
    auto tensor749 = vector<shared_ptr<Tensor>>{I1340, v2_};
    task749 = make_shared<Task749>(tensor749, pindex);
    task748->add_dep(task749);
    task749->add_dep(task728);
    sourceq->add_task(task749);
  }

  vector<IndexRange> I1342_index = {virt_, closed_, virt_, active_};
  auto I1342 = make_shared<Tensor>(I1342_index);
  auto tensor750 = vector<shared_ptr<Tensor>>{s, I1342};
  auto task750 = make_shared<Task750>(tensor750, pindex);
  task750->add_dep(task728);
  sourceq->add_task(task750);

  auto tensor751 = vector<shared_ptr<Tensor>>{I1342, Gamma27_(), v2_};
  auto task751 = make_shared<Task751>(tensor751, pindex);
  task750->add_dep(task751);
  task751->add_dep(task728);
  sourceq->add_task(task751);

  auto tensor752 = vector<shared_ptr<Tensor>>{I1342, v2_, Gamma27_()};
  auto task752 = make_shared<Task752>(tensor752, pindex);
  task750->add_dep(task752);
  task752->add_dep(task728);
  sourceq->add_task(task752);

  vector<IndexRange> I1346_index = {active_, active_, virt_, virt_};
  auto I1346 = make_shared<Tensor>(I1346_index);
  auto tensor753 = vector<shared_ptr<Tensor>>{s, I1346};
  auto task753 = make_shared<Task753>(tensor753, pindex);
  task753->add_dep(task728);
  sourceq->add_task(task753);

  auto tensor754 = vector<shared_ptr<Tensor>>{I1346, v2_, Gamma33_()};
  auto task754 = make_shared<Task754>(tensor754, pindex);
  task753->add_dep(task754);
  task754->add_dep(task728);
  sourceq->add_task(task754);

  return sourceq;
}


#endif
