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
#include <src/smith/relmrci/RelMRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelMRCI::RelMRCI::make_sourceq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  auto tensor731 = vector<shared_ptr<Tensor>>{s};
  auto task731 = make_shared<Task731>(tensor731, reset);
  sourceq->add_task(task731);

  vector<IndexRange> I1308_index = {closed_, active_, active_, active_};
  auto I1308 = make_shared<Tensor>(I1308_index);
  auto tensor732 = vector<shared_ptr<Tensor>>{s, I1308};
  auto task732 = make_shared<Task732>(tensor732, pindex);
  task732->add_dep(task731);
  sourceq->add_task(task732);

  auto tensor733 = vector<shared_ptr<Tensor>>{I1308, Gamma5_(), h1_};
  auto task733 = make_shared<Task733>(tensor733, pindex);
  task732->add_dep(task733);
  task733->add_dep(task731);
  sourceq->add_task(task733);

  auto tensor734 = vector<shared_ptr<Tensor>>{I1308, v2_, Gamma81_()};
  auto task734 = make_shared<Task734>(tensor734, pindex);
  task732->add_dep(task734);
  task734->add_dep(task731);
  sourceq->add_task(task734);

  auto tensor735 = vector<shared_ptr<Tensor>>{I1308, Gamma4_(), v2_};
  auto task735 = make_shared<Task735>(tensor735, pindex);
  task732->add_dep(task735);
  task735->add_dep(task731);
  sourceq->add_task(task735);

  vector<IndexRange> I1310_index = {active_, active_, closed_, virt_};
  auto I1310 = make_shared<Tensor>(I1310_index);
  auto tensor736 = vector<shared_ptr<Tensor>>{s, I1310};
  auto task736 = make_shared<Task736>(tensor736, pindex);
  task736->add_dep(task731);
  sourceq->add_task(task736);

  auto tensor737 = vector<shared_ptr<Tensor>>{I1310, h1_, Gamma27_()};
  auto task737 = make_shared<Task737>(tensor737, pindex);
  task736->add_dep(task737);
  task737->add_dep(task731);
  sourceq->add_task(task737);

  vector<IndexRange> I1325_index = {closed_, virt_, active_, active_};
  auto I1325 = make_shared<Tensor>(I1325_index);
  auto tensor738 = vector<shared_ptr<Tensor>>{I1310, Gamma24_(), I1325};
  auto task738 = make_shared<Task738>(tensor738, pindex);
  task736->add_dep(task738);
  task738->add_dep(task731);
  sourceq->add_task(task738);

  auto tensor739 = vector<shared_ptr<Tensor>>{I1325, v2_};
  auto task739 = make_shared<Task739>(tensor739, pindex);
  task738->add_dep(task739);
  task739->add_dep(task731);
  sourceq->add_task(task739);

  auto tensor740 = vector<shared_ptr<Tensor>>{I1310, v2_, Gamma5_()};
  auto task740 = make_shared<Task740>(tensor740, pindex);
  task736->add_dep(task740);
  task740->add_dep(task731);
  sourceq->add_task(task740);

  auto tensor741 = vector<shared_ptr<Tensor>>{I1310, v2_, Gamma24_()};
  auto task741 = make_shared<Task741>(tensor741, pindex);
  task736->add_dep(task741);
  task741->add_dep(task731);
  sourceq->add_task(task741);

  vector<IndexRange> I1312_index = {virt_, active_, active_, active_};
  auto I1312 = make_shared<Tensor>(I1312_index);
  auto tensor742 = vector<shared_ptr<Tensor>>{s, I1312};
  auto task742 = make_shared<Task742>(tensor742, pindex);
  task742->add_dep(task731);
  sourceq->add_task(task742);

  auto tensor743 = vector<shared_ptr<Tensor>>{I1312, Gamma33_(), h1_};
  auto task743 = make_shared<Task743>(tensor743, pindex);
  task742->add_dep(task743);
  task743->add_dep(task731);
  sourceq->add_task(task743);

  auto tensor744 = vector<shared_ptr<Tensor>>{I1312, v2_, Gamma32_()};
  auto task744 = make_shared<Task744>(tensor744, pindex);
  task742->add_dep(task744);
  task744->add_dep(task731);
  sourceq->add_task(task744);

  auto tensor745 = vector<shared_ptr<Tensor>>{I1312, v2_, Gamma31_()};
  auto task745 = make_shared<Task745>(tensor745, pindex);
  task742->add_dep(task745);
  task745->add_dep(task731);
  sourceq->add_task(task745);

  vector<IndexRange> I1314_index = {active_, active_, closed_, closed_};
  auto I1314 = make_shared<Tensor>(I1314_index);
  auto tensor746 = vector<shared_ptr<Tensor>>{s, I1314};
  auto task746 = make_shared<Task746>(tensor746, pindex);
  task746->add_dep(task731);
  sourceq->add_task(task746);

  auto tensor747 = vector<shared_ptr<Tensor>>{I1314, v2_, Gamma0_()};
  auto task747 = make_shared<Task747>(tensor747, pindex);
  task746->add_dep(task747);
  task747->add_dep(task731);
  sourceq->add_task(task747);

  vector<IndexRange> I1320_index = {active_, closed_, closed_, virt_};
  auto I1320 = make_shared<Tensor>(I1320_index);
  auto tensor748 = vector<shared_ptr<Tensor>>{s, I1320};
  auto task748 = make_shared<Task748>(tensor748, pindex);
  task748->add_dep(task731);
  sourceq->add_task(task748);

  auto tensor749 = vector<shared_ptr<Tensor>>{I1320, v2_, Gamma11_()};
  auto task749 = make_shared<Task749>(tensor749, pindex);
  task748->add_dep(task749);
  task749->add_dep(task731);
  sourceq->add_task(task749);

  auto tensor750 = vector<shared_ptr<Tensor>>{I1320, v2_, Gamma11_()};
  auto task750 = make_shared<Task750>(tensor750, pindex);
  task748->add_dep(task750);
  task750->add_dep(task731);
  sourceq->add_task(task750);

  shared_ptr<Tensor> I1336;
  if (diagonal) {
    vector<IndexRange> I1336_index = {closed_, virt_, closed_, virt_};
    I1336 = make_shared<Tensor>(I1336_index);
  }
  shared_ptr<Task751> task751;
  if (diagonal) {
    auto tensor751 = vector<shared_ptr<Tensor>>{s, I1336};
    task751 = make_shared<Task751>(tensor751, pindex);
    task751->add_dep(task731);
    sourceq->add_task(task751);
  }

  shared_ptr<Task752> task752;
  if (diagonal) {
    auto tensor752 = vector<shared_ptr<Tensor>>{I1336, v2_};
    task752 = make_shared<Task752>(tensor752, pindex);
    task751->add_dep(task752);
    task752->add_dep(task731);
    sourceq->add_task(task752);
  }

  vector<IndexRange> I1338_index = {active_, virt_, closed_, virt_};
  auto I1338 = make_shared<Tensor>(I1338_index);
  auto tensor753 = vector<shared_ptr<Tensor>>{s, I1338};
  auto task753 = make_shared<Task753>(tensor753, pindex);
  task753->add_dep(task731);
  sourceq->add_task(task753);

  auto tensor754 = vector<shared_ptr<Tensor>>{I1338, v2_, Gamma27_()};
  auto task754 = make_shared<Task754>(tensor754, pindex);
  task753->add_dep(task754);
  task754->add_dep(task731);
  sourceq->add_task(task754);

  auto tensor755 = vector<shared_ptr<Tensor>>{I1338, Gamma27_(), v2_};
  auto task755 = make_shared<Task755>(tensor755, pindex);
  task753->add_dep(task755);
  task755->add_dep(task731);
  sourceq->add_task(task755);

  vector<IndexRange> I1342_index = {active_, active_, virt_, virt_};
  auto I1342 = make_shared<Tensor>(I1342_index);
  auto tensor756 = vector<shared_ptr<Tensor>>{s, I1342};
  auto task756 = make_shared<Task756>(tensor756, pindex);
  task756->add_dep(task731);
  sourceq->add_task(task756);

  auto tensor757 = vector<shared_ptr<Tensor>>{I1342, v2_, Gamma33_()};
  auto task757 = make_shared<Task757>(tensor757, pindex);
  task756->add_dep(task757);
  task757->add_dep(task731);
  sourceq->add_task(task757);

  return sourceq;
}


#endif
