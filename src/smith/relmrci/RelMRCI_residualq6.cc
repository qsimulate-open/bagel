//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_residualq6.cc
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
#include <src/smith/relmrci/RelMRCI_tasks14.h>
#include <src/smith/relmrci/RelMRCI_tasks15.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void RelMRCI::RelMRCI::make_residualq6(shared_ptr<Queue> residualq, shared_ptr<Task> task83, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I173_index = {virt_, active_, active_, virt_};
  auto I173 = make_shared<Tensor>(I173_index);
  auto tensor675 = vector<shared_ptr<Tensor>>{r, I173};
  auto task675 = make_shared<Task675>(tensor675, pindex);
  task675->add_dep(task83);
  residualq->add_task(task675);

  vector<IndexRange> I174_index = {virt_, active_, active_, active_};
  auto I174 = make_shared<Tensor>(I174_index);
  auto tensor676 = vector<shared_ptr<Tensor>>{I173, h1_, I174};
  auto task676 = make_shared<Task676>(tensor676, pindex);
  task675->add_dep(task676);
  task676->add_dep(task83);
  residualq->add_task(task676);

  auto tensor677 = vector<shared_ptr<Tensor>>{I174, Gamma32_(), t2};
  auto task677 = make_shared<Task677>(tensor677, pindex);
  task676->add_dep(task677);
  task677->add_dep(task83);
  residualq->add_task(task677);

  vector<IndexRange> I177_index = {active_, active_, virt_, virt_};
  auto I177 = make_shared<Tensor>(I177_index);
  auto tensor678 = vector<shared_ptr<Tensor>>{I173, Gamma33_(), I177};
  auto task678 = make_shared<Task678>(tensor678, pindex);
  task675->add_dep(task678);
  task678->add_dep(task83);
  residualq->add_task(task678);

  auto tensor679 = vector<shared_ptr<Tensor>>{I177, t2, h1_};
  auto task679 = make_shared<Task679>(tensor679, pindex);
  task678->add_dep(task679);
  task679->add_dep(task83);
  residualq->add_task(task679);

  auto tensor680 = vector<shared_ptr<Tensor>>{I177, t2, h1_};
  auto task680 = make_shared<Task680>(tensor680, pindex);
  task678->add_dep(task680);
  task680->add_dep(task83);
  residualq->add_task(task680);

  auto tensor681 = vector<shared_ptr<Tensor>>{I177, t2, v2_};
  auto task681 = make_shared<Task681>(tensor681, pindex);
  task678->add_dep(task681);
  task681->add_dep(task83);
  residualq->add_task(task681);

  auto tensor682 = vector<shared_ptr<Tensor>>{I177, t2, v2_};
  auto task682 = make_shared<Task682>(tensor682, pindex);
  task678->add_dep(task682);
  task682->add_dep(task83);
  residualq->add_task(task682);

  auto tensor683 = vector<shared_ptr<Tensor>>{I177, t2, v2_};
  auto task683 = make_shared<Task683>(tensor683, pindex);
  task678->add_dep(task683);
  task683->add_dep(task83);
  residualq->add_task(task683);

  vector<IndexRange> I1217_index = {virt_, closed_, active_, active_, active_, active_};
  auto I1217 = make_shared<Tensor>(I1217_index);
  auto tensor684 = vector<shared_ptr<Tensor>>{I173, t2, I1217};
  auto task684 = make_shared<Task684>(tensor684, pindex);
  task675->add_dep(task684);
  task684->add_dep(task83);
  residualq->add_task(task684);

  auto tensor685 = vector<shared_ptr<Tensor>>{I1217, Gamma396_(), v2_};
  auto task685 = make_shared<Task685>(tensor685, pindex);
  task684->add_dep(task685);
  task685->add_dep(task83);
  residualq->add_task(task685);

  vector<IndexRange> I1220_index = {virt_, active_, active_, active_, active_, active_};
  auto I1220 = make_shared<Tensor>(I1220_index);
  auto tensor686 = vector<shared_ptr<Tensor>>{I173, t2, I1220};
  auto task686 = make_shared<Task686>(tensor686, pindex);
  task675->add_dep(task686);
  task686->add_dep(task83);
  residualq->add_task(task686);

  auto tensor687 = vector<shared_ptr<Tensor>>{I1220, Gamma397_(), v2_};
  auto task687 = make_shared<Task687>(tensor687, pindex);
  task686->add_dep(task687);
  task687->add_dep(task83);
  residualq->add_task(task687);

  auto tensor688 = vector<shared_ptr<Tensor>>{I1220, Gamma236_(), v2_};
  auto task688 = make_shared<Task688>(tensor688, pindex);
  task686->add_dep(task688);
  task688->add_dep(task83);
  residualq->add_task(task688);

  vector<IndexRange> I1226_index = {virt_, active_, active_, active_};
  auto I1226 = make_shared<Tensor>(I1226_index);
  auto tensor689 = vector<shared_ptr<Tensor>>{I173, v2_, I1226};
  auto task689 = make_shared<Task689>(tensor689, pindex);
  task675->add_dep(task689);
  task689->add_dep(task83);
  residualq->add_task(task689);

  auto tensor690 = vector<shared_ptr<Tensor>>{I1226, Gamma336_(), t2};
  auto task690 = make_shared<Task690>(tensor690, pindex);
  task689->add_dep(task690);
  task690->add_dep(task83);
  residualq->add_task(task690);

  vector<IndexRange> I1232_index = {closed_, active_, active_, active_};
  auto I1232 = make_shared<Tensor>(I1232_index);
  auto tensor691 = vector<shared_ptr<Tensor>>{I173, t2, I1232};
  auto task691 = make_shared<Task691>(tensor691, pindex);
  task675->add_dep(task691);
  task691->add_dep(task83);
  residualq->add_task(task691);

  auto tensor692 = vector<shared_ptr<Tensor>>{I1232, Gamma391_(), v2_};
  auto task692 = make_shared<Task692>(tensor692, pindex);
  task691->add_dep(task692);
  task692->add_dep(task83);
  residualq->add_task(task692);

  auto tensor693 = vector<shared_ptr<Tensor>>{I1232, Gamma32_(), v2_};
  auto task693 = make_shared<Task693>(tensor693, pindex);
  task691->add_dep(task693);
  task693->add_dep(task83);
  residualq->add_task(task693);

  vector<IndexRange> I1238_index = {active_, virt_, active_, virt_};
  auto I1238 = make_shared<Tensor>(I1238_index);
  auto tensor694 = vector<shared_ptr<Tensor>>{I173, Gamma368_(), I1238};
  auto task694 = make_shared<Task694>(tensor694, pindex);
  task675->add_dep(task694);
  task694->add_dep(task83);
  residualq->add_task(task694);

  auto tensor695 = vector<shared_ptr<Tensor>>{I1238, t2, v2_};
  auto task695 = make_shared<Task695>(tensor695, pindex);
  task694->add_dep(task695);
  task695->add_dep(task83);
  residualq->add_task(task695);

  vector<IndexRange> I1250_index = {virt_, active_, active_, active_, virt_, active_};
  auto I1250 = make_shared<Tensor>(I1250_index);
  auto tensor696 = vector<shared_ptr<Tensor>>{I173, Gamma391_(), I1250};
  auto task696 = make_shared<Task696>(tensor696, pindex);
  task675->add_dep(task696);
  task696->add_dep(task83);
  residualq->add_task(task696);

  vector<IndexRange> I1251_index = {virt_, virt_, active_, active_};
  auto I1251 = make_shared<Tensor>(I1251_index);
  auto tensor697 = vector<shared_ptr<Tensor>>{I1250, t2, I1251};
  auto task697 = make_shared<Task697>(tensor697, pindex);
  task696->add_dep(task697);
  task697->add_dep(task83);
  residualq->add_task(task697);

  auto tensor698 = vector<shared_ptr<Tensor>>{I1251, v2_};
  auto task698 = make_shared<Task698>(tensor698, pindex);
  task697->add_dep(task698);
  task698->add_dep(task83);
  residualq->add_task(task698);

  vector<IndexRange> I1253_index = {active_, active_, virt_, active_, virt_, active_};
  auto I1253 = make_shared<Tensor>(I1253_index);
  auto tensor699 = vector<shared_ptr<Tensor>>{I173, Gamma32_(), I1253};
  auto task699 = make_shared<Task699>(tensor699, pindex);
  task675->add_dep(task699);
  task699->add_dep(task83);
  residualq->add_task(task699);

  auto tensor700 = vector<shared_ptr<Tensor>>{I1253, t2, v2_};
  auto task700 = make_shared<Task700>(tensor700, pindex);
  task699->add_dep(task700);
  task700->add_dep(task83);
  residualq->add_task(task700);

  vector<IndexRange> I1256_index = {active_, virt_, active_, active_, virt_, active_};
  auto I1256 = make_shared<Tensor>(I1256_index);
  auto tensor701 = vector<shared_ptr<Tensor>>{I173, Gamma409_(), I1256};
  auto task701 = make_shared<Task701>(tensor701, pindex);
  task675->add_dep(task701);
  task701->add_dep(task83);
  residualq->add_task(task701);

  auto tensor702 = vector<shared_ptr<Tensor>>{I1256, t2, v2_};
  auto task702 = make_shared<Task702>(tensor702, pindex);
  task701->add_dep(task702);
  task702->add_dep(task83);
  residualq->add_task(task702);

  vector<IndexRange> I194_index = {closed_, closed_, active_, active_};
  auto I194 = make_shared<Tensor>(I194_index);
  auto tensor703 = vector<shared_ptr<Tensor>>{r, I194};
  auto task703 = make_shared<Task703>(tensor703, pindex);
  task703->add_dep(task83);
  residualq->add_task(task703);

  vector<IndexRange> I195_index = {closed_, closed_, active_, active_};
  auto I195 = make_shared<Tensor>(I195_index);
  auto tensor704 = vector<shared_ptr<Tensor>>{I194, Gamma2_(), I195};
  auto task704 = make_shared<Task704>(tensor704, pindex);
  task703->add_dep(task704);
  task704->add_dep(task83);
  residualq->add_task(task704);

  auto tensor705 = vector<shared_ptr<Tensor>>{I195, t2, v2_};
  auto task705 = make_shared<Task705>(tensor705, pindex);
  task704->add_dep(task705);
  task705->add_dep(task83);
  residualq->add_task(task705);

  auto tensor706 = vector<shared_ptr<Tensor>>{I195, t2, v2_};
  auto task706 = make_shared<Task706>(tensor706, pindex);
  task704->add_dep(task706);
  task706->add_dep(task83);
  residualq->add_task(task706);

  auto tensor707 = vector<shared_ptr<Tensor>>{I194, Gamma412_(), t2};
  auto task707 = make_shared<Task707>(tensor707, pindex);
  task703->add_dep(task707);
  task707->add_dep(task83);
  residualq->add_task(task707);

  auto tensor708 = vector<shared_ptr<Tensor>>{I194, Gamma413_(), t2};
  auto task708 = make_shared<Task708>(tensor708, pindex);
  task703->add_dep(task708);
  task708->add_dep(task83);
  residualq->add_task(task708);

  vector<IndexRange> I773_index = {closed_, closed_, virt_, virt_};
  auto I773 = make_shared<Tensor>(I773_index);
  auto tensor709 = vector<shared_ptr<Tensor>>{r, I773};
  auto task709 = make_shared<Task709>(tensor709, pindex);
  task709->add_dep(task83);
  residualq->add_task(task709);

  vector<IndexRange> I774_index = {closed_, closed_, active_, active_};
  auto I774 = make_shared<Tensor>(I774_index);
  auto tensor710 = vector<shared_ptr<Tensor>>{I773, v2_, I774};
  auto task710 = make_shared<Task710>(tensor710, pindex);
  task709->add_dep(task710);
  task710->add_dep(task83);
  residualq->add_task(task710);

  auto tensor711 = vector<shared_ptr<Tensor>>{I774, Gamma2_(), t2};
  auto task711 = make_shared<Task711>(tensor711, pindex);
  task710->add_dep(task711);
  task711->add_dep(task83);
  residualq->add_task(task711);

  shared_ptr<Task712> task712;
  if (diagonal) {
    auto tensor712 = vector<shared_ptr<Tensor>>{I773, t2, v2_};
    task712 = make_shared<Task712>(tensor712, pindex);
    task709->add_dep(task712);
    task712->add_dep(task83);
    residualq->add_task(task712);
  }

  shared_ptr<Task713> task713;
  if (diagonal) {
    auto tensor713 = vector<shared_ptr<Tensor>>{I773, t2, v2_};
    task713 = make_shared<Task713>(tensor713, pindex);
    task709->add_dep(task713);
    task713->add_dep(task83);
    residualq->add_task(task713);
  }

  vector<IndexRange> I977_index = {closed_, closed_, active_, active_};
  auto I977 = make_shared<Tensor>(I977_index);
  auto tensor714 = vector<shared_ptr<Tensor>>{I773, t2, I977};
  auto task714 = make_shared<Task714>(tensor714, pindex);
  task709->add_dep(task714);
  task714->add_dep(task83);
  residualq->add_task(task714);

  auto tensor715 = vector<shared_ptr<Tensor>>{I977, Gamma368_(), v2_};
  auto task715 = make_shared<Task715>(tensor715, pindex);
  task714->add_dep(task715);
  task715->add_dep(task83);
  residualq->add_task(task715);

  vector<IndexRange> I1289_index = {closed_, virt_, closed_, virt_};
  auto I1289 = make_shared<Tensor>(I1289_index);
  auto tensor716 = vector<shared_ptr<Tensor>>{I773, Gamma424_(), I1289};
  auto task716 = make_shared<Task716>(tensor716, pindex);
  task709->add_dep(task716);
  task716->add_dep(task83);
  residualq->add_task(task716);

  auto tensor717 = vector<shared_ptr<Tensor>>{I1289, t2};
  auto task717 = make_shared<Task717>(tensor717, pindex);
  task716->add_dep(task717);
  task717->add_dep(task83);
  residualq->add_task(task717);

  vector<IndexRange> I1293_index = {closed_, virt_, closed_, virt_};
  auto I1293 = make_shared<Tensor>(I1293_index);
  auto tensor718 = vector<shared_ptr<Tensor>>{I773, Gamma426_(), I1293};
  auto task718 = make_shared<Task718>(tensor718, pindex);
  task709->add_dep(task718);
  task718->add_dep(task83);
  residualq->add_task(task718);

  auto tensor719 = vector<shared_ptr<Tensor>>{I1293, t2};
  auto task719 = make_shared<Task719>(tensor719, pindex);
  task718->add_dep(task719);
  task719->add_dep(task83);
  residualq->add_task(task719);

  vector<IndexRange> I1228_index = {active_, active_, virt_, virt_};
  auto I1228 = make_shared<Tensor>(I1228_index);
  auto tensor720 = vector<shared_ptr<Tensor>>{r, I1228};
  auto task720 = make_shared<Task720>(tensor720, pindex);
  task720->add_dep(task83);
  residualq->add_task(task720);

  vector<IndexRange> I1229_index = {closed_, closed_, active_, active_};
  auto I1229 = make_shared<Tensor>(I1229_index);
  auto tensor721 = vector<shared_ptr<Tensor>>{I1228, t2, I1229};
  auto task721 = make_shared<Task721>(tensor721, pindex);
  task720->add_dep(task721);
  task721->add_dep(task83);
  residualq->add_task(task721);

  auto tensor722 = vector<shared_ptr<Tensor>>{I1229, Gamma368_(), v2_};
  auto task722 = make_shared<Task722>(tensor722, pindex);
  task721->add_dep(task722);
  task722->add_dep(task83);
  residualq->add_task(task722);

  vector<IndexRange> I1262_index = {virt_, virt_, active_, active_};
  auto I1262 = make_shared<Tensor>(I1262_index);
  auto tensor723 = vector<shared_ptr<Tensor>>{I1228, Gamma368_(), I1262};
  auto task723 = make_shared<Task723>(tensor723, pindex);
  task720->add_dep(task723);
  task723->add_dep(task83);
  residualq->add_task(task723);

  auto tensor724 = vector<shared_ptr<Tensor>>{I1262, t2, v2_};
  auto task724 = make_shared<Task724>(tensor724, pindex);
  task723->add_dep(task724);
  task724->add_dep(task83);
  residualq->add_task(task724);

  auto tensor725 = vector<shared_ptr<Tensor>>{I1228, Gamma432_(), t2};
  auto task725 = make_shared<Task725>(tensor725, pindex);
  task720->add_dep(task725);
  task725->add_dep(task83);
  residualq->add_task(task725);

  auto tensor726 = vector<shared_ptr<Tensor>>{I1228, Gamma433_(), t2};
  auto task726 = make_shared<Task726>(tensor726, pindex);
  task720->add_dep(task726);
  task726->add_dep(task83);
  residualq->add_task(task726);
}

#endif
