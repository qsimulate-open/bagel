//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_deci2qq.cc
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


#include <src/smith/caspt2/MSCASPT2.h>
#include <src/smith/caspt2/MSCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_deci2q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deci2q = make_shared<Queue>();
  auto tensor658 = vector<shared_ptr<Tensor>>{deci};
  auto task658 = make_shared<Task658>(tensor658, reset);
  deci2q->add_task(task658);

  vector<IndexRange> I1120_index = {ci_};
  auto I1120 = make_shared<Tensor>(I1120_index);
  auto tensor659 = vector<shared_ptr<Tensor>>{deci, I1120};
  auto task659 = make_shared<Task659>(tensor659, cindex);
  task659->add_dep(task658);
  deci2q->add_task(task659);

  vector<IndexRange> I1121_index = {active_, active_, active_, active_};
  auto I1121 = make_shared<Tensor>(I1121_index);
  auto tensor660 = vector<shared_ptr<Tensor>>{I1120, Gamma111_(), I1121};
  auto task660 = make_shared<Task660>(tensor660, cindex);
  task659->add_dep(task660);
  task660->add_dep(task658);
  deci2q->add_task(task660);

  auto tensor661 = vector<shared_ptr<Tensor>>{I1121, v2_, l2};
  auto task661 = make_shared<Task661>(tensor661, cindex);
  task660->add_dep(task661);
  task661->add_dep(task658);
  deci2q->add_task(task661);

  auto tensor662 = vector<shared_ptr<Tensor>>{I1121, v2_, l2};
  auto task662 = make_shared<Task662>(tensor662, cindex);
  task660->add_dep(task662);
  task662->add_dep(task658);
  deci2q->add_task(task662);

  vector<IndexRange> I1124_index = {active_, active_, active_, active_, active_, active_};
  auto I1124 = make_shared<Tensor>(I1124_index);
  auto tensor663 = vector<shared_ptr<Tensor>>{I1120, Gamma323_(), I1124};
  auto task663 = make_shared<Task663>(tensor663, cindex);
  task659->add_dep(task663);
  task663->add_dep(task658);
  deci2q->add_task(task663);

  auto tensor664 = vector<shared_ptr<Tensor>>{I1124, v2_, l2};
  auto task664 = make_shared<Task664>(tensor664, cindex);
  task663->add_dep(task664);
  task664->add_dep(task658);
  deci2q->add_task(task664);

  vector<IndexRange> I1127_index = {active_, active_, active_, active_, active_, active_};
  auto I1127 = make_shared<Tensor>(I1127_index);
  auto tensor665 = vector<shared_ptr<Tensor>>{I1120, Gamma116_(), I1127};
  auto task665 = make_shared<Task665>(tensor665, cindex);
  task659->add_dep(task665);
  task665->add_dep(task658);
  deci2q->add_task(task665);

  auto tensor666 = vector<shared_ptr<Tensor>>{I1127, v2_, l2};
  auto task666 = make_shared<Task666>(tensor666, cindex);
  task665->add_dep(task666);
  task666->add_dep(task658);
  deci2q->add_task(task666);

  auto tensor667 = vector<shared_ptr<Tensor>>{I1127, v2_, l2};
  auto task667 = make_shared<Task667>(tensor667, cindex);
  task665->add_dep(task667);
  task667->add_dep(task658);
  deci2q->add_task(task667);

  vector<IndexRange> I1130_index = {active_, active_};
  auto I1130 = make_shared<Tensor>(I1130_index);
  auto tensor668 = vector<shared_ptr<Tensor>>{I1120, Gamma126_(), I1130};
  auto task668 = make_shared<Task668>(tensor668, cindex);
  task659->add_dep(task668);
  task668->add_dep(task658);
  deci2q->add_task(task668);

  auto tensor669 = vector<shared_ptr<Tensor>>{I1130, v2_, l2};
  auto task669 = make_shared<Task669>(tensor669, cindex);
  task668->add_dep(task669);
  task669->add_dep(task658);
  deci2q->add_task(task669);

  auto tensor670 = vector<shared_ptr<Tensor>>{I1130, v2_, l2};
  auto task670 = make_shared<Task670>(tensor670, cindex);
  task668->add_dep(task670);
  task670->add_dep(task658);
  deci2q->add_task(task670);

  auto tensor671 = vector<shared_ptr<Tensor>>{I1130, v2_, l2};
  auto task671 = make_shared<Task671>(tensor671, cindex);
  task668->add_dep(task671);
  task671->add_dep(task658);
  deci2q->add_task(task671);

  auto tensor672 = vector<shared_ptr<Tensor>>{I1130, v2_, l2};
  auto task672 = make_shared<Task672>(tensor672, cindex);
  task668->add_dep(task672);
  task672->add_dep(task658);
  deci2q->add_task(task672);

  vector<IndexRange> I1136_index = {active_, active_, active_, active_};
  auto I1136 = make_shared<Tensor>(I1136_index);
  auto tensor673 = vector<shared_ptr<Tensor>>{I1120, Gamma145_(), I1136};
  auto task673 = make_shared<Task673>(tensor673, cindex);
  task659->add_dep(task673);
  task673->add_dep(task658);
  deci2q->add_task(task673);

  auto tensor674 = vector<shared_ptr<Tensor>>{I1136, v2_, l2};
  auto task674 = make_shared<Task674>(tensor674, cindex);
  task673->add_dep(task674);
  task674->add_dep(task658);
  deci2q->add_task(task674);

  auto tensor675 = vector<shared_ptr<Tensor>>{I1136, v2_, l2};
  auto task675 = make_shared<Task675>(tensor675, cindex);
  task673->add_dep(task675);
  task675->add_dep(task658);
  deci2q->add_task(task675);

  auto tensor676 = vector<shared_ptr<Tensor>>{I1136, v2_, l2};
  auto task676 = make_shared<Task676>(tensor676, cindex);
  task673->add_dep(task676);
  task676->add_dep(task658);
  deci2q->add_task(task676);

  auto tensor677 = vector<shared_ptr<Tensor>>{I1136, v2_, l2};
  auto task677 = make_shared<Task677>(tensor677, cindex);
  task673->add_dep(task677);
  task677->add_dep(task658);
  deci2q->add_task(task677);

  auto tensor678 = vector<shared_ptr<Tensor>>{I1136, v2_, l2};
  auto task678 = make_shared<Task678>(tensor678, cindex);
  task673->add_dep(task678);
  task678->add_dep(task658);
  deci2q->add_task(task678);

  auto tensor679 = vector<shared_ptr<Tensor>>{I1136, v2_, l2};
  auto task679 = make_shared<Task679>(tensor679, cindex);
  task673->add_dep(task679);
  task679->add_dep(task658);
  deci2q->add_task(task679);

  auto tensor680 = vector<shared_ptr<Tensor>>{I1136, v2_, l2};
  auto task680 = make_shared<Task680>(tensor680, cindex);
  task673->add_dep(task680);
  task680->add_dep(task658);
  deci2q->add_task(task680);

  auto tensor681 = vector<shared_ptr<Tensor>>{I1136, v2_, l2};
  auto task681 = make_shared<Task681>(tensor681, cindex);
  task673->add_dep(task681);
  task681->add_dep(task658);
  deci2q->add_task(task681);

  auto tensor682 = vector<shared_ptr<Tensor>>{I1136, v2_, l2};
  auto task682 = make_shared<Task682>(tensor682, cindex);
  task673->add_dep(task682);
  task682->add_dep(task658);
  deci2q->add_task(task682);

  auto tensor683 = vector<shared_ptr<Tensor>>{I1136, v2_, l2};
  auto task683 = make_shared<Task683>(tensor683, cindex);
  task673->add_dep(task683);
  task683->add_dep(task658);
  deci2q->add_task(task683);

  vector<IndexRange> I1139_index = {active_, active_, active_, active_};
  auto I1139 = make_shared<Tensor>(I1139_index);
  auto tensor684 = vector<shared_ptr<Tensor>>{I1120, Gamma139_(), I1139};
  auto task684 = make_shared<Task684>(tensor684, cindex);
  task659->add_dep(task684);
  task684->add_dep(task658);
  deci2q->add_task(task684);

  auto tensor685 = vector<shared_ptr<Tensor>>{I1139, v2_, l2};
  auto task685 = make_shared<Task685>(tensor685, cindex);
  task684->add_dep(task685);
  task685->add_dep(task658);
  deci2q->add_task(task685);

  vector<IndexRange> I1142_index = {active_, active_, active_, active_};
  auto I1142 = make_shared<Tensor>(I1142_index);
  auto tensor686 = vector<shared_ptr<Tensor>>{I1120, Gamma142_(), I1142};
  auto task686 = make_shared<Task686>(tensor686, cindex);
  task659->add_dep(task686);
  task686->add_dep(task658);
  deci2q->add_task(task686);

  auto tensor687 = vector<shared_ptr<Tensor>>{I1142, v2_, l2};
  auto task687 = make_shared<Task687>(tensor687, cindex);
  task686->add_dep(task687);
  task687->add_dep(task658);
  deci2q->add_task(task687);

  auto tensor688 = vector<shared_ptr<Tensor>>{I1142, v2_, l2};
  auto task688 = make_shared<Task688>(tensor688, cindex);
  task686->add_dep(task688);
  task688->add_dep(task658);
  deci2q->add_task(task688);

  vector<IndexRange> I1151_index = {active_, active_, active_, active_};
  auto I1151 = make_shared<Tensor>(I1151_index);
  auto tensor689 = vector<shared_ptr<Tensor>>{I1120, Gamma117_(), I1151};
  auto task689 = make_shared<Task689>(tensor689, cindex);
  task659->add_dep(task689);
  task689->add_dep(task658);
  deci2q->add_task(task689);

  auto tensor690 = vector<shared_ptr<Tensor>>{I1151, v2_, l2};
  auto task690 = make_shared<Task690>(tensor690, cindex);
  task689->add_dep(task690);
  task690->add_dep(task658);
  deci2q->add_task(task690);

  auto tensor691 = vector<shared_ptr<Tensor>>{I1151, h1_, l2};
  auto task691 = make_shared<Task691>(tensor691, cindex);
  task689->add_dep(task691);
  task691->add_dep(task658);
  deci2q->add_task(task691);

  vector<IndexRange> I1160_index = {active_, active_, active_, active_, active_, active_};
  auto I1160 = make_shared<Tensor>(I1160_index);
  auto tensor692 = vector<shared_ptr<Tensor>>{I1120, Gamma169_(), I1160};
  auto task692 = make_shared<Task692>(tensor692, cindex);
  task659->add_dep(task692);
  task692->add_dep(task658);
  deci2q->add_task(task692);

  auto tensor693 = vector<shared_ptr<Tensor>>{I1160, v2_, l2};
  auto task693 = make_shared<Task693>(tensor693, cindex);
  task692->add_dep(task693);
  task693->add_dep(task658);
  deci2q->add_task(task693);

  auto tensor694 = vector<shared_ptr<Tensor>>{I1160, v2_, l2};
  auto task694 = make_shared<Task694>(tensor694, cindex);
  task692->add_dep(task694);
  task694->add_dep(task658);
  deci2q->add_task(task694);

  vector<IndexRange> I1163_index = {active_, active_, active_, active_, active_, active_};
  auto I1163 = make_shared<Tensor>(I1163_index);
  auto tensor695 = vector<shared_ptr<Tensor>>{I1120, Gamma167_(), I1163};
  auto task695 = make_shared<Task695>(tensor695, cindex);
  task659->add_dep(task695);
  task695->add_dep(task658);
  deci2q->add_task(task695);

  auto tensor696 = vector<shared_ptr<Tensor>>{I1163, v2_, l2};
  auto task696 = make_shared<Task696>(tensor696, cindex);
  task695->add_dep(task696);
  task696->add_dep(task658);
  deci2q->add_task(task696);

  vector<IndexRange> I1166_index = {active_, active_};
  auto I1166 = make_shared<Tensor>(I1166_index);
  auto tensor697 = vector<shared_ptr<Tensor>>{I1120, Gamma148_(), I1166};
  auto task697 = make_shared<Task697>(tensor697, cindex);
  task659->add_dep(task697);
  task697->add_dep(task658);
  deci2q->add_task(task697);

  auto tensor698 = vector<shared_ptr<Tensor>>{I1166, v2_, l2};
  auto task698 = make_shared<Task698>(tensor698, cindex);
  task697->add_dep(task698);
  task698->add_dep(task658);
  deci2q->add_task(task698);

  auto tensor699 = vector<shared_ptr<Tensor>>{I1166, v2_, l2};
  auto task699 = make_shared<Task699>(tensor699, cindex);
  task697->add_dep(task699);
  task699->add_dep(task658);
  deci2q->add_task(task699);

  auto tensor700 = vector<shared_ptr<Tensor>>{I1166, v2_, l2};
  auto task700 = make_shared<Task700>(tensor700, cindex);
  task697->add_dep(task700);
  task700->add_dep(task658);
  deci2q->add_task(task700);

  auto tensor701 = vector<shared_ptr<Tensor>>{I1166, v2_, l2};
  auto task701 = make_shared<Task701>(tensor701, cindex);
  task697->add_dep(task701);
  task701->add_dep(task658);
  deci2q->add_task(task701);

  vector<IndexRange> I1233_index = {active_, closed_, virt_, active_};
  auto I1233 = make_shared<Tensor>(I1233_index);
  auto tensor702 = vector<shared_ptr<Tensor>>{I1166, h1_, I1233};
  auto task702 = make_shared<Task702>(tensor702, cindex);
  task697->add_dep(task702);
  task702->add_dep(task658);
  deci2q->add_task(task702);

  auto tensor703 = vector<shared_ptr<Tensor>>{I1233, l2};
  auto task703 = make_shared<Task703>(tensor703, cindex);
  task702->add_dep(task703);
  task703->add_dep(task658);
  deci2q->add_task(task703);

  vector<IndexRange> I1236_index = {active_, active_, virt_, closed_};
  auto I1236 = make_shared<Tensor>(I1236_index);
  auto tensor704 = vector<shared_ptr<Tensor>>{I1166, h1_, I1236};
  auto task704 = make_shared<Task704>(tensor704, cindex);
  task697->add_dep(task704);
  task704->add_dep(task658);
  deci2q->add_task(task704);

  auto tensor705 = vector<shared_ptr<Tensor>>{I1236, l2};
  auto task705 = make_shared<Task705>(tensor705, cindex);
  task704->add_dep(task705);
  task705->add_dep(task658);
  deci2q->add_task(task705);

  vector<IndexRange> I1172_index = {active_, active_, active_, active_};
  auto I1172 = make_shared<Tensor>(I1172_index);
  auto tensor706 = vector<shared_ptr<Tensor>>{I1120, Gamma170_(), I1172};
  auto task706 = make_shared<Task706>(tensor706, cindex);
  task659->add_dep(task706);
  task706->add_dep(task658);
  deci2q->add_task(task706);

  auto tensor707 = vector<shared_ptr<Tensor>>{I1172, v2_, l2};
  auto task707 = make_shared<Task707>(tensor707, cindex);
  task706->add_dep(task707);
  task707->add_dep(task658);
  deci2q->add_task(task707);

  auto tensor708 = vector<shared_ptr<Tensor>>{I1172, v2_, l2};
  auto task708 = make_shared<Task708>(tensor708, cindex);
  task706->add_dep(task708);
  task708->add_dep(task658);
  deci2q->add_task(task708);

  auto tensor709 = vector<shared_ptr<Tensor>>{I1172, h1_, l2};
  auto task709 = make_shared<Task709>(tensor709, cindex);
  task706->add_dep(task709);
  task709->add_dep(task658);
  deci2q->add_task(task709);

  auto tensor710 = vector<shared_ptr<Tensor>>{I1172, h1_, l2};
  auto task710 = make_shared<Task710>(tensor710, cindex);
  task706->add_dep(task710);
  task710->add_dep(task658);
  deci2q->add_task(task710);

  vector<IndexRange> I1178_index = {active_, active_, active_, active_, active_, active_};
  auto I1178 = make_shared<Tensor>(I1178_index);
  auto tensor711 = vector<shared_ptr<Tensor>>{I1120, Gamma341_(), I1178};
  auto task711 = make_shared<Task711>(tensor711, cindex);
  task659->add_dep(task711);
  task711->add_dep(task658);
  deci2q->add_task(task711);

  auto tensor712 = vector<shared_ptr<Tensor>>{I1178, v2_, l2};
  auto task712 = make_shared<Task712>(tensor712, cindex);
  task711->add_dep(task712);
  task712->add_dep(task658);
  deci2q->add_task(task712);

  vector<IndexRange> I1193_index = {active_, active_, active_, active_};
  auto I1193 = make_shared<Tensor>(I1193_index);
  auto tensor713 = vector<shared_ptr<Tensor>>{I1120, Gamma132_(), I1193};
  auto task713 = make_shared<Task713>(tensor713, cindex);
  task659->add_dep(task713);
  task713->add_dep(task658);
  deci2q->add_task(task713);

  auto tensor714 = vector<shared_ptr<Tensor>>{I1193, v2_, l2};
  auto task714 = make_shared<Task714>(tensor714, cindex);
  task713->add_dep(task714);
  task714->add_dep(task658);
  deci2q->add_task(task714);

  vector<IndexRange> I1205_index = {active_, active_, active_, active_};
  auto I1205 = make_shared<Tensor>(I1205_index);
  auto tensor715 = vector<shared_ptr<Tensor>>{I1120, Gamma122_(), I1205};
  auto task715 = make_shared<Task715>(tensor715, cindex);
  task659->add_dep(task715);
  task715->add_dep(task658);
  deci2q->add_task(task715);

  auto tensor716 = vector<shared_ptr<Tensor>>{I1205, v2_, l2};
  auto task716 = make_shared<Task716>(tensor716, cindex);
  task715->add_dep(task716);
  task716->add_dep(task658);
  deci2q->add_task(task716);

  auto tensor717 = vector<shared_ptr<Tensor>>{I1205, h1_, l2};
  auto task717 = make_shared<Task717>(tensor717, cindex);
  task715->add_dep(task717);
  task717->add_dep(task658);
  deci2q->add_task(task717);

  vector<IndexRange> I1217_index = {active_, active_, active_, active_, active_, active_};
  auto I1217 = make_shared<Tensor>(I1217_index);
  auto tensor718 = vector<shared_ptr<Tensor>>{I1120, Gamma161_(), I1217};
  auto task718 = make_shared<Task718>(tensor718, cindex);
  task659->add_dep(task718);
  task718->add_dep(task658);
  deci2q->add_task(task718);

  auto tensor719 = vector<shared_ptr<Tensor>>{I1217, v2_, l2};
  auto task719 = make_shared<Task719>(tensor719, cindex);
  task718->add_dep(task719);
  task719->add_dep(task658);
  deci2q->add_task(task719);

  return deci2q;
}


#endif
