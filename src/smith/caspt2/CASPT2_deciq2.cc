//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_deciq2.cc
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
#include <src/smith/caspt2/CASPT2_tasks13.h>
#include <src/smith/caspt2/CASPT2_tasks14.h>
#include <src/smith/caspt2/CASPT2_tasks15.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void CASPT2::CASPT2::make_deciq2(shared_ptr<Queue> deciq, shared_ptr<Task> task538, shared_ptr<Task> task539, const bool diagonal, shared_ptr<Tensor> I698) {
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  vector<IndexRange> I759_index = {active_, active_};
  auto I759 = make_shared<Tensor>(I759_index);
  auto tensor601 = vector<shared_ptr<Tensor>>{I698, Gamma264_(), I759};
  auto task601 = make_shared<Task601>(tensor601, cindex);
  task539->add_dep(task601);
  task601->add_dep(task538);
  deciq->add_task(task601);

  vector<IndexRange> I760_index = {active_, closed_, virt_, closed_};
  auto I760 = make_shared<Tensor>(I760_index);
  auto tensor602 = vector<shared_ptr<Tensor>>{I759, t2, I760};
  auto task602 = make_shared<Task602>(tensor602, cindex);
  task601->add_dep(task602);
  task602->add_dep(task538);
  deciq->add_task(task602);

  auto tensor603 = vector<shared_ptr<Tensor>>{I760, f1_, t2};
  auto task603 = make_shared<Task603>(tensor603, cindex);
  task602->add_dep(task603);
  task603->add_dep(task538);
  deciq->add_task(task603);

  vector<IndexRange> I764_index = {active_, closed_, virt_, closed_};
  auto I764 = make_shared<Tensor>(I764_index);
  auto tensor604 = vector<shared_ptr<Tensor>>{I759, t2, I764};
  auto task604 = make_shared<Task604>(tensor604, cindex);
  task601->add_dep(task604);
  task604->add_dep(task538);
  deciq->add_task(task604);

  auto tensor605 = vector<shared_ptr<Tensor>>{I764, f1_, t2};
  auto task605 = make_shared<Task605>(tensor605, cindex);
  task604->add_dep(task605);
  task605->add_dep(task538);
  deciq->add_task(task605);

  vector<IndexRange> I768_index = {active_, virt_, closed_, closed_};
  auto I768 = make_shared<Tensor>(I768_index);
  auto tensor606 = vector<shared_ptr<Tensor>>{I759, t2, I768};
  auto task606 = make_shared<Task606>(tensor606, cindex);
  task601->add_dep(task606);
  task606->add_dep(task538);
  deciq->add_task(task606);

  auto tensor607 = vector<shared_ptr<Tensor>>{I768, f1_, t2};
  auto task607 = make_shared<Task607>(tensor607, cindex);
  task606->add_dep(task607);
  task607->add_dep(task538);
  deciq->add_task(task607);

  vector<IndexRange> I772_index = {active_, closed_, closed_, virt_};
  auto I772 = make_shared<Tensor>(I772_index);
  auto tensor608 = vector<shared_ptr<Tensor>>{I759, t2, I772};
  auto task608 = make_shared<Task608>(tensor608, cindex);
  task601->add_dep(task608);
  task608->add_dep(task538);
  deciq->add_task(task608);

  auto tensor609 = vector<shared_ptr<Tensor>>{I772, f1_, t2};
  auto task609 = make_shared<Task609>(tensor609, cindex);
  task608->add_dep(task609);
  task609->add_dep(task538);
  deciq->add_task(task609);

  vector<IndexRange> I776_index = {active_, virt_, closed_, closed_};
  auto I776 = make_shared<Tensor>(I776_index);
  auto tensor610 = vector<shared_ptr<Tensor>>{I759, t2, I776};
  auto task610 = make_shared<Task610>(tensor610, cindex);
  task601->add_dep(task610);
  task610->add_dep(task538);
  deciq->add_task(task610);

  auto tensor611 = vector<shared_ptr<Tensor>>{I776, f1_, t2};
  auto task611 = make_shared<Task611>(tensor611, cindex);
  task610->add_dep(task611);
  task611->add_dep(task538);
  deciq->add_task(task611);

  vector<IndexRange> I780_index = {active_, closed_, closed_, virt_};
  auto I780 = make_shared<Tensor>(I780_index);
  auto tensor612 = vector<shared_ptr<Tensor>>{I759, t2, I780};
  auto task612 = make_shared<Task612>(tensor612, cindex);
  task601->add_dep(task612);
  task612->add_dep(task538);
  deciq->add_task(task612);

  auto tensor613 = vector<shared_ptr<Tensor>>{I780, f1_, t2};
  auto task613 = make_shared<Task613>(tensor613, cindex);
  task612->add_dep(task613);
  task613->add_dep(task538);
  deciq->add_task(task613);

  vector<IndexRange> I800_index = {active_, virt_};
  auto I800 = make_shared<Tensor>(I800_index);
  auto tensor614 = vector<shared_ptr<Tensor>>{I759, f1_, I800};
  auto task614 = make_shared<Task614>(tensor614, cindex);
  task601->add_dep(task614);
  task614->add_dep(task538);
  deciq->add_task(task614);

  auto tensor615 = vector<shared_ptr<Tensor>>{I800, t2};
  auto task615 = make_shared<Task615>(tensor615, cindex);
  task614->add_dep(task615);
  task615->add_dep(task538);
  deciq->add_task(task615);

  auto tensor616 = vector<shared_ptr<Tensor>>{I800, t2};
  auto task616 = make_shared<Task616>(tensor616, cindex);
  task614->add_dep(task616);
  task616->add_dep(task538);
  deciq->add_task(task616);

  vector<IndexRange> I943_index = {virt_, active_};
  auto I943 = make_shared<Tensor>(I943_index);
  auto tensor617 = vector<shared_ptr<Tensor>>{I759, f1_, I943};
  auto task617 = make_shared<Task617>(tensor617, cindex);
  task601->add_dep(task617);
  task617->add_dep(task538);
  deciq->add_task(task617);

  auto tensor618 = vector<shared_ptr<Tensor>>{I943, t2};
  auto task618 = make_shared<Task618>(tensor618, cindex);
  task617->add_dep(task618);
  task618->add_dep(task538);
  deciq->add_task(task618);

  vector<IndexRange> I947_index = {virt_, active_};
  auto I947 = make_shared<Tensor>(I947_index);
  auto tensor619 = vector<shared_ptr<Tensor>>{I759, f1_, I947};
  auto task619 = make_shared<Task619>(tensor619, cindex);
  task601->add_dep(task619);
  task619->add_dep(task538);
  deciq->add_task(task619);

  auto tensor620 = vector<shared_ptr<Tensor>>{I947, t2};
  auto task620 = make_shared<Task620>(tensor620, cindex);
  task619->add_dep(task620);
  task620->add_dep(task538);
  deciq->add_task(task620);

  auto tensor621 = vector<shared_ptr<Tensor>>{I759, t2};
  auto task621 = make_shared<Task621>(tensor621, cindex, this->e0_);
  task601->add_dep(task621);
  task621->add_dep(task538);
  deciq->add_task(task621);

  auto tensor622 = vector<shared_ptr<Tensor>>{I759, t2};
  auto task622 = make_shared<Task622>(tensor622, cindex, this->e0_);
  task601->add_dep(task622);
  task622->add_dep(task538);
  deciq->add_task(task622);

  auto tensor623 = vector<shared_ptr<Tensor>>{I759, v2_, t2};
  auto task623 = make_shared<Task623>(tensor623, cindex);
  task601->add_dep(task623);
  task623->add_dep(task538);
  deciq->add_task(task623);

  auto tensor624 = vector<shared_ptr<Tensor>>{I759, v2_, t2};
  auto task624 = make_shared<Task624>(tensor624, cindex);
  task601->add_dep(task624);
  task624->add_dep(task538);
  deciq->add_task(task624);

  auto tensor625 = vector<shared_ptr<Tensor>>{I759, v2_, t2};
  auto task625 = make_shared<Task625>(tensor625, cindex);
  task601->add_dep(task625);
  task625->add_dep(task538);
  deciq->add_task(task625);

  auto tensor626 = vector<shared_ptr<Tensor>>{I759, v2_, t2};
  auto task626 = make_shared<Task626>(tensor626, cindex);
  task601->add_dep(task626);
  task626->add_dep(task538);
  deciq->add_task(task626);

  vector<IndexRange> I783_index = {active_, active_, active_, active_};
  auto I783 = make_shared<Tensor>(I783_index);
  auto tensor627 = vector<shared_ptr<Tensor>>{I698, Gamma270_(), I783};
  auto task627 = make_shared<Task627>(tensor627, cindex);
  task539->add_dep(task627);
  task627->add_dep(task538);
  deciq->add_task(task627);

  vector<IndexRange> I784_index = {active_, closed_, virt_, active_};
  auto I784 = make_shared<Tensor>(I784_index);
  auto tensor628 = vector<shared_ptr<Tensor>>{I783, t2, I784};
  auto task628 = make_shared<Task628>(tensor628, cindex);
  task627->add_dep(task628);
  task628->add_dep(task538);
  deciq->add_task(task628);

  auto tensor629 = vector<shared_ptr<Tensor>>{I784, f1_, t2};
  auto task629 = make_shared<Task629>(tensor629, cindex);
  task628->add_dep(task629);
  task629->add_dep(task538);
  deciq->add_task(task629);

  auto tensor630 = vector<shared_ptr<Tensor>>{I783, v2_, t2};
  auto task630 = make_shared<Task630>(tensor630, cindex);
  task627->add_dep(task630);
  task630->add_dep(task538);
  deciq->add_task(task630);

  vector<IndexRange> I807_index = {active_, active_, active_, active_, active_, active_};
  auto I807 = make_shared<Tensor>(I807_index);
  auto tensor631 = vector<shared_ptr<Tensor>>{I698, Gamma276_(), I807};
  auto task631 = make_shared<Task631>(tensor631, cindex);
  task539->add_dep(task631);
  task631->add_dep(task538);
  deciq->add_task(task631);

  vector<IndexRange> I808_index = {active_, closed_, active_, active_};
  auto I808 = make_shared<Tensor>(I808_index);
  auto tensor632 = vector<shared_ptr<Tensor>>{I807, t2, I808};
  auto task632 = make_shared<Task632>(tensor632, cindex);
  task631->add_dep(task632);
  task632->add_dep(task538);
  deciq->add_task(task632);

  auto tensor633 = vector<shared_ptr<Tensor>>{I808, f1_, t2};
  auto task633 = make_shared<Task633>(tensor633, cindex);
  task632->add_dep(task633);
  task633->add_dep(task538);
  deciq->add_task(task633);

  vector<IndexRange> I811_index = {active_, active_, active_, active_};
  auto I811 = make_shared<Tensor>(I811_index);
  auto tensor634 = vector<shared_ptr<Tensor>>{I698, Gamma277_(), I811};
  auto task634 = make_shared<Task634>(tensor634, cindex);
  task539->add_dep(task634);
  task634->add_dep(task538);
  deciq->add_task(task634);

  vector<IndexRange> I812_index = {active_, virt_, closed_, active_};
  auto I812 = make_shared<Tensor>(I812_index);
  auto tensor635 = vector<shared_ptr<Tensor>>{I811, t2, I812};
  auto task635 = make_shared<Task635>(tensor635, cindex);
  task634->add_dep(task635);
  task635->add_dep(task538);
  deciq->add_task(task635);

  auto tensor636 = vector<shared_ptr<Tensor>>{I812, t2, f1_};
  auto task636 = make_shared<Task636>(tensor636, cindex);
  task635->add_dep(task636);
  task636->add_dep(task538);
  deciq->add_task(task636);

  auto tensor637 = vector<shared_ptr<Tensor>>{I811, v2_, t2};
  auto task637 = make_shared<Task637>(tensor637, cindex);
  task634->add_dep(task637);
  task637->add_dep(task538);
  deciq->add_task(task637);

  vector<IndexRange> I819_index = {active_, active_, active_, active_};
  auto I819 = make_shared<Tensor>(I819_index);
  auto tensor638 = vector<shared_ptr<Tensor>>{I698, Gamma279_(), I819};
  auto task638 = make_shared<Task638>(tensor638, cindex);
  task539->add_dep(task638);
  task638->add_dep(task538);
  deciq->add_task(task638);

  auto tensor639 = vector<shared_ptr<Tensor>>{I819, t2};
  auto task639 = make_shared<Task639>(tensor639, cindex);
  task638->add_dep(task639);
  task639->add_dep(task538);
  deciq->add_task(task639);

  vector<IndexRange> I822_index = {active_, active_, active_, active_};
  auto I822 = make_shared<Tensor>(I822_index);
  auto tensor640 = vector<shared_ptr<Tensor>>{I698, Gamma280_(), I822};
  auto task640 = make_shared<Task640>(tensor640, cindex);
  task539->add_dep(task640);
  task640->add_dep(task538);
  deciq->add_task(task640);

  vector<IndexRange> I823_index = {active_, virt_, active_, closed_};
  auto I823 = make_shared<Tensor>(I823_index);
  auto tensor641 = vector<shared_ptr<Tensor>>{I822, t2, I823};
  auto task641 = make_shared<Task641>(tensor641, cindex);
  task640->add_dep(task641);
  task641->add_dep(task538);
  deciq->add_task(task641);

  auto tensor642 = vector<shared_ptr<Tensor>>{I823, f1_, t2};
  auto task642 = make_shared<Task642>(tensor642, cindex);
  task641->add_dep(task642);
  task642->add_dep(task538);
  deciq->add_task(task642);

  vector<IndexRange> I827_index = {active_, closed_, active_, virt_};
  auto I827 = make_shared<Tensor>(I827_index);
  auto tensor643 = vector<shared_ptr<Tensor>>{I822, t2, I827};
  auto task643 = make_shared<Task643>(tensor643, cindex);
  task640->add_dep(task643);
  task643->add_dep(task538);
  deciq->add_task(task643);

  auto tensor644 = vector<shared_ptr<Tensor>>{I827, f1_, t2};
  auto task644 = make_shared<Task644>(tensor644, cindex);
  task643->add_dep(task644);
  task644->add_dep(task538);
  deciq->add_task(task644);

  vector<IndexRange> I858_index = {active_, active_, virt_, closed_};
  auto I858 = make_shared<Tensor>(I858_index);
  auto tensor645 = vector<shared_ptr<Tensor>>{I822, t2, I858};
  auto task645 = make_shared<Task645>(tensor645, cindex);
  task640->add_dep(task645);
  task645->add_dep(task538);
  deciq->add_task(task645);

  auto tensor646 = vector<shared_ptr<Tensor>>{I858, t2, f1_};
  auto task646 = make_shared<Task646>(tensor646, cindex);
  task645->add_dep(task646);
  task646->add_dep(task538);
  deciq->add_task(task646);

  vector<IndexRange> I985_index = {closed_, virt_, active_, active_};
  auto I985 = make_shared<Tensor>(I985_index);
  auto tensor647 = vector<shared_ptr<Tensor>>{I822, t2, I985};
  auto task647 = make_shared<Task647>(tensor647, cindex);
  task640->add_dep(task647);
  task647->add_dep(task538);
  deciq->add_task(task647);

  auto tensor648 = vector<shared_ptr<Tensor>>{I985, t2};
  auto task648 = make_shared<Task648>(tensor648, cindex, this->e0_);
  task647->add_dep(task648);
  task648->add_dep(task538);
  deciq->add_task(task648);

  auto tensor649 = vector<shared_ptr<Tensor>>{I985, f1_, t2};
  auto task649 = make_shared<Task649>(tensor649, cindex);
  task647->add_dep(task649);
  task649->add_dep(task538);
  deciq->add_task(task649);

  auto tensor650 = vector<shared_ptr<Tensor>>{I822, v2_, t2};
  auto task650 = make_shared<Task650>(tensor650, cindex);
  task640->add_dep(task650);
  task650->add_dep(task538);
  deciq->add_task(task650);

  auto tensor651 = vector<shared_ptr<Tensor>>{I822, v2_, t2};
  auto task651 = make_shared<Task651>(tensor651, cindex);
  task640->add_dep(task651);
  task651->add_dep(task538);
  deciq->add_task(task651);

  vector<IndexRange> I830_index = {active_, active_, active_, active_};
  auto I830 = make_shared<Tensor>(I830_index);
  auto tensor652 = vector<shared_ptr<Tensor>>{I698, Gamma282_(), I830};
  auto task652 = make_shared<Task652>(tensor652, cindex);
  task539->add_dep(task652);
  task652->add_dep(task538);
  deciq->add_task(task652);

  auto tensor653 = vector<shared_ptr<Tensor>>{I830, t2};
  auto task653 = make_shared<Task653>(tensor653, cindex);
  task652->add_dep(task653);
  task653->add_dep(task538);
  deciq->add_task(task653);

  auto tensor654 = vector<shared_ptr<Tensor>>{I830, t2};
  auto task654 = make_shared<Task654>(tensor654, cindex);
  task652->add_dep(task654);
  task654->add_dep(task538);
  deciq->add_task(task654);

  auto tensor655 = vector<shared_ptr<Tensor>>{I830, t2};
  auto task655 = make_shared<Task655>(tensor655, cindex);
  task652->add_dep(task655);
  task655->add_dep(task538);
  deciq->add_task(task655);

  vector<IndexRange> I833_index = {active_, active_, active_, active_};
  auto I833 = make_shared<Tensor>(I833_index);
  auto tensor656 = vector<shared_ptr<Tensor>>{I698, Gamma283_(), I833};
  auto task656 = make_shared<Task656>(tensor656, cindex);
  task539->add_dep(task656);
  task656->add_dep(task538);
  deciq->add_task(task656);

  vector<IndexRange> I834_index = {active_, virt_, active_, closed_};
  auto I834 = make_shared<Tensor>(I834_index);
  auto tensor657 = vector<shared_ptr<Tensor>>{I833, t2, I834};
  auto task657 = make_shared<Task657>(tensor657, cindex);
  task656->add_dep(task657);
  task657->add_dep(task538);
  deciq->add_task(task657);

  auto tensor658 = vector<shared_ptr<Tensor>>{I834, f1_, t2};
  auto task658 = make_shared<Task658>(tensor658, cindex);
  task657->add_dep(task658);
  task658->add_dep(task538);
  deciq->add_task(task658);

  vector<IndexRange> I838_index = {active_, closed_, active_, virt_};
  auto I838 = make_shared<Tensor>(I838_index);
  auto tensor659 = vector<shared_ptr<Tensor>>{I833, t2, I838};
  auto task659 = make_shared<Task659>(tensor659, cindex);
  task656->add_dep(task659);
  task659->add_dep(task538);
  deciq->add_task(task659);

  auto tensor660 = vector<shared_ptr<Tensor>>{I838, f1_, t2};
  auto task660 = make_shared<Task660>(tensor660, cindex);
  task659->add_dep(task660);
  task660->add_dep(task538);
  deciq->add_task(task660);

  auto tensor661 = vector<shared_ptr<Tensor>>{I838, f1_, t2};
  auto task661 = make_shared<Task661>(tensor661, cindex);
  task659->add_dep(task661);
  task661->add_dep(task538);
  deciq->add_task(task661);

  vector<IndexRange> I854_index = {active_, active_, closed_, virt_};
  auto I854 = make_shared<Tensor>(I854_index);
  auto tensor662 = vector<shared_ptr<Tensor>>{I833, t2, I854};
  auto task662 = make_shared<Task662>(tensor662, cindex);
  task656->add_dep(task662);
  task662->add_dep(task538);
  deciq->add_task(task662);

  auto tensor663 = vector<shared_ptr<Tensor>>{I854, t2, f1_};
  auto task663 = make_shared<Task663>(tensor663, cindex);
  task662->add_dep(task663);
  task663->add_dep(task538);
  deciq->add_task(task663);

  vector<IndexRange> I877_index = {active_, active_, virt_, closed_};
  auto I877 = make_shared<Tensor>(I877_index);
  auto tensor664 = vector<shared_ptr<Tensor>>{I833, t2, I877};
  auto task664 = make_shared<Task664>(tensor664, cindex);
  task656->add_dep(task664);
  task664->add_dep(task538);
  deciq->add_task(task664);

  auto tensor665 = vector<shared_ptr<Tensor>>{I877, f1_, t2};
  auto task665 = make_shared<Task665>(tensor665, cindex);
  task664->add_dep(task665);
  task665->add_dep(task538);
  deciq->add_task(task665);

  vector<IndexRange> I881_index = {active_, active_, closed_, virt_};
  auto I881 = make_shared<Tensor>(I881_index);
  auto tensor666 = vector<shared_ptr<Tensor>>{I833, t2, I881};
  auto task666 = make_shared<Task666>(tensor666, cindex);
  task656->add_dep(task666);
  task666->add_dep(task538);
  deciq->add_task(task666);

  auto tensor667 = vector<shared_ptr<Tensor>>{I881, f1_, t2};
  auto task667 = make_shared<Task667>(tensor667, cindex);
  task666->add_dep(task667);
  task667->add_dep(task538);
  deciq->add_task(task667);

  vector<IndexRange> I888_index = {active_, active_, virt_, closed_};
  auto I888 = make_shared<Tensor>(I888_index);
  auto tensor668 = vector<shared_ptr<Tensor>>{I833, t2, I888};
  auto task668 = make_shared<Task668>(tensor668, cindex);
  task656->add_dep(task668);
  task668->add_dep(task538);
  deciq->add_task(task668);

  auto tensor669 = vector<shared_ptr<Tensor>>{I888, f1_, t2};
  auto task669 = make_shared<Task669>(tensor669, cindex);
  task668->add_dep(task669);
  task669->add_dep(task538);
  deciq->add_task(task669);

  vector<IndexRange> I892_index = {active_, active_, closed_, virt_};
  auto I892 = make_shared<Tensor>(I892_index);
  auto tensor670 = vector<shared_ptr<Tensor>>{I833, t2, I892};
  auto task670 = make_shared<Task670>(tensor670, cindex);
  task656->add_dep(task670);
  task670->add_dep(task538);
  deciq->add_task(task670);

  auto tensor671 = vector<shared_ptr<Tensor>>{I892, f1_, t2};
  auto task671 = make_shared<Task671>(tensor671, cindex);
  task670->add_dep(task671);
  task671->add_dep(task538);
  deciq->add_task(task671);

  vector<IndexRange> I908_index = {active_, active_, closed_, virt_};
  auto I908 = make_shared<Tensor>(I908_index);
  auto tensor672 = vector<shared_ptr<Tensor>>{I833, t2, I908};
  auto task672 = make_shared<Task672>(tensor672, cindex);
  task656->add_dep(task672);
  task672->add_dep(task538);
  deciq->add_task(task672);

  auto tensor673 = vector<shared_ptr<Tensor>>{I908, t2, f1_};
  auto task673 = make_shared<Task673>(tensor673, cindex);
  task672->add_dep(task673);
  task673->add_dep(task538);
  deciq->add_task(task673);

  auto tensor674 = vector<shared_ptr<Tensor>>{I908, t2, f1_};
  auto task674 = make_shared<Task674>(tensor674, cindex);
  task672->add_dep(task674);
  task674->add_dep(task538);
  deciq->add_task(task674);

  vector<IndexRange> I981_index = {virt_, closed_, active_, active_};
  auto I981 = make_shared<Tensor>(I981_index);
  auto tensor675 = vector<shared_ptr<Tensor>>{I833, t2, I981};
  auto task675 = make_shared<Task675>(tensor675, cindex);
  task656->add_dep(task675);
  task675->add_dep(task538);
  deciq->add_task(task675);

  auto tensor676 = vector<shared_ptr<Tensor>>{I981, f1_, t2};
  auto task676 = make_shared<Task676>(tensor676, cindex);
  task675->add_dep(task676);
  task676->add_dep(task538);
  deciq->add_task(task676);

  vector<IndexRange> I993_index = {closed_, virt_, active_, active_};
  auto I993 = make_shared<Tensor>(I993_index);
  auto tensor677 = vector<shared_ptr<Tensor>>{I833, t2, I993};
  auto task677 = make_shared<Task677>(tensor677, cindex);
  task656->add_dep(task677);
  task677->add_dep(task538);
  deciq->add_task(task677);

  auto tensor678 = vector<shared_ptr<Tensor>>{I993, t2};
  auto task678 = make_shared<Task678>(tensor678, cindex, this->e0_);
  task677->add_dep(task678);
  task678->add_dep(task538);
  deciq->add_task(task678);

  auto tensor679 = vector<shared_ptr<Tensor>>{I993, f1_, t2};
  auto task679 = make_shared<Task679>(tensor679, cindex);
  task677->add_dep(task679);
  task679->add_dep(task538);
  deciq->add_task(task679);

  auto tensor680 = vector<shared_ptr<Tensor>>{I833, t2};
  auto task680 = make_shared<Task680>(tensor680, cindex, this->e0_);
  task656->add_dep(task680);
  task680->add_dep(task538);
  deciq->add_task(task680);

  auto tensor681 = vector<shared_ptr<Tensor>>{I833, t2};
  auto task681 = make_shared<Task681>(tensor681, cindex, this->e0_);
  task656->add_dep(task681);
  task681->add_dep(task538);
  deciq->add_task(task681);

  auto tensor682 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task682 = make_shared<Task682>(tensor682, cindex);
  task656->add_dep(task682);
  task682->add_dep(task538);
  deciq->add_task(task682);

  auto tensor683 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task683 = make_shared<Task683>(tensor683, cindex);
  task656->add_dep(task683);
  task683->add_dep(task538);
  deciq->add_task(task683);

  auto tensor684 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task684 = make_shared<Task684>(tensor684, cindex);
  task656->add_dep(task684);
  task684->add_dep(task538);
  deciq->add_task(task684);

  auto tensor685 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task685 = make_shared<Task685>(tensor685, cindex);
  task656->add_dep(task685);
  task685->add_dep(task538);
  deciq->add_task(task685);

  auto tensor686 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task686 = make_shared<Task686>(tensor686, cindex);
  task656->add_dep(task686);
  task686->add_dep(task538);
  deciq->add_task(task686);

  auto tensor687 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task687 = make_shared<Task687>(tensor687, cindex);
  task656->add_dep(task687);
  task687->add_dep(task538);
  deciq->add_task(task687);

  auto tensor688 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task688 = make_shared<Task688>(tensor688, cindex);
  task656->add_dep(task688);
  task688->add_dep(task538);
  deciq->add_task(task688);

  auto tensor689 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task689 = make_shared<Task689>(tensor689, cindex);
  task656->add_dep(task689);
  task689->add_dep(task538);
  deciq->add_task(task689);

  auto tensor690 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task690 = make_shared<Task690>(tensor690, cindex);
  task656->add_dep(task690);
  task690->add_dep(task538);
  deciq->add_task(task690);

  auto tensor691 = vector<shared_ptr<Tensor>>{I833, v2_, t2};
  auto task691 = make_shared<Task691>(tensor691, cindex);
  task656->add_dep(task691);
  task691->add_dep(task538);
  deciq->add_task(task691);

  vector<IndexRange> I841_index = {active_, active_, active_, active_, active_, active_};
  auto I841 = make_shared<Tensor>(I841_index);
  auto tensor692 = vector<shared_ptr<Tensor>>{I698, Gamma285_(), I841};
  auto task692 = make_shared<Task692>(tensor692, cindex);
  task539->add_dep(task692);
  task692->add_dep(task538);
  deciq->add_task(task692);

  vector<IndexRange> I842_index = {active_, virt_, active_, active_};
  auto I842 = make_shared<Tensor>(I842_index);
  auto tensor693 = vector<shared_ptr<Tensor>>{I841, t2, I842};
  auto task693 = make_shared<Task693>(tensor693, cindex);
  task692->add_dep(task693);
  task693->add_dep(task538);
  deciq->add_task(task693);

  auto tensor694 = vector<shared_ptr<Tensor>>{I842, f1_, t2};
  auto task694 = make_shared<Task694>(tensor694, cindex);
  task693->add_dep(task694);
  task694->add_dep(task538);
  deciq->add_task(task694);

  vector<IndexRange> I845_index = {active_, active_};
  auto I845 = make_shared<Tensor>(I845_index);
  auto tensor695 = vector<shared_ptr<Tensor>>{I698, Gamma286_(), I845};
  auto task695 = make_shared<Task695>(tensor695, cindex);
  task539->add_dep(task695);
  task695->add_dep(task538);
  deciq->add_task(task695);

  vector<IndexRange> I846_index = {closed_, virt_};
  auto I846 = make_shared<Tensor>(I846_index);
  auto tensor696 = vector<shared_ptr<Tensor>>{I845, t2, I846};
  auto task696 = make_shared<Task696>(tensor696, cindex);
  task695->add_dep(task696);
  task696->add_dep(task538);
  deciq->add_task(task696);

  auto tensor697 = vector<shared_ptr<Tensor>>{I846, t2, f1_};
  auto task697 = make_shared<Task697>(tensor697, cindex);
  task696->add_dep(task697);
  task697->add_dep(task538);
  deciq->add_task(task697);

  auto tensor698 = vector<shared_ptr<Tensor>>{I846, t2, f1_};
  auto task698 = make_shared<Task698>(tensor698, cindex);
  task696->add_dep(task698);
  task698->add_dep(task538);
  deciq->add_task(task698);

  vector<IndexRange> I900_index = {closed_, virt_};
  auto I900 = make_shared<Tensor>(I900_index);
  auto tensor699 = vector<shared_ptr<Tensor>>{I845, t2, I900};
  auto task699 = make_shared<Task699>(tensor699, cindex);
  task695->add_dep(task699);
  task699->add_dep(task538);
  deciq->add_task(task699);

  auto tensor700 = vector<shared_ptr<Tensor>>{I900, t2, f1_};
  auto task700 = make_shared<Task700>(tensor700, cindex);
  task699->add_dep(task700);
  task700->add_dep(task538);
  deciq->add_task(task700);

  auto tensor701 = vector<shared_ptr<Tensor>>{I900, t2, f1_};
  auto task701 = make_shared<Task701>(tensor701, cindex);
  task699->add_dep(task701);
  task701->add_dep(task538);
  deciq->add_task(task701);

  vector<IndexRange> I951_index = {virt_, closed_};
  auto I951 = make_shared<Tensor>(I951_index);
  auto tensor702 = vector<shared_ptr<Tensor>>{I845, t2, I951};
  auto task702 = make_shared<Task702>(tensor702, cindex);
  task695->add_dep(task702);
  task702->add_dep(task538);
  deciq->add_task(task702);

  auto tensor703 = vector<shared_ptr<Tensor>>{I951, f1_, t2};
  auto task703 = make_shared<Task703>(tensor703, cindex);
  task702->add_dep(task703);
  task703->add_dep(task538);
  deciq->add_task(task703);

  vector<IndexRange> I955_index = {virt_, closed_};
  auto I955 = make_shared<Tensor>(I955_index);
  auto tensor704 = vector<shared_ptr<Tensor>>{I845, t2, I955};
  auto task704 = make_shared<Task704>(tensor704, cindex);
  task695->add_dep(task704);
  task704->add_dep(task538);
  deciq->add_task(task704);

  auto tensor705 = vector<shared_ptr<Tensor>>{I955, f1_, t2};
  auto task705 = make_shared<Task705>(tensor705, cindex);
  task704->add_dep(task705);
  task705->add_dep(task538);
  deciq->add_task(task705);

  vector<IndexRange> I959_index = {virt_, closed_};
  auto I959 = make_shared<Tensor>(I959_index);
  auto tensor706 = vector<shared_ptr<Tensor>>{I845, t2, I959};
  auto task706 = make_shared<Task706>(tensor706, cindex);
  task695->add_dep(task706);
  task706->add_dep(task538);
  deciq->add_task(task706);

  auto tensor707 = vector<shared_ptr<Tensor>>{I959, f1_, t2};
  auto task707 = make_shared<Task707>(tensor707, cindex);
  task706->add_dep(task707);
  task707->add_dep(task538);
  deciq->add_task(task707);

  vector<IndexRange> I963_index = {virt_, closed_};
  auto I963 = make_shared<Tensor>(I963_index);
  auto tensor708 = vector<shared_ptr<Tensor>>{I845, t2, I963};
  auto task708 = make_shared<Task708>(tensor708, cindex);
  task695->add_dep(task708);
  task708->add_dep(task538);
  deciq->add_task(task708);

  auto tensor709 = vector<shared_ptr<Tensor>>{I963, f1_, t2};
  auto task709 = make_shared<Task709>(tensor709, cindex);
  task708->add_dep(task709);
  task709->add_dep(task538);
  deciq->add_task(task709);

  vector<IndexRange> I973_index = {closed_, active_};
  auto I973 = make_shared<Tensor>(I973_index);
  auto tensor710 = vector<shared_ptr<Tensor>>{I845, f1_, I973};
  auto task710 = make_shared<Task710>(tensor710, cindex);
  task695->add_dep(task710);
  task710->add_dep(task538);
  deciq->add_task(task710);

  auto tensor711 = vector<shared_ptr<Tensor>>{I973, t2};
  auto task711 = make_shared<Task711>(tensor711, cindex);
  task710->add_dep(task711);
  task711->add_dep(task538);
  deciq->add_task(task711);

  auto tensor712 = vector<shared_ptr<Tensor>>{I973, t2};
  auto task712 = make_shared<Task712>(tensor712, cindex);
  task710->add_dep(task712);
  task712->add_dep(task538);
  deciq->add_task(task712);

  vector<IndexRange> I1005_index = {active_, closed_};
  auto I1005 = make_shared<Tensor>(I1005_index);
  auto tensor713 = vector<shared_ptr<Tensor>>{I845, f1_, I1005};
  auto task713 = make_shared<Task713>(tensor713, cindex);
  task695->add_dep(task713);
  task713->add_dep(task538);
  deciq->add_task(task713);

  auto tensor714 = vector<shared_ptr<Tensor>>{I1005, t2};
  auto task714 = make_shared<Task714>(tensor714, cindex);
  task713->add_dep(task714);
  task714->add_dep(task538);
  deciq->add_task(task714);

  auto tensor715 = vector<shared_ptr<Tensor>>{I1005, t2};
  auto task715 = make_shared<Task715>(tensor715, cindex);
  task713->add_dep(task715);
  task715->add_dep(task538);
  deciq->add_task(task715);

  vector<IndexRange> I1019_index = {virt_, virt_, active_, closed_};
  auto I1019 = make_shared<Tensor>(I1019_index);
  auto tensor716 = vector<shared_ptr<Tensor>>{I845, t2, I1019};
  auto task716 = make_shared<Task716>(tensor716, cindex);
  task695->add_dep(task716);
  task716->add_dep(task538);
  deciq->add_task(task716);

  auto tensor717 = vector<shared_ptr<Tensor>>{I1019, f1_, t2};
  auto task717 = make_shared<Task717>(tensor717, cindex);
  task716->add_dep(task717);
  task717->add_dep(task538);
  deciq->add_task(task717);

  vector<IndexRange> I1023_index = {virt_, virt_, active_, closed_};
  auto I1023 = make_shared<Tensor>(I1023_index);
  auto tensor718 = vector<shared_ptr<Tensor>>{I845, t2, I1023};
  auto task718 = make_shared<Task718>(tensor718, cindex);
  task695->add_dep(task718);
  task718->add_dep(task538);
  deciq->add_task(task718);

  auto tensor719 = vector<shared_ptr<Tensor>>{I1023, f1_, t2};
  auto task719 = make_shared<Task719>(tensor719, cindex);
  task718->add_dep(task719);
  task719->add_dep(task538);
  deciq->add_task(task719);

  vector<IndexRange> I1027_index = {virt_, closed_, active_, virt_};
  auto I1027 = make_shared<Tensor>(I1027_index);
  auto tensor720 = vector<shared_ptr<Tensor>>{I845, t2, I1027};
  auto task720 = make_shared<Task720>(tensor720, cindex);
  task695->add_dep(task720);
  task720->add_dep(task538);
  deciq->add_task(task720);

  auto tensor721 = vector<shared_ptr<Tensor>>{I1027, f1_, t2};
  auto task721 = make_shared<Task721>(tensor721, cindex);
  task720->add_dep(task721);
  task721->add_dep(task538);
  deciq->add_task(task721);

  vector<IndexRange> I1031_index = {virt_, closed_, active_, virt_};
  auto I1031 = make_shared<Tensor>(I1031_index);
  auto tensor722 = vector<shared_ptr<Tensor>>{I845, t2, I1031};
  auto task722 = make_shared<Task722>(tensor722, cindex);
  task695->add_dep(task722);
  task722->add_dep(task538);
  deciq->add_task(task722);

  auto tensor723 = vector<shared_ptr<Tensor>>{I1031, f1_, t2};
  auto task723 = make_shared<Task723>(tensor723, cindex);
  task722->add_dep(task723);
  task723->add_dep(task538);
  deciq->add_task(task723);

  vector<IndexRange> I1035_index = {closed_, virt_, active_, virt_};
  auto I1035 = make_shared<Tensor>(I1035_index);
  auto tensor724 = vector<shared_ptr<Tensor>>{I845, t2, I1035};
  auto task724 = make_shared<Task724>(tensor724, cindex);
  task695->add_dep(task724);
  task724->add_dep(task538);
  deciq->add_task(task724);

  auto tensor725 = vector<shared_ptr<Tensor>>{I1035, f1_, t2};
  auto task725 = make_shared<Task725>(tensor725, cindex);
  task724->add_dep(task725);
  task725->add_dep(task538);
  deciq->add_task(task725);

  vector<IndexRange> I1039_index = {closed_, virt_, active_, virt_};
  auto I1039 = make_shared<Tensor>(I1039_index);
  auto tensor726 = vector<shared_ptr<Tensor>>{I845, t2, I1039};
  auto task726 = make_shared<Task726>(tensor726, cindex);
  task695->add_dep(task726);
  task726->add_dep(task538);
  deciq->add_task(task726);

  auto tensor727 = vector<shared_ptr<Tensor>>{I1039, f1_, t2};
  auto task727 = make_shared<Task727>(tensor727, cindex);
  task726->add_dep(task727);
  task727->add_dep(task538);
  deciq->add_task(task727);

  auto tensor728 = vector<shared_ptr<Tensor>>{I845, t2};
  auto task728 = make_shared<Task728>(tensor728, cindex, this->e0_);
  task695->add_dep(task728);
  task728->add_dep(task538);
  deciq->add_task(task728);

  auto tensor729 = vector<shared_ptr<Tensor>>{I845, t2};
  auto task729 = make_shared<Task729>(tensor729, cindex, this->e0_);
  task695->add_dep(task729);
  task729->add_dep(task538);
  deciq->add_task(task729);

  auto tensor730 = vector<shared_ptr<Tensor>>{I845, v2_, t2};
  auto task730 = make_shared<Task730>(tensor730, cindex);
  task695->add_dep(task730);
  task730->add_dep(task538);
  deciq->add_task(task730);

  auto tensor731 = vector<shared_ptr<Tensor>>{I845, v2_, t2};
  auto task731 = make_shared<Task731>(tensor731, cindex);
  task695->add_dep(task731);
  task731->add_dep(task538);
  deciq->add_task(task731);

  auto tensor732 = vector<shared_ptr<Tensor>>{I845, v2_, t2};
  auto task732 = make_shared<Task732>(tensor732, cindex);
  task695->add_dep(task732);
  task732->add_dep(task538);
  deciq->add_task(task732);

  auto tensor733 = vector<shared_ptr<Tensor>>{I845, v2_, t2};
  auto task733 = make_shared<Task733>(tensor733, cindex);
  task695->add_dep(task733);
  task733->add_dep(task538);
  deciq->add_task(task733);

  vector<IndexRange> I1209_index = {active_, closed_, virt_, active_};
  auto I1209 = make_shared<Tensor>(I1209_index);
  auto tensor734 = vector<shared_ptr<Tensor>>{I845, h1_, I1209};
  auto task734 = make_shared<Task734>(tensor734, cindex);
  task695->add_dep(task734);
  task734->add_dep(task538);
  deciq->add_task(task734);

  auto tensor735 = vector<shared_ptr<Tensor>>{I1209, t2};
  auto task735 = make_shared<Task735>(tensor735, cindex);
  task734->add_dep(task735);
  task735->add_dep(task538);
  deciq->add_task(task735);

  vector<IndexRange> I1212_index = {active_, active_, virt_, closed_};
  auto I1212 = make_shared<Tensor>(I1212_index);
  auto tensor736 = vector<shared_ptr<Tensor>>{I845, h1_, I1212};
  auto task736 = make_shared<Task736>(tensor736, cindex);
  task695->add_dep(task736);
  task736->add_dep(task538);
  deciq->add_task(task736);

  auto tensor737 = vector<shared_ptr<Tensor>>{I1212, t2};
  auto task737 = make_shared<Task737>(tensor737, cindex);
  task736->add_dep(task737);
  task737->add_dep(task538);
  deciq->add_task(task737);
}

#endif
