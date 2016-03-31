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

void CASPT2::CASPT2::make_deciq2(shared_ptr<Queue> deciq, shared_ptr<Task> task539, shared_ptr<Task> task540, const bool diagonal, shared_ptr<Tensor> I684) {
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  vector<IndexRange> I745_index = {active_, active_};
  auto I745 = make_shared<Tensor>(I745_index);
  auto tensor602 = vector<shared_ptr<Tensor>>{I684, Gamma264_(), I745};
  auto task602 = make_shared<Task602>(tensor602, cindex);
  task540->add_dep(task602);
  task602->add_dep(task539);
  deciq->add_task(task602);

  vector<IndexRange> I746_index = {active_, closed_, virt_, closed_};
  auto I746 = make_shared<Tensor>(I746_index);
  auto tensor603 = vector<shared_ptr<Tensor>>{I745, t2, I746};
  auto task603 = make_shared<Task603>(tensor603, cindex);
  task602->add_dep(task603);
  task603->add_dep(task539);
  deciq->add_task(task603);

  auto tensor604 = vector<shared_ptr<Tensor>>{I746, f1_, t2};
  auto task604 = make_shared<Task604>(tensor604, cindex);
  task603->add_dep(task604);
  task604->add_dep(task539);
  deciq->add_task(task604);

  vector<IndexRange> I750_index = {active_, closed_, virt_, closed_};
  auto I750 = make_shared<Tensor>(I750_index);
  auto tensor605 = vector<shared_ptr<Tensor>>{I745, t2, I750};
  auto task605 = make_shared<Task605>(tensor605, cindex);
  task602->add_dep(task605);
  task605->add_dep(task539);
  deciq->add_task(task605);

  auto tensor606 = vector<shared_ptr<Tensor>>{I750, f1_, t2};
  auto task606 = make_shared<Task606>(tensor606, cindex);
  task605->add_dep(task606);
  task606->add_dep(task539);
  deciq->add_task(task606);

  vector<IndexRange> I754_index = {active_, virt_, closed_, closed_};
  auto I754 = make_shared<Tensor>(I754_index);
  auto tensor607 = vector<shared_ptr<Tensor>>{I745, t2, I754};
  auto task607 = make_shared<Task607>(tensor607, cindex);
  task602->add_dep(task607);
  task607->add_dep(task539);
  deciq->add_task(task607);

  auto tensor608 = vector<shared_ptr<Tensor>>{I754, f1_, t2};
  auto task608 = make_shared<Task608>(tensor608, cindex);
  task607->add_dep(task608);
  task608->add_dep(task539);
  deciq->add_task(task608);

  vector<IndexRange> I758_index = {active_, closed_, closed_, virt_};
  auto I758 = make_shared<Tensor>(I758_index);
  auto tensor609 = vector<shared_ptr<Tensor>>{I745, t2, I758};
  auto task609 = make_shared<Task609>(tensor609, cindex);
  task602->add_dep(task609);
  task609->add_dep(task539);
  deciq->add_task(task609);

  auto tensor610 = vector<shared_ptr<Tensor>>{I758, f1_, t2};
  auto task610 = make_shared<Task610>(tensor610, cindex);
  task609->add_dep(task610);
  task610->add_dep(task539);
  deciq->add_task(task610);

  vector<IndexRange> I762_index = {active_, virt_, closed_, closed_};
  auto I762 = make_shared<Tensor>(I762_index);
  auto tensor611 = vector<shared_ptr<Tensor>>{I745, t2, I762};
  auto task611 = make_shared<Task611>(tensor611, cindex);
  task602->add_dep(task611);
  task611->add_dep(task539);
  deciq->add_task(task611);

  auto tensor612 = vector<shared_ptr<Tensor>>{I762, f1_, t2};
  auto task612 = make_shared<Task612>(tensor612, cindex);
  task611->add_dep(task612);
  task612->add_dep(task539);
  deciq->add_task(task612);

  vector<IndexRange> I766_index = {active_, closed_, closed_, virt_};
  auto I766 = make_shared<Tensor>(I766_index);
  auto tensor613 = vector<shared_ptr<Tensor>>{I745, t2, I766};
  auto task613 = make_shared<Task613>(tensor613, cindex);
  task602->add_dep(task613);
  task613->add_dep(task539);
  deciq->add_task(task613);

  auto tensor614 = vector<shared_ptr<Tensor>>{I766, f1_, t2};
  auto task614 = make_shared<Task614>(tensor614, cindex);
  task613->add_dep(task614);
  task614->add_dep(task539);
  deciq->add_task(task614);

  vector<IndexRange> I786_index = {active_, virt_};
  auto I786 = make_shared<Tensor>(I786_index);
  auto tensor615 = vector<shared_ptr<Tensor>>{I745, f1_, I786};
  auto task615 = make_shared<Task615>(tensor615, cindex);
  task602->add_dep(task615);
  task615->add_dep(task539);
  deciq->add_task(task615);

  auto tensor616 = vector<shared_ptr<Tensor>>{I786, t2};
  auto task616 = make_shared<Task616>(tensor616, cindex);
  task615->add_dep(task616);
  task616->add_dep(task539);
  deciq->add_task(task616);

  auto tensor617 = vector<shared_ptr<Tensor>>{I786, t2};
  auto task617 = make_shared<Task617>(tensor617, cindex);
  task615->add_dep(task617);
  task617->add_dep(task539);
  deciq->add_task(task617);

  vector<IndexRange> I929_index = {virt_, active_};
  auto I929 = make_shared<Tensor>(I929_index);
  auto tensor618 = vector<shared_ptr<Tensor>>{I745, f1_, I929};
  auto task618 = make_shared<Task618>(tensor618, cindex);
  task602->add_dep(task618);
  task618->add_dep(task539);
  deciq->add_task(task618);

  auto tensor619 = vector<shared_ptr<Tensor>>{I929, t2};
  auto task619 = make_shared<Task619>(tensor619, cindex);
  task618->add_dep(task619);
  task619->add_dep(task539);
  deciq->add_task(task619);

  vector<IndexRange> I933_index = {virt_, active_};
  auto I933 = make_shared<Tensor>(I933_index);
  auto tensor620 = vector<shared_ptr<Tensor>>{I745, f1_, I933};
  auto task620 = make_shared<Task620>(tensor620, cindex);
  task602->add_dep(task620);
  task620->add_dep(task539);
  deciq->add_task(task620);

  auto tensor621 = vector<shared_ptr<Tensor>>{I933, t2};
  auto task621 = make_shared<Task621>(tensor621, cindex);
  task620->add_dep(task621);
  task621->add_dep(task539);
  deciq->add_task(task621);

  auto tensor622 = vector<shared_ptr<Tensor>>{I745, t2};
  auto task622 = make_shared<Task622>(tensor622, cindex, this->e0_);
  task602->add_dep(task622);
  task622->add_dep(task539);
  deciq->add_task(task622);

  auto tensor623 = vector<shared_ptr<Tensor>>{I745, t2};
  auto task623 = make_shared<Task623>(tensor623, cindex, this->e0_);
  task602->add_dep(task623);
  task623->add_dep(task539);
  deciq->add_task(task623);

  auto tensor624 = vector<shared_ptr<Tensor>>{I745, v2_, t2};
  auto task624 = make_shared<Task624>(tensor624, cindex);
  task602->add_dep(task624);
  task624->add_dep(task539);
  deciq->add_task(task624);

  auto tensor625 = vector<shared_ptr<Tensor>>{I745, v2_, t2};
  auto task625 = make_shared<Task625>(tensor625, cindex);
  task602->add_dep(task625);
  task625->add_dep(task539);
  deciq->add_task(task625);

  auto tensor626 = vector<shared_ptr<Tensor>>{I745, v2_, t2};
  auto task626 = make_shared<Task626>(tensor626, cindex);
  task602->add_dep(task626);
  task626->add_dep(task539);
  deciq->add_task(task626);

  auto tensor627 = vector<shared_ptr<Tensor>>{I745, v2_, t2};
  auto task627 = make_shared<Task627>(tensor627, cindex);
  task602->add_dep(task627);
  task627->add_dep(task539);
  deciq->add_task(task627);

  vector<IndexRange> I769_index = {active_, active_, active_, active_};
  auto I769 = make_shared<Tensor>(I769_index);
  auto tensor628 = vector<shared_ptr<Tensor>>{I684, Gamma270_(), I769};
  auto task628 = make_shared<Task628>(tensor628, cindex);
  task540->add_dep(task628);
  task628->add_dep(task539);
  deciq->add_task(task628);

  vector<IndexRange> I770_index = {active_, closed_, virt_, active_};
  auto I770 = make_shared<Tensor>(I770_index);
  auto tensor629 = vector<shared_ptr<Tensor>>{I769, t2, I770};
  auto task629 = make_shared<Task629>(tensor629, cindex);
  task628->add_dep(task629);
  task629->add_dep(task539);
  deciq->add_task(task629);

  auto tensor630 = vector<shared_ptr<Tensor>>{I770, f1_, t2};
  auto task630 = make_shared<Task630>(tensor630, cindex);
  task629->add_dep(task630);
  task630->add_dep(task539);
  deciq->add_task(task630);

  auto tensor631 = vector<shared_ptr<Tensor>>{I769, v2_, t2};
  auto task631 = make_shared<Task631>(tensor631, cindex);
  task628->add_dep(task631);
  task631->add_dep(task539);
  deciq->add_task(task631);

  vector<IndexRange> I793_index = {active_, active_, active_, active_, active_, active_};
  auto I793 = make_shared<Tensor>(I793_index);
  auto tensor632 = vector<shared_ptr<Tensor>>{I684, Gamma276_(), I793};
  auto task632 = make_shared<Task632>(tensor632, cindex);
  task540->add_dep(task632);
  task632->add_dep(task539);
  deciq->add_task(task632);

  vector<IndexRange> I794_index = {active_, closed_, active_, active_};
  auto I794 = make_shared<Tensor>(I794_index);
  auto tensor633 = vector<shared_ptr<Tensor>>{I793, t2, I794};
  auto task633 = make_shared<Task633>(tensor633, cindex);
  task632->add_dep(task633);
  task633->add_dep(task539);
  deciq->add_task(task633);

  auto tensor634 = vector<shared_ptr<Tensor>>{I794, f1_, t2};
  auto task634 = make_shared<Task634>(tensor634, cindex);
  task633->add_dep(task634);
  task634->add_dep(task539);
  deciq->add_task(task634);

  vector<IndexRange> I797_index = {active_, active_, active_, active_};
  auto I797 = make_shared<Tensor>(I797_index);
  auto tensor635 = vector<shared_ptr<Tensor>>{I684, Gamma277_(), I797};
  auto task635 = make_shared<Task635>(tensor635, cindex);
  task540->add_dep(task635);
  task635->add_dep(task539);
  deciq->add_task(task635);

  vector<IndexRange> I798_index = {active_, virt_, closed_, active_};
  auto I798 = make_shared<Tensor>(I798_index);
  auto tensor636 = vector<shared_ptr<Tensor>>{I797, t2, I798};
  auto task636 = make_shared<Task636>(tensor636, cindex);
  task635->add_dep(task636);
  task636->add_dep(task539);
  deciq->add_task(task636);

  auto tensor637 = vector<shared_ptr<Tensor>>{I798, t2, f1_};
  auto task637 = make_shared<Task637>(tensor637, cindex);
  task636->add_dep(task637);
  task637->add_dep(task539);
  deciq->add_task(task637);

  auto tensor638 = vector<shared_ptr<Tensor>>{I797, v2_, t2};
  auto task638 = make_shared<Task638>(tensor638, cindex);
  task635->add_dep(task638);
  task638->add_dep(task539);
  deciq->add_task(task638);

  vector<IndexRange> I805_index = {active_, active_, active_, active_};
  auto I805 = make_shared<Tensor>(I805_index);
  auto tensor639 = vector<shared_ptr<Tensor>>{I684, Gamma279_(), I805};
  auto task639 = make_shared<Task639>(tensor639, cindex);
  task540->add_dep(task639);
  task639->add_dep(task539);
  deciq->add_task(task639);

  auto tensor640 = vector<shared_ptr<Tensor>>{I805, t2};
  auto task640 = make_shared<Task640>(tensor640, cindex);
  task639->add_dep(task640);
  task640->add_dep(task539);
  deciq->add_task(task640);

  vector<IndexRange> I808_index = {active_, active_, active_, active_};
  auto I808 = make_shared<Tensor>(I808_index);
  auto tensor641 = vector<shared_ptr<Tensor>>{I684, Gamma280_(), I808};
  auto task641 = make_shared<Task641>(tensor641, cindex);
  task540->add_dep(task641);
  task641->add_dep(task539);
  deciq->add_task(task641);

  vector<IndexRange> I809_index = {active_, virt_, active_, closed_};
  auto I809 = make_shared<Tensor>(I809_index);
  auto tensor642 = vector<shared_ptr<Tensor>>{I808, t2, I809};
  auto task642 = make_shared<Task642>(tensor642, cindex);
  task641->add_dep(task642);
  task642->add_dep(task539);
  deciq->add_task(task642);

  auto tensor643 = vector<shared_ptr<Tensor>>{I809, f1_, t2};
  auto task643 = make_shared<Task643>(tensor643, cindex);
  task642->add_dep(task643);
  task643->add_dep(task539);
  deciq->add_task(task643);

  vector<IndexRange> I813_index = {active_, closed_, active_, virt_};
  auto I813 = make_shared<Tensor>(I813_index);
  auto tensor644 = vector<shared_ptr<Tensor>>{I808, t2, I813};
  auto task644 = make_shared<Task644>(tensor644, cindex);
  task641->add_dep(task644);
  task644->add_dep(task539);
  deciq->add_task(task644);

  auto tensor645 = vector<shared_ptr<Tensor>>{I813, f1_, t2};
  auto task645 = make_shared<Task645>(tensor645, cindex);
  task644->add_dep(task645);
  task645->add_dep(task539);
  deciq->add_task(task645);

  vector<IndexRange> I844_index = {active_, active_, virt_, closed_};
  auto I844 = make_shared<Tensor>(I844_index);
  auto tensor646 = vector<shared_ptr<Tensor>>{I808, t2, I844};
  auto task646 = make_shared<Task646>(tensor646, cindex);
  task641->add_dep(task646);
  task646->add_dep(task539);
  deciq->add_task(task646);

  auto tensor647 = vector<shared_ptr<Tensor>>{I844, t2, f1_};
  auto task647 = make_shared<Task647>(tensor647, cindex);
  task646->add_dep(task647);
  task647->add_dep(task539);
  deciq->add_task(task647);

  vector<IndexRange> I971_index = {closed_, virt_, active_, active_};
  auto I971 = make_shared<Tensor>(I971_index);
  auto tensor648 = vector<shared_ptr<Tensor>>{I808, t2, I971};
  auto task648 = make_shared<Task648>(tensor648, cindex);
  task641->add_dep(task648);
  task648->add_dep(task539);
  deciq->add_task(task648);

  auto tensor649 = vector<shared_ptr<Tensor>>{I971, t2};
  auto task649 = make_shared<Task649>(tensor649, cindex, this->e0_);
  task648->add_dep(task649);
  task649->add_dep(task539);
  deciq->add_task(task649);

  auto tensor650 = vector<shared_ptr<Tensor>>{I971, f1_, t2};
  auto task650 = make_shared<Task650>(tensor650, cindex);
  task648->add_dep(task650);
  task650->add_dep(task539);
  deciq->add_task(task650);

  auto tensor651 = vector<shared_ptr<Tensor>>{I808, v2_, t2};
  auto task651 = make_shared<Task651>(tensor651, cindex);
  task641->add_dep(task651);
  task651->add_dep(task539);
  deciq->add_task(task651);

  auto tensor652 = vector<shared_ptr<Tensor>>{I808, v2_, t2};
  auto task652 = make_shared<Task652>(tensor652, cindex);
  task641->add_dep(task652);
  task652->add_dep(task539);
  deciq->add_task(task652);

  vector<IndexRange> I816_index = {active_, active_, active_, active_};
  auto I816 = make_shared<Tensor>(I816_index);
  auto tensor653 = vector<shared_ptr<Tensor>>{I684, Gamma282_(), I816};
  auto task653 = make_shared<Task653>(tensor653, cindex);
  task540->add_dep(task653);
  task653->add_dep(task539);
  deciq->add_task(task653);

  auto tensor654 = vector<shared_ptr<Tensor>>{I816, t2};
  auto task654 = make_shared<Task654>(tensor654, cindex);
  task653->add_dep(task654);
  task654->add_dep(task539);
  deciq->add_task(task654);

  auto tensor655 = vector<shared_ptr<Tensor>>{I816, t2};
  auto task655 = make_shared<Task655>(tensor655, cindex);
  task653->add_dep(task655);
  task655->add_dep(task539);
  deciq->add_task(task655);

  auto tensor656 = vector<shared_ptr<Tensor>>{I816, t2};
  auto task656 = make_shared<Task656>(tensor656, cindex);
  task653->add_dep(task656);
  task656->add_dep(task539);
  deciq->add_task(task656);

  vector<IndexRange> I819_index = {active_, active_, active_, active_};
  auto I819 = make_shared<Tensor>(I819_index);
  auto tensor657 = vector<shared_ptr<Tensor>>{I684, Gamma283_(), I819};
  auto task657 = make_shared<Task657>(tensor657, cindex);
  task540->add_dep(task657);
  task657->add_dep(task539);
  deciq->add_task(task657);

  vector<IndexRange> I820_index = {active_, virt_, active_, closed_};
  auto I820 = make_shared<Tensor>(I820_index);
  auto tensor658 = vector<shared_ptr<Tensor>>{I819, t2, I820};
  auto task658 = make_shared<Task658>(tensor658, cindex);
  task657->add_dep(task658);
  task658->add_dep(task539);
  deciq->add_task(task658);

  auto tensor659 = vector<shared_ptr<Tensor>>{I820, f1_, t2};
  auto task659 = make_shared<Task659>(tensor659, cindex);
  task658->add_dep(task659);
  task659->add_dep(task539);
  deciq->add_task(task659);

  vector<IndexRange> I824_index = {active_, closed_, active_, virt_};
  auto I824 = make_shared<Tensor>(I824_index);
  auto tensor660 = vector<shared_ptr<Tensor>>{I819, t2, I824};
  auto task660 = make_shared<Task660>(tensor660, cindex);
  task657->add_dep(task660);
  task660->add_dep(task539);
  deciq->add_task(task660);

  auto tensor661 = vector<shared_ptr<Tensor>>{I824, f1_, t2};
  auto task661 = make_shared<Task661>(tensor661, cindex);
  task660->add_dep(task661);
  task661->add_dep(task539);
  deciq->add_task(task661);

  auto tensor662 = vector<shared_ptr<Tensor>>{I824, f1_, t2};
  auto task662 = make_shared<Task662>(tensor662, cindex);
  task660->add_dep(task662);
  task662->add_dep(task539);
  deciq->add_task(task662);

  vector<IndexRange> I840_index = {active_, active_, closed_, virt_};
  auto I840 = make_shared<Tensor>(I840_index);
  auto tensor663 = vector<shared_ptr<Tensor>>{I819, t2, I840};
  auto task663 = make_shared<Task663>(tensor663, cindex);
  task657->add_dep(task663);
  task663->add_dep(task539);
  deciq->add_task(task663);

  auto tensor664 = vector<shared_ptr<Tensor>>{I840, t2, f1_};
  auto task664 = make_shared<Task664>(tensor664, cindex);
  task663->add_dep(task664);
  task664->add_dep(task539);
  deciq->add_task(task664);

  vector<IndexRange> I863_index = {active_, active_, virt_, closed_};
  auto I863 = make_shared<Tensor>(I863_index);
  auto tensor665 = vector<shared_ptr<Tensor>>{I819, t2, I863};
  auto task665 = make_shared<Task665>(tensor665, cindex);
  task657->add_dep(task665);
  task665->add_dep(task539);
  deciq->add_task(task665);

  auto tensor666 = vector<shared_ptr<Tensor>>{I863, f1_, t2};
  auto task666 = make_shared<Task666>(tensor666, cindex);
  task665->add_dep(task666);
  task666->add_dep(task539);
  deciq->add_task(task666);

  vector<IndexRange> I867_index = {active_, active_, closed_, virt_};
  auto I867 = make_shared<Tensor>(I867_index);
  auto tensor667 = vector<shared_ptr<Tensor>>{I819, t2, I867};
  auto task667 = make_shared<Task667>(tensor667, cindex);
  task657->add_dep(task667);
  task667->add_dep(task539);
  deciq->add_task(task667);

  auto tensor668 = vector<shared_ptr<Tensor>>{I867, f1_, t2};
  auto task668 = make_shared<Task668>(tensor668, cindex);
  task667->add_dep(task668);
  task668->add_dep(task539);
  deciq->add_task(task668);

  vector<IndexRange> I874_index = {active_, active_, virt_, closed_};
  auto I874 = make_shared<Tensor>(I874_index);
  auto tensor669 = vector<shared_ptr<Tensor>>{I819, t2, I874};
  auto task669 = make_shared<Task669>(tensor669, cindex);
  task657->add_dep(task669);
  task669->add_dep(task539);
  deciq->add_task(task669);

  auto tensor670 = vector<shared_ptr<Tensor>>{I874, f1_, t2};
  auto task670 = make_shared<Task670>(tensor670, cindex);
  task669->add_dep(task670);
  task670->add_dep(task539);
  deciq->add_task(task670);

  vector<IndexRange> I878_index = {active_, active_, closed_, virt_};
  auto I878 = make_shared<Tensor>(I878_index);
  auto tensor671 = vector<shared_ptr<Tensor>>{I819, t2, I878};
  auto task671 = make_shared<Task671>(tensor671, cindex);
  task657->add_dep(task671);
  task671->add_dep(task539);
  deciq->add_task(task671);

  auto tensor672 = vector<shared_ptr<Tensor>>{I878, f1_, t2};
  auto task672 = make_shared<Task672>(tensor672, cindex);
  task671->add_dep(task672);
  task672->add_dep(task539);
  deciq->add_task(task672);

  vector<IndexRange> I894_index = {active_, active_, closed_, virt_};
  auto I894 = make_shared<Tensor>(I894_index);
  auto tensor673 = vector<shared_ptr<Tensor>>{I819, t2, I894};
  auto task673 = make_shared<Task673>(tensor673, cindex);
  task657->add_dep(task673);
  task673->add_dep(task539);
  deciq->add_task(task673);

  auto tensor674 = vector<shared_ptr<Tensor>>{I894, t2, f1_};
  auto task674 = make_shared<Task674>(tensor674, cindex);
  task673->add_dep(task674);
  task674->add_dep(task539);
  deciq->add_task(task674);

  auto tensor675 = vector<shared_ptr<Tensor>>{I894, t2, f1_};
  auto task675 = make_shared<Task675>(tensor675, cindex);
  task673->add_dep(task675);
  task675->add_dep(task539);
  deciq->add_task(task675);

  vector<IndexRange> I967_index = {virt_, closed_, active_, active_};
  auto I967 = make_shared<Tensor>(I967_index);
  auto tensor676 = vector<shared_ptr<Tensor>>{I819, t2, I967};
  auto task676 = make_shared<Task676>(tensor676, cindex);
  task657->add_dep(task676);
  task676->add_dep(task539);
  deciq->add_task(task676);

  auto tensor677 = vector<shared_ptr<Tensor>>{I967, f1_, t2};
  auto task677 = make_shared<Task677>(tensor677, cindex);
  task676->add_dep(task677);
  task677->add_dep(task539);
  deciq->add_task(task677);

  vector<IndexRange> I979_index = {closed_, virt_, active_, active_};
  auto I979 = make_shared<Tensor>(I979_index);
  auto tensor678 = vector<shared_ptr<Tensor>>{I819, t2, I979};
  auto task678 = make_shared<Task678>(tensor678, cindex);
  task657->add_dep(task678);
  task678->add_dep(task539);
  deciq->add_task(task678);

  auto tensor679 = vector<shared_ptr<Tensor>>{I979, t2};
  auto task679 = make_shared<Task679>(tensor679, cindex, this->e0_);
  task678->add_dep(task679);
  task679->add_dep(task539);
  deciq->add_task(task679);

  auto tensor680 = vector<shared_ptr<Tensor>>{I979, f1_, t2};
  auto task680 = make_shared<Task680>(tensor680, cindex);
  task678->add_dep(task680);
  task680->add_dep(task539);
  deciq->add_task(task680);

  auto tensor681 = vector<shared_ptr<Tensor>>{I819, t2};
  auto task681 = make_shared<Task681>(tensor681, cindex, this->e0_);
  task657->add_dep(task681);
  task681->add_dep(task539);
  deciq->add_task(task681);

  auto tensor682 = vector<shared_ptr<Tensor>>{I819, t2};
  auto task682 = make_shared<Task682>(tensor682, cindex, this->e0_);
  task657->add_dep(task682);
  task682->add_dep(task539);
  deciq->add_task(task682);

  auto tensor683 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task683 = make_shared<Task683>(tensor683, cindex);
  task657->add_dep(task683);
  task683->add_dep(task539);
  deciq->add_task(task683);

  auto tensor684 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task684 = make_shared<Task684>(tensor684, cindex);
  task657->add_dep(task684);
  task684->add_dep(task539);
  deciq->add_task(task684);

  auto tensor685 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task685 = make_shared<Task685>(tensor685, cindex);
  task657->add_dep(task685);
  task685->add_dep(task539);
  deciq->add_task(task685);

  auto tensor686 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task686 = make_shared<Task686>(tensor686, cindex);
  task657->add_dep(task686);
  task686->add_dep(task539);
  deciq->add_task(task686);

  auto tensor687 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task687 = make_shared<Task687>(tensor687, cindex);
  task657->add_dep(task687);
  task687->add_dep(task539);
  deciq->add_task(task687);

  auto tensor688 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task688 = make_shared<Task688>(tensor688, cindex);
  task657->add_dep(task688);
  task688->add_dep(task539);
  deciq->add_task(task688);

  auto tensor689 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task689 = make_shared<Task689>(tensor689, cindex);
  task657->add_dep(task689);
  task689->add_dep(task539);
  deciq->add_task(task689);

  auto tensor690 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task690 = make_shared<Task690>(tensor690, cindex);
  task657->add_dep(task690);
  task690->add_dep(task539);
  deciq->add_task(task690);

  auto tensor691 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task691 = make_shared<Task691>(tensor691, cindex);
  task657->add_dep(task691);
  task691->add_dep(task539);
  deciq->add_task(task691);

  auto tensor692 = vector<shared_ptr<Tensor>>{I819, v2_, t2};
  auto task692 = make_shared<Task692>(tensor692, cindex);
  task657->add_dep(task692);
  task692->add_dep(task539);
  deciq->add_task(task692);

  vector<IndexRange> I827_index = {active_, active_, active_, active_, active_, active_};
  auto I827 = make_shared<Tensor>(I827_index);
  auto tensor693 = vector<shared_ptr<Tensor>>{I684, Gamma285_(), I827};
  auto task693 = make_shared<Task693>(tensor693, cindex);
  task540->add_dep(task693);
  task693->add_dep(task539);
  deciq->add_task(task693);

  vector<IndexRange> I828_index = {active_, virt_, active_, active_};
  auto I828 = make_shared<Tensor>(I828_index);
  auto tensor694 = vector<shared_ptr<Tensor>>{I827, t2, I828};
  auto task694 = make_shared<Task694>(tensor694, cindex);
  task693->add_dep(task694);
  task694->add_dep(task539);
  deciq->add_task(task694);

  auto tensor695 = vector<shared_ptr<Tensor>>{I828, f1_, t2};
  auto task695 = make_shared<Task695>(tensor695, cindex);
  task694->add_dep(task695);
  task695->add_dep(task539);
  deciq->add_task(task695);

  vector<IndexRange> I831_index = {active_, active_};
  auto I831 = make_shared<Tensor>(I831_index);
  auto tensor696 = vector<shared_ptr<Tensor>>{I684, Gamma286_(), I831};
  auto task696 = make_shared<Task696>(tensor696, cindex);
  task540->add_dep(task696);
  task696->add_dep(task539);
  deciq->add_task(task696);

  vector<IndexRange> I832_index = {closed_, virt_};
  auto I832 = make_shared<Tensor>(I832_index);
  auto tensor697 = vector<shared_ptr<Tensor>>{I831, t2, I832};
  auto task697 = make_shared<Task697>(tensor697, cindex);
  task696->add_dep(task697);
  task697->add_dep(task539);
  deciq->add_task(task697);

  auto tensor698 = vector<shared_ptr<Tensor>>{I832, t2, f1_};
  auto task698 = make_shared<Task698>(tensor698, cindex);
  task697->add_dep(task698);
  task698->add_dep(task539);
  deciq->add_task(task698);

  auto tensor699 = vector<shared_ptr<Tensor>>{I832, t2, f1_};
  auto task699 = make_shared<Task699>(tensor699, cindex);
  task697->add_dep(task699);
  task699->add_dep(task539);
  deciq->add_task(task699);

  vector<IndexRange> I886_index = {closed_, virt_};
  auto I886 = make_shared<Tensor>(I886_index);
  auto tensor700 = vector<shared_ptr<Tensor>>{I831, t2, I886};
  auto task700 = make_shared<Task700>(tensor700, cindex);
  task696->add_dep(task700);
  task700->add_dep(task539);
  deciq->add_task(task700);

  auto tensor701 = vector<shared_ptr<Tensor>>{I886, t2, f1_};
  auto task701 = make_shared<Task701>(tensor701, cindex);
  task700->add_dep(task701);
  task701->add_dep(task539);
  deciq->add_task(task701);

  auto tensor702 = vector<shared_ptr<Tensor>>{I886, t2, f1_};
  auto task702 = make_shared<Task702>(tensor702, cindex);
  task700->add_dep(task702);
  task702->add_dep(task539);
  deciq->add_task(task702);

  vector<IndexRange> I937_index = {virt_, closed_};
  auto I937 = make_shared<Tensor>(I937_index);
  auto tensor703 = vector<shared_ptr<Tensor>>{I831, t2, I937};
  auto task703 = make_shared<Task703>(tensor703, cindex);
  task696->add_dep(task703);
  task703->add_dep(task539);
  deciq->add_task(task703);

  auto tensor704 = vector<shared_ptr<Tensor>>{I937, f1_, t2};
  auto task704 = make_shared<Task704>(tensor704, cindex);
  task703->add_dep(task704);
  task704->add_dep(task539);
  deciq->add_task(task704);

  vector<IndexRange> I941_index = {virt_, closed_};
  auto I941 = make_shared<Tensor>(I941_index);
  auto tensor705 = vector<shared_ptr<Tensor>>{I831, t2, I941};
  auto task705 = make_shared<Task705>(tensor705, cindex);
  task696->add_dep(task705);
  task705->add_dep(task539);
  deciq->add_task(task705);

  auto tensor706 = vector<shared_ptr<Tensor>>{I941, f1_, t2};
  auto task706 = make_shared<Task706>(tensor706, cindex);
  task705->add_dep(task706);
  task706->add_dep(task539);
  deciq->add_task(task706);

  vector<IndexRange> I945_index = {virt_, closed_};
  auto I945 = make_shared<Tensor>(I945_index);
  auto tensor707 = vector<shared_ptr<Tensor>>{I831, t2, I945};
  auto task707 = make_shared<Task707>(tensor707, cindex);
  task696->add_dep(task707);
  task707->add_dep(task539);
  deciq->add_task(task707);

  auto tensor708 = vector<shared_ptr<Tensor>>{I945, f1_, t2};
  auto task708 = make_shared<Task708>(tensor708, cindex);
  task707->add_dep(task708);
  task708->add_dep(task539);
  deciq->add_task(task708);

  vector<IndexRange> I949_index = {virt_, closed_};
  auto I949 = make_shared<Tensor>(I949_index);
  auto tensor709 = vector<shared_ptr<Tensor>>{I831, t2, I949};
  auto task709 = make_shared<Task709>(tensor709, cindex);
  task696->add_dep(task709);
  task709->add_dep(task539);
  deciq->add_task(task709);

  auto tensor710 = vector<shared_ptr<Tensor>>{I949, f1_, t2};
  auto task710 = make_shared<Task710>(tensor710, cindex);
  task709->add_dep(task710);
  task710->add_dep(task539);
  deciq->add_task(task710);

  vector<IndexRange> I959_index = {closed_, active_};
  auto I959 = make_shared<Tensor>(I959_index);
  auto tensor711 = vector<shared_ptr<Tensor>>{I831, f1_, I959};
  auto task711 = make_shared<Task711>(tensor711, cindex);
  task696->add_dep(task711);
  task711->add_dep(task539);
  deciq->add_task(task711);

  auto tensor712 = vector<shared_ptr<Tensor>>{I959, t2};
  auto task712 = make_shared<Task712>(tensor712, cindex);
  task711->add_dep(task712);
  task712->add_dep(task539);
  deciq->add_task(task712);

  auto tensor713 = vector<shared_ptr<Tensor>>{I959, t2};
  auto task713 = make_shared<Task713>(tensor713, cindex);
  task711->add_dep(task713);
  task713->add_dep(task539);
  deciq->add_task(task713);

  vector<IndexRange> I991_index = {active_, closed_};
  auto I991 = make_shared<Tensor>(I991_index);
  auto tensor714 = vector<shared_ptr<Tensor>>{I831, f1_, I991};
  auto task714 = make_shared<Task714>(tensor714, cindex);
  task696->add_dep(task714);
  task714->add_dep(task539);
  deciq->add_task(task714);

  auto tensor715 = vector<shared_ptr<Tensor>>{I991, t2};
  auto task715 = make_shared<Task715>(tensor715, cindex);
  task714->add_dep(task715);
  task715->add_dep(task539);
  deciq->add_task(task715);

  auto tensor716 = vector<shared_ptr<Tensor>>{I991, t2};
  auto task716 = make_shared<Task716>(tensor716, cindex);
  task714->add_dep(task716);
  task716->add_dep(task539);
  deciq->add_task(task716);

  vector<IndexRange> I1005_index = {virt_, virt_, active_, closed_};
  auto I1005 = make_shared<Tensor>(I1005_index);
  auto tensor717 = vector<shared_ptr<Tensor>>{I831, t2, I1005};
  auto task717 = make_shared<Task717>(tensor717, cindex);
  task696->add_dep(task717);
  task717->add_dep(task539);
  deciq->add_task(task717);

  auto tensor718 = vector<shared_ptr<Tensor>>{I1005, f1_, t2};
  auto task718 = make_shared<Task718>(tensor718, cindex);
  task717->add_dep(task718);
  task718->add_dep(task539);
  deciq->add_task(task718);

  vector<IndexRange> I1009_index = {virt_, virt_, active_, closed_};
  auto I1009 = make_shared<Tensor>(I1009_index);
  auto tensor719 = vector<shared_ptr<Tensor>>{I831, t2, I1009};
  auto task719 = make_shared<Task719>(tensor719, cindex);
  task696->add_dep(task719);
  task719->add_dep(task539);
  deciq->add_task(task719);

  auto tensor720 = vector<shared_ptr<Tensor>>{I1009, f1_, t2};
  auto task720 = make_shared<Task720>(tensor720, cindex);
  task719->add_dep(task720);
  task720->add_dep(task539);
  deciq->add_task(task720);

  vector<IndexRange> I1013_index = {virt_, closed_, active_, virt_};
  auto I1013 = make_shared<Tensor>(I1013_index);
  auto tensor721 = vector<shared_ptr<Tensor>>{I831, t2, I1013};
  auto task721 = make_shared<Task721>(tensor721, cindex);
  task696->add_dep(task721);
  task721->add_dep(task539);
  deciq->add_task(task721);

  auto tensor722 = vector<shared_ptr<Tensor>>{I1013, f1_, t2};
  auto task722 = make_shared<Task722>(tensor722, cindex);
  task721->add_dep(task722);
  task722->add_dep(task539);
  deciq->add_task(task722);

  vector<IndexRange> I1017_index = {virt_, closed_, active_, virt_};
  auto I1017 = make_shared<Tensor>(I1017_index);
  auto tensor723 = vector<shared_ptr<Tensor>>{I831, t2, I1017};
  auto task723 = make_shared<Task723>(tensor723, cindex);
  task696->add_dep(task723);
  task723->add_dep(task539);
  deciq->add_task(task723);

  auto tensor724 = vector<shared_ptr<Tensor>>{I1017, f1_, t2};
  auto task724 = make_shared<Task724>(tensor724, cindex);
  task723->add_dep(task724);
  task724->add_dep(task539);
  deciq->add_task(task724);

  vector<IndexRange> I1021_index = {closed_, virt_, active_, virt_};
  auto I1021 = make_shared<Tensor>(I1021_index);
  auto tensor725 = vector<shared_ptr<Tensor>>{I831, t2, I1021};
  auto task725 = make_shared<Task725>(tensor725, cindex);
  task696->add_dep(task725);
  task725->add_dep(task539);
  deciq->add_task(task725);

  auto tensor726 = vector<shared_ptr<Tensor>>{I1021, f1_, t2};
  auto task726 = make_shared<Task726>(tensor726, cindex);
  task725->add_dep(task726);
  task726->add_dep(task539);
  deciq->add_task(task726);

  vector<IndexRange> I1025_index = {closed_, virt_, active_, virt_};
  auto I1025 = make_shared<Tensor>(I1025_index);
  auto tensor727 = vector<shared_ptr<Tensor>>{I831, t2, I1025};
  auto task727 = make_shared<Task727>(tensor727, cindex);
  task696->add_dep(task727);
  task727->add_dep(task539);
  deciq->add_task(task727);

  auto tensor728 = vector<shared_ptr<Tensor>>{I1025, f1_, t2};
  auto task728 = make_shared<Task728>(tensor728, cindex);
  task727->add_dep(task728);
  task728->add_dep(task539);
  deciq->add_task(task728);

  auto tensor729 = vector<shared_ptr<Tensor>>{I831, t2};
  auto task729 = make_shared<Task729>(tensor729, cindex, this->e0_);
  task696->add_dep(task729);
  task729->add_dep(task539);
  deciq->add_task(task729);

  auto tensor730 = vector<shared_ptr<Tensor>>{I831, t2};
  auto task730 = make_shared<Task730>(tensor730, cindex, this->e0_);
  task696->add_dep(task730);
  task730->add_dep(task539);
  deciq->add_task(task730);

  auto tensor731 = vector<shared_ptr<Tensor>>{I831, v2_, t2};
  auto task731 = make_shared<Task731>(tensor731, cindex);
  task696->add_dep(task731);
  task731->add_dep(task539);
  deciq->add_task(task731);

  auto tensor732 = vector<shared_ptr<Tensor>>{I831, v2_, t2};
  auto task732 = make_shared<Task732>(tensor732, cindex);
  task696->add_dep(task732);
  task732->add_dep(task539);
  deciq->add_task(task732);

  auto tensor733 = vector<shared_ptr<Tensor>>{I831, v2_, t2};
  auto task733 = make_shared<Task733>(tensor733, cindex);
  task696->add_dep(task733);
  task733->add_dep(task539);
  deciq->add_task(task733);

  auto tensor734 = vector<shared_ptr<Tensor>>{I831, v2_, t2};
  auto task734 = make_shared<Task734>(tensor734, cindex);
  task696->add_dep(task734);
  task734->add_dep(task539);
  deciq->add_task(task734);

  vector<IndexRange> I1195_index = {active_, closed_, virt_, active_};
  auto I1195 = make_shared<Tensor>(I1195_index);
  auto tensor735 = vector<shared_ptr<Tensor>>{I831, h1_, I1195};
  auto task735 = make_shared<Task735>(tensor735, cindex);
  task696->add_dep(task735);
  task735->add_dep(task539);
  deciq->add_task(task735);

  auto tensor736 = vector<shared_ptr<Tensor>>{I1195, t2};
  auto task736 = make_shared<Task736>(tensor736, cindex);
  task735->add_dep(task736);
  task736->add_dep(task539);
  deciq->add_task(task736);

  vector<IndexRange> I1198_index = {active_, active_, virt_, closed_};
  auto I1198 = make_shared<Tensor>(I1198_index);
  auto tensor737 = vector<shared_ptr<Tensor>>{I831, h1_, I1198};
  auto task737 = make_shared<Task737>(tensor737, cindex);
  task696->add_dep(task737);
  task737->add_dep(task539);
  deciq->add_task(task737);

  auto tensor738 = vector<shared_ptr<Tensor>>{I1198, t2};
  auto task738 = make_shared<Task738>(tensor738, cindex);
  task737->add_dep(task738);
  task738->add_dep(task539);
  deciq->add_task(task738);
}

#endif
