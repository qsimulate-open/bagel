//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_residualqq.cc
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


#include <src/smith/mrci/MRCI.h>
#include <src/smith/mrci/MRCI_tasks13.h>
#include <src/smith/mrci/MRCI_tasks14.h>
#include <src/smith/mrci/MRCI_tasks15.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void MRCI::MRCI::make_residualq7(shared_ptr<Queue> residualq, shared_ptr<Task> task108, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I162_index = {virt_, closed_, virt_, closed_};
  auto I162 = make_shared<Tensor>(I162_index);
  auto tensor603 = vector<shared_ptr<Tensor>>{r, I162};
  auto task603 = make_shared<Task603>(tensor603, pindex);
  task603->add_dep(task108);
  residualq->add_task(task603);

  vector<IndexRange> I163_index = {virt_, active_};
  auto I163 = make_shared<Tensor>(I163_index);
  auto tensor604 = vector<shared_ptr<Tensor>>{I162, t2, I163};
  auto task604 = make_shared<Task604>(tensor604, pindex);
  task603->add_dep(task604);
  task604->add_dep(task108);
  residualq->add_task(task604);

  auto tensor605 = vector<shared_ptr<Tensor>>{I163, Gamma12_(), h1_};
  auto task605 = make_shared<Task605>(tensor605, pindex);
  task604->add_dep(task605);
  task605->add_dep(task108);
  residualq->add_task(task605);

  vector<IndexRange> I166_index = {virt_, active_};
  auto I166 = make_shared<Tensor>(I166_index);
  auto tensor606 = vector<shared_ptr<Tensor>>{I162, t2, I166};
  auto task606 = make_shared<Task606>(tensor606, pindex);
  task603->add_dep(task606);
  task606->add_dep(task108);
  residualq->add_task(task606);

  auto tensor607 = vector<shared_ptr<Tensor>>{I166, Gamma12_(), h1_};
  auto task607 = make_shared<Task607>(tensor607, pindex);
  task606->add_dep(task607);
  task607->add_dep(task108);
  residualq->add_task(task607);

  vector<IndexRange> I169_index = {virt_, closed_};
  auto I169 = make_shared<Tensor>(I169_index);
  auto tensor608 = vector<shared_ptr<Tensor>>{I162, h1_, I169};
  auto task608 = make_shared<Task608>(tensor608, pindex);
  task603->add_dep(task608);
  task608->add_dep(task108);
  residualq->add_task(task608);

  vector<IndexRange> I170_index = {active_, virt_, closed_, active_};
  auto I170 = make_shared<Tensor>(I170_index);
  auto tensor609 = vector<shared_ptr<Tensor>>{I169, Gamma32_(), I170};
  auto task609 = make_shared<Task609>(tensor609, pindex);
  task608->add_dep(task609);
  task609->add_dep(task108);
  residualq->add_task(task609);

  auto tensor610 = vector<shared_ptr<Tensor>>{I170, t2};
  auto task610 = make_shared<Task610>(tensor610, pindex);
  task609->add_dep(task610);
  task610->add_dep(task108);
  residualq->add_task(task610);

  vector<IndexRange> I172_index = {virt_, closed_};
  auto I172 = make_shared<Tensor>(I172_index);
  auto tensor611 = vector<shared_ptr<Tensor>>{I162, h1_, I172};
  auto task611 = make_shared<Task611>(tensor611, pindex);
  task603->add_dep(task611);
  task611->add_dep(task108);
  residualq->add_task(task611);

  vector<IndexRange> I173_index = {active_, virt_, closed_, active_};
  auto I173 = make_shared<Tensor>(I173_index);
  auto tensor612 = vector<shared_ptr<Tensor>>{I172, Gamma32_(), I173};
  auto task612 = make_shared<Task612>(tensor612, pindex);
  task611->add_dep(task612);
  task612->add_dep(task108);
  residualq->add_task(task612);

  auto tensor613 = vector<shared_ptr<Tensor>>{I173, t2};
  auto task613 = make_shared<Task613>(tensor613, pindex);
  task612->add_dep(task613);
  task613->add_dep(task108);
  residualq->add_task(task613);

  vector<IndexRange> I181_index = {closed_, closed_};
  auto I181 = make_shared<Tensor>(I181_index);
  auto tensor614 = vector<shared_ptr<Tensor>>{I162, t2, I181};
  auto task614 = make_shared<Task614>(tensor614, pindex);
  task603->add_dep(task614);
  task614->add_dep(task108);
  residualq->add_task(task614);

  shared_ptr<Task615> task615;
  if (diagonal) {
    auto tensor615 = vector<shared_ptr<Tensor>>{I181, h1_};
    task615 = make_shared<Task615>(tensor615, pindex);
    task614->add_dep(task615);
    task615->add_dep(task108);
    residualq->add_task(task615);
  }

  vector<IndexRange> I1243_index = {closed_, closed_, active_, active_};
  auto I1243 = make_shared<Tensor>(I1243_index);
  auto tensor616 = vector<shared_ptr<Tensor>>{I181, Gamma32_(), I1243};
  auto task616 = make_shared<Task616>(tensor616, pindex);
  task614->add_dep(task616);
  task616->add_dep(task108);
  residualq->add_task(task616);

  auto tensor617 = vector<shared_ptr<Tensor>>{I1243, v2_};
  auto task617 = make_shared<Task617>(tensor617, pindex);
  task616->add_dep(task617);
  task617->add_dep(task108);
  residualq->add_task(task617);

  auto tensor618 = vector<shared_ptr<Tensor>>{I181, Gamma12_(), v2_};
  auto task618 = make_shared<Task618>(tensor618, pindex);
  task614->add_dep(task618);
  task618->add_dep(task108);
  residualq->add_task(task618);

  vector<IndexRange> I183_index = {closed_, closed_};
  auto I183 = make_shared<Tensor>(I183_index);
  auto tensor619 = vector<shared_ptr<Tensor>>{I162, t2, I183};
  auto task619 = make_shared<Task619>(tensor619, pindex);
  task603->add_dep(task619);
  task619->add_dep(task108);
  residualq->add_task(task619);

  shared_ptr<Task620> task620;
  if (diagonal) {
    auto tensor620 = vector<shared_ptr<Tensor>>{I183, h1_};
    task620 = make_shared<Task620>(tensor620, pindex);
    task619->add_dep(task620);
    task620->add_dep(task108);
    residualq->add_task(task620);
  }

  vector<IndexRange> I1246_index = {closed_, closed_, active_, active_};
  auto I1246 = make_shared<Tensor>(I1246_index);
  auto tensor621 = vector<shared_ptr<Tensor>>{I183, Gamma32_(), I1246};
  auto task621 = make_shared<Task621>(tensor621, pindex);
  task619->add_dep(task621);
  task621->add_dep(task108);
  residualq->add_task(task621);

  auto tensor622 = vector<shared_ptr<Tensor>>{I1246, v2_};
  auto task622 = make_shared<Task622>(tensor622, pindex);
  task621->add_dep(task622);
  task622->add_dep(task108);
  residualq->add_task(task622);

  auto tensor623 = vector<shared_ptr<Tensor>>{I183, Gamma12_(), v2_};
  auto task623 = make_shared<Task623>(tensor623, pindex);
  task619->add_dep(task623);
  task623->add_dep(task108);
  residualq->add_task(task623);

  vector<IndexRange> I185_index = {virt_, virt_};
  auto I185 = make_shared<Tensor>(I185_index);
  auto tensor624 = vector<shared_ptr<Tensor>>{I162, t2, I185};
  auto task624 = make_shared<Task624>(tensor624, pindex);
  task603->add_dep(task624);
  task624->add_dep(task108);
  residualq->add_task(task624);

  shared_ptr<Task625> task625;
  if (diagonal) {
    auto tensor625 = vector<shared_ptr<Tensor>>{I185, h1_};
    task625 = make_shared<Task625>(tensor625, pindex);
    task624->add_dep(task625);
    task625->add_dep(task108);
    residualq->add_task(task625);
  }

  vector<IndexRange> I1249_index = {virt_, virt_, active_, active_};
  auto I1249 = make_shared<Tensor>(I1249_index);
  auto tensor626 = vector<shared_ptr<Tensor>>{I185, Gamma32_(), I1249};
  auto task626 = make_shared<Task626>(tensor626, pindex);
  task624->add_dep(task626);
  task626->add_dep(task108);
  residualq->add_task(task626);

  auto tensor627 = vector<shared_ptr<Tensor>>{I1249, v2_};
  auto task627 = make_shared<Task627>(tensor627, pindex);
  task626->add_dep(task627);
  task627->add_dep(task108);
  residualq->add_task(task627);

  auto tensor628 = vector<shared_ptr<Tensor>>{I185, Gamma12_(), v2_};
  auto task628 = make_shared<Task628>(tensor628, pindex);
  task624->add_dep(task628);
  task628->add_dep(task108);
  residualq->add_task(task628);

  vector<IndexRange> I187_index = {virt_, virt_};
  auto I187 = make_shared<Tensor>(I187_index);
  auto tensor629 = vector<shared_ptr<Tensor>>{I162, t2, I187};
  auto task629 = make_shared<Task629>(tensor629, pindex);
  task603->add_dep(task629);
  task629->add_dep(task108);
  residualq->add_task(task629);

  shared_ptr<Task630> task630;
  if (diagonal) {
    auto tensor630 = vector<shared_ptr<Tensor>>{I187, h1_};
    task630 = make_shared<Task630>(tensor630, pindex);
    task629->add_dep(task630);
    task630->add_dep(task108);
    residualq->add_task(task630);
  }

  vector<IndexRange> I1252_index = {virt_, virt_, active_, active_};
  auto I1252 = make_shared<Tensor>(I1252_index);
  auto tensor631 = vector<shared_ptr<Tensor>>{I187, Gamma32_(), I1252};
  auto task631 = make_shared<Task631>(tensor631, pindex);
  task629->add_dep(task631);
  task631->add_dep(task108);
  residualq->add_task(task631);

  auto tensor632 = vector<shared_ptr<Tensor>>{I1252, v2_};
  auto task632 = make_shared<Task632>(tensor632, pindex);
  task631->add_dep(task632);
  task632->add_dep(task108);
  residualq->add_task(task632);

  auto tensor633 = vector<shared_ptr<Tensor>>{I187, Gamma12_(), v2_};
  auto task633 = make_shared<Task633>(tensor633, pindex);
  task629->add_dep(task633);
  task633->add_dep(task108);
  residualq->add_task(task633);

  vector<IndexRange> I189_index = {closed_, active_};
  auto I189 = make_shared<Tensor>(I189_index);
  auto tensor634 = vector<shared_ptr<Tensor>>{I162, t2, I189};
  auto task634 = make_shared<Task634>(tensor634, pindex);
  task603->add_dep(task634);
  task634->add_dep(task108);
  residualq->add_task(task634);

  auto tensor635 = vector<shared_ptr<Tensor>>{I189, Gamma32_(), h1_};
  auto task635 = make_shared<Task635>(tensor635, pindex);
  task634->add_dep(task635);
  task635->add_dep(task108);
  residualq->add_task(task635);

  vector<IndexRange> I192_index = {closed_, active_};
  auto I192 = make_shared<Tensor>(I192_index);
  auto tensor636 = vector<shared_ptr<Tensor>>{I162, t2, I192};
  auto task636 = make_shared<Task636>(tensor636, pindex);
  task603->add_dep(task636);
  task636->add_dep(task108);
  residualq->add_task(task636);

  auto tensor637 = vector<shared_ptr<Tensor>>{I192, Gamma32_(), h1_};
  auto task637 = make_shared<Task637>(tensor637, pindex);
  task636->add_dep(task637);
  task637->add_dep(task108);
  residualq->add_task(task637);

  vector<IndexRange> I1116_index = {closed_, active_};
  auto I1116 = make_shared<Tensor>(I1116_index);
  auto tensor638 = vector<shared_ptr<Tensor>>{I162, v2_, I1116};
  auto task638 = make_shared<Task638>(tensor638, pindex);
  task603->add_dep(task638);
  task638->add_dep(task108);
  residualq->add_task(task638);

  auto tensor639 = vector<shared_ptr<Tensor>>{I1116, Gamma10_(), t2};
  auto task639 = make_shared<Task639>(tensor639, pindex);
  task638->add_dep(task639);
  task639->add_dep(task108);
  residualq->add_task(task639);

  vector<IndexRange> I1119_index = {closed_, active_};
  auto I1119 = make_shared<Tensor>(I1119_index);
  auto tensor640 = vector<shared_ptr<Tensor>>{I162, v2_, I1119};
  auto task640 = make_shared<Task640>(tensor640, pindex);
  task603->add_dep(task640);
  task640->add_dep(task108);
  residualq->add_task(task640);

  auto tensor641 = vector<shared_ptr<Tensor>>{I1119, Gamma10_(), t2};
  auto task641 = make_shared<Task641>(tensor641, pindex);
  task640->add_dep(task641);
  task641->add_dep(task108);
  residualq->add_task(task641);

  vector<IndexRange> I1122_index = {virt_, active_};
  auto I1122 = make_shared<Tensor>(I1122_index);
  auto tensor642 = vector<shared_ptr<Tensor>>{I162, t2, I1122};
  auto task642 = make_shared<Task642>(tensor642, pindex);
  task603->add_dep(task642);
  task642->add_dep(task108);
  residualq->add_task(task642);

  auto tensor643 = vector<shared_ptr<Tensor>>{I1122, Gamma5_(), v2_};
  auto task643 = make_shared<Task643>(tensor643, pindex);
  task642->add_dep(task643);
  task643->add_dep(task108);
  residualq->add_task(task643);

  auto tensor644 = vector<shared_ptr<Tensor>>{I1122, Gamma197_(), v2_};
  auto task644 = make_shared<Task644>(tensor644, pindex);
  task642->add_dep(task644);
  task644->add_dep(task108);
  residualq->add_task(task644);

  vector<IndexRange> I1125_index = {virt_, active_};
  auto I1125 = make_shared<Tensor>(I1125_index);
  auto tensor645 = vector<shared_ptr<Tensor>>{I162, t2, I1125};
  auto task645 = make_shared<Task645>(tensor645, pindex);
  task603->add_dep(task645);
  task645->add_dep(task108);
  residualq->add_task(task645);

  auto tensor646 = vector<shared_ptr<Tensor>>{I1125, Gamma5_(), v2_};
  auto task646 = make_shared<Task646>(tensor646, pindex);
  task645->add_dep(task646);
  task646->add_dep(task108);
  residualq->add_task(task646);

  auto tensor647 = vector<shared_ptr<Tensor>>{I1125, Gamma197_(), v2_};
  auto task647 = make_shared<Task647>(tensor647, pindex);
  task645->add_dep(task647);
  task647->add_dep(task108);
  residualq->add_task(task647);

  vector<IndexRange> I1134_index = {closed_, closed_, virt_, active_};
  auto I1134 = make_shared<Tensor>(I1134_index);
  auto tensor648 = vector<shared_ptr<Tensor>>{I162, t2, I1134};
  auto task648 = make_shared<Task648>(tensor648, pindex);
  task603->add_dep(task648);
  task648->add_dep(task108);
  residualq->add_task(task648);

  vector<IndexRange> I1135_index = {active_, closed_, closed_, virt_};
  auto I1135 = make_shared<Tensor>(I1135_index);
  auto tensor649 = vector<shared_ptr<Tensor>>{I1134, Gamma12_(), I1135};
  auto task649 = make_shared<Task649>(tensor649, pindex);
  task648->add_dep(task649);
  task649->add_dep(task108);
  residualq->add_task(task649);

  auto tensor650 = vector<shared_ptr<Tensor>>{I1135, v2_};
  auto task650 = make_shared<Task650>(tensor650, pindex);
  task649->add_dep(task650);
  task650->add_dep(task108);
  residualq->add_task(task650);

  vector<IndexRange> I1137_index = {closed_, closed_, virt_, active_};
  auto I1137 = make_shared<Tensor>(I1137_index);
  auto tensor651 = vector<shared_ptr<Tensor>>{I162, t2, I1137};
  auto task651 = make_shared<Task651>(tensor651, pindex);
  task603->add_dep(task651);
  task651->add_dep(task108);
  residualq->add_task(task651);

  vector<IndexRange> I1138_index = {active_, closed_, closed_, virt_};
  auto I1138 = make_shared<Tensor>(I1138_index);
  auto tensor652 = vector<shared_ptr<Tensor>>{I1137, Gamma12_(), I1138};
  auto task652 = make_shared<Task652>(tensor652, pindex);
  task651->add_dep(task652);
  task652->add_dep(task108);
  residualq->add_task(task652);

  auto tensor653 = vector<shared_ptr<Tensor>>{I1138, v2_};
  auto task653 = make_shared<Task653>(tensor653, pindex);
  task652->add_dep(task653);
  task653->add_dep(task108);
  residualq->add_task(task653);

  vector<IndexRange> I1146_index = {closed_, closed_, virt_, active_};
  auto I1146 = make_shared<Tensor>(I1146_index);
  auto tensor654 = vector<shared_ptr<Tensor>>{I162, t2, I1146};
  auto task654 = make_shared<Task654>(tensor654, pindex);
  task603->add_dep(task654);
  task654->add_dep(task108);
  residualq->add_task(task654);

  vector<IndexRange> I1147_index = {active_, closed_, closed_, virt_};
  auto I1147 = make_shared<Tensor>(I1147_index);
  auto tensor655 = vector<shared_ptr<Tensor>>{I1146, Gamma12_(), I1147};
  auto task655 = make_shared<Task655>(tensor655, pindex);
  task654->add_dep(task655);
  task655->add_dep(task108);
  residualq->add_task(task655);

  auto tensor656 = vector<shared_ptr<Tensor>>{I1147, v2_};
  auto task656 = make_shared<Task656>(tensor656, pindex);
  task655->add_dep(task656);
  task656->add_dep(task108);
  residualq->add_task(task656);

  vector<IndexRange> I1149_index = {closed_, closed_, virt_, active_};
  auto I1149 = make_shared<Tensor>(I1149_index);
  auto tensor657 = vector<shared_ptr<Tensor>>{I162, t2, I1149};
  auto task657 = make_shared<Task657>(tensor657, pindex);
  task603->add_dep(task657);
  task657->add_dep(task108);
  residualq->add_task(task657);

  vector<IndexRange> I1150_index = {active_, closed_, closed_, virt_};
  auto I1150 = make_shared<Tensor>(I1150_index);
  auto tensor658 = vector<shared_ptr<Tensor>>{I1149, Gamma12_(), I1150};
  auto task658 = make_shared<Task658>(tensor658, pindex);
  task657->add_dep(task658);
  task658->add_dep(task108);
  residualq->add_task(task658);

  auto tensor659 = vector<shared_ptr<Tensor>>{I1150, v2_};
  auto task659 = make_shared<Task659>(tensor659, pindex);
  task658->add_dep(task659);
  task659->add_dep(task108);
  residualq->add_task(task659);

  vector<IndexRange> I1158_index = {closed_, virt_, closed_, active_};
  auto I1158 = make_shared<Tensor>(I1158_index);
  auto tensor660 = vector<shared_ptr<Tensor>>{I162, v2_, I1158};
  auto task660 = make_shared<Task660>(tensor660, pindex);
  task603->add_dep(task660);
  task660->add_dep(task108);
  residualq->add_task(task660);

  auto tensor661 = vector<shared_ptr<Tensor>>{I1158, Gamma12_(), t2};
  auto task661 = make_shared<Task661>(tensor661, pindex);
  task660->add_dep(task661);
  task661->add_dep(task108);
  residualq->add_task(task661);

  vector<IndexRange> I1161_index = {closed_, virt_, closed_, active_};
  auto I1161 = make_shared<Tensor>(I1161_index);
  auto tensor662 = vector<shared_ptr<Tensor>>{I162, v2_, I1161};
  auto task662 = make_shared<Task662>(tensor662, pindex);
  task603->add_dep(task662);
  task662->add_dep(task108);
  residualq->add_task(task662);

  auto tensor663 = vector<shared_ptr<Tensor>>{I1161, Gamma12_(), t2};
  auto task663 = make_shared<Task663>(tensor663, pindex);
  task662->add_dep(task663);
  task663->add_dep(task108);
  residualq->add_task(task663);

  vector<IndexRange> I1164_index = {closed_, virt_, active_, active_};
  auto I1164 = make_shared<Tensor>(I1164_index);
  auto tensor664 = vector<shared_ptr<Tensor>>{I162, t2, I1164};
  auto task664 = make_shared<Task664>(tensor664, pindex);
  task603->add_dep(task664);
  task664->add_dep(task108);
  residualq->add_task(task664);

  vector<IndexRange> I1165_index = {closed_, virt_, active_, active_};
  auto I1165 = make_shared<Tensor>(I1165_index);
  auto tensor665 = vector<shared_ptr<Tensor>>{I1164, Gamma29_(), I1165};
  auto task665 = make_shared<Task665>(tensor665, pindex);
  task664->add_dep(task665);
  task665->add_dep(task108);
  residualq->add_task(task665);

  auto tensor666 = vector<shared_ptr<Tensor>>{I1165, v2_};
  auto task666 = make_shared<Task666>(tensor666, pindex);
  task665->add_dep(task666);
  task666->add_dep(task108);
  residualq->add_task(task666);

  auto tensor667 = vector<shared_ptr<Tensor>>{I1164, Gamma18_(), v2_};
  auto task667 = make_shared<Task667>(tensor667, pindex);
  task664->add_dep(task667);
  task667->add_dep(task108);
  residualq->add_task(task667);

  auto tensor668 = vector<shared_ptr<Tensor>>{I1164, Gamma27_(), v2_};
  auto task668 = make_shared<Task668>(tensor668, pindex);
  task664->add_dep(task668);
  task668->add_dep(task108);
  residualq->add_task(task668);

  vector<IndexRange> I1167_index = {closed_, virt_, active_, active_};
  auto I1167 = make_shared<Tensor>(I1167_index);
  auto tensor669 = vector<shared_ptr<Tensor>>{I162, t2, I1167};
  auto task669 = make_shared<Task669>(tensor669, pindex);
  task603->add_dep(task669);
  task669->add_dep(task108);
  residualq->add_task(task669);

  vector<IndexRange> I1168_index = {closed_, virt_, active_, active_};
  auto I1168 = make_shared<Tensor>(I1168_index);
  auto tensor670 = vector<shared_ptr<Tensor>>{I1167, Gamma29_(), I1168};
  auto task670 = make_shared<Task670>(tensor670, pindex);
  task669->add_dep(task670);
  task670->add_dep(task108);
  residualq->add_task(task670);

  auto tensor671 = vector<shared_ptr<Tensor>>{I1168, v2_};
  auto task671 = make_shared<Task671>(tensor671, pindex);
  task670->add_dep(task671);
  task671->add_dep(task108);
  residualq->add_task(task671);

  auto tensor672 = vector<shared_ptr<Tensor>>{I1167, Gamma10_(), v2_};
  auto task672 = make_shared<Task672>(tensor672, pindex);
  task669->add_dep(task672);
  task672->add_dep(task108);
  residualq->add_task(task672);

  vector<IndexRange> I1188_index = {virt_, closed_};
  auto I1188 = make_shared<Tensor>(I1188_index);
  auto tensor673 = vector<shared_ptr<Tensor>>{I162, v2_, I1188};
  auto task673 = make_shared<Task673>(tensor673, pindex);
  task603->add_dep(task673);
  task673->add_dep(task108);
  residualq->add_task(task673);

  vector<IndexRange> I1189_index = {active_, virt_, closed_, active_};
  auto I1189 = make_shared<Tensor>(I1189_index);
  auto tensor674 = vector<shared_ptr<Tensor>>{I1188, Gamma32_(), I1189};
  auto task674 = make_shared<Task674>(tensor674, pindex);
  task673->add_dep(task674);
  task674->add_dep(task108);
  residualq->add_task(task674);

  auto tensor675 = vector<shared_ptr<Tensor>>{I1189, t2};
  auto task675 = make_shared<Task675>(tensor675, pindex);
  task674->add_dep(task675);
  task675->add_dep(task108);
  residualq->add_task(task675);

  vector<IndexRange> I1191_index = {virt_, closed_};
  auto I1191 = make_shared<Tensor>(I1191_index);
  auto tensor676 = vector<shared_ptr<Tensor>>{I162, v2_, I1191};
  auto task676 = make_shared<Task676>(tensor676, pindex);
  task603->add_dep(task676);
  task676->add_dep(task108);
  residualq->add_task(task676);

  vector<IndexRange> I1192_index = {active_, virt_, closed_, active_};
  auto I1192 = make_shared<Tensor>(I1192_index);
  auto tensor677 = vector<shared_ptr<Tensor>>{I1191, Gamma32_(), I1192};
  auto task677 = make_shared<Task677>(tensor677, pindex);
  task676->add_dep(task677);
  task677->add_dep(task108);
  residualq->add_task(task677);

  auto tensor678 = vector<shared_ptr<Tensor>>{I1192, t2};
  auto task678 = make_shared<Task678>(tensor678, pindex);
  task677->add_dep(task678);
  task678->add_dep(task108);
  residualq->add_task(task678);

  vector<IndexRange> I1194_index = {virt_, closed_};
  auto I1194 = make_shared<Tensor>(I1194_index);
  auto tensor679 = vector<shared_ptr<Tensor>>{I162, v2_, I1194};
  auto task679 = make_shared<Task679>(tensor679, pindex);
  task603->add_dep(task679);
  task679->add_dep(task108);
  residualq->add_task(task679);

  vector<IndexRange> I1195_index = {active_, virt_, closed_, active_};
  auto I1195 = make_shared<Tensor>(I1195_index);
  auto tensor680 = vector<shared_ptr<Tensor>>{I1194, Gamma32_(), I1195};
  auto task680 = make_shared<Task680>(tensor680, pindex);
  task679->add_dep(task680);
  task680->add_dep(task108);
  residualq->add_task(task680);

  auto tensor681 = vector<shared_ptr<Tensor>>{I1195, t2};
  auto task681 = make_shared<Task681>(tensor681, pindex);
  task680->add_dep(task681);
  task681->add_dep(task108);
  residualq->add_task(task681);

  vector<IndexRange> I1197_index = {virt_, closed_};
  auto I1197 = make_shared<Tensor>(I1197_index);
  auto tensor682 = vector<shared_ptr<Tensor>>{I162, v2_, I1197};
  auto task682 = make_shared<Task682>(tensor682, pindex);
  task603->add_dep(task682);
  task682->add_dep(task108);
  residualq->add_task(task682);

  vector<IndexRange> I1198_index = {active_, virt_, closed_, active_};
  auto I1198 = make_shared<Tensor>(I1198_index);
  auto tensor683 = vector<shared_ptr<Tensor>>{I1197, Gamma32_(), I1198};
  auto task683 = make_shared<Task683>(tensor683, pindex);
  task682->add_dep(task683);
  task683->add_dep(task108);
  residualq->add_task(task683);

  auto tensor684 = vector<shared_ptr<Tensor>>{I1198, t2};
  auto task684 = make_shared<Task684>(tensor684, pindex);
  task683->add_dep(task684);
  task684->add_dep(task108);
  residualq->add_task(task684);

  vector<IndexRange> I1200_index = {closed_, virt_, active_, active_};
  auto I1200 = make_shared<Tensor>(I1200_index);
  auto tensor685 = vector<shared_ptr<Tensor>>{I162, t2, I1200};
  auto task685 = make_shared<Task685>(tensor685, pindex);
  task603->add_dep(task685);
  task685->add_dep(task108);
  residualq->add_task(task685);

  vector<IndexRange> I1201_index = {closed_, virt_, active_, active_};
  auto I1201 = make_shared<Tensor>(I1201_index);
  auto tensor686 = vector<shared_ptr<Tensor>>{I1200, Gamma29_(), I1201};
  auto task686 = make_shared<Task686>(tensor686, pindex);
  task685->add_dep(task686);
  task686->add_dep(task108);
  residualq->add_task(task686);

  auto tensor687 = vector<shared_ptr<Tensor>>{I1201, v2_};
  auto task687 = make_shared<Task687>(tensor687, pindex);
  task686->add_dep(task687);
  task687->add_dep(task108);
  residualq->add_task(task687);

  auto tensor688 = vector<shared_ptr<Tensor>>{I1200, Gamma10_(), v2_};
  auto task688 = make_shared<Task688>(tensor688, pindex);
  task685->add_dep(task688);
  task688->add_dep(task108);
  residualq->add_task(task688);

  vector<IndexRange> I1203_index = {closed_, virt_, active_, active_};
  auto I1203 = make_shared<Tensor>(I1203_index);
  auto tensor689 = vector<shared_ptr<Tensor>>{I162, t2, I1203};
  auto task689 = make_shared<Task689>(tensor689, pindex);
  task603->add_dep(task689);
  task689->add_dep(task108);
  residualq->add_task(task689);

  vector<IndexRange> I1204_index = {closed_, virt_, active_, active_};
  auto I1204 = make_shared<Tensor>(I1204_index);
  auto tensor690 = vector<shared_ptr<Tensor>>{I1203, Gamma29_(), I1204};
  auto task690 = make_shared<Task690>(tensor690, pindex);
  task689->add_dep(task690);
  task690->add_dep(task108);
  residualq->add_task(task690);

  auto tensor691 = vector<shared_ptr<Tensor>>{I1204, v2_};
  auto task691 = make_shared<Task691>(tensor691, pindex);
  task690->add_dep(task691);
  task691->add_dep(task108);
  residualq->add_task(task691);

  auto tensor692 = vector<shared_ptr<Tensor>>{I1203, Gamma10_(), v2_};
  auto task692 = make_shared<Task692>(tensor692, pindex);
  task689->add_dep(task692);
  task692->add_dep(task108);
  residualq->add_task(task692);

  vector<IndexRange> I1236_index = {virt_, active_};
  auto I1236 = make_shared<Tensor>(I1236_index);
  auto tensor693 = vector<shared_ptr<Tensor>>{I162, v2_, I1236};
  auto task693 = make_shared<Task693>(tensor693, pindex);
  task603->add_dep(task693);
  task693->add_dep(task108);
  residualq->add_task(task693);

  auto tensor694 = vector<shared_ptr<Tensor>>{I1236, Gamma51_(), t2};
  auto task694 = make_shared<Task694>(tensor694, pindex);
  task693->add_dep(task694);
  task694->add_dep(task108);
  residualq->add_task(task694);

  vector<IndexRange> I1239_index = {virt_, active_};
  auto I1239 = make_shared<Tensor>(I1239_index);
  auto tensor695 = vector<shared_ptr<Tensor>>{I162, v2_, I1239};
  auto task695 = make_shared<Task695>(tensor695, pindex);
  task603->add_dep(task695);
  task695->add_dep(task108);
  residualq->add_task(task695);

  auto tensor696 = vector<shared_ptr<Tensor>>{I1239, Gamma51_(), t2};
  auto task696 = make_shared<Task696>(tensor696, pindex);
  task695->add_dep(task696);
  task696->add_dep(task108);
  residualq->add_task(task696);

  shared_ptr<Task697> task697;
  if (diagonal) {
    auto tensor697 = vector<shared_ptr<Tensor>>{I162, t2, v2_};
    task697 = make_shared<Task697>(tensor697, pindex);
    task603->add_dep(task697);
    task697->add_dep(task108);
    residualq->add_task(task697);
  }

  shared_ptr<Task698> task698;
  if (diagonal) {
    auto tensor698 = vector<shared_ptr<Tensor>>{I162, t2, v2_};
    task698 = make_shared<Task698>(tensor698, pindex);
    task603->add_dep(task698);
    task698->add_dep(task108);
    residualq->add_task(task698);
  }

  shared_ptr<Task699> task699;
  if (diagonal) {
    auto tensor699 = vector<shared_ptr<Tensor>>{I162, t2, v2_};
    task699 = make_shared<Task699>(tensor699, pindex);
    task603->add_dep(task699);
    task699->add_dep(task108);
    residualq->add_task(task699);
  }

  shared_ptr<Task700> task700;
  if (diagonal) {
    auto tensor700 = vector<shared_ptr<Tensor>>{I162, t2, v2_};
    task700 = make_shared<Task700>(tensor700, pindex);
    task603->add_dep(task700);
    task700->add_dep(task108);
    residualq->add_task(task700);
  }

  shared_ptr<Task701> task701;
  if (diagonal) {
    auto tensor701 = vector<shared_ptr<Tensor>>{I162, t2, v2_};
    task701 = make_shared<Task701>(tensor701, pindex);
    task603->add_dep(task701);
    task701->add_dep(task108);
    residualq->add_task(task701);
  }

  shared_ptr<Task702> task702;
  if (diagonal) {
    auto tensor702 = vector<shared_ptr<Tensor>>{I162, t2, v2_};
    task702 = make_shared<Task702>(tensor702, pindex);
    task603->add_dep(task702);
    task702->add_dep(task108);
    residualq->add_task(task702);
  }

  shared_ptr<Task703> task703;
  if (diagonal) {
    auto tensor703 = vector<shared_ptr<Tensor>>{I162, t2, v2_};
    task703 = make_shared<Task703>(tensor703, pindex);
    task603->add_dep(task703);
    task703->add_dep(task108);
    residualq->add_task(task703);
  }

  shared_ptr<Task704> task704;
  if (diagonal) {
    auto tensor704 = vector<shared_ptr<Tensor>>{I162, t2, v2_};
    task704 = make_shared<Task704>(tensor704, pindex);
    task603->add_dep(task704);
    task704->add_dep(task108);
    residualq->add_task(task704);
  }

  vector<IndexRange> I1310_index = {closed_, active_};
  auto I1310 = make_shared<Tensor>(I1310_index);
  auto tensor705 = vector<shared_ptr<Tensor>>{I162, t2, I1310};
  auto task705 = make_shared<Task705>(tensor705, pindex);
  task603->add_dep(task705);
  task705->add_dep(task108);
  residualq->add_task(task705);

  auto tensor706 = vector<shared_ptr<Tensor>>{I1310, Gamma29_(), v2_};
  auto task706 = make_shared<Task706>(tensor706, pindex);
  task705->add_dep(task706);
  task706->add_dep(task108);
  residualq->add_task(task706);

  auto tensor707 = vector<shared_ptr<Tensor>>{I1310, Gamma51_(), v2_};
  auto task707 = make_shared<Task707>(tensor707, pindex);
  task705->add_dep(task707);
  task707->add_dep(task108);
  residualq->add_task(task707);

  vector<IndexRange> I1313_index = {closed_, active_};
  auto I1313 = make_shared<Tensor>(I1313_index);
  auto tensor708 = vector<shared_ptr<Tensor>>{I162, t2, I1313};
  auto task708 = make_shared<Task708>(tensor708, pindex);
  task603->add_dep(task708);
  task708->add_dep(task108);
  residualq->add_task(task708);

  auto tensor709 = vector<shared_ptr<Tensor>>{I1313, Gamma29_(), v2_};
  auto task709 = make_shared<Task709>(tensor709, pindex);
  task708->add_dep(task709);
  task709->add_dep(task108);
  residualq->add_task(task709);

  auto tensor710 = vector<shared_ptr<Tensor>>{I1313, Gamma51_(), v2_};
  auto task710 = make_shared<Task710>(tensor710, pindex);
  task708->add_dep(task710);
  task710->add_dep(task108);
  residualq->add_task(task710);

  vector<IndexRange> I1322_index = {closed_, closed_, closed_, active_};
  auto I1322 = make_shared<Tensor>(I1322_index);
  auto tensor711 = vector<shared_ptr<Tensor>>{I162, t2, I1322};
  auto task711 = make_shared<Task711>(tensor711, pindex);
  task603->add_dep(task711);
  task711->add_dep(task108);
  residualq->add_task(task711);

  auto tensor712 = vector<shared_ptr<Tensor>>{I1322, Gamma32_(), v2_};
  auto task712 = make_shared<Task712>(tensor712, pindex);
  task711->add_dep(task712);
  task712->add_dep(task108);
  residualq->add_task(task712);

  vector<IndexRange> I1325_index = {closed_, closed_, closed_, active_};
  auto I1325 = make_shared<Tensor>(I1325_index);
  auto tensor713 = vector<shared_ptr<Tensor>>{I162, t2, I1325};
  auto task713 = make_shared<Task713>(tensor713, pindex);
  task603->add_dep(task713);
  task713->add_dep(task108);
  residualq->add_task(task713);

  auto tensor714 = vector<shared_ptr<Tensor>>{I1325, Gamma32_(), v2_};
  auto task714 = make_shared<Task714>(tensor714, pindex);
  task713->add_dep(task714);
  task714->add_dep(task108);
  residualq->add_task(task714);

  vector<IndexRange> I1328_index = {virt_, closed_, virt_, active_};
  auto I1328 = make_shared<Tensor>(I1328_index);
  auto tensor715 = vector<shared_ptr<Tensor>>{I162, t2, I1328};
  auto task715 = make_shared<Task715>(tensor715, pindex);
  task603->add_dep(task715);
  task715->add_dep(task108);
  residualq->add_task(task715);

  vector<IndexRange> I1329_index = {virt_, active_, closed_, virt_};
  auto I1329 = make_shared<Tensor>(I1329_index);
  auto tensor716 = vector<shared_ptr<Tensor>>{I1328, Gamma32_(), I1329};
  auto task716 = make_shared<Task716>(tensor716, pindex);
  task715->add_dep(task716);
  task716->add_dep(task108);
  residualq->add_task(task716);

  auto tensor717 = vector<shared_ptr<Tensor>>{I1329, v2_};
  auto task717 = make_shared<Task717>(tensor717, pindex);
  task716->add_dep(task717);
  task717->add_dep(task108);
  residualq->add_task(task717);

  vector<IndexRange> I1331_index = {virt_, closed_, virt_, active_};
  auto I1331 = make_shared<Tensor>(I1331_index);
  auto tensor718 = vector<shared_ptr<Tensor>>{I162, t2, I1331};
  auto task718 = make_shared<Task718>(tensor718, pindex);
  task603->add_dep(task718);
  task718->add_dep(task108);
  residualq->add_task(task718);

  vector<IndexRange> I1332_index = {virt_, active_, closed_, virt_};
  auto I1332 = make_shared<Tensor>(I1332_index);
  auto tensor719 = vector<shared_ptr<Tensor>>{I1331, Gamma32_(), I1332};
  auto task719 = make_shared<Task719>(tensor719, pindex);
  task718->add_dep(task719);
  task719->add_dep(task108);
  residualq->add_task(task719);

  auto tensor720 = vector<shared_ptr<Tensor>>{I1332, v2_};
  auto task720 = make_shared<Task720>(tensor720, pindex);
  task719->add_dep(task720);
  task720->add_dep(task108);
  residualq->add_task(task720);

  vector<IndexRange> I1334_index = {virt_, closed_, virt_, active_};
  auto I1334 = make_shared<Tensor>(I1334_index);
  auto tensor721 = vector<shared_ptr<Tensor>>{I162, t2, I1334};
  auto task721 = make_shared<Task721>(tensor721, pindex);
  task603->add_dep(task721);
  task721->add_dep(task108);
  residualq->add_task(task721);

  vector<IndexRange> I1335_index = {virt_, active_, closed_, virt_};
  auto I1335 = make_shared<Tensor>(I1335_index);
  auto tensor722 = vector<shared_ptr<Tensor>>{I1334, Gamma32_(), I1335};
  auto task722 = make_shared<Task722>(tensor722, pindex);
  task721->add_dep(task722);
  task722->add_dep(task108);
  residualq->add_task(task722);

  auto tensor723 = vector<shared_ptr<Tensor>>{I1335, v2_};
  auto task723 = make_shared<Task723>(tensor723, pindex);
  task722->add_dep(task723);
  task723->add_dep(task108);
  residualq->add_task(task723);

  vector<IndexRange> I1337_index = {virt_, closed_, virt_, active_};
  auto I1337 = make_shared<Tensor>(I1337_index);
  auto tensor724 = vector<shared_ptr<Tensor>>{I162, t2, I1337};
  auto task724 = make_shared<Task724>(tensor724, pindex);
  task603->add_dep(task724);
  task724->add_dep(task108);
  residualq->add_task(task724);

  vector<IndexRange> I1338_index = {virt_, active_, closed_, virt_};
  auto I1338 = make_shared<Tensor>(I1338_index);
  auto tensor725 = vector<shared_ptr<Tensor>>{I1337, Gamma32_(), I1338};
  auto task725 = make_shared<Task725>(tensor725, pindex);
  task724->add_dep(task725);
  task725->add_dep(task108);
  residualq->add_task(task725);

  auto tensor726 = vector<shared_ptr<Tensor>>{I1338, v2_};
  auto task726 = make_shared<Task726>(tensor726, pindex);
  task725->add_dep(task726);
  task726->add_dep(task108);
  residualq->add_task(task726);
}

#endif
