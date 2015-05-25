//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_gamma.cc
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
using namespace bagel::SMITH::RelMRCI;

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma0_() {
  vector<IndexRange> Gamma0_index = {active_, active_, active_, active_};
  auto Gamma0 = make_shared<Tensor>(Gamma0_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor0 = vector<shared_ptr<Tensor>>{Gamma0, rdm0_, rdm1_, rdm2_};
  auto task0 = make_shared<Task0>(tensor0, pindex);
  return make_shared<FutureTensor>(*Gamma0, task0);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma1_() {
  vector<IndexRange> Gamma1_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma1 = make_shared<Tensor>(Gamma1_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor1 = vector<shared_ptr<Tensor>>{Gamma1, rdm1_, rdm2_, rdm3_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  return make_shared<FutureTensor>(*Gamma1, task1);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma2_() {
  vector<IndexRange> Gamma2_index = {active_, active_, active_, active_};
  auto Gamma2 = make_shared<Tensor>(Gamma2_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor2 = vector<shared_ptr<Tensor>>{Gamma2, rdm0_, rdm1_, rdm2_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  return make_shared<FutureTensor>(*Gamma2, task2);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma84_() {
  vector<IndexRange> Gamma84_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma84 = make_shared<Tensor>(Gamma84_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor3 = vector<shared_ptr<Tensor>>{Gamma84, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  return make_shared<FutureTensor>(*Gamma84, task3);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma85_() {
  vector<IndexRange> Gamma85_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma85 = make_shared<Tensor>(Gamma85_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor4 = vector<shared_ptr<Tensor>>{Gamma85, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  return make_shared<FutureTensor>(*Gamma85, task4);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma86_() {
  vector<IndexRange> Gamma86_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma86 = make_shared<Tensor>(Gamma86_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor5 = vector<shared_ptr<Tensor>>{Gamma86, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  return make_shared<FutureTensor>(*Gamma86, task5);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma89_() {
  vector<IndexRange> Gamma89_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma89 = make_shared<Tensor>(Gamma89_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor6 = vector<shared_ptr<Tensor>>{Gamma89, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  return make_shared<FutureTensor>(*Gamma89, task6);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma90_() {
  vector<IndexRange> Gamma90_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma90 = make_shared<Tensor>(Gamma90_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor7 = vector<shared_ptr<Tensor>>{Gamma90, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task7 = make_shared<Task7>(tensor7, pindex);
  return make_shared<FutureTensor>(*Gamma90, task7);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma91_() {
  vector<IndexRange> Gamma91_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma91 = make_shared<Tensor>(Gamma91_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor8 = vector<shared_ptr<Tensor>>{Gamma91, rdm1_, rdm2_, rdm3_};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  return make_shared<FutureTensor>(*Gamma91, task8);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma92_() {
  vector<IndexRange> Gamma92_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma92 = make_shared<Tensor>(Gamma92_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor9 = vector<shared_ptr<Tensor>>{Gamma92, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  return make_shared<FutureTensor>(*Gamma92, task9);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma93_() {
  vector<IndexRange> Gamma93_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma93 = make_shared<Tensor>(Gamma93_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor10 = vector<shared_ptr<Tensor>>{Gamma93, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  return make_shared<FutureTensor>(*Gamma93, task10);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma98_() {
  vector<IndexRange> Gamma98_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma98 = make_shared<Tensor>(Gamma98_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor11 = vector<shared_ptr<Tensor>>{Gamma98, rdm1_, rdm2_, rdm3_};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  return make_shared<FutureTensor>(*Gamma98, task11);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma3_() {
  vector<IndexRange> Gamma3_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma3 = make_shared<Tensor>(Gamma3_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor12 = vector<shared_ptr<Tensor>>{Gamma3, rdm1_, rdm2_, rdm3_};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  return make_shared<FutureTensor>(*Gamma3, task12);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma4_() {
  vector<IndexRange> Gamma4_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma4 = make_shared<Tensor>(Gamma4_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor13 = vector<shared_ptr<Tensor>>{Gamma4, rdm1_, rdm2_, rdm3_};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  return make_shared<FutureTensor>(*Gamma4, task13);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma5_() {
  vector<IndexRange> Gamma5_index = {active_, active_, active_, active_};
  auto Gamma5 = make_shared<Tensor>(Gamma5_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor14 = vector<shared_ptr<Tensor>>{Gamma5, rdm1_, rdm2_};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  return make_shared<FutureTensor>(*Gamma5, task14);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma7_() {
  vector<IndexRange> Gamma7_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma7 = make_shared<Tensor>(Gamma7_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor15 = vector<shared_ptr<Tensor>>{Gamma7, rdm1_, rdm2_, rdm3_};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  return make_shared<FutureTensor>(*Gamma7, task15);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma101_() {
  vector<IndexRange> Gamma101_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma101 = make_shared<Tensor>(Gamma101_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor16 = vector<shared_ptr<Tensor>>{Gamma101, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  return make_shared<FutureTensor>(*Gamma101, task16);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma102_() {
  vector<IndexRange> Gamma102_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma102 = make_shared<Tensor>(Gamma102_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor17 = vector<shared_ptr<Tensor>>{Gamma102, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  return make_shared<FutureTensor>(*Gamma102, task17);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma104_() {
  vector<IndexRange> Gamma104_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma104 = make_shared<Tensor>(Gamma104_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor18 = vector<shared_ptr<Tensor>>{Gamma104, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  return make_shared<FutureTensor>(*Gamma104, task18);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma105_() {
  vector<IndexRange> Gamma105_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma105 = make_shared<Tensor>(Gamma105_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor19 = vector<shared_ptr<Tensor>>{Gamma105, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  return make_shared<FutureTensor>(*Gamma105, task19);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma106_() {
  vector<IndexRange> Gamma106_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma106 = make_shared<Tensor>(Gamma106_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor20 = vector<shared_ptr<Tensor>>{Gamma106, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  return make_shared<FutureTensor>(*Gamma106, task20);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma108_() {
  vector<IndexRange> Gamma108_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma108 = make_shared<Tensor>(Gamma108_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor21 = vector<shared_ptr<Tensor>>{Gamma108, rdm1_, rdm2_, rdm3_};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  return make_shared<FutureTensor>(*Gamma108, task21);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma111_() {
  vector<IndexRange> Gamma111_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma111 = make_shared<Tensor>(Gamma111_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor22 = vector<shared_ptr<Tensor>>{Gamma111, rdm1_, rdm2_, rdm3_};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  return make_shared<FutureTensor>(*Gamma111, task22);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma113_() {
  vector<IndexRange> Gamma113_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma113 = make_shared<Tensor>(Gamma113_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor23 = vector<shared_ptr<Tensor>>{Gamma113, rdm1_, rdm2_, rdm3_};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  return make_shared<FutureTensor>(*Gamma113, task23);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma118_() {
  vector<IndexRange> Gamma118_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma118 = make_shared<Tensor>(Gamma118_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor24 = vector<shared_ptr<Tensor>>{Gamma118, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  return make_shared<FutureTensor>(*Gamma118, task24);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma119_() {
  vector<IndexRange> Gamma119_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma119 = make_shared<Tensor>(Gamma119_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor25 = vector<shared_ptr<Tensor>>{Gamma119, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  return make_shared<FutureTensor>(*Gamma119, task25);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma123_() {
  vector<IndexRange> Gamma123_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma123 = make_shared<Tensor>(Gamma123_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor26 = vector<shared_ptr<Tensor>>{Gamma123, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  return make_shared<FutureTensor>(*Gamma123, task26);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma126_() {
  vector<IndexRange> Gamma126_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma126 = make_shared<Tensor>(Gamma126_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor27 = vector<shared_ptr<Tensor>>{Gamma126, rdm2_, rdm3_, rdm4_};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  return make_shared<FutureTensor>(*Gamma126, task27);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma558_() {
  vector<IndexRange> Gamma558_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma558 = make_shared<Tensor>(Gamma558_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor28 = vector<shared_ptr<Tensor>>{Gamma558, rdm1_, rdm2_, rdm3_, rdm4_, h1_};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  return make_shared<FutureTensor>(*Gamma558, task28);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma559_() {
  vector<IndexRange> Gamma559_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma559 = make_shared<Tensor>(Gamma559_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor29 = vector<shared_ptr<Tensor>>{Gamma559, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  return make_shared<FutureTensor>(*Gamma559, task29);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma10_() {
  vector<IndexRange> Gamma10_index = {active_, active_, active_, active_};
  auto Gamma10 = make_shared<Tensor>(Gamma10_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor30 = vector<shared_ptr<Tensor>>{Gamma10, rdm1_, rdm2_};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  return make_shared<FutureTensor>(*Gamma10, task30);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma12_() {
  vector<IndexRange> Gamma12_index = {active_, active_};
  auto Gamma12 = make_shared<Tensor>(Gamma12_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor31 = vector<shared_ptr<Tensor>>{Gamma12, rdm0_, rdm1_};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  return make_shared<FutureTensor>(*Gamma12, task31);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma18_() {
  vector<IndexRange> Gamma18_index = {active_, active_, active_, active_};
  auto Gamma18 = make_shared<Tensor>(Gamma18_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor32 = vector<shared_ptr<Tensor>>{Gamma18, rdm1_, rdm2_};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  return make_shared<FutureTensor>(*Gamma18, task32);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma201_() {
  vector<IndexRange> Gamma201_index = {active_, active_, active_, active_};
  auto Gamma201 = make_shared<Tensor>(Gamma201_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor33 = vector<shared_ptr<Tensor>>{Gamma201, rdm0_, rdm1_, rdm2_};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  return make_shared<FutureTensor>(*Gamma201, task33);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma130_() {
  vector<IndexRange> Gamma130_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma130 = make_shared<Tensor>(Gamma130_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor34 = vector<shared_ptr<Tensor>>{Gamma130, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  return make_shared<FutureTensor>(*Gamma130, task34);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma136_() {
  vector<IndexRange> Gamma136_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma136 = make_shared<Tensor>(Gamma136_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor35 = vector<shared_ptr<Tensor>>{Gamma136, rdm1_, rdm2_, rdm3_};
  auto task35 = make_shared<Task35>(tensor35, pindex);
  return make_shared<FutureTensor>(*Gamma136, task35);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma141_() {
  vector<IndexRange> Gamma141_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma141 = make_shared<Tensor>(Gamma141_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor36 = vector<shared_ptr<Tensor>>{Gamma141, rdm1_, rdm2_, rdm3_};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  return make_shared<FutureTensor>(*Gamma141, task36);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma159_() {
  vector<IndexRange> Gamma159_index = {active_, active_, active_, active_};
  auto Gamma159 = make_shared<Tensor>(Gamma159_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor37 = vector<shared_ptr<Tensor>>{Gamma159, rdm0_, rdm1_, rdm2_};
  auto task37 = make_shared<Task37>(tensor37, pindex);
  return make_shared<FutureTensor>(*Gamma159, task37);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma180_() {
  vector<IndexRange> Gamma180_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma180 = make_shared<Tensor>(Gamma180_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor38 = vector<shared_ptr<Tensor>>{Gamma180, rdm1_, rdm2_, rdm3_};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  return make_shared<FutureTensor>(*Gamma180, task38);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma182_() {
  vector<IndexRange> Gamma182_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma182 = make_shared<Tensor>(Gamma182_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor39 = vector<shared_ptr<Tensor>>{Gamma182, rdm1_, rdm2_, rdm3_};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  return make_shared<FutureTensor>(*Gamma182, task39);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma183_() {
  vector<IndexRange> Gamma183_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma183 = make_shared<Tensor>(Gamma183_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor40 = vector<shared_ptr<Tensor>>{Gamma183, rdm1_, rdm2_, rdm3_};
  auto task40 = make_shared<Task40>(tensor40, pindex);
  return make_shared<FutureTensor>(*Gamma183, task40);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma200_() {
  vector<IndexRange> Gamma200_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma200 = make_shared<Tensor>(Gamma200_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor41 = vector<shared_ptr<Tensor>>{Gamma200, rdm2_, rdm3_};
  auto task41 = make_shared<Task41>(tensor41, pindex);
  return make_shared<FutureTensor>(*Gamma200, task41);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma560_() {
  vector<IndexRange> Gamma560_index = {active_, active_};
  auto Gamma560 = make_shared<Tensor>(Gamma560_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor42 = vector<shared_ptr<Tensor>>{Gamma560, rdm0_, rdm1_, rdm2_, h1_};
  auto task42 = make_shared<Task42>(tensor42, pindex);
  return make_shared<FutureTensor>(*Gamma560, task42);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma562_() {
  vector<IndexRange> Gamma562_index = {active_, active_};
  auto Gamma562 = make_shared<Tensor>(Gamma562_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor43 = vector<shared_ptr<Tensor>>{Gamma562, rdm0_, rdm1_, rdm2_, rdm3_, v2_};
  auto task43 = make_shared<Task43>(tensor43, pindex);
  return make_shared<FutureTensor>(*Gamma562, task43);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma24_() {
  vector<IndexRange> Gamma24_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma24 = make_shared<Tensor>(Gamma24_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor44 = vector<shared_ptr<Tensor>>{Gamma24, rdm1_, rdm2_, rdm3_};
  auto task44 = make_shared<Task44>(tensor44, pindex);
  return make_shared<FutureTensor>(*Gamma24, task44);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma25_() {
  vector<IndexRange> Gamma25_index = {active_, active_, active_, active_};
  auto Gamma25 = make_shared<Tensor>(Gamma25_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor45 = vector<shared_ptr<Tensor>>{Gamma25, rdm1_, rdm2_};
  auto task45 = make_shared<Task45>(tensor45, pindex);
  return make_shared<FutureTensor>(*Gamma25, task45);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma27_() {
  vector<IndexRange> Gamma27_index = {active_, active_, active_, active_};
  auto Gamma27 = make_shared<Tensor>(Gamma27_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor46 = vector<shared_ptr<Tensor>>{Gamma27, rdm1_, rdm2_, rdm3_, h1_};
  auto task46 = make_shared<Task46>(tensor46, pindex);
  return make_shared<FutureTensor>(*Gamma27, task46);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma28_() {
  vector<IndexRange> Gamma28_index = {active_, active_, active_, active_};
  auto Gamma28 = make_shared<Tensor>(Gamma28_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor47 = vector<shared_ptr<Tensor>>{Gamma28, rdm1_, rdm2_};
  auto task47 = make_shared<Task47>(tensor47, pindex);
  return make_shared<FutureTensor>(*Gamma28, task47);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma30_() {
  vector<IndexRange> Gamma30_index = {active_, active_, active_, active_};
  auto Gamma30 = make_shared<Tensor>(Gamma30_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor48 = vector<shared_ptr<Tensor>>{Gamma30, rdm1_, rdm2_, rdm3_, h1_};
  auto task48 = make_shared<Task48>(tensor48, pindex);
  return make_shared<FutureTensor>(*Gamma30, task48);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma31_() {
  vector<IndexRange> Gamma31_index = {active_, active_, active_, active_};
  auto Gamma31 = make_shared<Tensor>(Gamma31_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor49 = vector<shared_ptr<Tensor>>{Gamma31, rdm1_, rdm2_};
  auto task49 = make_shared<Task49>(tensor49, pindex);
  return make_shared<FutureTensor>(*Gamma31, task49);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma33_() {
  vector<IndexRange> Gamma33_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma33 = make_shared<Tensor>(Gamma33_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor50 = vector<shared_ptr<Tensor>>{Gamma33, rdm2_, rdm3_};
  auto task50 = make_shared<Task50>(tensor50, pindex);
  return make_shared<FutureTensor>(*Gamma33, task50);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma34_() {
  vector<IndexRange> Gamma34_index = {active_, active_};
  auto Gamma34 = make_shared<Tensor>(Gamma34_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor51 = vector<shared_ptr<Tensor>>{Gamma34, rdm1_};
  auto task51 = make_shared<Task51>(tensor51, pindex);
  return make_shared<FutureTensor>(*Gamma34, task51);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma219_() {
  vector<IndexRange> Gamma219_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma219 = make_shared<Tensor>(Gamma219_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor52 = vector<shared_ptr<Tensor>>{Gamma219, rdm1_, rdm2_, rdm3_};
  auto task52 = make_shared<Task52>(tensor52, pindex);
  return make_shared<FutureTensor>(*Gamma219, task52);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma220_() {
  vector<IndexRange> Gamma220_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma220 = make_shared<Tensor>(Gamma220_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor53 = vector<shared_ptr<Tensor>>{Gamma220, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task53 = make_shared<Task53>(tensor53, pindex);
  return make_shared<FutureTensor>(*Gamma220, task53);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma221_() {
  vector<IndexRange> Gamma221_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma221 = make_shared<Tensor>(Gamma221_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor54 = vector<shared_ptr<Tensor>>{Gamma221, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task54 = make_shared<Task54>(tensor54, pindex);
  return make_shared<FutureTensor>(*Gamma221, task54);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma224_() {
  vector<IndexRange> Gamma224_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma224 = make_shared<Tensor>(Gamma224_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor55 = vector<shared_ptr<Tensor>>{Gamma224, rdm1_, rdm2_, rdm3_};
  auto task55 = make_shared<Task55>(tensor55, pindex);
  return make_shared<FutureTensor>(*Gamma224, task55);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma226_() {
  vector<IndexRange> Gamma226_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma226 = make_shared<Tensor>(Gamma226_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor56 = vector<shared_ptr<Tensor>>{Gamma226, rdm1_, rdm2_, rdm3_};
  auto task56 = make_shared<Task56>(tensor56, pindex);
  return make_shared<FutureTensor>(*Gamma226, task56);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma225_() {
  vector<IndexRange> Gamma225_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma225 = make_shared<Tensor>(Gamma225_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor57 = vector<shared_ptr<Tensor>>{Gamma225, rdm1_, rdm2_, rdm3_};
  auto task57 = make_shared<Task57>(tensor57, pindex);
  return make_shared<FutureTensor>(*Gamma225, task57);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma234_() {
  vector<IndexRange> Gamma234_index = {active_, active_, active_, active_};
  auto Gamma234 = make_shared<Tensor>(Gamma234_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor58 = vector<shared_ptr<Tensor>>{Gamma234, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task58 = make_shared<Task58>(tensor58, pindex);
  return make_shared<FutureTensor>(*Gamma234, task58);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma235_() {
  vector<IndexRange> Gamma235_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma235 = make_shared<Tensor>(Gamma235_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor59 = vector<shared_ptr<Tensor>>{Gamma235, rdm1_, rdm2_, rdm3_};
  auto task59 = make_shared<Task59>(tensor59, pindex);
  return make_shared<FutureTensor>(*Gamma235, task59);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma237_() {
  vector<IndexRange> Gamma237_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma237 = make_shared<Tensor>(Gamma237_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor60 = vector<shared_ptr<Tensor>>{Gamma237, rdm1_, rdm2_, rdm3_};
  auto task60 = make_shared<Task60>(tensor60, pindex);
  return make_shared<FutureTensor>(*Gamma237, task60);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma239_() {
  vector<IndexRange> Gamma239_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma239 = make_shared<Tensor>(Gamma239_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor61 = vector<shared_ptr<Tensor>>{Gamma239, rdm1_, rdm2_, rdm3_};
  auto task61 = make_shared<Task61>(tensor61, pindex);
  return make_shared<FutureTensor>(*Gamma239, task61);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma238_() {
  vector<IndexRange> Gamma238_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma238 = make_shared<Tensor>(Gamma238_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor62 = vector<shared_ptr<Tensor>>{Gamma238, rdm1_, rdm2_, rdm3_};
  auto task62 = make_shared<Task62>(tensor62, pindex);
  return make_shared<FutureTensor>(*Gamma238, task62);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma240_() {
  vector<IndexRange> Gamma240_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma240 = make_shared<Tensor>(Gamma240_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor63 = vector<shared_ptr<Tensor>>{Gamma240, rdm1_, rdm2_, rdm3_};
  auto task63 = make_shared<Task63>(tensor63, pindex);
  return make_shared<FutureTensor>(*Gamma240, task63);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma245_() {
  vector<IndexRange> Gamma245_index = {active_, active_, active_, active_};
  auto Gamma245 = make_shared<Tensor>(Gamma245_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor64 = vector<shared_ptr<Tensor>>{Gamma245, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task64 = make_shared<Task64>(tensor64, pindex);
  return make_shared<FutureTensor>(*Gamma245, task64);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma246_() {
  vector<IndexRange> Gamma246_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma246 = make_shared<Tensor>(Gamma246_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor65 = vector<shared_ptr<Tensor>>{Gamma246, rdm1_, rdm2_, rdm3_};
  auto task65 = make_shared<Task65>(tensor65, pindex);
  return make_shared<FutureTensor>(*Gamma246, task65);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma250_() {
  vector<IndexRange> Gamma250_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma250 = make_shared<Tensor>(Gamma250_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor66 = vector<shared_ptr<Tensor>>{Gamma250, rdm1_, rdm2_, rdm3_};
  auto task66 = make_shared<Task66>(tensor66, pindex);
  return make_shared<FutureTensor>(*Gamma250, task66);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma256_() {
  vector<IndexRange> Gamma256_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma256 = make_shared<Tensor>(Gamma256_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor67 = vector<shared_ptr<Tensor>>{Gamma256, rdm2_, rdm3_, rdm4_};
  auto task67 = make_shared<Task67>(tensor67, pindex);
  return make_shared<FutureTensor>(*Gamma256, task67);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma257_() {
  vector<IndexRange> Gamma257_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma257 = make_shared<Tensor>(Gamma257_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor68 = vector<shared_ptr<Tensor>>{Gamma257, rdm2_, rdm3_, rdm4_};
  auto task68 = make_shared<Task68>(tensor68, pindex);
  return make_shared<FutureTensor>(*Gamma257, task68);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma258_() {
  vector<IndexRange> Gamma258_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma258 = make_shared<Tensor>(Gamma258_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor69 = vector<shared_ptr<Tensor>>{Gamma258, rdm2_, rdm3_};
  auto task69 = make_shared<Task69>(tensor69, pindex);
  return make_shared<FutureTensor>(*Gamma258, task69);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma282_() {
  vector<IndexRange> Gamma282_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma282 = make_shared<Tensor>(Gamma282_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor70 = vector<shared_ptr<Tensor>>{Gamma282, rdm2_, rdm3_};
  auto task70 = make_shared<Task70>(tensor70, pindex);
  return make_shared<FutureTensor>(*Gamma282, task70);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma284_() {
  vector<IndexRange> Gamma284_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma284 = make_shared<Tensor>(Gamma284_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor71 = vector<shared_ptr<Tensor>>{Gamma284, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task71 = make_shared<Task71>(tensor71, pindex);
  return make_shared<FutureTensor>(*Gamma284, task71);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma303_() {
  vector<IndexRange> Gamma303_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma303 = make_shared<Tensor>(Gamma303_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor72 = vector<shared_ptr<Tensor>>{Gamma303, rdm1_, rdm2_, rdm3_};
  auto task72 = make_shared<Task72>(tensor72, pindex);
  return make_shared<FutureTensor>(*Gamma303, task72);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma320_() {
  vector<IndexRange> Gamma320_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma320 = make_shared<Tensor>(Gamma320_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor73 = vector<shared_ptr<Tensor>>{Gamma320, rdm2_, rdm3_, rdm4_};
  auto task73 = make_shared<Task73>(tensor73, pindex);
  return make_shared<FutureTensor>(*Gamma320, task73);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma321_() {
  vector<IndexRange> Gamma321_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma321 = make_shared<Tensor>(Gamma321_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor74 = vector<shared_ptr<Tensor>>{Gamma321, rdm2_, rdm3_, rdm4_};
  auto task74 = make_shared<Task74>(tensor74, pindex);
  return make_shared<FutureTensor>(*Gamma321, task74);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma346_() {
  vector<IndexRange> Gamma346_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma346 = make_shared<Tensor>(Gamma346_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor75 = vector<shared_ptr<Tensor>>{Gamma346, rdm2_, rdm3_};
  auto task75 = make_shared<Task75>(tensor75, pindex);
  return make_shared<FutureTensor>(*Gamma346, task75);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma52_() {
  vector<IndexRange> Gamma52_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma52 = make_shared<Tensor>(Gamma52_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor76 = vector<shared_ptr<Tensor>>{Gamma52, rdm2_, rdm3_};
  auto task76 = make_shared<Task76>(tensor76, pindex);
  return make_shared<FutureTensor>(*Gamma52, task76);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma53_() {
  vector<IndexRange> Gamma53_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma53 = make_shared<Tensor>(Gamma53_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor77 = vector<shared_ptr<Tensor>>{Gamma53, rdm2_, rdm3_};
  auto task77 = make_shared<Task77>(tensor77, pindex);
  return make_shared<FutureTensor>(*Gamma53, task77);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma54_() {
  vector<IndexRange> Gamma54_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma54 = make_shared<Tensor>(Gamma54_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor78 = vector<shared_ptr<Tensor>>{Gamma54, rdm2_, rdm3_};
  auto task78 = make_shared<Task78>(tensor78, pindex);
  return make_shared<FutureTensor>(*Gamma54, task78);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma55_() {
  vector<IndexRange> Gamma55_index = {active_, active_, active_, active_};
  auto Gamma55 = make_shared<Tensor>(Gamma55_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor79 = vector<shared_ptr<Tensor>>{Gamma55, rdm2_};
  auto task79 = make_shared<Task79>(tensor79, pindex);
  return make_shared<FutureTensor>(*Gamma55, task79);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma347_() {
  vector<IndexRange> Gamma347_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma347 = make_shared<Tensor>(Gamma347_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor80 = vector<shared_ptr<Tensor>>{Gamma347, rdm2_, rdm3_, rdm4_};
  auto task80 = make_shared<Task80>(tensor80, pindex);
  return make_shared<FutureTensor>(*Gamma347, task80);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma348_() {
  vector<IndexRange> Gamma348_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma348 = make_shared<Tensor>(Gamma348_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor81 = vector<shared_ptr<Tensor>>{Gamma348, rdm2_, rdm3_};
  auto task81 = make_shared<Task81>(tensor81, pindex);
  return make_shared<FutureTensor>(*Gamma348, task81);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma349_() {
  vector<IndexRange> Gamma349_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma349 = make_shared<Tensor>(Gamma349_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor82 = vector<shared_ptr<Tensor>>{Gamma349, rdm2_, rdm3_, rdm4_};
  auto task82 = make_shared<Task82>(tensor82, pindex);
  return make_shared<FutureTensor>(*Gamma349, task82);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma350_() {
  vector<IndexRange> Gamma350_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma350 = make_shared<Tensor>(Gamma350_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor83 = vector<shared_ptr<Tensor>>{Gamma350, rdm2_, rdm3_, rdm4_};
  auto task83 = make_shared<Task83>(tensor83, pindex);
  return make_shared<FutureTensor>(*Gamma350, task83);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma353_() {
  vector<IndexRange> Gamma353_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma353 = make_shared<Tensor>(Gamma353_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor84 = vector<shared_ptr<Tensor>>{Gamma353, rdm2_, rdm3_, rdm4_};
  auto task84 = make_shared<Task84>(tensor84, pindex);
  return make_shared<FutureTensor>(*Gamma353, task84);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma354_() {
  vector<IndexRange> Gamma354_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma354 = make_shared<Tensor>(Gamma354_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor85 = vector<shared_ptr<Tensor>>{Gamma354, rdm2_, rdm3_, rdm4_};
  auto task85 = make_shared<Task85>(tensor85, pindex);
  return make_shared<FutureTensor>(*Gamma354, task85);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma357_() {
  vector<IndexRange> Gamma357_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma357 = make_shared<Tensor>(Gamma357_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor86 = vector<shared_ptr<Tensor>>{Gamma357, rdm2_, rdm3_, rdm4_};
  auto task86 = make_shared<Task86>(tensor86, pindex);
  return make_shared<FutureTensor>(*Gamma357, task86);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma358_() {
  vector<IndexRange> Gamma358_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma358 = make_shared<Tensor>(Gamma358_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor87 = vector<shared_ptr<Tensor>>{Gamma358, rdm2_, rdm3_, rdm4_};
  auto task87 = make_shared<Task87>(tensor87, pindex);
  return make_shared<FutureTensor>(*Gamma358, task87);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma359_() {
  vector<IndexRange> Gamma359_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma359 = make_shared<Tensor>(Gamma359_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor88 = vector<shared_ptr<Tensor>>{Gamma359, rdm2_, rdm3_, rdm4_};
  auto task88 = make_shared<Task88>(tensor88, pindex);
  return make_shared<FutureTensor>(*Gamma359, task88);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma367_() {
  vector<IndexRange> Gamma367_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma367 = make_shared<Tensor>(Gamma367_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor89 = vector<shared_ptr<Tensor>>{Gamma367, rdm2_, rdm3_};
  auto task89 = make_shared<Task89>(tensor89, pindex);
  return make_shared<FutureTensor>(*Gamma367, task89);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma374_() {
  vector<IndexRange> Gamma374_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma374 = make_shared<Tensor>(Gamma374_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor90 = vector<shared_ptr<Tensor>>{Gamma374, rdm3_, rdm4_};
  auto task90 = make_shared<Task90>(tensor90, pindex);
  return make_shared<FutureTensor>(*Gamma374, task90);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma564_() {
  vector<IndexRange> Gamma564_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma564 = make_shared<Tensor>(Gamma564_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor91 = vector<shared_ptr<Tensor>>{Gamma564, rdm2_, rdm3_, rdm4_, h1_};
  auto task91 = make_shared<Task91>(tensor91, pindex);
  return make_shared<FutureTensor>(*Gamma564, task91);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma565_() {
  vector<IndexRange> Gamma565_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma565 = make_shared<Tensor>(Gamma565_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor92 = vector<shared_ptr<Tensor>>{Gamma565, rdm2_, rdm3_, rdm4_, v2_};
  auto task92 = make_shared<Task92>(tensor92, pindex);
  return make_shared<FutureTensor>(*Gamma565, task92);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma479_() {
  vector<IndexRange> Gamma479_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma479 = make_shared<Tensor>(Gamma479_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor93 = vector<shared_ptr<Tensor>>{Gamma479, rdm2_, rdm3_};
  auto task93 = make_shared<Task93>(tensor93, pindex);
  return make_shared<FutureTensor>(*Gamma479, task93);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma511_() {
  vector<IndexRange> Gamma511_index = {active_, active_, active_, active_};
  auto Gamma511 = make_shared<Tensor>(Gamma511_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor94 = vector<shared_ptr<Tensor>>{Gamma511, rdm2_};
  auto task94 = make_shared<Task94>(tensor94, pindex);
  return make_shared<FutureTensor>(*Gamma511, task94);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma534_() {
  vector<IndexRange> Gamma534_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma534 = make_shared<Tensor>(Gamma534_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor95 = vector<shared_ptr<Tensor>>{Gamma534, rdm3_};
  auto task95 = make_shared<Task95>(tensor95, pindex);
  return make_shared<FutureTensor>(*Gamma534, task95);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma570_() {
  vector<IndexRange> Gamma570_index = {active_, active_};
  auto Gamma570 = make_shared<Tensor>(Gamma570_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor96 = vector<shared_ptr<Tensor>>{Gamma570, rdm1_, rdm2_, h1_};
  auto task96 = make_shared<Task96>(tensor96, pindex);
  return make_shared<FutureTensor>(*Gamma570, task96);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma572_() {
  vector<IndexRange> Gamma572_index = {active_, active_};
  auto Gamma572 = make_shared<Tensor>(Gamma572_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor97 = vector<shared_ptr<Tensor>>{Gamma572, rdm1_, rdm2_, rdm3_, v2_};
  auto task97 = make_shared<Task97>(tensor97, pindex);
  return make_shared<FutureTensor>(*Gamma572, task97);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma539_() {
  vector<IndexRange> Gamma539_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma539 = make_shared<Tensor>(Gamma539_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor98 = vector<shared_ptr<Tensor>>{Gamma539, rdm2_, rdm3_};
  auto task98 = make_shared<Task98>(tensor98, pindex);
  return make_shared<FutureTensor>(*Gamma539, task98);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma540_() {
  vector<IndexRange> Gamma540_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma540 = make_shared<Tensor>(Gamma540_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor99 = vector<shared_ptr<Tensor>>{Gamma540, rdm2_, rdm3_};
  auto task99 = make_shared<Task99>(tensor99, pindex);
  return make_shared<FutureTensor>(*Gamma540, task99);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma541_() {
  vector<IndexRange> Gamma541_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma541 = make_shared<Tensor>(Gamma541_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor100 = vector<shared_ptr<Tensor>>{Gamma541, rdm3_, rdm4_};
  auto task100 = make_shared<Task100>(tensor100, pindex);
  return make_shared<FutureTensor>(*Gamma541, task100);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma553_() {
  vector<IndexRange> Gamma553_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma553 = make_shared<Tensor>(Gamma553_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor101 = vector<shared_ptr<Tensor>>{Gamma553, rdm3_};
  auto task101 = make_shared<Task101>(tensor101, pindex);
  return make_shared<FutureTensor>(*Gamma553, task101);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma556_() {
  vector<IndexRange> Gamma556_index = {active_, active_, active_, active_};
  auto Gamma556 = make_shared<Tensor>(Gamma556_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor102 = vector<shared_ptr<Tensor>>{Gamma556, rdm0_, rdm1_, rdm2_, rdm3_, h1_};
  auto task102 = make_shared<Task102>(tensor102, pindex);
  return make_shared<FutureTensor>(*Gamma556, task102);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma557_() {
  vector<IndexRange> Gamma557_index = {active_, active_, active_, active_};
  auto Gamma557 = make_shared<Tensor>(Gamma557_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor103 = vector<shared_ptr<Tensor>>{Gamma557, rdm0_, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task103 = make_shared<Task103>(tensor103, pindex);
  return make_shared<FutureTensor>(*Gamma557, task103);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma566_() {
  vector<IndexRange> Gamma566_index;
  auto Gamma566 = make_shared<Tensor>(Gamma566_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor104 = vector<shared_ptr<Tensor>>{Gamma566, rdm1_, h1_};
  auto task104 = make_shared<Task104>(tensor104, pindex);
  return make_shared<FutureTensor>(*Gamma566, task104);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma568_() {
  vector<IndexRange> Gamma568_index;
  auto Gamma568 = make_shared<Tensor>(Gamma568_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor105 = vector<shared_ptr<Tensor>>{Gamma568, rdm1_, rdm2_, v2_};
  auto task105 = make_shared<Task105>(tensor105, pindex);
  return make_shared<FutureTensor>(*Gamma568, task105);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma574_() {
  vector<IndexRange> Gamma574_index = {active_, active_, active_, active_};
  auto Gamma574 = make_shared<Tensor>(Gamma574_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor106 = vector<shared_ptr<Tensor>>{Gamma574, rdm2_, rdm3_, h1_};
  auto task106 = make_shared<Task106>(tensor106, pindex);
  return make_shared<FutureTensor>(*Gamma574, task106);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma575_() {
  vector<IndexRange> Gamma575_index = {active_, active_, active_, active_};
  auto Gamma575 = make_shared<Tensor>(Gamma575_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor107 = vector<shared_ptr<Tensor>>{Gamma575, rdm2_, rdm3_, rdm4_, v2_};
  auto task107 = make_shared<Task107>(tensor107, pindex);
  return make_shared<FutureTensor>(*Gamma575, task107);
}

#endif
