//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_gamma.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software; you can redistribute it and/or modify
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
#include <src/smith/relmrci/RelMRCI_tasks1.h>
#include <src/smith/relmrci/RelMRCI_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using bagel::SMITH::RelMRCI::FutureTensor;

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

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma58_() {
  vector<IndexRange> Gamma58_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma58 = make_shared<Tensor>(Gamma58_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor3 = vector<shared_ptr<Tensor>>{Gamma58, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  return make_shared<FutureTensor>(*Gamma58, task3);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma59_() {
  vector<IndexRange> Gamma59_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma59 = make_shared<Tensor>(Gamma59_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor4 = vector<shared_ptr<Tensor>>{Gamma59, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  return make_shared<FutureTensor>(*Gamma59, task4);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma60_() {
  vector<IndexRange> Gamma60_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma60 = make_shared<Tensor>(Gamma60_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor5 = vector<shared_ptr<Tensor>>{Gamma60, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  return make_shared<FutureTensor>(*Gamma60, task5);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma63_() {
  vector<IndexRange> Gamma63_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma63 = make_shared<Tensor>(Gamma63_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor6 = vector<shared_ptr<Tensor>>{Gamma63, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  return make_shared<FutureTensor>(*Gamma63, task6);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma64_() {
  vector<IndexRange> Gamma64_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma64 = make_shared<Tensor>(Gamma64_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor7 = vector<shared_ptr<Tensor>>{Gamma64, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task7 = make_shared<Task7>(tensor7, pindex);
  return make_shared<FutureTensor>(*Gamma64, task7);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma65_() {
  vector<IndexRange> Gamma65_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma65 = make_shared<Tensor>(Gamma65_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor8 = vector<shared_ptr<Tensor>>{Gamma65, rdm1_, rdm2_, rdm3_};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  return make_shared<FutureTensor>(*Gamma65, task8);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma66_() {
  vector<IndexRange> Gamma66_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma66 = make_shared<Tensor>(Gamma66_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor9 = vector<shared_ptr<Tensor>>{Gamma66, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  return make_shared<FutureTensor>(*Gamma66, task9);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma67_() {
  vector<IndexRange> Gamma67_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma67 = make_shared<Tensor>(Gamma67_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor10 = vector<shared_ptr<Tensor>>{Gamma67, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  return make_shared<FutureTensor>(*Gamma67, task10);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma3_() {
  vector<IndexRange> Gamma3_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma3 = make_shared<Tensor>(Gamma3_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor11 = vector<shared_ptr<Tensor>>{Gamma3, rdm1_, rdm2_, rdm3_};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  return make_shared<FutureTensor>(*Gamma3, task11);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma4_() {
  vector<IndexRange> Gamma4_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma4 = make_shared<Tensor>(Gamma4_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor12 = vector<shared_ptr<Tensor>>{Gamma4, rdm1_, rdm2_, rdm3_};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  return make_shared<FutureTensor>(*Gamma4, task12);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma5_() {
  vector<IndexRange> Gamma5_index = {active_, active_, active_, active_};
  auto Gamma5 = make_shared<Tensor>(Gamma5_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor13 = vector<shared_ptr<Tensor>>{Gamma5, rdm1_, rdm2_};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  return make_shared<FutureTensor>(*Gamma5, task13);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma74_() {
  vector<IndexRange> Gamma74_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma74 = make_shared<Tensor>(Gamma74_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor14 = vector<shared_ptr<Tensor>>{Gamma74, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  return make_shared<FutureTensor>(*Gamma74, task14);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma75_() {
  vector<IndexRange> Gamma75_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma75 = make_shared<Tensor>(Gamma75_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor15 = vector<shared_ptr<Tensor>>{Gamma75, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  return make_shared<FutureTensor>(*Gamma75, task15);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma77_() {
  vector<IndexRange> Gamma77_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma77 = make_shared<Tensor>(Gamma77_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor16 = vector<shared_ptr<Tensor>>{Gamma77, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  return make_shared<FutureTensor>(*Gamma77, task16);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma78_() {
  vector<IndexRange> Gamma78_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma78 = make_shared<Tensor>(Gamma78_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor17 = vector<shared_ptr<Tensor>>{Gamma78, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  return make_shared<FutureTensor>(*Gamma78, task17);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma79_() {
  vector<IndexRange> Gamma79_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma79 = make_shared<Tensor>(Gamma79_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor18 = vector<shared_ptr<Tensor>>{Gamma79, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  return make_shared<FutureTensor>(*Gamma79, task18);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma81_() {
  vector<IndexRange> Gamma81_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma81 = make_shared<Tensor>(Gamma81_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor19 = vector<shared_ptr<Tensor>>{Gamma81, rdm1_, rdm2_, rdm3_};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  return make_shared<FutureTensor>(*Gamma81, task19);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma84_() {
  vector<IndexRange> Gamma84_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma84 = make_shared<Tensor>(Gamma84_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor20 = vector<shared_ptr<Tensor>>{Gamma84, rdm1_, rdm2_, rdm3_};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  return make_shared<FutureTensor>(*Gamma84, task20);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma86_() {
  vector<IndexRange> Gamma86_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma86 = make_shared<Tensor>(Gamma86_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor21 = vector<shared_ptr<Tensor>>{Gamma86, rdm1_, rdm2_, rdm3_};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  return make_shared<FutureTensor>(*Gamma86, task21);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma92_() {
  vector<IndexRange> Gamma92_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma92 = make_shared<Tensor>(Gamma92_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor22 = vector<shared_ptr<Tensor>>{Gamma92, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  return make_shared<FutureTensor>(*Gamma92, task22);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma95_() {
  vector<IndexRange> Gamma95_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma95 = make_shared<Tensor>(Gamma95_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor23 = vector<shared_ptr<Tensor>>{Gamma95, rdm2_, rdm3_, rdm4_};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  return make_shared<FutureTensor>(*Gamma95, task23);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma98_() {
  vector<IndexRange> Gamma98_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma98 = make_shared<Tensor>(Gamma98_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor24 = vector<shared_ptr<Tensor>>{Gamma98, rdm1_, rdm2_, rdm3_};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  return make_shared<FutureTensor>(*Gamma98, task24);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma414_() {
  vector<IndexRange> Gamma414_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma414 = make_shared<Tensor>(Gamma414_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor25 = vector<shared_ptr<Tensor>>{Gamma414, rdm1_, rdm2_, rdm3_, rdm4_, h1_};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  return make_shared<FutureTensor>(*Gamma414, task25);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma415_() {
  vector<IndexRange> Gamma415_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma415 = make_shared<Tensor>(Gamma415_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor26 = vector<shared_ptr<Tensor>>{Gamma415, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  return make_shared<FutureTensor>(*Gamma415, task26);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma9_() {
  vector<IndexRange> Gamma9_index = {active_, active_, active_, active_};
  auto Gamma9 = make_shared<Tensor>(Gamma9_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor27 = vector<shared_ptr<Tensor>>{Gamma9, rdm1_, rdm2_};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  return make_shared<FutureTensor>(*Gamma9, task27);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma11_() {
  vector<IndexRange> Gamma11_index = {active_, active_};
  auto Gamma11 = make_shared<Tensor>(Gamma11_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor28 = vector<shared_ptr<Tensor>>{Gamma11, rdm0_, rdm1_};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  return make_shared<FutureTensor>(*Gamma11, task28);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma160_() {
  vector<IndexRange> Gamma160_index = {active_, active_, active_, active_};
  auto Gamma160 = make_shared<Tensor>(Gamma160_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor29 = vector<shared_ptr<Tensor>>{Gamma160, rdm0_, rdm1_, rdm2_};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  return make_shared<FutureTensor>(*Gamma160, task29);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma99_() {
  vector<IndexRange> Gamma99_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma99 = make_shared<Tensor>(Gamma99_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor30 = vector<shared_ptr<Tensor>>{Gamma99, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  return make_shared<FutureTensor>(*Gamma99, task30);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma105_() {
  vector<IndexRange> Gamma105_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma105 = make_shared<Tensor>(Gamma105_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor31 = vector<shared_ptr<Tensor>>{Gamma105, rdm1_, rdm2_, rdm3_};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  return make_shared<FutureTensor>(*Gamma105, task31);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma110_() {
  vector<IndexRange> Gamma110_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma110 = make_shared<Tensor>(Gamma110_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor32 = vector<shared_ptr<Tensor>>{Gamma110, rdm1_, rdm2_, rdm3_};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  return make_shared<FutureTensor>(*Gamma110, task32);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma128_() {
  vector<IndexRange> Gamma128_index = {active_, active_, active_, active_};
  auto Gamma128 = make_shared<Tensor>(Gamma128_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor33 = vector<shared_ptr<Tensor>>{Gamma128, rdm0_, rdm1_, rdm2_};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  return make_shared<FutureTensor>(*Gamma128, task33);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma151_() {
  vector<IndexRange> Gamma151_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma151 = make_shared<Tensor>(Gamma151_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor34 = vector<shared_ptr<Tensor>>{Gamma151, rdm1_, rdm2_, rdm3_};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  return make_shared<FutureTensor>(*Gamma151, task34);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma159_() {
  vector<IndexRange> Gamma159_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma159 = make_shared<Tensor>(Gamma159_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor35 = vector<shared_ptr<Tensor>>{Gamma159, rdm2_, rdm3_};
  auto task35 = make_shared<Task35>(tensor35, pindex);
  return make_shared<FutureTensor>(*Gamma159, task35);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma174_() {
  vector<IndexRange> Gamma174_index = {active_, active_, active_, active_};
  auto Gamma174 = make_shared<Tensor>(Gamma174_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor36 = vector<shared_ptr<Tensor>>{Gamma174, rdm1_, rdm2_};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  return make_shared<FutureTensor>(*Gamma174, task36);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma416_() {
  vector<IndexRange> Gamma416_index = {active_, active_};
  auto Gamma416 = make_shared<Tensor>(Gamma416_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor37 = vector<shared_ptr<Tensor>>{Gamma416, rdm0_, rdm1_, rdm2_, h1_};
  auto task37 = make_shared<Task37>(tensor37, pindex);
  return make_shared<FutureTensor>(*Gamma416, task37);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma418_() {
  vector<IndexRange> Gamma418_index = {active_, active_};
  auto Gamma418 = make_shared<Tensor>(Gamma418_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor38 = vector<shared_ptr<Tensor>>{Gamma418, rdm0_, rdm1_, rdm2_, rdm3_, v2_};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  return make_shared<FutureTensor>(*Gamma418, task38);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma24_() {
  vector<IndexRange> Gamma24_index = {active_, active_, active_, active_};
  auto Gamma24 = make_shared<Tensor>(Gamma24_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor39 = vector<shared_ptr<Tensor>>{Gamma24, rdm1_, rdm2_};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  return make_shared<FutureTensor>(*Gamma24, task39);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma26_() {
  vector<IndexRange> Gamma26_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma26 = make_shared<Tensor>(Gamma26_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor40 = vector<shared_ptr<Tensor>>{Gamma26, rdm2_, rdm3_};
  auto task40 = make_shared<Task40>(tensor40, pindex);
  return make_shared<FutureTensor>(*Gamma26, task40);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma27_() {
  vector<IndexRange> Gamma27_index = {active_, active_};
  auto Gamma27 = make_shared<Tensor>(Gamma27_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor41 = vector<shared_ptr<Tensor>>{Gamma27, rdm1_};
  auto task41 = make_shared<Task41>(tensor41, pindex);
  return make_shared<FutureTensor>(*Gamma27, task41);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma179_() {
  vector<IndexRange> Gamma179_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma179 = make_shared<Tensor>(Gamma179_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor42 = vector<shared_ptr<Tensor>>{Gamma179, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task42 = make_shared<Task42>(tensor42, pindex);
  return make_shared<FutureTensor>(*Gamma179, task42);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma183_() {
  vector<IndexRange> Gamma183_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma183 = make_shared<Tensor>(Gamma183_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor43 = vector<shared_ptr<Tensor>>{Gamma183, rdm1_, rdm2_, rdm3_};
  auto task43 = make_shared<Task43>(tensor43, pindex);
  return make_shared<FutureTensor>(*Gamma183, task43);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma193_() {
  vector<IndexRange> Gamma193_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma193 = make_shared<Tensor>(Gamma193_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor44 = vector<shared_ptr<Tensor>>{Gamma193, rdm1_, rdm2_, rdm3_};
  auto task44 = make_shared<Task44>(tensor44, pindex);
  return make_shared<FutureTensor>(*Gamma193, task44);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma203_() {
  vector<IndexRange> Gamma203_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma203 = make_shared<Tensor>(Gamma203_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor45 = vector<shared_ptr<Tensor>>{Gamma203, rdm2_, rdm3_, rdm4_};
  auto task45 = make_shared<Task45>(tensor45, pindex);
  return make_shared<FutureTensor>(*Gamma203, task45);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma204_() {
  vector<IndexRange> Gamma204_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma204 = make_shared<Tensor>(Gamma204_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor46 = vector<shared_ptr<Tensor>>{Gamma204, rdm2_, rdm3_, rdm4_};
  auto task46 = make_shared<Task46>(tensor46, pindex);
  return make_shared<FutureTensor>(*Gamma204, task46);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma229_() {
  vector<IndexRange> Gamma229_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma229 = make_shared<Tensor>(Gamma229_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor47 = vector<shared_ptr<Tensor>>{Gamma229, rdm2_, rdm3_};
  auto task47 = make_shared<Task47>(tensor47, pindex);
  return make_shared<FutureTensor>(*Gamma229, task47);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma420_() {
  vector<IndexRange> Gamma420_index = {active_, active_, active_, active_};
  auto Gamma420 = make_shared<Tensor>(Gamma420_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor48 = vector<shared_ptr<Tensor>>{Gamma420, rdm1_, rdm2_, rdm3_, h1_};
  auto task48 = make_shared<Task48>(tensor48, pindex);
  return make_shared<FutureTensor>(*Gamma420, task48);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma421_() {
  vector<IndexRange> Gamma421_index = {active_, active_, active_, active_};
  auto Gamma421 = make_shared<Tensor>(Gamma421_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor49 = vector<shared_ptr<Tensor>>{Gamma421, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task49 = make_shared<Task49>(tensor49, pindex);
  return make_shared<FutureTensor>(*Gamma421, task49);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma31_() {
  vector<IndexRange> Gamma31_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma31 = make_shared<Tensor>(Gamma31_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor50 = vector<shared_ptr<Tensor>>{Gamma31, rdm2_, rdm3_};
  auto task50 = make_shared<Task50>(tensor50, pindex);
  return make_shared<FutureTensor>(*Gamma31, task50);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma32_() {
  vector<IndexRange> Gamma32_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma32 = make_shared<Tensor>(Gamma32_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor51 = vector<shared_ptr<Tensor>>{Gamma32, rdm2_, rdm3_};
  auto task51 = make_shared<Task51>(tensor51, pindex);
  return make_shared<FutureTensor>(*Gamma32, task51);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma33_() {
  vector<IndexRange> Gamma33_index = {active_, active_, active_, active_};
  auto Gamma33 = make_shared<Tensor>(Gamma33_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor52 = vector<shared_ptr<Tensor>>{Gamma33, rdm2_};
  auto task52 = make_shared<Task52>(tensor52, pindex);
  return make_shared<FutureTensor>(*Gamma33, task52);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma230_() {
  vector<IndexRange> Gamma230_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma230 = make_shared<Tensor>(Gamma230_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor53 = vector<shared_ptr<Tensor>>{Gamma230, rdm2_, rdm3_, rdm4_};
  auto task53 = make_shared<Task53>(tensor53, pindex);
  return make_shared<FutureTensor>(*Gamma230, task53);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma231_() {
  vector<IndexRange> Gamma231_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma231 = make_shared<Tensor>(Gamma231_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor54 = vector<shared_ptr<Tensor>>{Gamma231, rdm2_, rdm3_};
  auto task54 = make_shared<Task54>(tensor54, pindex);
  return make_shared<FutureTensor>(*Gamma231, task54);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma232_() {
  vector<IndexRange> Gamma232_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma232 = make_shared<Tensor>(Gamma232_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor55 = vector<shared_ptr<Tensor>>{Gamma232, rdm2_, rdm3_, rdm4_};
  auto task55 = make_shared<Task55>(tensor55, pindex);
  return make_shared<FutureTensor>(*Gamma232, task55);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma233_() {
  vector<IndexRange> Gamma233_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma233 = make_shared<Tensor>(Gamma233_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor56 = vector<shared_ptr<Tensor>>{Gamma233, rdm2_, rdm3_, rdm4_};
  auto task56 = make_shared<Task56>(tensor56, pindex);
  return make_shared<FutureTensor>(*Gamma233, task56);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma236_() {
  vector<IndexRange> Gamma236_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma236 = make_shared<Tensor>(Gamma236_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor57 = vector<shared_ptr<Tensor>>{Gamma236, rdm2_, rdm3_, rdm4_};
  auto task57 = make_shared<Task57>(tensor57, pindex);
  return make_shared<FutureTensor>(*Gamma236, task57);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma237_() {
  vector<IndexRange> Gamma237_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma237 = make_shared<Tensor>(Gamma237_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor58 = vector<shared_ptr<Tensor>>{Gamma237, rdm2_, rdm3_, rdm4_};
  auto task58 = make_shared<Task58>(tensor58, pindex);
  return make_shared<FutureTensor>(*Gamma237, task58);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma238_() {
  vector<IndexRange> Gamma238_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma238 = make_shared<Tensor>(Gamma238_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor59 = vector<shared_ptr<Tensor>>{Gamma238, rdm2_, rdm3_, rdm4_};
  auto task59 = make_shared<Task59>(tensor59, pindex);
  return make_shared<FutureTensor>(*Gamma238, task59);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma245_() {
  vector<IndexRange> Gamma245_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma245 = make_shared<Tensor>(Gamma245_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor60 = vector<shared_ptr<Tensor>>{Gamma245, rdm2_, rdm3_};
  auto task60 = make_shared<Task60>(tensor60, pindex);
  return make_shared<FutureTensor>(*Gamma245, task60);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma246_() {
  vector<IndexRange> Gamma246_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma246 = make_shared<Tensor>(Gamma246_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor61 = vector<shared_ptr<Tensor>>{Gamma246, rdm2_, rdm3_};
  auto task61 = make_shared<Task61>(tensor61, pindex);
  return make_shared<FutureTensor>(*Gamma246, task61);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma253_() {
  vector<IndexRange> Gamma253_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma253 = make_shared<Tensor>(Gamma253_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor62 = vector<shared_ptr<Tensor>>{Gamma253, rdm3_, rdm4_};
  auto task62 = make_shared<Task62>(tensor62, pindex);
  return make_shared<FutureTensor>(*Gamma253, task62);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma422_() {
  vector<IndexRange> Gamma422_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma422 = make_shared<Tensor>(Gamma422_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor63 = vector<shared_ptr<Tensor>>{Gamma422, rdm2_, rdm3_, rdm4_, h1_};
  auto task63 = make_shared<Task63>(tensor63, pindex);
  return make_shared<FutureTensor>(*Gamma422, task63);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma423_() {
  vector<IndexRange> Gamma423_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma423 = make_shared<Tensor>(Gamma423_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor64 = vector<shared_ptr<Tensor>>{Gamma423, rdm2_, rdm3_, rdm4_, v2_};
  auto task64 = make_shared<Task64>(tensor64, pindex);
  return make_shared<FutureTensor>(*Gamma423, task64);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma317_() {
  vector<IndexRange> Gamma317_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma317 = make_shared<Tensor>(Gamma317_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor65 = vector<shared_ptr<Tensor>>{Gamma317, rdm1_, rdm2_, rdm3_};
  auto task65 = make_shared<Task65>(tensor65, pindex);
  return make_shared<FutureTensor>(*Gamma317, task65);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma318_() {
  vector<IndexRange> Gamma318_index = {active_, active_, active_, active_};
  auto Gamma318 = make_shared<Tensor>(Gamma318_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor66 = vector<shared_ptr<Tensor>>{Gamma318, rdm1_, rdm2_};
  auto task66 = make_shared<Task66>(tensor66, pindex);
  return make_shared<FutureTensor>(*Gamma318, task66);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma335_() {
  vector<IndexRange> Gamma335_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma335 = make_shared<Tensor>(Gamma335_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor67 = vector<shared_ptr<Tensor>>{Gamma335, rdm2_, rdm3_};
  auto task67 = make_shared<Task67>(tensor67, pindex);
  return make_shared<FutureTensor>(*Gamma335, task67);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma336_() {
  vector<IndexRange> Gamma336_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma336 = make_shared<Tensor>(Gamma336_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor68 = vector<shared_ptr<Tensor>>{Gamma336, rdm2_, rdm3_};
  auto task68 = make_shared<Task68>(tensor68, pindex);
  return make_shared<FutureTensor>(*Gamma336, task68);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma368_() {
  vector<IndexRange> Gamma368_index = {active_, active_, active_, active_};
  auto Gamma368 = make_shared<Tensor>(Gamma368_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor69 = vector<shared_ptr<Tensor>>{Gamma368, rdm2_};
  auto task69 = make_shared<Task69>(tensor69, pindex);
  return make_shared<FutureTensor>(*Gamma368, task69);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma363_() {
  vector<IndexRange> Gamma363_index = {active_, active_, active_, active_};
  auto Gamma363 = make_shared<Tensor>(Gamma363_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor70 = vector<shared_ptr<Tensor>>{Gamma363, rdm1_, rdm2_};
  auto task70 = make_shared<Task70>(tensor70, pindex);
  return make_shared<FutureTensor>(*Gamma363, task70);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma391_() {
  vector<IndexRange> Gamma391_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma391 = make_shared<Tensor>(Gamma391_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor71 = vector<shared_ptr<Tensor>>{Gamma391, rdm3_};
  auto task71 = make_shared<Task71>(tensor71, pindex);
  return make_shared<FutureTensor>(*Gamma391, task71);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma428_() {
  vector<IndexRange> Gamma428_index = {active_, active_};
  auto Gamma428 = make_shared<Tensor>(Gamma428_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor72 = vector<shared_ptr<Tensor>>{Gamma428, rdm1_, rdm2_, h1_};
  auto task72 = make_shared<Task72>(tensor72, pindex);
  return make_shared<FutureTensor>(*Gamma428, task72);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma430_() {
  vector<IndexRange> Gamma430_index = {active_, active_};
  auto Gamma430 = make_shared<Tensor>(Gamma430_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor73 = vector<shared_ptr<Tensor>>{Gamma430, rdm1_, rdm2_, rdm3_, v2_};
  auto task73 = make_shared<Task73>(tensor73, pindex);
  return make_shared<FutureTensor>(*Gamma430, task73);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma396_() {
  vector<IndexRange> Gamma396_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma396 = make_shared<Tensor>(Gamma396_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor74 = vector<shared_ptr<Tensor>>{Gamma396, rdm2_, rdm3_};
  auto task74 = make_shared<Task74>(tensor74, pindex);
  return make_shared<FutureTensor>(*Gamma396, task74);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma397_() {
  vector<IndexRange> Gamma397_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma397 = make_shared<Tensor>(Gamma397_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor75 = vector<shared_ptr<Tensor>>{Gamma397, rdm3_, rdm4_};
  auto task75 = make_shared<Task75>(tensor75, pindex);
  return make_shared<FutureTensor>(*Gamma397, task75);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma409_() {
  vector<IndexRange> Gamma409_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma409 = make_shared<Tensor>(Gamma409_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor76 = vector<shared_ptr<Tensor>>{Gamma409, rdm3_};
  auto task76 = make_shared<Task76>(tensor76, pindex);
  return make_shared<FutureTensor>(*Gamma409, task76);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma412_() {
  vector<IndexRange> Gamma412_index = {active_, active_, active_, active_};
  auto Gamma412 = make_shared<Tensor>(Gamma412_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor77 = vector<shared_ptr<Tensor>>{Gamma412, rdm0_, rdm1_, rdm2_, rdm3_, h1_};
  auto task77 = make_shared<Task77>(tensor77, pindex);
  return make_shared<FutureTensor>(*Gamma412, task77);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma413_() {
  vector<IndexRange> Gamma413_index = {active_, active_, active_, active_};
  auto Gamma413 = make_shared<Tensor>(Gamma413_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor78 = vector<shared_ptr<Tensor>>{Gamma413, rdm0_, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task78 = make_shared<Task78>(tensor78, pindex);
  return make_shared<FutureTensor>(*Gamma413, task78);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma424_() {
  vector<IndexRange> Gamma424_index;
  auto Gamma424 = make_shared<Tensor>(Gamma424_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor79 = vector<shared_ptr<Tensor>>{Gamma424, rdm1_, h1_};
  auto task79 = make_shared<Task79>(tensor79, pindex);
  return make_shared<FutureTensor>(*Gamma424, task79);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma426_() {
  vector<IndexRange> Gamma426_index;
  auto Gamma426 = make_shared<Tensor>(Gamma426_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor80 = vector<shared_ptr<Tensor>>{Gamma426, rdm1_, rdm2_, v2_};
  auto task80 = make_shared<Task80>(tensor80, pindex);
  return make_shared<FutureTensor>(*Gamma426, task80);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma432_() {
  vector<IndexRange> Gamma432_index = {active_, active_, active_, active_};
  auto Gamma432 = make_shared<Tensor>(Gamma432_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor81 = vector<shared_ptr<Tensor>>{Gamma432, rdm2_, rdm3_, h1_};
  auto task81 = make_shared<Task81>(tensor81, pindex);
  return make_shared<FutureTensor>(*Gamma432, task81);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma433_() {
  vector<IndexRange> Gamma433_index = {active_, active_, active_, active_};
  auto Gamma433 = make_shared<Tensor>(Gamma433_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor82 = vector<shared_ptr<Tensor>>{Gamma433, rdm2_, rdm3_, rdm4_, v2_};
  auto task82 = make_shared<Task82>(tensor82, pindex);
  return make_shared<FutureTensor>(*Gamma433, task82);
}

#endif
