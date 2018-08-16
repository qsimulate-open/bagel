//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: SPCASPT2_gamma.cc
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
#include <src/smith/caspt2/SPCASPT2.h>
#include <src/smith/caspt2/SPCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using bagel::SMITH::SPCASPT2::FutureTensor;

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma0_() {
  vector<IndexRange> Gamma0_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma0 = make_shared<Tensor>(Gamma0_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor0 = vector<shared_ptr<Tensor>>{Gamma0, ardm1_, ardm2_, ardm3_, rdm0_, rdm1_};
  auto task0 = make_shared<Task0>(tensor0, pindex);
  return make_shared<FutureTensor>(*Gamma0, task0);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma31_() {
  vector<IndexRange> Gamma31_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma31 = make_shared<Tensor>(Gamma31_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor1 = vector<shared_ptr<Tensor>>{Gamma31, ardm2_, ardm3_, rdm1_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  return make_shared<FutureTensor>(*Gamma31, task1);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma34_() {
  vector<IndexRange> Gamma34_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma34 = make_shared<Tensor>(Gamma34_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor2 = vector<shared_ptr<Tensor>>{Gamma34, ardm1_, ardm2_, ardm3_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  return make_shared<FutureTensor>(*Gamma34, task2);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma92_() {
  vector<IndexRange> Gamma92_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma92 = make_shared<Tensor>(Gamma92_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor3 = vector<shared_ptr<Tensor>>{Gamma92, ardm3_};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  return make_shared<FutureTensor>(*Gamma92, task3);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma1_() {
  vector<IndexRange> Gamma1_index = {active_, active_, active_, active_};
  auto Gamma1 = make_shared<Tensor>(Gamma1_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor4 = vector<shared_ptr<Tensor>>{Gamma1, ardm1_, ardm2_, rdm0_, rdm1_};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  return make_shared<FutureTensor>(*Gamma1, task4);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma32_() {
  vector<IndexRange> Gamma32_index = {active_, active_, active_, active_};
  auto Gamma32 = make_shared<Tensor>(Gamma32_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor5 = vector<shared_ptr<Tensor>>{Gamma32, ardm2_, rdm1_};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  return make_shared<FutureTensor>(*Gamma32, task5);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma35_() {
  vector<IndexRange> Gamma35_index = {active_, active_, active_, active_};
  auto Gamma35 = make_shared<Tensor>(Gamma35_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor6 = vector<shared_ptr<Tensor>>{Gamma35, ardm1_, ardm2_};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  return make_shared<FutureTensor>(*Gamma35, task6);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma2_() {
  vector<IndexRange> Gamma2_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma2 = make_shared<Tensor>(Gamma2_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor7 = vector<shared_ptr<Tensor>>{Gamma2, ardm1_, ardm2_, ardm3_, rdm1_, rdm2_};
  auto task7 = make_shared<Task7>(tensor7, pindex);
  return make_shared<FutureTensor>(*Gamma2, task7);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma37_() {
  vector<IndexRange> Gamma37_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma37 = make_shared<Tensor>(Gamma37_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor8 = vector<shared_ptr<Tensor>>{Gamma37, ardm2_, ardm3_, rdm2_};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  return make_shared<FutureTensor>(*Gamma37, task8);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma3_() {
  vector<IndexRange> Gamma3_index = {active_, active_, active_, active_};
  auto Gamma3 = make_shared<Tensor>(Gamma3_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor9 = vector<shared_ptr<Tensor>>{Gamma3, ardm1_, ardm2_, rdm0_, rdm1_};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  return make_shared<FutureTensor>(*Gamma3, task9);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma40_() {
  vector<IndexRange> Gamma40_index = {active_, active_, active_, active_};
  auto Gamma40 = make_shared<Tensor>(Gamma40_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor10 = vector<shared_ptr<Tensor>>{Gamma40, ardm1_, ardm2_};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  return make_shared<FutureTensor>(*Gamma40, task10);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma4_() {
  vector<IndexRange> Gamma4_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma4 = make_shared<Tensor>(Gamma4_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor11 = vector<shared_ptr<Tensor>>{Gamma4, ardm1_, ardm2_, ardm3_, rdm1_, rdm2_};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  return make_shared<FutureTensor>(*Gamma4, task11);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma56_() {
  vector<IndexRange> Gamma56_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma56 = make_shared<Tensor>(Gamma56_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor12 = vector<shared_ptr<Tensor>>{Gamma56, ardm2_, ardm3_, rdm2_};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  return make_shared<FutureTensor>(*Gamma56, task12);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma57_() {
  vector<IndexRange> Gamma57_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma57 = make_shared<Tensor>(Gamma57_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor13 = vector<shared_ptr<Tensor>>{Gamma57, ardm2_, ardm3_};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  return make_shared<FutureTensor>(*Gamma57, task13);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma5_() {
  vector<IndexRange> Gamma5_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma5 = make_shared<Tensor>(Gamma5_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor14 = vector<shared_ptr<Tensor>>{Gamma5, ardm1_, ardm2_, ardm3_, ardm4_, rdm1_, rdm2_};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  return make_shared<FutureTensor>(*Gamma5, task14);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma58_() {
  vector<IndexRange> Gamma58_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma58 = make_shared<Tensor>(Gamma58_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor15 = vector<shared_ptr<Tensor>>{Gamma58, ardm2_, ardm3_, ardm4_};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  return make_shared<FutureTensor>(*Gamma58, task15);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma6_() {
  vector<IndexRange> Gamma6_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma6 = make_shared<Tensor>(Gamma6_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor16 = vector<shared_ptr<Tensor>>{Gamma6, ardm1_, ardm2_, ardm3_, rdm1_, rdm2_};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  return make_shared<FutureTensor>(*Gamma6, task16);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma7_() {
  vector<IndexRange> Gamma7_index = {active_, active_, active_, active_};
  auto Gamma7 = make_shared<Tensor>(Gamma7_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor17 = vector<shared_ptr<Tensor>>{Gamma7, rdm1_, rdm2_};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  return make_shared<FutureTensor>(*Gamma7, task17);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma8_() {
  vector<IndexRange> Gamma8_index = {active_, active_, active_, active_};
  auto Gamma8 = make_shared<Tensor>(Gamma8_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor18 = vector<shared_ptr<Tensor>>{Gamma8, ardm1_, ardm2_, rdm1_};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  return make_shared<FutureTensor>(*Gamma8, task18);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma60_() {
  vector<IndexRange> Gamma60_index = {active_, active_, active_, active_};
  auto Gamma60 = make_shared<Tensor>(Gamma60_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor19 = vector<shared_ptr<Tensor>>{Gamma60, ardm2_};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  return make_shared<FutureTensor>(*Gamma60, task19);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma61_() {
  vector<IndexRange> Gamma61_index = {active_, active_, active_, active_};
  auto Gamma61 = make_shared<Tensor>(Gamma61_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor20 = vector<shared_ptr<Tensor>>{Gamma61, rdm2_};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  return make_shared<FutureTensor>(*Gamma61, task20);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma9_() {
  vector<IndexRange> Gamma9_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma9 = make_shared<Tensor>(Gamma9_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor21 = vector<shared_ptr<Tensor>>{Gamma9, ardm1_, ardm2_, ardm3_};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  return make_shared<FutureTensor>(*Gamma9, task21);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma62_() {
  vector<IndexRange> Gamma62_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma62 = make_shared<Tensor>(Gamma62_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor22 = vector<shared_ptr<Tensor>>{Gamma62, ardm2_, ardm3_};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  return make_shared<FutureTensor>(*Gamma62, task22);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma11_() {
  vector<IndexRange> Gamma11_index = {active_, active_, active_, active_};
  auto Gamma11 = make_shared<Tensor>(Gamma11_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor23 = vector<shared_ptr<Tensor>>{Gamma11, ardm1_, ardm2_, rdm0_, rdm1_};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  return make_shared<FutureTensor>(*Gamma11, task23);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma12_() {
  vector<IndexRange> Gamma12_index = {active_, active_, active_, active_};
  auto Gamma12 = make_shared<Tensor>(Gamma12_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor24 = vector<shared_ptr<Tensor>>{Gamma12, rdm1_, rdm2_};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  return make_shared<FutureTensor>(*Gamma12, task24);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma13_() {
  vector<IndexRange> Gamma13_index = {active_, active_, active_, active_};
  auto Gamma13 = make_shared<Tensor>(Gamma13_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor25 = vector<shared_ptr<Tensor>>{Gamma13, ardm1_, ardm2_, rdm1_};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  return make_shared<FutureTensor>(*Gamma13, task25);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma65_() {
  vector<IndexRange> Gamma65_index = {active_, active_};
  auto Gamma65 = make_shared<Tensor>(Gamma65_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor26 = vector<shared_ptr<Tensor>>{Gamma65, ardm1_};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  return make_shared<FutureTensor>(*Gamma65, task26);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma67_() {
  vector<IndexRange> Gamma67_index = {active_, active_};
  auto Gamma67 = make_shared<Tensor>(Gamma67_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor27 = vector<shared_ptr<Tensor>>{Gamma67, rdm1_};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  return make_shared<FutureTensor>(*Gamma67, task27);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma14_() {
  vector<IndexRange> Gamma14_index = {active_, active_, active_, active_};
  auto Gamma14 = make_shared<Tensor>(Gamma14_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor28 = vector<shared_ptr<Tensor>>{Gamma14, ardm1_, ardm2_, rdm0_};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  return make_shared<FutureTensor>(*Gamma14, task28);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma81_() {
  vector<IndexRange> Gamma81_index = {active_, active_, active_, active_};
  auto Gamma81 = make_shared<Tensor>(Gamma81_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor29 = vector<shared_ptr<Tensor>>{Gamma81, ardm2_};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  return make_shared<FutureTensor>(*Gamma81, task29);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma16_() {
  vector<IndexRange> Gamma16_index = {active_, active_};
  auto Gamma16 = make_shared<Tensor>(Gamma16_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor30 = vector<shared_ptr<Tensor>>{Gamma16, rdm0_, rdm1_};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  return make_shared<FutureTensor>(*Gamma16, task30);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma17_() {
  vector<IndexRange> Gamma17_index = {active_, active_};
  auto Gamma17 = make_shared<Tensor>(Gamma17_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor31 = vector<shared_ptr<Tensor>>{Gamma17, ardm1_, rdm0_};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  return make_shared<FutureTensor>(*Gamma17, task31);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma22_() {
  vector<IndexRange> Gamma22_index = {active_, active_, active_, active_};
  auto Gamma22 = make_shared<Tensor>(Gamma22_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor32 = vector<shared_ptr<Tensor>>{Gamma22, ardm1_, ardm2_};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  return make_shared<FutureTensor>(*Gamma22, task32);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma28_() {
  vector<IndexRange> Gamma28_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma28 = make_shared<Tensor>(Gamma28_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor33 = vector<shared_ptr<Tensor>>{Gamma28, ardm1_, ardm2_, ardm3_};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  return make_shared<FutureTensor>(*Gamma28, task33);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma29_() {
  vector<IndexRange> Gamma29_index = {active_, active_, active_, active_};
  auto Gamma29 = make_shared<Tensor>(Gamma29_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor34 = vector<shared_ptr<Tensor>>{Gamma29, ardm1_, ardm2_};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  return make_shared<FutureTensor>(*Gamma29, task34);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma33_() {
  vector<IndexRange> Gamma33_index = {active_, active_, active_, active_};
  auto Gamma33 = make_shared<Tensor>(Gamma33_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor35 = vector<shared_ptr<Tensor>>{Gamma33, ardm1_, ardm2_};
  auto task35 = make_shared<Task35>(tensor35, pindex);
  return make_shared<FutureTensor>(*Gamma33, task35);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma90_() {
  vector<IndexRange> Gamma90_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma90 = make_shared<Tensor>(Gamma90_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor36 = vector<shared_ptr<Tensor>>{Gamma90, ardm2_, ardm3_};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  return make_shared<FutureTensor>(*Gamma90, task36);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma49_() {
  vector<IndexRange> Gamma49_index = {active_, active_, active_, active_};
  auto Gamma49 = make_shared<Tensor>(Gamma49_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor37 = vector<shared_ptr<Tensor>>{Gamma49, rdm1_, rdm2_};
  auto task37 = make_shared<Task37>(tensor37, pindex);
  return make_shared<FutureTensor>(*Gamma49, task37);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma51_() {
  vector<IndexRange> Gamma51_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma51 = make_shared<Tensor>(Gamma51_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor38 = vector<shared_ptr<Tensor>>{Gamma51, ardm2_, ardm3_};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  return make_shared<FutureTensor>(*Gamma51, task38);
}

shared_ptr<FutureTensor> SPCASPT2::SPCASPT2::Gamma59_() {
  vector<IndexRange> Gamma59_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma59 = make_shared<Tensor>(Gamma59_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor39 = vector<shared_ptr<Tensor>>{Gamma59, ardm2_, ardm3_};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  return make_shared<FutureTensor>(*Gamma59, task39);
}

#endif
