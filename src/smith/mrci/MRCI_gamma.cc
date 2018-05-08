//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_gamma.cc
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
#include <src/smith/mrci/MRCI.h>
#include <src/smith/mrci/MRCI_tasks1.h>
#include <src/smith/mrci/MRCI_tasks2.h>
#include <src/smith/mrci/MRCI_tasks3.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using bagel::SMITH::MRCI::FutureTensor;

shared_ptr<FutureTensor> MRCI::MRCI::Gamma0_() {
  vector<IndexRange> Gamma0_index = {active_, active_, active_, active_};
  auto Gamma0 = make_shared<Tensor>(Gamma0_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor0 = vector<shared_ptr<Tensor>>{Gamma0, rdm0_, rdm1_, rdm2_};
  auto task0 = make_shared<Task0>(tensor0, pindex);
  return make_shared<FutureTensor>(*Gamma0, task0);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma1_() {
  vector<IndexRange> Gamma1_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma1 = make_shared<Tensor>(Gamma1_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor1 = vector<shared_ptr<Tensor>>{Gamma1, rdm1_, rdm2_, rdm3_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  return make_shared<FutureTensor>(*Gamma1, task1);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma2_() {
  vector<IndexRange> Gamma2_index = {active_, active_, active_, active_};
  auto Gamma2 = make_shared<Tensor>(Gamma2_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor2 = vector<shared_ptr<Tensor>>{Gamma2, rdm0_, rdm1_, rdm2_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  return make_shared<FutureTensor>(*Gamma2, task2);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma80_() {
  vector<IndexRange> Gamma80_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma80 = make_shared<Tensor>(Gamma80_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor3 = vector<shared_ptr<Tensor>>{Gamma80, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  return make_shared<FutureTensor>(*Gamma80, task3);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma81_() {
  vector<IndexRange> Gamma81_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma81 = make_shared<Tensor>(Gamma81_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor4 = vector<shared_ptr<Tensor>>{Gamma81, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  return make_shared<FutureTensor>(*Gamma81, task4);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma82_() {
  vector<IndexRange> Gamma82_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma82 = make_shared<Tensor>(Gamma82_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor5 = vector<shared_ptr<Tensor>>{Gamma82, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  return make_shared<FutureTensor>(*Gamma82, task5);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma85_() {
  vector<IndexRange> Gamma85_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma85 = make_shared<Tensor>(Gamma85_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor6 = vector<shared_ptr<Tensor>>{Gamma85, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  return make_shared<FutureTensor>(*Gamma85, task6);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma86_() {
  vector<IndexRange> Gamma86_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma86 = make_shared<Tensor>(Gamma86_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor7 = vector<shared_ptr<Tensor>>{Gamma86, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task7 = make_shared<Task7>(tensor7, pindex);
  return make_shared<FutureTensor>(*Gamma86, task7);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma87_() {
  vector<IndexRange> Gamma87_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma87 = make_shared<Tensor>(Gamma87_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor8 = vector<shared_ptr<Tensor>>{Gamma87, rdm1_, rdm2_, rdm3_};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  return make_shared<FutureTensor>(*Gamma87, task8);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma88_() {
  vector<IndexRange> Gamma88_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma88 = make_shared<Tensor>(Gamma88_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor9 = vector<shared_ptr<Tensor>>{Gamma88, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  return make_shared<FutureTensor>(*Gamma88, task9);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma89_() {
  vector<IndexRange> Gamma89_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma89 = make_shared<Tensor>(Gamma89_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor10 = vector<shared_ptr<Tensor>>{Gamma89, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  return make_shared<FutureTensor>(*Gamma89, task10);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma94_() {
  vector<IndexRange> Gamma94_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma94 = make_shared<Tensor>(Gamma94_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor11 = vector<shared_ptr<Tensor>>{Gamma94, rdm1_, rdm2_, rdm3_};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  return make_shared<FutureTensor>(*Gamma94, task11);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma3_() {
  vector<IndexRange> Gamma3_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma3 = make_shared<Tensor>(Gamma3_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor12 = vector<shared_ptr<Tensor>>{Gamma3, rdm1_, rdm2_, rdm3_};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  return make_shared<FutureTensor>(*Gamma3, task12);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma4_() {
  vector<IndexRange> Gamma4_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma4 = make_shared<Tensor>(Gamma4_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor13 = vector<shared_ptr<Tensor>>{Gamma4, rdm1_, rdm2_, rdm3_};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  return make_shared<FutureTensor>(*Gamma4, task13);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma5_() {
  vector<IndexRange> Gamma5_index = {active_, active_, active_, active_};
  auto Gamma5 = make_shared<Tensor>(Gamma5_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor14 = vector<shared_ptr<Tensor>>{Gamma5, rdm1_, rdm2_};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  return make_shared<FutureTensor>(*Gamma5, task14);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma7_() {
  vector<IndexRange> Gamma7_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma7 = make_shared<Tensor>(Gamma7_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor15 = vector<shared_ptr<Tensor>>{Gamma7, rdm1_, rdm2_, rdm3_};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  return make_shared<FutureTensor>(*Gamma7, task15);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma97_() {
  vector<IndexRange> Gamma97_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma97 = make_shared<Tensor>(Gamma97_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor16 = vector<shared_ptr<Tensor>>{Gamma97, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  return make_shared<FutureTensor>(*Gamma97, task16);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma98_() {
  vector<IndexRange> Gamma98_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma98 = make_shared<Tensor>(Gamma98_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor17 = vector<shared_ptr<Tensor>>{Gamma98, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  return make_shared<FutureTensor>(*Gamma98, task17);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma100_() {
  vector<IndexRange> Gamma100_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma100 = make_shared<Tensor>(Gamma100_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor18 = vector<shared_ptr<Tensor>>{Gamma100, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  return make_shared<FutureTensor>(*Gamma100, task18);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma101_() {
  vector<IndexRange> Gamma101_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma101 = make_shared<Tensor>(Gamma101_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor19 = vector<shared_ptr<Tensor>>{Gamma101, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  return make_shared<FutureTensor>(*Gamma101, task19);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma102_() {
  vector<IndexRange> Gamma102_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma102 = make_shared<Tensor>(Gamma102_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor20 = vector<shared_ptr<Tensor>>{Gamma102, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  return make_shared<FutureTensor>(*Gamma102, task20);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma104_() {
  vector<IndexRange> Gamma104_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma104 = make_shared<Tensor>(Gamma104_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor21 = vector<shared_ptr<Tensor>>{Gamma104, rdm1_, rdm2_, rdm3_};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  return make_shared<FutureTensor>(*Gamma104, task21);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma107_() {
  vector<IndexRange> Gamma107_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma107 = make_shared<Tensor>(Gamma107_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor22 = vector<shared_ptr<Tensor>>{Gamma107, rdm1_, rdm2_, rdm3_};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  return make_shared<FutureTensor>(*Gamma107, task22);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma109_() {
  vector<IndexRange> Gamma109_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma109 = make_shared<Tensor>(Gamma109_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor23 = vector<shared_ptr<Tensor>>{Gamma109, rdm1_, rdm2_, rdm3_};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  return make_shared<FutureTensor>(*Gamma109, task23);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma114_() {
  vector<IndexRange> Gamma114_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma114 = make_shared<Tensor>(Gamma114_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor24 = vector<shared_ptr<Tensor>>{Gamma114, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  return make_shared<FutureTensor>(*Gamma114, task24);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma115_() {
  vector<IndexRange> Gamma115_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma115 = make_shared<Tensor>(Gamma115_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor25 = vector<shared_ptr<Tensor>>{Gamma115, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  return make_shared<FutureTensor>(*Gamma115, task25);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma119_() {
  vector<IndexRange> Gamma119_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma119 = make_shared<Tensor>(Gamma119_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor26 = vector<shared_ptr<Tensor>>{Gamma119, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  return make_shared<FutureTensor>(*Gamma119, task26);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma122_() {
  vector<IndexRange> Gamma122_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma122 = make_shared<Tensor>(Gamma122_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor27 = vector<shared_ptr<Tensor>>{Gamma122, rdm2_, rdm3_, rdm4_};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  return make_shared<FutureTensor>(*Gamma122, task27);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma550_() {
  vector<IndexRange> Gamma550_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma550 = make_shared<Tensor>(Gamma550_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor28 = vector<shared_ptr<Tensor>>{Gamma550, rdm1_, rdm2_, rdm3_, rdm4_, h1_};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  return make_shared<FutureTensor>(*Gamma550, task28);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma551_() {
  vector<IndexRange> Gamma551_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma551 = make_shared<Tensor>(Gamma551_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor29 = vector<shared_ptr<Tensor>>{Gamma551, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  return make_shared<FutureTensor>(*Gamma551, task29);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma10_() {
  vector<IndexRange> Gamma10_index = {active_, active_, active_, active_};
  auto Gamma10 = make_shared<Tensor>(Gamma10_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor30 = vector<shared_ptr<Tensor>>{Gamma10, rdm1_, rdm2_};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  return make_shared<FutureTensor>(*Gamma10, task30);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma12_() {
  vector<IndexRange> Gamma12_index = {active_, active_};
  auto Gamma12 = make_shared<Tensor>(Gamma12_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor31 = vector<shared_ptr<Tensor>>{Gamma12, rdm0_, rdm1_};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  return make_shared<FutureTensor>(*Gamma12, task31);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma18_() {
  vector<IndexRange> Gamma18_index = {active_, active_, active_, active_};
  auto Gamma18 = make_shared<Tensor>(Gamma18_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor32 = vector<shared_ptr<Tensor>>{Gamma18, rdm1_, rdm2_};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  return make_shared<FutureTensor>(*Gamma18, task32);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma197_() {
  vector<IndexRange> Gamma197_index = {active_, active_, active_, active_};
  auto Gamma197 = make_shared<Tensor>(Gamma197_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor33 = vector<shared_ptr<Tensor>>{Gamma197, rdm0_, rdm1_, rdm2_};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  return make_shared<FutureTensor>(*Gamma197, task33);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma126_() {
  vector<IndexRange> Gamma126_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma126 = make_shared<Tensor>(Gamma126_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor34 = vector<shared_ptr<Tensor>>{Gamma126, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  return make_shared<FutureTensor>(*Gamma126, task34);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma132_() {
  vector<IndexRange> Gamma132_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma132 = make_shared<Tensor>(Gamma132_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor35 = vector<shared_ptr<Tensor>>{Gamma132, rdm1_, rdm2_, rdm3_};
  auto task35 = make_shared<Task35>(tensor35, pindex);
  return make_shared<FutureTensor>(*Gamma132, task35);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma137_() {
  vector<IndexRange> Gamma137_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma137 = make_shared<Tensor>(Gamma137_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor36 = vector<shared_ptr<Tensor>>{Gamma137, rdm1_, rdm2_, rdm3_};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  return make_shared<FutureTensor>(*Gamma137, task36);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma155_() {
  vector<IndexRange> Gamma155_index = {active_, active_, active_, active_};
  auto Gamma155 = make_shared<Tensor>(Gamma155_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor37 = vector<shared_ptr<Tensor>>{Gamma155, rdm0_, rdm1_, rdm2_};
  auto task37 = make_shared<Task37>(tensor37, pindex);
  return make_shared<FutureTensor>(*Gamma155, task37);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma176_() {
  vector<IndexRange> Gamma176_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma176 = make_shared<Tensor>(Gamma176_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor38 = vector<shared_ptr<Tensor>>{Gamma176, rdm1_, rdm2_, rdm3_};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  return make_shared<FutureTensor>(*Gamma176, task38);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma178_() {
  vector<IndexRange> Gamma178_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma178 = make_shared<Tensor>(Gamma178_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor39 = vector<shared_ptr<Tensor>>{Gamma178, rdm1_, rdm2_, rdm3_};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  return make_shared<FutureTensor>(*Gamma178, task39);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma179_() {
  vector<IndexRange> Gamma179_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma179 = make_shared<Tensor>(Gamma179_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor40 = vector<shared_ptr<Tensor>>{Gamma179, rdm1_, rdm2_, rdm3_};
  auto task40 = make_shared<Task40>(tensor40, pindex);
  return make_shared<FutureTensor>(*Gamma179, task40);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma196_() {
  vector<IndexRange> Gamma196_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma196 = make_shared<Tensor>(Gamma196_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor41 = vector<shared_ptr<Tensor>>{Gamma196, rdm2_, rdm3_};
  auto task41 = make_shared<Task41>(tensor41, pindex);
  return make_shared<FutureTensor>(*Gamma196, task41);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma552_() {
  vector<IndexRange> Gamma552_index = {active_, active_};
  auto Gamma552 = make_shared<Tensor>(Gamma552_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor42 = vector<shared_ptr<Tensor>>{Gamma552, rdm0_, rdm1_, rdm2_, h1_};
  auto task42 = make_shared<Task42>(tensor42, pindex);
  return make_shared<FutureTensor>(*Gamma552, task42);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma554_() {
  vector<IndexRange> Gamma554_index = {active_, active_};
  auto Gamma554 = make_shared<Tensor>(Gamma554_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor43 = vector<shared_ptr<Tensor>>{Gamma554, rdm0_, rdm1_, rdm2_, rdm3_, v2_};
  auto task43 = make_shared<Task43>(tensor43, pindex);
  return make_shared<FutureTensor>(*Gamma554, task43);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma24_() {
  vector<IndexRange> Gamma24_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma24 = make_shared<Tensor>(Gamma24_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor44 = vector<shared_ptr<Tensor>>{Gamma24, rdm1_, rdm2_, rdm3_};
  auto task44 = make_shared<Task44>(tensor44, pindex);
  return make_shared<FutureTensor>(*Gamma24, task44);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma25_() {
  vector<IndexRange> Gamma25_index = {active_, active_, active_, active_};
  auto Gamma25 = make_shared<Tensor>(Gamma25_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor45 = vector<shared_ptr<Tensor>>{Gamma25, rdm1_, rdm2_};
  auto task45 = make_shared<Task45>(tensor45, pindex);
  return make_shared<FutureTensor>(*Gamma25, task45);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma27_() {
  vector<IndexRange> Gamma27_index = {active_, active_, active_, active_};
  auto Gamma27 = make_shared<Tensor>(Gamma27_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor46 = vector<shared_ptr<Tensor>>{Gamma27, rdm1_, rdm2_};
  auto task46 = make_shared<Task46>(tensor46, pindex);
  return make_shared<FutureTensor>(*Gamma27, task46);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma29_() {
  vector<IndexRange> Gamma29_index = {active_, active_, active_, active_};
  auto Gamma29 = make_shared<Tensor>(Gamma29_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor47 = vector<shared_ptr<Tensor>>{Gamma29, rdm1_, rdm2_};
  auto task47 = make_shared<Task47>(tensor47, pindex);
  return make_shared<FutureTensor>(*Gamma29, task47);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma31_() {
  vector<IndexRange> Gamma31_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma31 = make_shared<Tensor>(Gamma31_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor48 = vector<shared_ptr<Tensor>>{Gamma31, rdm2_, rdm3_};
  auto task48 = make_shared<Task48>(tensor48, pindex);
  return make_shared<FutureTensor>(*Gamma31, task48);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma32_() {
  vector<IndexRange> Gamma32_index = {active_, active_};
  auto Gamma32 = make_shared<Tensor>(Gamma32_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor49 = vector<shared_ptr<Tensor>>{Gamma32, rdm1_};
  auto task49 = make_shared<Task49>(tensor49, pindex);
  return make_shared<FutureTensor>(*Gamma32, task49);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma215_() {
  vector<IndexRange> Gamma215_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma215 = make_shared<Tensor>(Gamma215_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor50 = vector<shared_ptr<Tensor>>{Gamma215, rdm1_, rdm2_, rdm3_};
  auto task50 = make_shared<Task50>(tensor50, pindex);
  return make_shared<FutureTensor>(*Gamma215, task50);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma216_() {
  vector<IndexRange> Gamma216_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma216 = make_shared<Tensor>(Gamma216_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor51 = vector<shared_ptr<Tensor>>{Gamma216, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task51 = make_shared<Task51>(tensor51, pindex);
  return make_shared<FutureTensor>(*Gamma216, task51);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma217_() {
  vector<IndexRange> Gamma217_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma217 = make_shared<Tensor>(Gamma217_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor52 = vector<shared_ptr<Tensor>>{Gamma217, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task52 = make_shared<Task52>(tensor52, pindex);
  return make_shared<FutureTensor>(*Gamma217, task52);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma220_() {
  vector<IndexRange> Gamma220_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma220 = make_shared<Tensor>(Gamma220_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor53 = vector<shared_ptr<Tensor>>{Gamma220, rdm1_, rdm2_, rdm3_};
  auto task53 = make_shared<Task53>(tensor53, pindex);
  return make_shared<FutureTensor>(*Gamma220, task53);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma222_() {
  vector<IndexRange> Gamma222_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma222 = make_shared<Tensor>(Gamma222_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor54 = vector<shared_ptr<Tensor>>{Gamma222, rdm1_, rdm2_, rdm3_};
  auto task54 = make_shared<Task54>(tensor54, pindex);
  return make_shared<FutureTensor>(*Gamma222, task54);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma221_() {
  vector<IndexRange> Gamma221_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma221 = make_shared<Tensor>(Gamma221_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor55 = vector<shared_ptr<Tensor>>{Gamma221, rdm1_, rdm2_, rdm3_};
  auto task55 = make_shared<Task55>(tensor55, pindex);
  return make_shared<FutureTensor>(*Gamma221, task55);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma230_() {
  vector<IndexRange> Gamma230_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma230 = make_shared<Tensor>(Gamma230_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor56 = vector<shared_ptr<Tensor>>{Gamma230, rdm1_, rdm2_, rdm3_};
  auto task56 = make_shared<Task56>(tensor56, pindex);
  return make_shared<FutureTensor>(*Gamma230, task56);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma232_() {
  vector<IndexRange> Gamma232_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma232 = make_shared<Tensor>(Gamma232_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor57 = vector<shared_ptr<Tensor>>{Gamma232, rdm1_, rdm2_, rdm3_};
  auto task57 = make_shared<Task57>(tensor57, pindex);
  return make_shared<FutureTensor>(*Gamma232, task57);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma234_() {
  vector<IndexRange> Gamma234_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma234 = make_shared<Tensor>(Gamma234_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor58 = vector<shared_ptr<Tensor>>{Gamma234, rdm1_, rdm2_, rdm3_};
  auto task58 = make_shared<Task58>(tensor58, pindex);
  return make_shared<FutureTensor>(*Gamma234, task58);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma233_() {
  vector<IndexRange> Gamma233_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma233 = make_shared<Tensor>(Gamma233_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor59 = vector<shared_ptr<Tensor>>{Gamma233, rdm1_, rdm2_, rdm3_};
  auto task59 = make_shared<Task59>(tensor59, pindex);
  return make_shared<FutureTensor>(*Gamma233, task59);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma235_() {
  vector<IndexRange> Gamma235_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma235 = make_shared<Tensor>(Gamma235_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor60 = vector<shared_ptr<Tensor>>{Gamma235, rdm1_, rdm2_, rdm3_};
  auto task60 = make_shared<Task60>(tensor60, pindex);
  return make_shared<FutureTensor>(*Gamma235, task60);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma240_() {
  vector<IndexRange> Gamma240_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma240 = make_shared<Tensor>(Gamma240_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor61 = vector<shared_ptr<Tensor>>{Gamma240, rdm1_, rdm2_, rdm3_};
  auto task61 = make_shared<Task61>(tensor61, pindex);
  return make_shared<FutureTensor>(*Gamma240, task61);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma244_() {
  vector<IndexRange> Gamma244_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma244 = make_shared<Tensor>(Gamma244_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor62 = vector<shared_ptr<Tensor>>{Gamma244, rdm1_, rdm2_, rdm3_};
  auto task62 = make_shared<Task62>(tensor62, pindex);
  return make_shared<FutureTensor>(*Gamma244, task62);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma250_() {
  vector<IndexRange> Gamma250_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma250 = make_shared<Tensor>(Gamma250_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor63 = vector<shared_ptr<Tensor>>{Gamma250, rdm2_, rdm3_, rdm4_};
  auto task63 = make_shared<Task63>(tensor63, pindex);
  return make_shared<FutureTensor>(*Gamma250, task63);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma251_() {
  vector<IndexRange> Gamma251_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma251 = make_shared<Tensor>(Gamma251_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor64 = vector<shared_ptr<Tensor>>{Gamma251, rdm2_, rdm3_, rdm4_};
  auto task64 = make_shared<Task64>(tensor64, pindex);
  return make_shared<FutureTensor>(*Gamma251, task64);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma252_() {
  vector<IndexRange> Gamma252_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma252 = make_shared<Tensor>(Gamma252_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor65 = vector<shared_ptr<Tensor>>{Gamma252, rdm2_, rdm3_};
  auto task65 = make_shared<Task65>(tensor65, pindex);
  return make_shared<FutureTensor>(*Gamma252, task65);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma276_() {
  vector<IndexRange> Gamma276_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma276 = make_shared<Tensor>(Gamma276_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor66 = vector<shared_ptr<Tensor>>{Gamma276, rdm2_, rdm3_};
  auto task66 = make_shared<Task66>(tensor66, pindex);
  return make_shared<FutureTensor>(*Gamma276, task66);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma568_() {
  vector<IndexRange> Gamma568_index = {active_, active_, active_, active_};
  auto Gamma568 = make_shared<Tensor>(Gamma568_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor67 = vector<shared_ptr<Tensor>>{Gamma568, rdm1_, rdm2_, rdm3_, h1_};
  auto task67 = make_shared<Task67>(tensor67, pindex);
  return make_shared<FutureTensor>(*Gamma568, task67);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma569_() {
  vector<IndexRange> Gamma569_index = {active_, active_, active_, active_};
  auto Gamma569 = make_shared<Tensor>(Gamma569_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor68 = vector<shared_ptr<Tensor>>{Gamma569, rdm1_, rdm2_, rdm3_, h1_};
  auto task68 = make_shared<Task68>(tensor68, pindex);
  return make_shared<FutureTensor>(*Gamma569, task68);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma572_() {
  vector<IndexRange> Gamma572_index = {active_, active_, active_, active_};
  auto Gamma572 = make_shared<Tensor>(Gamma572_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor69 = vector<shared_ptr<Tensor>>{Gamma572, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task69 = make_shared<Task69>(tensor69, pindex);
  return make_shared<FutureTensor>(*Gamma572, task69);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma573_() {
  vector<IndexRange> Gamma573_index = {active_, active_, active_, active_};
  auto Gamma573 = make_shared<Tensor>(Gamma573_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor70 = vector<shared_ptr<Tensor>>{Gamma573, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task70 = make_shared<Task70>(tensor70, pindex);
  return make_shared<FutureTensor>(*Gamma573, task70);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma278_() {
  vector<IndexRange> Gamma278_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma278 = make_shared<Tensor>(Gamma278_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor71 = vector<shared_ptr<Tensor>>{Gamma278, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task71 = make_shared<Task71>(tensor71, pindex);
  return make_shared<FutureTensor>(*Gamma278, task71);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma296_() {
  vector<IndexRange> Gamma296_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma296 = make_shared<Tensor>(Gamma296_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor72 = vector<shared_ptr<Tensor>>{Gamma296, rdm1_, rdm2_, rdm3_};
  auto task72 = make_shared<Task72>(tensor72, pindex);
  return make_shared<FutureTensor>(*Gamma296, task72);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma312_() {
  vector<IndexRange> Gamma312_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma312 = make_shared<Tensor>(Gamma312_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor73 = vector<shared_ptr<Tensor>>{Gamma312, rdm2_, rdm3_, rdm4_};
  auto task73 = make_shared<Task73>(tensor73, pindex);
  return make_shared<FutureTensor>(*Gamma312, task73);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma313_() {
  vector<IndexRange> Gamma313_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma313 = make_shared<Tensor>(Gamma313_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor74 = vector<shared_ptr<Tensor>>{Gamma313, rdm2_, rdm3_, rdm4_};
  auto task74 = make_shared<Task74>(tensor74, pindex);
  return make_shared<FutureTensor>(*Gamma313, task74);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma338_() {
  vector<IndexRange> Gamma338_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma338 = make_shared<Tensor>(Gamma338_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor75 = vector<shared_ptr<Tensor>>{Gamma338, rdm2_, rdm3_};
  auto task75 = make_shared<Task75>(tensor75, pindex);
  return make_shared<FutureTensor>(*Gamma338, task75);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma48_() {
  vector<IndexRange> Gamma48_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma48 = make_shared<Tensor>(Gamma48_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor76 = vector<shared_ptr<Tensor>>{Gamma48, rdm2_, rdm3_};
  auto task76 = make_shared<Task76>(tensor76, pindex);
  return make_shared<FutureTensor>(*Gamma48, task76);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma49_() {
  vector<IndexRange> Gamma49_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma49 = make_shared<Tensor>(Gamma49_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor77 = vector<shared_ptr<Tensor>>{Gamma49, rdm2_, rdm3_};
  auto task77 = make_shared<Task77>(tensor77, pindex);
  return make_shared<FutureTensor>(*Gamma49, task77);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma50_() {
  vector<IndexRange> Gamma50_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma50 = make_shared<Tensor>(Gamma50_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor78 = vector<shared_ptr<Tensor>>{Gamma50, rdm2_, rdm3_};
  auto task78 = make_shared<Task78>(tensor78, pindex);
  return make_shared<FutureTensor>(*Gamma50, task78);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma51_() {
  vector<IndexRange> Gamma51_index = {active_, active_, active_, active_};
  auto Gamma51 = make_shared<Tensor>(Gamma51_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor79 = vector<shared_ptr<Tensor>>{Gamma51, rdm2_};
  auto task79 = make_shared<Task79>(tensor79, pindex);
  return make_shared<FutureTensor>(*Gamma51, task79);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma339_() {
  vector<IndexRange> Gamma339_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma339 = make_shared<Tensor>(Gamma339_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor80 = vector<shared_ptr<Tensor>>{Gamma339, rdm2_, rdm3_, rdm4_};
  auto task80 = make_shared<Task80>(tensor80, pindex);
  return make_shared<FutureTensor>(*Gamma339, task80);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma340_() {
  vector<IndexRange> Gamma340_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma340 = make_shared<Tensor>(Gamma340_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor81 = vector<shared_ptr<Tensor>>{Gamma340, rdm2_, rdm3_};
  auto task81 = make_shared<Task81>(tensor81, pindex);
  return make_shared<FutureTensor>(*Gamma340, task81);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma341_() {
  vector<IndexRange> Gamma341_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma341 = make_shared<Tensor>(Gamma341_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor82 = vector<shared_ptr<Tensor>>{Gamma341, rdm2_, rdm3_, rdm4_};
  auto task82 = make_shared<Task82>(tensor82, pindex);
  return make_shared<FutureTensor>(*Gamma341, task82);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma342_() {
  vector<IndexRange> Gamma342_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma342 = make_shared<Tensor>(Gamma342_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor83 = vector<shared_ptr<Tensor>>{Gamma342, rdm2_, rdm3_, rdm4_};
  auto task83 = make_shared<Task83>(tensor83, pindex);
  return make_shared<FutureTensor>(*Gamma342, task83);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma345_() {
  vector<IndexRange> Gamma345_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma345 = make_shared<Tensor>(Gamma345_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor84 = vector<shared_ptr<Tensor>>{Gamma345, rdm2_, rdm3_, rdm4_};
  auto task84 = make_shared<Task84>(tensor84, pindex);
  return make_shared<FutureTensor>(*Gamma345, task84);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma346_() {
  vector<IndexRange> Gamma346_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma346 = make_shared<Tensor>(Gamma346_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor85 = vector<shared_ptr<Tensor>>{Gamma346, rdm2_, rdm3_, rdm4_};
  auto task85 = make_shared<Task85>(tensor85, pindex);
  return make_shared<FutureTensor>(*Gamma346, task85);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma349_() {
  vector<IndexRange> Gamma349_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma349 = make_shared<Tensor>(Gamma349_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor86 = vector<shared_ptr<Tensor>>{Gamma349, rdm2_, rdm3_, rdm4_};
  auto task86 = make_shared<Task86>(tensor86, pindex);
  return make_shared<FutureTensor>(*Gamma349, task86);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma350_() {
  vector<IndexRange> Gamma350_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma350 = make_shared<Tensor>(Gamma350_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor87 = vector<shared_ptr<Tensor>>{Gamma350, rdm2_, rdm3_, rdm4_};
  auto task87 = make_shared<Task87>(tensor87, pindex);
  return make_shared<FutureTensor>(*Gamma350, task87);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma351_() {
  vector<IndexRange> Gamma351_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma351 = make_shared<Tensor>(Gamma351_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor88 = vector<shared_ptr<Tensor>>{Gamma351, rdm2_, rdm3_, rdm4_};
  auto task88 = make_shared<Task88>(tensor88, pindex);
  return make_shared<FutureTensor>(*Gamma351, task88);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma359_() {
  vector<IndexRange> Gamma359_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma359 = make_shared<Tensor>(Gamma359_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor89 = vector<shared_ptr<Tensor>>{Gamma359, rdm2_, rdm3_};
  auto task89 = make_shared<Task89>(tensor89, pindex);
  return make_shared<FutureTensor>(*Gamma359, task89);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma366_() {
  vector<IndexRange> Gamma366_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma366 = make_shared<Tensor>(Gamma366_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor90 = vector<shared_ptr<Tensor>>{Gamma366, rdm3_, rdm4_};
  auto task90 = make_shared<Task90>(tensor90, pindex);
  return make_shared<FutureTensor>(*Gamma366, task90);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma556_() {
  vector<IndexRange> Gamma556_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma556 = make_shared<Tensor>(Gamma556_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor91 = vector<shared_ptr<Tensor>>{Gamma556, rdm2_, rdm3_, rdm4_, h1_};
  auto task91 = make_shared<Task91>(tensor91, pindex);
  return make_shared<FutureTensor>(*Gamma556, task91);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma557_() {
  vector<IndexRange> Gamma557_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma557 = make_shared<Tensor>(Gamma557_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor92 = vector<shared_ptr<Tensor>>{Gamma557, rdm2_, rdm3_, rdm4_, v2_};
  auto task92 = make_shared<Task92>(tensor92, pindex);
  return make_shared<FutureTensor>(*Gamma557, task92);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma471_() {
  vector<IndexRange> Gamma471_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma471 = make_shared<Tensor>(Gamma471_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor93 = vector<shared_ptr<Tensor>>{Gamma471, rdm2_, rdm3_};
  auto task93 = make_shared<Task93>(tensor93, pindex);
  return make_shared<FutureTensor>(*Gamma471, task93);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma503_() {
  vector<IndexRange> Gamma503_index = {active_, active_, active_, active_};
  auto Gamma503 = make_shared<Tensor>(Gamma503_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor94 = vector<shared_ptr<Tensor>>{Gamma503, rdm2_};
  auto task94 = make_shared<Task94>(tensor94, pindex);
  return make_shared<FutureTensor>(*Gamma503, task94);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma526_() {
  vector<IndexRange> Gamma526_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma526 = make_shared<Tensor>(Gamma526_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor95 = vector<shared_ptr<Tensor>>{Gamma526, rdm3_};
  auto task95 = make_shared<Task95>(tensor95, pindex);
  return make_shared<FutureTensor>(*Gamma526, task95);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma562_() {
  vector<IndexRange> Gamma562_index = {active_, active_};
  auto Gamma562 = make_shared<Tensor>(Gamma562_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor96 = vector<shared_ptr<Tensor>>{Gamma562, rdm1_, rdm2_, h1_};
  auto task96 = make_shared<Task96>(tensor96, pindex);
  return make_shared<FutureTensor>(*Gamma562, task96);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma564_() {
  vector<IndexRange> Gamma564_index = {active_, active_};
  auto Gamma564 = make_shared<Tensor>(Gamma564_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor97 = vector<shared_ptr<Tensor>>{Gamma564, rdm1_, rdm2_, rdm3_, v2_};
  auto task97 = make_shared<Task97>(tensor97, pindex);
  return make_shared<FutureTensor>(*Gamma564, task97);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma531_() {
  vector<IndexRange> Gamma531_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma531 = make_shared<Tensor>(Gamma531_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor98 = vector<shared_ptr<Tensor>>{Gamma531, rdm2_, rdm3_};
  auto task98 = make_shared<Task98>(tensor98, pindex);
  return make_shared<FutureTensor>(*Gamma531, task98);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma532_() {
  vector<IndexRange> Gamma532_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma532 = make_shared<Tensor>(Gamma532_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor99 = vector<shared_ptr<Tensor>>{Gamma532, rdm2_, rdm3_};
  auto task99 = make_shared<Task99>(tensor99, pindex);
  return make_shared<FutureTensor>(*Gamma532, task99);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma533_() {
  vector<IndexRange> Gamma533_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma533 = make_shared<Tensor>(Gamma533_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor100 = vector<shared_ptr<Tensor>>{Gamma533, rdm3_, rdm4_};
  auto task100 = make_shared<Task100>(tensor100, pindex);
  return make_shared<FutureTensor>(*Gamma533, task100);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma545_() {
  vector<IndexRange> Gamma545_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma545 = make_shared<Tensor>(Gamma545_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor101 = vector<shared_ptr<Tensor>>{Gamma545, rdm3_};
  auto task101 = make_shared<Task101>(tensor101, pindex);
  return make_shared<FutureTensor>(*Gamma545, task101);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma548_() {
  vector<IndexRange> Gamma548_index = {active_, active_, active_, active_};
  auto Gamma548 = make_shared<Tensor>(Gamma548_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor102 = vector<shared_ptr<Tensor>>{Gamma548, rdm0_, rdm1_, rdm2_, rdm3_, h1_};
  auto task102 = make_shared<Task102>(tensor102, pindex);
  return make_shared<FutureTensor>(*Gamma548, task102);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma549_() {
  vector<IndexRange> Gamma549_index = {active_, active_, active_, active_};
  auto Gamma549 = make_shared<Tensor>(Gamma549_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor103 = vector<shared_ptr<Tensor>>{Gamma549, rdm0_, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task103 = make_shared<Task103>(tensor103, pindex);
  return make_shared<FutureTensor>(*Gamma549, task103);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma558_() {
  vector<IndexRange> Gamma558_index;
  auto Gamma558 = make_shared<Tensor>(Gamma558_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor104 = vector<shared_ptr<Tensor>>{Gamma558, rdm1_, h1_};
  auto task104 = make_shared<Task104>(tensor104, pindex);
  return make_shared<FutureTensor>(*Gamma558, task104);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma560_() {
  vector<IndexRange> Gamma560_index;
  auto Gamma560 = make_shared<Tensor>(Gamma560_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor105 = vector<shared_ptr<Tensor>>{Gamma560, rdm1_, rdm2_, v2_};
  auto task105 = make_shared<Task105>(tensor105, pindex);
  return make_shared<FutureTensor>(*Gamma560, task105);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma566_() {
  vector<IndexRange> Gamma566_index = {active_, active_, active_, active_};
  auto Gamma566 = make_shared<Tensor>(Gamma566_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor106 = vector<shared_ptr<Tensor>>{Gamma566, rdm2_, rdm3_, h1_};
  auto task106 = make_shared<Task106>(tensor106, pindex);
  return make_shared<FutureTensor>(*Gamma566, task106);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma567_() {
  vector<IndexRange> Gamma567_index = {active_, active_, active_, active_};
  auto Gamma567 = make_shared<Tensor>(Gamma567_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor107 = vector<shared_ptr<Tensor>>{Gamma567, rdm2_, rdm3_, rdm4_, v2_};
  auto task107 = make_shared<Task107>(tensor107, pindex);
  return make_shared<FutureTensor>(*Gamma567, task107);
}

#endif
