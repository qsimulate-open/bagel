//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASA_gamma.cc
// Copyright (C) 2014 Shiozaki group
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
#include <src/smith/casa/CASA.h>
#include <src/smith/casa/CASA_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASA;

shared_ptr<FutureTensor> CASA::CASA::Gamma0_() {
  vector<IndexRange> Gamma0_index = {active_, active_, active_, active_};
  auto Gamma0 = make_shared<Tensor>(Gamma0_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor0 = vector<shared_ptr<Tensor>>{Gamma0, rdm0_, rdm1_, rdm2_};
  auto task0 = make_shared<Task0>(tensor0, pindex);
  return make_shared<FutureTensor>(*Gamma0, task0);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma1_() {
  vector<IndexRange> Gamma1_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma1 = make_shared<Tensor>(Gamma1_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor1 = vector<shared_ptr<Tensor>>{Gamma1, rdm1_, rdm2_, rdm3_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  return make_shared<FutureTensor>(*Gamma1, task1);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma2_() {
  vector<IndexRange> Gamma2_index = {active_, active_, active_, active_};
  auto Gamma2 = make_shared<Tensor>(Gamma2_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor2 = vector<shared_ptr<Tensor>>{Gamma2, rdm0_, rdm1_, rdm2_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  return make_shared<FutureTensor>(*Gamma2, task2);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma3_() {
  vector<IndexRange> Gamma3_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma3 = make_shared<Tensor>(Gamma3_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor3 = vector<shared_ptr<Tensor>>{Gamma3, rdm1_, rdm2_, rdm3_};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  return make_shared<FutureTensor>(*Gamma3, task3);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma4_() {
  vector<IndexRange> Gamma4_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma4 = make_shared<Tensor>(Gamma4_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor4 = vector<shared_ptr<Tensor>>{Gamma4, rdm1_, rdm2_, rdm3_};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  return make_shared<FutureTensor>(*Gamma4, task4);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma5_() {
  vector<IndexRange> Gamma5_index = {active_, active_, active_, active_};
  auto Gamma5 = make_shared<Tensor>(Gamma5_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor5 = vector<shared_ptr<Tensor>>{Gamma5, rdm1_, rdm2_};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  return make_shared<FutureTensor>(*Gamma5, task5);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma7_() {
  vector<IndexRange> Gamma7_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma7 = make_shared<Tensor>(Gamma7_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor6 = vector<shared_ptr<Tensor>>{Gamma7, rdm1_, rdm2_, rdm3_};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  return make_shared<FutureTensor>(*Gamma7, task6);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma82_() {
  vector<IndexRange> Gamma82_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma82 = make_shared<Tensor>(Gamma82_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor7 = vector<shared_ptr<Tensor>>{Gamma82, rdm1_, rdm2_, rdm3_, rdm4_, h1_};
  auto task7 = make_shared<Task7>(tensor7, pindex);
  return make_shared<FutureTensor>(*Gamma82, task7);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma83_() {
  vector<IndexRange> Gamma83_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma83 = make_shared<Tensor>(Gamma83_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor8 = vector<shared_ptr<Tensor>>{Gamma83, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  return make_shared<FutureTensor>(*Gamma83, task8);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma9_() {
  vector<IndexRange> Gamma9_index = {active_, active_};
  auto Gamma9 = make_shared<Tensor>(Gamma9_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor9 = vector<shared_ptr<Tensor>>{Gamma9, rdm0_, rdm1_};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  return make_shared<FutureTensor>(*Gamma9, task9);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma13_() {
  vector<IndexRange> Gamma13_index = {active_, active_, active_, active_};
  auto Gamma13 = make_shared<Tensor>(Gamma13_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor10 = vector<shared_ptr<Tensor>>{Gamma13, rdm1_, rdm2_};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  return make_shared<FutureTensor>(*Gamma13, task10);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma15_() {
  vector<IndexRange> Gamma15_index = {active_, active_, active_, active_};
  auto Gamma15 = make_shared<Tensor>(Gamma15_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor11 = vector<shared_ptr<Tensor>>{Gamma15, rdm1_, rdm2_};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  return make_shared<FutureTensor>(*Gamma15, task11);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma84_() {
  vector<IndexRange> Gamma84_index = {active_, active_};
  auto Gamma84 = make_shared<Tensor>(Gamma84_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor12 = vector<shared_ptr<Tensor>>{Gamma84, rdm0_, rdm1_, rdm2_, h1_};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  return make_shared<FutureTensor>(*Gamma84, task12);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma86_() {
  vector<IndexRange> Gamma86_index = {active_, active_};
  auto Gamma86 = make_shared<Tensor>(Gamma86_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor13 = vector<shared_ptr<Tensor>>{Gamma86, rdm0_, rdm1_, rdm2_, rdm3_, v2_};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  return make_shared<FutureTensor>(*Gamma86, task13);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma24_() {
  vector<IndexRange> Gamma24_index = {active_, active_, active_, active_};
  auto Gamma24 = make_shared<Tensor>(Gamma24_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor14 = vector<shared_ptr<Tensor>>{Gamma24, rdm1_, rdm2_};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  return make_shared<FutureTensor>(*Gamma24, task14);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma25_() {
  vector<IndexRange> Gamma25_index = {active_, active_, active_, active_};
  auto Gamma25 = make_shared<Tensor>(Gamma25_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor15 = vector<shared_ptr<Tensor>>{Gamma25, rdm1_, rdm2_};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  return make_shared<FutureTensor>(*Gamma25, task15);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma26_() {
  vector<IndexRange> Gamma26_index = {active_, active_, active_, active_};
  auto Gamma26 = make_shared<Tensor>(Gamma26_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor16 = vector<shared_ptr<Tensor>>{Gamma26, rdm1_, rdm2_};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  return make_shared<FutureTensor>(*Gamma26, task16);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma28_() {
  vector<IndexRange> Gamma28_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma28 = make_shared<Tensor>(Gamma28_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor17 = vector<shared_ptr<Tensor>>{Gamma28, rdm2_, rdm3_};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  return make_shared<FutureTensor>(*Gamma28, task17);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma29_() {
  vector<IndexRange> Gamma29_index = {active_, active_};
  auto Gamma29 = make_shared<Tensor>(Gamma29_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor18 = vector<shared_ptr<Tensor>>{Gamma29, rdm1_};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  return make_shared<FutureTensor>(*Gamma29, task18);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma33_() {
  vector<IndexRange> Gamma33_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma33 = make_shared<Tensor>(Gamma33_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor19 = vector<shared_ptr<Tensor>>{Gamma33, rdm1_, rdm2_, rdm3_};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  return make_shared<FutureTensor>(*Gamma33, task19);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma100_() {
  vector<IndexRange> Gamma100_index = {active_, active_, active_, active_};
  auto Gamma100 = make_shared<Tensor>(Gamma100_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor20 = vector<shared_ptr<Tensor>>{Gamma100, rdm1_, rdm2_, rdm3_, h1_};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  return make_shared<FutureTensor>(*Gamma100, task20);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma101_() {
  vector<IndexRange> Gamma101_index = {active_, active_, active_, active_};
  auto Gamma101 = make_shared<Tensor>(Gamma101_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor21 = vector<shared_ptr<Tensor>>{Gamma101, rdm1_, rdm2_, rdm3_, h1_};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  return make_shared<FutureTensor>(*Gamma101, task21);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma104_() {
  vector<IndexRange> Gamma104_index = {active_, active_, active_, active_};
  auto Gamma104 = make_shared<Tensor>(Gamma104_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor22 = vector<shared_ptr<Tensor>>{Gamma104, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  return make_shared<FutureTensor>(*Gamma104, task22);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma105_() {
  vector<IndexRange> Gamma105_index = {active_, active_, active_, active_};
  auto Gamma105 = make_shared<Tensor>(Gamma105_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor23 = vector<shared_ptr<Tensor>>{Gamma105, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  return make_shared<FutureTensor>(*Gamma105, task23);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma40_() {
  vector<IndexRange> Gamma40_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma40 = make_shared<Tensor>(Gamma40_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor24 = vector<shared_ptr<Tensor>>{Gamma40, rdm2_, rdm3_};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  return make_shared<FutureTensor>(*Gamma40, task24);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma48_() {
  vector<IndexRange> Gamma48_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma48 = make_shared<Tensor>(Gamma48_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor25 = vector<shared_ptr<Tensor>>{Gamma48, rdm2_, rdm3_};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  return make_shared<FutureTensor>(*Gamma48, task25);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma49_() {
  vector<IndexRange> Gamma49_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma49 = make_shared<Tensor>(Gamma49_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor26 = vector<shared_ptr<Tensor>>{Gamma49, rdm2_, rdm3_};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  return make_shared<FutureTensor>(*Gamma49, task26);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma50_() {
  vector<IndexRange> Gamma50_index = {active_, active_, active_, active_};
  auto Gamma50 = make_shared<Tensor>(Gamma50_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor27 = vector<shared_ptr<Tensor>>{Gamma50, rdm2_};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  return make_shared<FutureTensor>(*Gamma50, task27);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma52_() {
  vector<IndexRange> Gamma52_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma52 = make_shared<Tensor>(Gamma52_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor28 = vector<shared_ptr<Tensor>>{Gamma52, rdm2_, rdm3_};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  return make_shared<FutureTensor>(*Gamma52, task28);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma88_() {
  vector<IndexRange> Gamma88_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma88 = make_shared<Tensor>(Gamma88_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor29 = vector<shared_ptr<Tensor>>{Gamma88, rdm2_, rdm3_, rdm4_, h1_};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  return make_shared<FutureTensor>(*Gamma88, task29);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma89_() {
  vector<IndexRange> Gamma89_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma89 = make_shared<Tensor>(Gamma89_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor30 = vector<shared_ptr<Tensor>>{Gamma89, rdm2_, rdm3_, rdm4_, v2_};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  return make_shared<FutureTensor>(*Gamma89, task30);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma94_() {
  vector<IndexRange> Gamma94_index = {active_, active_};
  auto Gamma94 = make_shared<Tensor>(Gamma94_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor31 = vector<shared_ptr<Tensor>>{Gamma94, rdm1_, rdm2_, h1_};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  return make_shared<FutureTensor>(*Gamma94, task31);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma96_() {
  vector<IndexRange> Gamma96_index = {active_, active_};
  auto Gamma96 = make_shared<Tensor>(Gamma96_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor32 = vector<shared_ptr<Tensor>>{Gamma96, rdm1_, rdm2_, rdm3_, v2_};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  return make_shared<FutureTensor>(*Gamma96, task32);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma80_() {
  vector<IndexRange> Gamma80_index = {active_, active_, active_, active_};
  auto Gamma80 = make_shared<Tensor>(Gamma80_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor33 = vector<shared_ptr<Tensor>>{Gamma80, rdm0_, rdm1_, rdm2_, rdm3_, h1_};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  return make_shared<FutureTensor>(*Gamma80, task33);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma81_() {
  vector<IndexRange> Gamma81_index = {active_, active_, active_, active_};
  auto Gamma81 = make_shared<Tensor>(Gamma81_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor34 = vector<shared_ptr<Tensor>>{Gamma81, rdm0_, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  return make_shared<FutureTensor>(*Gamma81, task34);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma90_() {
  vector<IndexRange> Gamma90_index;
  auto Gamma90 = make_shared<Tensor>(Gamma90_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor35 = vector<shared_ptr<Tensor>>{Gamma90, rdm1_, h1_};
  auto task35 = make_shared<Task35>(tensor35, pindex);
  return make_shared<FutureTensor>(*Gamma90, task35);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma92_() {
  vector<IndexRange> Gamma92_index;
  auto Gamma92 = make_shared<Tensor>(Gamma92_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor36 = vector<shared_ptr<Tensor>>{Gamma92, rdm1_, rdm2_, v2_};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  return make_shared<FutureTensor>(*Gamma92, task36);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma98_() {
  vector<IndexRange> Gamma98_index = {active_, active_, active_, active_};
  auto Gamma98 = make_shared<Tensor>(Gamma98_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor37 = vector<shared_ptr<Tensor>>{Gamma98, rdm2_, rdm3_, h1_};
  auto task37 = make_shared<Task37>(tensor37, pindex);
  return make_shared<FutureTensor>(*Gamma98, task37);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma99_() {
  vector<IndexRange> Gamma99_index = {active_, active_, active_, active_};
  auto Gamma99 = make_shared<Tensor>(Gamma99_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor38 = vector<shared_ptr<Tensor>>{Gamma99, rdm2_, rdm3_, rdm4_, v2_};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  return make_shared<FutureTensor>(*Gamma99, task38);
}

shared_ptr<FutureTensor> CASA::CASA::Gamma109_() {
  vector<IndexRange> Gamma109_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma109 = make_shared<Tensor>(Gamma109_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor39 = vector<shared_ptr<Tensor>>{Gamma109, rdm1_, rdm2_, rdm3_};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  return make_shared<FutureTensor>(*Gamma109, task39);
}

#endif
