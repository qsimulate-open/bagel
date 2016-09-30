//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelCASA_gamma.cc
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
#include <src/smith/relcasa/RelCASA.h>
#include <src/smith/relcasa/RelCASA_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelCASA;

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma0_() {
  vector<IndexRange> Gamma0_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma0 = make_shared<Tensor>(Gamma0_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor0 = vector<shared_ptr<Tensor>>{Gamma0, rdm1_, rdm2_, rdm3_};
  auto task0 = make_shared<Task0>(tensor0, pindex);
  return make_shared<FutureTensor>(*Gamma0, task0);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma1_() {
  vector<IndexRange> Gamma1_index = {active_, active_, active_, active_};
  auto Gamma1 = make_shared<Tensor>(Gamma1_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor1 = vector<shared_ptr<Tensor>>{Gamma1, rdm0_, rdm1_, rdm2_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  return make_shared<FutureTensor>(*Gamma1, task1);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma2_() {
  vector<IndexRange> Gamma2_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma2 = make_shared<Tensor>(Gamma2_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor2 = vector<shared_ptr<Tensor>>{Gamma2, rdm1_, rdm2_, rdm3_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  return make_shared<FutureTensor>(*Gamma2, task2);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma3_() {
  vector<IndexRange> Gamma3_index = {active_, active_, active_, active_};
  auto Gamma3 = make_shared<Tensor>(Gamma3_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor3 = vector<shared_ptr<Tensor>>{Gamma3, rdm1_, rdm2_};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  return make_shared<FutureTensor>(*Gamma3, task3);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma5_() {
  vector<IndexRange> Gamma5_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma5 = make_shared<Tensor>(Gamma5_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor4 = vector<shared_ptr<Tensor>>{Gamma5, rdm1_, rdm2_, rdm3_};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  return make_shared<FutureTensor>(*Gamma5, task4);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma6_() {
  vector<IndexRange> Gamma6_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma6 = make_shared<Tensor>(Gamma6_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor5 = vector<shared_ptr<Tensor>>{Gamma6, rdm1_, rdm2_, rdm3_};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  return make_shared<FutureTensor>(*Gamma6, task5);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma70_() {
  vector<IndexRange> Gamma70_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma70 = make_shared<Tensor>(Gamma70_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor6 = vector<shared_ptr<Tensor>>{Gamma70, rdm1_, rdm2_, rdm3_, rdm4_, h1_};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  return make_shared<FutureTensor>(*Gamma70, task6);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma71_() {
  vector<IndexRange> Gamma71_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma71 = make_shared<Tensor>(Gamma71_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor7 = vector<shared_ptr<Tensor>>{Gamma71, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task7 = make_shared<Task7>(tensor7, pindex);
  return make_shared<FutureTensor>(*Gamma71, task7);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma7_() {
  vector<IndexRange> Gamma7_index = {active_, active_, active_, active_};
  auto Gamma7 = make_shared<Tensor>(Gamma7_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor8 = vector<shared_ptr<Tensor>>{Gamma7, rdm1_, rdm2_};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  return make_shared<FutureTensor>(*Gamma7, task8);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma9_() {
  vector<IndexRange> Gamma9_index = {active_, active_, active_, active_};
  auto Gamma9 = make_shared<Tensor>(Gamma9_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor9 = vector<shared_ptr<Tensor>>{Gamma9, rdm1_, rdm2_};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  return make_shared<FutureTensor>(*Gamma9, task9);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma13_() {
  vector<IndexRange> Gamma13_index = {active_, active_};
  auto Gamma13 = make_shared<Tensor>(Gamma13_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor10 = vector<shared_ptr<Tensor>>{Gamma13, rdm0_, rdm1_};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  return make_shared<FutureTensor>(*Gamma13, task10);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma72_() {
  vector<IndexRange> Gamma72_index = {active_, active_};
  auto Gamma72 = make_shared<Tensor>(Gamma72_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor11 = vector<shared_ptr<Tensor>>{Gamma72, rdm0_, rdm1_, rdm2_, h1_};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  return make_shared<FutureTensor>(*Gamma72, task11);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma74_() {
  vector<IndexRange> Gamma74_index = {active_, active_};
  auto Gamma74 = make_shared<Tensor>(Gamma74_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor12 = vector<shared_ptr<Tensor>>{Gamma74, rdm0_, rdm1_, rdm2_, rdm3_, v2_};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  return make_shared<FutureTensor>(*Gamma74, task12);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma18_() {
  vector<IndexRange> Gamma18_index = {active_, active_, active_, active_};
  auto Gamma18 = make_shared<Tensor>(Gamma18_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor13 = vector<shared_ptr<Tensor>>{Gamma18, rdm1_, rdm2_};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  return make_shared<FutureTensor>(*Gamma18, task13);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma20_() {
  vector<IndexRange> Gamma20_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma20 = make_shared<Tensor>(Gamma20_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor14 = vector<shared_ptr<Tensor>>{Gamma20, rdm2_, rdm3_};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  return make_shared<FutureTensor>(*Gamma20, task14);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma21_() {
  vector<IndexRange> Gamma21_index = {active_, active_};
  auto Gamma21 = make_shared<Tensor>(Gamma21_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor15 = vector<shared_ptr<Tensor>>{Gamma21, rdm1_};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  return make_shared<FutureTensor>(*Gamma21, task15);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma23_() {
  vector<IndexRange> Gamma23_index = {active_, active_, active_, active_};
  auto Gamma23 = make_shared<Tensor>(Gamma23_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor16 = vector<shared_ptr<Tensor>>{Gamma23, rdm1_, rdm2_};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  return make_shared<FutureTensor>(*Gamma23, task16);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma24_() {
  vector<IndexRange> Gamma24_index = {active_, active_, active_, active_};
  auto Gamma24 = make_shared<Tensor>(Gamma24_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor17 = vector<shared_ptr<Tensor>>{Gamma24, rdm1_, rdm2_};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  return make_shared<FutureTensor>(*Gamma24, task17);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma25_() {
  vector<IndexRange> Gamma25_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma25 = make_shared<Tensor>(Gamma25_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor18 = vector<shared_ptr<Tensor>>{Gamma25, rdm1_, rdm2_, rdm3_};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  return make_shared<FutureTensor>(*Gamma25, task18);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma88_() {
  vector<IndexRange> Gamma88_index = {active_, active_, active_, active_};
  auto Gamma88 = make_shared<Tensor>(Gamma88_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor19 = vector<shared_ptr<Tensor>>{Gamma88, rdm1_, rdm2_, rdm3_, h1_};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  return make_shared<FutureTensor>(*Gamma88, task19);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma89_() {
  vector<IndexRange> Gamma89_index = {active_, active_, active_, active_};
  auto Gamma89 = make_shared<Tensor>(Gamma89_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor20 = vector<shared_ptr<Tensor>>{Gamma89, rdm1_, rdm2_, rdm3_, h1_};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  return make_shared<FutureTensor>(*Gamma89, task20);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma92_() {
  vector<IndexRange> Gamma92_index = {active_, active_, active_, active_};
  auto Gamma92 = make_shared<Tensor>(Gamma92_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor21 = vector<shared_ptr<Tensor>>{Gamma92, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  return make_shared<FutureTensor>(*Gamma92, task21);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma93_() {
  vector<IndexRange> Gamma93_index = {active_, active_, active_, active_};
  auto Gamma93 = make_shared<Tensor>(Gamma93_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor22 = vector<shared_ptr<Tensor>>{Gamma93, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  return make_shared<FutureTensor>(*Gamma93, task22);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma30_() {
  vector<IndexRange> Gamma30_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma30 = make_shared<Tensor>(Gamma30_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor23 = vector<shared_ptr<Tensor>>{Gamma30, rdm2_, rdm3_};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  return make_shared<FutureTensor>(*Gamma30, task23);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma38_() {
  vector<IndexRange> Gamma38_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma38 = make_shared<Tensor>(Gamma38_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor24 = vector<shared_ptr<Tensor>>{Gamma38, rdm2_, rdm3_};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  return make_shared<FutureTensor>(*Gamma38, task24);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma39_() {
  vector<IndexRange> Gamma39_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma39 = make_shared<Tensor>(Gamma39_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor25 = vector<shared_ptr<Tensor>>{Gamma39, rdm2_, rdm3_};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  return make_shared<FutureTensor>(*Gamma39, task25);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma40_() {
  vector<IndexRange> Gamma40_index = {active_, active_, active_, active_};
  auto Gamma40 = make_shared<Tensor>(Gamma40_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor26 = vector<shared_ptr<Tensor>>{Gamma40, rdm2_};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  return make_shared<FutureTensor>(*Gamma40, task26);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma42_() {
  vector<IndexRange> Gamma42_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma42 = make_shared<Tensor>(Gamma42_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor27 = vector<shared_ptr<Tensor>>{Gamma42, rdm2_, rdm3_};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  return make_shared<FutureTensor>(*Gamma42, task27);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma76_() {
  vector<IndexRange> Gamma76_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma76 = make_shared<Tensor>(Gamma76_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor28 = vector<shared_ptr<Tensor>>{Gamma76, rdm2_, rdm3_, rdm4_, h1_};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  return make_shared<FutureTensor>(*Gamma76, task28);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma77_() {
  vector<IndexRange> Gamma77_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma77 = make_shared<Tensor>(Gamma77_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor29 = vector<shared_ptr<Tensor>>{Gamma77, rdm2_, rdm3_, rdm4_, v2_};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  return make_shared<FutureTensor>(*Gamma77, task29);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma82_() {
  vector<IndexRange> Gamma82_index = {active_, active_};
  auto Gamma82 = make_shared<Tensor>(Gamma82_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor30 = vector<shared_ptr<Tensor>>{Gamma82, rdm1_, rdm2_, h1_};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  return make_shared<FutureTensor>(*Gamma82, task30);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma84_() {
  vector<IndexRange> Gamma84_index = {active_, active_};
  auto Gamma84 = make_shared<Tensor>(Gamma84_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor31 = vector<shared_ptr<Tensor>>{Gamma84, rdm1_, rdm2_, rdm3_, v2_};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  return make_shared<FutureTensor>(*Gamma84, task31);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma68_() {
  vector<IndexRange> Gamma68_index = {active_, active_, active_, active_};
  auto Gamma68 = make_shared<Tensor>(Gamma68_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor32 = vector<shared_ptr<Tensor>>{Gamma68, rdm0_, rdm1_, rdm2_, rdm3_, h1_};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  return make_shared<FutureTensor>(*Gamma68, task32);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma69_() {
  vector<IndexRange> Gamma69_index = {active_, active_, active_, active_};
  auto Gamma69 = make_shared<Tensor>(Gamma69_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor33 = vector<shared_ptr<Tensor>>{Gamma69, rdm0_, rdm1_, rdm2_, rdm3_, rdm4_, v2_};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  return make_shared<FutureTensor>(*Gamma69, task33);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma78_() {
  vector<IndexRange> Gamma78_index;
  auto Gamma78 = make_shared<Tensor>(Gamma78_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor34 = vector<shared_ptr<Tensor>>{Gamma78, rdm1_, h1_};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  return make_shared<FutureTensor>(*Gamma78, task34);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma80_() {
  vector<IndexRange> Gamma80_index;
  auto Gamma80 = make_shared<Tensor>(Gamma80_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor35 = vector<shared_ptr<Tensor>>{Gamma80, rdm1_, rdm2_, v2_};
  auto task35 = make_shared<Task35>(tensor35, pindex);
  return make_shared<FutureTensor>(*Gamma80, task35);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma86_() {
  vector<IndexRange> Gamma86_index = {active_, active_, active_, active_};
  auto Gamma86 = make_shared<Tensor>(Gamma86_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor36 = vector<shared_ptr<Tensor>>{Gamma86, rdm2_, rdm3_, h1_};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  return make_shared<FutureTensor>(*Gamma86, task36);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma87_() {
  vector<IndexRange> Gamma87_index = {active_, active_, active_, active_};
  auto Gamma87 = make_shared<Tensor>(Gamma87_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor37 = vector<shared_ptr<Tensor>>{Gamma87, rdm2_, rdm3_, rdm4_, v2_};
  auto task37 = make_shared<Task37>(tensor37, pindex);
  return make_shared<FutureTensor>(*Gamma87, task37);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma96_() {
  vector<IndexRange> Gamma96_index = {active_, active_, active_, active_};
  auto Gamma96 = make_shared<Tensor>(Gamma96_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor38 = vector<shared_ptr<Tensor>>{Gamma96, rdm0_, rdm1_, rdm2_};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  return make_shared<FutureTensor>(*Gamma96, task38);
}

shared_ptr<FutureTensor> RelCASA::RelCASA::Gamma97_() {
  vector<IndexRange> Gamma97_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma97 = make_shared<Tensor>(Gamma97_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor39 = vector<shared_ptr<Tensor>>{Gamma97, rdm1_, rdm2_, rdm3_};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  return make_shared<FutureTensor>(*Gamma97, task39);
}

#endif
