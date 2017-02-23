//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_gamma.cc
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
#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma0_() {
  vector<IndexRange> Gamma0_index = {active_, active_, active_, active_};
  auto Gamma0 = make_shared<Tensor>(Gamma0_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor0 = vector<shared_ptr<Tensor>>{Gamma0, rdm0_, rdm1_, rdm2_, rdm3_, f1_};
  auto task0 = make_shared<Task0>(tensor0, pindex);
  return make_shared<FutureTensor>(*Gamma0, task0);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma92_() {
  vector<IndexRange> Gamma92_index = {active_, active_, active_, active_};
  auto Gamma92 = make_shared<Tensor>(Gamma92_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor1 = vector<shared_ptr<Tensor>>{Gamma92, rdm0_, rdm1_, rdm2_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  return make_shared<FutureTensor>(*Gamma92, task1);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma2_() {
  vector<IndexRange> Gamma2_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma2 = make_shared<Tensor>(Gamma2_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor2 = vector<shared_ptr<Tensor>>{Gamma2, rdm1_, rdm2_, rdm3_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  return make_shared<FutureTensor>(*Gamma2, task2);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma3_() {
  vector<IndexRange> Gamma3_index = {active_, active_, active_, active_};
  auto Gamma3 = make_shared<Tensor>(Gamma3_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor3 = vector<shared_ptr<Tensor>>{Gamma3, rdm0_, rdm1_, rdm2_};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  return make_shared<FutureTensor>(*Gamma3, task3);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma4_() {
  vector<IndexRange> Gamma4_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma4 = make_shared<Tensor>(Gamma4_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor4 = vector<shared_ptr<Tensor>>{Gamma4, rdm1_, rdm2_, rdm3_};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  return make_shared<FutureTensor>(*Gamma4, task4);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma5_() {
  vector<IndexRange> Gamma5_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma5 = make_shared<Tensor>(Gamma5_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor5 = vector<shared_ptr<Tensor>>{Gamma5, rdm1_, rdm2_, rdm3_, rdm4_, f1_};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  return make_shared<FutureTensor>(*Gamma5, task5);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma6_() {
  vector<IndexRange> Gamma6_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma6 = make_shared<Tensor>(Gamma6_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor6 = vector<shared_ptr<Tensor>>{Gamma6, rdm1_, rdm2_, rdm3_};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  return make_shared<FutureTensor>(*Gamma6, task6);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma7_() {
  vector<IndexRange> Gamma7_index = {active_, active_, active_, active_};
  auto Gamma7 = make_shared<Tensor>(Gamma7_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor7 = vector<shared_ptr<Tensor>>{Gamma7, rdm1_, rdm2_};
  auto task7 = make_shared<Task7>(tensor7, pindex);
  return make_shared<FutureTensor>(*Gamma7, task7);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma9_() {
  vector<IndexRange> Gamma9_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma9 = make_shared<Tensor>(Gamma9_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor8 = vector<shared_ptr<Tensor>>{Gamma9, rdm1_, rdm2_, rdm3_};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  return make_shared<FutureTensor>(*Gamma9, task8);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma12_() {
  vector<IndexRange> Gamma12_index = {active_, active_, active_, active_};
  auto Gamma12 = make_shared<Tensor>(Gamma12_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor9 = vector<shared_ptr<Tensor>>{Gamma12, rdm1_, rdm2_};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  return make_shared<FutureTensor>(*Gamma12, task9);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma14_() {
  vector<IndexRange> Gamma14_index = {active_, active_};
  auto Gamma14 = make_shared<Tensor>(Gamma14_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor10 = vector<shared_ptr<Tensor>>{Gamma14, rdm0_, rdm1_, rdm2_, f1_};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  return make_shared<FutureTensor>(*Gamma14, task10);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma16_() {
  vector<IndexRange> Gamma16_index = {active_, active_};
  auto Gamma16 = make_shared<Tensor>(Gamma16_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor11 = vector<shared_ptr<Tensor>>{Gamma16, rdm0_, rdm1_};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  return make_shared<FutureTensor>(*Gamma16, task11);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma22_() {
  vector<IndexRange> Gamma22_index = {active_, active_, active_, active_};
  auto Gamma22 = make_shared<Tensor>(Gamma22_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor12 = vector<shared_ptr<Tensor>>{Gamma22, rdm1_, rdm2_};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  return make_shared<FutureTensor>(*Gamma22, task12);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma28_() {
  vector<IndexRange> Gamma28_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma28 = make_shared<Tensor>(Gamma28_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor13 = vector<shared_ptr<Tensor>>{Gamma28, rdm1_, rdm2_, rdm3_};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  return make_shared<FutureTensor>(*Gamma28, task13);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma29_() {
  vector<IndexRange> Gamma29_index = {active_, active_, active_, active_};
  auto Gamma29 = make_shared<Tensor>(Gamma29_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor14 = vector<shared_ptr<Tensor>>{Gamma29, rdm1_, rdm2_};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  return make_shared<FutureTensor>(*Gamma29, task14);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma31_() {
  vector<IndexRange> Gamma31_index = {active_, active_, active_, active_};
  auto Gamma31 = make_shared<Tensor>(Gamma31_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor15 = vector<shared_ptr<Tensor>>{Gamma31, rdm1_, rdm2_, rdm3_, f1_};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  return make_shared<FutureTensor>(*Gamma31, task15);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma32_() {
  vector<IndexRange> Gamma32_index = {active_, active_, active_, active_};
  auto Gamma32 = make_shared<Tensor>(Gamma32_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor16 = vector<shared_ptr<Tensor>>{Gamma32, rdm1_, rdm2_};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  return make_shared<FutureTensor>(*Gamma32, task16);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma35_() {
  vector<IndexRange> Gamma35_index = {active_, active_, active_, active_};
  auto Gamma35 = make_shared<Tensor>(Gamma35_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor17 = vector<shared_ptr<Tensor>>{Gamma35, rdm1_, rdm2_};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  return make_shared<FutureTensor>(*Gamma35, task17);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma34_() {
  vector<IndexRange> Gamma34_index = {active_, active_, active_, active_};
  auto Gamma34 = make_shared<Tensor>(Gamma34_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor18 = vector<shared_ptr<Tensor>>{Gamma34, rdm1_, rdm2_, rdm3_, f1_};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  return make_shared<FutureTensor>(*Gamma34, task18);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma37_() {
  vector<IndexRange> Gamma37_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma37 = make_shared<Tensor>(Gamma37_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor19 = vector<shared_ptr<Tensor>>{Gamma37, rdm2_, rdm3_};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  return make_shared<FutureTensor>(*Gamma37, task19);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma38_() {
  vector<IndexRange> Gamma38_index = {active_, active_};
  auto Gamma38 = make_shared<Tensor>(Gamma38_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor20 = vector<shared_ptr<Tensor>>{Gamma38, rdm1_};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  return make_shared<FutureTensor>(*Gamma38, task20);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma51_() {
  vector<IndexRange> Gamma51_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma51 = make_shared<Tensor>(Gamma51_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor21 = vector<shared_ptr<Tensor>>{Gamma51, rdm2_, rdm3_};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  return make_shared<FutureTensor>(*Gamma51, task21);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma56_() {
  vector<IndexRange> Gamma56_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma56 = make_shared<Tensor>(Gamma56_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor22 = vector<shared_ptr<Tensor>>{Gamma56, rdm2_, rdm3_};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  return make_shared<FutureTensor>(*Gamma56, task22);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma57_() {
  vector<IndexRange> Gamma57_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma57 = make_shared<Tensor>(Gamma57_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor23 = vector<shared_ptr<Tensor>>{Gamma57, rdm2_, rdm3_};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  return make_shared<FutureTensor>(*Gamma57, task23);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma58_() {
  vector<IndexRange> Gamma58_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma58 = make_shared<Tensor>(Gamma58_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor24 = vector<shared_ptr<Tensor>>{Gamma58, rdm2_, rdm3_, rdm4_, f1_};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  return make_shared<FutureTensor>(*Gamma58, task24);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma59_() {
  vector<IndexRange> Gamma59_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma59 = make_shared<Tensor>(Gamma59_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor25 = vector<shared_ptr<Tensor>>{Gamma59, rdm2_, rdm3_};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  return make_shared<FutureTensor>(*Gamma59, task25);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma60_() {
  vector<IndexRange> Gamma60_index = {active_, active_, active_, active_};
  auto Gamma60 = make_shared<Tensor>(Gamma60_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor26 = vector<shared_ptr<Tensor>>{Gamma60, rdm2_};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  return make_shared<FutureTensor>(*Gamma60, task26);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma79_() {
  vector<IndexRange> Gamma79_index = {active_, active_};
  auto Gamma79 = make_shared<Tensor>(Gamma79_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor27 = vector<shared_ptr<Tensor>>{Gamma79, rdm2_, f1_};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  return make_shared<FutureTensor>(*Gamma79, task27);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma90_() {
  vector<IndexRange> Gamma90_index = {active_, active_, active_, active_};
  auto Gamma90 = make_shared<Tensor>(Gamma90_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor28 = vector<shared_ptr<Tensor>>{Gamma90, rdm3_, f1_};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  return make_shared<FutureTensor>(*Gamma90, task28);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma105_() {
  vector<IndexRange> Gamma105_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma105 = make_shared<Tensor>(Gamma105_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor29 = vector<shared_ptr<Tensor>>{Gamma105, rdm1_, rdm2_, rdm3_};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  return make_shared<FutureTensor>(*Gamma105, task29);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma138_() {
  vector<IndexRange> Gamma138_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma138 = make_shared<Tensor>(Gamma138_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor30 = vector<shared_ptr<Tensor>>{Gamma138, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  return make_shared<FutureTensor>(*Gamma138, task30);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma169_() {
  vector<IndexRange> Gamma169_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma169 = make_shared<Tensor>(Gamma169_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor31 = vector<shared_ptr<Tensor>>{Gamma169, rdm1_, rdm2_, rdm3_};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  return make_shared<FutureTensor>(*Gamma169, task31);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma172_() {
  vector<IndexRange> Gamma172_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma172 = make_shared<Tensor>(Gamma172_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor32 = vector<shared_ptr<Tensor>>{Gamma172, rdm1_, rdm2_, rdm3_};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  return make_shared<FutureTensor>(*Gamma172, task32);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma230_() {
  vector<IndexRange> Gamma230_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma230 = make_shared<Tensor>(Gamma230_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor33 = vector<shared_ptr<Tensor>>{Gamma230, rdm3_};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  return make_shared<FutureTensor>(*Gamma230, task33);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma143_() {
  vector<IndexRange> Gamma143_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma143 = make_shared<Tensor>(Gamma143_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor34 = vector<shared_ptr<Tensor>>{Gamma143, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  return make_shared<FutureTensor>(*Gamma143, task34);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma196_() {
  vector<IndexRange> Gamma196_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma196 = make_shared<Tensor>(Gamma196_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor35 = vector<shared_ptr<Tensor>>{Gamma196, rdm2_, rdm3_, rdm4_};
  auto task35 = make_shared<Task35>(tensor35, pindex);
  return make_shared<FutureTensor>(*Gamma196, task35);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma152_() {
  vector<IndexRange> Gamma152_index = {active_, active_, active_, active_};
  auto Gamma152 = make_shared<Tensor>(Gamma152_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor36 = vector<shared_ptr<Tensor>>{Gamma152, rdm0_, rdm1_, rdm2_};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  return make_shared<FutureTensor>(*Gamma152, task36);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma248_() {
  vector<IndexRange> Gamma248_index = {ci_, active_, active_, active_, active_};
  auto Gamma248 = make_shared<Tensor>(Gamma248_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor37 = vector<shared_ptr<Tensor>>{Gamma248, rdm0deriv_, rdm1deriv_, rdm2deriv_, rdm3deriv_, f1_};
  auto task37 = make_shared<Task37>(tensor37, cindex);
  return make_shared<FutureTensor>(*Gamma248, task37);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma249_() {
  vector<IndexRange> Gamma249_index = {ci_, active_, active_, active_, active_};
  auto Gamma249 = make_shared<Tensor>(Gamma249_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor38 = vector<shared_ptr<Tensor>>{Gamma249, rdm0deriv_, rdm1deriv_, rdm2deriv_};
  auto task38 = make_shared<Task38>(tensor38, cindex);
  return make_shared<FutureTensor>(*Gamma249, task38);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma250_() {
  vector<IndexRange> Gamma250_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma250 = make_shared<Tensor>(Gamma250_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor39 = vector<shared_ptr<Tensor>>{Gamma250, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task39 = make_shared<Task39>(tensor39, cindex);
  return make_shared<FutureTensor>(*Gamma250, task39);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma251_() {
  vector<IndexRange> Gamma251_index = {ci_, active_, active_, active_, active_};
  auto Gamma251 = make_shared<Tensor>(Gamma251_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor40 = vector<shared_ptr<Tensor>>{Gamma251, rdm0deriv_, rdm1deriv_, rdm2deriv_};
  auto task40 = make_shared<Task40>(tensor40, cindex);
  return make_shared<FutureTensor>(*Gamma251, task40);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma252_() {
  vector<IndexRange> Gamma252_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma252 = make_shared<Tensor>(Gamma252_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor41 = vector<shared_ptr<Tensor>>{Gamma252, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task41 = make_shared<Task41>(tensor41, cindex);
  return make_shared<FutureTensor>(*Gamma252, task41);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma253_() {
  vector<IndexRange> Gamma253_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma253 = make_shared<Tensor>(Gamma253_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor42 = vector<shared_ptr<Tensor>>{Gamma253, rdm1deriv_, rdm2deriv_, rdm3deriv_, rdm4deriv_, f1_};
  auto task42 = make_shared<Task42>(tensor42, cindex);
  return make_shared<FutureTensor>(*Gamma253, task42);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma254_() {
  vector<IndexRange> Gamma254_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma254 = make_shared<Tensor>(Gamma254_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor43 = vector<shared_ptr<Tensor>>{Gamma254, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task43 = make_shared<Task43>(tensor43, cindex);
  return make_shared<FutureTensor>(*Gamma254, task43);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma255_() {
  vector<IndexRange> Gamma255_index = {ci_, active_, active_, active_, active_};
  auto Gamma255 = make_shared<Tensor>(Gamma255_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor44 = vector<shared_ptr<Tensor>>{Gamma255, rdm1deriv_, rdm2deriv_};
  auto task44 = make_shared<Task44>(tensor44, cindex);
  return make_shared<FutureTensor>(*Gamma255, task44);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma257_() {
  vector<IndexRange> Gamma257_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma257 = make_shared<Tensor>(Gamma257_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor45 = vector<shared_ptr<Tensor>>{Gamma257, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task45 = make_shared<Task45>(tensor45, cindex);
  return make_shared<FutureTensor>(*Gamma257, task45);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma260_() {
  vector<IndexRange> Gamma260_index = {ci_, active_, active_, active_, active_};
  auto Gamma260 = make_shared<Tensor>(Gamma260_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor46 = vector<shared_ptr<Tensor>>{Gamma260, rdm1deriv_, rdm2deriv_};
  auto task46 = make_shared<Task46>(tensor46, cindex);
  return make_shared<FutureTensor>(*Gamma260, task46);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma262_() {
  vector<IndexRange> Gamma262_index = {ci_, active_, active_};
  auto Gamma262 = make_shared<Tensor>(Gamma262_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor47 = vector<shared_ptr<Tensor>>{Gamma262, rdm0deriv_, rdm1deriv_, rdm2deriv_, f1_};
  auto task47 = make_shared<Task47>(tensor47, cindex);
  return make_shared<FutureTensor>(*Gamma262, task47);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma264_() {
  vector<IndexRange> Gamma264_index = {ci_, active_, active_};
  auto Gamma264 = make_shared<Tensor>(Gamma264_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor48 = vector<shared_ptr<Tensor>>{Gamma264, rdm0deriv_, rdm1deriv_};
  auto task48 = make_shared<Task48>(tensor48, cindex);
  return make_shared<FutureTensor>(*Gamma264, task48);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma270_() {
  vector<IndexRange> Gamma270_index = {ci_, active_, active_, active_, active_};
  auto Gamma270 = make_shared<Tensor>(Gamma270_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor49 = vector<shared_ptr<Tensor>>{Gamma270, rdm1deriv_, rdm2deriv_};
  auto task49 = make_shared<Task49>(tensor49, cindex);
  return make_shared<FutureTensor>(*Gamma270, task49);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma276_() {
  vector<IndexRange> Gamma276_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma276 = make_shared<Tensor>(Gamma276_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor50 = vector<shared_ptr<Tensor>>{Gamma276, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task50 = make_shared<Task50>(tensor50, cindex);
  return make_shared<FutureTensor>(*Gamma276, task50);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma277_() {
  vector<IndexRange> Gamma277_index = {ci_, active_, active_, active_, active_};
  auto Gamma277 = make_shared<Tensor>(Gamma277_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor51 = vector<shared_ptr<Tensor>>{Gamma277, rdm1deriv_, rdm2deriv_};
  auto task51 = make_shared<Task51>(tensor51, cindex);
  return make_shared<FutureTensor>(*Gamma277, task51);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma279_() {
  vector<IndexRange> Gamma279_index = {ci_, active_, active_, active_, active_};
  auto Gamma279 = make_shared<Tensor>(Gamma279_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor52 = vector<shared_ptr<Tensor>>{Gamma279, rdm1deriv_, rdm2deriv_, rdm3deriv_, f1_};
  auto task52 = make_shared<Task52>(tensor52, cindex);
  return make_shared<FutureTensor>(*Gamma279, task52);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma280_() {
  vector<IndexRange> Gamma280_index = {ci_, active_, active_, active_, active_};
  auto Gamma280 = make_shared<Tensor>(Gamma280_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor53 = vector<shared_ptr<Tensor>>{Gamma280, rdm1deriv_, rdm2deriv_};
  auto task53 = make_shared<Task53>(tensor53, cindex);
  return make_shared<FutureTensor>(*Gamma280, task53);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma282_() {
  vector<IndexRange> Gamma282_index = {ci_, active_, active_, active_, active_};
  auto Gamma282 = make_shared<Tensor>(Gamma282_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor54 = vector<shared_ptr<Tensor>>{Gamma282, rdm1deriv_, rdm2deriv_, rdm3deriv_, f1_};
  auto task54 = make_shared<Task54>(tensor54, cindex);
  return make_shared<FutureTensor>(*Gamma282, task54);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma283_() {
  vector<IndexRange> Gamma283_index = {ci_, active_, active_, active_, active_};
  auto Gamma283 = make_shared<Tensor>(Gamma283_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor55 = vector<shared_ptr<Tensor>>{Gamma283, rdm1deriv_, rdm2deriv_};
  auto task55 = make_shared<Task55>(tensor55, cindex);
  return make_shared<FutureTensor>(*Gamma283, task55);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma285_() {
  vector<IndexRange> Gamma285_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma285 = make_shared<Tensor>(Gamma285_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor56 = vector<shared_ptr<Tensor>>{Gamma285, rdm2deriv_, rdm3deriv_};
  auto task56 = make_shared<Task56>(tensor56, cindex);
  return make_shared<FutureTensor>(*Gamma285, task56);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma286_() {
  vector<IndexRange> Gamma286_index = {ci_, active_, active_};
  auto Gamma286 = make_shared<Tensor>(Gamma286_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor57 = vector<shared_ptr<Tensor>>{Gamma286, rdm1deriv_};
  auto task57 = make_shared<Task57>(tensor57, cindex);
  return make_shared<FutureTensor>(*Gamma286, task57);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma299_() {
  vector<IndexRange> Gamma299_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma299 = make_shared<Tensor>(Gamma299_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor58 = vector<shared_ptr<Tensor>>{Gamma299, rdm2deriv_, rdm3deriv_};
  auto task58 = make_shared<Task58>(tensor58, cindex);
  return make_shared<FutureTensor>(*Gamma299, task58);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma304_() {
  vector<IndexRange> Gamma304_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma304 = make_shared<Tensor>(Gamma304_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor59 = vector<shared_ptr<Tensor>>{Gamma304, rdm2deriv_, rdm3deriv_};
  auto task59 = make_shared<Task59>(tensor59, cindex);
  return make_shared<FutureTensor>(*Gamma304, task59);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma305_() {
  vector<IndexRange> Gamma305_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma305 = make_shared<Tensor>(Gamma305_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor60 = vector<shared_ptr<Tensor>>{Gamma305, rdm2deriv_, rdm3deriv_};
  auto task60 = make_shared<Task60>(tensor60, cindex);
  return make_shared<FutureTensor>(*Gamma305, task60);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma306_() {
  vector<IndexRange> Gamma306_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma306 = make_shared<Tensor>(Gamma306_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor61 = vector<shared_ptr<Tensor>>{Gamma306, rdm2deriv_, rdm3deriv_, rdm4deriv_, f1_};
  auto task61 = make_shared<Task61>(tensor61, cindex);
  return make_shared<FutureTensor>(*Gamma306, task61);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma307_() {
  vector<IndexRange> Gamma307_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma307 = make_shared<Tensor>(Gamma307_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor62 = vector<shared_ptr<Tensor>>{Gamma307, rdm2deriv_, rdm3deriv_};
  auto task62 = make_shared<Task62>(tensor62, cindex);
  return make_shared<FutureTensor>(*Gamma307, task62);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma308_() {
  vector<IndexRange> Gamma308_index = {ci_, active_, active_, active_, active_};
  auto Gamma308 = make_shared<Tensor>(Gamma308_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor63 = vector<shared_ptr<Tensor>>{Gamma308, rdm2deriv_};
  auto task63 = make_shared<Task63>(tensor63, cindex);
  return make_shared<FutureTensor>(*Gamma308, task63);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma317_() {
  vector<IndexRange> Gamma317_index = {ci_};
  auto Gamma317 = make_shared<Tensor>(Gamma317_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor64 = vector<shared_ptr<Tensor>>{Gamma317, rdm1deriv_, f1_};
  auto task64 = make_shared<Task64>(tensor64, cindex);
  return make_shared<FutureTensor>(*Gamma317, task64);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma329_() {
  vector<IndexRange> Gamma329_index = {ci_, active_, active_};
  auto Gamma329 = make_shared<Tensor>(Gamma329_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor65 = vector<shared_ptr<Tensor>>{Gamma329, rdm2deriv_, f1_};
  auto task65 = make_shared<Task65>(tensor65, cindex);
  return make_shared<FutureTensor>(*Gamma329, task65);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma340_() {
  vector<IndexRange> Gamma340_index = {ci_, active_, active_, active_, active_};
  auto Gamma340 = make_shared<Tensor>(Gamma340_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor66 = vector<shared_ptr<Tensor>>{Gamma340, rdm3deriv_, f1_};
  auto task66 = make_shared<Task66>(tensor66, cindex);
  return make_shared<FutureTensor>(*Gamma340, task66);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma355_() {
  vector<IndexRange> Gamma355_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma355 = make_shared<Tensor>(Gamma355_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor67 = vector<shared_ptr<Tensor>>{Gamma355, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task67 = make_shared<Task67>(tensor67, cindex);
  return make_shared<FutureTensor>(*Gamma355, task67);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma373_() {
  vector<IndexRange> Gamma373_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma373 = make_shared<Tensor>(Gamma373_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor68 = vector<shared_ptr<Tensor>>{Gamma373, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task68 = make_shared<Task68>(tensor68, cindex);
  return make_shared<FutureTensor>(*Gamma373, task68);
}

#endif
