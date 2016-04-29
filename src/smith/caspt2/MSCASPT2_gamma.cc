//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_gamma.cc
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
#include <src/smith/caspt2/MSCASPT2.h>
#include <src/smith/caspt2/MSCASPT2_tasks1.h>
#include <src/smith/caspt2/MSCASPT2_tasks2.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MSCASPT2;

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma0_() {
  vector<IndexRange> Gamma0_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma0 = make_shared<Tensor>(Gamma0_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor0 = vector<shared_ptr<Tensor>>{Gamma0, rdm0_, rdm1_, rdm2_, rdm3_};
  auto task0 = make_shared<Task0>(tensor0, pindex);
  return make_shared<FutureTensor>(*Gamma0, task0);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma31_() {
  vector<IndexRange> Gamma31_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma31 = make_shared<Tensor>(Gamma31_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor1 = vector<shared_ptr<Tensor>>{Gamma31, rdm1_, rdm2_, rdm3_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  return make_shared<FutureTensor>(*Gamma31, task1);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma34_() {
  vector<IndexRange> Gamma34_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma34 = make_shared<Tensor>(Gamma34_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor2 = vector<shared_ptr<Tensor>>{Gamma34, rdm1_, rdm2_, rdm3_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  return make_shared<FutureTensor>(*Gamma34, task2);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma92_() {
  vector<IndexRange> Gamma92_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma92 = make_shared<Tensor>(Gamma92_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor3 = vector<shared_ptr<Tensor>>{Gamma92, rdm3_};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  return make_shared<FutureTensor>(*Gamma92, task3);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma1_() {
  vector<IndexRange> Gamma1_index = {active_, active_, active_, active_};
  auto Gamma1 = make_shared<Tensor>(Gamma1_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor4 = vector<shared_ptr<Tensor>>{Gamma1, rdm0_, rdm1_, rdm2_};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  return make_shared<FutureTensor>(*Gamma1, task4);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma32_() {
  vector<IndexRange> Gamma32_index = {active_, active_, active_, active_};
  auto Gamma32 = make_shared<Tensor>(Gamma32_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor5 = vector<shared_ptr<Tensor>>{Gamma32, rdm1_, rdm2_};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  return make_shared<FutureTensor>(*Gamma32, task5);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma35_() {
  vector<IndexRange> Gamma35_index = {active_, active_, active_, active_};
  auto Gamma35 = make_shared<Tensor>(Gamma35_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor6 = vector<shared_ptr<Tensor>>{Gamma35, rdm1_, rdm2_};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  return make_shared<FutureTensor>(*Gamma35, task6);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma2_() {
  vector<IndexRange> Gamma2_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma2 = make_shared<Tensor>(Gamma2_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor7 = vector<shared_ptr<Tensor>>{Gamma2, rdm1_, rdm2_, rdm3_};
  auto task7 = make_shared<Task7>(tensor7, pindex);
  return make_shared<FutureTensor>(*Gamma2, task7);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma37_() {
  vector<IndexRange> Gamma37_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma37 = make_shared<Tensor>(Gamma37_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor8 = vector<shared_ptr<Tensor>>{Gamma37, rdm2_, rdm3_};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  return make_shared<FutureTensor>(*Gamma37, task8);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma3_() {
  vector<IndexRange> Gamma3_index = {active_, active_, active_, active_};
  auto Gamma3 = make_shared<Tensor>(Gamma3_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor9 = vector<shared_ptr<Tensor>>{Gamma3, rdm0_, rdm1_, rdm2_};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  return make_shared<FutureTensor>(*Gamma3, task9);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma4_() {
  vector<IndexRange> Gamma4_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma4 = make_shared<Tensor>(Gamma4_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor10 = vector<shared_ptr<Tensor>>{Gamma4, rdm1_, rdm2_, rdm3_};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  return make_shared<FutureTensor>(*Gamma4, task10);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma56_() {
  vector<IndexRange> Gamma56_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma56 = make_shared<Tensor>(Gamma56_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor11 = vector<shared_ptr<Tensor>>{Gamma56, rdm2_, rdm3_};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  return make_shared<FutureTensor>(*Gamma56, task11);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma57_() {
  vector<IndexRange> Gamma57_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma57 = make_shared<Tensor>(Gamma57_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor12 = vector<shared_ptr<Tensor>>{Gamma57, rdm2_, rdm3_};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  return make_shared<FutureTensor>(*Gamma57, task12);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma5_() {
  vector<IndexRange> Gamma5_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma5 = make_shared<Tensor>(Gamma5_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor13 = vector<shared_ptr<Tensor>>{Gamma5, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  return make_shared<FutureTensor>(*Gamma5, task13);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma58_() {
  vector<IndexRange> Gamma58_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma58 = make_shared<Tensor>(Gamma58_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor14 = vector<shared_ptr<Tensor>>{Gamma58, rdm2_, rdm3_, rdm4_};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  return make_shared<FutureTensor>(*Gamma58, task14);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma6_() {
  vector<IndexRange> Gamma6_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma6 = make_shared<Tensor>(Gamma6_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor15 = vector<shared_ptr<Tensor>>{Gamma6, rdm1_, rdm2_, rdm3_};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  return make_shared<FutureTensor>(*Gamma6, task15);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma7_() {
  vector<IndexRange> Gamma7_index = {active_, active_, active_, active_};
  auto Gamma7 = make_shared<Tensor>(Gamma7_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor16 = vector<shared_ptr<Tensor>>{Gamma7, rdm1_, rdm2_};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  return make_shared<FutureTensor>(*Gamma7, task16);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma60_() {
  vector<IndexRange> Gamma60_index = {active_, active_, active_, active_};
  auto Gamma60 = make_shared<Tensor>(Gamma60_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor17 = vector<shared_ptr<Tensor>>{Gamma60, rdm2_};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  return make_shared<FutureTensor>(*Gamma60, task17);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma9_() {
  vector<IndexRange> Gamma9_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma9 = make_shared<Tensor>(Gamma9_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor18 = vector<shared_ptr<Tensor>>{Gamma9, rdm1_, rdm2_, rdm3_};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  return make_shared<FutureTensor>(*Gamma9, task18);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma62_() {
  vector<IndexRange> Gamma62_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma62 = make_shared<Tensor>(Gamma62_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor19 = vector<shared_ptr<Tensor>>{Gamma62, rdm2_, rdm3_};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  return make_shared<FutureTensor>(*Gamma62, task19);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma12_() {
  vector<IndexRange> Gamma12_index = {active_, active_, active_, active_};
  auto Gamma12 = make_shared<Tensor>(Gamma12_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor20 = vector<shared_ptr<Tensor>>{Gamma12, rdm1_, rdm2_};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  return make_shared<FutureTensor>(*Gamma12, task20);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma65_() {
  vector<IndexRange> Gamma65_index = {active_, active_};
  auto Gamma65 = make_shared<Tensor>(Gamma65_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor21 = vector<shared_ptr<Tensor>>{Gamma65, rdm1_};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  return make_shared<FutureTensor>(*Gamma65, task21);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma14_() {
  vector<IndexRange> Gamma14_index = {active_, active_, active_, active_};
  auto Gamma14 = make_shared<Tensor>(Gamma14_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor22 = vector<shared_ptr<Tensor>>{Gamma14, rdm0_, rdm1_, rdm2_};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  return make_shared<FutureTensor>(*Gamma14, task22);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma16_() {
  vector<IndexRange> Gamma16_index = {active_, active_};
  auto Gamma16 = make_shared<Tensor>(Gamma16_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor23 = vector<shared_ptr<Tensor>>{Gamma16, rdm0_, rdm1_};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  return make_shared<FutureTensor>(*Gamma16, task23);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma22_() {
  vector<IndexRange> Gamma22_index = {active_, active_, active_, active_};
  auto Gamma22 = make_shared<Tensor>(Gamma22_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor24 = vector<shared_ptr<Tensor>>{Gamma22, rdm1_, rdm2_};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  return make_shared<FutureTensor>(*Gamma22, task24);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma28_() {
  vector<IndexRange> Gamma28_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma28 = make_shared<Tensor>(Gamma28_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor25 = vector<shared_ptr<Tensor>>{Gamma28, rdm1_, rdm2_, rdm3_};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  return make_shared<FutureTensor>(*Gamma28, task25);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma29_() {
  vector<IndexRange> Gamma29_index = {active_, active_, active_, active_};
  auto Gamma29 = make_shared<Tensor>(Gamma29_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor26 = vector<shared_ptr<Tensor>>{Gamma29, rdm1_, rdm2_};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  return make_shared<FutureTensor>(*Gamma29, task26);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma51_() {
  vector<IndexRange> Gamma51_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma51 = make_shared<Tensor>(Gamma51_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor27 = vector<shared_ptr<Tensor>>{Gamma51, rdm2_, rdm3_};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  return make_shared<FutureTensor>(*Gamma51, task27);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma110_() {
  vector<IndexRange> Gamma110_index = {ci_, active_, active_, active_, active_};
  auto Gamma110 = make_shared<Tensor>(Gamma110_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor28 = vector<shared_ptr<Tensor>>{Gamma110, rdm0deriv_, rdm1deriv_, rdm2deriv_, rdm3deriv_, f1_};
  auto task28 = make_shared<Task28>(tensor28, cindex);
  return make_shared<FutureTensor>(*Gamma110, task28);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma111_() {
  vector<IndexRange> Gamma111_index = {ci_, active_, active_, active_, active_};
  auto Gamma111 = make_shared<Tensor>(Gamma111_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor29 = vector<shared_ptr<Tensor>>{Gamma111, rdm0deriv_, rdm1deriv_, rdm2deriv_};
  auto task29 = make_shared<Task29>(tensor29, cindex);
  return make_shared<FutureTensor>(*Gamma111, task29);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma112_() {
  vector<IndexRange> Gamma112_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma112 = make_shared<Tensor>(Gamma112_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor30 = vector<shared_ptr<Tensor>>{Gamma112, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task30 = make_shared<Task30>(tensor30, cindex);
  return make_shared<FutureTensor>(*Gamma112, task30);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma113_() {
  vector<IndexRange> Gamma113_index = {ci_, active_, active_, active_, active_};
  auto Gamma113 = make_shared<Tensor>(Gamma113_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor31 = vector<shared_ptr<Tensor>>{Gamma113, rdm0deriv_, rdm1deriv_, rdm2deriv_};
  auto task31 = make_shared<Task31>(tensor31, cindex);
  return make_shared<FutureTensor>(*Gamma113, task31);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma114_() {
  vector<IndexRange> Gamma114_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma114 = make_shared<Tensor>(Gamma114_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor32 = vector<shared_ptr<Tensor>>{Gamma114, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task32 = make_shared<Task32>(tensor32, cindex);
  return make_shared<FutureTensor>(*Gamma114, task32);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma115_() {
  vector<IndexRange> Gamma115_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma115 = make_shared<Tensor>(Gamma115_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor33 = vector<shared_ptr<Tensor>>{Gamma115, rdm1deriv_, rdm2deriv_, rdm3deriv_, rdm4deriv_, f1_};
  auto task33 = make_shared<Task33>(tensor33, cindex);
  return make_shared<FutureTensor>(*Gamma115, task33);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma116_() {
  vector<IndexRange> Gamma116_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma116 = make_shared<Tensor>(Gamma116_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor34 = vector<shared_ptr<Tensor>>{Gamma116, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task34 = make_shared<Task34>(tensor34, cindex);
  return make_shared<FutureTensor>(*Gamma116, task34);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma117_() {
  vector<IndexRange> Gamma117_index = {ci_, active_, active_, active_, active_};
  auto Gamma117 = make_shared<Tensor>(Gamma117_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor35 = vector<shared_ptr<Tensor>>{Gamma117, rdm1deriv_, rdm2deriv_};
  auto task35 = make_shared<Task35>(tensor35, cindex);
  return make_shared<FutureTensor>(*Gamma117, task35);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma119_() {
  vector<IndexRange> Gamma119_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma119 = make_shared<Tensor>(Gamma119_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor36 = vector<shared_ptr<Tensor>>{Gamma119, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task36 = make_shared<Task36>(tensor36, cindex);
  return make_shared<FutureTensor>(*Gamma119, task36);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma122_() {
  vector<IndexRange> Gamma122_index = {ci_, active_, active_, active_, active_};
  auto Gamma122 = make_shared<Tensor>(Gamma122_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor37 = vector<shared_ptr<Tensor>>{Gamma122, rdm1deriv_, rdm2deriv_};
  auto task37 = make_shared<Task37>(tensor37, cindex);
  return make_shared<FutureTensor>(*Gamma122, task37);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma124_() {
  vector<IndexRange> Gamma124_index = {ci_, active_, active_};
  auto Gamma124 = make_shared<Tensor>(Gamma124_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor38 = vector<shared_ptr<Tensor>>{Gamma124, rdm0deriv_, rdm1deriv_, rdm2deriv_, f1_};
  auto task38 = make_shared<Task38>(tensor38, cindex);
  return make_shared<FutureTensor>(*Gamma124, task38);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma126_() {
  vector<IndexRange> Gamma126_index = {ci_, active_, active_};
  auto Gamma126 = make_shared<Tensor>(Gamma126_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor39 = vector<shared_ptr<Tensor>>{Gamma126, rdm0deriv_, rdm1deriv_};
  auto task39 = make_shared<Task39>(tensor39, cindex);
  return make_shared<FutureTensor>(*Gamma126, task39);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma132_() {
  vector<IndexRange> Gamma132_index = {ci_, active_, active_, active_, active_};
  auto Gamma132 = make_shared<Tensor>(Gamma132_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor40 = vector<shared_ptr<Tensor>>{Gamma132, rdm1deriv_, rdm2deriv_};
  auto task40 = make_shared<Task40>(tensor40, cindex);
  return make_shared<FutureTensor>(*Gamma132, task40);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma138_() {
  vector<IndexRange> Gamma138_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma138 = make_shared<Tensor>(Gamma138_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor41 = vector<shared_ptr<Tensor>>{Gamma138, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task41 = make_shared<Task41>(tensor41, cindex);
  return make_shared<FutureTensor>(*Gamma138, task41);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma139_() {
  vector<IndexRange> Gamma139_index = {ci_, active_, active_, active_, active_};
  auto Gamma139 = make_shared<Tensor>(Gamma139_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor42 = vector<shared_ptr<Tensor>>{Gamma139, rdm1deriv_, rdm2deriv_};
  auto task42 = make_shared<Task42>(tensor42, cindex);
  return make_shared<FutureTensor>(*Gamma139, task42);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma141_() {
  vector<IndexRange> Gamma141_index = {ci_, active_, active_, active_, active_};
  auto Gamma141 = make_shared<Tensor>(Gamma141_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor43 = vector<shared_ptr<Tensor>>{Gamma141, rdm1deriv_, rdm2deriv_, rdm3deriv_, f1_};
  auto task43 = make_shared<Task43>(tensor43, cindex);
  return make_shared<FutureTensor>(*Gamma141, task43);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma142_() {
  vector<IndexRange> Gamma142_index = {ci_, active_, active_, active_, active_};
  auto Gamma142 = make_shared<Tensor>(Gamma142_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor44 = vector<shared_ptr<Tensor>>{Gamma142, rdm1deriv_, rdm2deriv_};
  auto task44 = make_shared<Task44>(tensor44, cindex);
  return make_shared<FutureTensor>(*Gamma142, task44);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma144_() {
  vector<IndexRange> Gamma144_index = {ci_, active_, active_, active_, active_};
  auto Gamma144 = make_shared<Tensor>(Gamma144_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor45 = vector<shared_ptr<Tensor>>{Gamma144, rdm1deriv_, rdm2deriv_, rdm3deriv_, f1_};
  auto task45 = make_shared<Task45>(tensor45, cindex);
  return make_shared<FutureTensor>(*Gamma144, task45);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma145_() {
  vector<IndexRange> Gamma145_index = {ci_, active_, active_, active_, active_};
  auto Gamma145 = make_shared<Tensor>(Gamma145_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor46 = vector<shared_ptr<Tensor>>{Gamma145, rdm1deriv_, rdm2deriv_};
  auto task46 = make_shared<Task46>(tensor46, cindex);
  return make_shared<FutureTensor>(*Gamma145, task46);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma147_() {
  vector<IndexRange> Gamma147_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma147 = make_shared<Tensor>(Gamma147_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor47 = vector<shared_ptr<Tensor>>{Gamma147, rdm2deriv_, rdm3deriv_};
  auto task47 = make_shared<Task47>(tensor47, cindex);
  return make_shared<FutureTensor>(*Gamma147, task47);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma148_() {
  vector<IndexRange> Gamma148_index = {ci_, active_, active_};
  auto Gamma148 = make_shared<Tensor>(Gamma148_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor48 = vector<shared_ptr<Tensor>>{Gamma148, rdm1deriv_};
  auto task48 = make_shared<Task48>(tensor48, cindex);
  return make_shared<FutureTensor>(*Gamma148, task48);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma161_() {
  vector<IndexRange> Gamma161_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma161 = make_shared<Tensor>(Gamma161_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor49 = vector<shared_ptr<Tensor>>{Gamma161, rdm2deriv_, rdm3deriv_};
  auto task49 = make_shared<Task49>(tensor49, cindex);
  return make_shared<FutureTensor>(*Gamma161, task49);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma166_() {
  vector<IndexRange> Gamma166_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma166 = make_shared<Tensor>(Gamma166_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor50 = vector<shared_ptr<Tensor>>{Gamma166, rdm2deriv_, rdm3deriv_};
  auto task50 = make_shared<Task50>(tensor50, cindex);
  return make_shared<FutureTensor>(*Gamma166, task50);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma167_() {
  vector<IndexRange> Gamma167_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma167 = make_shared<Tensor>(Gamma167_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor51 = vector<shared_ptr<Tensor>>{Gamma167, rdm2deriv_, rdm3deriv_};
  auto task51 = make_shared<Task51>(tensor51, cindex);
  return make_shared<FutureTensor>(*Gamma167, task51);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma168_() {
  vector<IndexRange> Gamma168_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma168 = make_shared<Tensor>(Gamma168_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor52 = vector<shared_ptr<Tensor>>{Gamma168, rdm2deriv_, rdm3deriv_, rdm4deriv_, f1_};
  auto task52 = make_shared<Task52>(tensor52, cindex);
  return make_shared<FutureTensor>(*Gamma168, task52);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma169_() {
  vector<IndexRange> Gamma169_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma169 = make_shared<Tensor>(Gamma169_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor53 = vector<shared_ptr<Tensor>>{Gamma169, rdm2deriv_, rdm3deriv_};
  auto task53 = make_shared<Task53>(tensor53, cindex);
  return make_shared<FutureTensor>(*Gamma169, task53);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma170_() {
  vector<IndexRange> Gamma170_index = {ci_, active_, active_, active_, active_};
  auto Gamma170 = make_shared<Tensor>(Gamma170_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor54 = vector<shared_ptr<Tensor>>{Gamma170, rdm2deriv_};
  auto task54 = make_shared<Task54>(tensor54, cindex);
  return make_shared<FutureTensor>(*Gamma170, task54);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma179_() {
  vector<IndexRange> Gamma179_index = {ci_};
  auto Gamma179 = make_shared<Tensor>(Gamma179_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor55 = vector<shared_ptr<Tensor>>{Gamma179, rdm1deriv_, f1_};
  auto task55 = make_shared<Task55>(tensor55, cindex);
  return make_shared<FutureTensor>(*Gamma179, task55);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma191_() {
  vector<IndexRange> Gamma191_index = {ci_, active_, active_};
  auto Gamma191 = make_shared<Tensor>(Gamma191_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor56 = vector<shared_ptr<Tensor>>{Gamma191, rdm2deriv_, f1_};
  auto task56 = make_shared<Task56>(tensor56, cindex);
  return make_shared<FutureTensor>(*Gamma191, task56);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma202_() {
  vector<IndexRange> Gamma202_index = {ci_, active_, active_, active_, active_};
  auto Gamma202 = make_shared<Tensor>(Gamma202_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor57 = vector<shared_ptr<Tensor>>{Gamma202, rdm3deriv_, f1_};
  auto task57 = make_shared<Task57>(tensor57, cindex);
  return make_shared<FutureTensor>(*Gamma202, task57);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma217_() {
  vector<IndexRange> Gamma217_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma217 = make_shared<Tensor>(Gamma217_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor58 = vector<shared_ptr<Tensor>>{Gamma217, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task58 = make_shared<Task58>(tensor58, cindex);
  return make_shared<FutureTensor>(*Gamma217, task58);
}

shared_ptr<FutureTensor> MSCASPT2::MSCASPT2::Gamma239_() {
  vector<IndexRange> Gamma239_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma239 = make_shared<Tensor>(Gamma239_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto tensor59 = vector<shared_ptr<Tensor>>{Gamma239, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task59 = make_shared<Task59>(tensor59, cindex);
  return make_shared<FutureTensor>(*Gamma239, task59);
}

#endif
