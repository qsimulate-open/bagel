//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2.cc
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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma0_() {
  vector<IndexRange> Gamma0_index = {active_, active_, active_, active_};
  auto Gamma0 = make_shared<Tensor>(Gamma0_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor0 = {Gamma0, rdm1_, rdm2_, rdm3_, f1_};
  auto task0 = make_shared<Task0>(tensor0, pindex);
  return make_shared<FutureTensor>(*Gamma0, task0);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma92_() {
  vector<IndexRange> Gamma92_index = {active_, active_, active_, active_};
  auto Gamma92 = make_shared<Tensor>(Gamma92_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor1 = {Gamma92, rdm1_, rdm2_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  return make_shared<FutureTensor>(*Gamma92, task1);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma2_() {
  vector<IndexRange> Gamma2_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma2 = make_shared<Tensor>(Gamma2_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor2 = {Gamma2, rdm1_, rdm2_, rdm3_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  return make_shared<FutureTensor>(*Gamma2, task2);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma3_() {
  vector<IndexRange> Gamma3_index = {active_, active_, active_, active_};
  auto Gamma3 = make_shared<Tensor>(Gamma3_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor3 = {Gamma3, rdm1_, rdm2_};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  return make_shared<FutureTensor>(*Gamma3, task3);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma4_() {
  vector<IndexRange> Gamma4_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma4 = make_shared<Tensor>(Gamma4_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor4 = {Gamma4, rdm1_, rdm2_, rdm3_};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  return make_shared<FutureTensor>(*Gamma4, task4);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma5_() {
  vector<IndexRange> Gamma5_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma5 = make_shared<Tensor>(Gamma5_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor5 = {Gamma5, rdm1_, rdm2_, rdm3_, rdm4_, f1_};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  return make_shared<FutureTensor>(*Gamma5, task5);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma6_() {
  vector<IndexRange> Gamma6_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma6 = make_shared<Tensor>(Gamma6_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor6 = {Gamma6, rdm1_, rdm2_, rdm3_};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  return make_shared<FutureTensor>(*Gamma6, task6);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma7_() {
  vector<IndexRange> Gamma7_index = {active_, active_, active_, active_};
  auto Gamma7 = make_shared<Tensor>(Gamma7_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor7 = {Gamma7, rdm1_, rdm2_};
  auto task7 = make_shared<Task7>(tensor7, pindex);
  return make_shared<FutureTensor>(*Gamma7, task7);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma9_() {
  vector<IndexRange> Gamma9_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma9 = make_shared<Tensor>(Gamma9_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor8 = {Gamma9, rdm1_, rdm2_, rdm3_};
  auto task8 = make_shared<Task8>(tensor8, pindex);
  return make_shared<FutureTensor>(*Gamma9, task8);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma105_() {
  vector<IndexRange> Gamma105_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma105 = make_shared<Tensor>(Gamma105_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor9 = {Gamma105, rdm1_, rdm2_, rdm3_};
  auto task9 = make_shared<Task9>(tensor9, pindex);
  return make_shared<FutureTensor>(*Gamma105, task9);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma12_() {
  vector<IndexRange> Gamma12_index = {active_, active_, active_, active_};
  auto Gamma12 = make_shared<Tensor>(Gamma12_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor10 = {Gamma12, rdm1_, rdm2_};
  auto task10 = make_shared<Task10>(tensor10, pindex);
  return make_shared<FutureTensor>(*Gamma12, task10);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma14_() {
  vector<IndexRange> Gamma14_index = {active_, active_};
  auto Gamma14 = make_shared<Tensor>(Gamma14_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor11 = {Gamma14, rdm1_, rdm2_, f1_};
  auto task11 = make_shared<Task11>(tensor11, pindex);
  return make_shared<FutureTensor>(*Gamma14, task11);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma16_() {
  vector<IndexRange> Gamma16_index = {active_, active_};
  auto Gamma16 = make_shared<Tensor>(Gamma16_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor12 = {Gamma16, rdm1_};
  auto task12 = make_shared<Task12>(tensor12, pindex);
  return make_shared<FutureTensor>(*Gamma16, task12);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma22_() {
  vector<IndexRange> Gamma22_index = {active_, active_, active_, active_};
  auto Gamma22 = make_shared<Tensor>(Gamma22_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor13 = {Gamma22, rdm1_, rdm2_};
  auto task13 = make_shared<Task13>(tensor13, pindex);
  return make_shared<FutureTensor>(*Gamma22, task13);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma28_() {
  vector<IndexRange> Gamma28_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma28 = make_shared<Tensor>(Gamma28_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor14 = {Gamma28, rdm1_, rdm2_, rdm3_};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  return make_shared<FutureTensor>(*Gamma28, task14);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma29_() {
  vector<IndexRange> Gamma29_index = {active_, active_, active_, active_};
  auto Gamma29 = make_shared<Tensor>(Gamma29_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor15 = {Gamma29, rdm1_, rdm2_};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  return make_shared<FutureTensor>(*Gamma29, task15);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma31_() {
  vector<IndexRange> Gamma31_index = {active_, active_, active_, active_};
  auto Gamma31 = make_shared<Tensor>(Gamma31_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor16 = {Gamma31, rdm1_, rdm2_, rdm3_, f1_};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  return make_shared<FutureTensor>(*Gamma31, task16);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma32_() {
  vector<IndexRange> Gamma32_index = {active_, active_, active_, active_};
  auto Gamma32 = make_shared<Tensor>(Gamma32_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor17 = {Gamma32, rdm1_, rdm2_};
  auto task17 = make_shared<Task17>(tensor17, pindex);
  return make_shared<FutureTensor>(*Gamma32, task17);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma34_() {
  vector<IndexRange> Gamma34_index = {active_, active_, active_, active_};
  auto Gamma34 = make_shared<Tensor>(Gamma34_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor18 = {Gamma34, rdm1_, rdm2_, rdm3_, f1_};
  auto task18 = make_shared<Task18>(tensor18, pindex);
  return make_shared<FutureTensor>(*Gamma34, task18);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma35_() {
  vector<IndexRange> Gamma35_index = {active_, active_, active_, active_};
  auto Gamma35 = make_shared<Tensor>(Gamma35_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor19 = {Gamma35, rdm1_, rdm2_};
  auto task19 = make_shared<Task19>(tensor19, pindex);
  return make_shared<FutureTensor>(*Gamma35, task19);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma37_() {
  vector<IndexRange> Gamma37_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma37 = make_shared<Tensor>(Gamma37_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor20 = {Gamma37, rdm2_, rdm3_};
  auto task20 = make_shared<Task20>(tensor20, pindex);
  return make_shared<FutureTensor>(*Gamma37, task20);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma38_() {
  vector<IndexRange> Gamma38_index = {active_, active_};
  auto Gamma38 = make_shared<Tensor>(Gamma38_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor21 = {Gamma38, rdm1_};
  auto task21 = make_shared<Task21>(tensor21, pindex);
  return make_shared<FutureTensor>(*Gamma38, task21);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma51_() {
  vector<IndexRange> Gamma51_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma51 = make_shared<Tensor>(Gamma51_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor22 = {Gamma51, rdm2_, rdm3_};
  auto task22 = make_shared<Task22>(tensor22, pindex);
  return make_shared<FutureTensor>(*Gamma51, task22);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma56_() {
  vector<IndexRange> Gamma56_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma56 = make_shared<Tensor>(Gamma56_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor23 = {Gamma56, rdm2_, rdm3_};
  auto task23 = make_shared<Task23>(tensor23, pindex);
  return make_shared<FutureTensor>(*Gamma56, task23);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma57_() {
  vector<IndexRange> Gamma57_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma57 = make_shared<Tensor>(Gamma57_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor24 = {Gamma57, rdm2_, rdm3_};
  auto task24 = make_shared<Task24>(tensor24, pindex);
  return make_shared<FutureTensor>(*Gamma57, task24);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma58_() {
  vector<IndexRange> Gamma58_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma58 = make_shared<Tensor>(Gamma58_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor25 = {Gamma58, rdm2_, rdm3_, rdm4_, f1_};
  auto task25 = make_shared<Task25>(tensor25, pindex);
  return make_shared<FutureTensor>(*Gamma58, task25);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma59_() {
  vector<IndexRange> Gamma59_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma59 = make_shared<Tensor>(Gamma59_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor26 = {Gamma59, rdm2_, rdm3_};
  auto task26 = make_shared<Task26>(tensor26, pindex);
  return make_shared<FutureTensor>(*Gamma59, task26);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma60_() {
  vector<IndexRange> Gamma60_index = {active_, active_, active_, active_};
  auto Gamma60 = make_shared<Tensor>(Gamma60_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor27 = {Gamma60, rdm2_};
  auto task27 = make_shared<Task27>(tensor27, pindex);
  return make_shared<FutureTensor>(*Gamma60, task27);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma79_() {
  vector<IndexRange> Gamma79_index = {active_, active_};
  auto Gamma79 = make_shared<Tensor>(Gamma79_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor28 = {Gamma79, rdm2_, f1_};
  auto task28 = make_shared<Task28>(tensor28, pindex);
  return make_shared<FutureTensor>(*Gamma79, task28);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma90_() {
  vector<IndexRange> Gamma90_index = {active_, active_, active_, active_};
  auto Gamma90 = make_shared<Tensor>(Gamma90_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor29 = {Gamma90, rdm3_, f1_};
  auto task29 = make_shared<Task29>(tensor29, pindex);
  return make_shared<FutureTensor>(*Gamma90, task29);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma160_() {
  vector<IndexRange> Gamma160_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma160 = make_shared<Tensor>(Gamma160_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor30 = {Gamma160, rdm1_, rdm2_, rdm3_};
  auto task30 = make_shared<Task30>(tensor30, pindex);
  return make_shared<FutureTensor>(*Gamma160, task30);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma191_() {
  vector<IndexRange> Gamma191_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma191 = make_shared<Tensor>(Gamma191_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor31 = {Gamma191, rdm1_, rdm2_, rdm3_};
  auto task31 = make_shared<Task31>(tensor31, pindex);
  return make_shared<FutureTensor>(*Gamma191, task31);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma194_() {
  vector<IndexRange> Gamma194_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma194 = make_shared<Tensor>(Gamma194_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor32 = {Gamma194, rdm1_, rdm2_, rdm3_};
  auto task32 = make_shared<Task32>(tensor32, pindex);
  return make_shared<FutureTensor>(*Gamma194, task32);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma252_() {
  vector<IndexRange> Gamma252_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma252 = make_shared<Tensor>(Gamma252_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor33 = {Gamma252, rdm3_};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  return make_shared<FutureTensor>(*Gamma252, task33);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma165_() {
  vector<IndexRange> Gamma165_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma165 = make_shared<Tensor>(Gamma165_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor34 = {Gamma165, rdm1_, rdm2_, rdm3_, rdm4_};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  return make_shared<FutureTensor>(*Gamma165, task34);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma218_() {
  vector<IndexRange> Gamma218_index = {active_, active_, active_, active_, active_, active_, active_, active_};
  auto Gamma218 = make_shared<Tensor>(Gamma218_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor35 = {Gamma218, rdm2_, rdm3_, rdm4_};
  auto task35 = make_shared<Task35>(tensor35, pindex);
  return make_shared<FutureTensor>(*Gamma218, task35);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma174_() {
  vector<IndexRange> Gamma174_index = {active_, active_, active_, active_};
  auto Gamma174 = make_shared<Tensor>(Gamma174_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor36 = {Gamma174, rdm1_, rdm2_};
  auto task36 = make_shared<Task36>(tensor36, pindex);
  return make_shared<FutureTensor>(*Gamma174, task36);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma270_() {
  vector<IndexRange> Gamma270_index = {ci_, active_, active_, active_, active_};
  auto Gamma270 = make_shared<Tensor>(Gamma270_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor37 = {Gamma270, rdm0deriv_, rdm1deriv_, rdm2deriv_, rdm3deriv_, f1_};
  auto task37 = make_shared<Task37>(tensor37, cindex);
  return make_shared<FutureTensor>(*Gamma270, task37);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma271_() {
  vector<IndexRange> Gamma271_index = {ci_, active_, active_, active_, active_};
  auto Gamma271 = make_shared<Tensor>(Gamma271_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor38 = {Gamma271, rdm0deriv_, rdm1deriv_, rdm2deriv_};
  auto task38 = make_shared<Task38>(tensor38, cindex);
  return make_shared<FutureTensor>(*Gamma271, task38);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma272_() {
  vector<IndexRange> Gamma272_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma272 = make_shared<Tensor>(Gamma272_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor39 = {Gamma272, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task39 = make_shared<Task39>(tensor39, cindex);
  return make_shared<FutureTensor>(*Gamma272, task39);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma273_() {
  vector<IndexRange> Gamma273_index = {ci_, active_, active_, active_, active_};
  auto Gamma273 = make_shared<Tensor>(Gamma273_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor40 = {Gamma273, rdm0deriv_, rdm1deriv_, rdm2deriv_};
  auto task40 = make_shared<Task40>(tensor40, cindex);
  return make_shared<FutureTensor>(*Gamma273, task40);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma274_() {
  vector<IndexRange> Gamma274_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma274 = make_shared<Tensor>(Gamma274_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor41 = {Gamma274, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task41 = make_shared<Task41>(tensor41, cindex);
  return make_shared<FutureTensor>(*Gamma274, task41);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma275_() {
  vector<IndexRange> Gamma275_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma275 = make_shared<Tensor>(Gamma275_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor42 = {Gamma275, rdm1deriv_, rdm2deriv_, rdm3deriv_, rdm4deriv_, f1_};
  auto task42 = make_shared<Task42>(tensor42, cindex);
  return make_shared<FutureTensor>(*Gamma275, task42);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma276_() {
  vector<IndexRange> Gamma276_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma276 = make_shared<Tensor>(Gamma276_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor43 = {Gamma276, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task43 = make_shared<Task43>(tensor43, cindex);
  return make_shared<FutureTensor>(*Gamma276, task43);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma277_() {
  vector<IndexRange> Gamma277_index = {ci_, active_, active_, active_, active_};
  auto Gamma277 = make_shared<Tensor>(Gamma277_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor44 = {Gamma277, rdm1deriv_, rdm2deriv_};
  auto task44 = make_shared<Task44>(tensor44, cindex);
  return make_shared<FutureTensor>(*Gamma277, task44);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma279_() {
  vector<IndexRange> Gamma279_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma279 = make_shared<Tensor>(Gamma279_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor45 = {Gamma279, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task45 = make_shared<Task45>(tensor45, cindex);
  return make_shared<FutureTensor>(*Gamma279, task45);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma282_() {
  vector<IndexRange> Gamma282_index = {ci_, active_, active_, active_, active_};
  auto Gamma282 = make_shared<Tensor>(Gamma282_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor46 = {Gamma282, rdm1deriv_, rdm2deriv_};
  auto task46 = make_shared<Task46>(tensor46, cindex);
  return make_shared<FutureTensor>(*Gamma282, task46);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma284_() {
  vector<IndexRange> Gamma284_index = {ci_, active_, active_};
  auto Gamma284 = make_shared<Tensor>(Gamma284_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor47 = {Gamma284, rdm0deriv_, rdm1deriv_, rdm2deriv_, f1_};
  auto task47 = make_shared<Task47>(tensor47, cindex);
  return make_shared<FutureTensor>(*Gamma284, task47);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma286_() {
  vector<IndexRange> Gamma286_index = {ci_, active_, active_};
  auto Gamma286 = make_shared<Tensor>(Gamma286_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor48 = {Gamma286, rdm0deriv_, rdm1deriv_};
  auto task48 = make_shared<Task48>(tensor48, cindex);
  return make_shared<FutureTensor>(*Gamma286, task48);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma292_() {
  vector<IndexRange> Gamma292_index = {ci_, active_, active_, active_, active_};
  auto Gamma292 = make_shared<Tensor>(Gamma292_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor49 = {Gamma292, rdm1deriv_, rdm2deriv_};
  auto task49 = make_shared<Task49>(tensor49, cindex);
  return make_shared<FutureTensor>(*Gamma292, task49);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma298_() {
  vector<IndexRange> Gamma298_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma298 = make_shared<Tensor>(Gamma298_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor50 = {Gamma298, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task50 = make_shared<Task50>(tensor50, cindex);
  return make_shared<FutureTensor>(*Gamma298, task50);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma299_() {
  vector<IndexRange> Gamma299_index = {ci_, active_, active_, active_, active_};
  auto Gamma299 = make_shared<Tensor>(Gamma299_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor51 = {Gamma299, rdm1deriv_, rdm2deriv_};
  auto task51 = make_shared<Task51>(tensor51, cindex);
  return make_shared<FutureTensor>(*Gamma299, task51);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma301_() {
  vector<IndexRange> Gamma301_index = {ci_, active_, active_, active_, active_};
  auto Gamma301 = make_shared<Tensor>(Gamma301_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor52 = {Gamma301, rdm1deriv_, rdm2deriv_, rdm3deriv_, f1_};
  auto task52 = make_shared<Task52>(tensor52, cindex);
  return make_shared<FutureTensor>(*Gamma301, task52);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma302_() {
  vector<IndexRange> Gamma302_index = {ci_, active_, active_, active_, active_};
  auto Gamma302 = make_shared<Tensor>(Gamma302_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor53 = {Gamma302, rdm1deriv_, rdm2deriv_};
  auto task53 = make_shared<Task53>(tensor53, cindex);
  return make_shared<FutureTensor>(*Gamma302, task53);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma304_() {
  vector<IndexRange> Gamma304_index = {ci_, active_, active_, active_, active_};
  auto Gamma304 = make_shared<Tensor>(Gamma304_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor54 = {Gamma304, rdm1deriv_, rdm2deriv_, rdm3deriv_, f1_};
  auto task54 = make_shared<Task54>(tensor54, cindex);
  return make_shared<FutureTensor>(*Gamma304, task54);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma305_() {
  vector<IndexRange> Gamma305_index = {ci_, active_, active_, active_, active_};
  auto Gamma305 = make_shared<Tensor>(Gamma305_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor55 = {Gamma305, rdm1deriv_, rdm2deriv_};
  auto task55 = make_shared<Task55>(tensor55, cindex);
  return make_shared<FutureTensor>(*Gamma305, task55);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma307_() {
  vector<IndexRange> Gamma307_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma307 = make_shared<Tensor>(Gamma307_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor56 = {Gamma307, rdm2deriv_, rdm3deriv_};
  auto task56 = make_shared<Task56>(tensor56, cindex);
  return make_shared<FutureTensor>(*Gamma307, task56);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma308_() {
  vector<IndexRange> Gamma308_index = {ci_, active_, active_};
  auto Gamma308 = make_shared<Tensor>(Gamma308_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor57 = {Gamma308, rdm1deriv_};
  auto task57 = make_shared<Task57>(tensor57, cindex);
  return make_shared<FutureTensor>(*Gamma308, task57);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma321_() {
  vector<IndexRange> Gamma321_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma321 = make_shared<Tensor>(Gamma321_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor58 = {Gamma321, rdm2deriv_, rdm3deriv_};
  auto task58 = make_shared<Task58>(tensor58, cindex);
  return make_shared<FutureTensor>(*Gamma321, task58);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma326_() {
  vector<IndexRange> Gamma326_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma326 = make_shared<Tensor>(Gamma326_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor59 = {Gamma326, rdm2deriv_, rdm3deriv_};
  auto task59 = make_shared<Task59>(tensor59, cindex);
  return make_shared<FutureTensor>(*Gamma326, task59);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma327_() {
  vector<IndexRange> Gamma327_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma327 = make_shared<Tensor>(Gamma327_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor60 = {Gamma327, rdm2deriv_, rdm3deriv_};
  auto task60 = make_shared<Task60>(tensor60, cindex);
  return make_shared<FutureTensor>(*Gamma327, task60);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma328_() {
  vector<IndexRange> Gamma328_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma328 = make_shared<Tensor>(Gamma328_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor61 = {Gamma328, rdm2deriv_, rdm3deriv_, rdm4deriv_, f1_};
  auto task61 = make_shared<Task61>(tensor61, cindex);
  return make_shared<FutureTensor>(*Gamma328, task61);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma329_() {
  vector<IndexRange> Gamma329_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma329 = make_shared<Tensor>(Gamma329_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor62 = {Gamma329, rdm2deriv_, rdm3deriv_};
  auto task62 = make_shared<Task62>(tensor62, cindex);
  return make_shared<FutureTensor>(*Gamma329, task62);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma330_() {
  vector<IndexRange> Gamma330_index = {ci_, active_, active_, active_, active_};
  auto Gamma330 = make_shared<Tensor>(Gamma330_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor63 = {Gamma330, rdm2deriv_};
  auto task63 = make_shared<Task63>(tensor63, cindex);
  return make_shared<FutureTensor>(*Gamma330, task63);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma339_() {
  vector<IndexRange> Gamma339_index = {ci_};
  auto Gamma339 = make_shared<Tensor>(Gamma339_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor64 = {Gamma339, rdm1deriv_, f1_};
  auto task64 = make_shared<Task64>(tensor64, cindex);
  return make_shared<FutureTensor>(*Gamma339, task64);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma351_() {
  vector<IndexRange> Gamma351_index = {ci_, active_, active_};
  auto Gamma351 = make_shared<Tensor>(Gamma351_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor65 = {Gamma351, rdm2deriv_, f1_};
  auto task65 = make_shared<Task65>(tensor65, cindex);
  return make_shared<FutureTensor>(*Gamma351, task65);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma362_() {
  vector<IndexRange> Gamma362_index = {ci_, active_, active_, active_, active_};
  auto Gamma362 = make_shared<Tensor>(Gamma362_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor66 = {Gamma362, rdm3deriv_, f1_};
  auto task66 = make_shared<Task66>(tensor66, cindex);
  return make_shared<FutureTensor>(*Gamma362, task66);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma377_() {
  vector<IndexRange> Gamma377_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma377 = make_shared<Tensor>(Gamma377_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor67 = {Gamma377, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task67 = make_shared<Task67>(tensor67, cindex);
  return make_shared<FutureTensor>(*Gamma377, task67);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma395_() {
  vector<IndexRange> Gamma395_index = {ci_, active_, active_, active_, active_, active_, active_};
  auto Gamma395 = make_shared<Tensor>(Gamma395_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor68 = {Gamma395, rdm1deriv_, rdm2deriv_, rdm3deriv_};
  auto task68 = make_shared<Task68>(tensor68, cindex);
  return make_shared<FutureTensor>(*Gamma395, task68);
}

