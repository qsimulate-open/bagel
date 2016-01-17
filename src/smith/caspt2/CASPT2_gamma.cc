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

#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::CASPT2;

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma0_() {
  auto Gamma0 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma0 = Gamma0->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task0 = make_shared<Task0>(array<shared_ptr<Tensor>,6>{{Gamma0, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma0, Gamma0, task0);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma92_() {
  auto Gamma92 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma92 = Gamma92->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task1 = make_shared<Task1>(array<shared_ptr<Tensor>,4>{{Gamma92, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma92, Gamma92, task1);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma2_() {
  auto Gamma2 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma2 = Gamma2->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task2 = make_shared<Task2>(array<shared_ptr<Tensor>,4>{{Gamma2, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma2, Gamma2, task2);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma3_() {
  auto Gamma3 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma3 = Gamma3->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task3 = make_shared<Task3>(array<shared_ptr<Tensor>,4>{{Gamma3, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma3, Gamma3, task3);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma4_() {
  auto Gamma4 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma4 = Gamma4->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task4 = make_shared<Task4>(array<shared_ptr<Tensor>,4>{{Gamma4, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma4, Gamma4, task4);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma5_() {
  auto Gamma5 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma5 = Gamma5->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task5 = make_shared<Task5>(array<shared_ptr<Tensor>,6>{{Gamma5, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma5, Gamma5, task5);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma6_() {
  auto Gamma6 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma6 = Gamma6->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task6 = make_shared<Task6>(array<shared_ptr<Tensor>,4>{{Gamma6, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma6, Gamma6, task6);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma7_() {
  auto Gamma7 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma7 = Gamma7->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task7 = make_shared<Task7>(array<shared_ptr<Tensor>,3>{{Gamma7, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma7, Gamma7, task7);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma9_() {
  auto Gamma9 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma9 = Gamma9->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task8 = make_shared<Task8>(array<shared_ptr<Tensor>,4>{{Gamma9, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma9, Gamma9, task8);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma12_() {
  auto Gamma12 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma12 = Gamma12->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task9 = make_shared<Task9>(array<shared_ptr<Tensor>,3>{{Gamma12, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma12, Gamma12, task9);
}

shared_ptr<FutureTATensor<double,2>> CASPT2::CASPT2::Gamma14_() {
  auto Gamma14 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma14 = Gamma14->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task10 = make_shared<Task10>(array<shared_ptr<Tensor>,5>{{Gamma14, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma14, Gamma14, task10);
}

shared_ptr<FutureTATensor<double,2>> CASPT2::CASPT2::Gamma16_() {
  auto Gamma16 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma16 = Gamma16->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task11 = make_shared<Task11>(array<shared_ptr<Tensor>,3>{{Gamma16, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma16, Gamma16, task11);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma22_() {
  auto Gamma22 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma22 = Gamma22->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task12 = make_shared<Task12>(array<shared_ptr<Tensor>,3>{{Gamma22, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma22, Gamma22, task12);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma28_() {
  auto Gamma28 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma28 = Gamma28->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task13 = make_shared<Task13>(array<shared_ptr<Tensor>,4>{{Gamma28, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma28, Gamma28, task13);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma29_() {
  auto Gamma29 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma29 = Gamma29->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task14 = make_shared<Task14>(array<shared_ptr<Tensor>,3>{{Gamma29, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma29, Gamma29, task14);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma31_() {
  auto Gamma31 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma31 = Gamma31->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task15 = make_shared<Task15>(array<shared_ptr<Tensor>,5>{{Gamma31, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma31, Gamma31, task15);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma32_() {
  auto Gamma32 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma32 = Gamma32->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task16 = make_shared<Task16>(array<shared_ptr<Tensor>,3>{{Gamma32, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma32, Gamma32, task16);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma34_() {
  auto Gamma34 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma34 = Gamma34->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task17 = make_shared<Task17>(array<shared_ptr<Tensor>,5>{{Gamma34, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma34, Gamma34, task17);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma35_() {
  auto Gamma35 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma35 = Gamma35->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task18 = make_shared<Task18>(array<shared_ptr<Tensor>,3>{{Gamma35, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma35, Gamma35, task18);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma37_() {
  auto Gamma37 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma37 = Gamma37->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task19 = make_shared<Task19>(array<shared_ptr<Tensor>,3>{{Gamma37, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma37, Gamma37, task19);
}

shared_ptr<FutureTATensor<double,2>> CASPT2::CASPT2::Gamma38_() {
  auto Gamma38 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma38 = Gamma38->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task20 = make_shared<Task20>(array<shared_ptr<Tensor>,2>{{Gamma38, make_shared<Tensor>(*rdm1_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma38, Gamma38, task20);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma51_() {
  auto Gamma51 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma51 = Gamma51->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task21 = make_shared<Task21>(array<shared_ptr<Tensor>,3>{{Gamma51, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma51, Gamma51, task21);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma56_() {
  auto Gamma56 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma56 = Gamma56->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task22 = make_shared<Task22>(array<shared_ptr<Tensor>,3>{{Gamma56, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma56, Gamma56, task22);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma57_() {
  auto Gamma57 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma57 = Gamma57->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task23 = make_shared<Task23>(array<shared_ptr<Tensor>,3>{{Gamma57, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma57, Gamma57, task23);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma58_() {
  auto Gamma58 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma58 = Gamma58->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task24 = make_shared<Task24>(array<shared_ptr<Tensor>,5>{{Gamma58, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma58, Gamma58, task24);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma59_() {
  auto Gamma59 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma59 = Gamma59->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task25 = make_shared<Task25>(array<shared_ptr<Tensor>,3>{{Gamma59, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma59, Gamma59, task25);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma60_() {
  auto Gamma60 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma60 = Gamma60->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task26 = make_shared<Task26>(array<shared_ptr<Tensor>,2>{{Gamma60, make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma60, Gamma60, task26);
}

shared_ptr<FutureTATensor<double,2>> CASPT2::CASPT2::Gamma79_() {
  auto Gamma79 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma79 = Gamma79->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task27 = make_shared<Task27>(array<shared_ptr<Tensor>,3>{{Gamma79, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma79, Gamma79, task27);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma90_() {
  auto Gamma90 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma90 = Gamma90->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task28 = make_shared<Task28>(array<shared_ptr<Tensor>,3>{{Gamma90, make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma90, Gamma90, task28);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma105_() {
  auto Gamma105 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma105 = Gamma105->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task29 = make_shared<Task29>(array<shared_ptr<Tensor>,4>{{Gamma105, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma105, Gamma105, task29);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma138_() {
  auto Gamma138 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma138 = Gamma138->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task30 = make_shared<Task30>(array<shared_ptr<Tensor>,5>{{Gamma138, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma138, Gamma138, task30);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma169_() {
  auto Gamma169 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma169 = Gamma169->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task31 = make_shared<Task31>(array<shared_ptr<Tensor>,4>{{Gamma169, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma169, Gamma169, task31);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma172_() {
  auto Gamma172 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma172 = Gamma172->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task32 = make_shared<Task32>(array<shared_ptr<Tensor>,4>{{Gamma172, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma172, Gamma172, task32);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma230_() {
  auto Gamma230 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma230 = Gamma230->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task33 = make_shared<Task33>(array<shared_ptr<Tensor>,2>{{Gamma230, make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma230, Gamma230, task33);
}

shared_ptr<FutureTATensor<double,8>> CASPT2::CASPT2::Gamma143_() {
  auto Gamma143 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma143 = Gamma143->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task34 = make_shared<Task34>(array<shared_ptr<Tensor>,5>{{Gamma143, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma143, Gamma143, task34);
}

shared_ptr<FutureTATensor<double,8>> CASPT2::CASPT2::Gamma196_() {
  auto Gamma196 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma196 = Gamma196->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task35 = make_shared<Task35>(array<shared_ptr<Tensor>,4>{{Gamma196, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma196, Gamma196, task35);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma152_() {
  auto Gamma152 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma152 = Gamma152->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task36 = make_shared<Task36>(array<shared_ptr<Tensor>,4>{{Gamma152, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma152, Gamma152, task36);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma248_() {
  auto Gamma248 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma248 = Gamma248->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task37 = make_shared<Task37>(array<shared_ptr<Tensor>,6>{{Gamma248, make_shared<Tensor>(*rdm0deriv_), make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma248, Gamma248, task37);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma249_() {
  auto Gamma249 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma249 = Gamma249->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task38 = make_shared<Task38>(array<shared_ptr<Tensor>,4>{{Gamma249, make_shared<Tensor>(*rdm0deriv_), make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma249, Gamma249, task38);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma250_() {
  auto Gamma250 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma250 = Gamma250->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task39 = make_shared<Task39>(array<shared_ptr<Tensor>,4>{{Gamma250, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma250, Gamma250, task39);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma251_() {
  auto Gamma251 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma251 = Gamma251->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task40 = make_shared<Task40>(array<shared_ptr<Tensor>,4>{{Gamma251, make_shared<Tensor>(*rdm0deriv_), make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma251, Gamma251, task40);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma252_() {
  auto Gamma252 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma252 = Gamma252->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task41 = make_shared<Task41>(array<shared_ptr<Tensor>,4>{{Gamma252, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma252, Gamma252, task41);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma253_() {
  auto Gamma253 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma253 = Gamma253->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task42 = make_shared<Task42>(array<shared_ptr<Tensor>,6>{{Gamma253, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_), make_shared<Tensor>(*rdm4deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma253, Gamma253, task42);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma254_() {
  auto Gamma254 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma254 = Gamma254->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task43 = make_shared<Task43>(array<shared_ptr<Tensor>,4>{{Gamma254, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma254, Gamma254, task43);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma255_() {
  auto Gamma255 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma255 = Gamma255->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task44 = make_shared<Task44>(array<shared_ptr<Tensor>,3>{{Gamma255, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma255, Gamma255, task44);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma257_() {
  auto Gamma257 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma257 = Gamma257->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task45 = make_shared<Task45>(array<shared_ptr<Tensor>,4>{{Gamma257, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma257, Gamma257, task45);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma260_() {
  auto Gamma260 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma260 = Gamma260->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task46 = make_shared<Task46>(array<shared_ptr<Tensor>,3>{{Gamma260, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma260, Gamma260, task46);
}

shared_ptr<FutureTATensor<double,3>> CASPT2::CASPT2::Gamma262_() {
  auto Gamma262 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_});
  auto TAGamma262 = Gamma262->tiledarray<3>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task47 = make_shared<Task47>(array<shared_ptr<Tensor>,5>{{Gamma262, make_shared<Tensor>(*rdm0deriv_), make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,3>>(*TAGamma262, Gamma262, task47);
}

shared_ptr<FutureTATensor<double,3>> CASPT2::CASPT2::Gamma264_() {
  auto Gamma264 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_});
  auto TAGamma264 = Gamma264->tiledarray<3>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task48 = make_shared<Task48>(array<shared_ptr<Tensor>,3>{{Gamma264, make_shared<Tensor>(*rdm0deriv_), make_shared<Tensor>(*rdm1deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,3>>(*TAGamma264, Gamma264, task48);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma270_() {
  auto Gamma270 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma270 = Gamma270->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task49 = make_shared<Task49>(array<shared_ptr<Tensor>,3>{{Gamma270, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma270, Gamma270, task49);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma276_() {
  auto Gamma276 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma276 = Gamma276->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task50 = make_shared<Task50>(array<shared_ptr<Tensor>,4>{{Gamma276, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma276, Gamma276, task50);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma277_() {
  auto Gamma277 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma277 = Gamma277->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task51 = make_shared<Task51>(array<shared_ptr<Tensor>,3>{{Gamma277, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma277, Gamma277, task51);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma279_() {
  auto Gamma279 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma279 = Gamma279->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task52 = make_shared<Task52>(array<shared_ptr<Tensor>,5>{{Gamma279, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma279, Gamma279, task52);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma280_() {
  auto Gamma280 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma280 = Gamma280->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task53 = make_shared<Task53>(array<shared_ptr<Tensor>,3>{{Gamma280, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma280, Gamma280, task53);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma282_() {
  auto Gamma282 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma282 = Gamma282->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task54 = make_shared<Task54>(array<shared_ptr<Tensor>,5>{{Gamma282, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma282, Gamma282, task54);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma283_() {
  auto Gamma283 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma283 = Gamma283->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task55 = make_shared<Task55>(array<shared_ptr<Tensor>,3>{{Gamma283, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma283, Gamma283, task55);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma285_() {
  auto Gamma285 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma285 = Gamma285->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task56 = make_shared<Task56>(array<shared_ptr<Tensor>,3>{{Gamma285, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma285, Gamma285, task56);
}

shared_ptr<FutureTATensor<double,3>> CASPT2::CASPT2::Gamma286_() {
  auto Gamma286 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_});
  auto TAGamma286 = Gamma286->tiledarray<3>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task57 = make_shared<Task57>(array<shared_ptr<Tensor>,2>{{Gamma286, make_shared<Tensor>(*rdm1deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,3>>(*TAGamma286, Gamma286, task57);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma299_() {
  auto Gamma299 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma299 = Gamma299->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task58 = make_shared<Task58>(array<shared_ptr<Tensor>,3>{{Gamma299, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma299, Gamma299, task58);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma304_() {
  auto Gamma304 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma304 = Gamma304->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task59 = make_shared<Task59>(array<shared_ptr<Tensor>,3>{{Gamma304, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma304, Gamma304, task59);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma305_() {
  auto Gamma305 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma305 = Gamma305->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task60 = make_shared<Task60>(array<shared_ptr<Tensor>,3>{{Gamma305, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma305, Gamma305, task60);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma306_() {
  auto Gamma306 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma306 = Gamma306->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task61 = make_shared<Task61>(array<shared_ptr<Tensor>,5>{{Gamma306, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_), make_shared<Tensor>(*rdm4deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma306, Gamma306, task61);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma307_() {
  auto Gamma307 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma307 = Gamma307->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task62 = make_shared<Task62>(array<shared_ptr<Tensor>,3>{{Gamma307, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma307, Gamma307, task62);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma308_() {
  auto Gamma308 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma308 = Gamma308->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task63 = make_shared<Task63>(array<shared_ptr<Tensor>,2>{{Gamma308, make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma308, Gamma308, task63);
}

shared_ptr<FutureTATensor<double,1>> CASPT2::CASPT2::Gamma317_() {
  auto Gamma317 = make_shared<Tensor>(std::vector<IndexRange>{ci_});
  auto TAGamma317 = Gamma317->tiledarray<1>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task64 = make_shared<Task64>(array<shared_ptr<Tensor>,3>{{Gamma317, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,1>>(*TAGamma317, Gamma317, task64);
}

shared_ptr<FutureTATensor<double,3>> CASPT2::CASPT2::Gamma329_() {
  auto Gamma329 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_});
  auto TAGamma329 = Gamma329->tiledarray<3>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task65 = make_shared<Task65>(array<shared_ptr<Tensor>,3>{{Gamma329, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,3>>(*TAGamma329, Gamma329, task65);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma340_() {
  auto Gamma340 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma340 = Gamma340->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task66 = make_shared<Task66>(array<shared_ptr<Tensor>,3>{{Gamma340, make_shared<Tensor>(*rdm3deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma340, Gamma340, task66);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma355_() {
  auto Gamma355 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma355 = Gamma355->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task67 = make_shared<Task67>(array<shared_ptr<Tensor>,4>{{Gamma355, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma355, Gamma355, task67);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma373_() {
  auto Gamma373 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma373 = Gamma373->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task68 = make_shared<Task68>(array<shared_ptr<Tensor>,4>{{Gamma373, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma373, Gamma373, task68);
}

#endif
