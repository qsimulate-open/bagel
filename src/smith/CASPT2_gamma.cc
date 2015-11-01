//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_gamma.cc
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

#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

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

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma105_() {
  auto Gamma105 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma105 = Gamma105->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task9 = make_shared<Task9>(array<shared_ptr<Tensor>,4>{{Gamma105, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma105, Gamma105, task9);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma12_() {
  auto Gamma12 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma12 = Gamma12->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task10 = make_shared<Task10>(array<shared_ptr<Tensor>,3>{{Gamma12, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma12, Gamma12, task10);
}

shared_ptr<FutureTATensor<double,2>> CASPT2::CASPT2::Gamma14_() {
  auto Gamma14 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma14 = Gamma14->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task11 = make_shared<Task11>(array<shared_ptr<Tensor>,5>{{Gamma14, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma14, Gamma14, task11);
}

shared_ptr<FutureTATensor<double,2>> CASPT2::CASPT2::Gamma16_() {
  auto Gamma16 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma16 = Gamma16->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task12 = make_shared<Task12>(array<shared_ptr<Tensor>,3>{{Gamma16, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma16, Gamma16, task12);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma22_() {
  auto Gamma22 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma22 = Gamma22->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task13 = make_shared<Task13>(array<shared_ptr<Tensor>,3>{{Gamma22, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma22, Gamma22, task13);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma28_() {
  auto Gamma28 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma28 = Gamma28->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task14 = make_shared<Task14>(array<shared_ptr<Tensor>,4>{{Gamma28, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma28, Gamma28, task14);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma29_() {
  auto Gamma29 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma29 = Gamma29->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task15 = make_shared<Task15>(array<shared_ptr<Tensor>,3>{{Gamma29, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma29, Gamma29, task15);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma31_() {
  auto Gamma31 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma31 = Gamma31->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task16 = make_shared<Task16>(array<shared_ptr<Tensor>,5>{{Gamma31, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma31, Gamma31, task16);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma32_() {
  auto Gamma32 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma32 = Gamma32->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task17 = make_shared<Task17>(array<shared_ptr<Tensor>,3>{{Gamma32, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma32, Gamma32, task17);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma34_() {
  auto Gamma34 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma34 = Gamma34->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task18 = make_shared<Task18>(array<shared_ptr<Tensor>,5>{{Gamma34, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma34, Gamma34, task18);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma35_() {
  auto Gamma35 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma35 = Gamma35->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task19 = make_shared<Task19>(array<shared_ptr<Tensor>,3>{{Gamma35, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma35, Gamma35, task19);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma37_() {
  auto Gamma37 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma37 = Gamma37->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task20 = make_shared<Task20>(array<shared_ptr<Tensor>,3>{{Gamma37, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma37, Gamma37, task20);
}

shared_ptr<FutureTATensor<double,2>> CASPT2::CASPT2::Gamma38_() {
  auto Gamma38 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma38 = Gamma38->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task21 = make_shared<Task21>(array<shared_ptr<Tensor>,2>{{Gamma38, make_shared<Tensor>(*rdm1_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma38, Gamma38, task21);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma51_() {
  auto Gamma51 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma51 = Gamma51->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task22 = make_shared<Task22>(array<shared_ptr<Tensor>,3>{{Gamma51, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma51, Gamma51, task22);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma56_() {
  auto Gamma56 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma56 = Gamma56->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task23 = make_shared<Task23>(array<shared_ptr<Tensor>,3>{{Gamma56, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma56, Gamma56, task23);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma57_() {
  auto Gamma57 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma57 = Gamma57->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task24 = make_shared<Task24>(array<shared_ptr<Tensor>,3>{{Gamma57, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma57, Gamma57, task24);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma58_() {
  auto Gamma58 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma58 = Gamma58->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task25 = make_shared<Task25>(array<shared_ptr<Tensor>,5>{{Gamma58, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma58, Gamma58, task25);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma59_() {
  auto Gamma59 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma59 = Gamma59->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task26 = make_shared<Task26>(array<shared_ptr<Tensor>,3>{{Gamma59, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma59, Gamma59, task26);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma60_() {
  auto Gamma60 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma60 = Gamma60->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task27 = make_shared<Task27>(array<shared_ptr<Tensor>,2>{{Gamma60, make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma60, Gamma60, task27);
}

shared_ptr<FutureTATensor<double,2>> CASPT2::CASPT2::Gamma79_() {
  auto Gamma79 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma79 = Gamma79->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task28 = make_shared<Task28>(array<shared_ptr<Tensor>,3>{{Gamma79, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma79, Gamma79, task28);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma90_() {
  auto Gamma90 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma90 = Gamma90->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task29 = make_shared<Task29>(array<shared_ptr<Tensor>,3>{{Gamma90, make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*f1_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma90, Gamma90, task29);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma160_() {
  auto Gamma160 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma160 = Gamma160->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task30 = make_shared<Task30>(array<shared_ptr<Tensor>,5>{{Gamma160, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma160, Gamma160, task30);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma191_() {
  auto Gamma191 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma191 = Gamma191->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task31 = make_shared<Task31>(array<shared_ptr<Tensor>,4>{{Gamma191, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma191, Gamma191, task31);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma194_() {
  auto Gamma194 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma194 = Gamma194->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task32 = make_shared<Task32>(array<shared_ptr<Tensor>,4>{{Gamma194, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma194, Gamma194, task32);
}

shared_ptr<FutureTATensor<double,6>> CASPT2::CASPT2::Gamma252_() {
  auto Gamma252 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma252 = Gamma252->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task33 = make_shared<Task33>(array<shared_ptr<Tensor>,2>{{Gamma252, make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma252, Gamma252, task33);
}

shared_ptr<FutureTATensor<double,8>> CASPT2::CASPT2::Gamma165_() {
  auto Gamma165 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma165 = Gamma165->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task34 = make_shared<Task34>(array<shared_ptr<Tensor>,5>{{Gamma165, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma165, Gamma165, task34);
}

shared_ptr<FutureTATensor<double,8>> CASPT2::CASPT2::Gamma218_() {
  auto Gamma218 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma218 = Gamma218->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task35 = make_shared<Task35>(array<shared_ptr<Tensor>,4>{{Gamma218, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma218, Gamma218, task35);
}

shared_ptr<FutureTATensor<double,4>> CASPT2::CASPT2::Gamma174_() {
  auto Gamma174 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma174 = Gamma174->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task36 = make_shared<Task36>(array<shared_ptr<Tensor>,4>{{Gamma174, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma174, Gamma174, task36);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma270_() {
  auto Gamma270 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma270 = Gamma270->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task37 = make_shared<Task37>(array<shared_ptr<Tensor>,6>{{Gamma270, make_shared<Tensor>(*rdm0deriv_), make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma270, Gamma270, task37);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma271_() {
  auto Gamma271 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma271 = Gamma271->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task38 = make_shared<Task38>(array<shared_ptr<Tensor>,4>{{Gamma271, make_shared<Tensor>(*rdm0deriv_), make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma271, Gamma271, task38);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma272_() {
  auto Gamma272 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma272 = Gamma272->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task39 = make_shared<Task39>(array<shared_ptr<Tensor>,4>{{Gamma272, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma272, Gamma272, task39);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma273_() {
  auto Gamma273 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma273 = Gamma273->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task40 = make_shared<Task40>(array<shared_ptr<Tensor>,4>{{Gamma273, make_shared<Tensor>(*rdm0deriv_), make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma273, Gamma273, task40);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma274_() {
  auto Gamma274 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma274 = Gamma274->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task41 = make_shared<Task41>(array<shared_ptr<Tensor>,4>{{Gamma274, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma274, Gamma274, task41);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma275_() {
  auto Gamma275 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma275 = Gamma275->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task42 = make_shared<Task42>(array<shared_ptr<Tensor>,6>{{Gamma275, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_), make_shared<Tensor>(*rdm4deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma275, Gamma275, task42);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma276_() {
  auto Gamma276 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma276 = Gamma276->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task43 = make_shared<Task43>(array<shared_ptr<Tensor>,4>{{Gamma276, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma276, Gamma276, task43);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma277_() {
  auto Gamma277 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma277 = Gamma277->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task44 = make_shared<Task44>(array<shared_ptr<Tensor>,3>{{Gamma277, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma277, Gamma277, task44);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma279_() {
  auto Gamma279 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma279 = Gamma279->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task45 = make_shared<Task45>(array<shared_ptr<Tensor>,4>{{Gamma279, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma279, Gamma279, task45);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma282_() {
  auto Gamma282 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma282 = Gamma282->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task46 = make_shared<Task46>(array<shared_ptr<Tensor>,3>{{Gamma282, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma282, Gamma282, task46);
}

shared_ptr<FutureTATensor<double,3>> CASPT2::CASPT2::Gamma284_() {
  auto Gamma284 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_});
  auto TAGamma284 = Gamma284->tiledarray<3>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task47 = make_shared<Task47>(array<shared_ptr<Tensor>,5>{{Gamma284, make_shared<Tensor>(*rdm0deriv_), make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,3>>(*TAGamma284, Gamma284, task47);
}

shared_ptr<FutureTATensor<double,3>> CASPT2::CASPT2::Gamma286_() {
  auto Gamma286 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_});
  auto TAGamma286 = Gamma286->tiledarray<3>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task48 = make_shared<Task48>(array<shared_ptr<Tensor>,3>{{Gamma286, make_shared<Tensor>(*rdm0deriv_), make_shared<Tensor>(*rdm1deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,3>>(*TAGamma286, Gamma286, task48);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma292_() {
  auto Gamma292 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma292 = Gamma292->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task49 = make_shared<Task49>(array<shared_ptr<Tensor>,3>{{Gamma292, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma292, Gamma292, task49);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma298_() {
  auto Gamma298 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma298 = Gamma298->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task50 = make_shared<Task50>(array<shared_ptr<Tensor>,4>{{Gamma298, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma298, Gamma298, task50);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma299_() {
  auto Gamma299 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma299 = Gamma299->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task51 = make_shared<Task51>(array<shared_ptr<Tensor>,3>{{Gamma299, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma299, Gamma299, task51);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma301_() {
  auto Gamma301 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma301 = Gamma301->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task52 = make_shared<Task52>(array<shared_ptr<Tensor>,5>{{Gamma301, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma301, Gamma301, task52);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma302_() {
  auto Gamma302 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma302 = Gamma302->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task53 = make_shared<Task53>(array<shared_ptr<Tensor>,3>{{Gamma302, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma302, Gamma302, task53);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma304_() {
  auto Gamma304 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma304 = Gamma304->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task54 = make_shared<Task54>(array<shared_ptr<Tensor>,5>{{Gamma304, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma304, Gamma304, task54);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma305_() {
  auto Gamma305 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma305 = Gamma305->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task55 = make_shared<Task55>(array<shared_ptr<Tensor>,3>{{Gamma305, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma305, Gamma305, task55);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma307_() {
  auto Gamma307 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma307 = Gamma307->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task56 = make_shared<Task56>(array<shared_ptr<Tensor>,3>{{Gamma307, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma307, Gamma307, task56);
}

shared_ptr<FutureTATensor<double,3>> CASPT2::CASPT2::Gamma308_() {
  auto Gamma308 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_});
  auto TAGamma308 = Gamma308->tiledarray<3>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task57 = make_shared<Task57>(array<shared_ptr<Tensor>,2>{{Gamma308, make_shared<Tensor>(*rdm1deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,3>>(*TAGamma308, Gamma308, task57);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma321_() {
  auto Gamma321 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma321 = Gamma321->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task58 = make_shared<Task58>(array<shared_ptr<Tensor>,3>{{Gamma321, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma321, Gamma321, task58);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma326_() {
  auto Gamma326 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma326 = Gamma326->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task59 = make_shared<Task59>(array<shared_ptr<Tensor>,3>{{Gamma326, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma326, Gamma326, task59);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma327_() {
  auto Gamma327 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma327 = Gamma327->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task60 = make_shared<Task60>(array<shared_ptr<Tensor>,3>{{Gamma327, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma327, Gamma327, task60);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma328_() {
  auto Gamma328 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma328 = Gamma328->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task61 = make_shared<Task61>(array<shared_ptr<Tensor>,5>{{Gamma328, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_), make_shared<Tensor>(*rdm4deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma328, Gamma328, task61);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma329_() {
  auto Gamma329 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma329 = Gamma329->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task62 = make_shared<Task62>(array<shared_ptr<Tensor>,3>{{Gamma329, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma329, Gamma329, task62);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma330_() {
  auto Gamma330 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma330 = Gamma330->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task63 = make_shared<Task63>(array<shared_ptr<Tensor>,2>{{Gamma330, make_shared<Tensor>(*rdm2deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma330, Gamma330, task63);
}

shared_ptr<FutureTATensor<double,1>> CASPT2::CASPT2::Gamma339_() {
  auto Gamma339 = make_shared<Tensor>(std::vector<IndexRange>{ci_});
  auto TAGamma339 = Gamma339->tiledarray<1>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task64 = make_shared<Task64>(array<shared_ptr<Tensor>,3>{{Gamma339, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,1>>(*TAGamma339, Gamma339, task64);
}

shared_ptr<FutureTATensor<double,3>> CASPT2::CASPT2::Gamma351_() {
  auto Gamma351 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_});
  auto TAGamma351 = Gamma351->tiledarray<3>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task65 = make_shared<Task65>(array<shared_ptr<Tensor>,3>{{Gamma351, make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,3>>(*TAGamma351, Gamma351, task65);
}

shared_ptr<FutureTATensor<double,5>> CASPT2::CASPT2::Gamma362_() {
  auto Gamma362 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_});
  auto TAGamma362 = Gamma362->tiledarray<5>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task66 = make_shared<Task66>(array<shared_ptr<Tensor>,3>{{Gamma362, make_shared<Tensor>(*rdm3deriv_), make_shared<Tensor>(*f1_)}}, cindex);
  return make_shared<FutureTATensor<double,5>>(*TAGamma362, Gamma362, task66);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma377_() {
  auto Gamma377 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma377 = Gamma377->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task67 = make_shared<Task67>(array<shared_ptr<Tensor>,4>{{Gamma377, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma377, Gamma377, task67);
}

shared_ptr<FutureTATensor<double,7>> CASPT2::CASPT2::Gamma395_() {
  auto Gamma395 = make_shared<Tensor>(std::vector<IndexRange>{ci_, active_, active_, active_, active_, active_, active_});
  auto TAGamma395 = Gamma395->tiledarray<7>();
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  auto task68 = make_shared<Task68>(array<shared_ptr<Tensor>,4>{{Gamma395, make_shared<Tensor>(*rdm1deriv_), make_shared<Tensor>(*rdm2deriv_), make_shared<Tensor>(*rdm3deriv_)}}, cindex);
  return make_shared<FutureTATensor<double,7>>(*TAGamma395, Gamma395, task68);
}

#endif
