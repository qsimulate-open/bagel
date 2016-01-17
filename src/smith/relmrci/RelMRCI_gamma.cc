//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: RelMRCI_gamma.cc
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
#include <src/smith/relmrci/RelMRCI.h>
#include <src/smith/relmrci/RelMRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma0_() {
  auto Gamma0 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma0 = Gamma0->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task0 = make_shared<Task0>(array<shared_ptr<Tensor>,4>{{Gamma0, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma0, Gamma0, task0);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma1_() {
  auto Gamma1 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma1 = Gamma1->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task1 = make_shared<Task1>(array<shared_ptr<Tensor>,4>{{Gamma1, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma1, Gamma1, task1);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma2_() {
  auto Gamma2 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma2 = Gamma2->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task2 = make_shared<Task2>(array<shared_ptr<Tensor>,4>{{Gamma2, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma2, Gamma2, task2);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma58_() {
  auto Gamma58 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma58 = Gamma58->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task3 = make_shared<Task3>(array<shared_ptr<Tensor>,5>{{Gamma58, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma58, Gamma58, task3);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma59_() {
  auto Gamma59 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma59 = Gamma59->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task4 = make_shared<Task4>(array<shared_ptr<Tensor>,5>{{Gamma59, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma59, Gamma59, task4);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma60_() {
  auto Gamma60 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma60 = Gamma60->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task5 = make_shared<Task5>(array<shared_ptr<Tensor>,5>{{Gamma60, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma60, Gamma60, task5);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma63_() {
  auto Gamma63 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma63 = Gamma63->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task6 = make_shared<Task6>(array<shared_ptr<Tensor>,5>{{Gamma63, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma63, Gamma63, task6);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma64_() {
  auto Gamma64 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma64 = Gamma64->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task7 = make_shared<Task7>(array<shared_ptr<Tensor>,5>{{Gamma64, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma64, Gamma64, task7);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma65_() {
  auto Gamma65 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma65 = Gamma65->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task8 = make_shared<Task8>(array<shared_ptr<Tensor>,4>{{Gamma65, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma65, Gamma65, task8);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma66_() {
  auto Gamma66 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma66 = Gamma66->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task9 = make_shared<Task9>(array<shared_ptr<Tensor>,5>{{Gamma66, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma66, Gamma66, task9);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma67_() {
  auto Gamma67 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma67 = Gamma67->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task10 = make_shared<Task10>(array<shared_ptr<Tensor>,5>{{Gamma67, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma67, Gamma67, task10);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma3_() {
  auto Gamma3 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma3 = Gamma3->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task11 = make_shared<Task11>(array<shared_ptr<Tensor>,4>{{Gamma3, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma3, Gamma3, task11);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma4_() {
  auto Gamma4 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma4 = Gamma4->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task12 = make_shared<Task12>(array<shared_ptr<Tensor>,4>{{Gamma4, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma4, Gamma4, task12);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma5_() {
  auto Gamma5 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma5 = Gamma5->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task13 = make_shared<Task13>(array<shared_ptr<Tensor>,3>{{Gamma5, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma5, Gamma5, task13);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma74_() {
  auto Gamma74 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma74 = Gamma74->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task14 = make_shared<Task14>(array<shared_ptr<Tensor>,5>{{Gamma74, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma74, Gamma74, task14);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma75_() {
  auto Gamma75 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma75 = Gamma75->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task15 = make_shared<Task15>(array<shared_ptr<Tensor>,5>{{Gamma75, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma75, Gamma75, task15);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma77_() {
  auto Gamma77 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma77 = Gamma77->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task16 = make_shared<Task16>(array<shared_ptr<Tensor>,5>{{Gamma77, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma77, Gamma77, task16);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma78_() {
  auto Gamma78 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma78 = Gamma78->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task17 = make_shared<Task17>(array<shared_ptr<Tensor>,5>{{Gamma78, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma78, Gamma78, task17);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma79_() {
  auto Gamma79 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma79 = Gamma79->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task18 = make_shared<Task18>(array<shared_ptr<Tensor>,5>{{Gamma79, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma79, Gamma79, task18);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma81_() {
  auto Gamma81 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma81 = Gamma81->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task19 = make_shared<Task19>(array<shared_ptr<Tensor>,4>{{Gamma81, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma81, Gamma81, task19);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma84_() {
  auto Gamma84 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma84 = Gamma84->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task20 = make_shared<Task20>(array<shared_ptr<Tensor>,4>{{Gamma84, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma84, Gamma84, task20);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma86_() {
  auto Gamma86 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma86 = Gamma86->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task21 = make_shared<Task21>(array<shared_ptr<Tensor>,4>{{Gamma86, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma86, Gamma86, task21);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma92_() {
  auto Gamma92 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma92 = Gamma92->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task22 = make_shared<Task22>(array<shared_ptr<Tensor>,5>{{Gamma92, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma92, Gamma92, task22);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma95_() {
  auto Gamma95 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma95 = Gamma95->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task23 = make_shared<Task23>(array<shared_ptr<Tensor>,4>{{Gamma95, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma95, Gamma95, task23);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma98_() {
  auto Gamma98 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma98 = Gamma98->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task24 = make_shared<Task24>(array<shared_ptr<Tensor>,4>{{Gamma98, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma98, Gamma98, task24);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma411_() {
  auto Gamma411 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma411 = Gamma411->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task25 = make_shared<Task25>(array<shared_ptr<Tensor>,6>{{Gamma411, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma411, Gamma411, task25);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma412_() {
  auto Gamma412 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma412 = Gamma412->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task26 = make_shared<Task26>(array<shared_ptr<Tensor>,6>{{Gamma412, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma412, Gamma412, task26);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma9_() {
  auto Gamma9 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma9 = Gamma9->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task27 = make_shared<Task27>(array<shared_ptr<Tensor>,3>{{Gamma9, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma9, Gamma9, task27);
}

shared_ptr<FutureTATensor<std::complex<double>,2>> RelMRCI::RelMRCI::Gamma11_() {
  auto Gamma11 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma11 = Gamma11->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task28 = make_shared<Task28>(array<shared_ptr<Tensor>,3>{{Gamma11, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,2>>(*TAGamma11, Gamma11, task28);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma160_() {
  auto Gamma160 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma160 = Gamma160->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task29 = make_shared<Task29>(array<shared_ptr<Tensor>,4>{{Gamma160, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma160, Gamma160, task29);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma99_() {
  auto Gamma99 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma99 = Gamma99->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task30 = make_shared<Task30>(array<shared_ptr<Tensor>,5>{{Gamma99, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma99, Gamma99, task30);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma105_() {
  auto Gamma105 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma105 = Gamma105->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task31 = make_shared<Task31>(array<shared_ptr<Tensor>,4>{{Gamma105, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma105, Gamma105, task31);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma110_() {
  auto Gamma110 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma110 = Gamma110->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task32 = make_shared<Task32>(array<shared_ptr<Tensor>,4>{{Gamma110, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma110, Gamma110, task32);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma128_() {
  auto Gamma128 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma128 = Gamma128->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task33 = make_shared<Task33>(array<shared_ptr<Tensor>,4>{{Gamma128, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma128, Gamma128, task33);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma151_() {
  auto Gamma151 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma151 = Gamma151->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task34 = make_shared<Task34>(array<shared_ptr<Tensor>,4>{{Gamma151, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma151, Gamma151, task34);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma159_() {
  auto Gamma159 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma159 = Gamma159->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task35 = make_shared<Task35>(array<shared_ptr<Tensor>,3>{{Gamma159, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma159, Gamma159, task35);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma174_() {
  auto Gamma174 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma174 = Gamma174->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task36 = make_shared<Task36>(array<shared_ptr<Tensor>,3>{{Gamma174, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma174, Gamma174, task36);
}

shared_ptr<FutureTATensor<std::complex<double>,2>> RelMRCI::RelMRCI::Gamma413_() {
  auto Gamma413 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma413 = Gamma413->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task37 = make_shared<Task37>(array<shared_ptr<Tensor>,5>{{Gamma413, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,2>>(*TAGamma413, Gamma413, task37);
}

shared_ptr<FutureTATensor<std::complex<double>,2>> RelMRCI::RelMRCI::Gamma415_() {
  auto Gamma415 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma415 = Gamma415->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task38 = make_shared<Task38>(array<shared_ptr<Tensor>,6>{{Gamma415, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,2>>(*TAGamma415, Gamma415, task38);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma24_() {
  auto Gamma24 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma24 = Gamma24->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task39 = make_shared<Task39>(array<shared_ptr<Tensor>,3>{{Gamma24, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma24, Gamma24, task39);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma26_() {
  auto Gamma26 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma26 = Gamma26->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task40 = make_shared<Task40>(array<shared_ptr<Tensor>,3>{{Gamma26, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma26, Gamma26, task40);
}

shared_ptr<FutureTATensor<std::complex<double>,2>> RelMRCI::RelMRCI::Gamma27_() {
  auto Gamma27 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma27 = Gamma27->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task41 = make_shared<Task41>(array<shared_ptr<Tensor>,2>{{Gamma27, make_shared<Tensor>(*rdm1_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,2>>(*TAGamma27, Gamma27, task41);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma179_() {
  auto Gamma179 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma179 = Gamma179->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task42 = make_shared<Task42>(array<shared_ptr<Tensor>,5>{{Gamma179, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma179, Gamma179, task42);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma183_() {
  auto Gamma183 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma183 = Gamma183->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task43 = make_shared<Task43>(array<shared_ptr<Tensor>,4>{{Gamma183, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma183, Gamma183, task43);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma193_() {
  auto Gamma193 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma193 = Gamma193->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task44 = make_shared<Task44>(array<shared_ptr<Tensor>,4>{{Gamma193, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma193, Gamma193, task44);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma203_() {
  auto Gamma203 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma203 = Gamma203->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task45 = make_shared<Task45>(array<shared_ptr<Tensor>,4>{{Gamma203, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma203, Gamma203, task45);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma204_() {
  auto Gamma204 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma204 = Gamma204->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task46 = make_shared<Task46>(array<shared_ptr<Tensor>,4>{{Gamma204, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma204, Gamma204, task46);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma229_() {
  auto Gamma229 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma229 = Gamma229->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task47 = make_shared<Task47>(array<shared_ptr<Tensor>,3>{{Gamma229, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma229, Gamma229, task47);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma417_() {
  auto Gamma417 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma417 = Gamma417->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task48 = make_shared<Task48>(array<shared_ptr<Tensor>,5>{{Gamma417, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma417, Gamma417, task48);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma418_() {
  auto Gamma418 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma418 = Gamma418->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task49 = make_shared<Task49>(array<shared_ptr<Tensor>,6>{{Gamma418, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma418, Gamma418, task49);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma31_() {
  auto Gamma31 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma31 = Gamma31->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task50 = make_shared<Task50>(array<shared_ptr<Tensor>,3>{{Gamma31, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma31, Gamma31, task50);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma32_() {
  auto Gamma32 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma32 = Gamma32->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task51 = make_shared<Task51>(array<shared_ptr<Tensor>,3>{{Gamma32, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma32, Gamma32, task51);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma33_() {
  auto Gamma33 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma33 = Gamma33->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task52 = make_shared<Task52>(array<shared_ptr<Tensor>,2>{{Gamma33, make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma33, Gamma33, task52);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma230_() {
  auto Gamma230 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma230 = Gamma230->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task53 = make_shared<Task53>(array<shared_ptr<Tensor>,4>{{Gamma230, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma230, Gamma230, task53);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma231_() {
  auto Gamma231 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma231 = Gamma231->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task54 = make_shared<Task54>(array<shared_ptr<Tensor>,3>{{Gamma231, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma231, Gamma231, task54);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma232_() {
  auto Gamma232 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma232 = Gamma232->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task55 = make_shared<Task55>(array<shared_ptr<Tensor>,4>{{Gamma232, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma232, Gamma232, task55);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma233_() {
  auto Gamma233 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma233 = Gamma233->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task56 = make_shared<Task56>(array<shared_ptr<Tensor>,4>{{Gamma233, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma233, Gamma233, task56);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma236_() {
  auto Gamma236 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma236 = Gamma236->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task57 = make_shared<Task57>(array<shared_ptr<Tensor>,4>{{Gamma236, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma236, Gamma236, task57);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma237_() {
  auto Gamma237 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma237 = Gamma237->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task58 = make_shared<Task58>(array<shared_ptr<Tensor>,4>{{Gamma237, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma237, Gamma237, task58);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma238_() {
  auto Gamma238 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma238 = Gamma238->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task59 = make_shared<Task59>(array<shared_ptr<Tensor>,4>{{Gamma238, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma238, Gamma238, task59);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma245_() {
  auto Gamma245 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma245 = Gamma245->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task60 = make_shared<Task60>(array<shared_ptr<Tensor>,3>{{Gamma245, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma245, Gamma245, task60);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma246_() {
  auto Gamma246 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma246 = Gamma246->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task61 = make_shared<Task61>(array<shared_ptr<Tensor>,3>{{Gamma246, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma246, Gamma246, task61);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma253_() {
  auto Gamma253 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma253 = Gamma253->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task62 = make_shared<Task62>(array<shared_ptr<Tensor>,3>{{Gamma253, make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma253, Gamma253, task62);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma419_() {
  auto Gamma419 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma419 = Gamma419->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task63 = make_shared<Task63>(array<shared_ptr<Tensor>,5>{{Gamma419, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma419, Gamma419, task63);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma420_() {
  auto Gamma420 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma420 = Gamma420->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task64 = make_shared<Task64>(array<shared_ptr<Tensor>,5>{{Gamma420, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma420, Gamma420, task64);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma317_() {
  auto Gamma317 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma317 = Gamma317->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task65 = make_shared<Task65>(array<shared_ptr<Tensor>,4>{{Gamma317, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma317, Gamma317, task65);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma318_() {
  auto Gamma318 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma318 = Gamma318->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task66 = make_shared<Task66>(array<shared_ptr<Tensor>,3>{{Gamma318, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma318, Gamma318, task66);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma335_() {
  auto Gamma335 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma335 = Gamma335->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task67 = make_shared<Task67>(array<shared_ptr<Tensor>,3>{{Gamma335, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma335, Gamma335, task67);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma336_() {
  auto Gamma336 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma336 = Gamma336->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task68 = make_shared<Task68>(array<shared_ptr<Tensor>,3>{{Gamma336, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma336, Gamma336, task68);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma368_() {
  auto Gamma368 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma368 = Gamma368->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task69 = make_shared<Task69>(array<shared_ptr<Tensor>,2>{{Gamma368, make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma368, Gamma368, task69);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma363_() {
  auto Gamma363 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma363 = Gamma363->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task70 = make_shared<Task70>(array<shared_ptr<Tensor>,3>{{Gamma363, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma363, Gamma363, task70);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma389_() {
  auto Gamma389 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma389 = Gamma389->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task71 = make_shared<Task71>(array<shared_ptr<Tensor>,2>{{Gamma389, make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma389, Gamma389, task71);
}

shared_ptr<FutureTATensor<std::complex<double>,2>> RelMRCI::RelMRCI::Gamma425_() {
  auto Gamma425 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma425 = Gamma425->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task72 = make_shared<Task72>(array<shared_ptr<Tensor>,4>{{Gamma425, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,2>>(*TAGamma425, Gamma425, task72);
}

shared_ptr<FutureTATensor<std::complex<double>,2>> RelMRCI::RelMRCI::Gamma427_() {
  auto Gamma427 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma427 = Gamma427->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task73 = make_shared<Task73>(array<shared_ptr<Tensor>,5>{{Gamma427, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,2>>(*TAGamma427, Gamma427, task73);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma394_() {
  auto Gamma394 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma394 = Gamma394->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task74 = make_shared<Task74>(array<shared_ptr<Tensor>,3>{{Gamma394, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma394, Gamma394, task74);
}

shared_ptr<FutureTATensor<std::complex<double>,8>> RelMRCI::RelMRCI::Gamma395_() {
  auto Gamma395 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma395 = Gamma395->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task75 = make_shared<Task75>(array<shared_ptr<Tensor>,3>{{Gamma395, make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,8>>(*TAGamma395, Gamma395, task75);
}

shared_ptr<FutureTATensor<std::complex<double>,6>> RelMRCI::RelMRCI::Gamma407_() {
  auto Gamma407 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma407 = Gamma407->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task76 = make_shared<Task76>(array<shared_ptr<Tensor>,2>{{Gamma407, make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,6>>(*TAGamma407, Gamma407, task76);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma409_() {
  auto Gamma409 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma409 = Gamma409->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task77 = make_shared<Task77>(array<shared_ptr<Tensor>,6>{{Gamma409, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma409, Gamma409, task77);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma410_() {
  auto Gamma410 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma410 = Gamma410->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task78 = make_shared<Task78>(array<shared_ptr<Tensor>,7>{{Gamma410, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma410, Gamma410, task78);
}

shared_ptr<FutureTATensor<std::complex<double>,0>> RelMRCI::RelMRCI::Gamma421_() {
  auto Gamma421 = make_shared<Tensor>(std::vector<IndexRange>{});
  auto TAGamma421 = Gamma421->tiledarray<0>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task79 = make_shared<Task79>(array<shared_ptr<Tensor>,3>{{Gamma421, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,0>>(*TAGamma421, Gamma421, task79);
}

shared_ptr<FutureTATensor<std::complex<double>,0>> RelMRCI::RelMRCI::Gamma423_() {
  auto Gamma423 = make_shared<Tensor>(std::vector<IndexRange>{});
  auto TAGamma423 = Gamma423->tiledarray<0>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task80 = make_shared<Task80>(array<shared_ptr<Tensor>,4>{{Gamma423, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,0>>(*TAGamma423, Gamma423, task80);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma429_() {
  auto Gamma429 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma429 = Gamma429->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task81 = make_shared<Task81>(array<shared_ptr<Tensor>,4>{{Gamma429, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma429, Gamma429, task81);
}

shared_ptr<FutureTATensor<std::complex<double>,4>> RelMRCI::RelMRCI::Gamma430_() {
  auto Gamma430 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma430 = Gamma430->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task82 = make_shared<Task82>(array<shared_ptr<Tensor>,5>{{Gamma430, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<std::complex<double>,4>>(*TAGamma430, Gamma430, task82);
}

#endif
