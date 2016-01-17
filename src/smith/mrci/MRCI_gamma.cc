//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MRCI_gamma.cc
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

#include <src/smith/mrci/MRCI.h>
#include <src/smith/mrci/MRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::MRCI;

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma0_() {
  auto Gamma0 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma0 = Gamma0->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task0 = make_shared<Task0>(array<shared_ptr<Tensor>,4>{{Gamma0, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma0, Gamma0, task0);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma1_() {
  auto Gamma1 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma1 = Gamma1->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task1 = make_shared<Task1>(array<shared_ptr<Tensor>,4>{{Gamma1, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma1, Gamma1, task1);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma2_() {
  auto Gamma2 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma2 = Gamma2->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task2 = make_shared<Task2>(array<shared_ptr<Tensor>,4>{{Gamma2, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma2, Gamma2, task2);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma80_() {
  auto Gamma80 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma80 = Gamma80->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task3 = make_shared<Task3>(array<shared_ptr<Tensor>,5>{{Gamma80, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma80, Gamma80, task3);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma81_() {
  auto Gamma81 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma81 = Gamma81->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task4 = make_shared<Task4>(array<shared_ptr<Tensor>,5>{{Gamma81, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma81, Gamma81, task4);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma82_() {
  auto Gamma82 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma82 = Gamma82->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task5 = make_shared<Task5>(array<shared_ptr<Tensor>,5>{{Gamma82, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma82, Gamma82, task5);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma85_() {
  auto Gamma85 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma85 = Gamma85->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task6 = make_shared<Task6>(array<shared_ptr<Tensor>,5>{{Gamma85, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma85, Gamma85, task6);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma86_() {
  auto Gamma86 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma86 = Gamma86->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task7 = make_shared<Task7>(array<shared_ptr<Tensor>,5>{{Gamma86, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma86, Gamma86, task7);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma87_() {
  auto Gamma87 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma87 = Gamma87->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task8 = make_shared<Task8>(array<shared_ptr<Tensor>,4>{{Gamma87, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma87, Gamma87, task8);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma88_() {
  auto Gamma88 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma88 = Gamma88->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task9 = make_shared<Task9>(array<shared_ptr<Tensor>,5>{{Gamma88, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma88, Gamma88, task9);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma89_() {
  auto Gamma89 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma89 = Gamma89->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task10 = make_shared<Task10>(array<shared_ptr<Tensor>,5>{{Gamma89, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma89, Gamma89, task10);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma94_() {
  auto Gamma94 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma94 = Gamma94->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task11 = make_shared<Task11>(array<shared_ptr<Tensor>,4>{{Gamma94, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma94, Gamma94, task11);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma3_() {
  auto Gamma3 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma3 = Gamma3->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task12 = make_shared<Task12>(array<shared_ptr<Tensor>,4>{{Gamma3, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma3, Gamma3, task12);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma4_() {
  auto Gamma4 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma4 = Gamma4->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task13 = make_shared<Task13>(array<shared_ptr<Tensor>,4>{{Gamma4, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma4, Gamma4, task13);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma5_() {
  auto Gamma5 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma5 = Gamma5->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task14 = make_shared<Task14>(array<shared_ptr<Tensor>,3>{{Gamma5, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma5, Gamma5, task14);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma7_() {
  auto Gamma7 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma7 = Gamma7->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task15 = make_shared<Task15>(array<shared_ptr<Tensor>,4>{{Gamma7, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma7, Gamma7, task15);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma97_() {
  auto Gamma97 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma97 = Gamma97->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task16 = make_shared<Task16>(array<shared_ptr<Tensor>,5>{{Gamma97, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma97, Gamma97, task16);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma98_() {
  auto Gamma98 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma98 = Gamma98->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task17 = make_shared<Task17>(array<shared_ptr<Tensor>,5>{{Gamma98, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma98, Gamma98, task17);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma100_() {
  auto Gamma100 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma100 = Gamma100->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task18 = make_shared<Task18>(array<shared_ptr<Tensor>,5>{{Gamma100, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma100, Gamma100, task18);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma101_() {
  auto Gamma101 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma101 = Gamma101->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task19 = make_shared<Task19>(array<shared_ptr<Tensor>,5>{{Gamma101, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma101, Gamma101, task19);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma102_() {
  auto Gamma102 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma102 = Gamma102->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task20 = make_shared<Task20>(array<shared_ptr<Tensor>,5>{{Gamma102, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma102, Gamma102, task20);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma104_() {
  auto Gamma104 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma104 = Gamma104->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task21 = make_shared<Task21>(array<shared_ptr<Tensor>,4>{{Gamma104, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma104, Gamma104, task21);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma107_() {
  auto Gamma107 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma107 = Gamma107->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task22 = make_shared<Task22>(array<shared_ptr<Tensor>,4>{{Gamma107, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma107, Gamma107, task22);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma109_() {
  auto Gamma109 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma109 = Gamma109->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task23 = make_shared<Task23>(array<shared_ptr<Tensor>,4>{{Gamma109, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma109, Gamma109, task23);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma114_() {
  auto Gamma114 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma114 = Gamma114->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task24 = make_shared<Task24>(array<shared_ptr<Tensor>,5>{{Gamma114, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma114, Gamma114, task24);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma115_() {
  auto Gamma115 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma115 = Gamma115->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task25 = make_shared<Task25>(array<shared_ptr<Tensor>,5>{{Gamma115, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma115, Gamma115, task25);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma119_() {
  auto Gamma119 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma119 = Gamma119->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task26 = make_shared<Task26>(array<shared_ptr<Tensor>,5>{{Gamma119, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma119, Gamma119, task26);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma122_() {
  auto Gamma122 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma122 = Gamma122->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task27 = make_shared<Task27>(array<shared_ptr<Tensor>,4>{{Gamma122, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma122, Gamma122, task27);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma547_() {
  auto Gamma547 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma547 = Gamma547->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task28 = make_shared<Task28>(array<shared_ptr<Tensor>,6>{{Gamma547, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma547, Gamma547, task28);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma548_() {
  auto Gamma548 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma548 = Gamma548->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task29 = make_shared<Task29>(array<shared_ptr<Tensor>,6>{{Gamma548, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma548, Gamma548, task29);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma10_() {
  auto Gamma10 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma10 = Gamma10->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task30 = make_shared<Task30>(array<shared_ptr<Tensor>,3>{{Gamma10, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma10, Gamma10, task30);
}

shared_ptr<FutureTATensor<double,2>> MRCI::MRCI::Gamma12_() {
  auto Gamma12 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma12 = Gamma12->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task31 = make_shared<Task31>(array<shared_ptr<Tensor>,3>{{Gamma12, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma12, Gamma12, task31);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma18_() {
  auto Gamma18 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma18 = Gamma18->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task32 = make_shared<Task32>(array<shared_ptr<Tensor>,3>{{Gamma18, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma18, Gamma18, task32);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma197_() {
  auto Gamma197 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma197 = Gamma197->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task33 = make_shared<Task33>(array<shared_ptr<Tensor>,4>{{Gamma197, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma197, Gamma197, task33);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma126_() {
  auto Gamma126 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma126 = Gamma126->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task34 = make_shared<Task34>(array<shared_ptr<Tensor>,5>{{Gamma126, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma126, Gamma126, task34);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma132_() {
  auto Gamma132 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma132 = Gamma132->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task35 = make_shared<Task35>(array<shared_ptr<Tensor>,4>{{Gamma132, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma132, Gamma132, task35);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma137_() {
  auto Gamma137 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma137 = Gamma137->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task36 = make_shared<Task36>(array<shared_ptr<Tensor>,4>{{Gamma137, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma137, Gamma137, task36);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma155_() {
  auto Gamma155 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma155 = Gamma155->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task37 = make_shared<Task37>(array<shared_ptr<Tensor>,4>{{Gamma155, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma155, Gamma155, task37);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma176_() {
  auto Gamma176 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma176 = Gamma176->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task38 = make_shared<Task38>(array<shared_ptr<Tensor>,4>{{Gamma176, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma176, Gamma176, task38);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma178_() {
  auto Gamma178 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma178 = Gamma178->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task39 = make_shared<Task39>(array<shared_ptr<Tensor>,4>{{Gamma178, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma178, Gamma178, task39);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma179_() {
  auto Gamma179 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma179 = Gamma179->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task40 = make_shared<Task40>(array<shared_ptr<Tensor>,4>{{Gamma179, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma179, Gamma179, task40);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma196_() {
  auto Gamma196 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma196 = Gamma196->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task41 = make_shared<Task41>(array<shared_ptr<Tensor>,3>{{Gamma196, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma196, Gamma196, task41);
}

shared_ptr<FutureTATensor<double,2>> MRCI::MRCI::Gamma549_() {
  auto Gamma549 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma549 = Gamma549->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task42 = make_shared<Task42>(array<shared_ptr<Tensor>,5>{{Gamma549, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma549, Gamma549, task42);
}

shared_ptr<FutureTATensor<double,2>> MRCI::MRCI::Gamma551_() {
  auto Gamma551 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma551 = Gamma551->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task43 = make_shared<Task43>(array<shared_ptr<Tensor>,6>{{Gamma551, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma551, Gamma551, task43);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma24_() {
  auto Gamma24 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma24 = Gamma24->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task44 = make_shared<Task44>(array<shared_ptr<Tensor>,4>{{Gamma24, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma24, Gamma24, task44);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma25_() {
  auto Gamma25 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma25 = Gamma25->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task45 = make_shared<Task45>(array<shared_ptr<Tensor>,3>{{Gamma25, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma25, Gamma25, task45);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma27_() {
  auto Gamma27 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma27 = Gamma27->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task46 = make_shared<Task46>(array<shared_ptr<Tensor>,3>{{Gamma27, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma27, Gamma27, task46);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma29_() {
  auto Gamma29 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma29 = Gamma29->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task47 = make_shared<Task47>(array<shared_ptr<Tensor>,3>{{Gamma29, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma29, Gamma29, task47);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma31_() {
  auto Gamma31 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma31 = Gamma31->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task48 = make_shared<Task48>(array<shared_ptr<Tensor>,3>{{Gamma31, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma31, Gamma31, task48);
}

shared_ptr<FutureTATensor<double,2>> MRCI::MRCI::Gamma32_() {
  auto Gamma32 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma32 = Gamma32->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task49 = make_shared<Task49>(array<shared_ptr<Tensor>,2>{{Gamma32, make_shared<Tensor>(*rdm1_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma32, Gamma32, task49);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma215_() {
  auto Gamma215 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma215 = Gamma215->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task50 = make_shared<Task50>(array<shared_ptr<Tensor>,4>{{Gamma215, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma215, Gamma215, task50);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma216_() {
  auto Gamma216 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma216 = Gamma216->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task51 = make_shared<Task51>(array<shared_ptr<Tensor>,5>{{Gamma216, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma216, Gamma216, task51);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma217_() {
  auto Gamma217 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma217 = Gamma217->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task52 = make_shared<Task52>(array<shared_ptr<Tensor>,5>{{Gamma217, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma217, Gamma217, task52);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma220_() {
  auto Gamma220 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma220 = Gamma220->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task53 = make_shared<Task53>(array<shared_ptr<Tensor>,4>{{Gamma220, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma220, Gamma220, task53);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma222_() {
  auto Gamma222 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma222 = Gamma222->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task54 = make_shared<Task54>(array<shared_ptr<Tensor>,4>{{Gamma222, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma222, Gamma222, task54);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma221_() {
  auto Gamma221 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma221 = Gamma221->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task55 = make_shared<Task55>(array<shared_ptr<Tensor>,4>{{Gamma221, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma221, Gamma221, task55);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma230_() {
  auto Gamma230 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma230 = Gamma230->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task56 = make_shared<Task56>(array<shared_ptr<Tensor>,4>{{Gamma230, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma230, Gamma230, task56);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma232_() {
  auto Gamma232 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma232 = Gamma232->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task57 = make_shared<Task57>(array<shared_ptr<Tensor>,4>{{Gamma232, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma232, Gamma232, task57);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma234_() {
  auto Gamma234 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma234 = Gamma234->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task58 = make_shared<Task58>(array<shared_ptr<Tensor>,4>{{Gamma234, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma234, Gamma234, task58);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma233_() {
  auto Gamma233 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma233 = Gamma233->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task59 = make_shared<Task59>(array<shared_ptr<Tensor>,4>{{Gamma233, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma233, Gamma233, task59);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma235_() {
  auto Gamma235 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma235 = Gamma235->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task60 = make_shared<Task60>(array<shared_ptr<Tensor>,4>{{Gamma235, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma235, Gamma235, task60);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma240_() {
  auto Gamma240 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma240 = Gamma240->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task61 = make_shared<Task61>(array<shared_ptr<Tensor>,4>{{Gamma240, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma240, Gamma240, task61);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma244_() {
  auto Gamma244 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma244 = Gamma244->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task62 = make_shared<Task62>(array<shared_ptr<Tensor>,4>{{Gamma244, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma244, Gamma244, task62);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma250_() {
  auto Gamma250 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma250 = Gamma250->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task63 = make_shared<Task63>(array<shared_ptr<Tensor>,4>{{Gamma250, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma250, Gamma250, task63);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma251_() {
  auto Gamma251 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma251 = Gamma251->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task64 = make_shared<Task64>(array<shared_ptr<Tensor>,4>{{Gamma251, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma251, Gamma251, task64);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma252_() {
  auto Gamma252 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma252 = Gamma252->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task65 = make_shared<Task65>(array<shared_ptr<Tensor>,3>{{Gamma252, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma252, Gamma252, task65);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma276_() {
  auto Gamma276 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma276 = Gamma276->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task66 = make_shared<Task66>(array<shared_ptr<Tensor>,3>{{Gamma276, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma276, Gamma276, task66);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma565_() {
  auto Gamma565 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma565 = Gamma565->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task67 = make_shared<Task67>(array<shared_ptr<Tensor>,5>{{Gamma565, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma565, Gamma565, task67);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma566_() {
  auto Gamma566 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma566 = Gamma566->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task68 = make_shared<Task68>(array<shared_ptr<Tensor>,5>{{Gamma566, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma566, Gamma566, task68);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma569_() {
  auto Gamma569 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma569 = Gamma569->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task69 = make_shared<Task69>(array<shared_ptr<Tensor>,6>{{Gamma569, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma569, Gamma569, task69);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma570_() {
  auto Gamma570 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma570 = Gamma570->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task70 = make_shared<Task70>(array<shared_ptr<Tensor>,6>{{Gamma570, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma570, Gamma570, task70);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma278_() {
  auto Gamma278 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma278 = Gamma278->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task71 = make_shared<Task71>(array<shared_ptr<Tensor>,5>{{Gamma278, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma278, Gamma278, task71);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma296_() {
  auto Gamma296 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma296 = Gamma296->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task72 = make_shared<Task72>(array<shared_ptr<Tensor>,4>{{Gamma296, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma296, Gamma296, task72);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma312_() {
  auto Gamma312 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma312 = Gamma312->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task73 = make_shared<Task73>(array<shared_ptr<Tensor>,4>{{Gamma312, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma312, Gamma312, task73);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma313_() {
  auto Gamma313 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma313 = Gamma313->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task74 = make_shared<Task74>(array<shared_ptr<Tensor>,4>{{Gamma313, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma313, Gamma313, task74);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma338_() {
  auto Gamma338 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma338 = Gamma338->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task75 = make_shared<Task75>(array<shared_ptr<Tensor>,3>{{Gamma338, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma338, Gamma338, task75);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma48_() {
  auto Gamma48 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma48 = Gamma48->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task76 = make_shared<Task76>(array<shared_ptr<Tensor>,3>{{Gamma48, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma48, Gamma48, task76);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma49_() {
  auto Gamma49 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma49 = Gamma49->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task77 = make_shared<Task77>(array<shared_ptr<Tensor>,3>{{Gamma49, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma49, Gamma49, task77);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma50_() {
  auto Gamma50 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma50 = Gamma50->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task78 = make_shared<Task78>(array<shared_ptr<Tensor>,3>{{Gamma50, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma50, Gamma50, task78);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma51_() {
  auto Gamma51 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma51 = Gamma51->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task79 = make_shared<Task79>(array<shared_ptr<Tensor>,2>{{Gamma51, make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma51, Gamma51, task79);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma339_() {
  auto Gamma339 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma339 = Gamma339->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task80 = make_shared<Task80>(array<shared_ptr<Tensor>,4>{{Gamma339, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma339, Gamma339, task80);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma340_() {
  auto Gamma340 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma340 = Gamma340->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task81 = make_shared<Task81>(array<shared_ptr<Tensor>,3>{{Gamma340, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma340, Gamma340, task81);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma341_() {
  auto Gamma341 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma341 = Gamma341->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task82 = make_shared<Task82>(array<shared_ptr<Tensor>,4>{{Gamma341, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma341, Gamma341, task82);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma342_() {
  auto Gamma342 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma342 = Gamma342->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task83 = make_shared<Task83>(array<shared_ptr<Tensor>,4>{{Gamma342, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma342, Gamma342, task83);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma345_() {
  auto Gamma345 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma345 = Gamma345->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task84 = make_shared<Task84>(array<shared_ptr<Tensor>,4>{{Gamma345, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma345, Gamma345, task84);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma346_() {
  auto Gamma346 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma346 = Gamma346->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task85 = make_shared<Task85>(array<shared_ptr<Tensor>,4>{{Gamma346, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma346, Gamma346, task85);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma349_() {
  auto Gamma349 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma349 = Gamma349->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task86 = make_shared<Task86>(array<shared_ptr<Tensor>,4>{{Gamma349, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma349, Gamma349, task86);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma350_() {
  auto Gamma350 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma350 = Gamma350->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task87 = make_shared<Task87>(array<shared_ptr<Tensor>,4>{{Gamma350, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma350, Gamma350, task87);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma351_() {
  auto Gamma351 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma351 = Gamma351->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task88 = make_shared<Task88>(array<shared_ptr<Tensor>,4>{{Gamma351, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma351, Gamma351, task88);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma359_() {
  auto Gamma359 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma359 = Gamma359->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task89 = make_shared<Task89>(array<shared_ptr<Tensor>,3>{{Gamma359, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma359, Gamma359, task89);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma366_() {
  auto Gamma366 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma366 = Gamma366->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task90 = make_shared<Task90>(array<shared_ptr<Tensor>,3>{{Gamma366, make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma366, Gamma366, task90);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma553_() {
  auto Gamma553 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma553 = Gamma553->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task91 = make_shared<Task91>(array<shared_ptr<Tensor>,5>{{Gamma553, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma553, Gamma553, task91);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma554_() {
  auto Gamma554 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma554 = Gamma554->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task92 = make_shared<Task92>(array<shared_ptr<Tensor>,5>{{Gamma554, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma554, Gamma554, task92);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma471_() {
  auto Gamma471 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma471 = Gamma471->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task93 = make_shared<Task93>(array<shared_ptr<Tensor>,3>{{Gamma471, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma471, Gamma471, task93);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma503_() {
  auto Gamma503 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma503 = Gamma503->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task94 = make_shared<Task94>(array<shared_ptr<Tensor>,2>{{Gamma503, make_shared<Tensor>(*rdm2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma503, Gamma503, task94);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma524_() {
  auto Gamma524 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma524 = Gamma524->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task95 = make_shared<Task95>(array<shared_ptr<Tensor>,2>{{Gamma524, make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma524, Gamma524, task95);
}

shared_ptr<FutureTATensor<double,2>> MRCI::MRCI::Gamma559_() {
  auto Gamma559 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma559 = Gamma559->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task96 = make_shared<Task96>(array<shared_ptr<Tensor>,4>{{Gamma559, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma559, Gamma559, task96);
}

shared_ptr<FutureTATensor<double,2>> MRCI::MRCI::Gamma561_() {
  auto Gamma561 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_});
  auto TAGamma561 = Gamma561->tiledarray<2>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task97 = make_shared<Task97>(array<shared_ptr<Tensor>,5>{{Gamma561, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<double,2>>(*TAGamma561, Gamma561, task97);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma529_() {
  auto Gamma529 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma529 = Gamma529->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task98 = make_shared<Task98>(array<shared_ptr<Tensor>,3>{{Gamma529, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma529, Gamma529, task98);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma530_() {
  auto Gamma530 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma530 = Gamma530->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task99 = make_shared<Task99>(array<shared_ptr<Tensor>,3>{{Gamma530, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma530, Gamma530, task99);
}

shared_ptr<FutureTATensor<double,8>> MRCI::MRCI::Gamma531_() {
  auto Gamma531 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_, active_, active_});
  auto TAGamma531 = Gamma531->tiledarray<8>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task100 = make_shared<Task100>(array<shared_ptr<Tensor>,3>{{Gamma531, make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_)}}, pindex);
  return make_shared<FutureTATensor<double,8>>(*TAGamma531, Gamma531, task100);
}

shared_ptr<FutureTATensor<double,6>> MRCI::MRCI::Gamma543_() {
  auto Gamma543 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_, active_, active_});
  auto TAGamma543 = Gamma543->tiledarray<6>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task101 = make_shared<Task101>(array<shared_ptr<Tensor>,2>{{Gamma543, make_shared<Tensor>(*rdm3_)}}, pindex);
  return make_shared<FutureTATensor<double,6>>(*TAGamma543, Gamma543, task101);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma545_() {
  auto Gamma545 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma545 = Gamma545->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task102 = make_shared<Task102>(array<shared_ptr<Tensor>,6>{{Gamma545, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma545, Gamma545, task102);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma546_() {
  auto Gamma546 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma546 = Gamma546->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task103 = make_shared<Task103>(array<shared_ptr<Tensor>,7>{{Gamma546, make_shared<Tensor>(*rdm0_), make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma546, Gamma546, task103);
}

shared_ptr<FutureTATensor<double,0>> MRCI::MRCI::Gamma555_() {
  auto Gamma555 = make_shared<Tensor>(std::vector<IndexRange>{});
  auto TAGamma555 = Gamma555->tiledarray<0>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task104 = make_shared<Task104>(array<shared_ptr<Tensor>,3>{{Gamma555, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<double,0>>(*TAGamma555, Gamma555, task104);
}

shared_ptr<FutureTATensor<double,0>> MRCI::MRCI::Gamma557_() {
  auto Gamma557 = make_shared<Tensor>(std::vector<IndexRange>{});
  auto TAGamma557 = Gamma557->tiledarray<0>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task105 = make_shared<Task105>(array<shared_ptr<Tensor>,4>{{Gamma557, make_shared<Tensor>(*rdm1_), make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<double,0>>(*TAGamma557, Gamma557, task105);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma563_() {
  auto Gamma563 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma563 = Gamma563->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task106 = make_shared<Task106>(array<shared_ptr<Tensor>,4>{{Gamma563, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*h1_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma563, Gamma563, task106);
}

shared_ptr<FutureTATensor<double,4>> MRCI::MRCI::Gamma564_() {
  auto Gamma564 = make_shared<Tensor>(std::vector<IndexRange>{active_, active_, active_, active_});
  auto TAGamma564 = Gamma564->tiledarray<4>();
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto task107 = make_shared<Task107>(array<shared_ptr<Tensor>,5>{{Gamma564, make_shared<Tensor>(*rdm2_), make_shared<Tensor>(*rdm3_), make_shared<Tensor>(*rdm4_), make_shared<Tensor>(*v2_)}}, pindex);
  return make_shared<FutureTATensor<double,4>>(*TAGamma564, Gamma564, task107);
}

#endif
