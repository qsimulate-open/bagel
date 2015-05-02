//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_gamma.cc
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
#include <src/smith/RelMRCI.h>
#include <src/smith/RelMRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;
using namespace bagel::SMITH::RelMRCI;

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma0_() {
  vector<IndexRange> Gamma0_index = {active_, active_};
  auto Gamma0 = make_shared<Tensor>(Gamma0_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor0 = vector<shared_ptr<Tensor>>{Gamma0, rdm1_};
  auto task0 = make_shared<Task0>(tensor0, pindex);
  return make_shared<FutureTensor>(*Gamma0, task0);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma4_() {
  vector<IndexRange> Gamma4_index = {active_, active_};
  auto Gamma4 = make_shared<Tensor>(Gamma4_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor1 = vector<shared_ptr<Tensor>>{Gamma4, rdm0_, rdm1_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  return make_shared<FutureTensor>(*Gamma4, task1);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma16_() {
  vector<IndexRange> Gamma16_index;
  auto Gamma16 = make_shared<Tensor>(Gamma16_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor2 = vector<shared_ptr<Tensor>>{Gamma16, rdm1_, h1_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  return make_shared<FutureTensor>(*Gamma16, task2);
}

shared_ptr<FutureTensor> RelMRCI::RelMRCI::Gamma18_() {
  vector<IndexRange> Gamma18_index;
  auto Gamma18 = make_shared<Tensor>(Gamma18_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto tensor3 = vector<shared_ptr<Tensor>>{Gamma18, rdm1_, rdm2_, v2_};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  return make_shared<FutureTensor>(*Gamma18, task3);
}

#endif
