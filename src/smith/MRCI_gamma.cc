//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_gamma.cc
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

#include <src/smith/MRCI.h>
#include <src/smith/MRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<FutureTensor> MRCI::MRCI::Gamma0_() {
  vector<IndexRange> Gamma0_index = {active_, active_, active_, active_};
  auto Gamma0 = make_shared<Tensor>(Gamma0_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor0 = {Gamma0, rdm3_, h1_};
  auto task0 = make_shared<Task0>(tensor0, pindex);
  return make_shared<FutureTensor>(*Gamma0, task0);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma2_() {
  vector<IndexRange> Gamma2_index = {active_, active_, active_, active_};
  auto Gamma2 = make_shared<Tensor>(Gamma2_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor1 = {Gamma2, rdm3_, rdm4_, v2_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  return make_shared<FutureTensor>(*Gamma2, task1);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma7_() {
  vector<IndexRange> Gamma7_index = {active_, active_, active_, active_};
  auto Gamma7 = make_shared<Tensor>(Gamma7_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor2 = {Gamma7, rdm2_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  return make_shared<FutureTensor>(*Gamma7, task2);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma8_() {
  vector<IndexRange> Gamma8_index = {active_, active_, active_, active_};
  auto Gamma8 = make_shared<Tensor>(Gamma8_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor3 = {Gamma8, rdm2_};
  auto task3 = make_shared<Task3>(tensor3, pindex);
  return make_shared<FutureTensor>(*Gamma8, task3);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma3_() {
  vector<IndexRange> Gamma3_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma3 = make_shared<Tensor>(Gamma3_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor4 = {Gamma3, rdm3_};
  auto task4 = make_shared<Task4>(tensor4, pindex);
  return make_shared<FutureTensor>(*Gamma3, task4);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma4_() {
  vector<IndexRange> Gamma4_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma4 = make_shared<Tensor>(Gamma4_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor5 = {Gamma4, rdm2_, rdm3_};
  auto task5 = make_shared<Task5>(tensor5, pindex);
  return make_shared<FutureTensor>(*Gamma4, task5);
}

shared_ptr<FutureTensor> MRCI::MRCI::Gamma5_() {
  vector<IndexRange> Gamma5_index = {active_, active_, active_, active_, active_, active_};
  auto Gamma5 = make_shared<Tensor>(Gamma5_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor6 = {Gamma5, rdm3_};
  auto task6 = make_shared<Task6>(tensor6, pindex);
  return make_shared<FutureTensor>(*Gamma5, task6);
}

#endif
