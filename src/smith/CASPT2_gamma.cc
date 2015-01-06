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
  vector<IndexRange> Gamma0_index = {active_, active_};
  auto Gamma0 = make_shared<Tensor>(Gamma0_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor0 = {Gamma0, rdm1_};
  auto task0 = make_shared<Task0>(tensor0, pindex);
  return make_shared<FutureTensor>(*Gamma0, task0);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma4_() {
  vector<IndexRange> Gamma4_index = {active_, active_};
  auto Gamma4 = make_shared<Tensor>(Gamma4_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor1 = {Gamma4, rdm2_, f1_};
  auto task1 = make_shared<Task1>(tensor1, pindex);
  return make_shared<FutureTensor>(*Gamma4, task1);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma26_() {
  vector<IndexRange> Gamma26_index = {active_, active_, active_, active_};
  auto Gamma26 = make_shared<Tensor>(Gamma26_index);
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  vector<shared_ptr<Tensor>> tensor2 = {Gamma26, rdm2_};
  auto task2 = make_shared<Task2>(tensor2, pindex);
  return make_shared<FutureTensor>(*Gamma26, task2);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma36_() {
  vector<IndexRange> Gamma36_index = {ci_};
  auto Gamma36 = make_shared<Tensor>(Gamma36_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor3 = {Gamma36, rdm1deriv_, f1_};
  auto task3 = make_shared<Task3>(tensor3, cindex);
  return make_shared<FutureTensor>(*Gamma36, task3);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma38_() {
  vector<IndexRange> Gamma38_index = {ci_, active_, active_};
  auto Gamma38 = make_shared<Tensor>(Gamma38_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor4 = {Gamma38, rdm1deriv_};
  auto task4 = make_shared<Task4>(tensor4, cindex);
  return make_shared<FutureTensor>(*Gamma38, task4);
}

shared_ptr<FutureTensor> CASPT2::CASPT2::Gamma42_() {
  vector<IndexRange> Gamma42_index = {ci_, active_, active_};
  auto Gamma42 = make_shared<Tensor>(Gamma42_index);
  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};
  vector<shared_ptr<Tensor>> tensor5 = {Gamma42, rdm2deriv_, f1_};
  auto task5 = make_shared<Task5>(tensor5, cindex);
  return make_shared<FutureTensor>(*Gamma42, task5);
}

