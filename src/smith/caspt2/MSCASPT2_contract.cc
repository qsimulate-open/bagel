//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_contract.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
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


#include <src/smith/caspt2/MSCASPT2.h>
#include <src/smith/caspt2/MSCASPT2_contract_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MSCASPT2::MSCASPT2::contract_rdm_deriv(shared_ptr<const CIWfn> ciwfn, shared_ptr<VectorB> bdata, int offset, int size, const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto contract = make_shared<Queue>();
  auto tensor900 = vector<shared_ptr<Tensor>>{deci};
  auto task900 = make_shared<Task900>(tensor900, reset);
  contract->add_task(task900);

  auto tensor901 = vector<shared_ptr<Tensor>>{deci, rdm0deriv_, den0cit};
  auto task901 = make_shared<Task901>(tensor901, cindex);
  task901->add_dep(task900);
  contract->add_task(task901);

  auto tensor902 = vector<shared_ptr<Tensor>>{deci, rdm1deriv_, den1cit};
  auto task902 = make_shared<Task902>(tensor902, cindex);
  task902->add_dep(task900);
  contract->add_task(task902);

  auto tensor903 = vector<shared_ptr<Tensor>>{deci, rdm2deriv_, den2cit};
  auto task903 = make_shared<Task903>(tensor903, cindex);
  task903->add_dep(task900);
  contract->add_task(task903);

  vector<IndexRange> I900_index = {ci_, active_, active_};
  auto I900 = make_shared<Tensor>(I900_index);
  auto tensor914 = vector<shared_ptr<Tensor>>{deci, I900};
  auto task914 = make_shared<Task914>(tensor914, cindex, ciwfn, bdata, offset, size);
  task914->add_dep(task900);
  contract->add_task(task914);

  vector<IndexRange> I901_index = {ci_, active_, active_};
  auto I901 = make_shared<Tensor>(I901_index);

  vector<IndexRange> I903_index = {active_, active_, active_, active_, active_, active_};
  auto I903 = make_shared<Tensor>(I903_index);
  auto tensor915 = vector<shared_ptr<Tensor>>{I901, rdm2deriv_, den3cit, I903, deci};
  auto task915 = make_shared<Task915>(tensor915, cindex);
  task914->add_dep(task915);
  task915->add_dep(task900);
  contract->add_task(task915);

  auto tensor921 = vector<shared_ptr<Tensor>>{I903, den4cit, f1_};
  auto task921 = make_shared<Task921>(tensor921, cindex);
  task915->add_dep(task921);
  task921->add_dep(task900);
  contract->add_task(task921);

  auto tensor916 = vector<shared_ptr<Tensor>>{deci, rdm2deriv_, den3cit, I903};
  auto task916 = make_shared<Task916>(tensor916, cindex);
  task916->add_dep(task921);
  task916->add_dep(task900);
  contract->add_task(task916);

  vector<IndexRange> I902_index = {ci_, active_, active_};
  auto I902 = make_shared<Tensor>(I902_index);
  auto tensor917 = vector<shared_ptr<Tensor>>{I900, I901, I902, deci};
  auto task917 = make_shared<Task917>(tensor917, cindex);
  task917->add_dep(task900);
  task917->add_dep(task916);
  task917->add_dep(task915);
  task914->add_dep(task917);
  contract->add_task(task917);

  auto tensor918 = vector<shared_ptr<Tensor>>{I902, rdm3fderiv_, den4cit, deci};
  auto task918 = make_shared<Task918>(tensor918, cindex);
  task917->add_dep(task918);
  task918->add_dep(task900);
  contract->add_task(task918);

  auto tensor923 = vector<shared_ptr<Tensor>>{deci, rdm3fderiv_, den4cit};
  auto task923 = make_shared<Task923>(tensor923, cindex);
  task923->add_dep(task900);
  contract->add_task(task923);

  return contract;
}

#endif
