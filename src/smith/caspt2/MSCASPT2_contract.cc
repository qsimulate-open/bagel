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

shared_ptr<Queue> MSCASPT2::MSCASPT2::contract_rdm_deriv(int number, const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto contract = make_shared<Queue>();
  auto tensor900 = vector<shared_ptr<Tensor>>{deci};
  auto task900 = make_shared<Task900>(tensor900, reset);
  contract->add_task(task900);

  auto tensor901 = vector<shared_ptr<Tensor>>{deci, rdm0deriv_, den0ci};
  auto task901 = make_shared<Task901>(tensor901, cindex);
  task901->add_dep(task900);
  contract->add_task(task901);

  if (number >= 1) {
    auto tensor902 = vector<shared_ptr<Tensor>>{deci, rdm1deriv_, den1ci};
    auto task902 = make_shared<Task902>(tensor902, cindex);
    task902->add_dep(task900);
    contract->add_task(task902);
  }

  if (number >= 2) {
    auto tensor903 = vector<shared_ptr<Tensor>>{deci, rdm2deriv_, den2ci};
    auto task903 = make_shared<Task903>(tensor903, cindex);
    task903->add_dep(task900);
    contract->add_task(task903);
  }

  if (number >= 3) {
    auto tensor904 = vector<shared_ptr<Tensor>>{deci, rdm3deriv_, den3ci};
    auto task904 = make_shared<Task904>(tensor904, cindex);
    task904->add_dep(task900);
    contract->add_task(task904);
  }

  if (number >= 4) {
    auto tensor905 = vector<shared_ptr<Tensor>>{deci, rdm4deriv_, den4ci};
    auto task905 = make_shared<Task905>(tensor905, cindex);
    task905->add_dep(task900);
    contract->add_task(task905);
  }

  return contract;
}

#endif
