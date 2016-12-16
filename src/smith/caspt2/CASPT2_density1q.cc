//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_density1qq.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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
#include <src/smith/caspt2/CASPT2_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_density1q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density1q = make_shared<Queue>();
  auto tensor508 = vector<shared_ptr<Tensor>>{den1};
  auto task508 = make_shared<Task508>(tensor508, reset);
  density1q->add_task(task508);

  vector<IndexRange> I650_index = {closed_, active_};
  auto I650 = make_shared<Tensor>(I650_index);
  auto tensor509 = vector<shared_ptr<Tensor>>{den1, I650};
  auto task509 = make_shared<Task509>(tensor509, pindex);
  task509->add_dep(task508);
  density1q->add_task(task509);

  auto tensor510 = vector<shared_ptr<Tensor>>{I650, Gamma12_(), t2};
  auto task510 = make_shared<Task510>(tensor510, pindex);
  task509->add_dep(task510);
  task510->add_dep(task508);
  density1q->add_task(task510);

  vector<IndexRange> I652_index = {virt_, closed_};
  auto I652 = make_shared<Tensor>(I652_index);
  auto tensor511 = vector<shared_ptr<Tensor>>{den1, I652};
  auto task511 = make_shared<Task511>(tensor511, pindex);
  task511->add_dep(task508);
  density1q->add_task(task511);

  auto tensor512 = vector<shared_ptr<Tensor>>{I652, t2, Gamma38_()};
  auto task512 = make_shared<Task512>(tensor512, pindex);
  task511->add_dep(task512);
  task512->add_dep(task508);
  density1q->add_task(task512);

  auto tensor513 = vector<shared_ptr<Tensor>>{I652, t2, Gamma38_()};
  auto task513 = make_shared<Task513>(tensor513, pindex);
  task511->add_dep(task513);
  task513->add_dep(task508);
  density1q->add_task(task513);

  vector<IndexRange> I656_index = {active_, virt_};
  auto I656 = make_shared<Tensor>(I656_index);
  auto tensor514 = vector<shared_ptr<Tensor>>{den1, I656};
  auto task514 = make_shared<Task514>(tensor514, pindex);
  task514->add_dep(task508);
  density1q->add_task(task514);

  auto tensor515 = vector<shared_ptr<Tensor>>{I656, t2, Gamma60_()};
  auto task515 = make_shared<Task515>(tensor515, pindex);
  task514->add_dep(task515);
  task515->add_dep(task508);
  density1q->add_task(task515);

  return density1q;
}


#endif
