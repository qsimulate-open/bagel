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
#include <src/smith/caspt2/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_density1q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto density1q = make_shared<Queue>();
  auto tensor507 = vector<shared_ptr<Tensor>>{den1};
  auto task507 = make_shared<Task507>(tensor507, reset);
  density1q->add_task(task507);

  vector<IndexRange> I664_index = {closed_, active_};
  auto I664 = make_shared<Tensor>(I664_index);
  auto tensor508 = vector<shared_ptr<Tensor>>{den1, I664};
  auto task508 = make_shared<Task508>(tensor508, pindex);
  task508->add_dep(task507);
  density1q->add_task(task508);

  auto tensor509 = vector<shared_ptr<Tensor>>{I664, Gamma12_(), t2};
  auto task509 = make_shared<Task509>(tensor509, pindex);
  task508->add_dep(task509);
  task509->add_dep(task507);
  density1q->add_task(task509);

  vector<IndexRange> I666_index = {virt_, closed_};
  auto I666 = make_shared<Tensor>(I666_index);
  auto tensor510 = vector<shared_ptr<Tensor>>{den1, I666};
  auto task510 = make_shared<Task510>(tensor510, pindex);
  task510->add_dep(task507);
  density1q->add_task(task510);

  auto tensor511 = vector<shared_ptr<Tensor>>{I666, t2, Gamma38_()};
  auto task511 = make_shared<Task511>(tensor511, pindex);
  task510->add_dep(task511);
  task511->add_dep(task507);
  density1q->add_task(task511);

  auto tensor512 = vector<shared_ptr<Tensor>>{I666, t2, Gamma38_()};
  auto task512 = make_shared<Task512>(tensor512, pindex);
  task510->add_dep(task512);
  task512->add_dep(task507);
  density1q->add_task(task512);

  vector<IndexRange> I670_index = {active_, virt_};
  auto I670 = make_shared<Tensor>(I670_index);
  auto tensor513 = vector<shared_ptr<Tensor>>{den1, I670};
  auto task513 = make_shared<Task513>(tensor513, pindex);
  task513->add_dep(task507);
  density1q->add_task(task513);

  auto tensor514 = vector<shared_ptr<Tensor>>{I670, t2, Gamma60_()};
  auto task514 = make_shared<Task514>(tensor514, pindex);
  task513->add_dep(task514);
  task514->add_dep(task507);
  density1q->add_task(task514);

  return density1q;
}


#endif
