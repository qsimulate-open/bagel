//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_deci2q.cc
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


#include <src/smith/caspt2/MSCASPT2.h>
#include <src/smith/caspt2/MSCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_deci2q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deci2q = make_shared<Queue>();
  auto tensor502 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci};
  auto task502 = make_shared<Task502>(tensor502, reset);
  deci2q->add_task(task502);

  vector<IndexRange> I703_index = {active_, active_, active_, active_};
  auto I703 = make_shared<Tensor>(I703_index);
  auto tensor504 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I703};
  auto task504 = make_shared<Task504>(tensor504, cindex);
  task504->add_dep(task502);
  deci2q->add_task(task504);

  auto tensor505 = vector<shared_ptr<Tensor>>{I703, t2, l2};
  auto task505 = make_shared<Task505>(tensor505, cindex, this->e0_);
  task504->add_dep(task505);
  task505->add_dep(task502);
  deci2q->add_task(task505);

  vector<IndexRange> I706_index = {active_, active_, active_, active_, active_, active_};
  auto I706 = make_shared<Tensor>(I706_index);
  auto tensor506 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I706};
  auto task506 = make_shared<Task506>(tensor506, cindex);
  task506->add_dep(task502);
  deci2q->add_task(task506);

  auto tensor507 = vector<shared_ptr<Tensor>>{I706, t2, l2};
  auto task507 = make_shared<Task507>(tensor507, cindex, this->e0_);
  task506->add_dep(task507);
  task507->add_dep(task502);
  deci2q->add_task(task507);

  vector<IndexRange> I709_index = {active_, active_};
  auto I709 = make_shared<Tensor>(I709_index);
  auto tensor508 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I709};
  auto task508 = make_shared<Task508>(tensor508, cindex);
  task508->add_dep(task502);
  deci2q->add_task(task508);

  auto tensor509 = vector<shared_ptr<Tensor>>{I709, t2, l2};
  auto task509 = make_shared<Task509>(tensor509, cindex, this->e0_);
  task508->add_dep(task509);
  task509->add_dep(task502);
  deci2q->add_task(task509);

  auto tensor510 = vector<shared_ptr<Tensor>>{I709, t2, l2};
  auto task510 = make_shared<Task510>(tensor510, cindex, this->e0_);
  task508->add_dep(task510);
  task510->add_dep(task502);
  deci2q->add_task(task510);

  vector<IndexRange> I715_index = {active_, active_, active_, active_};
  auto I715 = make_shared<Tensor>(I715_index);
  auto tensor511 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I715};
  auto task511 = make_shared<Task511>(tensor511, cindex);
  task511->add_dep(task502);
  deci2q->add_task(task511);

  auto tensor512 = vector<shared_ptr<Tensor>>{I715, t2, l2};
  auto task512 = make_shared<Task512>(tensor512, cindex, this->e0_);
  task511->add_dep(task512);
  task512->add_dep(task502);
  deci2q->add_task(task512);

  vector<IndexRange> I718_index = {active_, active_, active_, active_};
  auto I718 = make_shared<Tensor>(I718_index);
  auto tensor513 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I718};
  auto task513 = make_shared<Task513>(tensor513, cindex);
  task513->add_dep(task502);
  deci2q->add_task(task513);

  auto tensor514 = vector<shared_ptr<Tensor>>{I718, t2, l2};
  auto task514 = make_shared<Task514>(tensor514, cindex, this->e0_);
  task513->add_dep(task514);
  task514->add_dep(task502);
  deci2q->add_task(task514);

  auto tensor515 = vector<shared_ptr<Tensor>>{I718, t2, l2};
  auto task515 = make_shared<Task515>(tensor515, cindex, this->e0_);
  task513->add_dep(task515);
  task515->add_dep(task502);
  deci2q->add_task(task515);

  auto tensor516 = vector<shared_ptr<Tensor>>{I718, t2, l2};
  auto task516 = make_shared<Task516>(tensor516, cindex, this->e0_);
  task513->add_dep(task516);
  task516->add_dep(task502);
  deci2q->add_task(task516);

  vector<IndexRange> I727_index = {active_, active_, active_, active_, active_, active_};
  auto I727 = make_shared<Tensor>(I727_index);
  auto tensor517 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I727};
  auto task517 = make_shared<Task517>(tensor517, cindex);
  task517->add_dep(task502);
  deci2q->add_task(task517);

  auto tensor518 = vector<shared_ptr<Tensor>>{I727, t2, l2};
  auto task518 = make_shared<Task518>(tensor518, cindex, this->e0_);
  task517->add_dep(task518);
  task518->add_dep(task502);
  deci2q->add_task(task518);

  shared_ptr<Tensor> I730;
  if (diagonal) {
    vector<IndexRange> I730_index;
    I730 = make_shared<Tensor>(I730_index);
  }
  shared_ptr<Task519> task519;
  if (diagonal) {
    auto tensor519 = vector<shared_ptr<Tensor>>{den0ci, I730};
    task519 = make_shared<Task519>(tensor519, cindex);
    task519->add_dep(task502);
    deci2q->add_task(task519);
  }

  shared_ptr<Task520> task520;
  if (diagonal) {
    auto tensor520 = vector<shared_ptr<Tensor>>{I730, t2, l2};
    task520 = make_shared<Task520>(tensor520, cindex, this->e0_);
    task519->add_dep(task520);
    task520->add_dep(task502);
    deci2q->add_task(task520);
  }

  shared_ptr<Tensor> I733;
  if (diagonal) {
    vector<IndexRange> I733_index;
    I733 = make_shared<Tensor>(I733_index);
  }
  shared_ptr<Task521> task521;
  if (diagonal) {
    auto tensor521 = vector<shared_ptr<Tensor>>{den0ci, I733};
    task521 = make_shared<Task521>(tensor521, cindex);
    task521->add_dep(task502);
    deci2q->add_task(task521);
  }

  shared_ptr<Task522> task522;
  if (diagonal) {
    auto tensor522 = vector<shared_ptr<Tensor>>{I733, t2, l2};
    task522 = make_shared<Task522>(tensor522, cindex, this->e0_);
    task521->add_dep(task522);
    task522->add_dep(task502);
    deci2q->add_task(task522);
  }

  vector<IndexRange> I736_index = {active_, active_};
  auto I736 = make_shared<Tensor>(I736_index);
  auto tensor523 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I736};
  auto task523 = make_shared<Task523>(tensor523, cindex);
  task523->add_dep(task502);
  deci2q->add_task(task523);

  auto tensor524 = vector<shared_ptr<Tensor>>{I736, t2, l2};
  auto task524 = make_shared<Task524>(tensor524, cindex, this->e0_);
  task523->add_dep(task524);
  task524->add_dep(task502);
  deci2q->add_task(task524);

  auto tensor525 = vector<shared_ptr<Tensor>>{I736, t2, l2};
  auto task525 = make_shared<Task525>(tensor525, cindex, this->e0_);
  task523->add_dep(task525);
  task525->add_dep(task502);
  deci2q->add_task(task525);

  vector<IndexRange> I742_index = {active_, active_, active_, active_};
  auto I742 = make_shared<Tensor>(I742_index);
  auto tensor526 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I742};
  auto task526 = make_shared<Task526>(tensor526, cindex);
  task526->add_dep(task502);
  deci2q->add_task(task526);

  auto tensor527 = vector<shared_ptr<Tensor>>{I742, t2, l2};
  auto task527 = make_shared<Task527>(tensor527, cindex, this->e0_);
  task526->add_dep(task527);
  task527->add_dep(task502);
  deci2q->add_task(task527);

  return deci2q;
}


#endif
