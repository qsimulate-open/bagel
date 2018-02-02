//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_deci3q.cc
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

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_deci3q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto deci3q = make_shared<Queue>();
  auto tensor496 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci};
  auto task496 = make_shared<Task496>(tensor496, reset);
  deci3q->add_task(task496);

  vector<IndexRange> I745_index = {active_, active_, active_, active_};
  auto I745 = make_shared<Tensor>(I745_index);
  auto tensor498 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I745};
  auto task498 = make_shared<Task498>(tensor498, pindex);
  task498->add_dep(task496);
  deci3q->add_task(task498);

  auto tensor499 = vector<shared_ptr<Tensor>>{I745, v2_, l2};
  auto task499 = make_shared<Task499>(tensor499, pindex);
  task498->add_dep(task499);
  task499->add_dep(task496);
  deci3q->add_task(task499);

  vector<IndexRange> I748_index = {active_, active_, active_, active_, active_, active_};
  auto I748 = make_shared<Tensor>(I748_index);
  auto tensor500 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I748};
  auto task500 = make_shared<Task500>(tensor500, pindex);
  task500->add_dep(task496);
  deci3q->add_task(task500);

  auto tensor501 = vector<shared_ptr<Tensor>>{I748, v2_, l2};
  auto task501 = make_shared<Task501>(tensor501, pindex);
  task500->add_dep(task501);
  task501->add_dep(task496);
  deci3q->add_task(task501);

  vector<IndexRange> I751_index = {active_, active_, active_, active_, active_, active_};
  auto I751 = make_shared<Tensor>(I751_index);
  auto tensor502 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I751};
  auto task502 = make_shared<Task502>(tensor502, pindex);
  task502->add_dep(task496);
  deci3q->add_task(task502);

  auto tensor503 = vector<shared_ptr<Tensor>>{I751, v2_, l2};
  auto task503 = make_shared<Task503>(tensor503, pindex);
  task502->add_dep(task503);
  task503->add_dep(task496);
  deci3q->add_task(task503);

  vector<IndexRange> I754_index = {active_, active_};
  auto I754 = make_shared<Tensor>(I754_index);
  auto tensor504 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I754};
  auto task504 = make_shared<Task504>(tensor504, pindex);
  task504->add_dep(task496);
  deci3q->add_task(task504);

  auto tensor505 = vector<shared_ptr<Tensor>>{I754, v2_, l2};
  auto task505 = make_shared<Task505>(tensor505, pindex);
  task504->add_dep(task505);
  task505->add_dep(task496);
  deci3q->add_task(task505);

  auto tensor506 = vector<shared_ptr<Tensor>>{I754, v2_, l2};
  auto task506 = make_shared<Task506>(tensor506, pindex);
  task504->add_dep(task506);
  task506->add_dep(task496);
  deci3q->add_task(task506);

  vector<IndexRange> I760_index = {active_, active_, active_, active_};
  auto I760 = make_shared<Tensor>(I760_index);
  auto tensor507 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I760};
  auto task507 = make_shared<Task507>(tensor507, pindex);
  task507->add_dep(task496);
  deci3q->add_task(task507);

  auto tensor508 = vector<shared_ptr<Tensor>>{I760, v2_, l2};
  auto task508 = make_shared<Task508>(tensor508, pindex);
  task507->add_dep(task508);
  task508->add_dep(task496);
  deci3q->add_task(task508);

  auto tensor509 = vector<shared_ptr<Tensor>>{I760, v2_, l2};
  auto task509 = make_shared<Task509>(tensor509, pindex);
  task507->add_dep(task509);
  task509->add_dep(task496);
  deci3q->add_task(task509);

  auto tensor510 = vector<shared_ptr<Tensor>>{I760, v2_, l2};
  auto task510 = make_shared<Task510>(tensor510, pindex);
  task507->add_dep(task510);
  task510->add_dep(task496);
  deci3q->add_task(task510);

  auto tensor511 = vector<shared_ptr<Tensor>>{I760, v2_, l2};
  auto task511 = make_shared<Task511>(tensor511, pindex);
  task507->add_dep(task511);
  task511->add_dep(task496);
  deci3q->add_task(task511);

  auto tensor512 = vector<shared_ptr<Tensor>>{I760, v2_, l2};
  auto task512 = make_shared<Task512>(tensor512, pindex);
  task507->add_dep(task512);
  task512->add_dep(task496);
  deci3q->add_task(task512);

  vector<IndexRange> I763_index = {active_, active_, active_, active_};
  auto I763 = make_shared<Tensor>(I763_index);
  auto tensor513 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I763};
  auto task513 = make_shared<Task513>(tensor513, pindex);
  task513->add_dep(task496);
  deci3q->add_task(task513);

  auto tensor514 = vector<shared_ptr<Tensor>>{I763, v2_, l2};
  auto task514 = make_shared<Task514>(tensor514, pindex);
  task513->add_dep(task514);
  task514->add_dep(task496);
  deci3q->add_task(task514);

  vector<IndexRange> I766_index = {active_, active_, active_, active_};
  auto I766 = make_shared<Tensor>(I766_index);
  auto tensor515 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I766};
  auto task515 = make_shared<Task515>(tensor515, pindex);
  task515->add_dep(task496);
  deci3q->add_task(task515);

  auto tensor516 = vector<shared_ptr<Tensor>>{I766, v2_, l2};
  auto task516 = make_shared<Task516>(tensor516, pindex);
  task515->add_dep(task516);
  task516->add_dep(task496);
  deci3q->add_task(task516);

  vector<IndexRange> I775_index = {active_, active_, active_, active_};
  auto I775 = make_shared<Tensor>(I775_index);
  auto tensor517 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I775};
  auto task517 = make_shared<Task517>(tensor517, pindex);
  task517->add_dep(task496);
  deci3q->add_task(task517);

  auto tensor518 = vector<shared_ptr<Tensor>>{I775, v2_, l2};
  auto task518 = make_shared<Task518>(tensor518, pindex);
  task517->add_dep(task518);
  task518->add_dep(task496);
  deci3q->add_task(task518);

  auto tensor519 = vector<shared_ptr<Tensor>>{I775, h1_, l2};
  auto task519 = make_shared<Task519>(tensor519, pindex);
  task517->add_dep(task519);
  task519->add_dep(task496);
  deci3q->add_task(task519);

  vector<IndexRange> I784_index = {active_, active_, active_, active_, active_, active_};
  auto I784 = make_shared<Tensor>(I784_index);
  auto tensor520 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I784};
  auto task520 = make_shared<Task520>(tensor520, pindex);
  task520->add_dep(task496);
  deci3q->add_task(task520);

  auto tensor521 = vector<shared_ptr<Tensor>>{I784, v2_, l2};
  auto task521 = make_shared<Task521>(tensor521, pindex);
  task520->add_dep(task521);
  task521->add_dep(task496);
  deci3q->add_task(task521);

  vector<IndexRange> I787_index = {active_, active_, active_, active_, active_, active_};
  auto I787 = make_shared<Tensor>(I787_index);
  auto tensor522 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I787};
  auto task522 = make_shared<Task522>(tensor522, pindex);
  task522->add_dep(task496);
  deci3q->add_task(task522);

  auto tensor523 = vector<shared_ptr<Tensor>>{I787, v2_, l2};
  auto task523 = make_shared<Task523>(tensor523, pindex);
  task522->add_dep(task523);
  task523->add_dep(task496);
  deci3q->add_task(task523);

  shared_ptr<Tensor> I790;
  if (diagonal) {
    vector<IndexRange> I790_index;
    I790 = make_shared<Tensor>(I790_index);
  }
  shared_ptr<Task524> task524;
  if (diagonal) {
    auto tensor524 = vector<shared_ptr<Tensor>>{den0ci, I790};
    task524 = make_shared<Task524>(tensor524, pindex);
    task524->add_dep(task496);
    deci3q->add_task(task524);
  }

  shared_ptr<Task525> task525;
  if (diagonal) {
    auto tensor525 = vector<shared_ptr<Tensor>>{I790, v2_, l2};
    task525 = make_shared<Task525>(tensor525, pindex);
    task524->add_dep(task525);
    task525->add_dep(task496);
    deci3q->add_task(task525);
  }

  shared_ptr<Tensor> I793;
  if (diagonal) {
    vector<IndexRange> I793_index;
    I793 = make_shared<Tensor>(I793_index);
  }
  shared_ptr<Task526> task526;
  if (diagonal) {
    auto tensor526 = vector<shared_ptr<Tensor>>{den0ci, I793};
    task526 = make_shared<Task526>(tensor526, pindex);
    task526->add_dep(task496);
    deci3q->add_task(task526);
  }

  shared_ptr<Task527> task527;
  if (diagonal) {
    auto tensor527 = vector<shared_ptr<Tensor>>{I793, v2_, l2};
    task527 = make_shared<Task527>(tensor527, pindex);
    task526->add_dep(task527);
    task527->add_dep(task496);
    deci3q->add_task(task527);
  }

  vector<IndexRange> I796_index = {active_, active_};
  auto I796 = make_shared<Tensor>(I796_index);
  auto tensor528 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I796};
  auto task528 = make_shared<Task528>(tensor528, pindex);
  task528->add_dep(task496);
  deci3q->add_task(task528);

  auto tensor529 = vector<shared_ptr<Tensor>>{I796, v2_, l2};
  auto task529 = make_shared<Task529>(tensor529, pindex);
  task528->add_dep(task529);
  task529->add_dep(task496);
  deci3q->add_task(task529);

  auto tensor530 = vector<shared_ptr<Tensor>>{I796, v2_, l2};
  auto task530 = make_shared<Task530>(tensor530, pindex);
  task528->add_dep(task530);
  task530->add_dep(task496);
  deci3q->add_task(task530);

  auto tensor531 = vector<shared_ptr<Tensor>>{I796, h1_, l2};
  auto task531 = make_shared<Task531>(tensor531, pindex);
  task528->add_dep(task531);
  task531->add_dep(task496);
  deci3q->add_task(task531);

  auto tensor532 = vector<shared_ptr<Tensor>>{I796, h1_, l2};
  auto task532 = make_shared<Task532>(tensor532, pindex);
  task528->add_dep(task532);
  task532->add_dep(task496);
  deci3q->add_task(task532);

  vector<IndexRange> I802_index = {active_, active_, active_, active_};
  auto I802 = make_shared<Tensor>(I802_index);
  auto tensor533 = vector<shared_ptr<Tensor>>{den0ci, den1ci, den2ci, den3ci, den4ci, I802};
  auto task533 = make_shared<Task533>(tensor533, pindex);
  task533->add_dep(task496);
  deci3q->add_task(task533);

  auto tensor534 = vector<shared_ptr<Tensor>>{I802, v2_, l2};
  auto task534 = make_shared<Task534>(tensor534, pindex);
  task533->add_dep(task534);
  task534->add_dep(task496);
  deci3q->add_task(task534);

  auto tensor535 = vector<shared_ptr<Tensor>>{I802, h1_, l2};
  auto task535 = make_shared<Task535>(tensor535, pindex);
  task533->add_dep(task535);
  task535->add_dep(task496);
  deci3q->add_task(task535);

  return deci3q;
}


#endif
