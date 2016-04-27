//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: MSCASPT2_deci2qq.cc
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
#include <src/smith/caspt2/MSCASPT2_tasks11.h>
#include <src/smith/caspt2/MSCASPT2_tasks12.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MSCASPT2::MSCASPT2::make_deci2q(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deci2q = make_shared<Queue>();
  auto tensor517 = vector<shared_ptr<Tensor>>{deci};
  auto task517 = make_shared<Task517>(tensor517, reset);
  deci2q->add_task(task517);

  vector<IndexRange> I744_index = {ci_};
  auto I744 = make_shared<Tensor>(I744_index);
  auto tensor518 = vector<shared_ptr<Tensor>>{deci, I744};
  auto task518 = make_shared<Task518>(tensor518, cindex);
  task518->add_dep(task517);
  deci2q->add_task(task518);

  vector<IndexRange> I745_index = {active_, active_, active_, active_};
  auto I745 = make_shared<Tensor>(I745_index);
  auto tensor519 = vector<shared_ptr<Tensor>>{I744, Gamma111_(), I745};
  auto task519 = make_shared<Task519>(tensor519, cindex);
  task518->add_dep(task519);
  task519->add_dep(task517);
  deci2q->add_task(task519);

  auto tensor520 = vector<shared_ptr<Tensor>>{I745, v2_, l2};
  auto task520 = make_shared<Task520>(tensor520, cindex);
  task519->add_dep(task520);
  task520->add_dep(task517);
  deci2q->add_task(task520);

  vector<IndexRange> I748_index = {active_, active_, active_, active_, active_, active_};
  auto I748 = make_shared<Tensor>(I748_index);
  auto tensor521 = vector<shared_ptr<Tensor>>{I744, Gamma217_(), I748};
  auto task521 = make_shared<Task521>(tensor521, cindex);
  task518->add_dep(task521);
  task521->add_dep(task517);
  deci2q->add_task(task521);

  auto tensor522 = vector<shared_ptr<Tensor>>{I748, v2_, l2};
  auto task522 = make_shared<Task522>(tensor522, cindex);
  task521->add_dep(task522);
  task522->add_dep(task517);
  deci2q->add_task(task522);

  vector<IndexRange> I751_index = {active_, active_, active_, active_, active_, active_};
  auto I751 = make_shared<Tensor>(I751_index);
  auto tensor523 = vector<shared_ptr<Tensor>>{I744, Gamma116_(), I751};
  auto task523 = make_shared<Task523>(tensor523, cindex);
  task518->add_dep(task523);
  task523->add_dep(task517);
  deci2q->add_task(task523);

  auto tensor524 = vector<shared_ptr<Tensor>>{I751, v2_, l2};
  auto task524 = make_shared<Task524>(tensor524, cindex);
  task523->add_dep(task524);
  task524->add_dep(task517);
  deci2q->add_task(task524);

  vector<IndexRange> I754_index = {active_, active_};
  auto I754 = make_shared<Tensor>(I754_index);
  auto tensor525 = vector<shared_ptr<Tensor>>{I744, Gamma126_(), I754};
  auto task525 = make_shared<Task525>(tensor525, cindex);
  task518->add_dep(task525);
  task525->add_dep(task517);
  deci2q->add_task(task525);

  auto tensor526 = vector<shared_ptr<Tensor>>{I754, v2_, l2};
  auto task526 = make_shared<Task526>(tensor526, cindex);
  task525->add_dep(task526);
  task526->add_dep(task517);
  deci2q->add_task(task526);

  auto tensor527 = vector<shared_ptr<Tensor>>{I754, v2_, l2};
  auto task527 = make_shared<Task527>(tensor527, cindex);
  task525->add_dep(task527);
  task527->add_dep(task517);
  deci2q->add_task(task527);

  vector<IndexRange> I760_index = {active_, active_, active_, active_};
  auto I760 = make_shared<Tensor>(I760_index);
  auto tensor528 = vector<shared_ptr<Tensor>>{I744, Gamma145_(), I760};
  auto task528 = make_shared<Task528>(tensor528, cindex);
  task518->add_dep(task528);
  task528->add_dep(task517);
  deci2q->add_task(task528);

  auto tensor529 = vector<shared_ptr<Tensor>>{I760, v2_, l2};
  auto task529 = make_shared<Task529>(tensor529, cindex);
  task528->add_dep(task529);
  task529->add_dep(task517);
  deci2q->add_task(task529);

  auto tensor530 = vector<shared_ptr<Tensor>>{I760, v2_, l2};
  auto task530 = make_shared<Task530>(tensor530, cindex);
  task528->add_dep(task530);
  task530->add_dep(task517);
  deci2q->add_task(task530);

  auto tensor531 = vector<shared_ptr<Tensor>>{I760, v2_, l2};
  auto task531 = make_shared<Task531>(tensor531, cindex);
  task528->add_dep(task531);
  task531->add_dep(task517);
  deci2q->add_task(task531);

  auto tensor532 = vector<shared_ptr<Tensor>>{I760, v2_, l2};
  auto task532 = make_shared<Task532>(tensor532, cindex);
  task528->add_dep(task532);
  task532->add_dep(task517);
  deci2q->add_task(task532);

  auto tensor533 = vector<shared_ptr<Tensor>>{I760, v2_, l2};
  auto task533 = make_shared<Task533>(tensor533, cindex);
  task528->add_dep(task533);
  task533->add_dep(task517);
  deci2q->add_task(task533);

  vector<IndexRange> I763_index = {active_, active_, active_, active_};
  auto I763 = make_shared<Tensor>(I763_index);
  auto tensor534 = vector<shared_ptr<Tensor>>{I744, Gamma139_(), I763};
  auto task534 = make_shared<Task534>(tensor534, cindex);
  task518->add_dep(task534);
  task534->add_dep(task517);
  deci2q->add_task(task534);

  auto tensor535 = vector<shared_ptr<Tensor>>{I763, v2_, l2};
  auto task535 = make_shared<Task535>(tensor535, cindex);
  task534->add_dep(task535);
  task535->add_dep(task517);
  deci2q->add_task(task535);

  vector<IndexRange> I766_index = {active_, active_, active_, active_};
  auto I766 = make_shared<Tensor>(I766_index);
  auto tensor536 = vector<shared_ptr<Tensor>>{I744, Gamma142_(), I766};
  auto task536 = make_shared<Task536>(tensor536, cindex);
  task518->add_dep(task536);
  task536->add_dep(task517);
  deci2q->add_task(task536);

  auto tensor537 = vector<shared_ptr<Tensor>>{I766, v2_, l2};
  auto task537 = make_shared<Task537>(tensor537, cindex);
  task536->add_dep(task537);
  task537->add_dep(task517);
  deci2q->add_task(task537);

  vector<IndexRange> I775_index = {active_, active_, active_, active_};
  auto I775 = make_shared<Tensor>(I775_index);
  auto tensor538 = vector<shared_ptr<Tensor>>{I744, Gamma117_(), I775};
  auto task538 = make_shared<Task538>(tensor538, cindex);
  task518->add_dep(task538);
  task538->add_dep(task517);
  deci2q->add_task(task538);

  auto tensor539 = vector<shared_ptr<Tensor>>{I775, v2_, l2};
  auto task539 = make_shared<Task539>(tensor539, cindex);
  task538->add_dep(task539);
  task539->add_dep(task517);
  deci2q->add_task(task539);

  auto tensor540 = vector<shared_ptr<Tensor>>{I775, h1_, l2};
  auto task540 = make_shared<Task540>(tensor540, cindex);
  task538->add_dep(task540);
  task540->add_dep(task517);
  deci2q->add_task(task540);

  vector<IndexRange> I784_index = {active_, active_, active_, active_, active_, active_};
  auto I784 = make_shared<Tensor>(I784_index);
  auto tensor541 = vector<shared_ptr<Tensor>>{I744, Gamma169_(), I784};
  auto task541 = make_shared<Task541>(tensor541, cindex);
  task518->add_dep(task541);
  task541->add_dep(task517);
  deci2q->add_task(task541);

  auto tensor542 = vector<shared_ptr<Tensor>>{I784, v2_, l2};
  auto task542 = make_shared<Task542>(tensor542, cindex);
  task541->add_dep(task542);
  task542->add_dep(task517);
  deci2q->add_task(task542);

  vector<IndexRange> I787_index = {active_, active_, active_, active_, active_, active_};
  auto I787 = make_shared<Tensor>(I787_index);
  auto tensor543 = vector<shared_ptr<Tensor>>{I744, Gamma167_(), I787};
  auto task543 = make_shared<Task543>(tensor543, cindex);
  task518->add_dep(task543);
  task543->add_dep(task517);
  deci2q->add_task(task543);

  auto tensor544 = vector<shared_ptr<Tensor>>{I787, v2_, l2};
  auto task544 = make_shared<Task544>(tensor544, cindex);
  task543->add_dep(task544);
  task544->add_dep(task517);
  deci2q->add_task(task544);

  shared_ptr<Tensor> I790;
  if (diagonal) {
    vector<IndexRange> I790_index;
    I790 = make_shared<Tensor>(I790_index);
  }
  shared_ptr<Task545> task545;
  if (diagonal) {
    auto tensor545 = vector<shared_ptr<Tensor>>{I744, rdm0deriv_, I790};
    task545 = make_shared<Task545>(tensor545, cindex);
    task518->add_dep(task545);
    task545->add_dep(task517);
    deci2q->add_task(task545);
  }

  shared_ptr<Task546> task546;
  if (diagonal) {
    auto tensor546 = vector<shared_ptr<Tensor>>{I790, v2_, l2};
    task546 = make_shared<Task546>(tensor546, cindex);
    task545->add_dep(task546);
    task546->add_dep(task517);
    deci2q->add_task(task546);
  }

  shared_ptr<Tensor> I793;
  if (diagonal) {
    vector<IndexRange> I793_index;
    I793 = make_shared<Tensor>(I793_index);
  }
  shared_ptr<Task547> task547;
  if (diagonal) {
    auto tensor547 = vector<shared_ptr<Tensor>>{I744, rdm0deriv_, I793};
    task547 = make_shared<Task547>(tensor547, cindex);
    task518->add_dep(task547);
    task547->add_dep(task517);
    deci2q->add_task(task547);
  }

  shared_ptr<Task548> task548;
  if (diagonal) {
    auto tensor548 = vector<shared_ptr<Tensor>>{I793, v2_, l2};
    task548 = make_shared<Task548>(tensor548, cindex);
    task547->add_dep(task548);
    task548->add_dep(task517);
    deci2q->add_task(task548);
  }

  vector<IndexRange> I796_index = {active_, active_};
  auto I796 = make_shared<Tensor>(I796_index);
  auto tensor549 = vector<shared_ptr<Tensor>>{I744, Gamma148_(), I796};
  auto task549 = make_shared<Task549>(tensor549, cindex);
  task518->add_dep(task549);
  task549->add_dep(task517);
  deci2q->add_task(task549);

  auto tensor550 = vector<shared_ptr<Tensor>>{I796, v2_, l2};
  auto task550 = make_shared<Task550>(tensor550, cindex);
  task549->add_dep(task550);
  task550->add_dep(task517);
  deci2q->add_task(task550);

  auto tensor551 = vector<shared_ptr<Tensor>>{I796, v2_, l2};
  auto task551 = make_shared<Task551>(tensor551, cindex);
  task549->add_dep(task551);
  task551->add_dep(task517);
  deci2q->add_task(task551);

  auto tensor552 = vector<shared_ptr<Tensor>>{I796, h1_, l2};
  auto task552 = make_shared<Task552>(tensor552, cindex);
  task549->add_dep(task552);
  task552->add_dep(task517);
  deci2q->add_task(task552);

  auto tensor553 = vector<shared_ptr<Tensor>>{I796, h1_, l2};
  auto task553 = make_shared<Task553>(tensor553, cindex);
  task549->add_dep(task553);
  task553->add_dep(task517);
  deci2q->add_task(task553);

  vector<IndexRange> I802_index = {active_, active_, active_, active_};
  auto I802 = make_shared<Tensor>(I802_index);
  auto tensor554 = vector<shared_ptr<Tensor>>{I744, Gamma170_(), I802};
  auto task554 = make_shared<Task554>(tensor554, cindex);
  task518->add_dep(task554);
  task554->add_dep(task517);
  deci2q->add_task(task554);

  auto tensor555 = vector<shared_ptr<Tensor>>{I802, v2_, l2};
  auto task555 = make_shared<Task555>(tensor555, cindex);
  task554->add_dep(task555);
  task555->add_dep(task517);
  deci2q->add_task(task555);

  auto tensor556 = vector<shared_ptr<Tensor>>{I802, h1_, l2};
  auto task556 = make_shared<Task556>(tensor556, cindex);
  task554->add_dep(task556);
  task556->add_dep(task517);
  deci2q->add_task(task556);

  return deci2q;
}


#endif
