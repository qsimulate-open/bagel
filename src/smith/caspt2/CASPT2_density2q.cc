//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_density2qq.cc
// Copyright (C) 2014 Shiozaki group
//
// Author: Shiozaki group <shiozaki@northwestern.edu>
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

shared_ptr<Queue> CASPT2::CASPT2::make_density2q(const bool reset, const bool diagonal) {

  auto density2q = make_shared<Queue>();
  auto task514 = make_shared<Task514>(Den1, reset);
  density2q->add_task(task514);

  auto I672 = make_shared<TATensor<double,4>>({closed_, closed_, active_, active_});
  auto task515 = make_shared<Task515>(Den1, I672);
  task515->add_dep(task514);
  density2q->add_task(task515);

  auto task516 = make_shared<Task516>(I672, Gamma92_(), t2);
  task515->add_dep(task516);
  task516->add_dep(task514);
  density2q->add_task(task516);

  auto I674 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task517 = make_shared<Task517>(Den1, I674);
  task517->add_dep(task514);
  density2q->add_task(task517);

  auto task518 = make_shared<Task518>(I674, Gamma6_(), t2);
  task517->add_dep(task518);
  task518->add_dep(task514);
  density2q->add_task(task518);

  auto I676 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, active_});
  auto task519 = make_shared<Task519>(Den1, I676);
  task519->add_dep(task514);
  density2q->add_task(task519);

  auto I677 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, active_});
  auto task520 = make_shared<Task520>(I676, Gamma16_(), I677);
  task519->add_dep(task520);
  task520->add_dep(task514);
  density2q->add_task(task520);

  auto task521 = make_shared<Task521>(I677, t2);
  task520->add_dep(task521);
  task521->add_dep(task514);
  density2q->add_task(task521);

  auto I680 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task522 = make_shared<Task522>(Den1, I680);
  task522->add_dep(task514);
  density2q->add_task(task522);

  auto task523 = make_shared<Task523>(I680, Gamma32_(), t2);
  task522->add_dep(task523);
  task523->add_dep(task514);
  density2q->add_task(task523);

  auto task524 = make_shared<Task524>(I680, Gamma35_(), t2);
  task522->add_dep(task524);
  task524->add_dep(task514);
  density2q->add_task(task524);

  auto I684 = make_shared<TATensor<double,4>>({virt_, closed_, active_, active_});
  auto task525 = make_shared<Task525>(Den1, I684);
  task525->add_dep(task514);
  density2q->add_task(task525);

  auto I685 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task526 = make_shared<Task526>(I684, Gamma35_(), I685);
  task525->add_dep(task526);
  task526->add_dep(task514);
  density2q->add_task(task526);

  auto task527 = make_shared<Task527>(I685, t2);
  task526->add_dep(task527);
  task527->add_dep(task514);
  density2q->add_task(task527);

  auto I688 = make_shared<TATensor<double,4>>({virt_, active_, active_, active_});
  auto task528 = make_shared<Task528>(Den1, I688);
  task528->add_dep(task514);
  density2q->add_task(task528);

  auto task529 = make_shared<Task529>(I688, Gamma59_(), t2);
  task528->add_dep(task529);
  task529->add_dep(task514);
  density2q->add_task(task529);

  shared_ptr<TATensor<double,4>> I690;
  if (diagonal) {
    I690 = make_shared<TATensor<double,4>>({closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task530> task530;
  if (diagonal) {
    task530 = make_shared<Task530>(Den1, I690);
    task530->add_dep(task514);
    density2q->add_task(task530);
  }

  shared_ptr<Task531> task531;
  if (diagonal) {
    task531 = make_shared<Task531>(I690, t2);
    task530->add_dep(task531);
    task531->add_dep(task514);
    density2q->add_task(task531);
  }

  auto I692 = make_shared<TATensor<double,4>>({virt_, closed_, virt_, active_});
  auto task532 = make_shared<Task532>(Den1, I692);
  task532->add_dep(task514);
  density2q->add_task(task532);

  auto I693 = make_shared<TATensor<double,4>>({active_, virt_, closed_, virt_});
  auto task533 = make_shared<Task533>(I692, Gamma38_(), I693);
  task532->add_dep(task533);
  task533->add_dep(task514);
  density2q->add_task(task533);

  auto task534 = make_shared<Task534>(I693, t2);
  task533->add_dep(task534);
  task534->add_dep(task514);
  density2q->add_task(task534);

  auto I696 = make_shared<TATensor<double,4>>({virt_, virt_, active_, active_});
  auto task535 = make_shared<Task535>(Den1, I696);
  task535->add_dep(task514);
  density2q->add_task(task535);

  auto task536 = make_shared<Task536>(I696, Gamma60_(), t2);
  task535->add_dep(task536);
  task536->add_dep(task514);
  density2q->add_task(task536);

  return density2q;
}


#endif
