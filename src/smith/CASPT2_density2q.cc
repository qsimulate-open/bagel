//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_density2qq.cc
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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_density2q(const bool reset, const bool diagonal) {

  auto density2q = make_shared<Queue>();
  auto task531 = make_shared<Task531>(Den1, reset);
  density2q->add_task(task531);

  auto I742 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, closed_, active_, active_});
  auto task532 = make_shared<Task532>(Den1, I742);
  task532->add_dep(task531);
  density2q->add_task(task532);

  auto task533 = make_shared<Task533>(I742, Gamma92_(), t2);
  task532->add_dep(task533);
  task533->add_dep(task531);
  density2q->add_task(task533);

  auto I744 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, active_, active_, active_});
  auto task534 = make_shared<Task534>(Den1, I744);
  task534->add_dep(task531);
  density2q->add_task(task534);

  auto task535 = make_shared<Task535>(I744, Gamma6_(), t2);
  task534->add_dep(task535);
  task535->add_dep(task531);
  density2q->add_task(task535);

  auto I746 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task536 = make_shared<Task536>(Den1, I746);
  task536->add_dep(task531);
  density2q->add_task(task536);

  auto I747 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, active_});
  auto task537 = make_shared<Task537>(I746, Gamma16_(), I747);
  task536->add_dep(task537);
  task537->add_dep(task531);
  density2q->add_task(task537);

  auto task538 = make_shared<Task538>(I747, t2);
  task537->add_dep(task538);
  task538->add_dep(task531);
  density2q->add_task(task538);

  auto I750 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task539 = make_shared<Task539>(Den1, I750);
  task539->add_dep(task531);
  density2q->add_task(task539);

  auto task540 = make_shared<Task540>(I750, Gamma32_(), t2);
  task539->add_dep(task540);
  task540->add_dep(task531);
  density2q->add_task(task540);

  auto task541 = make_shared<Task541>(I750, Gamma35_(), t2);
  task539->add_dep(task541);
  task541->add_dep(task531);
  density2q->add_task(task541);

  auto I754 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, active_, active_});
  auto task542 = make_shared<Task542>(Den1, I754);
  task542->add_dep(task531);
  density2q->add_task(task542);

  auto I755 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task543 = make_shared<Task543>(I754, Gamma35_(), I755);
  task542->add_dep(task543);
  task543->add_dep(task531);
  density2q->add_task(task543);

  auto task544 = make_shared<Task544>(I755, t2);
  task543->add_dep(task544);
  task544->add_dep(task531);
  density2q->add_task(task544);

  auto I758 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, active_, active_, active_});
  auto task545 = make_shared<Task545>(Den1, I758);
  task545->add_dep(task531);
  density2q->add_task(task545);

  auto task546 = make_shared<Task546>(I758, Gamma59_(), t2);
  task545->add_dep(task546);
  task546->add_dep(task531);
  density2q->add_task(task546);

  shared_ptr<TATensor<double,4>> I760;
  if (diagonal) {
    I760 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{closed_, virt_, closed_, virt_});
  }
  shared_ptr<Task547> task547;
  if (diagonal) {
    task547 = make_shared<Task547>(Den1, I760);
    task547->add_dep(task531);
    density2q->add_task(task547);
  }

  shared_ptr<Task548> task548;
  if (diagonal) {
    task548 = make_shared<Task548>(I760, t2);
    task547->add_dep(task548);
    task548->add_dep(task531);
    density2q->add_task(task548);
  }

  auto I762 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, closed_, virt_, active_});
  auto task549 = make_shared<Task549>(Den1, I762);
  task549->add_dep(task531);
  density2q->add_task(task549);

  auto I763 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, virt_});
  auto task550 = make_shared<Task550>(I762, Gamma38_(), I763);
  task549->add_dep(task550);
  task550->add_dep(task531);
  density2q->add_task(task550);

  auto task551 = make_shared<Task551>(I763, t2);
  task550->add_dep(task551);
  task551->add_dep(task531);
  density2q->add_task(task551);

  auto I766 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{virt_, virt_, active_, active_});
  auto task552 = make_shared<Task552>(Den1, I766);
  task552->add_dep(task531);
  density2q->add_task(task552);

  auto task553 = make_shared<Task553>(I766, Gamma60_(), t2);
  task552->add_dep(task553);
  task553->add_dep(task531);
  density2q->add_task(task553);

  return density2q;
}


#endif
