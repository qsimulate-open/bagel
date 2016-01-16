//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_corrqq.cc
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


#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_corrq(const bool reset, const bool diagonal) {

  auto corrq = make_shared<Queue>();
  auto I335 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task261 = make_shared<Task261>(Gamma92_(), I335);
  corrq->add_task(task261);

  auto task262 = make_shared<Task262>(I335, t2);
  task261->add_dep(task262);
  corrq->add_task(task262);

  auto I338 = make_shared<TATensor<double,4>>({closed_, active_, active_, active_});
  auto task263 = make_shared<Task263>(t2, I338);
  task261->add_dep(task263);
  corrq->add_task(task263);

  auto task264 = make_shared<Task264>(I338, Gamma6_(), t2);
  task263->add_dep(task264);
  corrq->add_task(task264);

  auto I341 = make_shared<TATensor<double,2>>({active_, active_});
  auto task265 = make_shared<Task265>(Gamma16_(), I341);
  task261->add_dep(task265);
  corrq->add_task(task265);

  auto task266 = make_shared<Task266>(I341, t2);
  task265->add_dep(task266);
  corrq->add_task(task266);

  auto task267 = make_shared<Task267>(I341, t2);
  task265->add_dep(task267);
  corrq->add_task(task267);

  auto I347 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task268 = make_shared<Task268>(Gamma32_(), I347);
  task261->add_dep(task268);
  corrq->add_task(task268);

  auto task269 = make_shared<Task269>(I347, t2);
  task268->add_dep(task269);
  corrq->add_task(task269);

  auto I350 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task270 = make_shared<Task270>(Gamma35_(), I350);
  task261->add_dep(task270);
  corrq->add_task(task270);

  auto task271 = make_shared<Task271>(I350, t2);
  task270->add_dep(task271);
  corrq->add_task(task271);

  auto task272 = make_shared<Task272>(I350, t2);
  task270->add_dep(task272);
  corrq->add_task(task272);

  auto task273 = make_shared<Task273>(I350, t2);
  task270->add_dep(task273);
  corrq->add_task(task273);

  auto I359 = make_shared<TATensor<double,6>>({active_, active_, active_, active_, active_, active_});
  auto task274 = make_shared<Task274>(Gamma59_(), I359);
  task261->add_dep(task274);
  corrq->add_task(task274);

  auto task275 = make_shared<Task275>(I359, t2);
  task274->add_dep(task275);
  corrq->add_task(task275);

  shared_ptr<Task276> task276;
  if (diagonal) {
    task276 = make_shared<Task276>(t2);
    task261->add_dep(task276);
    corrq->add_task(task276);
  }

  shared_ptr<Task277> task277;
  if (diagonal) {
    task277 = make_shared<Task277>(t2);
    task261->add_dep(task277);
    corrq->add_task(task277);
  }

  auto I366 = make_shared<TATensor<double,2>>({active_, active_});
  auto task278 = make_shared<Task278>(Gamma38_(), I366);
  task261->add_dep(task278);
  corrq->add_task(task278);

  auto task279 = make_shared<Task279>(I366, t2);
  task278->add_dep(task279);
  corrq->add_task(task279);

  auto task280 = make_shared<Task280>(I366, t2);
  task278->add_dep(task280);
  corrq->add_task(task280);

  auto I372 = make_shared<TATensor<double,4>>({active_, active_, active_, active_});
  auto task281 = make_shared<Task281>(Gamma60_(), I372);
  task261->add_dep(task281);
  corrq->add_task(task281);

  auto task282 = make_shared<Task282>(I372, t2);
  task281->add_dep(task282);
  corrq->add_task(task282);

  return corrq;
}


#endif
