//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_density1qq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_density1q(const bool reset, const bool diagonal) {

  auto density1q = make_shared<Queue>();
  auto task523 = make_shared<Task523>(den1, reset);
  density1q->add_task(task523);

  auto I734 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{closed_, active_});
  auto task524 = make_shared<Task524>(den1, I734);
  task524->add_dep(task523);
  density1q->add_task(task524);

  auto task525 = make_shared<Task525>(I734, Gamma12_(), t2);
  task524->add_dep(task525);
  task525->add_dep(task523);
  density1q->add_task(task525);

  auto I736 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, closed_});
  auto task526 = make_shared<Task526>(den1, I736);
  task526->add_dep(task523);
  density1q->add_task(task526);

  auto I737 = make_shared<TATensor<double,4>>(std::vector<IndexRange>{active_, virt_, closed_, active_});
  auto task527 = make_shared<Task527>(I736, Gamma38_(), I737);
  task526->add_dep(task527);
  task527->add_dep(task523);
  density1q->add_task(task527);

  auto task528 = make_shared<Task528>(I737, t2);
  task527->add_dep(task528);
  task528->add_dep(task523);
  density1q->add_task(task528);

  auto I740 = make_shared<TATensor<double,2>>(std::vector<IndexRange>{virt_, active_});
  auto task529 = make_shared<Task529>(den1, I740);
  task529->add_dep(task523);
  density1q->add_task(task529);

  auto task530 = make_shared<Task530>(I740, Gamma60_(), t2);
  task529->add_dep(task530);
  task530->add_dep(task523);
  density1q->add_task(task530);

  return density1q;
}


#endif
