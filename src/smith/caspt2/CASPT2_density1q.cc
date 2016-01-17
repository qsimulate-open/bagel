//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_density1qq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_density1q(const bool reset, const bool diagonal) {

  auto density1q = make_shared<Queue>();
  auto task506 = make_shared<Task506>(den1, reset);
  density1q->add_task(task506);

  auto I664 = make_shared<TATensor<double,2>>({closed_, active_});
  auto task507 = make_shared<Task507>(den1, I664);
  task507->add_dep(task506);
  density1q->add_task(task507);

  auto task508 = make_shared<Task508>(I664, Gamma12_(), t2);
  task507->add_dep(task508);
  task508->add_dep(task506);
  density1q->add_task(task508);

  auto I666 = make_shared<TATensor<double,2>>({virt_, closed_});
  auto task509 = make_shared<Task509>(den1, I666);
  task509->add_dep(task506);
  density1q->add_task(task509);

  auto I667 = make_shared<TATensor<double,4>>({active_, virt_, closed_, active_});
  auto task510 = make_shared<Task510>(I666, Gamma38_(), I667);
  task509->add_dep(task510);
  task510->add_dep(task506);
  density1q->add_task(task510);

  auto task511 = make_shared<Task511>(I667, t2);
  task510->add_dep(task511);
  task511->add_dep(task506);
  density1q->add_task(task511);

  auto I670 = make_shared<TATensor<double,2>>({virt_, active_});
  auto task512 = make_shared<Task512>(den1, I670);
  task512->add_dep(task506);
  density1q->add_task(task512);

  auto task513 = make_shared<Task513>(I670, Gamma60_(), t2);
  task512->add_dep(task513);
  task513->add_dep(task506);
  density1q->add_task(task513);

  return density1q;
}


#endif
