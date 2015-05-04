//
// BAGEL - Parallel electron correlation program.
// Filename: RelCASPT2_energyqq.cc
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


#include <src/smith/RelCASPT2.h>
#include <src/smith/RelCASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelCASPT2::RelCASPT2::make_energyq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto energyq = make_shared<Queue>();
  vector<IndexRange> I30_index;
  auto I30 = make_shared<Tensor>(I30_index);
  vector<IndexRange> I31_index = {active_, active_};
  auto I31 = make_shared<Tensor>(I31_index);
  auto tensor14 = vector<shared_ptr<Tensor>>{I30, Gamma2_(), I31};
  auto task14 = make_shared<Task14>(tensor14, pindex);
  energyq->add_task(task14);

  auto tensor15 = vector<shared_ptr<Tensor>>{I31, v2_, t2};
  auto task15 = make_shared<Task15>(tensor15, pindex);
  task14->add_dep(task15);
  energyq->add_task(task15);

  auto tensor16 = vector<shared_ptr<Tensor>>{I31, t2, v2_};
  auto task16 = make_shared<Task16>(tensor16, pindex);
  task14->add_dep(task16);
  energyq->add_task(task16);

  return energyq;
}


#endif
