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


#include <src/smith/CASPT2.h>
#include <src/smith/CASPT2_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> CASPT2::CASPT2::make_corrq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto corrq = make_shared<Queue>();
  vector<IndexRange> I42_index;
  auto I42 = make_shared<Tensor>(I42_index);
  vector<IndexRange> I43_index = {closed_, virt_, closed_, virt_};
  auto I43 = make_shared<Tensor>(I43_index);
  vector<shared_ptr<Tensor>> tensor40 = {I42, t2, I43};
  auto task40 = make_shared<Task40>(tensor40, pindex);
  corrq->add_task(task40);

  vector<shared_ptr<Tensor>> tensor41 = {I43, t2};
  auto task41 = make_shared<Task41>(tensor41, pindex);
  task40->add_dep(task41);
  corrq->add_task(task41);

  return corrq;
}


