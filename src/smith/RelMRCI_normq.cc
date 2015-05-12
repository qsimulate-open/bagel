//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_normqq.cc
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


#include <src/smith/RelMRCI.h>
#include <src/smith/RelMRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> RelMRCI::RelMRCI::make_normq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto normq = make_shared<Queue>();
  auto tensor49 = vector<shared_ptr<Tensor>>{n};
  auto task49 = make_shared<Task49>(tensor49, reset);
  normq->add_task(task49);

  shared_ptr<Tensor> I90;
  if (diagonal) {
    vector<IndexRange> I90_index = {closed_, virt_, closed_, virt_};
    I90 = make_shared<Tensor>(I90_index);
  }
  shared_ptr<Task50> task50;
  if (diagonal) {
    auto tensor50 = vector<shared_ptr<Tensor>>{n, I90};
    task50 = make_shared<Task50>(tensor50, pindex);
    task50->add_dep(task49);
    normq->add_task(task50);
  }

  shared_ptr<Task51> task51;
  if (diagonal) {
    auto tensor51 = vector<shared_ptr<Tensor>>{I90, t2};
    task51 = make_shared<Task51>(tensor51, pindex);
    task50->add_dep(task51);
    task51->add_dep(task49);
    normq->add_task(task51);
  }

  return normq;
}


#endif
