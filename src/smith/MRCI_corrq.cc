//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_corrqq.cc
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


#include <src/smith/MRCI.h>
#include <src/smith/MRCI_tasks.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Queue> MRCI::MRCI::make_corrq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto corrq = make_shared<Queue>();
  vector<IndexRange> I26_index;
  auto I26 = make_shared<Tensor>(I26_index);
  vector<IndexRange> I27_index = {active_, active_, active_, active_};
  auto I27 = make_shared<Tensor>(I27_index);
  vector<shared_ptr<Tensor>> tensor38 = {I26, Gamma8_(), I27};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  corrq->add_task(task38);

  vector<IndexRange> I28_index = {virt_, active_, virt_, active_};
  auto I28 = make_shared<Tensor>(I28_index);
  vector<shared_ptr<Tensor>> tensor39 = {I27, t2, I28};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  task38->add_dep(task39);
  corrq->add_task(task39);

  vector<shared_ptr<Tensor>> tensor40 = {I28, t2};
  auto task40 = make_shared<Task40>(tensor40, pindex);
  task39->add_dep(task40);
  corrq->add_task(task40);

  return corrq;
}


#endif
