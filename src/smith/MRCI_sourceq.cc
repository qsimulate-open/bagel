//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_sourceqq.cc
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

shared_ptr<Queue> MRCI::MRCI::make_sourceq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor32 = {s};
  auto task32 = make_shared<Task32>(tensor32);
  sourceq->add_task(task32);

  vector<IndexRange> I22_index = {virt_, virt_, active_, active_};
  auto I22 = make_shared<Tensor>(I22_index);
  vector<shared_ptr<Tensor>> tensor33 = {s, I22};
  auto task33 = make_shared<Task33>(tensor33, pindex);
  task33->add_dep(task32);
  sourceq->add_task(task33);

  vector<IndexRange> I23_index = {active_, virt_, active_, virt_};
  auto I23 = make_shared<Tensor>(I23_index);
  vector<shared_ptr<Tensor>> tensor34 = {I22, Gamma1_(), I23};
  auto task34 = make_shared<Task34>(tensor34, pindex);
  task33->add_dep(task34);
  task34->add_dep(task32);
  sourceq->add_task(task34);

  vector<shared_ptr<Tensor>> tensor35 = {I23, v2_};
  auto task35 = make_shared<Task35>(tensor35, pindex);
  task34->add_dep(task35);
  task35->add_dep(task32);
  sourceq->add_task(task35);

  return sourceq;
}


#endif
