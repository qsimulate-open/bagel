//
// BAGEL - Parallel electron correlation program.
// Filename: MRCI_normqq.cc
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

shared_ptr<Queue> MRCI::MRCI::make_normq() {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto normq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor36 = {n};
  auto task36 = make_shared<Task36>(tensor36);
  normq->add_task(task36);

  vector<IndexRange> I24_index = {virt_, virt_, active_, active_};
  auto I24 = make_shared<Tensor>(I24_index);
  vector<shared_ptr<Tensor>> tensor37 = {n, I24};
  auto task37 = make_shared<Task37>(tensor37, pindex);
  task37->add_dep(task36);
  normq->add_task(task37);

  vector<IndexRange> I25_index = {active_, virt_, active_, virt_};
  auto I25 = make_shared<Tensor>(I25_index);
  vector<shared_ptr<Tensor>> tensor38 = {I24, Gamma1_(), I25};
  auto task38 = make_shared<Task38>(tensor38, pindex);
  task37->add_dep(task38);
  task38->add_dep(task36);
  normq->add_task(task38);

  vector<shared_ptr<Tensor>> tensor39 = {I25, t2};
  auto task39 = make_shared<Task39>(tensor39, pindex);
  task38->add_dep(task39);
  task39->add_dep(task36);
  normq->add_task(task39);

  return normq;
}


#endif
