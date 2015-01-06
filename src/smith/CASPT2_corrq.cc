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
  vector<IndexRange> I54_index;
  auto I54 = make_shared<Tensor>(I54_index);
  vector<IndexRange> I55_index = {virt_, closed_, virt_, closed_};
  auto I55 = make_shared<Tensor>(I55_index);
  vector<shared_ptr<Tensor>> tensor48 = {I54, t2, I55};
  auto task48 = make_shared<Task48>(tensor48, pindex);
  corrq->add_task(task48);

  vector<shared_ptr<Tensor>> tensor49 = {I55, t2};
  auto task49 = make_shared<Task49>(tensor49, pindex);
  task48->add_dep(task49);
  corrq->add_task(task49);

  vector<IndexRange> I57_index = {closed_, virt_, closed_, virt_};
  auto I57 = make_shared<Tensor>(I57_index);
  vector<shared_ptr<Tensor>> tensor50 = {I54, t2, I57};
  auto task50 = make_shared<Task50>(tensor50, pindex);
  task48->add_dep(task50);
  corrq->add_task(task50);

  vector<shared_ptr<Tensor>> tensor51 = {I57, t2};
  auto task51 = make_shared<Task51>(tensor51, pindex);
  task50->add_dep(task51);
  corrq->add_task(task51);

  vector<IndexRange> I59_index = {active_, active_};
  auto I59 = make_shared<Tensor>(I59_index);
  vector<shared_ptr<Tensor>> tensor52 = {I54, Gamma0_(), I59};
  auto task52 = make_shared<Task52>(tensor52, pindex);
  task48->add_dep(task52);
  corrq->add_task(task52);

  vector<IndexRange> I60_index = {virt_, closed_, virt_, active_};
  auto I60 = make_shared<Tensor>(I60_index);
  vector<shared_ptr<Tensor>> tensor53 = {I59, t2, I60};
  auto task53 = make_shared<Task53>(tensor53, pindex);
  task52->add_dep(task53);
  corrq->add_task(task53);

  vector<shared_ptr<Tensor>> tensor54 = {I60, t2};
  auto task54 = make_shared<Task54>(tensor54, pindex);
  task53->add_dep(task54);
  corrq->add_task(task54);

  vector<IndexRange> I63_index = {virt_, closed_, virt_, active_};
  auto I63 = make_shared<Tensor>(I63_index);
  vector<shared_ptr<Tensor>> tensor55 = {I59, t2, I63};
  auto task55 = make_shared<Task55>(tensor55, pindex);
  task52->add_dep(task55);
  corrq->add_task(task55);

  vector<shared_ptr<Tensor>> tensor56 = {I63, t2};
  auto task56 = make_shared<Task56>(tensor56, pindex);
  task55->add_dep(task56);
  corrq->add_task(task56);

  return corrq;
}


