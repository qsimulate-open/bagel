//
// BAGEL - Parallel electron correlation program.
// Filename: CASPT2_deciqq.cc
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

shared_ptr<Queue> CASPT2::CASPT2::make_deciq() {

  array<shared_ptr<const IndexRange>,4> cindex = {{rclosed_, ractive_, rvirt_, rci_}};

  auto deciq = make_shared<Queue>();
  vector<shared_ptr<Tensor>> tensor45 = {deci};
  auto task45 = make_shared<Task45>(tensor45);
  deciq->add_task(task45);

  vector<IndexRange> I40_index = {ci_};
  auto I40 = make_shared<Tensor>(I40_index);
  vector<shared_ptr<Tensor>> tensor46 = {deci, I40};
  auto task46 = make_shared<Task46>(tensor46, cindex);
  task46->add_dep(task45);
  deciq->add_task(task46);

  vector<IndexRange> I41_index;
  auto I41 = make_shared<Tensor>(I41_index);
  vector<shared_ptr<Tensor>> tensor47 = {I40, Gamma4_(), I41};
  auto task47 = make_shared<Task47>(tensor47, cindex);
  task46->add_dep(task47);
  task47->add_dep(task45);
  deciq->add_task(task47);

  vector<IndexRange> I42_index = {closed_, virt_, closed_, virt_};
  auto I42 = make_shared<Tensor>(I42_index);
  vector<shared_ptr<Tensor>> tensor48 = {I41, t2, I42};
  auto task48 = make_shared<Task48>(tensor48, cindex);
  task47->add_dep(task48);
  task48->add_dep(task45);
  deciq->add_task(task48);

  vector<shared_ptr<Tensor>> tensor49 = {I42, t2};
  auto task49 = make_shared<Task49>(tensor49, cindex);
  task48->add_dep(task49);
  task49->add_dep(task45);
  deciq->add_task(task49);

  vector<IndexRange> I45_index = {virt_, closed_, virt_, closed_};
  auto I45 = make_shared<Tensor>(I45_index);
  vector<shared_ptr<Tensor>> tensor50 = {I41, t2, I45};
  auto task50 = make_shared<Task50>(tensor50, cindex);
  task47->add_dep(task50);
  task50->add_dep(task45);
  deciq->add_task(task50);

  vector<shared_ptr<Tensor>> tensor51 = {I45, t2};
  auto task51 = make_shared<Task51>(tensor51, cindex);
  task50->add_dep(task51);
  task51->add_dep(task45);
  deciq->add_task(task51);

  return deciq;
}


