//
// BAGEL - Parallel electron correlation program.
// Filename: RelMRCI_sourceqq.cc
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

shared_ptr<Queue> RelMRCI::RelMRCI::make_sourceq(const bool reset, const bool diagonal) {

  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};
  auto sourceq = make_shared<Queue>();
  auto tensor46 = vector<shared_ptr<Tensor>>{s};
  auto task46 = make_shared<Task46>(tensor46, reset);
  sourceq->add_task(task46);

  shared_ptr<Tensor> I88;
  if (diagonal) {
    vector<IndexRange> I88_index = {closed_, virt_, closed_, virt_};
    I88 = make_shared<Tensor>(I88_index);
  }
  shared_ptr<Task47> task47;
  if (diagonal) {
    auto tensor47 = vector<shared_ptr<Tensor>>{s, I88};
    task47 = make_shared<Task47>(tensor47, pindex);
    task47->add_dep(task46);
    sourceq->add_task(task47);
  }

  shared_ptr<Task48> task48;
  if (diagonal) {
    auto tensor48 = vector<shared_ptr<Tensor>>{I88, v2_};
    task48 = make_shared<Task48>(tensor48, pindex);
    task47->add_dep(task48);
    task48->add_dep(task46);
    sourceq->add_task(task48);
  }

  return sourceq;
}


#endif
