//
// BAGEL - Parallel electron correlation program.
// Filename: queue.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <ga.h>
#include <src/smith/queue.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

std::shared_ptr<Task> Queue::next_compute() {
  auto i = tasklist_.begin();
  for ( ; i != tasklist_.end(); ++i)
    if ((*i)->ready()) break;

  assert(i != tasklist_.end());
  std::shared_ptr<Task> out = *i;
  // execute
  out->compute();

  // synchronize
  GA_Sync();

  // delete dependency (to remove intermediate storages)
  for (auto& j : tasklist_) j->delete_dep(out);
  // delete this task from the queue
  tasklist_.erase(i);
  return out;
}

#endif
