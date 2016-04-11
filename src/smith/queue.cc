//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: queue.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH

#include <src/smith/queue.h>
#include <src/util/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

shared_ptr<Task> Queue::next_compute() {
  auto i = tasklist_.begin();
  for ( ; i != tasklist_.end(); ++i)
    if ((*i)->ready()) break;

  assert(i != tasklist_.end());
  shared_ptr<Task> out = *i;
  // execute
  out->compute();

  // synchronize. This only works because add_block is local...
  mpi__->barrier();

  // delete dependency (to remove intermediate storages)
  for (auto& j : tasklist_) j->delete_dep(out);
  // delete this task from the queue
  tasklist_.erase(i);
  return out;
}

#endif
