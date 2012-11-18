//
// BAGEL - Parallel electron correlation program.
// Filename: dfblock.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <numeric>
#include <src/df/dfinttask.h>
#include <src/util/taskqueue.h>
#include <src/df/dfblock.h>
#include <src/rysint/libint.h>
#include <src/rysint/eribatch.h>

using namespace bagel;
using namespace std;

// protected functions
void DFBlock::ao_init() {
  // allocation of the data area
  data_ = unique_ptr<double[]>(new double[asize_*b1size_*b2size_]);

  const shared_ptr<const Shell> i3(new Shell(b1_.front()->spherical()));

  // making a task list
  vector<DFIntTask> tasks;
  tasks.reserve(b1_.size()*b2_.size()*aux_.size());

  // TODO this is not general, but for the time being I plan to have full basis functions (and limited aux basis functions); 
  assert(b1_ == b2_);
  
  auto j2 = b2off_.begin();
  for (auto& i2 : b2_) { 
    auto j1 = b1off_.begin();
    for (auto& i1 : b1_) { 
      // TODO using symmetry. This assumes that swap(i1, i2) integrals are also located in this block, which might not be the case in general.
      if (*j1 <= *j2) {
        auto j0 = aoff_.begin();
        for (auto& i0 : aux_) { 
          tasks.push_back(DFIntTask(array<shared_ptr<const Shell>,4>{{i3, i0, i1, i2}}, vector<int>{*j2, *j1, *j0}, this));
          ++j0;
        }
      }
      ++j1;
    }
    ++j2;
  }

  TaskQueue<DFIntTask> tq(tasks);
  tq.compute(resources__->max_num_threads());
}


pair<const double*, shared_ptr<RysInt> > DFBlock::compute_batch(array<shared_ptr<const Shell>,4>& input) {
#ifdef LIBINT_INTERFACE
  shared_ptr<Libint> eribatch(new Libint(input));
#else
  shared_ptr<ERIBatch> eribatch(new ERIBatch(input, 2.0));
#endif
  eribatch->compute();
  return make_pair(eribatch->data(), eribatch);
}
