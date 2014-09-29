//
// BAGEL - Parallel electron correlation program.
// Filename: pmatrix1e.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <algorithm>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <cassert>
#include <cmath>
#include <src/parallel/mpi_interface.h>
#include <src/parallel/resources.h>
#include <src/util/taskqueue.h>
#include <src/periodic/pmatrix1e.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(PMatrix1e)

PMatrix1e::PMatrix1e(const shared_ptr<const Lattice> lattice) : PData(lattice->primitive_cell()->nbasis(), lattice->num_lattice_vectors()) { }

namespace bagel {
class PMatrix1eTask {
  protected:
    PMatrix1e* parent_;
    size_t ob0, ob1;
    array<shared_ptr<const Shell>,2> bas;
    shared_ptr<const Lattice> lattice;
    int block;
  public:
    PMatrix1eTask(array<shared_ptr<const Shell>,2> sh, size_t b0, size_t b1, shared_ptr<const Lattice> l, PMatrix1e* p, int b)
      : parent_(p), ob0(b0), ob1(b1), bas(sh), lattice(l), block(b) { }
    void compute() const { parent_->computebatch(bas, ob0, ob1, lattice, block); }
};
}

void PMatrix1e::init(shared_ptr<const Lattice> lattice) {

  shared_ptr<const Geometry> cell0 = lattice->primitive_cell();
  const size_t nshell = accumulate(cell0->atoms().begin(), cell0->atoms().end(), 0, [](int r, shared_ptr<const Atom> p) { return r+p->nshell(); });
  TaskQueue<PMatrix1eTask> task(nshell * nshell * lattice->ncell());

  int g = 0;
  int u = 0;
  /* loop over all cells in direct space */
  for (auto& disp : lattice->lattice_vectors()) {
    auto cell = make_shared<const Geometry>(*(lattice->primitive_cell()), disp);

    size_t oa0 = 0;
    for (auto a0 = cell0->atoms().begin(); a0 != cell0->atoms().end(); ++a0) { /* cell 0 */
      size_t oa1 = 0;
      for (auto a1 = cell->atoms().begin(); a1 != cell->atoms().end(); ++a1) { /* cell g */

        size_t ob0 = oa0;
        for (auto& b0 : (*a0)->shells()) {
          size_t ob1 = oa1;
          for (auto& b1 : (*a1)->shells()) {
            if (u++ % mpi__->size() == mpi__->rank())
              task.emplace_back(array<shared_ptr<const Shell>,2>{{b1, b0}}, ob0, ob1, lattice, this, g);
            ob1 += b1->nbasis();
          }
          ob0 += b0->nbasis();
        }
        oa1 += (*a1)->nbasis();
      }
      oa0 += (*a0)->nbasis();
    }
    ++g;
  }

  task.compute();
  allreduce();
}
