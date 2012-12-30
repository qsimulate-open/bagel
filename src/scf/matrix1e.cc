//
// BAGEL - Parallel electron correlation program.
// Filename: matrix1e.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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


#include <src/scf/matrix1e.h>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <cassert>
#include <cmath>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;


Matrix1e::Matrix1e(const shared_ptr<const Geometry> geom) : Matrix(geom->nbasis(), geom->nbasis()), geom_(geom) {
  zero();
}


Matrix1e::Matrix1e(const shared_ptr<const Geometry> geom, const int n, const int m) : Matrix(n,m), geom_(geom) {
  zero();
}


Matrix1e::Matrix1e(const Matrix1e& o) : Matrix(o.ndim_, o.mdim_), geom_(o.geom_) {
  copy_n(o.data(), ndim_*mdim_, data()); 
}


void Matrix1e::init() {

  // only lower half will be stored

  auto o0 = geom_->offsets().begin();
  int u = 0;
  for (auto a0 = geom_->atoms().begin(); a0 != geom_->atoms().end(); ++a0, ++o0) {
    // iatom1 = iatom1;
    auto offset0 = o0->begin();
    for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++offset0) {
      auto offset1 = o0->begin();
      for (auto b1 = (*a0)->shells().begin(); b1 != (*a0)->shells().end(); ++b1, ++offset1) {
        if (u++ % mpi__->size() == mpi__->rank()) {
          array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
          computebatch(input, *offset0, *offset1);
        }
      }
    }

    auto o1 = o0+1;
    for (auto a1 = a0+1; a1 != geom_->atoms().end(); ++a1, ++o1) {
      auto offset0 = o0->begin();
      for (auto b0 = (*a0)->shells().begin(); b0 != (*a0)->shells().end(); ++b0, ++offset0) {
        auto offset1 = o1->begin();
        for (auto b1 = (*a1)->shells().begin(); b1 != (*a1)->shells().end(); ++b1, ++offset1) {
          if (u++ % mpi__->size() == mpi__->rank()) {
            array<shared_ptr<const Shell>,2> input = {{*b1, *b0}};
            computebatch(input, *offset0, *offset1);
          }
        }
      }
    }
  }
  mpi__->allreduce(data_.get(), size());

}


