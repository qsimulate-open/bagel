//
// BAGEL - Parallel electron correlation program.
// Filename: matrix1earray.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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


Matrix1eArray::Matrix1eArray(const shared_ptr<const Geometry> geom) : Matrix(geom->nbasis(), geom->nbasis()), geom_(geom) {
  
  zero();
}


Matrix1eArray::Matrix1eArray(const shared_ptr<const Geometry> geom, const int n, const int m) : Matrix(n,m), geom_(geom) {
  zero();
}


Matrix1eArray::Matrix1eArray(const Matrix1eArray& o) : Matrix(o.ndim_, o.mdim_), geom_(o.geom_) {
  copy_n(o.data(), ndim_*mdim_, data()); 
}
