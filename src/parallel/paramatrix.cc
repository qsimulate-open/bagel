//
// BAGEL - Parallel electron correlation program.
// Filename: paramatrix.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <src/parallel/paramatrix.h>
#include <config.h>
#ifdef HAVE_MPI_H
  #include <mpi.h>
#endif

using namespace std;
using namespace bagel;

void ParaMatrix::sum_reduce() {
#ifdef HAVE_MPI_H
  MPI::COMM_WORLD.Allreduce(MPI_IN_PLACE, static_cast<void*>(data_.get()), size(), MPI_DOUBLE, MPI::SUM);
#endif
}
