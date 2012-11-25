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


#include <stdexcept>
#include <src/parallel/paramatrix.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

void ParaMatrix::allreduce() {
  mpi__->allreduce(data_.get(), size());
}


void ParaMatrix::broadcast(const int root) {
  mpi__->broadcast(data_.get(), size(), root);
}

#ifdef HAVE_SCALAPACK
void ParaMatrix::diagonalize(double* eig) {
assert(false);
  if (ndim_ != mdim_) throw logic_error("illegal call of ParaMatrix::diagonalize");
  const int n = ndim_;
  int localrow, localcol;
  tie(localrow, localcol) = mpi__->numroc(n, n);
  unique_ptr<double[]> local(new double[localrow*localcol]);
  unique_ptr<double[]> coeff(new double[localrow*localcol]);

  unique_ptr<int[]> desc = mpi__->descinit(n, n);

  for (int i = 0; i != n; ++i)
    for (int j = 0; j != n; ++j)
      pdelset_(local.get(), j+1, i+1, desc.get(), data_[j+n*i]); 

  int info;
  int lwork = n*n;
  unique_ptr<double[]> work(new double[lwork]);
  pdsyev_("V", "U", n, local.get(), desc.get(), eig, coeff.get(), desc.get(), work.get(), lwork, info);

  if (info) throw runtime_error("pdsyev failed in paramatrix");
}
#endif
