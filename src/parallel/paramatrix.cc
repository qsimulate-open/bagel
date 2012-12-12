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


ParaMatrix::ParaMatrix(const int n, const int m) : Matrix(n,m) {
}


void ParaMatrix::allreduce() {
  mpi__->allreduce(data_.get(), size());
}


void ParaMatrix::broadcast(const int root) {
  mpi__->broadcast(data_.get(), size(), root);
}

#ifdef HAVE_SCALAPACK
void ParaMatrix::diagonalize(double* eig) {

  if (ndim_ != mdim_) throw logic_error("illegal call of ParaMatrix::diagonalize");
  const int n = ndim_;
  int localrow, localcol;
  tie(localrow, localcol) = mpi__->numroc(n, n);

  unique_ptr<int[]> desc = mpi__->descinit(n, n);
  unique_ptr<double[]> coeff(new double[localrow*localcol]);

  unique_ptr<double[]> local = getlocal_(desc);

  int info;
  // first compute worksize
  double wsize;
  pdsyev_("V", "U", n, local.get(), desc.get(), eig, coeff.get(), desc.get(), &wsize, -1, info);

  const int lwork = round(wsize);
  unique_ptr<double[]> work(new double[lwork]);
  pdsyev_("V", "U", n, local.get(), desc.get(), eig, coeff.get(), desc.get(), work.get(), lwork, info);
  if (info) throw runtime_error("pdsyev failed in paramatrix");

  setlocal_(coeff, desc);
}


unique_ptr<double[]> ParaMatrix::getlocal_(const unique_ptr<int[]>& desc) const {
  int localrow, localcol;
  tie(localrow, localcol) = mpi__->numroc(ndim_, mdim_);
  unique_ptr<double[]> local(new double[localrow*localcol]);

  const int nblock = localrow/blocksize__;
  const int mblock = localcol/blocksize__;
  const size_t nstride = blocksize__*mpi__->nprow();
  const size_t mstride = blocksize__*mpi__->npcol();
  const int myprow = mpi__->myprow()*blocksize__;
  const int mypcol = mpi__->mypcol()*blocksize__;

  for (int i = 0; i != mblock; ++i)
    for (int j = 0; j != nblock; ++j)
      for (int id = 0; id != blocksize__; ++id) 
        copy_n(element_ptr(myprow+j*nstride, mypcol+i*mstride+id), blocksize__, &local[j*blocksize__+localrow*(i*blocksize__+id)]); 

  for (int id = 0; id != localcol % blocksize__; ++id) {
    for (int j = 0; j != nblock; ++j)
      copy_n(element_ptr(myprow+j*nstride, mypcol+mblock*mstride+id), blocksize__, &local[j*blocksize__+localrow*(mblock*blocksize__+id)]); 
    for (int jd = 0; jd != localrow % blocksize__; ++jd)
      local[nblock*blocksize__+jd+localrow*(mblock*blocksize__+id)] = element(myprow+nblock*nstride+jd, mypcol+mblock*mstride+id);
  }
  for (int i = 0; i != mblock; ++i)
    for (int id = 0; id != blocksize__; ++id) 
      for (int jd = 0; jd != localrow % blocksize__; ++jd)
        local[nblock*blocksize__+jd+localrow*(i*blocksize__+id)] = element(myprow+nblock*nstride+jd, mypcol+i*mstride+id);

  return local;
}


void ParaMatrix::setlocal_(const unique_ptr<double[]>& local, const unique_ptr<int[]>& desc) {
  zero();

  int localrow, localcol;
  tie(localrow, localcol) = mpi__->numroc(ndim_, mdim_);
  const int nblock = localrow/blocksize__;
  const int mblock = localcol/blocksize__;
  const size_t nstride = blocksize__*mpi__->nprow();
  const size_t mstride = blocksize__*mpi__->npcol();
  const int myprow = mpi__->myprow()*blocksize__;
  const int mypcol = mpi__->mypcol()*blocksize__;

  for (int i = 0; i != mblock; ++i)
    for (int j = 0; j != nblock; ++j)
      for (int id = 0; id != blocksize__; ++id) 
        copy_n(&local[j*blocksize__+localrow*(i*blocksize__+id)], blocksize__, element_ptr(myprow+j*nstride, mypcol+i*mstride+id)); 

  for (int id = 0; id != localcol % blocksize__; ++id) {
    for (int j = 0; j != nblock; ++j)
      copy_n(&local[j*blocksize__+localrow*(mblock*blocksize__+id)], blocksize__, element_ptr(myprow+j*nstride, mypcol+mblock*mstride+id)); 
    for (int jd = 0; jd != localrow % blocksize__; ++jd)
      element(myprow+nblock*nstride+jd, mypcol+mblock*mstride+id) = local[nblock*blocksize__+jd+localrow*(mblock*blocksize__+id)];
  }
  for (int i = 0; i != mblock; ++i)
    for (int id = 0; id != blocksize__; ++id) 
      for (int jd = 0; jd != localrow % blocksize__; ++jd)
        element(myprow+nblock*nstride+jd, mypcol+i*mstride+id) = local[nblock*blocksize__+jd+localrow*(i*blocksize__+id)];

  // syncronize (this can be improved, but...)
  allreduce();
}
#endif
