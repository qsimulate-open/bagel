//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: matview.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/util/math/matrix_base.h>
#include <src/util/math/matview.h>

using namespace std;
using namespace bagel;


template<typename DataType>
void MatView_<DataType>::init() {
#ifdef HAVE_SCALAPACK
  if (!localized_) {
    desc_ = mpi__->descinit(ndim(), mdim());
    localsize_ = mpi__->numroc(ndim(), mdim());
  }
#endif
}


#ifdef HAVE_SCALAPACK
template<typename DataType>
void MatView_<DataType>::setlocal_(const std::unique_ptr<DataType[]>& local) {
  zero();

  const int localrow = std::get<0>(localsize_);
  const int localcol = std::get<1>(localsize_);

  const int nblock = localrow/blocksize__;
  const int mblock = localcol/blocksize__;
  const size_t nstride = blocksize__*mpi__->nprow();
  const size_t mstride = blocksize__*mpi__->npcol();
  const int myprow = mpi__->myprow()*blocksize__;
  const int mypcol = mpi__->mypcol()*blocksize__;

  for (int i = 0; i != mblock; ++i)
    for (int j = 0; j != nblock; ++j)
      for (int id = 0; id != blocksize__; ++id)
        std::copy_n(&local[j*blocksize__+localrow*(i*blocksize__+id)], blocksize__, element_ptr(myprow+j*nstride, mypcol+i*mstride+id));

  for (int id = 0; id != localcol % blocksize__; ++id) {
    for (int j = 0; j != nblock; ++j)
      std::copy_n(&local[j*blocksize__+localrow*(mblock*blocksize__+id)], blocksize__, element_ptr(myprow+j*nstride, mypcol+mblock*mstride+id));
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


template<typename DataType>
std::unique_ptr<DataType[]> MatView_<DataType>::getlocal() const {
  const int localrow = std::get<0>(localsize_);
  const int localcol = std::get<1>(localsize_);

  std::unique_ptr<DataType[]> local(new DataType[localrow*localcol]);

  const int nblock = localrow/blocksize__;
  const int mblock = localcol/blocksize__;
  const size_t nstride = blocksize__*mpi__->nprow();
  const size_t mstride = blocksize__*mpi__->npcol();
  const int myprow = mpi__->myprow()*blocksize__;
  const int mypcol = mpi__->mypcol()*blocksize__;

  for (int i = 0; i != mblock; ++i)
    for (int j = 0; j != nblock; ++j)
      for (int id = 0; id != blocksize__; ++id)
        std::copy_n(element_ptr(myprow+j*nstride, mypcol+i*mstride+id), blocksize__, &local[j*blocksize__+localrow*(i*blocksize__+id)]);

  for (int id = 0; id != localcol % blocksize__; ++id) {
    for (int j = 0; j != nblock; ++j)
      std::copy_n(element_ptr(myprow+j*nstride, mypcol+mblock*mstride+id), blocksize__, &local[j*blocksize__+localrow*(mblock*blocksize__+id)]);
    for (int jd = 0; jd != localrow % blocksize__; ++jd)
      local[nblock*blocksize__+jd+localrow*(mblock*blocksize__+id)] = element(myprow+nblock*nstride+jd, mypcol+mblock*mstride+id);
  }
  for (int i = 0; i != mblock; ++i)
    for (int id = 0; id != blocksize__; ++id)
      for (int jd = 0; jd != localrow % blocksize__; ++jd)
        local[nblock*blocksize__+jd+localrow*(i*blocksize__+id)] = element(myprow+nblock*nstride+jd, mypcol+i*mstride+id);
  return local;
}
#endif


template class bagel::MatView_<double>;
template class bagel::MatView_<complex<double>>;
