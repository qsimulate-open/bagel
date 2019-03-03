//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: distmatrix_base.cc
// Copyright (C) 2009 Toru Shiozaki
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

#include <src/util/math/distmatrix_base.h>

using namespace std;
using namespace bagel;

#ifdef HAVE_SCALAPACK

template<typename DataType>
DistMatrix_base<DataType>::DistMatrix_base(const int n, const int m) : ndim_(n), mdim_(m), desc_(mpi__->descinit(ndim_, mdim_)), localsize_(mpi__->numroc(ndim_, mdim_)) {
  local_ = unique_ptr<DataType[]>(new DataType[size()]);
  zero();
}


template<typename DataType>
DistMatrix_base<DataType>::DistMatrix_base(const DistMatrix_base& o) : ndim_(o.ndim_), mdim_(o.mdim_), desc_(mpi__->descinit(ndim_, mdim_)), localsize_(o.localsize_) {
  local_ = unique_ptr<DataType[]>(new DataType[size()]);
  copy_n(o.local_.get(), size(), local_.get());
}


template<typename DataType>
DistMatrix_base<DataType>::DistMatrix_base(DistMatrix_base&& o) : ndim_(o.ndim_), mdim_(o.mdim_), local_(move(o.local_)), desc_(move(o.desc_)), localsize_(o.localsize_) {
}


template<typename DataType>
void DistMatrix_base<DataType>::add_diag(const DataType& a, const size_t start, const size_t fence) {
  assert(ndim_ == mdim_ && start <= fence);
  const int localrow = get<0>(localsize_);
  const int localcol = get<1>(localsize_);

  const int nblock = localrow/blocksize__;
  const int mblock = localcol/blocksize__;
  const size_t nstride = blocksize__*mpi__->nprow();
  const size_t mstride = blocksize__*mpi__->npcol();
  const int myprow = mpi__->myprow()*blocksize__;
  const int mypcol = mpi__->mypcol()*blocksize__;

  for (int i = 0; i != mblock; ++i)
    for (int j = 0; j != nblock; ++j)
      for (int id = 0; id != blocksize__; ++id)
        for (int jd = 0; jd != blocksize__; ++jd) {
          const size_t m = mypcol+i*mstride+id;
          const size_t n = myprow+j*nstride+jd;
          if (m == n && start <= n && n < fence)
            local_[j*blocksize__+jd+localrow*(i*blocksize__+id)] += a;
        }
  for (int id = 0; id != localcol % blocksize__; ++id)
    for (int jd = 0; jd != localrow % blocksize__; ++jd) {
      const size_t m = mblock*blocksize__+id;
      const size_t n = nblock*blocksize__+jd;
      if (m == n && start <= n && n < fence)
        local_[nblock*blocksize__+jd+localrow*(mblock*blocksize__+id)] += a;
    }
}


template<typename DataType>
void DistMatrix_base<DataType>::scale(const double* vec) {
  const int localrow = get<0>(localsize_);
  const int localcol = get<1>(localsize_);

  const int nblock = localrow/blocksize__;
  const int mblock = localcol/blocksize__;
  const size_t mstride = blocksize__*mpi__->npcol();
  const int mypcol = mpi__->mypcol()*blocksize__;
  for (int i = 0; i != mblock; ++i)
    for (int j = 0; j != nblock; ++j)
      for (int id = 0; id != blocksize__; ++id) {
        DataType* c = local_.get()+j*blocksize__+localrow*(i*blocksize__+id);
        blas::scale_n(vec[mypcol+i*mstride+id], c, blocksize__);
      }

  for (int id = 0; id != localcol % blocksize__; ++id) {
    for (int j = 0; j != nblock; ++j) {
      DataType* c = local_.get()+j*blocksize__+localrow*(mblock*blocksize__+id);
      blas::scale_n(vec[mypcol+mblock*mstride+id], c, blocksize__);
    }
    DataType* c = local_.get()+nblock*blocksize__+localrow*(mblock*blocksize__+id);
    blas::scale_n(vec[mypcol+mblock*mstride+id], c, localrow%blocksize__);
  }
  for (int i = 0; i != mblock; ++i)
    for (int id = 0; id != blocksize__; ++id) {
      DataType* c = local_.get()+nblock*blocksize__+localrow*(i*blocksize__+id);
      blas::scale_n(vec[mypcol+i*mstride+id], c, localrow%blocksize__);
    }
}


template<typename DataType>
pair<int, int> DistMatrix_base<DataType>::locate_row(const int i) { // Returns prow and local row offset for ith row
  const int rowstride = mpi__->nprow() * blocksize__;
  const int istride = i/rowstride;

  const int prow = (i - rowstride * istride) / blocksize__;
  const int off = i - rowstride * istride - prow * blocksize__ + istride * blocksize__;

  return {prow, off};
}


template<typename DataType>
pair<int, int> DistMatrix_base<DataType>::locate_column(const int j) { // Returns pcol and local col offset for jth col
  const int colstride = mpi__->npcol() * blocksize__;
  const int jstride = j/colstride;

  const int pcol = (j - colstride * jstride) / blocksize__;
  const int off = j - colstride * jstride - pcol * blocksize__ + jstride * blocksize__;

  return {pcol, off};
}


template class bagel::DistMatrix_base<double>;
template class bagel::DistMatrix_base<complex<double>>;

#endif
