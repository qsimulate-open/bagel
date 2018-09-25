//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: matrix_base.cc
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

using namespace std;
using namespace bagel;


template<typename DataType>
Matrix_base<DataType>::Matrix_base(const size_t n, const size_t m, const bool local) : btas::Tensor2<DataType>(n, m), localized_(local) {
#ifdef HAVE_SCALAPACK
  if (!localized_) {
    desc_ = mpi__->descinit(ndim(), mdim());
    localsize_ = mpi__->numroc(ndim(), mdim());
  }
#endif
  zero();
}


template<typename DataType>
Matrix_base<DataType>::Matrix_base(const Matrix_base& o) : btas::Tensor2<DataType>(o.ndim(), o.mdim()), localized_(o.localized_) {
#ifdef HAVE_SCALAPACK
  if (!localized_) {
    desc_ = mpi__->descinit(ndim(), mdim());
    localsize_ = mpi__->numroc(ndim(), mdim());
  }
#endif
  copy_n(o.data(), size(), data());
}


template<typename DataType>
Matrix_base<DataType>::Matrix_base(const MatView_<DataType>& o) : btas::Tensor2<DataType>(o.ndim(), o.mdim()), localized_(o.localized()) {
  copy_n(o.data(), o.size(), data());
#ifdef HAVE_SCALAPACK
  if (!localized_) {
    desc_ = mpi__->descinit(ndim(), mdim());
    localsize_ = mpi__->numroc(ndim(), mdim());
  }
#endif
}


template<typename DataType>
Matrix_base<DataType>::Matrix_base(Matrix_base&& o) : btas::Tensor2<DataType>(forward<Matrix_base<DataType>>(o)), localized_(o.localized_) {
#ifdef HAVE_SCALAPACK
  if (!localized_) {
    desc_ = mpi__->descinit(ndim(), mdim());
    localsize_ = mpi__->numroc(ndim(), mdim());
  }
#endif
}


template<typename DataType>
Matrix_base<DataType>& Matrix_base<DataType>::operator=(const Matrix_base<DataType>& o) {
  btas::Tensor2<DataType>::operator=(o);
  localized_ = o.localized_;
#ifdef HAVE_SCALAPACK
  if (!localized_) {
    desc_ = o.desc_;
    localsize_ = o.localsize_;
  }
#endif
  return *this;
}


template<typename DataType>
Matrix_base<DataType>& Matrix_base<DataType>::operator=(Matrix_base<DataType>&& o) {
  btas::Tensor2<DataType>::operator=(move(o));
  localized_ = o.localized_;
#ifdef HAVE_SCALAPACK
  if (!localized_) {
    desc_ = o.desc_;
    localsize_ = o.localsize_;
  }
#endif
  return *this;
}


template<typename DataType>
void Matrix_base<DataType>::fill_lower() {
  assert(ndim() == mdim());
  for (size_t i = 0; i != mdim(); ++i)
    for (size_t j = i+1; j != ndim(); ++j)
      element(j, i) = element(i, j);
}


template<typename DataType>
void Matrix_base<DataType>::fill_upper() {
  assert(ndim() == mdim());
  for (size_t i = 0; i != mdim(); ++i)
    for (size_t j = i+1; j != ndim(); ++j)
      element(i, j) = element(j, i);
}


template<typename DataType>
void Matrix_base<DataType>::fill_upper_conjg() {
  assert(ndim() == mdim());
  for (size_t i = 0; i != mdim(); ++i)
    for (size_t j = i+1; j != ndim(); ++j)
      element(i, j) = detail::conj(element(j, i));
}


template<typename DataType>
void Matrix_base<DataType>::fill_upper_negative() {
  assert(ndim() == mdim());
  for (size_t i = 0; i != mdim(); ++i) {
    assert(abs(element(i, i)) < 1e-15);
    for (size_t j = i+1; j != ndim(); ++j)
      element(i, j) = -element(j, i);
  }
}


template<typename DataType>
void Matrix_base<DataType>::symmetrize() {
  assert(ndim() == mdim());
  const size_t n = mdim();
  for (size_t i = 0; i != n; ++i)
    for (size_t j = i+1; j != n; ++j)
      element(i, j) = element(j, i) = 0.5*(element(i, j)+element(j, i));
}


template<typename DataType>
void Matrix_base<DataType>::antisymmetrize() {
  assert(ndim() == mdim());
  const size_t n = mdim();
  for (size_t i = 0; i != n; ++i) {
    for (size_t j = i; j != n; ++j) {
      element(i, j) = 0.5*(element(i, j)-element(j, i));
      element(j, i) = -element(i, j);
    }
  }
}


template<typename DataType>
void Matrix_base<DataType>::hermite() {
  assert(ndim() == mdim());
  const size_t n = mdim();
  for (size_t i = 0; i != n; ++i) {
    for (size_t j = i; j != n; ++j) {
      element(i, j) = 0.5*(element(i, j)+detail::conj(element(j, i)));
      element(j, i) = detail::conj(element(i, j));
    }
  }
}


template<typename DataType>
void Matrix_base<DataType>::copy_block(const int nstart, const int mstart, const int nsize, const int msize, const DataType* o) {
  assert(nstart >=0 && mstart >=0 && nstart + nsize <= ndim() && mstart + msize <= mdim());
  for (size_t i = mstart, j = 0; i != mstart + msize; ++i, ++j)
    copy_n(o + j*nsize, nsize, data() + nstart + i*ndim());
}


template<typename DataType>
void Matrix_base<DataType>::copy_block(const int nstart, const int mstart, const int nsize, const int msize, const btas::TensorView2<DataType> o) {
  assert(nsize == o.extent(0) && msize == o.extent(1) && o.range().ordinal().contiguous());
  copy_block(nstart, mstart, nsize, msize, &*o.begin());
}


template<typename DataType>
void Matrix_base<DataType>::add_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize, const DataType* o) {
  assert(nstart >=0 && mstart >=0 && nstart + nsize <= ndim() && mstart + msize <= mdim());
  for (size_t i = mstart, j = 0; i != mstart + msize ; ++i, ++j)
    blas::ax_plus_y_n(a, o+j*nsize, nsize, element_ptr(nstart, i));
}


template<typename DataType>
void Matrix_base<DataType>::add_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize, const btas::TensorView2<DataType> o) {
  assert(nsize == o.extent(0) && msize == o.extent(1) && o.range().ordinal().contiguous());
  add_block(a, nstart, mstart, nsize, msize, &*o.begin());
}


template<typename DataType>
void Matrix_base<DataType>::delocalize() {
  localized_ = false;
#ifdef HAVE_SCALAPACK
  desc_ = mpi__->descinit(ndim(), mdim());
  localsize_ = mpi__->numroc(ndim(), mdim());
#endif
}


template<typename DataType>
DataType Matrix_base<DataType>::trace() const {
  DataType out(0.0);
  assert(ndim() == mdim());
  for (int i = 0; i != ndim(); ++i)
    out += element(i, i);
  return out;
}


template<typename DataType>
Vector_<DataType> Matrix_base<DataType>::diag() const {
  if (ndim() != mdim()) throw logic_error("illegal call of Matrix::diag()");
  Vector_<DataType> out(ndim());
  for (int i = 0; i != ndim(); ++i)
    out[i] = element(i,i);
  return out;
}


template<typename DataType>
void Matrix_base<DataType>::add_diag(const DataType& a, const int i, const int j) {
  assert(ndim() == mdim());
  for (int ii = i; ii != j; ++ii) element(ii,ii) += a;
}


#ifdef HAVE_SCALAPACK
template<typename DataType>
void Matrix_base<DataType>::setlocal_(const unique_ptr<DataType[]>& local) {
  zero();

  const size_t localrow = get<0>(localsize_);
  const size_t localcol = get<1>(localsize_);

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


template<typename DataType>
unique_ptr<DataType[]> Matrix_base<DataType>::getlocal() const {
  const size_t localrow = get<0>(localsize_);
  const size_t localcol = get<1>(localsize_);

  unique_ptr<DataType[]> local(new DataType[localrow*localcol]);

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

#endif


template class bagel::Matrix_base<double>;
template class bagel::Matrix_base<complex<double>>;

BOOST_CLASS_EXPORT_IMPLEMENT(Matrix_base<double>)
BOOST_CLASS_EXPORT_IMPLEMENT(Matrix_base<complex<double>>)
