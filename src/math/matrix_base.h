//
// BAGEL - Parallel electron correlation program.
// Filename: matrix_base.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_UTIL_MATRIX_BASE_H
#define __SRC_UTIL_MATRIX_BASE_H

#include <cassert>
#include <numeric>
#include <algorithm>
#include <src/math/algo.h>
#include <src/parallel/scalapack.h>
#include <src/parallel/mpi_interface.h>

namespace bagel {

template<typename DataType>
class Matrix_base {
  protected:
    std::unique_ptr<DataType[]> data_;
    const size_t ndim_;
    const size_t mdim_;

    // if this matrix is used within node
    bool localized_;

    // for Scalapack BLAS3 operation
#ifdef HAVE_SCALAPACK
    std::unique_ptr<int[]> desc_;
    std::tuple<int, int> localsize_;

    void setlocal_(const std::unique_ptr<DataType[]>& local) {
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
#endif

    // some functions for implementation in derived classes
    template<class T>
    std::shared_ptr<T> get_submatrix_impl(const int nstart, const int mstart, const int nsize, const int msize) const {
      assert(nstart >= 0 && mstart >= 0 && nsize >= 0 && msize >= 0 && nstart+nsize <= ndim_ && mstart+msize <= mdim_);
      auto out = std::make_shared<T>(nsize, msize, localized_);
      for (int i = mstart, j = 0; i != mstart + msize ; ++i, ++j)
        std::copy_n(element_ptr(nstart, i), nsize, out->element_ptr(0, j));
      return out;
    }
    template<class T>
    std::shared_ptr<T> resize_impl(const int n, const int m) const {
    assert(n >= ndim_ && m >= mdim_);
      auto out = std::make_shared<T>(n, m, localized_);
      for (int i = 0; i != mdim_; ++i)
        std::copy_n(data()+i*ndim_, ndim_, out->data()+i*n);
      return out;
    }
    template<class T>
    std::shared_ptr<T> merge_impl(const std::shared_ptr<const T> o) const {
      assert(ndim_ == o->ndim_ && localized_ == o->localized_);
      auto out = std::make_shared<T>(ndim_, mdim_ + o->mdim_, localized_);
      std::copy_n(data_.get(), ndim_*mdim_, out->data_.get());
      std::copy_n(o->data_.get(), o->ndim_*o->mdim_, out->data_.get()+ndim_*mdim_);
      return out;
    }

    template<class T>
    void ax_plus_y_impl(const DataType& a, const T& o) {
      std::transform(o.data(), o.data()+size(), data(), data(), [&a](DataType p, DataType q){ return a*p+q; }); 
    }
    template<class T>
    DataType dot_product_impl(const T& o) const {
      return std::inner_product(data(), data()+size(), o.data(), DataType(0.0), std::plus<DataType>(), [](DataType p, DataType q){ return detail::conj(p)*q; });
    }
    template<class T>
    double orthog_impl(const std::list<std::shared_ptr<const T>> o) { 
      for (auto& it : o) {
        const DataType m = detail::conj(this->dot_product(it));
        ax_plus_y(-m, it);
      }
      const double n = norm();
      scale(1.0/n);
      return n;
    }
    template<class T>
    std::shared_ptr<T> diagonalize_blocks_impl(double* eig, std::vector<int> blocks) {
      if (!((ndim_ == mdim_) && (ndim_ == std::accumulate(blocks.begin(), blocks.end(), 0))))
        throw std::logic_error("illegal call of Matrix::diagonalize_blocks");
      auto out = std::make_shared<T>(ndim_,ndim_);
      int location = 0;
      for (auto& block_size : blocks) {
        if (block_size == 0) continue;
        auto submat = get_submatrix_impl<T>(location, location, block_size, block_size);
        submat->diagonalize(eig + location);
        out->copy_block(location, location, block_size, block_size, submat);
        location += block_size;
      }
      return out;
    }


  public:
    Matrix_base(const size_t n, const size_t m, const bool local = false) : data_(new DataType[n*m]), ndim_(n), mdim_(m), localized_(local) {
#ifdef HAVE_SCALAPACK
      if (!localized_) {
        desc_ = mpi__->descinit(ndim_, mdim_);
        localsize_ = mpi__->numroc(ndim_, mdim_);
      }
#endif
      zero();
    }

    Matrix_base(const Matrix_base& o) : data_(new DataType[o.ndim_*o.mdim_]), ndim_(o.ndim_), mdim_(o.mdim_), localized_(o.localized_) {
#ifdef HAVE_SCALAPACK
      if (!localized_) {
        desc_ = mpi__->descinit(ndim_, mdim_);
        localsize_ = mpi__->numroc(ndim_, mdim_);
      }
#endif
      std::copy_n(o.data_.get(), size(), data_.get());
    }

    Matrix_base(Matrix_base&& o) : data_(std::move(o.data_)), ndim_(o.ndim_), mdim_(o.mdim_), localized_(o.localized_) {
#ifdef HAVE_SCALAPACK
      if (!localized_) {
        desc_ = mpi__->descinit(ndim_, mdim_);
        localsize_ = mpi__->numroc(ndim_, mdim_);
      }
#endif
    }

    size_t size() const { return ndim_*mdim_; }
    int ndim() const { return ndim_; }
    int mdim() const { return mdim_; }

    virtual void fill_upper() {
      assert(ndim_ == mdim_);
      for (size_t i = 0; i != mdim_; ++i)
        for (size_t j = i+1; j != ndim_; ++j)
          data_[i+j*ndim_] = data_[j+i*ndim_];
    }

    void symmetrize() {
      assert(ndim_ == mdim_);
      const size_t n = mdim_;
      for (size_t i = 0; i != n; ++i)
        for (size_t j = i+1; j != n; ++j)
          data_[i+j*n] = data_[j+i*n] = 0.5*(data_[i+j*n]+data_[j+i*n]);
    }

    virtual void diagonalize(double* vec) = 0;

    void zero() { DataType z(0.0); fill(z); }
    void fill(const DataType a) { std::fill_n(data(), size(), a); }
    void unit() { zero(); for (int i = 0; i != ndim_; ++i) element(i,i) = DataType(1.0); assert(ndim_ == mdim_);}

    void copy_block(const int nstart, const int mstart, const int nsize, const int msize, const DataType* o) {
      for (size_t i = mstart, j = 0; i != mstart + msize; ++i, ++j)
        std::copy_n(o + j*nsize, nsize, data_.get() + nstart + i*ndim_);
    }
    void copy_block(const int nstart, const int mstart, const int nsize, const int msize, const std::shared_ptr<const Matrix_base<DataType>> o) {
      assert(nsize == o->ndim() && msize == o->mdim()); 
      copy_block(nstart, mstart, nsize, msize, o->data());
    }
    void copy_block(const int nstart, const int mstart, const int nsize, const int msize, const std::unique_ptr<DataType[]>& o) {
      copy_block(nstart, mstart, nsize, msize, o.get());
    }

    void add_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize, const DataType* o) {
      for (size_t i = mstart, j = 0; i != mstart + msize ; ++i, ++j)
        std::transform(o+j*nsize, o+(j+1)*nsize, data_.get()+nstart+i*ndim_, data_.get()+nstart+i*ndim_,
                       [&a](DataType p, DataType q){ return q + a*p; });
    }
    void add_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize, const std::shared_ptr<const Matrix_base<DataType>> o) {
      assert(nsize == o->ndim() && msize == o->mdim());
      add_block(a, nstart, mstart, nsize, msize, o->data());
    }
    void add_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize, const std::unique_ptr<DataType[]>& o) {
      add_block(a, nstart, mstart, nsize, msize, o.get());
    }

    void add_strided_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize,
                            const int ld, const DataType* o) {
      for (size_t i = mstart, j = 0; i != mstart + msize; ++i, ++j)
        std::transform(o + j*ld, o + (j+1)*ld, data_.get() + nstart + i*ndim_, data_.get() + nstart + i*ndim_,
                        [&a] (const DataType& p, const DataType& q) { return q + a*p; });
    }

    std::unique_ptr<DataType[]> get_block(const int nstart, const int mstart, const int nsize, const int msize) const {
      std::unique_ptr<DataType[]> out(new DataType[nsize*msize]);
      for (size_t i = mstart, j = 0; i != mstart + msize ; ++i, ++j)
        std::copy_n(data_.get() + nstart + i*ndim_, nsize, out.get() + j*nsize);
      return out;
    }

    DataType& operator()(const size_t& i, const size_t& j) { return data_[i+j*ndim_]; }
    const DataType& operator()(const size_t& i, const size_t& j) const { return data_[i+j*ndim_]; }

    DataType* data() { return data_.get(); }
    const DataType* data() const { return data_.get(); }

    DataType* begin() { return data_.get(); }
    DataType* end() { return data_.get() + size(); }
    const DataType* cbegin() const { return data_.get(); }
    const DataType* cend() const { return data_.get() + size(); }

    DataType& element(size_t i, size_t j) { return *element_ptr(i, j); }
    DataType* element_ptr(size_t i, size_t j) { return data()+i+j*ndim_; }
    const DataType& element(size_t i, size_t j) const { return *element_ptr(i, j); }
    const DataType* element_ptr(size_t i, size_t j) const { return data()+i+j*ndim_; }

    void ax_plus_y(const DataType a, const std::shared_ptr<const Matrix_base<DataType>> o) { ax_plus_y_impl(a, *o); } 
    DataType dot_product(const std::shared_ptr<const Matrix_base<DataType>> o) const { return dot_product_impl(*o); }

    double norm() const { return std::sqrt(detail::real(dot_product_impl(*this))); }
    double variance() const { return detail::real(dot_product_impl(*this)) / (ndim_ * mdim_); }
    double rms() const { return std::sqrt(variance()); }

    DataType trace() const {
      DataType out(0.0);
      assert(ndim_ == mdim_);
      for (int i = 0; i != ndim_; ++i)
        out += data_[i * ndim_ + i];
      return out;
    }

    void scale(const DataType& a) { std::for_each(data(), data()+size(), [&a](DataType& p){ p *= a; }); }

    void allreduce() {
      mpi__->allreduce(data_.get(), size());
    }

    void broadcast(const int root = 0) {
      mpi__->broadcast(data_.get(), size(), root);
    }

    // if we use this matrix within node, or in parallel
    void localize() { localized_ = true; }
    bool localized() const { return localized_; }

    void add_diag(const DataType& a, const int i, const int j) {
      assert(ndim_ == mdim_);
      for (int ii = i; ii != j; ++ii) element(ii,ii) += a;
    }

    void add_diag(const DataType& a) { add_diag(a,0,ndim_); }

    // returns diagonal elements
    std::unique_ptr<DataType[]> diag() const {
      if (ndim_ != mdim_) throw std::logic_error("illegal call of Matrix::diag()");
      std::unique_ptr<DataType[]> out(new DataType[ndim_]);
      for (int i = 0; i != ndim_; ++i) {
        out[i] = element(i,i);
      }
      return move(out);
    }



#ifdef HAVE_SCALAPACK
    std::unique_ptr<DataType[]> getlocal() const {
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
};

}
#endif
