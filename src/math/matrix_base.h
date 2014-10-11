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
#include <src/util/serialization.h>
#include <src/math/vectorb.h>

#define MATRIX_BASE
#include <src/math/matview.h>
#undef MATRIX_BASE

namespace bagel {

template<typename DataType>
class Matrix_base : public btas::Tensor2<DataType> {
  public:
    using data_type = DataType;
    using btas::Tensor2<DataType>::data;
    using btas::Tensor2<DataType>::begin;
    using btas::Tensor2<DataType>::cbegin;
    using btas::Tensor2<DataType>::end;
    using btas::Tensor2<DataType>::cend;

  protected:
    // if this matrix is used within node
    bool localized_;

    // for Scalapack BLAS3 operation
#ifdef HAVE_SCALAPACK
    std::vector<int> desc_;
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
      assert(nstart >= 0 && mstart >= 0 && nsize >= 0 && msize >= 0 && nstart+nsize <= ndim() && mstart+msize <= mdim());
      auto out = std::make_shared<T>(nsize, msize, localized_);
      for (int i = mstart, j = 0; i != mstart + msize ; ++i, ++j)
        std::copy_n(element_ptr(nstart, i), nsize, out->element_ptr(0, j));
      return out;
    }
    template<class T>
    std::shared_ptr<T> resize_impl(const int n, const int m) const {
    assert(n >= ndim() && m >= mdim());
      auto out = std::make_shared<T>(n, m, localized_);
      for (int i = 0; i != mdim(); ++i)
        std::copy_n(data()+i*ndim(), ndim(), out->data()+i*n);
      return out;
    }
    template<class T>
    std::shared_ptr<T> merge_impl(const std::shared_ptr<const T> o) const {
      assert(ndim() == o->ndim() && localized_ == o->localized_);
      auto out = std::make_shared<T>(ndim(), mdim() + o->mdim(), localized_);
      std::copy_n(data(), ndim()*mdim(), out->data());
      std::copy_n(o->data(), o->ndim()*o->mdim(), out->data()+ndim()*mdim());
      return out;
    }

    template<class T>
    void ax_plus_y_impl(const DataType& a, const T& o) { btas::axpy(a, o, *this); }

    template<class T>
    DataType dot_product_impl(const T& o) const {
      return blas::dot_product(data(), size(), o.data());
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
    std::shared_ptr<T> diagonalize_blocks_impl(VectorB& eig, std::vector<int> blocks) {
      if (!((ndim() == mdim()) && (ndim() == std::accumulate(blocks.begin(), blocks.end(), 0))))
        throw std::logic_error("illegal call of Matrix::diagonalize_blocks");
      assert(eig.size() >= ndim());
      auto out = std::make_shared<T>(ndim(),ndim());
      int location = 0;
      for (auto& block_size : blocks) {
        if (block_size == 0) continue;
        auto submat = get_submatrix_impl<T>(location, location, block_size, block_size);
        submat->diagonalize(eig.slice(location, location+block_size));
        out->copy_block(location, location, block_size, block_size, submat);
        location += block_size;
      }
      return out;
    }

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      ar & boost::serialization::base_object<btas::Tensor2<DataType>>(*this) & localized_;
#ifdef HAVE_SCALAPACK
      ar & desc_ & localsize_;
#endif
    }

  public:
    Matrix_base(const size_t n, const size_t m, const bool local = false) : btas::Tensor2<DataType>(n, m), localized_(local) {
#ifdef HAVE_SCALAPACK
      if (!localized_) {
        desc_ = mpi__->descinit(ndim(), mdim());
        localsize_ = mpi__->numroc(ndim(), mdim());
      }
#endif
      zero();
    }

    Matrix_base(const Matrix_base& o) : btas::Tensor2<DataType>(o.ndim(), o.mdim()), localized_(o.localized_) {
#ifdef HAVE_SCALAPACK
      if (!localized_) {
        desc_ = mpi__->descinit(ndim(), mdim());
        localsize_ = mpi__->numroc(ndim(), mdim());
      }
#endif
      std::copy_n(o.data(), size(), data());
    }

    Matrix_base(const MatView_<DataType>& o) : btas::Tensor2<DataType>(o.ndim(), o.mdim()), localized_(o.localized()) {
      std::copy_n(o.data(), o.size(), data());
#ifdef HAVE_SCALAPACK
      if (!localized_) {
        desc_ = mpi__->descinit(ndim(), mdim());
        localsize_ = mpi__->numroc(ndim(), mdim());
      }
#endif
    }

    Matrix_base(Matrix_base&& o) : btas::Tensor2<DataType>(std::forward<Matrix_base<DataType>>(o)), localized_(o.localized_) {
#ifdef HAVE_SCALAPACK
      if (!localized_) {
        desc_ = mpi__->descinit(ndim(), mdim());
        localsize_ = mpi__->numroc(ndim(), mdim());
      }
#endif
    }

    Matrix_base() : localized_(true) { }

    virtual ~Matrix_base() { }

    Matrix_base<DataType>& operator=(const Matrix_base<DataType>& o) {
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

    Matrix_base<DataType>& operator=(Matrix_base<DataType>&& o) {
      btas::Tensor2<DataType>::operator=(std::move(o));
      localized_ = o.localized_;
#ifdef HAVE_SCALAPACK
      if (!localized_) {
        desc_ = o.desc_;
        localsize_ = o.localsize_;
      }
#endif
      return *this;
    }

    size_t size() const { return ndim()*mdim(); }
    int ndim() const { return this->extent(0); }
    int mdim() const { return this->extent(1); }

    void fill_upper() {
      assert(ndim() == mdim());
      for (size_t i = 0; i != mdim(); ++i)
        for (size_t j = i+1; j != ndim(); ++j)
          element(i, j) = element(j, i);
    }

    void fill_upper_conjg() {
      assert(ndim() == mdim());
      for (size_t i = 0; i != mdim(); ++i)
        for (size_t j = i+1; j != ndim(); ++j)
          element(i, j) = detail::conj(element(j, i));
    }

    virtual void fill_upper_negative() {
      assert(ndim() == mdim());
      for (size_t i = 0; i != mdim(); ++i) {
        assert(abs(element(i, i)) < 1e-15);
        for (size_t j = i+1; j != ndim(); ++j)
          element(i, j) = -element(j, i);
      }
    }

    // Three functions to average upper and lower halves, enforcing a certain symmetry
    // It is recommended that you assert is_symmetric() (etc.) if using these to suppress numerical noise
    void symmetrize() {
      assert(ndim() == mdim());
      const size_t n = mdim();
      for (size_t i = 0; i != n; ++i)
        for (size_t j = i+1; j != n; ++j)
          element(i, j) = element(j, i) = 0.5*(element(i, j)+element(j, i));
    }

    void antisymmetrize() {
      assert(ndim() == mdim());
      const size_t n = mdim();
      for (size_t i = 0; i != n; ++i) {
        for (size_t j = i; j != n; ++j) {
          element(i, j) = 0.5*(element(i, j)-element(j, i));
          element(j, i) = -element(i, j);
        }
      }
    }

    void hermite() {
      assert(ndim() == mdim());
      const size_t n = mdim();
      for (size_t i = 0; i != n; ++i) {
        for (size_t j = i; j != n; ++j) {
          element(i, j) = 0.5*(element(i, j)+detail::conj(element(j, i)));
          element(j, i) = detail::conj(element(i, j));
        }
      }
    }

    // check the symmetry of a matrix
    virtual bool is_symmetric(const double thresh = 1.0e-8) const = 0;
    virtual bool is_antisymmetric(const double thresh = 1.0e-8) const = 0;
    virtual bool is_hermitian(const double thresh = 1.0e-8) const = 0;
    virtual bool is_identity(const double thresh = 1.0e-8) const = 0;

    virtual void diagonalize(VecView vec) = 0;

    void zero() { DataType z(0.0); fill(z); }
    void fill(const DataType a) { std::fill_n(data(), size(), a); }
    void unit() { zero(); for (int i = 0; i != ndim(); ++i) element(i,i) = DataType(1.0); assert(ndim() == mdim());}

    void copy_block(const int nstart, const int mstart, const int nsize, const int msize, const DataType* o) {
      for (size_t i = mstart, j = 0; i != mstart + msize; ++i, ++j)
        std::copy_n(o + j*nsize, nsize, data() + nstart + i*ndim());
    }
    void copy_block(const int nstart, const int mstart, const int nsize, const int msize, const btas::TensorView2<DataType> o) {
      assert(nsize == o.extent(0) && msize == o.extent(1) && o.range().ordinal().contiguous());
      copy_block(nstart, mstart, nsize, msize, &*o.begin());
    }
    template <typename T>
    void copy_block(const int nstart, const int mstart, const int nsize, const int msize, std::shared_ptr<T> o) {
      copy_block(nstart, mstart, nsize, msize, *o);
    }

    void add_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize, const DataType* o) {
      for (size_t i = mstart, j = 0; i != mstart + msize ; ++i, ++j)
        blas::ax_plus_y_n(a, o+j*nsize, nsize, element_ptr(nstart, i));
    }
    void add_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize, const btas::TensorView2<DataType> o) {
      assert(nsize == o.extent(0) && msize == o.extent(1) && o.range().ordinal().contiguous());
      add_block(a, nstart, mstart, nsize, msize, &*o.begin());
    }
    template <typename T>
    void add_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize, std::shared_ptr<T> o) {
      add_block(a, nstart, mstart, nsize, msize, *o);
    }

    void add_strided_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize,
                            const int ld, const DataType* o) {
      for (size_t i = mstart, j = 0; i != mstart + msize; ++i, ++j)
        blas::ax_plus_y_n(a, o+j*ld, nsize, element_ptr(nstart, i));
    }

    DataType& element(size_t i, size_t j) { return (*this)(i, j); }
    DataType* element_ptr(size_t i, size_t j) { return data()+i+j*ndim(); }
    const DataType& element(size_t i, size_t j) const { return (*this)(i, j); }
    const DataType* element_ptr(size_t i, size_t j) const { return data()+i+j*ndim(); }

    void ax_plus_y(const DataType a, const std::shared_ptr<const Matrix_base<DataType>> o) { ax_plus_y_impl(a, *o); }
    DataType dot_product(const std::shared_ptr<const Matrix_base<DataType>> o) const { return dot_product_impl(*o); }

    double norm() const { return std::sqrt(detail::real(dot_product_impl(*this))); }
    double variance() const { return detail::real(dot_product_impl(*this)) / (ndim() * mdim()); }
    double rms() const { return std::sqrt(variance()); }

    DataType trace() const {
      DataType out(0.0);
      assert(ndim() == mdim());
      for (int i = 0; i != ndim(); ++i)
        out += element(i, i);
      return out;
    }

    void scale(const DataType& a) { std::for_each(data(), data()+size(), [&a](DataType& p){ p *= a; }); }

    void allreduce() {
      mpi__->allreduce(data(), size());
    }

    void broadcast(const int root = 0) {
      mpi__->broadcast(data(), size(), root);
    }

    void synchronize() {
      broadcast();
    }

    virtual void print(const std::string tag = "", const int size = 10) const { btas::print(*this, tag, size); }

    // if we use this matrix within node, or in parallel
    void delocalize() { localized_ = false;
#ifdef HAVE_SCALAPACK
      desc_ = mpi__->descinit(ndim(), mdim());
      localsize_ = mpi__->numroc(ndim(), mdim());
#endif
    }
    void localize() { localized_ = true; }
    bool localized() const { return localized_; }

    void add_diag(const DataType& a, const int i, const int j) {
      assert(ndim() == mdim());
      for (int ii = i; ii != j; ++ii) element(ii,ii) += a;
    }

    void add_diag(const DataType& a) { add_diag(a,0,ndim()); }

    // returns diagonal elements
    std::unique_ptr<DataType[]> diag() const {
      if (ndim() != mdim()) throw std::logic_error("illegal call of Matrix::diag()");
      std::unique_ptr<DataType[]> out(new DataType[ndim()]);
      for (int i = 0; i != ndim(); ++i) {
        out[i] = element(i,i);
      }
      return move(out);
    }



#ifdef HAVE_SCALAPACK
    const std::vector<int>& desc() const { return desc_; }
    void setlocal(const std::unique_ptr<DataType[]>& local) { setlocal_(local); }

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

extern template class bagel::Matrix_base<double>;
extern template class bagel::Matrix_base<std::complex<double>>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Matrix_base<double>)
BOOST_CLASS_EXPORT_KEY(bagel::Matrix_base<std::complex<double>>)

#endif
