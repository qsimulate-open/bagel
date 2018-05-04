//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: matrix_base.h
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


#ifndef __SRC_UTIL_MATRIX_BASE_H
#define __SRC_UTIL_MATRIX_BASE_H

#include <cassert>
#include <numeric>
#include <algorithm>
#include <src/util/math/algo.h>
#include <src/util/math/vectorb.h>
#include <src/util/parallel/scalapack.h>
#include <src/util/parallel/mpi_interface.h>
#include <src/util/serialization.h>

#define MATRIX_BASE
#include <src/util/math/matview.h>
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
    void setlocal_(const std::unique_ptr<DataType[]>& local);
#endif

  private:
    // serialization
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      boost::serialization::split_member(ar, *this, file_version);
    }

    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> boost::serialization::base_object<btas::Tensor2<DataType>>(*this) >> localized_;
#ifdef HAVE_SCALAPACK
      if (!localized_) {
        desc_ = mpi__->descinit(ndim(), mdim());
        localsize_ = mpi__->numroc(ndim(), mdim());
      }
#endif
    }

    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << boost::serialization::base_object<btas::Tensor2<DataType>>(*this) << localized_;
    }

  public:
    Matrix_base(const size_t n, const size_t m, const bool local = false);
    Matrix_base(const Matrix_base& o);
    Matrix_base(const MatView_<DataType>& o);
    Matrix_base(Matrix_base&& o);
    Matrix_base() : localized_(true) { }

    virtual ~Matrix_base() { }

    Matrix_base<DataType>& operator=(const Matrix_base<DataType>& o);
    Matrix_base<DataType>& operator=(Matrix_base<DataType>&& o);

    size_t size() const { return ndim()*mdim(); }
    size_t ndim() const { return this->extent(0); }
    size_t mdim() const { return this->extent(1); }

    void fill_lower();
    void fill_upper();
    void fill_upper_conjg();
    void fill_upper_negative();

    // Three functions to average upper and lower halves, enforcing a certain symmetry
    // It is recommended that you assert is_symmetric() (etc.) if using these to suppress numerical noise
    void symmetrize();
    void antisymmetrize();
    void hermite();

    // check the symmetry of a matrix
    virtual bool is_symmetric(const double thresh = 1.0e-8) const = 0;
    virtual bool is_antisymmetric(const double thresh = 1.0e-8) const = 0;
    virtual bool is_hermitian(const double thresh = 1.0e-8) const = 0;
    virtual bool is_identity(const double thresh = 1.0e-8) const = 0;

    virtual void diagonalize(VecView vec) = 0;

    void zero() { fill(0.0); }
    void fill(const DataType a) { std::fill_n(data(), size(), a); }
    void unit() { zero(); add_diag(1.0); }

    void copy_block(const int nstart, const int mstart, const int nsize, const int msize, const DataType* o);
    void copy_block(const int nstart, const int mstart, const int nsize, const int msize, const btas::TensorView2<DataType> o);
    template <typename T>
    void copy_block(const int nstart, const int mstart, const int nsize, const int msize, std::shared_ptr<T> o) { copy_block(nstart, mstart, nsize, msize, *o); }

    void add_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize, const DataType* o);
    void add_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize, const btas::TensorView2<DataType> o);
    template <typename T>
    void add_block(const DataType a, const int nstart, const int mstart, const int nsize, const int msize, std::shared_ptr<T> o) { add_block(a, nstart, mstart, nsize, msize, *o); }

    DataType& operator()(size_t i, size_t j) { return element(i, j); }
    DataType& element(size_t i, size_t j) { return *element_ptr(i, j); }
    DataType* element_ptr(size_t i, size_t j) { return data()+i+j*ndim(); }
    const DataType& operator()(size_t i, size_t j) const { return element(i, j); }
    const DataType& element(size_t i, size_t j) const { return *element_ptr(i, j); }
    const DataType* element_ptr(size_t i, size_t j) const { return data()+i+j*ndim(); }

    void ax_plus_y(const DataType a, const std::shared_ptr<const Matrix_base<DataType>> o) { ax_plus_y_impl(a, *o); }
    DataType dot_product(const std::shared_ptr<const Matrix_base<DataType>> o) const { return dot_product_impl(*o); }

    double norm() const { return std::sqrt(detail::real(dot_product_impl(*this))); }
    double variance() const { return detail::real(dot_product_impl(*this)) / (ndim() * mdim()); }
    double rms() const { return std::sqrt(variance()); }
    DataType trace() const;

    void scale(const DataType& a) { blas::scale_n(a, data(), size()); }

    void allreduce() { mpi__->allreduce(data(), size()); }
    void synchronize() { broadcast(); }
    void broadcast(const int root = 0) { mpi__->broadcast(data(), size(), root); }

    virtual void print(const std::string tag = "", const int size = 10) const { btas::print(*this, tag, size); }

    // if we use this matrix within node, or in parallel
    void delocalize();
    void localize() { localized_ = true; }
    bool localized() const { return localized_; }

    Vector_<DataType> diag() const;
    void add_diag(const DataType& a) { add_diag(a,0,ndim()); }
    void add_diag(const DataType& a, const int i, const int j);

#ifdef HAVE_SCALAPACK
    const std::vector<int>& desc() const { return desc_; }
    void setlocal(const std::unique_ptr<DataType[]>& local) { setlocal_(local); }
    std::unique_ptr<DataType[]> getlocal() const;
#endif

  protected:
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
      for (auto& it : o)
        ax_plus_y(-detail::conj(this->dot_product(it)), it);
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

    template<class T>
    std::shared_ptr<T> swap_columns_impl(const int i, const int iblock, const int j, const int jblock) const {
      assert(jblock >= 0 && iblock >= 0);
      assert(i >= 0 && j + jblock <= mdim());
      assert(i + iblock <= j);
      const int n = ndim();
      const int m = mdim();
      auto low0 = {0, 0}, low1 = {0, j}, low2 = {0, i+iblock}, low3 = {0, i}, low4 = {0, j+jblock};
      auto up0 = {n, i}, up1 = {n, j+jblock}, up2 = {n, j}, up3 = {n, i+iblock}, up4 = {n, m};
      auto out = std::make_shared<T>(ndim(), mdim(), localized_);
      out->copy_block(0,                   0, ndim(),                 i, btas::make_rwview(this->range().slice(low0, up0), this->storage()));
      out->copy_block(0,                   i, ndim(),            jblock, btas::make_rwview(this->range().slice(low1, up1), this->storage()));
      out->copy_block(0,          i + jblock, ndim(),    j - (i+iblock), btas::make_rwview(this->range().slice(low2, up2), this->storage()));
      out->copy_block(0, j + jblock - iblock, ndim(),            iblock, btas::make_rwview(this->range().slice(low3, up3), this->storage()));
      out->copy_block(0,          j + jblock, ndim(), mdim()-(j+jblock), btas::make_rwview(this->range().slice(low4, up4), this->storage()));
      return out;
    }

};

}

extern template class bagel::Matrix_base<double>;
extern template class bagel::Matrix_base<std::complex<double>>;

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::Matrix_base<double>)
BOOST_CLASS_EXPORT_KEY(bagel::Matrix_base<std::complex<double>>)

#endif
