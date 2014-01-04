//
// BAGEL - Parallel electron correlation program.
// Filename: distmatrix_base.h
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


#ifndef __SRC_UTIL_DISTMATRIX_BASE_H
#define __SRC_UTIL_DISTMATRIX_BASE_H

#include <cassert>
#include <algorithm>
#include <src/parallel/scalapack.h>
#include <src/parallel/mpi_interface.h>

namespace bagel {

#ifdef HAVE_SCALAPACK

template<typename DataType>
class DistMatrix_base {
  protected:
    // global dimension
    const int ndim_;
    const int mdim_;

    // distributed data
    std::unique_ptr<DataType[]> local_;

    // Scalapack specific
    std::unique_ptr<int[]> desc_;
    const std::tuple<int, int> localsize_;

    template<class T>
    void ax_plus_y_impl(const DataType a, const T& o) {
      assert(size() == o.size());
      std::transform(o.local_.get(), o.local_.get()+size(), local_.get(), local_.get(),
                     [&a](DataType p, DataType q) { return a*p+q; });
    }
    template<class T>
    DataType dot_product_impl(const T& o) const {
      assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
      DataType sum = size() ? blas::dot_product(local_.get(), size(), o.local_.get()) : 0.0;
      mpi__->allreduce(&sum, 1);
      return sum;
    }

    std::pair<int, int> locate_row(const int i) { // Returns prow and local row offset for ith row
      const int rowstride = mpi__->nprow() * blocksize__;
      const int istride = i/rowstride;

      const int prow = (i - rowstride * istride) / blocksize__;
      const int off = i - rowstride * istride - prow * blocksize__ + istride * blocksize__;

      return std::make_pair(prow, off);
    }

    std::pair<int, int> locate_column(const int j) { // Returns pcol and local col offset for jth col
      const int colstride = mpi__->npcol() * blocksize__;
      const int jstride = j/colstride;

      const int pcol = (j - colstride * jstride) / blocksize__;
      const int off = j - colstride * jstride - pcol * blocksize__ + jstride * blocksize__;

      return std::make_pair(pcol, off);
    }

  public:
    DistMatrix_base(const int n, const int m) : ndim_(n), mdim_(m), desc_(mpi__->descinit(ndim_, mdim_)), localsize_(mpi__->numroc(ndim_, mdim_)) {
      local_ = std::unique_ptr<DataType[]>(new DataType[size()]);
      zero();
    }

    DistMatrix_base(const DistMatrix_base& o) : ndim_(o.ndim_), mdim_(o.mdim_), desc_(mpi__->descinit(ndim_, mdim_)), localsize_(o.localsize_) {
      local_ = std::unique_ptr<DataType[]>(new DataType[size()]);
      std::copy_n(o.local_.get(), size(), local_.get());
    }

    DistMatrix_base(DistMatrix_base&& o) : ndim_(o.ndim_), mdim_(o.mdim_), local_(std::move(o.local_)), desc_(std::move(o.desc_)), localsize_(o.localsize_) {
    }

    const std::unique_ptr<DataType[]>& local() const { return local_; }
    const std::unique_ptr<int[]>& desc() const { return desc_; }

    size_t size() const { return std::get<0>(localsize_)*std::get<1>(localsize_); }
    int ndim() const { return ndim_; }
    int mdim() const { return mdim_; }

    virtual void diagonalize(double* vec) = 0;

    void fill(const DataType a) { std::fill_n(local_.get(), size(), a); }
    void zero() { const DataType zero(0.0); fill(zero); }

    void add_diag(const DataType& a, const size_t start, const size_t fence) {
      assert(ndim_ == mdim_ && start <= fence);
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

    void ax_plus_y(const double a, const std::shared_ptr<const DistMatrix_base<DataType>> o) { ax_plus_y_impl(a, *o); }

    DataType dot_product(const std::shared_ptr<const DistMatrix_base<DataType>> o) const { return dot_product_impl(*o); }

    double norm() const { return std::sqrt(detail::real(dot_product_impl(*this))); }
    double rms() const { return norm()/std::sqrt(ndim_*mdim_); }

    void scale(const DataType a) { std::for_each(local_.get(), local_.get()+size(), [&a](DataType& p) { p *= a; }); }

    void scale(const double* vec) {
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
          for (int id = 0; id != blocksize__; ++id) {
            DataType* c = local_.get()+j*blocksize__+localrow*(i*blocksize__+id);
            std::for_each(c, c+blocksize__, [&x = vec[mypcol+i*mstride+id]](DataType& a) { a*= x; });
          }

      for (int id = 0; id != localcol % blocksize__; ++id) {
        for (int j = 0; j != nblock; ++j) {
          DataType* c = local_.get()+j*blocksize__+localrow*(mblock*blocksize__+id);
          std::for_each(c, c+blocksize__, [&x = vec[mypcol+mblock*mstride+id]](DataType& a) { a*= x; });
        }

        DataType* c = local_.get()+nblock*blocksize__+localrow*(mblock*blocksize__+id);
        std::for_each(c, c+(localrow%blocksize__), [&x = vec[mypcol+mblock*mstride+id]](DataType& a) { a*= x; });
      }
      for (int i = 0; i != mblock; ++i)
        for (int id = 0; id != blocksize__; ++id) {
          DataType* c = local_.get()+nblock*blocksize__+localrow*(i*blocksize__+id);
          std::for_each(c, c+(localrow%blocksize__), [&x = vec[mypcol+i*mstride+id]](DataType& a) { a*= x; });
        }
    }


};
#endif

}

#endif
