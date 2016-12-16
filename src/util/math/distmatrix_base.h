//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: distmatrix_base.h
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


#ifndef __SRC_UTIL_DISTMATRIX_BASE_H
#define __SRC_UTIL_DISTMATRIX_BASE_H

#include <cassert>
#include <algorithm>
#include <src/util/math/vectorb.h>
#include <src/util/parallel/scalapack.h>
#include <src/util/parallel/mpi_interface.h>

namespace bagel {

#ifdef HAVE_SCALAPACK

template<typename DataType>
class DistMatrix_base {
  protected:
    // global dimension
    int ndim_;
    int mdim_;

    // distributed data
    std::unique_ptr<DataType[]> local_;

    // Scalapack specific
    std::vector<int> desc_;
    std::tuple<int, int> localsize_;

    template<class T>
    void ax_plus_y_impl(const DataType a, const T& o) {
      assert(size() == o.size());
      blas::ax_plus_y_n(a, o.local_.get(), size(), local_.get());
    }
    template<class T>
    DataType dot_product_impl(const T& o) const {
      assert(ndim_ == o.ndim_ && mdim_ == o.mdim_);
      DataType sum = size() ? blas::dot_product(local_.get(), size(), o.local_.get()) : 0.0;
      mpi__->allreduce(&sum, 1);
      return sum;
    }

    std::pair<int, int> locate_row(const int i);
    std::pair<int, int> locate_column(const int j);

  public:
    DistMatrix_base() { }
    DistMatrix_base(const int n, const int m);
    DistMatrix_base(const DistMatrix_base& o);
    DistMatrix_base(DistMatrix_base&& o);

    virtual ~DistMatrix_base() { }

    const std::unique_ptr<DataType[]>& local() const { return local_; }
    const std::vector<int>& desc() const { return desc_; }

    size_t size() const { return std::get<0>(localsize_)*std::get<1>(localsize_); }
    int ndim() const { return ndim_; }
    int mdim() const { return mdim_; }

    virtual void diagonalize(VecView vec) = 0;

    void fill(const DataType a) { std::fill_n(local_.get(), size(), a); }
    void zero() { const DataType zero(0.0); fill(zero); }

    void add_diag(const DataType& a, const size_t start, const size_t fence);
    void ax_plus_y(const double a, const std::shared_ptr<const DistMatrix_base<DataType>> o) { ax_plus_y_impl(a, *o); }

    DataType dot_product(const std::shared_ptr<const DistMatrix_base<DataType>> o) const { return dot_product_impl(*o); }

    double norm() const { return std::sqrt(detail::real(dot_product_impl(*this))); }
    double rms() const { return norm()/std::sqrt(ndim_*mdim_); }

    void scale(const DataType a) { blas::scale_n(a, local_.get(), size()); }
    void scale(const double* vec);

};

extern template class DistMatrix_base<double>;
extern template class DistMatrix_base<std::complex<double>>;

#endif

}

#endif
