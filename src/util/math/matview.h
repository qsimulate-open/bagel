//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: matview.h
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

#ifdef MATRIX_BASE
#ifndef __SRC_MATH_MATVIEW_H
#define __SRC_MATH_MATVIEW_H

#include <complex>
#include <src/util/math/btas_interface.h>
#include <src/util/math/matrix_base.h>

namespace bagel {

template <typename DataType>
class Matrix_base;

template <typename DataType>
class MatView_ : public btas::TensorView2<DataType> {
  protected:
    bool localized_;

#ifdef HAVE_SCALAPACK
    std::vector<int> desc_;
    std::tuple<int, int> localsize_;
    void setlocal_(const std::unique_ptr<DataType[]>& local);
#endif

    void init();

  public:
    MatView_(      MatView_<DataType>& o) : btas::TensorView2<DataType>(o), localized_(o.localized()) { init(); }
    MatView_(const MatView_<DataType>& o) : btas::TensorView2<DataType>(o), localized_(o.localized()) { init(); }
    MatView_(      MatView_<DataType>&& o) : btas::TensorView2<DataType>(std::move(o)), localized_(o.localized()) { init(); }
    MatView_(      btas::TensorView2<DataType>& o, const bool lo) : btas::TensorView2<DataType>(o), localized_(lo) { init(); }
    MatView_(const btas::TensorView2<DataType>& o, const bool lo) : btas::TensorView2<DataType>(o), localized_(lo) { init(); }
    MatView_(      btas::TensorView2<DataType>&& o, const bool lo) : btas::TensorView2<DataType>(std::move(o)), localized_(lo) { init(); }
    MatView_(      Matrix_base<DataType>& o) : btas::TensorView2<DataType>(o), localized_(o.localized()) { init(); }
    MatView_(const Matrix_base<DataType>& o) : btas::TensorView2<DataType>(o), localized_(o.localized()) { init(); }

    MatView_<DataType> operator=(const MatView_<DataType>& o) { localized_ = o.localized_; btas::TensorView2<DataType>::operator=(o); return *this; }
    MatView_<DataType> operator=(MatView_<DataType>&& o)      { localized_ = o.localized_; btas::TensorView2<DataType>::operator=(o); return *this; }

    int ndim() const { return this->extent(0); }
    int mdim() const { return this->extent(1); }
    size_t size() const { return ndim() * mdim(); }

    DataType* data()             { assert(contiguous()); return &*this->begin(); }
    const DataType* data() const { assert(contiguous()); return &*this->cbegin(); }

    DataType& operator()(const int i, const int j) { return element(i,j); }
    const DataType& operator()(const int i, const int j) const { return element(i,j); }

    bool localized() const { return localized_; }
    bool contiguous() const { return this->range().ordinal().contiguous(); }

    void zero() { std::fill_n(data(), size(), 0.0); }
    DataType& element(const int i, const int j) { return *element_ptr(i,j); }
    const DataType& element(const int i, const int j) const { return *element_ptr(i,j); }
    DataType* element_ptr(const int i, const int j) { return data()+i+ndim()*j; }
    const DataType* element_ptr(const int i, const int j) const { return data()+i+ndim()*j; }

    void allreduce() { assert(!localized_); mpi__->allreduce(data(), size()); }

#ifdef HAVE_SCALAPACK
    const std::vector<int>& desc() const { return desc_; }
    void setlocal(const std::unique_ptr<DataType[]>& local) { setlocal_(local); }
    std::unique_ptr<DataType[]> getlocal() const;
#endif
};

using MatView = MatView_<double>;
using ZMatView = MatView_<std::complex<double>>;

extern template class MatView_<double>;
extern template class MatView_<std::complex<double>>;

}

#endif
#endif
