//
// BAGEL - Parallel electron correlation program.
// Filename: vectorb.h
// Copyright (C) 2014 Toru Shiozaki
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


#ifndef __SRC_MATH_BVECTOR_H
#define __SRC_MATH_BVECTOR_H

#include <src/math/btas_interface.h>
#include <complex>

namespace bagel {

template<typename DataType>
class Vector_;

template<typename DataType>
class VecView_ : public btas::TensorView1<DataType> {
  public:
    using data_type = DataType;
    using btas::TensorView1<DataType>::begin;
    using btas::TensorView1<DataType>::cbegin;
    using btas::TensorView1<DataType>::end;
    using btas::TensorView1<DataType>::cend;

  public:
    VecView_(const VecView_<DataType>& o) : btas::TensorView1<DataType>(o) { }
    VecView_(const btas::TensorView1<DataType>& o) : btas::TensorView1<DataType>(o) { }
    VecView_(VecView_&& o) : btas::TensorView1<DataType>(std::move(o)) { }
    VecView_(Vector_<DataType>& o) : btas::TensorView1<DataType>(o) { }
    VecView_(const Vector_<DataType>& o) : btas::TensorView1<DataType>(o) { }
    VecView_() { }
    virtual ~VecView_() { }

    size_t size() const { return this->storage().size(); }
    double rms() const { return std::sqrt(detail::real(btas::dotc(*this, *this))/size()); }

    DataType* data() { assert(contiguous()); return &*begin(); }
    const DataType* data() const { assert(contiguous()); return &*cbegin(); }

    bool contiguous() const { return this->range().ordinal().contiguous(); }

    template<typename U = DataType, class = typename std::enable_if<std::is_same<std::complex<double>,U>::value>::type>
    std::shared_ptr<Vector_<double>> get_real_part() const {
      assert(contiguous());
      auto out = std::make_shared<Vector_<double>>(size());
      auto iter = out->begin();
      for (auto& i : *this) *iter++ = std::real(i);
      return out;
    }

    template<typename U = DataType, class = typename std::enable_if<std::is_same<std::complex<double>,U>::value>::type>
    std::shared_ptr<Vector_<double>> get_imag_part() const {
      assert(contiguous());
      auto out = std::make_shared<Vector_<double>>(size());
      auto iter = out->begin();
      for (auto& i : *this) *iter++ = std::imag(i);
      return out;
    }
};

using VecView = VecView_<double>;
using ZVecView = VecView_<std::complex<double>>;

template<typename DataType>
class Vector_ : public btas::Tensor1<DataType> {
  public:
    using data_type = DataType;
    using btas::Tensor1<DataType>::data;
    using btas::Tensor1<DataType>::begin;
    using btas::Tensor1<DataType>::cbegin;
    using btas::Tensor1<DataType>::end;
    using btas::Tensor1<DataType>::cend;

  public:
    Vector_(const size_t n) : btas::Tensor1<DataType>(n) { this->fill(0.0); }
    Vector_(const Vector_<DataType>& o) : btas::Tensor1<DataType>(o) { }
    Vector_(const btas::TensorView1<DataType>& o) : btas::Tensor1<DataType>(o) { }
    Vector_(btas::TensorView1<DataType>&& o) : btas::Tensor1<DataType>(std::move(o)) { }
    Vector_(Vector_&& o) : btas::Tensor1<DataType>(std::move(o)) { }
    Vector_() { }
    virtual ~Vector_() { }

    template<typename U = DataType, class = typename std::enable_if<std::is_same<std::complex<double>,U>::value>::type>
    Vector_(const Vector_<double>& r, const Vector_<double>& i) : btas::Tensor1<std::complex<double>>(r.size()) {
      assert(r.size() == i.size());
      auto riter = r.begin(); auto iiter = i.begin();
      for (auto& i : *this)
        i = *riter++ + *iiter++ * (std::complex<double>(0.0,1.0));
    }

    std::shared_ptr<Vector_<DataType>> clone() const { return std::make_shared<Vector_<DataType>>(size()); }
    std::shared_ptr<Vector_<DataType>> copy()  const { return std::make_shared<Vector_<DataType>>(*this); }

    std::shared_ptr<VecView_<DataType>> slice(const int mstart, const int mend) const {
      auto low = {mstart};
      auto up  = {mend};
      assert(mstart >= 0 && mend <= size());
      btas::TensorView1<DataType> tmp(this->range().slice(low, up), this->storage());
      return std::make_shared<VecView_<DataType>>(std::move(tmp));
    }

    size_t size() const { return this->storage().size(); }
    double rms() const { return std::sqrt(detail::real(btas::dotc(*this, *this))/size()); }

    Vector_<DataType>& operator=(const Vector_<DataType>& o) { btas::Tensor1<DataType>::operator=(o);            return *this; }
    Vector_<DataType>& operator=(Vector_<DataType>&& o)      { btas::Tensor1<DataType>::operator=(std::move(o)); return *this; }

    template<typename U = DataType, class = typename std::enable_if<std::is_same<std::complex<double>,U>::value>::type>
    std::shared_ptr<Vector_<double>> get_real_part() const {
      auto out = std::make_shared<Vector_<double>>(size());
      auto iter = out->begin();
      for (auto& i : *this) *iter++ = std::real(i);
      return out;
    }

    template<typename U = DataType, class = typename std::enable_if<std::is_same<std::complex<double>,U>::value>::type>
    std::shared_ptr<Vector_<double>> get_imag_part() const {
      auto out = std::make_shared<Vector_<double>>(size());
      auto iter = out->begin();
      for (auto& i : *this) *iter++ = std::imag(i);
      return out;
    }

    void allreduce() {
      mpi__->allreduce(data(), size());
    }

};

using VectorB = Vector_<double>;
using ZVectorB = Vector_<std::complex<double>>;

}

#endif
