//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: vectorb.h
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


#ifndef __SRC_MATH_BVECTOR_H
#define __SRC_MATH_BVECTOR_H

#include <src/util/math/algo.h>
#include <src/util/math/btas_interface.h>
#include <src/util/parallel/mpi_interface.h>

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
    VecView_(VecView_<DataType>& o) : btas::TensorView1<DataType>(o) { }
    VecView_(VecView_<DataType>&& o) : btas::TensorView1<DataType>(std::move(o)) { }
    VecView_(const VecView_<DataType>& o) : btas::TensorView1<DataType>(o) { }
    VecView_(btas::TensorView1<DataType>& o) : btas::TensorView1<DataType>(o) { }
    VecView_(btas::TensorView1<DataType>&& o) : btas::TensorView1<DataType>(std::move(o)) { }
    VecView_(const btas::TensorView1<DataType>& o) : btas::TensorView1<DataType>(o) { }
    VecView_(Vector_<DataType>& o) : btas::TensorView1<DataType>(o) { }
    VecView_(const Vector_<DataType>& o) : btas::TensorView1<DataType>(o) { }
    VecView_() { }
    virtual ~VecView_() { }

    size_t size() const { return this->extent(0); }
    double rms() const { return std::sqrt(detail::real(btas::dotc(*this, *this))/size()); }

    DataType* data() { /*assert(contiguous());*/ return &*begin(); }
    const DataType* data() const { /*assert(contiguous());*/ return &*cbegin(); }

    DataType& operator()(const int i) { return *(data()+i); }
    const DataType& operator()(const int i) const { return *(data()+i); }

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

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int file_version) {
      ar & boost::serialization::base_object<btas::Tensor1<DataType>>(*this);
    }

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

    DataType& operator()(size_t i) { return *(data()+i); }
    const DataType& operator()(size_t i) const { return *(data()+i); }

    VecView_<DataType> slice(const int mstart, const int mend) {
      auto low = {mstart};
      auto up  = {mend};
      assert(mstart >= 0 && mend <= size());
      return VecView_<DataType>(btas::make_rwview(this->range().slice(low, up), this->storage()));
    }

    const VecView_<DataType> slice(const int mstart, const int mend) const {
      auto low = {mstart};
      auto up  = {mend};
      assert(mstart >= 0 && mend <= size());
      return VecView_<DataType>(btas::make_rwview(this->range().slice(low, up), this->storage()));
    }

    size_t size() const { return this->extent(0); }
    double norm() const { return std::sqrt(detail::real(btas::dotc(*this, *this))); }
    double rms() const { return std::sqrt(detail::real(btas::dotc(*this, *this))/size()); }
    DataType dot_product(std::shared_ptr<const Vector_<DataType>> o) const { return blas::dot_product(data(), size(), o->data()); }
    DataType dot_product(const Vector_<DataType>& o) const { return blas::dot_product(data(), size(), o.data()); }
    void zero() { this->fill(0.0); }
    void scale(const DataType& a) { blas::scale_n(a, data(), size()); }
    void ax_plus_y(const DataType a, std::shared_ptr<const Vector_<DataType>> o) { btas::axpy(a, *o, *this); }
    void ax_plus_y(const DataType a, const Vector_<DataType>& o) { btas::axpy(a, o, *this); }

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

    void print(const std::string tag = "") const {
      if (tag != "")
        std::cout << std::endl << "  ++ " << tag << " ++" << std::endl << std::endl;

      for (int i = 0; i != size(); ++i)
        std::cout << std::setw(6) << i << std::setw(20) << std::setprecision(10) << *(data()+i) << std::endl;
    }

    void allreduce() {
      mpi__->allreduce(data(), size());
    }

};

using VectorB = Vector_<double>;
using ZVectorB = Vector_<std::complex<double>>;

}

#endif
