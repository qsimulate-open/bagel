//
// BAGEL - Parallel electron correlation program.
// Filename: multitensor.h
// Copyright (C) 2015 Toru Shiozaki
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


#ifndef __SRC_SMITH_MULTITENSOR_H
#define __SRC_SMITH_MULTITENSOR_H

#include <src/smith/tensor.h>

namespace bagel {
namespace SMITH {

template<typename DataType>
class MultiTensor_ {
  protected:
    using TType = Tensor_<DataType>;
    using Iterator = typename std::vector<std::shared_ptr<TType>>::iterator;
    using ConstIterator = typename std::vector<std::shared_ptr<TType>>::const_iterator;

    // reference coefficients
    std::vector<DataType> fac_;
    // external coefficients
    std::vector<std::shared_ptr<TType>> tensors_;

  public:
    MultiTensor_() { }
    MultiTensor_(const int n) : fac_(n, 0.0), tensors_(n) { }
    MultiTensor_(const std::vector<DataType>& f, const std::vector<std::shared_ptr<TType>>& o) : fac_(f), tensors_(o) { }

    MultiTensor_(const MultiTensor_<DataType>& o) : fac_(o.fac_.begin(), o.fac_.end()) {
      for (auto& i : o.tensors_)
        tensors_.push_back(i->copy());
    }

    std::shared_ptr<MultiTensor_<DataType>> copy() const {
      return std::make_shared<MultiTensor_<DataType>>(*this);
    }

    std::shared_ptr<MultiTensor_<DataType>> clone() const {
      auto out = std::make_shared<MultiTensor_<DataType>>(nref());
      auto oiter = out->tensors_.begin();
      for (auto& i : tensors_)
        *oiter++ = i->clone();
      return out;
    }

    void ax_plus_y(const DataType& a, const MultiTensor_<DataType>& o) {
      blas::ax_plus_y_n(a, o.fac_.data(), o.fac_.size(), fac_.data());
      auto oiter = o.tensors_.begin();
      for (auto& i : tensors_)
        i->ax_plus_y(a, *oiter++);
    }

    void ax_plus_y(const DataType& a, std::shared_ptr<const MultiTensor_<DataType>> o) { ax_plus_y(a, *o); }

    Iterator begin() { return tensors_.begin(); }
    Iterator end()   { return tensors_.end(); }
    ConstIterator begin() const { return tensors_.begin(); }
    ConstIterator end()   const { return tensors_.end(); }
    ConstIterator cbegin() const { return tensors_.cbegin(); }
    ConstIterator cend()   const { return tensors_.cend(); }

    std::shared_ptr<TType>& operator[](const size_t i) { assert(i < tensors_.size()); return tensors_[i]; }
    std::shared_ptr<TType>& at(const size_t i) { return tensors_.at(i); }
    std::shared_ptr<const TType> operator[](const size_t i) const { assert(i < tensors_.size()); return tensors_[i]; }
    std::shared_ptr<const TType> at(const size_t i) const { return tensors_.at(i); }

    std::vector<std::shared_ptr<TType>>& tensors() { return tensors(); }

    void scale(const DataType& a) {
      blas::scale_n(a, fac_.data(), fac_.size());
      for (auto& i : tensors_)
        i->scale(a);
    }

    void zero() { scale(0.0); }

    size_t nref() const { assert(fac_.size() == tensors_.size()); return fac_.size(); }
    size_t size = fac_.size();
 
    double norm() const {
      double out = 0.0;
      size_t size = fac_.size();
      for (auto& i : fac_)
        out += detail::real(detail::conj(i)*i);
      for (auto& i : tensors_) {
        out += std::pow(i->norm(),2);
        size += i->size_alloc();
      }
      return std::sqrt(out);
    }

     double rms() const { 
     size_t size = fac_.size();
     return norm() / std::sqrt(size); }

    DataType dot_product(const MultiTensor_<DataType>& o) const {
      DataType out = 0.0;
      assert(fac_.size() == o.fac_.size());
      for (int i = 0; i != fac_.size(); ++i) {
        out += detail::conj(fac_[i]) * o.fac_[i];
        out += tensors_[i]->dot_product(o.tensors_[i]);
      }
      return out;
    }

    DataType& fac(const int i) { return fac_[i]; }
    const DataType& fac(const int i) const { return fac_[i]; }
};

namespace CASPT2 { using MultiTensor = MultiTensor_<double>; }
namespace MRCI   { using MultiTensor = MultiTensor_<double>; }
namespace RelCASPT2 { using MultiTensor = MultiTensor_<std::complex<double>>; }
namespace RelMRCI   { using MultiTensor = MultiTensor_<std::complex<double>>; }

}
}

#endif

