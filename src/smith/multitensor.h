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

class MultiTensor {
  protected:
    // reference coefficients
    std::vector<double> fac_;
    // external coefficients
    std::vector<std::shared_ptr<Tensor>> tensors_;

  public:
    MultiTensor() { }
    MultiTensor(const int n) : fac_(n, 0.0), tensors_(n) { }
    MultiTensor(const std::vector<double>& f, const std::vector<std::shared_ptr<Tensor>>& o) : fac_(f), tensors_(o) { }

    MultiTensor(const MultiTensor& o) : fac_(o.fac_.begin(), o.fac_.end()) {
      for (auto& i : o.tensors_)
        tensors_.push_back(i->copy());
    }

    std::shared_ptr<MultiTensor> copy() const {
      return std::make_shared<MultiTensor>(*this);
    }

    std::shared_ptr<MultiTensor> clone() const {
      auto out = std::make_shared<MultiTensor>(nref());
      auto oiter = out->tensors_.begin();
      for (auto& i : tensors_)
        *oiter++ = i->clone();
      return out;
    }

    void ax_plus_y(const double& a, const MultiTensor& o) {
      blas::ax_plus_y_n(a, o.fac_.data(), o.fac_.size(), fac_.data());
      auto oiter = o.tensors_.begin();
      for (auto& i : tensors_)
        i->ax_plus_y(a, *oiter++);
    }

    void ax_plus_y(const double& a, std::shared_ptr<const MultiTensor> o) { ax_plus_y(a, *o); }

    std::vector<std::shared_ptr<Tensor>>::iterator begin() { return tensors_.begin(); } 
    std::vector<std::shared_ptr<Tensor>>::iterator end()   { return tensors_.end(); } 
    std::vector<std::shared_ptr<Tensor>>::const_iterator begin() const { return tensors_.begin(); } 
    std::vector<std::shared_ptr<Tensor>>::const_iterator end()   const { return tensors_.end(); } 
    std::vector<std::shared_ptr<Tensor>>::const_iterator cbegin() const { return tensors_.cbegin(); } 
    std::vector<std::shared_ptr<Tensor>>::const_iterator cend()   const { return tensors_.cend(); } 

    std::shared_ptr<Tensor>& operator[](const size_t i) { assert(i < tensors_.size()); return tensors_[i]; }
    std::shared_ptr<Tensor>& at(const size_t i) { return tensors_.at(i); }
    std::shared_ptr<const Tensor> operator[](const size_t i) const { assert(i < tensors_.size()); return tensors_[i]; }
    std::shared_ptr<const Tensor> at(const size_t i) const { return tensors_.at(i); }

    std::vector<std::shared_ptr<Tensor>>& tensors() { return tensors(); }

    void scale(const double a) {
      blas::scale_n(a, fac_.data(), fac_.size());
      for (auto& i : tensors_)
        i->scale(a);
    }

    void zero() { scale(0.0); }

    size_t nref() const { assert(fac_.size() == tensors_.size()); return fac_.size(); } 

    double rms() const {
      double out = 0.0;
      size_t size = fac_.size();
      for (auto& i : fac_)
        out += i*i;
      for (auto& i : tensors_) {
        out += std::pow(i->norm(),2);
        size += i->size_alloc();
      }
      return std::sqrt(out/size);
    }

    double& fac(const int i) { return fac_[i]; }
    const double& fac(const int i) const { return fac_[i]; }
};

}
}

#endif

