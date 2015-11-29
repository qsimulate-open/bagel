//
// BAGEL - Parallel electron correlation program.
// Filename: extrap.h
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


#ifndef __SRC_SMITH_AMPLITUDE_H
#define __SRC_SMITH_AMPLITUDE_H

#include <src/smith/multitensor.h>
#include <src/smith/spinfreebase.h>

namespace bagel {
namespace SMITH {

// DIIS utility classes
template<typename DataType>
class Amplitude;

template<typename DataType>
class Residual {
  friend class Amplitude<DataType>;
  protected:
    std::shared_ptr<MultiTATensor<DataType,4>> res_;
    SpinFreeMethod<DataType>* me_;

  public:
    Residual(std::shared_ptr<const MultiTATensor<DataType,4>> r, SpinFreeMethod<DataType>* m) : res_(r->copy()), me_(m) { }
    Residual(std::shared_ptr<MultiTATensor<DataType,4>>&& r, SpinFreeMethod<DataType>* m) : res_(std::move(r)), me_(m) { }

    std::shared_ptr<Residual<DataType>> clone() const {
      auto out = std::make_shared<Residual<DataType>>(res_, me_);
      out->res_->zero();
      return out;
    }
    void synchronize() { }

    void ax_plus_y(const DataType& a, const Residual<DataType>& o) { res_->ax_plus_y(a, o.res_); }
    void ax_plus_y(const DataType& a, std::shared_ptr<const Residual<DataType>> o) { ax_plus_y(a, *o); }

    void ax_plus_y(const DataType& a, const Amplitude<DataType>& o);
    void ax_plus_y(const DataType& a, std::shared_ptr<const Amplitude<DataType>> o);

    DataType dot_product(const Amplitude<DataType>& o) const;
    DataType dot_product(std::shared_ptr<const Amplitude<DataType>> o) const;

    std::shared_ptr<MultiTATensor<DataType,4>> tensor() { return res_; }
    std::shared_ptr<const MultiTATensor<DataType,4>> tensor() const { return res_; }
};

template<typename DataType>
class Amplitude {
  friend class Residual<DataType>;
  protected:
    std::shared_ptr<MultiTATensor<DataType,4>> amp_;
    // the result of <ab/ij|\phi> so that one can easily compute the overlap
    std::shared_ptr<MultiTATensor<DataType,4>> left_;

    SpinFreeMethod<DataType>* me_;

  public:
    Amplitude(std::shared_ptr<const MultiTATensor<DataType,4>> t, std::shared_ptr<const MultiTATensor<DataType,4>> l, SpinFreeMethod<DataType>* m)
     : amp_(t->copy()), left_(l->copy()), me_(m) { }

    std::shared_ptr<Amplitude<DataType>> clone() const {
      auto out = std::make_shared<Amplitude<DataType>>(amp_, left_, me_);
      out->amp_->zero();
      out->left_->zero();
      return out;
    }
    void synchronize() { }

    void ax_plus_y(const DataType& a, const Amplitude<DataType>& o) { amp_->ax_plus_y(a, o.amp_); left_->ax_plus_y(a, o.left_); }
    void ax_plus_y(const DataType& a, std::shared_ptr<const Amplitude<DataType>> o) { ax_plus_y(a, *o); }

    DataType dot_product(const Amplitude<DataType>& o) const { return detail::conj(me_->dot_product_transpose(o.left_, amp_)); }
    DataType dot_product(std::shared_ptr<const Amplitude<DataType>> o) const { return dot_product(*o); }

    DataType dot_product(const Residual<DataType>& o) const { return detail::conj(me_->dot_product_transpose(o.res_, amp_)); }
    DataType dot_product(std::shared_ptr<const Residual<DataType>> o) const { return dot_product(*o); }

    std::shared_ptr<const MultiTATensor<DataType,4>> tensor() const { return amp_; }
    std::shared_ptr<const MultiTATensor<DataType,4>> left() const { return left_; }
};

template<typename DataType>
void Residual<DataType>::ax_plus_y(const DataType& a, const Amplitude<DataType>& o) { res_->ax_plus_y(a, o.left_); }
template<typename DataType>
void Residual<DataType>::ax_plus_y(const DataType& a, std::shared_ptr<const Amplitude<DataType>> o) { ax_plus_y(a, *o); }
template<typename DataType>
DataType Residual<DataType>::dot_product(const Amplitude<DataType>& o) const { return me_->dot_product_transpose(res_, o.amp_); }
template<typename DataType>
DataType Residual<DataType>::dot_product(std::shared_ptr<const Amplitude<DataType>> o) const { return dot_product(*o); }

}
}

#endif
