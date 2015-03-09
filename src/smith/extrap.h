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

#include <src/smith/tensor.h>
#include <src/smith/spinfreebase.h>

namespace bagel {
namespace SMITH {

// DIIS utility classes
class Amplitude;

class Residual {
  friend class Amplitude;
  protected:
    std::vector<double> refcoeff_;
    std::vector<std::shared_ptr<Tensor>> res_;
    SpinFreeMethod* me_;

    void check_consistency(const Residual& o) const {
      assert(refcoeff_.size() == o.refcoeff_.size());
      assert(res_.size() == o.res_.size());
      assert(res_.size() == refcoeff_.size());
    }

    void check_consistency(const Amplitude& o) const {
      assert(res_.size() == refcoeff_.size());
      assert(refcoeff_.size() == o.refcoeff_.size());
      assert(res_.size() == o.amp_.size());
    }

  public:
    Residual(const double c, std::shared_ptr<const Tensor> r, SpinFreeMethod* m) : refcoeff_{c}, res_{r->copy()}, me_(m) { }
    Residual(const double c, std::shared_ptr<Tensor>&& r, SpinFreeMethod* m) : refcoeff_{c}, res_{std::move(r)}, me_(m) { }

    Residual(const std::vector<double>& c, const std::vector<std::shared_ptr<Tensor>>& r, SpinFreeMethod* m) : refcoeff_(c), me_(m) {
      for (auto& i : r) res_.push_back(i->copy());
    }

    void zero() {
      std::fill(refcoeff_.begin(), refcoeff_.end(), 0.0); 
      for (auto& i : res_)
        i->zero();
    }

    std::shared_ptr<Residual> clone() const {
      auto out = std::make_shared<Residual>(refcoeff_, res_, me_);
      out->zero();
      return out;
    }
    void synchronize() { }

    void ax_plus_y(const double a, const Residual& o) {
      check_consistency(o);
      auto j = o.refcoeff_.begin();
      for (auto i : refcoeff_)
        i += a * *j++;
      auto n = o.res_.begin();
      for (auto i : res_)
        i->ax_plus_y(a, *n++);
    }
    void ax_plus_y(const double a, std::shared_ptr<const Residual> o) { ax_plus_y(a, *o); }

    void ax_plus_y(const double a, const Amplitude& o);
    void ax_plus_y(const double a, std::shared_ptr<const Amplitude> o);

    double dot_product(const Amplitude& o) const;
    double dot_product(std::shared_ptr<const Amplitude> o) const;

    std::vector<std::shared_ptr<Tensor>> tensor() { return res_; }
};

class Amplitude {
  friend class Residual;
  protected:
    std::vector<double> refcoeff_;
    std::vector<std::shared_ptr<Tensor>> amp_;
    // the result of <ab/ij|\phi> so that one can easily compute the overlap
    std::vector<std::shared_ptr<Tensor>> left_;

    SpinFreeMethod* me_;

    void check_() const {
      assert(refcoeff_.size() == amp_.size());
      assert(refcoeff_.size() == left_.size());
    }

    void check_consistency(const Amplitude& o) const {
      check_();
      assert(refcoeff_.size() == o.refcoeff_.size()); 
      assert(amp_.size() == o.amp_.size()); 
      assert(left_.size() == o.left_.size()); 
    }

    void check_consistency(const Residual& o) const {
      check_();
      assert(refcoeff_.size() == o.refcoeff_.size()); 
      assert(amp_.size() == o.res_.size()); 
    }

  public:
    Amplitude(const double c, std::shared_ptr<const Tensor> t, std::shared_ptr<const Tensor> l, SpinFreeMethod* m)
     : refcoeff_{c}, amp_{t->copy()}, left_{l->copy()}, me_(m) { }

    Amplitude(const std::vector<double>& c, const std::vector<std::shared_ptr<Tensor>>& t, const std::vector<std::shared_ptr<Tensor>>& l, SpinFreeMethod* m)
     : refcoeff_(c), me_(m) {
      for (auto& i : t) amp_.push_back(i->copy());
      for (auto& i : l) left_.push_back(i->copy());
    }

    void zero() {
      std::fill(refcoeff_.begin(), refcoeff_.end(), 0.0);
      for (auto& i : amp_)
        i->zero();
      for (auto& i : left_)
        i->zero();
    }

    std::shared_ptr<Amplitude> clone() const {
      auto out = std::make_shared<Amplitude>(refcoeff_, amp_, left_, me_);
      out->zero();
      return out;
    }
    void synchronize() { }

    void ax_plus_y(const double a, const Amplitude& o) {
      check_consistency(o);
      auto j = o.refcoeff_.begin(); 
      for (auto& i : refcoeff_)
        i += a * *j++;
      auto n = o.amp_.begin();
      for (auto& i : amp_)
        i->ax_plus_y(a, *n++);
      n = o.left_.begin();
      for (auto& i : left_)
        i->ax_plus_y(a, *n++);
    }
    void ax_plus_y(const double a, std::shared_ptr<const Amplitude> o) { ax_plus_y(a, *o); }

    double dot_product(const Amplitude& o) const {
      check_consistency(o);
      double out = blas::dot_product(refcoeff_.data(), refcoeff_.size(), o.refcoeff_.data());
      auto j = o.left_.begin();
      for (auto& i : amp_)
        out += me_->dot_product_transpose(*j++, i);
      return out;  
    }
    double dot_product(std::shared_ptr<const Amplitude> o) const { return dot_product(*o); }

    double dot_product(const Residual& o) const {
      check_consistency(o);
      double out = blas::dot_product(refcoeff_.data(), refcoeff_.size(), o.refcoeff_.data());
      auto j = o.res_.begin();
      for (auto& i : amp_)
        out += me_->dot_product_transpose(*j++, i); 
      return out;
    }
    double dot_product(std::shared_ptr<const Residual> o) const { return dot_product(*o); }

    std::vector<std::shared_ptr<Tensor>> tensor() { return amp_; }
};

void Residual::ax_plus_y(const double a, const Amplitude& o) {
  check_consistency(o);
  auto j = o.refcoeff_.begin();
  for (auto& i : refcoeff_)
    i += *j++;
  auto n = o.left_.begin();
  for (auto& i : res_)
    i->ax_plus_y(a, *n++);
}
void Residual::ax_plus_y(const double a, std::shared_ptr<const Amplitude> o) { ax_plus_y(a, *o); }
double Residual::dot_product(const Amplitude& o) const {
  check_consistency(o);
  double out = blas::dot_product(refcoeff_.data(), refcoeff_.size(), o.refcoeff_.data());
  auto j = o.amp_.begin();
  for (auto& i : res_)
    out += me_->dot_product_transpose(i, *j++);
  return out;
}
double Residual::dot_product(std::shared_ptr<const Amplitude> o) const { return dot_product(*o); }

}
}

#endif
