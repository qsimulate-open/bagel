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
    double refcoeff_;
    std::shared_ptr<Tensor> amp_;
    SpinFreeMethod* me_;

  public:
    Residual(const double c, std::shared_ptr<const Tensor> r, SpinFreeMethod* m) : refcoeff_(c), amp_(r->copy()), me_(m) { }

    std::shared_ptr<Residual> clone() const { return std::make_shared<Residual>(refcoeff_, amp_, me_); }
    void synchronize() { }

    void ax_plus_y(const double a, const Residual& o) { refcoeff_ += a*o.refcoeff_; amp_->ax_plus_y(a, o.amp_); }
    void ax_plus_y(const double a, std::shared_ptr<const Residual> o) { ax_plus_y(a, *o); }

    void ax_plus_y(const double a, const Amplitude& o);
    void ax_plus_y(const double a, std::shared_ptr<const Amplitude> o);

    double dot_product(const Amplitude& o) const;
    double dot_product(std::shared_ptr<const Amplitude> o) const;

    std::shared_ptr<Tensor> tensor() { return amp_; }
};

class Amplitude {
  friend class Residual;
  protected:
    double refcoeff_;
    std::shared_ptr<Tensor> amp_;
    // the result of <ab/ij|\phi> so that one can easily compute the overlap
    std::shared_ptr<Tensor> left_;

    SpinFreeMethod* me_;

  public:
    Amplitude(const double c, std::shared_ptr<const Tensor> t, std::shared_ptr<const Tensor> l, SpinFreeMethod* m)
     : refcoeff_(c), amp_(t->copy()), left_(l->copy()), me_(m) { }

    std::shared_ptr<Amplitude> clone() const { return std::make_shared<Amplitude>(refcoeff_, amp_, left_, me_); }
    void synchronize() { }

    void ax_plus_y(const double a, const Amplitude& o) { refcoeff_ += a*o.refcoeff_; amp_->ax_plus_y(a, o.amp_); left_->ax_plus_y(a, o.left_); }
    void ax_plus_y(const double a, std::shared_ptr<const Amplitude> o) { ax_plus_y(a, *o); }

    double dot_product(const Amplitude& o) const { return refcoeff_*o.refcoeff_ + me_->dot_product_transpose(o.left_, amp_); }
    double dot_product(std::shared_ptr<const Amplitude> o) const { return dot_product(*o); }

    double dot_product(const Residual& o) const { return refcoeff_*o.refcoeff_ + me_->dot_product_transpose(o.amp_, amp_); }
    double dot_product(std::shared_ptr<const Residual> o) const { return dot_product(*o); }

    std::shared_ptr<Tensor> tensor() { return amp_; }
};

void Residual::ax_plus_y(const double a, const Amplitude& o) { refcoeff_ += a*o.refcoeff_; amp_->ax_plus_y(a, o.left_); }
void Residual::ax_plus_y(const double a, std::shared_ptr<const Amplitude> o) { ax_plus_y(a, *o); }
double Residual::dot_product(const Amplitude& o) const { return refcoeff_*o.refcoeff_ + me_->dot_product_transpose(amp_, o.amp_); }
double Residual::dot_product(std::shared_ptr<const Amplitude> o) const { return dot_product(*o); }

}
}

#endif
