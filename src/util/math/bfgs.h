//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: bfgs.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_UTIL_BFGS_H
#define __SRC_UTIL_BFGS_H

// implements BFGS based on Fischer & Almlof JPC 1992.
// T needs clone, ddot and daxpy, along with overloaded operators, a copy constructor
//         data(const size_t)

#include <stddef.h>
#include <vector>
#include <memory>

namespace bagel {

template<typename T>
class BFGS {
  protected:
    std::vector<std::shared_ptr<const T>> delta_;
    std::vector<std::shared_ptr<const T>> y_;
    std::vector<std::shared_ptr<const T>> D_;
    std::vector<std::shared_ptr<const T>> Y_;

    std::shared_ptr<const T> prev_grad;
    std::shared_ptr<const T> prev_value;

    // initial guess for hessian in a diagonal form
    const std::shared_ptr<const T> denom_;

    const bool debug_;

  public:
    BFGS(std::shared_ptr<const T> denom, bool debug = false) : denom_(denom), debug_(debug) {}

    std::shared_ptr<const T> denom() const { return denom_; }
    std::vector<std::shared_ptr<const T>> delta() { return delta_; }
    std::vector<std::shared_ptr<const T>> Y() const { return Y_; }
    std::vector<std::shared_ptr<const T>> D() const { return D_; }
    std::vector<std::shared_ptr<const T>> y() const { return y_; }

    // returns a displacement
    std::shared_ptr<T> extrapolate(std::shared_ptr<const T> _grad, std::shared_ptr<const T> _value) {
      // to make sure, inputs are copied.
      auto grad = std::make_shared<const T>(*_grad);
      auto value = std::make_shared<const T>(*_value);

      auto out = std::make_shared<T>(*grad);
      // (1)
      *out /= *denom_;

      if (prev_value != nullptr && !debug_) {
        // (3)
        std::shared_ptr<T> yy = grad->clone();
        {
          auto DD = std::make_shared<T>(*grad - *prev_grad);
          D_.push_back(DD);

          *yy = *DD / *denom_;

          auto vv = std::make_shared<T>(*value - *prev_value);
          delta_.push_back(vv);

        }
        const int n = delta_.size()-1;
        assert(delta_.size() == y_.size()+1 && y_.size()+1 == D_.size());

        // (4)
        for (int i = 0; i < n; ++i) {
          auto s1 = detail::real(1.0 / D_[i]->dot_product(delta_[i]));
          auto s2 = detail::real(1.0 / D_[i]->dot_product(y_[i]));
          auto s3 = delta_[i]->dot_product(grad);
          auto s4 = y_[i]->dot_product(grad);
          auto s5 = delta_[i]->dot_product(D_[n]);
          auto s6 = y_[i]->dot_product(D_[n]);
          auto t1 = (1.0 + s1/s2) * s1 * s3 - s1 * s4;
          auto t2 = s1 * s3;
          auto t3 = (1.0 + s1/s2) * s1 * s5 - s1 * s6;
          auto t4 = s1 * s5;
          out->ax_plus_y(t1, delta_[i]);
          out->ax_plus_y(-t2, y_[i]);
          yy->ax_plus_y(t3, delta_[i]);
          yy->ax_plus_y(-t4, y_[i]);
        }
        { // (5)
          auto s1 = detail::real(1.0 / D_[n]->dot_product(delta_[n]));
          auto s2 = detail::real(1.0 / D_[n]->dot_product(std::shared_ptr<const T>(yy)));
          auto s3 = delta_[n]->dot_product(grad);
          auto s4 = yy->dot_product(grad);
          auto t1 = (1.0 + s1/s2) * s1 * s3 - s1 * s4;
          auto t2 = s1 * s3;
          out->ax_plus_y(t1, delta_[n]);
          out->ax_plus_y(-t2, std::shared_ptr<const T>(yy));
        }
        y_.push_back(yy);
        assert(y_.size() == n+1);
      }
      prev_grad = grad;
      prev_value = value;
      return out;
    }

    // returns Hessian * value ; should be called after extrapolate
    std::shared_ptr<T> interpolate_hessian(std::shared_ptr<const T> _value, const bool update = false) {
      // to make sure, inputs are copied.
      auto value = std::make_shared<const T>(*_value);

      auto out = std::make_shared<T>(*value);
      // (1)
      *out *= *denom_;

      if (!delta_.empty() && !debug_) {
        // (3)
        std::shared_ptr<T> yy = value->clone(); // used to accumulate Y
        {
          // D and delta_ available from extrapolation
          auto DD = delta_.back();

          *yy = *DD * *denom_;

        }
        const int n = delta_.size()-1;
        assert(delta_.size() == Y_.size()+1 && Y_.size()+1 == D_.size());

        // (4)
        for (int i = 0; i < n; ++i) {
          auto s1 = detail::real(1.0 / D_[i]->dot_product(delta_[i]));
          auto s2 = detail::real(1.0 / Y_[i]->dot_product(delta_[i]));
          auto s3 = D_[i]->dot_product(delta_[n]);
          auto s4 = Y_[i]->dot_product(delta_[n]);
          auto s5 = D_[i]->dot_product(value);
          auto s6 = Y_[i]->dot_product(value);
          auto t1 = s1 * s5; // updating w
          auto t2 = s2 * s6; // updating w
          auto t3 = s1 * s3; // updating Y
          auto t4 = s2 * s4; // updating Y
          out->ax_plus_y(t1, D_[i]);
          out->ax_plus_y(-t2, Y_[i]);
          yy->ax_plus_y(t3, D_[i]);
          yy->ax_plus_y(-t4, Y_[i]);
        }
        { // (5)
          auto s1 = detail::real(1.0 / D_[n]->dot_product(delta_[n]));
          auto s2 = detail::real(1.0 / delta_[n]->dot_product(std::shared_ptr<const T>(yy)));
          auto s3 = D_[n]->dot_product(value);
          auto s4 =   yy->dot_product(value);
          auto t1 = s1 * s3;
          auto t2 = s2 * s4;
          out->ax_plus_y(t1, D_[n]);
          out->ax_plus_y(-t2, std::shared_ptr<const T>(yy));
        }
        if (update) {
          Y_.push_back(yy);
          assert(Y_.size() == n+1);
       }
      }
      return out;
    }

};

}

#endif
