//
// BAGEL - Parallel electron correlation program.
// Filename: bfgs.h
// Copyright (C) 2012 Toru Shiozaki
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
    std::vector<std::shared_ptr<const T>> delta;
    std::vector<std::shared_ptr<const T>> y;
    std::vector<std::shared_ptr<const T>> D;
    std::vector<std::shared_ptr<const T>> Y;

    std::shared_ptr<const T> prev_grad;
    std::shared_ptr<const T> prev_value;

    // initial guess for hessian in a diagonal form
    const std::shared_ptr<const T> denom_;

    const bool debug_;

  public:
    BFGS(std::shared_ptr<const T> denom, bool debug = false) : denom_(denom), debug_(debug) {}
    std::shared_ptr<const T> denom() const { return denom_; }

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
          D.push_back(DD);

          *yy = *DD / *denom_;

          auto vv = std::make_shared<T>(*value - *prev_value);
          delta.push_back(vv);

        }
        const int n = delta.size()-1;
        assert(delta.size() == y.size()+1 && y.size()+1 == D.size());

        // (4)
        for (int i = 0; i < n; ++i) {
          auto s1 = detail::real(1.0 / D[i]->dot_product(delta[i]));
          auto s2 = detail::real(1.0 / D[i]->dot_product(y[i]));
          auto s3 = delta[i]->dot_product(grad);
          auto s4 = y[i]->dot_product(grad);
          auto s5 = delta[i]->dot_product(D[n]);
          auto s6 = y[i]->dot_product(D[n]);
          auto t1 = (1.0 + s1/s2) * s1 * s3 - s1 * s4;
          auto t2 = s1 * s3;
          auto t3 = (1.0 + s1/s2) * s1 * s5 - s1 * s6;
          auto t4 = s1 * s5;
          out->ax_plus_y(t1, delta[i]);
          out->ax_plus_y(-t2, y[i]);
          yy->ax_plus_y(t3, delta[i]);
          yy->ax_plus_y(-t4, y[i]);
        }
        { // (5)
          auto s1 = detail::real(1.0 / D[n]->dot_product(delta[n]));
          auto s2 = detail::real(1.0 / D[n]->dot_product(std::shared_ptr<const T>(yy)));
          auto s3 = delta[n]->dot_product(grad);
          auto s4 =       yy->dot_product(grad);
          auto t1 = (1.0 + s1/s2) * s1 * s3 - s1 * s4;
          auto t2 = s1 * s3;
          out->ax_plus_y(t1, delta[n]);
          out->ax_plus_y(-t2, std::shared_ptr<const T>(yy));
        }
        y.push_back(yy);
        assert(y.size() == n+1);
      }
      prev_grad = grad;
      prev_value = value;
      return out;
    }

    // returns Hessian * kappa as part of the trust radius procedure ; should be called after extrapolate
    std::shared_ptr<T> interpolate(std::shared_ptr<const T> _value) {
      // to make sure, inputs are copied.
      auto value = std::make_shared<const T>(*_value);

      auto out = std::make_shared<T>(*value);
      // (1)
      *out *= *denom_;

      if (!delta.empty() && !debug_) {
        // (3)
        std::shared_ptr<T> yy = value->clone(); // used to accumulate Y
        {
          // D and delta available from extrapolation
          auto DD = delta.back();

          *yy = *DD * *denom_;

        }
        const int n = delta.size()-1;
        assert(delta.size() == Y.size()+1 && Y.size()+1 == D.size());

        // (4)
        for (int i = 0; i < n; ++i) {
          auto s1 = detail::real(1.0 / D[i]->dot_product(delta[i]));
          auto s2 = detail::real(1.0 / Y[i]->dot_product(delta[i]));
          auto s3 = D[i]->dot_product(delta[n]);
          auto s4 = Y[i]->dot_product(delta[n]);
          auto s5 = D[i]->dot_product(value);
          auto s6 = Y[i]->dot_product(value);
          auto t1 = s1 * s5; // updating w
          auto t2 = s2 * s6; // updating w
          auto t3 = s1 * s3; // updating Y
          auto t4 = s2 * s4; // updating Y
          out->ax_plus_y(t1, D[i]);
          out->ax_plus_y(-t2, Y[i]);
          yy->ax_plus_y(t3, D[i]);
          yy->ax_plus_y(-t4, Y[i]);
        }
        { // (5)
          auto s1 = detail::real(1.0 / D[n]->dot_product(delta[n]));
          auto s2 = detail::real(1.0 / delta[n]->dot_product(std::shared_ptr<const T>(yy)));
          auto s3 = D[n]->dot_product(value);
          auto s4 =   yy->dot_product(value);
          auto t1 = s1 * s3;
          auto t2 = s2 * s4;
          out->ax_plus_y(t1, D[n]);
          out->ax_plus_y(-t2, std::shared_ptr<const T>(yy));
        }
        Y.push_back(yy);
        assert(Y.size() == n+1);
      }
      return out;
    }
};

}

#endif
