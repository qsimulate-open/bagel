//
// BAGEL - Parallel electron correlation program.
// Filename: bfgs.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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

    std::shared_ptr<const T> prev_grad;
    std::shared_ptr<const T> prev_value;

    // initial guess for hessian in a diagonal form
    const std::shared_ptr<const T> denom_;

    const bool debug_;

  public:
    BFGS(std::shared_ptr<const T> denom, bool debug = false) : denom_(denom), debug_(debug) {}

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
          const double s1 = 1.0 / delta[i]->ddot(D[i]);
          const double s2 = 1.0 / D[i]->ddot(y[i]);
          const double s3 = delta[i]->ddot(grad);
          const double s4 =     y[i]->ddot(grad);
          const double s5 = delta[i]->ddot(D[n]);
          const double s6 =     y[i]->ddot(D[n]);
          const double t1 = (1.0 + s1/s2) * s1 * s3 - s1 * s4;
          const double t2 = s1 * s3;
          const double t3 = (1.0 + s1/s2) * s1 * s5 - s1 * s6;
          const double t4 = s1 * s5;
          out->daxpy(t1, delta[i]);
          out->daxpy(-t2, y[i]);
          yy->daxpy(t3, delta[i]);
          yy->daxpy(-t4, y[i]);
        }
        { // (5)
          const double s1 = 1.0 / delta[n]->ddot(D[n]);
          const double s2 = 1.0 /     D[n]->ddot(std::shared_ptr<const T>(yy));
          const double s3 = delta[n]->ddot(grad);
          const double s4 =       yy->ddot(grad);
          const double t1 = (1.0 + s1/s2) * s1 * s3 - s1 * s4;
          const double t2 = s1 * s3;
          out->daxpy(t1, delta[n]);
          out->daxpy(-t2, std::shared_ptr<const T>(yy));
        }
        y.push_back(yy);
        assert(y.size() == n+1);
      }
      prev_grad = grad;
      prev_value = value;
      return out;
    }

};

}

#endif
