//
// BAGEL - Parallel electron correlation program.
// Filename: step_restrict_bfgs.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Jefferson E Bates <jefferson.bates@northwestern.edu>
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


#ifndef __SRC_UTIL_SRBFGS_H
#define __SRC_UTIL_SRBFGS_H

// implements BFGS updates based on Fischer & Almlof JPC 1992.
// step restriction loosely based on Jorgensen, Swanstrom, & Yeager (JSY) JCP 78 347 (1983).
// TODO : implement level shift and trust radius procedures
// T needs clone, ddot and daxpy, along with overloaded operators, a copy constructor
//         data(const size_t)

#include <stddef.h>
#include <vector>
#include <memory>

namespace bagel {

template<typename T>
class SRBFGS {
  protected:
    std::vector<std::shared_ptr<const T>> delta_;
    std::vector<std::shared_ptr<const T>> y_;
    std::vector<std::shared_ptr<const T>> D_;
    std::vector<std::shared_ptr<const T>> Y_;

    std::shared_ptr<const T> prev_grad_;
    std::shared_ptr<const T> prev_value_;

    // initial guess for hessian in a diagonal form
    std::shared_ptr<const T> denom_;

    double trust_radius_;
    double rk_;
    double level_shift_;
    const bool debug_;
    const int hebden_iter_ = 10;
    const int maxiter_ = 300;
    // default convergence parameters taken from JSY
    const double rmin_ = 0.6;
    const double rgood_ = 0.85;
    const double alpha_ = 1.3;

  public:
    SRBFGS(std::shared_ptr<const T> denom, double level_shift = 0.0, bool debug = false) : denom_(denom), level_shift_(level_shift), debug_(debug) {}

    std::shared_ptr<const T> denom() const { return denom_; }
    std::vector<std::shared_ptr<const T>> delta() { return delta_; }
    std::vector<std::shared_ptr<const T>> Y() const { return Y_; }
    std::vector<std::shared_ptr<const T>> D() const { return D_; }
    std::vector<std::shared_ptr<const T>> y() const { return y_; }
    std::shared_ptr<const T> prev_grad() const { return prev_grad_; }
    std::shared_ptr<const T> prev_value() const { return prev_value_; }
    std::shared_ptr<const T> level_shift() const { return level_shift_; }

    // sets initial trust radius ; presently not used
    void initiate_trust_radius(std::shared_ptr<const T> _grad) {
      // to make sure, inputs are copied.
      auto grad = std::make_shared<const T>(*_grad);
      auto tmp1 = std::make_shared<const T>(*grad / *denom_);
      trust_radius_ = std::sqrt(detail::real(tmp1->dot_product(tmp1)));
    }


    // returns Hessian * value as part of the restricted step procedure ; should be called after extrapolate
    std::shared_ptr<T> interpolate_hessian(std::shared_ptr<const T> _value, std::shared_ptr<const T> _shift, const bool update = false) {
      // to make sure, inputs are copied.
      auto value = std::make_shared<const T>(*_value);
      auto shift = std::make_shared<const T>(*_shift);

      auto out = std::make_shared<T>(*value);
      // (1)
      *out *= (*denom_ + *shift);

      if (!delta_.empty() && !debug_) {
        // (3)
        std::shared_ptr<T> yy = value->clone(); // used to accumulate Y
        {
          // D and delta_ available from extrapolation
          auto DD = delta_.back();

          *yy = *DD * (*denom_ + *shift);

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


    // return search direction for given level shift
    std::shared_ptr<T> extrapolate_micro(std::shared_ptr<const T> _grad, std::shared_ptr<const T> _value, std::shared_ptr<const T> _shift, const bool update = false) {
      // to make sure, inputs are copied.
      auto grad = std::make_shared<const T>(*_grad);
      auto value = std::make_shared<const T>(*_value);
      auto shift = std::make_shared<const T>(*_shift);

      auto out = std::make_shared<T>(*grad);
      // (1)
      *out /= (*denom_ + *shift);

      if (prev_value_ != nullptr && !debug_) {
        // (3)
        std::shared_ptr<T> yy = grad->clone();
        {
          auto DD = std::make_shared<T>(*grad - *prev_grad_);
          D_.push_back(DD);

          *yy = *DD / (*denom_ + *shift);

          auto vv = std::make_shared<T>(*value - *prev_value_);
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
      if (update) {
        prev_grad_ = grad;
        prev_value_ = value;
      }
      return out;
    }


    // apply inverse hessian according to BFGS recursion without updating intermediate quantities
    std::shared_ptr<T> apply_inverse_hessian(std::shared_ptr<const T> _grad, std::shared_ptr<const T> _value, std::shared_ptr<const T> _vector, std::shared_ptr<const T> _shift) {
      // applies hessian inverse without updating intermediates
      // to make sure, inputs are copied.
      auto grad = std::make_shared<const T>(*_grad);
      auto value = std::make_shared<const T>(*_value);
      auto shift = std::make_shared<const T>(*_shift);
      auto vector = std::make_shared<const T>(*_vector);
      auto out = std::make_shared<T>(*vector);
      // (1)
      *out /= (*denom_ + *shift);
      if (prev_value_ != nullptr && !debug_) {
        // (3)
        std::shared_ptr<T> yy = grad->clone();

        auto DD = std::make_shared<T>(*grad - *prev_grad_);

        *yy = *DD / (*denom_ + *shift);

        auto vv = std::make_shared<T>(*value - *prev_value_);

        const int n = delta_.size()-1;
//        assert(delta_.size() == y_.size()+1 && y_.size()+1 == D_.size());
  
        // (4)
        for (int i = 0; i < n; ++i) {
          auto s1 = detail::real(1.0 / D_[i]->dot_product(delta_[i]));
          auto s2 = detail::real(1.0 / D_[i]->dot_product(y_[i]));
          auto s3 = delta_[i]->dot_product(vector);
          auto s4 = y_[i]->dot_product(vector);
          auto s5 = delta_[i]->dot_product(DD);
          auto s6 = y_[i]->dot_product(DD);
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
          auto s1 = detail::real(1.0 / DD->dot_product(vv));
          auto s2 = detail::real(1.0 / DD->dot_product(std::shared_ptr<const T>(yy)));
          auto s3 = vv->dot_product(vector);
          auto s4 = yy->dot_product(vector);
          auto t1 = (1.0 + s1/s2) * s1 * s3 - s1 * s4;
          auto t2 = s1 * s3;
          out->ax_plus_y(t1, vv);
          out->ax_plus_y(-t2, std::shared_ptr<const T>(yy));
        }
      }
      return out;
    }


    double taylor_series_validity_ratio(const std::vector<double> _f, std::shared_ptr<const T> _grad, std::shared_ptr<const T> _a, std::shared_ptr<const T> _Ha) const {
    // Following JSY and Jensen and Jorgensen (JCP 80 1204 1984)
    // Returns r_k = ( f(a_k) - f(a_(k-1)) )/( f^(2)(a_k) - f(a_(k-1)) ) (Eq 64) 
    //   where f^(2) is the second order Taylor expansion, f^(2)(a) = f(0) + g*a + 1/2 * a^T H a
    // TODO : remove dependency on iter, should be able to use back, etc to refernce the proper elements

      // to make sure, inputs are copied.
      auto f    = std::vector<double>(_f);
      auto a    = std::make_shared<const T>(*_a);
      auto Ha   = std::make_shared<const T>(*_Ha);
      auto grad = std::make_shared<const T>(*_grad);
      assert(!f.empty());
    
      auto f1     = ( f.size() > 1 ? *(f.end()-2) : 0.0 );
      auto DeltaE = f.back() - f1;
    
      auto e2 = 0.5 * (detail::real(a->dot_product(Ha)));
      auto e1 = detail::real(a->dot_product(grad));
      auto e0 = f.back() + e1 + e2;
      e0 -= f1;
      DeltaE /= e0;
      return DeltaE;
    }


    // returns restricted step displacement
    std::shared_ptr<T> step_restricted_extrapolate(const std::vector<double> _f, std::shared_ptr<const T> _grad, std::shared_ptr<const T> _value, std::shared_ptr<const T> _shift, const bool tight = false) {
      // to make sure, inputs are copied
      auto f     = std::vector<double>(_f);
      auto grad  = std::make_shared<const T>(*_grad);
      auto value = std::make_shared<const T>(*_value);
      auto shift = std::make_shared<const T>(*_shift);

      double alpha = alpha_;
      if (tight) {
        alpha = alpha + static_cast<double>(f.size()-1) * 0.25;
      }

      std::shared_ptr<T> acopy = extrapolate_micro(grad, value, shift, true);
#if 0
      std::shared_ptr<T> acopy = extrapolate_micro(grad, value, shift, false);
      for (int mi = 0; mi!= maxiter_; ++mi) {
        std::shared_ptr<T> V = interpolate_hessian(acopy, shift, false); // H_n * a
        double rk = taylor_series_validity_ratio(f, grad, acopy, V);
        std::cout << std::setprecision(4) << " Taylor expansion validity parameter  = " << rk << std::endl;
        if (!tight) {
          if (rk < 0.25 && rk > 0) {
            std::cout << " condition (i) satisfied in microiteration " << mi << std::endl;
            std::cout << " end microiteration " << mi << std::endl;
            if (delta().size() != 0)
              decrement_intermediates();
            auto tmp1 = extrapolate_micro(grad, value, shift, true);
            auto v    = interpolate_hessian(value, shift, true);
            trust_radius_ = trust_radius_ * 2.0/3.0;
            break;
          } else if ((rk > 0.25) && (rk < 0.75)) {
            std::cout << " condition (ii) satisfied in microiteration " << mi << std::endl;
            std::cout << " end microiteration " << mi << std::endl;
            trust_radius_ = trust_radius_;
            if (delta().size() != 0)
              decrement_intermediates();
            auto tmp1 = extrapolate_micro(grad, value, shift, true);
            auto v    = interpolate_hessian(value, shift, true);
            break;
          } else if (rk > 0.75) {
            std::cout << " condition (iii) satisfied in microiteration " << mi << std::endl;
            std::cout << " end microiteration " << mi << std::endl;
            trust_radius_ = std::min(1.2 * trust_radius_, 0.75);
            if (delta().size() != 0)
              decrement_intermediates();
            auto tmp1 = extrapolate_micro(grad, value, shift, true);
            auto v    = interpolate_hessian(value, shift, true);
            break;
          } else if (rk < 0) {
            std::cout << " step does not satisfy Taylor expansion criteria " << std::endl;
            std::cout << " scaling down the step vector " << std::endl;
            acopy->scale(1.0/alpha);
            trust_radius_ = trust_radius_ * 2.0/3.0;
            std::cout << " end microiteration " << mi << std::endl;
          }
        } else { // conditions for tighter optimization
          if (((rmin_ < rk) && (rk < rgood_)) || ((2.0 - rgood_ < rk) && (rk < 2.0 - rmin_))) {
            std::cout << " condition (i) satisfied in microiteration " << mi << std::endl;
            std::cout << " end microiteration " << mi << std::endl;
            if (delta().size() != 0)
              decrement_intermediates();
            auto tmp1 = extrapolate_micro(grad, value, shift, true);
            auto v    = interpolate_hessian(value, shift, true);
            trust_radius_ = trust_radius_;
            break;
          } else if ((rgood_ < rk) && (rk < 2 - rgood_) ) {
            std::cout << " condition (ii) satisfied in microiteration " << mi << std::endl;
            trust_radius_ = std::min(alpha * trust_radius_, 0.75);
            std::cout << " new trust radius   = " << trust_radius_ << std::endl;
            if (delta().size() != 0)
              decrement_intermediates();
            auto tmp1 = extrapolate_micro(grad, value, shift, true);
            auto v    = interpolate_hessian(value, shift, true);
            std::cout << " end microiteration " << mi << std::endl;
            break;
          } else if (rk < rmin_ || rk > 2 - rmin_ ) {
            std::cout << " step does not satisfy Taylor expansion criteria " << std::endl;
            std::cout << " scaling down the step vector " << std::endl;
            if (delta().size() == Y().size() && delta().size() != 0)
              decrement_intermediates();
            acopy->scale(1.0/alpha);
            trust_radius_ = trust_radius_ / alpha;
            std::cout << " pk * grad    = " << acopy->dot_product(grad) << std::endl;
            std::cout << " end microiteration " << mi << std::endl;
          }
        }
      }
#endif
      return acopy;
    }


   double newton_levelshift(std::shared_ptr<const T> _grad, std::shared_ptr<const T> _value) {
     // iteratively finds an appropriate level shift according to Newton-Raphson algorithm
     // that satisfies     fn(v)  = h_k^2 ; hk is the trust radius and fn(v) = gn^+ (Hn + vI)^-2 gn
     // exact gradient and hessian 
     // TODO: need a way to guess intial value from denom
     
     // to make sure, inputs are copied
     auto grad  = std::make_shared<const T>(*_grad);
     auto value = std::make_shared<const T>(*_value);
     if (delta().size() == 0) 
       initiate_trust_radius(grad);
     std::cout << " trust radius    = " << trust_radius_ << std::endl;
     
     double shift = fabs(level_shift_);
     std::cout << "shift = " << shift << std::endl;
     auto shift_vec = value->clone();
     shift_vec->fill(level_shift_);
     for (int k = 0; k != 20; ++k) {
       auto dl  = apply_inverse_hessian(grad, value, grad, shift_vec); // Hn^-1 * gn
       std::cout << std::setprecision(6) << " dl norm     = " << 
            std::sqrt(detail::real(dl->dot_product(dl))) << std::endl;
       auto dl2 = apply_inverse_hessian(grad, value, dl, shift_vec);   // Hn^-2 * gn
       auto dl3 = detail::real(dl2->dot_product(dl));                  // gn * Hn^-3 * gn
       auto dl4 = detail::real(dl2->dot_product(dl2));                 // gn * Hn^-4 * gn
       auto dshift = dl3 / dl4 / 3.0;                                  // factor comes from derivatives
       std::cout << " dshift  = " << dshift << std::endl;
       shift += dshift;
       shift_vec->fill(shift);
       dl       = apply_inverse_hessian(grad, value, grad, shift_vec); 
       auto dl_norm  = std::sqrt(detail::real(dl->dot_product(dl)));
       std::cout << " step norm with new shift    = " << 
          std::sqrt(detail::real(dl->dot_product(dl))) << std::endl;
       if (dl_norm <= trust_radius_) break;
       std::cout << "shift end loop = " << shift << std::endl;
     }
     level_shift_ = shift;
     return shift;
   }
    // decrement intermediates
    void decrement_y() { y_.pop_back(); }
    void decrement_delta() { delta_.pop_back(); }
    void decrement_D() { D_.pop_back(); }
    void decrement_Y() { Y_.pop_back(); }

    void decrement_intermediates() {
      decrement_y();
      decrement_delta();
      decrement_D();
    }

};

}

#endif
