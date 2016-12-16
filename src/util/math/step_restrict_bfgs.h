//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: step_restrict_bfgs.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Jefferson E Bates <jefferson.bates@northwestern.edu>
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


#ifndef __SRC_UTIL_SRBFGS_H
#define __SRC_UTIL_SRBFGS_H

// implements BFGS updates based on Fischer & Almlof JPC 1992.
// step restriction based on Jorgensen, Swanstrom, & Yeager (JSY) JCP 78 347 (1983).
// level shifted BFGS based on Erway, Marcia (EM) Linear Algebra and its Applications 473 333 (2012)
//      and Erway, Jain, Marcia, arXiv:1209.5141.
// Hebden algorithm : Hebden, AERE Harwell Report TP515 (1973).
// T needs clone, ddot and daxpy, along with overloaded operators, a copy constructor
//         data(const size_t), norm, fill

#include <stddef.h>
#include <vector>
#include <memory>

namespace bagel {

template<typename T>
class SRBFGS {
  protected:
    std::vector<std::shared_ptr<const T>> delta_;
    std::vector<std::shared_ptr<const T>> avec_;
    std::vector<std::shared_ptr<const T>> y_;
    std::vector<std::shared_ptr<const T>> D_;
    std::vector<std::shared_ptr<const T>> Y_;
    std::vector<double> rho_;

    std::shared_ptr<const T> prev_grad_;
    std::shared_ptr<const T> prev_value_;

    // initial guess for hessian in a diagonal form
    std::shared_ptr<const T> denom_;

    double trust_radius_ = 0.0;
    double rk_;
    double level_shift_;
    double prev_level_shift_;
    const int hebden_iter_ = 75;
    // default convergence parameters taken from JSY
    const double alpha_ = 1.3;
    const double rmin_ = 0.6;
    const double rgood_ = 0.85;
    const bool debug_;

  public:
    SRBFGS(std::shared_ptr<const T> _denom, const double rad = 0.4, const int _hebden_iter = 75, const double _alpha = 1.3,
           const bool _debug = false) : denom_(_denom), trust_radius_(rad), hebden_iter_(_hebden_iter), alpha_(_alpha), debug_(_debug) {}

    std::shared_ptr<const T> denom() const { return denom_; }
    std::vector<std::shared_ptr<const T>> delta() { return delta_; }
    std::vector<std::shared_ptr<const T>> Y() const { return Y_; }
    std::vector<std::shared_ptr<const T>> D() const { return D_; }
    std::vector<std::shared_ptr<const T>> y() const { return y_; }
    std::vector<std::shared_ptr<const T>> avec() const { return avec_; }
    std::vector<double> rho() const { return rho_; }
    std::shared_ptr<const T> prev_grad() const { return prev_grad_; }
    std::shared_ptr<const T> prev_value() const { return prev_value_; }
    double level_shift() const { return level_shift_; }
    double prev_level_shift() const { return prev_level_shift_; }

    // returns Hessian * value as part of the restricted step procedure ; should be called after extrapolate
    std::shared_ptr<T> interpolate_hessian(std::shared_ptr<const T> _value, std::shared_ptr<const T> _shift, const bool update = false) {
      // to make sure, inputs are copied.
      auto value = std::make_shared<const T>(*_value);
      auto shift = std::make_shared<const T>(*_shift);

      auto out = std::make_shared<T>(*value);
      // (1)
      *out *= (*denom_ + *shift);

      if (!delta().empty() && !debug_) {
        // (3)
        std::shared_ptr<T> yy = value->clone(); // used to accumulate Y
        {
          // D and delta_ available from extrapolation
          auto DD = delta().back();

          *yy = *DD * (*denom_ + *shift);

        }
        const int n = delta().size()-1;
        assert(delta().size() == Y().size()+1 && Y().size()+1 == D().size()); // NOT SURE ABOUT THIS

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


    // apply inverse hessian according to Fischer and Almlof BFGS recursion ; y_ is updated
    std::shared_ptr<T> apply_inverse_hessian(std::shared_ptr<const T> _vector) {
      // to make sure, inputs are copied.
      auto vector = std::make_shared<const T>(*_vector);
      auto out = std::make_shared<T>(*vector);
      // (1)
      *out /= *denom_;
      if (prev_value_ != nullptr && !debug_) {
        // (3)
        std::shared_ptr<T> yy = vector->clone();

        auto DD = D().back();

        *yy = *DD / *denom_;

        auto vv = delta().back();

        const int n = delta().size()-1;
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
        y_.push_back(yy);
        assert(y_.size() == n+1);
      }
      return out;
    }


    double taylor_series_validity_ratio(const std::vector<double> _f, std::shared_ptr<const T> _grad, std::shared_ptr<const T> _a) const {
    // Following JSY and Jensen and Jorgensen (JCP 80 1204 1984)
    // Returns r_k = ( f(a_k) - f(a_(k-1)) )/( f^(2)(a_k) - f(a_(k-1)) ) (Eq 64)
    //   where f^(2) is the second order Taylor expansion, f^(2)(a_k) = f(a_(k-1)) + g_k*a_k + 1/2 * a_k^T H_k a_k

      // to make sure, inputs are copied.
      auto f    = std::vector<double>(_f);
      auto a    = std::make_shared<const T>(*_a);
      auto grad = std::make_shared<const T>(*_grad);
      assert(!f.empty());

      auto f1     = ( f.size() > 1 ? *(f.end()-2) : 0.0 );
      auto DeltaE = f.back() - f1;

      auto dg = -1.0 * detail::real(a->dot_product(grad));
      auto dd = detail::real(a->dot_product(a));

      auto e1 = detail::real(a->dot_product(grad));
      auto e2 = 0.5 * ( dg - prev_level_shift() * dd );
      auto e0 = e1 + e2;
      DeltaE /= e0;
      return DeltaE;
    }


   // iteratively finds an appropriate level shift according to hebden's algorithm (JSY)
   double hebden_levelshift(std::shared_ptr<const T> _grad) {
     // No shift is preferred when steps are within the trust radius

     // to make sure, inputs are copied
     auto grad  = std::make_shared<const T>(*_grad);
     bool converged = false;

     double shift = 1e-12;
     auto shift_vec = grad->clone();
     shift_vec->fill(shift);
     double dl_norm = 0.0;
     for (int k = 0; k != hebden_iter_; ++k) {
       auto dl  = level_shift_inverse_hessian(grad, shift_vec); // Hn^-1 * gn
       dl_norm = dl->norm();//std::sqrt(detail::real(dl->dot_product(dl)));
       if (dl_norm <= trust_radius_ && k != 0) {
         std::cout << " Hebden algorithm converged in " << k << " iterations. " << std::endl;
         std::cout << " Level Shift = " << shift << std::endl;
         converged = true;
         break;
       }
       auto dl2 = level_shift_inverse_hessian(dl, shift_vec);   // Hn^-2 * gn
       auto dl2_norm = detail::real(dl->dot_product(dl2)) / dl_norm;
       auto t1  =  dl_norm / dl2_norm;
       auto t2  =  dl_norm / trust_radius_;
       auto dshift = (t2 - 1.0) * t1;
       shift += dshift;
       shift_vec->fill(shift);
       if (k == hebden_iter_ - 1) {
         std::cout << " Hebden algorithm did not converge to appropriate level shift within " << k << " iterations " << std::endl;
         std::cout << " step norm with shift   = " << dl_norm << std::endl;
         converged = false;
       }
     }

     if (!converged) {
       // If the step size has exploded, then we throw away the level shift
       auto temp = grad->clone();
       temp->fill(1e-12);
       auto dl2 = level_shift_inverse_hessian(grad, temp); // Hn^-1 * gn
       const double unshifted_norm = dl2->norm();
       if (unshifted_norm < dl_norm) {
         std::cout << " Level shift will be discarded." << std::endl;
         shift = 1e-12;
       }
     }
     level_shift_ = shift;
     return shift;
   }


   // returns a level-shifted displacement according to EM12
    std::shared_ptr<T> level_shift_inverse_hessian(std::shared_ptr<const T> _vector, std::shared_ptr<const T> _shift) {
      // to make sure, inputs are copied.
      auto vector = std::make_shared<const T>(*_vector);
      auto shift = std::make_shared<const T>(*_shift);
      auto out = std::make_shared<T>(*vector);

      *out /= (*denom_ + *shift);

      if (prev_value() != nullptr) {
        const int n = delta().size();
        assert(delta().size() == avec().size() && delta().size() == D().size());

        std::vector<std::shared_ptr<T>> pvec;
        std::vector<std::shared_ptr<T>> vvec;
        std::vector<double> vcoeff;
        for (int j = 0; j != 2*n; ++j) {
          auto ptmp  = vector->clone();
          auto utmp  = vector->clone();
          auto vtmp  = vector->clone();
          auto vtild = vector->clone();
          if (j%2 == 0) {
            vtmp = std::make_shared<T>(*avec().at(j/2));
            utmp = vtmp->copy();
            auto s1 = 1.0 / detail::real(utmp->dot_product(delta().at(j/2)));
            assert(fabs(s1) > 1e-20);
            utmp->scale(s1);
          } else {
            vtmp = std::make_shared<T>(*D().at((j-1)/2));
            utmp = vtmp->copy();
            auto s1 = rho().at((j-1)/2);
            assert(fabs(s1) > 1e-20);
            utmp->scale(s1);
          }
          *ptmp  = *utmp / (*denom_ + *shift);
          *vtild = *vtmp / (*denom_ + *shift);
          for (int i = 0; i != j; ++i) {
            auto s1 = detail::real(vvec.at(i)->dot_product(utmp));
            auto s2 = pow(-1.0, i%2) * vcoeff.at(i) * s1;
            ptmp->ax_plus_y(s2, pvec.at(i));
                 s1 = detail::real(vvec.at(i)->dot_product(vtmp));
                 s2 = pow(-1.0, i%2) * vcoeff.at(i) * s1;
            vtild->ax_plus_y(s2, pvec.at(i));
          }
          auto s1   = detail::real(vtmp->dot_product(ptmp));
          auto s2   = 1.0 + pow(-1.0, j%2+1) * s1;
          assert(fabs(s2) > 1e-20);
          auto tau  = 1.0 / s2;
          s1        = detail::real(vtild->dot_product(vector));
          s2        = pow(-1.0, j%2) * tau * s1;
          out->ax_plus_y(s2, ptmp);
          vcoeff.push_back(tau);
          pvec.push_back(ptmp);
          vvec.push_back(vtild);
        }
      }
      return out;
   }


   // decrement intermediates
   void decrement_y() { y_.pop_back(); }
   void decrement_delta() { delta_.pop_back(); }
   void decrement_D() { D_.pop_back(); }
   void decrement_Y() { Y_.pop_back(); }
   void decrement_avec() { avec_.pop_back(); }
   void decrement_rho() { rho_.pop_back(); }

   void decrement_intermediates() {
     decrement_avec();
     decrement_delta();
     decrement_D();
     decrement_rho();
   }

    std::shared_ptr<T> unrolled_hessian(std::shared_ptr<const T> _value, const bool update = true) {
      // TODO : reorganize loop structure ; first double loop can probably be replaced
      // to make sure, inputs are copied.
      auto value = std::make_shared<const T>(*_value);
      auto out = std::make_shared<T>(*_value);

      // apply H0 * v
      *out *= *denom_;

      // get size of intermediate
      const int nd = delta().empty() ? 0 : delta().size()-1;
      for (int i = 0; i != nd; ++i) {
        auto yy =  std::make_shared<T>(*delta().at(nd-1));
        *yy *= *denom_;
        for (int j = 0; j != i; ++j) {
          auto t1 = D().at(j)->copy();
          auto t2 = delta().at(j)->copy();
          auto t3 = Y().at(j)->copy();
          auto t4 = delta().at(i)->copy();
          auto s1 = rho().at(j);                       // Delta_j * delta_j
          auto s2 = detail::real(t3->dot_product(t2)); // y_j     * delta_j
          auto s3 = t1->dot_product(t4);               // Delta_j * delta_i
          auto s4 = t3->dot_product(t4);               // y_j     * delta_i
          yy->ax_plus_y( s3 * s1, t1);
          yy->ax_plus_y(-s4 / s2, t3);
        }
        if (i == nd-1)
          Y_.push_back(yy);
      }
      for (int i = 0; i != nd; ++i) {
        auto t1 = D().at(i)->copy();
        auto t2 = delta().at(i)->copy();
        auto t3 = Y().at(i)->copy();
        auto s1 = rho().at(i);                       // Delta_i * delta_i
        auto s2 = detail::real(t3->dot_product(t2)); // Y_i     * delta_i
        auto s3 = t1->dot_product(value);            // Delta_i * v
        auto s4 = t3->dot_product(value);            // Y_i     * v
        out->ax_plus_y( s3 * s1, t1);
        out->ax_plus_y(-s4 / s2, t3);
      }
      if (!update && !Y().empty()) {
        decrement_Y();
      }
      return out;
    }


    // returns displacement using two-loop formula : no y_ intermediate just requires delta and Delta and vector to apply on
    // complex generalization : Sorber, van Barel, de Lathauwer, SIAM J. Optim. 22 (3) 379 (2012)
    std::shared_ptr<T> two_loop_inverse_hessian(std::shared_ptr<const T> _vector) {
      // to make sure, inputs are copied.
      auto vector = std::make_shared<const T>(*_vector);
      auto out = std::make_shared<T>(*vector);

      // set dimension
      const int n = delta().size();
      // first loop
      std::vector<double> alpha;
      if (n > 0) {
        const int np = n-1;
        for (int i = np; i > -1; --i) {
          auto deltai = delta().at(i)->copy();
          auto yi     = D().at(i)->copy();
          auto s1     = detail::real(deltai->dot_product(out));
          auto s2     = rho().at(i);
          auto alphai = s1 * s2;
          out->ax_plus_y(-alphai, yi);
          alpha.insert(alpha.begin(), alphai);
        }
      }

      // apply h0^-1
      *out /= *denom();

      // second loop
      for (int i = 0; i != n; ++i) {
        auto deltai = delta().at(i)->copy();
        auto yi = D().at(i)->copy();
        auto s1 = yi->dot_product(out);
        auto s2 = detail::real(yi->dot_product(deltai));
        auto b  = alpha.at(i) - s1/s2;
        out->ax_plus_y(b, deltai);
      }

      return out;
    }


    // returns restricted step displacement
    std::shared_ptr<T> extrapolate(std::shared_ptr<const T> _grad, std::shared_ptr<const T> _value) {
      bool with_shift = false;
      // to make sure, inputs are copied
      auto grad  = std::make_shared<const T>(*_grad);
      auto value  = std::make_shared<const T>(*_value);
      assert(trust_radius_ != 0.0);

      // compute Newton step and compare norm to trust radius
      auto acopy = two_loop_inverse_hessian(grad);
      auto anorm = acopy->norm();
      std::cout << std::setprecision(6) << " Initial Step Length = " << anorm << std::endl;
      if (anorm  > trust_radius_) {
        std::cout << std::setprecision(6) << " step norm = " << anorm << " | trust_radius = " << trust_radius_ << std::endl;
        std::cout << " NEWTON STEP NORM EXCEEDS THE TRUST RADIUS : Level shifting will be used " << std::endl;
        with_shift = true;
      }
      if (with_shift) {
        std::cout << " +++ Hebden algorithm used to determine level shift +++ " << std::endl;
        auto tshift = hebden_levelshift(grad);
        auto shift = grad->clone();
        shift->fill(tshift);
        acopy = level_shift_inverse_hessian(grad, shift);
        std::cout << " Level Shift = " << tshift << std::endl;
        std::cout << " Step Length = " << acopy->norm() << " | trust radius = " << trust_radius_ << std::endl;
      }
      acopy->scale(-1.0);
      if (!with_shift) {
        level_shift_ = 0.0;
      }
      prev_level_shift_ = level_shift_;
      prev_grad_ = grad;
      prev_value_ = value;

      return acopy;
    }


    // More-Sorensen method to compute level shift and displacement as discussed in Erway arXiv:1212.1525
    std::shared_ptr<T> more_sorensen_extrapolate(std::shared_ptr<const T> _grad, std::shared_ptr<const T> _value) {
      // to be sure; inputs are copied
      auto grad  = std::make_shared<T>(*_grad);
      auto value  = std::make_shared<T>(*_value);
      assert(trust_radius_ != 0.0);

      auto p = two_loop_inverse_hessian(grad);
      p->scale(-1.0);
      auto phat = p->clone();
      double tshift = 0.0;

      std::cout << std::setprecision(6) << " Initial Step Length = " << p->norm() << std::endl;
      for (int k = 0; k != hebden_iter_; ++k) {
        if (p->norm() <= trust_radius_) {
          std::cout << " More-Sorensen converged in " << k << " iterations. " << std::endl;
          std::cout << " Level Shift = " << tshift << std::endl;
          std::cout << " Step Length = " << p->norm() << " | trust radius = " << trust_radius_ << std::endl;
          break;
        } else {
          auto phi = 1.0 / p->norm() - 1.0 / trust_radius_;
          auto shift = p->clone();
          if (tshift < 1.0e-8) {
            phat = two_loop_inverse_hessian(p);
            phat->scale(-1.0);
            tshift = 0.0;
          } else {
            shift->fill(tshift);
            phat = level_shift_inverse_hessian(p, shift);
            phat->scale(-1.0);
          }
          auto phiprime = detail::real(p->dot_product(phat));
          phiprime /= -1.0 * pow(p->norm(), 3);
          tshift = tshift - phi / phiprime;
          if (tshift < 1.0e-8) {
            p = two_loop_inverse_hessian(grad);
            p->scale(-1.0);
            tshift = 0.0;
          } else {
            shift->fill(tshift);
            p = level_shift_inverse_hessian(grad, shift);
            p->scale(-1.0);
          }
          if (k == hebden_iter_ - 1) {
            std::cout << "+++ More-Sorensen extrapolate did NOT converge +++ " << std::endl;
            std::cout << " step norm with level shift = " << p->norm() << std::endl;
            std::cout << " Trying Hebden's method " << std::endl;
            tshift = hebden_levelshift(grad);
            shift->fill(tshift);
            p = level_shift_inverse_hessian(grad, shift);
            p->scale(-1.0);
            if (tshift < 0.0) {
              std::cout << " Hebden found negative level shift = " << tshift << std::endl;
            }
          }
        }
      }
      level_shift_ = tshift;
      prev_level_shift_ = level_shift_;
      prev_value_ = value;
      prev_grad_ = grad;

      return p;
    }


    std::shared_ptr<T> conjugate_gradient(std::shared_ptr<const T> _grad, std::shared_ptr<const T> _value) {
      // to be sure; inputs are copied
      auto grad  = std::make_shared<T>(*_grad);
      auto value  = std::make_shared<T>(*_value);
      auto out = std::make_shared<T>(*_grad);
      out->scale(-1.0);

      prev_grad_ = grad;
      prev_value_ = _value;

      return out;
    }


    // check if a step satisfies the quadratic approximation ; if not return bool to tell user to reset to previous expansion point
    bool check_step(std::vector<double> _f, std::shared_ptr<T> _grad, std::shared_ptr<T> _value, const bool tight = false, const int limited_memory = 0) {
      // f, grad, and value may all be modified by this function if the quadratic approximation fails for the current expansion point (CEP)
      bool reset = false;

      // update intermediates
      if (prev_value() != nullptr) {
        if (limited_memory > 0 && delta().size() > limited_memory) {
          std::cout << " Limited Memory : keeping the " << limited_memory << " most recent vectors " << std::endl;
          assert(limited_memory >= 1);
          D_.erase(D_.begin());
          delta_.erase(delta_.begin());
          avec_.erase(avec_.begin());
          rho_.erase(rho_.begin());
        }
        auto DD = std::make_shared<T>(*_grad - *prev_grad());
        D_.push_back(DD);
        auto yy = std::make_shared<T>(*_value - *prev_value());
        delta_.push_back(yy);
        {
          auto xx = std::make_shared<T>(*prev_grad());
          xx->scale(-1.0);
          auto xx2 = std::make_shared<T>(*yy);
          xx2->scale(prev_level_shift());
          *xx -= *xx2;
          avec_.push_back(xx);
        }
        auto rr = 1.0 / detail::real(DD->dot_product(yy));
        rho_.push_back(rr);
      }

      if (_f.size() > 1) {
        if (_f.at(_f.size()-1) > _f.at(_f.size()-2))
          std::cout << " OBJECTIVE FUNCTION VALUE INCREASED BY PREVIOUS STEP " << std::endl;
      }

      if (prev_grad() != nullptr) {
        std::cout << std::setprecision(10) << " previous gradient norm   =  " << prev_grad()->rms() << std::endl;
        std::cout << std::setprecision(10) << " current  gradient norm   =  " << _grad->rms() << std::endl;
      }

      if (!delta().empty()) {
        if (tight)
          std::cout << " Tight optimization specified " << std::endl;
        double rk = taylor_series_validity_ratio(_f, prev_grad(), delta().back());
        std::cout << std::setprecision(4) << " Taylor expansion validity parameter  = " << rk << std::endl;
#if 1
        if (!tight) {
          if (rk < 0.25 && rk > 0) {
            std::cout << " condition (i) satisfied " << std::endl;
            std::cout << " scaling down the trust radius " << std::endl;
            trust_radius_ = trust_radius_ * 2.0/3.0;
            std::cout << std::setprecision(8) << " trust radius   = " << trust_radius_ << std::endl;
          } else if ((rk >= 0.25) && (rk < 0.75)) {
            std::cout << " condition (ii) satisfied  : trust radius kept from previous iteration" << std::endl;
            trust_radius_ = trust_radius_;
            std::cout << std::setprecision(8) << " trust radius   = " << trust_radius_ << std::endl;
          } else if (rk > 0.75) {
            std::cout << " condition (iii) satisfied : trust radius increased " << std::endl;
            trust_radius_ = std::min(alpha_ * trust_radius_, 0.75);
            std::cout << std::setprecision(8) << " trust radius   = " << trust_radius_ << std::endl;
          } else if (rk < 0) {
            std::cout << " step does not satisfy Taylor expansion criteria " << std::endl;
            std::cout << " scaling down the trust radius " << std::endl;
            trust_radius_ = trust_radius_ / alpha_;
            reset = true;
          }
        } else { // conditions for tighter optimization
          if (((rmin_ < rk) && (rk < rgood_)) || ((2.0 - rgood_ < rk) && (rk < 2.0 - rmin_))) {
            std::cout << " condition (i) satisfied " << std::endl;
            trust_radius_ = trust_radius_;
            std::cout << std::setprecision(8) << " trust radius   = " << trust_radius_ << std::endl;
          } else if ((rgood_ < rk) && (rk < 2.0 - rgood_) ) {
            std::cout << " condition (ii) satisfied " << std::endl;
            trust_radius_ = std::min(alpha_ * trust_radius_, 0.75);
            std::cout << std::setprecision(8) << " trust radius   = " << trust_radius_ << std::endl;
          } else if (rk < rmin_ || rk > 2.0 - rmin_ ) {
            std::cout << " step does not satisfy Taylor expansion criteria " << std::endl;
            std::cout << " scaling down the trust radius " << std::endl;
            trust_radius_ = trust_radius_ / alpha_;
            reset = true;
          }
        }
#endif
      }

      return reset;
    }
};

}

#endif
