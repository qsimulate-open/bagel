//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: step_restrict_bfgs.h
// Copyright (C) 2014 Toru Shiozaki
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


#ifndef __SRC_UTIL_SRBFGS_H
#define __SRC_UTIL_SRBFGS_H

// Wrapper to Alglib

#include <stddef.h>
#include <thread>
#include <vector>
#include <memory>
#include <src/alglib/optimization.h>

namespace bagel {

template<typename T>
class SRBFGS {
  protected:
    std::atomic_bool flag_;
    alglib::minlbfgsstate state_;
    alglib::real_1d_array denom_;

    std::shared_ptr<T> current_;
    std::shared_ptr<const T> gradient_;
    double energy_;

    void evaluate(const alglib::real_1d_array& x, double& en, alglib::real_1d_array& grad, void* ptr) {
      // set new x to the member 
      for (int i = 0; i != current_->size(); ++i)
        *(current_->data()+i) = x[i];
      flag_ = true;
      // then release the mutex and sleep till mutex2 is released
      while (true) {
        if (flag_) {
          std::this_thread::sleep_for(std::chrono::milliseconds(1));
          continue;
        }
        break;
      }
      // set gradient and energy_
      for (int i = 0; i != gradient_->size(); ++i)
        grad[i] = *(gradient_->data()+i);
      en = energy_;
    }

    std::shared_ptr<std::thread> server_;

    using eval_type = std::function<void(const alglib::real_1d_array&, double&, alglib::real_1d_array&, void*)>;

  public:
    SRBFGS(std::shared_ptr<const T> denom, const bool dummy = false) : current_(denom->clone()) {
      denom_.setcontent(denom->size(), denom->data());
    }

    ~SRBFGS() {
      
    }

    // _value is supposed to be identical to cunrret_
    std::shared_ptr<T> extrapolate(std::shared_ptr<const T> _grad, std::shared_ptr<const T> _value, const double en) {
    std::cout << "extrapolate is called" << std::endl;
      // start the server thread in the first iteration
      if (!server_) {
        alglib::real_1d_array x;
        x.setcontent(_value->size(), _value->data()); 
        eval_type eval = std::bind(&SRBFGS<T>::evaluate, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);

        alglib::minlbfgsreport rep;

        alglib::minlbfgscreate(1, x, state_); 
        alglib::minlbfgssetcond(state_, /*essentially zero*/1.0e-15, 0.0, 0.0, /*essentially infty*/1000);
        alglib::minlbfgssetstpmax(state_, /*maxstep*/ 0.0);
        alglib::minlbfgssetprecdiag(state_, denom_);

        flag_ = false;
        server_ = std::make_shared<std::thread>([&eval,this](){ alglib::minlbfgsoptimize(state_, eval); });
        // wait till the server thread lock the mutex_
        while (true) {
          if (!flag_) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            continue;
          }
          break;
        }
      }
      energy_ = en;
      gradient_ = _grad;
      flag_ = false;

      // wait till the server thread lock the mutex_; then a new vector is ready in current_
      while (true) {
        if (!flag_) {
          std::this_thread::sleep_for(std::chrono::milliseconds(1));
          continue;
        }
        break;
      }
      return current_->copy();
    }
};

}

#endif
