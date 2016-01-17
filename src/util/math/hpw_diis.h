//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: hpw_diis.h
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


// implements special class of DIIS
// Hampel, Peterson, Werner, Chem. Phys. Lett. 190, 1 (1992).
//
// In addition to the requirement in diis.h, T should have
// functions exp(), log(), unit() that return shared_ptr<T>.
// Also operator "*=" should be overloaded in T.
// Copy constructor is needed as well.

#ifndef __BAGEL_UTIL_HPW_DIIS_H
#define __BAGEL_UTIL_HPW_DIIS_H

#include <src/util/math/diis.h>
#include <src/util/math/bfgs.h>

namespace bagel {

template<class T, class Mat = Matrix,
         class = typename std::enable_if<std::is_same<typename Mat::data_type, typename T::data_type>::value>::type
        >
class HPW_DIIS  {
  using RefT = std::shared_ptr<const T>;
  protected:
    DIIS<T,Mat> diis_;
    RefT base_;
    const RefT orig_;
    const bool testing_;

  public:
    HPW_DIIS(const int n, RefT o, RefT init, const bool test = false) : diis_(n), base_(init), orig_(o), testing_(test) {
    }

    RefT extrapolate(const RefT rot, const RefT errin = RefT()) {
      // prev = log(base)
      RefT expo = (*base_**rot).log(100);
      RefT prev_ = base_->log(100);
      RefT err;
      if (errin) {
        err = errin;
      } else {
        err = std::make_shared<T>(prev_ != nullptr ? (*expo-*prev_) : (*expo));
      }

      std::shared_ptr<T> extrap = diis_.extrapolate({expo, err})->exp(100);
      // this is important
      extrap->purify_unitary();
      base_ = extrap;
      // returns unitary matrix with respect to the original matrix
      return std::make_shared<T>(*orig_* *extrap);
    }

    RefT extrap() const { return base_; }
    RefT start() const { return orig_; }

};

}

#endif
