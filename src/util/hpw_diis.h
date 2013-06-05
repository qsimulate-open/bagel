//
// BAGEL - Parallel electron correlation program.
// Filename: hpw_diis.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


// implements special class of DIIS
// Hampel, Peterson, Werner, Chem. Phys. Lett. 190, 1 (1992).
//
// In addition to the requirement in diis.h, T should have
// functions exp(), log(), unit() that return shared_ptr<T>.
// Also operator "*=" should be overloaded in T.
// Copy constructor is needed as well.

#ifndef __BAGEL_UTIL_HPW_DIIS_H
#define __BAGEL_UTIL_HPW_DIIS_H

#include <src/util/diis.h>
#include <src/util/bfgs.h>

namespace bagel {

template<class T>
class HPW_DIIS  {
  using RefT = std::shared_ptr<const T>;
  protected:
    DIIS<T> diis_;
    RefT base_;
    const RefT orig_;
    const bool testing_;

  public:
    HPW_DIIS(const int n, RefT o, const bool test = false) : diis_(n), orig_(o), testing_(test) {
      std::shared_ptr<T> b = o->clone();
      b->unit();
      base_ = b;
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

      std::shared_ptr<T> extrap = diis_.extrapolate(std::make_pair(expo, err))->exp(100);
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
