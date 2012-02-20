//
// Newint - Parallel electron correlation program.
// Filename: hpw_diis.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


// implements special class of DIIS
// Hampel, Peterson, Werner, Chem. Phys. Lett. 190, 1 (1992).
//
// In addition to the requirement in diis.h, T should have
// functions exp(), log(), unit() that return shared_ptr<T>.
// Also operator "*=" should be overloaded in T.
// Copy constructor is needed as well.

#ifndef __NEWINT_UTIL_HPW_DIIS_H
#define __NEWINT_UTIL_HPW_DIIS_H

#include <src/util/diis.h>

template<class T>
class HPW_DIIS  {
  typedef std::shared_ptr<T> RefT;
  protected:
    DIIS<T> diis_;
    RefT base_;
    const RefT orig_;

  public:
    HPW_DIIS(const int n, RefT o) : diis_(n), base_(o->clone()), orig_(o) { base_->unit(); };
    ~HPW_DIIS() {};

    RefT extrapolate(const RefT rot) {
      // prev = log(base)
      RefT expo = (*base_**rot).log();
      RefT prev_ = base_->log();
      RefT err(new T(prev_ ? (*expo-*prev_) : (*expo)));

      RefT extrap = diis_.extrapolate(std::make_pair(expo, err))->exp();
      // this is important
      extrap->purify_unitary();
      // returns unitary matrix with respect to the original matrix
      RefT out(new T(*orig_* *extrap));
      *base_ = *extrap;
      return out;
    };

};

#endif
