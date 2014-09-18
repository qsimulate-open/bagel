//
// BAGEL - Parallel electron correlation program.
// Filename: kronecker.h
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef __SRC_ASD_DMRG_KRONECKER_H
#define __SRC_ASD_DMRG_KRONECKER_H

#include <src/math/matop.h>

namespace bagel {

template <class T,
          class = typename std::enable_if<detail::is_mat<T>::value>::type
         >
T kronecker_product(const bool atrans, const T& A, const bool btrans, const T& B) {
  const int n = atrans ? A.mdim() : A.ndim();
  const int m = atrans ? A.ndim() : A.mdim();
  const int p = btrans ? B.mdim() : B.ndim();
  const int q = btrans ? B.ndim() : B.mdim();

  std::shared_ptr<const Matrix> BB = btrans ? B.transpose() : nullptr;
  const Matrix* Bptr = btrans ? BB.get() : &B;

  T out(n*p, m*q);
  for (int ja = 0; ja < m; ++ja) {
    for (int ia = 0; ia < n; ++ia) {
      const double aval = atrans ? A(ja, ia) : A(ia, ja);
      out.add_block(aval, ia*p, ja*q, p, q, Bptr->data());
    }
  }

  return out;
}

}

#endif
