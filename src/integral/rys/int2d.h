//
// BAGEL - Parallel electron correlation program.
// Filename: int2d.h
// Copyright (C) 2009 Toru Shiozaki
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

#ifndef __SRC_RYSINT_INT2D_H
#define __SRC_RYSINT_INT2D_H

#include <src/integral/rys/_vrr.h>

namespace bagel {

template <int a_, int c_, int rank_, typename DataType>
void int2d(const DataType& P, const DataType& Q, const DataType& A, const DataType& B, const DataType& C, const DataType& D,
           const DataType& xp, const DataType& xq, const DataType& one_2p, const DataType& one_2q, const DataType& one_pq, const DataType* roots, DataType* const data) {
  /// for recursion
  alignas(32) DataType C00_[rank_];
  alignas(32) DataType D00_[rank_];
  alignas(32) DataType B00_[rank_];
  alignas(32) DataType B10_[rank_];
  alignas(32) DataType B01_[rank_];

  const DataType xqopq = xq * one_pq;
  const DataType xpopq = xp * one_pq;

  const DataType c00i0 = P - A;
  const DataType c00i1 = (P - Q) * xqopq;
  const DataType d00i0 = Q - C;
  const DataType d00i1 = (P - Q) * xpopq;
  const DataType b00i0 = 0.5 * one_pq;
  const DataType b10i0 = xqopq * one_2p;
  const DataType b01i0 = xpopq * one_2q;

  for (int i = 0; i != rank_; ++i) {
    const DataType tsq = roots[i];
    C00_[i] = c00i0 - c00i1 * tsq;
    D00_[i] = d00i0 + d00i1 * tsq;
    B00_[i] = b00i0 * tsq;
    B10_[i] = one_2p - b10i0 * tsq;
    B01_[i] = one_2q - b01i0 * tsq;
  }

  vrr<a_,c_,rank_, DataType>(data, C00_, D00_, B00_, B01_, B10_);
}

}

#endif
