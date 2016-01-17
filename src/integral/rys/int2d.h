//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: int2d.h
// Copyright (C) 2009 Toru Shiozaki
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

#ifndef __SRC_RYSINT_INT2D_H
#define __SRC_RYSINT_INT2D_H

#include <src/integral/rys/_vrr.h>

namespace bagel {

template <int a_, int c_, int rank_, typename DataType>
void int2d(const DataType& P, const DataType& Q, const DataType& A, const DataType& B, const DataType& C, const DataType& D,
           const double& xp, const double& xq, const double& one_2p, const double& one_2q, const double& one_pq, const DataType* roots, DataType* const data) {
  /// for recursion
  alignas(32) DataType C00_[rank_];
  alignas(32) DataType D00_[rank_];
  alignas(32) DataType B00_[rank_];
  alignas(32) DataType B10_[rank_];
  alignas(32) DataType B01_[rank_];

  const double xqopq = xq * one_pq;
  const double xpopq = xp * one_pq;

  const DataType c00i0 = P - A;
  const DataType c00i1 = (P - Q) * xqopq;
  const DataType d00i0 = Q - C;
  const DataType d00i1 = (P - Q) * xpopq;
  const double b00i0 = 0.5 * one_pq;
  const double b10i0 = xqopq * one_2p;
  const double b01i0 = xpopq * one_2q;

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
