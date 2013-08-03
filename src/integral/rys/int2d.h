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

template <int a_, int c_, int rank_>
void int2d(const double& P, const double& Q, const double& A, const double& B, const double& C, const double& D,
           const double& xp, const double& xq, const double& one_2p, const double& one_2q, const double& one_pq, const double* roots, double* const data) {
  /// for recursion
#if __GNUC__ == 4 && __GNUC_MINOR__ <= 7
  double C00_[rank_]__attribute__((aligned(32))); // TODO deprecated
  double D00_[rank_]__attribute__((aligned(32)));
  double B00_[rank_]__attribute__((aligned(32)));
  double B10_[rank_]__attribute__((aligned(32)));
  double B01_[rank_]__attribute__((aligned(32)));
#else
  alignas(32) double C00_[rank_];
  alignas(32) double D00_[rank_];
  alignas(32) double B00_[rank_];
  alignas(32) double B10_[rank_];
  alignas(32) double B01_[rank_];
#endif

  const double xqopq = xq * one_pq;
  const double xpopq = xp * one_pq;

  const double c00i0 = P - A;
  const double c00i1 = (P - Q) * xqopq;
  const double d00i0 = Q - C;
  const double d00i1 = (P - Q) * xpopq;
  const double b00i0 = 0.5 * one_pq;
  const double b10i0 = xqopq * one_2p;
  const double b01i0 = xpopq * one_2q;

  for (int i = 0; i != rank_; ++i) {
    const double tsq = roots[i];
    C00_[i] = c00i0 - c00i1 * tsq;
    D00_[i] = d00i0 + d00i1 * tsq;
    B00_[i] = b00i0 * tsq;
    B10_[i] = one_2p - b10i0 * tsq;
    B01_[i] = one_2q - b01i0 * tsq;
  }

  vrr<a_,c_,rank_>(data, C00_, D00_, B00_, B01_, B10_);
}

}

#endif
