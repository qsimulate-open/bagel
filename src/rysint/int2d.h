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

#ifndef __src_rysint_int2d_h
#define __src_rysint_int2d_h

#include <vector>
#include <map>
#include <array>
#include <src/rysint/scalelist.h>

namespace bagel {

template <int RANK>
class Int2D {

  protected:
    /// for recursion
    double C00_[RANK];
    double D00_[RANK];
    double B00_[RANK];
    double B10_[RANK];
    double B01_[RANK];

    // two index intermediate
    double* const data_;
    const int datasize_;

    const ScaleList scale_;

  public:
    Int2D(const std::array<double,11>& dparam, const double* roots, const int datasize, double* const data_pointer,
       void (*vrrfunc)(double*, const double*, const double*, const double*, const double*, const double*))
     : data_(data_pointer), datasize_(datasize) {

      const double P = dparam[0];
      const double Q = dparam[1];
      const double A = dparam[2];
      const double B = dparam[3];
      const double C = dparam[4];
      const double D = dparam[5];
      const double xp = dparam[6];
      const double xq = dparam[7];

      const double one_2p = dparam[8];
      const double one_2q = dparam[9];

      const double one_pq = dparam[10];
      const double xqopq = xq * one_pq;
      const double xpopq = xp * one_pq;

      const double c00i0 = P - A;
      const double c00i1 = (P - Q) * xqopq;
      const double d00i0 = Q - C;
      const double d00i1 = (P - Q) * xpopq;
      const double b00i0 = 0.5 * one_pq;
      const double b10i0 = xqopq * one_2p;
      const double b01i0 = xpopq * one_2q;

      for (int i = 0; i != RANK; ++i) {
        const double tsq = roots[i];
        C00_[i] = c00i0 - c00i1 * tsq;
        D00_[i] = d00i0 + d00i1 * tsq;
        B00_[i] = b00i0 * tsq;
        B10_[i] = one_2p - b10i0 * tsq;
        B01_[i] = one_2q - b01i0 * tsq;
      }

      vrrfunc(data_, C00_, D00_, B00_, B01_, B10_);
    }



    const double* data() const { return data_; }
    int datasize() const { return datasize_; }

    void scale_data(const double* a, const double c) {
      scale_.scalefunc[RANK](data_, a, c, data_, datasize_);
    }

    void scale_data_t(double* target, const double* a, const double c) {
      scale_.scalefunc[RANK](target, a, c, data_, datasize_);
    }
};

}

#endif
