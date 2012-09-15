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
#include <src/rysint/macros.h>

namespace bagel {

class Int2D {

  protected:

    /// for recursion
    double C00_[RYS_MAX];
    double D00_[RYS_MAX];
    double B00_[RYS_MAX];
    double B10_[RYS_MAX];
    double B01_[RYS_MAX];

    // two index intermediate
    int rank_;
    double* data_;
    int datasize_;

    ScaleList scale_;

  public:

    Int2D(const std::array<double, 11>&, const double*, const int, const int, double*,
          void (*vrrfunc)(double*, const double*, const double*, const double*, const double*, const double*));
    Int2D() {};
    ~Int2D();

    const double* data() const { return data_; };
    int datasize() const { return datasize_; };

    void scale_data(const double* a, const double c) {
      scale_.scalefunc[rank_](data_, a, c, data_, datasize_);
    };

    void scale_data_t(double* target, const double* a, const double c) {
      scale_.scalefunc[rank_](target, a, c, data_, datasize_);
    };
};

}

#endif
