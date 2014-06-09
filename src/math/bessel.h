//
// BAGEL - Parallel electron correlation program.
// Filename: bessel.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

//

#ifndef __BAGEL_UTIL_BESSEL_H
#define __BAGEL_UTIL_BESSEL_H

#include <iostream>
#include <cmath>
#include <vector>
#include <src/math/factorial.h>

namespace bagel {

class MSphBesselI {

  private:
    int l_;

    const static Factorial f;

    double R_l(const double x) const {
      double sum = 0.0;
      for (int i = 0; i <= l_; ++i) sum += f(l_ + i) / (f(i) * f(l_ -  i) * std::pow(2.0 * x, i));
      return sum;
    }

  public:
    MSphBesselI(const int l) : l_(l) {}
    ~MSphBesselI() {}

    double compute(const double x) const {

      double bessel = 0.0;

      if (x < 1e-7) {
        bessel = (1.0 - x) * std::pow(2.0 * x, l_) * f(l_) / f(2 * l_);
      } else if (x > 16) {
        bessel = 0.5 * R_l(-x) / x;
      } else {
        const double halfxsq = 0.5 * x * x;
        double current = std::exp(-x);
        if (l_ != 0)
          for (int ll = 1; ll <= l_; ++ll) current *= x / (2*ll+1);

        int j = 0;
        do {
            bessel += current;
            ++j;
            current *= halfxsq / (j*(2*l_ + 2*j + 1));
        } while (std::fabs(current) > std::numeric_limits<double>::epsilon());
      }

      return bessel;
    }

};

}

#endif
