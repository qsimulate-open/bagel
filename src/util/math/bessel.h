//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: bessel.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

//

#ifndef __BAGEL_UTIL_MATH_BESSEL_H
#define __BAGEL_UTIL_MATH_BESSEL_H

#include <iostream>
#include <cmath>
#include <vector>
#include <src/util/math/factorial.h>

namespace bagel {

class MSphBesselI {

  private:
    const Factorial f;

    double R_l(const double x, const int l) const {
      double sum = 0.0;
      for (int i = 0; i <= l; ++i) sum += f(l + i) / (f(i) * f(l -  i) * std::pow(2.0 * x, i));
      return sum;
    }

  public:
    MSphBesselI() : f() {}
    ~MSphBesselI() {}

    double compute(const int l, const double x) const {

      double bessel = 0.0;

      if (x < 1e-7) {
        bessel = (1.0 - x) * std::pow(2.0 * x, l) * f(l) / f(2 * l);
      } else if (x > 16) {
        bessel = 0.5 * R_l(-x, l) / x;
      } else {
        const double halfxsq = 0.5 * x * x;
        double current = std::exp(-x);
        if (l != 0)
          for (int ll = 1; ll <= l; ++ll) current *= x / (2*ll+1);

        int j = 0;
        do {
            bessel += current;
            ++j;
            current *= halfxsq / (j*(2*l + 2*j + 1));
        } while (std::fabs(current) > std::numeric_limits<double>::epsilon());
      }

      return bessel;
    }

};

}

#endif
