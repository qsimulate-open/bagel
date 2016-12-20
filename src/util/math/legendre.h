//
// BAGEL - Parallel electron correlation program.
// Filename: legendre.h
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

#ifndef __BAGEL_MATH_LEGENDRE_H
#define __BAGEL_MATH_LEGENDRE_H

#include <stdexcept>
#include <src/util/constants.h>

namespace bagel {

class Legendre {

  private:

  public:
    Legendre() {}
    ~Legendre() {}

    double compute(const int l, const int am, const double x) const {

      if (am < 0 || am > l || fabs(x) > 1.0)
        throw std::runtime_error("SH: abs(m) must be in [0, l] and x in [-1, 1]");
      double pmm = 1.0;
      if (am > 0) {
        double somx2 = sqrt((1.0 - x)*(1.0 + x));
        double fact = 1.0;
        for (int i = 1; i <= am; ++i) {
          pmm *= -fact * somx2;
          fact += 2.0;
        }
      }
      if (l == am) {
        return pmm;
      } else {
        double pmmp1 = x * (2.0 * am + 1) * pmm;
        if (l == am+1) {
          return pmmp1;
        } else {
          double plm = 0.0;
          for (int i = am + 2; i <= l; ++i) {
            plm = (x * (2 * i -1) * pmmp1 - (i + am - 1) * pmm) / (i - am);
            pmm = pmmp1;
            pmmp1 = plm;
          }
          return plm;
        }
      }
    }
};


class Legendre_renorm {

  public:
    Legendre_renorm() { }
    ~Legendre_renorm() { }


    double compute(const int l, const int am, const double x) const {
      if (am < 0 || am > l || fabs(x) > 1.0)
        throw std::runtime_error("Legendre: m must be in [0, l] and x in [-1, 1]");
      double pmm = 1.0;
      if (am > 0) {
        double omx2 = (1.0 - x) * (1.0 + x);
        double fact = 1.0;
        for (int i = 1; i <= am; ++i) {
          pmm *= omx2 * fact / (fact + 1.0);
          fact += 2.0;
        }
      }
      pmm = sqrt((2 * am + 1) * pmm / (4.0 * pi__));
      if (am & 1)
        pmm = -pmm;
      if (l == am) {
        return pmm;
      } else {
        double pmmp1 = x * sqrt(2.0 * am + 3.0) * pmm;
        if (l == am + 1) {
          return pmmp1;
        } else {
          double oldfact = sqrt(2.0 * am + 3.0);
          double pll = 0.0;
          for (int ll = am + 2; ll <= l; ++ll) {
            double fact = sqrt((4.0 * ll * ll - 1.0) / (ll * ll - am * am));
            pll = (x * pmmp1 - pmm / oldfact) * fact;
            oldfact = fact;
            pmm = pmmp1;
            pmmp1 = pll;
          }
          return pll;
        }
      }
    }


    // should agree with Legendre's compute(l, am, x)
    double compute_plm(const int l, const int am, const double x) const {
      if (am < 0 || am > l || abs(x) > 1.0)
        throw std::runtime_error("Legendre: bad arguments");

      double coef = 1.0;
      for (int j=l-am+1; j<=l+am; ++j)
        coef *= j;

      return sqrt(4.0*pi__*coef/(2*l+1)) * compute(l,am,x);
    }
};

}

#endif
