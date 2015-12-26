//
// BAGEL - Parallel electron correlation program.
// Filename: gamma.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __BAGEL_UTIL_MATH_GAMMA_H
#define __BAGEL_UTIL_MATH_GAMMA_H

#include <iostream>
#include <cmath>
#include <vector>
#include <src/util/constants.h>
#include <src/util/math/factorial.h>
#include <boost/math/special_functions/gamma.hpp>

namespace bagel {

class Gamma { // special case for half integers G(l + 1/2)

  private:

  public:
    Gamma() { }
    ~Gamma() { }

    double compute(const int l) const {
      const double sqrtpi = sqrt(pi__);
      double out = 1.0;
      double ft = 1.0;
      for (int i = 1; i <= l; ++i) {
         out *= ft / 2.0;
         ft += 2.0;
      }

      return out * sqrtpi;
    }

    double operator()(const int l) const { return compute(l); }
};


class Gamma_scaled { // r^l/G(l + 1/2)

  private:

  public:
    Gamma_scaled() { }
    ~Gamma_scaled() { }

    double compute(const int l, const double r) const {
      const double sqrtpi = sqrt(pi__);
      double out = 1.0;
      double ft = 1.0;
      for (int i = 1; i <= l; ++i) {
         out *= 2.0 * r / ft;
         ft += 2.0;
      }

      return out / sqrtpi;
    }

    double operator()(const int l, const double r) const { return compute(l, r); }
};


class Gamma_upper { // Upper incomplete gamma function G(l+1/2, x)

  private:

  public:
    Gamma_upper() { }
    ~Gamma_upper() { }

    double compute(const int l, const double x) const {

      if (l > 1200 || x > 300)
        throw std::runtime_error("Failed to compute incomplete gamma function!");

      double gamma = sqrt(pi__) * erfc(sqrt(x));
      for (int i = 1; i <= l; ++i) {
        const double r = i + 0.5;
        gamma = (r-1.0) * gamma + pow(x, r-1.0) * exp(-x);
      }

      return gamma;
    }

    double operator()(const int l, const double x) const { assert(l < 1200 && x < 300); return compute(l, x); }

};


class Gamma_lower_scaled { // g(l+1/2, z^2)/x^{l+1} where z = beta^2 x^2

  private:
    double F_l(const int l, const double t) const {
      double out;
      if (t < 300) {
        const int istart = 1200;
        Gamma_upper gupper;
        const double glower = boost::math::tgamma(istart+0.5) - gupper(istart+0.5, t);
        double prev = glower / (2.0 * pow(t, istart+0.5));
        for (int i = istart-1; i <=l; --i) {
          out = (2.0 * t * prev + exp(-t)) / (2 * i + 1.0);
          prev = out;
        }
      } else {
        out = pow(pi__, 0.5) / (pow(2.0, l+1.0) * pow(t, l+0.5));
        double df = 1.0;
        for (int i = 1; i <= 2*l-1; i += 2) {
          out *= df;
          df += 2.0;
        }
      }

      return out;
    }

  public:
    Gamma_lower_scaled() { }
    ~Gamma_lower_scaled() { }

    double compute(const int l, const double z, const double beta) const {
      const double f = F_l(l, z);
      const double rsq = z / (beta * beta);
      const double out = 2.0 * pow(beta, 2*l+1.0) * pow(rsq, 0.5*l) * f;

      return out;
    }

    double operator() (const int l, const double z, const double beta) const {
      assert(z > 0);
      return compute(l, z, beta);
    }
};

}

#endif
