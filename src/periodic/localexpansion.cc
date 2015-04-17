//
// BAGEL - Parallel electron correlation program.
// Filename: localexpansion.cc
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


#include <src/periodic/localexpansion.h>

using namespace std;
using namespace bagel;

const static Legendre plm;
const static Factorial f;

LocalExpansion::LocalExpansion(const array<double, 3> c, vector<shared_ptr<const ZMatrix>> m, const int lmax)
 : centre_(c), moments_(m), lmax_(lmax) {

  nbasis1_ = m.front()->ndim();
  nbasis0_ = m.front()->mdim();
  num_multipoles_ = (lmax + 1) * (lmax + 1);
  assert(m.size() == num_multipoles_);
}


vector<shared_ptr<const ZMatrix>> LocalExpansion::compute_local_moments() {

  const double r = sqrt(centre_[0]*centre_[0] + centre_[1]*centre_[1] + centre_[2]*centre_[2]);
  const double ctheta = centre_[2]/r;
  const double phi = atan2(centre_[1], centre_[0]);

  vector<shared_ptr<const ZMatrix>> out(num_multipoles_);

  int i1 = 0;
  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m, ++i1) {

      ZMatrix local(nbasis1_, nbasis0_);
      int i2 = 0;
      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2 * j; ++k, ++i2) {

          const int a = l + j;
          const int b = m - l + k - j;

          const double prefactor = f(a - abs(b)) * plm.compute(a, abs(b), ctheta) / pow(r, a + 1);
          const double real = (b >= 0) ? (prefactor * cos(abs(b) * phi)) : (-1.0 * prefactor * cos(abs(b) * phi));
          const double imag = prefactor * sin(abs(b) * phi);
          const complex<double> coeff(real, imag);

          if (abs(coeff) > numerical_zero__)
            zaxpy_(nbasis0_ * nbasis1_, coeff, moments_[i2]->data(), 1, local.data(), 1);
        }
      }
      assert(i2 == num_multipoles_);
      out[i1] = make_shared<const ZMatrix>(local);
    }
  }

  return out;
}


/* given O(a) and centre (b-a) compute O(b) */
vector<shared_ptr<const ZMatrix>> LocalExpansion::compute_shifted_moments() {

  const double r = sqrt(centre_[0]*centre_[0] + centre_[1]*centre_[1] + centre_[2]*centre_[2]);
  const double ctheta = centre_[2]/r;
  const double phi = atan2(centre_[1], centre_[0]);

  vector<shared_ptr<const ZMatrix>> out(num_multipoles_);

  int i1 = 0;
  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m, ++i1) {

      ZMatrix shifted(nbasis1_, nbasis0_);
      int i2 = 0;
      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2 * j; ++k, ++i2) {

          const int a = l - j;
          const int b = m - l - k + j;
            const double prefactor = pow(r, a) * plm.compute(a, abs(b), ctheta) / f(a + abs(b));
          if (abs(b) <= a && a >= 0) {
            }
            const double real = (b >= 0) ? (prefactor * cos(abs(b) * phi)) : (-1.0 * prefactor * cos(abs(b) * phi));
            const double imag = prefactor * sin(abs(b) * phi);
            const complex<double> coeff(real, imag);

            if (abs(coeff) > numerical_zero__)
              zaxpy_(nbasis0_ * nbasis1_, coeff, moments_[i2]->data(), 1, shifted.data(), 1);
          }
        }
      }
      assert(i2 == num_multipoles_);
      out[i1] = make_shared<const ZMatrix>(shifted);
    }
  }

  return out;
}


/* given L(a) and centre (a-b) compute L(b) */
vector<shared_ptr<const ZMatrix>> LocalExpansion::compute_shifted_local_expansions() {

  const double r = sqrt(centre_[0]*centre_[0] + centre_[1]*centre_[1] + centre_[2]*centre_[2]);
  const double ctheta = centre_[2]/r;
  const double phi = atan2(centre_[1], centre_[0]);

  vector<shared_ptr<const ZMatrix>> out(num_multipoles_);

  int i1 = 0;
  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m, ++i1) {

      ZMatrix shifted(nbasis1_, nbasis0_);
      int i2 = 0;
      for (int j = 0; j <= lmax_; ++j) {
        for (int k = 0; k <= 2 * j; ++k, ++i2) {

          const int a = j - l;
          const int b = k - j - m + l;
            const double prefactor = pow(r, a) * plm.compute(a, abs(b), ctheta) / f(a + abs(b));
          if (abs(b) <= a && a >= 0) {
            }
            const double real = (b >= 0) ? (prefactor * cos(abs(b) * phi)) : (-1.0 * prefactor * cos(abs(b) * phi));
            const double imag = prefactor * sin(abs(b) * phi);
            const complex<double> coeff(real, imag);

            if (abs(coeff) > numerical_zero__)
              zaxpy_(nbasis0_ * nbasis1_, coeff, moments_[i2]->data(), 1, shifted.data(), 1);
          }
        }
      }
      assert(i2 == num_multipoles_);
      out[i1] = make_shared<const ZMatrix>(shifted);
    }
  }

  return out;
}
