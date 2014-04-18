//
// BAGEL - Parallel electron correlation program.
// Filename: sphharmonics.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <iostream>
#include <iomanip>
#include <src/integral/ecp/sphharmonics.h>

using namespace std;
using namespace bagel;

SphHarmonics::SphHarmonics(const array<int, 2> lm, const array<double, 3> c)
 : angular_momentum_(lm), centre_(c) {

  const double r = sqrt(centre_[0]*centre_[0] + centre_[1]*centre_[1] + centre_[2]*centre_[2]);
  theta_ = acos(centre_[2]/r);
  phi_ = atan2(centre_[1], centre_[0]);

}

SphHarmonics::SphHarmonics(const int l, const int m, const array<double, 3> c)
 : angular_momentum_{{l, m}}, centre_(c) {

  const double r = sqrt(centre_[0]*centre_[0] + centre_[1]*centre_[1] + centre_[2]*centre_[2]);
  theta_ = acos(centre_[2]/r);
  phi_ = atan2(centre_[1], centre_[0]);

}

SphHarmonics::SphHarmonics(const array<int, 2> lm)
 : angular_momentum_(lm), centre_{{0.0, 0.0, 0.0}} {

  theta_ = 0.0;
  phi_ = 0.0;
}

SphHarmonics::SphHarmonics(const int l, const int m)
 : angular_momentum_{{l, m}}, centre_{{0.0, 0.0, 0.0}} {

  theta_ = 0.0;
  phi_ = 0.0;
}

double SphHarmonics::LegendrePolynomial(const double x) const {

  const int l = angular_momentum_[0];
  const int m = fabs(angular_momentum_[1]);

  if (m < 0 || m > l || fabs(x) > 1.0) throw std::runtime_error("SH: m must be in [0, l] and x in [-1, 1]");
  double pmm = 1.0;
  if (m > 0) {
    double somx2 = std::sqrt((1.0 - x)*(1.0 + x));
    double fact = 1.0;
    for (int i = 1; i <= m; ++i) {
      pmm *= -fact * somx2;
      fact += 2.0;
    }
  }
  if (l == m) {
    return pmm;
  } else {
    double pmmp1 = x * (2.0 * m + 1) * pmm;
    if (l == m+1) {
      return pmmp1;
    } else {
      double plm = 0.0;
      for (int i = m + 2; i <= l; ++i) {
        plm = (x * (2 * i -1) * pmmp1 - (i + m - 1) * pmm) / (i - m);
        pmm = pmmp1;
        pmmp1 = plm;
      }
      return plm;
    }
  }

}

complex<double> SphHarmonics::ylm() const {

  const int l = angular_momentum_[0];
  const int m = angular_momentum_[1];

  const double cth = cos(theta_);
  if (fabs(m) > l) throw std::runtime_error ("SphHarmonics.ylm: |m| > l");

  const int am = fabs(m);
  const double plm = LegendrePolynomial(cth);
  double fact = 1.0;
  for (int i = 1; i <= 2*am; ++i) {
    fact *= l - am + i;
  }
  const double coef = sqrt((2*l+1) * 0.25/pi__ / fact);
  double real = coef * plm * cos(am * phi_);
  double imag = coef * plm * sin(am * phi_);
  if (m < 0) {
    real *= std::pow(-1, m);
    imag *= std::pow(-1, m+1);
  }

  return complex<double>(real, imag);

}

double SphHarmonics::zlm() const {

  const int l = angular_momentum_[0];
  const int m = angular_momentum_[1];

  const double cth = cos(theta_);
  if (fabs(m) > l) throw runtime_error ("SphHarmonics.zlm: |m| > l");

  const int am = fabs(m);
  const double plm = LegendrePolynomial(cth);
  double fact = 1.0;
  for (int i = 1; i <= 2*am; ++i) {
    fact *= l - am + i;
  }
  const double coef = sqrt((2*l+1) * 0.25/pi__ / fact);
  if (m == 0) {
    return coef * plm;
  } else if (m > 0) {
    return sqrt(2.0) * pow(-1.0, m) * coef * plm * cos(am*phi_);
  } else {
    return sqrt(2.0) * pow(-1.0, m) * coef * plm * sin(am*phi_);
  }

}

double SphHarmonics::zlm(const int l, const int m) const {

  const double cth = cos(theta_);
  if (fabs(m) > l) throw runtime_error ("SphHarmonics.zlm: |m| > l");

  const int am = fabs(m);
  const double plm = LegendrePolynomial(cth);
  double fact = 1.0;
  for (int i = 1; i <= 2*am; ++i) {
    fact *= l - am + i;
  }
  const double coef = sqrt((2*l+1) * 0.25/pi__ / fact);
  if (m == 0) {
    return coef * plm;
  } else if (m > 0) {
    return sqrt(2.0) * pow(-1.0, m) * coef * plm * cos(am*phi_);
  } else {
    return sqrt(2.0) * pow(-1.0, m) * coef * plm * sin(am*phi_);
  }

}

void SphHarmonics::print() {
  cout << "Computing Real Spherical Harmonics" << endl;
  cout << "Angular momentum (lm) = (" << angular_momentum_[0] << ", "
                                      << angular_momentum_[1] << ")" << endl;
  cout << "Centre = ";
  for (int i = 0; i != 3; ++i) cout << setw(17) << setprecision(9) << centre_[i] << ";   ";
  cout << endl;
}

double SphHarmonics::sph_to_USP(const int lx, const int ly) const {

  const int am = fabs(angular_momentum_[1]);
  const int j = (lx + ly - am) / 2;
  if ((lx + ly - am) % 2 != 0 || j < 0) {
    return 0.0;
  } else {
    Factorial fact;
    double flam = 1.0;
    for (int i = 1; i <= 2*am; ++i) {
      flam *= angular_momentum_[0] - am + i;
    }
    flam = 1.0 / flam;

    const double lf = fact(angular_momentum_[0]);
    const double prefactor = sqrt(0.5 * (2*angular_momentum_[0]+1) / flam / pi__) / pow(2, angular_momentum_[0]) / lf;
    const int parity = (am - lx) % 2;
    double factor;
    if (angular_momentum_[1] > 0 && parity == 0) {
      factor = 1.0;
    } else if (angular_momentum_[1] == 0 && lx % 2 == 0) {
      factor = 1.0 / sqrt(2.0);
    } else if (angular_momentum_[1] < 0 && parity != 0) {
      factor = 1.0;
    } else {
      factor = 0.0;
    }

    double si = 0.0;
    Comb comb;
    const int jp = (angular_momentum_[0] - am) / 2;
    for (int i = j; i <= jp; ++i) {
      si += comb(angular_momentum_[0], i) * comb(i, j) * pow(-1.0, i)
          * ((2 * angular_momentum_[0] - 2 * i) >= 0 ? fact(2 * angular_momentum_[0] - 2 * i) : 0.0)
          / ((angular_momentum_[0] - am - 2 * i) >= 0 ? fact(angular_momentum_[0] - am - 2 * i) : 0.0);
    }
    double sk = 0.0;
    for (int k = 0; k <= j; ++k) {
      if (lx - 2*k >= 0 && lx - 2*k <= am) {
        sk += (k <= j ? comb(j, k) : 0.0) * comb(am, lx - 2*k) * pow(-1.0, static_cast<int>((am - lx + 2*k) / 2));
      }
    }
    return prefactor * si * sk * factor;
  }
}

vector<pair<double, int>> SphHarmonics::sph_to_USPs_expansion() {

  vector<pair<double, int>> usp;

  int cnt = 0;
  for (int lz = 0; lz <= angular_momentum_[0]; ++lz) {
    for (int ly = 0; ly <= angular_momentum_[0] - lz; ++ly) {
      ++cnt;
      const int lx = angular_momentum_[0] - lz - ly;
      const double coeff = sph_to_USP(lx, ly);
      pair<double, int> c_usp(coeff, cnt);
      usp.push_back(c_usp);
    }
  }

  return usp;

}

void SphHarmonics::sph_to_USPs_expansion_print() {

  cout << "** Real spherical harmonics to unitary sphere polynomials **" << endl;
  cout << "Angular momentum (lm) = (" << angular_momentum_[0] << ", "
                                      << angular_momentum_[1] << ")" << endl;
  cout << "Y_lm = sum_i c_i * usp_i" << endl;
  vector<pair<double, int>> usp = this->sph_to_USPs_expansion();
  for (auto& it : usp) {
    cout << "(" << setw(17) << setprecision(9) << it.first << ", " << it.second << ")" << endl;
  }

}

