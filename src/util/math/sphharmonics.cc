//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: sphharmonics.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <iostream>
#include <iomanip>
#include <src/util/math/sphharmonics.h>

using namespace std;
using namespace bagel;

const static Legendre legendre;
const static Factorial fact;

SphHarmonics::SphHarmonics(const array<int, 2> lm, const array<double, 3> c)
 : angular_momentum_(lm), centre_(c) {

  const double r = sqrt(centre_[0]*centre_[0] + centre_[1]*centre_[1] + centre_[2]*centre_[2]);
  theta_ = acos(centre_[2]/r);
  phi_ = atan2(centre_[1], centre_[0]);

}


SphHarmonics::SphHarmonics(const int l, const int m, const array<double, 3> c)
 : angular_momentum_{{l, m}}, centre_(c) {

  if (abs(m) > l)
    throw runtime_error ("SphHarmonics.zlm: |m| > l");

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


complex<double> SphHarmonics::ylm() const {

  const int l = angular_momentum_[0];
  const int m = angular_momentum_[1];

  const double cth = cos(theta_);
  const int am = abs(m);
  if (am > l)
    throw runtime_error ("SphHarmonics.ylm: |m| > l");

  const double plm = legendre.compute(l, am, cth);
  double fact = 1.0;
  for (int i = 1; i <= 2*am; ++i)
    fact *= l - am + i;

  const double coef = sqrt((2*l+1)*0.25 * fact / pi__);
  double real = coef * plm * cos(am*phi_);
  double imag = coef * plm * sin(am*phi_);

  if (m < 0) {
    real *= pow(-1, m);
    imag *= pow(-1, m+1);
  }

  return complex<double>(real, imag);

}


double SphHarmonics::zlm() const {

  const int l = angular_momentum_[0];
  const int m = angular_momentum_[1];

  const double cth = cos(theta_);
  const int am = abs(m);

  const double plm = legendre.compute(l, am, cth);

  double coef0 = 1.0;
  for (unsigned int i = l + am; i > (l - am); i--) coef0 *= i;
  const double coef = sqrt((2*l+1)/(4.0*pi__*coef0));

  double out = coef * plm;
  if (m > 0) {
    out *= sqrt(2.0) * pow(-1.0, m) * cos(am*phi_);
  } else if (m < 0) {
    out *= sqrt(2.0) * pow(-1.0, m) * sin(am*phi_);
  }

  return out;
}


double SphHarmonics::zlm(const int l, const int m) const {

  const int am = abs(m);
  if (am > l)
    throw runtime_error ("SphHarmonics.zlm: |m| > l");
  const double cth = cos(theta_);

  const double plm = legendre.compute(l, am, cth);

  double coef0 = 1.0;
  for (unsigned int i = l + am; i > (l - am); i--) coef0 *= i;
  const double coef = sqrt((2*l+1)/(4.0*pi__*coef0));

  double out = coef * plm;
  if (m > 0) {
    out *= sqrt(2.0) * pow(-1.0, m) * cos(am*phi_);
  } else if (m < 0) {
    out *= sqrt(2.0) * pow(-1.0, m) * sin(am*phi_);
  }

  return out;
}

void SphHarmonics::print() const {
  cout << "Computing Real Spherical Harmonics" << endl;
  cout << "Angular momentum (lm) = (" << angular_momentum_[0] << ", "
                                      << angular_momentum_[1] << ")" << endl;
  cout << "Centre = ";
  for (int i = 0; i != 3; ++i) cout << setw(17) << setprecision(9) << centre_[i] << ";   ";
  cout << endl;
}

double SphHarmonics::sph_to_USP(const int lx, const int ly) const {

  const int m = angular_momentum_[1];
  const int am = abs(m);
  const int j = (lx + ly - am) / 2;
  if ((lx + ly - am) % 2 != 0 || j < 0) {
    return 0.0;
  } else {
    const int l = angular_momentum_[0];
    const double lmamf = fact(l - am);
    const double lpamf = fact(l + am);
    const double lf = fact(l);
    const double prefactor = sqrt(0.5 * (2 * l + 1) * lmamf / lpamf / pi__) / pow(2, l) / lf;

    const int parity = (am - lx) % 2;
    double factor;

    if (m > 0 && parity == 0) {
      factor = 1.0;
    } else if (m == 0 && lx % 2 == 0) {
      factor = 1.0 / sqrt(2.0);
    } else if (m < 0 && parity != 0) {
      factor = 1.0;
    } else {
      factor = 0.0;
    }

    Comb comb;
    const int jp = (l - am) / 2;
    double si = 0.0;
    for (int i = j; i <= jp; ++i) {
      si += comb(l, i) * comb(i, j) * pow(-1.0, i) * ((2 * l - 2 * i) >= 0 ? fact(2 * l - 2 * i) : 0.0)
          / ((l - am - 2 * i) >= 0 ? fact(l - am - 2 * i) : 0.0);
    }
    double sk = 0.0;
    for (int k = 0; k <= j; ++k)
      if (lx - 2*k >= 0 && lx - 2*k <= am)
        sk += (k <= j ? comb(j, k) : 0.0) * comb(am, lx - 2*k) * pow(-1.0, static_cast<int>((am - lx + 2*k) / 2));

    return prefactor * si * sk * factor;

  }

}


vector<pair<double, int>> SphHarmonics::sph_to_USPs_expansion() const {

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


void SphHarmonics::sph_to_USPs_expansion_print() const {

  cout << "** Real spherical harmonics to unitary sphere polynomials **" << endl;
  cout << "Angular momentum (lm) = (" << angular_momentum_[0] << ", "
                                      << angular_momentum_[1] << ")" << endl;
  cout << "Y_lm = sum_i c_i * usp_i" << endl;
  vector<pair<double, int>> usp = this->sph_to_USPs_expansion();
  for (auto& it : usp) {
    cout << "(" << setw(17) << setprecision(9) << it.first << ", " << it.second << ")" << endl;
  }

}

