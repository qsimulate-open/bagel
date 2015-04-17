//
// BAGEL - Parallel electron correlation program.
// Filename: sphmultipole.cc
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


#include <iomanip>
#include <src/periodic/sphmultipole.h>

using namespace std;
using namespace bagel;

const static Legendre plm;
const static Factorial f;

SphMultipole::SphMultipole(const array<double, 3> c, const bool do_complex, const int l) : centre_(c), do_complex_(do_complex), lmax_(l) {

  num_multipoles_ = (lmax_ + 1) * (lmax_ + 1);
  if (do_complex_) {
    compute_multipoles();
  } else {
    compute_real_multipoles();
  }
}


void SphMultipole::compute_multipoles() {

  const double r = sqrt(centre_[0]*centre_[0] + centre_[1]*centre_[1] + centre_[2]*centre_[2]);
  const double ctheta = (r > numerical_zero__) ? centre_[2]/r : 0.0;
  const double phi = atan2(centre_[1], centre_[0]);

  multipole_.resize(num_multipoles_);

  int count = 1;
  multipole_[0] = 1.0;
  for (int l = 1; l <= lmax_; ++l) {
    for (int mm = 0; mm <= 2 * l; ++mm, ++count) {
      const int m = mm - l;
      const int am = abs(m);

      double coeff = pow(r, l) * plm.compute(l, am, ctheta);
      double ft = 1.0;
      for (int i = 1; i <= l + am; ++i) {
        coeff /= ft;
        ++ft;
      }

      const double real = (m >=0) ? (coeff * cos(am * phi)) : (-1.0 * coeff * cos(am * phi));
      const double imag = coeff * sin(am * phi);
      multipole_[count] = complex<double>(real, imag);

    }
  }
}


complex<double> SphMultipole::multipole(const int l, const int m) const {
  assert (l <= lmax_ && abs(m) <= l); return multipole_[l * l + l + m];
}



vector<std::complex<double>> SphMultipole::multipoles(const int l) {
  assert (l <= lmax_);
  vector<std::complex<double>> out(2 * l + 1);
  const int i0 = (l + 1) * (l + 1);
  for (int i = 0; i != 2 * l + 1; ++i) out[i] = multipole_[i + i0];

  return out;
}


void SphMultipole::compute_real_multipoles() {

  const double r = sqrt(centre_[0]*centre_[0] + centre_[1]*centre_[1] + centre_[2]*centre_[2]);
  const double ctheta = centre_[2]/r;
  const double phi = atan2(centre_[1], centre_[0]);

  real_multipole_.resize(num_multipoles_);

  real_multipole_[0] = 1.0;
  const double sqrttwo = sqrt(2.0);
  for (int l = 1; l <= lmax_; ++l) {
    real_multipole_[l * l] = 0.5 * pow(r, l) * plm.compute(l, 0, ctheta);
    for (int m = 1; m <= l; ++m) {

      const double coeff = pow(r, l) * plm.compute(l, m, ctheta) * sqrt(f(l - m) / f(l + m));

      real_multipole_[l*l+l+m] = sqrttwo * coeff * cos(m * phi);
      real_multipole_[l*l+l-m] = sqrttwo * coeff * sin(m * phi);
    }
  }
}


void SphMultipole::print_multipoles() const {

  cout << "LMAX = " << lmax_ << endl;
  int cnt = 0;
  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m, ++cnt)
      if (do_complex_) {
        cout << setprecision(6) << multipole_[cnt] << "   ";
      } else {
        cout << setprecision(6) << real_multipole_[cnt] << "   ";
      }
    cout << endl;
  }
}
