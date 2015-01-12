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

SphMultipole::SphMultipole(const array<double, 3> c, const int l) : centre_(c), lmax_(l) {

  num_multipoles_ = (lmax_ + 1) * (lmax_ + 1);
  compute_multipoles();
}


void SphMultipole::compute_multipoles() {

  const double r = sqrt(centre_[0]*centre_[0] + centre_[1]*centre_[1] + centre_[2]*centre_[2]);
  const double ctheta = centre_[2]/r;
  const double phi = atan2(centre_[1], centre_[0]);

  multipole_.resize(num_multipoles_);

  int count = 0;
  for (int l = 0; l <= lmax_; ++l) {
    for (int mm = 0; mm <= 2 * l; ++mm) {
      const int m = mm - l;
      const double coeff = pow(r, l) * plm.compute(l, abs(m), ctheta) / f(l + abs(m));

      double real = coeff * cos(-m * phi);
      double imag = coeff * sin(-m * phi);
      multipole_[count] = complex<double>(real, imag);

      ++count;
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


void SphMultipole::print_multipoles() const {

  cout << "LMAX = " << lmax_ << endl;
  int cnt = 0;
  for (int l = 0; l <= lmax_; ++l) {
    for (int m = 0; m <= 2 * l; ++m, ++cnt)
      cout << setprecision(9) << multipole_[cnt] << "   ";
    cout << endl;
  }
}

