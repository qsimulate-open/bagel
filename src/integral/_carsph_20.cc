//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_20.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/integral/carsphlist.h>
#include <algorithm>

using namespace std;
using namespace bagel;


void CarSphList::carsph_20(const int nloop, const double* source, double* target) {
  const double c1 = 1.7320508075688772;
  const double c0 = 0.8660254037844386;
  const double c2 = 0.5;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 5, source += 6) {
    target[0] =  c0 * source[0] - c0 * source[2];
    target[1] =  c1 * source[1];
    target[2] =  c1 * source[3];
    target[3] =  c1 * source[4];
    target[4] =  source[5] - c2 * source[0] - c2 * source[2];
  }
}

void CCarSphList::carsph_20(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c1 = 1.7320508075688772;
  const double c0 = 0.8660254037844386;
  const double c2 = 0.5;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 5, source += 6) {
    target[0] =  c0 * source[0] - c0 * source[2];
    target[1] =  c1 * source[1];
    target[2] =  c1 * source[3];
    target[3] =  c1 * source[4];
    target[4] =  source[5] - c2 * source[0] - c2 * source[2];
  }
}

