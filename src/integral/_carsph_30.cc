//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_30.cc
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


void CarSphList::carsph_30(const int nloop, const double* source, double* target) {
  const double c3 = 3.872983346207417;
  const double c4 = 2.4494897427831779;
  const double c1 = 2.3717082451262845;
  const double c2 = 1.9364916731037085;
  const double c6 = 1.5;
  const double c0 = 0.79056941504209488;
  const double c5 = 0.61237243569579447;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 7, source += 10) {
    target[0] =  c0 * source[0] - c1 * source[2];
    target[1] =  c1 * source[1] - c0 * source[3];
    target[2] =  c2 * source[4] - c2 * source[6];
    target[3] =  c3 * source[5];
    target[4] =  c4 * source[7] - c5 * source[0] - c5 * source[2];
    target[5] =  c4 * source[8] - c5 * source[1] - c5 * source[3];
    target[6] =  source[9] - c6 * source[4] - c6 * source[6];
  }
}

void CCarSphList::carsph_30(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c3 = 3.872983346207417;
  const double c4 = 2.4494897427831779;
  const double c1 = 2.3717082451262845;
  const double c2 = 1.9364916731037085;
  const double c6 = 1.5;
  const double c0 = 0.79056941504209488;
  const double c5 = 0.61237243569579447;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 7, source += 10) {
    target[0] =  c0 * source[0] - c1 * source[2];
    target[1] =  c1 * source[1] - c0 * source[3];
    target[2] =  c2 * source[4] - c2 * source[6];
    target[3] =  c3 * source[5];
    target[4] =  c4 * source[7] - c5 * source[0] - c5 * source[2];
    target[5] =  c4 * source[8] - c5 * source[1] - c5 * source[3];
    target[6] =  source[9] - c6 * source[4] - c6 * source[6];
  }
}

