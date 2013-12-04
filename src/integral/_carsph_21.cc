//
// BAGEL - Parallel electron correlation program.
// Filename: _carsph_21.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/integral/carsphlist.h>
#include <algorithm>

using namespace std;
using namespace bagel;


void CarSphList::carsph_21(const int nloop, const double* source, double* target) {
  const double c1 = 1.7320508075688772;
  const double c0 = 0.8660254037844386;
  const double c2 = 0.5;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 15, source += 18) {
    target[0] =  c0 * source[0] - c0 * source[6];
    target[1] =  c0 * source[1] - c0 * source[7];
    target[2] =  c0 * source[2] - c0 * source[8];
    target[3] =  c1 * source[3];
    target[4] =  c1 * source[4];
    target[5] =  c1 * source[5];
    target[6] =  c1 * source[9];
    target[7] =  c1 * source[10];
    target[8] =  c1 * source[11];
    target[9] =  c1 * source[12];
    target[10] =  c1 * source[13];
    target[11] =  c1 * source[14];
    target[12] =  source[15] - c2 * source[0] - c2 * source[6];
    target[13] =  source[16] - c2 * source[1] - c2 * source[7];
    target[14] =  source[17] - c2 * source[2] - c2 * source[8];
  }
}

void CCarSphList::carsph_21(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c1 = 1.7320508075688772;
  const double c0 = 0.8660254037844386;
  const double c2 = 0.5;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 15, source += 18) {
    target[0] =  c0 * source[0] - c0 * source[6];
    target[1] =  c0 * source[1] - c0 * source[7];
    target[2] =  c0 * source[2] - c0 * source[8];
    target[3] =  c1 * source[3];
    target[4] =  c1 * source[4];
    target[5] =  c1 * source[5];
    target[6] =  c1 * source[9];
    target[7] =  c1 * source[10];
    target[8] =  c1 * source[11];
    target[9] =  c1 * source[12];
    target[10] =  c1 * source[13];
    target[11] =  c1 * source[14];
    target[12] =  source[15] - c2 * source[0] - c2 * source[6];
    target[13] =  source[16] - c2 * source[1] - c2 * source[7];
    target[14] =  source[17] - c2 * source[2] - c2 * source[8];
  }
}

