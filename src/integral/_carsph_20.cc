//
// BAGEL - Parallel electron correlation program.
// Filename: _carsph_20.cc
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

