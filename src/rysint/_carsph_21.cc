//
// Newint - Parallel electron correlation program.
// Filename: _carsph_21.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


// editied by hand a little

#include <src/rysint/carsphlist.h>
#include <src/util/f77.h>


void CarSphList::carsph_21(const int nloop, const double* source, double* target) {
  const double c1 = 1.7320508075688772;
  const double c0 = 0.8660254037844386;
  const double c2 = 0.5;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 15, source += 18) {
    for (int i = 0; i != 3; ++i)
      target[  i] = c0 * (source[i] - source[6+i]);
    for (int i = 0; i != 3; ++i)
      target[3+i] = c1 * source[3+i];
    for (int i = 0; i != 7; ++i)
      target[6+i] = c1 * source[9+i];
    for (int i = 0; i != 3; ++i)
      target[12+i] = source[15+i] - c2 * (source[i] + source[6+i]);
  }
}

