//
// BAGEL - Parallel electron correlation program.
// Filename: meh_spin.cc
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#include <src/meh/meh_spin.h>

using namespace bagel;
using namespace std;

shared_ptr<Matrix> MEHSpin::apply(const Matrix& o) const {
  return make_shared<Matrix>(*this * o);
}

void MEHSpin::filter(Matrix& o, const int desired_spin) const {
  const int n = o.ndim();
  assert( n == this->mdim() );
  const int m = o.mdim();

  for (int ispin = 0; ispin != max_spin_; ++ispin) {
    if (ispin == desired_spin) continue;

    Matrix S2 = *this * o;

    const double factor = -4.0/(static_cast<double>(ispin*(ispin+2)));
    o.ax_plus_y(factor, S2);
  }
}
