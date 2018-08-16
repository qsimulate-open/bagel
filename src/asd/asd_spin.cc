//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_spin.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#include <src/asd/asd_spin.h>

using namespace bagel;
using namespace std;

shared_ptr<Matrix> ASDSpin::apply(const Matrix& o) const {
  return make_shared<Matrix>(*this * o);
}

void ASDSpin::filter(Matrix& o, const int desired_spin) const {
  assert(o.ndim() == this->mdim());

  for (int ispin = 0; ispin != max_spin_; ++ispin) {
    if (ispin == desired_spin) continue;

    Matrix S2 = *this * o;

    if (ispin == 0) {
      o = S2;
    }
    else {
      const double factor = -4.0/(static_cast<double>(ispin*(ispin+2)));
      o.ax_plus_y(factor, S2);
    }
  }
}
