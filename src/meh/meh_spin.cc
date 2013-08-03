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
// the Free Software Foundation; either version 2, or (at your option)
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

  Matrix tmp(n, m);
  for (int ispin = 0; ispin != max_spin_; ++ispin) {
    if (ispin == desired_spin) continue;
    const double kk1 = 0.25 * static_cast<double>( ispin*(ispin+2) );

    copy_n(o.data(), n * m, tmp.data());
    // o = (SpinMatrix - eye*kk1) * o
    for (int j = 0; j < m; ++j) {
      double* target = o.element_ptr(0,j);
      const double* source = tmp.element_ptr(0,j);
      for (int i = 0; i < n; ++i) {
        target[i] = (diagonal(i) - kk1) * source[i];
      }
      for (auto& offdiag : offdiagonal_) {
        target[offdiag.i] += offdiag.value * source[offdiag.j];
        target[offdiag.j] += offdiag.value * source[offdiag.i];
      }
    }
  }
}
