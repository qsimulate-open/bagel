//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dvop.cc
// Copyright (C) 2017 Nils Strand
//
// Author: Nils Strand <nilsstrand2022@u.northwestern.edu>
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


#include <src/wfn/diagvec.h>

namespace bagel {

  Matrix operator*(const DiagVec &v, const Matrix &m) {
    assert(m.ndim() == v.size());
    Matrix out(v.size(), v.size());
    for (int i = 0; i < m.mdim(); i++) {
      for (int j = 0; j < m.ndim(); j++) {
        out(j, i) = m(j, i) * v(j);
      }
    }
    return out;
  }

  Matrix operator*(const Matrix &m, const DiagVec &v) {
    assert(m.mdim() == v.size());
    Matrix out = m;
    for (int i = 0; i != m.mdim(); ++i)
      blas::scale_n(v(i), out.element_ptr(0, i), m.ndim());
    return out;
  }

}
