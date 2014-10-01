//
// BAGEL - Parallel electron correlation program.
// Filename: kronecker.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/asd_dmrg/kronecker.h>
#include <src/math/matrix.h>
#include <src/util/taskqueue.h>

using namespace bagel;
using namespace std;

template <>
Matrix bagel::kronecker_product(const bool atrans, const Matrix& A, const bool btrans, const Matrix& B) {
  const int n = atrans ? A.mdim() : A.ndim();
  const int m = atrans ? A.ndim() : A.mdim();
  const int p = btrans ? B.mdim() : B.ndim();
  const int q = btrans ? B.ndim() : B.mdim();

  std::shared_ptr<const Matrix> BB = btrans ? B.transpose() : nullptr;
  const Matrix* Bptr = btrans ? BB.get() : &B;

  Matrix out(n*p, m*q);
  const int astride = atrans ? A.ndim() : 1;

  for (int jb = 0; jb < q; ++jb) {
    for (int ja = 0; ja < m; ++ja) {
      const double* Adata = atrans ? A.element_ptr(ja, 0) : A.element_ptr(0, ja);
      dger_(p, n, 1.0, Bptr->element_ptr(0, jb), 1, Adata, astride, out.element_ptr(0, jb + ja*q), p);
    }
  }

  return out;
}
