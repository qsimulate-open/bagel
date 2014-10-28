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

  Matrix out(n*p, m*q);
  kronecker_product(1.0, atrans, A.ndim(), A.mdim(), A.data(), A.ndim(), btrans, B.ndim(), B.mdim(), B.data(), B.ndim(), out.data(), out.ndim());

  return out;
}

void bagel::kronecker_product(const double fac, const bool atrans, const Matrix& A, const bool btrans, const Matrix& B, Matrix& C) {
  assert((atrans ? A.mdim() : A.ndim())*(btrans ? B.mdim() : B.ndim()) == C.ndim() &&
         (atrans ? A.ndim() : A.mdim())*(btrans ? B.ndim() : B.mdim()) == C.mdim());
  kronecker_product(fac, atrans, A.ndim(), A.mdim(), A.data(), A.ndim(), btrans, B.ndim(), B.mdim(), B.data(), B.ndim(), C.data(), C.ndim());
}

void bagel::kronecker_product(const double fac, const bool atrans, const int ndimA, const int mdimA, const double* A, const int ldA,
                                         const bool btrans, const int ndimB, const int mdimB, const double* B, const int ldB,
                                                                                                    double* C, const int ldC)
{
  const int n = atrans ? mdimA : ndimA;
  const int m = atrans ? ndimA : mdimA;
  const int p = btrans ? mdimB : ndimB;
  const int q = btrans ? ndimB : mdimB;

  assert(ldC >= n * p);

  const int astride = atrans ? ldA : 1;
  const int bstride = btrans ? ldB : 1;

  for (int ja = 0; ja < m; ++ja) {
    const double* Adata = atrans ? A + ja : A + ldA * ja;
    for (int jb = 0; jb < q; ++jb) {
      const double* Bdata = btrans ? B + jb : B + ldB * jb;
      double* Cdata = C + ldC * (jb + ja*q);
      dger_(p, n, fac, Bdata, bstride, Adata, astride, Cdata, p);
    }
  }
}
