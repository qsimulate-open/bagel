//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: kronecker.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/asd/dmrg/kronecker.h>
#include <src/util/math/matrix.h>
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

void bagel::kronecker_product_I_B(const double fac, const int Idim, const bool btrans, const int ndimB, const int mdimB, const double* B, const int ldB, double* C, const int ldC) {
  const int p = btrans ? mdimB : ndimB;
  const int q = btrans ? ndimB : mdimB;

  assert(ldC >= p*Idim);

  const int bstride = btrans ? ldB : 1;

  for (int ja = 0; ja < Idim; ++ja) {
    for (int jb = 0; jb < q; ++jb) {
      const double* bdata = btrans ? B + jb : B + ldB * jb;
      double* cdata = C + ja*p + (ja*q + jb)*ldC;
      daxpy_(p, fac, bdata, bstride, cdata, 1);
    }
  }
}

void bagel::kronecker_product_I_B(const double fac, const int Idim, const bool btrans, const Matrix& B, Matrix& C) {
  assert((btrans ? B.mdim() : B.ndim()) * Idim == C.ndim() && (btrans ? B.ndim() : B.mdim()) * Idim == C.mdim());
  kronecker_product_I_B(fac, Idim, btrans, B.ndim(), B.mdim(), B.data(), B.ndim(), C.data(), C.ndim());
}

void bagel::kronecker_product_A_I(const double fac, const bool atrans, const int ndimA, const int mdimA, const double* A, const int ldA, const int Idim, double* C, const int ldC) {
  const int n = atrans ? mdimA : ndimA;
  const int m = atrans ? ndimA : mdimA;

  assert(ldC >= n*Idim);

  const int astride = atrans ? ldA : 1;

  for (int ja = 0; ja < m; ++ja) {
    const double* adata = atrans ? A + ja : A + ldA * ja;
    for (int jb = 0; jb < Idim; ++jb) {
      double* cdata = C + jb + (jb + Idim*ja)*ldC;
      daxpy_(n, fac, adata, astride, cdata, Idim);
    }
  }
}

void bagel::kronecker_product_A_I(const double fac, const bool atrans, const Matrix& A, const int Idim,  Matrix& C) {
  assert((atrans ? A.mdim() : A.ndim()) * Idim == C.ndim() && (atrans ? A.ndim() : A.mdim()) * Idim == C.mdim());
  kronecker_product_A_I(fac, atrans, A.ndim(), A.mdim(), A.data(), A.ndim(), Idim, C.data(), C.ndim());
}
