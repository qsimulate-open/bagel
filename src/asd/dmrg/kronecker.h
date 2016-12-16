//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: kronecker.h
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

#ifndef __SRC_ASD_DMRG_KRONECKER_H
#define __SRC_ASD_DMRG_KRONECKER_H

#include <src/util/math/algo.h>
#include <src/util/math/matop.h>
#include <src/util/parallel/mpi_interface.h>

namespace bagel {

template <class T,
          class = typename std::enable_if<detail::is_mat<T>::value>::type
         >
T kronecker_product(const bool atrans, const T& A, const bool btrans, const T& B) {
  const int n = atrans ? A.mdim() : A.ndim();
  const int m = atrans ? A.ndim() : A.mdim();
  const int p = btrans ? B.mdim() : B.ndim();
  const int q = btrans ? B.ndim() : B.mdim();

  std::shared_ptr<const T> BB = btrans ? B.transpose() : nullptr;
  const T* Bptr = btrans ? BB.get() : &B;

  T out(n*p, m*q);
  for (int ja = 0; ja < m; ++ja) {
    for (int ia = 0; ia < n; ++ia) {
      const double aval = atrans ? A(ja, ia) : A(ia, ja);
      out.add_block(aval, ia*p, ja*q, p, q, *Bptr);
    }
  }

  return out;
}

/// Kronecker product between the identity and Matrix B
void kronecker_product_I_B(const double fac, const int Idim, const bool btrans, const int ndimB, const int mdimB, const double* B, const int ldB, double* C, const int ldC);
/// Kronecker product between the identity and Matrix B
void kronecker_product_I_B(const double fac, const int Idim, const bool btrans, const Matrix& B, Matrix& C);

/// Kronecker product between Matrix A and the identity
void kronecker_product_A_I(const double fac, const bool atrans, const int ndimA, const int mdimA, const double* A, const int ldA, const int Idim, double* C, const int ldC);
/// Kronecker product between Matrix A and the identity
void kronecker_product_A_I(const double fac, const bool atrans, const Matrix& A, const int Idim, Matrix& C);

template <>
Matrix kronecker_product(const bool atrans, const Matrix& A, const bool btrans, const Matrix& B);

void kronecker_product(const double fac, const bool atrans, const Matrix& A, const bool btrans, const Matrix& B, Matrix& C);

void kronecker_product(const double fac, const bool atrans, const int ndimA, const int mdimA, const double* A, const int ldA,
                                         const bool btrans, const int ndimB, const int mdimB, const double* B, const int ldB,
                                                                                                    double* C, const int ldC);

}

#endif
