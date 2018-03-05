//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mkl_ext.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __SRC_UTIL_MKL_EXT_H
#define __SRC_UTIL_MKL_EXT_H

#include <bagel_config.h>

#ifdef HAVE_MKL_H

#include <memory>

extern "C" {
  // mkl_sparse routines
  void mkl_dcsrmm_(const char* transa, const int* m, const int* n, const int* k, const double* alpha, const char* matdescra,
                   const double* val, const int* indx, const int* pntrb, const int* pntre,
                   const double* b, const int* ldb, const double* beta, double* c, const int* ldc);

  // MKL Hadamard product
  void vdmul_(const int* n, const double* a, const double* b, double* y); 
}

// All arguments passed in
static void mkl_dcsrmm_(const char* transa, const int m, const int n, const int k, const double alpha, const char* matdescra,
                        const double* val, const int* indx, const int* pntrb, const int* pntre,
                        const double* b, const int ldb, const double beta, double* c, const int ldc) {
  mkl_dcsrmm_(transa, &m, &n, &k, &alpha, matdescra, val, indx, pntrb, pntre, b, &ldb, &beta, c, &ldc);
}
// Special case
static void mkl_dcsrmm_(const char* transa, const int m, const int n, const int k, const double alpha,
                        const double* val, const int* cols, const int* rind,
                        const double* b, const int ldb, const double beta, double* c, const int ldc) {
  std::unique_ptr<char[]> matdescra(new char[6]);
  matdescra[0] = 'G'; // General
  matdescra[3] = 'F'; // Fortran 1-based indexing (this is required to use a column-major dense matrix)
  mkl_dcsrmm_(transa, m, n, k, alpha, matdescra.get(), val, cols, rind, rind+1, b, ldb, beta, c, ldc);
}

#endif

#endif
