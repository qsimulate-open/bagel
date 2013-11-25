//
// BAGEL - Parallel electron correlation program.
// Filename: mkl_sparse.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef __SRC_UTIL_MKL_SPARSE_H
#define __SRC_UTIL_MKL_SPARSE_H

#include <bagel_config.h>

#ifdef HAVE_MKL_H

//#include <mkl.h>
#include <memory>

extern "C" {
  // mkl_sparse routines
  //void mkl_dcsrmm (char *transa, MKL_INT *m, MKL_INT *n, MKL_INT *k, double *alpha, char *matdescra, double *val, MKL_INT *indx, MKL_INT *pntrb, MKL_INT *pntre, double *b, MKL_INT *ldb, double *beta, double *c, MKL_INT *ldc);
  void mkl_dcsrmm_(const char *transa, const int *m, const int *n, const int *k, const double *alpha, const char *matdescra, const double *val, const int *indx, const int *pntrb, const int *pntre,
                   const double *b, const int *ldb, const double *beta, double *c, const int *ldc);
}

// All arguments passed in
static void mkl_dcsrmm_(const char *transa, int m, int n, int k, double alpha, const char *matdescra, const double* val, const int* indx, const int* pntrb, const int* pntre, const double* b, int ldb, double beta, double* c, int ldc)
    { mkl_dcsrmm_(transa, &m, &n, &k, &alpha, matdescra, val, indx, pntrb, pntre, b, &ldb, &beta, c, &ldc); }

// Special case
static void mkl_dcsrmm_(const char *transa, int m, int n, int k, double alpha, const double* val, const int* cols, const int* rind, const double* b, int ldb, double beta, double* c, int ldc) {
  std::unique_ptr<char[]> matdescra(new char[6]);
  matdescra[0] = 'G'; // General
  matdescra[3] = 'F'; // Fortran 1-based indexing (this is required to use a column-major dense matrix)
  mkl_dcsrmm_(transa, m, n, k, alpha, matdescra.get(), val, cols, rind, rind+1, b, ldb, beta, c, ldc);
}

#endif

#endif
