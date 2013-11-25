//
// BAGEL - Parallel electron correlation program.
// Filename: sparsealgebra.cc
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

#include <src/math/algo.h>
#include <bagel_config.h>
#ifdef HAVE_MKL_H
#include <src/util/mkl_sparse.h>
#endif

using namespace std;

namespace bagel {
void dcsrmm_(const char *transa, const int m, const int n, const int k, const double alpha, const double* adata,
                                 const int* acols, const int* arind, const double* b, const int ldb, const double beta,
                                 double* c, const int ldc) {
#ifdef HAVE_MKL_H
  mkl_dcsrmm_(transa, m, n, k, alpha, adata, acols, arind, b, ldb, beta, c, ldc);
#else
  if (transa != "N") throw logic_error("Only \"N\" case implemented for dcsrmm_");
  for (int j = 0; j < n; ++j) {
    double* target = c + j*ldc;
    const double* source = b + j*ldb;

    for (int i = 0; i < m; ++i) {
      for (int rowdata = arind[i] - 1; rowdata < arind[i+1] - 1; ++rowdata) {
        target[i] += adata[rowdata] * source[acols[rowdata] - 1];
      }
    }
  }
#endif
}
}
