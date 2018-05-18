//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: algo.cc
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

#include <cstring>
#include <src/util/math/algo.h>
#include <bagel_config.h>
#ifdef HAVE_MKL_H
#include <src/util/mkl_ext.h>
#endif

using namespace std;

namespace bagel {

void dcsrmm_(const char *transa, const int m, const int n, const int k, const double alpha, const double* adata,
             const int* acols, const int* arind, const double* b, const int ldb, const double beta,
             double* c, const int ldc) {
#ifdef HAVE_MKL_H
  mkl_dcsrmm_(transa, m, n, k, alpha, adata, acols, arind, b, ldb, beta, c, ldc);
#else
  if (strcmp(transa, "N") != 0) throw logic_error("Only \"N\" case implemented for dcsrmm_");
  for (int j = 0; j < n; ++j) {
    double* target = c + j*ldc;
    const double* source = b + j*ldb;

    blas::scale_n(beta, target, m);

    for (int i = 0; i < m; ++i) {
      for (int rowdata = arind[i] - 1; rowdata < arind[i+1] - 1; ++rowdata) {
        target[i] += alpha * adata[rowdata] * source[acols[rowdata] - 1];
      }
    }
  }
#endif
}


void vdmul_(const int n, const double* a, const double* b, double* y) {
#ifdef HAVE_MKL_H
  ::vdmul_(&n, a, b, y);
#else
  for (size_t i = 0; i != n; ++i)
    y[i] = a[i] * b[i];
#endif
}

}
