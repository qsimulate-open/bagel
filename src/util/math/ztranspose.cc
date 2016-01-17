//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ztranspose.cc
// Copyright (C) 2008 Toru Shiozaki
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

#include <bagel_config.h>
#include <src/util/math/algo.h>

using namespace std;

#ifdef HAVE_MKL_H
extern "C" {
  void mkl_zomatcopy_(const char*, const char*, const int*, const int*, const complex<double>*, const complex<double>*, const int*, complex<double>*, const int*);
}
#endif

namespace bagel {
namespace blas {

template<>
void transpose(const complex<double>* h, const int m, const int n, complex<double>* vec, const double fac) {
  transpose(h, m, n, vec, static_cast<complex<double>>(fac));
}

template<>
void transpose(const complex<double>* h, const int m, const int n, complex<double>* vec, const complex<double> fac) {
#ifdef HAVE_MKL_H
  mkl_zomatcopy_("c", "t", &m, &n, &fac, h, &m, vec, &n);
#else
  const int mresidual = m % 10;
  const int nresidual = n % 10;
  const int mlim = m - mresidual;
  const int nlim = n - nresidual;
  for (int i = 0; i < nlim; i += 10) {
    for (int j = 0; j < mlim; j += 10) {
      vec[i   + n*(j  )] = fac*h[j   + m*(i  )];
      vec[i+1 + n*(j  )] = fac*h[j   + m*(i+1)];
      vec[i+2 + n*(j  )] = fac*h[j   + m*(i+2)];
      vec[i+3 + n*(j  )] = fac*h[j   + m*(i+3)];
      vec[i+4 + n*(j  )] = fac*h[j   + m*(i+4)];
      vec[i+5 + n*(j  )] = fac*h[j   + m*(i+5)];
      vec[i+6 + n*(j  )] = fac*h[j   + m*(i+6)];
      vec[i+7 + n*(j  )] = fac*h[j   + m*(i+7)];
      vec[i+8 + n*(j  )] = fac*h[j   + m*(i+8)];
      vec[i+9 + n*(j  )] = fac*h[j   + m*(i+9)];
      vec[i   + n*(j+1)] = fac*h[j+1 + m*(i  )];
      vec[i+1 + n*(j+1)] = fac*h[j+1 + m*(i+1)];
      vec[i+2 + n*(j+1)] = fac*h[j+1 + m*(i+2)];
      vec[i+3 + n*(j+1)] = fac*h[j+1 + m*(i+3)];
      vec[i+4 + n*(j+1)] = fac*h[j+1 + m*(i+4)];
      vec[i+5 + n*(j+1)] = fac*h[j+1 + m*(i+5)];
      vec[i+6 + n*(j+1)] = fac*h[j+1 + m*(i+6)];
      vec[i+7 + n*(j+1)] = fac*h[j+1 + m*(i+7)];
      vec[i+8 + n*(j+1)] = fac*h[j+1 + m*(i+8)];
      vec[i+9 + n*(j+1)] = fac*h[j+1 + m*(i+9)];
      vec[i   + n*(j+2)] = fac*h[j+2 + m*(i  )];
      vec[i+1 + n*(j+2)] = fac*h[j+2 + m*(i+1)];
      vec[i+2 + n*(j+2)] = fac*h[j+2 + m*(i+2)];
      vec[i+3 + n*(j+2)] = fac*h[j+2 + m*(i+3)];
      vec[i+4 + n*(j+2)] = fac*h[j+2 + m*(i+4)];
      vec[i+5 + n*(j+2)] = fac*h[j+2 + m*(i+5)];
      vec[i+6 + n*(j+2)] = fac*h[j+2 + m*(i+6)];
      vec[i+7 + n*(j+2)] = fac*h[j+2 + m*(i+7)];
      vec[i+8 + n*(j+2)] = fac*h[j+2 + m*(i+8)];
      vec[i+9 + n*(j+2)] = fac*h[j+2 + m*(i+9)];
      vec[i   + n*(j+3)] = fac*h[j+3 + m*(i  )];
      vec[i+1 + n*(j+3)] = fac*h[j+3 + m*(i+1)];
      vec[i+2 + n*(j+3)] = fac*h[j+3 + m*(i+2)];
      vec[i+3 + n*(j+3)] = fac*h[j+3 + m*(i+3)];
      vec[i+4 + n*(j+3)] = fac*h[j+3 + m*(i+4)];
      vec[i+5 + n*(j+3)] = fac*h[j+3 + m*(i+5)];
      vec[i+6 + n*(j+3)] = fac*h[j+3 + m*(i+6)];
      vec[i+7 + n*(j+3)] = fac*h[j+3 + m*(i+7)];
      vec[i+8 + n*(j+3)] = fac*h[j+3 + m*(i+8)];
      vec[i+9 + n*(j+3)] = fac*h[j+3 + m*(i+9)];
      vec[i   + n*(j+4)] = fac*h[j+4 + m*(i  )];
      vec[i+1 + n*(j+4)] = fac*h[j+4 + m*(i+1)];
      vec[i+2 + n*(j+4)] = fac*h[j+4 + m*(i+2)];
      vec[i+3 + n*(j+4)] = fac*h[j+4 + m*(i+3)];
      vec[i+4 + n*(j+4)] = fac*h[j+4 + m*(i+4)];
      vec[i+5 + n*(j+4)] = fac*h[j+4 + m*(i+5)];
      vec[i+6 + n*(j+4)] = fac*h[j+4 + m*(i+6)];
      vec[i+7 + n*(j+4)] = fac*h[j+4 + m*(i+7)];
      vec[i+8 + n*(j+4)] = fac*h[j+4 + m*(i+8)];
      vec[i+9 + n*(j+4)] = fac*h[j+4 + m*(i+9)];
      vec[i   + n*(j+5)] = fac*h[j+5 + m*(i  )];
      vec[i+1 + n*(j+5)] = fac*h[j+5 + m*(i+1)];
      vec[i+2 + n*(j+5)] = fac*h[j+5 + m*(i+2)];
      vec[i+3 + n*(j+5)] = fac*h[j+5 + m*(i+3)];
      vec[i+4 + n*(j+5)] = fac*h[j+5 + m*(i+4)];
      vec[i+5 + n*(j+5)] = fac*h[j+5 + m*(i+5)];
      vec[i+6 + n*(j+5)] = fac*h[j+5 + m*(i+6)];
      vec[i+7 + n*(j+5)] = fac*h[j+5 + m*(i+7)];
      vec[i+8 + n*(j+5)] = fac*h[j+5 + m*(i+8)];
      vec[i+9 + n*(j+5)] = fac*h[j+5 + m*(i+9)];
      vec[i   + n*(j+6)] = fac*h[j+6 + m*(i  )];
      vec[i+1 + n*(j+6)] = fac*h[j+6 + m*(i+1)];
      vec[i+2 + n*(j+6)] = fac*h[j+6 + m*(i+2)];
      vec[i+3 + n*(j+6)] = fac*h[j+6 + m*(i+3)];
      vec[i+4 + n*(j+6)] = fac*h[j+6 + m*(i+4)];
      vec[i+5 + n*(j+6)] = fac*h[j+6 + m*(i+5)];
      vec[i+6 + n*(j+6)] = fac*h[j+6 + m*(i+6)];
      vec[i+7 + n*(j+6)] = fac*h[j+6 + m*(i+7)];
      vec[i+8 + n*(j+6)] = fac*h[j+6 + m*(i+8)];
      vec[i+9 + n*(j+6)] = fac*h[j+6 + m*(i+9)];
      vec[i   + n*(j+7)] = fac*h[j+7 + m*(i  )];
      vec[i+1 + n*(j+7)] = fac*h[j+7 + m*(i+1)];
      vec[i+2 + n*(j+7)] = fac*h[j+7 + m*(i+2)];
      vec[i+3 + n*(j+7)] = fac*h[j+7 + m*(i+3)];
      vec[i+4 + n*(j+7)] = fac*h[j+7 + m*(i+4)];
      vec[i+5 + n*(j+7)] = fac*h[j+7 + m*(i+5)];
      vec[i+6 + n*(j+7)] = fac*h[j+7 + m*(i+6)];
      vec[i+7 + n*(j+7)] = fac*h[j+7 + m*(i+7)];
      vec[i+8 + n*(j+7)] = fac*h[j+7 + m*(i+8)];
      vec[i+9 + n*(j+7)] = fac*h[j+7 + m*(i+9)];
      vec[i   + n*(j+8)] = fac*h[j+8 + m*(i  )];
      vec[i+1 + n*(j+8)] = fac*h[j+8 + m*(i+1)];
      vec[i+2 + n*(j+8)] = fac*h[j+8 + m*(i+2)];
      vec[i+3 + n*(j+8)] = fac*h[j+8 + m*(i+3)];
      vec[i+4 + n*(j+8)] = fac*h[j+8 + m*(i+4)];
      vec[i+5 + n*(j+8)] = fac*h[j+8 + m*(i+5)];
      vec[i+6 + n*(j+8)] = fac*h[j+8 + m*(i+6)];
      vec[i+7 + n*(j+8)] = fac*h[j+8 + m*(i+7)];
      vec[i+8 + n*(j+8)] = fac*h[j+8 + m*(i+8)];
      vec[i+9 + n*(j+8)] = fac*h[j+8 + m*(i+9)];
      vec[i   + n*(j+9)] = fac*h[j+9 + m*(i  )];
      vec[i+1 + n*(j+9)] = fac*h[j+9 + m*(i+1)];
      vec[i+2 + n*(j+9)] = fac*h[j+9 + m*(i+2)];
      vec[i+3 + n*(j+9)] = fac*h[j+9 + m*(i+3)];
      vec[i+4 + n*(j+9)] = fac*h[j+9 + m*(i+4)];
      vec[i+5 + n*(j+9)] = fac*h[j+9 + m*(i+5)];
      vec[i+6 + n*(j+9)] = fac*h[j+9 + m*(i+6)];
      vec[i+7 + n*(j+9)] = fac*h[j+9 + m*(i+7)];
      vec[i+8 + n*(j+9)] = fac*h[j+9 + m*(i+8)];
      vec[i+9 + n*(j+9)] = fac*h[j+9 + m*(i+9)];
    }
    for (int j = mlim; j != m; ++j) {
      vec[i   + n*j] = fac*h[j  + m*(i  )];
      vec[i+1 + n*j] = fac*h[j  + m*(i+1)];
      vec[i+2 + n*j] = fac*h[j  + m*(i+2)];
      vec[i+3 + n*j] = fac*h[j  + m*(i+3)];
      vec[i+4 + n*j] = fac*h[j  + m*(i+4)];
      vec[i+5 + n*j] = fac*h[j  + m*(i+5)];
      vec[i+6 + n*j] = fac*h[j  + m*(i+6)];
      vec[i+7 + n*j] = fac*h[j  + m*(i+7)];
      vec[i+8 + n*j] = fac*h[j  + m*(i+8)];
      vec[i+9 + n*j] = fac*h[j  + m*(i+9)];
    }
  }
  for (int i = nlim; i != n; ++i) {
    for (int j = 0; j < mlim; j += 10) {
      vec[i   + n*(j  )] = fac*h[j   + m*(i  )];
      vec[i   + n*(j+1)] = fac*h[j+1 + m*(i  )];
      vec[i   + n*(j+2)] = fac*h[j+2 + m*(i  )];
      vec[i   + n*(j+3)] = fac*h[j+3 + m*(i  )];
      vec[i   + n*(j+4)] = fac*h[j+4 + m*(i  )];
      vec[i   + n*(j+5)] = fac*h[j+5 + m*(i  )];
      vec[i   + n*(j+6)] = fac*h[j+6 + m*(i  )];
      vec[i   + n*(j+7)] = fac*h[j+7 + m*(i  )];
      vec[i   + n*(j+8)] = fac*h[j+8 + m*(i  )];
      vec[i   + n*(j+9)] = fac*h[j+9 + m*(i  )];
    }
    for (int j = mlim; j !=  m; ++j) {
      vec[i   + n*(j  )] = fac*h[j   + m*(i  )];
    }
  }
#endif
}


template<>
void transpose_add(const complex<double>* h, const int m, const int n, complex<double>* vec, const double fac) {
  transpose_add(h, m, n, vec, static_cast<complex<double>>(fac));
}

template<>
void transpose_add(const complex<double>* h, const int m, const int n, complex<double>* vec, const complex<double> fac) {
  const int mresidual = m % 10;
  const int nresidual = n % 10;
  const int mlim = m - mresidual;
  const int nlim = n - nresidual;
  for (int i = 0; i < nlim; i += 10) {
    for (int j = 0; j < mlim; j += 10) {
      vec[i   + n*(j  )] += fac*h[j   + m*(i  )];
      vec[i+1 + n*(j  )] += fac*h[j   + m*(i+1)];
      vec[i+2 + n*(j  )] += fac*h[j   + m*(i+2)];
      vec[i+3 + n*(j  )] += fac*h[j   + m*(i+3)];
      vec[i+4 + n*(j  )] += fac*h[j   + m*(i+4)];
      vec[i+5 + n*(j  )] += fac*h[j   + m*(i+5)];
      vec[i+6 + n*(j  )] += fac*h[j   + m*(i+6)];
      vec[i+7 + n*(j  )] += fac*h[j   + m*(i+7)];
      vec[i+8 + n*(j  )] += fac*h[j   + m*(i+8)];
      vec[i+9 + n*(j  )] += fac*h[j   + m*(i+9)];
      vec[i   + n*(j+1)] += fac*h[j+1 + m*(i  )];
      vec[i+1 + n*(j+1)] += fac*h[j+1 + m*(i+1)];
      vec[i+2 + n*(j+1)] += fac*h[j+1 + m*(i+2)];
      vec[i+3 + n*(j+1)] += fac*h[j+1 + m*(i+3)];
      vec[i+4 + n*(j+1)] += fac*h[j+1 + m*(i+4)];
      vec[i+5 + n*(j+1)] += fac*h[j+1 + m*(i+5)];
      vec[i+6 + n*(j+1)] += fac*h[j+1 + m*(i+6)];
      vec[i+7 + n*(j+1)] += fac*h[j+1 + m*(i+7)];
      vec[i+8 + n*(j+1)] += fac*h[j+1 + m*(i+8)];
      vec[i+9 + n*(j+1)] += fac*h[j+1 + m*(i+9)];
      vec[i   + n*(j+2)] += fac*h[j+2 + m*(i  )];
      vec[i+1 + n*(j+2)] += fac*h[j+2 + m*(i+1)];
      vec[i+2 + n*(j+2)] += fac*h[j+2 + m*(i+2)];
      vec[i+3 + n*(j+2)] += fac*h[j+2 + m*(i+3)];
      vec[i+4 + n*(j+2)] += fac*h[j+2 + m*(i+4)];
      vec[i+5 + n*(j+2)] += fac*h[j+2 + m*(i+5)];
      vec[i+6 + n*(j+2)] += fac*h[j+2 + m*(i+6)];
      vec[i+7 + n*(j+2)] += fac*h[j+2 + m*(i+7)];
      vec[i+8 + n*(j+2)] += fac*h[j+2 + m*(i+8)];
      vec[i+9 + n*(j+2)] += fac*h[j+2 + m*(i+9)];
      vec[i   + n*(j+3)] += fac*h[j+3 + m*(i  )];
      vec[i+1 + n*(j+3)] += fac*h[j+3 + m*(i+1)];
      vec[i+2 + n*(j+3)] += fac*h[j+3 + m*(i+2)];
      vec[i+3 + n*(j+3)] += fac*h[j+3 + m*(i+3)];
      vec[i+4 + n*(j+3)] += fac*h[j+3 + m*(i+4)];
      vec[i+5 + n*(j+3)] += fac*h[j+3 + m*(i+5)];
      vec[i+6 + n*(j+3)] += fac*h[j+3 + m*(i+6)];
      vec[i+7 + n*(j+3)] += fac*h[j+3 + m*(i+7)];
      vec[i+8 + n*(j+3)] += fac*h[j+3 + m*(i+8)];
      vec[i+9 + n*(j+3)] += fac*h[j+3 + m*(i+9)];
      vec[i   + n*(j+4)] += fac*h[j+4 + m*(i  )];
      vec[i+1 + n*(j+4)] += fac*h[j+4 + m*(i+1)];
      vec[i+2 + n*(j+4)] += fac*h[j+4 + m*(i+2)];
      vec[i+3 + n*(j+4)] += fac*h[j+4 + m*(i+3)];
      vec[i+4 + n*(j+4)] += fac*h[j+4 + m*(i+4)];
      vec[i+5 + n*(j+4)] += fac*h[j+4 + m*(i+5)];
      vec[i+6 + n*(j+4)] += fac*h[j+4 + m*(i+6)];
      vec[i+7 + n*(j+4)] += fac*h[j+4 + m*(i+7)];
      vec[i+8 + n*(j+4)] += fac*h[j+4 + m*(i+8)];
      vec[i+9 + n*(j+4)] += fac*h[j+4 + m*(i+9)];
      vec[i   + n*(j+5)] += fac*h[j+5 + m*(i  )];
      vec[i+1 + n*(j+5)] += fac*h[j+5 + m*(i+1)];
      vec[i+2 + n*(j+5)] += fac*h[j+5 + m*(i+2)];
      vec[i+3 + n*(j+5)] += fac*h[j+5 + m*(i+3)];
      vec[i+4 + n*(j+5)] += fac*h[j+5 + m*(i+4)];
      vec[i+5 + n*(j+5)] += fac*h[j+5 + m*(i+5)];
      vec[i+6 + n*(j+5)] += fac*h[j+5 + m*(i+6)];
      vec[i+7 + n*(j+5)] += fac*h[j+5 + m*(i+7)];
      vec[i+8 + n*(j+5)] += fac*h[j+5 + m*(i+8)];
      vec[i+9 + n*(j+5)] += fac*h[j+5 + m*(i+9)];
      vec[i   + n*(j+6)] += fac*h[j+6 + m*(i  )];
      vec[i+1 + n*(j+6)] += fac*h[j+6 + m*(i+1)];
      vec[i+2 + n*(j+6)] += fac*h[j+6 + m*(i+2)];
      vec[i+3 + n*(j+6)] += fac*h[j+6 + m*(i+3)];
      vec[i+4 + n*(j+6)] += fac*h[j+6 + m*(i+4)];
      vec[i+5 + n*(j+6)] += fac*h[j+6 + m*(i+5)];
      vec[i+6 + n*(j+6)] += fac*h[j+6 + m*(i+6)];
      vec[i+7 + n*(j+6)] += fac*h[j+6 + m*(i+7)];
      vec[i+8 + n*(j+6)] += fac*h[j+6 + m*(i+8)];
      vec[i+9 + n*(j+6)] += fac*h[j+6 + m*(i+9)];
      vec[i   + n*(j+7)] += fac*h[j+7 + m*(i  )];
      vec[i+1 + n*(j+7)] += fac*h[j+7 + m*(i+1)];
      vec[i+2 + n*(j+7)] += fac*h[j+7 + m*(i+2)];
      vec[i+3 + n*(j+7)] += fac*h[j+7 + m*(i+3)];
      vec[i+4 + n*(j+7)] += fac*h[j+7 + m*(i+4)];
      vec[i+5 + n*(j+7)] += fac*h[j+7 + m*(i+5)];
      vec[i+6 + n*(j+7)] += fac*h[j+7 + m*(i+6)];
      vec[i+7 + n*(j+7)] += fac*h[j+7 + m*(i+7)];
      vec[i+8 + n*(j+7)] += fac*h[j+7 + m*(i+8)];
      vec[i+9 + n*(j+7)] += fac*h[j+7 + m*(i+9)];
      vec[i   + n*(j+8)] += fac*h[j+8 + m*(i  )];
      vec[i+1 + n*(j+8)] += fac*h[j+8 + m*(i+1)];
      vec[i+2 + n*(j+8)] += fac*h[j+8 + m*(i+2)];
      vec[i+3 + n*(j+8)] += fac*h[j+8 + m*(i+3)];
      vec[i+4 + n*(j+8)] += fac*h[j+8 + m*(i+4)];
      vec[i+5 + n*(j+8)] += fac*h[j+8 + m*(i+5)];
      vec[i+6 + n*(j+8)] += fac*h[j+8 + m*(i+6)];
      vec[i+7 + n*(j+8)] += fac*h[j+8 + m*(i+7)];
      vec[i+8 + n*(j+8)] += fac*h[j+8 + m*(i+8)];
      vec[i+9 + n*(j+8)] += fac*h[j+8 + m*(i+9)];
      vec[i   + n*(j+9)] += fac*h[j+9 + m*(i  )];
      vec[i+1 + n*(j+9)] += fac*h[j+9 + m*(i+1)];
      vec[i+2 + n*(j+9)] += fac*h[j+9 + m*(i+2)];
      vec[i+3 + n*(j+9)] += fac*h[j+9 + m*(i+3)];
      vec[i+4 + n*(j+9)] += fac*h[j+9 + m*(i+4)];
      vec[i+5 + n*(j+9)] += fac*h[j+9 + m*(i+5)];
      vec[i+6 + n*(j+9)] += fac*h[j+9 + m*(i+6)];
      vec[i+7 + n*(j+9)] += fac*h[j+9 + m*(i+7)];
      vec[i+8 + n*(j+9)] += fac*h[j+9 + m*(i+8)];
      vec[i+9 + n*(j+9)] += fac*h[j+9 + m*(i+9)];
    }
    for (int j = mlim; j != m; ++j) {
      vec[i   + n*j] += fac*h[j  + m*(i  )];
      vec[i+1 + n*j] += fac*h[j  + m*(i+1)];
      vec[i+2 + n*j] += fac*h[j  + m*(i+2)];
      vec[i+3 + n*j] += fac*h[j  + m*(i+3)];
      vec[i+4 + n*j] += fac*h[j  + m*(i+4)];
      vec[i+5 + n*j] += fac*h[j  + m*(i+5)];
      vec[i+6 + n*j] += fac*h[j  + m*(i+6)];
      vec[i+7 + n*j] += fac*h[j  + m*(i+7)];
      vec[i+8 + n*j] += fac*h[j  + m*(i+8)];
      vec[i+9 + n*j] += fac*h[j  + m*(i+9)];
    }
  }
  for (int i = nlim; i != n; ++i) {
    for (int j = 0; j < mlim; j += 10) {
      vec[i   + n*(j  )] += fac*h[j   + m*(i  )];
      vec[i   + n*(j+1)] += fac*h[j+1 + m*(i  )];
      vec[i   + n*(j+2)] += fac*h[j+2 + m*(i  )];
      vec[i   + n*(j+3)] += fac*h[j+3 + m*(i  )];
      vec[i   + n*(j+4)] += fac*h[j+4 + m*(i  )];
      vec[i   + n*(j+5)] += fac*h[j+5 + m*(i  )];
      vec[i   + n*(j+6)] += fac*h[j+6 + m*(i  )];
      vec[i   + n*(j+7)] += fac*h[j+7 + m*(i  )];
      vec[i   + n*(j+8)] += fac*h[j+8 + m*(i  )];
      vec[i   + n*(j+9)] += fac*h[j+9 + m*(i  )];
    }
    for (int j = mlim; j !=  m; ++j) {
      vec[i   + n*(j  )] += fac*h[j   + m*(i  )];
    }
  }
}


template<>
void transpose_conjg(const complex<double>* h, const int m, const int n, complex<double>* vec, const double fac) {
  transpose_conjg(h, m, n, vec, static_cast<complex<double>>(fac));
}

template<>
void transpose_conjg(const complex<double>* h, const int m, const int n, complex<double>* vec, const complex<double> fac) {
#ifdef HAVE_MKL_H
  mkl_zomatcopy_("c", "c", &m, &n, &fac, h, &m, vec, &n);
#else
  const int mresidual = m % 10;
  const int nresidual = n % 10;
  const int mlim = m - mresidual;
  const int nlim = n - nresidual;
  for (int i = 0; i < nlim; i += 10) {
    for (int j = 0; j < mlim; j += 10) {
      vec[i   + n*(j  )] = fac*conj(h[j   + m*(i  )]);
      vec[i+1 + n*(j  )] = fac*conj(h[j   + m*(i+1)]);
      vec[i+2 + n*(j  )] = fac*conj(h[j   + m*(i+2)]);
      vec[i+3 + n*(j  )] = fac*conj(h[j   + m*(i+3)]);
      vec[i+4 + n*(j  )] = fac*conj(h[j   + m*(i+4)]);
      vec[i+5 + n*(j  )] = fac*conj(h[j   + m*(i+5)]);
      vec[i+6 + n*(j  )] = fac*conj(h[j   + m*(i+6)]);
      vec[i+7 + n*(j  )] = fac*conj(h[j   + m*(i+7)]);
      vec[i+8 + n*(j  )] = fac*conj(h[j   + m*(i+8)]);
      vec[i+9 + n*(j  )] = fac*conj(h[j   + m*(i+9)]);
      vec[i   + n*(j+1)] = fac*conj(h[j+1 + m*(i  )]);
      vec[i+1 + n*(j+1)] = fac*conj(h[j+1 + m*(i+1)]);
      vec[i+2 + n*(j+1)] = fac*conj(h[j+1 + m*(i+2)]);
      vec[i+3 + n*(j+1)] = fac*conj(h[j+1 + m*(i+3)]);
      vec[i+4 + n*(j+1)] = fac*conj(h[j+1 + m*(i+4)]);
      vec[i+5 + n*(j+1)] = fac*conj(h[j+1 + m*(i+5)]);
      vec[i+6 + n*(j+1)] = fac*conj(h[j+1 + m*(i+6)]);
      vec[i+7 + n*(j+1)] = fac*conj(h[j+1 + m*(i+7)]);
      vec[i+8 + n*(j+1)] = fac*conj(h[j+1 + m*(i+8)]);
      vec[i+9 + n*(j+1)] = fac*conj(h[j+1 + m*(i+9)]);
      vec[i   + n*(j+2)] = fac*conj(h[j+2 + m*(i  )]);
      vec[i+1 + n*(j+2)] = fac*conj(h[j+2 + m*(i+1)]);
      vec[i+2 + n*(j+2)] = fac*conj(h[j+2 + m*(i+2)]);
      vec[i+3 + n*(j+2)] = fac*conj(h[j+2 + m*(i+3)]);
      vec[i+4 + n*(j+2)] = fac*conj(h[j+2 + m*(i+4)]);
      vec[i+5 + n*(j+2)] = fac*conj(h[j+2 + m*(i+5)]);
      vec[i+6 + n*(j+2)] = fac*conj(h[j+2 + m*(i+6)]);
      vec[i+7 + n*(j+2)] = fac*conj(h[j+2 + m*(i+7)]);
      vec[i+8 + n*(j+2)] = fac*conj(h[j+2 + m*(i+8)]);
      vec[i+9 + n*(j+2)] = fac*conj(h[j+2 + m*(i+9)]);
      vec[i   + n*(j+3)] = fac*conj(h[j+3 + m*(i  )]);
      vec[i+1 + n*(j+3)] = fac*conj(h[j+3 + m*(i+1)]);
      vec[i+2 + n*(j+3)] = fac*conj(h[j+3 + m*(i+2)]);
      vec[i+3 + n*(j+3)] = fac*conj(h[j+3 + m*(i+3)]);
      vec[i+4 + n*(j+3)] = fac*conj(h[j+3 + m*(i+4)]);
      vec[i+5 + n*(j+3)] = fac*conj(h[j+3 + m*(i+5)]);
      vec[i+6 + n*(j+3)] = fac*conj(h[j+3 + m*(i+6)]);
      vec[i+7 + n*(j+3)] = fac*conj(h[j+3 + m*(i+7)]);
      vec[i+8 + n*(j+3)] = fac*conj(h[j+3 + m*(i+8)]);
      vec[i+9 + n*(j+3)] = fac*conj(h[j+3 + m*(i+9)]);
      vec[i   + n*(j+4)] = fac*conj(h[j+4 + m*(i  )]);
      vec[i+1 + n*(j+4)] = fac*conj(h[j+4 + m*(i+1)]);
      vec[i+2 + n*(j+4)] = fac*conj(h[j+4 + m*(i+2)]);
      vec[i+3 + n*(j+4)] = fac*conj(h[j+4 + m*(i+3)]);
      vec[i+4 + n*(j+4)] = fac*conj(h[j+4 + m*(i+4)]);
      vec[i+5 + n*(j+4)] = fac*conj(h[j+4 + m*(i+5)]);
      vec[i+6 + n*(j+4)] = fac*conj(h[j+4 + m*(i+6)]);
      vec[i+7 + n*(j+4)] = fac*conj(h[j+4 + m*(i+7)]);
      vec[i+8 + n*(j+4)] = fac*conj(h[j+4 + m*(i+8)]);
      vec[i+9 + n*(j+4)] = fac*conj(h[j+4 + m*(i+9)]);
      vec[i   + n*(j+5)] = fac*conj(h[j+5 + m*(i  )]);
      vec[i+1 + n*(j+5)] = fac*conj(h[j+5 + m*(i+1)]);
      vec[i+2 + n*(j+5)] = fac*conj(h[j+5 + m*(i+2)]);
      vec[i+3 + n*(j+5)] = fac*conj(h[j+5 + m*(i+3)]);
      vec[i+4 + n*(j+5)] = fac*conj(h[j+5 + m*(i+4)]);
      vec[i+5 + n*(j+5)] = fac*conj(h[j+5 + m*(i+5)]);
      vec[i+6 + n*(j+5)] = fac*conj(h[j+5 + m*(i+6)]);
      vec[i+7 + n*(j+5)] = fac*conj(h[j+5 + m*(i+7)]);
      vec[i+8 + n*(j+5)] = fac*conj(h[j+5 + m*(i+8)]);
      vec[i+9 + n*(j+5)] = fac*conj(h[j+5 + m*(i+9)]);
      vec[i   + n*(j+6)] = fac*conj(h[j+6 + m*(i  )]);
      vec[i+1 + n*(j+6)] = fac*conj(h[j+6 + m*(i+1)]);
      vec[i+2 + n*(j+6)] = fac*conj(h[j+6 + m*(i+2)]);
      vec[i+3 + n*(j+6)] = fac*conj(h[j+6 + m*(i+3)]);
      vec[i+4 + n*(j+6)] = fac*conj(h[j+6 + m*(i+4)]);
      vec[i+5 + n*(j+6)] = fac*conj(h[j+6 + m*(i+5)]);
      vec[i+6 + n*(j+6)] = fac*conj(h[j+6 + m*(i+6)]);
      vec[i+7 + n*(j+6)] = fac*conj(h[j+6 + m*(i+7)]);
      vec[i+8 + n*(j+6)] = fac*conj(h[j+6 + m*(i+8)]);
      vec[i+9 + n*(j+6)] = fac*conj(h[j+6 + m*(i+9)]);
      vec[i   + n*(j+7)] = fac*conj(h[j+7 + m*(i  )]);
      vec[i+1 + n*(j+7)] = fac*conj(h[j+7 + m*(i+1)]);
      vec[i+2 + n*(j+7)] = fac*conj(h[j+7 + m*(i+2)]);
      vec[i+3 + n*(j+7)] = fac*conj(h[j+7 + m*(i+3)]);
      vec[i+4 + n*(j+7)] = fac*conj(h[j+7 + m*(i+4)]);
      vec[i+5 + n*(j+7)] = fac*conj(h[j+7 + m*(i+5)]);
      vec[i+6 + n*(j+7)] = fac*conj(h[j+7 + m*(i+6)]);
      vec[i+7 + n*(j+7)] = fac*conj(h[j+7 + m*(i+7)]);
      vec[i+8 + n*(j+7)] = fac*conj(h[j+7 + m*(i+8)]);
      vec[i+9 + n*(j+7)] = fac*conj(h[j+7 + m*(i+9)]);
      vec[i   + n*(j+8)] = fac*conj(h[j+8 + m*(i  )]);
      vec[i+1 + n*(j+8)] = fac*conj(h[j+8 + m*(i+1)]);
      vec[i+2 + n*(j+8)] = fac*conj(h[j+8 + m*(i+2)]);
      vec[i+3 + n*(j+8)] = fac*conj(h[j+8 + m*(i+3)]);
      vec[i+4 + n*(j+8)] = fac*conj(h[j+8 + m*(i+4)]);
      vec[i+5 + n*(j+8)] = fac*conj(h[j+8 + m*(i+5)]);
      vec[i+6 + n*(j+8)] = fac*conj(h[j+8 + m*(i+6)]);
      vec[i+7 + n*(j+8)] = fac*conj(h[j+8 + m*(i+7)]);
      vec[i+8 + n*(j+8)] = fac*conj(h[j+8 + m*(i+8)]);
      vec[i+9 + n*(j+8)] = fac*conj(h[j+8 + m*(i+9)]);
      vec[i   + n*(j+9)] = fac*conj(h[j+9 + m*(i  )]);
      vec[i+1 + n*(j+9)] = fac*conj(h[j+9 + m*(i+1)]);
      vec[i+2 + n*(j+9)] = fac*conj(h[j+9 + m*(i+2)]);
      vec[i+3 + n*(j+9)] = fac*conj(h[j+9 + m*(i+3)]);
      vec[i+4 + n*(j+9)] = fac*conj(h[j+9 + m*(i+4)]);
      vec[i+5 + n*(j+9)] = fac*conj(h[j+9 + m*(i+5)]);
      vec[i+6 + n*(j+9)] = fac*conj(h[j+9 + m*(i+6)]);
      vec[i+7 + n*(j+9)] = fac*conj(h[j+9 + m*(i+7)]);
      vec[i+8 + n*(j+9)] = fac*conj(h[j+9 + m*(i+8)]);
      vec[i+9 + n*(j+9)] = fac*conj(h[j+9 + m*(i+9)]);
    }
    for (int j = mlim; j != m; ++j) {
      vec[i   + n*j] = fac*conj(h[j  + m*(i  )]);
      vec[i+1 + n*j] = fac*conj(h[j  + m*(i+1)]);
      vec[i+2 + n*j] = fac*conj(h[j  + m*(i+2)]);
      vec[i+3 + n*j] = fac*conj(h[j  + m*(i+3)]);
      vec[i+4 + n*j] = fac*conj(h[j  + m*(i+4)]);
      vec[i+5 + n*j] = fac*conj(h[j  + m*(i+5)]);
      vec[i+6 + n*j] = fac*conj(h[j  + m*(i+6)]);
      vec[i+7 + n*j] = fac*conj(h[j  + m*(i+7)]);
      vec[i+8 + n*j] = fac*conj(h[j  + m*(i+8)]);
      vec[i+9 + n*j] = fac*conj(h[j  + m*(i+9)]);
    }
  }
  for (int i = nlim; i != n; ++i) {
    for (int j = 0; j < mlim; j += 10) {
      vec[i   + n*(j  )] = fac*conj(h[j   + m*(i  )]);
      vec[i   + n*(j+1)] = fac*conj(h[j+1 + m*(i  )]);
      vec[i   + n*(j+2)] = fac*conj(h[j+2 + m*(i  )]);
      vec[i   + n*(j+3)] = fac*conj(h[j+3 + m*(i  )]);
      vec[i   + n*(j+4)] = fac*conj(h[j+4 + m*(i  )]);
      vec[i   + n*(j+5)] = fac*conj(h[j+5 + m*(i  )]);
      vec[i   + n*(j+6)] = fac*conj(h[j+6 + m*(i  )]);
      vec[i   + n*(j+7)] = fac*conj(h[j+7 + m*(i  )]);
      vec[i   + n*(j+8)] = fac*conj(h[j+8 + m*(i  )]);
      vec[i   + n*(j+9)] = fac*conj(h[j+9 + m*(i  )]);
    }
    for (int j = mlim; j !=  m; ++j) {
      vec[i   + n*(j  )] = fac*conj(h[j   + m*(i  )]);
    }
  }
#endif
}

}}
