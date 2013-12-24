//
// BAGEL - Parallel electron correlation program.
// Filename: transpose.cc
// Copyright (C) 2008 Toru Shiozaki
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

#include <src/math/algo.h>
#include <bagel_config.h>
#ifdef HAVE_MKL_H
#define MKL_Complex16 std::complex<double>
#include <mkl.h>
#endif

using namespace std;

namespace bagel {
void mytranspose_(const double* h, const int m, const int n, double* vec, const double fac) {
#ifdef HAVE_MKL_H
  mkl_domatcopy('c', 't', m, n, fac, h, m, vec, n);
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
}
