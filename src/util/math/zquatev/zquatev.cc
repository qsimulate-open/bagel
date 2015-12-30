//
// ZQUATEV: Diagonalization of quaternionic matrices
// File   : zquatev.cc
// Copyright (c) 2013, Toru Shiozaki (shiozaki@northwestern.edu)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies,
// either expressed or implied, of the FreeBSD Project.
//

#include "zquatev.h"
#include "f77.h"
#include <cassert>
#include <algorithm>

using namespace std;
using namespace ts::impl;

namespace ts {

int zquatev(const int n2, complex<double>* const D, const int ld2, double* const eig) {
  assert(n2 % 2 == 0);
  const int n = n2/2;
  const int ld = ld2/2;

  // rearrange data
  complex<double>* const D0 = D;
  complex<double>* const D1 = D + n*ld;
  copy_n(D, ld2*n, D+ld2*n);
  for (int i = 0; i != n; ++i) {
    copy_n(D+ld2*n+i*ld2, n, D0+i*ld);
    copy_n(D+ld2*n+i*ld2+n, n, D1+i*ld);
  }

  // identity matrix of n2 dimension
  complex<double>* const Q0 = D + ld2*n;
  complex<double>* const Q1 = D + ld2*n + ld*n;
  fill_n(Q0, ld*n, 0.0);
  fill_n(Q1, ld*n, 0.0);
  for (int i = 0; i != n; ++i) Q0[i+ld*i] = 1.0;

  // TODO allocation will move out
  const int nb = 20;
  const size_t alloc_size = nb*5 + nb*nb*15 + (max(4,nb)+10*nb+3)*n;
  unique_ptr<complex<double>[]> tmp_mem(new complex<double>[alloc_size]);

  for (int p = 0; p < n; p += nb)
    if (n-p > nb+1)
      panel_update(n-p, nb, D0+p*ld+p, D1+p*ld+p, Q0+p*ld, Q1+p*ld, ld, n, tmp_mem.get());
    else
      unblocked_update(n-p, D0+p*ld+p, D1+p*ld+p, Q0+p*ld, Q1+p*ld, ld, n, tmp_mem.get());

  // diagonalize this tri-diagonal matrix (this step is much cheaper than
  // the Householder transformation above).
  auto work1 = tmp_mem.get();
  auto work2 = work1 + n;
  auto work3 = work2 + n;
  for (int i = 0; i != n; ++i) {
    work3[2*i] = D0[i+i*ld];
    work3[2*i+1] = D0[i+i*ld+1];
  }
  int info;
  zhbev_("V", "L", n, 1, work3, 2, eig, D0+ld, ld2, work1, reinterpret_cast<double*>(work2), info);

  // form the coefficient matrix in D
  zgemm3m_("N", "N", n, n, n, 1.0, Q0, ld, D+ld, ld2, 0.0, D, ld2);
  for (int i = 0; i != n; ++i)
    copy_n(D+ld+i*ld2, n, Q0+i*ld);
  zgemm3m_("N", "N", n, n, n, 1.0, Q1, ld, Q0, n, 0.0, D+ld, ld2);

  // eigen vectors using symmetry
  for (int i = 0; i != n; ++i) {
    for (int j = 0; j != n; ++j) {
       D[j+ld2*(i+n)] = -conj(D[j+ld+ld2*i]);
       D[j+ld+ld2*(i+n)] = conj(D[j+ld2*i]);
    }
  }
  return info;
}

}
