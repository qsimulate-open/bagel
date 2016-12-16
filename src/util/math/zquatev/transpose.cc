//
// ZQUATEV: Diagonalization of quaternionic matrices
// File:    transpose.cc
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

#include <complex>
#include <bagel_config.h>

using namespace std;

#ifdef HAVE_MKL_H
extern "C" {
 void mkl_zomatcopy_(const char*, const char*, const int*, const int*, const complex<double>*, const complex<double>*, const int*, complex<double>*, const int*);
}
#endif

namespace ts {
namespace impl {

void transpose(const int m, const int n, const complex<double>* h, const int ld, complex<double>* vec, const int ldt) {
#ifdef HAVE_MKL_H
  complex<double> one = 1.0;
  mkl_zomatcopy_("c", "t", &m, &n, &one, h, &ld, vec, &ldt);
#else
  const int mresidual = m % 10;
  const int nresidual = n % 10;
  const int mlim = m - mresidual;
  const int nlim = n - nresidual;
  for (int i = 0; i < nlim; i += 10) {
    for (int j = 0; j < mlim; j += 10) {
      vec[i   + ldt*(j  )] = h[j   + ld*(i  )];
      vec[i+1 + ldt*(j  )] = h[j   + ld*(i+1)];
      vec[i+2 + ldt*(j  )] = h[j   + ld*(i+2)];
      vec[i+3 + ldt*(j  )] = h[j   + ld*(i+3)];
      vec[i+4 + ldt*(j  )] = h[j   + ld*(i+4)];
      vec[i+5 + ldt*(j  )] = h[j   + ld*(i+5)];
      vec[i+6 + ldt*(j  )] = h[j   + ld*(i+6)];
      vec[i+7 + ldt*(j  )] = h[j   + ld*(i+7)];
      vec[i+8 + ldt*(j  )] = h[j   + ld*(i+8)];
      vec[i+9 + ldt*(j  )] = h[j   + ld*(i+9)];
      vec[i   + ldt*(j+1)] = h[j+1 + ld*(i  )];
      vec[i+1 + ldt*(j+1)] = h[j+1 + ld*(i+1)];
      vec[i+2 + ldt*(j+1)] = h[j+1 + ld*(i+2)];
      vec[i+3 + ldt*(j+1)] = h[j+1 + ld*(i+3)];
      vec[i+4 + ldt*(j+1)] = h[j+1 + ld*(i+4)];
      vec[i+5 + ldt*(j+1)] = h[j+1 + ld*(i+5)];
      vec[i+6 + ldt*(j+1)] = h[j+1 + ld*(i+6)];
      vec[i+7 + ldt*(j+1)] = h[j+1 + ld*(i+7)];
      vec[i+8 + ldt*(j+1)] = h[j+1 + ld*(i+8)];
      vec[i+9 + ldt*(j+1)] = h[j+1 + ld*(i+9)];
      vec[i   + ldt*(j+2)] = h[j+2 + ld*(i  )];
      vec[i+1 + ldt*(j+2)] = h[j+2 + ld*(i+1)];
      vec[i+2 + ldt*(j+2)] = h[j+2 + ld*(i+2)];
      vec[i+3 + ldt*(j+2)] = h[j+2 + ld*(i+3)];
      vec[i+4 + ldt*(j+2)] = h[j+2 + ld*(i+4)];
      vec[i+5 + ldt*(j+2)] = h[j+2 + ld*(i+5)];
      vec[i+6 + ldt*(j+2)] = h[j+2 + ld*(i+6)];
      vec[i+7 + ldt*(j+2)] = h[j+2 + ld*(i+7)];
      vec[i+8 + ldt*(j+2)] = h[j+2 + ld*(i+8)];
      vec[i+9 + ldt*(j+2)] = h[j+2 + ld*(i+9)];
      vec[i   + ldt*(j+3)] = h[j+3 + ld*(i  )];
      vec[i+1 + ldt*(j+3)] = h[j+3 + ld*(i+1)];
      vec[i+2 + ldt*(j+3)] = h[j+3 + ld*(i+2)];
      vec[i+3 + ldt*(j+3)] = h[j+3 + ld*(i+3)];
      vec[i+4 + ldt*(j+3)] = h[j+3 + ld*(i+4)];
      vec[i+5 + ldt*(j+3)] = h[j+3 + ld*(i+5)];
      vec[i+6 + ldt*(j+3)] = h[j+3 + ld*(i+6)];
      vec[i+7 + ldt*(j+3)] = h[j+3 + ld*(i+7)];
      vec[i+8 + ldt*(j+3)] = h[j+3 + ld*(i+8)];
      vec[i+9 + ldt*(j+3)] = h[j+3 + ld*(i+9)];
      vec[i   + ldt*(j+4)] = h[j+4 + ld*(i  )];
      vec[i+1 + ldt*(j+4)] = h[j+4 + ld*(i+1)];
      vec[i+2 + ldt*(j+4)] = h[j+4 + ld*(i+2)];
      vec[i+3 + ldt*(j+4)] = h[j+4 + ld*(i+3)];
      vec[i+4 + ldt*(j+4)] = h[j+4 + ld*(i+4)];
      vec[i+5 + ldt*(j+4)] = h[j+4 + ld*(i+5)];
      vec[i+6 + ldt*(j+4)] = h[j+4 + ld*(i+6)];
      vec[i+7 + ldt*(j+4)] = h[j+4 + ld*(i+7)];
      vec[i+8 + ldt*(j+4)] = h[j+4 + ld*(i+8)];
      vec[i+9 + ldt*(j+4)] = h[j+4 + ld*(i+9)];
      vec[i   + ldt*(j+5)] = h[j+5 + ld*(i  )];
      vec[i+1 + ldt*(j+5)] = h[j+5 + ld*(i+1)];
      vec[i+2 + ldt*(j+5)] = h[j+5 + ld*(i+2)];
      vec[i+3 + ldt*(j+5)] = h[j+5 + ld*(i+3)];
      vec[i+4 + ldt*(j+5)] = h[j+5 + ld*(i+4)];
      vec[i+5 + ldt*(j+5)] = h[j+5 + ld*(i+5)];
      vec[i+6 + ldt*(j+5)] = h[j+5 + ld*(i+6)];
      vec[i+7 + ldt*(j+5)] = h[j+5 + ld*(i+7)];
      vec[i+8 + ldt*(j+5)] = h[j+5 + ld*(i+8)];
      vec[i+9 + ldt*(j+5)] = h[j+5 + ld*(i+9)];
      vec[i   + ldt*(j+6)] = h[j+6 + ld*(i  )];
      vec[i+1 + ldt*(j+6)] = h[j+6 + ld*(i+1)];
      vec[i+2 + ldt*(j+6)] = h[j+6 + ld*(i+2)];
      vec[i+3 + ldt*(j+6)] = h[j+6 + ld*(i+3)];
      vec[i+4 + ldt*(j+6)] = h[j+6 + ld*(i+4)];
      vec[i+5 + ldt*(j+6)] = h[j+6 + ld*(i+5)];
      vec[i+6 + ldt*(j+6)] = h[j+6 + ld*(i+6)];
      vec[i+7 + ldt*(j+6)] = h[j+6 + ld*(i+7)];
      vec[i+8 + ldt*(j+6)] = h[j+6 + ld*(i+8)];
      vec[i+9 + ldt*(j+6)] = h[j+6 + ld*(i+9)];
      vec[i   + ldt*(j+7)] = h[j+7 + ld*(i  )];
      vec[i+1 + ldt*(j+7)] = h[j+7 + ld*(i+1)];
      vec[i+2 + ldt*(j+7)] = h[j+7 + ld*(i+2)];
      vec[i+3 + ldt*(j+7)] = h[j+7 + ld*(i+3)];
      vec[i+4 + ldt*(j+7)] = h[j+7 + ld*(i+4)];
      vec[i+5 + ldt*(j+7)] = h[j+7 + ld*(i+5)];
      vec[i+6 + ldt*(j+7)] = h[j+7 + ld*(i+6)];
      vec[i+7 + ldt*(j+7)] = h[j+7 + ld*(i+7)];
      vec[i+8 + ldt*(j+7)] = h[j+7 + ld*(i+8)];
      vec[i+9 + ldt*(j+7)] = h[j+7 + ld*(i+9)];
      vec[i   + ldt*(j+8)] = h[j+8 + ld*(i  )];
      vec[i+1 + ldt*(j+8)] = h[j+8 + ld*(i+1)];
      vec[i+2 + ldt*(j+8)] = h[j+8 + ld*(i+2)];
      vec[i+3 + ldt*(j+8)] = h[j+8 + ld*(i+3)];
      vec[i+4 + ldt*(j+8)] = h[j+8 + ld*(i+4)];
      vec[i+5 + ldt*(j+8)] = h[j+8 + ld*(i+5)];
      vec[i+6 + ldt*(j+8)] = h[j+8 + ld*(i+6)];
      vec[i+7 + ldt*(j+8)] = h[j+8 + ld*(i+7)];
      vec[i+8 + ldt*(j+8)] = h[j+8 + ld*(i+8)];
      vec[i+9 + ldt*(j+8)] = h[j+8 + ld*(i+9)];
      vec[i   + ldt*(j+9)] = h[j+9 + ld*(i  )];
      vec[i+1 + ldt*(j+9)] = h[j+9 + ld*(i+1)];
      vec[i+2 + ldt*(j+9)] = h[j+9 + ld*(i+2)];
      vec[i+3 + ldt*(j+9)] = h[j+9 + ld*(i+3)];
      vec[i+4 + ldt*(j+9)] = h[j+9 + ld*(i+4)];
      vec[i+5 + ldt*(j+9)] = h[j+9 + ld*(i+5)];
      vec[i+6 + ldt*(j+9)] = h[j+9 + ld*(i+6)];
      vec[i+7 + ldt*(j+9)] = h[j+9 + ld*(i+7)];
      vec[i+8 + ldt*(j+9)] = h[j+9 + ld*(i+8)];
      vec[i+9 + ldt*(j+9)] = h[j+9 + ld*(i+9)];
    }
    for (int j = mlim; j != m; ++j) {
      vec[i   + ldt*j] = h[j  + ld*(i  )];
      vec[i+1 + ldt*j] = h[j  + ld*(i+1)];
      vec[i+2 + ldt*j] = h[j  + ld*(i+2)];
      vec[i+3 + ldt*j] = h[j  + ld*(i+3)];
      vec[i+4 + ldt*j] = h[j  + ld*(i+4)];
      vec[i+5 + ldt*j] = h[j  + ld*(i+5)];
      vec[i+6 + ldt*j] = h[j  + ld*(i+6)];
      vec[i+7 + ldt*j] = h[j  + ld*(i+7)];
      vec[i+8 + ldt*j] = h[j  + ld*(i+8)];
      vec[i+9 + ldt*j] = h[j  + ld*(i+9)];
    }
  }
  for (int i = nlim; i != n; ++i) {
    for (int j = 0; j < mlim; j += 10) {
      vec[i   + ldt*(j  )] = h[j   + ld*(i  )];
      vec[i   + ldt*(j+1)] = h[j+1 + ld*(i  )];
      vec[i   + ldt*(j+2)] = h[j+2 + ld*(i  )];
      vec[i   + ldt*(j+3)] = h[j+3 + ld*(i  )];
      vec[i   + ldt*(j+4)] = h[j+4 + ld*(i  )];
      vec[i   + ldt*(j+5)] = h[j+5 + ld*(i  )];
      vec[i   + ldt*(j+6)] = h[j+6 + ld*(i  )];
      vec[i   + ldt*(j+7)] = h[j+7 + ld*(i  )];
      vec[i   + ldt*(j+8)] = h[j+8 + ld*(i  )];
      vec[i   + ldt*(j+9)] = h[j+9 + ld*(i  )];
    }
    for (int j = mlim; j !=  m; ++j) {
      vec[i   + ldt*(j  )] = h[j   + ld*(i  )];
    }
  }
#endif
}


void transpose_conj(const int m, const int n, const complex<double>* h, const int ld, complex<double>* vec, const int ldt) {
#ifdef HAVE_MKL_H
  complex<double> one = 1.0;
  mkl_zomatcopy_("c", "c", &m, &n, &one, h, &ld, vec, &ldt);
#else
  const int mresidual = m % 10;
  const int nresidual = n % 10;
  const int mlim = m - mresidual;
  const int nlim = n - nresidual;
  for (int i = 0; i < nlim; i += 10) {
    for (int j = 0; j < mlim; j += 10) {
      vec[i   + ldt*(j  )] = conj(h[j   + ld*(i  )]);
      vec[i+1 + ldt*(j  )] = conj(h[j   + ld*(i+1)]);
      vec[i+2 + ldt*(j  )] = conj(h[j   + ld*(i+2)]);
      vec[i+3 + ldt*(j  )] = conj(h[j   + ld*(i+3)]);
      vec[i+4 + ldt*(j  )] = conj(h[j   + ld*(i+4)]);
      vec[i+5 + ldt*(j  )] = conj(h[j   + ld*(i+5)]);
      vec[i+6 + ldt*(j  )] = conj(h[j   + ld*(i+6)]);
      vec[i+7 + ldt*(j  )] = conj(h[j   + ld*(i+7)]);
      vec[i+8 + ldt*(j  )] = conj(h[j   + ld*(i+8)]);
      vec[i+9 + ldt*(j  )] = conj(h[j   + ld*(i+9)]);
      vec[i   + ldt*(j+1)] = conj(h[j+1 + ld*(i  )]);
      vec[i+1 + ldt*(j+1)] = conj(h[j+1 + ld*(i+1)]);
      vec[i+2 + ldt*(j+1)] = conj(h[j+1 + ld*(i+2)]);
      vec[i+3 + ldt*(j+1)] = conj(h[j+1 + ld*(i+3)]);
      vec[i+4 + ldt*(j+1)] = conj(h[j+1 + ld*(i+4)]);
      vec[i+5 + ldt*(j+1)] = conj(h[j+1 + ld*(i+5)]);
      vec[i+6 + ldt*(j+1)] = conj(h[j+1 + ld*(i+6)]);
      vec[i+7 + ldt*(j+1)] = conj(h[j+1 + ld*(i+7)]);
      vec[i+8 + ldt*(j+1)] = conj(h[j+1 + ld*(i+8)]);
      vec[i+9 + ldt*(j+1)] = conj(h[j+1 + ld*(i+9)]);
      vec[i   + ldt*(j+2)] = conj(h[j+2 + ld*(i  )]);
      vec[i+1 + ldt*(j+2)] = conj(h[j+2 + ld*(i+1)]);
      vec[i+2 + ldt*(j+2)] = conj(h[j+2 + ld*(i+2)]);
      vec[i+3 + ldt*(j+2)] = conj(h[j+2 + ld*(i+3)]);
      vec[i+4 + ldt*(j+2)] = conj(h[j+2 + ld*(i+4)]);
      vec[i+5 + ldt*(j+2)] = conj(h[j+2 + ld*(i+5)]);
      vec[i+6 + ldt*(j+2)] = conj(h[j+2 + ld*(i+6)]);
      vec[i+7 + ldt*(j+2)] = conj(h[j+2 + ld*(i+7)]);
      vec[i+8 + ldt*(j+2)] = conj(h[j+2 + ld*(i+8)]);
      vec[i+9 + ldt*(j+2)] = conj(h[j+2 + ld*(i+9)]);
      vec[i   + ldt*(j+3)] = conj(h[j+3 + ld*(i  )]);
      vec[i+1 + ldt*(j+3)] = conj(h[j+3 + ld*(i+1)]);
      vec[i+2 + ldt*(j+3)] = conj(h[j+3 + ld*(i+2)]);
      vec[i+3 + ldt*(j+3)] = conj(h[j+3 + ld*(i+3)]);
      vec[i+4 + ldt*(j+3)] = conj(h[j+3 + ld*(i+4)]);
      vec[i+5 + ldt*(j+3)] = conj(h[j+3 + ld*(i+5)]);
      vec[i+6 + ldt*(j+3)] = conj(h[j+3 + ld*(i+6)]);
      vec[i+7 + ldt*(j+3)] = conj(h[j+3 + ld*(i+7)]);
      vec[i+8 + ldt*(j+3)] = conj(h[j+3 + ld*(i+8)]);
      vec[i+9 + ldt*(j+3)] = conj(h[j+3 + ld*(i+9)]);
      vec[i   + ldt*(j+4)] = conj(h[j+4 + ld*(i  )]);
      vec[i+1 + ldt*(j+4)] = conj(h[j+4 + ld*(i+1)]);
      vec[i+2 + ldt*(j+4)] = conj(h[j+4 + ld*(i+2)]);
      vec[i+3 + ldt*(j+4)] = conj(h[j+4 + ld*(i+3)]);
      vec[i+4 + ldt*(j+4)] = conj(h[j+4 + ld*(i+4)]);
      vec[i+5 + ldt*(j+4)] = conj(h[j+4 + ld*(i+5)]);
      vec[i+6 + ldt*(j+4)] = conj(h[j+4 + ld*(i+6)]);
      vec[i+7 + ldt*(j+4)] = conj(h[j+4 + ld*(i+7)]);
      vec[i+8 + ldt*(j+4)] = conj(h[j+4 + ld*(i+8)]);
      vec[i+9 + ldt*(j+4)] = conj(h[j+4 + ld*(i+9)]);
      vec[i   + ldt*(j+5)] = conj(h[j+5 + ld*(i  )]);
      vec[i+1 + ldt*(j+5)] = conj(h[j+5 + ld*(i+1)]);
      vec[i+2 + ldt*(j+5)] = conj(h[j+5 + ld*(i+2)]);
      vec[i+3 + ldt*(j+5)] = conj(h[j+5 + ld*(i+3)]);
      vec[i+4 + ldt*(j+5)] = conj(h[j+5 + ld*(i+4)]);
      vec[i+5 + ldt*(j+5)] = conj(h[j+5 + ld*(i+5)]);
      vec[i+6 + ldt*(j+5)] = conj(h[j+5 + ld*(i+6)]);
      vec[i+7 + ldt*(j+5)] = conj(h[j+5 + ld*(i+7)]);
      vec[i+8 + ldt*(j+5)] = conj(h[j+5 + ld*(i+8)]);
      vec[i+9 + ldt*(j+5)] = conj(h[j+5 + ld*(i+9)]);
      vec[i   + ldt*(j+6)] = conj(h[j+6 + ld*(i  )]);
      vec[i+1 + ldt*(j+6)] = conj(h[j+6 + ld*(i+1)]);
      vec[i+2 + ldt*(j+6)] = conj(h[j+6 + ld*(i+2)]);
      vec[i+3 + ldt*(j+6)] = conj(h[j+6 + ld*(i+3)]);
      vec[i+4 + ldt*(j+6)] = conj(h[j+6 + ld*(i+4)]);
      vec[i+5 + ldt*(j+6)] = conj(h[j+6 + ld*(i+5)]);
      vec[i+6 + ldt*(j+6)] = conj(h[j+6 + ld*(i+6)]);
      vec[i+7 + ldt*(j+6)] = conj(h[j+6 + ld*(i+7)]);
      vec[i+8 + ldt*(j+6)] = conj(h[j+6 + ld*(i+8)]);
      vec[i+9 + ldt*(j+6)] = conj(h[j+6 + ld*(i+9)]);
      vec[i   + ldt*(j+7)] = conj(h[j+7 + ld*(i  )]);
      vec[i+1 + ldt*(j+7)] = conj(h[j+7 + ld*(i+1)]);
      vec[i+2 + ldt*(j+7)] = conj(h[j+7 + ld*(i+2)]);
      vec[i+3 + ldt*(j+7)] = conj(h[j+7 + ld*(i+3)]);
      vec[i+4 + ldt*(j+7)] = conj(h[j+7 + ld*(i+4)]);
      vec[i+5 + ldt*(j+7)] = conj(h[j+7 + ld*(i+5)]);
      vec[i+6 + ldt*(j+7)] = conj(h[j+7 + ld*(i+6)]);
      vec[i+7 + ldt*(j+7)] = conj(h[j+7 + ld*(i+7)]);
      vec[i+8 + ldt*(j+7)] = conj(h[j+7 + ld*(i+8)]);
      vec[i+9 + ldt*(j+7)] = conj(h[j+7 + ld*(i+9)]);
      vec[i   + ldt*(j+8)] = conj(h[j+8 + ld*(i  )]);
      vec[i+1 + ldt*(j+8)] = conj(h[j+8 + ld*(i+1)]);
      vec[i+2 + ldt*(j+8)] = conj(h[j+8 + ld*(i+2)]);
      vec[i+3 + ldt*(j+8)] = conj(h[j+8 + ld*(i+3)]);
      vec[i+4 + ldt*(j+8)] = conj(h[j+8 + ld*(i+4)]);
      vec[i+5 + ldt*(j+8)] = conj(h[j+8 + ld*(i+5)]);
      vec[i+6 + ldt*(j+8)] = conj(h[j+8 + ld*(i+6)]);
      vec[i+7 + ldt*(j+8)] = conj(h[j+8 + ld*(i+7)]);
      vec[i+8 + ldt*(j+8)] = conj(h[j+8 + ld*(i+8)]);
      vec[i+9 + ldt*(j+8)] = conj(h[j+8 + ld*(i+9)]);
      vec[i   + ldt*(j+9)] = conj(h[j+9 + ld*(i  )]);
      vec[i+1 + ldt*(j+9)] = conj(h[j+9 + ld*(i+1)]);
      vec[i+2 + ldt*(j+9)] = conj(h[j+9 + ld*(i+2)]);
      vec[i+3 + ldt*(j+9)] = conj(h[j+9 + ld*(i+3)]);
      vec[i+4 + ldt*(j+9)] = conj(h[j+9 + ld*(i+4)]);
      vec[i+5 + ldt*(j+9)] = conj(h[j+9 + ld*(i+5)]);
      vec[i+6 + ldt*(j+9)] = conj(h[j+9 + ld*(i+6)]);
      vec[i+7 + ldt*(j+9)] = conj(h[j+9 + ld*(i+7)]);
      vec[i+8 + ldt*(j+9)] = conj(h[j+9 + ld*(i+8)]);
      vec[i+9 + ldt*(j+9)] = conj(h[j+9 + ld*(i+9)]);
    }
    for (int j = mlim; j != m; ++j) {
      vec[i   + ldt*j] = conj(h[j  + ld*(i  )]);
      vec[i+1 + ldt*j] = conj(h[j  + ld*(i+1)]);
      vec[i+2 + ldt*j] = conj(h[j  + ld*(i+2)]);
      vec[i+3 + ldt*j] = conj(h[j  + ld*(i+3)]);
      vec[i+4 + ldt*j] = conj(h[j  + ld*(i+4)]);
      vec[i+5 + ldt*j] = conj(h[j  + ld*(i+5)]);
      vec[i+6 + ldt*j] = conj(h[j  + ld*(i+6)]);
      vec[i+7 + ldt*j] = conj(h[j  + ld*(i+7)]);
      vec[i+8 + ldt*j] = conj(h[j  + ld*(i+8)]);
      vec[i+9 + ldt*j] = conj(h[j  + ld*(i+9)]);
    }
  }
  for (int i = nlim; i != n; ++i) {
    for (int j = 0; j < mlim; j += 10) {
      vec[i   + ldt*(j  )] = conj(h[j   + ld*(i  )]);
      vec[i   + ldt*(j+1)] = conj(h[j+1 + ld*(i  )]);
      vec[i   + ldt*(j+2)] = conj(h[j+2 + ld*(i  )]);
      vec[i   + ldt*(j+3)] = conj(h[j+3 + ld*(i  )]);
      vec[i   + ldt*(j+4)] = conj(h[j+4 + ld*(i  )]);
      vec[i   + ldt*(j+5)] = conj(h[j+5 + ld*(i  )]);
      vec[i   + ldt*(j+6)] = conj(h[j+6 + ld*(i  )]);
      vec[i   + ldt*(j+7)] = conj(h[j+7 + ld*(i  )]);
      vec[i   + ldt*(j+8)] = conj(h[j+8 + ld*(i  )]);
      vec[i   + ldt*(j+9)] = conj(h[j+9 + ld*(i  )]);
    }
    for (int j = mlim; j !=  m; ++j) {
      vec[i   + ldt*(j  )] = conj(h[j   + ld*(i  )]);
    }
  }
#endif
}

}}
