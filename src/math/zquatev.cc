//
// Filename: zquatev.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: TS
//
// You can redistribute this program and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//

#include <src/math/algo.h>
#include <src/util/f77.h>
#include <cassert>
#include <algorithm>

using namespace std;

namespace bagel {

// some local functions..
static auto givens = [](const complex<double> a, const complex<double> b) {
  const double absa = abs(a);
  const double c = absa == 0.0 ? 0.0 : absa / sqrt(absa*absa + norm(b));
  const complex<double> s = absa == 0.0 ? 1.0 : (a / absa * conj(b) / sqrt(absa*absa + norm(b)));
  return make_pair(c, s);
};

static auto householder = [](const complex<double>* const hin, complex<double>* const out, const int len) {
  const bool trivial = abs(real(hin[0])) == 0.0;
  complex<double> fac = 0.0;
  copy_n(hin, len, out);
  if (!trivial) {
    const double norm = sqrt(real(zdotc_(len, hin, 1, hin, 1)));
    const double sign = real(hin[0])/abs(real(hin[0]));
    out[0] = hin[0] + sign*norm;
    fac = conj(1.0 / (conj(out[0]) * (sign*norm)));
  }
  return fac;
};

// implementation...

void zquatev_(const int n2, complex<double>* const D, double* const eig) {
  assert(n2 % 2 == 0);
  const int n = n2/2;

  // rearrange data
  complex<double>* const D0 = D;
  complex<double>* const D1 = D + n*n;
  copy_n(D, n2*n, D+n2*n);
  for (int i = 0; i != n; ++i) {
    copy_n(D+n2*n+i*n2, n, D0+i*n);
    copy_n(D+n2*n+i*n2+n, n, D1+i*n);
  }

  // identity matrix of n2 dimension
  complex<double>* const Q0 = D + n2*n;
  complex<double>* const Q1 = D + n2*n + n*n;
  fill_n(Q0, n*n, 0.0);
  fill_n(Q1, n*n, 0.0);
  for (int i = 0; i != n; ++i) Q0[i+n*i] = 1.0;

  unique_ptr<complex<double>[]> buf(new complex<double>[n2]);
  unique_ptr<complex<double>[]> hout(new complex<double>[n2]);
  unique_ptr<complex<double>[]> choutf(new complex<double>[n2]);


  // Reference - arXiv:1203.6151v4
  for (int k = 0; k != n-1; ++k) {
    const int len = n-k-1;
    if (len > 1) {
      complex<double>* const hin = D1+n*k+k+1;
      complex<double> tau = householder(hin, hout.get(), len);

      for (int i = 0; i != len; ++i) choutf[i] = conj(hout[i]);

      // 00-1
      zgemv_("T", len, len+1, 1.0, D0+k+1+(k)*n, n, hout.get(), 1, 0.0, buf.get(), 1);
      zgeru_(len, len+1, -conj(tau), choutf.get(), 1, buf.get(), 1, D0+k+1+(k)*n, n);

      // 00-2
      zgemv_("N", len+1, len, 1.0, D0+(k+1)*n+(k), n, choutf.get(), 1, 0.0, buf.get(), 1);
      zgeru_(len+1, len, -tau, buf.get(), 1, hout.get(), 1, D+(k+1)*n+(k), n);

      // 10-1
      zgemv_("T", len, len+1, -1.0, D1+k+1+(k)*n, n, choutf.get(), 1, 0.0, buf.get(), 1);
      zgeru_(len, len+1, tau, hout.get(), 1, buf.get(), 1, D1+k+1+(k)*n, n);

      // 10-2
      zgemv_("N", len+1, len, 1.0, D1+k+(k+1)*n, n, choutf.get(), 1, 0.0, buf.get(), 1);
      zgeru_(len+1, len, -tau, buf.get(), 1, hout.get(), 1, D1+(k+1)*n+(k), n);

      // 00-2
      zgemv_("N", n, len, 1.0, Q0+(k+1)*n, n, choutf.get(), 1, 0.0, buf.get(), 1);
      zgeru_(n, len, -tau, buf.get(), 1, hout.get(), 1, Q0+(k+1)*n, n);

      // 10-2
      zgemv_("N", n, len, 1.0, Q1+(k+1)*n, n, choutf.get(), 1, 0.0, buf.get(), 1);
      zgeru_(n, len, -tau, buf.get(), 1, hout.get(), 1, Q1+(k+1)*n, n);

    }

    // symplectic Givens rotation to clear out D(k+n, k)
    pair<double,complex<double>> gr = givens(D0[k+1+k*n], D1[k+1+k*n]);

    zrot_(len+1, D0+k+1+k*n, n, D1+k+1+k*n, n, gr.first, gr.second);

    for (int i = 0; i != len+1; ++i)
      D1[(k+1)*n+k+i] = -conj(D1[(k+1)*n+k+i]);
    zrot_(len+1, D0+(k+1)*n+k, 1, D1+(k+1)*n+k, 1, gr.first, conj(gr.second));
    for (int i = 0; i != len+1; ++i)
      D1[(k+1)*n+k+i] = -conj(D1[(k+1)*n+k+i]);

    for (int i = 0; i != n; ++i)
      Q1[(k+1)*n+i] = -conj(Q1[(k+1)*n+i]);
    zrot_(n, Q0+(k+1)*n, 1, Q1+(k+1)*n, 1, gr.first, conj(gr.second));
    for (int i = 0; i != n; ++i)
      Q1[(k+1)*n+i] = -conj(Q1[(k+1)*n+i]);

    // Householder to fix top half in column k
    if (len > 1) {
      complex<double>* const hin = D0+n*k+k+1;
      complex<double> tau = householder(hin, hout.get(), len);

      for (int i = 0; i != len; ++i) choutf[i] = conj(hout[i]);

      // 00-1
      zgemv_("C", len, len+1, 1.0, D0+k+1+(k)*n, n, hout.get(), 1, 0.0, buf.get(), 1);
      zgerc_(len, len+1, -tau, hout.get(), 1, buf.get(), 1, D0+k+1+(k)*n, n);

      // 00-2
      zgemv_("N", len+1, len, 1.0, D0+(k+1)*n+(k), n, hout.get(), 1, 0.0, buf.get(), 1);
      zgerc_(len+1, len, -conj(tau), buf.get(), 1, hout.get(), 1, D0+(k+1)*n+(k), n);

      // 01-1
      zgemv_("T", len, len+1, 1.0, D1+k+1+(k)*n, n, hout.get(), 1, 0.0, buf.get(), 1);
      zgeru_(len, len+1, -conj(tau), choutf.get(), 1, buf.get(), 1, D1+k+1+(k)*n, n);

      // 01-2
      zgemv_("N", len+1, len, -1.0, D1+(k+1)*n+(k), n, hout.get(), 1, 0.0, buf.get(), 1);
      zgerc_(len+1, len, conj(tau), buf.get(), 1, hout.get(), 1, D1+(k+1)*n+(k), n);

      // 00-2
      zgemv_("N", n, len, 1.0, Q0+(k+1)*n, n, hout.get(), 1, 0.0, buf.get(), 1);
      zgerc_(n, len, -conj(tau), buf.get(), 1, hout.get(), 1, Q0+(k+1)*n, n);

      // 01-2
      zgemv_("N", n, len, -1.0, Q1+(k+1)*n, n, hout.get(), 1, 0.0, buf.get(), 1);
      zgerc_(n, len, conj(tau), buf.get(), 1, hout.get(), 1, Q1+(k+1)*n, n);

    }

  }

  // diagonalize this tri-diagonal matrix (this step is much cheaper than
  // the Householder transformation above).
  unique_ptr<complex<double>[]> Cmat(new complex<double>[n*n]);
  unique_ptr<complex<double>[]> Work(new complex<double>[n]);
  int info;
  unique_ptr<double[]> rwork(new double[n*3]);
  for (int i = 0; i != n; ++i)
    for (int j = 0; j <= i; ++j)
      D0[i-j+j*n] = D0[i+j*n];
  zhbev_("V", "L", n, 1, D0, n, eig, Cmat.get(), n, Work.get(), rwork.get(), info);
  if (info) throw runtime_error("zhbev failed in quaternion diagonalization");

  // form the coefficient matrix in D
  zgemm3m_("N", "N", n, n, n, 1.0, Q0, n, Cmat.get(), n, 0.0, D, n2);
  zgemm3m_("N", "N", n, n, n, 1.0, Q1, n, Cmat.get(), n, 0.0, D+n, n2);

  // eigen vectors using symmetry
  for (int i = 0; i != n; ++i) {
    for (int j = 0; j != n; ++j) {
       D[j+n2*(i+n)] = -conj(D[j+n+n2*i]);
       D[j+n+n2*(i+n)] = conj(D[j+n2*i]);
    }
  }

  // eigen values using symmetry
  for (int i = 0; i != n; ++i) {
    eig[n+i] = eig[i];
  }
}

}
