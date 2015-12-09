//
// Author: Hai-Anh Le
// Date  : August 2015
//

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include "mpreal.h"
#include "gmp_macros.h"
#include <algorithm>
#include <vector>
#include <boost/math/special_functions/gamma.hpp>
#include <src/util/math/legendre.h>

using namespace std;
using namespace mpfr;

mpreal boys(const mpreal ta, const int nrank) {

  mpfr::mpreal::set_default_prec(GMPPREC);
  mpreal fm[nrank+1];

  const mpreal T = ta;
  assert(T > 0);

  mpreal out;

  const mpreal zero = "0.0";
  const mpreal half = "0.5";
  const mpreal quater = "0.25";

  const mpreal sqrtt = mpfr::sqrt(T);
  // target Fm(T)
  const mpreal pi = GMPPI;
  const mpreal halfpT = half / T;
  fm[0] = mpfr::sqrt(pi) / sqrtt * half * mpfr::erf(sqrtt);
  for (int i = 1; i <= nrank; ++i) {
    fm[i] = halfpT * ((2*i-1) * fm[i - 1] - mpfr::exp(-T));
  }
  out = fm[nrank];

  const double boost_gamma = 0.5 * boost::math::tgamma_lower(nrank + 0.5, T.toDouble()) / std::pow(T.toDouble(), nrank+0.5);
  assert(abs(boost_gamma - out) < 1e-14);
  if (abs(boost_gamma - out) > 1e-16)
    cout << " *** Warning: " << nrank << "   " << setprecision(16) << ta.toDouble() << " * " << boost_gamma << "  " << out << endl;

  return out;
}


static bool sort_vector(std::array<int, 3> v1, std::array<int, 3> v2) {
  int rad1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
  int rad2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
  return rad1 < rad2;
}


mpreal dot(array<mpreal, 3> b, array<mpreal, 3> c) {
  mpfr::mpreal::set_default_prec(GMPPREC);
  return b[0]*c[0]+b[1]*c[1]+b[2]*c[2];
}


array<mpreal, 3> cross(array<mpreal, 3> b, array<mpreal, 3> c, mpreal s = "1.0") {
  mpfr::mpreal::set_default_prec(GMPPREC);
  array<mpreal, 3> out;
  out[0] = (b[1]*c[2] - b[2]*c[1]) * s;
  out[1] = (b[2]*c[0] - b[0]*c[2]) * s;
  out[2] = (b[0]*c[1] - b[1]*c[0]) * s;
  return out;
}

mpreal plm(const int l, const int am, const mpreal x) {

  mpfr::mpreal::set_default_prec(GMPPREC);
  const mpreal zero = "0.0";
  const mpreal one  = "1.0";
  const mpreal two  = "2.0";
  assert(am >= 0 && am <= l && abs(x) <= one);
  mpreal pmm = one;
  if (am > 0) {
    mpreal somx2 = mpfr::sqrt((one - x)*(one + x));
    mpreal fact = one;
    for (int i = 1; i <= am; ++i) {
      pmm *= -fact * somx2;
      fact += two;
    }
  }
  if (l == am) {
    return pmm;
  } else {
    mpreal pmmp1 = x * (two * am + one) * pmm;
    if (l == am + 1) {
      return pmmp1;
    } else {
      mpreal out = zero;
      for (int i = am + 2; i <= l; ++i) {
        out = (x * (2 * i -1) * pmmp1 - (i + am - 1) * pmm) / (i - am);
        pmm = pmmp1;
        pmmp1 = out;
      }
      return out;
    }
  }
}


mpreal gamma(const int l) {
  mpfr::mpreal::set_default_prec(GMPPREC);
  const mpreal one  = "1.0";
  const mpreal two  = "2.0";
  mpreal out = one;
  mpreal ft = one;
  for (int i = 1; i <= l; ++i) {
    out *= ft / two;
    ft += two;
  }

  return out * mpfr::sqrt(GMPPI);
}


int main() {

  mpfr::mpreal::set_default_prec(GMPPREC);
  const mpreal pi = GMPPI;
  const mpreal beta = GMPPISQRT;
  const mpreal zero = "0.0";
  const mpreal half = "0.5";
  const mpreal one  = "1.0";
  const mpreal two  = "2.0";

  vector<array<mpreal, 3>> primvecs(3);
  primvecs[0] = {{one, zero, zero}};
  primvecs[1] = {{zero, one, zero}};
  primvecs[2] = {{zero, zero, "0.0"}};
  const int a = 10;
  const int ndim = 3;
  const int nvec = std::pow(2*a+1, ndim);
  const int lmax = 20;
  const int ws = 1;
  vector<complex<mpreal>> mlm((lmax+1)*(lmax+1));

  vector<array<int, 3>> vidx(nvec);
  int cnt = 0;
  if (ndim == 3) {
    for (int i = -a; i <= a; ++i)
      for (int j = -a; j <= a; ++j)
        for (int k = -a; k <= a; ++k, ++cnt)
          vidx[cnt] = {{k, j, i}};
  } else if (ndim == 2) {
    for (int j = -a; j <= a; ++j)
      for (int k = -a; k <= a; ++k, ++cnt)
        vidx[cnt] = {{k, j, 0}};
  } else if (ndim == 1) {
    for (int k = -a; k <= a; ++k, ++cnt)
      vidx[cnt] = {{k, 0, 0}};
  }
  assert(cnt == nvec);
//  std::sort(vidx.begin(), vidx.end(), sort_vector);

#if 1
  for(int l = 0; l <= lmax; ++l) {
    mpreal g = gamma(l);
    const int rank = l + 1;

    for (int ivec = 0; ivec != nvec; ++ivec) {
      array<mpreal, 3> v;
      array<int, 3> idx = vidx[ivec];
      v[0] = idx[0] * primvecs[0][0] + idx[1] * primvecs[1][0] + idx[2] * primvecs[2][0];
      v[1] = idx[0] * primvecs[0][1] + idx[1] * primvecs[1][1] + idx[2] * primvecs[2][1];
      v[2] = idx[0] * primvecs[0][2] + idx[1] * primvecs[1][2] + idx[2] * primvecs[2][2];
      const mpreal rsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
      if (rsq > zero) {
        const mpreal r = mpfr::sqrt(rsq);
        const mpreal ctheta = (rsq > 0.0) ? v[2]/r : 0.0;
        const mpreal phi = mpfr::atan2(v[1], v[0]);
        const mpreal b2r2 = beta * beta * rsq;
        mpreal glower = boys(b2r2, l);

        const complex<mpreal> coeffl = std::pow(complex<mpreal>(zero, one), l) * mpfr::pow(pi, l-half);
        const mpreal coeff = two * mpfr::pow(beta, 2*l+1) * pow(r, l) / g;
        glower = glower * coeff;
        const mpreal gupper = 1.0/mpfr::pow(r, l+1.0) - glower;

        for (int m = -l; m <= l; ++m) {
          const int im = l * l + m + l;
          const int am = abs(m);
          mpreal plm_tilde = plm(l, am, ctheta);

#if 0
          const static bagel::Legendre plm;
          const double plmbagel = plm.compute(l, am, ctheta.toDouble());
          const double err = abs(plm_tilde.toDouble() - plmbagel);
          if (err > 1e-14)
            cout << "*** Warning (1): " << l << ", " << am << ", " << ctheta.toDouble() << "  " << setprecision(20) << plmbagel << "  " << plm_tilde.toDouble() << endl;
#endif

          mpreal ft = one;
          for (int i = 1; i <= l - am; ++i) {
            plm_tilde *= ft;
            ft += one;
          }
          const mpreal sign = (m >=0) ? (mpfr::cos(am * phi)) : (-one * mpfr::cos(am * phi));

          if (abs(idx[0]) > ws || abs(idx[1]) > ws || abs(idx[2]) > ws) {
            const mpreal real = gupper * sign * plm_tilde;
            const mpreal imag = gupper * mpfr::sin(am * phi) * plm_tilde;
            mlm[im] += complex<mpreal>(real, imag);
          }
          if (abs(idx[0]) <= ws && abs(idx[1]) <= ws && abs(idx[2]) <= ws) {
            const mpreal real = glower * sign * plm_tilde;
            const mpreal imag = glower * mpfr::sin(am * phi) * plm_tilde;
            mlm[im] -= complex<mpreal>(real, imag);
          }
        }
      }
    }
  }
#endif

#if 1
  vector<array<mpreal, 3>> primkvecs(3);
  switch (ndim) {
    case 1:
      {
        const mpreal a1sq = dot(primvecs[0], primvecs[0]);
        for (int i = 0; i != 3; ++i)
          primkvecs[0][i] = primvecs[0][i] / a1sq;
      }
    case 2:
      {
        array<mpreal, 3> a12 = cross(primvecs[0], primvecs[1]);
        const mpreal scale = one / dot(a12, a12);
        primkvecs[0] = cross(primvecs[1], a12, scale);
        primkvecs[1] = cross(a12, primvecs[0], scale);
      }
    case 3:
      {
        array<mpreal, 3> a23 = cross(primvecs[1], primvecs[2]);
        const mpreal scale = one / dot(primvecs[0], a23);
        primkvecs[0] = cross(primvecs[1], primvecs[2], scale);
        primkvecs[1] = cross(primvecs[2], primvecs[0], scale);
        primkvecs[2] = cross(primvecs[0], primvecs[1], scale);
      }
  }
  for (int i = ndim; i != 3; ++i)
    primkvecs[i] = {{zero, zero, zero}};

  for (int ivec = 0; ivec != nvec; ++ivec) {
    array<mpreal, 3> kvec;
    array<int, 3> idx = vidx[ivec];
    kvec[0] = idx[0] * primkvecs[0][0] + idx[1] * primkvecs[1][0] + idx[2] * primkvecs[2][0];
    kvec[1] = idx[0] * primkvecs[0][1] + idx[1] * primkvecs[1][1] + idx[2] * primkvecs[2][1];
    kvec[2] = idx[0] * primkvecs[0][2] + idx[1] * primkvecs[1][2] + idx[2] * primkvecs[2][2];
    const mpreal rsq = kvec[0]*kvec[0] + kvec[1]*kvec[1] + kvec[2]*kvec[2];
    const mpreal r = mpfr::sqrt(rsq);
    const mpreal ctheta = kvec[2]/r;
    const mpreal phi = mpfr::atan2(kvec[1], kvec[0]);

    if (rsq > zero) {
      for (int l = 0; l <= lmax; ++l) {
        const complex<mpreal> coeffl = std::pow(complex<mpreal>(zero, one), l) * mpfr::pow(pi, l-half);

        for (int m = -l; m <= l; ++m) {
          const int am = abs(m);
          const int im = l * l + m + l;

          mpreal plm_tilde = plm(l, am, ctheta);
#if 0
          const static bagel::Legendre plm;
          const double plmbagel = plm.compute(l, am, ctheta.toDouble());
          const double err = abs(plm_tilde.toDouble() - plmbagel);
          if (err > 1e-14)
            cout << "*** Warning (2): " << l << ", " << am << ", " << ctheta.toDouble() << "  " << setprecision(20) << plmbagel << "  " << plm_tilde.toDouble() << endl;
#endif

          mpreal ft = one;
          for (int i = 1; i <= l - am; ++i) {
            plm_tilde *= ft;
            ft += one;
          }
          const mpreal sign = (m >=0) ? (mpfr::cos(am * phi)) : (-one * mpfr::cos(am * phi));

          // smooth term
          const mpreal coeffm = plm_tilde * mpfr::pow(r, l-2)*  mpfr::exp(-rsq * pi) / gamma(l);
          mpreal real = coeffm * sign;
          mpreal imag = coeffm * mpfr::sin(am * phi);
          mlm[im] += coeffl * complex<mpreal>(real, imag);
        }
      }
    }
  }
#endif

  for(int l = 0; l <= lmax; ++l)
    for(int m = 0; m <= l; ++m) {
      const int i = l * l + m + l;
      if (l % 2 == 0 && m % 4 == 0)
        cout << l << "  " << m << "       " << scientific << setprecision(32) << mlm[i].real().toDouble() << "   " << mlm[i].imag().toDouble() << endl;
    }

  return 0;
}
