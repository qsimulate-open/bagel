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

  const mpreal sqrtt = sqrt(T);
  // target Fm(T)
  const mpreal pi = GMPPI;
  const mpreal halfpT = half / T;
  fm[0] = sqrt(pi) / sqrtt * half * erf(sqrtt);
  for (int i = 1; i <= nrank; ++i) {
    fm[i] = halfpT * ((2*i-1) * fm[i - 1] - exp(-T));
  }
  out = fm[nrank];

  //const double boost_gamma = 0.5 * boost::math::tgamma_lower(nrank + 0.5, T.toDouble()) / std::pow(T.toDouble(), nrank+0.5);
  //assert(abs(boost_gamma - out) < 1e-14);
  //if (abs(boost_gamma - out) < 1e-16)
  //  cout << nrank << "   " << setprecision(16) << ta.toDouble() << " * " << boost_gamma << "  " << out << endl;

  return out;
}


static bool sort_vector(std::array<int, 3> v1, std::array<int, 3> v2) {
  int rad1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
  int rad2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
  return rad1 < rad2;
}


mpreal plm(const int l, const int am, const mpreal x) {

  const mpreal zero = "0.0";
  const mpreal one  = "1.0";
  const mpreal two  = "2.0";
  assert(am >= 0 && am <= l && abs(x) <= one);
  mpreal pmm = one;
  if (am > 0) {
    mpreal somx2 = sqrt((one - x)*(one + x));
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
  const mpreal one  = "1.0";
  const mpreal two  = "2.0";
  mpreal out = one;
  mpreal ft = one;
  for (int i = 1; i <= l; ++i) {
    out *= ft / two;
    ft += two;
  }

  return out * sqrt(GMPPI);
}


int main() {

  const mpreal pi = GMPPI;
  const int a = 10;
  const int n = std::pow(2*a+1, 3);
  vector<array<int, 3>> idx(n);
  int cnt = 0;
  for (int i = -a; i <= a; ++i)
    for (int j = -a; j <= a; ++j)
      for (int k = -a; k <= a; ++k, ++cnt)
        idx[cnt] = {{i, j, k}};

  std::sort(idx.begin(), idx.end(), sort_vector);

  const mpreal zero = "0.0";
  const mpreal half = "0.5";
  const mpreal one  = "1.0";
  const mpreal two  = "2.0";
  const int lmax = 20;
  const int ws = 2;
  vector<complex<mpreal>> mlm((lmax+1)*(lmax+1));
  for (int ivec = 1; ivec != n; ++ivec) {
    array<mpreal, 3> v;
    v[0] = static_cast<mpreal>(idx[ivec][0]);
    v[1] = static_cast<mpreal>(idx[ivec][1]);
    v[2] = static_cast<mpreal>(idx[ivec][2]);
    const mpreal rsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    const mpreal r = sqrt(rsq);
    const mpreal ctheta = (rsq > 0.0) ? v[2]/r : 0.0;
    const mpreal phi = atan2(v[1], v[0]);
    const mpreal b2r2 = pi * rsq;

    for(int l = 0; l <= lmax; ++l) {
      mpreal g = gamma(l);
      mpreal glower = boys(b2r2, l);
      const complex<mpreal> coeffl = std::pow(complex<mpreal>(zero, one), l) * pow(pi, l-half);
      const mpreal coeff = two * pow(sqrt(pi), 2*l+1) * pow(r, l) / g;
      const mpreal gupper = pow(r, -l-1) - glower * coeff;
      for (int m = -l; m <= l; ++m) {
        const int im = l * l + m + l;
        const int am = abs(m);
        mpreal plm_tilde = plm(l, am, ctheta);
        mpreal ft = one;
        for (int i = 1; i <= l - am; ++i) {
          plm_tilde *= ft;
          ft += one;
        }
        const mpreal sign = (m >=0) ? (cos(am * phi)) : (-one * cos(am * phi));

        if (abs(v[0]) > ws || abs(v[1]) > ws || abs(v[2]) > ws) {
          const mpreal real = gupper * sign * plm_tilde;
          const mpreal imag = gupper * sin(am * phi) * plm_tilde;
          mlm[im] += complex<mpreal>(real, imag);
        }
        if (abs(v[0]) <= ws && abs(v[1]) <= ws && abs(v[2]) <= ws) {
          const mpreal real = coeff * glower * sign * plm_tilde;
          const mpreal imag = coeff * glower * sin(am * phi) * plm_tilde;
          mlm[im] -= complex<mpreal>(real, imag);
        }
        const mpreal coeffm = plm_tilde * pow(r, l-2)*  exp(-rsq * pi) / g;
        mpreal real = coeffm * sign;
        mpreal imag = coeffm * sin(am * phi);
        mlm[im] += coeffl * complex<mpreal>(real, imag);
      }
    }
  }

  for(int l = 0; l <= lmax; ++l)
    for(int m = 0; m <= l; ++m) {
      const int i = l * l + m + l;
      //if (abs((mlm[i].real()).toDouble()) > 1e-10)
      if (l % 2 == 0 && m % 4 == 0)
        cout << l << "  " << m << "       " << scientific << setprecision(32) << mlm[i] << endl;
    }

  return 0;
}
