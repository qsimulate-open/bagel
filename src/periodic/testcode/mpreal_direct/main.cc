//
// Author: Hai-Anh Le
// Date  : November 2015
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

mpreal plm(const int l, const int am, const mpreal x) {

  const mpreal zero = "0.0";
  const mpreal one  = "1.0";
  const mpreal two  = "2.0";
  //assert(am >= 0 && am <= l && abs(x) <= one);
  if (am < 0 || am > l || abs(x) > one)
    cout << "|m| = " << am << " l = " << l << " x = " << x << endl;
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


static bool sort_vector(std::array<int, 3> v1, std::array<int, 3> v2) {
  int rad1 = v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2];
  int rad2 = v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2];
  return rad1 < rad2;
}


int main() {

  const mpreal zero = "0.0";
  const mpreal half = "0.5";
  const mpreal one  = "1.0";
  const mpreal two  = "2.0";

  vector<array<mpreal, 3>> primvecs(3);
  primvecs[0] = {{one, zero, zero}};
  primvecs[1] = {{zero, one, zero}};
  primvecs[2] = {{zero, zero, one}};

  const int lmax = 10;
  const int max_rank = lmax*2 + 1;
  const int ws = 2;
  const int osize = (lmax+1)*(lmax+1);
  const int msize = (2*lmax+1)*(2*lmax+1);
  vector<complex<mpreal>> mstar(osize);

  // M* = sum of M in [-1, 1]
  const int a = 1;
  const int n0 = std::pow(2*a+1, 3);
  vector<array<int, 3>> vidx0(n0);
  int cnt = 0;
  for (int i = -a; i <= a; ++i)
    for (int j = -a; j <= a; ++j)
      for (int k = -a; k <= a; ++k, ++cnt)
        vidx0[cnt] = {{i, j, k}};
  std::sort(vidx0.begin(), vidx0.end(), sort_vector);

  for (int n = 0; n != n0; ++n) {
    array<int, 3> idx = vidx0[n];
    array<mpreal, 3> mvec;
    mvec[0] = idx[0] * primvecs[0][0] + idx[1] * primvecs[1][0] + idx[2] * primvecs[2][0];
    mvec[1] = idx[0] * primvecs[0][1] + idx[1] * primvecs[1][1] + idx[2] * primvecs[2][1];
    mvec[2] = idx[0] * primvecs[0][2] + idx[1] * primvecs[1][2] + idx[2] * primvecs[2][2];
    const mpreal rsq = mvec[0] * mvec[0] + mvec[1] * mvec[1] + mvec[2] * mvec[2];
    const mpreal r = sqrt(rsq);
    const mpreal ctheta = (r > zero) ? mvec[2]/r : zero;
    const mpreal phi = atan2(mvec[1], mvec[0]);

    for (int l = 0; l <= lmax; ++l) {
      for (int m = 0; m <= 2 * l; ++m) {
        const int am = abs(m - l);
        const int imul = l * l + m;

        mpreal plm_tilde = plm(l, am, ctheta) * pow(r, l);
        mpreal ft = one;
        for (int i = 1; i <= l + am; ++i) {
          plm_tilde /= ft;
          ft += one;
        }

        const mpreal sign = (m - l >= 0) ? (cos(am * phi)) : (-1.0 * cos(am * phi));
        const mpreal real = sign * plm_tilde;
        const mpreal imag = sin(am * phi) * plm_tilde;

        mstar[imul] += complex<mpreal>(real, imag);
      }
    }
  }

  cout << "*** Mstar ***" << endl;
  for(int l = 0; l <= lmax; ++l)
    for(int m = 0; m <= l; ++m) {
      const int i = l * l + m + l;
      if (l % 2 == 0 && m % 4 == 0)
        cout << l << "  " << m << "       " << scientific << setprecision(32) << (mstar[i].real()).toDouble() << endl;
    }
  cout << endl;


  // get L* = sum of L in FF'
  const int ws1 = 3 * ws + 1;
  const int n1 = std::pow(2*ws1+1, 3);

  vector<array<int, 3>> tmp(n1);
  cnt = 0;
  for (int i = -ws1; i <= ws1; ++i)
    for (int j = -ws1; j <= ws1; ++j)
      for (int k = -ws1; k <= ws1; ++k, ++cnt)
        tmp[cnt] = {{i, j, k}};
  std::sort(tmp.begin(), tmp.end(), sort_vector);

  vector<array<int, 3>> vidx1;
  for (int n = 0; n != n1; ++n) {
    array<int, 3> v = tmp[n];
    if (abs(v[0]) > ws || abs(v[1]) > ws || abs(v[2]) > ws)
      vidx1.push_back(v);
  }

  const int nvec = vidx1.size();
  vector<complex<mpreal>> lstar(msize);

  for (int ivec = 0; ivec != nvec; ++ivec) {
    array<int, 3> idx = vidx1[ivec];
    array<mpreal, 3> mvec;
    mvec[0] = idx[0] * primvecs[0][0] + idx[1] * primvecs[1][0] + idx[2] * primvecs[2][0];
    mvec[1] = idx[0] * primvecs[0][1] + idx[1] * primvecs[1][1] + idx[2] * primvecs[2][1];
    mvec[2] = idx[0] * primvecs[0][2] + idx[1] * primvecs[1][2] + idx[2] * primvecs[2][2];

    const mpreal rsq = mvec[0] * mvec[0] + mvec[1] * mvec[1] + mvec[2] * mvec[2];
    const mpreal r = sqrt(rsq);
    const mpreal ctheta = (r > zero) ? mvec[2]/r : zero;
    const mpreal phi = atan2(mvec[1], mvec[0]);

    for (int l = 0; l < max_rank; ++l) {
      for (int m = 0; m <= 2 * l; ++m) {
        const int am = abs(m - l);
        const int imul = l * l + m;

        mpreal plm_tilde = plm(l, am, ctheta) / pow(r, l+1);
        mpreal ft = one;
        for (int i = 1; i <= l - am; ++i) {
          plm_tilde *= ft;
          ft += one;
        }

        const mpreal sign = (m - l >=0) ? (cos(am * phi)) : (-1.0 * cos(am * phi));
        const mpreal real = sign * plm_tilde;
        const mpreal imag = sin(am * phi) * plm_tilde;

        lstar[imul] += complex<mpreal>(real, imag);
      }
    }
  }
  cout << "*** Lstar ***" << endl;
  for(int l = 0; l < max_rank; ++l)
    for(int m = 0; m <= l; ++m) {
      const int i = l * l + m + l;
      if (l % 2 == 0 && m % 4 == 0)
        cout << l << "  " << m << "       " << scientific << setprecision(32) << (lstar[i].real()).toDouble() << endl;
    }
  cout << endl;

  // Mlm(n)
  vector<pair<int, int>> lm_map;
  for (int l = 0; l < max_rank; ++l)
    for (int m = 0; m <= 2 * l; ++m)
      lm_map.push_back(make_pair(l, m-l));
  assert(lm_map.size() == msize);

  vector<complex<mpreal>> mlm = lstar; // iter 0
  const int max_iter = 15;
  for (int n = 1; n <= max_iter; ++n) {
    vector<complex<mpreal>> previous(msize);
    for (int i = 0; i != msize; ++i) {
      const int l = lm_map[i].first;
      previous[i] = mlm[i] / pow(3.0, l+1);
      mlm[i] = 0.0;
    }

    for (int l = 0; l < max_rank; ++l) {
      for (int m = 0; m <= 2*l; ++m) {
        const int im0 = l * l + m;

        for (int j = 0; j <= lmax - l; ++j) {
          for (int k = 0; k <= 2*j; ++k) {
            const int im1 = j * j + k;
            const int im = (l+j)*(l+j) + m + k;
            assert(l + j < max_rank);
            mlm[im0] += previous[im] * mstar[im1];
          }
        }
        mlm[im0] += lstar[im0];
      }
    }
  }


  cout << "*******************   Mlm ***********************" << endl;
  for(int l = 0; l < max_rank; ++l)
    for(int m = 0; m <= l; ++m) {
      const int i = l * l + m + l;
      if (l % 2 == 0 && m % 4 == 0)
        cout << l << "  " << m << "       " << scientific << setprecision(32) << (mlm[i].real()).toDouble() << endl;
    }

  return 0;
}
