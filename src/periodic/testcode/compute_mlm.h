//
// Author  : Hai-Anh Le <anh@u.northwestern.edu>
// Date    : July 2015
//


#include <boost/math/constants/constants.hpp>
#include <src/util/math/legendre.h>
#include "gmp_macros.h"
#include "mpreal.h"
#include "gamma.h"

using namespace std;
using namespace bagel;
using namespace mpfr;

const static Legendre plm;

mpreal dot(array<mpreal, 3> b, array<mpreal, 3> c) { return b[0] * c[0] + b[1] * c[1] + b[2] * c[2]; }


array<mpreal, 3> cross(array<mpreal, 3> b, array<mpreal, 3> c, mpreal s) {

  array<mpreal, 3> out;
  out[0] = (b[1] * c[2] - b[2] * c[1]) * s;
  out[1] = (b[2] * c[0] - b[0] * c[2]) * s;
  out[2] = (b[0] * c[1] - b[1] * c[0]) * s;

  return out;
};


bool is_in_cff(const int ws, const int n0, const int n1, const int n2) {

  const bool out = (abs(n0) > ws && abs(n1) > ws && abs(n2) > ws) ? true : false;
  return out;
};


#if 1
vector<complex<mpreal>> compute_mlm(const int ws, const int lmax, const int limit, const mpreal thresh) {

  mpfr::mpreal::set_default_prec(GMPPREC);
  const mpreal pi = GMPPI;
  const mpreal beta = GMPPISQRT; // convergence parameter
  const mpreal zero = "0.0";
  const mpreal one  = "1.0";
  const mpreal half = "0.5";
  const mpreal two  = "2.0";
  const mpreal mone = "-1.0";

  const size_t ndim = 3;
  const size_t num_multipoles = (lmax + 1) * (lmax + 1);
  vector<complex<mpreal>> out(num_multipoles);

  mpreal pibeta = (pi / beta) * (pi / beta);

  // generate lattice vectors - 3D for now
  vector<array<mpreal, 3>> primitive_vectors(3);
  primitive_vectors[0] = {{one, zero, zero}};
  primitive_vectors[1] = {{zero, one, zero}};
  primitive_vectors[2] = {{zero, zero, one}};

  const int nvec = std::pow(2*limit+1, ndim);
  vector<array<mpreal, 3>> rvec(nvec);
  vector<array<mpreal, 3>> kvec(nvec);
  vector<array<int, 3>> vindex(nvec);

  array<mpreal, 3> v;
  int cnt = 0;
  for (int n3 = -limit; n3 <= limit; ++n3) {
    for (int n2 = -limit; n2 <= limit; ++n2) {
      for (int n1 = -limit; n1 <= limit; ++n1, ++cnt) {
        v[0] = n1 * primitive_vectors[0][0] + n2 * primitive_vectors[1][0] + n3 * primitive_vectors[2][0];
        v[1] = n1 * primitive_vectors[0][1] + n2 * primitive_vectors[1][1] + n3 * primitive_vectors[2][1];
        v[2] = n1 * primitive_vectors[0][2] + n2 * primitive_vectors[1][2] + n3 * primitive_vectors[2][2];
        vindex[cnt] = {{n1, n2, n3}};
        rvec[cnt] = v;
      }
    }
  }

  int ivec = 0;
  for (auto& v : rvec) {
    array<int, 3> id = vindex[ivec];

    const mpreal rsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    const mpreal r = sqrt(rsq);
    const mpreal ctheta = (rsq > static_cast<mpreal>(numerical_zero__)) ? v[2]/r : zero;
    const mpreal phi = atan2(v[1], v[0]);
    const mpreal b2r2 = beta * beta * rsq;

    int cnt = 0;
    for (int l = 0; l <= lmax; ++l) {
      const mpreal lhalf = l + half;
      const mpreal gamma = compute_gamma(2*l+1);
      const mpreal glower = compute_gamma_lower_scaled(l, b2r2, beta);
      const mpreal gupper = gamma / pow(r, l+one) - gupper;

      for (int mm = 0; mm <= 2 * l; ++mm, ++cnt) {
        const int m = mm - l;
        const int am = abs(m);

        mpreal plm_tilde = static_cast<mpreal>(plm.compute(l, abs(m), ctheta.toDouble()));
        mpreal ft = one;
        for (int i = 1; i <= l - abs(m); ++i) {
          plm_tilde *= ft;
          ft += one;
        }
        const mpreal sign = (m >=0) ? (cos(am * phi)) : (mone * cos(am * phi));

        if (is_in_cff(ws, id[0], id[1], id[2])) {
          // real term
          const mpreal real = gupper * static_cast<mpreal>(sign) * plm_tilde / gamma;
          const mpreal imag = gupper * static_cast<mpreal>(sin(am * phi)) * plm_tilde;
          out[cnt] += complex<mpreal>(real, imag);
        }
      }
    }

    ++ivec;
  }

  array<mpreal, 3> a23 = cross(primitive_vectors[1], primitive_vectors[2], 1.0);
  const mpreal scale = one / dot(primitive_vectors[0], a23);
  //const mpreal scale = 2.0 * pi__ / dot(primitive_vectors[0], a23);
  vector<array<mpreal, 3>> primitive_kvectors(3);
  primitive_kvectors[0] = cross(primitive_vectors[1], primitive_vectors[2], scale);
  primitive_kvectors[1] = cross(primitive_vectors[2], primitive_vectors[0], scale);
  primitive_kvectors[2] = cross(primitive_vectors[0], primitive_vectors[1], scale);

  cnt = 0;
  for (int n3 = -limit; n3 <= limit; ++n3) {
    for (int n2 = -limit; n2 <= limit; ++n2) {
      for (int n1 = -limit; n1 <= limit; ++n1, ++cnt) {
        v[0] = n1 * primitive_kvectors[0][0] + n2 * primitive_kvectors[1][0] + n3 * primitive_kvectors[2][0];
        v[1] = n1 * primitive_kvectors[0][1] + n2 * primitive_kvectors[1][1] + n3 * primitive_kvectors[2][1];
        v[2] = n1 * primitive_kvectors[0][2] + n2 * primitive_kvectors[1][2] + n3 * primitive_kvectors[2][2];
        kvec[cnt] = v;
      }
    }
  }

  ivec = 0;
  for (auto& v : kvec) {
    array<int, 3> id = vindex[ivec];

    const mpreal rsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    const mpreal r = sqrt(rsq);
    const mpreal ctheta = (rsq > static_cast<mpreal>(numerical_zero__)) ? v[2]/r : zero;
    const mpreal phi = atan2(v[1], v[0]);
    const mpreal b2r2 = beta * beta * rsq;

    int cnt = 0;
    for (int l = 0; l <= lmax; ++l) {
      const complex<mpreal> coeffl = static_cast<complex<mpreal>>(std::pow(complex<double>(0.0, 1.0), l)) * pow(pi, l-half);
      const mpreal lhalf = l + half;
      const mpreal gamma  = compute_gamma(2*l+1);
      const mpreal glower = compute_gamma_lower_scaled(l, b2r2, beta);
      for (int mm = 0; mm <= 2 * l; ++mm, ++cnt) {
        const int m = mm - l;
        const int am = abs(m);

        mpreal plm_tilde = static_cast<mpreal>(plm.compute(l, abs(m), ctheta.toDouble()));
        mpreal ft = one;
        for (int i = 1; i <= l - abs(m); ++i) {
          plm_tilde *= ft;
          ft += one;
        }

        const mpreal sign = (m >=0) ? (cos(am * phi)) : (mone * cos(am * phi));

        if (!is_in_cff(ws, id[0], id[1], id[2])) {
          // substract smooth part within ws
          const mpreal real = glower * static_cast<mpreal>(sign) * plm_tilde / gamma;
          const mpreal imag = glower * static_cast<mpreal>(sin(am * phi)) * plm_tilde;
          out[cnt] -= complex<mpreal>(real, imag);
        }
        // smooth term
        const mpreal coeffm = plm_tilde * pow(r, l-2) *  exp(-rsq * pibeta);
        mpreal real = coeffm * static_cast<mpreal>(sign) / gamma;
        mpreal imag = coeffm * static_cast<mpreal>(sin(am * phi));
        out[cnt] += coeffl * complex<mpreal>(real, imag);
      }
    }

    ++ivec;
  }

  return out;
};
#endif
