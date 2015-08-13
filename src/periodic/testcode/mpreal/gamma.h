//
// Author  : Hai-Anh Le <anh@u.northwestern.edu>
// Date    : Aug 2015
//

#include <iostream>
#include <cmath>
#include <vector>
#include <src/util/constants.h>
#include <src/util/math/factorial.h>
#include "mpreal.h"
#include "gmp_macros.h"

using namespace mpfr;

mpreal compute_gamma(const int n) {

  mpfr::mpreal::set_default_prec(GMPPREC);
  assert(n > 0 && std::abs(n) % 2 == 1);
  const int l = (n-1)/2;

  const mpreal sqrtpi = GMPPISQRT;
  mpreal one = "1.0";
  mpreal two = "2.0";

  mpreal out = one;
  mpreal ft = one;
  for (int i = 1; i <= l; ++i) {
     out *= ft / two;
     ft += two;
  }

  return out * sqrtpi;
};

mpreal compute_gamma_upper(const int l, const mpreal x) {

  mpfr::mpreal::set_default_prec(GMPPREC);
  assert(l < 1200 && x < 300);
  const mpreal sqrtpi = GMPPISQRT;

  if (l > 1200 || x > 300)
    throw std::runtime_error("Failed to compute incomplete gamma function!");

  mpreal gamma = sqrtpi * erfc(sqrt(x));
  mpreal half = "0.5";
  mpreal one  = "1.0";
  for (int i = 1; i <= l; ++i) {
    const mpreal r = i + half;
    gamma = (r - one) * gamma + pow(x, r - one) * exp(-x);
  }

  return gamma;
};


mpreal compute_gamma_lower_scaled(const int l, const mpreal z, const mpreal beta) {
  mpfr::mpreal::set_default_prec(GMPPREC);
  mpreal sqrtpi = GMPPISQRT;
  mpreal half = "0.5";
  mpreal one  = "1.0";
  mpreal two  = "2.0";

  mpreal f;
  mpreal thresh = "300.0";

  if (z < thresh) {
    const int istart = 1200;
    mpreal gupper = compute_gamma_upper(istart, z);
    mpreal gamma = compute_gamma(2*istart+1);
    mpreal glower = gamma - gupper;
    mpreal prev = glower / (two * pow(z, istart + half));
    for (int i = istart-1; i <=l; --i) {
      f = (two * z * prev + exp(-z)) / (two * i + one);
      prev = f;
    }
  } else {
    f = sqrtpi / (pow(two, l+one) * pow(z, l+half));
    mpreal df = 1.0;
    for (int i = 1; i <= 2*l-1; i += 2) {
      f *= df;
      df += two;
    }
  }

  mpreal rsq = z / (beta * beta);
  mpreal out = 2.0 * pow(beta, 2*l+one) * pow(rsq, half*l) * f;

  return out;
};
