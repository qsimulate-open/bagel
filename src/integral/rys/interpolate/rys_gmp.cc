//
// Author: Toru Shiozaki
// Date  : April 2009
// Last Updated: May 2009
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
#include <gsl/gsl_sf_dawson.h>

using namespace std;
using namespace mpfr;

void rysroot_gmp(const vector<mpreal>& ta, vector<mpreal>& dx, vector<mpreal>& dw, const int nrank, const int nbatch) {

  mpfr::mpreal::set_default_prec(GMPPREC);
  mpreal mlp[150];
  mpreal sigma[150];
  mpreal fm[150];
  mpreal w[150];
  mpreal x[150];
  mpreal lp[150];

  for (int ibatch = 0; ibatch != nbatch; ++ibatch) {
    const int offset = ibatch * nrank;

    if (ta[ibatch] < 0) continue;

    const mpreal T = ta[ibatch];
    assert(T > 0);

    mpreal mone;

    {
      const mpreal zero = "0.0";
      const mpreal half = "0.5";
      const mpreal quater = "0.25";

      const mpreal sqrtt = sqrt(T);
      // target Fm(T)
      const mpreal pi = GMPPI;
      const mpreal halfpT = half / T;
#ifndef DAWSON
      fm[0] = sqrt(pi) / sqrtt * half * erf(sqrtt);
      for (int i = 1; i != 150; ++i) {
        fm[i] = halfpT * ((2*i-1) * fm[i - 1] - exp(-T));
      }
#ifdef BREIT
      // in case of Breit, we shift by one
      for (int i = 0; i != 39; ++i)
        fm[i] = fm[i+1];
#endif
#ifdef SPIN2
      // in case of Spin-Spin, we shift by two
      for (int i = 0; i != 38; ++i)
        fm[i] = fm[i+2];
#endif
#else
      const mpreal one = "1.0";
      double dawson_double = gsl_sf_dawson(sqrtt.toDouble());
      // based on taylor expansion
      for (int i = 0; i != 150; ++i) {
        fm[i] = 0;
        mpreal current = one/(2*i+1) * exp(-T);
        int j = 0;
        do {
          fm[i] += current;
          ++j;
          current *= mpreal(2*i+2*j-1)/mpreal(2*i+2*j+1) * T / j;
        } while (abs(current) > machine_epsilon(fm[i]));
      }
      assert(abs(gsl_sf_dawson(sqrtt.toDouble())/sqrtt.toDouble() - fm[0].toDouble()) < 1.0e-12);
#endif
      mone = fm[0];
    }


    {
      // Chebyshev algorithm; VERY unstable
      const mpreal mpone = "1.0";
      const int n = nrank;

      mlp[0] = mpone / fm[0];
      x[0] = fm[1] * mlp[0];
      w[0] = "0.0";

      for (int k = 0; k <= n - 2; k += 2) {
        for (int l = k; l <= 2 * n - k - 3; ++l) {
          sigma[l + 1] = fm[l + 2] - x[k] * fm[l + 1] - w[k] * sigma[l + 1] ;
        }
        mlp[k + 1] = mpone / sigma[k + 1];
        x[k + 1] = - fm[k + 1] * mlp[k] + sigma[k + 2] * mlp[k + 1];
        w[k + 1] = sigma[k + 1] * mlp[k];

        if(k != n - 2) {
          for (int l = k + 1; l <= 2 * n - k - 4; ++l)
            fm[l + 1] = sigma[l + 2] - x[k + 1] * sigma[l + 1] - w[k + 1] * fm[l + 1];
          mlp[k + 2] = mpone / fm[k + 2];
          x[k + 2] = - sigma[k + 2] * mlp[k + 1] + fm[k + 3] * mlp[k + 2];
          w[k + 2] = fm[k + 2] * mlp[k + 1];
        }
      }

      dx[offset + 0] = x[0];
      for (int i = 1; i != n; ++i) {
        dw[offset + i - 1] = sqrt(w[i]);
        dx[offset + i] = x[i];
      }
      dw[offset + n - 1] = "0.0";

      // solve tri-diagonal linear equation
      const mpreal zero = "0.0";
      const mpreal one = "1.0";

      lp[0] = one;
      for (int i = 1; i <= n + n - 2; ++i) lp[i] = zero;

      for (int l = 0; l <= n - 1; ++l) {
        int iter = 0;
line1:
        int mm;
        for (mm = l; mm <= n - 2; ++mm) {
          mpreal dd = abs(dx[offset + mm]) + abs(dx[offset + mm + 1]);
          if (fabs(dw[offset + mm]) + dd == dd) goto line2;
        }
        mm = n - 1;
line2:
        if(mm != l) {
          ++iter;
          if (iter == 150) throw logic_error("bad");
          mpreal g = (dx[offset + l + 1] - dx[offset + l]) / (dw[offset + l] * 2);
          const mpreal r = sqrt(g * g + one);
          g = dx[offset + mm] - dx[offset + l] + dw[offset + l] / (g + (g >= 0 ? fabs(r) : -fabs(r)));
          mpreal s = one;
          mpreal c = one;
          mpreal p = zero;
          for (int i = mm - 1; i >= l; --i) {
            mpreal f = s * dw[offset + i];
            mpreal bb = c * dw[offset + i];
            mpreal r = sqrt(f * f + g * g);
            dw[offset + i + 1] = r;
            if(r == zero) {
              dx[offset + i + 1] = dx[offset + i + 1] - p;
              dw[offset + mm] = zero;
              goto line1;
            }
            const mpreal pr = one / r;
            s = f * pr;
            c = g * pr;
            g = dx[offset + i + 1] - p;
            mpreal td = c * bb;
            r = (dx[offset + i] - g) * s + td * 2;
            p = s * r;
            dx[offset + i + 1] = g + p;
            g = c * r - bb;
            f= lp[i + 1];
            lp[i + 1] = s * lp[i] + c * f;
            lp[i] = c * lp[i] - s * f;
          }
          dx[offset + l] = dx[offset + l] - p;
          dw[offset + l] = g;
          dw[offset + mm] = zero;
          goto line1;
        }
      }
      for (int i = 0; i != n; ++i) dw[offset + i] = lp[i] * lp[i] * mone;
    }
  }
}

