//
// Newint - Parallel electron correlation program.
// Filename: rys_gmp.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

//

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <src/rysint/mpreal.h>
#include <src/rysint/gmp_macros.h>
#include <src/rysint/eribatch.h>
#include <algorithm>

using namespace std;
using namespace mpfr;

void ERIBatch::rysroot_gmp(const double* ta, double* dx, double* dw, const int nrank, const int nbatch) {

  mpfr::mpreal::set_default_prec(GMPPREC);
  mpreal mlp[40];
  mpreal sigma[40];
  mpreal g_tu[40];
  mpreal w[40]; 
  mpreal x[40]; 
  double lp[40];

  for (int ibatch = 0; ibatch != nbatch; ++ibatch) {

    const int offset = ibatch * nrank;
    if (ta[ibatch] < 0.0) continue;

    mpreal T = ta[ibatch];
    if (ta[ibatch] == 0.0) { const mpreal ttt = "0.0000000000000001"; T = ttt;}
    
    const mpreal U = "0.0";

    double mone;
  
    {
      const mpreal zero = "0.0";
      const mpreal half = "0.5";
      const mpreal quater = "0.25";

      // some parameter
      const mpreal sqrtt = sqrt(T);
      const mpreal sqrtu = sqrt(U);
      const mpreal kappa  = - sqrtt + sqrtu;
      const mpreal lambda = sqrtt + sqrtu;

      // target Gm(T, U)
      const mpreal expmt = exp(-T);
      const mpreal prefactor = expmt * quater * GMPPISQRT; 
      const mpreal expkk = exp(kappa * kappa);
      const mpreal expll = exp(lambda * lambda);
      const mpreal erfck = erfc(kappa) * expkk;
      const mpreal erfcl = erfc(lambda) * expll;
      const mpreal halfpT = half / T;
      const mpreal twoU = U + U;
      mpreal g_tu_1 = U == zero ? zero : prefactor / sqrtu * (erfck + erfcl); 
      g_tu[0] = prefactor / sqrtt * (erfck - erfcl); 
      mone = g_tu[0];
      g_tu[1] = halfpT * (g_tu[0] + twoU * g_tu_1 - expmt);
      const int gtuend = 2 * nrank + 1;
      for (int i = 2; i != gtuend; ++i) {
        g_tu[i] = halfpT * (static_cast<mpreal>(2 * i - 1) * g_tu[i - 1] + twoU * g_tu[i - 2] - expmt);
      }
    }  


    {
      // Chebyshev algorithm; VERY unstable
      const mpreal mpone = "1.0";
      const int n = nrank; 

      mlp[0] = mpone / g_tu[0];
      x[0] = g_tu[1] * mlp[0];
      w[0] = "0.0";
 
      for (int k = 0; k <= n - 2; k += 2) {
        for (int l = k; l <= 2 * n - k - 3; ++l) {
          sigma[l + 1] = g_tu[l + 2] - x[k] * g_tu[l + 1] - w[k] * sigma[l + 1] ;
        }
        mlp[k + 1] = mpone / sigma[k + 1];
        x[k + 1] = - g_tu[k + 1] * mlp[k] + sigma[k + 2] * mlp[k + 1];
        w[k + 1] = sigma[k + 1] * mlp[k];

        if(k != n - 2) { 
          for (int l = k + 1; l <= 2 * n - k - 2; ++l)  
            g_tu[l + 1] = sigma[l + 2] - x[k + 1] * sigma[l + 1] - w[k + 1] * g_tu[l + 1];
          mlp[k + 2] = mpone / g_tu[k + 2];
          x[k + 2] = - sigma[k + 2] * mlp[k + 1] + g_tu[k + 3] * mlp[k + 2];
          w[k + 2] = g_tu[k + 2] * mlp[k + 1];
        }
      }

      // switching to double as the diagonalization is stable

      double* cdx = &dx[offset];
      double* cdw = &dw[offset];
      cdx[0] = x[0];
      for (int i = 1; i != n; ++i) {
        cdw[i - 1] = std::sqrt(static_cast<double>(w[i]));
        cdx[i] = x[i]; 
      }
      cdw[n - 1] = 0.0;

      // solve tri-diagonal linear equation 
      const double zero = 0.0;
      const double one = 1.0;
    
      lp[0] = one;
      for (int i = 1; i <= n + n - 2; ++i) lp[i] = zero;

      for (int l = 0; l <= n - 1; ++l) {
        int iter = 0;
line1:
        int mm;
        for (mm = l; mm <= n - 2; ++mm) { 
          double dd = fabs(cdx[mm]) + fabs(cdx[mm + 1]);
          if (fabs(cdw[mm]) + dd == dd) goto line2;
        }
        mm = n - 1;
line2:
        if(mm != l) {
          ++iter;
          double g = (cdx[l + 1] - cdx[l]) / (cdw[l] * 2);
          const double r = std::sqrt(g * g + one);
          g = cdx[mm] - cdx[l] + cdw[l] / (g + (g >= 0 ? fabs(r) : -fabs(r)));
          double s = one;
          double c = one;
          double p = zero;
          for (int i = mm - 1; i >= l; --i) {
            double f = s * cdw[i];
            double bb = c * cdw[i];
            double r = std::sqrt(f * f + g * g);
            cdw[i + 1] = r;
            if(r == zero) {
              cdx[i + 1] = cdx[i + 1] - p;
              cdw[mm] = zero;
              goto line1;
            }
            const double pr = 1.0 / r;
            s = f * pr;
            c = g * pr;
            g = cdx[i + 1] - p;
            double td = c * bb;
            r = (cdx[i] - g) * s + td * 2;
            p = s * r;
            cdx[i + 1] = g + p;     
            g = c * r - bb;
            f= lp[i + 1];
            lp[i + 1] = s * lp[i] + c * f;
            lp[i] = c * lp[i] - s * f;
          }
          cdx[l] = cdx[l] - p;
          cdw[l] = g;
          cdw[mm] = zero;
          goto line1;
        }
      }
      for (int i = 0; i != n; ++i) cdw[i] = lp[i] * lp[i] * mone;
    }
  }
}

