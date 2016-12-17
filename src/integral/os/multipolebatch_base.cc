//
// BAGEL - Parallel electron correlation program.
// Filename: multipolebatch_base.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/integral/os/multipolebatch_base.h>
#include <src/util/math/legendre.h>

using namespace std;
using namespace bagel;

const static Legendre plm;

static const double pisqrt__ = sqrt(pi__);

MultipoleBatch_base::MultipoleBatch_base(const array<shared_ptr<const Shell>,2>& sh, const array<double, 3> c,
                                         const int lmax, shared_ptr<StackMem> stack)
 : OSIntegral(sh, stack), centre_(c), lmax_(lmax) {

  num_multipoles_ = (lmax_ + 1) * (lmax_ + 1); // spherical
  common_init();
}


void MultipoleBatch_base::allocate_arrays(const size_t ps) {

  multipole_.resize(ps * num_multipoles_);
}


void MultipoleBatch_base::compute_ss(const double thr) {

  const vector<double> exp0 = basisinfo_[0]->exponents();
  const vector<double> exp1 = basisinfo_[1]->exponents();

  int iprim = 0;
  for (auto e0 = exp0.begin(); e0 != exp0.end(); ++e0) {
    for (auto e1 = exp1.begin(); e1 != exp1.end(); ++e1, ++iprim) {

      const double cxp = *e0 + *e1;
      const double cxp_inv = 1.0 / cxp;

      array<double, 3> P, PQ;
      P[0] = (basisinfo_[0]->position(0) * *e0 + basisinfo_[1]->position(0) * *e1) * cxp_inv;
      P[1] = (basisinfo_[0]->position(1) * *e0 + basisinfo_[1]->position(1) * *e1) * cxp_inv;
      P[2] = (basisinfo_[0]->position(2) * *e0 + basisinfo_[1]->position(2) * *e1) * cxp_inv;

      PQ[0] = P[0] - centre_[0];
      PQ[1] = P[1] - centre_[1];
      PQ[2] = P[2] - centre_[2];

      vector<complex<double>> olm = compute_OlmPQ(PQ);
      const double coeff = pisqrt__ * pi__ * sqrt(cxp_inv) * cxp_inv;

      for (int i = 0; i != num_multipoles_; ++i) {
        const double rABsq = AB_[0] * AB_[0] + AB_[1] * AB_[1] + AB_[2] * AB_[2];
        const int index = iprim + i * prim0_ * prim1_;
        const double Sab = coeff * exp(- *e0 * *e1 * cxp_inv * rABsq);
        multipole_[index] = Sab * olm[i];
        if (swap01_)
          multipole_[index] = conj(multipole_[index]);
      }
    }
  }
}


vector<complex<double>> MultipoleBatch_base::compute_OlmPQ(const array<double, 3>& PQ) {

  const double r = sqrt(PQ[0]*PQ[0] + PQ[1]*PQ[1] + PQ[2]*PQ[2]);
  const double ctheta = (r > numerical_zero__) ? PQ[2]/r : 0.0;
  const double phi = atan2(PQ[1], PQ[0]);

  vector<complex<double>> out(num_multipoles_);

  out[0] = 1.0;
  for (int l = 1; l <= lmax_; ++l) {
    for (int mm = 0; mm <= 2 * l; ++mm) {
      const int m = mm - l;
      const int am = abs(m);

      double coeff = pow(r, l) * plm.compute(l, am, ctheta);
      double ft = 1.0;
      for (int i = 1; i <= l + am; ++i) {
        coeff /= ft;
        ft += 1.0;
      }

      const double real = (m >= 0) ? (coeff * cos(am * phi)) : (pow(-1.0, m) * coeff * cos(am * phi));
      const double imag = (m >= 0) ? (coeff * sin(am * phi)) : (pow(-1.0, m+1) * coeff * sin(am * phi));
      out[l*l+mm] = complex<double>(real, imag);
    }
  }

  return out;
}
