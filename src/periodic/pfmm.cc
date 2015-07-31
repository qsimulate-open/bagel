//
// BAGEL - Parallel electron correlation program.
// Filename: pfmm.cc
// Copyright (C) 2015 Toru Shiozaki
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


#include <boost/math/special_functions/gamma.hpp>
#include <src/util/math/legendre.h>
#include <src/periodic/pfmm.h>

using namespace std;
using namespace bagel;

const static Legendre plm;
const static Factorial f;

const static double beta__ = sqrt(pi__); // convergence parameter

PFMM::PFMM(shared_ptr<const SimulationCell> scell, const int lmax, const int ws, const double thresh)
  : scell_(scell), lmax_(lmax), ws_(ws) {

  assert(lmax_ <= ANG_HRR_END);
  ndim_ = scell->ndim();
  num_multipoles_ = (lmax_ + 1) * (lmax_ + 1);
  mlm_.resize(num_multipoles_);
}


bool PFMM::is_in_cff(array<double, 3> L) {

  const double extent = scell_->extent();

  const double rsq = L[0]*L[0] + L[1]*L[1] + L[2]*L[2];
  const bool out = (rsq > 2.0 * (1 + ws_) *  extent) ? true : false;

  return out;
}


void PFMM::compute_mlm(const int limit, const double thresh) {

  const double pibeta = pi__ * pi__ / (beta__ * beta__);

  array<double, 3> v = {{0.0, 0.0, 0.0}};
  for (int n0 = -limit; n0 != limit; ++n0) {
    array<double, 3> prim0 = scell_->primitive_vectors(0);
    v[0] += n0 * prim0[0];
    v[1] += n0 * prim0[1];
    v[2] += n0 * prim0[2];
    for (int n1 = -limit; n1 != limit; ++n1) {
      if (ndim_ > 1) {
        array<double, 3> prim1 = scell_->primitive_vectors(1);
        v[0] += n1 * prim1[0];
        v[1] += n1 * prim1[1];
        v[2] += n1 * prim1[2];
      }
      for (int n2 = -limit; n2 != limit; ++n2) {
        if (ndim_ > 2) {
          array<double, 3> prim2 = scell_->primitive_vectors(2);
          v[0] += n2 * prim2[0];
          v[1] += n2 * prim2[1];
          v[2] += n2 * prim2[2];
        }

        const double rsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        const double ctheta = (rsq > numerical_zero__) ? v[2]/sqrt(rsq) : 0.0;
        const double phi = atan2(v[1], v[0]);
        const double b2r2 = beta__ * beta__ * rsq;

        if (is_in_cff(v)) {
          // real term
          int count = 1;
          mlm_[0] += 1.0;
          for (int l = 1; l <= lmax_; ++l) {
            for (int mm = 0; mm <= 2 * l; ++mm, ++count) {
              const int m = mm - l;
              const int am = abs(m);

              double coeff = plm.compute(l, abs(m), ctheta) / pow(rsq, (l + 1)/2) * boost::math::gamma_q(l+0.5, b2r2);
              double ft = 1.0;
              for (int i = 1; i <= l - abs(m); ++i) {
                coeff *= ft;
                ++ft;
              }

              const double real = (m >=0) ? (coeff * cos(am * phi)) : (-1.0 * coeff * cos(am * phi));
              const double imag = coeff * sin(am * phi);
              mlm_[count] += complex<double>(real, imag);

            }
          }
        } else {
          // substract smooth part within ws
          int count = 1;
          for (int l = 1; l <= lmax_; ++l) {
            for (int mm = 0; mm <= 2 * l; ++mm, ++count) {
              const int m = mm - l;
              const int am = abs(m);

              double coeff = plm.compute(l, abs(m), ctheta) / pow(rsq, (l + 1)/2) * boost::math::gamma_p(l+0.5, b2r2);
              double ft = 1.0;
              for (int i = 1; i <= l - abs(m); ++i) {
                coeff *= ft;
                ++ft;
              }

              const double real = (m >=0) ? (coeff * cos(am * phi)) : (-1.0 * coeff * cos(am * phi));
              const double imag = coeff * sin(am * phi);
              mlm_[count] -= complex<double>(real, imag);
            }
          }
        }

        // smooth term
        int count = 1;
        for (int l = 1; l <= lmax_; ++l) {
          const complex<double> coeffl = pow(complex<double>(0.0, 1.0), l) * pow(pi__, l-0.5) / boost::math::tgamma(l+0.5);
          for (int mm = 0; mm <= 2 * l; ++mm, ++count) {
            const int m = mm - l;
            const int am = abs(m);

            double coeffm = plm.compute(l, abs(m), ctheta) * pow(rsq, (l - 2)/2) * exp(-rsq * pibeta);
            double ft = 1.0;
            for (int i = 1; i <= l - abs(m); ++i) {
              coeffm *= ft;
              ++ft;
            }

            const double real = (m >=0) ? (coeffm * cos(am * phi)) : (-1.0 * coeffm * cos(am * phi));
            const double imag = coeffm * sin(am * phi);
            mlm_[count] += coeffl * complex<double>(real, imag);
          }
        }

      } //n2
    } //n1
  } //n0
}


#if 0
void PFMM::compute_Sn(const double thresh, const int max_iter) { // S(n+1) = U_M[S(n)] O* + M*

  for (int iter = 0; iter != max_iter; ++i) {

    const double error = 0.0;
    if (error < thresh) {
      cout << "  * Sn converged." << endl << endl;
      break;
    } else if (iter == max_iter-1) {
      cout << "  * Max iteration reached when in compute_Sn." << endl << endl;
      break;
    }
  }
}


shared_ptr<PData> PFMM::compute_Jop(shared_ptr<const PData> density) {

}
#endif
