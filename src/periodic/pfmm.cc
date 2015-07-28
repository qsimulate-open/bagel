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
  num_multipoles_ = (lmax_ + 1) * (lmax_ + 1);
  mlm_.resize(num_multipoles_);
}


bool PFMM::is_in_cff(array<double, 3> L) {

  const double extent = scell_->extent();

  const double rsq = L[0]*L[0] + L[1]*L[1] + L[2]*L[2];
  const bool out = (rsq > 2.0 * (1 + ws_) *  extent) ? true : false;

  return out;
}


void PFMM::compute_mlm(vector<array<double, 3>> rvec) {

  // real term
  for (auto& v : rvec) {
    const double rsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    const double ctheta = (rsq > numerical_zero__) ? v[2]/sqrt(rsq) : 0.0;
    const double phi = atan2(v[1], v[0]);
    const double b2r2 = beta__ * beta__ * rsq;

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
  }
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


vector<const array<double, 3>> PFMM::get_intermediate_layer(const int imax) {

  int cnt = 0;
  for (auto& L : lattice_vectors_) {

  }
}


shared_ptr<PData> PFMM::compute_Jop(shared_ptr<const PData> density) {

}
#endif
