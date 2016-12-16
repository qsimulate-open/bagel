//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: slaterbatch.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <src/util/f77.h>
#include <src/util/constants.h>
#include <src/integral/rys/inline.h>
#include <src/integral/rys/sinline.h>
#include <src/integral/rys/slaterbatch.h>
#include <src/integral/rys/eribatch_base.h>
#include <src/integral/rys/coulombbatch_base.h>

#ifdef HAVE_LIBSLATER
#include <srootlist.h>

using namespace std;
using namespace bagel;

#define MIN_EXPONENT -200.0
#define GM1_THRESH 1.0e-15
#define GM1_THRESH_MAX 1.0e+15

static const double pitwohalf__ = pow(pi__, 2.5);

// TODO need to update
const static SRootList root;

SlaterBatch::SlaterBatch(const array<shared_ptr<const Shell>,4>& _info, const double max_density, const double gmm, const bool yukawa, shared_ptr<StackMem> stack)
:  RysIntegral<double>(_info, stack), gamma_(gmm), yukawa_(yukawa) {
  // ten-no == 1 means it is Slater/Yukawa int
  ++tenno_;

  const double integral_thresh = std::max(1.0e-30, PRIM_SCREEN_THRESH * 1.0e-10);// / max_density;

  set_swap_info(false);
  set_ab_cd();

  // set primsize_ and contsize_, as well as relevant members
  set_prim_contsizes();

  // set members for angular information
  int asize_final, csize_final, asize_final_sph, csize_final_sph;
  tie(asize_final, csize_final, asize_final_sph, csize_final_sph) = set_angular_info();

  allocate_data(asize_final, csize_final, asize_final_sph, csize_final_sph);

  allocate_arrays(primsize_);

  compute_ssss(integral_thresh);

  // obtain quadrature grid
  root_weight(primsize_);
}


SlaterBatch::~SlaterBatch() {
}


void SlaterBatch::root_weight(const int prim) {
  // determine the quadrature grid
  fill_n(weights_, primsize_, 0.0);
  if (rank_ == 1)
    root1_direct();
  else if (rank_ == 2)
    root2_direct();
  else {
    root.srootfunc_call(rank_, T_, U_, roots_, weights_, primsize_);
  }
}


void SlaterBatch::root1_direct() {
  const double prefac = sqrt(pi__) * 0.25;
  for (int i = 0; i != screening_size_; ++i) {
    const int ii = screening_[i];
    const double ct = T_[ii];
    const double cu = U_[ii];
    const double sqct = sqrt(ct);
    const double sqcu = sqrt(cu);
    const double lambda = sqct + sqcu;
    const double kappa = sqcu - sqct;
    double gm1, g0;

    if (fabs(ct) < 1.0e-10) {
      const double experfc_lambda = experfc(lambda);
      gm1 = sqrt(pi__)*0.5 / sqcu * experfc_lambda;
      g0  = 1 - (cu + cu) * gm1;
    } else {
      if (kappa < 0.0) {
        const double experfc_lambda = experfc(lambda);
        const double erfc_kappa = experfcm(kappa);
        const double factor1 = exp(-ct) * prefac;
        const double factor2 = exp(kappa * kappa - ct) * prefac;
        gm1 = (factor2 * erfc_kappa + factor1 * experfc_lambda) / sqcu;
        g0 =  (factor2 * erfc_kappa - factor1 * experfc_lambda) / sqct;
      } else {
        const double experfc_lambda = experfc(lambda);
        const double experfc_kappa  = experfc(kappa);
        const double factor = exp(-ct) * prefac;
        gm1 = factor / sqcu * (experfc_kappa + experfc_lambda);
        g0 =  factor / sqct * (experfc_kappa - experfc_lambda);
      }
    }

    if (gm1 > GM1_THRESH) {
      weights_[ii] = gm1;
      roots_[ii] = g0 / gm1;
    } else {
      weights_[ii] = 0.0;
      roots_[ii] = 0.0;
    }
  }
}


void SlaterBatch::root2_direct() {
  const double prefac = sqrt(pi__) * 0.25;
  for (int i = 0; i != screening_size_; ++i) {
    const int ii = screening_[i];
    const double ct = T_[ii];
    const double cu = U_[ii];
    const double sqct = sqrt(ct);
    const double sqcu = sqrt(cu);
    const double lambda = sqct + sqcu;
    const double kappa = sqcu - sqct;
    double gm1, g0, g1, g2;

    if (fabs(ct) < 1.0e-10) {
      const double experfc_lambda = experfc(lambda);
      const double cu2 = cu + cu;
      gm1 = sqrt(pi__)*0.5 / sqcu * experfc_lambda;
      g0  = 1.0 - cu2 * gm1;
      g1  = (1.0 - cu2 * g0) * 0.33333333333333333;
      g2  = (1.0 - cu2 * g1) * 0.2;
    } else {
      if (kappa < 0.0) {
        const double experfc_lambda = experfc(lambda);
        const double erfc_kappa = experfcm(kappa);
        const double expct = exp(-ct);
        const double factor1 = expct * prefac;
        const double factor2 = exp(kappa * kappa - ct) * prefac;
        gm1 = (factor2 * erfc_kappa + factor1 * experfc_lambda) / sqcu;
        g0 = (factor2 * erfc_kappa - factor1 * experfc_lambda) / sqct;
        const double cu2 = cu + cu;
        const double one2t = 0.5 / ct;
        g1 = (g0 + cu2 * gm1 - expct) * one2t;
        g2 = (3.0 * g1 + cu2 * g0 - expct) * one2t;
      } else {
        const double experfc_lambda = experfc(lambda);
        const double experfc_kappa  = experfc(kappa);
        const double expct = exp(-ct);
        const double factor = expct * prefac;
        gm1 = factor / sqcu * (experfc_kappa + experfc_lambda);
        g0 = factor / sqct * (experfc_kappa - experfc_lambda);
        const double cu2 = cu + cu;
        const double one2t = 0.5 / ct;
        g1 = (g0 + cu2 * gm1 - expct) * one2t;
        g2 = (3.0 * g1 + cu2 * g0 - expct) * one2t;
      }
    }

    if (gm1 < GM1_THRESH) {
      const int offset =  ii + ii;
      roots_[offset] = 0.0;
      roots_[offset + 1] = 0.0;
      weights_[offset] = 0.0;
      weights_[offset + 1] = 0.0;
    } else {
      double x[2];
      double w[1];
      x[0] = g0 / gm1;
      const double sigma11 = g1 - x[0] * g0;
      const double sigma12 = g2 - x[0] * g1;
      x[1] = sigma12 / sigma11 - x[0];
      w[0] = sigma11 / gm1;

      const double dd = x[0] + x[1];
      double g = x[0] - x[1];
      g = sqrt(g * g + 4.0 * w[0]);
      const double s = (dd - g) * 0.5;
      const double c = (dd + g) * 0.5;
      const double bb = w[0];
      const double p = x[0] - s;
      const double f = x[0] - c;
      const double lp1 = bb / (p * p + bb);
      const double lp2 = bb / (f * f + bb);

      const int offset =  ii + ii;
      roots_[offset] = s;
      roots_[offset + 1] = c;
      weights_[offset] = lp1 * gm1;
      weights_[offset + 1] = lp2 * gm1;
    }
  }
}


void SlaterBatch::compute_ssss(const double integral_thresh) {
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);

  const vector<double> exp0 = basisinfo_[0]->exponents();
  const vector<double> exp1 = basisinfo_[1]->exponents();
  const vector<double> exp2 = basisinfo_[2]->exponents();
  const vector<double> exp3 = basisinfo_[3]->exponents();

  vector<double> Ecd_save(prim2size_ * prim3size_);
  vector<double> qx_save(prim2size_ * prim3size_);
  vector<double> qy_save(prim2size_ * prim3size_);
  vector<double> qz_save(prim2size_ * prim3size_);

// TODO no primitive screening yet
#if 0
  const double minexp0 = *min_element(exp0.begin(), exp0.end());
  const double minexp1 = *min_element(exp1.begin(), exp1.end());
  const double minexp2 = *min_element(exp2.begin(), exp2.end());
  const double minexp3 = *min_element(exp3.begin(), exp3.end());
  const double min_ab = minexp0 * minexp1;
  const double min_cd = minexp2 * minexp3;
#endif

  indexpair23_.reserve(prim2size_ * prim3size_);

  const double r01_sq = AB_[0] * AB_[0] + AB_[1] * AB_[1] + AB_[2] * AB_[2];
  const double r23_sq = CD_[0] * CD_[0] + CD_[1] * CD_[1] + CD_[2] * CD_[2];

  {
#if 0
    const double cxp_min = minexp0 + minexp1;
    const double cxp_inv_min = 1.0 / cxp_min;
    const double min_Eab = exp(-r01_sq * min_ab * cxp_inv_min);
#endif
    int index23 = 0;
    for (auto expi2 = exp2.begin(); expi2 != exp2.end(); ++expi2) {
      for (auto expi3 = exp3.begin(); expi3 != exp3.end(); ++expi3, ++index23) {
        const double cxq = *expi2 + *expi3;
        const double cd = *expi2 * *expi3;
        const double cxq_inv = 1.0 / cxq;

        if (-r23_sq * (cd * cxq_inv) < MIN_EXPONENT) continue;
        Ecd_save[index23] = exp(-r23_sq * (cd * cxq_inv) );
        qx_save[index23] = (cx * *expi2 + dx * *expi3) * cxq_inv;
        qy_save[index23] = (cy * *expi2 + dy * *expi3) * cxq_inv;
        qz_save[index23] = (cz * *expi2 + dz * *expi3) * cxq_inv;

        indexpair23_.push_back(make_tuple(index23, *expi2, *expi3));
      }
    }
  }

  int index01 = 0;
  fill_n(coeff_,  primsize_, 0.0);
  fill_n(coeffy_, primsize_, 0.0);
  fill_n(T_, primsize_, -1.0);
  fill_n(U_, primsize_, 1.0e-100);

  // this should be allocated in RysInt
  screening_size_ = 0;

  const double twogamma = 2.0 / gamma_;
  for (auto expi0 = exp0.begin(); expi0 != exp0.end(); ++expi0) {
    for (auto expi1 = exp1.begin(); expi1 != exp1.end(); ++expi1, ++index01) {
      const double cxp = *expi0 + *expi1;
      const double ab = *expi0 * *expi1;
      const double cxp_inv = 1.0 / cxp;
      if (-r01_sq * (ab * cxp_inv) < MIN_EXPONENT) continue;
      const double Eab = exp(-r01_sq * (ab * cxp_inv) );
      const double coeff_half = 2 * Eab * pitwohalf__;
      const double px = (ax * *expi0 + bx * *expi1) * cxp_inv;
      const double py = (ay * *expi0 + by * *expi1) * cxp_inv;
      const double pz = (az * *expi0 + bz * *expi1) * cxp_inv;

      const int index_base = prim2size_ * prim3size_ * index01;

      for (auto expi23 = indexpair23_.begin(); expi23 != indexpair23_.end(); ++expi23) {
          const int index23 = get<0>(*expi23);
          const int index = index_base + index23;
          const double exp2value = get<1>(*expi23);
          const double exp3value = get<2>(*expi23);
          const double cxq = exp2value + exp3value;
          xp_[index] = cxp;
          xq_[index] = cxq;
          const double cxpxq = cxp * cxq;
          const double onepqp_q = 1.0 / (::sqrt(cxp + cxq) * cxpxq);
          coeffy_[index] = Ecd_save[index23] * coeff_half * onepqp_q;
          const double rho = cxpxq / (cxp + cxq);
          const double xpq = qx_save[index23] - px;
          const double ypq = qy_save[index23] - py;
          const double zpq = qz_save[index23] - pz;
          const double pq_sq = xpq * xpq + ypq * ypq + zpq * zpq;
          const double U = 0.25 * gamma_ * gamma_ / rho;
          if (U - gamma_ * sqrt(pq_sq) < MIN_EXPONENT) continue;
          if (coeffy_[index] * exp(U - gamma_ * sqrt(pq_sq)) < integral_thresh) continue;

          const double T = rho * pq_sq;
          coeff_[index] = coeffy_[index] * twogamma * U;
          const int index3 = index * 3;
          P_[index3] = px;
          P_[index3 + 1] = py;
          P_[index3 + 2] = pz;
          Q_[index3] =     qx_save[index23];
          Q_[index3 + 1] = qy_save[index23];
          Q_[index3 + 2] = qz_save[index23];
          T_[index] = T;
          U_[index] = U;
          screening_[screening_size_] = index;
          ++screening_size_;
      }
    }
  }
}

void SlaterBatch::allocate_data(const int asize_final, const int csize_final, const int asize_final_sph, const int csize_final_sph) {
  size_final_ = asize_final_sph * csize_final_sph * contsize_;

  const unsigned int size_start = asize_ * csize_ * primsize_;
  const unsigned int size_intermediate = asize_final * csize_ * contsize_;
  const unsigned int size_intermediate2 = asize_final_sph * csize_final * contsize_;
  size_block_ = std::max(size_start, std::max(size_intermediate, size_intermediate2));
  size_alloc_ = size_block_;

  stack_save_ = stack_->get(size_alloc_);
  stack_save2_ = nullptr;

  // if Slater/Yukawa integrals
  if (tenno_)
    stack_save2_ = stack_->get(size_alloc_);

  data_ = stack_save_;
  data2_ = stack_save2_;
}

#endif
