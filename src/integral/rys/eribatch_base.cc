//
// BAGEL - Parallel electron correlation program.
// Filename: eribatch_base.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/integral/rys/inline.h>
#include <src/integral/rys/eribatch_base.h>
#include <src/integral/rys/naibatch_base.h>
#include <src/integral/rys/erirootlist.h>
#include <src/integral/rys/breitrootlist.h>
#include <src/integral/rys/spin2rootlist.h>
// #include <src/integral/comprys/comperirootlist.h>
#include <src/integral/rys/breitbatch.h>
#include <src/util/constants.h>
#include <algorithm>
#include <cmath>

// same function for NAI below.
using namespace std;
using namespace bagel;


const static ERIRootList eri;
// const static ComplexERIRootList comperi;
const static BreitRootList br;
const static Spin2RootList s2;

static const double pitwohalf__ = pow(pi__, 2.5);
static const double pimhalf__ = 1.0/sqrt(pi__);
static const double T_thresh__ = 1.0e-8;

ERIBatch_base::ERIBatch_base(const array<shared_ptr<const Shell>,4>& o, const double max_density, const int deriv, const int breit, shared_ptr<StackMem> stack)
 : RysInt(o, stack) {

  breit_ = breit;

  deriv_rank_ = deriv;

  // determines if we want to swap shells
  set_swap_info(true);

  // stores AB and CD
  set_ab_cd();

  // set primsize_ and contsize_, as well as relevant members
  set_prim_contsizes();

  // sets angular info
  int asize_final, csize_final, asize_final_sph, csize_final_sph;
  tie(asize_final, csize_final, asize_final_sph, csize_final_sph) = set_angular_info();

  // allocate
  allocate_data(asize_final, csize_final, asize_final_sph, csize_final_sph);
  allocate_arrays(primsize_);

  const double integral_thresh = (max_density != 0.0) ? (PRIM_SCREEN_THRESH / max_density) : 0.0;
  compute_ssss(integral_thresh);

  root_weight(primsize_);

}


void ERIBatch_base::root_weight(const int ps) {
  if (breit_ == 0) {
    if (london_) {
//      comperi.root(rank_, T_, roots_, weights_, ps);   // This line goes on the other branch
    } else if (amax_ + cmax_ == 0) {
      for (int j = 0; j != screening_size_; ++j) {
        int i = screening_[j];
        if (T_[i] < T_thresh__) {
          weights_[i] = 1.0;
        } else {
          const double sqrtt = sqrt(T_[i]);
          const double erfsqt = inline_erf(sqrtt);
          weights_[i] = erfsqt * sqrt(pi__) * 0.5 / sqrtt;
        }
      }
    } else {
      eri.root(rank_, T_, roots_, weights_, ps);
    }
  } else if (breit_ == 1) {
    br.root(rank_, T_, roots_, weights_, ps);
  } else if (breit_ == 2) {
    s2.root(rank_, T_, roots_, weights_, ps);
  } else {
    assert(false);
  }
}

void NAIBatch_base::root_weight(const int ps) {
  if (amax_ + cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int i = screening_[j];
      if (T_[i] < T_thresh__) {
        weights_[i] = 1.0;
      } else {
        const double sqrtt = sqrt(T_[i]);
        const double erfsqt = inline_erf(sqrtt);
        weights_[i] = erfsqt * sqrt(pi__) * 0.5 / sqrtt;
      }
    }
  } else {
    eri.root(rank_, T_, roots_, weights_, ps);
  }
}


// a hack for screening of three-center integrals
static double rnd(const double& a) { return (a > 0.0) ? a : 1.0; };

// this sets T_ (+ U_), p_, q_, xp_, xq_, coeff_, and screening_size_
// for ERI evaulation. Other than that, we need to overload this function in a derived class
void ERIBatch_base::compute_ssss(const double integral_thresh) {
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

  const double* exp0 = (basisinfo_[0]->exponents_pointer());
  const double* exp1 = (basisinfo_[1]->exponents_pointer());
  const double* exp2 = (basisinfo_[2]->exponents_pointer());
  const double* exp3 = (basisinfo_[3]->exponents_pointer());
  const int nexp0 = basisinfo_[0]->num_primitive();
  const int nexp1 = basisinfo_[1]->num_primitive();
  const int nexp2 = basisinfo_[2]->num_primitive();
  const int nexp3 = basisinfo_[3]->num_primitive();

  double* const Ecd_save = stack_->get(prim2size_ * prim3size_);
  double* const qx_save = stack_->get(prim2size_ * prim3size_);
  double* const qy_save = stack_->get(prim2size_ * prim3size_);
  double* const qz_save = stack_->get(prim2size_ * prim3size_);

  const double minexp0 = *min_element(exp0, exp0+nexp0);
  const double minexp1 = *min_element(exp1, exp1+nexp1);
  const double minexp2 = *min_element(exp2, exp2+nexp2);
  const double minexp3 = *min_element(exp3, exp3+nexp3);
  const double min_ab = rnd(minexp0) * rnd(minexp1);
  const double min_cd = rnd(minexp2) * rnd(minexp3);
  const double min_abp = minexp0 * minexp1;
  const double min_cdp = minexp2 * minexp3;

  // minimum distance between two lines (AB and CD)
  const double x_ab_cd = AB_[1] * CD_[2] - AB_[2] * CD_[1];
  const double y_ab_cd = AB_[2] * CD_[0] - AB_[0] * CD_[2];
  const double z_ab_cd = AB_[0] * CD_[1] - AB_[1] * CD_[0];
  const double x_n_ac = x_ab_cd * (ax - cx);
  const double y_n_ac = y_ab_cd * (ay - cy);
  const double z_n_ac = z_ab_cd * (az - cz);
  const double innerproduct = x_n_ac + y_n_ac + z_n_ac;
  const double norm_ab_cd_sq = x_ab_cd * x_ab_cd + y_ab_cd * y_ab_cd + z_ab_cd * z_ab_cd;
  const double min_pq_sq = norm_ab_cd_sq == 0.0 ? 0.0 : innerproduct * innerproduct / norm_ab_cd_sq;

  const double r01_sq = AB_[0] * AB_[0] + AB_[1] * AB_[1] + AB_[2] * AB_[2];
  const double r23_sq = CD_[0] * CD_[0] + CD_[1] * CD_[1] + CD_[2] * CD_[2];

  unsigned int tuple_length = 0u;
  double* const tuple_field = stack_->get(nexp2*nexp3*3);
  int* tuple_index = (int*)(tuple_field+nexp2*nexp3*2);
  {
    const double cxp_min = minexp0 + minexp1;
    const double cxp_inv_min = 1.0 / cxp_min;
    const double min_Eab = exp(-r01_sq * min_abp * cxp_inv_min);
    int index23 = 0;
    for (const double* expi2 = exp2; expi2 != exp2+nexp2; ++expi2) {
      for (const double* expi3 = exp3; expi3 != exp3+nexp3; ++expi3, ++index23) {
        const double cxq = *expi2 + *expi3;
        const double cdp = *expi2 * *expi3;
        const double cd = rnd(*expi2) * rnd(*expi3);
        const double cxq_inv = 1.0 / cxq;
        Ecd_save[index23] = exp(-r23_sq * (cdp * cxq_inv) );
        qx_save[index23] = (cx * *expi2 + dx * *expi3) * cxq_inv;
        qy_save[index23] = (cy * *expi2 + dy * *expi3) * cxq_inv;
        qz_save[index23] = (cz * *expi2 + dz * *expi3) * cxq_inv;
        if (integral_thresh != 0.0) {
          const double rho = cxp_min * cxq / (cxp_min + cxq);
          const double T = rho * min_pq_sq;
          const double onepqp_q = 1.0 / (sqrt(cxp_min + cxq) * cxp_min * cxq);
          const double abcd_sc = min_ab * cd;
          const double abcd_sc_3 = abcd_sc * abcd_sc * abcd_sc;
          const double abcd_sc_3_4 = sqrt(sqrt(abcd_sc_3));
          const double tsqrt = sqrt(T);
          const double ssss = 16.0 * Ecd_save[index23] * min_Eab * abcd_sc_3_4 * onepqp_q
                              * (T > T_thresh__ ? inline_erf(tsqrt) * 0.5 / tsqrt : pimhalf__);
          if (ssss > integral_thresh) {
            tuple_field[tuple_length*2  ] = *expi2;
            tuple_field[tuple_length*2+1] = *expi3;
            tuple_index[tuple_length] = index23;
            ++tuple_length;
          }
        } else {
          tuple_field[tuple_length*2  ] = *expi2;
          tuple_field[tuple_length*2+1] = *expi3;
          tuple_index[tuple_length] = index23;
          ++tuple_length;
        }
      }
    }
  }

  int index = 0;
  int index01 = 0;
  fill_n(T_, primsize_, -1.0);
  screening_size_ = 0;

  const double cxq_min = minexp2 + minexp3;
  const double cxq_inv_min = 1.0 / cxq_min;
  const double min_Ecd = exp(-r23_sq * min_cdp * cxq_inv_min);
  for (const double* expi0 = exp0; expi0 != exp0+nexp0; ++expi0) {
    for (const double* expi1 = exp1; expi1 != exp1+nexp1; ++expi1, ++index01) {
      const double cxp = *expi0 + *expi1;
      const double abp = *expi0 * *expi1;
      const double ab = rnd(*expi0) * rnd(*expi1);
      const double cxp_inv = 1.0 / cxp;
      const double Eab = exp(-r01_sq * (abp * cxp_inv) );
      const double coeff_half = 2 * Eab * pitwohalf__;
      const double px = (ax * *expi0 + bx * *expi1) * cxp_inv;
      const double py = (ay * *expi0 + by * *expi1) * cxp_inv;
      const double pz = (az * *expi0 + bz * *expi1) * cxp_inv;
      if (integral_thresh != 0.0) {
        const double rho_sc = cxp * cxq_min / (cxp + cxq_min);
        const double T_sc = rho_sc * min_pq_sq;
        const double onepqp_q_sc = 1.0 / (sqrt(cxp + cxq_min) * cxp * cxq_min);
        const double tsqrt = sqrt(T_sc);
        const double abcd_sc = ab * min_cd;
        const double abcd_sc_3 = abcd_sc * abcd_sc * abcd_sc;
        const double abcd_sc_3_4 = sqrt(sqrt(abcd_sc_3));
        const double ssss = 16.0 * min_Ecd * Eab * abcd_sc_3_4 * onepqp_q_sc
                          * (T_sc > T_thresh__ ? inline_erf(tsqrt) * 0.5 / tsqrt : pimhalf__);
        if (ssss < integral_thresh) continue;
      }

      const int index_base = prim2size_ * prim3size_ * index01;

      for (unsigned int i = 0; i != tuple_length; ++i) {
          const int index23 = tuple_index[i];
          const int index = index_base + index23;
          const double exp2value = tuple_field[2*i];
          const double exp3value = tuple_field[2*i+1];
          const double cxq = exp2value + exp3value;
          xp_[index] = cxp;
          xq_[index] = cxq;
          const double cxpxq = cxp * cxq;
          const double onepqp_q = 1.0 / (sqrt(cxp + cxq) * cxpxq);
          coeff_[index] = Ecd_save[index23] * coeff_half * onepqp_q;
          const double rho = cxpxq / (cxp + cxq);
          const double xpq = qx_save[index23] - px;
          const double ypq = qy_save[index23] - py;
          const double zpq = qz_save[index23] - pz;
          const double T = rho * (xpq * xpq + ypq * ypq + zpq * zpq);
          const int index3 = index * 3;
          p_[index3] = px;
          p_[index3 + 1] = py;
          p_[index3 + 2] = pz;
          q_[index3] =     qx_save[index23];
          q_[index3 + 1] = qy_save[index23];
          q_[index3 + 2] = qz_save[index23];
          T_[index] = T;
          screening_[screening_size_] = index;
          ++screening_size_;
      }
    }
  }
  stack_->release(nexp2*nexp3*3, tuple_field);
  stack_->release(prim2size_*prim3size_, qz_save);
  stack_->release(prim2size_*prim3size_, qy_save);
  stack_->release(prim2size_*prim3size_, qx_save);
  stack_->release(prim2size_*prim3size_, Ecd_save);

}

// TODO this is not a good design. Should refactor at certain point using virtual functions...
void ERIBatch_base::allocate_data(const int asize_final, const int csize_final, const int asize_final_sph, const int csize_final_sph) {
  size_final_ = asize_final_sph * csize_final_sph * contsize_;
  if (deriv_rank_ == 0) {
    const unsigned int size_start = asize_ * csize_ * primsize_;
    const unsigned int size_intermediate = asize_final * csize_ * contsize_;
    const unsigned int size_intermediate2 = asize_final_sph * csize_final * contsize_;
    size_block_ = std::max(size_start, std::max(size_intermediate, size_intermediate2));
    size_alloc_ = size_block_;

    // if this is a two-electron Breit integral
    if (breit_)
      size_alloc_ = 6 * size_block_;

    stack_save_ = stack_->get(size_alloc_);
    stack_save2_ = nullptr;

    // if Slater/Yukawa integrals
    if (tenno_)
      stack_save2_ = stack_->get(size_alloc_);

  // derivative integrals
  } else if (deriv_rank_ == 1) {
    size_block_ = asize_final * csize_final * primsize_;
///*
    // if this is a two-electron gradient integral
    if (dynamic_cast<ERIBatch_base*>(this)) {
      size_alloc_ = 12 * size_block_;
    // if this is an NAI gradient integral
    } else if (dynamic_cast<NAIBatch_base*>(this)) {
      // in this case, we store everything
      size_alloc_ = (dynamic_cast<NAIBatch_base*>(this)->mol()->natom()) * 3.0 * size_block_;
      assert(csize_final == 1);
    } else {
      throw std::logic_error("something is strange in RysInt::allocate_data");
    }
//*/
    stack_save_ = stack_->get(size_alloc_);
    stack_save2_ = nullptr;
  }
  data_ = stack_save_;
  data2_ = stack_save2_;
}

