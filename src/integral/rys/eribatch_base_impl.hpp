//
// BAGEL - Parallel electron correlation program.
// Filename: rysintegral_impl.hpp
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Ryan Reynolds <ryanreynolds2018@u.northwestern.edu>
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

#ifdef ERIBATCH_BASE_HEADERS

#ifndef __SRC_INTEGRAL_RYS_ERIBATCH_BASE_HPP
#define __SRC_INTEGRAL_RYS_ERIBATCH_BASE_HPP

namespace bagel {

// this sets T_ (+ U_), P_, Q_, xp_, xq_, coeff_, and screening_size_
// for ERI evaulation. Other than that, we need to overload this function in a derived class
template <typename DataType, Int_t IntType>
void ERIBatch_Base<DataType, IntType>::compute_ssss(const double integral_thresh) {
  static_assert(IntType != Int_t::London || std::is_same<DataType, std::complex<double>>::value, "London-orbital integrals should be complex");

  // set atomic coordinates
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

  // set Gaussian exponents
  const double* exp0 = basisinfo_[0]->exponents_pointer();
  const double* exp1 = basisinfo_[1]->exponents_pointer();
  const double* exp2 = basisinfo_[2]->exponents_pointer();
  const double* exp3 = basisinfo_[3]->exponents_pointer();

  // set number of primitive Gaussians for each center
  const int nexp0 = basisinfo_[0]->num_primitive();
  const int nexp1 = basisinfo_[1]->num_primitive();
  const int nexp2 = basisinfo_[2]->num_primitive();
  const int nexp3 = basisinfo_[3]->num_primitive();

  // allocate temporary Q storage
  double* const Ecd_save = stack_->template get<double>(prim2size_*prim3size_);
  DataType* const qx_save = stack_->template get<DataType>(prim2size_*prim3size_);
  DataType* const qy_save = stack_->template get<DataType>(prim2size_*prim3size_);
  DataType* const qz_save = stack_->template get<DataType>(prim2size_*prim3size_);
  DataType* factor_cd_save;
  if (IntType == Int_t::London) {
    factor_cd_save = stack_->template get<DataType>(prim2size_*prim3size_);
  }

  // determine smallest a, b, c, d, p, q
  const double minexp0 = *std::min_element(exp0, exp0+nexp0);
  const double minexp1 = *std::min_element(exp1, exp1+nexp1);
  const double minexp2 = *std::min_element(exp2, exp2+nexp2);
  const double minexp3 = *std::min_element(exp3, exp3+nexp3);
  const double min_ab = rnd(minexp0) * rnd(minexp1);
  const double min_cd = rnd(minexp2) * rnd(minexp3);
  const double min_abp = minexp0 * minexp1;
  const double min_cdp = minexp2 * minexp3;

  // minimum distance between two lines (AB and CD) - used only for integral screening
  const double x_ab_cd = AB_[1] * CD_[2] - AB_[2] * CD_[1];
  const double y_ab_cd = AB_[2] * CD_[0] - AB_[0] * CD_[2];
  const double z_ab_cd = AB_[0] * CD_[1] - AB_[1] * CD_[0];
  const double x_n_ac = x_ab_cd * (ax - cx);
  const double y_n_ac = y_ab_cd * (ay - cy);
  const double z_n_ac = z_ab_cd * (az - cz);
  const double innerproduct = x_n_ac + y_n_ac + z_n_ac;
  const double norm_ab_cd_sq = x_ab_cd * x_ab_cd + y_ab_cd * y_ab_cd + z_ab_cd * z_ab_cd;
  const double min_pq_sq = norm_ab_cd_sq == 0.0 ? 0.0 : innerproduct * innerproduct / norm_ab_cd_sq;

  // used for prefactors Eab and Ecd (and integral screening)
  const double r01_sq = AB_[0] * AB_[0] + AB_[1] * AB_[1] + AB_[2] * AB_[2];
  const double r23_sq = CD_[0] * CD_[0] + CD_[1] * CD_[1] + CD_[2] * CD_[2];

  unsigned int tuple_length = 0u;
  double* const tuple_field = stack_->template get<double>(nexp2*nexp3*3);
  int* tuple_index = (int*)(tuple_field+nexp2*nexp3*2);

  // compute Q values and save them for later
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
        qx_save[index23] = get_PQ(cx, dx, *expi2, *expi3, cxq_inv, 2, 0, swap23_);
        qy_save[index23] = get_PQ(cy, dy, *expi2, *expi3, cxq_inv, 2, 1, swap23_);
        qz_save[index23] = get_PQ(cz, dz, *expi2, *expi3, cxq_inv, 2, 2, swap23_);

        if (IntType == Int_t::London) {
          const double A_DC_x = (basisinfo_[3]->vector_potential(0) - basisinfo_[2]->vector_potential(0));
          const double A_DC_y = (basisinfo_[3]->vector_potential(1) - basisinfo_[2]->vector_potential(1));
          const double A_DC_z = (basisinfo_[3]->vector_potential(2) - basisinfo_[2]->vector_potential(2));
          const double DC_DC_innerproduct = A_DC_x * A_DC_x + A_DC_y * A_DC_y + A_DC_z * A_DC_z;
          const double Q_DC_innerproduct = std::real(qx_save[index23]) * A_DC_x + std::real(qy_save[index23]) * A_DC_y + std::real(qz_save[index23]) * A_DC_z;
          if (swap23_) detail::make_complex(-0.25*cxq_inv*DC_DC_innerproduct, Q_DC_innerproduct, factor_cd_save[index23]);
          else detail::make_complex(-0.25*cxq_inv*DC_DC_innerproduct, -1*Q_DC_innerproduct, factor_cd_save[index23]);
        }

        // integral screening using Q
        if (integral_thresh != 0.0) {
          const double rho = cxp_min * cxq / (cxp_min + cxq);
          const double T = rho * min_pq_sq;
          const double onepqp_q = 1.0 / (std::sqrt(cxp_min + cxq) * cxp_min * cxq);
          const double abcd_sc = min_ab * cd;
          const double abcd_sc_3 = abcd_sc * abcd_sc * abcd_sc;
          const double abcd_sc_3_4 = std::sqrt(std::sqrt(abcd_sc_3));
          const double tsqrt = std::sqrt(T);
          const double ssss = 16.0 * Ecd_save[index23] * min_Eab * abcd_sc_3_4 * onepqp_q
                              * (T > 1.0e-8 ? inline_erf(tsqrt) * 0.5 / tsqrt : 1.0/std::sqrt(std::atan(1.0)*4.0) );
          if (std::abs(ssss) > integral_thresh) {
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

  // compute P values
  int index = 0;
  int index01 = 0;
  std::fill_n(T_, primsize_, -1.0);
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
      const double Eab = std::exp(-r01_sq * (abp * cxp_inv) );
      const double coeff_half = 2 * Eab * std::pow(std::atan(1.0)*4.0, 2.5);
      const DataType px = get_PQ(ax, bx, *expi0, *expi1, cxp_inv, 0, 0, swap01_);
      const DataType py = get_PQ(ay, by, *expi0, *expi1, cxp_inv, 0, 1, swap01_);
      const DataType pz = get_PQ(az, bz, *expi0, *expi1, cxp_inv, 0, 2, swap01_);

      // integral screening using P
      if (integral_thresh != 0.0) {
        const double rho_sc = cxp * cxq_min / (cxp + cxq_min);
        const double T_sc = rho_sc * min_pq_sq;
        const double onepqp_q_sc = 1.0 / (std::sqrt(cxp + cxq_min) * cxp * cxq_min);
        const double tsqrt = std::sqrt(T_sc);
        const double abcd_sc = ab * min_cd;
        const double abcd_sc_3 = abcd_sc * abcd_sc * abcd_sc;
        const double abcd_sc_3_4 = std::sqrt(std::sqrt(abcd_sc_3));
        const double ssss = 16.0 * min_Ecd * Eab * abcd_sc_3_4 * onepqp_q_sc
                          * (T_sc > 1.0e-8 ? inline_erf(tsqrt) * 0.5 / tsqrt : 1.0/std::sqrt(std::atan(1.0)*4.0) );
        if (std::abs(ssss) < std::abs(integral_thresh)) continue;
      }

      const int index_base = prim2size_ * prim3size_ * index01;

      DataType factor_ab;
      if (IntType == Int_t::London) {
        const double A_BA_x = (basisinfo_[1]->vector_potential(0) - basisinfo_[0]->vector_potential(0));
        const double A_BA_y = (basisinfo_[1]->vector_potential(1) - basisinfo_[0]->vector_potential(1));
        const double A_BA_z = (basisinfo_[1]->vector_potential(2) - basisinfo_[0]->vector_potential(2));
        const double BA_BA_innerproduct = A_BA_x * A_BA_x + A_BA_y * A_BA_y + A_BA_z * A_BA_z;
        const double P_BA_innerproduct = std::real(px) * A_BA_x + std::real(py) * A_BA_y + std::real(pz) * A_BA_z;
        if (swap01_) detail::make_complex(-0.25*cxp_inv*BA_BA_innerproduct , P_BA_innerproduct, factor_ab);
        else detail::make_complex(-0.25*cxp_inv*BA_BA_innerproduct , -1*P_BA_innerproduct, factor_ab);
      }

      // store calculated values as members of RysIntegral
      for (unsigned int i = 0; i != tuple_length; ++i) {
        const int index23 = tuple_index[i];
        const int index = index_base + index23;
        const double exp2value = tuple_field[2*i];
        const double exp3value = tuple_field[2*i+1];
        const double cxq = exp2value + exp3value;
        xp_[index] = cxp;
        xq_[index] = cxq;
        const double cxpxq = cxp * cxq;
        const double onepqp_q = 1.0 / (std::sqrt(cxp + cxq) * cxpxq);

        coeff_[index] = Ecd_save[index23] * coeff_half * onepqp_q;

        // only used for London orbitals
        // It should be generic enough to let us use complex functions other than London orbitals
        if (IntType == Int_t::London) {
          coeff_[index] *= std::exp(factor_ab + factor_cd_save[index23]);
        }

        const double rho = cxpxq / (cxp + cxq);
        const DataType xpq = qx_save[index23] - px;
        const DataType ypq = qy_save[index23] - py;
        const DataType zpq = qz_save[index23] - pz;
        const DataType T = rho * (xpq * xpq + ypq * ypq + zpq * zpq);
        const int index3 = index * 3;
        P_[index3] = px;
        P_[index3 + 1] = py;
        P_[index3 + 2] = pz;
        Q_[index3] =     qx_save[index23];
        Q_[index3 + 1] = qy_save[index23];
        Q_[index3 + 2] = qz_save[index23];
        T_[index] = T;
        screening_[screening_size_] = index;
        ++screening_size_;
/*
        if (IntType == Int_t::London && ax != bx && bx != cx && cx != dx && ax != cx && ax != dx && bx != dx) {
          std::cout << std::setprecision(14) << "A = " << ax << ", " << ay << ", " << az << std::endl;
          std::cout << "B = " << bx << ", " << by << ", " << bz << std::endl;
          std::cout << "C = " << cx << ", " << cy << ", " << cz << std::endl;
          std::cout << "D = " << dx << ", " << dy << ", " << dz << std::endl;
          std::cout << "A_A = " << basisinfo_[0]->vector_potential(0) << ", " << basisinfo_[0]->vector_potential(1) << ", " << basisinfo_[0]->vector_potential(2) << std::endl;
          std::cout << "A_B = " << basisinfo_[1]->vector_potential(0) << ", " << basisinfo_[1]->vector_potential(1) << ", " << basisinfo_[1]->vector_potential(2) << std::endl;
          std::cout << "A_C = " << basisinfo_[2]->vector_potential(0) << ", " << basisinfo_[2]->vector_potential(1) << ", " << basisinfo_[2]->vector_potential(2) << std::endl;
          std::cout << "A_D = " << basisinfo_[3]->vector_potential(0) << ", " << basisinfo_[3]->vector_potential(1) << ", " << basisinfo_[3]->vector_potential(2) << std::endl;
          std::cout << "exponents = " << *expi0 << ", " << *expi1 << ", " << exp2value << ", " << exp3value << std::endl;
          std::cout << "p = " << cxp << ", q = " << cxq << ", rho = " << rho << std::endl;
          std::cout << "P = " << P_[index3] << ", " << P_[index3+1] << ", " << P_[index3+2] << std::endl;
          std::cout << "Q = " << Q_[index3] << ", " << Q_[index3+1] << ", " << Q_[index3+2] << std::endl;
          std::cout << "T = " << T_[index] << std::endl;
          std::cout << "Eab (real) = " << Eab << ", Ecd (real) = " << Ecd_save[index23] << std::endl;
          std::cout << "factor_ab = " << factor_ab << ", factor_cd = " << factor_cd_save[index23] << std::endl;
          std::cout << "coefficient = " << coeff_[index] << std::endl << std::endl;
        }
*/
      }
    }
  }

  // release temporary Q storage
  stack_->release(nexp2*nexp3*3, tuple_field);
  if (IntType == Int_t::London) {
    stack_->release(prim2size_*prim3size_, factor_cd_save);
  }
  stack_->release(prim2size_*prim3size_, qz_save);
  stack_->release(prim2size_*prim3size_, qy_save);
  stack_->release(prim2size_*prim3size_, qx_save);
  stack_->release(prim2size_*prim3size_, Ecd_save);
}


template <typename DataType, Int_t IntType>
void ERIBatch_Base<DataType, IntType>::allocate_data(const int asize_final, const int csize_final, const int asize_final_sph, const int csize_final_sph) {
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

    stack_save_ = stack_->template get<DataType>(size_alloc_);
    stack_save2_ = nullptr;

    // if Slater/Yukawa integrals
    if (tenno_)
      stack_save2_ = stack_->get(size_alloc_);

  // derivative integrals
  } else if (deriv_rank_ == 1) {
    size_block_ = asize_final * csize_final * primsize_;
    size_alloc_ = 12 * size_block_;
    if (breit_ || tenno_)
      throw std::logic_error("Gradient integrals for the Breit and Slater operators not implemented");
    stack_save_ = stack_->template get<DataType>(size_alloc_);
    stack_save2_ = nullptr;
  }
  data_ = stack_save_;
  data2_ = stack_save2_;
}

}

#endif
#endif
