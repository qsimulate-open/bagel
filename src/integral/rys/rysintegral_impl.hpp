//
// BAGEL - Parallel electron correlation program.
// Filename: rysintegral_impl.hpp
// Copyright (C) 2013 Toru Shiozaki
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

#ifdef RYSINTEGRAL_HEADERS

#ifndef __SRC_INTEGRAL_RYS_RYSINTEGRAL_IMPL_HPP
#define __SRC_INTEGRAL_RYS_RYSINTEGRAL_IMPL_HPP

namespace bagel {

template<typename DataType>
void RysIntegral<DataType>::set_ab_cd() {
  AB_[0] = basisinfo_[0]->position(0) - basisinfo_[1]->position(0);
  AB_[1] = basisinfo_[0]->position(1) - basisinfo_[1]->position(1);
  AB_[2] = basisinfo_[0]->position(2) - basisinfo_[1]->position(2);
  CD_[0] = basisinfo_[2]->position(0) - basisinfo_[3]->position(0);
  CD_[1] = basisinfo_[2]->position(1) - basisinfo_[3]->position(1);
  CD_[2] = basisinfo_[2]->position(2) - basisinfo_[3]->position(2);
}


template<typename DataType>
void RysIntegral<DataType>::set_prim_contsizes() {
  prim0size_ = basisinfo_[0]->num_primitive();
  prim1size_ = basisinfo_[1]->num_primitive();
  prim2size_ = basisinfo_[2]->num_primitive();
  prim3size_ = basisinfo_[3]->num_primitive();
  primsize_ = prim0size_ * prim1size_ * prim2size_ * prim3size_;
  cont0size_ = basisinfo_[0]->num_contracted();
  cont1size_ = basisinfo_[1]->num_contracted();
  cont2size_ = basisinfo_[2]->num_contracted();
  cont3size_ = basisinfo_[3]->num_contracted();
  contsize_ = cont0size_ * cont1size_ * cont2size_ * cont3size_;
}


template<typename DataType>
std::tuple<int, int, int, int> RysIntegral<DataType>::set_angular_info() {
  const int ang0 = basisinfo_[0]->angular_number();
  const int ang1 = basisinfo_[1]->angular_number();
  const int ang2 = basisinfo_[2]->angular_number();
  const int ang3 = basisinfo_[3]->angular_number();
  rank_ = ceil(0.5 * (ang0 + ang1 + ang2 + ang3 + 1 + deriv_rank_ + tenno_ + breit_));
  assert(2 * rank_ >= ang0 + ang1 + ang2 + ang3 + 1 + deriv_rank_ + tenno_ + breit_);

  amax_ = ang0 + ang1 + deriv_rank_;
  cmax_ = ang2 + ang3 + deriv_rank_;
  amin_ = std::max(ang0 - deriv_rank_, 0);
  cmin_ = std::max(ang2 - deriv_rank_, 0);
  amax1_ = amax_ + 1;
  cmax1_ = cmax_ + 1;

  asize_ = 0;
  csize_ = 0;
  for (int i = amin_; i <= amax_; ++i) asize_ += (i + 1) * (i + 2) / 2;
  for (int i = cmin_; i <= cmax_; ++i) csize_ += (i + 1) * (i + 2) / 2;

  const int asize_final = (ang0 + 1) * (ang0 + 2) * (ang1 + 1) * (ang1 + 2) / 4;
  const int csize_final = (ang2 + 1) * (ang2 + 2) * (ang3 + 1) * (ang3 + 2) / 4;

  const int asize_final_sph = spherical1_ ? (2 * ang0 + 1) * (2 * ang1 + 1) : asize_final;
  const int csize_final_sph = spherical2_ ? (2 * ang2 + 1) * (2 * ang3 + 1) : csize_final;

  int cnt = 0;
  for (int i = cmin_; i <= cmax_; ++i) {
    for (int iz = 0; iz <= i; ++iz) {
      for (int iy = 0; iy <= i - iz; ++iy) {
        const int ix = i - iy - iz;
        if (ix >= 0)
          cmapping_[ix + cmax1_ * (iy + cmax1_ * iz)] = cnt++;
      }
    }
  }
  cnt = 0;
  for (int j = amin_; j <= amax_; ++j) {
    for (int jz = 0; jz <= j; ++jz) {
      for (int jy = 0; jy <= j - jz; ++jy) {
        const int jx = j - jy - jz;
        if (jx >= 0)
          amapping_[jx + amax1_ * (jy + amax1_ * jz)] = cnt++;
      }
    }
  }
  return std::make_tuple(asize_final, csize_final, asize_final_sph, csize_final_sph);
}


template<typename DataType>
void RysIntegral<DataType>::set_swap_info(const bool swap_bra_ket) {
  // swap 01 indices when needed: Larger angular momentum function comes first
  if (basisinfo_[0]->angular_number() < basisinfo_[1]->angular_number() || basisinfo_[0]->angular_number() == 0) {
    std::swap(basisinfo_[0], basisinfo_[1]);
    swap01_ = true;
  } else {
    swap01_ = false;
  }
  // swap 23 indices when needed
  if (basisinfo_[2]->angular_number() < basisinfo_[3]->angular_number() || basisinfo_[2]->angular_number() == 0) {
    std::swap(basisinfo_[2], basisinfo_[3]);
    swap23_ = true;
  } else {
    swap23_ = false;
  }

  swap0123_ = false;
  if (swap_bra_ket) {
    if (!basisinfo_[0]->angular_number() && !basisinfo_[2]->angular_number()) {
      swap0123_ = true;
      std::tie(basisinfo_[0], basisinfo_[1], basisinfo_[2], basisinfo_[3], swap01_, swap23_)
        = std::make_tuple(basisinfo_[2], basisinfo_[3], basisinfo_[0], basisinfo_[1], swap23_, swap01_);
      std::swap(spherical1_, spherical2_);
    }
  }
}


template<typename DataType>
void RysIntegral<DataType>::allocate_arrays(const size_t ps) {
  size_allocated_ = tenno_ > 0 ? ((rank_ * 2 + 13) * ps) : ((rank_ * 2 + 11) * ps);

  buff_ = stack_->get<DataType>(size_allocated_);  // stack_->get(size_alloc_) stack_->get((rank_ * 2 + 10) * ps)
  DataType* pointer = buff_;
  screening_ = (int*)pointer;
  pointer += ps;
  P_ = pointer;     pointer += ps * 3;
  Q_ = pointer;     pointer += ps * 3;
  xp_ = reinterpret_cast<double*>(pointer);    pointer += ps;
  xq_ = reinterpret_cast<double*>(pointer);    pointer += ps;
  coeff_ = pointer; pointer += ps;
  T_ = pointer;     pointer += ps;
  roots_ = pointer; pointer += rank_ * ps;
  weights_ = pointer; pointer += rank_ * ps;
  if (tenno_) {
    coeffy_ = pointer;pointer += ps;
    U_ = pointer;     pointer += ps;
  }
}


template<typename DataType>
void RysIntegral<DataType>::perform_contraction_new_outer(const int nsize, const DataType* prim, const int pdim0, const int pdim1, DataType* cont,
                     const std::vector<std::vector<double>>& coeff0, const std::vector<int>& upper0, const std::vector<int>& lower0, const int cdim0,
                     const std::vector<std::vector<double>>& coeff1, const std::vector<int>& upper1, const std::vector<int>& lower1, const int cdim1) {
  const int worksize = nsize * pdim1;
  DataType* const work = stack_->get<DataType>(worksize);
  DataType* current_cont = cont;

  for (int i = 0; i != cdim0; ++i) {
    const int begin0 = lower0[i];
    const int end0   = upper0[i];
    std::fill_n(work, worksize, DataType(0.0));
    for (int j = begin0; j != end0; ++j)
      blas::ax_plus_y_n(coeff0[i][j], prim+j*worksize, worksize, work);

    for (int k = 0; k != cdim1; ++k, current_cont += nsize) {
      const int begin1 = lower1[k];
      const int end1   = upper1[k];
      std::fill_n(current_cont, nsize, DataType(0.0));
      for (int j = begin1; j != end1; ++j)
        blas::ax_plus_y_n(coeff1[k][j], work+j*nsize, nsize, current_cont);
    }
  }

  stack_->release(worksize, work);
}


template<typename DataType>
void RysIntegral<DataType>::perform_contraction_new_inner(const int nsize, const int ac, const DataType* prim, const int pdim0, const int pdim1, DataType* cont,
                 const std::vector<std::vector<double>>& coeff0, const std::vector<int>& upper0, const std::vector<int>& lower0, const int cdim0,
                 const std::vector<std::vector<double>>& coeff1, const std::vector<int>& upper1, const std::vector<int>& lower1, const int cdim1) {
  const int worksize = pdim1 * ac;
  DataType* const work = stack_->get<DataType>(worksize);
  DataType* current_cont = cont;

  for (int n = 0; n != nsize; ++n) { // loop of cdim * cdim
    const DataType* current_prim = &prim[ac * pdim1 * pdim0 * n];

    for (int i = 0; i != cdim0; ++i) {

      const int begin0 = lower0[i];
      const int end0   = upper0[i];
      std::fill_n(work, worksize,  DataType(0.0));
      for (int j = begin0; j != end0; ++j)
        blas::ax_plus_y_n(coeff0[i][j], current_prim+j*worksize, worksize, work);

      for (int k = 0; k != cdim1; ++k, current_cont += ac) {
        const int begin1 = lower1[k];
        const int end1   = upper1[k];
        std::fill_n(current_cont, ac, DataType(0.0));
        for (int j = begin1; j != end1; ++j)
          blas::ax_plus_y_n(coeff1[k][j], work+j*ac, ac, current_cont);
      }
    }
  }
  stack_->release(worksize, work);
}


template<typename DataType>
void RysIntegral<DataType>::perform_contraction(const int asize, const DataType* prim, const int pdim0, const int pdim1, DataType* cont,
                           const std::vector<std::vector<double>>& coeff0, const std::vector<std::pair<int, int>>& ranges0, const int cdim0,
                           const std::vector<std::vector<double>>& coeff1, const std::vector<std::pair<int, int>>& ranges1, const int cdim1) {
  // transformation of index1
  const int worksize = pdim1 * asize;
  DataType* const work = stack_->get<DataType>(worksize);

  for (int i = 0; i != cdim0; ++i) {
    const int begin0 = ranges0[i].first;
    const int end0   = ranges0[i].second;
    std::fill_n(work, worksize, DataType(0.0));
    for (int j = begin0; j != end0; ++j)
      blas::ax_plus_y_n(coeff0[i][j], prim+j*worksize, worksize, work);

    for (int k = 0; k != cdim1; ++k, cont += asize) {
      const int begin1 = ranges1[k].first;
      const int end1   = ranges1[k].second;
      std::fill_n(cont, asize, DataType(0.0));
      for (int j = begin1; j != end1; ++j)
        blas::ax_plus_y_n(coeff1[k][j], work+j*asize, asize, cont);
    }
  }
  stack_->release(worksize, work);
}

}

#endif
#endif
