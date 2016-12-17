//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: osintegral.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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

#include <src/integral/os/osintegral.h>

using namespace std;
using namespace bagel;

static const double pisqrt__ = sqrt(pi__);

template <typename DataType, Int_t IntType>
OSIntegral<DataType, IntType>::OSIntegral(const array<shared_ptr<const Shell>,2>& basis, shared_ptr<StackMem> stack)
 : basisinfo_(basis), spherical_(basis.front()->spherical()) {

  if (stack == nullptr) {
    stack_ = resources__->get();
    allocated_here_ = true;
  } else {
    stack_ = stack;
    allocated_here_ = false;
  }
}


template <typename DataType, Int_t IntType>
void OSIntegral<DataType, IntType>::common_init() {
  static_assert(IntType != Int_t::London || is_same<DataType, complex<double>>::value, "London-orbital integrals should be complex");

  assert(basisinfo_.size() == 2);

  ang0_ = basisinfo_[0]->angular_number();
  ang1_ = basisinfo_[1]->angular_number();

  if (ang0_ < ang1_) {
    swap(basisinfo_[0], basisinfo_[1]);
    swap(ang0_, ang1_);
    swap01_ = true;
  } else {
    swap01_ = false;
  }

  prim0_ = basisinfo_[0]->num_primitive();
  cont0_ = basisinfo_[0]->num_contracted();
  const vector<double> exponents0 = basisinfo_[0]->exponents();
  prim1_ = basisinfo_[1]->num_primitive();
  cont1_ = basisinfo_[1]->num_contracted();
  const vector<double> exponents1 = basisinfo_[1]->exponents();

  AB_[0] = basisinfo_[0]->position(0) - basisinfo_[1]->position(0);
  AB_[1] = basisinfo_[0]->position(1) - basisinfo_[1]->position(1);
  AB_[2] = basisinfo_[0]->position(2) - basisinfo_[1]->position(2);

  vector<double>::const_iterator expi0, expi1;
  P_.reserve(3 * prim0_ * prim1_);
  xa_.reserve(prim0_ * prim1_);
  xb_.reserve(prim0_ * prim1_);
  xp_.reserve(prim0_ * prim1_);
  coeffsx_.reserve(prim0_ * prim1_);
  coeffsy_.reserve(prim0_ * prim1_);
  coeffsz_.reserve(prim0_ * prim1_);
  coefftx_.reserve(prim0_ * prim1_);
  coeffty_.reserve(prim0_ * prim1_);
  coefftz_.reserve(prim0_ * prim1_);
  for (expi0 = exponents0.begin(); expi0 != exponents0.end(); ++expi0) {
    for (expi1 = exponents1.begin(); expi1 != exponents1.end(); ++expi1) {
      xa_.push_back(*expi0);
      xb_.push_back(*expi1);
      const double cxp = *expi0 + *expi1;
      const double cxp_inv = 1.0 / cxp;

      const DataType px = get_P(basisinfo_[0]->position(0), basisinfo_[1]->position(0), *expi0, *expi1, cxp_inv, 0, swap01_);
      const DataType py = get_P(basisinfo_[0]->position(1), basisinfo_[1]->position(1), *expi0, *expi1, cxp_inv, 1, swap01_);
      const DataType pz = get_P(basisinfo_[0]->position(2), basisinfo_[1]->position(2), *expi0, *expi1, cxp_inv, 2, swap01_);

      xp_.push_back(cxp);
      P_.push_back(px);
      P_.push_back(py);
      P_.push_back(pz);
      const double tmp = pisqrt__ * sqrt(cxp_inv);
      DataType coeffsx = tmp * exp(- *expi0 * *expi1 * cxp_inv * (AB_[0] * AB_[0]));
      DataType coeffsy = tmp * exp(- *expi0 * *expi1 * cxp_inv * (AB_[1] * AB_[1]));
      DataType coeffsz = tmp * exp(- *expi0 * *expi1 * cxp_inv * (AB_[2] * AB_[2]));
      if (IntType == Int_t::London) {
        DataType factor_x;
        DataType factor_y;
        DataType factor_z;
        const double A_BA_x = (basisinfo_[1]->vector_potential(0) - basisinfo_[0]->vector_potential(0));
        const double A_BA_y = (basisinfo_[1]->vector_potential(1) - basisinfo_[0]->vector_potential(1));
        const double A_BA_z = (basisinfo_[1]->vector_potential(2) - basisinfo_[0]->vector_potential(2));
        const double BA_xsq = A_BA_x * A_BA_x;
        const double BA_ysq = A_BA_y * A_BA_y;
        const double BA_zsq = A_BA_z * A_BA_z;
        const double P_BA_x = real(px) * A_BA_x;
        const double P_BA_y = real(py) * A_BA_y;
        const double P_BA_z = real(pz) * A_BA_z;
        if (swap01_) {
          detail::make_complex(-0.25*cxp_inv*BA_xsq , P_BA_x, factor_x);
          detail::make_complex(-0.25*cxp_inv*BA_ysq , P_BA_y, factor_y);
          detail::make_complex(-0.25*cxp_inv*BA_zsq , P_BA_z, factor_z);
        } else {
          detail::make_complex(-0.25*cxp_inv*BA_xsq , -1*P_BA_x, factor_x);
          detail::make_complex(-0.25*cxp_inv*BA_ysq , -1*P_BA_y, factor_y);
          detail::make_complex(-0.25*cxp_inv*BA_zsq , -1*P_BA_z, factor_z);
        }
        coeffsx *= exp(factor_x);
        coeffsy *= exp(factor_y);
        coeffsz *= exp(factor_z);
      }
      //ss overlap integral
      coeffsx_.push_back(coeffsx);
      coeffsy_.push_back(coeffsy);
      coeffsz_.push_back(coeffsz);

      const DataType xpa = px - basisinfo_[0]->position(0);
      const DataType ypa = py - basisinfo_[0]->position(1);
      const DataType zpa = pz - basisinfo_[0]->position(2);

      // ss kinetic energy integral (not used for London orbitals)
      coefftx_.push_back(coeffsx_.back() * (*expi0 - 2 * *expi0 * *expi0 * (xpa * xpa + 0.5 * cxp_inv)));
      coeffty_.push_back(coeffsy_.back() * (*expi0 - 2 * *expi0 * *expi0 * (ypa * ypa + 0.5 * cxp_inv)));
      coefftz_.push_back(coeffsz_.back() * (*expi0 - 2 * *expi0 * *expi0 * (zpa * zpa + 0.5 * cxp_inv)));
    }
  }

  assert(P_.size() == 3 * prim0_ * prim1_);
  assert(xp_.size() == prim0_ * prim1_);
  assert(coeffsx_.size() == prim0_ * prim1_);
  assert(coefftx_.size() == prim0_ * prim1_);

  amax_ = ang0_ + ang1_ + nrank();
  amax1_ = amax_ + 1;
  amin_ = ang0_;

  asize_ = 0;
  for (int i = amin_; i != amax1_; ++i) asize_ += (i+1)*(i+2) / 2;
  asize_intermediate_ = (ang0_+1) * (ang0_+2) * (ang1_+1) * (ang1_+2) / 4;

  // note: for relativistic cases, mixed spherical and cartesian basis is considered
  asize_final_ = (basisinfo_[0]->spherical() ? (2*ang0_+1) : (ang0_+1)*(ang0_+2)/2)
               * (basisinfo_[1]->spherical() ? (2*ang1_+1) : (ang1_+1)*(ang1_+2)/2);

  // note: this is larger than what is needed. just to make the code simpler..
  size_block_ = prim0_ * prim1_ * max(asize_intermediate_, asize_);
  size_alloc_ = size_block_ * nblocks();

  stack_save_ = stack_->template get<DataType>(size_alloc_);
  data_ = stack_save_;

  amapping_.resize(amax1_ * amax1_ * amax1_);
  int cnt = 0;
  for (int i = amin_; i <= amax_; ++i) {
    for (int iz = 0; iz <= i; ++iz) {
      for (int iy = 0; iy <= i - iz; ++iy) {
        const int ix = i - iy - iz;
        amapping_[ix + amax1_ * (iy + amax1_ * iz)] = cnt;
        ++cnt;
      }
    }
  }
}


template <typename DataType, Int_t IntType>
OSIntegral<DataType, IntType>::~OSIntegral() {
  stack_->release(size_alloc_, stack_save_);
  if (allocated_here_) resources__->release(stack_);
}


template <typename DataType, Int_t IntType>
shared_ptr<GradFile> OSIntegral<DataType, IntType>::compute_gradient(shared_ptr<const Matrix> d, const int iatom0, const int iatom1, const int natom) const {
  if (nblocks() != 6) throw logic_error("OSIntegral::contract_density called unexpectedly");
  auto out = make_shared<GradFile>(natom);
  const int jatom0 = swap01() ? iatom1 : iatom0;
  const int jatom1 = swap01() ? iatom0 : iatom1;

  // TODO This is used to avoid compiler errors in ComplexOverlapBatch.  Probably better to define compute_gradient in derived classes.
  const double* data = reinterpret_cast<double*> (data_);
  if (IntType == Int_t::London) throw runtime_error("Gradient computation has not been set up for London orbitals");

  for (int k = 0; k != 3; ++k) {
    out->element(k, jatom1) += ddot_(d->size(), d->data(), 1, data+size_block_*k, 1);
    out->element(k, jatom0) += ddot_(d->size(), d->data(), 1, data+size_block_*(k+3), 1);
  }
  return out;
}


template <typename DataType, Int_t IntType>
void OSIntegral<DataType, IntType>::perform_contraction(const int asize, const DataType* prim, const int pdim0, const int pdim1, DataType* cont,
                                                        const vector<vector<double>>& coeff0, const vector<pair<int, int>>& ranges0, const int cdim0,
                                                        const vector<vector<double>>& coeff1, const vector<pair<int, int>>& ranges1, const int cdim1) {
  // transformation of index1
  const int worksize = pdim1 * asize;
  DataType* const work = stack_->template get<DataType>(worksize);
  fill_n(cont, asize*cdim0*cdim1, 0.0);

  for (int i = 0; i != cdim0; ++i) {
    const int begin0 = ranges0[i].first;
    const int end0   = ranges0[i].second;
    fill_n(work, worksize, 0.0);
    for (int j = begin0; j != end0; ++j)
      for (int n=0; n!=worksize; n++) work[n] += coeff0[i][j]*prim[j*worksize+n];
    for (int k = 0; k != cdim1; ++k, cont += asize) {
      const int begin1 = ranges1[k].first;
      const int end1   = ranges1[k].second;
      for (int j = begin1; j != end1; ++j) {
        blas::ax_plus_y_n(coeff1[k][j], work+j*asize, asize, cont);
      }
    }
  }
  stack_->release(worksize, work);
}


template class bagel::OSIntegral<double,bagel::Int_t::Standard>;
template class bagel::OSIntegral<complex<double>,bagel::Int_t::London>;
template class bagel::OSIntegral<complex<double>,bagel::Int_t::Standard>;
