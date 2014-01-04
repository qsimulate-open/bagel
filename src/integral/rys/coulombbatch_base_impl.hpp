//
// BAGEL - Parallel electron correlation program.
// Filename: coulombbatch_base.cc
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


#ifdef COULOMBBATCH_BASE_HEADERS

#ifndef __SRC_INTEGRAL_RYS_COULOMBBATCH_BASE_HPP
#define __SRC_INTEGRAL_RYS_COULOMBBATCH_BASE_HPP


namespace bagel {

static constexpr double T_thresh__ = 1.0e-8;

template <typename DataType>
CoulombBatch_Base<DataType>::CoulombBatch_Base(const std::array<std::shared_ptr<const Shell>,2>& _info, const std::shared_ptr<const Molecule> mol, const int deriv,
                             std::shared_ptr<StackMem> stack, const int L, const double A)
 : RysIntegral<DataType>(_info, stack), mol_(mol), L_(L), A_(A) {

  deriv_rank_ = deriv;
  const double integral_thresh = PRIM_SCREEN_THRESH;

  if (_info.size() != 2) throw std::logic_error("CoulombBatch_base should be called with shell pairs");

  // natom_
  natom_ = mol_->natom() * (2 * L + 1);

  this->set_swap_info();
  this->set_ab_cd();
  this->set_prim_contsizes();

  int asize_intermediate, asize_final, dum0, dum1;
  std::tie(asize_intermediate, dum0, asize_final, dum1) = this->set_angular_info();

  allocate_data(asize_intermediate, 1, asize_final, 1);

  this->allocate_arrays(primsize_*natom_);

}

template <typename DataType>
void CoulombBatch_Base<DataType>::compute_ssss(const double integral_thresh) {
  screening_size_ = 0;

  const std::vector<double> exp0 = basisinfo_[0]->exponents();
  const std::vector<double> exp1 = basisinfo_[1]->exponents();

  int index = 0;
  std::vector<std::shared_ptr<const Atom>> atoms = mol_->atoms();

  const double onepi2 = 1.0 / (pi__ * pi__);
  const double sqrtpi = std::sqrt(pi__);
  for (auto expi0 = exp0.begin(); expi0 != exp0.end(); ++expi0) {
    for (auto expi1 = exp1.begin(); expi1 != exp1.end(); ++expi1) {
      for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++index) {
        double Z = (*aiter)->atom_charge();
        const double cxp = *expi0 + *expi1;
        xp_[index] = cxp;
        const double ab = *expi0 * *expi1;
        const double cxp_inv = 1.0 / cxp;
        P_[index * 3    ] = (basisinfo_[0]->position(0) * *expi0 + basisinfo_[1]->position(0) * *expi1) * cxp_inv;
        P_[index * 3 + 1] = (basisinfo_[0]->position(1) * *expi0 + basisinfo_[1]->position(1) * *expi1) * cxp_inv;
        P_[index * 3 + 2] = (basisinfo_[0]->position(2) * *expi0 + basisinfo_[1]->position(2) * *expi1) * cxp_inv;
        const double Eab = exp(-(AB_[0] * AB_[0] + AB_[1] * AB_[1] + AB_[2] * AB_[2]) * (ab * cxp_inv) );
        coeff_[index] = - 2 * Z * pi__ * cxp_inv * Eab;
        const double PCx = P_[index * 3    ] - (*aiter)->position(0);
        const double PCy = P_[index * 3 + 1] - (*aiter)->position(1);
        const double PCz = P_[index * 3 + 2] - (*aiter)->position(2);
        const double T = cxp * (PCx * PCx + PCy * PCy + PCz * PCz);
        const double sqrtt = sqrt(T);
        const double ss = - coeff_[index] * pow(4.0 * ab * onepi2, 0.75) * (T > 1.0e-15 ? sqrtpi * erf(sqrtt) / sqrtt * 0.5 : 1.0);
        if (ss > integral_thresh && !(*aiter)->finite_nucleus()) { // TODO reorganize loops
          T_[index] = T;
          screening_[screening_size_] = index;
          ++screening_size_;
        } else {
          T_[index] = -1.0;
          coeff_[index] = 0.0;
        }
      }
    }
  }
}

template <typename DataType>
void CoulombBatch_Base<DataType>::allocate_data(const int asize_final, const int csize_final, const int asize_final_sph, const int csize_final_sph) {
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
    // if this is an NAI gradient integral
    if (dynamic_cast<CoulombBatch_base*>(this)) {
      // in this case, we store everything
      size_alloc_ = (dynamic_cast<CoulombBatch_base*>(this)->mol()->natom()) * 3.0 * size_block_;
      assert(csize_final == 1);
    } else {
      throw std::logic_error("something is strange in CoulombBatch_base::allocate_data");
    }
    stack_save_ = stack_->get(size_alloc_);
    stack_save2_ = nullptr;
  }
  data_ = stack_save_;
  data2_ = stack_save2_;
}

template <typename DataType>
void CoulombBatch_Base<DataType>::root_weight(const int ps) {
  if (amax_ + cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int i = screening_[j];
      if (T_[i] < T_thresh__) {
        weights_[i] = 1.0;
      } else {
        const double sqrtt = std::sqrt(T_[i]);
        const double erfsqt = inline_erf(sqrtt);
        weights_[i] = erfsqt * std::sqrt(pi__) * 0.5 / sqrtt;
      }
    }
  } else {
    eriroot__.root(rank_, T_, roots_, weights_, ps);
  }
}

}

#endif
#endif
