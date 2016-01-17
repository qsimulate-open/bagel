//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: coulombbatch_base.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/integral/rys/coulombbatch_base.h>

using namespace std;
using namespace bagel;


template <typename DataType, Int_t IntType>
CoulombBatch_Base<DataType, IntType>::CoulombBatch_Base(const array<shared_ptr<const Shell>,2>& _info, shared_ptr<const Molecule> mol, const int deriv, const int breit,
                                                        shared_ptr<StackMem> stack)
 : RysIntegral<DataType, IntType>(_info, stack), mol_(mol) {

  breit_ = breit;
  deriv_rank_ = deriv;

  if (_info.size() != 2) throw logic_error("CoulombBatch_base should be called with shell pairs");

  // natom_
  natom_ = mol_->natom();

  this->set_swap_info();
  this->set_ab_cd();
  this->set_prim_contsizes();

  int asize_intermediate, asize_final, dum0, dum1;
  tie(asize_intermediate, dum0, asize_final, dum1) = this->set_angular_info();

  allocate_data(asize_intermediate, 1, asize_final, 1);

}

template <typename DataType, Int_t IntType>
void CoulombBatch_Base<DataType, IntType>::compute_ssss(const double integral_thresh) {
  static_assert(IntType != Int_t::London || is_same<DataType, complex<double>>::value, "London-orbital integrals should be complex");
  static_assert(IntType != Int_t::Standard || is_same<DataType, double>::value, "Standard Guassian-orbital integrals should be real");
  screening_size_ = 0;

  const vector<double> exp0 = basisinfo_[0]->exponents();
  const vector<double> exp1 = basisinfo_[1]->exponents();

  int index = 0;
  vector<shared_ptr<const Atom>> atoms = mol_->atoms();

  const double onepi2 = 1.0 / (pi__ * pi__);
  const double sqrtpi = sqrt(pi__);
  for (auto expi0 = exp0.begin(); expi0 != exp0.end(); ++expi0) {
    for (auto expi1 = exp1.begin(); expi1 != exp1.end(); ++expi1) {
      const double cxp = *expi0 + *expi1;
      const double ab = *expi0 * *expi1;
      const double cxp_inv = 1.0 / cxp;
      const DataType px = get_PQ(basisinfo_[0]->position(0), basisinfo_[1]->position(0), *expi0, *expi1, cxp_inv, 0, 0, swap01_);
      const DataType py = get_PQ(basisinfo_[0]->position(1), basisinfo_[1]->position(1), *expi0, *expi1, cxp_inv, 0, 1, swap01_);
      const DataType pz = get_PQ(basisinfo_[0]->position(2), basisinfo_[1]->position(2), *expi0, *expi1, cxp_inv, 0, 2, swap01_);
      // For London orbitals, calculate the correction needed for the pre-integral coefficient
      DataType factor_ab;
      if (IntType == Int_t::London) {
        const double A_BA_x = (basisinfo_[1]->vector_potential(0) - basisinfo_[0]->vector_potential(0));
        const double A_BA_y = (basisinfo_[1]->vector_potential(1) - basisinfo_[0]->vector_potential(1));
        const double A_BA_z = (basisinfo_[1]->vector_potential(2) - basisinfo_[0]->vector_potential(2));
        const double BA_BA_innerproduct = A_BA_x * A_BA_x + A_BA_y * A_BA_y + A_BA_z * A_BA_z;
        const double P_BA_innerproduct = real(px) * A_BA_x + real(py) * A_BA_y + real(pz) * A_BA_z;
        if (swap01_) detail::make_complex(-0.25*cxp_inv*BA_BA_innerproduct , P_BA_innerproduct, factor_ab);
        else detail::make_complex(-0.25*cxp_inv*BA_BA_innerproduct , -1*P_BA_innerproduct, factor_ab);
      }
      // Now loop over all nuclei
      for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter, ++index) {
        double Z = (*aiter)->atom_charge();
        if ((*aiter)->use_ecp_basis()) Z -= (*aiter)->ecp_parameters()->ecp_ncore();
        xp_[index] = cxp;
        P_[index * 3    ] = px;
        P_[index * 3 + 1] = py;
        P_[index * 3 + 2] = pz;
        const double Eab = exp(-(AB_[0] * AB_[0] + AB_[1] * AB_[1] + AB_[2] * AB_[2]) * (ab * cxp_inv) );
        const double coeff_real = - 2 * Z * pi__ * cxp_inv * Eab;
        coeff_[index] = coeff_real;
        if (IntType == Int_t::London) coeff_[index] *= exp(factor_ab);
        const DataType PCx = P_[index * 3    ] - (*aiter)->position(0);
        const DataType PCy = P_[index * 3 + 1] - (*aiter)->position(1);
        const DataType PCz = P_[index * 3 + 2] - (*aiter)->position(2);
        const DataType T = cxp * (PCx * PCx + PCy * PCy + PCz * PCz);
        const double sqrtt = sqrt(abs(T));
        const double ss = - coeff_real * pow(4.0 * ab * onepi2, 0.75) * (abs(T) > 1.0e-15 ? sqrtpi * erf(sqrtt) / sqrtt * 0.5 : 1.0);
        if (abs(ss) > integral_thresh && !(*aiter)->finite_nucleus()) { // TODO reorganize loops
          T_[index] = T;
          screening_[screening_size_] = index;
          ++screening_size_;
        } else {
          T_[index] = nan("");
          coeff_[index] = 0.0;
        }
      }
    }
  }
}

template <typename DataType, Int_t IntType>
void CoulombBatch_Base<DataType, IntType>::allocate_data(const int asize_final, const int csize_final, const int asize_final_sph, const int csize_final_sph) {
  size_final_ = asize_final_sph * csize_final_sph * contsize_;
  if (deriv_rank_ == 0) {
    const unsigned int size_start = asize_ * csize_ * primsize_;
    const unsigned int size_intermediate = asize_final * csize_ * contsize_;
    const unsigned int size_intermediate2 = asize_final_sph * csize_final * contsize_;
    size_block_ = max(size_start, max(size_intermediate, size_intermediate2));
    size_alloc_ = size_block_;

    // if this is a one-electron Breit-type integral
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
    // if this is an NAI gradient integral
    if (is_same<DataType,double>::value) {
      assert(IntType == Int_t::Standard);
      // in this case, we store everything
      size_alloc_ = mol()->natom() * 3 * size_block_;
      assert(csize_final == 1);
    } else {
      throw logic_error("something is strange in CoulombBatch_base::allocate_data");
    }
    stack_save_ = stack_->template get<DataType>(size_alloc_);
    stack_save2_ = nullptr;
  }
  data_ = stack_save_;
  data2_ = stack_save2_;
}


template class bagel::CoulombBatch_Base<double>;
template class bagel::CoulombBatch_Base<complex<double>,Int_t::London>;
