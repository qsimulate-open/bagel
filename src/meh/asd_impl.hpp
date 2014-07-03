//
// BAGEL - Parallel electron correlation program.
// Filename: asd_impl.hpp
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifdef MEH_HEADERS

#ifndef BAGEL_MEH_ASD_IMPL_H
#define BAGEL_MEH_ASD_IMPL_H

template <>
template <class VecType>
std::shared_ptr<Matrix> asd::ASD_impl<true>::compute_diagonal_block(MultiExcitonHamiltonian<VecType>* me, DSb<VecType>& subspace) {
  const double core = me->dimer_->sref()->geom()->nuclear_repulsion() + me->jop_->core_energy();

  auto out = me->compute_intra(subspace, me->jop_, core);
  *out += *compute_inter_2e(me, subspace, subspace);

  return out;
}

// This term will couple off-diagonal blocks since it has no delta functions involved
template <>
template <class VecType>
std::shared_ptr<Matrix> asd::ASD_impl<true>::compute_inter_2e(MultiExcitonHamiltonian<VecType>* me, DSb<VecType>& AB, DSb<VecType>& ApBp) {
  auto out = std::make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());

  // alpha-alpha
  Matrix gamma_AA_alpha = *me->gammaforest_->template get<0>(AB.template tag<0>(), ApBp.template tag<0>(), {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha});
  Matrix gamma_BB_alpha = *me->gammaforest_->template get<1>(AB.template tag<1>(), ApBp.template tag<1>(), {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha});

  // beta-beta
  Matrix gamma_AA_beta = *me->gammaforest_->template get<0>(AB.template tag<0>(), ApBp.template tag<0>(), {GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta});
  Matrix gamma_BB_beta = *me->gammaforest_->template get<1>(AB.template tag<1>(), ApBp.template tag<1>(), {GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta});

  // build J and K matrices
  std::shared_ptr<const Matrix> Jmatrix = me->jop_->template coulomb_matrix<0,1,0,1>();
  std::shared_ptr<const Matrix> Kmatrix = me->jop_->template coulomb_matrix<0,1,1,0>();

  Matrix tmp((gamma_AA_alpha + gamma_AA_beta) * (*Jmatrix) ^ (gamma_BB_alpha + gamma_BB_beta));

  tmp -= gamma_AA_alpha * (*Kmatrix) ^ gamma_BB_alpha;
  tmp -= gamma_AA_beta * (*Kmatrix) ^ gamma_BB_beta;

  // sort: (A',A,B',B) --> (A,B,A',B') + block(A,B,A',B')
  SMITH::sort_indices<1,3,0,2,0,1,1,1>(tmp.data(), out->data(), ApBp.template nstates<0>(), AB.template nstates<0>(), ApBp.template nstates<1>(), AB.template nstates<1>());
  return out;
}

#endif
#endif
