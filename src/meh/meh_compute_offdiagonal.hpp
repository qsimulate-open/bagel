//
// BAGEL - Parallel electron correlation program.
// Filename: meh_compute_offdiagonal.hpp
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef BAGEL_MEH_COMPUTE_OFFDIAGONAL_H
#define BAGEL_MEH_COMPUTE_OFFDIAGONAL_H

template <class VecType>
std::shared_ptr<Matrix> MultiExcitonHamiltonian<VecType>::compute_offdiagonal_1e(const DSubSpace& AB, const DSubSpace& ApBp, std::shared_ptr<const Matrix> hAB) const {
  Coupling term_type = coupling_type(AB,ApBp);

  GammaSQ operatorA;
  GammaSQ operatorB;
  int neleA = AB.template ci<0>()->det()->nelea() + AB.template ci<0>()->det()->neleb();

  switch(term_type) {
    case Coupling::aET :
      operatorA = GammaSQ::CreateAlpha;
      operatorB = GammaSQ::AnnihilateAlpha;
      break;
    case Coupling::inv_aET :
      operatorA = GammaSQ::AnnihilateAlpha;
      operatorB = GammaSQ::CreateAlpha;
      --neleA;
      break;
    case Coupling::bET :
      operatorA = GammaSQ::CreateBeta;
      operatorB = GammaSQ::AnnihilateBeta;
      break;
    case Coupling::inv_bET :
      operatorA = GammaSQ::AnnihilateBeta;
      operatorB = GammaSQ::CreateBeta;
      --neleA;
      break;
    default :
      return std::make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());
  }

  Matrix gamma_A = *gammaforest_->template get<0>(AB.offset(), ApBp.offset(), operatorA);
  Matrix gamma_B = *gammaforest_->template get<1>(AB.offset(), ApBp.offset(), operatorB);
  Matrix tmp = gamma_A * (*hAB) ^ gamma_B;

  auto out = std::make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());

  if ((neleA % 2) == 1) {
    // sort: (A',A,B',B) --> -1.0 * (A,B,A',B')
    SMITH::sort_indices<1,3,0,2,0,1,-1,1>(tmp.data(), out->data(), ApBp.template nstates<0>(), AB.template nstates<0>(), ApBp.template nstates<1>(), AB.template nstates<1>());
  }
  else {
    // sort: (A',A,B',B) --> (A,B,A',B')
    SMITH::sort_indices<1,3,0,2,0,1,1,1>(tmp.data(), out->data(), ApBp.template nstates<0>(), AB.template nstates<0>(), ApBp.template nstates<1>(), AB.template nstates<1>());
  }

  return out;
}


namespace {
  template<typename T>
  void transpose_call(std::shared_ptr<T>& o) { assert(false); }
  template<>
  void transpose_call(std::shared_ptr<Matrix>& o) { o = o->transpose(); }
  template<>
  void transpose_call(std::shared_ptr<RDM<2>>& o) { /* doing nothing */ }
}

template <class VecType>
template <bool _N, typename return_type>
std::shared_ptr<return_type> MultiExcitonHamiltonian<VecType>::couple_blocks(DSubSpace& AB, DSubSpace& ApBp) {

  Coupling term_type = coupling_type(AB, ApBp);

  DSubSpace* space1 = &AB;
  DSubSpace* space2 = &ApBp;

  bool flip = (static_cast<int>(term_type) < 0);

  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    std::swap(space1,space2);
  }

  std::shared_ptr<return_type> out;

  switch(term_type) {
    case Coupling::none :
      out = nullptr; break;
    case Coupling::diagonal :
      out = std::make_shared<return_type>(AB.dimerstates(), ApBp.dimerstates());
      compute_inter_2e<_N>(*out, *space1, *space2); break;
    case Coupling::aET :
      out = compute_aET<_N>(*space1, *space2); break;
    case Coupling::bET :
      out = compute_bET<_N>(*space1, *space2); break;
    case Coupling::abFlip :
      out = compute_abFlip<_N>(*space1, *space2); break;
    case Coupling::abET :
      out = compute_abET<_N>(*space1, *space2); break;
    case Coupling::aaET :
      out = compute_aaET<_N>(*space1, *space2); break;
    case Coupling::bbET :
      out = compute_bbET<_N>(*space1, *space2); break;
    default :
      throw std::logic_error("Asking for a coupling type that has not been written.");
  }

  /* if we are computing the Hamiltonian and flip = true, then we tranpose the output (see above) */
  if (flip) transpose_call(out);

  return out;
}


template <>
template <class VecType>
std::shared_ptr<Matrix> asd::ASD_impl<true>::compute_aET(MultiExcitonHamiltonian<VecType>* me, DSb<VecType>& AB, DSb<VecType>& ApBp) {
  auto out = std::make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());
  Matrix tmp(AB.template nstates<0>() * ApBp.template nstates<0>(), AB.template nstates<1>() * ApBp.template nstates<1>());

  // One-body aET
  {
    Matrix gamma_A = *me->gammaforest_->template get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateAlpha);
    Matrix gamma_B = *me->gammaforest_->template get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha);

    std::shared_ptr<const Matrix> Fmatrix = me->jop_->cross_mo1e();

    tmp += gamma_A * (*Fmatrix) ^ gamma_B;
  }

  //Two-body aET, type 1
  {
    Matrix gamma_A = *me->gammaforest_->template get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateAlpha);
    Matrix gamma_B1 = *me->gammaforest_->template get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha);
    Matrix gamma_B2 = *me->gammaforest_->template get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta);

    std::shared_ptr<const Matrix> Jmatrix = me->jop_->template coulomb_matrix<0,1,1,1>();

    tmp -= gamma_A * (*Jmatrix) ^ (gamma_B1 + gamma_B2);
  }

  //Two-body aET, type 2
  {
    Matrix gamma_A1 = *me->gammaforest_->template get<0>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha);
    Matrix gamma_A2 = *me->gammaforest_->template get<0>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha);
    Matrix gamma_B = *me->gammaforest_->template get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha);

    std::shared_ptr<const Matrix> Jmatrix = me->jop_->template coulomb_matrix<0,0,1,0>();

    tmp += (gamma_A1 + gamma_A2) * (*Jmatrix) ^ gamma_B;
  }

  const int neleA = AB.template ci<0>()->det()->nelea() + AB.template ci<0>()->det()->neleb();
  if ((neleA % 2) == 1) {
    // sort: (A',A,B',B) --> -1.0 * (A,B,A',B')
    SMITH::sort_indices<1,3,0,2,0,1,-1,1>(tmp.data(), out->data(), ApBp.template nstates<0>(), AB.template nstates<0>(), ApBp.template nstates<1>(), AB.template nstates<1>());
  }
  else {
    // sort: (A',A,B',B) --> (A,B,A',B')
    SMITH::sort_indices<1,3,0,2,0,1,1,1>(tmp.data(), out->data(), ApBp.template nstates<0>(), AB.template nstates<0>(), ApBp.template nstates<1>(), AB.template nstates<1>());
  }

  return out;
}


template <>
template <class VecType>
std::shared_ptr<Matrix> asd::ASD_impl<true>::compute_bET(MultiExcitonHamiltonian<VecType>* me, DSb<VecType>& AB, DSb<VecType>& ApBp) {
  auto out = std::make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());
  Matrix tmp(AB.template nstates<0>() * ApBp.template nstates<0>(), AB.template nstates<1>() * ApBp.template nstates<1>());

  // One-body bET
  {
    Matrix gamma_A = *me->gammaforest_->template get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateBeta);
    Matrix gamma_B = *me->gammaforest_->template get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta);

    std::shared_ptr<const Matrix> Fmatrix = me->jop_->cross_mo1e();

    tmp += gamma_A * (*Fmatrix) ^ gamma_B;
  }


  //Two-body bET, type 1
  {
    Matrix gamma_A = *me->gammaforest_->template get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateBeta);
    Matrix gamma_B1 = *me->gammaforest_->template get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::CreateAlpha);
    Matrix gamma_B2 = *me->gammaforest_->template get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta);

    std::shared_ptr<const Matrix> Jmatrix = me->jop_->template coulomb_matrix<0,1,1,1>();

    tmp -= gamma_A * (*Jmatrix) ^ (gamma_B1 + gamma_B2);
  }

  //Two-body aET, type 2
  {
    Matrix gamma_A1 = *me->gammaforest_->template get<0>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta);
    Matrix gamma_A2 = *me->gammaforest_->template get<0>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta);
    Matrix gamma_B = *me->gammaforest_->template get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta);

    std::shared_ptr<const Matrix> Jmatrix = me->jop_->template coulomb_matrix<0,0,1,0>();

    tmp += (gamma_A1 + gamma_A2) * (*Jmatrix) ^ gamma_B;
  }

  const int neleA = AB.template ci<0>()->det()->nelea() + AB.template ci<0>()->det()->neleb();
  if ((neleA % 2) == 1) {
    // sort: (A',A,B',B) --> -1.0 * (A,B,A',B')
    SMITH::sort_indices<1,3,0,2,0,1,-1,1>(tmp.data(), out->data(), ApBp.template nstates<0>(), AB.template nstates<0>(), ApBp.template nstates<1>(), AB.template nstates<1>());
  }
  else {
    // sort: (A',A,B',B) --> (A,B,A',B')
    SMITH::sort_indices<1,3,0,2,0,1,1,1>(tmp.data(), out->data(), ApBp.template nstates<0>(), AB.template nstates<0>(), ApBp.template nstates<1>(), AB.template nstates<1>());
  }

  return out;
}


template <>
template <class VecType>
std::shared_ptr<Matrix> asd::ASD_impl<true>::compute_abFlip(MultiExcitonHamiltonian<VecType>* me, DSb<VecType>& AB, DSb<VecType>& ApBp) {
  auto out = std::make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());

  Matrix gamma_A = *me->gammaforest_->template get<0>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta);
  Matrix gamma_B = *me->gammaforest_->template get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::CreateAlpha);

  std::shared_ptr<const Matrix> Kmatrix = me->jop_->template coulomb_matrix<0,1,1,0>();

  Matrix tmp = gamma_A * (*Kmatrix) ^ gamma_B;

  // sort: (A',A,B',B) --> -1.0 * (A,B,A',B')
  SMITH::sort_indices<1,3,0,2,0,1,-1,1>(tmp.data(), out->data(), ApBp.template nstates<0>(), AB.template nstates<0>(), ApBp.template nstates<1>(), AB.template nstates<1>());

  return out;
}


template <>
template <class VecType>
std::shared_ptr<Matrix> asd::ASD_impl<true>::compute_abET(MultiExcitonHamiltonian<VecType>* me, DSb<VecType>& AB, DSb<VecType>& ApBp) {
  auto out = std::make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());

  Matrix gamma_A = *me->gammaforest_->template get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateBeta, GammaSQ::CreateAlpha);
  Matrix gamma_B = *me->gammaforest_->template get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha);

  std::shared_ptr<const Matrix> Jmatrix = me->jop_->template coulomb_matrix<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) ^ gamma_B;

  // sort: (A',A,B',B) --> -1.0 * (A,B,A',B')
  SMITH::sort_indices<1,3,0,2,0,1,-1,1>(tmp.data(), out->data(), ApBp.template nstates<0>(), AB.template nstates<0>(), ApBp.template nstates<1>(), AB.template nstates<1>());

  return out;
}


template <>
template <class VecType>
std::shared_ptr<Matrix> asd::ASD_impl<true>::compute_aaET(MultiExcitonHamiltonian<VecType>* me, DSb<VecType>& AB, DSb<VecType>& ApBp) {
  auto out = std::make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());

  Matrix gamma_A = *me->gammaforest_->template get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateAlpha, GammaSQ::CreateAlpha);
  Matrix gamma_B = *me->gammaforest_->template get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha);

  std::shared_ptr<const Matrix> Jmatrix = me->jop_->template coulomb_matrix<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) ^ gamma_B;

  // sort: (A',A,B',B) --> -0.5 * (A,B,A',B')
  SMITH::sort_indices<1,3,0,2,0,1,-1,2>(tmp.data(), out->data(), ApBp.template nstates<0>(), AB.template nstates<0>(), ApBp.template nstates<1>(), AB.template nstates<1>());

  return out;
}


template <>
template <class VecType>
std::shared_ptr<Matrix> asd::ASD_impl<true>::compute_bbET(MultiExcitonHamiltonian<VecType>* me, DSb<VecType>& AB, DSb<VecType>& ApBp) {
  auto out = std::make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());

  Matrix gamma_A = *me->gammaforest_->template get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateBeta, GammaSQ::CreateBeta);
  Matrix gamma_B = *me->gammaforest_->template get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta);

  std::shared_ptr<const Matrix> Jmatrix = me->jop_->template coulomb_matrix<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) ^ gamma_B;

  // sort: (A',A,B',B) --> -0.5 * (A,B,A',B')
  SMITH::sort_indices<1,3,0,2,0,1,-1,2>(tmp.data(), out->data(), ApBp.template nstates<0>(), AB.template nstates<0>(), ApBp.template nstates<1>(), AB.template nstates<1>());

  return out;
}

#endif

#endif
