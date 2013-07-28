//
// BAGEL - Parallel electron correlation program.
// Filename: meh_offdiagonal.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#include <src/meh/meh.h>

using namespace std;
using namespace bagel;

shared_ptr<Matrix> MultiExcitonHamiltonian::couple_blocks(DimerSubspace& AB, DimerSubspace& ApBp) {
  Coupling term_type = coupling_type(AB, ApBp);

  DimerSubspace* space1 = &AB;
  DimerSubspace* space2 = &ApBp;

  bool flip = (static_cast<int>(term_type) < 0);

  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    swap(space1,space2);
  }

  shared_ptr<Matrix> out;

  switch(term_type) {
    case Coupling::none :
      out = make_shared<Matrix>(space1->dimerstates(), space2->dimerstates()); break;
    case Coupling::diagonal :
      out = compute_inter_2e(*space1, *space2); break;
    case Coupling::aET :
      out = compute_aET(*space1, *space2); break;
    case Coupling::bET :
      out = compute_bET(*space1, *space2); break;
    case Coupling::abFlip :
      out = compute_abFlip(*space1, *space2); break;
    case Coupling::abET :
      out = compute_abET(*space1, *space2); break;
    case Coupling::aaET :
      out = compute_aaET(*space1, *space2); break;
    case Coupling::bbET :
      out = compute_bbET(*space1, *space2); break;
    default :
      throw logic_error("Asking for a coupling type that has not been written.");
  }

  if (flip) out = out->transpose();

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_aET(DimerSubspace& AB, DimerSubspace& ApBp) {
  auto out = make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());
  Matrix tmp(AB.nstates<0>() * ApBp.nstates<0>(), AB.nstates<1>() * ApBp.nstates<1>());

  // One-body aET
  {
    Matrix gamma_A = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateAlpha);
    Matrix gamma_B = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha);

    shared_ptr<Matrix> Fmatrix = jop_->cross_mo1e();

    tmp += gamma_A * (*Fmatrix) ^ gamma_B;
  }

  //Two-body aET, type 1
  {
    Matrix gamma_A = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateAlpha);
    Matrix gamma_B1 = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha);
    Matrix gamma_B2 = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta);

    shared_ptr<Matrix> Jmatrix = form_coulomb_matrix<0,1,1,1>();

    tmp -= gamma_A * (*Jmatrix) ^ (gamma_B1 + gamma_B2);
  }

  //Two-body aET, type 2
  {
    Matrix gamma_A1 = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha);
    Matrix gamma_A2 = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha);
    Matrix gamma_B = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha);

    shared_ptr<Matrix> Jmatrix = form_coulomb_matrix<0,0,1,0>();

    tmp += (gamma_A1 + gamma_A2) * (*Jmatrix) ^ gamma_B;
  }

  reorder_matrix(tmp.data(), out->data(), AB.nstates<0>(), ApBp.nstates<0>(), AB.nstates<1>(), ApBp.nstates<1>());

  const int neleA = AB.ci<0>()->det()->nelea() + AB.ci<0>()->det()->neleb();
  if ((neleA % 2) == 1) out->scale(-1.0);

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_bET(DimerSubspace& AB, DimerSubspace& ApBp) {
  auto out = make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());
  Matrix tmp(AB.nstates<0>() * ApBp.nstates<0>(), AB.nstates<1>() * ApBp.nstates<1>());

  // One-body bET
  {
    Matrix gamma_A = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateBeta);
    Matrix gamma_B = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta);

    shared_ptr<Matrix> Fmatrix = jop_->cross_mo1e();

    tmp += gamma_A * (*Fmatrix) ^ gamma_B;
  }


  //Two-body bET, type 1
  {
    Matrix gamma_A = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateBeta);
    Matrix gamma_B1 = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::CreateAlpha);
    Matrix gamma_B2 = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta);

    shared_ptr<Matrix> Jmatrix = form_coulomb_matrix<0,1,1,1>();

    tmp -= gamma_A * (*Jmatrix) ^ (gamma_B1 + gamma_B2);
  }

  //Two-body aET, type 2
  {
    Matrix gamma_A1 = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta);
    Matrix gamma_A2 = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta);
    Matrix gamma_B = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta);

    shared_ptr<Matrix> Jmatrix = form_coulomb_matrix<0,0,1,0>();

    tmp += (gamma_A1 + gamma_A2) * (*Jmatrix) ^ gamma_B;
  }

  reorder_matrix(tmp.data(), out->data(), AB.nstates<0>(), ApBp.nstates<0>(), AB.nstates<1>(), ApBp.nstates<1>());

  const int neleA = AB.ci<0>()->det()->nelea() + AB.ci<0>()->det()->neleb();
  if ((neleA % 2) == 1) out->scale(-1.0);

  return out;
}

// Currently defined as an alpha->beta flip in A and a beta->alpha flip in B
shared_ptr<Matrix> MultiExcitonHamiltonian::compute_abFlip(DimerSubspace& AB, DimerSubspace& ApBp) {
  auto out = make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());

  const int nstatesA = AB.nstates<0>();
  const int nstatesAp = ApBp.nstates<0>();
  const int nstatesB = AB.nstates<1>();
  const int nstatesBp = ApBp.nstates<1>();

  Matrix gamma_A = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta);
  Matrix gamma_B = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::CreateAlpha);

  shared_ptr<Matrix> Kmatrix = form_coulomb_matrix<0,1,1,0>();

  Matrix tmp = gamma_A * (*Kmatrix) ^ gamma_B;
  tmp *= -1.0;

  reorder_matrix(tmp.data(), out->data(), nstatesA, nstatesAp, nstatesB, nstatesBp);

  return out;
}


shared_ptr<Matrix> MultiExcitonHamiltonian::compute_abET(DimerSubspace& AB, DimerSubspace& ApBp) {
  auto out = make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());

  Matrix gamma_A = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateBeta, GammaSQ::CreateAlpha);
  Matrix gamma_B = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha);

  shared_ptr<Matrix> Jmatrix = form_coulomb_matrix<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) ^ gamma_B;
  tmp *= -1.0;

  reorder_matrix(tmp.data(), out->data(), AB.nstates<0>(), ApBp.nstates<0>(), AB.nstates<1>(), ApBp.nstates<1>());

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_aaET(DimerSubspace& AB, DimerSubspace& ApBp) {
  auto out = make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());

  Matrix gamma_A = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateAlpha, GammaSQ::CreateAlpha);
  Matrix gamma_B = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha);

  shared_ptr<Matrix> Jmatrix = form_coulomb_matrix<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) ^ gamma_B;
  tmp *= -0.5;

  reorder_matrix(tmp.data(), out->data(), AB.nstates<0>(), ApBp.nstates<0>(), AB.nstates<1>(), ApBp.nstates<1>());

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_bbET(DimerSubspace& AB, DimerSubspace& ApBp) {
  auto out = make_shared<Matrix>(AB.dimerstates(), ApBp.dimerstates());

  Matrix gamma_A = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::CreateBeta, GammaSQ::CreateBeta);
  Matrix gamma_B = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta);

  shared_ptr<Matrix> Jmatrix = form_coulomb_matrix<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) ^ gamma_B;
  tmp *= -0.5;

  reorder_matrix(tmp.data(), out->data(), AB.nstates<0>(), ApBp.nstates<0>(), AB.nstates<1>(), ApBp.nstates<1>());

  return out;
}
