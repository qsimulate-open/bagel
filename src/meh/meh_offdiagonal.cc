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

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <vector>

#include <src/util/matrix.h>
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

  if (term_type == Coupling::none) { out = shared_ptr<Matrix>(new Matrix(space1->dimerstates(), space2->dimerstates())); }
  else if (term_type == Coupling::diagonal) { throw logic_error("Diagonal blocks are being called from the offdiagonal code. Why?"); }
  else if (term_type == Coupling::aET)      { out = compute_aET(*space1, *space2); }
  else if (term_type == Coupling::bET)      { out = compute_bET(*space1, *space2); }
  else if (term_type == Coupling::abFlip)   { out = compute_abFlip(*space1, *space2); }
  else if (term_type == Coupling::abET)     { out = compute_abET(*space1, *space2); }
  else if (term_type == Coupling::aaET)     { out = compute_aaET(*space1, *space2); }
  else if (term_type == Coupling::bbET)     { out = compute_bbET(*space1, *space2); }
  else { throw logic_error("Blocks with an unknown relationship are being coupled"); }

  if (flip) out = out->transpose();

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_aET(DimerSubspace& AB, DimerSubspace& ApBp) {
  shared_ptr<Matrix> out(new Matrix(AB.dimerstates(), ApBp.dimerstates()));

  // One-body aET
  {

  }

  Matrix twobody(AB.nstates<0>() * ApBp.nstates<0>(), AB.nstates<1>() * ApBp.nstates<1>());

  //Two-body aET, type 1
  {
    shared_ptr<Quantization> one(new OneBody<SQ::CreateAlpha>());
    shared_ptr<Quantization> three_alpha(new ThreeBody<SQ::CreateAlpha,SQ::AnnihilateAlpha,SQ::AnnihilateAlpha>());
    shared_ptr<Quantization> three_beta(new ThreeBody<SQ::CreateBeta,SQ::AnnihilateBeta,SQ::AnnihilateAlpha>());

    Matrix gamma_A = *form_gamma(AB.ci<0>(), ApBp.ci<0>(), one);
    Matrix gamma_B = (*form_gamma(AB.ci<1>(), ApBp.ci<1>(), three_alpha)) + (*form_gamma(AB.ci<1>(), ApBp.ci<1>(), three_beta));

    shared_ptr<Matrix> Jmatrix, Kmatrix;
    tie(Jmatrix,Kmatrix) = form_JKmatrices<0,1,1,1>();

    twobody += gamma_A * (*Jmatrix) ^ gamma_B;
  }

  //Two-body aET, type 2
  {
    shared_ptr<Quantization> one(new OneBody<SQ::AnnihilateAlpha>());
    shared_ptr<Quantization> three_alpha(new ThreeBody<SQ::CreateAlpha,SQ::CreateAlpha,SQ::AnnihilateAlpha>());
    shared_ptr<Quantization> three_beta(new ThreeBody<SQ::CreateAlpha,SQ::CreateBeta,SQ::AnnihilateBeta>());

    Matrix gamma_A = (*form_gamma(AB.ci<0>(), ApBp.ci<0>(), three_alpha)) + (*form_gamma(AB.ci<0>(), ApBp.ci<0>(), three_beta));
    Matrix gamma_B = *form_gamma(AB.ci<1>(), ApBp.ci<1>(), one);

    shared_ptr<Matrix> Jmatrix, Kmatrix;
    tie(Jmatrix, Kmatrix) = form_JKmatrices<0,0,1,0>();

    twobody += gamma_A * (*Jmatrix) ^ gamma_B;
  }

  Matrix tmp_twobody(AB.dimerstates(), ApBp.dimerstates());
  reorder_matrix(twobody.data(), tmp_twobody.data(), AB.nstates<0>(), ApBp.nstates<0>(), AB.nstates<1>(), ApBp.nstates<1>());
  *out += tmp_twobody;

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_bET(DimerSubspace& AB, DimerSubspace& ApBp) {
  shared_ptr<Matrix> out(new Matrix(AB.dimerstates(), ApBp.dimerstates()));

  // One-body bET
  {

  }

  Matrix twobody(AB.nstates<0>() * ApBp.nstates<0>(), AB.nstates<1>() * ApBp.nstates<1>());

  //Two-body bET, type 1
  {
    shared_ptr<Quantization> one(new OneBody<SQ::CreateBeta>());
    shared_ptr<Quantization> three_alpha(new ThreeBody<SQ::CreateAlpha,SQ::AnnihilateAlpha,SQ::AnnihilateBeta>());
    shared_ptr<Quantization> three_beta(new ThreeBody<SQ::CreateBeta,SQ::AnnihilateBeta,SQ::AnnihilateBeta>());

    Matrix gamma_A = *form_gamma(AB.ci<0>(), ApBp.ci<0>(), one);
    Matrix gamma_B = (*form_gamma(AB.ci<1>(), ApBp.ci<1>(), three_alpha)) + (*form_gamma(AB.ci<1>(), ApBp.ci<1>(), three_beta));

    shared_ptr<Matrix> Jmatrix, Kmatrix;
    tie(Jmatrix,Kmatrix) = form_JKmatrices<0,1,1,1>();

    twobody += gamma_A * (*Jmatrix) ^ gamma_B;
  }

  //Two-body aET, type 2
  {
    shared_ptr<Quantization> one(new OneBody<SQ::AnnihilateBeta>());
    shared_ptr<Quantization> three_alpha(new ThreeBody<SQ::CreateBeta,SQ::CreateAlpha,SQ::AnnihilateAlpha>());
    shared_ptr<Quantization> three_beta(new ThreeBody<SQ::CreateBeta,SQ::CreateBeta,SQ::AnnihilateBeta>());

    Matrix gamma_A = (*form_gamma(AB.ci<0>(), ApBp.ci<0>(), three_alpha)) + (*form_gamma(AB.ci<0>(), ApBp.ci<0>(), three_beta));
    Matrix gamma_B = *form_gamma(AB.ci<1>(), ApBp.ci<1>(), one);

    shared_ptr<Matrix> Jmatrix, Kmatrix;
    tie(Jmatrix, Kmatrix) = form_JKmatrices<0,0,1,0>();

    twobody += gamma_A * (*Jmatrix) ^ gamma_B;
  }

  Matrix tmp_twobody(AB.dimerstates(), ApBp.dimerstates());
  reorder_matrix(twobody.data(), tmp_twobody.data(), AB.nstates<0>(), ApBp.nstates<0>(), AB.nstates<1>(), ApBp.nstates<1>());
  *out += tmp_twobody;

  return out;
}

// Currently defined as an alpha->beta flip in A and a beta->alpha flip in B
shared_ptr<Matrix> MultiExcitonHamiltonian::compute_abFlip(DimerSubspace& AB, DimerSubspace& ApBp) {
  shared_ptr<Matrix> out(new Matrix(AB.dimerstates(), ApBp.dimerstates()));

  const int nstatesA = AB.nstates<0>();
  const int nstatesAp = ApBp.nstates<0>();
  const int nstatesB = AB.nstates<1>();
  const int nstatesBp = ApBp.nstates<1>();

  shared_ptr<Quantization> ab_oper(new TwoBody<SQ::CreateBeta,SQ::AnnihilateAlpha>());
  shared_ptr<Quantization> ba_oper(new TwoBody<SQ::CreateAlpha,SQ::AnnihilateBeta>());

  Matrix gamma_A = *form_gamma(AB.ci<0>(), ApBp.ci<0>(), ab_oper);
  Matrix gamma_B = *form_gamma(AB.ci<1>(), ApBp.ci<1>(), ba_oper);
  shared_ptr<Matrix> Jmatrix, Kmatrix;
  tie(Jmatrix,Kmatrix) = form_JKmatrices<0,1,0,1>();

  Matrix tmp = gamma_A * (*Kmatrix) ^ gamma_B;
  tmp *= -1.0;

  reorder_matrix(tmp.data(), out->data(), nstatesA, nstatesAp, nstatesB, nstatesBp);

  return out;
}


shared_ptr<Matrix> MultiExcitonHamiltonian::compute_abET(DimerSubspace& AB, DimerSubspace& ApBp) {
  shared_ptr<Matrix> out(new Matrix(AB.dimerstates(), ApBp.dimerstates()));

  shared_ptr<Quantization> creation(new TwoBody<SQ::CreateAlpha,SQ::CreateBeta>());
  shared_ptr<Quantization> annihilation(new TwoBody<SQ::AnnihilateBeta,SQ::AnnihilateAlpha>());

  Matrix gamma_A = *form_gamma(AB.ci<0>(), ApBp.ci<0>(), creation);
  Matrix gamma_B = *form_gamma(AB.ci<1>(), ApBp.ci<1>(), annihilation);
  shared_ptr<Matrix> Jmatrix, Kmatrix;
  tie(Jmatrix,Kmatrix) = form_JKmatrices<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) * gamma_B;

  reorder_matrix(tmp.data(), out->data(), AB.nstates<0>(), ApBp.nstates<0>(), AB.nstates<1>(), ApBp.nstates<1>());

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_aaET(DimerSubspace& AB, DimerSubspace& ApBp) {
  shared_ptr<Matrix> out(new Matrix(AB.dimerstates(), ApBp.dimerstates()));

  shared_ptr<Quantization> creation(new TwoBody<SQ::CreateAlpha,SQ::CreateAlpha>());
  shared_ptr<Quantization> annihilation(new TwoBody<SQ::AnnihilateAlpha,SQ::AnnihilateAlpha>());

  Matrix gamma_A = *form_gamma(AB.ci<0>(), ApBp.ci<0>(), creation);
  Matrix gamma_B = *form_gamma(AB.ci<1>(), ApBp.ci<1>(), annihilation);
  shared_ptr<Matrix> Jmatrix, Kmatrix;
  tie(Jmatrix,Kmatrix) = form_JKmatrices<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) * gamma_B;
  tmp *= 0.5;

  reorder_matrix(tmp.data(), out->data(), AB.nstates<0>(), ApBp.nstates<0>(), AB.nstates<1>(), ApBp.nstates<1>());

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_bbET(DimerSubspace& AB, DimerSubspace& ApBp) {
  shared_ptr<Matrix> out(new Matrix(AB.dimerstates(), ApBp.dimerstates()));

  shared_ptr<Quantization> creation(new TwoBody<SQ::CreateBeta,SQ::CreateBeta>());
  shared_ptr<Quantization> annihilation(new TwoBody<SQ::AnnihilateBeta,SQ::AnnihilateBeta>());

  Matrix gamma_A = *form_gamma(AB.ci<0>(), ApBp.ci<0>(), creation);
  Matrix gamma_B = *form_gamma(AB.ci<1>(), ApBp.ci<1>(), annihilation);
  shared_ptr<Matrix> Jmatrix, Kmatrix;
  tie(Jmatrix,Kmatrix) = form_JKmatrices<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) * gamma_B;
  tmp *= 0.5;

  reorder_matrix(tmp.data(), out->data(), AB.nstates<0>(), ApBp.nstates<0>(), AB.nstates<1>(), ApBp.nstates<1>());

  return out;
}
