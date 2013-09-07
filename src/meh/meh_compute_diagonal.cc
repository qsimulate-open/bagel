//
// BAGEL - Parallel electron correlation program.
// Filename: meh_compute_diagonal.cc
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

#include <src/meh/meh.h>

using namespace std;
using namespace bagel;

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_intra(const DimerSubspace& AB, shared_ptr<const DimerJop> jop, const double diag) {
  shared_ptr<const Dvec> ccvecA = AB.ci<0>();
  shared_ptr<const Dvec> ccvecB = AB.ci<1>();
  shared_ptr<const Dvec> sigmavecA = form_sigma(ccvecA, jop->mo1e_first(), jop->mo2e_first());
  shared_ptr<const Dvec> sigmavecB = form_sigma(ccvecB, jop->mo1e_second(), jop->mo2e_second());

  const int nstatesA = AB.nstates<0>();
  const int nstatesB = AB.nstates<1>();
  const int dimerstates = AB.dimerstates();

  auto out = make_shared<Matrix>(dimerstates, dimerstates);

  // first H^{AA}_{AA}
  for(int stateA = 0; stateA < nstatesA; ++stateA) {
    shared_ptr<const Civec> isigma = sigmavecA->data(stateA);
    for(int stateAp = 0; stateAp < stateA; ++stateAp) {
      shared_ptr<const Civec> icc = ccvecA->data(stateAp);
      const double dotproduct = isigma->dot_product(*icc);
      for(int stateB = 0; stateB < nstatesB; ++stateB) {
        const int stateApB = AB.dimerindex(stateAp, stateB);
        const int stateAB = AB.dimerindex(stateA, stateB);
        out->element(stateAB, stateApB) += dotproduct;
        out->element(stateApB, stateAB) += dotproduct;
      }
    }
    const double dotproduct = ccvecA->data(stateA)->dot_product(*isigma);
    for(int stateB = 0; stateB < nstatesB; ++stateB) {
      const int stateAB = AB.dimerindex(stateA, stateB);
      out->element(stateAB,stateAB) += dotproduct;
    }
  }

  // H^{BB}_{BB}
  for(int stateB = 0; stateB < nstatesB; ++stateB) {
    shared_ptr<const Civec> isigma = sigmavecB->data(stateB);
    for(int stateBp = 0; stateBp < stateB; ++stateBp) {
      shared_ptr<const Civec> icc = ccvecB->data(stateBp);
      const double dotproduct = isigma->dot_product(*icc);
      for(int stateA = 0; stateA < nstatesA; ++stateA) {
        const int stateAB = AB.dimerindex(stateA, stateB);
        const int stateABp = AB.dimerindex(stateA, stateBp);
        out->element(stateAB, stateABp) += dotproduct;
        out->element(stateABp, stateAB) += dotproduct;
      }
    }
    const double dotproduct = ccvecB->data(stateB)->dot_product(*isigma);
    for(int stateA = 0; stateA < nstatesA; ++stateA) {
      const int stateAB = AB.dimerindex(stateA, stateB);
      out->element(stateAB,stateAB) += dotproduct;
    }
  }

  out->add_diag(diag);

  return out;
}

shared_ptr<Matrix> MultiExcitonHamiltonian::compute_diagonal_1e(const DimerSubspace& AB, const double* hAA, const double* hBB, const double diag) const {
  shared_ptr<const Dvec> ccvecA = AB.ci<0>();
  shared_ptr<const Dvec> ccvecB = AB.ci<1>();
  shared_ptr<const Dvec> sigmavecA = form_sigma_1e(ccvecA, hAA);
  shared_ptr<const Dvec> sigmavecB = form_sigma_1e(ccvecB, hBB);

  const int nstatesA = AB.nstates<0>();
  const int nstatesB = AB.nstates<1>();
  const int dimerstates = AB.dimerstates();

  auto out = make_shared<Matrix>(dimerstates, dimerstates);

  // first H^{AA}_{AA}
  for(int stateA = 0; stateA < nstatesA; ++stateA) {
    shared_ptr<const Civec> isigma = sigmavecA->data(stateA);
    for(int stateAp = 0; stateAp < stateA; ++stateAp) {
      shared_ptr<const Civec> icc = ccvecA->data(stateAp);
      const double dotproduct = isigma->dot_product(*icc);
      for(int stateB = 0; stateB < nstatesB; ++stateB) {
        const int stateApB = AB.dimerindex(stateAp, stateB);
        const int stateAB = AB.dimerindex(stateA, stateB);
        out->element(stateAB, stateApB) += dotproduct;
        out->element(stateApB, stateAB) += dotproduct;
      }
    }
    const double dotproduct = isigma->dot_product(*ccvecA->data(stateA));
    for(int stateB = 0; stateB < nstatesB; ++stateB) {
      const int stateAB = AB.dimerindex(stateA, stateB);
      out->element(stateAB,stateAB) += dotproduct;
    }
  }

  // H^{BB}_{BB}
  for(int stateB = 0; stateB < nstatesB; ++stateB) {
    shared_ptr<const Civec> isigma = sigmavecB->data(stateB);
    for(int stateBp = 0; stateBp < stateB; ++stateBp) {
      shared_ptr<const Civec> icc = ccvecB->data(stateBp);
      const double dotproduct = isigma->dot_product(*icc);
      for(int stateA = 0; stateA < nstatesA; ++stateA) {
        const int stateAB = AB.dimerindex(stateA, stateB);
        const int stateABp = AB.dimerindex(stateA, stateBp);
        out->element(stateAB, stateABp) += dotproduct;
        out->element(stateABp, stateAB) += dotproduct;
      }
    }
    const double dotproduct = ccvecB->data(stateB)->dot_product(*isigma);
    for(int stateA = 0; stateA < nstatesA; ++stateA) {
      const int stateAB = AB.dimerindex(stateA, stateB);
      out->element(stateAB,stateAB) += dotproduct;
    }
  }

  out->add_diag(diag);

  return out;
}


shared_ptr<Matrix> MultiExcitonHamiltonian::compute_diagonal_block(DimerSubspace& subspace) {
  const double core = ref_->geom()->nuclear_repulsion() + jop_->core_energy();

  // Would be better to allocate here and then send to subprocesses
  shared_ptr<Matrix> out = compute_intra(subspace, jop_, core);
  *out += *compute_inter_2e(subspace, subspace);

  return out;
}

// This term will couple off-diagonal blocks since it has no delta functions involved
shared_ptr<Matrix> MultiExcitonHamiltonian::compute_inter_2e(DimerSubspace& AB, DimerSubspace& ApBp) {
  shared_ptr<const Dvec> ccvecA = AB.ci<0>();
  shared_ptr<const Dvec> ccvecB = AB.ci<1>();
  shared_ptr<const Dvec> ccvecAp = ApBp.ci<0>();
  shared_ptr<const Dvec> ccvecBp = ApBp.ci<1>();

  const int nstatesA = AB.nstates<0>();
  const int nstatesB = AB.nstates<1>();
  const int nstates = nstatesA * nstatesB;

  const int nstatesAp = ApBp.nstates<0>();
  const int nstatesBp = ApBp.nstates<1>();
  const int nstatesp = nstatesAp * nstatesBp;

  // alpha-alpha
  Matrix gamma_AA_alpha = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha);
  Matrix gamma_BB_alpha = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha);

  // beta-beta
  Matrix gamma_AA_beta = *gammaforest_->get<0>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta);
  Matrix gamma_BB_beta = *gammaforest_->get<1>(AB.offset(), ApBp.offset(), GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta);

  // build J and K matrices
  shared_ptr<const Matrix> Jmatrix = jop_->coulomb_matrix<0,1,0,1>();
  shared_ptr<const Matrix> Kmatrix = jop_->coulomb_matrix<0,1,1,0>();

  Matrix tmp(nstatesA*nstatesAp, nstatesB*nstatesBp);

  tmp += (gamma_AA_alpha + gamma_AA_beta) * (*Jmatrix) ^ (gamma_BB_alpha + gamma_BB_beta);

  tmp -= gamma_AA_alpha * (*Kmatrix) ^ gamma_BB_alpha;
  tmp -= gamma_AA_beta * (*Kmatrix) ^ gamma_BB_beta;

  auto out = make_shared<Matrix>(nstates, nstatesp);
  reorder_matrix(tmp.data(), out->data(), nstatesA, nstatesAp, nstatesB, nstatesBp);

  return out;
}
