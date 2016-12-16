//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_compute_diagonal.hpp
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifdef ASD_HEADERS

#ifndef BAGEL_ASD_COMPUTE_DIAGONAL_H
#define BAGEL_ASD_COMPUTE_DIAGONAL_H

template <class VecType>
void ASD<VecType>::compute_pure_terms(DSubSpace& AB, std::shared_ptr<const DimerJop> jop) {
  {
    std::shared_ptr<const VecType> ccvecA = AB.template ci<0>();
    std::shared_ptr<const VecType> sigmavecA = form_sigma(ccvecA, jop->monomer_jop<0>());

    const int nstatesA = AB.template nstates<0>();

    auto Asigma = std::make_shared<CSymMatrix>(nstatesA, true);
    for (int i = 0; i < nstatesA; ++i) {
      for (int j = 0; j < i; ++j) {
        Asigma->element(j,i) = ccvecA->data(i)->dot_product(*sigmavecA->data(j));
      }
      Asigma->element(i,i) = ccvecA->data(i)->dot_product(*sigmavecA->data(i));
    }
    AB.template set_sigma<0>(Asigma);
  }

  {
    std::shared_ptr<const VecType> ccvecB = AB.template ci<1>();
    std::shared_ptr<const VecType> sigmavecB = form_sigma(ccvecB, jop->monomer_jop<1>());

    const int nstatesB = AB.template nstates<1>();

    auto Bsigma = std::make_shared<CSymMatrix>(nstatesB, true);
    for (int i = 0; i < nstatesB; ++i) {
      for (int j = 0; j < i; ++j) {
        Bsigma->element(j,i) = ccvecB->data(i)->dot_product(*sigmavecB->data(j));
      }
        Bsigma->element(i,i) = ccvecB->data(i)->dot_product(*sigmavecB->data(i));
    }
    AB.template set_sigma<1>(Bsigma);
  }
}


template <class VecType>
std::shared_ptr<Matrix> ASD<VecType>::compute_diagonal_1e(const DSubSpace& AB, const double* hAA, const double* hBB, const double diag) const {
  std::shared_ptr<const VecType> ccvecA = AB.template ci<0>();
  std::shared_ptr<const VecType> ccvecB = AB.template ci<1>();
  std::shared_ptr<const VecType> sigmavecA = form_sigma_1e(ccvecA, hAA);
  std::shared_ptr<const VecType> sigmavecB = form_sigma_1e(ccvecB, hBB);

  const int nstatesA = AB.template nstates<0>();
  const int nstatesB = AB.template nstates<1>();
  const int dimerstates = AB.dimerstates();

  auto out = std::make_shared<Matrix>(dimerstates, dimerstates);

  // first H^{AA}_{AA}
  for(int stateA = 0; stateA < nstatesA; ++stateA) {
    std::shared_ptr<const CiType> isigma = sigmavecA->data(stateA);
    for(int stateAp = 0; stateAp < stateA; ++stateAp) {
      std::shared_ptr<const CiType> icc = ccvecA->data(stateAp);
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
    std::shared_ptr<const CiType> isigma = sigmavecB->data(stateB);
    for(int stateBp = 0; stateBp < stateB; ++stateBp) {
      std::shared_ptr<const CiType> icc = ccvecB->data(stateBp);
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


#endif

#endif

