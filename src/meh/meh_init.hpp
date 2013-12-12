//
// BAGEL - Parallel electron correlation program.
// Filename: meh_init.h
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

#ifndef BAGEL_MEH_INIT_H
#define BAGEL_MEH_INIT_H

template <class VecType>
MultiExcitonHamiltonian<VecType>::MultiExcitonHamiltonian(const std::shared_ptr<const PTree> input, std::shared_ptr<Dimer> dimer, std::shared_ptr<DimerCISpace_base<VecType>> cispace) :
  ref_(dimer->sref()), coeff_(dimer->scoeff()), cispace_(cispace),
  dimerbasis_(dimer->dimerbasis()), dimerclosed_(dimer->sref()->nclosed()), dimeractive_(dimer->sref()->nact()),
  nact_(dimer->nact()), nbasis_(dimer->nbasis())
{
  nstates_ = input->get<int>("nstates", 10);
  max_iter_ = input->get<int>("max_iter", 50);
  davidsonceiling_ = input->get<int>("davidsonceiling", 10);
  dipoles_ = input->get<bool>("dipoles", false);
  thresh_ = input->get<double>("thresh", 1.0e-12);
  print_thresh_ = input->get<double>("print_thresh", 0.05);
  store_matrix_ = input->get<bool>("store_matrix", false);
  nspin_ = 0; // TODO hardcoded to singlets for now

  jop_ = std::make_shared<DimerJop>(ref_, dimerclosed_, dimerclosed_ + nact_.first, dimerclosed_ + dimeractive_, coeff_);

  cispace_->complete();

  dimerstates_ = 0;

  int maxspin = 0;
  // Process DimerCISpace to form and organize needed Civecs and figure out maxspin
  for ( auto& aiter : cispace_->template cispace<0>() ) {
    SpaceKey akey = aiter.first;
    SpaceKey bkey( akey.S, -akey.m_s, -akey.q );
    std::shared_ptr<const VecType> bspace = cispace_->template ccvec<1>(bkey);
    if ( bspace ) {
      subspaces_.emplace_back(dimerstates_, akey, bkey, make_pair(aiter.second, bspace));
      maxspin = std::max(aiter.second->det()->nspin(), maxspin);
    }
  }
  maxspin = 2 * maxspin + 1;
  max_spin_ = maxspin;

  energies_ = std::vector<double>(nstates_, 0.0);

  gammaforest_ = std::make_shared<GammaForest<VecType, 2>>();
}

template <class VecType>
const Coupling MultiExcitonHamiltonian<VecType>::coupling_type(const DSubSpace& AB, const DSubSpace& ApBp) const {
  std::pair<int,int> neleaAB = std::make_pair(AB.template ci<0>()->det()->nelea(), AB.template ci<1>()->det()->nelea());
  std::pair<int,int> nelebAB = std::make_pair(AB.template ci<0>()->det()->neleb(), AB.template ci<1>()->det()->neleb());

  std::pair<int,int> neleaApBp = std::make_pair(ApBp.template ci<0>()->det()->nelea(), ApBp.template ci<1>()->det()->nelea());
  std::pair<int,int> nelebApBp = std::make_pair(ApBp.template ci<0>()->det()->neleb(), ApBp.template ci<1>()->det()->neleb());

  // AlphaTransfer and BetaTransfer
  std::pair<int,int> AT = std::make_pair(neleaApBp.first - neleaAB.first, neleaApBp.second - neleaAB.second);
  std::pair<int,int> BT = std::make_pair(nelebApBp.first - nelebAB.first, nelebApBp.second - nelebAB.second);

  const int stride = 8; // Should be sufficient
  auto coupling_index = [&stride] (const int a, const int b, const int c, const int d) { return a + b * stride + stride*stride * (c + d * stride); };

  /************************************************************
  *  BT\AT  | ( 0, 0) | (+1,-1) | (-1,+1) | (+2,-2) | (-2,+2) *
  *-----------------------------------------------------------*
  * ( 0, 0) |  diag   |  aET    |  -aET   |  aaET   | -aaET   *
  * (+1,-1) |  bET    |  dABT   |  ABflp  |         |         *
  * (-1,+1) | -bET    | BAflp   | -dABT   |         |         *
  * (+2,-2) |  bbET   |         |         |         |         *
  * (-2,+2) | -bbET   |         |         |         |         *
  ************************************************************/

  const int icouple = coupling_index(AT.first, AT.second, BT.first, BT.second);

  if      ( icouple == coupling_index( 0, 0, 0, 0) ) return Coupling::diagonal;
  else if ( icouple == coupling_index( 0, 0,+1,-1) ) return Coupling::bET;
  else if ( icouple == coupling_index( 0, 0,-1,+1) ) return Coupling::inv_bET;
  else if ( icouple == coupling_index(+1,-1, 0, 0) ) return Coupling::aET;
  else if ( icouple == coupling_index(+1,-1,+1,-1) ) return Coupling::abET;
  else if ( icouple == coupling_index(+1,-1,-1,+1) ) return Coupling::baFlip;
  else if ( icouple == coupling_index(-1,+1, 0, 0) ) return Coupling::inv_aET;
  else if ( icouple == coupling_index(-1,+1,+1,-1) ) return Coupling::abFlip;
  else if ( icouple == coupling_index(-1,+1,-1,+1) ) return Coupling::inv_abET;
  else if ( icouple == coupling_index(+2,-2, 0, 0) ) return Coupling::aaET;
  else if ( icouple == coupling_index(-2,+2, 0, 0) ) return Coupling::inv_aaET;
  else if ( icouple == coupling_index( 0, 0,+2,-2) ) return Coupling::bbET;
  else if ( icouple == coupling_index( 0, 0,-2,+2) ) return Coupling::inv_bbET;
  else                                               return Coupling::none;
}

template <class VecType>
std::shared_ptr<const Matrix> MultiExcitonHamiltonian<VecType>::apply_hamiltonian(const Matrix& o) {
  if (store_matrix_) {
    return std::make_shared<const Matrix>(*hamiltonian_ * o);
  }
  else {
    auto out = std::make_shared<Matrix>(o.ndim(), o.mdim());
    for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
      const int ioff = iAB->offset();
      for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
        const int joff = jAB->offset();

        std::shared_ptr<const Matrix> block = couple_blocks(*iAB, *jAB);

        dgemm_("N", "N", block->ndim(), o.mdim(), block->mdim(), 1.0, block->data(), block->ndim(), o.element_ptr(joff, 0), dimerstates_, 1.0, out->element_ptr(ioff, 0), o.ndim());
        dgemm_("T", "N", block->mdim(), o.mdim(), block->ndim(), 1.0, block->data(), block->ndim(), o.element_ptr(ioff, 0), dimerstates_, 1.0, out->element_ptr(joff, 0), o.ndim());
      }

      std::shared_ptr<const Matrix> block = compute_diagonal_block(*iAB);
      dgemm_("N", "N", block->ndim(), o.mdim(), block->mdim(), 1.0, block->data(), block->ndim(), o.element_ptr(ioff, 0), dimerstates_, 1.0, out->element_ptr(ioff, 0), o.ndim());
    }

    return out;
  }
}


template <class VecType>
std::shared_ptr<Matrix> MultiExcitonHamiltonian<VecType>::compute_1e_prop(std::shared_ptr<const Matrix> hAA, std::shared_ptr<const Matrix> hBB, std::shared_ptr<const Matrix> hAB, const double core) const {

  auto out = std::make_shared<Matrix>(dimerstates_, dimerstates_);

  for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
    const int ioff = iAB->offset();
    for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
      const int joff = jAB->offset();

      std::shared_ptr<Matrix> out_block = compute_offdiagonal_1e(*iAB, *jAB, hAB);

      out->add_block(1.0, ioff, joff, out_block->ndim(), out_block->mdim(), out_block);
      out->add_block(1.0, joff, ioff, out_block->mdim(), out_block->ndim(), out_block->transpose());
    }
    std::shared_ptr<const Matrix> tmp = compute_diagonal_1e(*iAB, hAA->data(), hBB->data(), core);
    out->add_block(1.0, ioff, ioff, tmp->ndim(), tmp->mdim(), tmp);
  }

  return out;
}


template <class VecType>
void MultiExcitonHamiltonian<VecType>::reorder_matrix(const double* source, double* target,
  const int nA, const int nAp, const int nB, const int nBp) const
{
  const int nstatesAB = nA * nB;
  const int nstatesAA = nA * nAp;

  for (int Bp = 0; Bp < nBp; ++Bp) {
    for (int Ap = 0; Ap < nAp; ++Ap) {
      const int ABp = Ap + nAp * Bp;
      for (int B = 0; B < nB; ++B) {
        const int BBp = Bp + nBp * B;
        for (int A = 0; A < nA; ++A) {
          const int AAp = Ap + nAp * A;
          const int AB = A + nA * B;
          target[AB + nstatesAB * ABp] = source[AAp + nstatesAA * BBp];
        }
      }
    }
  }
}


template <class VecType>
void MultiExcitonHamiltonian<VecType>::print_hamiltonian(const std::string title, const int nstates) const {
  hamiltonian_->print(title, nstates);
}


template <class VecType>
void MultiExcitonHamiltonian<VecType>::print_adiabats(const double thresh, const std::string title, const int nstates) const {
  const int end = std::min(nstates, adiabats_->mdim());
  std::shared_ptr<Matrix> spn = spin_->apply(*adiabats_);
  std::cout << std::endl << " ===== " << title << " =====" << std::endl;
  for (int istate = 0; istate < end; ++istate) {
    std::cout << "   state  " << std::setw(3) << istate << ": "
         << std::setprecision(8) << std::setw(17) << std::fixed << energies_.at(istate)
         << "   <S^2> = " << std::setw(4) << std::setprecision(4) << std::fixed << ddot_(dimerstates_, spn->element_ptr(0,istate), 1, adiabats_->element_ptr(0,istate), 1) << std::endl;
    double *eigendata = adiabats_->element_ptr(0,istate);
    for (auto& subspace : subspaces_) {
      const int nA = subspace.template nstates<0>();
      const int nB = subspace.template nstates<1>();
      for (int i = 0; i < nA; ++i) {
        for (int j = 0; j < nB; ++j, ++eigendata) {
          if ( (*eigendata)*(*eigendata) > thresh ) {
            std::cout << "      " << subspace.string(i,j) << std::setprecision(12) << std::setw(20) << *eigendata << std::endl;
          }
        }
      }
    }
    std::cout << std::endl;
  }
}

template <class VecType>
void MultiExcitonHamiltonian<VecType>::print_property(const std::string label, std::shared_ptr<const Matrix> property , const int nstates) const {
  const std::string indent("   ");
  const int nprint = std::min(nstates, property->ndim());

  std::cout << indent << " " << label << "    |0>";
  for (int istate = 1; istate < nprint; ++istate) std::cout << "         |" << istate << ">";
  std::cout << std::endl;
  for (int istate = 0; istate < nprint; ++istate) {
    std::cout << indent << "<" << istate << "|";
    for (int jstate = 0; jstate < nprint; ++jstate) {
      std::cout << std::setw(12) << std::setprecision(6) << property->element(jstate, istate);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template <class VecType>
void MultiExcitonHamiltonian<VecType>::print(const int nstates, const double thresh) const {
  print_adiabats(thresh, "Adiabatic States", nstates);
  if (dipoles_) {for (auto& prop : properties_) print_property(prop.first, prop.second, nstates); }
}

#endif

#endif
