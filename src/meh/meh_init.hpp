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
  MEH_base(input, dimer), cispace_(cispace)
{
  Timer timer;

  cispace_->complete();
  std::cout << "  o completing CI spin space: " << timer.tick() << std::endl;

  // Organize subspaces
  dimerstates_ = 0;
  int maxspin = 0;
  // Process DimerCISpace to form and organize needed Civecs and figure out maxspin
  for ( auto& aiter : cispace_->template cispace<0>() ) {
    SpaceKey akey = aiter.first;
    // values of S_b that could give the proper spin:
    //   S_a-nspin_, S_a-nspin+2,...,S_a+nspin
    for (int Sb = std::abs(akey.S-std::abs(nspin_)); Sb <= akey.S+std::abs(nspin_); Sb+=2) {
      SpaceKey bkey( Sb, nspin_ - akey.m_s, charge_ - akey.q );
      std::shared_ptr<const VecType> bspace = cispace_->template ccvec<1>(bkey);
      if ( bspace ) {
        subspaces_.emplace_back(dimerstates_, akey, bkey, make_pair(aiter.second, bspace));
        maxspin = std::max(akey.S+Sb, maxspin);
      }
    }
  }
  max_spin_ = maxspin + 1;

}

template <class VecType>
Coupling MultiExcitonHamiltonian<VecType>::coupling_type(const DSubSpace& AB, const DSubSpace& ApBp) const {
  std::array<MonomerKey,4> keys {{ AB.template monomerkey<0>(), AB.template monomerkey<1>(), ApBp.template monomerkey<0>(), ApBp.template monomerkey<1>()}};
  return MEH_base::coupling_type(keys);
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
void MultiExcitonHamiltonian<VecType>::print_hamiltonian(const std::string title, const int nstates) const {
  hamiltonian_->print(title, nstates);
}


template <class VecType>
void MultiExcitonHamiltonian<VecType>::print_states(const Matrix& cc, const std::vector<double>& energies, const double thresh, const std::string title) const {
  const int nstates = cc.mdim();
  std::shared_ptr<Matrix> spn = spin_->apply(cc);
  std::cout << std::endl << " ===== " << title << " =====" << std::endl;
  for (int istate = 0; istate < nstates; ++istate) {
    std::cout << "   state  " << std::setw(3) << istate << ": "
         << std::setprecision(8) << std::setw(17) << std::fixed << energies.at(istate)
         << "   <S^2> = " << std::setw(4) << std::setprecision(4) << std::fixed << ddot_(dimerstates_, spn->element_ptr(0,istate), 1, cc.element_ptr(0,istate), 1) << std::endl;
    const double *eigendata = cc.element_ptr(0,istate);
    double printed = 0.0;
    for (auto& subspace : subspaces_) {
      const int nA = subspace.template nstates<0>();
      const int nB = subspace.template nstates<1>();
      for (int i = 0; i < nA; ++i) {
        for (int j = 0; j < nB; ++j, ++eigendata) {
          if ( (*eigendata)*(*eigendata) > thresh ) {
            std::cout << "      " << subspace.string(i,j) << std::setprecision(12) << std::setw(20) << *eigendata << std::endl;
            printed += (*eigendata)*(*eigendata);
          }
        }
      }
    }
    std::cout << "    total weight of printed elements: " << std::setprecision(12) << std::setw(20) << printed << std::endl << std::endl;
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
void MultiExcitonHamiltonian<VecType>::print(const double thresh) const {
  print_states(*adiabats_, energies_, thresh, "Adiabatic States");
  if (dipoles_) {for (auto& prop : properties_) print_property(prop.first, prop.second, nstates_); }
}

#endif

#endif
