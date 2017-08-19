//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_compute.hpp
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

#ifndef BAGEL_ASD_COMPUTE_H
#define BAGEL_ASD_COMPUTE_H

template <class VecType>
void ASD<VecType>::compute() {
  Timer asdtime;
  std::cout << std::endl << " ===== Starting construction of dimer Hamiltonian " << std::endl;
  std::cout << "   o Dimer space:" << std::endl;
  std::cout << "     -  spin: " << nspin_ << std::endl;
  std::cout << "     -  charge: " << charge_ << std::endl;
  std::cout << "     -  dimer states: " << dimerstates_ << std::endl << std::endl;

  if (!fix_ci_) {
    auto gammaforest = std::make_shared<GammaForest<VecType, 2>>();
    {
      auto spinmap = std::make_shared<ASDSpinMap<VecType>>();
      for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
        for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
          gammaforest->couple_blocks(*jAB, *iAB);
          spinmap->couple_blocks(*jAB, *iAB);
        }
        gammaforest->couple_blocks(*iAB, *iAB);
        spinmap->diagonal_block(*iAB);
      }
      spin_ = std::make_shared<ASDSpin>(dimerstates_, *spinmap, max_spin_);
    }

    std::cout << "  o Preparing Gamma trees and building spin operator - " << std::setw(9) << std::fixed << std::setprecision(2) << asdtime.tick() << std::endl;
    std::cout << "    - spin elements: " << spin_->size() << std::endl;

    gammaforest->compute();
    gammatensor_ = { std::make_shared<GammaTensor>(asd::Wrap<GammaForest<VecType,2>,0>(gammaforest), subspaces_),
                     std::make_shared<GammaTensor>(asd::Wrap<GammaForest<VecType,2>,1>(gammaforest), subspaces_) };

    std::cout << "  o Computing Gamma trees - " << std::setw(9) << std::fixed << std::setprecision(2) << asdtime.tick() << std::endl;
  } else {
    std::cout << "  o Monomer CI coefficients are fixed. Gamma trees from previous calculation will be used." << std::endl;
  }

  if (store_matrix_) hamiltonian_ = std::make_shared<Matrix>(dimerstates_, dimerstates_);

  denom_ = std::unique_ptr<double[]>(new double[dimerstates_]);

  for (auto& subspace : subspaces_) {
    compute_pure_terms(subspace, jop_);
    std::shared_ptr<Matrix> block = compute_diagonal_block(subspace);
    if (store_matrix_)
      hamiltonian_->add_block(1.0, subspace.offset(), subspace.offset(), block->ndim(), block->mdim(), block);
    const int n = block->ndim();
    for (int i = 0; i < n; ++i) denom_[subspace.offset() + i] = block->element(i,i);
  }
  std::cout << "  o Computing diagonal blocks and building denominator - time " << std::setw(9) << std::fixed << std::setprecision(2) << asdtime.tick() << std::endl;

  if (store_matrix_) {
    for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
      const int ioff = iAB->offset();
      for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
        const int joff = jAB->offset();

// TODO remove this comment once the gammaforst issue has been fixed (bra and ket have been exchanged)
        std::shared_ptr<Matrix> block = couple_blocks(*jAB, *iAB);

        if (block) {
          hamiltonian_->add_block(1.0, joff, ioff, block->ndim(), block->mdim(), block);
          hamiltonian_->add_block(1.0, ioff, joff, block->mdim(), block->ndim(), block->transpose());
        }
      }
    }
    std::cout << "  o Computing off-diagonal blocks - time " << std::setw(9) << std::fixed << std::setprecision(2) << asdtime.tick() << std::endl;
  }

  std::cout << "  o Diagonalizing ASD Hamiltonian with a Davidson procedure" << std::endl;
  auto cc = std::make_shared<Matrix>(dimerstates_, nstates_);
  generate_initial_guess(cc, subspaces_base(), nstates_);
  std::cout << "    - initial guess time " << std::setw(9) << std::fixed << std::setprecision(2) << asdtime.tick() << std::endl << std::endl;

  energies_ = diagonalize(cc, subspaces_base());

  adiabats_ = cc->copy();

//adiabats_->print("ADIABATS", adiabats_->ndim());

  if (compute_rdm_) compute_rdm12();

  if (dipoles_) { // TODO Redo to make better use of memory
#if 0
    std::cout << "  o Computing properties" << std::endl;
    std::shared_ptr<const Reference> dimerref = dimer_->sref();
    DimerDipole dipole(dimerref, dimerref->nclosed(), dimerref->nclosed() + dimer_->active_refs().first->nact(), dimerref->nclosed() + dimerref->nact(), dimerref->coeff());
    std::array<std::string,3> mu_labels = {{"x", "y", "z"}};
    for (int i = 0; i < 3; ++i) {
      std::string label("mu_");
      label += mu_labels[i];
      std::shared_ptr<Matrix> tmp = compute_1e_prop(dipole.template dipoles<0>(i), dipole.template dipoles<1>(i), dipole.cross_dipole(i), dipole.core_dipole(i));
      auto prop = std::make_shared<Matrix>( (*adiabats_) % (*tmp) * (*adiabats_) );
      properties_.emplace_back(label, prop);
    }
#else
    throw std::logic_error("Dipole moments should be computed from density matrices");
#endif
  }

  print(print_thresh_);

  modelize();
}

#endif

#endif
