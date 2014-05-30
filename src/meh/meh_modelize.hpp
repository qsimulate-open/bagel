//
// BAGEL - Parallel electron correlation program.
// Filename: meh_modelize.hpp
// Copyright (C) 2014 Toru Shiozaki
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

#ifndef BAGEL_MEH_MODELIZE_H
#define BAGEL_MEH_MODELIZE_H

template <class VecType>
void MultiExcitonHamiltonian<VecType>::modelize() {
  // makes, prints, and stores multiple model Hamiltonians
  for(auto& model : models_to_form_) {
    std::vector<std::shared_ptr<Matrix>> modelstates;
    std::cout << "Asking to build a model Hamiltonian from:" << std::endl;
    for(auto& block : model) {
      int S_A, q_A, S_B, q_B;
      std::tie(S_A, S_B) = block.S_;
      std::tie(q_A, q_B) = block.charge_;
      const int nstates = block.nstates_;
      std::cout << "  -";
      std::cout << " (" << S_A << ", " << q_A << ", " << nstates << ")";
      std::cout << " (" << S_B << ", " << q_B << ", " << nstates << ")";
      std::cout << std::endl;

      std::vector<DSubSpace> blocks_in_subspace;
      for (auto& sp : subspaces_)
        if ( block.S_ == sp.S() && block.charge_ == sp.charge() )
          blocks_in_subspace.push_back(sp);

      // Diagonalize in subspace

      // initial guess
      auto cc = std::make_shared<Matrix>(dimerstates_, nstates);
      generate_initial_guess(cc, blocks_in_subspace, nstates);

      // Davidson diagonalize

      // append to model space
      modelstates.push_back(cc);
    }

    const int modelsize = std::accumulate(modelstates.begin(), modelstates.end(), 0, [] (int x, std::shared_ptr<const Matrix> m) { return x + m->mdim(); });
    auto modelcc = std::make_shared<Matrix>(dimerstates_, modelsize);
    int offset = 0;
    for(auto& m : modelstates) {
      modelcc->copy_block(0, offset, dimerstates_, m->mdim(), *m);
      offset += m->mdim();
    }
    std::shared_ptr<Matrix> modelsigma = apply_hamiltonian(*modelcc);
    auto model_hamiltonian = std::make_shared<Matrix>(*modelsigma % *modelcc);
    model_hamiltonian->print("Model Hamiltonian", modelsize);
  }
}

#endif

#endif
