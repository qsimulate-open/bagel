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
    std::vector<double> diagonals;
    std::cout << "Building model Hamiltonian from model states:" << std::endl;

    std::vector<int> space_bounds;

    for(auto& block : model) {
      int S_A, q_A, S_B, q_B;
      std::tie(S_A, S_B) = block.S_;
      std::tie(q_A, q_B) = block.charge_;
      const int nstates = block.nstates_;
      std::cout << "  o S: (" << S_A << ", " << S_B << "), Q: (" << q_A << ", " << q_B << "), nstates: " << nstates << std::endl;

      std::vector<DSubSpace> blocks_in_subspace;
      for (auto& sp : subspaces_)
        if ( block.S_ == sp.S() && block.charge_ == sp.charge() )
          blocks_in_subspace.push_back(sp);
      std::for_each(blocks_in_subspace.begin(), blocks_in_subspace.end(), [&space_bounds] (const DSubSpace& s) { space_bounds.push_back(s.offset()); });

      // diagonalize in subspace
      auto cc = std::make_shared<Matrix>(dimerstates_, nstates);
      generate_initial_guess(cc, blocks_in_subspace, nstates);
      std::vector<double> eigenvalues = diagonalize(cc, blocks_in_subspace, /*mute=*/true);

      // append to model space
      modelstates.push_back(cc);
      diagonals.insert(diagonals.end(), eigenvalues.begin(), eigenvalues.end());
    }

    std::cout << std::endl;

    const int modelsize = std::accumulate(modelstates.begin(), modelstates.end(), 0, [] (int x, std::shared_ptr<const Matrix> m) { return x + m->mdim(); });
    auto modelcc = std::make_shared<Matrix>(dimerstates_, modelsize);
    int offset = 0;
    for(auto& m : modelstates) {
      modelcc->copy_block(0, offset, dimerstates_, m->mdim(), *m);
      offset += m->mdim();
    }
    print_states(*modelcc, diagonals, print_thresh_, "Model states");
    std::shared_ptr<Matrix> modelsigma = apply_hamiltonian(*modelcc, subspaces_);
    auto model_hamiltonian = std::make_shared<Matrix>(*modelsigma % *modelcc);

    const double E0 = model_hamiltonian->element(0,0);
    std::cout << "E_0 = " << E0 << " H" << std::endl;
    model_hamiltonian->add_diag(-E0);
    model_hamiltonian->scale(au2eV__);
    model_hamiltonian->print("Model Hamiltonian (eV)", modelsize);

    // Compute perturbative correction
    std::shared_ptr<Matrix> perturbation = model_hamiltonian->clone();
    for (int i = 0; i < modelsize; ++i) {
      const double Ei = diagonals[i];
      for (int j = 0; j <= i; ++j) {
        const double Ej = diagonals[j];
        double perturb = 0;
        for (auto& space : subspaces_) {
          if (std::count(space_bounds.begin(), space_bounds.end(), space.offset()) == 1) continue;
          for (int k = 0; k < space.dimerstates(); ++k) {
            const int kk = k + space.offset();
            const double fac = 1.0/(denom_[kk] - Ei) + 1.0/(denom_[kk] - Ej);
            const double V = modelsigma->element(kk,i) * modelsigma->element(kk,j);
            perturb -= 0.5*V*fac;
          }
        }
        perturbation->element(i,j) = perturbation->element(j,i) = perturb;
      }
    }
    perturbation->scale(au2eV__);
    std::cout << std::endl;
    perturbation->print("Perturbative correction (eV)", modelsize);

    models_.emplace_back(model_hamiltonian, perturbation);
  }
}

#endif

#endif
