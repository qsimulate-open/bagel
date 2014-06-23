//
// BAGEL - Parallel electron correlation program.
// Filename: meh_compute.hpp
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

#ifndef BAGEL_MEH_COMPUTE_H
#define BAGEL_MEH_COMPUTE_H

template <class VecType>
void MultiExcitonHamiltonian<VecType>::generate_initial_guess(std::shared_ptr<Matrix> cc, std::vector<DSubSpace>& subspaces, int nstates) {
  int trialsize = 0;
  int nguess = nguess_;

  const int subspace_states = std::accumulate(subspaces.begin(), subspaces.end(), 0,
                                              [] (int x, const DSubSpace& s) { return s.dimerstates()+x; });

  std::map<double, std::vector<int>> seeds;
  for(auto& ispace : subspaces) {
    for(int state = ispace.offset(); state < ispace.offset()+ispace.dimerstates(); ++state)
      seeds[denom_[state]].push_back(state);
  }

  while (trialsize < nstates) {
    std::vector<int> b;
    b.reserve(nguess);

    for (auto& i : seeds) {
      b.insert(b.end(), i.second.begin(), i.second.end());
      if (b.size() >= nguess) break;
    }

    // build matrix
    auto basis = std::make_shared<Matrix>(dimerstates_, b.size());
    for (int i = 0; i < b.size(); ++i)
      basis->element(b[i], i) = 1.0;

    // build spin operator
    std::shared_ptr<Matrix> spn = spin_->apply(*basis);
    spn = std::make_shared<Matrix>( *spn % *basis );
    std::vector<double> spin_values(b.size(), 0.0);
    spn->diagonalize(spin_values.data());
    const double expected_spin = 0.25 * static_cast<double>(nspin_ * (nspin_ + 2));
    int start, end;
    for (start = 0; start < nguess; ++start)
      if (std::fabs(spin_values[start] - expected_spin) < 1.0e-4) break;
    for (end = start; end < nguess; ++end)
      if (std::fabs(spin_values[end] - expected_spin) > 1.0e-4) break;

    trialsize = end - start;

    if (trialsize >= nstates) {
      basis = (*basis * *spn).slice_copy(start, end);

      std::shared_ptr<const Matrix> sigma = apply_hamiltonian(*basis, subspaces_);
      auto H = std::make_shared<Matrix>(*sigma % *basis);
      std::vector<double> energies(trialsize, 0.0);
      H->diagonalize(energies.data());

      basis = std::make_shared<Matrix>(*basis * *H);
      for (int i = 0; i < nstates; ++i)
        std::copy_n(basis->element_ptr(0, i), basis->ndim(), cc->element_ptr(0, i));
    }
    else if (nguess >= subspace_states) {
      throw std::runtime_error("Requesting more spin allowed states than exist in MEH space");
    }
    else {
      nguess *= 2;
    }
  }
}

template <class VecType>
void MultiExcitonHamiltonian<VecType>::compute() {
  Timer mehtime;
  std::cout << std::endl << " ===== Starting construction of dimer Hamiltonian " << std::endl;
  std::cout << "   o Dimer space:" << std::endl;
  std::cout << "     -  spin: " << nspin_ << std::endl;
  std::cout << "     -  charge: " << charge_ << std::endl;
  std::cout << "     -  dimer states: " << dimerstates_ << std::endl << std::endl;

  {
    std::map<std::pair<int,int>, double> spinmap;
    for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
      for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
        gamma_couple_blocks(*iAB, *jAB);
        spin_couple_blocks(*iAB, *jAB, spinmap);
      }
      gamma_couple_blocks(*iAB, *iAB);
      compute_diagonal_spin_block(*iAB, spinmap);
    }
    spin_ = std::make_shared<MEHSpin>(dimerstates_, spinmap, max_spin_);
  }

  std::cout << "  o Preparing Gamma trees and building spin operator - " << std::setw(9) << std::fixed << std::setprecision(2) << mehtime.tick() << std::endl;
  std::cout << "    - spin elements: " << spin_->size() << std::endl;

  gammaforest_->compute();

  std::cout << "  o Computing Gamma trees - " << std::setw(9) << std::fixed << std::setprecision(2) << mehtime.tick() << std::endl;

  if (store_matrix_) hamiltonian_ = std::make_shared<Matrix>(dimerstates_, dimerstates_);

  denom_ = std::unique_ptr<double[]>(new double[dimerstates_]);

  for (auto& subspace : subspaces_) {
    compute_pure_terms(subspace, jop_);
    std::shared_ptr<Matrix> block = compute_diagonal_block(subspace);
    if (store_matrix_)
      hamiltonian_->add_block(1.0, subspace.offset(), subspace.offset(), block->ndim(), block->mdim(), block);
    const int n = block->ndim();
    for ( int i = 0; i < n; ++i ) denom_[subspace.offset() + i] = block->element(i,i);
  }
  std::cout << "  o Computing diagonal blocks and building denominator - time " << std::setw(9) << std::fixed << std::setprecision(2) << mehtime.tick() << std::endl;

  if (store_matrix_) {
    for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
      const int ioff = iAB->offset();
      for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
        const int joff = jAB->offset();

        std::shared_ptr<Matrix> block = couple_blocks(*iAB, *jAB);

        if (block) {
          hamiltonian_->add_block(1.0, ioff, joff, block->ndim(), block->mdim(), block);
          hamiltonian_->add_block(1.0, joff, ioff, block->mdim(), block->ndim(), block->transpose());
        }
      }
    }
    std::cout << "  o Computing off-diagonal blocks - time " << std::setw(9) << std::fixed << std::setprecision(2) << mehtime.tick() << std::endl;
  }

  std::cout << "  o Diagonalizing ME Hamiltonian with a Davidson procedure" << std::endl;
  auto cc = std::make_shared<Matrix>(dimerstates_, nstates_);
  generate_initial_guess(cc, subspaces_, nstates_);
  std::cout << "    - initial guess time " << std::setw(9) << std::fixed << std::setprecision(2) << mehtime.tick() << std::endl << std::endl;

  energies_ = diagonalize(cc, subspaces_);

  adiabats_ = cc->copy();

  if ( dipoles_ ) { // TODO Redo to make better use of memory
    std::cout << "  o Computing properties" << std::endl;
    std::shared_ptr<const Reference> dimerref = dimer_->sref();
    DimerDipole dipole = DimerDipole(dimerref, dimerref->nclosed(), dimerref->nclosed() + dimer_->nact().first, dimerref->nclosed() + dimerref->nact(), dimerref->coeff());
    std::array<std::string,3> mu_labels = {{"x", "y", "z"}};
    for (int i = 0; i < 3; ++i) {
      std::string label("mu_");
      label += mu_labels[i];
      std::shared_ptr<Matrix> tmp = compute_1e_prop(dipole.template dipoles<0>(i), dipole.template dipoles<1>(i), dipole.cross_dipole(i), dipole.core_dipole(i));
      auto prop = std::make_shared<Matrix>( (*adiabats_) % (*tmp) * (*adiabats_) );
      properties_.emplace_back(label, prop);
    }
  }

  print(print_thresh_);

  modelize();
}

#endif

#endif
