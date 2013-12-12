//
// BAGEL - Parallel electron correlation program.
// Filename: meh_compute.h
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

#ifndef BAGEL_MEH_H
#define BAGEL_MEH_H

template <class VecType>
void MultiExcitonHamiltonian<VecType>::compute() {
  Timer mehtime;
  std::cout << std::endl << " ===== Starting construction of dimer Hamiltonian with " << dimerstates_ << " states ===== " << std::endl;

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

  std::cout << "  o Preparing Gamma trees and building spin operator - " << std::setw(9) << std::fixed << std::setprecision(2) << mehtime.tick() << std::endl;
  std::cout << "    - spin elements: " << spin_->size() << std::endl;

  gammaforest_->compute();

  std::cout << "  o Computing Gamma trees - " << std::setw(9) << std::fixed << std::setprecision(2) << mehtime.tick() << std::endl;

  if (store_matrix_) hamiltonian_ = std::make_shared<Matrix>(dimerstates_, dimerstates_);

  denom_ = std::unique_ptr<double[]>(new double[dimerstates_]);

  for (auto& subspace : subspaces_) {
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

        hamiltonian_->add_block(1.0, ioff, joff, block->ndim(), block->mdim(), block);
        hamiltonian_->add_block(1.0, joff, ioff, block->mdim(), block->ndim(), block->transpose());
      }
    }
    std::cout << "  o Computing off-diagonal blocks - time " << std::setw(9) << std::fixed << std::setprecision(2) << mehtime.tick() << std::endl;
  }

  std::cout << "  o Diagonalizing ME Hamiltonian with a Davidson procedure" << std::endl;
  DavidsonDiag<Matrix> davidson(nstates_, davidsonceiling_);

  // Generate initial guesses
  auto cc = std::make_shared<Matrix>(dimerstates_, nstates_);
  int trialsize = std::min(2 * spin_->max() * nstates_, dimerstates_);
  {
    int multiple = 1;
    std::multimap<double,int> seeds;
    while (seeds.size() < trialsize) {
      seeds.clear();
      for ( auto& iAB : subspaces_ ) {
        const int ioff = iAB.offset();
        const int nstatesA = iAB.template nstates<0>();
        const int nstatesB = iAB.template nstates<1>();
        for (int mult = 0; mult != multiple; ++mult) {
          for (int i = 0; i != nstatesA; ++i) {
            const int dindex = iAB.dimerindex(i,mult) + ioff;
            seeds.insert(std::make_pair(denom_[dindex], dindex));
          }
          for (int i = mult + 1; i != nstatesB; ++i) {
            const int dindex = iAB.dimerindex(mult,i) + ioff;
            seeds.insert(std::make_pair(denom_[dindex], dindex));
          }
        }
      }
      ++multiple;
    }

    auto iseed = seeds.begin();
    auto initial = std::make_shared<Matrix>(dimerstates_,trialsize);
    for (int istate = 0; istate != trialsize; ++istate, ++iseed)
      initial->element(iseed->second, istate) = 1.0;
    spin_->filter(*initial, nspin_);

    // Symmetric orthogonalization to get rid of linear dependencies
    Matrix overlap = *initial % *initial;
    overlap.inverse_half();
    *initial = *initial * overlap;

    // Average diagonal elements to pick the best options
    std::multimap<double, int> tmpmap;
    for (int j = 0; j < trialsize; ++j) {
      double energy = 0;
      for (int i = 0; i < dimerstates_; ++i)
        energy += initial->element(i,j) * initial->element(i,j) * denom_[i];
      tmpmap.emplace(energy,j);
    }

    auto imap = tmpmap.begin();
    for (int istate = 0; istate < nstates_; ++istate, ++imap) {
      int ii = imap->second;
      std::copy_n(initial->element_ptr(0, ii), dimerstates_, cc->element_ptr(0, istate));
    }
  }
  std::cout << "    - initial guess time " << std::setw(9) << std::fixed << std::setprecision(2) << mehtime.tick() << std::endl << std::endl;

  std::vector<int> conv(nstates_, static_cast<int>(false));

  for (int iter = 0; iter != max_iter_; ++iter) {
    std::shared_ptr<const Matrix> sigma = apply_hamiltonian(*cc);

    std::vector<std::shared_ptr<const Matrix>> sigman;
    std::vector<std::shared_ptr<const Matrix>> ccn;
    for (int i = 0; i != nstates_; ++i) {
      if (!conv[i]) {
        sigman.push_back(sigma->slice(i,i+1));
        ccn.push_back(cc->slice(i,i+1));
      }
    }
    const std::vector<double> energies = davidson.compute(ccn, sigman);

    // residual
    std::vector<std::shared_ptr<Matrix>> errvec = davidson.residual();
    std::vector<double> errors;
    for (int i = 0; i != nstates_; ++i) {
      errors.push_back(errvec.at(i)->variance());
      conv.at(i) = static_cast<int>(errors.at(i) < thresh_);
    }

    if (!*min_element(conv.begin(), conv.end())) {
      for (int ist = 0; ist != nstates_; ++ist) {
        if (conv.at(ist)) continue;
        const int size = dimerstates_;
        auto tmp_cc = std::make_shared<Matrix>(dimerstates_, 1);
        double* target_array = tmp_cc->data();
        double* source_array = errvec.at(ist)->data();
        const double en = energies.at(ist);
        for (int i = 0; i != size; ++i) {
          target_array[i] = source_array[i] / std::min(en - denom_[i], -0.1);
        }
        davidson.orthog(tmp_cc);
        std::list<std::shared_ptr<const Matrix>> tmp;
        for (int jst = 0; jst != ist; ++jst) tmp.push_back(cc->slice(jst, jst+1));
        tmp_cc->orthog(tmp);
        spin_->filter(*tmp_cc, nspin_);
        double nrm = tmp_cc->norm();
        double scal = (nrm > 1.0e-15 ? 1.0/nrm : 0.0);
        tmp_cc->scale(scal);
        std::copy_n(tmp_cc->data(), dimerstates_, cc->element_ptr(0, ist));
      }
    }

    if (nstates_ != 1 && iter) std::cout << std::endl;
    for (int i = 0; i != nstates_; ++i) {
      std::cout << std::setw(7) << iter << std::setw(3) << i << std::setw(2) << (conv[i] ? "*" : " ")
                              << std::setw(17) << std::fixed << std::setprecision(8) << energies[i] << "   "
                              << std::setw(10) << std::scientific << std::setprecision(2) << errors[i] << "   "
                              << std::fixed << std::setw(10) << std::setprecision(2) << mehtime.tick() << std::endl;
      energies_.at(i) = energies[i];
    }
    if (*min_element(conv.begin(), conv.end())) break;
  }

  adiabats_ = std::make_shared<Matrix>(dimerstates_, nstates_);
  std::vector<std::shared_ptr<Matrix>> advec = davidson.civec();
  for (int i = 0; i != nstates_; ++i) {
    std::copy_n(advec.at(i)->data(), dimerstates_, adiabats_->element_ptr(0, i));
  }

  if ( dipoles_ ) { // TODO Redo to make better use of memory
    std::cout << "  o Computing properties" << std::endl;
    DimerDipole dipole = DimerDipole(ref_, dimerclosed_, dimerclosed_ + nact_.first, dimerclosed_ + dimeractive_, coeff_);
    std::array<std::string,3> mu_labels = {{"x", "y", "z"}};
    for (int i = 0; i < 3; ++i) {
      std::string label("mu_");
      label += mu_labels[i];
      std::shared_ptr<Matrix> tmp = compute_1e_prop(dipole.template dipoles<0>(i), dipole.template dipoles<1>(i), dipole.cross_dipole(i), dipole.core_dipole(i));
      auto prop = std::make_shared<Matrix>( (*adiabats_) % (*tmp) * (*adiabats_) );
      properties_.emplace_back(label, prop);
    }
  }

  print(nstates_, print_thresh_);
}

#endif

#endif
