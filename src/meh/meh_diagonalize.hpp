//
// BAGEL - Parallel electron correlation program.
// Filename: meh_diagonalize.h
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

#ifndef BAGEL_MEH_DIAGONALIZE_H
#define BAGEL_MEH_DIAGONALIZE_H

template <class VecType>
std::shared_ptr<Matrix> MultiExcitonHamiltonian<VecType>::apply_hamiltonian(const Matrix& o, std::vector<DSubSpace>& subspaces) {
  const int nstates = o.mdim();

  std::shared_ptr<Matrix> out = o.clone();
  for (auto iAB = subspaces.begin(); iAB != subspaces.end(); ++iAB) {
    const int ioff = iAB->offset();
    for (auto jAB = subspaces.begin(); jAB != iAB; ++jAB) {
      const int joff = jAB->offset();

      if (store_matrix_) {
        dgemm_("N", "N", iAB->dimerstates(), nstates, jAB->dimerstates(), 1.0, hamiltonian_->element_ptr(ioff, joff), hamiltonian_->ndim(),
                                                                               o.element_ptr(joff, 0), o.ndim(),
                                                                          1.0, out->element_ptr(ioff, 0), out->ndim());
        dgemm_("T", "N", jAB->dimerstates(), nstates, iAB->dimerstates(), 1.0, hamiltonian_->element_ptr(ioff, joff), hamiltonian_->ndim(),
                                                                               o.element_ptr(ioff, 0), o.ndim(),
                                                                          1.0, out->element_ptr(joff, 0), out->ndim());
      }
      else {
        std::shared_ptr<const Matrix> block = couple_blocks<true>(*iAB, *jAB);

        if (block) {
          dgemm_("N", "N", block->ndim(), nstates, block->mdim(), 1.0, block->data(), block->ndim(), o.element_ptr(joff, 0), dimerstates_, 1.0, out->element_ptr(ioff, 0), o.ndim());
          dgemm_("T", "N", block->mdim(), nstates, block->ndim(), 1.0, block->data(), block->ndim(), o.element_ptr(ioff, 0), dimerstates_, 1.0, out->element_ptr(joff, 0), o.ndim());
        }
      }
    }

    if (store_matrix_) {
      dgemm_("N", "N", iAB->dimerstates(), nstates, iAB->dimerstates(), 1.0, hamiltonian_->element_ptr(ioff, ioff), hamiltonian_->ndim(),
                                                                             o.element_ptr(ioff, 0), o.ndim(),
                                                                        1.0, out->element_ptr(ioff, 0), out->ndim());
    }
    else {
      std::shared_ptr<const Matrix> block = compute_diagonal_block(*iAB);
      dgemm_("N", "N", block->ndim(), nstates, block->mdim(), 1.0, block->data(), block->ndim(), o.element_ptr(ioff, 0), dimerstates_, 1.0, out->element_ptr(ioff, 0), out->ndim());
    }
  }

  return out;
}

// Davidson diagonalize
//  - cc is initial guess on input, eigenvectors on exit
//  - subspaces is the subspaces over which to apply the Hamiltonian
//  - mute is whether to print convergence info
template <class VecType>
std::vector<double> MultiExcitonHamiltonian<VecType>::diagonalize(std::shared_ptr<Matrix>& cc, std::vector<DSubSpace>& subspaces, const bool mute) {
  Timer mehtime;
  const int nstates = cc->mdim();

  DavidsonDiag<Matrix> davidson(nstates, davidson_subspace_);

  std::vector<bool> conv(nstates, false);
  std::vector<double> out(nstates, 0.0);

  for (int iter = 0; iter != max_iter_; ++iter) {
    std::shared_ptr<const Matrix> sigma = apply_hamiltonian(*cc, subspaces);

    std::vector<std::shared_ptr<const Matrix>> sigman;
    std::vector<std::shared_ptr<const Matrix>> ccn;
    for (int i = 0; i != nstates; ++i) {
      if (!conv[i]) {
        sigman.push_back(sigma->slice_copy(i,i+1));
        ccn.push_back(cc->slice_copy(i,i+1));
      }
      else {
        sigman.push_back(std::shared_ptr<const Matrix>());
        ccn.push_back(std::shared_ptr<const Matrix>());
      }
    }
    const std::vector<double> energies = davidson.compute(ccn, sigman);

    // residual
    std::vector<std::shared_ptr<Matrix>> errvec = davidson.residual();
    std::vector<double> errors;
    for (int i = 0; i != nstates; ++i) {
      errors.push_back(errvec.at(i)->rms());
      conv.at(i) = errors.at(i) < thresh_;
    }

    if (std::any_of(conv.begin(), conv.end(), [] (const bool t) { return (!t); })) {
      for (int ist = 0; ist != nstates; ++ist) {
        if (conv.at(ist)) continue;
        auto tmp_cc = std::make_shared<Matrix>(dimerstates_, 1);
        double* target_array = tmp_cc->data();
        double* source_array = errvec.at(ist)->data();
        const double en = energies.at(ist);
        for (auto& space : subspaces) {
          for (int i = space.offset(); i != space.offset()+space.dimerstates(); ++i) {
            target_array[i] = source_array[i] / std::min(en - denom_[i], -0.1);
          }
        }
        std::list<std::shared_ptr<const Matrix>> tmp;
        spin_->filter(*tmp_cc, nspin_);
        double nrm = tmp_cc->norm();
        double scal = (nrm > 1.0e-15 ? 1.0/nrm : 0.0);
        tmp_cc->scale(scal);
        std::copy_n(tmp_cc->data(), dimerstates_, cc->element_ptr(0, ist));
      }
      cc->broadcast();
    }

    if (!mute) {
      if (nstates != 1 && iter) std::cout << std::endl;
      for (int i = 0; i != nstates; ++i) {
        std::cout << std::setw(7) << iter << std::setw(3) << i << std::setw(2) << (conv[i] ? "*" : " ")
                                << std::setw(17) << std::fixed << std::setprecision(8) << energies[i] << "   "
                                << std::setw(10) << std::scientific << std::setprecision(2) << errors[i] << "   "
                                << std::fixed << std::setw(10) << std::setprecision(2) << mehtime.tick() << std::endl;
      }
    }
    std::copy(energies.begin(), energies.end(), out.begin());
    if (std::all_of(conv.begin(), conv.end(), [] (const bool t) { return (t); })) break;
  }

  std::vector<std::shared_ptr<Matrix>> eigenstates = davidson.civec();
  for (int i = 0; i != nstates; ++i)
    std::copy_n(eigenstates.at(i)->data(), dimerstates_, cc->element_ptr(0, i));

  return out;
}

#endif

#endif
