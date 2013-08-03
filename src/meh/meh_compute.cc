//
// BAGEL - Parallel electron correlation program.
// Filename: meh.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/dimer/dimer_prop.h>
#include <src/math/davidson.h>
#include <src/meh/meh.h>

using namespace std;
using namespace bagel;

void MultiExcitonHamiltonian::compute() {
  Timer mehtime;
  cout << endl << " ===== Starting construction of dimer Hamiltonian with " << dimerstates_ << " states ===== " << endl;

  for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
    for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
      gamma_couple_blocks(*iAB, *jAB);
      spin_couple_blocks(*iAB, *jAB);
    }
    gamma_couple_blocks(*iAB, *iAB);
    compute_diagonal_spin_block(*iAB);
  }

  cout << "  o Preparing Gamma trees and building spin operator - " << setw(9) << fixed << setprecision(2) << mehtime.tick() << endl;
  cout << "    - offdiagonal spin elements: " << spin_->offdiagonal().size() << endl;

  gammaforest_->compute();

  cout << "  o Computing Gamma trees - " << setw(9) << fixed << setprecision(2) << mehtime.tick() << endl;

  if (store_matrix_) hamiltonian_ = make_shared<Matrix>(dimerstates_, dimerstates_);

  denom_ = unique_ptr<double[]>(new double[dimerstates_]);

  for (auto& subspace : subspaces_) {
    shared_ptr<Matrix> block = compute_diagonal_block(subspace);
    if (store_matrix_)
      hamiltonian_->add_block(subspace.offset(), subspace.offset(), block);
    const int n = block->ndim();
    for ( int i = 0; i < n; ++i ) denom_[subspace.offset() + i] = block->element(i,i);
  }
  cout << "  o Computing diagonal blocks and building denominator - time " << setw(9) << fixed << setprecision(2) << mehtime.tick() << endl;

  if (store_matrix_) {
    for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
      const int ioff = iAB->offset();
      for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
        const int joff = jAB->offset();

        shared_ptr<Matrix> block = couple_blocks(*iAB, *jAB);

        hamiltonian_->add_block(ioff, joff, block);
        hamiltonian_->add_block(joff, ioff, block->transpose());
      }
    }
    cout << "  o Computing off-diagonal blocks - time " << setw(9) << fixed << setprecision(2) << mehtime.tick() << endl;
  }

  cout << "  o Diagonalizing ME Hamiltonian with a Davidson procedure" << endl;
  DavidsonDiag<Matrix> davidson(nstates_, max_iter_);

  // Generate initial guesses
  auto cc = make_shared<Matrix>(dimerstates_, nstates_);
  int trialsize = min(2 * spin_->max() * nstates_, dimerstates_);
  {
    int multiple = 1;
    multimap<double,int> seeds;
    while (seeds.size() < trialsize) {
      seeds.clear();
      for ( auto& iAB : subspaces_ ) {
        const int ioff = iAB.offset();
        const int nstatesA = iAB.nstates<0>();
        const int nstatesB = iAB.nstates<1>();
        for (int mult = 0; mult != multiple; ++mult) {
          for (int i = 0; i != nstatesA; ++i) {
            const int dindex = iAB.dimerindex(i,mult) + ioff;
            seeds.insert(make_pair(denom_[dindex], dindex));
          }
          for (int i = mult + 1; i != nstatesB; ++i) {
            const int dindex = iAB.dimerindex(mult,i) + ioff;
            seeds.insert(make_pair(denom_[dindex], dindex));
          }
        }
      }
      ++multiple;
    }

    auto iseed = seeds.begin();
    auto initial = make_shared<Matrix>(dimerstates_,trialsize);
    for (int istate = 0; istate != trialsize; ++istate, ++iseed)
      initial->element(iseed->second, istate) = 1.0;
    spin_->filter(*initial, nspin_);

    // Symmetric orthogonalization to get rid of linear dependencies
    Matrix overlap = *initial % *initial;
    overlap.inverse_half();
    *initial = *initial * overlap;

    // Average diagonal elements to pick the best options
    multimap<double, int> tmpmap;
    for (int j = 0; j < trialsize; ++j) {
      double energy = 0;
      for (int i = 0; i < dimerstates_; ++i)
        energy += initial->element(i,j) * initial->element(i,j) * denom_[i];
      tmpmap.insert(make_pair(energy,j));
    }

    auto imap = tmpmap.begin();
    for (int istate = 0; istate < nstates_; ++istate, ++imap) {
      int ii = imap->second;
      copy_n(initial->element_ptr(0, ii), dimerstates_, cc->element_ptr(0, istate));
    }
  }
  cout << "    - initial guess time " << setw(9) << fixed << setprecision(2) << mehtime.tick() << endl << endl;

  vector<int> conv(nstates_, static_cast<int>(false));

  for (int iter = 0; iter != max_iter_; ++iter) {
    shared_ptr<const Matrix> sigma = apply_hamiltonian(*cc);

    vector<shared_ptr<const Matrix>> sigman;
    vector<shared_ptr<const Matrix>> ccn;
    for (int i = 0; i != nstates_; ++i) {
      if (!conv[i]) {
        sigman.push_back(sigma->slice(i,i+1));
        ccn.push_back(cc->slice(i,i+1));
      }
    }
    const vector<double> energies = davidson.compute(ccn, sigman);

    // residual
    vector<shared_ptr<Matrix>> errvec = davidson.residual();
    vector<double> errors;
    for (int i = 0; i != nstates_; ++i) {
      errors.push_back(errvec.at(i)->variance());
      conv.at(i) = static_cast<int>(errors.at(i) < thresh_);
    }

    if (!*min_element(conv.begin(), conv.end())) {
      for (int ist = 0; ist != nstates_; ++ist) {
        if (conv.at(ist)) continue;
        const int size = dimerstates_;
        auto tmp_cc = make_shared<Matrix>(dimerstates_, 1);
        double* target_array = tmp_cc->data();
        double* source_array = errvec.at(ist)->data();
        const double en = energies.at(ist);
        for (int i = 0; i != size; ++i) {
          target_array[i] = source_array[i] / min(en - denom_[i], -0.1);
        }
        davidson.orthog(tmp_cc);
        list<shared_ptr<const Matrix>> tmp;
        for (int jst = 0; jst != ist; ++jst) tmp.push_back(cc->slice(jst, jst+1));
        tmp_cc->orthog(tmp);
        spin_->filter(*tmp_cc, nspin_);
        double nrm = tmp_cc->norm();
        double scal = (nrm > 1.0e-15 ? 1.0/nrm : 0.0);
        tmp_cc->scale(scal);
        copy_n(tmp_cc->data(), dimerstates_, cc->element_ptr(0, ist));
      }
    }

    if (nstates_ != 1 && iter) cout << endl;
    for (int i = 0; i != nstates_; ++i) {
      cout << setw(7) << iter << setw(3) << i << setw(2) << (conv[i] ? "*" : " ")
                              << setw(17) << fixed << setprecision(8) << energies[i] << "   "
                              << setw(10) << scientific << setprecision(2) << errors[i] << "   "
                              << fixed << setw(10) << setprecision(2) << mehtime.tick() << endl;
      energies_.at(i) = energies[i];
    }
    if (*min_element(conv.begin(), conv.end())) break;
  }

  adiabats_ = make_shared<Matrix>(dimerstates_, nstates_);
  vector<shared_ptr<Matrix>> advec = davidson.civec();
  for (int i = 0; i != nstates_; ++i) {
    copy_n(advec.at(i)->data(), dimerstates_, adiabats_->element_ptr(0, i));
  }

  if ( dipoles_ ) { // TODO Redo to make better use of memory
    cout << "  o Computing properties" << endl;
    DimerDipole dipole = DimerDipole(ref_, dimerclosed_, dimerclosed_ + nact_.first, dimerclosed_ + dimeractive_, coeff_);
    array<string,3> mu_labels = {{"x", "y", "z"}};
    for (int i = 0; i < 3; ++i) {
      string label("mu_");
      label += mu_labels[i];
      shared_ptr<Matrix> tmp = compute_1e_prop(dipole.dipoles<0>(i), dipole.dipoles<1>(i), dipole.cross_dipole(i), dipole.core_dipole(i));
      shared_ptr<Matrix> prop(new Matrix( (*adiabats_) % (*tmp) * (*adiabats_) ));
      properties_.push_back(make_pair(label, prop));
    }
  }

  print(nstates_, print_thresh_);
}
