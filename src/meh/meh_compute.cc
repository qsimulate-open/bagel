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
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

  hamiltonian_ = make_shared<Matrix>(dimerstates_, dimerstates_);

  cout << "  o Computing diagonal blocks" << endl;
  for (auto& subspace : subspaces_) {
    hamiltonian_->add_block(subspace.offset(), subspace.offset(), compute_diagonal_block(subspace));
  }
  cout << "    - time " << setw(9) << fixed << setprecision(2) << mehtime.tick() << endl;

  cout << "  o Computing off-diagonal blocks" << endl;
  for (auto iAB = subspaces_.begin(); iAB != subspaces_.end(); ++iAB) {
    const int ioff = iAB->offset();
    for (auto jAB = subspaces_.begin(); jAB != iAB; ++jAB) {
      const int joff = jAB->offset();

      shared_ptr<Matrix> block = couple_blocks(*iAB, *jAB);

      hamiltonian_->add_block(ioff, joff, block);
      hamiltonian_->add_block(joff, ioff, block->transpose());
    }
  }
  cout << "    - time " << setw(9) << fixed << setprecision(2) << mehtime.tick() << endl;

  cout << "  o Diagonalizing ME Hamiltonian with a Davidson procedure" << endl;
  DavidsonDiag<Matrix> davidson(nstates_, max_iter_);

  // Generate initial guesses
  vector<shared_ptr<Matrix>> cc;
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
            seeds.insert(make_pair(hamiltonian_->element(dindex,dindex), dindex));
          }
          for (int i = mult + 1; i != nstatesB; ++i) {
            const int dindex = iAB.dimerindex(mult,i) + ioff;
            seeds.insert(make_pair(hamiltonian_->element(dindex,dindex), dindex));
          }
        }
      }
      ++multiple;
    }

    auto iseed = seeds.begin();
    auto initial = make_shared<Matrix>(dimerstates_,trialsize);
    for (int istate = 0; istate != trialsize; ++istate, ++iseed) {
      initial->element(iseed->second, istate) = 1.0;
    }
    spin_->filter(*initial, nspin_);

    Matrix overlap = *initial % *initial;
    overlap.inverse_half();
    *initial = *initial * overlap;
    Matrix initialham = *initial % *hamiltonian_ * *initial;
    multimap<double, int> tmpmap;
    for (int i = 0; i < trialsize; ++i) {
      tmpmap.insert(make_pair(initialham(i,i),i));
    }

    auto imap = tmpmap.begin();
    for (int istate = 0; istate < nstates_; ++istate, ++imap) {
      int ii = imap->second;
      cc.push_back(initial->slice(ii,ii+1));
    }
  }

  cout << "    - initial guess time " << setw(9) << fixed << setprecision(2) << mehtime.tick() << endl << endl;

  vector<int> conv(nstates_, static_cast<int>(false));

  for (int iter = 0; iter != max_iter_; ++iter) {
    vector<shared_ptr<const Matrix>> sigma;
    vector<shared_ptr<const Matrix>> ccn;
    for (int i = 0; i != nstates_; ++i) {
      if (!conv[i]) {
        sigma.push_back(make_shared<const Matrix>(*hamiltonian_ * *cc.at(i)));
        ccn.push_back(make_shared<const Matrix>(*cc.at(i)));
      }
    }
    const vector<double> energies = davidson.compute(ccn, sigma);

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
        double* target_array = cc.at(ist)->data();
        double* source_array = errvec.at(ist)->data();
        const double en = energies.at(ist);
        for (int i = 0; i != size; ++i) {
          target_array[i] = source_array[i] / min(en - hamiltonian_->element(i,i), -0.1);
        }
        davidson.orthog(cc.at(ist));
        list<shared_ptr<const Matrix>> tmp;
        for (int jst = 0; jst != ist; ++jst) tmp.push_back(cc.at(jst));
        cc.at(ist)->orthog(tmp);
        spin_->filter(*cc.at(ist), nspin_);
        double nrm = cc.at(ist)->norm();
        double scal = (nrm > 1.0e-15 ? 1.0/nrm : 0.0);
        cc.at(ist)->scale(scal);
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

  if ( dipoles_ ) {
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
