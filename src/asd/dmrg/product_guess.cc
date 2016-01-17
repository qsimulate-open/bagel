//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_dmrg/product_guess.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
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

#include <src/asd/dmrg/product_modelci.h>
#include <src/asd/dmrg/product_rasci.h>
#include <src/asd/dmrg/form_sigma.h>
#include <src/util/math/davidson.h>

using namespace std;
using namespace bagel;

void ProductRASCI::model_guess(vector<shared_ptr<ProductRASCivec>>& out) {
  // First, order elements of denom_ by energy
  map<double, vector<PCI::Basis>> ordered_elements;
  for (auto& sector : denom_->sectors()) {
    const int nstate = sector.second->nstates();
    for (int ist = 0; ist < nstate; ++ist) {
      const RASCivecView denomview = sector.second->civec(ist);
      for (auto& block : denomview.blocks()) {
        if (!block) continue;
        const double* d = block->data();
        for (auto& abit : *block->stringsa())
          for (auto& bbit : *block->stringsb())
            ordered_elements[*d++].emplace_back(sector.first, ist, abit, bbit);
      }
    }
  }

  // Now, pick the lowest lying elements to form a basis
  vector<PCI::Basis> basis;

  for (auto& p : ordered_elements) {
    basis.insert(basis.end(), p.second.begin(), p.second.end());
    if (basis.size() >= max(nguess_, nstate_))
      break;
  }
  const int guess_size = basis.size();

  // In the future, I may want to do some sort of spin decontaminate
  shared_ptr<Matrix> hamiltonian = make_shared<ProductCIHamiltonian>(basis, blockops_, jop_);
  VectorB eigs(guess_size);
  hamiltonian->diagonalize(eigs);

#if 0
    cout << "energies from the initial guess:" << endl;
    const double nuc_core = ref_->geom()->nuclear_repulsion() + jop_->core_energy();
    for (int i = 0; i < guess_size; ++i)
      cout << setw(12) << setprecision(8) << eigs(i) + nuc_core << endl;
#endif

  for (int ibasis = 0; ibasis < guess_size; ++ibasis) {
    const PCI::Basis b = basis[ibasis];
    for (int state = 0; state < nstate_; ++state)
      out[state]->sector(b.key())->civec(b.state).element(b.beta, b.alpha) = hamiltonian->element(ibasis,state);
  }

// run a few cycles using only "diagonal" terms
  if (preconverge_) {
    vector<shared_ptr<ProductRASCivec>> guessvecs = out;
    const double nuc_core = ref_->geom()->nuclear_repulsion() + jop_->core_energy();
    DavidsonDiag<ProductRASCivec> davidson(nstate_, davidson_subspace_);
    FormSigmaProdRAS form_sigma(batchsize_);
    cout << "  === Preconverging with particle number preserving terms ===" << endl;

    vector<bool> converged(nstate_, false);

    for (int iter = 0; iter != preconv_iter_; ++iter) {
      Timer calctime(2);
      vector<shared_ptr<ProductRASCivec>> sigma = form_sigma.diagonal(guessvecs, blockops_, jop_, converged);
      calctime.tick_print("diagonal sigma formation");

      // feeding sigma vectors into Davidson
      vector<shared_ptr<const ProductRASCivec>> ccn, sigman;
      for (int i = 0; i < nstate_; ++i) {
        if (!converged[i]) {
          ccn.push_back(make_shared<const ProductRASCivec>(*guessvecs.at(i))); // copy in civec
          sigman.push_back(sigma.at(i)); // place in sigma
        }
        else {
          ccn.push_back(nullptr);
          sigman.push_back(nullptr);
        }
      }
      const vector<double> energies = davidson.compute(ccn, sigman);
      // get residual and new vectors
      vector<shared_ptr<ProductRASCivec>> errvec = davidson.residual();
      calctime.tick_print("davidson");

      // compute errors
      vector<double> errors;
      for (int i = 0; i != nstate_; ++i) {
        errors.push_back(errvec.at(i)->rms());
        converged.at(i) = errors.at(i) < preconv_thresh_;
      }
      calctime.tick_print("error");

      if (!all_of(converged.begin(), converged.end(), [] (const bool b) { return b; })) {
        // denominator scaling
        for (int ist = 0; ist != nstate_; ++ist) {
          if (converged[ist]) continue;
          for (auto& denomsector : denom_->sectors()) {
            auto dbegin = denomsector.second->begin();
            auto dend = denomsector.second->end();
            auto targ_iter = guessvecs.at(ist)->sector(denomsector.first)->begin();
            auto source_iter = errvec.at(ist)->sector(denomsector.first)->begin();
            const double en = energies.at(ist);
            transform(dbegin, dend, source_iter, targ_iter, [&en] (const double den, const double cc) { return cc / min(en - den, -0.1); });
          }
          guessvecs.at(ist)->normalize();
        }
      }

      calctime.tick_print("denominator");

      // printing out
      if (nstate_ != 1 && iter > 0) cout << endl;
      for (int i = 0; i != nstate_; ++i) {
        cout << setw(7) << iter << setw(3) << i << setw(2) << (converged[i] ? "*" : " ")
                                << setw(17) << fixed << setprecision(8) << energies[i]+nuc_core << "   "
                                << setw(10) << scientific << setprecision(2) << errors[i] << fixed << setw(10) << setprecision(2)
                                << calctime.tick() << endl;
      }
      if (all_of(converged.begin(), converged.end(), [] (const bool b) { return b; })) break;
    }

    out = davidson.civec();

  }
}
