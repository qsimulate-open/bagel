//
// BAGEL - Parallel electron correlation program.
// Filename: asd_dmrg/product_guess.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
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

#include <src/asd_dmrg/product_rasci.h>

using namespace std;
using namespace bagel;

namespace bagel {
  namespace PRAS {
    struct Basis : public BlockKey {
      bitset<nbit__> alpha;
      bitset<nbit__> beta;
      int state;

      Basis(const BlockKey K, const int i, const bitset<nbit__> a, const bitset<nbit__> b) : BlockKey(K), alpha(a), beta(b), state(i) {}
      BlockKey key() const { return BlockKey(this->nelea, this->neleb); }
    };
  }
}

void ProductRASCI::model_guess(vector<shared_ptr<ProductRASCivec>>& out) {
  // First, order elements of denom_ by energy
  map<double, vector<PRAS::Basis>> ordered_elements;
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
  vector<PRAS::Basis> basis;

  for (auto& p : ordered_elements) {
    basis.insert(basis.end(), p.second.begin(), p.second.end());
    if (basis.size() >= nguess_)
      break;
  }
  const int nguess = basis.size();

  // In the future, I will want to do some sort of spin decontaminate
#if 0
  shared_ptr<Matrix> hamiltonian = make_shared<ProductCIHamiltonian>(basis, left_, jop_);
#else
  shared_ptr<Matrix> hamiltonian = make_shared<Matrix>(basis.size(), basis.size());
#endif
  VectorB eigs(nguess);
  hamiltonian->diagonalize(eigs);

  #if 1
    const double nuc_core = ref_->geom()->nuclear_repulsion() + jop_->core_energy();
    for (int i = 0; i < nguess; ++i)
      cout << setw(12) << setprecision(8) << eigs(i) + nuc_core << endl;
  #endif

  for (int ibasis = 0; ibasis < nguess; ++ibasis) {
    const PRAS::Basis b = basis[ibasis];
    for (int state = 0; state < nguess; ++state)
      out[state]->sector(b.key())->civec(b.state).element(b.beta, b.alpha) = hamiltonian->element(ibasis,state);
  }
}
