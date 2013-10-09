//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_forest.cc
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#include <src/meh/gamma_forest.h>

using namespace std;
using namespace bagel;

template <>
void GammaForest<DistDvec, 2>::compute() {
  const int nops = 4;
  using DistBranch = GammaBranch<DistDvec>;
  using DistTree = GammaTree<DistDvec>;

  // Allocate memory while counting tasks
  for (auto& iforest : forests_) {
    for (auto& itreemap : iforest) {
      shared_ptr<DistTree> itree = itreemap.second;
      const int nA = itree->ket()->ij();
      const int norb = itree->ket()->det()->norb();

      // Allocation sweep
      for (int i = 0; i < nops; ++i) {
        shared_ptr<DistBranch> first = itree->base()->branch(i);
        if (!first->active()) continue;
        for (auto& ibra : first->bras()) {
          const int nAp = ibra.second->ij();
          const int nstates = nA * nAp;
          first->gammas().emplace(ibra.first, make_shared<Matrix>(nstates, norb));
        }

        for (int j = 0; j < nops; ++j) {
          shared_ptr<DistBranch> second = first->branch(j);
          if (!second->active()) continue;
          for (auto& jbra : second->bras()) {
            const int nAp = jbra.second->ij();
            const int nstates = nA * nAp;
            second->gammas().emplace(jbra.first, make_shared<Matrix>(nstates, norb * norb));
          }

          for (int k = 0; k < nops; ++k) {
            shared_ptr<DistBranch> third = second->branch(k);
            if (!third->active()) continue;
            for (auto& kbra : third->bras()) {
              const int nAp = kbra.second->ij();
              const int nstates = nA * nAp;
              third->gammas().emplace(kbra.first, make_shared<Matrix>(nstates, norb * norb * norb));
            }
          }
        }
      }
    }
  }

  // Compute tasks
  for (auto& iforest : forests_) {
    for (auto& itreemap : iforest) {
      shared_ptr<DistTree> itree = itreemap.second;

      const int norb = itree->ket()->det()->norb();
      for (int i = 0; i < nops; ++i) {
        shared_ptr<DistBranch> first = itree->base()->branch(i);
        if (!first->active()) continue;
        for (int a = 0; a < norb; ++a) {
          GammaTask<DistDvec> task(itree, GammaSQ(i), a);
          task.compute();
        }
      }   
    }   
  }   
}
