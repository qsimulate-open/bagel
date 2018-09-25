//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gamma_forest.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#include <src/asd/gamma_forest.h>
#include <src/ci/ras/apply_block.h>
#include <src/util/taskqueue.h>
#include <src/asd/coupling.h>
#include <src/asd/gamma_task.h>

using namespace std;
using namespace bagel;

template<typename VecType, int N>
int GammaForest<VecType,N>::allocate_and_count() {
  constexpr int nops = 4;
  int ntasks = 0;
  // Allocate memory while counting tasks
  for (auto& iforest : forests_) {
    for (auto& itreemap : iforest) {
      shared_ptr<GammaTree<VecType>> itree = itreemap.second;
      const int nA = itree->ket()->ij();
      const int norb = itree->norb();
      for (auto& basebra : itree->base()->bras())
        itree->base()->gammas().emplace(basebra.first, make_shared<Matrix>(nA*basebra.second->ij(), 1));

      // Allocation sweep
      for (int i = 0; i < nops; ++i) {
        shared_ptr<GammaBranch<VecType>> first = itree->base()->branch(i);
        if (!first->active()) continue;
        ++ntasks;
        for (auto& ibra : first->bras()) {
          const int nAp = ibra.second->ij();
          const int nstates = nA * nAp;
          first->gammas().emplace(ibra.first, make_shared<Matrix>(nstates, norb));
        }

        for (int j = 0; j < nops; ++j) {
          shared_ptr<GammaBranch<VecType>> second = first->branch(j);
          if (!second->active()) continue;
          for (auto& jbra : second->bras()) {
            const int nAp = jbra.second->ij();
            const int nstates = nA * nAp;
            second->gammas().emplace(jbra.first, make_shared<Matrix>(nstates, norb * norb));
          }

          for (int k = 0; k < nops; ++k) {
            shared_ptr<GammaBranch<VecType>> third = second->branch(k);
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
  return ntasks;
}


template<typename VecType, int N>
void GammaForest<VecType,N>::compute() {
  constexpr int nops = 4;
  const int ntasks = allocate_and_count();
  TaskQueue<GammaTask<VecType>> tasks(ntasks);

  // Add tasks
  for (auto& iforest : forests_) {
    for (auto& itreemap : iforest) {
      shared_ptr<GammaTree<VecType>> itree = itreemap.second;
      const int nket = itree->ket()->ij();
      for (auto& brapair : itree->base()->bras()) {
        double* target = itree->base()->gammas().at(brapair.first)->data();
        const int nbra = brapair.second->ij();
        for (int k = 0; k < nket; ++k)
          for (int b = 0; b < nbra; ++b, ++target)
            *target = brapair.second->data(b)->dot_product(*itree->ket()->data(k));
      }

      const int norb = itree->norb();
      for (int i = 0; i < nops; ++i) {
        shared_ptr<GammaBranch<VecType>> first = itree->base()->branch(i);
        if (!first->active()) continue;
        for (int a = 0; a < norb; ++a) tasks.emplace_back(itree, GammaSQ(i), a);
      }
    }
  }
  tasks.compute();
}


template <class VecType, int N>
void GammaForest<VecType,N>::couple_blocks(const DimerSubspace<VecType>& AB, const DimerSubspace<VecType>& ApBp) {
  if (N != 2)
    throw logic_error("GammaForest<VecType,N>::couple_blocks assumes N == 2");

  Coupling term_type = coupling_type(AB, ApBp);

  auto* space1 = &AB;   // bra
  auto* space2 = &ApBp; // ket

  bool flip = (static_cast<int>(term_type) < 0);

  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    swap(space1,space2);
  }

  // Throughout, I'm going to consider space1 to be the bra state, space2 to be the ket state
  if (term_type != Coupling::none) {
    const array<int,2> bra_tags {{space1->template tag<0>(), space1->template tag<1>()}};
    const array<int,2> ket_tags {{space2->template tag<0>(), space2->template tag<1>()}};

    shared_ptr<const VecType> bra_A = space1->template ci<0>();
    shared_ptr<const VecType> bra_B = space1->template ci<1>();

    shared_ptr<const VecType> ket_A = space2->template ci<0>();
    shared_ptr<const VecType> ket_B = space2->template ci<1>();

    switch(term_type) {
      case Coupling::none :
        assert(false); // Control should never be able to reach here
        break;
      case Coupling::diagonal :
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});
        break;
      case Coupling::aET : // Alpha ET
        // One-body aET
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateAlpha});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::AnnihilateAlpha});

        //Two-body aET, type 1
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateAlpha});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});

        //Two-body aET, type 2
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::AnnihilateAlpha});
        break;
      case Coupling::bET : // Beta ET
        // One-body bET
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateBeta});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::AnnihilateBeta});
        //Two-body bET, type 1
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateBeta});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});

        //Two-body aET, type 2
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::AnnihilateBeta});
        break;
      case Coupling::abFlip :
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta});
        break;
      case Coupling::abET :
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateAlpha, GammaSQ::CreateBeta});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});
        break;
      case Coupling::aaET :
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});
        break;
      case Coupling::bbET :
        this->insert<0>(bra_A, bra_tags[0], ket_A, ket_tags[0], {GammaSQ::CreateBeta, GammaSQ::CreateBeta});
        this->insert<1>(bra_B, bra_tags[1], ket_B, ket_tags[1], {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});
        break;
      default :
        assert(false); break; // Control should never reach here
    }
  }
}


template class bagel::GammaForest<CASDvec,1>;
template class bagel::GammaForest<CASDvec,2>;
template class bagel::GammaForest<CASDvec,3>;
template class bagel::GammaForest<CASDvec,4>;
template class bagel::GammaForest<CASDvec,5>;
template class bagel::GammaForest<CASDvec,6>;
template class bagel::GammaForest<CASDvec,7>;
template class bagel::GammaForest<RASDvec,1>;
template class bagel::GammaForest<RASDvec,2>;
template class bagel::GammaForest<RASDvec,3>;
template class bagel::GammaForest<RASDvec,4>;
template class bagel::GammaForest<RASDvec,5>;
template class bagel::GammaForest<RASDvec,6>;
template class bagel::GammaForest<RASDvec,7>;
