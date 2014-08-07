//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_forest_one.h
// Copyright (C) 2014 Shane Parker
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

#ifndef __ASD_GAMMA_FOREST_ASD_H
#define __ASD_GAMMA_FOREST_ASD_H

#include <vector>
#include <list>

#include <src/asd/gamma_forest.h>
#include <src/asd_dmrg/block_key.h>

namespace bagel {

template <class VecType>
class GammaForestASD : public GammaForest<VecType, 1> {
  protected:
    using SparseList = std::list<std::tuple<std::list<GammaSQ>, BlockInfo, BlockInfo>>;
    SparseList sparselist_;

  public:
    GammaForestASD(std::map<BlockKey, std::shared_ptr<const VecType>> monomer_states) {
      std::vector<std::list<GammaSQ>> possible_couplings = {
        {GammaSQ::CreateAlpha},
        {GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha},
        {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},
        {GammaSQ::AnnihilateAlpha,  GammaSQ::AnnihilateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha},
        {GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha},
        {GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateAlpha},
      };

      for (auto& i : monomer_states) {
        for (auto& coupling : possible_couplings) {
          int new_nelea = i.first.nelea;
          int new_neleb = i.first.neleb;
          for (GammaSQ& operation : coupling) {
            switch (operation) {
              case GammaSQ::AnnihilateAlpha:
                --new_nelea; break;
              case GammaSQ::CreateAlpha:
                ++new_nelea; break;
              case GammaSQ::AnnihilateBeta:
                --new_neleb; break;
              case GammaSQ::CreateBeta:
                ++new_neleb; break;
            }
          }

          for (auto& j : monomer_states) {
            if (BlockKey(new_nelea, new_neleb)==j.first) {
#ifdef DEBUG
              std::cout << "inserting: <" << j.first.nelea << ", " << j.first.neleb << "|";
              for (auto opiter = coupling.rbegin(); opiter != coupling.rend(); ++opiter) {
                GammaSQ operation = *opiter;
                if ( operation == GammaSQ::AnnihilateAlpha || operation == GammaSQ::CreateAlpha )
                  std::cout << "(A)";
                else
                  std::cout << "(B)";
                if ( operation == GammaSQ::CreateAlpha || operation == GammaSQ::CreateBeta )
                  std::cout << "^t";
              }
              std::cout << "|" << i.first.nelea << ", " << i.first.neleb << ">" << std::endl;
#endif
              sparselist_.emplace_back(coupling, BlockInfo(j.first.nelea, j.first.neleb, j.second->ij()), BlockInfo(i.first.nelea, i.first.neleb, i.second->ij()));
              this->template insert<0>(i.second, block_tag(i.first), j.second, block_tag(j.first), coupling);
            }
          }
        }
      }

      this->compute();

#ifdef DEBUG
      this->for_each_branch([] (std::shared_ptr<GammaBranch<VecType>> b) { for (auto& i : b->gammas()) i.second->print(); });
#endif
    }

    SparseList sparselist() const { return sparselist_; }

    int block_tag(BlockKey b) const { return b.nelea + (b.neleb << 10); }

};

}

#endif
