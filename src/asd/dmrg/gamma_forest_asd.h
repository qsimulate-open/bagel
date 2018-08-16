//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gamma_forest_asd.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#ifndef __ASD_DMRG_GAMMA_FOREST_ASD_H
#define __ASD_DMRG_GAMMA_FOREST_ASD_H

#include <vector>
#include <list>

#include <src/asd/gamma_forest.h>
#include <src/asd/dmrg/block_key.h>

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
        {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha},
        {GammaSQ::CreateAlpha,     GammaSQ::AnnihilateAlpha},
        {GammaSQ::CreateBeta,      GammaSQ::AnnihilateBeta},
        {GammaSQ::CreateBeta,      GammaSQ::AnnihilateAlpha},
        {GammaSQ::CreateBeta,      GammaSQ::CreateAlpha,     GammaSQ::AnnihilateAlpha},
        {GammaSQ::CreateAlpha,     GammaSQ::CreateAlpha,     GammaSQ::AnnihilateAlpha},
        {GammaSQ::CreateBeta,      GammaSQ::CreateBeta,      GammaSQ::AnnihilateBeta},
        {GammaSQ::CreateAlpha,     GammaSQ::CreateBeta,      GammaSQ::AnnihilateBeta},
      };

      for (auto& ket : monomer_states) {
        for (auto& coupling : possible_couplings) {
          int new_nelea = ket.first.nelea;
          int new_neleb = ket.first.neleb;
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

          for (auto& bra : monomer_states) {
            if (BlockKey(new_nelea, new_neleb)==bra.first) {
#ifdef DEBUG
              std::cout << "inserting: <" << bra.first.nelea << ", " << bra.first.neleb << "|";
              for (auto opiter = coupling.begin(); opiter != coupling.end(); ++opiter)
                std::cout << is_alpha(*opiter) ? "(A)" : "(B)" << is_creation(*opiter) ? "^t" : "";
              std::cout << "|" << ket.first.nelea << ", " << ket.first.neleb << ">" << std::endl;
#endif
              sparselist_.emplace_back(coupling, BlockInfo(bra.first.nelea, bra.first.neleb, bra.second->ij()), BlockInfo(ket.first.nelea, ket.first.neleb, ket.second->ij()));
              this->template insert<0>(bra.second, block_tag(bra.first), ket.second, block_tag(ket.first), coupling);
            }
          }
        }
      }
    }

    SparseList sparselist() const { return sparselist_; }

    size_t block_tag(const BlockKey& b) const { return b.nelea + (b.neleb << 10); }

    bool exist(const BlockKey& bra, const BlockKey& ket, const std::list<GammaSQ>& gammalist) const {
      const size_t bratag = block_tag(bra);
      const size_t kettag = block_tag(ket);
      return GammaForest<VecType,1>::template exist<0>(bratag, kettag, gammalist);
    }

    std::shared_ptr<const Matrix> get(const BlockKey& bra, const BlockKey& ket, const std::list<GammaSQ>& gammalist) const {
      const size_t bratag = block_tag(bra);
      const size_t kettag = block_tag(ket);
      return GammaForest<VecType,1>::template get<0>(bratag, kettag, gammalist);
    }

};

}

#endif
