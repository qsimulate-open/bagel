//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gamma_forest_prod_asd.h
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

#ifndef __BAGEL_ASD_DMRG_GAMMA_FOREST_PROD_ASD_H
#define __BAGEL_ASD_DMRG_GAMMA_FOREST_PROD_ASD_H

#include <vector>
#include <list>

#include <src/asd/gamma_forest.h>
#include <src/asd/dmrg/block_key.h>

namespace bagel {

struct ProductState {
  BlockKey block;
  BlockKey ci;
  int state;
  ProductState(BlockKey b, BlockKey c, int s) : block(b), ci(c), state(s) {}
  BlockKey key() const { return BlockKey(block.nelea+ci.nelea, block.neleb+ci.neleb); }
  bool operator<(const ProductState& o) const { return std::make_tuple(block, ci, state) < std::make_tuple(o.block, o.ci, o.state); }
};

class ProductRASCivec;

class GammaForestProdASD {
  protected:
    std::set<std::tuple<std::list<GammaSQ>, BlockKey, BlockKey>> sparseset_;
    using SparseList = std::list<std::tuple<std::list<GammaSQ>, BlockInfo, BlockInfo>>;
    SparseList sparselist_;
    std::map<BlockKey, std::vector<std::shared_ptr<ProductRASCivec>>> block_states_;
    std::vector<std::list<GammaSQ>> possible_couplings_;

    std::shared_ptr<GammaForest<RASDvec, 1>> forest_;
    std::map<std::tuple<std::list<GammaSQ>, BlockKey, BlockKey>, std::shared_ptr<Matrix>> gammas_;

    size_t nele_;

  public:
    GammaForestProdASD(std::map<BlockKey, std::vector<std::shared_ptr<ProductRASCivec>>> block_states);

    void compute();

    bool exist(const BlockKey& bra, const BlockKey& ket, const std::list<GammaSQ>& gammalist) const {
      auto iter = std::find_if(sparselist_.begin(), sparselist_.end(), [&bra, &ket, &gammalist] (std::tuple<std::list<GammaSQ>, BlockInfo, BlockInfo> sparse) {
        return std::make_tuple(bra, ket, gammalist) == std::make_tuple(std::get<1>(sparse).key(), std::get<2>(sparse).key(), std::get<0>(sparse)); });
      return iter != sparselist_.end();
    }
    std::shared_ptr<const Matrix> get(const BlockKey& bra, const BlockKey& ket, const std::list<GammaSQ>& gammalist) const { return gammas_.at(make_tuple(gammalist,bra,ket)); }
    SparseList sparselist() const { return sparselist_; }

  private:
    BlockKey apply_key(BlockKey bk, const std::list<GammaSQ>& op) const {
      for (const GammaSQ& operation : op) {
        switch (operation) {
          case GammaSQ::AnnihilateAlpha:
            --bk.nelea; break;
          case GammaSQ::CreateAlpha:
            ++bk.nelea; break;
          case GammaSQ::AnnihilateBeta:
            --bk.neleb; break;
          case GammaSQ::CreateBeta:
            ++bk.neleb; break;
        }
      }
      return bk;
    }

    size_t state_tag(const ProductState& b) const {
      return static_cast<size_t>(b.block.nelea) + static_cast<size_t>(b.block.neleb)*nele_
             + (static_cast<size_t>(b.ci.nelea) + static_cast<size_t>(b.ci.neleb)*nele_ + static_cast<size_t>(b.state)*nele_*nele_)*nele_*nele_;
    }

    std::tuple</*conj*/bool, /*rev*/bool, std::list<GammaSQ>> try_permutations(const std::list<GammaSQ>& gammalist) const;
    std::tuple</*full_ijk*/size_t, /*block_index*/size_t, /*ci_index*/size_t> get_indices(std::bitset<3> bit, int size, size_t ijk_local, int lnorb, bool block_is_reversed, int rnorb, bool ci_is_reversed) const;
};

}

#endif
