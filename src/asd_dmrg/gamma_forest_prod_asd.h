//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_forest_prod_asd.h
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

#ifndef __BAGEL_ASD_DMRG_GAMMA_FOREST_PROD_ASD_H
#define __BAGEL_ASD_DMRG_GAMMA_FOREST_PROD_ASD_H

#include <vector>
#include <list>

#include <src/asd/gamma_forest.h>
#include <src/asd_dmrg/block_key.h>

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

    bool exist(BlockKey bra, BlockKey ket, std::list<GammaSQ> gammalist) const {
      auto iter = std::find_if(sparselist_.begin(), sparselist_.end(), [&bra, &ket, &gammalist] (std::tuple<std::list<GammaSQ>, BlockInfo, BlockInfo> sparse) {
        return std::make_tuple(bra, ket, gammalist) == std::make_tuple(std::get<1>(sparse).key(), std::get<2>(sparse).key(), std::get<0>(sparse)); });
      return iter != sparselist_.end();
    }
    std::shared_ptr<const Matrix> get(BlockKey bra, BlockKey ket, std::list<GammaSQ> gammalist) const { return gammas_.at(make_tuple(gammalist,bra,ket)); }
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
      return b.block.nelea + b.block.neleb*nele_ + (b.ci.nelea + b.ci.neleb*nele_ + b.state*nele_*nele_)*nele_*nele_;
    }

    std::tuple</*conj*/bool, /*rev*/bool, std::list<GammaSQ>> try_permutations(const std::list<GammaSQ>& gammalist) const;
    std::tuple</*full_ijk*/int, /*block_index*/int, /*ci_index*/int> get_indices(std::bitset<3> bit, int size, int ijk_local, int lnorb, bool block_is_reversed, int rnorb, bool ci_is_reversed) const;
};

}

#endif
