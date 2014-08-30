//
// BAGEL - Parallel electron correlation program.
// Filename: dmrg_block.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef BAGEL_ASD_DMRG_DMRG_BLOCK_H
#define BAGEL_ASD_DMRG_DMRG_BLOCK_H

#include <algorithm>

#include <src/asd_dmrg/gamma_forest_asd.h>
#include <src/asd_dmrg/gamma_forest_prod_asd.h>
#include <src/ras/civector.h>

namespace bagel {

struct CouplingBlock {
  std::pair<BlockInfo, BlockInfo> states;
  std::shared_ptr<btas::Tensor3<double>> data;
  CouplingBlock(BlockInfo ikey, BlockInfo jkey, std::shared_ptr<btas::Tensor3<double>> d) : states({ikey,jkey}), data(d) {}
  std::pair<BlockKey, BlockKey> key() const { return {states.first.key(), states.second.key()}; }
};

/// Store matrix representations of all of the operators needed to apply the Hamiltonian to a tensor-product state
class DMRG_Block {
  protected:
    using SparseMap = std::map<std::list<GammaSQ>, std::map<std::pair<BlockKey, BlockKey>, CouplingBlock>>;

    SparseMap sparse_;
    std::map<BlockKey, std::shared_ptr<const Matrix>> H2e_;
    std::map<BlockKey, std::shared_ptr<const Matrix>> spin_;
    std::set<BlockInfo> blocks_;

    std::shared_ptr<const Matrix> coeff_; ///< Coefficients for orbitals stored in DMRG_Block

  public:
    /// default constructor
    DMRG_Block() { }

    /// constructor that takes an Rvalue reference to a GammaForestASD
    DMRG_Block(GammaForestASD<RASDvec>&& forest, const std::map<BlockKey, std::shared_ptr<const Matrix>> h2e, const std::map<BlockKey,
               std::shared_ptr<const Matrix>> spin, std::shared_ptr<const Matrix> coeff);
    DMRG_Block(GammaForestProdASD&& forest, const std::map<BlockKey, std::shared_ptr<const Matrix>> h2e, const std::map<BlockKey,
               std::shared_ptr<const Matrix>> spin, std::shared_ptr<const Matrix> coeff);

    bool contains(const BlockKey& k) const {
      auto iter = std::find_if(blocks_.begin(), blocks_.end(), [&k] (const BlockInfo& b) { return (b.nelea==k.nelea && b.neleb==k.neleb); });
      return iter!=blocks_.end();
    }

    int norb() const { return coeff_->mdim(); }
    int nstates() const { return std::accumulate(blocks_.begin(), blocks_.end(), 0, [] (const int o, const BlockInfo& b) { return o + b.nstates; }); }

    std::set<BlockInfo>& blocks() { return blocks_; }
    const std::set<BlockInfo>& blocks() const { return blocks_; }

    BlockInfo blockinfo(const BlockKey k) const {
      auto iter = std::find_if(blocks_.begin(), blocks_.end(), [&k] (const BlockInfo& b) { return (b.nelea==k.nelea && b.neleb==k.neleb); });
      assert(iter!=blocks_.end());
      return *iter;
    }

    const std::map<std::pair<BlockKey, BlockKey>, CouplingBlock>& coupling(const std::list<GammaSQ>& l) const { return sparse_.at(l); }
    std::shared_ptr<const Matrix> h2e(const BlockKey& b) const { return H2e_.at(b); }
    std::shared_ptr<const Matrix> spin(const BlockKey& b) const { return spin_.at(b); }

    std::shared_ptr<const Matrix> coeff() const { return coeff_; }
};

}

#endif
