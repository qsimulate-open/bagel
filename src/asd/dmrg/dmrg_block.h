//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dmrg_block.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#ifndef BAGEL_ASD_DMRG_DMRG_BLOCK_H
#define BAGEL_ASD_DMRG_DMRG_BLOCK_H

#include <algorithm>

#include <src/asd/dmrg/gamma_forest_asd.h>
#include <src/asd/dmrg/gamma_forest_prod_asd.h>
#include <src/asd/dmrg/block_operators.h>
#include <src/asd/dmrg/kronecker.h>

namespace bagel {

struct CouplingBlock {
  std::pair<BlockInfo, BlockInfo> states;
  std::shared_ptr<btas::Tensor3<double>> data;
  CouplingBlock(BlockInfo ikey, BlockInfo jkey, std::shared_ptr<btas::Tensor3<double>> d) : states({ikey,jkey}), data(d) {}
  std::pair<BlockKey, BlockKey> key() const { return {states.first.key(), states.second.key()}; }
};

class DMRG_Block {
  protected:
    std::set<BlockInfo> blocks_;
    std::shared_ptr<const Matrix> coeff_; ///< Coefficients for orbitals stored in DMRG_Block

  public:
    DMRG_Block() {}
    DMRG_Block(std::shared_ptr<const Matrix> c) : coeff_(c) {}

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

    virtual std::string block_info_to_string(const BlockKey k, const int state) const = 0;

    std::shared_ptr<const Matrix> coeff() const { return coeff_; }

    virtual std::shared_ptr<Matrix> spin(const BlockKey k) const = 0;
    virtual std::shared_ptr<Matrix> spin_lower(const BlockKey k) const = 0;
    virtual std::shared_ptr<Matrix> spin_raise(const BlockKey k) const = 0;

    virtual std::shared_ptr<const BlockOperators> compute_block_ops(std::shared_ptr<DimerJop> jop) const = 0;
};

/// Store matrix representations of all of the operators needed to apply the Hamiltonian to a tensor-product state
class DMRG_Block1 : public std::enable_shared_from_this<DMRG_Block1>, public DMRG_Block {
  protected:
    using SparseMap = std::map<std::list<GammaSQ>, std::map<std::pair<BlockKey, BlockKey>, CouplingBlock>>;

    SparseMap sparse_;
    std::map<BlockKey, std::shared_ptr<const Matrix>> H2e_;
    std::map<BlockKey, std::shared_ptr<const Matrix>> spin_;

  public:
    /// default constructor
    DMRG_Block1() { }

    /// constructor that takes an Rvalue reference to a GammaForestASD
    template <typename GammaForestType>
    DMRG_Block1(GammaForestType&& forest, const std::map<BlockKey, std::shared_ptr<const Matrix>> h2e,
                const std::map<BlockKey, std::shared_ptr<const Matrix>> spin, std::shared_ptr<const Matrix> coeff)
     : DMRG_Block(coeff), H2e_(h2e), spin_(spin) {
      Timer dmrgtime(2);
      // Build set of blocks
      for (auto& i : h2e) {
        assert(i.second->ndim() == i.second->mdim());
        blocks_.emplace(i.first.nelea, i.first.neleb, i.second->ndim());
      }
      dmrgtime.tick_print("prepare blocks");

      // initialize operator space
      for (auto o : forest.sparselist()) {
        std::list<GammaSQ> gammalist = std::get<0>(o);

        BlockInfo brakey = std::get<1>(o);
        BlockInfo ketkey = std::get<2>(o);
        assert(blocks_.count(brakey) && blocks_.count(ketkey) && forest.exist(brakey, ketkey, gammalist));

        std::shared_ptr<const Matrix> mat = forest.get(brakey, ketkey, gammalist);
        btas::CRange<3> range(brakey.nstates, ketkey.nstates, lrint(pow(norb(), gammalist.size())));
        // checking dimensions
        assert(mat->ndim() == range.extent(0)*range.extent(1));
        assert(mat->mdim() == range.extent(2));
        assert(mat->storage().size() == range.area());
        // convert this matrix to 3-tensor
        auto tensor = std::make_shared<btas::Tensor3<double>>(range, std::move(mat->storage()));
#ifdef HAVE_MPI_H
        mpi__->broadcast(tensor->data(), tensor->size(), 0);
        dmrgtime.tick_print("broadcast");
#endif
        // add matrix
        CouplingBlock cb(brakey, ketkey, tensor);
        sparse_[gammalist].emplace(cb.key(), cb);
      }
    }

    std::string block_info_to_string(const BlockKey k, const int state) const override;

    const std::map<std::pair<BlockKey, BlockKey>, CouplingBlock>& coupling(const std::list<GammaSQ>& l) const { return sparse_.at(l); }
    std::shared_ptr<const Matrix> h2e(const BlockKey& b) const { return H2e_.at(b); }

    std::shared_ptr<Matrix> spin(const BlockKey b) const override { return spin_.at(b)->copy(); }
    std::shared_ptr<Matrix> spin_lower(const BlockKey b) const override;
    std::shared_ptr<Matrix> spin_raise(const BlockKey b) const override;

    std::shared_ptr<const BlockOperators> compute_block_ops(std::shared_ptr<DimerJop> jop) const override;
};

namespace DMRG {
struct BlockPair {
  BlockInfo left;
  BlockInfo right;
  int offset;

  BlockPair(BlockInfo l, BlockInfo r, int o) : left(l), right(r), offset(o) {}

  int nstates() const { return left.nstates*right.nstates; }
};
}

/// Combines two DMRG_Blocks to make one double block
class DMRG_Block2 : public std::enable_shared_from_this<DMRG_Block2>, public DMRG_Block {
  protected:
    std::shared_ptr<const DMRG_Block1> left_block_;
    std::shared_ptr<const DMRG_Block1> right_block_;

    std::map<BlockKey, std::vector<DMRG::BlockPair>> pairmap_;

  public:
    ///default constructor
    DMRG_Block2() { }

    /// constructor taking a pair of DMRG_Block1 objects
    DMRG_Block2(std::shared_ptr<const DMRG_Block1> lb, std::shared_ptr<const DMRG_Block1> rb);

    std::string block_info_to_string(const BlockKey k, const int state) const override;

    const std::vector<DMRG::BlockPair>& blockpairs(BlockKey bk) const { return pairmap_.at(bk); }

    std::shared_ptr<const DMRG_Block1> left_block() const { return left_block_; }
    std::shared_ptr<const DMRG_Block1> right_block() const { return right_block_; }

    std::shared_ptr<Matrix> spin(const BlockKey b) const override;
    std::shared_ptr<Matrix> spin_lower(const BlockKey b) const override;
    std::shared_ptr<Matrix> spin_raise(const BlockKey b) const override;

    std::shared_ptr<const BlockOperators> compute_block_ops(std::shared_ptr<DimerJop> jop) const override;
};

}

#endif
