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

#include <src/meh/gamma_tensor.h>
#include <src/asd_dmrg/gamma_forest_asd.h>

namespace bagel {

/// Store matrix representations of all of the operators needed to apply the Hamiltonian to a tensor-product state
class DMRG_Block {
  protected:
    struct CouplingBlock {
      std::pair<BlockInfo, BlockInfo> states;
      std::shared_ptr<btas::Tensor3<double>> data;
      CouplingBlock(BlockInfo ikey, BlockInfo jkey, std::shared_ptr<btas::Tensor3<double>> d) : states({ikey,jkey}), data(d) {}
    };
    using SparseMap = std::map<std::list<GammaSQ>, std::vector<CouplingBlock>>;

    SparseMap sparse_;
    std::map<BlockKey, std::shared_ptr<const Matrix>> H2e_;
    std::set<BlockInfo> blocks_;

    int norb_;

  public:
    /// default constructor
    DMRG_Block() { }

    /// constructor that takes an Rvalue reference to a GammaForestASD
    template <typename T>
    DMRG_Block(GammaForestASD<T>&& forest, const std::map<BlockKey, std::shared_ptr<const Matrix>> h2e, const int norb) : H2e_(h2e), norb_(norb) {
      // Build set of blocks
      for (auto& i : h2e) {
        assert(i.second->ndim() == i.second->mdim());
        blocks_.emplace(i.first.nelea, i.first.neleb, i.second->ndim());
      }

      // initialize operator space
      for (auto o : forest.sparselist()) {
        std::list<GammaSQ> gammalist = std::get<0>(o);

        BlockInfo ikey = std::get<1>(o);
        BlockInfo jkey = std::get<2>(o);

        assert(blocks_.count(ikey));
        assert(blocks_.count(jkey));

        const int itag = forest.block_tag(ikey);
        const int jtag = forest.block_tag(jkey);

        assert(forest.template exist<0>(itag, jtag, gammalist));
        std::shared_ptr<const Matrix> mat = forest.template get<0>(itag, jtag, gammalist);
        btas::CRange<3> range(ikey.nstates, jkey.nstates, std::lrint(std::pow(norb_, gammalist.size())));
        // checking ndim
        assert(mat->ndim() == range.extent(0)*range.extent(1));
        // checking mdim
        assert(mat->mdim() == range.extent(2));
        // internal check
        assert(mat->storage().size() == range.area());
        // convert this matrix to 3-tensor
        auto tensor = std::make_shared<btas::Tensor3<double>>(range, std::move(mat->storage()));
        // add matrix
        sparse_[gammalist].emplace_back(ikey, jkey, tensor);
      }
    }
};

}

#endif
