//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_tensor.h
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

#ifndef BAGEL_MEH_GAMMA_TENSOR_H
#define BAGEL_MEH_GAMMA_TENSOR_H

#include <src/meh/gamma_forest.h>

namespace bagel {

// Store Gamma Tensor as a block-sparse tensor.
class GammaTensor {
  protected:
    using OperatorClass = std::vector<GammaSQ>;

    using BlockSize     = std::map<SpaceKey, int>;
    using OperatorSize  = std::map<OperatorClass, int>;
    // if blocks exist of not. TODO has to be generalized if there are more than 2 active spaces
    using SparseMap     = std::map<std::tuple<OperatorClass, SpaceKey, SpaceKey>, std::shared_ptr<Matrix>>;

    SparseMap    sparse_;
    OperatorSize opsize_;
    BlockSize    blocksize_;

  public:
    // constructor that takes GammaForst
    template <typename T, int N>
    GammaTensor(const GammaForest<T,N>& f) {

    }

    // constructor that takes sparsity information
    GammaTensor(const SparseMap& o, const OperatorSize& op, const BlockSize& b) : opsize_(op), blocksize_(b) {
      for (auto& i : o)
        sparse_.emplace(i.first, i.second->clone());
    }

    // copy constructor
    GammaTensor(const GammaTensor& o) : opsize_(o.opsize_), blocksize_(o.blocksize_) {
      for (auto& i : o.sparse_)
        sparse_.emplace(i.first, i.second->copy());
    }

    std::shared_ptr<GammaTensor> clone() const { return std::make_shared<GammaTensor>(sparse_, opsize_, blocksize_); }
    std::shared_ptr<GammaTensor> copy() const { return std::make_shared<GammaTensor>(*this); }

};

}

#endif
