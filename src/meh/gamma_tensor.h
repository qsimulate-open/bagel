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
#include <src/meh/dimersubspace.h>

namespace bagel {

namespace asd {
  template <typename T, int N>
  struct Wrap {
    std::shared_ptr<T> f;
    Wrap(std::shared_ptr<T> p) : f(p) { }
  };
}

// Store Gamma Tensor as a block-sparse tensor.
class GammaTensor {
  protected:
    using OperatorClass = std::list<GammaSQ>;

    using BlockSize     = std::map<SpaceKey, int>;
    using OperatorSize  = std::map<OperatorClass, int>;
    // if blocks exist of not. TODO has to be generalized if there are more than 2 active spaces
    using SparseMap     = std::map<std::tuple<OperatorClass, SpaceKey, SpaceKey>, std::shared_ptr<const Matrix>>;

    SparseMap    sparse_;
    std::list<std::list<GammaSQ>> oplist_;

  public:
    // constructor that takes GammaForst
    template <typename T, int N, int M>
    GammaTensor(asd::Wrap<GammaForest<T,N>,M> forest, const std::vector<DimerSubspace_base<T>>& sp) {
      std::shared_ptr<const GammaForest<T,N>> f = forest.f;
      // operator
      oplist_ = {
        {GammaSQ::CreateAlpha},
        {GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateAlpha},
        {GammaSQ::AnnihilateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha},
        {GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateBeta,  GammaSQ::CreateAlpha},
        {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha},
        {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},
        {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha},
        {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha},
        {GammaSQ::CreateBeta,  GammaSQ::CreateBeta},
        {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha},
        {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha},
        {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::CreateAlpha},
        {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta},
        {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha},
        {GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta}
      };
      // initialize operator space
      for (auto o : oplist_) {
        int dim = 0;
        for (auto& i : sp)
          for (auto& j : sp)
            if (f->template exist<M>(i.template tag<M>(), j.template tag<M>(), o)) {
              std::shared_ptr<const Matrix> mat = f->template get<M>(i.template tag<M>(), j.template tag<M>(), o);
              if (dim != 0) assert(mat->ndim() == dim);
              sparse_.emplace(std::make_tuple(o, i.template spacekey<M>(), j.template spacekey<M>()), mat);
            }
      }
    }

#if 0
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

    std::shared_ptr<const Matrix> get_block(const SpaceKey& i, const SpaceKey& j, const std::initializer_list<GammaSQ>& o) const { return sparse_[std::make_tuple(std::vector<GammaSQ>(o), i, j)]; }
#endif

};

}

#endif
