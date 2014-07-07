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
    struct listGammaSQ {
      std::list<GammaSQ> list;
      int size;
      listGammaSQ(std::list<GammaSQ>& l, int s) : list(l), size(s) { }
      listGammaSQ(std::list<GammaSQ>&& l, int s) : list(l), size(s) { }

      bool operator<(const listGammaSQ& o) const { return std::make_pair(list,size) < std::make_pair(o.list,o.size); }
      bool operator==(const listGammaSQ& o) const { return std::make_pair(list,size) == std::make_pair(o.list,o.size); }

      int calculate_norb() const {
        const double o = std::pow(size, 1.0/list.size());
        const long int n = std::lrint(o);
        assert(std::fabs(o - n) < 1.0e-15);
        return n;
      }
    };

    using SparseMap     = std::map<std::tuple<listGammaSQ, MonomerKey, MonomerKey>, std::shared_ptr<Matrix>>;
    SparseMap sparse_;
    int norb_;

    static const std::list<std::list<GammaSQ>> oplist_;

  public:
    // default constructor
    GammaTensor() { }

    // constructor that takes GammaForst
    template <typename T, int N, int M>
    GammaTensor(asd::Wrap<GammaForest<T,N>,M> forest, const std::vector<DimerSubspace<T>>& sp) {
      std::shared_ptr<GammaForest<T,N>> f = forest.f;

      norb_ = f->norb();
      // initialize operator space
      for (auto o : oplist_) {
        for (auto& i : sp)
          for (auto& j : sp)
            if (f->template exist<M>(i.template tag<M>(), j.template tag<M>(), o)) {
              std::shared_ptr<Matrix> mat = f->template get<M>(i.template tag<M>(), j.template tag<M>(), o);
              // checking ndim
              assert(mat->ndim() == i.template monomerkey<M>().nstates()*j.template monomerkey<M>().nstates());
              // checking mdim
              assert(mat->mdim() == std::lrint(std::pow(norb_, o.size())));
              // add matrix
              sparse_.emplace(std::make_tuple(listGammaSQ(o, mat->mdim()), i.template monomerkey<M>(), j.template monomerkey<M>()), mat);
            }
      }
    }

    // constructor that takes sparsity information
    GammaTensor(const std::list<std::tuple<listGammaSQ, MonomerKey, MonomerKey>>& o) {
      norb_ = std::get<0>(o.front()).calculate_norb();
      for (auto& i : o) {
        assert(norb_ == std::get<0>(i).calculate_norb());
        sparse_.emplace(i, std::make_shared<Matrix>(std::get<0>(i).size, std::get<1>(i).nstates()*std::get<2>(i).nstates()));
      }
    }

    // copy constructor
    GammaTensor(const GammaTensor& o) : norb_(o.norb_) {
      for (auto& i : o.sparse_)
        sparse_.emplace(i.first, i.second->copy());
    }
    // move constructor
    GammaTensor(GammaTensor&& o) : norb_(o.norb_) {
      for (auto& i : o.sparse_)
        sparse_.emplace(i.first, i.second);
    }

    // copy and move assignment operators
    GammaTensor& operator=(const GammaTensor& o);
    GammaTensor& operator=(GammaTensor&& o);

    std::shared_ptr<GammaTensor> clone() const { return std::make_shared<GammaTensor>(blocks()); }
    std::shared_ptr<GammaTensor> copy() const { return std::make_shared<GammaTensor>(*this); }

    auto begin() -> decltype(sparse_.begin()) { return sparse_.begin(); }
    auto end() -> decltype(sparse_.end()) { return sparse_.end(); }
    auto begin() const -> decltype(sparse_.cbegin()) { return sparse_.cbegin(); }
    auto end() const -> decltype(sparse_.cend()) { return sparse_.cend(); }
    auto cbegin() const -> decltype(sparse_.cbegin()) { return sparse_.cbegin(); }
    auto cend() const -> decltype(sparse_.cend()) { return sparse_.cend(); }

    int nblocks() const { return sparse_.size(); }

    template <typename... Args>
    auto emplace(Args... args) -> decltype(sparse_.emplace(std::forward<Args>(args)...)) { return sparse_.emplace(std::forward<Args>(args)...); }

    std::shared_ptr<const Matrix> get_block(const MonomerKey& i, const MonomerKey& j, const std::initializer_list<GammaSQ>& o) const {
      return sparse_.at(std::make_tuple(listGammaSQ(std::list<GammaSQ>(o), std::lrint(std::pow(norb_, o.size()))), i, j));
    }

    std::list<std::tuple<listGammaSQ, MonomerKey, MonomerKey>> blocks() const {
      std::list<std::tuple<listGammaSQ, MonomerKey, MonomerKey>> out;
      for (auto& i : sparse_) out.push_back(i.first);
      return out;
    }
};

}

#endif
