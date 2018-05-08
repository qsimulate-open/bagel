//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gamma_forest.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
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

#ifndef __ASD_GAMMA_FOREST_H
#define __ASD_GAMMA_FOREST_H

#include <src/asd/gamma_tree.h>
#include <src/asd/dimersubspace.h>

namespace bagel {

template <typename VecType, int N>
class GammaForest {
  protected:
    std::array<std::map<size_t, std::shared_ptr<GammaTree<VecType>>>, N> forests_;

  public:
    GammaForest() {}

    template <int unit>
    void insert(std::shared_ptr<const VecType> bra, const size_t bra_tag, std::shared_ptr<const VecType> ket, const size_t ket_tag, const std::list<GammaSQ>& ops) {
      std::shared_ptr<GammaTree<VecType>> gtree = tree<unit>(ket, ket_tag);
      gtree->insert(bra, bra_tag, ops);
    }

    template <int unit>
    std::shared_ptr<Matrix> get(const size_t bra_tag, const size_t ket_tag, const std::list<GammaSQ>& ops) {
      auto itree = forests_[unit].find(ket_tag); assert(itree!=forests_[unit].end());
      return itree->second->search(bra_tag, ops);
    }

    template <int unit>
    std::shared_ptr<const Matrix> get(const size_t bra_tag, const size_t ket_tag, const std::list<GammaSQ>& ops) const {
      auto itree = forests_[unit].find(ket_tag); assert(itree!=forests_[unit].end());
      return itree->second->search(bra_tag, ops);
    }

    template <int unit>
    bool exist(const size_t bra_tag, const size_t ket_tag, const std::list<GammaSQ>& ops) const {
      auto itree = forests_[unit].find(ket_tag);
      return itree == forests_[unit].end() ? false : itree->second->exist(bra_tag, ops);
    }

    int allocate_and_count();
    void compute();
    void couple_blocks(const DimerSubspace<VecType>& AB, const DimerSubspace<VecType>& ABp);

  private:
    template <int unit>
    std::shared_ptr<GammaTree<VecType>> tree(std::shared_ptr<const VecType> ket, const size_t ket_tag) {
      typename std::map<size_t, std::shared_ptr<GammaTree<VecType>>>::iterator itree = forests_[unit].find(ket_tag);
      if (itree == forests_[unit].end()) {
        forests_[unit].emplace(ket_tag, std::make_shared<GammaTree<VecType>>(ket));
        itree = forests_[unit].find(ket_tag);
      }
      return itree->second;
    }
};

extern template class GammaForest<CASDvec,1>;
extern template class GammaForest<CASDvec,2>;
extern template class GammaForest<CASDvec,3>;
extern template class GammaForest<CASDvec,4>;
extern template class GammaForest<CASDvec,5>;
extern template class GammaForest<CASDvec,6>;
extern template class GammaForest<CASDvec,7>;
extern template class GammaForest<RASDvec,1>;
extern template class GammaForest<RASDvec,2>;
extern template class GammaForest<RASDvec,3>;
extern template class GammaForest<RASDvec,4>;
extern template class GammaForest<RASDvec,5>;
extern template class GammaForest<RASDvec,6>;
extern template class GammaForest<RASDvec,7>;

}

#endif
