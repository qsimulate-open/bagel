//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gamma_tree.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Toru Shiozaki
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

#ifndef __ASD_GAMMA_TREE_H
#define __ASD_GAMMA_TREE_H

#include <set>
#include <src/ci/fci/civec.h>
#include <src/ci/ras/civector.h>
#include <src/util/math/matrix.h>
#include <src/asd/gamma_sq.h>

namespace bagel {

template <typename VecType>
class GammaBranch {
  protected:
    std::array<std::shared_ptr<GammaBranch<VecType>>, 4> branches_;
    std::map<size_t, std::shared_ptr<const VecType>> bras_; // use tags as unique identifier
    std::map<size_t, std::shared_ptr<Matrix>> gammas_;

    bool active_;

  public:
    GammaBranch() : active_(false) {}

    void insert(std::shared_ptr<const VecType> bra, const size_t bra_tag, const std::list<GammaSQ>& gsq);

    std::shared_ptr<Matrix> search(const size_t tag, const std::list<GammaSQ>& gsq);
    std::shared_ptr<const Matrix> search(const size_t tag, const std::list<GammaSQ>& gsq) const;

    bool exist(const size_t tag, const std::list<GammaSQ>& gsq) const;

    void activate() { active_ = true; }
    bool active() const { return active_; }

    std::shared_ptr<GammaBranch<VecType>>& branch(const int i) { return branches_[i]; }
    const std::shared_ptr<GammaBranch<VecType>>& branch(const int i) const { return branches_[i]; }
    std::shared_ptr<GammaBranch<VecType>>& branch(const GammaSQ sq) { return branch(static_cast<int>(sq)); }
    const std::shared_ptr<GammaBranch<VecType>>& branch(const GammaSQ sq) const { return branch(static_cast<int>(sq)); }

    const std::map<size_t, std::shared_ptr<const VecType>>& bras() const { return bras_; }
    std::map<size_t, std::shared_ptr<Matrix>>& gammas() { return gammas_; }

    bool if_contributes(std::set<int> needed);

    // distterminal means there are no alpha operators higher in the tree
    bool is_distterminal();
};

template <typename VecType>
class GammaTree {
  protected:
    std::shared_ptr<const VecType> ket_;
    std::shared_ptr<GammaBranch<VecType>> base_;

  public:
    GammaTree(std::shared_ptr<const VecType> ket);

    std::shared_ptr<GammaBranch<VecType>> base() { return base_; }

    void insert(std::shared_ptr<const VecType> bra, const size_t tag, const std::list<GammaSQ>& ops) { base_->insert(bra, tag, ops); }

    std::shared_ptr<Matrix> search(const size_t tag, const std::list<GammaSQ>& address) { return base_->search(tag, address); }
    std::shared_ptr<const Matrix> search(const size_t tag, const std::list<GammaSQ>& address) const { return base_->search(tag, address); }

    bool exist(const size_t tag, const std::list<GammaSQ>& address) const { return base_->exist(tag, address); }

    std::shared_ptr<const VecType> ket() const { return ket_; }

    int norb() const { return ket()->det()->norb(); }
};

extern template class GammaBranch<CASDvec>;
extern template class GammaBranch<RASDvec>;
extern template class GammaTree<CASDvec>;
extern template class GammaTree<RASDvec>;

}

#endif
