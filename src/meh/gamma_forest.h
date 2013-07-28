//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_forest.h
// Copyright (C) 2013 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: NU theory
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#ifndef __meh_gamma_forest_h
#define __meh_gamma_forest_h

#include <utility>

#include <src/fci/dvec.h>
#include <src/math/matrix.h>

namespace bagel {

enum class GammaSQ {
  CreateAlpha = 0,
  AnnihilateAlpha = 1,
  CreateBeta = 2,
  AnnihilateBeta = 3
};

class GammaBranch {
  protected:
    std::array<std::shared_ptr<GammaBranch>, 4> branches_;
    std::map<int, std::shared_ptr<const Dvec>> bras_; // use offsets in the hamiltonian as a unique identifier
    std::map<int, std::shared_ptr<Matrix>> gammas_;

    bool active_;

  public:
    GammaBranch() : active_(false) {}

    template <class ...GSQs>
    void insert(std::shared_ptr<const Dvec> bra, const int offset, const GammaSQ first, const GSQs... rest) {
      std::shared_ptr<GammaBranch> target = branches_[static_cast<int>(first)];

      target->activate();
      target->insert(bra, offset, rest...);
    }
    void insert(std::shared_ptr<const Dvec> bra, const int offset) { bras_.insert(std::make_pair(offset,bra)); }

    template <class ...GSQs>
    std::shared_ptr<const Matrix> search(const int offset, GammaSQ first, GSQs... rest) const { return branch(first)->search(offset, rest...); }
    std::shared_ptr<const Matrix> search(const int offset) const {
      assert(gammas_.find(offset)!=gammas_.end()); return gammas_.find(offset)->second;
    }

    void activate() { active_ = true; }
    bool active() const { return active_; }

    std::shared_ptr<GammaBranch>& branch(const int i) { return branches_[i]; }
    const std::shared_ptr<GammaBranch>& branch(const int i) const { return branches_[i]; }
    std::shared_ptr<GammaBranch>& branch(const GammaSQ sq) { return branch(static_cast<int>(sq)); }
    const std::shared_ptr<GammaBranch>& branch(const GammaSQ sq) const { return branch(static_cast<int>(sq)); }

    const std::map<int, std::shared_ptr<const Dvec>>& bras() const { return bras_; }
    std::map<int, std::shared_ptr<Matrix>>& gammas() { return gammas_; }
};

class GammaTree {
  protected:
    std::shared_ptr<const Dvec> ket_;
    std::shared_ptr<GammaBranch> base_;

  public:
    GammaTree(std::shared_ptr<const Dvec> ket);

    std::shared_ptr<GammaBranch> base() { return base_; }

    template <class ...GSQs>
    void insert(std::shared_ptr<const Dvec> bra, const int offset, const GSQs... ops) { base_->insert(bra, offset, ops...); }

    template <class ...GSQs>
    std::shared_ptr<const Matrix> search(const int offset, GSQs... address) const { return base_->search(offset, address...); }

    void compute();

  private:
    static const int Alpha = 0;
    static const int Beta = 1;
    static const int Create = 0;
    static const int Annihilate = 1;

    std::shared_ptr<const Dvec> apply(std::shared_ptr<const Dvec> kets, const GammaSQ operation, const int orbital) const;
    std::unique_ptr<double[]> ddot(std::shared_ptr<const Dvec> bras, std::shared_ptr<const Dvec> kets) const;
};

template <int N> class GammaForest {
  protected:
    std::array<std::map<int, std::shared_ptr<GammaTree>>, N> forests_;

  public:
    GammaForest() {};

    void compute() {
      for (auto& iforest : forests_) {
        for (auto& itree : iforest) {
          itree.second->compute();
        }
      }
    }

    template <int unit, class ...operations>
    void insert(std::shared_ptr<const Dvec> ket, const int ioffset, std::shared_ptr<const Dvec> bra, const int joffset, const operations... ops) {
      std::shared_ptr<GammaTree> gtree = tree<unit>(ket, ioffset);
      gtree->insert(bra, joffset, ops...);
    }

    template <int unit, class ...operations>
    std::shared_ptr<const Matrix> get(const int ioffset, const int joffset, const operations... ops) const {
      auto itree = forests_[unit].find(ioffset); assert(itree!=forests_[unit].end());
      return itree->second->search(joffset, ops...);
    }

  private:
    template <int unit>
    std::shared_ptr<GammaTree> tree(std::shared_ptr<const Dvec> ket, const int ioffset) {
      auto itree = forests_[unit].find(ioffset);
      if (itree == forests_[unit].end()) {
        forests_[unit].insert(std::make_pair(ioffset, std::make_shared<GammaTree>(ket)));
        itree = forests_[unit].find(ioffset);
      }
      return itree->second;
    }
};

}

#endif
