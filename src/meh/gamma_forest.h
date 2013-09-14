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

#ifndef __meh_gamma_forest_h
#define __meh_gamma_forest_h

#include <array>

#include <src/fci/dvec.h>
#include <src/math/matrix.h>
#include <src/util/taskqueue.h>

namespace bagel {

enum class GammaSQ {
  CreateAlpha = 0,
  AnnihilateAlpha = 1,
  CreateBeta = 2,
  AnnihilateBeta = 3
};

template <typename VecType>
class GammaBranch {
  protected:
    std::array<std::shared_ptr<GammaBranch<VecType>>, 4> branches_;
    std::map<int, std::shared_ptr<const VecType>> bras_; // use offsets in the hamiltonian as a unique identifier
    std::map<int, std::shared_ptr<Matrix>> gammas_;

    bool active_;

  public:
    GammaBranch() : active_(false) {}

    template <class ...GSQs>
    void insert(std::shared_ptr<const VecType> bra, const int offset, const GammaSQ first, const GSQs... rest) {
      std::shared_ptr<GammaBranch<VecType>> target = branches_[static_cast<int>(first)];

      target->activate();
      target->insert(bra, offset, rest...);
    }
    void insert(std::shared_ptr<const VecType> bra, const int offset) { bras_.emplace(offset,bra); }

    template <class ...GSQs>
    std::shared_ptr<const Matrix> search(const int offset, GammaSQ first, GSQs... rest) const { return branch(first)->search(offset, rest...); }
    std::shared_ptr<const Matrix> search(const int offset) const {
      assert(gammas_.find(offset)!=gammas_.end()); return gammas_.find(offset)->second;
    }

    void activate() { active_ = true; }
    bool active() const { return active_; }

    std::shared_ptr<GammaBranch<VecType>>& branch(const int i) { return branches_[i]; }
    const std::shared_ptr<GammaBranch<VecType>>& branch(const int i) const { return branches_[i]; }
    std::shared_ptr<GammaBranch<VecType>>& branch(const GammaSQ sq) { return branch(static_cast<int>(sq)); }
    const std::shared_ptr<GammaBranch<VecType>>& branch(const GammaSQ sq) const { return branch(static_cast<int>(sq)); }

    const std::map<int, std::shared_ptr<const VecType>>& bras() const { return bras_; }
    std::map<int, std::shared_ptr<Matrix>>& gammas() { return gammas_; }
};

template <typename VecType>
class GammaTree {
  protected:
    std::shared_ptr<const VecType> ket_;
    std::shared_ptr<GammaBranch<VecType>> base_;

  public:
    GammaTree(std::shared_ptr<const VecType> ket) : ket_(ket) {
      base_ = std::make_shared<GammaBranch<VecType>>();
      const int nops = 4;

      for (int i = 0; i < nops; ++i) {
        base_->branch(i) = std::make_shared<GammaBranch<VecType>>();
        for (int j = 0; j < nops; ++j) {
          base_->branch(i)->branch(j) = std::make_shared<GammaBranch<VecType>>();
          for (int k = 0; k < nops; ++k) {
            base_->branch(i)->branch(j)->branch(k) = std::make_shared<GammaBranch<VecType>>();
          }
        }
      }
    }

    std::shared_ptr<GammaBranch<VecType>> base() { return base_; }

    template <class ...GSQs>
    void insert(std::shared_ptr<const VecType> bra, const int offset, const GSQs... ops) { base_->insert(bra, offset, ops...); }

    template <class ...GSQs>
    std::shared_ptr<const Matrix> search(const int offset, GSQs... address) const { return base_->search(offset, address...); }

    std::shared_ptr<const VecType> ket() const { return ket_; }

  private:
    static const int Alpha = 0;
    static const int Beta = 1;
    static const int Create = 0;
    static const int Annihilate = 1;

    std::unique_ptr<double[]> dot_product(std::shared_ptr<const VecType> bras, std::shared_ptr<const VecType> kets) const;
};

template <typename VecType>
class GammaTask; // Forward declaration of GammaTask

template <typename VecType, int N>
class GammaForest {
  protected:
    std::array<std::map<int, std::shared_ptr<GammaTree<VecType>>>, N> forests_;

  public:
    GammaForest() {}

    template <int unit, class ...operations>
    void insert(std::shared_ptr<const VecType> ket, const int ioffset, std::shared_ptr<const VecType> bra, const int joffset, const operations... ops) {
      std::shared_ptr<GammaTree<VecType>> gtree = tree<unit>(ket, ioffset);
      gtree->insert(bra, joffset, ops...);
    }

    template <int unit, class ...operations>
    std::shared_ptr<const Matrix> get(const int ioffset, const int joffset, const operations... ops) const {
      auto itree = forests_[unit].find(ioffset); assert(itree!=forests_[unit].end());
      return itree->second->search(joffset, ops...);
    }

    void compute() {
      const int nops = 4;

      int ntasks = 0;
      // Allocate memory while counting tasks
      for (auto& iforest : forests_) {
        for (auto& itreemap : iforest) {
          std::shared_ptr<GammaTree<VecType>> itree = itreemap.second;
          const int nA = itree->ket()->ij();
          const int norb = itree->ket()->det()->norb();

          // Allocation sweep
          for (int i = 0; i < nops; ++i) {
            std::shared_ptr<GammaBranch<VecType>> first = itree->base()->branch(i);
            if (!first->active()) continue;
            ++ntasks;
            for (auto& ibra : first->bras()) {
              const int nAp = ibra.second->ij();
              const int nstates = nA * nAp;
              first->gammas().emplace(ibra.first, std::make_shared<Matrix>(nstates, norb));
            }

            for (int j = 0; j < nops; ++j) {
              std::shared_ptr<GammaBranch<VecType>> second = first->branch(j);
              if (!second->active()) continue;
              for (auto& jbra : second->bras()) {
                const int nAp = jbra.second->ij();
                const int nstates = nA * nAp;
                second->gammas().emplace(jbra.first, std::make_shared<Matrix>(nstates, norb * norb));
              }

              for (int k = 0; k < nops; ++k) {
                std::shared_ptr<GammaBranch<VecType>> third = second->branch(k);
                if (!third->active()) continue;
                for (auto& kbra : third->bras()) {
                  const int nAp = kbra.second->ij();
                  const int nstates = nA * nAp;
                  third->gammas().emplace(kbra.first, std::make_shared<Matrix>(nstates, norb * norb * norb));
                }
              }
            }
          }
        }
      }

      TaskQueue<GammaTask<VecType>> tasks(ntasks);

      // Add tasks
      for (auto& iforest : forests_) {
        for (auto& itreemap : iforest) {
          std::shared_ptr<GammaTree<VecType>> itree = itreemap.second;

          const int norb = itree->ket()->det()->norb();
          for (int i = 0; i < nops; ++i) {
            std::shared_ptr<GammaBranch<VecType>> first = itree->base()->branch(i);
            if (!first->active()) continue;
            for (int a = 0; a < norb; ++a) tasks.emplace_back(itree, GammaSQ(i), a);
          }
        }
      }

      tasks.compute();
    }



  private:
    template <int unit>
    std::shared_ptr<GammaTree<VecType>> tree(std::shared_ptr<const VecType> ket, const int ioffset) {
      typename std::map<int, std::shared_ptr<GammaTree<VecType>>>::iterator itree = forests_[unit].find(ioffset);
      if (itree == forests_[unit].end()) {
        forests_[unit].emplace(ioffset, std::make_shared<GammaTree<VecType>>(ket));
        itree = forests_[unit].find(ioffset);
      }
      return itree->second;
    }
};

template <typename VecType>
class GammaTask {
  protected:
    const int a_;                            // Orbital
    const GammaSQ  operation_;               // Which operation
    const std::shared_ptr<GammaTree<VecType>> tree_;  // destination

  public:
    GammaTask(const std::shared_ptr<GammaTree<VecType>> tree, const GammaSQ operation, const int a) : a_(a), operation_(operation), tree_(tree) {}

    void compute() {
      const int nops = 4;
      const int norb = tree_->ket()->det()->norb();

      auto action = [] (const int op) { return (GammaSQ(op)==GammaSQ::CreateAlpha || GammaSQ(op)==GammaSQ::CreateBeta); };
      auto spin = [] (const int op) { return (GammaSQ(op)==GammaSQ::CreateAlpha || GammaSQ(op)==GammaSQ::AnnihilateAlpha); };

      std::shared_ptr<GammaBranch<VecType>> first = tree_->base()->branch(operation_);
      assert(first->active()); // This should have been checked before sending it to the TaskQueue

      std::shared_ptr<const VecType> avec = tree_->ket()->apply(a_, action(static_cast<int>(operation_)), spin(static_cast<int>(operation_)));
      for (auto& ibra : first->bras()) {
        const int nstates = avec->ij() * ibra.second->ij();
        std::unique_ptr<double[]> tmp = dot_product(ibra.second, avec);
        double* target = first->gammas().find(ibra.first)->second->element_ptr(0,a_);
        std::copy_n(tmp.get(), nstates, target);
      }

      for (int j = 0; j < nops; ++j) {
        auto second = first->branch(j);
        if (!second->active()) continue;

        for (int b = 0; b < norb; ++b) {
          std::shared_ptr<const VecType> bvec = avec->apply(b, action(j), spin(j));
          for (auto& jbra : second->bras()) {
            const int nstates = bvec->ij() * jbra.second->ij();
            std::unique_ptr<double[]> tmp = dot_product(jbra.second, bvec);
            double* target = second->gammas().find(jbra.first)->second->element_ptr(0, a_ + norb*b);
            std::copy_n(tmp.get(), nstates, target);
          }

          for (int k = 0; k < nops; ++k) {
            std::shared_ptr<GammaBranch<VecType>> third = second->branch(k);
            if (!third->active()) continue;

            for (int c = 0; c < norb; ++c) {
              std::shared_ptr<const VecType> cvec = bvec->apply(c, action(k), spin(k));
              for (auto& kbra : third->bras()) {
                const int nstates = cvec->ij() * kbra.second->ij();
                std::unique_ptr<double[]> tmp = dot_product(kbra.second, cvec);
                double* target = third->gammas().find(kbra.first)->second->element_ptr(0, a_ + norb * b + norb * norb * c);
                std::copy_n(tmp.get(), nstates, target);
              }
            }
          }
        }
      }
    }

    private:
      std::unique_ptr<double[]> dot_product(std::shared_ptr<const VecType> bras, std::shared_ptr<const VecType> kets) const {
        const int nbras = bras->ij();
        const int nkets = kets->ij();

        std::unique_ptr<double[]> out(new double[nbras*nkets]);
        double* odata = out.get();

        for (int iket = 0; iket < nkets; ++iket) {
          for (int jbra = 0; jbra < nbras; ++jbra, ++odata) {
            *odata = bras->data(jbra)->dot_product(*kets->data(iket));
          }
        }

        return out;
      }
};

}

#endif
