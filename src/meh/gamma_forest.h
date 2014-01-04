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
#include <src/ras/civector.h>
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
      for (auto& ibra : first->bras())
        dot_product(ibra.second, avec, first->gammas().find(ibra.first)->second->element_ptr(0,a_));

      for (int j = 0; j < nops; ++j) {
        auto second = first->branch(j);
        if (!second->active()) continue;

        for (int b = 0; b < norb; ++b) {
          if (b==a_ && j==static_cast<int>(operation_)) continue;
          std::shared_ptr<const VecType> bvec = avec->apply(b, action(j), spin(j));
          for (auto& jbra : second->bras())
            dot_product(jbra.second, bvec, second->gammas().find(jbra.first)->second->element_ptr(0, a_ + norb*b));

          for (int k = 0; k < nops; ++k) {
            std::shared_ptr<GammaBranch<VecType>> third = second->branch(k);
            if (!third->active()) continue;

            for (int c = 0; c < norb; ++c) {
              if (b==c && k==j) continue;
              std::shared_ptr<const VecType> cvec = bvec->apply(c, action(k), spin(k));
              for (auto& kbra : third->bras())
                dot_product(kbra.second, cvec, third->gammas().find(kbra.first)->second->element_ptr(0, a_ + norb * b + norb * norb * c));
            }
          }
        }
      }
    }

    private:
      void dot_product(std::shared_ptr<const VecType> bras, std::shared_ptr<const VecType> kets, double* target) const {
        const int nbras = bras->ij();
        const int nkets = kets->ij();

        for (int iket = 0; iket < nkets; ++iket) {
          for (int jbra = 0; jbra < nbras; ++jbra, ++target) {
            *target = bras->data(jbra)->dot_product(*kets->data(iket));
          }
        }
      }
};


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

template <>
void GammaForest<DistDvec, 2>::compute();

template <>
class GammaTask<RASDvec> {
  protected:
    const int a_;                            // Orbital
    const GammaSQ  operation_;               // Which operation
    const std::shared_ptr<GammaTree<RASDvec>> tree_;  // destination

    // to avoid rebuilding the stringspaces repeatedly
    std::map<std::tuple<int, int, int, int, int, int>, std::shared_ptr<const StringSpace>> stringspaces_;

  public:
    GammaTask(const std::shared_ptr<GammaTree<RASDvec>> tree, const GammaSQ operation, const int a) : a_(a), operation_(operation), tree_(tree) {}

    void compute() {
      const int nops = 4;
      const int norb = tree_->ket()->det()->norb();

      auto action = [] (const int op) { return (GammaSQ(op)==GammaSQ::CreateAlpha || GammaSQ(op)==GammaSQ::CreateBeta); };
      auto spin = [] (const int op) { return (GammaSQ(op)==GammaSQ::CreateAlpha || GammaSQ(op)==GammaSQ::AnnihilateAlpha); };

      const bool base_action = action(static_cast<int>(operation_));
      const bool base_spin = spin(static_cast<int>(operation_));

      std::shared_ptr<GammaBranch<RASDvec>> first = tree_->base()->branch(operation_);
      assert(first->active()); // This should have been checked before sending it to the TaskQueue

      std::shared_ptr<const RASDeterminants> base_det = tree_->ket()->det();

      const int nkets = tree_->ket()->ij();
      for (int iket = 0; iket < nkets; ++iket) {
        std::shared_ptr<const RASCivec> ketvec = tree_->ket()->data(iket);
        for (auto& ketblock : ketvec->blocks()) {
          if (!ketblock) continue;
          std::shared_ptr<const RASBlock<double>> ablock
            = next_block(ketblock, a_, action(static_cast<int>(operation_)), spin(static_cast<int>(operation_)));
          if (!ablock) continue;

          for (auto& ibra : first->bras())
            dot_product(ibra.second, ablock, first->gammas().find(ibra.first)->second->element_ptr(iket*ibra.second->ij(), a_));

          for (int j = 0; j < nops; ++j) {
            auto second = first->branch(j);
            if (!second->active()) continue;

            for (int b = 0; b < norb; ++b) {
              if (b==a_ && j==static_cast<int>(operation_)) continue;
              std::shared_ptr<const RASBlock<double>> bblock = next_block(ablock, b, action(j), spin(j));
              if (!bblock) continue;

              for (auto& jbra : second->bras())
                dot_product(jbra.second, bblock, second->gammas().find(jbra.first)->second->element_ptr(iket*jbra.second->ij(), a_ + norb*b));

              for (int k = 0; k < nops; ++k) {
                std::shared_ptr<GammaBranch<RASDvec>> third = second->branch(k);
                if (!third->active()) continue;

                for (int c = 0; c < norb; ++c) {
                  if (b==c && k==j) continue;
                  std::shared_ptr<const RASBlock<double>> cblock = next_block(bblock, c, action(k), spin(k));
                  if (!cblock) continue;
                  for (auto& kbra : third->bras())
                    dot_product(kbra.second, cblock, third->gammas().find(kbra.first)->second->element_ptr(iket*kbra.second->ij(), a_+norb*b+norb*norb*c));
                }
              }
            }
          }
        }
      }
    }

    private:
      void dot_product(std::shared_ptr<const RASDvec> bras, std::shared_ptr<const RASBlock<double>> ketblock, double* target) const {
        const int nbras = bras->ij();

        if (bras->det()->allowed(ketblock->stringb(), ketblock->stringa())) {
          for (int jbra = 0; jbra < nbras; ++jbra, ++target) {
            std::shared_ptr<const RASBlock<double>> brablock = bras->data(jbra)->block(ketblock->stringb(), ketblock->stringa());
            if (brablock)
              *target += blas::dot_product(brablock->data(), brablock->size(), ketblock->data());
          }
        }
      }

      std::shared_ptr<const StringSpace> stringspace(const int a, const int b, const int c, const int d, const int e, const int f) {
        auto iter = stringspaces_.find(std::make_tuple(a,b,c,d,e,f));
        if (iter != stringspaces_.end()) {
          return iter->second;
        }
        else {
          stringspaces_.emplace(std::make_tuple(a,b,c,d,e,f), std::make_shared<StringSpace>(a,b,c,d,e,f));
          return stringspaces_[std::make_tuple(a,b,c,d,e,f)];
        }
      }

      std::shared_ptr<RASBlock<double>> next_block(std::shared_ptr<const RASBlock<double>> base_block, const int& orbital, const bool& action, const bool& spin) {
        std::shared_ptr<const StringSpace> sa = base_block->stringa();
        std::shared_ptr<const StringSpace> sb = base_block->stringb();

        std::array<int, 3> ras{{sa->ras<0>().second, sa->ras<1>().second, sa->ras<2>().second}};

        const int ras_space = ( orbital >= ras[0] ) + ( orbital >= (ras[0]+ras[1]) );
        std::array<int, 6> info{{ras[0] - sa->nholes(), ras[0] - sb->nholes(), sa->nele2(), sb->nele2(), sa->nparticles(), sb->nparticles()}};

        const int mod = action ? +1 : -1;
        info[2*ras_space]   += spin ? mod : 0;
        info[2*ras_space+1] += spin ? 0 : mod;

        // make sure it is a valid result
        for (int i = 0; i < 6; ++i)
          if (info[i] < 0 || info[i] > ras[i/2]) return std::shared_ptr<RASBlock<double>>();

        std::shared_ptr<const StringSpace> ta = spin ? stringspace(info[0], ras[0], info[2], ras[1], info[4], ras[2]) : sa;
        std::shared_ptr<const StringSpace> tb = spin ? sb : stringspace(info[1], ras[0], info[3], ras[1], info[5], ras[2]);

        auto out = std::make_shared<RASBlock<double>>(ta,tb);

        std::shared_ptr<RAS::apply_block_base<double>> apply_block;
        switch ( 2*static_cast<int>(action) + static_cast<int>(spin) ) {
          case 0:
            apply_block = std::make_shared<RAS::apply_block_impl<double, false, false>>(orbital); break;
          case 1:
            apply_block = std::make_shared<RAS::apply_block_impl<double, false, true>>(orbital);  break;
          case 2:
            apply_block = std::make_shared<RAS::apply_block_impl<double, true, false>>(orbital);  break;
          case 3:
            apply_block = std::make_shared<RAS::apply_block_impl<double, true, true>>(orbital);   break;
          default:
            assert(false);
        }
        (*apply_block)(base_block, out);

        return out;
      }
};

}

#endif
