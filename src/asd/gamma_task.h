//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gamma_task.h
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

#ifndef __ASD_GAMMA_TASK_H
#define __ASD_GAMMA_TASK_H

#include <src/asd/gamma_tree.h>
#include <src/ci/ras/apply_block.h>

namespace bagel {

template <typename VecType>
class GammaTask {
  protected:
    const int a_;                            // Orbital
    const GammaSQ  operation_;               // Which operation
    const std::shared_ptr<GammaTree<VecType>> tree_;  // destination

  public:
    GammaTask(const std::shared_ptr<GammaTree<VecType>> tree, const GammaSQ operation, const int a) : a_(a), operation_(operation), tree_(tree) {}

    void compute() {
      constexpr int nops = 4;
      const int norb = tree_->norb();

      auto action = [] (const int op) { return is_creation(GammaSQ(op)); };
      auto spin = [] (const int op) { return is_alpha(GammaSQ(op)); };

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
            dot_product(jbra.second, bvec, second->gammas().find(jbra.first)->second->element_ptr(0, a_*norb + b));

          for (int k = 0; k < nops; ++k) {
            std::shared_ptr<GammaBranch<VecType>> third = second->branch(k);
            if (!third->active()) continue;

            for (int c = 0; c < norb; ++c) {
              if (b==c && k==j) continue;
              std::shared_ptr<const VecType> cvec = bvec->apply(c, action(k), spin(k));
              for (auto& kbra : third->bras())
                dot_product(kbra.second, cvec, third->gammas().find(kbra.first)->second->element_ptr(0, a_*norb*norb + b*norb + c));
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

template <>
void GammaTask<CASDvec>::compute();


template <class Branch>
class RASTask {
  protected:
    const int holes_;
    const int particles_;

  public:
    RASTask(const int h, const int p) : holes_(h), particles_(p) {}
    virtual std::shared_ptr<const RASString> stringspace(const int, const int, const int, const int , const int, const int) = 0;

    std::shared_ptr<RASBlock<double>> next_block(std::shared_ptr<Branch> branch, std::shared_ptr<const RASBlock<double>> base_block,
                                                 const int& orbital, const bool& action, const bool& spin) {
      std::shared_ptr<const RASString> sa = base_block->stringsa();
      std::shared_ptr<const RASString> sb = base_block->stringsb();

      std::array<int, 3> ras{{sa->ras<0>().second, sa->ras<1>().second, sa->ras<2>().second}};

      const int ras_space = ( orbital >= ras[0] ) + ( orbital >= (ras[0]+ras[1]) );
      std::array<int, 6> info{{ras[0] - sa->nholes(), ras[0] - sb->nholes(), sa->nele2(), sb->nele2(), sa->nparticles(), sb->nparticles()}};

      const int mod = action ? +1 : -1;
      info[2*ras_space]   += spin ? mod : 0;
      info[2*ras_space+1] += spin ? 0 : mod;

      // make sure it is a valid result
      for (int i = 0; i < 6; ++i)
        if (info[i] < 0 || info[i] > ras[i/2]) return nullptr;

      // is it out of space?
      const int nholes = 2*ras[0] - (info[0] + info[1]);
      const int nparts = info[4] + info[5];
      if (nholes > holes_ || nparts > particles_) {
        // peek ahead along the branch to see whether it will come back into the space

        // set of operations that would bring state back
        std::set<int> needed;
        if (nholes == (holes_+1)) {
          if (info[0] == ras[0]) {
            // no alpha holes, only beta holes
            needed.insert(static_cast<int>(GammaSQ::CreateBeta));
          } else if (info[1] == ras[0]) {
            // no beta holes, only alpha holes
            needed.insert(static_cast<int>(GammaSQ::CreateAlpha));
          } else {
            // both types of holes present
            needed.insert(static_cast<int>(GammaSQ::CreateAlpha));
            needed.insert(static_cast<int>(GammaSQ::CreateBeta));
          }
        } else if (nparts == (particles_+1)) {
          if (info[4] == 0) {
            // no alpha particles, only beta particles
            needed.insert(static_cast<int>(GammaSQ::AnnihilateBeta));
          } else if (info[5] == 0) {
            // no beta particles, only alpha particles
            needed.insert(static_cast<int>(GammaSQ::AnnihilateAlpha));
          } else {
            // both present
            needed.insert(static_cast<int>(GammaSQ::AnnihilateBeta));
            needed.insert(static_cast<int>(GammaSQ::AnnihilateAlpha));
          }
        } else {
          // impossible to contribute
          return nullptr;
        }

        // search down branch for active branches with needed operations
        if (!branch->if_contributes(needed))
          return nullptr;
      }

      std::shared_ptr<const RASString> ta = spin ? stringspace(info[0], ras[0], info[2], ras[1], info[4], ras[2]) : sa;
      std::shared_ptr<const RASString> tb = spin ? sb : stringspace(info[1], ras[0], info[3], ras[1], info[5], ras[2]);

      auto out = std::make_shared<RASBlock_alloc<double>>(ta,tb);

      RAS::Apply_block apply_block(orbital, action, spin);
      apply_block(base_block, out, branch->is_distterminal());

      return out;
    }
};

template <>
class GammaTask<RASDvec> : public RASTask<GammaBranch<RASDvec>> {
  protected:
    const int a_;                                     // Orbital
    const GammaSQ  operation_;                        // Which operation
    const std::shared_ptr<GammaTree<RASDvec>> tree_;  // destination

    // to avoid rebuilding the stringspaces repeatedly
    std::map<std::tuple<int, int, int, int, int, int>, std::shared_ptr<const RASString>> stringspaces_;

  public:
    GammaTask(const std::shared_ptr<GammaTree<RASDvec>> tree, const GammaSQ operation, const int a)
              : RASTask<GammaBranch<RASDvec>>(tree->ket()->det()->max_holes(), tree->ket()->det()->max_particles()),
                a_(a), operation_(operation), tree_(tree) {}

    void compute() {
      constexpr int nops = 4;
      const int norb = tree_->norb();

      auto action = [] (const int op) { return is_creation(GammaSQ(op)); };
      auto spin = [] (const int op) { return is_alpha(GammaSQ(op)); };

      std::shared_ptr<GammaBranch<RASDvec>> first = tree_->base()->branch(operation_);
      assert(first->active()); // This should have been checked before sending it to the TaskQueue

      std::shared_ptr<const RASDeterminants> base_det = tree_->ket()->det();

      const int nkets = tree_->ket()->ij();
      for (int iket = 0; iket < nkets; ++iket) {
        std::shared_ptr<const RASCivec> ketvec = tree_->ket()->data(iket);
        for (auto& ketblock : ketvec->blocks()) {
          if (!ketblock) continue;
          std::shared_ptr<const RASBlock<double>> ablock
            = next_block(first, ketblock, a_, action(static_cast<int>(operation_)), spin(static_cast<int>(operation_)));
          if (!ablock) continue;

          for (auto& ibra : first->bras())
            dot_product(ibra.second, ablock, first->gammas().find(ibra.first)->second->element_ptr(iket*ibra.second->ij(), a_));

          for (int j = 0; j < nops; ++j) {
            auto second = first->branch(j);
            if (!second->active()) continue;

            for (int b = 0; b < norb; ++b) {
              if (b==a_ && j==static_cast<int>(operation_)) continue;
              std::shared_ptr<const RASBlock<double>> bblock = next_block(second, ablock, b, action(j), spin(j));
              if (!bblock) continue;

              for (auto& jbra : second->bras())
                dot_product(jbra.second, bblock, second->gammas().find(jbra.first)->second->element_ptr(iket*jbra.second->ij(), a_*norb + b));

              for (int k = 0; k < nops; ++k) {
                std::shared_ptr<GammaBranch<RASDvec>> third = second->branch(k);
                if (!third->active()) continue;

                for (int c = 0; c < norb; ++c) {
                  if (b==c && k==j) continue;
                  std::shared_ptr<const RASBlock<double>> cblock = next_block(third, bblock, c, action(k), spin(k));
                  if (!cblock) continue;
                  for (auto& kbra : third->bras())
                    dot_product(kbra.second, cblock, third->gammas().find(kbra.first)->second->element_ptr(iket*kbra.second->ij(), a_*norb*norb + b*norb + c));
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

        if (bras->det()->allowed(ketblock->stringsb(), ketblock->stringsa())) {
          for (int jbra = 0; jbra < nbras; ++jbra, ++target) {
            std::shared_ptr<const RASBlock<double>> brablock = bras->data(jbra)->block(ketblock->stringsb(), ketblock->stringsa());
            if (brablock)
              *target += blas::dot_product(brablock->data(), brablock->size(), ketblock->data());
          }
        }
      }

      std::shared_ptr<const RASString> stringspace(const int a, const int b, const int c, const int d, const int e, const int f) final {
        auto iter = stringspaces_.find(std::make_tuple(a,b,c,d,e,f));
        if (iter != stringspaces_.end()) {
          return iter->second;
        }
        else {
          stringspaces_.emplace(std::make_tuple(a,b,c,d,e,f), std::make_shared<RASString>(a,b,c,d,e,f));
          return stringspaces_[std::make_tuple(a,b,c,d,e,f)];
        }
      }
};

}

#endif
