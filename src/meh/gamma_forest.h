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
#include <memory>

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

    std::shared_ptr<const Dvec> ket() const { return ket_; }

  private:
    static const int Alpha = 0;
    static const int Beta = 1;
    static const int Create = 0;
    static const int Annihilate = 1;

    std::shared_ptr<const Dvec> apply(std::shared_ptr<const Dvec> kets, const GammaSQ operation, const int orbital) const;
    std::unique_ptr<double[]> ddot(std::shared_ptr<const Dvec> bras, std::shared_ptr<const Dvec> kets) const;
};

class GammaTask; // Forward declaration of GammaTask

template <int N> class GammaForest {
  protected:
    std::array<std::map<int, std::shared_ptr<GammaTree>>, N> forests_;

  public:
    GammaForest() {};

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

    void compute() {
      const int nops = 4;

      int ntasks = 0;
      // Allocate memory while counting tasks
      for (auto& iforest : forests_) {
        for (auto& itreemap : iforest) {
          std::shared_ptr<GammaTree> itree = itreemap.second;
          const int nA = itree->ket()->ij();
          const int norb = itree->ket()->det()->norb();

          // Allocation sweep
          for (int i = 0; i < nops; ++i) {
            auto first = itree->base()->branch(i);
            if (!first->active()) continue;
            ++ntasks;
            for (auto& ibra : first->bras()) {
              const int nAp = ibra.second->ij();
              const int nstates = nA * nAp;
              first->gammas().insert(std::make_pair(ibra.first, std::make_shared<Matrix>(nstates, norb)));
            }

            for (int j = 0; j < nops; ++j) {
              auto second = first->branch(j);
              if (!second->active()) continue;
              for (auto& jbra : second->bras()) {
                const int nAp = jbra.second->ij();
                const int nstates = nA * nAp;
                second->gammas().insert(std::make_pair(jbra.first, std::make_shared<Matrix>(nstates, norb * norb)));
              }

              for (int k = 0; k < nops; ++k) {
                auto third = second->branch(k);
                if (!third->active()) continue;
                for (auto& kbra : third->bras()) {
                  const int nAp = kbra.second->ij();
                  const int nstates = nA * nAp;
                  third->gammas().insert(std::make_pair(kbra.first, std::make_shared<Matrix>(nstates, norb * norb * norb)));
                }
              }
            }
          }
        }
      }

      std::vector<GammaTask> tasks;
      tasks.reserve(ntasks);

      // Add tasks
      for (auto& iforest : forests_) {
        for (auto& itreemap : iforest) {
          std::shared_ptr<GammaTree> itree = itreemap.second;

          const int norb = itree->ket()->det()->norb();
          for (int i = 0; i < nops; ++i) {
            auto first = itree->base()->branch(i);
            if (!first->active()) continue;
            for (int a = 0; a < norb; ++a) tasks.emplace_back(itree, GammaSQ(i), a);
          }
        }
      }

      TaskQueue<GammaTask> tq(tasks);
      tq.compute(resources__->max_num_threads());
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

class GammaTask {
  protected:
    const int a_;                            // Orbital
    const GammaSQ  operation_;               // Which operation
    const std::shared_ptr<GammaTree> tree_;  // destination

  public:
    GammaTask(const std::shared_ptr<GammaTree> tree, const GammaSQ operation, const int a) : a_(a), operation_(operation), tree_(tree) {}

    void compute() {
      const int nops = 4;
      const int norb = tree_->ket()->det()->norb();

      std::shared_ptr<GammaBranch> first = tree_->base()->branch(operation_);
      assert(first->active()); // This should have been checked before sending it to the TaskQueue

      std::shared_ptr<const Dvec> avec = apply(tree_->ket(), operation_, a_);
      for (auto& ibra : first->bras()) {
        const int nstates = avec->ij() * ibra.second->ij();
        std::unique_ptr<double[]> tmp = ddot(ibra.second, avec);
        double* target = first->gammas().find(ibra.first)->second->element_ptr(0,a_);
        std::copy_n(tmp.get(), nstates, target);
      }

      for (int j = 0; j < nops; ++j) {
        auto second = first->branch(j);
        if (!second->active()) continue;

        for (int b = 0; b < norb; ++b) {
          std::shared_ptr<const Dvec> bvec = apply(avec, GammaSQ(j), b);
          for (auto& jbra : second->bras()) {
            const int nstates = bvec->ij() * jbra.second->ij();
            std::unique_ptr<double[]> tmp = ddot(jbra.second, bvec);
            double* target = second->gammas().find(jbra.first)->second->element_ptr(0, a_ + norb*b);
            std::copy_n(tmp.get(), nstates, target);
          }

          for (int k = 0; k < nops; ++k) {
            auto third = second->branch(k);
            if (!third->active()) continue;

            for (int c = 0; c < norb; ++c) {
              std::shared_ptr<const Dvec> cvec = apply(bvec, GammaSQ(k), c);
              for (auto& kbra : third->bras()) {
                const int nstates = cvec->ij() * kbra.second->ij();
                std::unique_ptr<double[]> tmp = ddot(kbra.second, cvec);
                double* target = third->gammas().find(kbra.first)->second->element_ptr(0, a_ + norb * b + norb * norb * c);
                std::copy_n(tmp.get(), nstates, target);
              }
            }
          }
        }
      }
    }

    private:
      std::shared_ptr<const Dvec> apply(std::shared_ptr<const Dvec> kets, const GammaSQ operation, const int orbital) const {
        const int Alpha = 0;
        const int Beta = 1;
        const int Create = 0;
        const int Annihilate = 1;

        const int action = ( (operation==GammaSQ::CreateAlpha || operation==GammaSQ::CreateBeta) ? Create : Annihilate);
        const int spin = ( (operation==GammaSQ::CreateAlpha || operation==GammaSQ::AnnihilateAlpha) ? Alpha : Beta );

        std::shared_ptr<const Determinants> source_det = kets->det();
        const int norb = source_det->norb();
        const int nstates = kets->ij();

        const int source_lena = source_det->lena();
        const int source_lenb = source_det->lenb();

        if (spin == Alpha) {
          std::shared_ptr<const Determinants> target_det = ( action == Annihilate ? source_det->remalpha() : source_det->addalpha() );

          auto out = std::make_shared<Dvec>(target_det, nstates);

          const int target_lena = target_det->lena();
          const int target_lenb = target_det->lenb();

          for (int i = 0; i < nstates; ++i) {
            double* target_base = out->data(i)->data();
            const double* source_base = kets->data(i)->data();

            for (auto& iter : (action == Annihilate ? source_det->phidowna(orbital) : source_det->phiupa(orbital)) ) {
              const double sign = static_cast<double>(iter.sign);
              double* target = target_base + target_lenb * iter.target;
              const double* source = source_base + source_lenb * iter.source;
              daxpy_(target_lenb, sign, source, 1, target, 1);
            }
          }

          return out;
        }
        else {
          std::shared_ptr<const Determinants> target_det = (action == Annihilate ? source_det->rembeta() : source_det->addbeta());

          const int target_lena = target_det->lena();
          const int target_lenb = target_det->lenb();

          auto out = std::make_shared<Dvec>(target_det, nstates);

          for (int i = 0; i < target_lena; ++i) {
            for (int istate = 0; istate < nstates; ++istate) {
              double* target_base = out->data(istate)->element_ptr(0,i);
              const double* source_base = kets->data(istate)->element_ptr(0,i);
              for (auto& iter : (action == Annihilate ? source_det->phidownb(orbital) : source_det->phiupb(orbital))) {
                const double sign = static_cast<double>(iter.sign);
                target_base[iter.target] += sign * source_base[iter.source];
              }
            }
          }

          return out;
        }
      }

      std::unique_ptr<double[]> ddot(std::shared_ptr<const Dvec> bras, std::shared_ptr<const Dvec> kets) const {
        const int nbras = bras->ij();
        const int nkets = kets->ij();

        std::unique_ptr<double[]> out(new double[nbras*nkets]);
        double* odata = out.get();

        for (int iket = 0; iket < nkets; ++iket) {
          for (int jbra = 0; jbra < nbras; ++jbra, ++odata) {
            *odata = bras->data(jbra)->ddot(*kets->data(iket));
          }
        }

        return out;
      }
};

}

#endif
