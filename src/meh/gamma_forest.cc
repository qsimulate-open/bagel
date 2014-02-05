//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_forest.cc
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

#include <src/meh/gamma_forest.h>

using namespace std;
using namespace bagel;

// TODO this needs to be improved
template <>
void GammaForest<DistDvec, 2>::compute() {
  using DistBranch = GammaBranch<DistDvec>;
  using DistTree = GammaTree<DistDvec>;

  constexpr int nops = 4;

  allocate_and_count();

  // Compute tasks
  for (auto& iforest : forests_) {
    for (auto& itreemap : iforest) {
      shared_ptr<DistTree> itree = itreemap.second;

      const int norb = itree->ket()->det()->norb();
      for (int i = 0; i < nops; ++i) {
        shared_ptr<DistBranch> first = itree->base()->branch(i);
        if (!first->active()) continue;
        for (int a = 0; a < norb; ++a) {
          GammaTask<DistDvec> task(itree, GammaSQ(i), a);
          task.compute();
        }
      }
    }
  }
}

class GammaDistRASTask : public RASTask<GammaBranch<DistRASDvec>> {
  protected:
    shared_ptr<GammaTree<DistRASDvec>> tree_;
    shared_ptr<RASBlock<double>> starting_block_;
    const int istate_;
    const int operation_;
    const int a_;
    mutex* mutex_;
    map<tuple<int,int,int,int,int,int>, shared_ptr<const RASString>>* stringspaces_;
    mutex* ssmutex_;

  public:
    GammaDistRASTask(shared_ptr<GammaTree<DistRASDvec>> t, shared_ptr<RASBlock<double>> b, const int ist, const GammaSQ op, const int a,
      mutex* m, map<tuple<int,int,int,int,int,int>, shared_ptr<const RASString>>* ss, mutex* ssmutex) :
        RASTask<GammaBranch<DistRASDvec>>(t->ket()->det()->max_holes(), t->ket()->det()->max_particles()),
        tree_(t), starting_block_(b), istate_(ist), operation_(static_cast<int>(op)), a_(a), mutex_(m), stringspaces_(ss), ssmutex_(ssmutex) {}

    void compute() {
      constexpr int nops = 4;
      const int norb = tree_->ket()->det()->norb();

      auto action = [] (const int op) { return (GammaSQ(op)==GammaSQ::CreateAlpha || GammaSQ(op)==GammaSQ::CreateBeta); };
      auto spin = [] (const int op) { return (GammaSQ(op)==GammaSQ::CreateAlpha || GammaSQ(op)==GammaSQ::AnnihilateAlpha); };

      const bool base_action = action(operation_);
      const bool base_spin = spin(operation_);

      shared_ptr<GammaBranch<DistRASDvec>> first = tree_->base()->branch(operation_);
      assert(first->active()); // This should have been checked before sending it to the TaskQueue

      shared_ptr<const RASBlock<double>> ablock = next_block(first, starting_block_, a_, base_action, base_spin);

      if (ablock) {
        for (auto& ibra : first->bras())
          dot_product(ibra.second, ablock, first->gammas()[ibra.first]->element_ptr(istate_*ibra.second->ij(), a_));

        for (int j = 0; j < nops; ++j) {
          auto second = first->branch(j);
          if (!second->active()) continue;

          for (int b = 0; b < norb; ++b) {
            if (b==a_ && j==operation_) continue;
            shared_ptr<const RASBlock<double>> bblock = next_block(second, ablock, b, action(j), spin(j));
            if (!bblock) continue;

            for (auto& jbra : second->bras())
              dot_product(jbra.second, bblock, second->gammas()[jbra.first]->element_ptr(istate_*jbra.second->ij(), a_ + norb*b));

            for (int k = 0; k < nops; ++k) {
              shared_ptr<GammaBranch<DistRASDvec>> third = second->branch(k);
              if (!third->active()) continue;

              for (int c = 0; c < norb; ++c) {
                if (b==c && k==j) continue;
                shared_ptr<const RASBlock<double>> cblock = next_block(third, bblock, c, action(k), spin(k));
                if (!cblock) continue;
                for (auto& kbra : third->bras())
                  dot_product(kbra.second, cblock, third->gammas()[kbra.first]->element_ptr(istate_*kbra.second->ij(), a_+norb*b+norb*norb*c));
              }
            }
          }
        }
      }
    }

  private:
    void dot_product(shared_ptr<const DistRASDvec> bras, shared_ptr<const RASBlock<double>> ketblock, double* target) const {
      const int nbras = bras->ij();

      if (bras->det()->allowed(ketblock->stringb(), ketblock->stringa())) {
        vector<double> values(nbras, 0.0);
        for (int jbra = 0; jbra < nbras; ++jbra) {
          shared_ptr<const DistCIBlock<double>> brablock = bras->data(jbra)->block(ketblock->stringb(), ketblock->stringa());
          if (brablock) {
            const double val = blas::dot_product(brablock->local(), brablock->size(), ketblock->data() + brablock->astart()*brablock->lenb());
            values[jbra] = val;
          }
        }

        lock_guard<mutex> lk(*mutex_);
        for (int j = 0; j < nbras; ++j)
          target[j] += values[j];
      }
    }

    shared_ptr<const RASString> stringspace(const int a, const int b, const int c, const int d, const int e, const int f) final {
      lock_guard<mutex> lk(*ssmutex_);
      auto iter = stringspaces_->find(make_tuple(a,b,c,d,e,f));
      if (iter != stringspaces_->end()) {
        return iter->second;
      }
      else {
        auto out = make_shared<RASString>(a,b,c,d,e,f);
        stringspaces_->emplace(make_tuple(a,b,c,d,e,f), out);
        return out;
      }
    }
};

template <>
void GammaForest<DistRASDvec, 2>::compute() {
  using DistBranch = GammaBranch<DistRASDvec>;
  using DistTree = GammaTree<DistRASDvec>;

  constexpr int nops = 4;

  allocate_and_count();

  map<tuple<int,int,int,int,int,int>, shared_ptr<const RASString>> ssmap;
  mutex ssmut;

  // Compute tasks
  for (auto& iforest : forests_) {
    for (auto& itreemap : iforest) {
      shared_ptr<DistTree> itree = itreemap.second;

      const int norb = itree->ket()->det()->norb();
      const int nstates = itree->ket()->ij();

      vector<shared_ptr<RASCivec>> localvecs;
      for(auto& i : itree->ket()->dvec()) localvecs.push_back(make_shared<RASCivec>(*i));

      TaskQueue<GammaDistRASTask> tasks;
      vector<mutex> mutexes(nstates*nops*norb);

      for(int istate = 0; istate < nstates; ++istate) {
        for (int i = 0; i < nops; ++i) {
          shared_ptr<DistBranch> first = itree->base()->branch(i);
          if (first->active()) {
            for (int a = 0; a < norb; ++a) {
              for (auto& block : localvecs[istate]->blocks()) {
                if (block) tasks.emplace_back(itree, block, istate, GammaSQ(i), a, &mutexes[istate + nstates*(a + i*norb)], &ssmap, &ssmut);
              }
            }
          }
        }
      }

      tasks.compute();
    }
  }
  for_each_branch([] (shared_ptr<DistBranch> b) {
    for (auto& g : b->gammas())
      g.second->allreduce();
  });
}
