//
// BAGEL - Parallel electron correlation program.
// Filename: gradeval_base.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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


#ifndef __SRC_GRAD_GRADEVAL_BASE_H
#define __SRC_GRAD_GRADEVAL_BASE_H

#include <src/grad/gradfile.h>
#include <mutex>

namespace bagel {

class GradEval_base;

class GradTask {
  protected:
    std::array<std::shared_ptr<const Shell>,4> shell_;
    std::array<std::shared_ptr<const Shell>,2> shell2_;
    std::array<int,4> atomindex_;
    std::array<int,4> offset_;
    std::shared_ptr<const DFDist> den_;
    std::shared_ptr<const Matrix> den2_;
    std::shared_ptr<const Matrix> den3_;
    std::shared_ptr<const Matrix> eden_;
    GradEval_base* ge_;
    int rank_;

  public:
    GradTask(const std::array<std::shared_ptr<const Shell>,4>& s, const std::vector<int>& a, const std::vector<int>& o,
             const std::shared_ptr<const DFDist> d, GradEval_base* p)
      : shell_(s), den_(d), ge_(p), rank_(3) { common_init(a,o); }

    GradTask(const std::array<std::shared_ptr<const Shell>,4>& s, const std::vector<int>& a, const std::vector<int>& o,
             const std::shared_ptr<const Matrix> d, GradEval_base* p)
      : shell_(s), den2_(d), ge_(p), rank_(2) { common_init(a,o); }

    GradTask(const std::array<std::shared_ptr<const Shell>,2>& s, const std::vector<int>& a, const std::vector<int>& o,
             const std::shared_ptr<const Matrix> nmat, const std::shared_ptr<const Matrix> kmat, const std::shared_ptr<const Matrix> omat, GradEval_base* p)
      : shell2_(s), den2_(nmat), den3_(kmat), eden_(omat), ge_(p), rank_(1) { common_init(a,o); }

    void common_init(const std::vector<int>& a, const std::vector<int>& o) {
      assert(a.size() == o.size());
      int k = 0;
      for (auto i = a.begin(), j = o.begin(); i != a.end(); ++i, ++j, ++k) {
        atomindex_[k] = *i;
        offset_[k] = *j;
      }
    }

    std::shared_ptr<GradFile> compute_nai(std::shared_ptr<const Matrix> den) const;

    template<typename TBatch>
    std::shared_ptr<GradFile> compute_os(std::shared_ptr<const Matrix> den) const;

    void compute();
};



// base class for gradient evaluations
class GradEval_base {
  friend class GradTask;
  protected:
    const std::shared_ptr<const Geometry> geom_;

    /// contract 1-electron gradient integrals with density matrix "d" and energy weighted density matrix (or equivalent) "w"
    std::vector<GradTask> contract_grad1e(const std::shared_ptr<const Matrix> d, const std::shared_ptr<const Matrix> w);
    /// same as above, but one can specify density matrices to each integral kernel
    std::vector<GradTask> contract_grad1e(const std::shared_ptr<const Matrix> n, const std::shared_ptr<const Matrix> k, const std::shared_ptr<const Matrix> o);

    /// contract 3-index 2-electron gradient integrals with density matrix "o".
    std::vector<GradTask> contract_grad2e(const std::shared_ptr<const DFDist> o);

    /// contract 2-index 2-electron gradient integrals with density matrix "o".
    std::vector<GradTask> contract_grad2e_2index(const std::shared_ptr<const Matrix> o);

    // the results will be stored in grad_
    std::shared_ptr<GradFile> grad_;
    std::vector<std::mutex> mutex_;

  public:
    GradEval_base(const std::shared_ptr<const Geometry> g) : geom_(g), grad_(std::make_shared<GradFile>(g->natom())), mutex_(g->natom()) { };
    ~GradEval_base() {};

    /// compute gradient given density matrices
    std::shared_ptr<GradFile> contract_gradient(const std::shared_ptr<const Matrix> d, const std::shared_ptr<const Matrix> w,
                                                const std::shared_ptr<const DFDist> o, const std::shared_ptr<const Matrix> o2);
    virtual std::shared_ptr<GradFile> compute() { assert(false); return std::shared_ptr<GradFile>(); };

};

template<typename TBatch>
std::shared_ptr<GradFile> GradTask::compute_os(std::shared_ptr<const Matrix> den) const {
  const int iatom0 = atomindex_[0];
  const int iatom1 = atomindex_[1];
  const int nbasis = ge_->geom_->nbasis();
  auto grad_local = std::make_shared<GradFile>(ge_->geom_->natom());
  const int dimb1 = shell2_[0]->nbasis();
  const int dimb0 = shell2_[1]->nbasis();

  TBatch batch(shell2_);
  batch.compute();
  const double* data = batch.data();
  const size_t s = batch.size_block();
  for (int i = offset_[0], cnt = 0; i != dimb0 + offset_[0]; ++i) {
    for (int j = offset_[1]; j != dimb1 + offset_[1]; ++j, ++cnt) {
      int jatom0 = batch.swap01() ? iatom1 : iatom0;
      int jatom1 = batch.swap01() ? iatom0 : iatom1;
      for (int k = 0; k != 3; ++k) {
        grad_local->data(k, jatom1) += data[cnt+s*k    ] * den->data(i*nbasis+j);
        grad_local->data(k, jatom0) += data[cnt+s*(k+3)] * den->data(i*nbasis+j);
      }
    }
  }
  return grad_local;
}

}

#endif
