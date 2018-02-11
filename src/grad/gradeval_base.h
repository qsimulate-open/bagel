//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gradeval_base.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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

#ifndef __SRC_GRAD_GRADEVAL_BASE_H
#define __SRC_GRAD_GRADEVAL_BASE_H

#include <mutex>
#include <src/util/math/xyzfile.h>
#include <src/wfn/geometry.h>

namespace bagel {

// GradTask is a helper class of GradEval_base and declared as friend
class GradEval_base;

#define GRADTASK_INCLUDE
#include <src/grad/gradtask.hpp>
#undef  GRADTASK_INCLUDE

// base class for gradient evaluations
class GradEval_base {
  protected:
    std::shared_ptr<const Geometry> geom_;

    /// contract 1-electron gradient integrals with density matrix "d" and energy weighted density matrix (or equivalent) "w"
    template<typename TaskType>
    std::vector<std::shared_ptr<GradTask>> contract_grad1e(const std::shared_ptr<const Matrix> d, const std::shared_ptr<const Matrix> w);
    /// same as above, but one can specify density matrices to each integral kernel
    template<typename TaskType>
    std::vector<std::shared_ptr<GradTask>> contract_grad1e(const std::shared_ptr<const Matrix> n, const std::shared_ptr<const Matrix> k, const std::shared_ptr<const Matrix> o);
    /// contract 1-electron gradient integrals used for DKH
    std::vector<std::shared_ptr<GradTask>> contract_graddkh1e(std::array<std::shared_ptr<const Matrix>, 4>);
    /// contract small NAI gradient integrals with an array of densities
    std::vector<std::shared_ptr<GradTask>> contract_gradsmall1e(std::array<std::shared_ptr<const Matrix>,6>);
    /// contract finite-nucleus NAI gradient
    std::vector<std::shared_ptr<GradTask>> contract_grad1e_fnai(const std::shared_ptr<const Matrix> n);
    std::vector<std::shared_ptr<GradTask>> contract_grad1e_fnai(std::array<std::shared_ptr<const Matrix>,6> n,  const std::shared_ptr<const Geometry> geom = nullptr);

    /// contract 3-index 2-electron gradient integrals with density matrix "o".
    std::vector<std::shared_ptr<GradTask>> contract_grad2e(const std::shared_ptr<const DFDist> o,               const std::shared_ptr<const Geometry> geom = nullptr);
    std::vector<std::shared_ptr<GradTask>> contract_grad2e(const std::array<std::shared_ptr<const DFDist>,6> o, const std::shared_ptr<const Geometry> geom = nullptr);

    /// contract 2-index 2-electron gradient integrals with density matrix "o".
    std::vector<std::shared_ptr<GradTask>> contract_grad2e_2index(const std::shared_ptr<const Matrix> o,        const std::shared_ptr<const Geometry> geom = nullptr);

    // the results will be stored in grad_
    std::shared_ptr<GradFile> grad_;
    std::vector<std::mutex> mutex_;

  public:
    GradEval_base(const std::shared_ptr<const Geometry> g) : geom_(g), grad_(std::make_shared<GradFile>(g->natom())), mutex_(g->natom()) { }

    /// compute gradient given density matrices
    std::shared_ptr<GradFile> contract_gradient(const std::shared_ptr<const Matrix> d, const std::shared_ptr<const Matrix> w,
                                                const std::shared_ptr<const DFDist> o, const std::shared_ptr<const Matrix> o2,
                                                const std::shared_ptr<const Matrix> v = nullptr, const bool numerical = false,
                                                const std::shared_ptr<const Geometry> g2 = nullptr,
                                                const std::shared_ptr<const DFDist> g2o = nullptr,
                                                const std::shared_ptr<const Matrix> g2o2 = nullptr);
    virtual std::shared_ptr<GradFile> compute() { assert(false); return nullptr; }

    friend class GradTask1;
    friend class GradTask1s;
    friend class GradTask2;
    friend class GradTask3;
    friend class GradTask1r;
    friend class GradTask1f;
    friend class GradTask3r;
    friend class GradTask1rf;
    friend class GradTask1d;
};

template<typename TBatch>
std::shared_ptr<GradFile> GradTask1::compute_os(std::shared_ptr<const Matrix> den) const {
  const int dimb1 = shell_[0]->nbasis();
  const int dimb0 = shell_[1]->nbasis();
  std::shared_ptr<const Matrix> cden = den->get_submatrix(offset_[1], offset_[0], dimb1, dimb0);
  TBatch batch(shell_);
  batch.compute();
  return batch.compute_gradient(cden, atomindex_[0], atomindex_[1], ge_->geom_->natom());
}

template<typename TBatch>
std::shared_ptr<GradFile> GradTask1s::compute_os(std::shared_ptr<const Matrix> den) const {
  const int dimb1 = shell_[0]->nbasis();
  const int dimb0 = shell_[1]->nbasis();
  std::shared_ptr<const Matrix> cden = den->get_submatrix(offset_[1], offset_[0], dimb1, dimb0);
  TBatch batch(shell_);
  batch.compute();
  return batch.compute_gradient(cden, atomindex_[0], atomindex_[1], ge_->geom_->natom());
}

template<typename TBatch>
std::shared_ptr<GradFile> GradTask1d::compute_os(std::shared_ptr<const Matrix> den) const {
  const int dimb1 = shell_[0]->nbasis();
  const int dimb0 = shell_[1]->nbasis();
  std::shared_ptr<const Matrix> cden = den->get_submatrix(offset_[1], offset_[0], dimb1, dimb0);
  TBatch batch(shell_);
  batch.compute();
  return batch.compute_gradient(cden, atomindex_[0], atomindex_[1], ge_->geom_->natom());
}


extern template
std::vector<std::shared_ptr<GradTask>> GradEval_base::contract_grad1e<GradTask1>(const std::shared_ptr<const Matrix> d, const std::shared_ptr<const Matrix> w);
extern template
std::vector<std::shared_ptr<GradTask>> GradEval_base::contract_grad1e<GradTask1s>(const std::shared_ptr<const Matrix> d, const std::shared_ptr<const Matrix> w);
extern template
std::vector<std::shared_ptr<GradTask>> GradEval_base::contract_grad1e<GradTask1>(const std::shared_ptr<const Matrix> n, const std::shared_ptr<const Matrix> k,
                                                                                 const std::shared_ptr<const Matrix> o);
extern template
std::vector<std::shared_ptr<GradTask>> GradEval_base::contract_grad1e<GradTask1s>(const std::shared_ptr<const Matrix> n, const std::shared_ptr<const Matrix> k,
                                                                                  const std::shared_ptr<const Matrix> o);

}

#endif
