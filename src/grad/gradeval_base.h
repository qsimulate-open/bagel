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

#ifndef __SRC_GRAD_GRADEVAL_BASE_H
#define __SRC_GRAD_GRADEVAL_BASE_H

#include <mutex>
#include <src/math/xyzfile.h>
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
    std::vector<std::shared_ptr<GradTask>> contract_grad1e(const std::shared_ptr<const Matrix> d, const std::shared_ptr<const Matrix> w);
    /// same as above, but one can specify density matrices to each integral kernel
    std::vector<std::shared_ptr<GradTask>> contract_grad1e(const std::shared_ptr<const Matrix> n, const std::shared_ptr<const Matrix> k, const std::shared_ptr<const Matrix> o);
    /// contract small NAI gradient integrals with an array of densities
    std::vector<std::shared_ptr<GradTask>> contract_gradsmall1e(std::array<std::shared_ptr<const Matrix>,6>);
    /// contract finite-nucleus NAI gradient
    std::vector<std::shared_ptr<GradTask>> contract_grad1e_fnai(const std::shared_ptr<const Matrix> n);
    std::vector<std::shared_ptr<GradTask>> contract_grad1e_fnai(std::array<std::shared_ptr<const Matrix>,6> n, const std::shared_ptr<const Geometry> geom = std::shared_ptr<const Geometry>());

    /// contract 3-index 2-electron gradient integrals with density matrix "o".
    std::vector<std::shared_ptr<GradTask>> contract_grad2e(const std::shared_ptr<const DFDist> o,
                                                           const std::shared_ptr<const Geometry> geom = std::shared_ptr<const Geometry>());
    std::vector<std::shared_ptr<GradTask>> contract_grad2e(const std::array<std::shared_ptr<const DFDist>,6> o,
                                                           const std::shared_ptr<const Geometry> geom = std::shared_ptr<const Geometry>());

    /// contract 2-index 2-electron gradient integrals with density matrix "o".
    std::vector<std::shared_ptr<GradTask>> contract_grad2e_2index(const std::shared_ptr<const Matrix> o,
                                                                  const std::shared_ptr<const Geometry> geom = std::shared_ptr<const Geometry>());

    // the results will be stored in grad_
    std::shared_ptr<GradFile> grad_;
    std::vector<std::mutex> mutex_;

  public:
    GradEval_base(const std::shared_ptr<const Geometry> g) : geom_(g), grad_(std::make_shared<GradFile>(g->natom())), mutex_(g->natom()) { }

    /// compute gradient given density matrices
    std::shared_ptr<GradFile> contract_gradient(const std::shared_ptr<const Matrix> d, const std::shared_ptr<const Matrix> w,
                                                const std::shared_ptr<const DFDist> o, const std::shared_ptr<const Matrix> o2,
                                                const std::shared_ptr<const Geometry> g2 = std::shared_ptr<const Geometry>(),
                                                const std::shared_ptr<const DFDist> g2o = std::shared_ptr<const DFDist>(),
                                                const std::shared_ptr<const Matrix> g2o2 = std::shared_ptr<const Matrix>());
    virtual std::shared_ptr<GradFile> compute() { assert(false); return std::shared_ptr<GradFile>(); }

  friend class GradTask1;
  friend class GradTask2;
  friend class GradTask3;
  friend class GradTask1r;
  friend class GradTask1f;
  friend class GradTask3r;
  friend class GradTask1rf;
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



}

#endif
