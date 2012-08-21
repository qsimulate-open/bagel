//
// Newint - Parallel electron correlation program.
// Filename: gradeval_base.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_GRAD_GRADEVAL_BASE_H
#define __SRC_GRAD_GRADEVAL_BASE_H

#include <memory>
#include <tuple>
#include <src/grad/gradfile.h>
#include <src/scf/geometry.h>
#include <src/scf/matrix1e.h>
#include <boost/thread/mutex.hpp>



class GradEval_base;

class GradTask {
  protected:
    std::array<std::shared_ptr<const Shell>,4> shell_;
    std::array<std::shared_ptr<const Shell>,2> shell2_;
    std::vector<int> atomindex_;
    std::vector<int> offset_;
    std::shared_ptr<const DF_AO> den_;
    std::shared_ptr<const Matrix1e> den2_;
    std::shared_ptr<const Matrix1e> eden_;
    GradEval_base* ge_;
    int rank_;

  public:
    GradTask(const std::array<std::shared_ptr<const Shell>,4>& s, const std::vector<int>& a, const std::vector<int>& o, const std::shared_ptr<const DF_AO> d, GradEval_base* p)
      : shell_(s), atomindex_(a), offset_(o), den_(d), ge_(p), rank_(3) {};
    GradTask(const std::array<std::shared_ptr<const Shell>,4>& s, const std::vector<int>& a, const std::vector<int>& o, const std::shared_ptr<const Matrix1e> d, GradEval_base* p)
      : shell_(s), atomindex_(a), offset_(o), den2_(d), ge_(p), rank_(3) {};
    GradTask(const std::array<std::shared_ptr<const Shell>,2>& s, const std::vector<int>& a, const std::vector<int>& o, const std::shared_ptr<const Matrix1e> d, const std::shared_ptr<const Matrix1e> w, GradEval_base* p)
      : shell2_(s), atomindex_(a), offset_(o), den2_(d), eden_(w), ge_(p), rank_(2) {
    };
    ~GradTask() {};

    void compute();
};



// base class for gradient evaluations
class GradEval_base {
  friend class GradTask;
  protected:
    const std::shared_ptr<const Geometry> geom_;

    /// contract 1-electron gradient integrals with density matrix "d" and energy weighted density matrix (or equivalent) "w"
    std::vector<GradTask> contract_grad1e(const std::shared_ptr<const Matrix1e> d, const std::shared_ptr<const Matrix1e> w);

    /// contract 3-index 2-electron gradient integrals with density matrix "o".
    std::vector<GradTask> contract_grad2e(const std::shared_ptr<const DF_AO> o);

    /// contract 2-index 2-electron gradient integrals with density matrix "o".
    std::vector<GradTask> contract_grad2e_2index(const std::unique_ptr<double[]>& o);

    // the results will be stored in grad_ and grad2_ (seperate area for multi-threading);
    std::shared_ptr<GradFile> grad_;
    std::vector<boost::mutex> mutex_;

  public:
    GradEval_base(const std::shared_ptr<const Geometry> g) : geom_(g), grad_(new GradFile(g->natom())), mutex_(g->natom()) { };
    ~GradEval_base() {};

    /// compute gradient given density matrices
    std::shared_ptr<GradFile> contract_gradient(const std::shared_ptr<const Matrix1e> d, const std::shared_ptr<const Matrix1e> w,
                                                const std::shared_ptr<const DF_AO> o, const std::unique_ptr<double[]>& o2);
    virtual std::shared_ptr<GradFile> compute() { assert(false); return std::shared_ptr<GradFile>(); };

};

#endif
