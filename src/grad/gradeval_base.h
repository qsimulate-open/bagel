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

class Grad1eFile {
  protected:
    const std::shared_ptr<const Geometry> geom_;
    std::vector<std::shared_ptr<Matrix1e> > data_;

  public:
    Grad1eFile(const std::shared_ptr<const Geometry> g) : geom_(g) {
      for (int i = 0; i != g->natom()*3; ++i)
        data_.push_back(std::shared_ptr<Matrix1e>(new Matrix1e(geom_)));
    };
    ~Grad1eFile() {};

    std::shared_ptr<Matrix1e> data(const int i) { return data_[i]; };
    std::shared_ptr<const Matrix1e> data(const int i) const { return data_[i]; };

    std::shared_ptr<Matrix1e>& operator[](const int i) { return data_[i]; };

};


class GradEval_base;

class GradTask {
  protected:
    std::vector<std::shared_ptr<const Shell> > shell_;
    std::vector<int> atomindex_;
    std::vector<int> offset_;
    std::shared_ptr<const DF_AO> den_;
    std::shared_ptr<const Matrix1e> den2_;
    GradEval_base* ge_;

  public:
    GradTask(const std::vector<std::shared_ptr<const Shell> >& s, const std::vector<int>& a, const std::vector<int>& o, const std::shared_ptr<const DF_AO> d, GradEval_base* p)
      : shell_(s), atomindex_(a), offset_(o), den_(d), ge_(p) {};
    GradTask(const std::vector<std::shared_ptr<const Shell> >& s, const std::vector<int>& a, const std::vector<int>& o, const std::shared_ptr<const Matrix1e> d, GradEval_base* p)
      : shell_(s), atomindex_(a), offset_(o), den2_(d), ge_(p) {};
    ~GradTask() {};

    void compute();
};



// base class for gradient evaluations
class GradEval_base {
  friend class GradTask;
  protected:
    const std::shared_ptr<const Geometry> geom_;

    void compute_grad1e_integrals(std::shared_ptr<Grad1eFile>, std::shared_ptr<Grad1eFile>) const;

    /// contract 1-electron gradient integrals with density matrices "d" and "w" (the latter to be contracted with overlap derivatives) 
    std::vector<double> contract_grad1e(const std::shared_ptr<const Matrix1e> d, const std::shared_ptr<const Matrix1e> w,
                                        const std::shared_ptr<const Grad1eFile> g1, const std::shared_ptr<const Grad1eFile> go) const;

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
