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
#include <src/scf/geometry.h>
#include <src/scf/matrix1e.h>

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

// base class for gradient evaluations
class GradEval_base {
  protected:
    const std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<Grad1eFile> grad1e_;
    std::shared_ptr<Grad1eFile> grado_;
    void compute_grad1e();

  public:
    GradEval_base(const std::shared_ptr<const Geometry> g) : geom_(g), grad1e_(new Grad1eFile(g)), grado_(new Grad1eFile(g)) {
      compute_grad1e();
    };
    ~GradEval_base() {};

    std::shared_ptr<Matrix1e> grad1e(const int i) { return grad1e_->data(i); };
    std::shared_ptr<Matrix1e> grado(const int i) { return grado_->data(i); };

    std::tuple<std::shared_ptr<Grad1eFile>, std::shared_ptr<Grad1eFile> > data() const { return std::make_tuple(grad1e_, grado_); };

};

#endif
