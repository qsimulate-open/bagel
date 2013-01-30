//
// BAGEL - Parallel electron correlation program.
// Filename: gradfile.h
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


#ifndef __SRC_GRAD_GRADFILE_H
#define __SRC_GRAD_GRADFILE_H

// a class for using the BFGS solver, which requires
// "clone, ddot and daxpy, along with overloaded operators and a copy constructor"

#include <vector>
#include <memory>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <src/util/matrix.h>
#include <src/wfn/geometry.h>
#include <src/util/f77.h>

namespace bagel {

class GradFile {
  protected:
    // matrix of 3*natom
    std::shared_ptr<Matrix> data_;

  public:
    GradFile(const size_t natom, const double a = 0.0) : data_(new Matrix(3, natom)) { data_->fill(a); }
    GradFile(const std::shared_ptr<const Matrix> a) : data_(new Matrix(*a)) { }
    GradFile(const GradFile& o) : data_(new Matrix(*o.data_)) {}
    ~GradFile() {}

    double ddot(const GradFile& o) const { return data_->ddot(o.data_); }
    double ddot(const std::shared_ptr<const GradFile> o) const { return ddot(*o); }

    std::shared_ptr<Matrix> data() { return data_; }
    std::shared_ptr<const Matrix> data() const { return data_; }
    size_t size() const { return data_->size(); }

    double& data(const size_t i) { return *(data_->data()+i); }
    const double& data(const size_t i) const { return *(data_->data()+i); }

    void daxpy(const double a, const GradFile& o) { data_->daxpy(a, o.data_); }
    void daxpy(const double a, const std::shared_ptr<const GradFile> o) { daxpy(a, *o); }

    GradFile operator+(const GradFile& o) const;
    GradFile operator-(const GradFile& o) const;
    GradFile& operator+=(const GradFile& o) { daxpy( 1.0, o); return *this; }
    GradFile& operator-=(const GradFile& o) { daxpy(-1.0, o); return *this; }
    GradFile& operator/=(const GradFile& o) { *data_ /= *o.data_; return *this; }
    GradFile operator/(const GradFile& o) const { GradFile out(*this); out /= o; return out; }

    std::shared_ptr<GradFile> clone() const;

    const std::shared_ptr<const Matrix> xyz() const { return data_; }
    double& data(int i, int j) { return data_->element(i, j); }
    const double& data(int i, int j) const { return data_->element(i, j); }

    void scale(const double a) { data_->scale(a); }
    double norm() const { return std::sqrt(ddot(*this)); }
    void print() const;

    // this function assumes that double[] has data_.size()*data_size() elements.
    std::shared_ptr<GradFile> transform(const std::unique_ptr<double[]>&, const bool transpose) const;

};

}

#endif
