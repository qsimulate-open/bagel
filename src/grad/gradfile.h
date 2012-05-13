//
// Newint - Parallel electron correlation program.
// Filename: gradfile.h
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


#ifndef __SRC_GRAD_GRADFILE_H
#define __SRC_GRAD_GRADFILE_H

// a class for using the BFGS solver, which requires
// "clone, ddot and daxpy, along with overloaded operators and a copy constructor"

#include <vector>
#include <memory>
#include <cassert>
#include <src/util/f77.h>

class GradFile {
  protected:
    std::vector<double> data_;

  public:
    GradFile(const size_t i, const double a = 0.0) : data_(i, a) {};
    GradFile(std::vector<double> a) : data_(a.begin(), a.end()) { assert(data_.size()%3 == 0); };
    GradFile(const GradFile& o) : data_(o.data_.begin(), o.data_.end()) {};
    ~GradFile() {};

    double ddot(const GradFile& o) const { return ddot_(size(), data(), 1, o.data(), 1); };
    double ddot(const std::shared_ptr<const GradFile> o) const { return ddot(*o); };

    double* data() { return &data_[0]; };
    const double* data() const { return &data_[0]; };
    size_t size() const { return data_.size(); }; 

    void daxpy(const double a, const GradFile& o) { daxpy_(size(), a, o.data(), 1, data(), 1); }; 
    void daxpy(const double a, const std::shared_ptr<const GradFile> o) { daxpy(a, *o); };

    GradFile operator-(const GradFile& o) const {
      GradFile out(*this);
      out.daxpy(-1.0, o); 
      return out;
    };

    GradFile& operator-=(const GradFile& o) { daxpy(-1.0, o); return *this; }; 

    std::shared_ptr<GradFile> clone() const { return std::shared_ptr<GradFile>(new GradFile(data_.size())); };

    double& data(int i, int j) { return data_.at(3*i+j); };
    const double& data(int i, int j) const { return data_.at(3*i+j); };

    double& data(size_t i) { return data_[i]; };
    const double& data(size_t i) const { return data_[i]; };

    void scale(const double a) { dscal_(size(), a, data(), 1); };

    double norm() const { return std::sqrt(ddot(*this)); };

};

#endif
