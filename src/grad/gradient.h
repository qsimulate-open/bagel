//
// Newint - Parallel electron correlation program.
// Filename: gradient.h
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


#ifndef __SRC_GRAD_GRADIENT_H
#define __SRC_GRAD_GRADIENT_H

// a class for using the BFGS solver, which requires
// "clone, ddot and daxpy, along with overloaded operators and a copy constructor"

#include <vector>
#include <memory>
#include <cassert>
#include <src/util/f77.h>

class Gradient {
  protected:
    std::vector<double> data_;

  public:
    Gradient(const size_t i, const double a = 0.0) : data_(i, a) {};
    Gradient(std::vector<double> a) : data_(a.begin(), a.end()) { assert(data_.size()%3 == 0); };
    Gradient(const Gradient& o) : data_(o.data_.begin(), o.data_.end()) {};
    ~Gradient() {};

    double ddot(const Gradient& o) const { return ddot_(size(), data(), 1, o.data(), 1); };
    double ddot(const std::shared_ptr<Gradient> o) const { return ddot(*o); };

    double* data() { return &data_[0]; };
    const double* data() const { return &data_[0]; };
    size_t size() const { return data_.size(); }; 

    void daxpy(const double a, const Gradient& o) { daxpy_(size(), a, o.data(), 1, data(), 1); }; 
    void daxpy(const double a, const std::shared_ptr<Gradient> o) { daxpy(a, *o); };

    Gradient operator-(const Gradient& o) {
      Gradient out(*this);
      out.daxpy(-1.0, o); 
      return out;
    };

    double data(int i, int j) { return data_.at(3*i+j); };

};

#endif
