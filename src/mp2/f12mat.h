//
// Newint - Parallel electron correlation program.
// Filename: f12mat.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


#ifndef __SRC_MP2_F12MAT_H
#define __SRC_MP2_F12MAT_H

#include <memory>
#include <algorithm>
#include <src/util/f77.h>

class F12Mat {
  protected:
    std::unique_ptr<double[]> data_;
    const size_t nocc_;

  public:
    F12Mat(const size_t i) : data_(new double[i*i*i*i]), nocc_(i) { };
    F12Mat(const size_t i, std::unique_ptr<double[]> b) : data_(std::move(b)), nocc_(i) { };
    F12Mat(const F12Mat& o) : data_(new double[o.size()]), nocc_(o.nocc()) {
      std::copy(o.data(), o.data()+o.size(), data_.get());
    }
    ~F12Mat() {};

    F12Mat operator*(const double a) const {
      F12Mat f(*this);
      dscal_(size(), a, f.data(), 1); 
      return f;
    };
    F12Mat operator-(const F12Mat& o) const {
      F12Mat f(*this);
      daxpy_(size(), -1.0, o.data(), 1, f.data(), 1);
      return f;
    };
    F12Mat& operator-=(const F12Mat& o) {
      daxpy_(size(), -1.0, o.data(), 1, data(), 1);
      return *this;
    };
    F12Mat operator+(const F12Mat& o) const {
      F12Mat f(*this);
      daxpy_(size(), 1.0, o.data(), 1, f.data(), 1);
      return f;
    };
    F12Mat& operator+=(const F12Mat& o) {
      daxpy_(size(), 1.0, o.data(), 1, data(), 1);
      return *this;
    };

    double* data() { return data_.get(); };
    const double* data() const { return data_.get(); };
    size_t nocc() const { return nocc_; };
    size_t size() const { return nocc_*nocc_*nocc_*nocc_; };

    double& data(const size_t i) { return data_[i]; };
    const double& data(const size_t i) const { return data_[i]; };

    double& data(int i, int j, int k, int l) { return data_[i+nocc_*(j+nocc_*(k+nocc_*l))]; };
    const double& data(int i, int j, int k, int l) const { return data_[i+nocc_*(j+nocc_*(k+nocc_*l))]; };

    std::shared_ptr<F12Mat> clone() const {
      std::shared_ptr<F12Mat> out(new F12Mat(nocc_));
      return out;
    };
    std::shared_ptr<F12Mat> copy() const {
      std::shared_ptr<F12Mat> out(new F12Mat(*this));
      return out;
    };

};

#endif


