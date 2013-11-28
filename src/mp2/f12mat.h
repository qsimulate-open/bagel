//
// BAGEL - Parallel electron correlation program.
// Filename: f12mat.h
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


#ifndef __SRC_MP2_F12MAT_H
#define __SRC_MP2_F12MAT_H

#include <stddef.h>
#include <iostream>
#include <iomanip>
#include <memory>
#include <algorithm>
#include <src/util/f77.h>
#include <stdexcept>
#include <src/math/matrix.h>

namespace bagel {

class F12Mat;

class F12Ten {
  protected:
    const size_t nocc_;
    const size_t dim0_;
    const size_t dim1_;
    std::unique_ptr<double[]> data_;

  public:
    F12Ten(const size_t i, const size_t d0, const size_t d1) : nocc_(i), dim0_(d0), dim1_(d1), data_(new double[i*i*d0*d1]) {}
    F12Ten(const size_t i, const size_t d0, const size_t d1, std::unique_ptr<double[]> b)
     : nocc_(i), dim0_(d0), dim1_(d1), data_(std::move(b)) {}
    F12Ten(const size_t i, const size_t d0, const size_t d1, std::shared_ptr<const Matrix> b)
     : nocc_(i), dim0_(d0), dim1_(d1), data_(new double[b->size()]) { std::copy_n(b->data(), b->size(), data_.get()); }
    F12Ten(const F12Ten& o) : nocc_(o.nocc_), dim0_(o.dim0_), dim1_(o.dim1_), data_(new double[o.size()]) {
      std::copy_n(o.data(), o.size(), data_.get());
    }

    F12Ten operator*(const double a) const {
      F12Ten f(*this);
      std::for_each(f.data(), f.data()+size(), [&a](double& b) { b *= a; });
      return f;
    }
    F12Ten operator-(const F12Ten& o) const {
      F12Ten f(*this);
      daxpy_(size(), -1.0, o.data(), 1, f.data(), 1);
      return f;
    }
    F12Ten& operator-=(const F12Ten& o) {
      daxpy_(size(), -1.0, o.data(), 1, data(), 1);
      return *this;
    }
    F12Ten operator+(const F12Ten& o) const {
      F12Ten f(*this);
      daxpy_(size(), 1.0, o.data(), 1, f.data(), 1);
      return f;
    }
    F12Ten& operator+=(const F12Ten& o) {
      daxpy_(size(), 1.0, o.data(), 1, data(), 1);
      return *this;
    }

    double* data() { return data_.get(); }
    const double* data() const { return data_.get(); }

    size_t nocc() const { return nocc_; }
    size_t dim0() const { return dim0_; }
    size_t dim1() const { return dim1_; }

    size_t size() const { return nocc_*nocc_*dim0_*dim1_; }

    std::shared_ptr<F12Mat> contract(const std::shared_ptr<const F12Ten> o) const;
};

class F12Mat : public F12Ten {
  protected:

  public:
    F12Mat(const size_t i) : F12Ten(i,i,i) { }
    F12Mat(const size_t i, std::unique_ptr<double[]> b) : F12Ten(i, i, i, std::move(b)) { }
    F12Mat(const size_t i, std::shared_ptr<const Matrix> b) : F12Ten(i, i, i, b) { }
    F12Mat(const F12Mat& o) : F12Ten(o) { }

    F12Mat operator*(const double a) const {
      F12Mat f(*this);
      std::for_each(f.data(), f.data()+size(), [&a](double& b) { b *= a; });
      return f;
    }
    F12Mat operator-(const F12Mat& o) const {
      F12Mat f(*this);
      daxpy_(size(), -1.0, o.data(), 1, f.data(), 1);
      return f;
    }
    F12Mat& operator-=(const F12Mat& o) {
      daxpy_(size(), -1.0, o.data(), 1, data(), 1);
      return *this;
    }
    F12Mat operator+(const F12Mat& o) const {
      F12Mat f(*this);
      daxpy_(size(), 1.0, o.data(), 1, f.data(), 1);
      return f;
    }
    F12Mat& operator+=(const F12Mat& o) {
      daxpy_(size(), 1.0, o.data(), 1, data(), 1);
      return *this;
    }

    double* data() { return data_.get(); }
    const double* data() const { return data_.get(); }

    double& data(const size_t i) { return data_[i]; }
    const double& data(const size_t i) const { return data_[i]; }

    double& data(int i, int j, int k, int l) { return data_[i+nocc_*(j+nocc_*(k+nocc_*l))]; }
    const double& data(int i, int j, int k, int l) const { return data_[i+nocc_*(j+nocc_*(k+nocc_*l))]; }

    std::shared_ptr<F12Mat> clone() const {
      return std::make_shared<F12Mat>(nocc_);
    }
    std::shared_ptr<F12Mat> copy() const {
      return std::make_shared<F12Mat>(*this);
    }

    #define OUTSIZE 4
    void print() const {
      std::cout << "++++++++" << std::endl;
      for (int i = 0; i != OUTSIZE; ++i) {
        for (int j = 0; j != OUTSIZE; ++j) {
          for (int k = 0; k != OUTSIZE; ++k) {
            for (int l = 0; l != OUTSIZE; ++l) {
              std::cout << std::fixed << std::setw(9) << std::setprecision(6) << data(l,k,j,i)  << " ";
            }
          }
          std::cout << std::endl;
        }
      }
    }

    void symmetrize(const bool braket = true);
};

}

#endif


