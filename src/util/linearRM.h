//
// BAGEL - Parallel electron correlation program.
// Filename: linearRM.h
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


#ifndef __NEWINT_SRC_UTIL_LinearRM_H
#define __NEWINT_SRC_UTIL_LinearRM_H

#include <memory>
#include <list>
#include <stdexcept>
#include <src/util/f77.h>

// Solves a linear equation using residual minimization.
// Compared to Linear, this is robust because we don't assume
// positive definiteness of the A matrix here.

namespace bagel {

template<typename T>
class LinearRM {

  protected:
    std::list<std::shared_ptr<const T> > c_;
    std::list<std::shared_ptr<const T> > sigma_;

    const int max_;
    int size_;
    const std::shared_ptr<const T> grad_;

    // contains
    std::unique_ptr<double[]> mat_;
    std::unique_ptr<double[]> scr_;
    std::unique_ptr<double[]> vec_;
    std::unique_ptr<double[]> prod_;
    std::unique_ptr<int[]> ipiv_;
    int info;

    // for convenience below
    double& mat(int i, int j) { return mat_[i+j*max_]; };
    double& scr(int i, int j) { return scr_[i+j*max_]; };


  public:
    LinearRM(const int ndim, const std::shared_ptr<const T> grad) : max_(ndim), size_(0), grad_(grad),
      mat_(new double[max_*max_]),
      scr_(new double[max_*max_]),
      vec_(new double[max_]),
      prod_(new double[max_]),
      ipiv_(new int[max_*2]) {
    };
    ~LinearRM() {};

    std::shared_ptr<T> compute_residual(const std::shared_ptr<const T> c, const std::shared_ptr<const T> s) {

      if (size_ == max_) throw std::runtime_error("max size reached in Linear");
      // register new vectors
      c_.push_back(c);
      sigma_.push_back(s);
      // first set mat (=x(i)A^dag Ax(j)) and prod (= x(i)A^dag *y)
      ++size_;
      auto citer = sigma_.begin();
      for (int i = 0; i != size_; ++i, ++citer) {
        mat(i,size_-1) = mat(size_-1,i) = s->ddot(**citer);
      }
      // NOTE THE MINUS SIGN HERE!!
      prod_[size_-1] = - s->ddot(*grad_);

      // set to scr_
      std::copy(mat_.get(), mat_.get()+max_*max_, scr_.get());
      std::copy(prod_.get(), prod_.get()+max_, vec_.get());
      dgesv_(size_, 1, scr_, max_, ipiv_, vec_, size_, info);
      if (info) throw std::runtime_error("dgesv failed in Linear");

      std::shared_ptr<T> out(new T(*grad_));
      int cnt = 0;
      for (auto j = sigma_.begin(); j != sigma_.end(); ++j, ++cnt)
        out->daxpy(vec_[cnt], *j);
      assert(cnt == size_);
      return out;
    };

    std::shared_ptr<T> civec() const {
      std::shared_ptr<T> out = c_.front()->clone();
      int cnt = 0;
      for (auto i = c_.begin(); i != c_.end(); ++i, ++cnt) {
        out->daxpy(vec_[cnt], *i);
      }
      return out;
    };

    // make cc orthogonal to cc_ vectors
    double orthog(std::shared_ptr<T>& cc) const { return cc->orthog(c_); }

};

}

#endif
