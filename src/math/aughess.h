//
// BAGEL - Parallel electron correlation program.
// Filename: aughess.h
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


#ifndef __BAGEL_SRC_UTIL_AUGHESS_H
#define __BAGEL_SRC_UTIL_AUGHESS_H

#include <memory>
#include <list>
#include <stdexcept>
#include <src/util/f77.h>

namespace bagel {

template<typename T>
class AugHess {

  protected:
    std::list<std::shared_ptr<const T>> c_;
    std::list<std::shared_ptr<const T>> sigma_;

    const int max_;
    int size_;
    const std::shared_ptr<const T> grad_;

    // contains
    std::shared_ptr<Matrix> mat_;
    std::unique_ptr<double[]> prod_;
    // scratch area for diagonalization
    std::shared_ptr<Matrix> scr_;
    std::unique_ptr<double[]> vec_;
    // an eigenvector
    std::unique_ptr<double[]> eig_;

    // for convenience below
    double& mat(int i, int j) { return mat_->element(i,j); }
    double& scr(int i, int j) { return scr_->element(i,j); }


  public:
    AugHess(const int ndim, const std::shared_ptr<const T> grad) : max_(ndim), size_(0), grad_(grad),
      mat_(std::make_shared<Matrix>(ndim,ndim,true)),
      prod_(new double[ndim]),
      scr_(std::make_shared<Matrix>(ndim,ndim,true)),
      vec_(new double[ndim]),
      eig_(new double[ndim]) {
    }

    std::shared_ptr<T> compute_residual(const std::shared_ptr<const T> c, const std::shared_ptr<const T> s) {

      if (size_+1 == max_) throw std::runtime_error("max size reached in AugHess");
      // register new vectors
      assert(std::fabs(c->norm()-1.0) < 1.0e-8);
      c_.push_back(c);
      sigma_.push_back(s);
      // first set mat (=x(i)Ax(j)) and prod (= x(i)*y)
      ++size_;
      auto citer = c_.begin();
      for (int i = 0; i != size_; ++i, ++citer) {
        mat(i,size_-1) = mat(size_-1,i) = s->ddot(**citer);
      }
      prod_[size_-1] = c->ddot(*grad_);

      // set to scr_
      *scr_ = *mat_;
      // adding (1,0) vector as an additional basis function
      for (int i = 0; i != size_; ++i) {
        scr(size_, i) = scr(i, size_) = prod_[i];
      }
      scr(size_, size_) = 0.0;
      scr_->diagonalize(eig_.get());

      // scale eigenfunction
      for (int i = 0; i != size_; ++i)
        vec_[i] = scr_->element(i,0) / scr_->element(size_,0);

      auto out = std::make_shared<T>(*grad_);
      int cnt = 0;
      for (auto i = c_.begin(), j = sigma_.begin(); i != c_.end(); ++i, ++j, ++cnt) {
        out->daxpy(vec_[cnt], *j);
        out->daxpy(-vec_[cnt]*eig_[0], *i);
      }
      assert(cnt == size_);
      return out;
    }

    double eig() const { return eig_[0]; }

    std::shared_ptr<T> civec() const {
      std::shared_ptr<T> out = c_.front()->clone();
      int cnt = 0;
      for (auto i = c_.begin(); i != c_.end(); ++i, ++cnt) {
        out->daxpy(vec_[cnt], *i);
      }
      return out;
    }

    // make cc orthogonal to cc_ vectors
    double orthog(std::shared_ptr<T>& cc) { return cc->orthog(c_); }

};

}

#endif
