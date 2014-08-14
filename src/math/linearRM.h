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


#ifndef __BAGEL_SRC_MATH_LINEARRM_H
#define __BAGEL_SRC_MATH_LINEARRM_H

#include <src/math/matrix.h>
#include <src/math/matop.h>

// Solves a linear equation using residual minimization.
// Compared to Linear, this is robust because we don't assume
// positive definiteness of the A matrix here.

namespace bagel {

template<typename T>
class LinearRM {

  protected:
    std::list<std::shared_ptr<const T>> c_;
    std::list<std::shared_ptr<const T>> sigma_;

    const int max_;
    int size_;
    const std::shared_ptr<const T> grad_;

    // contains
    std::shared_ptr<Matrix> mat_;
    std::shared_ptr<Matrix> vec_;
    std::shared_ptr<Matrix> prod_;

  public:
    LinearRM(const int ndim, const std::shared_ptr<const T> grad) : max_(ndim), size_(0), grad_(grad) {
      mat_ = std::make_shared<Matrix>(max_, max_);
      prod_ = std::make_shared<Matrix>(max_, 1);
    }

    std::shared_ptr<T> compute_residual(const std::shared_ptr<const T> c, const std::shared_ptr<const T> s) {

      if (size_ == max_) throw std::runtime_error("max size reached in Linear");
      // register new vectors
      c_.push_back(c);
      sigma_.push_back(s);
      // first set mat (=x(i)A^dag Ax(j)) and prod (= x(i)A^dag *y)
      ++size_;
      auto citer = sigma_.begin();
      for (int i = 0; i != size_; ++i) {
        mat_->element(i,size_-1) = mat_->element(size_-1,i) = s->dot_product(**citer++);
      }
      // NOTE THE MINUS SIGN HERE!!
      prod_->element(size_-1,0) = - s->dot_product(*grad_);

      // set to scr_
      vec_ = prod_->solve(mat_, size_);

      auto out = std::make_shared<T>(*grad_);
      int cnt = 0;
      for (auto& j : sigma_)
        out->ax_plus_y(vec_->element(cnt++, 0), j);
      assert(cnt == size_);
      return out;
    }

    std::shared_ptr<T> civec() const {
      std::shared_ptr<T> out = c_.front()->clone();
      int cnt = 0;
      for (auto& i : c_)
        out->ax_plus_y(vec_->element(cnt++, 0), i);
      return out;
    }

    // make cc orthogonal to cc_ vectors
    double orthog(std::shared_ptr<T>& cc) const { return cc->orthog(c_); }

};

}

#endif
