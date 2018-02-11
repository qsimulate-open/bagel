//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: linearRM.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#ifndef __BAGEL_SRC_MATH_LINEARRM_H
#define __BAGEL_SRC_MATH_LINEARRM_H

#include <src/util/math/matrix.h>
#include <src/util/math/matop.h>

// Solves a linear equation using residual minimization.
// Compared to Linear, this is robust because we don't assume
// positive definiteness of the A matrix here.

namespace bagel {

template<typename T, typename MatType = Matrix>
class LinearRM {

  protected:
    std::list<std::shared_ptr<const T>> c_;
    std::list<std::shared_ptr<const T>> sigma_;

    const int max_;
    int size_;
    const std::shared_ptr<const T> grad_;

    // contains
    std::shared_ptr<MatType> mat_;
    std::shared_ptr<MatType> overlap_;
    std::shared_ptr<MatType> vec_;
    std::shared_ptr<MatType> prod_;

  public:
    LinearRM(const int max, const std::shared_ptr<const T> grad) : max_(max), size_(0), grad_(grad) {
      if (max_ < 3)
        throw std::runtime_error("LinearRM works only if max >= 3");
      mat_     = std::make_shared<MatType>(max_, max_);
      overlap_ = std::make_shared<MatType>(max_, max_);
      prod_ = std::make_shared<MatType>(max_, 1);
    }

    std::shared_ptr<T> compute_residual(const std::shared_ptr<const T> c, const std::shared_ptr<const T> s) {

      // drop the oldest vector (second element)
      if (size_ == max_) {
        c_.erase(++c_.begin());
        sigma_.erase(++sigma_.begin());

        MatType trans(size_, size_-1);
        trans(0,0) = 1.0;
        for (int i = 1; i != size_-1; ++i)
          trans(i+1,i) = 1.0;

        mat_->copy_block(0, 0, size_-1, size_-1, trans % *mat_->get_submatrix(0, 0, size_, size_) * trans);
        overlap_->copy_block(0, 0, size_-1, size_-1, trans % *overlap_->get_submatrix(0, 0, size_, size_) * trans);
        prod_->copy_block(0, 0, size_-1, 1, trans % *prod_->get_submatrix(0, 0, size_, 1));
        --size_;
      }

      // register new vectors
      c_.push_back(c);
      sigma_.push_back(s);
      // first set mat (=x(i)A^dag Ax(j)) and prod (= x(i)A^dag *y)
      ++size_;
      auto citer = c_.begin();
      auto siter = sigma_.begin();
      for (int i = 0; i != size_; ++i) {
        mat_->element(size_-1,i) = s->dot_product(**siter++);
        mat_->element(i,size_-1) = detail::conj(mat_->element(size_-1,i));
        overlap_->element(size_-1,i) = c->dot_product(**citer++);
        overlap_->element(i,size_-1) = detail::conj(overlap_->element(size_-1,i));
      }
      // NOTE THE MINUS SIGN HERE!!
      prod_->element(size_-1,0) = - s->dot_product(*grad_);

      // temp areas
      std::shared_ptr<MatType> mat = mat_->get_submatrix(0,0,size_,size_);
      std::shared_ptr<MatType> ov = overlap_->get_submatrix(0,0,size_,size_);
      std::shared_ptr<MatType> prod = prod_->get_submatrix(0,0,size_,1);

      // canonical orthogonalization
      std::shared_ptr<const MatType> ovlp_scr = ov->tildex();

      // diagonalize matrix to get
      mat = std::make_shared<MatType>(*ovlp_scr % *mat * *ovlp_scr);
      prod = std::make_shared<MatType>(*ovlp_scr % *prod);

      // set to scr_
      vec_ = prod->solve(mat, prod->ndim());
      vec_ = std::make_shared<MatType>(*ovlp_scr * *vec_);

      // overwrite the first vector with the optimal vector
      {
        auto overwrite = [this](std::list<std::shared_ptr<const T>>& list) {
          int cnt = 0;
          std::shared_ptr<T> cn = list.front()->clone();
          for (auto& i : list)
            cn->ax_plus_y(vec_->element(cnt++, 0), i);
          list.front() = cn;
        };
        overwrite(c_);
        overwrite(sigma_);
        // matrix elements should be updated as well.
        MatType trans(size_, size_);
        trans.unit();
        std::copy_n(vec_->element_ptr(0,0), size_, trans.element_ptr(0,0));
        mat_->copy_block(0, 0, size_, size_, trans % *mat_->get_submatrix(0, 0, size_, size_) * trans);
        overlap_->copy_block(0, 0, size_, size_, trans % *overlap_->get_submatrix(0, 0, size_, size_) * trans);
        prod_->copy_block(0, 0, size_, 1, trans % *prod_->get_submatrix(0, 0, size_, 1));
      }

      std::shared_ptr<T> out = grad_->copy();
      out->ax_plus_y(1.0, sigma_.front());
      return out;
    }

    std::shared_ptr<T> civec() const { return c_.front()->copy(); }

};

}

#endif
