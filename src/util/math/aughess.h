//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: aughess.h
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


#ifndef __BAGEL_SRC_UTIL_AUGHESS_H
#define __BAGEL_SRC_UTIL_AUGHESS_H

#include <memory>
#include <list>
#include <stdexcept>
#include <src/util/f77.h>
#include <src/util/math/matrix.h>
#include <src/util/math/matop.h>

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
    VectorB prod_;
    VectorB vec_;
    // an eigenvector
    VectorB eig_;

    // for convenience below
    double& mat(int i, int j) { return mat_->element(i,j); }

    const double maxstepsize = 1.0;

    // carbon copy from ORZ ((c) Yanatech)
    std::tuple<double,double> compute_lambda(const Matrix& mat1, const Matrix& mat2) const {
      const int nlast = mat1.ndim();
      VectorB v(nlast);
      assert(nlast > 1);
      double lambda_test = 1.0, lambda_lasttest, stepsize_lasttest, stepsize;
      int iok = 0;
      for (int i = 0; i < 10; ++i) {
        Matrix scr = mat1 + mat2 * (1.0/lambda_test);
        scr.diagonalize(v);

        if (std::fabs(scr(nlast-1,0)) < 0.1)
          std::cout << "*** [aughess.cc] warning : you may have convergence issue due to very small : " << scr(nlast-1,0) << std::endl;

        blas::scale_n(1.0/scr(nlast-1,0), scr.data(), nlast-1);
        std::shared_ptr<T> x = c_.front()->clone();
        auto citer = c_.begin();
        for (int ii = 0; ii != nlast-1; ++ii)
          x->ax_plus_y(scr(ii,0), *citer++);

        stepsize = x->norm() / std::fabs(lambda_test);

        if (i == 0) {
          if (stepsize <= maxstepsize) break;
          lambda_lasttest = lambda_test;
          lambda_test = stepsize/maxstepsize;
        } else {
          if (std::fabs(stepsize-maxstepsize)/maxstepsize < 0.01) break;

          if (stepsize > maxstepsize) {
            lambda_lasttest = lambda_test;
            lambda_test *= stepsize/maxstepsize;
          } else {
            if (iok++ > 2) break;
            const double d1 = maxstepsize - stepsize;
            const double d2 = stepsize_lasttest - maxstepsize;
            if (d1 == 0.0 || d1 == -d2) break;
            const double lambda_lasttest_ = lambda_lasttest;
            lambda_lasttest = lambda_test;
            lambda_test = (d1/(d1+d2)) * lambda_lasttest_ + (d2/(d1+d2)) * lambda_test;
          }
        }
        if (lambda_test < 1.0) lambda_test = 1.0;
        stepsize_lasttest = stepsize;
      }
      return std::make_tuple(lambda_test, stepsize);
    }

  public:
    AugHess(const int ndim, const std::shared_ptr<const T> grad) : max_(ndim), size_(0), grad_(grad),
      mat_(std::make_shared<Matrix>(ndim+1,ndim+1)), prod_(ndim), vec_(ndim), eig_(ndim) {
    }

    std::tuple<std::shared_ptr<T>,double,double,double> compute_residual(const std::shared_ptr<const T> c, const std::shared_ptr<const T> s) {

      if (size_+1 == max_) throw std::runtime_error("max size reached in AugHess");
      // register new vectors
      assert(std::fabs(c->norm()-1.0) < 1.0e-8);
      c_.push_back(c);
      sigma_.push_back(s);
      // first set mat (=x(i)Ax(j)) and prod (= x(i)*y)
      ++size_;
      auto citer = c_.begin();
      auto siter = sigma_.begin();
      for (int i = 0; i != size_; ++i, ++citer, ++siter) {
        mat(i,size_-1) = mat(size_-1,i) = 0.5*(s->dot_product(**citer) + c->dot_product(**siter));
      }
      prod_(size_-1) = c->dot_product(*grad_);

      // set to scr1
      Matrix scr1(size_+1, size_+1);
      const Matrix scr2 = *mat_->get_submatrix(0, 0, size_+1, size_+1);
      // adding (1,0) vector as an additional basis function
      for (int i = 0; i != size_; ++i)
        scr1(size_, i) = scr1(i, size_) = prod_(i);

      double lambda, stepsize;
      std::tie(lambda, stepsize) = compute_lambda(scr1, scr2);
      Matrix scr = scr1 + scr2 * (1.0/lambda);
      scr.diagonalize(eig_);

      // find the best vector ((c) Yanatech)
      int ivec = -1;
      for (int i = 0; i != size_+1; ++i)
        if (std::fabs(scr.element(size_,i)) <= 1.1 && std::fabs(scr.element(size_,i)) > 0.1) {
          ivec = i;
          break;
        }
      if (ivec < 0)
       throw std::logic_error("logical error in AugHess");

      // scale eigenfunction
      for (int i = 0; i != size_; ++i)
        vec_(i) = scr.element(i,ivec) / (lambda*scr.element(size_,ivec));

      auto out = std::make_shared<T>(*grad_);
      int cnt = 0;
      for (auto i = c_.begin(), j = sigma_.begin(); i != c_.end(); ++i, ++j, ++cnt) {
        out->ax_plus_y(vec_(cnt), *j);
        out->ax_plus_y(-vec_(cnt)*lambda*eig_(ivec), *i);
      }
      assert(cnt == size_);
      return std::make_tuple(out, lambda, eig_(ivec), stepsize);
    }

    std::shared_ptr<T> civec() const {
      std::shared_ptr<T> out = c_.front()->clone();
      int cnt = 0;
      for (auto i = c_.begin(); i != c_.end(); ++i, ++cnt) {
        out->ax_plus_y(vec_(cnt), *i);
      }
      return out;
    }

    // make cc orthogonal to cc_ vectors
    double orthog(std::shared_ptr<T>& cc) { return cc->orthog(c_); }

};

}

#endif
