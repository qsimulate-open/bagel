//
// Newint - Parallel electron correlation program.
// Filename: aughess.h
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


#ifndef __NEWINT_SRC_UTIL_AUGHESS_H
#define __NEWINT_SRC_UTIL_AUGHESS_H

#include <memory>
#include <list>
#include <stdexcept>
#include <src/util/f77.h>

template<typename T>
class AugHess {
  typedef std::shared_ptr<T> RefT;
  typedef std::list<std::pair<RefT, RefT> > Container_type_;
  typedef typename Container_type_::iterator iterator;

  protected:
    std::list<std::shared_ptr<T> > c_;
    std::list<std::shared_ptr<T> > sigma_;

    const int max_;    
    const std::shared_ptr<T> grad_;

    // contains 
    std::unique_ptr<double[]> mat_;
    std::unique_ptr<double[]> prod_;
    // scratch area for diagonalization
    std::unique_ptr<double[]> scr_;
    std::unique_ptr<double[]> vec_; 
    // an eigenvector 
    std::unique_ptr<double[]> eig_;
    // work area in a lapack routine
    std::unique_ptr<double[]> work_;
    int lwork_;
    int info;

    int size_;
    // for convenience below
    double& mat(int i, int j) { return mat_[i+j*max_]; };
    double& scr(int i, int j) { return scr_[i+j*max_]; };


  public:
    AugHess(const int ndim, std::shared_ptr<T> grad) : max_(ndim), size_(0), grad_(grad), 
      mat_(new double[ndim*ndim]),
      scr_(new double[ndim*ndim]),
      vec_(new double[ndim]),
      prod_(new double[ndim]),
      work_(new double[ndim*5]),
      eig_(new double[ndim]),
      lwork_(ndim*5) {
    };
    ~AugHess() {};

    std::shared_ptr<T> compute_residual(std::shared_ptr<T> c, std::shared_ptr<T> s) {
      if (size_+1 == max_) throw std::runtime_error("max size reached in AugHess");
      // register new vectors
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
      std::copy(mat_.get(), mat_.get()+max_*max_, scr_.get());
      for (int i = 0; i != size_; ++i) {
        scr(size_, i) = scr(i, size_) = prod_[i];
      }
      scr(size_, size_) = 0.0;
      dsyev_("V", "U", size_+1, scr_, max_, eig_, work_, lwork_, info); 
      if (info) throw std::runtime_error("dsyev failed in AugHess");

      // scale eigenfunction
      for (int i = 0; i != size_; ++i) vec_[i] = scr_[i] / scr_[size_];
      
      std::shared_ptr<T> out(new T(*grad_)); 
      int cnt = 0;
      for (auto i = c_.begin(), j = sigma_.begin(); i != c_.end(); ++i, ++j, ++cnt) {
        out->daxpy(vec_[cnt], *j);
        out->daxpy(-vec_[cnt]*eig_[0], *i);
      }
      assert(cnt == size_);
      return out;
    };

    double eig() const { return eig_[0]; };

    std::shared_ptr<T> civec() const {
      std::shared_ptr<T> out = c_.front()->clone();
      int cnt = 0;
      for (auto i = c_.begin(); i != c_.end(); ++i, ++cnt) {
        out->daxpy(vec_[cnt], *i); 
      } 
      return out;
    };

    // make cc orthogonal to cc_ vectors
    double orthog(std::shared_ptr<T> cc) { return cc->orthog(c_); }

};

#endif
