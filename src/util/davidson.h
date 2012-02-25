//
// Newint - Parallel electron correlation program.
// Filename: davidson.h
// Copyright (C) 2011 Toru Shiozaki
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

// T should have
//  - double ddot(const T&)
//  - void daxpy(double, const T&) // added to self
//  - copy constructor (values are not used, though).
//  - void orthog(list<shared_ptr<T> >) 
//

#ifndef __NEWINT_UTIL_DAVIDSON
#define __NEWINT_UTIL_DAVIDSON

#include <list>
#include <memory>
#include <stdexcept>
#include <src/util/f77.h>

template <typename T>
class DavidsonDiag {
  protected:
    std::list<std::shared_ptr<T> > c_;
    std::list<std::shared_ptr<T> > sigma_;
    const int max_;
    const int nstate_;
    int size_;

    // contains 
    std::vector<double> mat_;
    // scratch area for diagonalization
    std::vector<double> scr_;
    std::vector<double> vec_; 
    // an eigenvector 
    std::vector<double> eig_;
    // work area in a lapack routine
    std::vector<double> work_;
    int lwork_;

    // for convenience below
    double& mat(int i, int j) { return mat_[i+j*max_]; };

  public:
    DavidsonDiag(int n, int m) : nstate_(n), max_(m), size_(0) {
      mat_.resize(max_*max_);
      scr_.resize(max_*max_);
      vec_.resize(max_);
      work_.resize(max_*5);
      eig_.resize(max_*nstate_);
      lwork_ = max_*5;
    };
    ~DavidsonDiag(){};

    double compute(std::shared_ptr<T> cc, std::shared_ptr<T> cs) {
      assert(nstate_ == 1);
      std::vector<std::shared_ptr<T> > c; c.push_back(cc);
      std::vector<std::shared_ptr<T> > s; s.push_back(cs);
      return compute(c,s).front();
    };

    std::vector<double> compute(std::vector<std::shared_ptr<T> > cc, std::vector<std::shared_ptr<T> > cs) {
      if (size_ == max_) throw std::runtime_error("max size reached in Davidson");
      // add entry
      for (auto iter = cc.begin(); iter != cc.end(); ++iter) c_.push_back(*iter);
      for (auto iter = cs.begin(); iter != cs.end(); ++iter) sigma_.push_back(*iter);
      // adding new matrix elements 
      for (auto siter = cs.begin(); siter != cs.end(); ++siter) {
        auto iter = c_.begin();
        ++size_;
        for (int i = 0; i != size_; ++iter, ++i) {
          mat(i,size_-1) = mat(size_-1,i) = (*iter)->ddot(**siter);
        }
      }
      // diagonalize matrix to get 
      std::copy(mat_.begin(), mat_.end(), scr_.begin());
      int info = 0;
      dsyev_("V", "U", &size_, &(scr_[0]), &max_, &(vec_[0]), &(work_[0]), &lwork_, &info); 
      if (info != 0) throw std::runtime_error("dsyev in davidson failed");
      // copy energies and eigen functions
      std::copy(scr_.begin(), scr_.begin()+nstate_*max_, eig_.begin());
      std::vector<double> out(nstate_);
      std::copy(vec_.begin(), vec_.begin()+nstate_, out.begin());
      return out;
    };

    // perhaps can be cleaner.
    std::vector<std::shared_ptr<T> > residual() {
      std::vector<std::shared_ptr<T> > out;
      for (int i = 0; i != nstate_; ++i) {
        std::shared_ptr<T> tmp(new T(*c_.front()));
        tmp->zero(); // <- waste of time 
        int k = 0;
        for (auto iter = c_.begin(); iter != c_.end(); ++iter, ++k) {
          tmp->daxpy(-vec_[i]*eig_[i*max_+k], **iter);
        }
        k = 0;
        for (auto iter = sigma_.begin(); iter != sigma_.end(); ++iter, ++k) {
          tmp->daxpy(eig_[i*max_+k], **iter);
        }
        out.push_back(tmp);
      }
      return out;
    };

    // returns ci vector
    std::vector<std::shared_ptr<T> > civec() {
      std::vector<std::shared_ptr<T> > out;
      for (int i = 0; i != nstate_; ++i) {
        std::shared_ptr<T> tmp(new T(*c_.front()));
        tmp->zero(); // <- waste of time 
        int k = 0;
        for (auto iter = c_.begin(); iter != c_.end(); ++iter, ++k) {
          tmp->daxpy(eig_[i*max_+k], **iter);
        }
        out.push_back(tmp);
      }
      return out;
    };

    // make cc orthogonal to cc_ vectors
    double orthog(std::shared_ptr<T> cc) { return cc->orthog(c_); }

}; 

#endif
