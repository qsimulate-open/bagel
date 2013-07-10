//
// BAGEL - Parallel electron correlation program.
// Filename: zdavidson.h
// Copyright (C) 2013 Michael Caldwell
//
// Author: Michael Caldwell <caldwell@u.northwestern.edu>
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

// T should have
//  - double ddot(const T&)
//  - void daxpy(double, const T&) // added to self
//  - copy constructor (values are not used, though).
//  - void orthog(list<shared_ptr<T>>)
//

#ifndef __BAGEL_UTIL_ZDAVIDSON
#define __BAGEL_UTIL_ZDAVIDSON

#include <list>
#include <memory>
#include <stdexcept>
#include <src/util/f77.h>
#include <src/math/zmatrix.h>

namespace bagel {

template <typename T>
class ZDavidsonDiag {
  protected:
    const int nstate_;
    const int max_;
    int size_;
    bool orthogonalize_; // has linear dependence been encountered yet?

    std::list<std::shared_ptr<const T>> c_;
    std::list<std::shared_ptr<const T>> sigma_;

    // contains
    std::shared_ptr<ZMatrix> mat_;
    // scratch area for diagonalization
    std::shared_ptr<ZMatrix> scr_;
    std::unique_ptr<double[]> vec_;
    // an eigenvector
    std::shared_ptr<ZMatrix> eig_;
    // overlap matrix
    std::shared_ptr<ZMatrix> overlap_;
    std::shared_ptr<ZMatrix> ovlp_scr_;

    // for convenience below
    std::complex<double>& mat(int i, int j) { return mat_->element(i,j); }

  public:
    ZDavidsonDiag(int n, int m) : nstate_(n), max_(m*n), size_(0), orthogonalize_(false), mat_(std::make_shared<ZMatrix>(max_,max_,true)),
                                 scr_(std::make_shared<ZMatrix>(max_,max_,true)), vec_(new double[max_]), overlap_(std::make_shared<ZMatrix>(max_,max_,true)),
                                 ovlp_scr_(std::make_shared<ZMatrix>(max_,max_,true)) {
    }
    ~ZDavidsonDiag(){}

    double compute(std::shared_ptr<const T> cc, std::shared_ptr<const T> cs) {
      assert(nstate_ == 1);
      return compute(std::vector<std::shared_ptr<const T>>{cc},
                     std::vector<std::shared_ptr<const T>>{cs}).front();
    }

    std::vector<double> compute(std::vector<std::shared_ptr<const T>> cc, std::vector<std::shared_ptr<const T>> cs) {
      if (size_ == max_) throw std::runtime_error("max size reached in Davidson");
      // add entry
      for (auto& it : cc) c_.push_back(it);
      for (auto& it : cs) sigma_.push_back(it);

      // adding new matrix elements
      auto icivec = cc.begin();
      for (auto isigma = cs.begin(); isigma != cs.end(); ++isigma, ++icivec) {
        ++size_;
        double overlap_row = 0.0;
        auto cciter = c_.begin();
        for (int i = 0; i != size_; ++i, ++cciter) {
          mat(i, size_-1) = (**cciter).zdotc(*isigma);
          mat(size_-1, i) = conj(mat(i,size_-1));

          overlap_->element(i, size_-1) = (**cciter).zdotc(*icivec);
          overlap_->element(size_-1, i) = conj(overlap_->element(i, size_-1));

          if (!orthogonalize_) {
            overlap_row += abs(overlap_->element(i, size_-1));
          }
        }
        if ( fabs(overlap_row - 1.0) > 1.0e-8 ) {
          orthogonalize_ = true;
        }
      }

      if (orthogonalize_) {
       std::shared_ptr<ZMatrix> tmp = overlap_->get_submatrix(0, 0, size_, size_);
       tmp->inverse_half();

        ovlp_scr_->copy_block(0, 0, tmp);
      }

      // diagonalize matrix to get
      *scr_ = orthogonalize_ ? (*ovlp_scr_) % (*mat_ * *ovlp_scr_) : *mat_;
      std::shared_ptr<ZMatrix> tmp = scr_->get_submatrix(0, 0, size_, size_);
      tmp->diagonalize(vec_.get());
      scr_->copy_block(0, 0, tmp);
      if ( orthogonalize_ ) *scr_ = *ovlp_scr_ * *scr_;
      eig_ = scr_->slice(0,nstate_);

      return std::vector<double>(vec_.get(), vec_.get()+nstate_);
    }

    // perhaps can be cleaner.
    std::vector<std::shared_ptr<T>> residual() {
      std::vector<std::shared_ptr<T>> out;
      for (int i = 0; i != nstate_; ++i) {
        auto tmp = std::make_shared<T>(*c_.front());
        tmp->zero(); // <- waste of time
        int k = 0;
        for (auto iter = c_.begin(); iter != c_.end(); ++iter, ++k) {
          tmp->zaxpy(-vec_[i]*eig_->element(k,i), **iter);
        }
        k = 0;
        for (auto iter = sigma_.begin(); iter != sigma_.end(); ++iter, ++k) {
          tmp->zaxpy(eig_->element(k,i), **iter);
        }
        out.push_back(tmp);
      }
      return out;
    }

    // returns ci vector
    std::vector<std::shared_ptr<T>> civec() {
      std::vector<std::shared_ptr<T>> out;
      for (int i = 0; i != nstate_; ++i) {
        auto tmp = std::make_shared<T>(*c_.front());
        tmp->zero(); // <- waste of time
        int k = 0;
        for (auto iter = c_.begin(); iter != c_.end(); ++iter, ++k) {
          tmp->zaxpy((*eig_)[i*max_+k], **iter);
        }
        out.push_back(tmp);
      }
      return out;
    }

    // make cc orthogonal to cc_ vectors
    std::complex<double> orthog(std::shared_ptr<T>& cc) { return cc->orthog(c_); }

};

}

#endif
