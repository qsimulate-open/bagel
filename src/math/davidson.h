//
// BAGEL - Parallel electron correlation program.
// Filename: davidson.h
// Copyright (C) 2011 Toru Shiozaki
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

// T should have
//  - double dot_product(const T&)
//  - void ax_plus_y(double, const T&) // added to self
//  - copy constructor (values are not used, though).
//  - void orthog(list<shared_ptr<T>>)
//

#ifndef __BAGEL_UTIL_DAVIDSON
#define __BAGEL_UTIL_DAVIDSON

#include <list>
#include <src/math/algo.h>
#include <src/util/f77.h>

namespace bagel {

template <typename T, class MatType = Matrix>
class DavidsonDiag {
  protected:
    const int nstate_;
    const int max_;
    int size_;

    std::list<std::shared_ptr<const T>> c_;
    std::list<std::shared_ptr<const T>> sigma_;

    // contains
    std::shared_ptr<MatType> mat_;
    // scratch area for diagonalization
    std::shared_ptr<MatType> scr_;
    std::unique_ptr<double[]> vec_;
    // an eigenvector
    std::shared_ptr<MatType> eig_;
    // overlap matrix
    std::shared_ptr<MatType> overlap_;
    std::shared_ptr<MatType> ovlp_scr_;

  public:
    // Davidson with periodic collapse of the subspace
    DavidsonDiag(int n, int max) : nstate_(n), max_(max*n), size_(0), mat_(std::make_shared<MatType>(max_,max_,true)),
                                 vec_(new double[max_]), overlap_(std::make_shared<MatType>(max_,max_,true)) {
    }

    double compute(std::shared_ptr<const T> cc, std::shared_ptr<const T> cs) {
      assert(nstate_ == 1);
      return compute(std::vector<std::shared_ptr<const T>>{cc},
                     std::vector<std::shared_ptr<const T>>{cs}).front();
    }

    std::vector<double> compute(std::vector<std::shared_ptr<const T>> cc, std::vector<std::shared_ptr<const T>> cs) {
      if (size_ + cc.size() > max_)
        collapse_subspace();

      // add entry
      for (auto& it : cc) c_.push_back(it);
      for (auto& it : cs) sigma_.push_back(it);

      // adding new matrix elements
      auto icivec = cc.begin();
      for (auto isigma = cs.begin(); isigma != cs.end(); ++isigma, ++icivec) {
        ++size_;
        auto cciter = c_.begin();
        for (int i = 0; i != size_; ++i, ++cciter) {
          mat_->element(i, size_-1) = (*cciter)->dot_product(**isigma);
          mat_->element(size_-1, i) = detail::conj(mat_->element(i, size_-1));

          overlap_->element(i, size_-1) = (*cciter)->dot_product(**icivec);
          overlap_->element(size_-1, i) = detail::conj(overlap_->element(i, size_-1));
        }
      }

      if ( std::fabs(static_cast<double>(size_) - overlap_->dot_product(*overlap_)) > 1.0e-6 )
        ovlp_scr_ = overlap_->get_submatrix(0, 0, size_, size_)->tildex();

      // diagonalize matrix to get
      scr_ = mat_->get_submatrix(0, 0, size_, size_);
      if (ovlp_scr_) scr_ = std::make_shared<MatType>(*ovlp_scr_ % *scr_ * *ovlp_scr_);
      scr_->diagonalize(vec_.get());
      if (ovlp_scr_) scr_ = std::make_shared<MatType>(*ovlp_scr_ * *scr_);

      // orthogonalize ci vectors
      if (ovlp_scr_)
        orthogonalize_subspace();

      eig_ = scr_->slice(0,nstate_);

      return std::vector<double>(vec_.get(), vec_.get()+nstate_);
    }

    // perhaps can be cleaner.
    std::vector<std::shared_ptr<T>> residual() {
      std::vector<std::shared_ptr<T>> out;
      for (int i = 0; i != nstate_; ++i) {
        auto tmp = c_.front()->clone();
        int k = 0;
        for (auto& iv : c_) {
          tmp->ax_plus_y(-vec_[i]*eig_->element(k++,i), iv);
        }
        k = 0;
        for (auto& iv : sigma_) {
          tmp->ax_plus_y(eig_->element(k++,i), iv);
        }
        out.push_back(tmp);
      }
      return out;
    }

    // returns ci vector
    std::vector<std::shared_ptr<T>> civec() {
      std::vector<std::shared_ptr<T>> out;
      for (int i = 0; i != nstate_; ++i) {
        auto tmp = c_.front()->clone();
        int k = 0;
        for (auto& iv : c_) {
          tmp->ax_plus_y(eig_->element(k++,i), iv);
        }
        out.push_back(tmp);
      }
      return out;
    }

    // make cc orthogonal to cc_ vectors
    double orthog(std::shared_ptr<T>& cc) { return cc->orthog(c_); }

    void collapse_subspace() {
      mat_->zero();
      for (int i = 0; i < nstate_; ++i)
        mat_->element(i,i) = vec_[i];

      overlap_->zero();
      for (int i = 0; i < nstate_; ++i)
        overlap_->element(i,i) = 1.0;

      {
        std::list<std::shared_ptr<const T>> collapsed_c;
        for (int i = 0; i < nstate_; ++i) {
          auto tmp_c = c_.front()->clone();
          int k = 0;
          for (auto ic = c_.begin(); ic != c_.end(); ++ic, ++k)
            tmp_c->ax_plus_y(eig_->element(k, i), *ic);
          collapsed_c.push_back(tmp_c);
        }
        c_ = std::move(collapsed_c);
      }
      {
        std::list<std::shared_ptr<const T>> collapsed_s;
        for (int i = 0; i < nstate_; ++i) {
          auto tmp_s = sigma_.front()->clone();
          int k = 0;
          for (auto is = sigma_.begin(); is != sigma_.end(); ++is, ++k)
            tmp_s->ax_plus_y(eig_->element(k, i), *is);
          collapsed_s.push_back(tmp_s);
        }
        sigma_ = std::move(collapsed_s);
      }

      std::fill(vec_.get() + nstate_, vec_.get() + max_, 0.0);
      size_ = nstate_;
    }

    void orthogonalize_subspace() {
      const int size = scr_->mdim();
      {
        std::list<std::shared_ptr<const T>> collapsed_c;
        for (int i = 0; i < size; ++i) {
          auto tmp_c = c_.front()->clone();
          int k = 0;
          for (auto ic = c_.begin(); ic != c_.end(); ++ic, ++k)
            tmp_c->ax_plus_y(scr_->element(k, i), *ic);
          collapsed_c.push_back(tmp_c);
        }
        c_ = std::move(collapsed_c);
      }
      {
        std::list<std::shared_ptr<const T>> collapsed_s;
        for (int i = 0; i < size; ++i) {
          auto tmp_s = sigma_.front()->clone();
          int k = 0;
          for (auto is = sigma_.begin(); is != sigma_.end(); ++is, ++k)
            tmp_s->ax_plus_y(scr_->element(k, i), *is);
          collapsed_s.push_back(tmp_s);
        }
        sigma_ = std::move(collapsed_s);
      }

      // update mat_ which should be diagonal
      mat_->zero();
      for (int i = 0; i < size; ++i)
        mat_->element(i,i) = vec_[i];

      overlap_->zero();
      for (int i = 0; i < size; ++i)
        overlap_->element(i,i) = 1.0;

      scr_ = std::make_shared<MatType>(size, size); scr_->unit();
      ovlp_scr_.reset();

      std::fill(vec_.get() + size, vec_.get() + max_, 0.0);
      size_ = size;
    }

};

}

#endif
