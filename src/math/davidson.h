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

#include <vector>
#include <src/math/algo.h>
#include <src/util/f77.h>

namespace bagel {

template <typename T, class MatType = Matrix>
class DavidsonDiag {
  protected:
    struct BasisPair {
      std::shared_ptr<const T> cc;
      std::shared_ptr<const T> sigma;
      BasisPair(std::shared_ptr<const T> a, std::shared_ptr<const T> b) : cc(a), sigma(b) { }
    };

    const int nstate_;
    const int max_;
    int size_;

    std::vector<std::shared_ptr<BasisPair>> basis_;

    // Hamiltonian
    std::shared_ptr<MatType> mat_;
    // eivenvalues
    std::unique_ptr<double[]> vec_;
    // an eigenvector
    std::shared_ptr<MatType> eig_;
    // overlap matrix
    std::shared_ptr<MatType> overlap_;

    std::vector<bool> converged_;

  public:
    // Davidson with periodic collapse of the subspace
    DavidsonDiag(int n, int max) : nstate_(n), max_((max+1)*n), size_(0), vec_(new double[max_]), converged_(n,false) {
      if (max < 2) throw std::runtime_error("Davidson diagonalization requires at least two trial vectors per root.");
    }

    double compute(std::shared_ptr<const T> cc, std::shared_ptr<const T> cs) {
      assert(nstate_ == 1);
      return compute(std::vector<std::shared_ptr<const T>>{cc},
                     std::vector<std::shared_ptr<const T>>{cs}).front();
    }

    std::vector<double> compute(std::vector<std::shared_ptr<const T>> cc, std::vector<std::shared_ptr<const T>> cs) {
      // reset the convergence flags
      std::fill(converged_.begin(), converged_.end(), false);

      // new pairs
      std::vector<std::shared_ptr<BasisPair>> newbasis;
      assert(cc.size() == nstate_ && cs.size() == nstate_);
      for (int ic = 0; ic < nstate_; ++ic) {
        assert(!cc[ic] == !cs[ic]);
        if (cc[ic] && cs[ic])
          newbasis.push_back(std::make_shared<BasisPair>(cc[ic], cs[ic]));
        else
          converged_[ic] = true;
      }

      // adding new matrix elements
      {
        const int n = newbasis.size();
        mat_ = mat_ ? mat_->resize(size_+n, size_+n) : std::make_shared<MatType>(n, n);
        overlap_ = overlap_ ? overlap_->resize(size_+n, size_+n) : std::make_shared<MatType>(n, n);
      }

      basis_.insert(basis_.end(), newbasis.begin(), newbasis.end());
      for (auto& ib : newbasis) {
        ++size_;
        int i = 0;
        for (auto& b : basis_) {
          if (i > size_-1) break;
          mat_->element(i, size_-1) = b->cc->dot_product(ib->sigma);
          mat_->element(size_-1, i) = detail::conj(mat_->element(i, size_-1));

          overlap_->element(i, size_-1) = b->cc->dot_product(ib->cc);
          overlap_->element(size_-1, i) = detail::conj(overlap_->element(i, size_-1));
          ++i;
        }
      }

      // canonical orthogonalization
      std::shared_ptr<const MatType> ovlp_scr = overlap_->tildex();

      // diagonalize matrix to get
      eig_ = std::make_shared<MatType>(*ovlp_scr % *mat_ * *ovlp_scr);
      eig_->diagonalize(vec_.get());
      eig_ = std::make_shared<MatType>(*ovlp_scr * *eig_);
      eig_ = eig_->slice(0,nstate_);

      // first basis vector is always the current best guess
      std::vector<std::shared_ptr<T>> cv = civec();
      std::vector<std::shared_ptr<T>> sv = sigmavec();
      for (int i = 0; i != nstate_; ++i)
        basis_[i] = std::make_shared<BasisPair>(cv[i], sv[i]);

      // due to this, we need to transform mat_ and overlap_
      auto trans = eig_->resize(eig_->ndim(), eig_->ndim());
      for (int i = nstate_; i != eig_->ndim(); ++i)
        trans->element(i, i) = 1.0;
      mat_ = std::make_shared<MatType>(*trans % *mat_ * *trans);
      overlap_ = std::make_shared<MatType>(*trans % *overlap_ * *trans);

      eig_->zero();
      for (int i = 0; i != nstate_; ++i)
        eig_->element(i, i) = 1.0;

      // possibly reduce the dimension
      assert(size_ == basis_.size());
      if (size_ > max_-nstate_) {
        std::map<int, int> remove;
        const int soff = size_ - newbasis.size();
        for (int i = 0; i != nstate_; ++i) {
          if (converged_[i]) continue;
          // a vector with largest weight will be removed.
          int n = 0;
          double abs = 1.0e10;
          for (int j = nstate_; j < soff; ++j) {
            if (std::abs(trans->element(j, i)) < abs) {
              if (remove.find(j) != remove.end()) continue;
              abs = std::abs(trans->element(j, i));
              n = j;
            }
          }
          remove.insert(std::make_pair(n, soff+remove.size()));
        }
        assert(newbasis.size() == remove.size());
        std::cout << "    ** throwing out " << remove.size() << " trial vectors **" << std::endl;
        for (auto m : remove) {
          basis_[m.first] = basis_[m.second];
          mat_->copy_block(0, m.first, size_, 1, mat_->get_submatrix(0, m.second, size_, 1));
          mat_->copy_block(m.first, 0, 1, size_, mat_->get_submatrix(m.second, 0, 1, size_));
          overlap_->copy_block(0, m.first, size_, 1, overlap_->get_submatrix(0, m.second, size_, 1));
          overlap_->copy_block(m.first, 0, 1, size_, overlap_->get_submatrix(m.second, 0, 1, size_));

          trans->copy_block(m.first, 0, 1, nstate_, trans->get_submatrix(m.second, 0, 1, nstate_));
        }
        basis_ = std::vector<std::shared_ptr<BasisPair>>(basis_.begin(), basis_.end()-remove.size());
        size_ = basis_.size();
        mat_ = mat_->get_submatrix(0, 0, size_, size_);
        overlap_ = overlap_->get_submatrix(0, 0, size_, size_);
      }

      return std::vector<double>(vec_.get(), vec_.get()+nstate_);
    }

    // perhaps can be cleaner.
    std::vector<std::shared_ptr<T>> residual() {
      std::vector<std::shared_ptr<T>> out;
      for (int i = 0; i != nstate_; ++i) {
        auto tmp = basis_.front()->cc->clone();
        int k = 0;
        for (auto& iv : basis_) {
          tmp->ax_plus_y(-vec_[i]*eig_->element(k++,i), iv->cc);
        }
        k = 0;
        for (auto& iv : basis_) {
          tmp->ax_plus_y(eig_->element(k++,i), iv->sigma);
        }
        out.push_back(tmp);
      }
      return out;
    }

    // returns ci vector
    std::vector<std::shared_ptr<T>> civec() {
      std::vector<std::shared_ptr<T>> out;
      for (int i = 0; i != nstate_; ++i) {
        auto tmp = basis_.front()->cc->clone();
        int k = 0;
        for (auto& iv : basis_) {
          tmp->ax_plus_y(eig_->element(k++,i), iv->cc);
        }
        out.push_back(tmp);
      }
      return out;
    }

    // return sigma vector
    std::vector<std::shared_ptr<T>> sigmavec() {
      std::vector<std::shared_ptr<T>> out;
      for (int i = 0; i != nstate_; ++i) {
        auto tmp = basis_.front()->sigma->clone();
        int k = 0;
        for (auto& iv : basis_) {
          tmp->ax_plus_y(eig_->element(k++,i), iv->sigma);
        }
        out.push_back(tmp);
      }
      return out;
    }

};

}

#endif
