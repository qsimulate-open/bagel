//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: davidson.h
// Copyright (C) 2011 Toru Shiozaki
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

// T should have
//  - double dot_product(const T&)
//  - void ax_plus_y(double, const T&) // added to self

#ifndef __BAGEL_UTIL_DAVIDSON
#define __BAGEL_UTIL_DAVIDSON

#include <vector>
#include <src/util/math/algo.h>
#include <src/util/math/matrix.h>
#include <src/util/f77.h>
#include <src/util/serialization.h>

namespace bagel {

template <typename T, typename U, class MatType = Matrix>
class DavidsonDiag_ {
  protected:
    struct BasisPair {
      public:
        std::shared_ptr<const T> cc;
        std::shared_ptr<const U> sigma;
        BasisPair() { }
        BasisPair(std::shared_ptr<const T> a, std::shared_ptr<const U> b) : cc(a), sigma(b) { }
      private:
        // serialization
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive& ar, const unsigned int) { ar & cc & sigma; }
    };

    int nstate_;
    int max_;
    int size_;

    std::vector<std::shared_ptr<BasisPair>> basis_;

    // Hamiltonian
    std::shared_ptr<MatType> mat_;
    // eivenvalues
    VectorB vec_;
    // an eigenvector
    std::shared_ptr<MatType> eig_;
    // overlap matrix
    std::shared_ptr<MatType> overlap_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & nstate_ & max_ & size_ & basis_ & mat_ & vec_ & eig_ & overlap_;
    }

  public:
    // Davidson with periodic collapse of the subspace
    DavidsonDiag_() { }
    DavidsonDiag_(int n, int max) : nstate_(n), max_((max+1)*n), size_(0), vec_(max_) {
      if (max < 2) throw std::runtime_error("Davidson diagonalization requires at least two trial vectors per root.");
    }

    double compute(std::shared_ptr<const T> cc, std::shared_ptr<const U> cs) {
      assert(nstate_ == 1);
      return compute(std::vector<std::shared_ptr<const T>>{cc},
                     std::vector<std::shared_ptr<const U>>{cs}).front();
    }

    std::vector<double> compute(std::vector<std::shared_ptr<const T>> cc, std::vector<std::shared_ptr<const U>> cs) {
      // reset the convergence flags
      std::vector<bool> converged(nstate_, false);

      // new pairs
      std::vector<std::shared_ptr<BasisPair>> newbasis;
      assert(cc.size() == nstate_ && cs.size() == nstate_);
      for (int ic = 0; ic < nstate_; ++ic) {
        assert(!cc[ic] == !cs[ic]);
        if (cc[ic] && cs[ic])
          newbasis.push_back(std::make_shared<BasisPair>(cc[ic], cs[ic]));
        else
          converged[ic] = true;
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

      mat_->synchronize();
      overlap_->synchronize();

      // canonical orthogonalization
      std::shared_ptr<const MatType> ovlp_scr = overlap_->tildex();
      if (ovlp_scr->mdim() < nstate_)
        throw std::runtime_error("Too much linear dependency in guess vectors provided to DavidsonDiag; cannot obtain the requested number of states.");

      // diagonalize matrix to get
      eig_ = std::make_shared<MatType>(*ovlp_scr % *mat_ * *ovlp_scr);
      eig_->diagonalize(vec_);
      eig_ = std::make_shared<MatType>(*ovlp_scr * eig_->slice(0,nstate_));
      eig_->synchronize();

      // first basis vector is always the current best guess
      std::vector<std::shared_ptr<T>> cv = civec();
      std::vector<std::shared_ptr<U>> sv = sigmavec();
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
          if (converged[i]) continue;
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
          remove.emplace(n, soff+remove.size());
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

      mat_->synchronize();
      overlap_->synchronize();

      return std::vector<double>(vec_.begin(), vec_.begin()+nstate_);
    }

    // perhaps can be cleaner.
    std::vector<std::shared_ptr<U>> residual() {
      std::vector<std::shared_ptr<U>> out;
      for (int i = 0; i != nstate_; ++i) {
        auto tmp = basis_.front()->sigma->clone();
        int k = 0;
        for (auto& iv : basis_) {
          if (std::abs(eig_->element(k++,i)) > 1.0e-16)
            tmp->ax_plus_y(-vec_(i)*eig_->element(k-1,i), iv->cc);
        }
        k = 0;
        for (auto& iv : basis_) {
          if (std::abs(eig_->element(k++,i)) > 1.0e-16)
            tmp->ax_plus_y(eig_->element(k-1,i), iv->sigma);
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
        tmp->synchronize();
        out.push_back(tmp);
      }
      return out;
    }

    // return sigma vector
    std::vector<std::shared_ptr<U>> sigmavec() {
      std::vector<std::shared_ptr<U>> out;
      for (int i = 0; i != nstate_; ++i) {
        auto tmp = basis_.front()->sigma->clone();
        int k = 0;
        for (auto& iv : basis_) {
          tmp->ax_plus_y(eig_->element(k++,i), iv->sigma);
        }
        tmp->synchronize();
        out.push_back(tmp);
      }
      return out;
    }

};

template <typename T, class MatType = Matrix>
using DavidsonDiag = DavidsonDiag_<T, T, MatType>;

}

#endif
