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
#include <list>

#include <src/math/algo.h>
#include <src/util/f77.h>

namespace bagel {

template <typename T, class MatType = Matrix>
class DavidsonDiag {
  protected:
    struct RitzPair {
      std::shared_ptr<const T> cc;
      std::shared_ptr<const T> sigma;
      int position;

      RitzPair() : position(-1) {}
      RitzPair(std::shared_ptr<const T> c, std::shared_ptr<const T> s, const int p) : cc(c), sigma(s), position(p) {}
    };

  protected:
    const int nstate_;
    const int max_;
    int size_;

    std::vector<std::pair<RitzPair, std::list<RitzPair>>> trials_;

    // projected hamiltonian
    std::shared_ptr<MatType> Hsub_;
    // eigenvalues
    std::unique_ptr<double[]> vec_;
    // overlap matrix
    std::shared_ptr<MatType> overlap_;

  public:
    // Davidson with periodic collapse of the subspace
    DavidsonDiag(const int nstate, const int max) :
                  nstate_(nstate), max_(max), size_(0), trials_(nstate_),
                  Hsub_(std::make_shared<MatType>((max_+1)*nstate_, (max_+1)*nstate_, true)),
                  vec_(new double[nstate_]), overlap_(std::make_shared<MatType>((max_+1)*nstate_,(max_+1)*nstate_,true))
    { }

    double compute(std::shared_ptr<const T> cc, std::shared_ptr<const T> cs) {
      assert(nstate_ == 1);
      return compute(std::vector<std::shared_ptr<const T>>{cc},
                     std::vector<std::shared_ptr<const T>>{cs}).front();
    }

    // expecting blank shared_ptrs for converged roots
    std::vector<double> compute(std::vector<std::shared_ptr<const T>> cc, std::vector<std::shared_ptr<const T>> cs) {
      std::vector<RitzPair> corrections;

      // adding new matrix elements
      assert( cc.size() == nstate_ && cs.size() == nstate_);
      for (int ic = 0; ic < nstate_; ++ic) {
        assert(!cc[ic] == !cs[ic]);
        if (cc[ic] && cs[ic]) corrections.emplace_back(cc[ic], cs[ic], ic);
      }

      auto expandedH = std::make_shared<MatType>(size_ + corrections.size(), size_ + corrections.size(), true);
      expandedH->copy_block(0, 0, size_, size_, Hsub_->get_submatrix(0, 0, size_, size_));

      auto expandedS = std::make_shared<MatType>(size_ + corrections.size(), size_ + corrections.size(), true);
      expandedS->copy_block(0, 0, size_, size_, overlap_->get_submatrix(0, 0, size_, size_));

      int current = size_;
      for (auto& pair : corrections) {
        for (auto& state : trials_) {
          for (auto& i : state.second) {
            expandedH->element(current, i.position) = pair.cc->dot_product(*i.sigma);
            expandedH->element(i.position, current) = detail::conj(expandedH->element(current, i.position));

            expandedS->element(current, i.position) = pair.cc->dot_product(*i.cc);
            expandedS->element(i.position, current) = detail::conj(expandedS->element(current, i.position));
          }
          if (state.first.position >= 0) {
            expandedH->element(current, state.first.position) = pair.cc->dot_product(*state.first.sigma);
            expandedH->element(state.first.position, current) = detail::conj(expandedH->element(current, state.first.position));

            expandedS->element(current, state.first.position) = pair.cc->dot_product(*state.first.cc);
            expandedS->element(state.first.position, current) = detail::conj(expandedS->element(current, state.first.position));
          }
        }
        for (int i = size_; i < current; ++i) {
          expandedH->element(current, i) = pair.cc->dot_product(*corrections[i-size_].sigma);
          expandedH->element(i, current) = detail::conj(expandedH->element(current, i));

          expandedS->element(current, i) = pair.cc->dot_product(*corrections[i-size_].cc);
          expandedS->element(i, current) = detail::conj(expandedS->element(current, i));
        }

        expandedH->element(current, current) = pair.cc->dot_product(*pair.sigma);
        expandedS->element(current, current) = pair.cc->dot_product(*pair.cc);

        ++current;
      }

      // solve generalized eigenvalue problem
      std::shared_ptr<MatType> Sinvhalf = expandedS->tildex();

      std::vector<double> eigenvalues(current, 0.0);
      expandedH = std::make_shared<MatType>(*Sinvhalf % *expandedH * *Sinvhalf);
      expandedH->diagonalize(eigenvalues.data());
      std::copy_n(eigenvalues.data(), nstate_, vec_.get());
      expandedH = std::make_shared<MatType>(*Sinvhalf * *expandedH);

      // form new trial vectors for input vectors
      int current_size = size_;
      for (auto& pair : corrections) {
        std::shared_ptr<T> nccvec = pair.cc->clone();
        std::shared_ptr<T> nsgvec = pair.sigma->clone();
        for (auto& state : trials_) {
          for (auto& i : state.second) {
            nccvec->ax_plus_y(expandedH->element(i.position, pair.position), *i.cc);
            nsgvec->ax_plus_y(expandedH->element(i.position, pair.position), *i.sigma);
          }
          if (state.first.position >= 0) {
            nccvec->ax_plus_y(expandedH->element(state.first.position, pair.position), *state.first.cc);
            nsgvec->ax_plus_y(expandedH->element(state.first.position, pair.position), *state.first.sigma);
          }
        }
        for (int i = 0; i < corrections.size(); ++i) {
          nccvec->ax_plus_y(expandedH->element(i+size_, pair.position), *corrections[i].cc);
          nsgvec->ax_plus_y(expandedH->element(i+size_, pair.position), *corrections[i].sigma);
        }

        std::list<RitzPair>& correctionlist = trials_[pair.position].second;
        RitzPair& old_trial = trials_[pair.position].first;
        if (correctionlist.size() == max_) {
          const int trial_position = (old_trial.position < 0) ? current_size : old_trial.position;
          trials_[pair.position].first = RitzPair(nccvec, nsgvec, trial_position);
          insert_ritz_pair(trials_[pair.position].first);
          if ( trial_position == current_size ) ++current_size;
        }
        else {
          trials_[pair.position].first = RitzPair(nccvec, nsgvec, -1);
        }

        const int corr_position = (correctionlist.size() == max_) ? correctionlist.back().position : current_size;
        correctionlist.emplace_front(pair.cc, pair.sigma, corr_position);
        if (correctionlist.size() > max_)
          correctionlist.pop_back();
        else
          ++current_size;
        insert_ritz_pair(correctionlist.front());
      }
      size_ = current_size;

      return std::vector<double>(vec_.get(), vec_.get()+nstate_);
    }

    // perhaps can be cleaner.
    std::vector<std::shared_ptr<T>> residual() {
      std::vector<std::shared_ptr<T>> out;
      for (int i = 0; i != nstate_; ++i) {
        auto tmp = std::make_shared<T>(*trials_[i].first.sigma);
        tmp->ax_plus_y(-vec_[i], trials_[i].first.cc);

        out.push_back(tmp);
      }
      return out;
    }

    // returns ci vector
    std::vector<std::shared_ptr<T>> civec() {
      std::vector<std::shared_ptr<T>> out;
      for (auto& state : trials_)
        out.push_back(std::make_shared<T>(*state.first.cc));
      return out;
    }

    // make cc orthogonal to cc_ vectors
    double orthog(std::shared_ptr<T>& cc) {
#if 1 // toggle orthogonalization with the current set of best guess solutions
      std::list<std::shared_ptr<const T>> orthog_list;
      for (auto& state : trials_)
        orthog_list.push_back(state.first.cc);

      return cc->orthog(orthog_list);
#else
      return 1.0;
#endif
    }

  private:
    void insert_ritz_pair(RitzPair& p) {
      const int position = p.position;
      for (auto& state : trials_) {
        for (auto& i : state.second) {
          Hsub_->element(position, i.position) = p.cc->dot_product(*i.sigma);
          overlap_->element(position, i.position) = p.cc->dot_product(*i.cc);

          if (i.position != position) {
            Hsub_->element(i.position, position) = detail::conj(Hsub_->element(position, i.position));
            overlap_->element(i.position, position) = detail::conj(overlap_->element(position, i.position));
          }
        }
        if (state.first.position >= 0) {
          Hsub_->element(position, state.first.position) = p.cc->dot_product(*state.first.sigma);
          overlap_->element(position, state.first.position) = p.cc->dot_product(*state.first.cc);

          if (state.first.position != position) {
            Hsub_->element(state.first.position, position) = detail::conj(Hsub_->element(position, state.first.position));
            overlap_->element(state.first.position, position) = detail::conj(overlap_->element(position, state.first.position));
          }
        }
      }
    }

};

}

#endif
