//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gamma_tree.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Toru Shiozaki
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

#include <src/asd/gamma_tree.h>

using namespace std;
using namespace bagel;

template<typename VecType>
void GammaBranch<VecType>::insert(shared_ptr<const VecType> bra, const size_t bra_tag, const std::list<GammaSQ>& gsq) {
  if (gsq.empty()) {
    bras_.emplace(bra_tag, bra);
  } else {
    auto first = gsq.back();
    auto rest = gsq; rest.pop_back();
    shared_ptr<GammaBranch<VecType>> target = branches_[static_cast<int>(first)];

    target->activate();
    target->insert(bra, bra_tag, rest);
  }
}


template<typename VecType>
shared_ptr<const Matrix> GammaBranch<VecType>::search(const size_t tag, const std::list<GammaSQ>& gsq) const {
  if (gsq.empty()) {
    assert(gammas_.find(tag)!=gammas_.end()); return gammas_.find(tag)->second;
  } else {
    auto first = gsq.back();
    auto rest = gsq; rest.pop_back();
    return branch(first)->search(tag, rest);
  }
}


template<typename VecType>
shared_ptr<Matrix> GammaBranch<VecType>::search(const size_t tag, const std::list<GammaSQ>& gsq) {
  if (gsq.empty()) {
    assert(gammas_.find(tag)!=gammas_.end()); return gammas_.find(tag)->second;
  } else {
    auto first = gsq.back();
    auto rest = gsq; rest.pop_back();
    return branch(first)->search(tag, rest);
  }
}


template<typename VecType>
bool GammaBranch<VecType>::exist(const size_t tag, const std::list<GammaSQ>& gsq) const {
  if (gsq.empty())
    return gammas_.find(tag) != gammas_.end();
  else {
    auto first = gsq.back();
    auto rest = gsq; rest.pop_back();
    return branch(first) ? branch(first)->exist(tag, rest) : false;
  }
}


template<typename VecType>
bool GammaBranch<VecType>::if_contributes(std::set<int> needed) {
  bool contributes = false;
  for (const int& b : needed) {
    if (branch(b)) contributes |= branch(b)->active();
  }
  if (!contributes) {
    for (int i = 0; i < 4; ++i) {
      if (branch(i)) {
        if (branch(i)->active()) contributes |= branch(i)->if_contributes(needed);
      }
    }
  }
  return contributes;
}


template<typename VecType>
bool GammaBranch<VecType>::is_distterminal() {
  return false;
#if 0
  bool not_terminal = false;
  if (branch(GammaSQ::CreateAlpha))
    not_terminal |= branch(GammaSQ::CreateAlpha)->active();
  if (branch(GammaSQ::AnnihilateAlpha))
    not_terminal |= branch(GammaSQ::AnnihilateAlpha)->active();
  if (!not_terminal) {
    if (branch(GammaSQ::CreateBeta) && branch(GammaSQ::CreateBeta)->active())
      not_terminal |= branch(GammaSQ::CreateBeta)->is_distterminal();
    if (branch(GammaSQ::AnnihilateBeta) && branch(GammaSQ::AnnihilateBeta)->active())
      not_terminal |= branch(GammaSQ::AnnihilateBeta)->is_distterminal();
  }
  return !not_terminal;
#endif
}

template<typename VecType>
GammaTree<VecType>::GammaTree(shared_ptr<const VecType> ket) : ket_(ket) {
  base_ = make_shared<GammaBranch<VecType>>();
  constexpr int nops = 4;

  for (int i = 0; i < nops; ++i) {
    base_->branch(i) = make_shared<GammaBranch<VecType>>();
    for (int j = 0; j < nops; ++j) {
      base_->branch(i)->branch(j) = make_shared<GammaBranch<VecType>>();
      for (int k = 0; k < nops; ++k)
        base_->branch(i)->branch(j)->branch(k) = make_shared<GammaBranch<VecType>>();
    }
  }
}


template class bagel::GammaBranch<CASDvec>;
template class bagel::GammaBranch<RASDvec>;
template class bagel::GammaTree<CASDvec>;
template class bagel::GammaTree<RASDvec>;
