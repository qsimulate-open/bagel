//
// BAGEL - Parallel electron correlation program.
// Filename: gamma_tensor.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
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

#include <src/asd/gamma_tensor.h>

using namespace bagel;

namespace {
struct init {
  static std::list<std::list<GammaSQ>> call() { return {
    {GammaSQ::CreateAlpha},
    {GammaSQ::CreateBeta},
    {GammaSQ::AnnihilateAlpha},
    {GammaSQ::AnnihilateBeta},
    {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha},
    {GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta},
    {GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta},
    {GammaSQ::AnnihilateBeta,  GammaSQ::CreateAlpha},
    {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha},
    {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},
    {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha},
    {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha},
    {GammaSQ::CreateBeta,  GammaSQ::CreateBeta},
    {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha},
    {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha},
    {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta},
    {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::CreateAlpha},
    {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta},
    {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta},
    {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha},
    {GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateAlpha},
    {GammaSQ::AnnihilateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta}
  }; }
};
}

const std::list<std::list<GammaSQ>> GammaTensor::oplist_(init::call());


GammaTensor& GammaTensor::operator=(const GammaTensor& o) {
  assert(sparse_.size() == o.sparse_.size() && norb_ == o.norb_);
  auto ptr = o.sparse_.begin();
  for (auto& i : sparse_) {
    assert(i.first == ptr->first);
    *i.second = *ptr->second; // data copy
    ++ptr;
  }
  return *this;
}


GammaTensor& GammaTensor::operator=(GammaTensor&& o) {
  assert(sparse_.size() == o.sparse_.size() && norb_ == o.norb_);
  auto ptr = o.sparse_.begin();
  for (auto& i : sparse_) {
    assert(i.first == ptr->first);
    i.second = ptr->second; // pointer copy
    ++ptr;
  }
  return *this;
}
