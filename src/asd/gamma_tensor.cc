//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gamma_tensor.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
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

#include <src/asd/gamma_tensor.h>

using namespace bagel;

namespace bagel {
namespace {
struct init {
  static std::list<std::list<GammaSQ>> call() { return {
    {GammaSQ::CreateAlpha},
    {GammaSQ::CreateBeta},
    {GammaSQ::AnnihilateAlpha},
    {GammaSQ::AnnihilateBeta},
    {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},
    {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},
    {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha},
    {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta},
    {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha},
    {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},
    {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta},
    {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha},
    {GammaSQ::CreateBeta,  GammaSQ::CreateBeta},
    {GammaSQ::CreateAlpha, GammaSQ::CreateBeta},
    {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha},
    {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},
    {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha},
    {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta},
    {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},
    {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha},
    {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},
    {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta}
  }; }
};
}
}

const std::list<std::list<GammaSQ>> GammaTensor::oplist_(init::call());


GammaTensor& GammaTensor::operator=(const GammaTensor& o) {
  auto ptr = o.sparse_.begin();
  for (auto& i : sparse_) {
    assert(i.first == ptr->first);
    *i.second = *ptr->second; // data copy
    ++ptr;
  }
  return *this;
}


GammaTensor& GammaTensor::operator=(GammaTensor&& o) {
  auto ptr = o.sparse_.begin();
  for (auto& i : sparse_) {
    assert(i.first == ptr->first);
    i.second = ptr->second; // pointer copy
    ++ptr;
  }
  return *this;
}
