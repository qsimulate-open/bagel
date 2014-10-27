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
    {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, //a' aa
    {GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  //b' bb
    {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateAlpha}, //a' ba
    {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta},  //b' ab
    {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, // b'a' a
    {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, // a'a' a
    {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  // a'b' b
    {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  // b'b' b
//ADDED
    {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  //a' bb
    {GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, //b' aa

    {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta},  //a'a' b
    {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha}, //b'b' a

    {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}, //a'a'a'
    {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta},  //a'a'b'
    {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateBeta},  //a'b'b'
    {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta},  //b'b'b'

    {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, //aaa 
    {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta},  //aab 
    {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  //abb 
    {GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  //bbb 

    {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, //a'a'aa
    {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta},  //a'b'ab
    {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  //b'b'bb

    {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, //a'b'aa
    {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta},  //b'b'ab

    {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, //a'a'a'a
    {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  //a'a'b'b
    {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, //a'b'a'a
    {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  //b'a'b'b
    {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta},  //b'b'b'b
    {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, //b'b'a'a

    {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, //a'a'a'aa
    {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta},  //a'b'a'ab
    {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta},  //a'b'b'bb
    {GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha}, //b'a'a'aa
    {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta},  //b'b'a'ab
    {GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta,  GammaSQ::AnnihilateBeta}  //b'b'b'bb

//END ADDED
  }; }
};
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
