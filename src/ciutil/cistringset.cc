//
// BAGEL - Parallel electron correlation program.
// Filename: cistringset.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@u.northwestern.edu>
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

#include <src/ciutil/cistringset.h>

template class bagel::CIStringSet<bagel::FCIString>;
template class bagel::CIStringSet<bagel::RASString>;

using namespace std;
using namespace bagel;

// FCI version
template<>
void CIStringSet<FCIString>::construct_phi() {
  phi_ = make_shared<StringMap>(norb()*(norb()+1)/2);
  phi_->reserve(size());

  uncompressed_phi_ = make_shared<StringMap>(norb()*norb());
  uncompressed_phi_->reserve(size());

  for (auto& istring : strings_) {
    for (unsigned int i = 0; i != norb(); ++i) { // annihilation
      // compress_ means that we store info only for i <= j
      if (istring[i]) {
        const unsigned int source = lexical_zero(istring);
        bitset<nbit__> nbit = istring; nbit.reset(i); // annihilated.
        for (unsigned int j = 0; j != norb(); ++j) { // creation
          if (!nbit[j]) {
            bitset<nbit__> mbit = nbit;
            mbit.set(j);
            int minij, maxij;
            tie(minij, maxij) = minmax(i,j);
            auto detmap = DetMap(lexical_zero(mbit), sign(mbit, i, j), source, i+norb()*j);
            (*phi_)[minij+((maxij*(maxij+1))>>1)].push_back(detmap);
            (*uncompressed_phi_)[i + j*norb()].push_back(detmap);
          }
        }
      }
    }
  }
}

// RAS version
template<>
void CIStringSet<RASString>::construct_phi() {
  phi_ = make_shared<StringMap>(size_);
  phi_->reserve(norb_*norb_);

  unordered_map<size_t, size_t> lexmap;
  for (size_t i = 0; i < size_; ++i)
    lexmap[strings_[i].to_ullong()] = i;

  size_t tindex = 0;
  for (auto& istring : strings_) {
    for (int j = 0; j < norb_; ++j) {
      if (!istring[j]) continue;
      bitset<nbit__> intermediatebit = istring; intermediatebit.reset(j);
      for (int i = 0; i < norb_; ++i) {
        if (intermediatebit[i]) continue;
        bitset<nbit__> sourcebit = intermediatebit; sourcebit.set(i);
        if (allowed(sourcebit)) {
          assert(lexmap.find(sourcebit.to_ullong()) != lexmap.end());
          (*phi_)[tindex].emplace_back(tindex, sign(istring, i, j), lexmap[sourcebit.to_ullong()], j+i*norb_);
        }
      }
    }
    (*phi_)[tindex++].shrink_to_fit();
  }
}
