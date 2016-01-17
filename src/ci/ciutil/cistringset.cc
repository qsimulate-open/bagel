//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cistringset.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@u.northwestern.edu>
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

#include <src/ci/ciutil/cistringset.h>

template class bagel::CIStringSet<bagel::FCIString>;
template class bagel::CIStringSet<bagel::RASString>;

BOOST_CLASS_EXPORT_IMPLEMENT(bagel::CIStringSet<bagel::FCIString>)
BOOST_CLASS_EXPORT_IMPLEMENT(bagel::CIStringSet<bagel::RASString>)

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

  unordered_map<bitset<nbit__>, size_t> lexmap;
  for (size_t i = 0; i < size_; ++i)
    lexmap[strings_[i]] = i;

  size_t tindex = 0;
  for (auto& istring : strings_) {
    for (int j = 0; j < norb_; ++j) {
      if (!istring[j]) continue;
      bitset<nbit__> intermediatebit = istring; intermediatebit.reset(j);
      for (int i = 0; i < norb_; ++i) {
        if (intermediatebit[i]) continue;
        bitset<nbit__> sourcebit = intermediatebit; sourcebit.set(i);
        if (allowed(sourcebit)) {
          assert(lexmap.find(sourcebit) != lexmap.end());
          (*phi_)[tindex].emplace_back(tindex, sign(istring, i, j), lexmap[sourcebit], j+i*norb_);
        }
      }
    }
    (*phi_)[tindex++].shrink_to_fit();
  }
}
