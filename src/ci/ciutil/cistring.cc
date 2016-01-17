//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cistring.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <cassert>
#include <src/ci/ciutil/bitutil.h>
#include <src/ci/ciutil/cistring.h>
#include <src/util/combination.hpp>

// explicit instantiation
template class bagel::CIString_base_impl<1,bagel::FCIString>;
template class bagel::CIString_base_impl<3,bagel::RASString>;

// instantiation of serialization code
BOOST_CLASS_EXPORT_IMPLEMENT(bagel::CIGraph)
BOOST_CLASS_EXPORT_IMPLEMENT(bagel::CIString_base)
BOOST_CLASS_EXPORT_IMPLEMENT(bagel::RASString_base)
BOOST_CLASS_EXPORT_IMPLEMENT(bagel::FCIString_base)
BOOST_CLASS_EXPORT_IMPLEMENT(bagel::RASString)
BOOST_CLASS_EXPORT_IMPLEMENT(bagel::FCIString)

using namespace std;
using namespace bagel;

CIGraph::CIGraph(const size_t nele, const size_t norb) : nele_(nele), norb_(norb), size_(1) {
  // Maximum active space size is hard-coded as nbit__ in src/util/constants.h
  if (norb_ > nbit__)
    throw runtime_error("Active space is too large.");

  if (nele*norb != 0) {
    weights_ = vector<size_t>(nele * norb, 0ull);

    const size_t nholes = norb - nele;
    for (size_t orb = 0; orb <= nholes; ++orb) {
      weight(orb, 0) = orb;
    }

    for (size_t ele = 1; ele < nele; ++ele) {
      for (size_t orb = ele; orb <= nholes + ele; ++orb) {
        weight(orb, ele) = weight(orb-1, ele-1) + weight(orb-1, ele);
      }
    }

    // compute size by "filling in" bottom row of graph
    size_ = nholes+1;

    for (size_t ele = 1; ele < nele; ++ele) {
      size_ += weight(ele+nholes, ele);
    }
  }
}


RASString::RASString(const size_t nele1, const size_t norb1, const size_t nele2, const size_t norb2, const size_t nele3, const size_t norb3, const size_t offset)
 : RASString_base{nele1, norb1, nele2, norb2, nele3, norb3, offset} {

  init();
}


void RASString::compute_strings_impl() {
  const size_t size = graphs_[0]->size()*graphs_[1]->size()*graphs_[2]->size();
  // Lexical ordering done, now fill in all the strings
  strings_ = vector<bitset<nbit__>>(size, bitset<nbit__>(0ul));

  const int nele1 = subspace_[0].first;
  const int nele2 = subspace_[1].first;
  const int nele3 = subspace_[2].first;
  const int norb1 = subspace_[0].second;
  const int norb2 = subspace_[1].second;
  const int norb3 = subspace_[2].second;

  size_t cnt = 0;
  vector<int> holes(norb1);
  iota(holes.begin(), holes.end(), 0);
  do {
    vector<int> active(norb2);
    iota(active.begin(), active.end(), norb1);
    do {
      vector<int> particles(norb3);
      iota(particles.begin(), particles.end(), norb1 + norb2);
      do {
        bitset<nbit__> bit(0ul);
        for (int i = 0; i != nele1; ++i) bit.set(holes[i]);
        for (int i = 0; i != nele2; ++i) bit.set(active[i]);
        for (int i = 0; i != nele3; ++i) bit.set(particles[i]);
        strings_[lexical_zero(bit)] = bit;

        ++cnt;
      } while (boost::next_combination(particles.begin(), particles.begin() + nele3, particles.end()));
    } while (boost::next_combination(active.begin(), active.begin() + nele2, active.end()));
  } while (boost::next_combination(holes.begin(), holes.begin() + nele1, holes.end()));
  assert(cnt == size);
}


FCIString::FCIString(const size_t nele1, const size_t norb1, const size_t offset)
 : FCIString_base{nele1, norb1, offset} {

  init();
}


void FCIString::compute_strings_impl() {
  vector<int> data(norb_);
  iota(data.begin(), data.end(), 0);
  // Lexical ordering done, now fill in all the strings
  strings_ = vector<bitset<nbit__>>(graphs_[0]->size(), bitset<nbit__>(0ul));
  size_t cnt = 0;
  do {
    bitset<nbit__> bit(0lu);
    for (int i=0; i!=nele_; ++i) bit.set(data[i]);
    strings_[lexical_zero(bit)] = bit;
    ++cnt;
  } while (boost::next_combination(data.begin(), data.begin()+nele_, data.end()));
  assert(cnt == graphs_[0]->size());
}


