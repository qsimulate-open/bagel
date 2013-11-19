//
// BAGEL - Parallel electron correlation program.
// Filename: ras/stringspace.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <cassert>

#include <src/ras/stringspace.h>
#include <src/util/combination.hpp>

using namespace std;
using namespace bagel;

RASGraph::RASGraph(const size_t nele, const size_t norb) : nele_(nele), norb_(norb), size_(1) {
  if ( nele*norb != 0 ) {
    weights_ = unique_ptr<size_t[]>(new size_t[nele * norb]);
    fill_n(weights_.get(), nele * norb, 0ull);

    Comb comb;

    size_ = comb.c(norb, nele);

    const size_t nholes = norb - nele;
    for(size_t k = 1; k <= nele; ++k) {
      for (size_t l = k; l < nholes + k; ++l) {
        size_t node_val = comb.c(l, k);
        weight(l, k-1) = node_val;
      }
    }
  }
}

StringSpace::StringSpace(const int nele1, const int norb1, const int nele2, const int norb2, const int nele3, const int norb3, const size_t offset) :
  ras_{{make_pair(nele1, norb1), make_pair(nele2, norb2), make_pair(nele3, norb3)}},
    graphs_{{ make_shared<RASGraph>(nele1, norb1), make_shared<RASGraph>(nele2, norb2), make_shared<RASGraph>(nele3, norb3) }},
    dist_(graphs_[0]->size()*graphs_[1]->size()*graphs_[2]->size(), mpi__->size()),
    norb_(norb1 + norb2 + norb3), nele_( nele1 + nele2 + nele3 ), offset_(offset)
{
  const size_t size = graphs_[0]->size()*graphs_[1]->size()*graphs_[2]->size();

  // Lexical ordering done, now fill in all the strings
  strings_ = vector<bitset<nbit__>>(size, bitset<nbit__>(0ul));

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
        strings_[lexical<0>(bit)] = bit;

        ++cnt;
      } while (boost::next_combination(particles.begin(), particles.begin() + nele3, particles.end()));
    } while (boost::next_combination(active.begin(), active.begin() + nele2, active.end()));
  } while (boost::next_combination(holes.begin(), holes.begin() + nele1, holes.end()));

  assert(cnt == size);
}
