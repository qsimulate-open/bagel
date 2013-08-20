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

#include <src/ras/stringspace.h>

using namespace std;
using namespace bagel;

StringSpace::StringSpace(const int nele1, const int norb1, const int nele2, const int norb2, const int nele3, const int norb3, const int offset) :
  ras_{ make_pair(nele1, norb1), make_pair(nele2, norb2), make_pair(nele3, norb3) }, nele_( nele1 + nele2 + nele3 ), norb_(norb1 + norb2 + norb3), offset_(offset)
{
  RASGraph graph(norb_ + 1, nele_ + 1);

  auto fill_ras_graph = [&graph](const int istart, const int jstart, const int norb, const int nele) {
    const int nholes = norb - nele;

    for (int i = 0; i <= nholes; ++i) {
      for (int j = 0; j <= nele; ++j) {
        if (i + j + istart == 0) graph(i+j+istart,j+jstart) = 1;
        else if (j + jstart == 0) graph(i+j+istart,j+jstart) = graph(i+j+istart-1,j+jstart);
        else graph(i+j+istart,j+jstart) = max(0, graph(i+j+istart-1,j+jstart-1)) + max(0,graph(i+j+istart-1,j+jstart));
      }
    }
  };

  fill_ras_graph(0, 0, norb1, nele1);
  fill_ras_graph(norb1, nele1, norb2, nele2);
  fill_ras_graph(norb1 + norb2, nele1 + nele2, norb3, nele3);

  size_ = graph.max();

  weights_.reserve( (norb1 - nele1)*nele1 + (norb2 - nele2)*nele2 + (norb3 - nele3)*nele3 );
  offsets_.reserve( nele_ );

  for (int j = 0; j < nele_; ++j) {
    int i = 0;
    while ( (graph(i+1,j+1) < 0) || (graph(i,j) < 0) ) ++i;

    offsets_.push_back( weights_.size() - i );

    for ( ; i < norb_; ++i) {
      if (graph(i+1,j+1) < 0 || graph(i,j) < 0) break;
      weights_.push_back( max(0, graph(i,j+1)) );
    }
  }

  // Lexical ordering done, now fill in all the strings
  strings_ = vector<bitset<nbiti__>>(size_, bitset<nbit__>(0ul));
  auto istring = strings_.begin();

  vector<int> holes(norb1);
  iota(holes.begin(), holes.end(), 0);
  do {
    vector<int> active(norb2);
    iota(active.begin(), active.end(), norb1);
    do {
      vector<int> particles(norb3);
      iota(particles.begin(), particles.end(), norb1 + norb2);
      do {
        for (int i = 0; i != nele1; ++i) istring->set(holes[i]);
        for (int i = 0; i != nele2; ++i) istring->set(active[i]);
        for (int i = 0; i != nele3; ++i) istring->set(particles[i]);

        ++istring;
      } while (boost::next_combination(particles.begin(), particles.begin() + nele3, particles.end()));
    } while (boost::next_combination(active.begin(), active.begin() + nele2, active.end()));
  } while (boost::next_combination(holes.begin(), holes.begin() + nele1, holes.end()));
}
