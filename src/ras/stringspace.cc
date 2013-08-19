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

RASStringSpace(const int nele, const int max_holes, const int max_particles, const int norb1, const int norb2, const int norb3) :
  norb_{norb1, norb2, norb3}, nele_(nele), max_holes_(max_holes), max_particles_(max_particles)
{
  lexical_.reserve( (max_holes_ + 1) * (max_particles_ + 1) );

  for (int nparticles = 0; nparticles <= max_particles_; ++nparticles) {
    for (int nholes = 0; nholes <= max_holes_; ++nholes) {
      const int nele1 = norb1 - nholes;
      const int nele2 = nele_ - (nparticles - nholes);
      const int nele3 = nparticles;

      lexical_.push_back(make_shared<RASLexical>(nele1, norb_[0], nele2, norb_[1], nele3, norb_[2], strings_.size()));

      // Generate all possible strings
      vector<bitset<nbit__>> tmpstrings(strings->size(), bitset<nbit__>(0));

      vector<int> holes(norb_[0]);
      iota(holes.begin(), holes.end(), 0);

      auto istring = tmpstrings.begin();
      do {
        vector<int> active(norb_[1]);
        iota(active.begin(), active.end(), norb_[0]);
        do {
          vector<int> particles(norb_[2]);
          iota(particles.begin(), particles.end(), norb_[0] + norb_[1]);
          do {
            for (int i = 0; i != nele1; ++i) istring->set(holes[i]);
            for (int i = 0; i != nele2; ++i) istring->set(active[i]);
            for (int i = 0; i != nele3; ++i) istring->set(particles[i]);

            ++istring;
          } while (boost::next_combination(particles.begin(), particles.begin() + nele3, particles.end()));
        } while (boost::next_combination(active.begin(), active.begin() + nele2, active.end()));
      } while (boost::next_combination(holes.begin(), holes.begin() + nele1, holes.end()));
      // boost::next_combination should produce combinations in the same ordering as the lexical ordering,
      // but it could be worthwhile to double check this at some point

      strings_.insert(strings_.end(), tmpstrings.begin(), tmpstrings.end());
    }
  }
}

