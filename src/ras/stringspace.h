//
// BAGEL - Parallel electron correlation program.
// Filename: ras/stringspace.h
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


#ifndef BAGEL_RAS_STRINGSPACE_H
#define BAGEL_RAS_STRINGSPACE_H

#include <src/ras/lexical.h>

namespace bagel {

// Different RAS strings are internally stored in a matrix form
//  where the location (i,j) is for i holes in RAS-I and j particles in RAS-III
class RASStringSpace {
  protected:
    std::vector<std::shared_ptr<RASLexical>> lexical_;
    std::vector<std::bitset<nbit__>> strings_;

    const std::array<const int, 3> norb_;

    const int nele_;

    const int max_holes_;
    const int max_particles_;

    const int subspace(const int nholes, const int nparticles) const { return nholes + nparticles * max_holes_; }
    const int subspace(const std::bitset<nbit__> string) const { return nholes(string) + nparticles(string) * max_holes_; }

    const int nholes(const std::bitset<nbit__> string) const { return ( (string & std::bitset<nbit__>((1ul << norb_[0]) - 1)).count() ); }
    const int nparticles(const std::bitset<nbit__> string) const { return ( (string & std::bitset<nbit__>(((1ul << norb_[2]) - 1) << (norb_[0] + norb_[1]))).count() ); }

  public:
    RASStringSpace(const int nele, const int max_holes, const int max_particles, const int norb1, const int norb2, const int norb3);

    size_t lexical(const int nholes, const int nparticles, const std::bitset<nbit__> string) const { return lexical_[subspace(nholes,nparticles)]->address(string); }
    size_t lexical(const std::bitset<nbit__> string) const { return lexical_[subspace(string)]->address(string); }
};

}

#endif
