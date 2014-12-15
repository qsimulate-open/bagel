//
// BAGEL - Parallel electron correlation program.
// Filename: block_key.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#ifndef BAGEL_ASD_DMRG_BLOCK_KEY_H
#define BAGEL_ASD_DMRG_BLOCK_KEY_H

namespace bagel {

/// Convenient key for organizing blocks in ASD. Only holds nelea and neleb.
/// Useful in maps.
struct BlockKey {
  int nelea;
  int neleb;

  BlockKey(const int na, const int nb) : nelea(na), neleb(nb) {}

  bool operator<(const BlockKey& o) const {
    if (nelea+neleb==o.nelea+o.neleb)
      return std::make_pair(nelea,neleb) < std::make_pair(o.nelea, o.neleb);
    else
      return nelea+neleb < o.nelea+o.neleb;
  }
  bool operator==(const BlockKey& o) const { return std::make_pair(nelea, neleb) == std::make_pair(o.nelea, o.neleb); }
};

/// Extends BlockKey to also include the number of states.
/// Use primarily as a store of information, NOT to organize data such as in a map
struct BlockInfo : public BlockKey {
  int nstates;
  BlockInfo(const int na, const int nb, const int ns) : BlockKey(na,nb), nstates(ns) {}
  BlockKey key() const { return BlockKey(nelea,neleb); }

  bool operator==(const BlockInfo& o) const { return (key()==o.key() && nstates==o.nstates); }
};

}

namespace std {
template <> struct hash<bagel::BlockKey> {
  typedef bagel::BlockKey argument_type;
  typedef std::size_t result_type;

  result_type operator()(const argument_type& k) const { return k.nelea + (k.neleb << 16); }
};

}

#endif
