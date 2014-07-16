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

// Only orders based on nelea and neleb. nstates is just extra info
struct BlockKey {
  int nelea;
  int neleb;
  int nstates;

  BlockKey(const int na, const int nb, const int ns = 0) : nelea(na), neleb(nb), nstates(ns) {}

  bool operator<(const BlockKey& o) const { return std::make_pair(nelea, neleb) < std::make_pair(o.nelea, o.neleb); }
  bool operator==(const BlockKey& o) const { return std::make_pair(nelea, neleb) == std::make_pair(o.nelea, o.neleb); }
};

}

#endif
