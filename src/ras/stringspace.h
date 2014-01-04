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

#include <vector>
#include <array>
#include <memory>
#include <bitset>
#include <algorithm>

#include <src/util/constants.h>
#include <src/math/comb.h>
#include <src/parallel/staticdist.h>
#include <src/parallel/mpi_interface.h>

namespace bagel {

// Contains all the strings and information for lexical ordering for one particular graph (set of strings)
//   comprised of three subgraphs (one each for RASI, RASII, RASIII)
class RASGraph {
  protected:
    const size_t nele_;
    const size_t norb_;

    size_t size_;

    std::unique_ptr<size_t[]> weights_;

  public:
    RASGraph(const size_t nele, const size_t norb);

    size_t& weight(const size_t i, const size_t j) { assert(nele_*norb_ > 0); return weights_[i + j*norb_]; }
    const size_t& weight(const size_t i, const size_t j) const { assert(nele_*norb_ > 0); return weights_[i + j*norb_]; }

    const size_t size() const { return size_; }

    const size_t lexical(const int& start, const int& fence, const std::bitset<nbit__>& abit) const {
      size_t out = 0;

      int k = 0;
      for (int i = start; i < fence; ++i)
        if (abit[i]) out += weight(i-start,k++);
      return out;
    }
};

class StringSpace {
  protected:
    std::array<const std::pair<const int, const int>, 3> ras_;

    std::vector<std::bitset<nbit__>> strings_;

    std::array<std::shared_ptr<RASGraph>, 3> graphs_;

    StaticDist dist_;

    const int norb_;
    const int nele_;

    const size_t offset_;

  public:
    StringSpace(const int nele1, const int norb1, const int nele2, const int norb2, const int nele3, const int norb3, const size_t offset = 0);

    const int nele() const { return nele_; }
    const int norb() const { return norb_; }

    template <int subspace> const std::pair<const int, const int> ras() const { return std::get<subspace>(ras_); }

    const int nholes() const { return ras_[0].second - ras_[0].first; }
    const int nele2() const { return nele_ - ras_[0].first - ras_[2].first; }
    const int nparticles() const { return ras_[2].first; }

    const size_t size() const { return strings_.size(); }
    const size_t offset() const { return offset_; }

    const size_t size1() const { return graphs_[0]->size(); }
    const size_t size2() const { return graphs_[1]->size(); }
    const size_t size3() const { return graphs_[2]->size(); }

    const std::vector<std::bitset<nbit__>>& strings() const { return strings_; }
    const std::bitset<nbit__> strings(const size_t i) const { return strings_[i]; }

    std::vector<std::bitset<nbit__>>::iterator begin() { return strings_.begin(); }
    std::vector<std::bitset<nbit__>>::iterator end() { return strings_.end(); }
    std::vector<std::bitset<nbit__>>::const_iterator begin() const { return strings_.cbegin(); }
    std::vector<std::bitset<nbit__>>::const_iterator end() const { return strings_.cend(); }

    const StaticDist& dist() const { return dist_; }

    // Assumes bit is within this graph
    template <int off = 1>
    size_t lexical(const std::bitset<nbit__>& bit) const {
      size_t out = ( off == 1 ? offset_ : 0 );

      const size_t r1 = ras_[0].second;
      const size_t r2 = ras_[1].second;
      const size_t r3 = ras_[2].second;

      const size_t n2 = graphs_[1]->size();
      const size_t n1 = graphs_[0]->size();

      out += graphs_[1]->lexical(r1, r1+r2, bit);
      out += n2 * graphs_[0]->lexical(0, r1, bit);
      out += n2 * n1 * graphs_[2]->lexical(r1+r2, r1+r2+r3, bit);

      return out;
    }
};

}

#endif
