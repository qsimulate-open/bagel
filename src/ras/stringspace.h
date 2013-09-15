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

namespace bagel {

// Contains all the strings and information for lexical ordering for one particular graph (set of strings)
class StringSpace {
  protected:
    std::array<const std::pair<const int, const int>, 3> ras_;

    std::vector<std::bitset<nbit__>> strings_;

    // Lexical ordering
    std::vector<size_t> weights_;
    std::vector<size_t> offsets_;

    const int norb_;
    const int nele_;

    const size_t offset_;

    struct RASGraph {
      const int ndim_;
      const int mdim_;
      std::unique_ptr<size_t[]> data_;

      RASGraph(const size_t n, const size_t m) : ndim_(n), mdim_(m) {
        data_ = std::unique_ptr<size_t[]>(new size_t[n*m]);
        std::fill_n(data_.get(), ndim_*mdim_, 0);
      }

      size_t& operator()(const int i, const int j) { return data_[j*ndim_ + i]; }
      const size_t max() const { return *std::max_element(data_.get(), data_.get() + ndim_ * mdim_); }
    };

  public:
    StringSpace(const int nele1, const int norb1, const int nele2, const int norb2, const int nele3, const int norb3, const size_t offset = 0);

    const int nele() const { return nele_; }
    const int norb() const { return norb_; }

    template <int subspace> const std::pair<const int, const int> ras() { return std::get<subspace>(ras_); }

    const int nholes() const { return ras_[0].second - ras_[0].first; }
    const int nele2() const { return nele_ - ras_[0].first - ras_[2].first; }
    const int nparticles() const { return ras_[2].first; }

    const size_t size() const { return strings_.size(); }
    const size_t offset() const { return offset_; }

    const std::vector<std::bitset<nbit__>>& strings() const { return strings_; }
    const std::bitset<nbit__> strings(const size_t i) const { return strings_[i]; }

    std::vector<std::bitset<nbit__>>::iterator begin() { return strings_.begin(); }
    std::vector<std::bitset<nbit__>>::iterator end() { return strings_.end(); }
    std::vector<std::bitset<nbit__>>::const_iterator begin() const { return strings_.cbegin(); }
    std::vector<std::bitset<nbit__>>::const_iterator end() const { return strings_.cend(); }

    // Assumes bit is within this graph
    template <int off = 1>
    size_t lexical(const std::bitset<nbit__>& bit) const {
      size_t out = ( off == 1 ? offset_ : 0 );
      int nele = 0;
      for (int i = 0; i != norb_; ++i)
        if (bit[i]) { out += weights_[offsets_[nele++] + i]; }
      return out;
    }
};

}

#endif
