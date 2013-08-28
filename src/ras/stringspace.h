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
    std::vector<double> weights_;
    std::vector<int> offsets_;

    const int norb_;
    const int nele_;

    const int offset_;
    int size_;

    struct RASGraph {
      const int ndim_;
      const int mdim_;
      std::unique_ptr<int[]> data_;

      RASGraph(const int n, const int m) : ndim_(n), mdim_(m) {
        data_ = std::unique_ptr<int[]>(new int[n*m]);
        std::fill_n(data_.get(), ndim_*mdim_, -1);
      }

      int& operator()(const int i, const int j) { return data_[j*ndim_ + i]; }
      const int max() const { return *std::max_element(data_.get(), data_.get() + ndim_ * mdim_); }
    };

  public:
    StringSpace(const int nele1, const int norb1, const int nele2, const int norb2, const int nele3, const int norb3, const int offset = 0);

    const int nele() const { return nele_; }
    const int norb() const { return norb_; }

    template <int subspace> const std::pair<const int, const int> ras() { return std::get<subspace>(ras_); }

    const int nholes() const { return ras_[0].second - ras_[0].first; }
    const int nparticles() const { return ras_[2].first; }

    const int size() const { return size_; }
    const int offset() const { return offset_; }

    const std::vector<std::bitset<nbit__>>& strings() const { return strings_; }
    const std::bitset<nbit__> strings(const int i) const { return strings_[i]; }

    std::vector<std::bitset<nbit__>>::iterator begin() { return strings_.begin(); }
    std::vector<std::bitset<nbit__>>::iterator end() { return strings_.end(); }
    std::vector<std::bitset<nbit__>>::const_iterator cbegin() const { return strings_.cbegin(); }
    std::vector<std::bitset<nbit__>>::const_iterator cend() const { return strings_.cend(); }

    // Assumes bit is within this graph
    template <int off = 1>
    unsigned int lexical(std::bitset<nbit__> bit) const {
      unsigned int out = 0;
      int nele = 0;
      for (int i = 0; i != norb_; ++i)
        if (bit[i]) { out += weights_[offsets_[nele] + i]; ++nele; }
      return out + off * offset_;
    }
};

}

#endif
