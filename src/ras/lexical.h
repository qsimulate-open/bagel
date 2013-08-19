//
// BAGEL - Parallel electron correlation program.
// Filename: ras/lexical.h
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


#ifndef BAGEL_RAS_LEXICAL_H
#define BAGEL_RAS_LEXICAL_H

namespace bagel {

class RASLexical {
  protected:
    std::array<const std::pair<const int, const int>, 3> ras_;
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
      const int max() const { return *std::max_element(data_.get(), data_.get() * ndim_ * mdim_); }
    };

  public:
    RASLexical(const int nele1, const int norb1, const int nele2, const int norb2, const int nele3, const int norb3, const int offset = 0);

    const int size() const { return size_; }

    // Assumes bit is within this graph
    unsigned int address(std::bitset<nbit__> bit) const {
      unsigned int out = 0;
      int nele = 0;
      for (int i = 0; i != norb_; ++i)
        if (bit[i]) { out += weights_[offsets_[nele] + i]; ++nele; }
      return out + offset_;
    }
};

}

#endif
