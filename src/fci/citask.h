//
// BAGEL - Parallel electron correlation program.
// Filename: citask.h
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

#ifndef __SRC_FCI_CITASK_H
#define __SRC_FCI_CITASK_H

namespace bagel {

// ColumnTask computes two columns for the sake of load balancing
class CITask {
  protected:
    using SD = std::pair<std::bitset<nbit__>, std::bitset<nbit__>>;

  protected:
    std::vector<SD>* basis_;
    const int norb_;
    const size_t col1_;
    const size_t col2_;
    double* const dest1_;
    double* const dest2_;


    CITask(std::vector<SD>* b, const int norb, const size_t c1, double* d1, const size_t c2, double* d2) :
      basis_(b), norb_(norb), col1_(c1), col2_(c2), dest1_(d1), dest2_(d2) { }

    std::vector<int> to_vector(const std::bitset<nbit__> bit) const {
      std::vector<int> out;
      for(int i = 0; i < norb_; ++i)
        if (bit[i]) out.push_back(i);
      return out;
    }

    int sign(std::bitset<nbit__> bit, int i, int j) {
      // masking irrelevant bits
      int min, max;
      std::tie(min,max) = std::minmax(i,j);
      bit &= ~((1ull << (min+1)) - 1ull);
      bit &= (1ull << max) - 1ull;
      return 1 - ((bit.count() & 1) << 1);
    }

    virtual double matrix_element(const SD& bra, const SD& ket) = 0;

  public:
    void compute() {
      double* odata = dest1_;
      const SD bra1 = basis_->at(col1_);
      for_each(basis_->begin() + col1_, basis_->end(), [this, &bra1, &odata] (const SD& ket) { *odata = matrix_element(bra1, ket); ++odata; });

      if (col1_ != col2_) {
        odata = dest2_;
        const SD bra2 = basis_->at(col2_);
        for_each(basis_->begin() + col2_, basis_->end(), [this, &bra2, &odata] (const SD& ket) { *odata = matrix_element(bra2, ket); ++odata; });
      }
    }
};

}

#endif
