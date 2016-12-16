//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: citask.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __SRC_FCI_CITASK_H
#define __SRC_FCI_CITASK_H

namespace bagel {

// CITask computes two columns for the sake of load balancing
template <class Derived>
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

    double matrix_element(const SD& bra, const SD& ket) { return static_cast<Derived*>(this)->matrix_element_impl(bra,ket); }

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
