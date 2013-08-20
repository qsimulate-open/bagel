//
// BAGEL - Parallel electron correlation program.
// Filename: ras/block.h
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


#ifndef BAGEL_RAS_BLOCK_H
#define BAGEL_RAS_BLOCK_H

namespace bagel {

// Contains all the data for a sub block of the CI coefficient matrix
class RASBlock {
  protected:
    std::shared_ptr<const StringSpace> astrings_;
    std::shared_ptr<const StringSpace> bstrings_;

    std::unique_ptr<double[]> data_;

    const int size_;

  public:
    RASBlock(std::shared_ptr<const StringSpace> astrings, std::shared_ptr<const StringSpace> bstrings) :
      astrings_(astrings), bstrings_(bstrings), size_(astrings->size() * bstrings->size())
    {
      data_ = std::unique_ptr<double[]>(new double[size_]);
      std::fill_n(data_.get(), size_, 0.0);
    }

    const int size() const { return size_; }

    double* data() { return data_.get(); }
    const double* data() const { return data_.get(); }

    template <int spin> std::shared_ptr<const StringSpace> strings() const { return ( spin == 0 ? astrings_ : bstrings_ ); }
};

}

#endif
