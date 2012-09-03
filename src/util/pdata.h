//
// BAGEL - Parallel electron correlation program.
// Filename: pdata.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __src_util_pdata_h
#define __src_util_pdata_h

#include <cstddef>
#include <complex>
#include <algorithm>

namespace bagel {

class PData {
  protected:
    std::complex<double>* data_;
    const int length_;

  public:
    PData(int i) : length_(i) { data_ = new std::complex<double>[i];
                                const std::complex<double> zero(0.0, 0.0);
                                std::fill(data_, data_ + i, zero);
                              };
    ~PData() { delete[] data_; };

    std::complex<double>& operator[] (int i) { assert(i < length_ && i >= 0); return data_[i]; };

    std::complex<double>* front() { return data_; };
    const std::complex<double>* cfront() const { return data_; };
    std::complex<double>* pointer(size_t j) { return data_ + j; };

};

}

#endif
