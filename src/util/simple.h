//
// BAGEL - Parallel electron correlation program.
// Filename: simple.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __SRC_UTIL_SIMPLE_H
#define __SRC_UTIL_SIMPLE_H

#include <algorithm>

namespace bagel {

struct CopyBlockTask {
  private:
    const double* const a_;
    const size_t astride_;
    double* const b_;
    const size_t bstride_;
    const size_t n_;
    const size_t m_;
  public:
    CopyBlockTask(const double* const a, const size_t& ast, double* const b, const size_t& bst, const size_t& n, const size_t& m)
      : a_(a), astride_(ast), b_(b), bstride_(bst), n_(n), m_(m) {}

    void compute() {
      for (size_t j = 0; j != m_; ++j)
        std::copy_n(a_+j*astride_, n_, b_+j*bstride_);
    }
};

}

#endif
