//
// BAGEL - Parallel electron correlation program.
// Filename: oscontraction_impl.hpp
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifdef OSINTEGRAL_HEADERS

#ifndef __SRC_INTEGRAL_OS_OSCONTRACTION_IMPL_HPP
#define __SRC_INTEGRAL_OS_OSCONTRACTION_IMPL_HPP

namespace bagel {

template <typename DataType, Int_t IntType>
void OSIntegral<DataType, IntType>::perform_contraction(const int asize, const DataType* prim, const int pdim0, const int pdim1, DataType* cont,
                                const std::vector<std::vector<double>>& coeff0, const std::vector<std::pair<int, int>>& ranges0, const int cdim0,
                                const std::vector<std::vector<double>>& coeff1, const std::vector<std::pair<int, int>>& ranges1, const int cdim1) {
  // transformation of index1
  const int worksize = pdim1 * asize;
  DataType* const work = stack_->template get<DataType>(worksize);
  std::fill_n(cont, asize*cdim0*cdim1, 0.0);

  for (int i = 0; i != cdim0; ++i) {
    const int begin0 = ranges0[i].first;
    const int end0   = ranges0[i].second;
    std::fill_n(work, worksize, 0.0);
    for (int j = begin0; j != end0; ++j)
      for (int n=0; n!=worksize; n++) work[n] += coeff0[i][j]*prim[j*worksize+n];
    for (int k = 0; k != cdim1; ++k, cont += asize) {
      const int begin1 = ranges1[k].first;
      const int end1   = ranges1[k].second;
      for (int j = begin1; j != end1; ++j) {
        for (int n=0; n!=asize; n++) cont[n] += coeff1[k][j]*work[j*asize+n];
      }
    }
  }
  stack_->release(worksize, work);
}

}

#endif
#endif
