//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: scaledata.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_RYSINT_SCALEDATA_H
#define __SRC_RYSINT_SCALEDATA_H

#include <cassert>

namespace bagel {

template<int rank_, int worksize, typename DataType>
void scaledata(DataType* out, const DataType* a, const DataType c, const DataType* in) {
  static_assert(worksize % rank_ == 0, "worksize and rank_ inconsistent");
  DataType ca[rank_];
  for (int i = 0; i != rank_; ++i)
    ca[i] = c * a[i];

  for (int of = 0; of != worksize; of += rank_)
    for (int i = 0; i != rank_; ++i)
      out[of+i] = in[of+i] * ca[i];
}

}

#endif

