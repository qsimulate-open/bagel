//
// BAGEL - Parallel electron correlation program.
// Filename: scaledata.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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


#ifndef __SRC_RYSINT_SCALEDATA_H
#define __SRC_RYSINT_SCALEDATA_H

#include <cassert>

namespace bagel {

template<int rank_, int worksize>
void scaledata(double* out, const double* a, const double c, const double* in) {
  static_assert(worksize % rank_ == 0, "worksize and rank_ inconsistent");
  double ca[rank_];
  for (int i = 0; i != rank_; ++i)
    ca[i] = c * a[i];

  for (int of = 0; of != worksize; of += rank_)
    for (int i = 0; i != rank_; ++i)
      out[of+i] = in[of+i] * ca[i];
}

}

#endif

