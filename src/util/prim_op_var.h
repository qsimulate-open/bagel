//
// BAGEL - Parallel electron correlation program.
// Filename: prim_op_var.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __SRC_UTIL_PRIM_OP_VAR_H
#define __SRC_UTIL_PRIM_OP_VAR_H

#include <src/util/prim_op.h>
#include <src/util/constants.h>

namespace bagel {

template<typename DataType>
void sort_indices(const std::array<int,4>& o, const double a, const double b, const DataType* in, DataType* out, const std::array<int,4>& d) {
  if (o[0] == 0 && o[1] == 1 && o[2] == 2 && o[3] == 3 && std::abs(a-1.0) < numerical_zero__ && std::abs(b) < numerical_zero__)
    sort_indices<0,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3]);
  else if (o[0] == 1 && o[1] == 0 && o[2] == 2 && o[3] == 3 && std::abs(a-1.0) < numerical_zero__ && std::abs(b) < numerical_zero__)
    sort_indices<1,0,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3]);
  else if (o[0] == 0 && o[1] == 1 && o[2] == 3 && o[3] == 2 && std::abs(a-1.0) < numerical_zero__ && std::abs(b) < numerical_zero__)
    sort_indices<0,1,3,2,0,1,1,1>(in, out, d[0], d[1], d[2], d[3]);
  else if (o[0] == 1 && o[1] == 0 && o[2] == 3 && o[3] == 2 && std::abs(a-1.0) < numerical_zero__ && std::abs(b) < numerical_zero__)
    sort_indices<1,0,3,2,0,1,1,1>(in, out, d[0], d[1], d[2], d[3]);
  else if (o[0] == 0 && o[1] == 1 && o[2] == 2 && o[3] == 3 && std::abs(a+1.0) < numerical_zero__ && std::abs(b) < numerical_zero__)
    sort_indices<0,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3]);
  else if (o[0] == 1 && o[1] == 0 && o[2] == 2 && o[3] == 3 && std::abs(a+1.0) < numerical_zero__ && std::abs(b) < numerical_zero__)
    sort_indices<1,0,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3]);
  else if (o[0] == 0 && o[1] == 1 && o[2] == 3 && o[3] == 2 && std::abs(a+1.0) < numerical_zero__ && std::abs(b) < numerical_zero__)
    sort_indices<0,1,3,2,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3]);
  else if (o[0] == 1 && o[1] == 0 && o[2] == 3 && o[3] == 2 && std::abs(a+1.0) < numerical_zero__ && std::abs(b) < numerical_zero__)
    sort_indices<1,0,3,2,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3]);
  else
    throw std::logic_error("This case has not been implemented in prim_op_var.h");
}

}

#endif
