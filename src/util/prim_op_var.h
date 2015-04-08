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


template<typename DataType>
void sort_indices(const std::array<int,6>& o, const double a, const double b, const DataType* in, DataType* out, const std::array<int,6>& d) {
  size_t tag = 0ull;
  for (int i = 0; i != 6; ++i)
    tag = (tag << 4) + o[i];

  if (std::abs(a-1.0) < numerical_zero__ && std::abs(b) < numerical_zero__) {
    switch (tag) {
    case 74565ull :
      sort_indices<0,1,2,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 136005ull :
      sort_indices<0,2,1,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1057605ull :
      sort_indices<1,0,2,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1180485ull :
      sort_indices<1,2,0,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2102085ull :
      sort_indices<2,0,1,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2163525ull :
      sort_indices<2,1,0,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 74580ull :
      sort_indices<0,1,2,3,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 136020ull :
      sort_indices<0,2,1,3,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1057620ull :
      sort_indices<1,0,2,3,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1180500ull :
      sort_indices<1,2,0,3,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2102100ull :
      sort_indices<2,0,1,3,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2163540ull :
      sort_indices<2,1,0,3,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 74805ull :
      sort_indices<0,1,2,4,3,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 136245ull :
      sort_indices<0,2,1,4,3,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1057845ull :
      sort_indices<1,0,2,4,3,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1180725ull :
      sort_indices<1,2,0,4,3,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2102325ull :
      sort_indices<2,0,1,4,3,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2163765ull :
      sort_indices<2,1,0,4,3,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 74835ull :
      sort_indices<0,1,2,4,5,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 136275ull :
      sort_indices<0,2,1,4,5,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1057875ull :
      sort_indices<1,0,2,4,5,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1180755ull :
      sort_indices<1,2,0,4,5,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2102355ull :
      sort_indices<2,0,1,4,5,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2163795ull :
      sort_indices<2,1,0,4,5,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 75060ull :
      sort_indices<0,1,2,5,3,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 136500ull :
      sort_indices<0,2,1,5,3,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1058100ull :
      sort_indices<1,0,2,5,3,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1180980ull :
      sort_indices<1,2,0,5,3,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2102580ull :
      sort_indices<2,0,1,5,3,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2164020ull :
      sort_indices<2,1,0,5,3,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 75075ull :
      sort_indices<0,1,2,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 136515ull :
      sort_indices<0,2,1,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1058115ull :
      sort_indices<1,0,2,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1180995ull :
      sort_indices<1,2,0,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2102595ull :
      sort_indices<2,0,1,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2164035ull :
      sort_indices<2,1,0,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    }
  } else if (std::abs(a+1.0) < numerical_zero__ && std::abs(b) < numerical_zero__) {
    switch (tag) {
    case 74565ull :
      sort_indices<0,1,2,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 136005ull :
      sort_indices<0,2,1,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1057605ull :
      sort_indices<1,0,2,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1180485ull :
      sort_indices<1,2,0,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2102085ull :
      sort_indices<2,0,1,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2163525ull :
      sort_indices<2,1,0,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 74580ull :
      sort_indices<0,1,2,3,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 136020ull :
      sort_indices<0,2,1,3,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1057620ull :
      sort_indices<1,0,2,3,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1180500ull :
      sort_indices<1,2,0,3,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2102100ull :
      sort_indices<2,0,1,3,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2163540ull :
      sort_indices<2,1,0,3,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 74805ull :
      sort_indices<0,1,2,4,3,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 136245ull :
      sort_indices<0,2,1,4,3,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1057845ull :
      sort_indices<1,0,2,4,3,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1180725ull :
      sort_indices<1,2,0,4,3,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2102325ull :
      sort_indices<2,0,1,4,3,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2163765ull :
      sort_indices<2,1,0,4,3,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 74835ull :
      sort_indices<0,1,2,4,5,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 136275ull :
      sort_indices<0,2,1,4,5,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1057875ull :
      sort_indices<1,0,2,4,5,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1180755ull :
      sort_indices<1,2,0,4,5,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2102355ull :
      sort_indices<2,0,1,4,5,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2163795ull :
      sort_indices<2,1,0,4,5,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 75060ull :
      sort_indices<0,1,2,5,3,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 136500ull :
      sort_indices<0,2,1,5,3,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1058100ull :
      sort_indices<1,0,2,5,3,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1180980ull :
      sort_indices<1,2,0,5,3,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2102580ull :
      sort_indices<2,0,1,5,3,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2164020ull :
      sort_indices<2,1,0,5,3,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 75075ull :
      sort_indices<0,1,2,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 136515ull :
      sort_indices<0,2,1,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1058115ull :
      sort_indices<1,0,2,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 1180995ull :
      sort_indices<1,2,0,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2102595ull :
      sort_indices<2,0,1,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    case 2164035ull :
      sort_indices<2,1,0,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
    }
  } else {
    throw std::logic_error("This case has not been implemented in prim_op_var.h");
  }
}


template<typename DataType>
void sort_indices(const std::array<int,8>& o, const double a, const double b, const DataType* in, DataType* out, const std::array<int,8>& d) {
  size_t tag = 0lu;
  for (int i = 0; i != 8; ++i)
    tag = (tag << 4) + o[i];

  if (std::abs(a-1.0) < numerical_zero__ && std::abs(b) < numerical_zero__) {
    switch (tag) {
    case 19088743ull :
      sort_indices<0,1,2,3,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20071783ull :
      sort_indices<0,1,3,2,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34817383ull :
      sort_indices<0,2,1,3,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36783463ull :
      sort_indices<0,2,3,1,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51529063ull :
      sort_indices<0,3,1,2,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52512103ull :
      sort_indices<0,3,2,1,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270746983ull :
      sort_indices<1,0,2,3,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271730023ull :
      sort_indices<1,0,3,2,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302204263ull :
      sort_indices<1,2,0,3,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305153383ull :
      sort_indices<1,2,3,0,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318915943ull :
      sort_indices<1,3,0,2,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320882023ull :
      sort_indices<1,3,2,0,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538133863ull :
      sort_indices<2,0,1,3,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540099943ull :
      sort_indices<2,0,3,1,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553862503ull :
      sort_indices<2,1,0,3,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556811623ull :
      sort_indices<2,1,3,0,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587285863ull :
      sort_indices<2,3,0,1,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588268903ull :
      sort_indices<2,3,1,0,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806503783ull :
      sort_indices<3,0,1,2,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807486823ull :
      sort_indices<3,0,2,1,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822232423ull :
      sort_indices<3,1,0,2,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824198503ull :
      sort_indices<3,1,2,0,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838944103ull :
      sort_indices<3,2,0,1,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839927143ull :
      sort_indices<3,2,1,0,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19088758ull :
      sort_indices<0,1,2,3,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20071798ull :
      sort_indices<0,1,3,2,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34817398ull :
      sort_indices<0,2,1,3,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36783478ull :
      sort_indices<0,2,3,1,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51529078ull :
      sort_indices<0,3,1,2,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52512118ull :
      sort_indices<0,3,2,1,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270746998ull :
      sort_indices<1,0,2,3,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271730038ull :
      sort_indices<1,0,3,2,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302204278ull :
      sort_indices<1,2,0,3,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305153398ull :
      sort_indices<1,2,3,0,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318915958ull :
      sort_indices<1,3,0,2,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320882038ull :
      sort_indices<1,3,2,0,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538133878ull :
      sort_indices<2,0,1,3,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540099958ull :
      sort_indices<2,0,3,1,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553862518ull :
      sort_indices<2,1,0,3,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556811638ull :
      sort_indices<2,1,3,0,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587285878ull :
      sort_indices<2,3,0,1,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588268918ull :
      sort_indices<2,3,1,0,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806503798ull :
      sort_indices<3,0,1,2,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807486838ull :
      sort_indices<3,0,2,1,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822232438ull :
      sort_indices<3,1,0,2,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824198518ull :
      sort_indices<3,1,2,0,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838944118ull :
      sort_indices<3,2,0,1,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839927158ull :
      sort_indices<3,2,1,0,4,5,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19088983ull :
      sort_indices<0,1,2,3,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20072023ull :
      sort_indices<0,1,3,2,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34817623ull :
      sort_indices<0,2,1,3,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36783703ull :
      sort_indices<0,2,3,1,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51529303ull :
      sort_indices<0,3,1,2,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52512343ull :
      sort_indices<0,3,2,1,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270747223ull :
      sort_indices<1,0,2,3,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271730263ull :
      sort_indices<1,0,3,2,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302204503ull :
      sort_indices<1,2,0,3,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305153623ull :
      sort_indices<1,2,3,0,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318916183ull :
      sort_indices<1,3,0,2,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320882263ull :
      sort_indices<1,3,2,0,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538134103ull :
      sort_indices<2,0,1,3,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540100183ull :
      sort_indices<2,0,3,1,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553862743ull :
      sort_indices<2,1,0,3,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556811863ull :
      sort_indices<2,1,3,0,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587286103ull :
      sort_indices<2,3,0,1,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588269143ull :
      sort_indices<2,3,1,0,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806504023ull :
      sort_indices<3,0,1,2,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807487063ull :
      sort_indices<3,0,2,1,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822232663ull :
      sort_indices<3,1,0,2,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824198743ull :
      sort_indices<3,1,2,0,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838944343ull :
      sort_indices<3,2,0,1,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839927383ull :
      sort_indices<3,2,1,0,4,6,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19089013ull :
      sort_indices<0,1,2,3,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20072053ull :
      sort_indices<0,1,3,2,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34817653ull :
      sort_indices<0,2,1,3,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36783733ull :
      sort_indices<0,2,3,1,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51529333ull :
      sort_indices<0,3,1,2,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52512373ull :
      sort_indices<0,3,2,1,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270747253ull :
      sort_indices<1,0,2,3,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271730293ull :
      sort_indices<1,0,3,2,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302204533ull :
      sort_indices<1,2,0,3,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305153653ull :
      sort_indices<1,2,3,0,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318916213ull :
      sort_indices<1,3,0,2,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320882293ull :
      sort_indices<1,3,2,0,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538134133ull :
      sort_indices<2,0,1,3,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540100213ull :
      sort_indices<2,0,3,1,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553862773ull :
      sort_indices<2,1,0,3,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556811893ull :
      sort_indices<2,1,3,0,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587286133ull :
      sort_indices<2,3,0,1,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588269173ull :
      sort_indices<2,3,1,0,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806504053ull :
      sort_indices<3,0,1,2,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807487093ull :
      sort_indices<3,0,2,1,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822232693ull :
      sort_indices<3,1,0,2,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824198773ull :
      sort_indices<3,1,2,0,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838944373ull :
      sort_indices<3,2,0,1,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839927413ull :
      sort_indices<3,2,1,0,4,6,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19089238ull :
      sort_indices<0,1,2,3,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20072278ull :
      sort_indices<0,1,3,2,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34817878ull :
      sort_indices<0,2,1,3,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36783958ull :
      sort_indices<0,2,3,1,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51529558ull :
      sort_indices<0,3,1,2,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52512598ull :
      sort_indices<0,3,2,1,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270747478ull :
      sort_indices<1,0,2,3,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271730518ull :
      sort_indices<1,0,3,2,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302204758ull :
      sort_indices<1,2,0,3,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305153878ull :
      sort_indices<1,2,3,0,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318916438ull :
      sort_indices<1,3,0,2,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320882518ull :
      sort_indices<1,3,2,0,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538134358ull :
      sort_indices<2,0,1,3,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540100438ull :
      sort_indices<2,0,3,1,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553862998ull :
      sort_indices<2,1,0,3,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556812118ull :
      sort_indices<2,1,3,0,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587286358ull :
      sort_indices<2,3,0,1,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588269398ull :
      sort_indices<2,3,1,0,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806504278ull :
      sort_indices<3,0,1,2,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807487318ull :
      sort_indices<3,0,2,1,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822232918ull :
      sort_indices<3,1,0,2,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824198998ull :
      sort_indices<3,1,2,0,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838944598ull :
      sort_indices<3,2,0,1,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839927638ull :
      sort_indices<3,2,1,0,4,7,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19089253ull :
      sort_indices<0,1,2,3,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20072293ull :
      sort_indices<0,1,3,2,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34817893ull :
      sort_indices<0,2,1,3,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36783973ull :
      sort_indices<0,2,3,1,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51529573ull :
      sort_indices<0,3,1,2,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52512613ull :
      sort_indices<0,3,2,1,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270747493ull :
      sort_indices<1,0,2,3,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271730533ull :
      sort_indices<1,0,3,2,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302204773ull :
      sort_indices<1,2,0,3,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305153893ull :
      sort_indices<1,2,3,0,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318916453ull :
      sort_indices<1,3,0,2,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320882533ull :
      sort_indices<1,3,2,0,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538134373ull :
      sort_indices<2,0,1,3,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540100453ull :
      sort_indices<2,0,3,1,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553863013ull :
      sort_indices<2,1,0,3,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556812133ull :
      sort_indices<2,1,3,0,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587286373ull :
      sort_indices<2,3,0,1,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588269413ull :
      sort_indices<2,3,1,0,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806504293ull :
      sort_indices<3,0,1,2,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807487333ull :
      sort_indices<3,0,2,1,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822232933ull :
      sort_indices<3,1,0,2,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824199013ull :
      sort_indices<3,1,2,0,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838944613ull :
      sort_indices<3,2,0,1,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839927653ull :
      sort_indices<3,2,1,0,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19092583ull :
      sort_indices<0,1,2,3,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20075623ull :
      sort_indices<0,1,3,2,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34821223ull :
      sort_indices<0,2,1,3,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36787303ull :
      sort_indices<0,2,3,1,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51532903ull :
      sort_indices<0,3,1,2,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52515943ull :
      sort_indices<0,3,2,1,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270750823ull :
      sort_indices<1,0,2,3,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271733863ull :
      sort_indices<1,0,3,2,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302208103ull :
      sort_indices<1,2,0,3,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305157223ull :
      sort_indices<1,2,3,0,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318919783ull :
      sort_indices<1,3,0,2,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320885863ull :
      sort_indices<1,3,2,0,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538137703ull :
      sort_indices<2,0,1,3,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540103783ull :
      sort_indices<2,0,3,1,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553866343ull :
      sort_indices<2,1,0,3,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556815463ull :
      sort_indices<2,1,3,0,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587289703ull :
      sort_indices<2,3,0,1,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588272743ull :
      sort_indices<2,3,1,0,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806507623ull :
      sort_indices<3,0,1,2,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807490663ull :
      sort_indices<3,0,2,1,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822236263ull :
      sort_indices<3,1,0,2,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824202343ull :
      sort_indices<3,1,2,0,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838947943ull :
      sort_indices<3,2,0,1,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839930983ull :
      sort_indices<3,2,1,0,5,4,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19092598ull :
      sort_indices<0,1,2,3,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20075638ull :
      sort_indices<0,1,3,2,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34821238ull :
      sort_indices<0,2,1,3,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36787318ull :
      sort_indices<0,2,3,1,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51532918ull :
      sort_indices<0,3,1,2,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52515958ull :
      sort_indices<0,3,2,1,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270750838ull :
      sort_indices<1,0,2,3,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271733878ull :
      sort_indices<1,0,3,2,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302208118ull :
      sort_indices<1,2,0,3,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305157238ull :
      sort_indices<1,2,3,0,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318919798ull :
      sort_indices<1,3,0,2,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320885878ull :
      sort_indices<1,3,2,0,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538137718ull :
      sort_indices<2,0,1,3,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540103798ull :
      sort_indices<2,0,3,1,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553866358ull :
      sort_indices<2,1,0,3,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556815478ull :
      sort_indices<2,1,3,0,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587289718ull :
      sort_indices<2,3,0,1,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588272758ull :
      sort_indices<2,3,1,0,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806507638ull :
      sort_indices<3,0,1,2,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807490678ull :
      sort_indices<3,0,2,1,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822236278ull :
      sort_indices<3,1,0,2,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824202358ull :
      sort_indices<3,1,2,0,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838947958ull :
      sort_indices<3,2,0,1,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839930998ull :
      sort_indices<3,2,1,0,5,4,7,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19093063ull :
      sort_indices<0,1,2,3,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20076103ull :
      sort_indices<0,1,3,2,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34821703ull :
      sort_indices<0,2,1,3,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36787783ull :
      sort_indices<0,2,3,1,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51533383ull :
      sort_indices<0,3,1,2,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52516423ull :
      sort_indices<0,3,2,1,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270751303ull :
      sort_indices<1,0,2,3,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271734343ull :
      sort_indices<1,0,3,2,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302208583ull :
      sort_indices<1,2,0,3,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305157703ull :
      sort_indices<1,2,3,0,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318920263ull :
      sort_indices<1,3,0,2,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320886343ull :
      sort_indices<1,3,2,0,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538138183ull :
      sort_indices<2,0,1,3,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540104263ull :
      sort_indices<2,0,3,1,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553866823ull :
      sort_indices<2,1,0,3,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556815943ull :
      sort_indices<2,1,3,0,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587290183ull :
      sort_indices<2,3,0,1,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588273223ull :
      sort_indices<2,3,1,0,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806508103ull :
      sort_indices<3,0,1,2,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807491143ull :
      sort_indices<3,0,2,1,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822236743ull :
      sort_indices<3,1,0,2,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824202823ull :
      sort_indices<3,1,2,0,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838948423ull :
      sort_indices<3,2,0,1,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839931463ull :
      sort_indices<3,2,1,0,5,6,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19093108ull :
      sort_indices<0,1,2,3,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20076148ull :
      sort_indices<0,1,3,2,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34821748ull :
      sort_indices<0,2,1,3,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36787828ull :
      sort_indices<0,2,3,1,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51533428ull :
      sort_indices<0,3,1,2,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52516468ull :
      sort_indices<0,3,2,1,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270751348ull :
      sort_indices<1,0,2,3,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271734388ull :
      sort_indices<1,0,3,2,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302208628ull :
      sort_indices<1,2,0,3,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305157748ull :
      sort_indices<1,2,3,0,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318920308ull :
      sort_indices<1,3,0,2,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320886388ull :
      sort_indices<1,3,2,0,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538138228ull :
      sort_indices<2,0,1,3,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540104308ull :
      sort_indices<2,0,3,1,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553866868ull :
      sort_indices<2,1,0,3,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556815988ull :
      sort_indices<2,1,3,0,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587290228ull :
      sort_indices<2,3,0,1,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588273268ull :
      sort_indices<2,3,1,0,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806508148ull :
      sort_indices<3,0,1,2,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807491188ull :
      sort_indices<3,0,2,1,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822236788ull :
      sort_indices<3,1,0,2,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824202868ull :
      sort_indices<3,1,2,0,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838948468ull :
      sort_indices<3,2,0,1,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839931508ull :
      sort_indices<3,2,1,0,5,6,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19093318ull :
      sort_indices<0,1,2,3,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20076358ull :
      sort_indices<0,1,3,2,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34821958ull :
      sort_indices<0,2,1,3,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36788038ull :
      sort_indices<0,2,3,1,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51533638ull :
      sort_indices<0,3,1,2,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52516678ull :
      sort_indices<0,3,2,1,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270751558ull :
      sort_indices<1,0,2,3,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271734598ull :
      sort_indices<1,0,3,2,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302208838ull :
      sort_indices<1,2,0,3,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305157958ull :
      sort_indices<1,2,3,0,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318920518ull :
      sort_indices<1,3,0,2,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320886598ull :
      sort_indices<1,3,2,0,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538138438ull :
      sort_indices<2,0,1,3,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540104518ull :
      sort_indices<2,0,3,1,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553867078ull :
      sort_indices<2,1,0,3,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556816198ull :
      sort_indices<2,1,3,0,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587290438ull :
      sort_indices<2,3,0,1,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588273478ull :
      sort_indices<2,3,1,0,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806508358ull :
      sort_indices<3,0,1,2,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807491398ull :
      sort_indices<3,0,2,1,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822236998ull :
      sort_indices<3,1,0,2,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824203078ull :
      sort_indices<3,1,2,0,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838948678ull :
      sort_indices<3,2,0,1,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839931718ull :
      sort_indices<3,2,1,0,5,7,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19093348ull :
      sort_indices<0,1,2,3,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20076388ull :
      sort_indices<0,1,3,2,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34821988ull :
      sort_indices<0,2,1,3,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36788068ull :
      sort_indices<0,2,3,1,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51533668ull :
      sort_indices<0,3,1,2,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52516708ull :
      sort_indices<0,3,2,1,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270751588ull :
      sort_indices<1,0,2,3,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271734628ull :
      sort_indices<1,0,3,2,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302208868ull :
      sort_indices<1,2,0,3,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305157988ull :
      sort_indices<1,2,3,0,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318920548ull :
      sort_indices<1,3,0,2,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320886628ull :
      sort_indices<1,3,2,0,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538138468ull :
      sort_indices<2,0,1,3,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540104548ull :
      sort_indices<2,0,3,1,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553867108ull :
      sort_indices<2,1,0,3,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556816228ull :
      sort_indices<2,1,3,0,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587290468ull :
      sort_indices<2,3,0,1,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588273508ull :
      sort_indices<2,3,1,0,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806508388ull :
      sort_indices<3,0,1,2,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807491428ull :
      sort_indices<3,0,2,1,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822237028ull :
      sort_indices<3,1,0,2,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824203108ull :
      sort_indices<3,1,2,0,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838948708ull :
      sort_indices<3,2,0,1,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839931748ull :
      sort_indices<3,2,1,0,5,7,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19096663ull :
      sort_indices<0,1,2,3,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20079703ull :
      sort_indices<0,1,3,2,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34825303ull :
      sort_indices<0,2,1,3,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36791383ull :
      sort_indices<0,2,3,1,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51536983ull :
      sort_indices<0,3,1,2,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52520023ull :
      sort_indices<0,3,2,1,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270754903ull :
      sort_indices<1,0,2,3,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271737943ull :
      sort_indices<1,0,3,2,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302212183ull :
      sort_indices<1,2,0,3,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305161303ull :
      sort_indices<1,2,3,0,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318923863ull :
      sort_indices<1,3,0,2,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320889943ull :
      sort_indices<1,3,2,0,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538141783ull :
      sort_indices<2,0,1,3,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540107863ull :
      sort_indices<2,0,3,1,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553870423ull :
      sort_indices<2,1,0,3,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556819543ull :
      sort_indices<2,1,3,0,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587293783ull :
      sort_indices<2,3,0,1,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588276823ull :
      sort_indices<2,3,1,0,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806511703ull :
      sort_indices<3,0,1,2,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807494743ull :
      sort_indices<3,0,2,1,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822240343ull :
      sort_indices<3,1,0,2,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824206423ull :
      sort_indices<3,1,2,0,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838952023ull :
      sort_indices<3,2,0,1,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839935063ull :
      sort_indices<3,2,1,0,6,4,5,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19096693ull :
      sort_indices<0,1,2,3,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20079733ull :
      sort_indices<0,1,3,2,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34825333ull :
      sort_indices<0,2,1,3,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36791413ull :
      sort_indices<0,2,3,1,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51537013ull :
      sort_indices<0,3,1,2,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52520053ull :
      sort_indices<0,3,2,1,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270754933ull :
      sort_indices<1,0,2,3,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271737973ull :
      sort_indices<1,0,3,2,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302212213ull :
      sort_indices<1,2,0,3,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305161333ull :
      sort_indices<1,2,3,0,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318923893ull :
      sort_indices<1,3,0,2,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320889973ull :
      sort_indices<1,3,2,0,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538141813ull :
      sort_indices<2,0,1,3,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540107893ull :
      sort_indices<2,0,3,1,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553870453ull :
      sort_indices<2,1,0,3,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556819573ull :
      sort_indices<2,1,3,0,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587293813ull :
      sort_indices<2,3,0,1,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588276853ull :
      sort_indices<2,3,1,0,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806511733ull :
      sort_indices<3,0,1,2,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807494773ull :
      sort_indices<3,0,2,1,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822240373ull :
      sort_indices<3,1,0,2,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824206453ull :
      sort_indices<3,1,2,0,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838952053ull :
      sort_indices<3,2,0,1,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839935093ull :
      sort_indices<3,2,1,0,6,4,7,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19096903ull :
      sort_indices<0,1,2,3,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20079943ull :
      sort_indices<0,1,3,2,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34825543ull :
      sort_indices<0,2,1,3,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36791623ull :
      sort_indices<0,2,3,1,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51537223ull :
      sort_indices<0,3,1,2,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52520263ull :
      sort_indices<0,3,2,1,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270755143ull :
      sort_indices<1,0,2,3,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271738183ull :
      sort_indices<1,0,3,2,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302212423ull :
      sort_indices<1,2,0,3,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305161543ull :
      sort_indices<1,2,3,0,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318924103ull :
      sort_indices<1,3,0,2,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320890183ull :
      sort_indices<1,3,2,0,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538142023ull :
      sort_indices<2,0,1,3,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540108103ull :
      sort_indices<2,0,3,1,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553870663ull :
      sort_indices<2,1,0,3,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556819783ull :
      sort_indices<2,1,3,0,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587294023ull :
      sort_indices<2,3,0,1,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588277063ull :
      sort_indices<2,3,1,0,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806511943ull :
      sort_indices<3,0,1,2,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807494983ull :
      sort_indices<3,0,2,1,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822240583ull :
      sort_indices<3,1,0,2,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824206663ull :
      sort_indices<3,1,2,0,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838952263ull :
      sort_indices<3,2,0,1,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839935303ull :
      sort_indices<3,2,1,0,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19096948ull :
      sort_indices<0,1,2,3,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20079988ull :
      sort_indices<0,1,3,2,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34825588ull :
      sort_indices<0,2,1,3,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36791668ull :
      sort_indices<0,2,3,1,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51537268ull :
      sort_indices<0,3,1,2,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52520308ull :
      sort_indices<0,3,2,1,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270755188ull :
      sort_indices<1,0,2,3,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271738228ull :
      sort_indices<1,0,3,2,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302212468ull :
      sort_indices<1,2,0,3,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305161588ull :
      sort_indices<1,2,3,0,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318924148ull :
      sort_indices<1,3,0,2,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320890228ull :
      sort_indices<1,3,2,0,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538142068ull :
      sort_indices<2,0,1,3,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540108148ull :
      sort_indices<2,0,3,1,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553870708ull :
      sort_indices<2,1,0,3,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556819828ull :
      sort_indices<2,1,3,0,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587294068ull :
      sort_indices<2,3,0,1,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588277108ull :
      sort_indices<2,3,1,0,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806511988ull :
      sort_indices<3,0,1,2,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807495028ull :
      sort_indices<3,0,2,1,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822240628ull :
      sort_indices<3,1,0,2,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824206708ull :
      sort_indices<3,1,2,0,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838952308ull :
      sort_indices<3,2,0,1,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839935348ull :
      sort_indices<3,2,1,0,6,5,7,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19097413ull :
      sort_indices<0,1,2,3,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20080453ull :
      sort_indices<0,1,3,2,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34826053ull :
      sort_indices<0,2,1,3,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36792133ull :
      sort_indices<0,2,3,1,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51537733ull :
      sort_indices<0,3,1,2,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52520773ull :
      sort_indices<0,3,2,1,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270755653ull :
      sort_indices<1,0,2,3,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271738693ull :
      sort_indices<1,0,3,2,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302212933ull :
      sort_indices<1,2,0,3,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305162053ull :
      sort_indices<1,2,3,0,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318924613ull :
      sort_indices<1,3,0,2,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320890693ull :
      sort_indices<1,3,2,0,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538142533ull :
      sort_indices<2,0,1,3,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540108613ull :
      sort_indices<2,0,3,1,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553871173ull :
      sort_indices<2,1,0,3,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556820293ull :
      sort_indices<2,1,3,0,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587294533ull :
      sort_indices<2,3,0,1,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588277573ull :
      sort_indices<2,3,1,0,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806512453ull :
      sort_indices<3,0,1,2,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807495493ull :
      sort_indices<3,0,2,1,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822241093ull :
      sort_indices<3,1,0,2,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824207173ull :
      sort_indices<3,1,2,0,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838952773ull :
      sort_indices<3,2,0,1,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839935813ull :
      sort_indices<3,2,1,0,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19097428ull :
      sort_indices<0,1,2,3,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20080468ull :
      sort_indices<0,1,3,2,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34826068ull :
      sort_indices<0,2,1,3,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36792148ull :
      sort_indices<0,2,3,1,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51537748ull :
      sort_indices<0,3,1,2,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52520788ull :
      sort_indices<0,3,2,1,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270755668ull :
      sort_indices<1,0,2,3,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271738708ull :
      sort_indices<1,0,3,2,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302212948ull :
      sort_indices<1,2,0,3,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305162068ull :
      sort_indices<1,2,3,0,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318924628ull :
      sort_indices<1,3,0,2,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320890708ull :
      sort_indices<1,3,2,0,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538142548ull :
      sort_indices<2,0,1,3,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540108628ull :
      sort_indices<2,0,3,1,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553871188ull :
      sort_indices<2,1,0,3,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556820308ull :
      sort_indices<2,1,3,0,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587294548ull :
      sort_indices<2,3,0,1,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588277588ull :
      sort_indices<2,3,1,0,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806512468ull :
      sort_indices<3,0,1,2,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807495508ull :
      sort_indices<3,0,2,1,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822241108ull :
      sort_indices<3,1,0,2,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824207188ull :
      sort_indices<3,1,2,0,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838952788ull :
      sort_indices<3,2,0,1,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839935828ull :
      sort_indices<3,2,1,0,6,7,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19100758ull :
      sort_indices<0,1,2,3,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20083798ull :
      sort_indices<0,1,3,2,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34829398ull :
      sort_indices<0,2,1,3,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36795478ull :
      sort_indices<0,2,3,1,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51541078ull :
      sort_indices<0,3,1,2,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52524118ull :
      sort_indices<0,3,2,1,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270758998ull :
      sort_indices<1,0,2,3,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271742038ull :
      sort_indices<1,0,3,2,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302216278ull :
      sort_indices<1,2,0,3,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305165398ull :
      sort_indices<1,2,3,0,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318927958ull :
      sort_indices<1,3,0,2,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320894038ull :
      sort_indices<1,3,2,0,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538145878ull :
      sort_indices<2,0,1,3,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540111958ull :
      sort_indices<2,0,3,1,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553874518ull :
      sort_indices<2,1,0,3,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556823638ull :
      sort_indices<2,1,3,0,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587297878ull :
      sort_indices<2,3,0,1,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588280918ull :
      sort_indices<2,3,1,0,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806515798ull :
      sort_indices<3,0,1,2,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807498838ull :
      sort_indices<3,0,2,1,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822244438ull :
      sort_indices<3,1,0,2,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824210518ull :
      sort_indices<3,1,2,0,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838956118ull :
      sort_indices<3,2,0,1,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839939158ull :
      sort_indices<3,2,1,0,7,4,5,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19100773ull :
      sort_indices<0,1,2,3,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20083813ull :
      sort_indices<0,1,3,2,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34829413ull :
      sort_indices<0,2,1,3,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36795493ull :
      sort_indices<0,2,3,1,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51541093ull :
      sort_indices<0,3,1,2,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52524133ull :
      sort_indices<0,3,2,1,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270759013ull :
      sort_indices<1,0,2,3,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271742053ull :
      sort_indices<1,0,3,2,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302216293ull :
      sort_indices<1,2,0,3,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305165413ull :
      sort_indices<1,2,3,0,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318927973ull :
      sort_indices<1,3,0,2,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320894053ull :
      sort_indices<1,3,2,0,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538145893ull :
      sort_indices<2,0,1,3,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540111973ull :
      sort_indices<2,0,3,1,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553874533ull :
      sort_indices<2,1,0,3,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556823653ull :
      sort_indices<2,1,3,0,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587297893ull :
      sort_indices<2,3,0,1,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588280933ull :
      sort_indices<2,3,1,0,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806515813ull :
      sort_indices<3,0,1,2,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807498853ull :
      sort_indices<3,0,2,1,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822244453ull :
      sort_indices<3,1,0,2,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824210533ull :
      sort_indices<3,1,2,0,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838956133ull :
      sort_indices<3,2,0,1,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839939173ull :
      sort_indices<3,2,1,0,7,4,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19100998ull :
      sort_indices<0,1,2,3,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20084038ull :
      sort_indices<0,1,3,2,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34829638ull :
      sort_indices<0,2,1,3,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36795718ull :
      sort_indices<0,2,3,1,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51541318ull :
      sort_indices<0,3,1,2,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52524358ull :
      sort_indices<0,3,2,1,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270759238ull :
      sort_indices<1,0,2,3,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271742278ull :
      sort_indices<1,0,3,2,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302216518ull :
      sort_indices<1,2,0,3,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305165638ull :
      sort_indices<1,2,3,0,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318928198ull :
      sort_indices<1,3,0,2,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320894278ull :
      sort_indices<1,3,2,0,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538146118ull :
      sort_indices<2,0,1,3,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540112198ull :
      sort_indices<2,0,3,1,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553874758ull :
      sort_indices<2,1,0,3,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556823878ull :
      sort_indices<2,1,3,0,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587298118ull :
      sort_indices<2,3,0,1,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588281158ull :
      sort_indices<2,3,1,0,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806516038ull :
      sort_indices<3,0,1,2,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807499078ull :
      sort_indices<3,0,2,1,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822244678ull :
      sort_indices<3,1,0,2,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824210758ull :
      sort_indices<3,1,2,0,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838956358ull :
      sort_indices<3,2,0,1,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839939398ull :
      sort_indices<3,2,1,0,7,5,4,6,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19101028ull :
      sort_indices<0,1,2,3,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20084068ull :
      sort_indices<0,1,3,2,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34829668ull :
      sort_indices<0,2,1,3,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36795748ull :
      sort_indices<0,2,3,1,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51541348ull :
      sort_indices<0,3,1,2,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52524388ull :
      sort_indices<0,3,2,1,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270759268ull :
      sort_indices<1,0,2,3,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271742308ull :
      sort_indices<1,0,3,2,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302216548ull :
      sort_indices<1,2,0,3,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305165668ull :
      sort_indices<1,2,3,0,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318928228ull :
      sort_indices<1,3,0,2,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320894308ull :
      sort_indices<1,3,2,0,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538146148ull :
      sort_indices<2,0,1,3,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540112228ull :
      sort_indices<2,0,3,1,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553874788ull :
      sort_indices<2,1,0,3,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556823908ull :
      sort_indices<2,1,3,0,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587298148ull :
      sort_indices<2,3,0,1,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588281188ull :
      sort_indices<2,3,1,0,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806516068ull :
      sort_indices<3,0,1,2,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807499108ull :
      sort_indices<3,0,2,1,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822244708ull :
      sort_indices<3,1,0,2,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824210788ull :
      sort_indices<3,1,2,0,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838956388ull :
      sort_indices<3,2,0,1,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839939428ull :
      sort_indices<3,2,1,0,7,5,6,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19101253ull :
      sort_indices<0,1,2,3,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20084293ull :
      sort_indices<0,1,3,2,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34829893ull :
      sort_indices<0,2,1,3,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36795973ull :
      sort_indices<0,2,3,1,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51541573ull :
      sort_indices<0,3,1,2,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52524613ull :
      sort_indices<0,3,2,1,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270759493ull :
      sort_indices<1,0,2,3,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271742533ull :
      sort_indices<1,0,3,2,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302216773ull :
      sort_indices<1,2,0,3,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305165893ull :
      sort_indices<1,2,3,0,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318928453ull :
      sort_indices<1,3,0,2,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320894533ull :
      sort_indices<1,3,2,0,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538146373ull :
      sort_indices<2,0,1,3,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540112453ull :
      sort_indices<2,0,3,1,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553875013ull :
      sort_indices<2,1,0,3,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556824133ull :
      sort_indices<2,1,3,0,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587298373ull :
      sort_indices<2,3,0,1,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588281413ull :
      sort_indices<2,3,1,0,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806516293ull :
      sort_indices<3,0,1,2,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807499333ull :
      sort_indices<3,0,2,1,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822244933ull :
      sort_indices<3,1,0,2,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824211013ull :
      sort_indices<3,1,2,0,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838956613ull :
      sort_indices<3,2,0,1,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839939653ull :
      sort_indices<3,2,1,0,7,6,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19101268ull :
      sort_indices<0,1,2,3,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20084308ull :
      sort_indices<0,1,3,2,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34829908ull :
      sort_indices<0,2,1,3,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36795988ull :
      sort_indices<0,2,3,1,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51541588ull :
      sort_indices<0,3,1,2,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52524628ull :
      sort_indices<0,3,2,1,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270759508ull :
      sort_indices<1,0,2,3,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271742548ull :
      sort_indices<1,0,3,2,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302216788ull :
      sort_indices<1,2,0,3,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305165908ull :
      sort_indices<1,2,3,0,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318928468ull :
      sort_indices<1,3,0,2,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320894548ull :
      sort_indices<1,3,2,0,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538146388ull :
      sort_indices<2,0,1,3,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540112468ull :
      sort_indices<2,0,3,1,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553875028ull :
      sort_indices<2,1,0,3,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556824148ull :
      sort_indices<2,1,3,0,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587298388ull :
      sort_indices<2,3,0,1,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588281428ull :
      sort_indices<2,3,1,0,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806516308ull :
      sort_indices<3,0,1,2,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807499348ull :
      sort_indices<3,0,2,1,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822244948ull :
      sort_indices<3,1,0,2,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824211028ull :
      sort_indices<3,1,2,0,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838956628ull :
      sort_indices<3,2,0,1,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839939668ull :
      sort_indices<3,2,1,0,7,6,5,4,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    }
  } else if (std::abs(a+1.0) < numerical_zero__ && std::abs(b) < numerical_zero__) {
    switch (tag) {
    case 19088743ull :
      sort_indices<0,1,2,3,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20071783ull :
      sort_indices<0,1,3,2,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34817383ull :
      sort_indices<0,2,1,3,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36783463ull :
      sort_indices<0,2,3,1,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51529063ull :
      sort_indices<0,3,1,2,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52512103ull :
      sort_indices<0,3,2,1,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270746983ull :
      sort_indices<1,0,2,3,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271730023ull :
      sort_indices<1,0,3,2,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302204263ull :
      sort_indices<1,2,0,3,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305153383ull :
      sort_indices<1,2,3,0,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318915943ull :
      sort_indices<1,3,0,2,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320882023ull :
      sort_indices<1,3,2,0,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538133863ull :
      sort_indices<2,0,1,3,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540099943ull :
      sort_indices<2,0,3,1,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553862503ull :
      sort_indices<2,1,0,3,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556811623ull :
      sort_indices<2,1,3,0,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587285863ull :
      sort_indices<2,3,0,1,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588268903ull :
      sort_indices<2,3,1,0,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806503783ull :
      sort_indices<3,0,1,2,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807486823ull :
      sort_indices<3,0,2,1,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822232423ull :
      sort_indices<3,1,0,2,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824198503ull :
      sort_indices<3,1,2,0,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838944103ull :
      sort_indices<3,2,0,1,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839927143ull :
      sort_indices<3,2,1,0,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19088758ull :
      sort_indices<0,1,2,3,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20071798ull :
      sort_indices<0,1,3,2,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34817398ull :
      sort_indices<0,2,1,3,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36783478ull :
      sort_indices<0,2,3,1,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51529078ull :
      sort_indices<0,3,1,2,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52512118ull :
      sort_indices<0,3,2,1,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270746998ull :
      sort_indices<1,0,2,3,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271730038ull :
      sort_indices<1,0,3,2,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302204278ull :
      sort_indices<1,2,0,3,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305153398ull :
      sort_indices<1,2,3,0,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318915958ull :
      sort_indices<1,3,0,2,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320882038ull :
      sort_indices<1,3,2,0,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538133878ull :
      sort_indices<2,0,1,3,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540099958ull :
      sort_indices<2,0,3,1,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553862518ull :
      sort_indices<2,1,0,3,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556811638ull :
      sort_indices<2,1,3,0,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587285878ull :
      sort_indices<2,3,0,1,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588268918ull :
      sort_indices<2,3,1,0,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806503798ull :
      sort_indices<3,0,1,2,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807486838ull :
      sort_indices<3,0,2,1,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822232438ull :
      sort_indices<3,1,0,2,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824198518ull :
      sort_indices<3,1,2,0,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838944118ull :
      sort_indices<3,2,0,1,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839927158ull :
      sort_indices<3,2,1,0,4,5,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19088983ull :
      sort_indices<0,1,2,3,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20072023ull :
      sort_indices<0,1,3,2,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34817623ull :
      sort_indices<0,2,1,3,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36783703ull :
      sort_indices<0,2,3,1,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51529303ull :
      sort_indices<0,3,1,2,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52512343ull :
      sort_indices<0,3,2,1,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270747223ull :
      sort_indices<1,0,2,3,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271730263ull :
      sort_indices<1,0,3,2,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302204503ull :
      sort_indices<1,2,0,3,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305153623ull :
      sort_indices<1,2,3,0,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318916183ull :
      sort_indices<1,3,0,2,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320882263ull :
      sort_indices<1,3,2,0,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538134103ull :
      sort_indices<2,0,1,3,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540100183ull :
      sort_indices<2,0,3,1,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553862743ull :
      sort_indices<2,1,0,3,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556811863ull :
      sort_indices<2,1,3,0,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587286103ull :
      sort_indices<2,3,0,1,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588269143ull :
      sort_indices<2,3,1,0,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806504023ull :
      sort_indices<3,0,1,2,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807487063ull :
      sort_indices<3,0,2,1,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822232663ull :
      sort_indices<3,1,0,2,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824198743ull :
      sort_indices<3,1,2,0,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838944343ull :
      sort_indices<3,2,0,1,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839927383ull :
      sort_indices<3,2,1,0,4,6,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19089013ull :
      sort_indices<0,1,2,3,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20072053ull :
      sort_indices<0,1,3,2,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34817653ull :
      sort_indices<0,2,1,3,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36783733ull :
      sort_indices<0,2,3,1,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51529333ull :
      sort_indices<0,3,1,2,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52512373ull :
      sort_indices<0,3,2,1,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270747253ull :
      sort_indices<1,0,2,3,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271730293ull :
      sort_indices<1,0,3,2,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302204533ull :
      sort_indices<1,2,0,3,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305153653ull :
      sort_indices<1,2,3,0,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318916213ull :
      sort_indices<1,3,0,2,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320882293ull :
      sort_indices<1,3,2,0,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538134133ull :
      sort_indices<2,0,1,3,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540100213ull :
      sort_indices<2,0,3,1,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553862773ull :
      sort_indices<2,1,0,3,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556811893ull :
      sort_indices<2,1,3,0,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587286133ull :
      sort_indices<2,3,0,1,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588269173ull :
      sort_indices<2,3,1,0,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806504053ull :
      sort_indices<3,0,1,2,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807487093ull :
      sort_indices<3,0,2,1,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822232693ull :
      sort_indices<3,1,0,2,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824198773ull :
      sort_indices<3,1,2,0,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838944373ull :
      sort_indices<3,2,0,1,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839927413ull :
      sort_indices<3,2,1,0,4,6,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19089238ull :
      sort_indices<0,1,2,3,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20072278ull :
      sort_indices<0,1,3,2,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34817878ull :
      sort_indices<0,2,1,3,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36783958ull :
      sort_indices<0,2,3,1,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51529558ull :
      sort_indices<0,3,1,2,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52512598ull :
      sort_indices<0,3,2,1,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270747478ull :
      sort_indices<1,0,2,3,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271730518ull :
      sort_indices<1,0,3,2,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302204758ull :
      sort_indices<1,2,0,3,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305153878ull :
      sort_indices<1,2,3,0,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318916438ull :
      sort_indices<1,3,0,2,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320882518ull :
      sort_indices<1,3,2,0,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538134358ull :
      sort_indices<2,0,1,3,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540100438ull :
      sort_indices<2,0,3,1,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553862998ull :
      sort_indices<2,1,0,3,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556812118ull :
      sort_indices<2,1,3,0,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587286358ull :
      sort_indices<2,3,0,1,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588269398ull :
      sort_indices<2,3,1,0,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806504278ull :
      sort_indices<3,0,1,2,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807487318ull :
      sort_indices<3,0,2,1,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822232918ull :
      sort_indices<3,1,0,2,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824198998ull :
      sort_indices<3,1,2,0,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838944598ull :
      sort_indices<3,2,0,1,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839927638ull :
      sort_indices<3,2,1,0,4,7,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19089253ull :
      sort_indices<0,1,2,3,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20072293ull :
      sort_indices<0,1,3,2,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34817893ull :
      sort_indices<0,2,1,3,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36783973ull :
      sort_indices<0,2,3,1,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51529573ull :
      sort_indices<0,3,1,2,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52512613ull :
      sort_indices<0,3,2,1,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270747493ull :
      sort_indices<1,0,2,3,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271730533ull :
      sort_indices<1,0,3,2,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302204773ull :
      sort_indices<1,2,0,3,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305153893ull :
      sort_indices<1,2,3,0,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318916453ull :
      sort_indices<1,3,0,2,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320882533ull :
      sort_indices<1,3,2,0,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538134373ull :
      sort_indices<2,0,1,3,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540100453ull :
      sort_indices<2,0,3,1,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553863013ull :
      sort_indices<2,1,0,3,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556812133ull :
      sort_indices<2,1,3,0,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587286373ull :
      sort_indices<2,3,0,1,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588269413ull :
      sort_indices<2,3,1,0,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806504293ull :
      sort_indices<3,0,1,2,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807487333ull :
      sort_indices<3,0,2,1,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822232933ull :
      sort_indices<3,1,0,2,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824199013ull :
      sort_indices<3,1,2,0,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838944613ull :
      sort_indices<3,2,0,1,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839927653ull :
      sort_indices<3,2,1,0,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19092583ull :
      sort_indices<0,1,2,3,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20075623ull :
      sort_indices<0,1,3,2,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34821223ull :
      sort_indices<0,2,1,3,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36787303ull :
      sort_indices<0,2,3,1,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51532903ull :
      sort_indices<0,3,1,2,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52515943ull :
      sort_indices<0,3,2,1,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270750823ull :
      sort_indices<1,0,2,3,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271733863ull :
      sort_indices<1,0,3,2,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302208103ull :
      sort_indices<1,2,0,3,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305157223ull :
      sort_indices<1,2,3,0,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318919783ull :
      sort_indices<1,3,0,2,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320885863ull :
      sort_indices<1,3,2,0,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538137703ull :
      sort_indices<2,0,1,3,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540103783ull :
      sort_indices<2,0,3,1,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553866343ull :
      sort_indices<2,1,0,3,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556815463ull :
      sort_indices<2,1,3,0,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587289703ull :
      sort_indices<2,3,0,1,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588272743ull :
      sort_indices<2,3,1,0,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806507623ull :
      sort_indices<3,0,1,2,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807490663ull :
      sort_indices<3,0,2,1,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822236263ull :
      sort_indices<3,1,0,2,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824202343ull :
      sort_indices<3,1,2,0,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838947943ull :
      sort_indices<3,2,0,1,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839930983ull :
      sort_indices<3,2,1,0,5,4,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19092598ull :
      sort_indices<0,1,2,3,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20075638ull :
      sort_indices<0,1,3,2,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34821238ull :
      sort_indices<0,2,1,3,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36787318ull :
      sort_indices<0,2,3,1,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51532918ull :
      sort_indices<0,3,1,2,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52515958ull :
      sort_indices<0,3,2,1,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270750838ull :
      sort_indices<1,0,2,3,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271733878ull :
      sort_indices<1,0,3,2,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302208118ull :
      sort_indices<1,2,0,3,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305157238ull :
      sort_indices<1,2,3,0,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318919798ull :
      sort_indices<1,3,0,2,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320885878ull :
      sort_indices<1,3,2,0,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538137718ull :
      sort_indices<2,0,1,3,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540103798ull :
      sort_indices<2,0,3,1,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553866358ull :
      sort_indices<2,1,0,3,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556815478ull :
      sort_indices<2,1,3,0,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587289718ull :
      sort_indices<2,3,0,1,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588272758ull :
      sort_indices<2,3,1,0,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806507638ull :
      sort_indices<3,0,1,2,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807490678ull :
      sort_indices<3,0,2,1,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822236278ull :
      sort_indices<3,1,0,2,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824202358ull :
      sort_indices<3,1,2,0,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838947958ull :
      sort_indices<3,2,0,1,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839930998ull :
      sort_indices<3,2,1,0,5,4,7,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19093063ull :
      sort_indices<0,1,2,3,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20076103ull :
      sort_indices<0,1,3,2,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34821703ull :
      sort_indices<0,2,1,3,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36787783ull :
      sort_indices<0,2,3,1,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51533383ull :
      sort_indices<0,3,1,2,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52516423ull :
      sort_indices<0,3,2,1,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270751303ull :
      sort_indices<1,0,2,3,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271734343ull :
      sort_indices<1,0,3,2,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302208583ull :
      sort_indices<1,2,0,3,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305157703ull :
      sort_indices<1,2,3,0,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318920263ull :
      sort_indices<1,3,0,2,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320886343ull :
      sort_indices<1,3,2,0,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538138183ull :
      sort_indices<2,0,1,3,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540104263ull :
      sort_indices<2,0,3,1,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553866823ull :
      sort_indices<2,1,0,3,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556815943ull :
      sort_indices<2,1,3,0,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587290183ull :
      sort_indices<2,3,0,1,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588273223ull :
      sort_indices<2,3,1,0,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806508103ull :
      sort_indices<3,0,1,2,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807491143ull :
      sort_indices<3,0,2,1,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822236743ull :
      sort_indices<3,1,0,2,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824202823ull :
      sort_indices<3,1,2,0,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838948423ull :
      sort_indices<3,2,0,1,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839931463ull :
      sort_indices<3,2,1,0,5,6,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19093108ull :
      sort_indices<0,1,2,3,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20076148ull :
      sort_indices<0,1,3,2,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34821748ull :
      sort_indices<0,2,1,3,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36787828ull :
      sort_indices<0,2,3,1,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51533428ull :
      sort_indices<0,3,1,2,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52516468ull :
      sort_indices<0,3,2,1,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270751348ull :
      sort_indices<1,0,2,3,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271734388ull :
      sort_indices<1,0,3,2,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302208628ull :
      sort_indices<1,2,0,3,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305157748ull :
      sort_indices<1,2,3,0,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318920308ull :
      sort_indices<1,3,0,2,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320886388ull :
      sort_indices<1,3,2,0,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538138228ull :
      sort_indices<2,0,1,3,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540104308ull :
      sort_indices<2,0,3,1,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553866868ull :
      sort_indices<2,1,0,3,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556815988ull :
      sort_indices<2,1,3,0,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587290228ull :
      sort_indices<2,3,0,1,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588273268ull :
      sort_indices<2,3,1,0,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806508148ull :
      sort_indices<3,0,1,2,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807491188ull :
      sort_indices<3,0,2,1,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822236788ull :
      sort_indices<3,1,0,2,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824202868ull :
      sort_indices<3,1,2,0,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838948468ull :
      sort_indices<3,2,0,1,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839931508ull :
      sort_indices<3,2,1,0,5,6,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19093318ull :
      sort_indices<0,1,2,3,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20076358ull :
      sort_indices<0,1,3,2,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34821958ull :
      sort_indices<0,2,1,3,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36788038ull :
      sort_indices<0,2,3,1,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51533638ull :
      sort_indices<0,3,1,2,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52516678ull :
      sort_indices<0,3,2,1,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270751558ull :
      sort_indices<1,0,2,3,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271734598ull :
      sort_indices<1,0,3,2,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302208838ull :
      sort_indices<1,2,0,3,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305157958ull :
      sort_indices<1,2,3,0,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318920518ull :
      sort_indices<1,3,0,2,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320886598ull :
      sort_indices<1,3,2,0,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538138438ull :
      sort_indices<2,0,1,3,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540104518ull :
      sort_indices<2,0,3,1,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553867078ull :
      sort_indices<2,1,0,3,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556816198ull :
      sort_indices<2,1,3,0,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587290438ull :
      sort_indices<2,3,0,1,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588273478ull :
      sort_indices<2,3,1,0,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806508358ull :
      sort_indices<3,0,1,2,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807491398ull :
      sort_indices<3,0,2,1,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822236998ull :
      sort_indices<3,1,0,2,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824203078ull :
      sort_indices<3,1,2,0,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838948678ull :
      sort_indices<3,2,0,1,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839931718ull :
      sort_indices<3,2,1,0,5,7,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19093348ull :
      sort_indices<0,1,2,3,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20076388ull :
      sort_indices<0,1,3,2,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34821988ull :
      sort_indices<0,2,1,3,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36788068ull :
      sort_indices<0,2,3,1,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51533668ull :
      sort_indices<0,3,1,2,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52516708ull :
      sort_indices<0,3,2,1,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270751588ull :
      sort_indices<1,0,2,3,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271734628ull :
      sort_indices<1,0,3,2,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302208868ull :
      sort_indices<1,2,0,3,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305157988ull :
      sort_indices<1,2,3,0,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318920548ull :
      sort_indices<1,3,0,2,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320886628ull :
      sort_indices<1,3,2,0,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538138468ull :
      sort_indices<2,0,1,3,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540104548ull :
      sort_indices<2,0,3,1,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553867108ull :
      sort_indices<2,1,0,3,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556816228ull :
      sort_indices<2,1,3,0,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587290468ull :
      sort_indices<2,3,0,1,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588273508ull :
      sort_indices<2,3,1,0,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806508388ull :
      sort_indices<3,0,1,2,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807491428ull :
      sort_indices<3,0,2,1,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822237028ull :
      sort_indices<3,1,0,2,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824203108ull :
      sort_indices<3,1,2,0,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838948708ull :
      sort_indices<3,2,0,1,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839931748ull :
      sort_indices<3,2,1,0,5,7,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19096663ull :
      sort_indices<0,1,2,3,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20079703ull :
      sort_indices<0,1,3,2,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34825303ull :
      sort_indices<0,2,1,3,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36791383ull :
      sort_indices<0,2,3,1,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51536983ull :
      sort_indices<0,3,1,2,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52520023ull :
      sort_indices<0,3,2,1,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270754903ull :
      sort_indices<1,0,2,3,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271737943ull :
      sort_indices<1,0,3,2,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302212183ull :
      sort_indices<1,2,0,3,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305161303ull :
      sort_indices<1,2,3,0,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318923863ull :
      sort_indices<1,3,0,2,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320889943ull :
      sort_indices<1,3,2,0,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538141783ull :
      sort_indices<2,0,1,3,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540107863ull :
      sort_indices<2,0,3,1,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553870423ull :
      sort_indices<2,1,0,3,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556819543ull :
      sort_indices<2,1,3,0,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587293783ull :
      sort_indices<2,3,0,1,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588276823ull :
      sort_indices<2,3,1,0,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806511703ull :
      sort_indices<3,0,1,2,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807494743ull :
      sort_indices<3,0,2,1,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822240343ull :
      sort_indices<3,1,0,2,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824206423ull :
      sort_indices<3,1,2,0,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838952023ull :
      sort_indices<3,2,0,1,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839935063ull :
      sort_indices<3,2,1,0,6,4,5,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19096693ull :
      sort_indices<0,1,2,3,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20079733ull :
      sort_indices<0,1,3,2,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34825333ull :
      sort_indices<0,2,1,3,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36791413ull :
      sort_indices<0,2,3,1,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51537013ull :
      sort_indices<0,3,1,2,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52520053ull :
      sort_indices<0,3,2,1,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270754933ull :
      sort_indices<1,0,2,3,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271737973ull :
      sort_indices<1,0,3,2,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302212213ull :
      sort_indices<1,2,0,3,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305161333ull :
      sort_indices<1,2,3,0,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318923893ull :
      sort_indices<1,3,0,2,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320889973ull :
      sort_indices<1,3,2,0,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538141813ull :
      sort_indices<2,0,1,3,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540107893ull :
      sort_indices<2,0,3,1,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553870453ull :
      sort_indices<2,1,0,3,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556819573ull :
      sort_indices<2,1,3,0,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587293813ull :
      sort_indices<2,3,0,1,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588276853ull :
      sort_indices<2,3,1,0,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806511733ull :
      sort_indices<3,0,1,2,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807494773ull :
      sort_indices<3,0,2,1,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822240373ull :
      sort_indices<3,1,0,2,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824206453ull :
      sort_indices<3,1,2,0,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838952053ull :
      sort_indices<3,2,0,1,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839935093ull :
      sort_indices<3,2,1,0,6,4,7,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19096903ull :
      sort_indices<0,1,2,3,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20079943ull :
      sort_indices<0,1,3,2,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34825543ull :
      sort_indices<0,2,1,3,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36791623ull :
      sort_indices<0,2,3,1,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51537223ull :
      sort_indices<0,3,1,2,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52520263ull :
      sort_indices<0,3,2,1,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270755143ull :
      sort_indices<1,0,2,3,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271738183ull :
      sort_indices<1,0,3,2,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302212423ull :
      sort_indices<1,2,0,3,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305161543ull :
      sort_indices<1,2,3,0,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318924103ull :
      sort_indices<1,3,0,2,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320890183ull :
      sort_indices<1,3,2,0,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538142023ull :
      sort_indices<2,0,1,3,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540108103ull :
      sort_indices<2,0,3,1,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553870663ull :
      sort_indices<2,1,0,3,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556819783ull :
      sort_indices<2,1,3,0,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587294023ull :
      sort_indices<2,3,0,1,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588277063ull :
      sort_indices<2,3,1,0,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806511943ull :
      sort_indices<3,0,1,2,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807494983ull :
      sort_indices<3,0,2,1,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822240583ull :
      sort_indices<3,1,0,2,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824206663ull :
      sort_indices<3,1,2,0,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838952263ull :
      sort_indices<3,2,0,1,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839935303ull :
      sort_indices<3,2,1,0,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19096948ull :
      sort_indices<0,1,2,3,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20079988ull :
      sort_indices<0,1,3,2,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34825588ull :
      sort_indices<0,2,1,3,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36791668ull :
      sort_indices<0,2,3,1,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51537268ull :
      sort_indices<0,3,1,2,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52520308ull :
      sort_indices<0,3,2,1,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270755188ull :
      sort_indices<1,0,2,3,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271738228ull :
      sort_indices<1,0,3,2,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302212468ull :
      sort_indices<1,2,0,3,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305161588ull :
      sort_indices<1,2,3,0,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318924148ull :
      sort_indices<1,3,0,2,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320890228ull :
      sort_indices<1,3,2,0,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538142068ull :
      sort_indices<2,0,1,3,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540108148ull :
      sort_indices<2,0,3,1,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553870708ull :
      sort_indices<2,1,0,3,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556819828ull :
      sort_indices<2,1,3,0,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587294068ull :
      sort_indices<2,3,0,1,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588277108ull :
      sort_indices<2,3,1,0,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806511988ull :
      sort_indices<3,0,1,2,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807495028ull :
      sort_indices<3,0,2,1,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822240628ull :
      sort_indices<3,1,0,2,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824206708ull :
      sort_indices<3,1,2,0,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838952308ull :
      sort_indices<3,2,0,1,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839935348ull :
      sort_indices<3,2,1,0,6,5,7,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19097413ull :
      sort_indices<0,1,2,3,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20080453ull :
      sort_indices<0,1,3,2,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34826053ull :
      sort_indices<0,2,1,3,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36792133ull :
      sort_indices<0,2,3,1,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51537733ull :
      sort_indices<0,3,1,2,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52520773ull :
      sort_indices<0,3,2,1,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270755653ull :
      sort_indices<1,0,2,3,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271738693ull :
      sort_indices<1,0,3,2,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302212933ull :
      sort_indices<1,2,0,3,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305162053ull :
      sort_indices<1,2,3,0,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318924613ull :
      sort_indices<1,3,0,2,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320890693ull :
      sort_indices<1,3,2,0,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538142533ull :
      sort_indices<2,0,1,3,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540108613ull :
      sort_indices<2,0,3,1,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553871173ull :
      sort_indices<2,1,0,3,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556820293ull :
      sort_indices<2,1,3,0,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587294533ull :
      sort_indices<2,3,0,1,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588277573ull :
      sort_indices<2,3,1,0,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806512453ull :
      sort_indices<3,0,1,2,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807495493ull :
      sort_indices<3,0,2,1,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822241093ull :
      sort_indices<3,1,0,2,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824207173ull :
      sort_indices<3,1,2,0,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838952773ull :
      sort_indices<3,2,0,1,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839935813ull :
      sort_indices<3,2,1,0,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19097428ull :
      sort_indices<0,1,2,3,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20080468ull :
      sort_indices<0,1,3,2,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34826068ull :
      sort_indices<0,2,1,3,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36792148ull :
      sort_indices<0,2,3,1,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51537748ull :
      sort_indices<0,3,1,2,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52520788ull :
      sort_indices<0,3,2,1,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270755668ull :
      sort_indices<1,0,2,3,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271738708ull :
      sort_indices<1,0,3,2,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302212948ull :
      sort_indices<1,2,0,3,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305162068ull :
      sort_indices<1,2,3,0,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318924628ull :
      sort_indices<1,3,0,2,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320890708ull :
      sort_indices<1,3,2,0,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538142548ull :
      sort_indices<2,0,1,3,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540108628ull :
      sort_indices<2,0,3,1,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553871188ull :
      sort_indices<2,1,0,3,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556820308ull :
      sort_indices<2,1,3,0,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587294548ull :
      sort_indices<2,3,0,1,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588277588ull :
      sort_indices<2,3,1,0,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806512468ull :
      sort_indices<3,0,1,2,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807495508ull :
      sort_indices<3,0,2,1,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822241108ull :
      sort_indices<3,1,0,2,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824207188ull :
      sort_indices<3,1,2,0,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838952788ull :
      sort_indices<3,2,0,1,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839935828ull :
      sort_indices<3,2,1,0,6,7,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19100758ull :
      sort_indices<0,1,2,3,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20083798ull :
      sort_indices<0,1,3,2,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34829398ull :
      sort_indices<0,2,1,3,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36795478ull :
      sort_indices<0,2,3,1,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51541078ull :
      sort_indices<0,3,1,2,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52524118ull :
      sort_indices<0,3,2,1,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270758998ull :
      sort_indices<1,0,2,3,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271742038ull :
      sort_indices<1,0,3,2,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302216278ull :
      sort_indices<1,2,0,3,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305165398ull :
      sort_indices<1,2,3,0,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318927958ull :
      sort_indices<1,3,0,2,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320894038ull :
      sort_indices<1,3,2,0,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538145878ull :
      sort_indices<2,0,1,3,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540111958ull :
      sort_indices<2,0,3,1,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553874518ull :
      sort_indices<2,1,0,3,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556823638ull :
      sort_indices<2,1,3,0,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587297878ull :
      sort_indices<2,3,0,1,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588280918ull :
      sort_indices<2,3,1,0,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806515798ull :
      sort_indices<3,0,1,2,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807498838ull :
      sort_indices<3,0,2,1,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822244438ull :
      sort_indices<3,1,0,2,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824210518ull :
      sort_indices<3,1,2,0,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838956118ull :
      sort_indices<3,2,0,1,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839939158ull :
      sort_indices<3,2,1,0,7,4,5,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19100773ull :
      sort_indices<0,1,2,3,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20083813ull :
      sort_indices<0,1,3,2,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34829413ull :
      sort_indices<0,2,1,3,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36795493ull :
      sort_indices<0,2,3,1,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51541093ull :
      sort_indices<0,3,1,2,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52524133ull :
      sort_indices<0,3,2,1,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270759013ull :
      sort_indices<1,0,2,3,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271742053ull :
      sort_indices<1,0,3,2,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302216293ull :
      sort_indices<1,2,0,3,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305165413ull :
      sort_indices<1,2,3,0,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318927973ull :
      sort_indices<1,3,0,2,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320894053ull :
      sort_indices<1,3,2,0,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538145893ull :
      sort_indices<2,0,1,3,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540111973ull :
      sort_indices<2,0,3,1,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553874533ull :
      sort_indices<2,1,0,3,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556823653ull :
      sort_indices<2,1,3,0,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587297893ull :
      sort_indices<2,3,0,1,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588280933ull :
      sort_indices<2,3,1,0,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806515813ull :
      sort_indices<3,0,1,2,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807498853ull :
      sort_indices<3,0,2,1,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822244453ull :
      sort_indices<3,1,0,2,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824210533ull :
      sort_indices<3,1,2,0,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838956133ull :
      sort_indices<3,2,0,1,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839939173ull :
      sort_indices<3,2,1,0,7,4,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19100998ull :
      sort_indices<0,1,2,3,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20084038ull :
      sort_indices<0,1,3,2,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34829638ull :
      sort_indices<0,2,1,3,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36795718ull :
      sort_indices<0,2,3,1,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51541318ull :
      sort_indices<0,3,1,2,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52524358ull :
      sort_indices<0,3,2,1,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270759238ull :
      sort_indices<1,0,2,3,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271742278ull :
      sort_indices<1,0,3,2,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302216518ull :
      sort_indices<1,2,0,3,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305165638ull :
      sort_indices<1,2,3,0,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318928198ull :
      sort_indices<1,3,0,2,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320894278ull :
      sort_indices<1,3,2,0,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538146118ull :
      sort_indices<2,0,1,3,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540112198ull :
      sort_indices<2,0,3,1,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553874758ull :
      sort_indices<2,1,0,3,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556823878ull :
      sort_indices<2,1,3,0,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587298118ull :
      sort_indices<2,3,0,1,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588281158ull :
      sort_indices<2,3,1,0,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806516038ull :
      sort_indices<3,0,1,2,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807499078ull :
      sort_indices<3,0,2,1,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822244678ull :
      sort_indices<3,1,0,2,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824210758ull :
      sort_indices<3,1,2,0,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838956358ull :
      sort_indices<3,2,0,1,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839939398ull :
      sort_indices<3,2,1,0,7,5,4,6,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19101028ull :
      sort_indices<0,1,2,3,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20084068ull :
      sort_indices<0,1,3,2,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34829668ull :
      sort_indices<0,2,1,3,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36795748ull :
      sort_indices<0,2,3,1,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51541348ull :
      sort_indices<0,3,1,2,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52524388ull :
      sort_indices<0,3,2,1,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270759268ull :
      sort_indices<1,0,2,3,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271742308ull :
      sort_indices<1,0,3,2,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302216548ull :
      sort_indices<1,2,0,3,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305165668ull :
      sort_indices<1,2,3,0,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318928228ull :
      sort_indices<1,3,0,2,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320894308ull :
      sort_indices<1,3,2,0,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538146148ull :
      sort_indices<2,0,1,3,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540112228ull :
      sort_indices<2,0,3,1,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553874788ull :
      sort_indices<2,1,0,3,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556823908ull :
      sort_indices<2,1,3,0,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587298148ull :
      sort_indices<2,3,0,1,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588281188ull :
      sort_indices<2,3,1,0,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806516068ull :
      sort_indices<3,0,1,2,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807499108ull :
      sort_indices<3,0,2,1,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822244708ull :
      sort_indices<3,1,0,2,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824210788ull :
      sort_indices<3,1,2,0,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838956388ull :
      sort_indices<3,2,0,1,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839939428ull :
      sort_indices<3,2,1,0,7,5,6,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19101253ull :
      sort_indices<0,1,2,3,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20084293ull :
      sort_indices<0,1,3,2,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34829893ull :
      sort_indices<0,2,1,3,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36795973ull :
      sort_indices<0,2,3,1,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51541573ull :
      sort_indices<0,3,1,2,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52524613ull :
      sort_indices<0,3,2,1,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270759493ull :
      sort_indices<1,0,2,3,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271742533ull :
      sort_indices<1,0,3,2,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302216773ull :
      sort_indices<1,2,0,3,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305165893ull :
      sort_indices<1,2,3,0,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318928453ull :
      sort_indices<1,3,0,2,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320894533ull :
      sort_indices<1,3,2,0,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538146373ull :
      sort_indices<2,0,1,3,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540112453ull :
      sort_indices<2,0,3,1,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553875013ull :
      sort_indices<2,1,0,3,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556824133ull :
      sort_indices<2,1,3,0,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587298373ull :
      sort_indices<2,3,0,1,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588281413ull :
      sort_indices<2,3,1,0,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806516293ull :
      sort_indices<3,0,1,2,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807499333ull :
      sort_indices<3,0,2,1,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822244933ull :
      sort_indices<3,1,0,2,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824211013ull :
      sort_indices<3,1,2,0,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838956613ull :
      sort_indices<3,2,0,1,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839939653ull :
      sort_indices<3,2,1,0,7,6,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 19101268ull :
      sort_indices<0,1,2,3,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 20084308ull :
      sort_indices<0,1,3,2,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 34829908ull :
      sort_indices<0,2,1,3,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 36795988ull :
      sort_indices<0,2,3,1,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 51541588ull :
      sort_indices<0,3,1,2,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 52524628ull :
      sort_indices<0,3,2,1,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 270759508ull :
      sort_indices<1,0,2,3,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 271742548ull :
      sort_indices<1,0,3,2,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 302216788ull :
      sort_indices<1,2,0,3,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 305165908ull :
      sort_indices<1,2,3,0,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 318928468ull :
      sort_indices<1,3,0,2,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 320894548ull :
      sort_indices<1,3,2,0,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 538146388ull :
      sort_indices<2,0,1,3,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 540112468ull :
      sort_indices<2,0,3,1,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 553875028ull :
      sort_indices<2,1,0,3,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 556824148ull :
      sort_indices<2,1,3,0,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 587298388ull :
      sort_indices<2,3,0,1,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 588281428ull :
      sort_indices<2,3,1,0,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 806516308ull :
      sort_indices<3,0,1,2,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 807499348ull :
      sort_indices<3,0,2,1,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 822244948ull :
      sort_indices<3,1,0,2,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 824211028ull :
      sort_indices<3,1,2,0,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 838956628ull :
      sort_indices<3,2,0,1,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    case 839939668ull :
      sort_indices<3,2,1,0,7,6,5,4,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
    }
  } else {
    throw std::logic_error("This case has not been implemented in prim_op_var.h");
  }
}

}

#endif
