//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: prim_op_var.h
// Copyright (C) 2015 Toru Shiozaki
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

#ifndef __SRC_UTIL_PRIM_OP_VAR_H
#define __SRC_UTIL_PRIM_OP_VAR_H

#include <src/util/prim_op.h>
#include <src/util/constants.h>

namespace bagel {

template<typename DataType>
void sort_indices(const std::array<int,1>& o, const double a, const double b, const DataType* in, DataType* out, const std::array<int,1>& d) {
  throw std::logic_error("This case has not been implemented in prim_op_var.h");
}

template<typename DataType>
void sort_indices(const std::array<int,2>& o, const double a, const double b, const DataType* in, DataType* out, const std::array<int,2>& d) {
  throw std::logic_error("This case has not been implemented in prim_op_var.h");
}

template<typename DataType>
void sort_indices(const std::array<int,3>& o, const double a, const double b, const DataType* in, DataType* out, const std::array<int,3>& d) {
  throw std::logic_error("This case has not been implemented in prim_op_var.h");
}

template<typename DataType>
void sort_indices(const std::array<int,5>& o, const double a, const double b, const DataType* in, DataType* out, const std::array<int,5>& d) {
  throw std::logic_error("This case has not been implemented in prim_op_var.h");
}

template<typename DataType>
void sort_indices(const std::array<int,7>& o, const double a, const double b, const DataType* in, DataType* out, const std::array<int,7>& d) {
  throw std::logic_error("This case has not been implemented in prim_op_var.h");
}
template<typename DataType>
void sort_indices(const std::array<int,4>& o, const double a, const double b, const DataType* in, DataType* out, const std::array<int,4>& d) {
  size_t tag = 0ull;
  for (int i = 0; i != 4; ++i)
    tag = (tag << 4) + o[i];

  if (std::abs(1.0-a) < numerical_zero__ && std::abs(b) < numerical_zero__) {
    switch (tag) {
        case 291ull :
          sort_indices<0,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 801ull :
          sort_indices<0,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 4146ull :
          sort_indices<1,0,3,2,0,1,1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 4656ull :
          sort_indices<1,2,3,0,0,1,1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 8451ull :
          sort_indices<2,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 8961ull :
          sort_indices<2,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 12306ull :
          sort_indices<3,0,1,2,0,1,1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 12816ull :
          sort_indices<3,2,1,0,0,1,1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        default:
          throw std::logic_error("This case has not been implemented in prim_op_var.h");
    }
  }
  else if (std::abs(-1.0-a) < numerical_zero__ && std::abs(b) < numerical_zero__) {
    switch (tag) {
        case 291ull :
          sort_indices<0,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 801ull :
          sort_indices<0,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 4146ull :
          sort_indices<1,0,3,2,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 4656ull :
          sort_indices<1,2,3,0,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 8451ull :
          sort_indices<2,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 8961ull :
          sort_indices<2,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 12306ull :
          sort_indices<3,0,1,2,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        case 12816ull :
          sort_indices<3,2,1,0,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3]); break;
        default:
          throw std::logic_error("This case has not been implemented in prim_op_var.h");
    }
  }
  else {
    throw std::logic_error("This case has not been implemented in prim_op_var.h");
  }
}

template<typename DataType>
void sort_indices(const std::array<int,6>& o, const double a, const double b, const DataType* in, DataType* out, const std::array<int,6>& d) {
  size_t tag = 0ull;
  for (int i = 0; i != 6; ++i)
    tag = (tag << 4) + o[i];

  if (std::abs(1.0-a) < numerical_zero__ && std::abs(b) < numerical_zero__) {
    switch (tag) {
        case 74565ull :
          sort_indices<0,1,2,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 75075ull :
          sort_indices<0,1,2,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 205125ull :
          sort_indices<0,3,2,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 206145ull :
          sort_indices<0,3,2,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 336195ull :
          sort_indices<0,5,2,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 336705ull :
          sort_indices<0,5,2,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 82725ull :
          sort_indices<0,1,4,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 83235ull :
          sort_indices<0,1,4,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 213285ull :
          sort_indices<0,3,4,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 214305ull :
          sort_indices<0,3,4,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 344355ull :
          sort_indices<0,5,4,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 344865ull :
          sort_indices<0,5,4,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2163525ull :
          sort_indices<2,1,0,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2164035ull :
          sort_indices<2,1,0,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2294085ull :
          sort_indices<2,3,0,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2295105ull :
          sort_indices<2,3,0,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2425155ull :
          sort_indices<2,5,0,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2425665ull :
          sort_indices<2,5,0,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2179845ull :
          sort_indices<2,1,4,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2180355ull :
          sort_indices<2,1,4,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2310405ull :
          sort_indices<2,3,4,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2311425ull :
          sort_indices<2,3,4,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2441475ull :
          sort_indices<2,5,4,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2441985ull :
          sort_indices<2,5,4,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4260645ull :
          sort_indices<4,1,0,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4261155ull :
          sort_indices<4,1,0,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4391205ull :
          sort_indices<4,3,0,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4392225ull :
          sort_indices<4,3,0,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4522275ull :
          sort_indices<4,5,0,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4522785ull :
          sort_indices<4,5,0,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4268805ull :
          sort_indices<4,1,2,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4269315ull :
          sort_indices<4,1,2,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4399365ull :
          sort_indices<4,3,2,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4400385ull :
          sort_indices<4,3,2,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4530435ull :
          sort_indices<4,5,2,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4530945ull :
          sort_indices<4,5,2,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        default:
          throw std::logic_error("This case has not been implemented in prim_op_var.h");
    }
  }
  else if (std::abs(-1.0-a) < numerical_zero__ && std::abs(b) < numerical_zero__) {
    switch (tag) {
        case 74565ull :
          sort_indices<0,1,2,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 75075ull :
          sort_indices<0,1,2,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 205125ull :
          sort_indices<0,3,2,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 206145ull :
          sort_indices<0,3,2,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 336195ull :
          sort_indices<0,5,2,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 336705ull :
          sort_indices<0,5,2,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 82725ull :
          sort_indices<0,1,4,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 83235ull :
          sort_indices<0,1,4,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 213285ull :
          sort_indices<0,3,4,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 214305ull :
          sort_indices<0,3,4,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 344355ull :
          sort_indices<0,5,4,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 344865ull :
          sort_indices<0,5,4,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2163525ull :
          sort_indices<2,1,0,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2164035ull :
          sort_indices<2,1,0,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2294085ull :
          sort_indices<2,3,0,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2295105ull :
          sort_indices<2,3,0,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2425155ull :
          sort_indices<2,5,0,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2425665ull :
          sort_indices<2,5,0,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2179845ull :
          sort_indices<2,1,4,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2180355ull :
          sort_indices<2,1,4,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2310405ull :
          sort_indices<2,3,4,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2311425ull :
          sort_indices<2,3,4,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2441475ull :
          sort_indices<2,5,4,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 2441985ull :
          sort_indices<2,5,4,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4260645ull :
          sort_indices<4,1,0,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4261155ull :
          sort_indices<4,1,0,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4391205ull :
          sort_indices<4,3,0,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4392225ull :
          sort_indices<4,3,0,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4522275ull :
          sort_indices<4,5,0,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4522785ull :
          sort_indices<4,5,0,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4268805ull :
          sort_indices<4,1,2,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4269315ull :
          sort_indices<4,1,2,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4399365ull :
          sort_indices<4,3,2,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4400385ull :
          sort_indices<4,3,2,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4530435ull :
          sort_indices<4,5,2,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        case 4530945ull :
          sort_indices<4,5,2,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5]); break;
        default:
          throw std::logic_error("This case has not been implemented in prim_op_var.h");
    }
  }
  else {
    throw std::logic_error("This case has not been implemented in prim_op_var.h");
  }
}

template<typename DataType>
void sort_indices(const std::array<int,8>& o, const double a, const double b, const DataType* in, DataType* out, const std::array<int,8>& d) {
  size_t tag = 0ull;
  for (int i = 0; i != 8; ++i)
    tag = (tag << 4) + o[i];

  if (std::abs(1.0-a) < numerical_zero__ && std::abs(b) < numerical_zero__) {
    switch (tag) {
        case 19088743ull :
          sort_indices<0,1,2,3,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19089253ull :
          sort_indices<0,1,2,3,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19219303ull :
          sort_indices<0,1,2,5,4,3,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19220323ull :
          sort_indices<0,1,2,5,4,7,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19350373ull :
          sort_indices<0,1,2,7,4,3,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19350883ull :
          sort_indices<0,1,2,7,4,5,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52512103ull :
          sort_indices<0,3,2,1,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52512613ull :
          sort_indices<0,3,2,1,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52773223ull :
          sort_indices<0,3,2,5,4,1,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52774753ull :
          sort_indices<0,3,2,5,4,7,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52904293ull :
          sort_indices<0,3,2,7,4,1,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52905313ull :
          sort_indices<0,3,2,7,4,5,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86066023ull :
          sort_indices<0,5,2,1,4,3,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86067043ull :
          sort_indices<0,5,2,1,4,7,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86196583ull :
          sort_indices<0,5,2,3,4,1,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86198113ull :
          sort_indices<0,5,2,3,4,7,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86458723ull :
          sort_indices<0,5,2,7,4,1,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86459233ull :
          sort_indices<0,5,2,7,4,3,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119620453ull :
          sort_indices<0,7,2,1,4,3,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119620963ull :
          sort_indices<0,7,2,1,4,5,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119751013ull :
          sort_indices<0,7,2,3,4,1,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119752033ull :
          sort_indices<0,7,2,3,4,5,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119882083ull :
          sort_indices<0,7,2,5,4,1,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119882593ull :
          sort_indices<0,7,2,5,4,3,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19096903ull :
          sort_indices<0,1,2,3,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19097413ull :
          sort_indices<0,1,2,3,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19227463ull :
          sort_indices<0,1,2,5,6,3,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19228483ull :
          sort_indices<0,1,2,5,6,7,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19358533ull :
          sort_indices<0,1,2,7,6,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19359043ull :
          sort_indices<0,1,2,7,6,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52520263ull :
          sort_indices<0,3,2,1,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52520773ull :
          sort_indices<0,3,2,1,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52781383ull :
          sort_indices<0,3,2,5,6,1,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52782913ull :
          sort_indices<0,3,2,5,6,7,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52912453ull :
          sort_indices<0,3,2,7,6,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52913473ull :
          sort_indices<0,3,2,7,6,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86074183ull :
          sort_indices<0,5,2,1,6,3,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86075203ull :
          sort_indices<0,5,2,1,6,7,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86204743ull :
          sort_indices<0,5,2,3,6,1,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86206273ull :
          sort_indices<0,5,2,3,6,7,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86466883ull :
          sort_indices<0,5,2,7,6,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86467393ull :
          sort_indices<0,5,2,7,6,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119628613ull :
          sort_indices<0,7,2,1,6,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119629123ull :
          sort_indices<0,7,2,1,6,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119759173ull :
          sort_indices<0,7,2,3,6,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119760193ull :
          sort_indices<0,7,2,3,6,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119890243ull :
          sort_indices<0,7,2,5,6,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119890753ull :
          sort_indices<0,7,2,5,6,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21177703ull :
          sort_indices<0,1,4,3,2,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21178213ull :
          sort_indices<0,1,4,3,2,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21308263ull :
          sort_indices<0,1,4,5,2,3,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21309283ull :
          sort_indices<0,1,4,5,2,7,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21439333ull :
          sort_indices<0,1,4,7,2,3,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21439843ull :
          sort_indices<0,1,4,7,2,5,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54601063ull :
          sort_indices<0,3,4,1,2,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54601573ull :
          sort_indices<0,3,4,1,2,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54862183ull :
          sort_indices<0,3,4,5,2,1,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54863713ull :
          sort_indices<0,3,4,5,2,7,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54993253ull :
          sort_indices<0,3,4,7,2,1,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54994273ull :
          sort_indices<0,3,4,7,2,5,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88154983ull :
          sort_indices<0,5,4,1,2,3,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88156003ull :
          sort_indices<0,5,4,1,2,7,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88285543ull :
          sort_indices<0,5,4,3,2,1,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88287073ull :
          sort_indices<0,5,4,3,2,7,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88547683ull :
          sort_indices<0,5,4,7,2,1,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88548193ull :
          sort_indices<0,5,4,7,2,3,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121709413ull :
          sort_indices<0,7,4,1,2,3,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121709923ull :
          sort_indices<0,7,4,1,2,5,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121839973ull :
          sort_indices<0,7,4,3,2,1,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121840993ull :
          sort_indices<0,7,4,3,2,5,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121971043ull :
          sort_indices<0,7,4,5,2,1,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121971553ull :
          sort_indices<0,7,4,5,2,3,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21194023ull :
          sort_indices<0,1,4,3,6,5,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21194533ull :
          sort_indices<0,1,4,3,6,7,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21324583ull :
          sort_indices<0,1,4,5,6,3,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21325603ull :
          sort_indices<0,1,4,5,6,7,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21455653ull :
          sort_indices<0,1,4,7,6,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21456163ull :
          sort_indices<0,1,4,7,6,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54617383ull :
          sort_indices<0,3,4,1,6,5,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54617893ull :
          sort_indices<0,3,4,1,6,7,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54878503ull :
          sort_indices<0,3,4,5,6,1,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54880033ull :
          sort_indices<0,3,4,5,6,7,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 55009573ull :
          sort_indices<0,3,4,7,6,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 55010593ull :
          sort_indices<0,3,4,7,6,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88171303ull :
          sort_indices<0,5,4,1,6,3,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88172323ull :
          sort_indices<0,5,4,1,6,7,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88301863ull :
          sort_indices<0,5,4,3,6,1,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88303393ull :
          sort_indices<0,5,4,3,6,7,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88564003ull :
          sort_indices<0,5,4,7,6,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88564513ull :
          sort_indices<0,5,4,7,6,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121725733ull :
          sort_indices<0,7,4,1,6,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121726243ull :
          sort_indices<0,7,4,1,6,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121856293ull :
          sort_indices<0,7,4,3,6,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121857313ull :
          sort_indices<0,7,4,3,6,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121987363ull :
          sort_indices<0,7,4,5,6,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121987873ull :
          sort_indices<0,7,4,5,6,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23274823ull :
          sort_indices<0,1,6,3,2,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23275333ull :
          sort_indices<0,1,6,3,2,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23405383ull :
          sort_indices<0,1,6,5,2,3,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23406403ull :
          sort_indices<0,1,6,5,2,7,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23536453ull :
          sort_indices<0,1,6,7,2,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23536963ull :
          sort_indices<0,1,6,7,2,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56698183ull :
          sort_indices<0,3,6,1,2,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56698693ull :
          sort_indices<0,3,6,1,2,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56959303ull :
          sort_indices<0,3,6,5,2,1,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56960833ull :
          sort_indices<0,3,6,5,2,7,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 57090373ull :
          sort_indices<0,3,6,7,2,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 57091393ull :
          sort_indices<0,3,6,7,2,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90252103ull :
          sort_indices<0,5,6,1,2,3,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90253123ull :
          sort_indices<0,5,6,1,2,7,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90382663ull :
          sort_indices<0,5,6,3,2,1,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90384193ull :
          sort_indices<0,5,6,3,2,7,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90644803ull :
          sort_indices<0,5,6,7,2,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90645313ull :
          sort_indices<0,5,6,7,2,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123806533ull :
          sort_indices<0,7,6,1,2,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123807043ull :
          sort_indices<0,7,6,1,2,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123937093ull :
          sort_indices<0,7,6,3,2,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123938113ull :
          sort_indices<0,7,6,3,2,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 124068163ull :
          sort_indices<0,7,6,5,2,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 124068673ull :
          sort_indices<0,7,6,5,2,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23282983ull :
          sort_indices<0,1,6,3,4,5,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23283493ull :
          sort_indices<0,1,6,3,4,7,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23413543ull :
          sort_indices<0,1,6,5,4,3,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23414563ull :
          sort_indices<0,1,6,5,4,7,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23544613ull :
          sort_indices<0,1,6,7,4,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23545123ull :
          sort_indices<0,1,6,7,4,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56706343ull :
          sort_indices<0,3,6,1,4,5,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56706853ull :
          sort_indices<0,3,6,1,4,7,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56967463ull :
          sort_indices<0,3,6,5,4,1,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56968993ull :
          sort_indices<0,3,6,5,4,7,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 57098533ull :
          sort_indices<0,3,6,7,4,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 57099553ull :
          sort_indices<0,3,6,7,4,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90260263ull :
          sort_indices<0,5,6,1,4,3,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90261283ull :
          sort_indices<0,5,6,1,4,7,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90390823ull :
          sort_indices<0,5,6,3,4,1,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90392353ull :
          sort_indices<0,5,6,3,4,7,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90652963ull :
          sort_indices<0,5,6,7,4,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90653473ull :
          sort_indices<0,5,6,7,4,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123814693ull :
          sort_indices<0,7,6,1,4,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123815203ull :
          sort_indices<0,7,6,1,4,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123945253ull :
          sort_indices<0,7,6,3,4,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123946273ull :
          sort_indices<0,7,6,3,4,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 124076323ull :
          sort_indices<0,7,6,5,4,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 124076833ull :
          sort_indices<0,7,6,5,4,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 553862503ull :
          sort_indices<2,1,0,3,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 553863013ull :
          sort_indices<2,1,0,3,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 553993063ull :
          sort_indices<2,1,0,5,4,3,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 553994083ull :
          sort_indices<2,1,0,5,4,7,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 554124133ull :
          sort_indices<2,1,0,7,4,3,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 554124643ull :
          sort_indices<2,1,0,7,4,5,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587285863ull :
          sort_indices<2,3,0,1,4,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587286373ull :
          sort_indices<2,3,0,1,4,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587546983ull :
          sort_indices<2,3,0,5,4,1,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587548513ull :
          sort_indices<2,3,0,5,4,7,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587678053ull :
          sort_indices<2,3,0,7,4,1,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587679073ull :
          sort_indices<2,3,0,7,4,5,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620839783ull :
          sort_indices<2,5,0,1,4,3,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620840803ull :
          sort_indices<2,5,0,1,4,7,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620970343ull :
          sort_indices<2,5,0,3,4,1,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620971873ull :
          sort_indices<2,5,0,3,4,7,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 621232483ull :
          sort_indices<2,5,0,7,4,1,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 621232993ull :
          sort_indices<2,5,0,7,4,3,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654394213ull :
          sort_indices<2,7,0,1,4,3,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654394723ull :
          sort_indices<2,7,0,1,4,5,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654524773ull :
          sort_indices<2,7,0,3,4,1,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654525793ull :
          sort_indices<2,7,0,3,4,5,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654655843ull :
          sort_indices<2,7,0,5,4,1,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654656353ull :
          sort_indices<2,7,0,5,4,3,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 553870663ull :
          sort_indices<2,1,0,3,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 553871173ull :
          sort_indices<2,1,0,3,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 554001223ull :
          sort_indices<2,1,0,5,6,3,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 554002243ull :
          sort_indices<2,1,0,5,6,7,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 554132293ull :
          sort_indices<2,1,0,7,6,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 554132803ull :
          sort_indices<2,1,0,7,6,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587294023ull :
          sort_indices<2,3,0,1,6,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587294533ull :
          sort_indices<2,3,0,1,6,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587555143ull :
          sort_indices<2,3,0,5,6,1,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587556673ull :
          sort_indices<2,3,0,5,6,7,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587686213ull :
          sort_indices<2,3,0,7,6,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587687233ull :
          sort_indices<2,3,0,7,6,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620847943ull :
          sort_indices<2,5,0,1,6,3,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620848963ull :
          sort_indices<2,5,0,1,6,7,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620978503ull :
          sort_indices<2,5,0,3,6,1,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620980033ull :
          sort_indices<2,5,0,3,6,7,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 621240643ull :
          sort_indices<2,5,0,7,6,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 621241153ull :
          sort_indices<2,5,0,7,6,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654402373ull :
          sort_indices<2,7,0,1,6,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654402883ull :
          sort_indices<2,7,0,1,6,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654532933ull :
          sort_indices<2,7,0,3,6,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654533953ull :
          sort_indices<2,7,0,3,6,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654664003ull :
          sort_indices<2,7,0,5,6,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654664513ull :
          sort_indices<2,7,0,5,6,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558040423ull :
          sort_indices<2,1,4,3,0,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558040933ull :
          sort_indices<2,1,4,3,0,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558170983ull :
          sort_indices<2,1,4,5,0,3,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558172003ull :
          sort_indices<2,1,4,5,0,7,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558302053ull :
          sort_indices<2,1,4,7,0,3,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558302563ull :
          sort_indices<2,1,4,7,0,5,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591463783ull :
          sort_indices<2,3,4,1,0,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591464293ull :
          sort_indices<2,3,4,1,0,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591724903ull :
          sort_indices<2,3,4,5,0,1,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591726433ull :
          sort_indices<2,3,4,5,0,7,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591855973ull :
          sort_indices<2,3,4,7,0,1,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591856993ull :
          sort_indices<2,3,4,7,0,5,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625017703ull :
          sort_indices<2,5,4,1,0,3,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625018723ull :
          sort_indices<2,5,4,1,0,7,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625148263ull :
          sort_indices<2,5,4,3,0,1,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625149793ull :
          sort_indices<2,5,4,3,0,7,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625410403ull :
          sort_indices<2,5,4,7,0,1,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625410913ull :
          sort_indices<2,5,4,7,0,3,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658572133ull :
          sort_indices<2,7,4,1,0,3,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658572643ull :
          sort_indices<2,7,4,1,0,5,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658702693ull :
          sort_indices<2,7,4,3,0,1,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658703713ull :
          sort_indices<2,7,4,3,0,5,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658833763ull :
          sort_indices<2,7,4,5,0,1,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658834273ull :
          sort_indices<2,7,4,5,0,3,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558064903ull :
          sort_indices<2,1,4,3,6,5,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558065413ull :
          sort_indices<2,1,4,3,6,7,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558195463ull :
          sort_indices<2,1,4,5,6,3,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558196483ull :
          sort_indices<2,1,4,5,6,7,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558326533ull :
          sort_indices<2,1,4,7,6,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558327043ull :
          sort_indices<2,1,4,7,6,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591488263ull :
          sort_indices<2,3,4,1,6,5,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591488773ull :
          sort_indices<2,3,4,1,6,7,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591749383ull :
          sort_indices<2,3,4,5,6,1,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591750913ull :
          sort_indices<2,3,4,5,6,7,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591880453ull :
          sort_indices<2,3,4,7,6,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591881473ull :
          sort_indices<2,3,4,7,6,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625042183ull :
          sort_indices<2,5,4,1,6,3,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625043203ull :
          sort_indices<2,5,4,1,6,7,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625172743ull :
          sort_indices<2,5,4,3,6,1,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625174273ull :
          sort_indices<2,5,4,3,6,7,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625434883ull :
          sort_indices<2,5,4,7,6,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625435393ull :
          sort_indices<2,5,4,7,6,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658596613ull :
          sort_indices<2,7,4,1,6,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658597123ull :
          sort_indices<2,7,4,1,6,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658727173ull :
          sort_indices<2,7,4,3,6,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658728193ull :
          sort_indices<2,7,4,3,6,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658858243ull :
          sort_indices<2,7,4,5,6,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658858753ull :
          sort_indices<2,7,4,5,6,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560137543ull :
          sort_indices<2,1,6,3,0,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560138053ull :
          sort_indices<2,1,6,3,0,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560268103ull :
          sort_indices<2,1,6,5,0,3,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560269123ull :
          sort_indices<2,1,6,5,0,7,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560399173ull :
          sort_indices<2,1,6,7,0,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560399683ull :
          sort_indices<2,1,6,7,0,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593560903ull :
          sort_indices<2,3,6,1,0,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593561413ull :
          sort_indices<2,3,6,1,0,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593822023ull :
          sort_indices<2,3,6,5,0,1,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593823553ull :
          sort_indices<2,3,6,5,0,7,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593953093ull :
          sort_indices<2,3,6,7,0,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593954113ull :
          sort_indices<2,3,6,7,0,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627114823ull :
          sort_indices<2,5,6,1,0,3,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627115843ull :
          sort_indices<2,5,6,1,0,7,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627245383ull :
          sort_indices<2,5,6,3,0,1,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627246913ull :
          sort_indices<2,5,6,3,0,7,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627507523ull :
          sort_indices<2,5,6,7,0,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627508033ull :
          sort_indices<2,5,6,7,0,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660669253ull :
          sort_indices<2,7,6,1,0,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660669763ull :
          sort_indices<2,7,6,1,0,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660799813ull :
          sort_indices<2,7,6,3,0,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660800833ull :
          sort_indices<2,7,6,3,0,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660930883ull :
          sort_indices<2,7,6,5,0,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660931393ull :
          sort_indices<2,7,6,5,0,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560153863ull :
          sort_indices<2,1,6,3,4,5,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560154373ull :
          sort_indices<2,1,6,3,4,7,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560284423ull :
          sort_indices<2,1,6,5,4,3,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560285443ull :
          sort_indices<2,1,6,5,4,7,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560415493ull :
          sort_indices<2,1,6,7,4,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560416003ull :
          sort_indices<2,1,6,7,4,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593577223ull :
          sort_indices<2,3,6,1,4,5,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593577733ull :
          sort_indices<2,3,6,1,4,7,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593838343ull :
          sort_indices<2,3,6,5,4,1,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593839873ull :
          sort_indices<2,3,6,5,4,7,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593969413ull :
          sort_indices<2,3,6,7,4,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593970433ull :
          sort_indices<2,3,6,7,4,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627131143ull :
          sort_indices<2,5,6,1,4,3,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627132163ull :
          sort_indices<2,5,6,1,4,7,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627261703ull :
          sort_indices<2,5,6,3,4,1,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627263233ull :
          sort_indices<2,5,6,3,4,7,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627523843ull :
          sort_indices<2,5,6,7,4,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627524353ull :
          sort_indices<2,5,6,7,4,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660685573ull :
          sort_indices<2,7,6,1,4,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660686083ull :
          sort_indices<2,7,6,1,4,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660816133ull :
          sort_indices<2,7,6,3,4,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660817153ull :
          sort_indices<2,7,6,3,4,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660947203ull :
          sort_indices<2,7,6,5,4,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660947713ull :
          sort_indices<2,7,6,5,4,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090725223ull :
          sort_indices<4,1,0,3,2,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090725733ull :
          sort_indices<4,1,0,3,2,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090855783ull :
          sort_indices<4,1,0,5,2,3,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090856803ull :
          sort_indices<4,1,0,5,2,7,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090986853ull :
          sort_indices<4,1,0,7,2,3,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090987363ull :
          sort_indices<4,1,0,7,2,5,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124148583ull :
          sort_indices<4,3,0,1,2,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124149093ull :
          sort_indices<4,3,0,1,2,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124409703ull :
          sort_indices<4,3,0,5,2,1,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124411233ull :
          sort_indices<4,3,0,5,2,7,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124540773ull :
          sort_indices<4,3,0,7,2,1,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124541793ull :
          sort_indices<4,3,0,7,2,5,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157702503ull :
          sort_indices<4,5,0,1,2,3,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157703523ull :
          sort_indices<4,5,0,1,2,7,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157833063ull :
          sort_indices<4,5,0,3,2,1,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157834593ull :
          sort_indices<4,5,0,3,2,7,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1158095203ull :
          sort_indices<4,5,0,7,2,1,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1158095713ull :
          sort_indices<4,5,0,7,2,3,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191256933ull :
          sort_indices<4,7,0,1,2,3,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191257443ull :
          sort_indices<4,7,0,1,2,5,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191387493ull :
          sort_indices<4,7,0,3,2,1,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191388513ull :
          sort_indices<4,7,0,3,2,5,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191518563ull :
          sort_indices<4,7,0,5,2,1,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191519073ull :
          sort_indices<4,7,0,5,2,3,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090741543ull :
          sort_indices<4,1,0,3,6,5,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090742053ull :
          sort_indices<4,1,0,3,6,7,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090872103ull :
          sort_indices<4,1,0,5,6,3,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090873123ull :
          sort_indices<4,1,0,5,6,7,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1091003173ull :
          sort_indices<4,1,0,7,6,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1091003683ull :
          sort_indices<4,1,0,7,6,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124164903ull :
          sort_indices<4,3,0,1,6,5,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124165413ull :
          sort_indices<4,3,0,1,6,7,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124426023ull :
          sort_indices<4,3,0,5,6,1,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124427553ull :
          sort_indices<4,3,0,5,6,7,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124557093ull :
          sort_indices<4,3,0,7,6,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124558113ull :
          sort_indices<4,3,0,7,6,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157718823ull :
          sort_indices<4,5,0,1,6,3,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157719843ull :
          sort_indices<4,5,0,1,6,7,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157849383ull :
          sort_indices<4,5,0,3,6,1,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157850913ull :
          sort_indices<4,5,0,3,6,7,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1158111523ull :
          sort_indices<4,5,0,7,6,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1158112033ull :
          sort_indices<4,5,0,7,6,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191273253ull :
          sort_indices<4,7,0,1,6,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191273763ull :
          sort_indices<4,7,0,1,6,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191403813ull :
          sort_indices<4,7,0,3,6,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191404833ull :
          sort_indices<4,7,0,3,6,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191534883ull :
          sort_indices<4,7,0,5,6,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191535393ull :
          sort_indices<4,7,0,5,6,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092814183ull :
          sort_indices<4,1,2,3,0,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092814693ull :
          sort_indices<4,1,2,3,0,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092944743ull :
          sort_indices<4,1,2,5,0,3,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092945763ull :
          sort_indices<4,1,2,5,0,7,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1093075813ull :
          sort_indices<4,1,2,7,0,3,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1093076323ull :
          sort_indices<4,1,2,7,0,5,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126237543ull :
          sort_indices<4,3,2,1,0,5,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126238053ull :
          sort_indices<4,3,2,1,0,7,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126498663ull :
          sort_indices<4,3,2,5,0,1,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126500193ull :
          sort_indices<4,3,2,5,0,7,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126629733ull :
          sort_indices<4,3,2,7,0,1,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126630753ull :
          sort_indices<4,3,2,7,0,5,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159791463ull :
          sort_indices<4,5,2,1,0,3,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159792483ull :
          sort_indices<4,5,2,1,0,7,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159922023ull :
          sort_indices<4,5,2,3,0,1,6,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159923553ull :
          sort_indices<4,5,2,3,0,7,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1160184163ull :
          sort_indices<4,5,2,7,0,1,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1160184673ull :
          sort_indices<4,5,2,7,0,3,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193345893ull :
          sort_indices<4,7,2,1,0,3,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193346403ull :
          sort_indices<4,7,2,1,0,5,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193476453ull :
          sort_indices<4,7,2,3,0,1,6,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193477473ull :
          sort_indices<4,7,2,3,0,5,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193607523ull :
          sort_indices<4,7,2,5,0,1,6,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193608033ull :
          sort_indices<4,7,2,5,0,3,6,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092838663ull :
          sort_indices<4,1,2,3,6,5,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092839173ull :
          sort_indices<4,1,2,3,6,7,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092969223ull :
          sort_indices<4,1,2,5,6,3,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092970243ull :
          sort_indices<4,1,2,5,6,7,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1093100293ull :
          sort_indices<4,1,2,7,6,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1093100803ull :
          sort_indices<4,1,2,7,6,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126262023ull :
          sort_indices<4,3,2,1,6,5,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126262533ull :
          sort_indices<4,3,2,1,6,7,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126523143ull :
          sort_indices<4,3,2,5,6,1,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126524673ull :
          sort_indices<4,3,2,5,6,7,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126654213ull :
          sort_indices<4,3,2,7,6,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126655233ull :
          sort_indices<4,3,2,7,6,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159815943ull :
          sort_indices<4,5,2,1,6,3,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159816963ull :
          sort_indices<4,5,2,1,6,7,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159946503ull :
          sort_indices<4,5,2,3,6,1,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159948033ull :
          sort_indices<4,5,2,3,6,7,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1160208643ull :
          sort_indices<4,5,2,7,6,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1160209153ull :
          sort_indices<4,5,2,7,6,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193370373ull :
          sort_indices<4,7,2,1,6,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193370883ull :
          sort_indices<4,7,2,1,6,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193500933ull :
          sort_indices<4,7,2,3,6,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193501953ull :
          sort_indices<4,7,2,3,6,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193632003ull :
          sort_indices<4,7,2,5,6,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193632513ull :
          sort_indices<4,7,2,5,6,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097008423ull :
          sort_indices<4,1,6,3,0,5,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097008933ull :
          sort_indices<4,1,6,3,0,7,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097138983ull :
          sort_indices<4,1,6,5,0,3,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097140003ull :
          sort_indices<4,1,6,5,0,7,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097270053ull :
          sort_indices<4,1,6,7,0,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097270563ull :
          sort_indices<4,1,6,7,0,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130431783ull :
          sort_indices<4,3,6,1,0,5,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130432293ull :
          sort_indices<4,3,6,1,0,7,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130692903ull :
          sort_indices<4,3,6,5,0,1,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130694433ull :
          sort_indices<4,3,6,5,0,7,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130823973ull :
          sort_indices<4,3,6,7,0,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130824993ull :
          sort_indices<4,3,6,7,0,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1163985703ull :
          sort_indices<4,5,6,1,0,3,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1163986723ull :
          sort_indices<4,5,6,1,0,7,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164116263ull :
          sort_indices<4,5,6,3,0,1,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164117793ull :
          sort_indices<4,5,6,3,0,7,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164378403ull :
          sort_indices<4,5,6,7,0,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164378913ull :
          sort_indices<4,5,6,7,0,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197540133ull :
          sort_indices<4,7,6,1,0,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197540643ull :
          sort_indices<4,7,6,1,0,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197670693ull :
          sort_indices<4,7,6,3,0,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197671713ull :
          sort_indices<4,7,6,3,0,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197801763ull :
          sort_indices<4,7,6,5,0,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197802273ull :
          sort_indices<4,7,6,5,0,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097016583ull :
          sort_indices<4,1,6,3,2,5,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097017093ull :
          sort_indices<4,1,6,3,2,7,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097147143ull :
          sort_indices<4,1,6,5,2,3,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097148163ull :
          sort_indices<4,1,6,5,2,7,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097278213ull :
          sort_indices<4,1,6,7,2,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097278723ull :
          sort_indices<4,1,6,7,2,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130439943ull :
          sort_indices<4,3,6,1,2,5,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130440453ull :
          sort_indices<4,3,6,1,2,7,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130701063ull :
          sort_indices<4,3,6,5,2,1,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130702593ull :
          sort_indices<4,3,6,5,2,7,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130832133ull :
          sort_indices<4,3,6,7,2,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130833153ull :
          sort_indices<4,3,6,7,2,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1163993863ull :
          sort_indices<4,5,6,1,2,3,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1163994883ull :
          sort_indices<4,5,6,1,2,7,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164124423ull :
          sort_indices<4,5,6,3,2,1,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164125953ull :
          sort_indices<4,5,6,3,2,7,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164386563ull :
          sort_indices<4,5,6,7,2,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164387073ull :
          sort_indices<4,5,6,7,2,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197548293ull :
          sort_indices<4,7,6,1,2,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197548803ull :
          sort_indices<4,7,6,1,2,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197678853ull :
          sort_indices<4,7,6,3,2,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197679873ull :
          sort_indices<4,7,6,3,2,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197809923ull :
          sort_indices<4,7,6,5,2,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197810433ull :
          sort_indices<4,7,6,5,2,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627596103ull :
          sort_indices<6,1,0,3,2,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627596613ull :
          sort_indices<6,1,0,3,2,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627726663ull :
          sort_indices<6,1,0,5,2,3,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627727683ull :
          sort_indices<6,1,0,5,2,7,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627857733ull :
          sort_indices<6,1,0,7,2,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627858243ull :
          sort_indices<6,1,0,7,2,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661019463ull :
          sort_indices<6,3,0,1,2,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661019973ull :
          sort_indices<6,3,0,1,2,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661280583ull :
          sort_indices<6,3,0,5,2,1,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661282113ull :
          sort_indices<6,3,0,5,2,7,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661411653ull :
          sort_indices<6,3,0,7,2,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661412673ull :
          sort_indices<6,3,0,7,2,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694573383ull :
          sort_indices<6,5,0,1,2,3,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694574403ull :
          sort_indices<6,5,0,1,2,7,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694703943ull :
          sort_indices<6,5,0,3,2,1,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694705473ull :
          sort_indices<6,5,0,3,2,7,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694966083ull :
          sort_indices<6,5,0,7,2,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694966593ull :
          sort_indices<6,5,0,7,2,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728127813ull :
          sort_indices<6,7,0,1,2,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728128323ull :
          sort_indices<6,7,0,1,2,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728258373ull :
          sort_indices<6,7,0,3,2,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728259393ull :
          sort_indices<6,7,0,3,2,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728389443ull :
          sort_indices<6,7,0,5,2,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728389953ull :
          sort_indices<6,7,0,5,2,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627604263ull :
          sort_indices<6,1,0,3,4,5,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627604773ull :
          sort_indices<6,1,0,3,4,7,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627734823ull :
          sort_indices<6,1,0,5,4,3,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627735843ull :
          sort_indices<6,1,0,5,4,7,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627865893ull :
          sort_indices<6,1,0,7,4,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627866403ull :
          sort_indices<6,1,0,7,4,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661027623ull :
          sort_indices<6,3,0,1,4,5,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661028133ull :
          sort_indices<6,3,0,1,4,7,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661288743ull :
          sort_indices<6,3,0,5,4,1,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661290273ull :
          sort_indices<6,3,0,5,4,7,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661419813ull :
          sort_indices<6,3,0,7,4,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661420833ull :
          sort_indices<6,3,0,7,4,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694581543ull :
          sort_indices<6,5,0,1,4,3,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694582563ull :
          sort_indices<6,5,0,1,4,7,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694712103ull :
          sort_indices<6,5,0,3,4,1,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694713633ull :
          sort_indices<6,5,0,3,4,7,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694974243ull :
          sort_indices<6,5,0,7,4,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694974753ull :
          sort_indices<6,5,0,7,4,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728135973ull :
          sort_indices<6,7,0,1,4,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728136483ull :
          sort_indices<6,7,0,1,4,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728266533ull :
          sort_indices<6,7,0,3,4,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728267553ull :
          sort_indices<6,7,0,3,4,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728397603ull :
          sort_indices<6,7,0,5,4,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728398113ull :
          sort_indices<6,7,0,5,4,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629685063ull :
          sort_indices<6,1,2,3,0,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629685573ull :
          sort_indices<6,1,2,3,0,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629815623ull :
          sort_indices<6,1,2,5,0,3,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629816643ull :
          sort_indices<6,1,2,5,0,7,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629946693ull :
          sort_indices<6,1,2,7,0,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629947203ull :
          sort_indices<6,1,2,7,0,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663108423ull :
          sort_indices<6,3,2,1,0,5,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663108933ull :
          sort_indices<6,3,2,1,0,7,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663369543ull :
          sort_indices<6,3,2,5,0,1,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663371073ull :
          sort_indices<6,3,2,5,0,7,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663500613ull :
          sort_indices<6,3,2,7,0,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663501633ull :
          sort_indices<6,3,2,7,0,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696662343ull :
          sort_indices<6,5,2,1,0,3,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696663363ull :
          sort_indices<6,5,2,1,0,7,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696792903ull :
          sort_indices<6,5,2,3,0,1,4,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696794433ull :
          sort_indices<6,5,2,3,0,7,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1697055043ull :
          sort_indices<6,5,2,7,0,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1697055553ull :
          sort_indices<6,5,2,7,0,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730216773ull :
          sort_indices<6,7,2,1,0,3,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730217283ull :
          sort_indices<6,7,2,1,0,5,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730347333ull :
          sort_indices<6,7,2,3,0,1,4,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730348353ull :
          sort_indices<6,7,2,3,0,5,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730478403ull :
          sort_indices<6,7,2,5,0,1,4,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730478913ull :
          sort_indices<6,7,2,5,0,3,4,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629701383ull :
          sort_indices<6,1,2,3,4,5,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629701893ull :
          sort_indices<6,1,2,3,4,7,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629831943ull :
          sort_indices<6,1,2,5,4,3,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629832963ull :
          sort_indices<6,1,2,5,4,7,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629963013ull :
          sort_indices<6,1,2,7,4,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629963523ull :
          sort_indices<6,1,2,7,4,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663124743ull :
          sort_indices<6,3,2,1,4,5,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663125253ull :
          sort_indices<6,3,2,1,4,7,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663385863ull :
          sort_indices<6,3,2,5,4,1,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663387393ull :
          sort_indices<6,3,2,5,4,7,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663516933ull :
          sort_indices<6,3,2,7,4,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663517953ull :
          sort_indices<6,3,2,7,4,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696678663ull :
          sort_indices<6,5,2,1,4,3,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696679683ull :
          sort_indices<6,5,2,1,4,7,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696809223ull :
          sort_indices<6,5,2,3,4,1,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696810753ull :
          sort_indices<6,5,2,3,4,7,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1697071363ull :
          sort_indices<6,5,2,7,4,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1697071873ull :
          sort_indices<6,5,2,7,4,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730233093ull :
          sort_indices<6,7,2,1,4,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730233603ull :
          sort_indices<6,7,2,1,4,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730363653ull :
          sort_indices<6,7,2,3,4,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730364673ull :
          sort_indices<6,7,2,3,4,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730494723ull :
          sort_indices<6,7,2,5,4,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730495233ull :
          sort_indices<6,7,2,5,4,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631782183ull :
          sort_indices<6,1,4,3,0,5,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631782693ull :
          sort_indices<6,1,4,3,0,7,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631912743ull :
          sort_indices<6,1,4,5,0,3,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631913763ull :
          sort_indices<6,1,4,5,0,7,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1632043813ull :
          sort_indices<6,1,4,7,0,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1632044323ull :
          sort_indices<6,1,4,7,0,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665205543ull :
          sort_indices<6,3,4,1,0,5,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665206053ull :
          sort_indices<6,3,4,1,0,7,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665466663ull :
          sort_indices<6,3,4,5,0,1,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665468193ull :
          sort_indices<6,3,4,5,0,7,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665597733ull :
          sort_indices<6,3,4,7,0,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665598753ull :
          sort_indices<6,3,4,7,0,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698759463ull :
          sort_indices<6,5,4,1,0,3,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698760483ull :
          sort_indices<6,5,4,1,0,7,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698890023ull :
          sort_indices<6,5,4,3,0,1,2,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698891553ull :
          sort_indices<6,5,4,3,0,7,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1699152163ull :
          sort_indices<6,5,4,7,0,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1699152673ull :
          sort_indices<6,5,4,7,0,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732313893ull :
          sort_indices<6,7,4,1,0,3,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732314403ull :
          sort_indices<6,7,4,1,0,5,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732444453ull :
          sort_indices<6,7,4,3,0,1,2,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732445473ull :
          sort_indices<6,7,4,3,0,5,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732575523ull :
          sort_indices<6,7,4,5,0,1,2,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732576033ull :
          sort_indices<6,7,4,5,0,3,2,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631790343ull :
          sort_indices<6,1,4,3,2,5,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631790853ull :
          sort_indices<6,1,4,3,2,7,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631920903ull :
          sort_indices<6,1,4,5,2,3,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631921923ull :
          sort_indices<6,1,4,5,2,7,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1632051973ull :
          sort_indices<6,1,4,7,2,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1632052483ull :
          sort_indices<6,1,4,7,2,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665213703ull :
          sort_indices<6,3,4,1,2,5,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665214213ull :
          sort_indices<6,3,4,1,2,7,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665474823ull :
          sort_indices<6,3,4,5,2,1,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665476353ull :
          sort_indices<6,3,4,5,2,7,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665605893ull :
          sort_indices<6,3,4,7,2,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665606913ull :
          sort_indices<6,3,4,7,2,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698767623ull :
          sort_indices<6,5,4,1,2,3,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698768643ull :
          sort_indices<6,5,4,1,2,7,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698898183ull :
          sort_indices<6,5,4,3,2,1,0,7,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698899713ull :
          sort_indices<6,5,4,3,2,7,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1699160323ull :
          sort_indices<6,5,4,7,2,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1699160833ull :
          sort_indices<6,5,4,7,2,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732322053ull :
          sort_indices<6,7,4,1,2,3,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732322563ull :
          sort_indices<6,7,4,1,2,5,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732452613ull :
          sort_indices<6,7,4,3,2,1,0,5,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732453633ull :
          sort_indices<6,7,4,3,2,5,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732583683ull :
          sort_indices<6,7,4,5,2,1,0,3,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732584193ull :
          sort_indices<6,7,4,5,2,3,0,1,0,1,1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        default:
          throw std::logic_error("This case has not been implemented in prim_op_var.h");
    }
  }
  else if (std::abs(-1.0-a) < numerical_zero__ && std::abs(b) < numerical_zero__) {
    switch (tag) {
        case 19088743ull :
          sort_indices<0,1,2,3,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19089253ull :
          sort_indices<0,1,2,3,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19219303ull :
          sort_indices<0,1,2,5,4,3,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19220323ull :
          sort_indices<0,1,2,5,4,7,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19350373ull :
          sort_indices<0,1,2,7,4,3,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19350883ull :
          sort_indices<0,1,2,7,4,5,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52512103ull :
          sort_indices<0,3,2,1,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52512613ull :
          sort_indices<0,3,2,1,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52773223ull :
          sort_indices<0,3,2,5,4,1,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52774753ull :
          sort_indices<0,3,2,5,4,7,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52904293ull :
          sort_indices<0,3,2,7,4,1,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52905313ull :
          sort_indices<0,3,2,7,4,5,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86066023ull :
          sort_indices<0,5,2,1,4,3,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86067043ull :
          sort_indices<0,5,2,1,4,7,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86196583ull :
          sort_indices<0,5,2,3,4,1,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86198113ull :
          sort_indices<0,5,2,3,4,7,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86458723ull :
          sort_indices<0,5,2,7,4,1,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86459233ull :
          sort_indices<0,5,2,7,4,3,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119620453ull :
          sort_indices<0,7,2,1,4,3,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119620963ull :
          sort_indices<0,7,2,1,4,5,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119751013ull :
          sort_indices<0,7,2,3,4,1,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119752033ull :
          sort_indices<0,7,2,3,4,5,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119882083ull :
          sort_indices<0,7,2,5,4,1,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119882593ull :
          sort_indices<0,7,2,5,4,3,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19096903ull :
          sort_indices<0,1,2,3,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19097413ull :
          sort_indices<0,1,2,3,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19227463ull :
          sort_indices<0,1,2,5,6,3,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19228483ull :
          sort_indices<0,1,2,5,6,7,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19358533ull :
          sort_indices<0,1,2,7,6,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 19359043ull :
          sort_indices<0,1,2,7,6,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52520263ull :
          sort_indices<0,3,2,1,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52520773ull :
          sort_indices<0,3,2,1,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52781383ull :
          sort_indices<0,3,2,5,6,1,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52782913ull :
          sort_indices<0,3,2,5,6,7,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52912453ull :
          sort_indices<0,3,2,7,6,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 52913473ull :
          sort_indices<0,3,2,7,6,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86074183ull :
          sort_indices<0,5,2,1,6,3,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86075203ull :
          sort_indices<0,5,2,1,6,7,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86204743ull :
          sort_indices<0,5,2,3,6,1,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86206273ull :
          sort_indices<0,5,2,3,6,7,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86466883ull :
          sort_indices<0,5,2,7,6,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 86467393ull :
          sort_indices<0,5,2,7,6,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119628613ull :
          sort_indices<0,7,2,1,6,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119629123ull :
          sort_indices<0,7,2,1,6,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119759173ull :
          sort_indices<0,7,2,3,6,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119760193ull :
          sort_indices<0,7,2,3,6,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119890243ull :
          sort_indices<0,7,2,5,6,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 119890753ull :
          sort_indices<0,7,2,5,6,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21177703ull :
          sort_indices<0,1,4,3,2,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21178213ull :
          sort_indices<0,1,4,3,2,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21308263ull :
          sort_indices<0,1,4,5,2,3,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21309283ull :
          sort_indices<0,1,4,5,2,7,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21439333ull :
          sort_indices<0,1,4,7,2,3,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21439843ull :
          sort_indices<0,1,4,7,2,5,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54601063ull :
          sort_indices<0,3,4,1,2,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54601573ull :
          sort_indices<0,3,4,1,2,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54862183ull :
          sort_indices<0,3,4,5,2,1,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54863713ull :
          sort_indices<0,3,4,5,2,7,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54993253ull :
          sort_indices<0,3,4,7,2,1,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54994273ull :
          sort_indices<0,3,4,7,2,5,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88154983ull :
          sort_indices<0,5,4,1,2,3,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88156003ull :
          sort_indices<0,5,4,1,2,7,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88285543ull :
          sort_indices<0,5,4,3,2,1,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88287073ull :
          sort_indices<0,5,4,3,2,7,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88547683ull :
          sort_indices<0,5,4,7,2,1,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88548193ull :
          sort_indices<0,5,4,7,2,3,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121709413ull :
          sort_indices<0,7,4,1,2,3,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121709923ull :
          sort_indices<0,7,4,1,2,5,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121839973ull :
          sort_indices<0,7,4,3,2,1,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121840993ull :
          sort_indices<0,7,4,3,2,5,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121971043ull :
          sort_indices<0,7,4,5,2,1,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121971553ull :
          sort_indices<0,7,4,5,2,3,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21194023ull :
          sort_indices<0,1,4,3,6,5,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21194533ull :
          sort_indices<0,1,4,3,6,7,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21324583ull :
          sort_indices<0,1,4,5,6,3,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21325603ull :
          sort_indices<0,1,4,5,6,7,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21455653ull :
          sort_indices<0,1,4,7,6,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 21456163ull :
          sort_indices<0,1,4,7,6,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54617383ull :
          sort_indices<0,3,4,1,6,5,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54617893ull :
          sort_indices<0,3,4,1,6,7,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54878503ull :
          sort_indices<0,3,4,5,6,1,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 54880033ull :
          sort_indices<0,3,4,5,6,7,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 55009573ull :
          sort_indices<0,3,4,7,6,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 55010593ull :
          sort_indices<0,3,4,7,6,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88171303ull :
          sort_indices<0,5,4,1,6,3,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88172323ull :
          sort_indices<0,5,4,1,6,7,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88301863ull :
          sort_indices<0,5,4,3,6,1,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88303393ull :
          sort_indices<0,5,4,3,6,7,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88564003ull :
          sort_indices<0,5,4,7,6,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 88564513ull :
          sort_indices<0,5,4,7,6,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121725733ull :
          sort_indices<0,7,4,1,6,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121726243ull :
          sort_indices<0,7,4,1,6,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121856293ull :
          sort_indices<0,7,4,3,6,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121857313ull :
          sort_indices<0,7,4,3,6,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121987363ull :
          sort_indices<0,7,4,5,6,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 121987873ull :
          sort_indices<0,7,4,5,6,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23274823ull :
          sort_indices<0,1,6,3,2,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23275333ull :
          sort_indices<0,1,6,3,2,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23405383ull :
          sort_indices<0,1,6,5,2,3,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23406403ull :
          sort_indices<0,1,6,5,2,7,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23536453ull :
          sort_indices<0,1,6,7,2,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23536963ull :
          sort_indices<0,1,6,7,2,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56698183ull :
          sort_indices<0,3,6,1,2,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56698693ull :
          sort_indices<0,3,6,1,2,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56959303ull :
          sort_indices<0,3,6,5,2,1,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56960833ull :
          sort_indices<0,3,6,5,2,7,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 57090373ull :
          sort_indices<0,3,6,7,2,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 57091393ull :
          sort_indices<0,3,6,7,2,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90252103ull :
          sort_indices<0,5,6,1,2,3,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90253123ull :
          sort_indices<0,5,6,1,2,7,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90382663ull :
          sort_indices<0,5,6,3,2,1,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90384193ull :
          sort_indices<0,5,6,3,2,7,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90644803ull :
          sort_indices<0,5,6,7,2,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90645313ull :
          sort_indices<0,5,6,7,2,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123806533ull :
          sort_indices<0,7,6,1,2,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123807043ull :
          sort_indices<0,7,6,1,2,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123937093ull :
          sort_indices<0,7,6,3,2,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123938113ull :
          sort_indices<0,7,6,3,2,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 124068163ull :
          sort_indices<0,7,6,5,2,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 124068673ull :
          sort_indices<0,7,6,5,2,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23282983ull :
          sort_indices<0,1,6,3,4,5,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23283493ull :
          sort_indices<0,1,6,3,4,7,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23413543ull :
          sort_indices<0,1,6,5,4,3,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23414563ull :
          sort_indices<0,1,6,5,4,7,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23544613ull :
          sort_indices<0,1,6,7,4,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 23545123ull :
          sort_indices<0,1,6,7,4,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56706343ull :
          sort_indices<0,3,6,1,4,5,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56706853ull :
          sort_indices<0,3,6,1,4,7,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56967463ull :
          sort_indices<0,3,6,5,4,1,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 56968993ull :
          sort_indices<0,3,6,5,4,7,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 57098533ull :
          sort_indices<0,3,6,7,4,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 57099553ull :
          sort_indices<0,3,6,7,4,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90260263ull :
          sort_indices<0,5,6,1,4,3,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90261283ull :
          sort_indices<0,5,6,1,4,7,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90390823ull :
          sort_indices<0,5,6,3,4,1,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90392353ull :
          sort_indices<0,5,6,3,4,7,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90652963ull :
          sort_indices<0,5,6,7,4,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 90653473ull :
          sort_indices<0,5,6,7,4,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123814693ull :
          sort_indices<0,7,6,1,4,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123815203ull :
          sort_indices<0,7,6,1,4,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123945253ull :
          sort_indices<0,7,6,3,4,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 123946273ull :
          sort_indices<0,7,6,3,4,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 124076323ull :
          sort_indices<0,7,6,5,4,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 124076833ull :
          sort_indices<0,7,6,5,4,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 553862503ull :
          sort_indices<2,1,0,3,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 553863013ull :
          sort_indices<2,1,0,3,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 553993063ull :
          sort_indices<2,1,0,5,4,3,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 553994083ull :
          sort_indices<2,1,0,5,4,7,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 554124133ull :
          sort_indices<2,1,0,7,4,3,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 554124643ull :
          sort_indices<2,1,0,7,4,5,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587285863ull :
          sort_indices<2,3,0,1,4,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587286373ull :
          sort_indices<2,3,0,1,4,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587546983ull :
          sort_indices<2,3,0,5,4,1,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587548513ull :
          sort_indices<2,3,0,5,4,7,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587678053ull :
          sort_indices<2,3,0,7,4,1,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587679073ull :
          sort_indices<2,3,0,7,4,5,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620839783ull :
          sort_indices<2,5,0,1,4,3,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620840803ull :
          sort_indices<2,5,0,1,4,7,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620970343ull :
          sort_indices<2,5,0,3,4,1,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620971873ull :
          sort_indices<2,5,0,3,4,7,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 621232483ull :
          sort_indices<2,5,0,7,4,1,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 621232993ull :
          sort_indices<2,5,0,7,4,3,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654394213ull :
          sort_indices<2,7,0,1,4,3,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654394723ull :
          sort_indices<2,7,0,1,4,5,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654524773ull :
          sort_indices<2,7,0,3,4,1,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654525793ull :
          sort_indices<2,7,0,3,4,5,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654655843ull :
          sort_indices<2,7,0,5,4,1,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654656353ull :
          sort_indices<2,7,0,5,4,3,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 553870663ull :
          sort_indices<2,1,0,3,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 553871173ull :
          sort_indices<2,1,0,3,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 554001223ull :
          sort_indices<2,1,0,5,6,3,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 554002243ull :
          sort_indices<2,1,0,5,6,7,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 554132293ull :
          sort_indices<2,1,0,7,6,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 554132803ull :
          sort_indices<2,1,0,7,6,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587294023ull :
          sort_indices<2,3,0,1,6,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587294533ull :
          sort_indices<2,3,0,1,6,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587555143ull :
          sort_indices<2,3,0,5,6,1,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587556673ull :
          sort_indices<2,3,0,5,6,7,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587686213ull :
          sort_indices<2,3,0,7,6,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 587687233ull :
          sort_indices<2,3,0,7,6,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620847943ull :
          sort_indices<2,5,0,1,6,3,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620848963ull :
          sort_indices<2,5,0,1,6,7,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620978503ull :
          sort_indices<2,5,0,3,6,1,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 620980033ull :
          sort_indices<2,5,0,3,6,7,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 621240643ull :
          sort_indices<2,5,0,7,6,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 621241153ull :
          sort_indices<2,5,0,7,6,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654402373ull :
          sort_indices<2,7,0,1,6,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654402883ull :
          sort_indices<2,7,0,1,6,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654532933ull :
          sort_indices<2,7,0,3,6,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654533953ull :
          sort_indices<2,7,0,3,6,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654664003ull :
          sort_indices<2,7,0,5,6,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 654664513ull :
          sort_indices<2,7,0,5,6,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558040423ull :
          sort_indices<2,1,4,3,0,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558040933ull :
          sort_indices<2,1,4,3,0,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558170983ull :
          sort_indices<2,1,4,5,0,3,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558172003ull :
          sort_indices<2,1,4,5,0,7,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558302053ull :
          sort_indices<2,1,4,7,0,3,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558302563ull :
          sort_indices<2,1,4,7,0,5,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591463783ull :
          sort_indices<2,3,4,1,0,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591464293ull :
          sort_indices<2,3,4,1,0,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591724903ull :
          sort_indices<2,3,4,5,0,1,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591726433ull :
          sort_indices<2,3,4,5,0,7,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591855973ull :
          sort_indices<2,3,4,7,0,1,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591856993ull :
          sort_indices<2,3,4,7,0,5,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625017703ull :
          sort_indices<2,5,4,1,0,3,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625018723ull :
          sort_indices<2,5,4,1,0,7,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625148263ull :
          sort_indices<2,5,4,3,0,1,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625149793ull :
          sort_indices<2,5,4,3,0,7,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625410403ull :
          sort_indices<2,5,4,7,0,1,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625410913ull :
          sort_indices<2,5,4,7,0,3,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658572133ull :
          sort_indices<2,7,4,1,0,3,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658572643ull :
          sort_indices<2,7,4,1,0,5,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658702693ull :
          sort_indices<2,7,4,3,0,1,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658703713ull :
          sort_indices<2,7,4,3,0,5,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658833763ull :
          sort_indices<2,7,4,5,0,1,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658834273ull :
          sort_indices<2,7,4,5,0,3,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558064903ull :
          sort_indices<2,1,4,3,6,5,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558065413ull :
          sort_indices<2,1,4,3,6,7,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558195463ull :
          sort_indices<2,1,4,5,6,3,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558196483ull :
          sort_indices<2,1,4,5,6,7,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558326533ull :
          sort_indices<2,1,4,7,6,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 558327043ull :
          sort_indices<2,1,4,7,6,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591488263ull :
          sort_indices<2,3,4,1,6,5,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591488773ull :
          sort_indices<2,3,4,1,6,7,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591749383ull :
          sort_indices<2,3,4,5,6,1,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591750913ull :
          sort_indices<2,3,4,5,6,7,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591880453ull :
          sort_indices<2,3,4,7,6,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 591881473ull :
          sort_indices<2,3,4,7,6,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625042183ull :
          sort_indices<2,5,4,1,6,3,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625043203ull :
          sort_indices<2,5,4,1,6,7,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625172743ull :
          sort_indices<2,5,4,3,6,1,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625174273ull :
          sort_indices<2,5,4,3,6,7,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625434883ull :
          sort_indices<2,5,4,7,6,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 625435393ull :
          sort_indices<2,5,4,7,6,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658596613ull :
          sort_indices<2,7,4,1,6,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658597123ull :
          sort_indices<2,7,4,1,6,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658727173ull :
          sort_indices<2,7,4,3,6,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658728193ull :
          sort_indices<2,7,4,3,6,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658858243ull :
          sort_indices<2,7,4,5,6,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 658858753ull :
          sort_indices<2,7,4,5,6,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560137543ull :
          sort_indices<2,1,6,3,0,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560138053ull :
          sort_indices<2,1,6,3,0,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560268103ull :
          sort_indices<2,1,6,5,0,3,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560269123ull :
          sort_indices<2,1,6,5,0,7,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560399173ull :
          sort_indices<2,1,6,7,0,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560399683ull :
          sort_indices<2,1,6,7,0,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593560903ull :
          sort_indices<2,3,6,1,0,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593561413ull :
          sort_indices<2,3,6,1,0,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593822023ull :
          sort_indices<2,3,6,5,0,1,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593823553ull :
          sort_indices<2,3,6,5,0,7,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593953093ull :
          sort_indices<2,3,6,7,0,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593954113ull :
          sort_indices<2,3,6,7,0,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627114823ull :
          sort_indices<2,5,6,1,0,3,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627115843ull :
          sort_indices<2,5,6,1,0,7,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627245383ull :
          sort_indices<2,5,6,3,0,1,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627246913ull :
          sort_indices<2,5,6,3,0,7,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627507523ull :
          sort_indices<2,5,6,7,0,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627508033ull :
          sort_indices<2,5,6,7,0,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660669253ull :
          sort_indices<2,7,6,1,0,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660669763ull :
          sort_indices<2,7,6,1,0,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660799813ull :
          sort_indices<2,7,6,3,0,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660800833ull :
          sort_indices<2,7,6,3,0,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660930883ull :
          sort_indices<2,7,6,5,0,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660931393ull :
          sort_indices<2,7,6,5,0,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560153863ull :
          sort_indices<2,1,6,3,4,5,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560154373ull :
          sort_indices<2,1,6,3,4,7,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560284423ull :
          sort_indices<2,1,6,5,4,3,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560285443ull :
          sort_indices<2,1,6,5,4,7,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560415493ull :
          sort_indices<2,1,6,7,4,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 560416003ull :
          sort_indices<2,1,6,7,4,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593577223ull :
          sort_indices<2,3,6,1,4,5,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593577733ull :
          sort_indices<2,3,6,1,4,7,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593838343ull :
          sort_indices<2,3,6,5,4,1,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593839873ull :
          sort_indices<2,3,6,5,4,7,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593969413ull :
          sort_indices<2,3,6,7,4,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 593970433ull :
          sort_indices<2,3,6,7,4,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627131143ull :
          sort_indices<2,5,6,1,4,3,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627132163ull :
          sort_indices<2,5,6,1,4,7,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627261703ull :
          sort_indices<2,5,6,3,4,1,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627263233ull :
          sort_indices<2,5,6,3,4,7,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627523843ull :
          sort_indices<2,5,6,7,4,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 627524353ull :
          sort_indices<2,5,6,7,4,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660685573ull :
          sort_indices<2,7,6,1,4,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660686083ull :
          sort_indices<2,7,6,1,4,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660816133ull :
          sort_indices<2,7,6,3,4,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660817153ull :
          sort_indices<2,7,6,3,4,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660947203ull :
          sort_indices<2,7,6,5,4,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 660947713ull :
          sort_indices<2,7,6,5,4,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090725223ull :
          sort_indices<4,1,0,3,2,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090725733ull :
          sort_indices<4,1,0,3,2,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090855783ull :
          sort_indices<4,1,0,5,2,3,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090856803ull :
          sort_indices<4,1,0,5,2,7,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090986853ull :
          sort_indices<4,1,0,7,2,3,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090987363ull :
          sort_indices<4,1,0,7,2,5,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124148583ull :
          sort_indices<4,3,0,1,2,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124149093ull :
          sort_indices<4,3,0,1,2,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124409703ull :
          sort_indices<4,3,0,5,2,1,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124411233ull :
          sort_indices<4,3,0,5,2,7,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124540773ull :
          sort_indices<4,3,0,7,2,1,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124541793ull :
          sort_indices<4,3,0,7,2,5,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157702503ull :
          sort_indices<4,5,0,1,2,3,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157703523ull :
          sort_indices<4,5,0,1,2,7,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157833063ull :
          sort_indices<4,5,0,3,2,1,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157834593ull :
          sort_indices<4,5,0,3,2,7,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1158095203ull :
          sort_indices<4,5,0,7,2,1,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1158095713ull :
          sort_indices<4,5,0,7,2,3,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191256933ull :
          sort_indices<4,7,0,1,2,3,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191257443ull :
          sort_indices<4,7,0,1,2,5,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191387493ull :
          sort_indices<4,7,0,3,2,1,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191388513ull :
          sort_indices<4,7,0,3,2,5,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191518563ull :
          sort_indices<4,7,0,5,2,1,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191519073ull :
          sort_indices<4,7,0,5,2,3,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090741543ull :
          sort_indices<4,1,0,3,6,5,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090742053ull :
          sort_indices<4,1,0,3,6,7,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090872103ull :
          sort_indices<4,1,0,5,6,3,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1090873123ull :
          sort_indices<4,1,0,5,6,7,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1091003173ull :
          sort_indices<4,1,0,7,6,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1091003683ull :
          sort_indices<4,1,0,7,6,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124164903ull :
          sort_indices<4,3,0,1,6,5,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124165413ull :
          sort_indices<4,3,0,1,6,7,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124426023ull :
          sort_indices<4,3,0,5,6,1,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124427553ull :
          sort_indices<4,3,0,5,6,7,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124557093ull :
          sort_indices<4,3,0,7,6,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1124558113ull :
          sort_indices<4,3,0,7,6,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157718823ull :
          sort_indices<4,5,0,1,6,3,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157719843ull :
          sort_indices<4,5,0,1,6,7,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157849383ull :
          sort_indices<4,5,0,3,6,1,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1157850913ull :
          sort_indices<4,5,0,3,6,7,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1158111523ull :
          sort_indices<4,5,0,7,6,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1158112033ull :
          sort_indices<4,5,0,7,6,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191273253ull :
          sort_indices<4,7,0,1,6,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191273763ull :
          sort_indices<4,7,0,1,6,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191403813ull :
          sort_indices<4,7,0,3,6,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191404833ull :
          sort_indices<4,7,0,3,6,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191534883ull :
          sort_indices<4,7,0,5,6,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1191535393ull :
          sort_indices<4,7,0,5,6,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092814183ull :
          sort_indices<4,1,2,3,0,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092814693ull :
          sort_indices<4,1,2,3,0,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092944743ull :
          sort_indices<4,1,2,5,0,3,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092945763ull :
          sort_indices<4,1,2,5,0,7,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1093075813ull :
          sort_indices<4,1,2,7,0,3,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1093076323ull :
          sort_indices<4,1,2,7,0,5,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126237543ull :
          sort_indices<4,3,2,1,0,5,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126238053ull :
          sort_indices<4,3,2,1,0,7,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126498663ull :
          sort_indices<4,3,2,5,0,1,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126500193ull :
          sort_indices<4,3,2,5,0,7,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126629733ull :
          sort_indices<4,3,2,7,0,1,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126630753ull :
          sort_indices<4,3,2,7,0,5,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159791463ull :
          sort_indices<4,5,2,1,0,3,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159792483ull :
          sort_indices<4,5,2,1,0,7,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159922023ull :
          sort_indices<4,5,2,3,0,1,6,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159923553ull :
          sort_indices<4,5,2,3,0,7,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1160184163ull :
          sort_indices<4,5,2,7,0,1,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1160184673ull :
          sort_indices<4,5,2,7,0,3,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193345893ull :
          sort_indices<4,7,2,1,0,3,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193346403ull :
          sort_indices<4,7,2,1,0,5,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193476453ull :
          sort_indices<4,7,2,3,0,1,6,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193477473ull :
          sort_indices<4,7,2,3,0,5,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193607523ull :
          sort_indices<4,7,2,5,0,1,6,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193608033ull :
          sort_indices<4,7,2,5,0,3,6,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092838663ull :
          sort_indices<4,1,2,3,6,5,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092839173ull :
          sort_indices<4,1,2,3,6,7,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092969223ull :
          sort_indices<4,1,2,5,6,3,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1092970243ull :
          sort_indices<4,1,2,5,6,7,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1093100293ull :
          sort_indices<4,1,2,7,6,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1093100803ull :
          sort_indices<4,1,2,7,6,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126262023ull :
          sort_indices<4,3,2,1,6,5,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126262533ull :
          sort_indices<4,3,2,1,6,7,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126523143ull :
          sort_indices<4,3,2,5,6,1,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126524673ull :
          sort_indices<4,3,2,5,6,7,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126654213ull :
          sort_indices<4,3,2,7,6,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1126655233ull :
          sort_indices<4,3,2,7,6,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159815943ull :
          sort_indices<4,5,2,1,6,3,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159816963ull :
          sort_indices<4,5,2,1,6,7,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159946503ull :
          sort_indices<4,5,2,3,6,1,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1159948033ull :
          sort_indices<4,5,2,3,6,7,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1160208643ull :
          sort_indices<4,5,2,7,6,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1160209153ull :
          sort_indices<4,5,2,7,6,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193370373ull :
          sort_indices<4,7,2,1,6,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193370883ull :
          sort_indices<4,7,2,1,6,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193500933ull :
          sort_indices<4,7,2,3,6,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193501953ull :
          sort_indices<4,7,2,3,6,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193632003ull :
          sort_indices<4,7,2,5,6,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1193632513ull :
          sort_indices<4,7,2,5,6,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097008423ull :
          sort_indices<4,1,6,3,0,5,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097008933ull :
          sort_indices<4,1,6,3,0,7,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097138983ull :
          sort_indices<4,1,6,5,0,3,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097140003ull :
          sort_indices<4,1,6,5,0,7,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097270053ull :
          sort_indices<4,1,6,7,0,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097270563ull :
          sort_indices<4,1,6,7,0,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130431783ull :
          sort_indices<4,3,6,1,0,5,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130432293ull :
          sort_indices<4,3,6,1,0,7,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130692903ull :
          sort_indices<4,3,6,5,0,1,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130694433ull :
          sort_indices<4,3,6,5,0,7,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130823973ull :
          sort_indices<4,3,6,7,0,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130824993ull :
          sort_indices<4,3,6,7,0,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1163985703ull :
          sort_indices<4,5,6,1,0,3,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1163986723ull :
          sort_indices<4,5,6,1,0,7,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164116263ull :
          sort_indices<4,5,6,3,0,1,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164117793ull :
          sort_indices<4,5,6,3,0,7,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164378403ull :
          sort_indices<4,5,6,7,0,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164378913ull :
          sort_indices<4,5,6,7,0,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197540133ull :
          sort_indices<4,7,6,1,0,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197540643ull :
          sort_indices<4,7,6,1,0,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197670693ull :
          sort_indices<4,7,6,3,0,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197671713ull :
          sort_indices<4,7,6,3,0,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197801763ull :
          sort_indices<4,7,6,5,0,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197802273ull :
          sort_indices<4,7,6,5,0,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097016583ull :
          sort_indices<4,1,6,3,2,5,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097017093ull :
          sort_indices<4,1,6,3,2,7,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097147143ull :
          sort_indices<4,1,6,5,2,3,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097148163ull :
          sort_indices<4,1,6,5,2,7,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097278213ull :
          sort_indices<4,1,6,7,2,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1097278723ull :
          sort_indices<4,1,6,7,2,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130439943ull :
          sort_indices<4,3,6,1,2,5,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130440453ull :
          sort_indices<4,3,6,1,2,7,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130701063ull :
          sort_indices<4,3,6,5,2,1,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130702593ull :
          sort_indices<4,3,6,5,2,7,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130832133ull :
          sort_indices<4,3,6,7,2,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1130833153ull :
          sort_indices<4,3,6,7,2,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1163993863ull :
          sort_indices<4,5,6,1,2,3,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1163994883ull :
          sort_indices<4,5,6,1,2,7,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164124423ull :
          sort_indices<4,5,6,3,2,1,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164125953ull :
          sort_indices<4,5,6,3,2,7,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164386563ull :
          sort_indices<4,5,6,7,2,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1164387073ull :
          sort_indices<4,5,6,7,2,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197548293ull :
          sort_indices<4,7,6,1,2,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197548803ull :
          sort_indices<4,7,6,1,2,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197678853ull :
          sort_indices<4,7,6,3,2,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197679873ull :
          sort_indices<4,7,6,3,2,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197809923ull :
          sort_indices<4,7,6,5,2,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1197810433ull :
          sort_indices<4,7,6,5,2,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627596103ull :
          sort_indices<6,1,0,3,2,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627596613ull :
          sort_indices<6,1,0,3,2,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627726663ull :
          sort_indices<6,1,0,5,2,3,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627727683ull :
          sort_indices<6,1,0,5,2,7,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627857733ull :
          sort_indices<6,1,0,7,2,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627858243ull :
          sort_indices<6,1,0,7,2,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661019463ull :
          sort_indices<6,3,0,1,2,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661019973ull :
          sort_indices<6,3,0,1,2,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661280583ull :
          sort_indices<6,3,0,5,2,1,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661282113ull :
          sort_indices<6,3,0,5,2,7,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661411653ull :
          sort_indices<6,3,0,7,2,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661412673ull :
          sort_indices<6,3,0,7,2,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694573383ull :
          sort_indices<6,5,0,1,2,3,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694574403ull :
          sort_indices<6,5,0,1,2,7,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694703943ull :
          sort_indices<6,5,0,3,2,1,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694705473ull :
          sort_indices<6,5,0,3,2,7,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694966083ull :
          sort_indices<6,5,0,7,2,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694966593ull :
          sort_indices<6,5,0,7,2,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728127813ull :
          sort_indices<6,7,0,1,2,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728128323ull :
          sort_indices<6,7,0,1,2,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728258373ull :
          sort_indices<6,7,0,3,2,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728259393ull :
          sort_indices<6,7,0,3,2,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728389443ull :
          sort_indices<6,7,0,5,2,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728389953ull :
          sort_indices<6,7,0,5,2,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627604263ull :
          sort_indices<6,1,0,3,4,5,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627604773ull :
          sort_indices<6,1,0,3,4,7,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627734823ull :
          sort_indices<6,1,0,5,4,3,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627735843ull :
          sort_indices<6,1,0,5,4,7,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627865893ull :
          sort_indices<6,1,0,7,4,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1627866403ull :
          sort_indices<6,1,0,7,4,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661027623ull :
          sort_indices<6,3,0,1,4,5,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661028133ull :
          sort_indices<6,3,0,1,4,7,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661288743ull :
          sort_indices<6,3,0,5,4,1,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661290273ull :
          sort_indices<6,3,0,5,4,7,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661419813ull :
          sort_indices<6,3,0,7,4,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1661420833ull :
          sort_indices<6,3,0,7,4,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694581543ull :
          sort_indices<6,5,0,1,4,3,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694582563ull :
          sort_indices<6,5,0,1,4,7,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694712103ull :
          sort_indices<6,5,0,3,4,1,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694713633ull :
          sort_indices<6,5,0,3,4,7,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694974243ull :
          sort_indices<6,5,0,7,4,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1694974753ull :
          sort_indices<6,5,0,7,4,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728135973ull :
          sort_indices<6,7,0,1,4,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728136483ull :
          sort_indices<6,7,0,1,4,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728266533ull :
          sort_indices<6,7,0,3,4,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728267553ull :
          sort_indices<6,7,0,3,4,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728397603ull :
          sort_indices<6,7,0,5,4,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1728398113ull :
          sort_indices<6,7,0,5,4,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629685063ull :
          sort_indices<6,1,2,3,0,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629685573ull :
          sort_indices<6,1,2,3,0,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629815623ull :
          sort_indices<6,1,2,5,0,3,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629816643ull :
          sort_indices<6,1,2,5,0,7,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629946693ull :
          sort_indices<6,1,2,7,0,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629947203ull :
          sort_indices<6,1,2,7,0,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663108423ull :
          sort_indices<6,3,2,1,0,5,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663108933ull :
          sort_indices<6,3,2,1,0,7,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663369543ull :
          sort_indices<6,3,2,5,0,1,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663371073ull :
          sort_indices<6,3,2,5,0,7,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663500613ull :
          sort_indices<6,3,2,7,0,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663501633ull :
          sort_indices<6,3,2,7,0,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696662343ull :
          sort_indices<6,5,2,1,0,3,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696663363ull :
          sort_indices<6,5,2,1,0,7,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696792903ull :
          sort_indices<6,5,2,3,0,1,4,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696794433ull :
          sort_indices<6,5,2,3,0,7,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1697055043ull :
          sort_indices<6,5,2,7,0,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1697055553ull :
          sort_indices<6,5,2,7,0,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730216773ull :
          sort_indices<6,7,2,1,0,3,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730217283ull :
          sort_indices<6,7,2,1,0,5,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730347333ull :
          sort_indices<6,7,2,3,0,1,4,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730348353ull :
          sort_indices<6,7,2,3,0,5,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730478403ull :
          sort_indices<6,7,2,5,0,1,4,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730478913ull :
          sort_indices<6,7,2,5,0,3,4,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629701383ull :
          sort_indices<6,1,2,3,4,5,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629701893ull :
          sort_indices<6,1,2,3,4,7,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629831943ull :
          sort_indices<6,1,2,5,4,3,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629832963ull :
          sort_indices<6,1,2,5,4,7,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629963013ull :
          sort_indices<6,1,2,7,4,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1629963523ull :
          sort_indices<6,1,2,7,4,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663124743ull :
          sort_indices<6,3,2,1,4,5,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663125253ull :
          sort_indices<6,3,2,1,4,7,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663385863ull :
          sort_indices<6,3,2,5,4,1,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663387393ull :
          sort_indices<6,3,2,5,4,7,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663516933ull :
          sort_indices<6,3,2,7,4,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1663517953ull :
          sort_indices<6,3,2,7,4,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696678663ull :
          sort_indices<6,5,2,1,4,3,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696679683ull :
          sort_indices<6,5,2,1,4,7,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696809223ull :
          sort_indices<6,5,2,3,4,1,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1696810753ull :
          sort_indices<6,5,2,3,4,7,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1697071363ull :
          sort_indices<6,5,2,7,4,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1697071873ull :
          sort_indices<6,5,2,7,4,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730233093ull :
          sort_indices<6,7,2,1,4,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730233603ull :
          sort_indices<6,7,2,1,4,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730363653ull :
          sort_indices<6,7,2,3,4,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730364673ull :
          sort_indices<6,7,2,3,4,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730494723ull :
          sort_indices<6,7,2,5,4,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1730495233ull :
          sort_indices<6,7,2,5,4,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631782183ull :
          sort_indices<6,1,4,3,0,5,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631782693ull :
          sort_indices<6,1,4,3,0,7,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631912743ull :
          sort_indices<6,1,4,5,0,3,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631913763ull :
          sort_indices<6,1,4,5,0,7,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1632043813ull :
          sort_indices<6,1,4,7,0,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1632044323ull :
          sort_indices<6,1,4,7,0,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665205543ull :
          sort_indices<6,3,4,1,0,5,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665206053ull :
          sort_indices<6,3,4,1,0,7,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665466663ull :
          sort_indices<6,3,4,5,0,1,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665468193ull :
          sort_indices<6,3,4,5,0,7,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665597733ull :
          sort_indices<6,3,4,7,0,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665598753ull :
          sort_indices<6,3,4,7,0,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698759463ull :
          sort_indices<6,5,4,1,0,3,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698760483ull :
          sort_indices<6,5,4,1,0,7,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698890023ull :
          sort_indices<6,5,4,3,0,1,2,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698891553ull :
          sort_indices<6,5,4,3,0,7,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1699152163ull :
          sort_indices<6,5,4,7,0,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1699152673ull :
          sort_indices<6,5,4,7,0,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732313893ull :
          sort_indices<6,7,4,1,0,3,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732314403ull :
          sort_indices<6,7,4,1,0,5,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732444453ull :
          sort_indices<6,7,4,3,0,1,2,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732445473ull :
          sort_indices<6,7,4,3,0,5,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732575523ull :
          sort_indices<6,7,4,5,0,1,2,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732576033ull :
          sort_indices<6,7,4,5,0,3,2,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631790343ull :
          sort_indices<6,1,4,3,2,5,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631790853ull :
          sort_indices<6,1,4,3,2,7,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631920903ull :
          sort_indices<6,1,4,5,2,3,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1631921923ull :
          sort_indices<6,1,4,5,2,7,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1632051973ull :
          sort_indices<6,1,4,7,2,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1632052483ull :
          sort_indices<6,1,4,7,2,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665213703ull :
          sort_indices<6,3,4,1,2,5,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665214213ull :
          sort_indices<6,3,4,1,2,7,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665474823ull :
          sort_indices<6,3,4,5,2,1,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665476353ull :
          sort_indices<6,3,4,5,2,7,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665605893ull :
          sort_indices<6,3,4,7,2,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1665606913ull :
          sort_indices<6,3,4,7,2,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698767623ull :
          sort_indices<6,5,4,1,2,3,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698768643ull :
          sort_indices<6,5,4,1,2,7,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698898183ull :
          sort_indices<6,5,4,3,2,1,0,7,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1698899713ull :
          sort_indices<6,5,4,3,2,7,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1699160323ull :
          sort_indices<6,5,4,7,2,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1699160833ull :
          sort_indices<6,5,4,7,2,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732322053ull :
          sort_indices<6,7,4,1,2,3,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732322563ull :
          sort_indices<6,7,4,1,2,5,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732452613ull :
          sort_indices<6,7,4,3,2,1,0,5,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732453633ull :
          sort_indices<6,7,4,3,2,5,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732583683ull :
          sort_indices<6,7,4,5,2,1,0,3,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        case 1732584193ull :
          sort_indices<6,7,4,5,2,3,0,1,0,1,-1,1>(in, out, d[0], d[1], d[2], d[3], d[4], d[5], d[6], d[7]); break;
        default:
          throw std::logic_error("This case has not been implemented in prim_op_var.h");
    }
  }
  else {
    throw std::logic_error("This case has not been implemented in prim_op_var.h");
  }
}

}
#endif

