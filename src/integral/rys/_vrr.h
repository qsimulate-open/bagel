//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _vrr.h
// Copyright (C) 2012 Toru Shiozaki
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


// replaces generated codes _vrr_xxxx.cc etc

#ifndef __SRC_RYSINT____VRR_H
#define __SRC_RYSINT____VRR_H

#include <algorithm>
#include <stdexcept>

namespace bagel {

template<int a_, int c_, int rank_, typename DataType>
void vrr(DataType* data_, const DataType* C00, const DataType* D00, const DataType* B00, const DataType* B01, const DataType* B10) {

  static_assert(a_>=0 && c_>=0 && rank_ >= 1, "parameter(s) wrong in vrr");

  alignas(32) DataType C00_[rank_];
  alignas(32) DataType D00_[rank_];
  alignas(32) DataType B00_[rank_];
  alignas(32) DataType B01_[rank_];
  alignas(32) DataType B10_[rank_];
  std::copy_n(C00, rank_, C00_);
  std::copy_n(D00, rank_, D00_);
  std::copy_n(B00, rank_, B00_);
  std::copy_n(B01, rank_, B01_);
  std::copy_n(B10, rank_, B10_);

  if ((a_ > 1 && c_ > 1) || (a_ == 1 && c_ >  1)) {
    // c == 0
    for (int t = 0; t != rank_; ++t)
      data_[rank_*0+t] = 1.0;

    for (int t = 0; t != rank_; ++t)
      data_[rank_*1+t] = C00_[t];

    if (a_ == 2) {
      for (int t = 0; t != rank_; ++t)
        data_[rank_*2+t] = C00_[t] * data_[rank_*1+t] + B10_[t];
    } else if (a_ > 1) {
      alignas(32) DataType B10_current[rank_];
      // a == 2
      for (int t = 0; t != rank_; ++t)
        B10_current[t] = B10_[t];
      for (int t = 0; t != rank_; ++t)
        data_[rank_*2+t] = C00_[t] * data_[rank_*1+t] + B10_current[t];
      // a > 2
      for (int a = 3; a != a_+1; ++a) {
        for (int t = 0; t != rank_; ++t)
          B10_current[t] += B10_[t];
        for (int t = 0; t != rank_; ++t)
          data_[rank_*a+t] = C00_[t] * data_[rank_*(a-1)+t] + B10_current[t] * data_[rank_*(a-2)+t];
      }
    }
    // c == 1
    for (int t = 0; t != rank_; ++t)
      data_[rank_*(a_+1)+t] = D00_[t];

    alignas(32) DataType cB00_current[rank_];
    for (int t = 0; t != rank_; ++t)
      cB00_current[t] = B00_[t];

    for (int t = 0; t != rank_; ++t)
      data_[rank_*(a_+2)+t] = C00_[t] * data_[rank_*(a_+1)+t] + cB00_current[t];

    if (a_ == 2) {
      for (int t = 0; t != rank_; ++t)
        data_[rank_*(a_+3)+t] = C00_[t] * data_[rank_*(a_+2)+t] + B10_[t] * data_[rank_*(a_+1)+t] + cB00_current[t] * data_[rank_+t];
    } else if (a_ > 1) {

      alignas(32) DataType B10_current[rank_];
      for (int t = 0; t != rank_; ++t)
        B10_current[t] = B10_[t];

      for (int t = 0; t != rank_; ++t)
        data_[rank_*(a_+3)+t] = C00_[t] * data_[rank_*(a_+2)+t] + B10_current[t] * data_[rank_*(a_+1)+t] + cB00_current[t] * data_[rank_+t];

      for (int a = 3; a != a_+1; ++a) {
        for (int t = 0; t != rank_; ++t)
          B10_current[t] += B10_[t];

        for (int t = 0; t != rank_; ++t)
         data_[rank_*(a_+1+a)+t] = C00_[t] * data_[rank_*(a_+a)+t] + B10_current[t] * data_[rank_*(a_-1+a)+t] + cB00_current[t] * data_[rank_*(a-1)+t];
      }
    }
    // c > 1
    alignas(32) DataType B01_current[rank_] = {0.0};
    for (int c = 2; c != c_+1; ++c) {
      for (int t = 0; t != rank_; ++t)
        B01_current[t] += B01_[t];

      for (int t = 0; t != rank_; ++t)
        data_[rank_*(a_+1)*c+t] = D00_[t] * data_[rank_*(a_+1)*(c-1)+t] + B01_current[t] * data_[rank_*(a_+1)*(c-2)+t];

      for (int t = 0; t != rank_; ++t)
        cB00_current[t] += B00_[t];

      for (int t = 0; t != rank_; ++t)
        data_[rank_*((a_+1)*c+1)+t] = C00_[t] * data_[rank_*((a_+1)*c)+t] + cB00_current[t] * data_[rank_*((a_+1)*(c-1))+t];

      if (a_ > 1) {
        alignas(32) DataType B10_current[rank_];
        for (int t = 0; t != rank_; ++t)
          B10_current[t] = B10_[t];

        for (int t = 0; t != rank_; ++t)
          data_[rank_*((a_+1)*c+2)+t] = C00_[t] * data_[rank_*((a_+1)*c+1)+t] + B10_current[t] * data_[rank_*(a_+1)*c+t] + cB00_current[t] * data_[rank_*((a_+1)*(c-1)+1)+t];

        for (int a = 3; a != a_+1; ++a) {
          for (int t = 0; t != rank_; ++t)
            B10_current[t] += B10_[t];

          for (int t = 0; t != rank_; ++t)
            data_[rank_*((a_+1)*c+a)+t] = C00_[t] * data_[rank_*((a_+1)*c+a-1)+t] + B10_current[t] * data_[rank_*((a_+1)*c+a-2)+t] + cB00_current[t] * data_[rank_*((a_+1)*(c-1)+a-1)+t];
        }
      }
    }

  } else if (a_ == 0 && c_ >  0) {
    for (int t = 0; t != rank_; ++t)
      data_[rank_*0+t] = 1.0;

    for (int t = 0; t != rank_; ++t)
      data_[rank_*1+t] = D00_[t];

    if (c_ == 2) {
      for (int t = 0; t != rank_; ++t)
        data_[rank_*2+t] = D00_[t] * data_[rank_*1+t] + B01_[t];
    } else if (c_ > 2) {
      alignas(32) DataType B01_current[rank_];
      // c == 2
      for (int t = 0; t != rank_; ++t)
        B01_current[t] = B01_[t];
      for (int t = 0; t != rank_; ++t)
        data_[rank_*2+t] = D00_[t] * data_[rank_+t] + B01_current[t];
      // c > 2
      for (int c = 3; c != c_+1; ++c) {
        for (int t = 0; t != rank_; ++t)
          B01_current[t] += B01_[t];

        for (int t = 0; t != rank_; ++t)
          data_[rank_*c+t] = D00_[t] * data_[rank_*(c-1)+t] + B01_current[t] * data_[rank_*(c-2)+t];
      }
    }
  } else if (a_ >  1 && c_ == 1) {
    for (int t = 0; t != rank_; ++t)
      data_[rank_*0+t] = 1.0;

    for (int t = 0; t != rank_; ++t)
      data_[rank_*1+t] = C00_[t];

    if (a_ == 2) {
      for (int t = 0; t != rank_; ++t)
        data_[rank_*2+t] = C00_[t] * data_[rank_*1+t] + B10_[t];

      for (int t = 0; t != rank_; ++t)
        data_[rank_*(a_+1)+t] = D00_[t];

      for (int t = 0; t != rank_; ++t)
        data_[rank_*(a_+2)+t] = C00_[t] * data_[rank_*(a_+1)+t] + B00_[t];

      for (int t = 0; t != rank_; ++t)
        data_[rank_*(a_+3)+t] = C00_[t] * data_[rank_*(a_+2)+t] + B10_[t] * data_[rank_*(a_+1)+t] + B00_[t] * data_[rank_+t];
    } else {
      alignas(32) DataType B10_current[rank_];
      // a == 2
      for (int t = 0; t != rank_; ++t)
        B10_current[t] = B10_[t];
      for (int t = 0; t != rank_; ++t)
        data_[rank_*2+t] = C00_[t] * data_[rank_*1+t] + B10_current[t];
      // a > 2
      for (int a = 3; a != a_+1; ++a) {
        for (int t = 0; t != rank_; ++t)
          B10_current[t] += B10_[t];
        for (int t = 0; t != rank_; ++t)
          data_[rank_*a+t] = C00_[t] * data_[rank_*(a-1)+t] + B10_current[t] * data_[rank_*(a-2)+t];
      }

      for (int t = 0; t != rank_; ++t)
        data_[rank_*(a_+1)+t] = D00_[t];

      for (int t = 0; t != rank_; ++t)
        data_[rank_*(a_+2)+t] = C00_[t] * data_[rank_*(a_+1)+t] + B00_[t];

      // a == 2
      for (int t = 0; t != rank_; ++t)
        B10_current[t] = B10_[t];
      for (int t = 0; t != rank_; ++t)
        data_[rank_*(a_+3)+t] = C00_[t] * data_[rank_*(a_+2)+t] + B10_current[t] * data_[rank_*(a_+1)+t] + B00_[t] * data_[rank_+t];

      for (int a = 3; a != a_+1; ++a) {
        for (int t = 0; t != rank_; ++t)
          B10_current[t] += B10_[t];
        for (int t = 0; t != rank_; ++t)
          data_[rank_*(a_+1+a)+t] = C00_[t] * data_[rank_*(a_+a)+t] + B10_current[t] * data_[rank_*(a_+a-1)+t] + B00_[t] * data_[rank_*(a-1)+t];
      }
    }

  } else if (a_ == 1 && c_ == 1) {
    for (int t = 0; t != rank_; ++t)
      data_[rank_*0+t] = 1.0;

    for (int t = 0; t != rank_; ++t)
      data_[rank_*1+t] = C00_[t];

    for (int t = 0; t != rank_; ++t)
      data_[rank_*2+t] = D00_[t];

    for (int t = 0; t != rank_; ++t)
      data_[rank_*3+t] = C00_[t] * data_[rank_*2+t] + B00_[t];

  } else if (a_ >  0 && c_ == 0) {
    for (int t = 0; t != rank_; ++t)
      data_[rank_*0+t] = 1.0;

    for (int t = 0; t != rank_; ++t)
      data_[rank_*1+t] = C00_[t];

    if (a_ == 2) {
      for (int t = 0; t != rank_; ++t)
        data_[rank_*2+t] = C00_[t] * data_[rank_*1+t] + B10_[t];
    } else if (a_ > 2) {
      alignas(32) DataType B10_current[rank_];
      // a == 2
      for (int t = 0; t != rank_; ++t)
        B10_current[t] = B10_[t];
      for (int t = 0; t != rank_; ++t)
        data_[rank_*2+t] = C00_[t] * data_[rank_*1+t] + B10_current[t];
      // a > 2
      for (int a = 3; a != a_+1; ++a) {
        for (int t = 0; t != rank_; ++t)
          B10_current[t] += B10_[t];
        for (int t = 0; t != rank_; ++t)
          data_[rank_*a+t] = C00_[t] * data_[rank_*(a-1)+t] + B10_current[t] * data_[rank_*(a-2)+t];
      }
    }
  } else if (a_ == 0 && c_ == 0) {
    for (int t = 0; t != rank_; ++t)
      data_[t] = 1.0;
  } else {
    throw std::runtime_error("strange");
  }
}


}
#endif
