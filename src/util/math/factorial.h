//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: factorial.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

//

#ifndef __BAGEL_UTIL_MATH_FACTORIAL_H
#define __BAGEL_UTIL_MATH_FACTORIAL_H

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib> // for std::abs

namespace bagel {

class Factorial {
    constexpr static int max_ = 21;
    std::array<size_t, max_> f_;
  public:
    Factorial() {
      std::fill(f_.begin(), f_.end(), 0);
      f_[ 0] =                    1ull; f_[ 1] =                    1ull; f_[ 2] =                    2ull; f_[ 3] =                    6ull;
      f_[ 4] =                   24ull; f_[ 5] =                  120ull; f_[ 6] =                  720ull; f_[ 7] =                 5040ull;
      f_[ 8] =                40320ull; f_[ 9] =               362880ull; f_[10] =              3628800ull; f_[11] =             39916800ull;
      f_[12] =            479001600ull; f_[13] =           6227020800ull; f_[14] =          87178291200ull; f_[15] =        1307674368000ull;
      f_[16] =       20922789888000ull; f_[17] =      355687428096000ull; f_[18] =     6402373705728000ull; f_[19] =   121645100408832000ull;
      f_[20] =  2432902008176640000ull;
    }
    size_t operator()(const int i) const { assert(i >= 0 && i < max_); return f_[i]; }
    int max() const { return max_; }
};

class DoubleFactorial {
    constexpr static int max_ = 18;
    std::array<size_t, max_> df_;
  public:
    DoubleFactorial() {
      std::fill(df_.begin(), df_.end(), 0);
      df_[ 0] =                    1ull; df_[ 1] =                    1ull; df_[ 2] =                    3ull; df_[ 3] =                   15ull;
      df_[ 4] =                  105ull; df_[ 5] =                  945ull; df_[ 6] =                10395ull; df_[ 7] =               135135ull;
      df_[ 8] =              2027025ull; df_[ 9] =             34459425ull; df_[10] =            654729075ull; df_[11] =          13749310575ull;
      df_[12] =         316234143225ull; df_[13] =        7905853580625ull; df_[14] =      213458046676875ull; df_[15] =     6190283353629375ull;
      df_[16] =   191898783962510625ull; df_[17] =  6332659870762850625ull;
    }
    size_t operator()(const int i) const {
      assert(i >= -1 && i < 2*max_-1 && std::abs(i) % 2 == 1);
      const int j = (i+1)/2;
      return df_[j];
    }
};

}

#endif

#if 0
// Code to generate this table
#include <iostream>
#include <iomanip>
#include <cassert>
#include "mpreal.h"
using namespace mpfr;
using namespace std;

mpreal fac(const mpreal i) {
  if (i < 1.1) { return 1; }
  else { return i * fac(i-1); }
}

mpreal dfac(const mpreal i) {
  if (i < 1.1) { return 1; }
  else { return i * dfac(i-2); }
}

int main() {
  mpreal::set_default_prec(10000);
  for (int i = 0; i != 21; ++i) {
    mpreal out = fac(i);
    // out should be an integer
    if (abs(out-mpreal(out.toULong())) > 1.0e-10) assert(false);
    cout << "f_[" << setw(2) << i << "] = " << setw(20) << out.toULong() << "ull;";
    if (i % 4 == 3) cout << endl;
    else cout << " ";
  }
  cout << endl;

  for (int i = 1; i != 18; ++i) {
    mpreal out = dfac(2*i-1);
    // out should be an integer
    if (abs(out-mpreal(out.toULong())) > 1.0e-10) assert(false);
    cout << "df_[" << setw(2) << i << "] = " << setw(20) << out.toULong() << "ull;";
    if (i % 4 == 3) cout << endl;
    else cout << " ";
  }
  cout << endl;

  return 0;
}
#endif
