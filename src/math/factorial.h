//
// BAGEL - Parallel electron correlation program.
// Filename: factorial.h
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

//

#ifndef __BAGEL_UTIL_FACTORIAL_H
#define __BAGEL_UTIL_FACTORIAL_H

#include <algorithm>
#include <array>
#include <cassert>

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
  return 0;
}
#endif
