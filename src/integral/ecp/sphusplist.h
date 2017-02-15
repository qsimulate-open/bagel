//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: sphusplist.h
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


#ifndef __SRC_INTEGRAL_ECP_SPHUSPLIST_H
#define __SRC_INTEGRAL_ECP_SPHUSPLIST_H

#include <src/util/constants.h>
#include <vector>
#include <functional>

namespace bagel {

// USP = Unitary Sphere Polynomial = x^i y^j z^k
struct SphUSPList {
  private:

    std::function<std::vector<double>(const int)> sphuspfunc[11];

    static std::vector<double> sphusp_0(const int);
    static std::vector<double> sphusp_1(const int);
    static std::vector<double> sphusp_2(const int);
    static std::vector<double> sphusp_3(const int);
    static std::vector<double> sphusp_4(const int);
    static std::vector<double> sphusp_5(const int);
    static std::vector<double> sphusp_6(const int);
    static std::vector<double> sphusp_7(const int);
    static std::vector<double> sphusp_8(const int);
    static std::vector<double> sphusp_9(const int);
    static std::vector<double> sphusp_10(const int);

  public:
    SphUSPList() {
      sphuspfunc[0] = &sphusp_0;
      sphuspfunc[1] = &sphusp_1;
      sphuspfunc[2] = &sphusp_2;
      sphuspfunc[3] = &sphusp_3;
      sphuspfunc[4] = &sphusp_4;
      sphuspfunc[5] = &sphusp_5;
      sphuspfunc[6] = &sphusp_6;
      sphuspfunc[7] = &sphusp_7;
      sphuspfunc[8] = &sphusp_8;
      sphuspfunc[9] = &sphusp_9;
      sphuspfunc[10] = &sphusp_10;
    }

    std::vector<double> sphuspfunc_call(const int l, const int m) const { assert(l <= 10); return sphuspfunc[l](m); }

};

}

#endif
