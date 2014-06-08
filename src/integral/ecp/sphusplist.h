//
// BAGEL - Parallel electron correlation program.
// Filename: sphusplist.h
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


#ifndef __SRC_INTEGRAL_ECP_SPHUSPLIST_H
#define __SRC_INTEGRAL_ECP_SPHUSPLIST_H

#include <src/util/constants.h>
#include <functional>

namespace bagel {

struct SphUSPList {
  private:

    std::function<void (const int, std::vector<double>&)> sphuspfunc[ANG_HRR_END+1];

    static void sphusp_0(const int, std::vector<double>&);
    static void sphusp_1(const int, std::vector<double>&);
    static void sphusp_2(const int, std::vector<double>&);
    static void sphusp_3(const int, std::vector<double>&);
    static void sphusp_4(const int, std::vector<double>&);
    static void sphusp_5(const int, std::vector<double>&);
    static void sphusp_6(const int, std::vector<double>&);
    static void sphusp_7(const int, std::vector<double>&);

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
    }

    void sphuspfunc_call(const int l, const int m, std::vector<double>& c) const { sphuspfunc[l](m, c); }

};

}

#endif
