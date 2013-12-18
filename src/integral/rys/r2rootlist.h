//
// BAGEL - Parallel electron correlation program.
// Filename: r2rootlist.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __SRC_INTEGRAL_RYS_R2ROOTLIST_H
#define __SRC_INTEGRAL_RYS_R2ROOTLIST_H

#include <functional>
#include <src/util/constants.h>

namespace bagel {

struct R2RootList  {
  private:
    std::function<void (const double*, double*, double*, const int)> rfunc[RYS_MAX + 1];

    static void r2root1(const double*, double*, double*, const int);
    static void r2root2(const double*, double*, double*, const int);
    static void r2root3(const double*, double*, double*, const int);
    static void r2root4(const double*, double*, double*, const int);
    static void r2root5(const double*, double*, double*, const int);
    static void r2root6(const double*, double*, double*, const int);
    static void r2root7(const double*, double*, double*, const int);
    static void r2root8(const double*, double*, double*, const int);
    static void r2root9(const double*, double*, double*, const int);
    static void r2root10(const double*, double*, double*, const int);
    static void r2root11(const double*, double*, double*, const int);
    static void r2root12(const double*, double*, double*, const int);
    static void r2root13(const double*, double*, double*, const int);


  public:
    R2RootList() {
      rfunc[1] = &r2root1;
      rfunc[2] = &r2root2;
      rfunc[3] = &r2root3;
      rfunc[4] = &r2root4;
      rfunc[5] = &r2root5;
      rfunc[6] = &r2root6;
      rfunc[7] = &r2root7;
      rfunc[8] = &r2root8;
      rfunc[9] = &r2root9;
      rfunc[10] = &r2root10;
      rfunc[11] = &r2root11;
      rfunc[12] = &r2root12;
      rfunc[13] = &r2root13;
    }

    void root(const int i, const double* a1, double* a2, double* a3, const int a4) const { rfunc[i](a1, a2, a3, a4); }

};

const static R2RootList r2root__;

}

#endif

