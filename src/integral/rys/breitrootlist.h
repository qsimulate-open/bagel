//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: breitrootlist.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __SRC_RYSINT_BREITROOTLIST_H
#define __SRC_RYSINT_BREITROOTLIST_H

#include <functional>
#include <src/util/constants.h>

namespace bagel {

struct BreitRootList  {
  private:
    std::function<void (const double*, double*, double*, const int)> rfunc[RYS_MAX + 1];

    static void breitroot1(const double*, double*, double*, const int);
    static void breitroot2(const double*, double*, double*, const int);
    static void breitroot3(const double*, double*, double*, const int);
    static void breitroot4(const double*, double*, double*, const int);
    static void breitroot5(const double*, double*, double*, const int);
    static void breitroot6(const double*, double*, double*, const int);
    static void breitroot7(const double*, double*, double*, const int);
    static void breitroot8(const double*, double*, double*, const int);
    static void breitroot9(const double*, double*, double*, const int);
    static void breitroot10(const double*, double*, double*, const int);
    static void breitroot11(const double*, double*, double*, const int);
    static void breitroot12(const double*, double*, double*, const int);
    static void breitroot13(const double*, double*, double*, const int);

  public:
    BreitRootList() {
      rfunc[1] = &breitroot1;
      rfunc[2] = &breitroot2;
      rfunc[3] = &breitroot3;
      rfunc[4] = &breitroot4;
      rfunc[5] = &breitroot5;
      rfunc[6] = &breitroot6;
      rfunc[7] = &breitroot7;
      rfunc[8] = &breitroot8;
      rfunc[9] = &breitroot9;
      rfunc[10] = &breitroot10;
      rfunc[11] = &breitroot11;
      rfunc[12] = &breitroot12;
      rfunc[13] = &breitroot13;
    }

    void root(const int i, const double* a1, double* a2, double* a3, const int a4) const { rfunc[i](a1, a2, a3, a4); }

};

const static BreitRootList breitroot__;

}

#endif

