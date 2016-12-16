//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: spin2rootlist.h
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


#ifndef __SRC_RYSINT_SPIN2ROOTLIST_H
#define __SRC_RYSINT_SPIN2ROOTLIST_H

#include <functional>
#include <src/util/constants.h>

namespace bagel {

struct Spin2RootList  {
  private:
    std::function<void (const double*, double*, double*, const int)> rfunc[RYS_MAX + 1];

    static void spin2root1(const double*, double*, double*, const int);
    static void spin2root2(const double*, double*, double*, const int);
    static void spin2root3(const double*, double*, double*, const int);
    static void spin2root4(const double*, double*, double*, const int);
    static void spin2root5(const double*, double*, double*, const int);
    static void spin2root6(const double*, double*, double*, const int);
    static void spin2root7(const double*, double*, double*, const int);
    static void spin2root8(const double*, double*, double*, const int);
    static void spin2root9(const double*, double*, double*, const int);
    static void spin2root10(const double*, double*, double*, const int);
    static void spin2root11(const double*, double*, double*, const int);
    static void spin2root12(const double*, double*, double*, const int);
    static void spin2root13(const double*, double*, double*, const int);

  public:
    Spin2RootList() {
      rfunc[1] = &spin2root1;
      rfunc[2] = &spin2root2;
      rfunc[3] = &spin2root3;
      rfunc[4] = &spin2root4;
      rfunc[5] = &spin2root5;
      rfunc[6] = &spin2root6;
      rfunc[7] = &spin2root7;
      rfunc[8] = &spin2root8;
      rfunc[9] = &spin2root9;
      rfunc[10] = &spin2root10;
      rfunc[11] = &spin2root11;
      rfunc[12] = &spin2root12;
      rfunc[13] = &spin2root13;
    }

    void root(const int i, const double* a1, double* a2, double* a3, const int a4) const { rfunc[i](a1, a2, a3, a4); }

};

const static Spin2RootList spin2root__;

}

#endif

