//
// BAGEL - Parallel electron correlation program.
// Filename: erirootlist.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#ifndef __SRC_RYSINT_ERIROOTLIST_H
#define __SRC_RYSINT_ERIROOTLIST_H

#include <functional>
#include <src/rysint/intparam.h>

namespace bagel {

struct ERIRootList  {
  private:
    std::function<void (const double*, double*, double*, const int)> rfunc[RYS_MAX + 1];

    static void eriroot1(const double*, double*, double*, const int);
    static void eriroot2(const double*, double*, double*, const int);
    static void eriroot3(const double*, double*, double*, const int);
    static void eriroot4(const double*, double*, double*, const int);
    static void eriroot5(const double*, double*, double*, const int);
    static void eriroot6(const double*, double*, double*, const int);
    static void eriroot7(const double*, double*, double*, const int);
    static void eriroot8(const double*, double*, double*, const int);
    static void eriroot9(const double*, double*, double*, const int);
    static void eriroot10(const double*, double*, double*, const int);
    static void eriroot11(const double*, double*, double*, const int);
    static void eriroot12(const double*, double*, double*, const int);
    static void eriroot13(const double*, double*, double*, const int);


  public:
    ERIRootList() {
      rfunc[1] = &eriroot1;
      rfunc[2] = &eriroot2;
      rfunc[3] = &eriroot3;
      rfunc[4] = &eriroot4;
      rfunc[5] = &eriroot5;
      rfunc[6] = &eriroot6;
      rfunc[7] = &eriroot7;
      rfunc[8] = &eriroot8;
      rfunc[9] = &eriroot9;
      rfunc[10] = &eriroot10;
      rfunc[11] = &eriroot11;
      rfunc[12] = &eriroot12;
      rfunc[13] = &eriroot13;
    }

    void root(const int i, const double* a1, double* a2, double* a3, const int a4) const { rfunc[i](a1, a2, a3, a4); }

};

}

#endif

