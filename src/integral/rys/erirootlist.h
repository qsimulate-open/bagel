//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: erirootlist.h
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


#ifndef __SRC_INTEGRAL_RYS_ERIROOTLIST_H
#define __SRC_INTEGRAL_RYS_ERIROOTLIST_H

#include <functional>
#include <src/util/constants.h>

namespace bagel {

struct ERIRootList  {
  private:
    std::function<void (const double*, double*, double*, const int)> rfunc[51];

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
    static void eriroot14(const double*, double*, double*, const int);
    static void eriroot15(const double*, double*, double*, const int);
    static void eriroot16(const double*, double*, double*, const int);
    static void eriroot17(const double*, double*, double*, const int);
    static void eriroot18(const double*, double*, double*, const int);
    static void eriroot19(const double*, double*, double*, const int);
    static void eriroot20(const double*, double*, double*, const int);
    static void eriroot21(const double*, double*, double*, const int);
    static void eriroot22(const double*, double*, double*, const int);
    static void eriroot23(const double*, double*, double*, const int);
    static void eriroot24(const double*, double*, double*, const int);
    static void eriroot25(const double*, double*, double*, const int);
    static void eriroot26(const double*, double*, double*, const int);
    static void eriroot27(const double*, double*, double*, const int);
    static void eriroot28(const double*, double*, double*, const int);
    static void eriroot29(const double*, double*, double*, const int);
    static void eriroot30(const double*, double*, double*, const int);
    static void eriroot31(const double*, double*, double*, const int);
    static void eriroot32(const double*, double*, double*, const int);
    static void eriroot33(const double*, double*, double*, const int);
    static void eriroot34(const double*, double*, double*, const int);
    static void eriroot35(const double*, double*, double*, const int);
    static void eriroot36(const double*, double*, double*, const int);
    static void eriroot37(const double*, double*, double*, const int);
    static void eriroot38(const double*, double*, double*, const int);
    static void eriroot39(const double*, double*, double*, const int);
    static void eriroot40(const double*, double*, double*, const int);
    static void eriroot41(const double*, double*, double*, const int);
    static void eriroot42(const double*, double*, double*, const int);
    static void eriroot43(const double*, double*, double*, const int);
    static void eriroot44(const double*, double*, double*, const int);
    static void eriroot45(const double*, double*, double*, const int);
    static void eriroot46(const double*, double*, double*, const int);
    static void eriroot47(const double*, double*, double*, const int);
    static void eriroot48(const double*, double*, double*, const int);
    static void eriroot49(const double*, double*, double*, const int);
    static void eriroot50(const double*, double*, double*, const int);


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
      rfunc[14] = &eriroot14;
      rfunc[15] = &eriroot15;
      rfunc[16] = &eriroot16;
      rfunc[17] = &eriroot17;
      rfunc[18] = &eriroot18;
      rfunc[19] = &eriroot19;
      rfunc[20] = &eriroot20;
      rfunc[21] = &eriroot21;
      rfunc[22] = &eriroot22;
      rfunc[23] = &eriroot23;
      rfunc[24] = &eriroot24;
      rfunc[25] = &eriroot25;
      rfunc[26] = &eriroot26;
      rfunc[27] = &eriroot27;
      rfunc[28] = &eriroot28;
      rfunc[29] = &eriroot29;
      rfunc[30] = &eriroot30;
      rfunc[31] = &eriroot31;
      rfunc[32] = &eriroot32;
      rfunc[33] = &eriroot33;
      rfunc[34] = &eriroot34;
      rfunc[35] = &eriroot35;
      rfunc[36] = &eriroot36;
      rfunc[37] = &eriroot37;
      rfunc[38] = &eriroot38;
      rfunc[39] = &eriroot39;
      rfunc[40] = &eriroot40;
      rfunc[41] = &eriroot41;
      rfunc[42] = &eriroot42;
      rfunc[43] = &eriroot43;
      rfunc[44] = &eriroot44;
      rfunc[45] = &eriroot45;
      rfunc[46] = &eriroot46;
      rfunc[47] = &eriroot47;
      rfunc[48] = &eriroot48;
      rfunc[49] = &eriroot49;
      rfunc[50] = &eriroot50;
    }

    void root(const int i, const double* a1, double* a2, double* a3, const int a4) const { rfunc[i](a1, a2, a3, a4); }

};

const static ERIRootList eriroot__;

}

#endif

