//
// BAGEL - Parallel electron correlation program.
// Filename: lebedevlist.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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


#ifndef __SRC_KS_LEBEDEVLIST_H
#define __SRC_KS_LEBEDEVLIST_H

#include <functional>
#include <map>
#include <cassert>

namespace bagel {

struct LebedevList  {
  private:
    std::map<int, int> map_;
    std::function<void (double*, double*, double*, double*)> rfunc[32];

    static void ld0006(double* x, double* y, double* z, double* weight);
    static void ld0014(double* x, double* y, double* z, double* weight);
    static void ld0026(double* x, double* y, double* z, double* weight);
    static void ld0038(double* x, double* y, double* z, double* weight);
    static void ld0050(double* x, double* y, double* z, double* weight);
    static void ld0074(double* x, double* y, double* z, double* weight);
    static void ld0086(double* x, double* y, double* z, double* weight);
    static void ld0110(double* x, double* y, double* z, double* weight);
    static void ld0146(double* x, double* y, double* z, double* weight);
    static void ld0170(double* x, double* y, double* z, double* weight);
    static void ld0194(double* x, double* y, double* z, double* weight);
    static void ld0230(double* x, double* y, double* z, double* weight);
    static void ld0266(double* x, double* y, double* z, double* weight);
    static void ld0302(double* x, double* y, double* z, double* weight);
    static void ld0350(double* x, double* y, double* z, double* weight);
    static void ld0434(double* x, double* y, double* z, double* weight);
    static void ld0590(double* x, double* y, double* z, double* weight);
    static void ld0770(double* x, double* y, double* z, double* weight);
    static void ld0974(double* x, double* y, double* z, double* weight);
    static void ld1202(double* x, double* y, double* z, double* weight);
    static void ld1454(double* x, double* y, double* z, double* weight);
    static void ld1730(double* x, double* y, double* z, double* weight);
    static void ld2030(double* x, double* y, double* z, double* weight);
    static void ld2354(double* x, double* y, double* z, double* weight);
    static void ld2702(double* x, double* y, double* z, double* weight);
    static void ld3074(double* x, double* y, double* z, double* weight);
    static void ld3470(double* x, double* y, double* z, double* weight);
    static void ld3890(double* x, double* y, double* z, double* weight);
    static void ld4334(double* x, double* y, double* z, double* weight);
    static void ld4802(double* x, double* y, double* z, double* weight);
    static void ld5294(double* x, double* y, double* z, double* weight);
    static void ld5810(double* x, double* y, double* z, double* weight);

  public:
    LebedevList() {
      rfunc[0] = &ld0006; map_.insert(std::make_pair(6, 0));
      rfunc[1] = &ld0014; map_.insert(std::make_pair(14, 1));
      rfunc[2] = &ld0026; map_.insert(std::make_pair(26, 2));
      rfunc[3] = &ld0038; map_.insert(std::make_pair(38, 3));
      rfunc[4] = &ld0050; map_.insert(std::make_pair(50, 4));
      rfunc[5] = &ld0074; map_.insert(std::make_pair(74, 5));
      rfunc[6] = &ld0086; map_.insert(std::make_pair(86, 6));
      rfunc[7] = &ld0110; map_.insert(std::make_pair(110, 7));
      rfunc[8] = &ld0146; map_.insert(std::make_pair(146, 8));
      rfunc[9] = &ld0170; map_.insert(std::make_pair(170, 9));
      rfunc[10] = &ld0194; map_.insert(std::make_pair(194, 10));
      rfunc[11] = &ld0230; map_.insert(std::make_pair(230, 11));
      rfunc[12] = &ld0266; map_.insert(std::make_pair(266, 12));
      rfunc[13] = &ld0302; map_.insert(std::make_pair(302, 13));
      rfunc[14] = &ld0350; map_.insert(std::make_pair(350, 14));
      rfunc[15] = &ld0434; map_.insert(std::make_pair(434, 15));
      rfunc[16] = &ld0590; map_.insert(std::make_pair(590, 16));
      rfunc[17] = &ld0770; map_.insert(std::make_pair(770, 17));
      rfunc[18] = &ld0974; map_.insert(std::make_pair(974, 18));
      rfunc[19] = &ld1202; map_.insert(std::make_pair(1202, 19));
      rfunc[20] = &ld1454; map_.insert(std::make_pair(1454, 20));
      rfunc[21] = &ld1730; map_.insert(std::make_pair(1730, 21));
      rfunc[22] = &ld2030; map_.insert(std::make_pair(2030, 22));
      rfunc[23] = &ld2354; map_.insert(std::make_pair(2354, 23));
      rfunc[24] = &ld2702; map_.insert(std::make_pair(2702, 24));
      rfunc[25] = &ld3074; map_.insert(std::make_pair(3074, 25));
      rfunc[26] = &ld3470; map_.insert(std::make_pair(3470, 26));
      rfunc[27] = &ld3890; map_.insert(std::make_pair(3890, 27));
      rfunc[28] = &ld4334; map_.insert(std::make_pair(4334, 28));
      rfunc[29] = &ld4802; map_.insert(std::make_pair(4802, 29));
      rfunc[30] = &ld5294; map_.insert(std::make_pair(5294, 30));
      rfunc[31] = &ld5810; map_.insert(std::make_pair(5810, 31));
    }

    void root(const int ngrid, double* a1, double* a2, double* a3, double* a4) const {
      auto iter = map_.find(ngrid);
      assert(iter != map_.end());
      rfunc[iter->second](a1, a2, a3, a4);
    }

};

}

#endif

