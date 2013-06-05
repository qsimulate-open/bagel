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


#ifndef __SRC_KS_LEBEDEVLIST_H
#define __SRC_KS_LEBEDEVLIST_H

#include <functional>
#include <map>
#include <cassert>

extern "C" {
  void ld0006_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0014_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0026_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0038_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0050_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0074_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0086_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0110_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0146_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0170_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0194_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0230_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0266_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0302_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0350_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0434_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0590_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0770_(double* x, double* y, double* z, double* weight, const int* n);
  void ld0974_(double* x, double* y, double* z, double* weight, const int* n);
  void ld1202_(double* x, double* y, double* z, double* weight, const int* n);
  void ld1454_(double* x, double* y, double* z, double* weight, const int* n);
  void ld1730_(double* x, double* y, double* z, double* weight, const int* n);
  void ld2030_(double* x, double* y, double* z, double* weight, const int* n);
  void ld2354_(double* x, double* y, double* z, double* weight, const int* n);
  void ld2702_(double* x, double* y, double* z, double* weight, const int* n);
  void ld3074_(double* x, double* y, double* z, double* weight, const int* n);
  void ld3470_(double* x, double* y, double* z, double* weight, const int* n);
  void ld3890_(double* x, double* y, double* z, double* weight, const int* n);
  void ld4334_(double* x, double* y, double* z, double* weight, const int* n);
  void ld4802_(double* x, double* y, double* z, double* weight, const int* n);
  void ld5294_(double* x, double* y, double* z, double* weight, const int* n);
  void ld5810_(double* x, double* y, double* z, double* weight, const int* n);
}

namespace bagel {

struct LebedevList  {
  private:
    std::map<int, int> map_;
    std::function<void (double*, double*, double*, double*, const int*)> rfunc[32];

  public:
    LebedevList() {
      rfunc[0] = &ld0006_; map_.insert(std::make_pair(6, 0));
      rfunc[1] = &ld0014_; map_.insert(std::make_pair(14, 1));
      rfunc[2] = &ld0026_; map_.insert(std::make_pair(26, 2));
      rfunc[3] = &ld0038_; map_.insert(std::make_pair(38, 3));
      rfunc[4] = &ld0050_; map_.insert(std::make_pair(50, 4));
      rfunc[5] = &ld0074_; map_.insert(std::make_pair(74, 5));
      rfunc[6] = &ld0086_; map_.insert(std::make_pair(86, 6));
      rfunc[7] = &ld0110_; map_.insert(std::make_pair(110, 7));
      rfunc[8] = &ld0146_; map_.insert(std::make_pair(146, 8));
      rfunc[9] = &ld0170_; map_.insert(std::make_pair(170, 9));
      rfunc[10] = &ld0194_; map_.insert(std::make_pair(194, 10));
      rfunc[11] = &ld0230_; map_.insert(std::make_pair(230, 11));
      rfunc[12] = &ld0266_; map_.insert(std::make_pair(266, 12));
      rfunc[13] = &ld0302_; map_.insert(std::make_pair(302, 13));
      rfunc[14] = &ld0350_; map_.insert(std::make_pair(350, 14));
      rfunc[15] = &ld0434_; map_.insert(std::make_pair(434, 15));
      rfunc[16] = &ld0590_; map_.insert(std::make_pair(590, 16));
      rfunc[17] = &ld0770_; map_.insert(std::make_pair(770, 17));
      rfunc[18] = &ld0974_; map_.insert(std::make_pair(974, 18));
      rfunc[19] = &ld1202_; map_.insert(std::make_pair(1202, 19));
      rfunc[20] = &ld1454_; map_.insert(std::make_pair(1454, 20));
      rfunc[21] = &ld1730_; map_.insert(std::make_pair(1730, 21));
      rfunc[22] = &ld2030_; map_.insert(std::make_pair(2030, 22));
      rfunc[23] = &ld2354_; map_.insert(std::make_pair(2354, 23));
      rfunc[24] = &ld2702_; map_.insert(std::make_pair(2702, 24));
      rfunc[25] = &ld3074_; map_.insert(std::make_pair(3074, 25));
      rfunc[26] = &ld3470_; map_.insert(std::make_pair(3470, 26));
      rfunc[27] = &ld3890_; map_.insert(std::make_pair(3890, 27));
      rfunc[28] = &ld4334_; map_.insert(std::make_pair(4334, 28));
      rfunc[29] = &ld4802_; map_.insert(std::make_pair(4802, 29));
      rfunc[30] = &ld5294_; map_.insert(std::make_pair(5294, 30));
      rfunc[31] = &ld5810_; map_.insert(std::make_pair(5810, 31));
    }

    void root(const int ngrid, double* a1, double* a2, double* a3, double* a4) const {
      auto iter = map_.find(ngrid);
      assert(iter != map_.end());
      int n;
      rfunc[iter->second](a1, a2, a3, a4, &n);
      assert(ngrid == n);
    }

};

}

#endif

