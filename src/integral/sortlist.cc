//
// BAGEL - Parallel electron correlation program.
// Filename: sortlist.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/integral/sortlist.h>

using namespace bagel;

SortList::SortList(const bool spherical) {
  if (spherical) {
    sortfunc[ANG_HRR_END * 0 + 0] = &sort_indices_00_sph;
    sortfunc[ANG_HRR_END * 0 + 1] = &sort_indices_01_sph;
    sortfunc[ANG_HRR_END * 0 + 2] = &sort_indices_02_sph;
    sortfunc[ANG_HRR_END * 0 + 3] = &sort_indices_03_sph;
    sortfunc[ANG_HRR_END * 0 + 4] = &sort_indices_04_sph;
    sortfunc[ANG_HRR_END * 0 + 5] = &sort_indices_05_sph;
    sortfunc[ANG_HRR_END * 0 + 6] = &sort_indices_06_sph;
    sortfunc[ANG_HRR_END * 1 + 1] = &sort_indices_11_sph;
    sortfunc[ANG_HRR_END * 1 + 2] = &sort_indices_12_sph;
    sortfunc[ANG_HRR_END * 1 + 3] = &sort_indices_13_sph;
    sortfunc[ANG_HRR_END * 1 + 4] = &sort_indices_14_sph;
    sortfunc[ANG_HRR_END * 1 + 5] = &sort_indices_15_sph;
    sortfunc[ANG_HRR_END * 1 + 6] = &sort_indices_16_sph;
    sortfunc[ANG_HRR_END * 2 + 2] = &sort_indices_22_sph;
    sortfunc[ANG_HRR_END * 2 + 3] = &sort_indices_23_sph;
    sortfunc[ANG_HRR_END * 2 + 4] = &sort_indices_24_sph;
    sortfunc[ANG_HRR_END * 2 + 5] = &sort_indices_25_sph;
    sortfunc[ANG_HRR_END * 2 + 6] = &sort_indices_26_sph;
    sortfunc[ANG_HRR_END * 3 + 3] = &sort_indices_33_sph;
    sortfunc[ANG_HRR_END * 3 + 4] = &sort_indices_34_sph;
    sortfunc[ANG_HRR_END * 3 + 5] = &sort_indices_35_sph;
    sortfunc[ANG_HRR_END * 3 + 6] = &sort_indices_36_sph;
    sortfunc[ANG_HRR_END * 4 + 4] = &sort_indices_44_sph;
    sortfunc[ANG_HRR_END * 4 + 5] = &sort_indices_45_sph;
    sortfunc[ANG_HRR_END * 4 + 6] = &sort_indices_46_sph;
    sortfunc[ANG_HRR_END * 5 + 5] = &sort_indices_55_sph;
    sortfunc[ANG_HRR_END * 5 + 6] = &sort_indices_56_sph;
    sortfunc[ANG_HRR_END * 6 + 6] = &sort_indices_66_sph;
  } else {
    sortfunc[ANG_HRR_END * 0 + 0] = &sort_indices_00;
    sortfunc[ANG_HRR_END * 0 + 1] = &sort_indices_01;
    sortfunc[ANG_HRR_END * 0 + 2] = &sort_indices_02;
    sortfunc[ANG_HRR_END * 0 + 3] = &sort_indices_03;
    sortfunc[ANG_HRR_END * 0 + 4] = &sort_indices_04;
    sortfunc[ANG_HRR_END * 0 + 5] = &sort_indices_05;
    sortfunc[ANG_HRR_END * 0 + 6] = &sort_indices_06;
    sortfunc[ANG_HRR_END * 1 + 1] = &sort_indices_11;
    sortfunc[ANG_HRR_END * 1 + 2] = &sort_indices_12;
    sortfunc[ANG_HRR_END * 1 + 3] = &sort_indices_13;
    sortfunc[ANG_HRR_END * 1 + 4] = &sort_indices_14;
    sortfunc[ANG_HRR_END * 1 + 5] = &sort_indices_15;
    sortfunc[ANG_HRR_END * 1 + 6] = &sort_indices_16;
    sortfunc[ANG_HRR_END * 2 + 2] = &sort_indices_22;
    sortfunc[ANG_HRR_END * 2 + 3] = &sort_indices_23;
    sortfunc[ANG_HRR_END * 2 + 4] = &sort_indices_24;
    sortfunc[ANG_HRR_END * 2 + 5] = &sort_indices_25;
    sortfunc[ANG_HRR_END * 2 + 6] = &sort_indices_26;
    sortfunc[ANG_HRR_END * 3 + 3] = &sort_indices_33;
    sortfunc[ANG_HRR_END * 3 + 4] = &sort_indices_34;
    sortfunc[ANG_HRR_END * 3 + 5] = &sort_indices_35;
    sortfunc[ANG_HRR_END * 3 + 6] = &sort_indices_36;
    sortfunc[ANG_HRR_END * 4 + 4] = &sort_indices_44;
    sortfunc[ANG_HRR_END * 4 + 5] = &sort_indices_45;
    sortfunc[ANG_HRR_END * 4 + 6] = &sort_indices_46;
    sortfunc[ANG_HRR_END * 5 + 5] = &sort_indices_55;
    sortfunc[ANG_HRR_END * 5 + 6] = &sort_indices_56;
    sortfunc[ANG_HRR_END * 6 + 6] = &sort_indices_66;
  }
}


CSortList::CSortList(const bool spherical) {
  if (spherical) {
    sortfunc[ANG_HRR_END * 0 + 0] = &sort_indices_00_sph;
    sortfunc[ANG_HRR_END * 0 + 1] = &sort_indices_01_sph;
    sortfunc[ANG_HRR_END * 0 + 2] = &sort_indices_02_sph;
    sortfunc[ANG_HRR_END * 0 + 3] = &sort_indices_03_sph;
    sortfunc[ANG_HRR_END * 0 + 4] = &sort_indices_04_sph;
    sortfunc[ANG_HRR_END * 0 + 5] = &sort_indices_05_sph;
    sortfunc[ANG_HRR_END * 0 + 6] = &sort_indices_06_sph;
    sortfunc[ANG_HRR_END * 1 + 1] = &sort_indices_11_sph;
    sortfunc[ANG_HRR_END * 1 + 2] = &sort_indices_12_sph;
    sortfunc[ANG_HRR_END * 1 + 3] = &sort_indices_13_sph;
    sortfunc[ANG_HRR_END * 1 + 4] = &sort_indices_14_sph;
    sortfunc[ANG_HRR_END * 1 + 5] = &sort_indices_15_sph;
    sortfunc[ANG_HRR_END * 1 + 6] = &sort_indices_16_sph;
    sortfunc[ANG_HRR_END * 2 + 2] = &sort_indices_22_sph;
    sortfunc[ANG_HRR_END * 2 + 3] = &sort_indices_23_sph;
    sortfunc[ANG_HRR_END * 2 + 4] = &sort_indices_24_sph;
    sortfunc[ANG_HRR_END * 2 + 5] = &sort_indices_25_sph;
    sortfunc[ANG_HRR_END * 2 + 6] = &sort_indices_26_sph;
    sortfunc[ANG_HRR_END * 3 + 3] = &sort_indices_33_sph;
    sortfunc[ANG_HRR_END * 3 + 4] = &sort_indices_34_sph;
    sortfunc[ANG_HRR_END * 3 + 5] = &sort_indices_35_sph;
    sortfunc[ANG_HRR_END * 3 + 6] = &sort_indices_36_sph;
    sortfunc[ANG_HRR_END * 4 + 4] = &sort_indices_44_sph;
    sortfunc[ANG_HRR_END * 4 + 5] = &sort_indices_45_sph;
    sortfunc[ANG_HRR_END * 4 + 6] = &sort_indices_46_sph;
    sortfunc[ANG_HRR_END * 5 + 5] = &sort_indices_55_sph;
    sortfunc[ANG_HRR_END * 5 + 6] = &sort_indices_56_sph;
    sortfunc[ANG_HRR_END * 6 + 6] = &sort_indices_66_sph;
  } else {
    sortfunc[ANG_HRR_END * 0 + 0] = &sort_indices_00;
    sortfunc[ANG_HRR_END * 0 + 1] = &sort_indices_01;
    sortfunc[ANG_HRR_END * 0 + 2] = &sort_indices_02;
    sortfunc[ANG_HRR_END * 0 + 3] = &sort_indices_03;
    sortfunc[ANG_HRR_END * 0 + 4] = &sort_indices_04;
    sortfunc[ANG_HRR_END * 0 + 5] = &sort_indices_05;
    sortfunc[ANG_HRR_END * 0 + 6] = &sort_indices_06;
    sortfunc[ANG_HRR_END * 1 + 1] = &sort_indices_11;
    sortfunc[ANG_HRR_END * 1 + 2] = &sort_indices_12;
    sortfunc[ANG_HRR_END * 1 + 3] = &sort_indices_13;
    sortfunc[ANG_HRR_END * 1 + 4] = &sort_indices_14;
    sortfunc[ANG_HRR_END * 1 + 5] = &sort_indices_15;
    sortfunc[ANG_HRR_END * 1 + 6] = &sort_indices_16;
    sortfunc[ANG_HRR_END * 2 + 2] = &sort_indices_22;
    sortfunc[ANG_HRR_END * 2 + 3] = &sort_indices_23;
    sortfunc[ANG_HRR_END * 2 + 4] = &sort_indices_24;
    sortfunc[ANG_HRR_END * 2 + 5] = &sort_indices_25;
    sortfunc[ANG_HRR_END * 2 + 6] = &sort_indices_26;
    sortfunc[ANG_HRR_END * 3 + 3] = &sort_indices_33;
    sortfunc[ANG_HRR_END * 3 + 4] = &sort_indices_34;
    sortfunc[ANG_HRR_END * 3 + 5] = &sort_indices_35;
    sortfunc[ANG_HRR_END * 3 + 6] = &sort_indices_36;
    sortfunc[ANG_HRR_END * 4 + 4] = &sort_indices_44;
    sortfunc[ANG_HRR_END * 4 + 5] = &sort_indices_45;
    sortfunc[ANG_HRR_END * 4 + 6] = &sort_indices_46;
    sortfunc[ANG_HRR_END * 5 + 5] = &sort_indices_55;
    sortfunc[ANG_HRR_END * 5 + 6] = &sort_indices_56;
    sortfunc[ANG_HRR_END * 6 + 6] = &sort_indices_66;
  }
}
