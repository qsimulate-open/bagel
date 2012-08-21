//
// Newint - Parallel electron correlation program.
// Filename: _hrr_40_31.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include "hrrlist.h"
#include <array>
#include <algorithm>

using namespace std;

void HRRList::perform_HRR_40_31(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    const double* current_data = &data_start[c * 25];
    double* current_out = &data_out[c * 30];
   {
     //current index a: xxx
      const double a0_0 = current_data[0];
      const double ax_0 = current_data[10];
      const double ay_0 = current_data[11];
      const double az_0 = current_data[15];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxy
      const double a0_0 = current_data[1];
      const double ax_0 = current_data[11];
      const double ay_0 = current_data[12];
      const double az_0 = current_data[16];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyy
      const double a0_0 = current_data[2];
      const double ax_0 = current_data[12];
      const double ay_0 = current_data[13];
      const double az_0 = current_data[17];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyy
      const double a0_0 = current_data[3];
      const double ax_0 = current_data[13];
      const double ay_0 = current_data[14];
      const double az_0 = current_data[18];

      current_out[9] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[10] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[11] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxz
      const double a0_0 = current_data[4];
      const double ax_0 = current_data[15];
      const double ay_0 = current_data[16];
      const double az_0 = current_data[19];

      current_out[12] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[13] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[14] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyz
      const double a0_0 = current_data[5];
      const double ax_0 = current_data[16];
      const double ay_0 = current_data[17];
      const double az_0 = current_data[20];

      current_out[15] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[16] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[17] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyz
      const double a0_0 = current_data[6];
      const double ax_0 = current_data[17];
      const double ay_0 = current_data[18];
      const double az_0 = current_data[21];

      current_out[18] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[19] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[20] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xzz
      const double a0_0 = current_data[7];
      const double ax_0 = current_data[19];
      const double ay_0 = current_data[20];
      const double az_0 = current_data[22];

      current_out[21] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[22] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[23] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yzz
      const double a0_0 = current_data[8];
      const double ax_0 = current_data[20];
      const double ay_0 = current_data[21];
      const double az_0 = current_data[23];

      current_out[24] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[25] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[26] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: zzz
      const double a0_0 = current_data[9];
      const double ax_0 = current_data[22];
      const double ay_0 = current_data[23];
      const double az_0 = current_data[24];

      current_out[27] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[28] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[29] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}

