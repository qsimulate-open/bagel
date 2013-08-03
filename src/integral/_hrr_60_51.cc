//
// BAGEL - Parallel electron correlation program.
// Filename: _hrr_60_51.cc
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

#include <src/integral/hrrlist.h>
#include <array>
#include <algorithm>

using namespace std;
using namespace bagel;

void HRRList::perform_HRR_60_51(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    const double* current_data = &data_start[c * 49];
    double* current_out = &data_out[c * 63];
   {
     //current index a: xxxxx
      const double a0_0 = current_data[0];
      const double ax_0 = current_data[21];
      const double ay_0 = current_data[22];
      const double az_0 = current_data[28];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxy
      const double a0_0 = current_data[1];
      const double ax_0 = current_data[22];
      const double ay_0 = current_data[23];
      const double az_0 = current_data[29];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyy
      const double a0_0 = current_data[2];
      const double ax_0 = current_data[23];
      const double ay_0 = current_data[24];
      const double az_0 = current_data[30];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyy
      const double a0_0 = current_data[3];
      const double ax_0 = current_data[24];
      const double ay_0 = current_data[25];
      const double az_0 = current_data[31];

      current_out[9] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[10] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[11] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyy
      const double a0_0 = current_data[4];
      const double ax_0 = current_data[25];
      const double ay_0 = current_data[26];
      const double az_0 = current_data[32];

      current_out[12] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[13] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[14] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyy
      const double a0_0 = current_data[5];
      const double ax_0 = current_data[26];
      const double ay_0 = current_data[27];
      const double az_0 = current_data[33];

      current_out[15] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[16] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[17] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxz
      const double a0_0 = current_data[6];
      const double ax_0 = current_data[28];
      const double ay_0 = current_data[29];
      const double az_0 = current_data[34];

      current_out[18] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[19] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[20] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyz
      const double a0_0 = current_data[7];
      const double ax_0 = current_data[29];
      const double ay_0 = current_data[30];
      const double az_0 = current_data[35];

      current_out[21] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[22] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[23] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyz
      const double a0_0 = current_data[8];
      const double ax_0 = current_data[30];
      const double ay_0 = current_data[31];
      const double az_0 = current_data[36];

      current_out[24] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[25] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[26] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyz
      const double a0_0 = current_data[9];
      const double ax_0 = current_data[31];
      const double ay_0 = current_data[32];
      const double az_0 = current_data[37];

      current_out[27] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[28] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[29] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyz
      const double a0_0 = current_data[10];
      const double ax_0 = current_data[32];
      const double ay_0 = current_data[33];
      const double az_0 = current_data[38];

      current_out[30] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[31] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[32] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxzz
      const double a0_0 = current_data[11];
      const double ax_0 = current_data[34];
      const double ay_0 = current_data[35];
      const double az_0 = current_data[39];

      current_out[33] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[34] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[35] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyzz
      const double a0_0 = current_data[12];
      const double ax_0 = current_data[35];
      const double ay_0 = current_data[36];
      const double az_0 = current_data[40];

      current_out[36] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[37] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[38] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyzz
      const double a0_0 = current_data[13];
      const double ax_0 = current_data[36];
      const double ay_0 = current_data[37];
      const double az_0 = current_data[41];

      current_out[39] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[40] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[41] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyzz
      const double a0_0 = current_data[14];
      const double ax_0 = current_data[37];
      const double ay_0 = current_data[38];
      const double az_0 = current_data[42];

      current_out[42] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[43] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[44] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxzzz
      const double a0_0 = current_data[15];
      const double ax_0 = current_data[39];
      const double ay_0 = current_data[40];
      const double az_0 = current_data[43];

      current_out[45] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[46] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[47] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyzzz
      const double a0_0 = current_data[16];
      const double ax_0 = current_data[40];
      const double ay_0 = current_data[41];
      const double az_0 = current_data[44];

      current_out[48] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[49] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[50] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyzzz
      const double a0_0 = current_data[17];
      const double ax_0 = current_data[41];
      const double ay_0 = current_data[42];
      const double az_0 = current_data[45];

      current_out[51] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[52] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[53] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xzzzz
      const double a0_0 = current_data[18];
      const double ax_0 = current_data[43];
      const double ay_0 = current_data[44];
      const double az_0 = current_data[46];

      current_out[54] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[55] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[56] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yzzzz
      const double a0_0 = current_data[19];
      const double ax_0 = current_data[44];
      const double ay_0 = current_data[45];
      const double az_0 = current_data[47];

      current_out[57] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[58] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[59] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: zzzzz
      const double a0_0 = current_data[20];
      const double ax_0 = current_data[46];
      const double ay_0 = current_data[47];
      const double az_0 = current_data[48];

      current_out[60] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[61] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[62] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}

