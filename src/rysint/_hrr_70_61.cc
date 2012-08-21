//
// Newint - Parallel electron correlation program.
// Filename: _hrr_70_61.cc
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

void HRRList::perform_HRR_70_61(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    const double* current_data = &data_start[c * 64];
    double* current_out = &data_out[c * 84];
   {
     //current index a: xxxxxx
      const double a0_0 = current_data[0];
      const double ax_0 = current_data[28];
      const double ay_0 = current_data[29];
      const double az_0 = current_data[36];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxxy
      const double a0_0 = current_data[1];
      const double ax_0 = current_data[29];
      const double ay_0 = current_data[30];
      const double az_0 = current_data[37];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxyy
      const double a0_0 = current_data[2];
      const double ax_0 = current_data[30];
      const double ay_0 = current_data[31];
      const double az_0 = current_data[38];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyyy
      const double a0_0 = current_data[3];
      const double ax_0 = current_data[31];
      const double ay_0 = current_data[32];
      const double az_0 = current_data[39];

      current_out[9] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[10] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[11] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyyy
      const double a0_0 = current_data[4];
      const double ax_0 = current_data[32];
      const double ay_0 = current_data[33];
      const double az_0 = current_data[40];

      current_out[12] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[13] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[14] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyyy
      const double a0_0 = current_data[5];
      const double ax_0 = current_data[33];
      const double ay_0 = current_data[34];
      const double az_0 = current_data[41];

      current_out[15] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[16] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[17] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyyy
      const double a0_0 = current_data[6];
      const double ax_0 = current_data[34];
      const double ay_0 = current_data[35];
      const double az_0 = current_data[42];

      current_out[18] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[19] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[20] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxxz
      const double a0_0 = current_data[7];
      const double ax_0 = current_data[36];
      const double ay_0 = current_data[37];
      const double az_0 = current_data[43];

      current_out[21] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[22] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[23] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxyz
      const double a0_0 = current_data[8];
      const double ax_0 = current_data[37];
      const double ay_0 = current_data[38];
      const double az_0 = current_data[44];

      current_out[24] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[25] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[26] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyyz
      const double a0_0 = current_data[9];
      const double ax_0 = current_data[38];
      const double ay_0 = current_data[39];
      const double az_0 = current_data[45];

      current_out[27] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[28] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[29] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyyz
      const double a0_0 = current_data[10];
      const double ax_0 = current_data[39];
      const double ay_0 = current_data[40];
      const double az_0 = current_data[46];

      current_out[30] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[31] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[32] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyyz
      const double a0_0 = current_data[11];
      const double ax_0 = current_data[40];
      const double ay_0 = current_data[41];
      const double az_0 = current_data[47];

      current_out[33] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[34] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[35] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyyz
      const double a0_0 = current_data[12];
      const double ax_0 = current_data[41];
      const double ay_0 = current_data[42];
      const double az_0 = current_data[48];

      current_out[36] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[37] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[38] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxzz
      const double a0_0 = current_data[13];
      const double ax_0 = current_data[43];
      const double ay_0 = current_data[44];
      const double az_0 = current_data[49];

      current_out[39] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[40] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[41] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyzz
      const double a0_0 = current_data[14];
      const double ax_0 = current_data[44];
      const double ay_0 = current_data[45];
      const double az_0 = current_data[50];

      current_out[42] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[43] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[44] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyzz
      const double a0_0 = current_data[15];
      const double ax_0 = current_data[45];
      const double ay_0 = current_data[46];
      const double az_0 = current_data[51];

      current_out[45] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[46] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[47] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyzz
      const double a0_0 = current_data[16];
      const double ax_0 = current_data[46];
      const double ay_0 = current_data[47];
      const double az_0 = current_data[52];

      current_out[48] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[49] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[50] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyzz
      const double a0_0 = current_data[17];
      const double ax_0 = current_data[47];
      const double ay_0 = current_data[48];
      const double az_0 = current_data[53];

      current_out[51] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[52] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[53] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxzzz
      const double a0_0 = current_data[18];
      const double ax_0 = current_data[49];
      const double ay_0 = current_data[50];
      const double az_0 = current_data[54];

      current_out[54] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[55] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[56] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyzzz
      const double a0_0 = current_data[19];
      const double ax_0 = current_data[50];
      const double ay_0 = current_data[51];
      const double az_0 = current_data[55];

      current_out[57] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[58] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[59] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyzzz
      const double a0_0 = current_data[20];
      const double ax_0 = current_data[51];
      const double ay_0 = current_data[52];
      const double az_0 = current_data[56];

      current_out[60] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[61] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[62] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyzzz
      const double a0_0 = current_data[21];
      const double ax_0 = current_data[52];
      const double ay_0 = current_data[53];
      const double az_0 = current_data[57];

      current_out[63] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[64] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[65] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxzzzz
      const double a0_0 = current_data[22];
      const double ax_0 = current_data[54];
      const double ay_0 = current_data[55];
      const double az_0 = current_data[58];

      current_out[66] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[67] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[68] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyzzzz
      const double a0_0 = current_data[23];
      const double ax_0 = current_data[55];
      const double ay_0 = current_data[56];
      const double az_0 = current_data[59];

      current_out[69] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[70] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[71] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyzzzz
      const double a0_0 = current_data[24];
      const double ax_0 = current_data[56];
      const double ay_0 = current_data[57];
      const double az_0 = current_data[60];

      current_out[72] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[73] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[74] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xzzzzz
      const double a0_0 = current_data[25];
      const double ax_0 = current_data[58];
      const double ay_0 = current_data[59];
      const double az_0 = current_data[61];

      current_out[75] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[76] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[77] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yzzzzz
      const double a0_0 = current_data[26];
      const double ax_0 = current_data[59];
      const double ay_0 = current_data[60];
      const double az_0 = current_data[62];

      current_out[78] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[79] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[80] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: zzzzzz
      const double a0_0 = current_data[27];
      const double ax_0 = current_data[61];
      const double ay_0 = current_data[62];
      const double az_0 = current_data[63];

      current_out[81] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[82] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[83] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}

