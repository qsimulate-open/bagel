//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _hrr_60_51.cc
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

#include <src/integral/hrrlist.h>
#include <array>
#include <algorithm>

using namespace std;
using namespace bagel;

void HRRList::perform_HRR_60_51(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 49];
    auto current_out = &data_out[c * 63];
   {
     //current index a: xxxxx
      const auto a0_0 = current_data[0];
      const auto ax_0 = current_data[21];
      const auto ay_0 = current_data[22];
      const auto az_0 = current_data[28];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxy
      const auto a0_0 = current_data[1];
      const auto ax_0 = current_data[22];
      const auto ay_0 = current_data[23];
      const auto az_0 = current_data[29];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyy
      const auto a0_0 = current_data[2];
      const auto ax_0 = current_data[23];
      const auto ay_0 = current_data[24];
      const auto az_0 = current_data[30];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyy
      const auto a0_0 = current_data[3];
      const auto ax_0 = current_data[24];
      const auto ay_0 = current_data[25];
      const auto az_0 = current_data[31];

      current_out[9] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[10] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[11] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyy
      const auto a0_0 = current_data[4];
      const auto ax_0 = current_data[25];
      const auto ay_0 = current_data[26];
      const auto az_0 = current_data[32];

      current_out[12] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[13] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[14] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyy
      const auto a0_0 = current_data[5];
      const auto ax_0 = current_data[26];
      const auto ay_0 = current_data[27];
      const auto az_0 = current_data[33];

      current_out[15] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[16] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[17] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxz
      const auto a0_0 = current_data[6];
      const auto ax_0 = current_data[28];
      const auto ay_0 = current_data[29];
      const auto az_0 = current_data[34];

      current_out[18] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[19] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[20] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyz
      const auto a0_0 = current_data[7];
      const auto ax_0 = current_data[29];
      const auto ay_0 = current_data[30];
      const auto az_0 = current_data[35];

      current_out[21] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[22] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[23] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyz
      const auto a0_0 = current_data[8];
      const auto ax_0 = current_data[30];
      const auto ay_0 = current_data[31];
      const auto az_0 = current_data[36];

      current_out[24] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[25] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[26] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyz
      const auto a0_0 = current_data[9];
      const auto ax_0 = current_data[31];
      const auto ay_0 = current_data[32];
      const auto az_0 = current_data[37];

      current_out[27] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[28] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[29] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyz
      const auto a0_0 = current_data[10];
      const auto ax_0 = current_data[32];
      const auto ay_0 = current_data[33];
      const auto az_0 = current_data[38];

      current_out[30] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[31] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[32] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxzz
      const auto a0_0 = current_data[11];
      const auto ax_0 = current_data[34];
      const auto ay_0 = current_data[35];
      const auto az_0 = current_data[39];

      current_out[33] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[34] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[35] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyzz
      const auto a0_0 = current_data[12];
      const auto ax_0 = current_data[35];
      const auto ay_0 = current_data[36];
      const auto az_0 = current_data[40];

      current_out[36] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[37] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[38] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyzz
      const auto a0_0 = current_data[13];
      const auto ax_0 = current_data[36];
      const auto ay_0 = current_data[37];
      const auto az_0 = current_data[41];

      current_out[39] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[40] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[41] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyzz
      const auto a0_0 = current_data[14];
      const auto ax_0 = current_data[37];
      const auto ay_0 = current_data[38];
      const auto az_0 = current_data[42];

      current_out[42] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[43] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[44] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxzzz
      const auto a0_0 = current_data[15];
      const auto ax_0 = current_data[39];
      const auto ay_0 = current_data[40];
      const auto az_0 = current_data[43];

      current_out[45] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[46] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[47] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyzzz
      const auto a0_0 = current_data[16];
      const auto ax_0 = current_data[40];
      const auto ay_0 = current_data[41];
      const auto az_0 = current_data[44];

      current_out[48] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[49] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[50] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyzzz
      const auto a0_0 = current_data[17];
      const auto ax_0 = current_data[41];
      const auto ay_0 = current_data[42];
      const auto az_0 = current_data[45];

      current_out[51] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[52] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[53] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xzzzz
      const auto a0_0 = current_data[18];
      const auto ax_0 = current_data[43];
      const auto ay_0 = current_data[44];
      const auto az_0 = current_data[46];

      current_out[54] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[55] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[56] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yzzzz
      const auto a0_0 = current_data[19];
      const auto ax_0 = current_data[44];
      const auto ay_0 = current_data[45];
      const auto az_0 = current_data[47];

      current_out[57] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[58] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[59] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: zzzzz
      const auto a0_0 = current_data[20];
      const auto ax_0 = current_data[46];
      const auto ay_0 = current_data[47];
      const auto az_0 = current_data[48];

      current_out[60] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[61] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[62] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}


void CHRRList::perform_HRR_60_51(const int nloop, const complex<double>* data_start, const array<double,3>& AB, complex<double>* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 49];
    auto current_out = &data_out[c * 63];
   {
     //current index a: xxxxx
      const auto a0_0 = current_data[0];
      const auto ax_0 = current_data[21];
      const auto ay_0 = current_data[22];
      const auto az_0 = current_data[28];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxy
      const auto a0_0 = current_data[1];
      const auto ax_0 = current_data[22];
      const auto ay_0 = current_data[23];
      const auto az_0 = current_data[29];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyy
      const auto a0_0 = current_data[2];
      const auto ax_0 = current_data[23];
      const auto ay_0 = current_data[24];
      const auto az_0 = current_data[30];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyy
      const auto a0_0 = current_data[3];
      const auto ax_0 = current_data[24];
      const auto ay_0 = current_data[25];
      const auto az_0 = current_data[31];

      current_out[9] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[10] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[11] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyy
      const auto a0_0 = current_data[4];
      const auto ax_0 = current_data[25];
      const auto ay_0 = current_data[26];
      const auto az_0 = current_data[32];

      current_out[12] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[13] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[14] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyy
      const auto a0_0 = current_data[5];
      const auto ax_0 = current_data[26];
      const auto ay_0 = current_data[27];
      const auto az_0 = current_data[33];

      current_out[15] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[16] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[17] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxz
      const auto a0_0 = current_data[6];
      const auto ax_0 = current_data[28];
      const auto ay_0 = current_data[29];
      const auto az_0 = current_data[34];

      current_out[18] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[19] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[20] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyz
      const auto a0_0 = current_data[7];
      const auto ax_0 = current_data[29];
      const auto ay_0 = current_data[30];
      const auto az_0 = current_data[35];

      current_out[21] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[22] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[23] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyz
      const auto a0_0 = current_data[8];
      const auto ax_0 = current_data[30];
      const auto ay_0 = current_data[31];
      const auto az_0 = current_data[36];

      current_out[24] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[25] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[26] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyz
      const auto a0_0 = current_data[9];
      const auto ax_0 = current_data[31];
      const auto ay_0 = current_data[32];
      const auto az_0 = current_data[37];

      current_out[27] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[28] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[29] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyz
      const auto a0_0 = current_data[10];
      const auto ax_0 = current_data[32];
      const auto ay_0 = current_data[33];
      const auto az_0 = current_data[38];

      current_out[30] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[31] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[32] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxzz
      const auto a0_0 = current_data[11];
      const auto ax_0 = current_data[34];
      const auto ay_0 = current_data[35];
      const auto az_0 = current_data[39];

      current_out[33] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[34] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[35] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyzz
      const auto a0_0 = current_data[12];
      const auto ax_0 = current_data[35];
      const auto ay_0 = current_data[36];
      const auto az_0 = current_data[40];

      current_out[36] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[37] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[38] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyzz
      const auto a0_0 = current_data[13];
      const auto ax_0 = current_data[36];
      const auto ay_0 = current_data[37];
      const auto az_0 = current_data[41];

      current_out[39] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[40] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[41] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyzz
      const auto a0_0 = current_data[14];
      const auto ax_0 = current_data[37];
      const auto ay_0 = current_data[38];
      const auto az_0 = current_data[42];

      current_out[42] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[43] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[44] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxzzz
      const auto a0_0 = current_data[15];
      const auto ax_0 = current_data[39];
      const auto ay_0 = current_data[40];
      const auto az_0 = current_data[43];

      current_out[45] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[46] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[47] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyzzz
      const auto a0_0 = current_data[16];
      const auto ax_0 = current_data[40];
      const auto ay_0 = current_data[41];
      const auto az_0 = current_data[44];

      current_out[48] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[49] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[50] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyzzz
      const auto a0_0 = current_data[17];
      const auto ax_0 = current_data[41];
      const auto ay_0 = current_data[42];
      const auto az_0 = current_data[45];

      current_out[51] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[52] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[53] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xzzzz
      const auto a0_0 = current_data[18];
      const auto ax_0 = current_data[43];
      const auto ay_0 = current_data[44];
      const auto az_0 = current_data[46];

      current_out[54] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[55] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[56] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yzzzz
      const auto a0_0 = current_data[19];
      const auto ax_0 = current_data[44];
      const auto ay_0 = current_data[45];
      const auto az_0 = current_data[47];

      current_out[57] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[58] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[59] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: zzzzz
      const auto a0_0 = current_data[20];
      const auto ax_0 = current_data[46];
      const auto ay_0 = current_data[47];
      const auto az_0 = current_data[48];

      current_out[60] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[61] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[62] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}

