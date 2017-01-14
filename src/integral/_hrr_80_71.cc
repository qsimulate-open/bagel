//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _hrr_80_71.cc
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

#ifdef COMPILE_J_ORB
#include <src/integral/hrrlist.h>
#include <array>
#include <algorithm>

using namespace std;
using namespace bagel;

void HRRList::perform_HRR_80_71(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 81];
    auto current_out = &data_out[c * 108];
   {
     //current index a: xxxxxxx
      const auto a0_0 = current_data[0];
      const auto ax_0 = current_data[36];
      const auto ay_0 = current_data[37];
      const auto az_0 = current_data[45];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxxxy
      const auto a0_0 = current_data[1];
      const auto ax_0 = current_data[37];
      const auto ay_0 = current_data[38];
      const auto az_0 = current_data[46];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxxyy
      const auto a0_0 = current_data[2];
      const auto ax_0 = current_data[38];
      const auto ay_0 = current_data[39];
      const auto az_0 = current_data[47];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxyyy
      const auto a0_0 = current_data[3];
      const auto ax_0 = current_data[39];
      const auto ay_0 = current_data[40];
      const auto az_0 = current_data[48];

      current_out[9] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[10] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[11] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyyyy
      const auto a0_0 = current_data[4];
      const auto ax_0 = current_data[40];
      const auto ay_0 = current_data[41];
      const auto az_0 = current_data[49];

      current_out[12] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[13] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[14] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyyyy
      const auto a0_0 = current_data[5];
      const auto ax_0 = current_data[41];
      const auto ay_0 = current_data[42];
      const auto az_0 = current_data[50];

      current_out[15] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[16] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[17] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyyyy
      const auto a0_0 = current_data[6];
      const auto ax_0 = current_data[42];
      const auto ay_0 = current_data[43];
      const auto az_0 = current_data[51];

      current_out[18] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[19] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[20] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyyyy
      const auto a0_0 = current_data[7];
      const auto ax_0 = current_data[43];
      const auto ay_0 = current_data[44];
      const auto az_0 = current_data[52];

      current_out[21] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[22] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[23] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxxxz
      const auto a0_0 = current_data[8];
      const auto ax_0 = current_data[45];
      const auto ay_0 = current_data[46];
      const auto az_0 = current_data[53];

      current_out[24] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[25] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[26] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxxyz
      const auto a0_0 = current_data[9];
      const auto ax_0 = current_data[46];
      const auto ay_0 = current_data[47];
      const auto az_0 = current_data[54];

      current_out[27] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[28] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[29] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxyyz
      const auto a0_0 = current_data[10];
      const auto ax_0 = current_data[47];
      const auto ay_0 = current_data[48];
      const auto az_0 = current_data[55];

      current_out[30] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[31] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[32] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyyyz
      const auto a0_0 = current_data[11];
      const auto ax_0 = current_data[48];
      const auto ay_0 = current_data[49];
      const auto az_0 = current_data[56];

      current_out[33] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[34] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[35] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyyyz
      const auto a0_0 = current_data[12];
      const auto ax_0 = current_data[49];
      const auto ay_0 = current_data[50];
      const auto az_0 = current_data[57];

      current_out[36] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[37] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[38] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyyyz
      const auto a0_0 = current_data[13];
      const auto ax_0 = current_data[50];
      const auto ay_0 = current_data[51];
      const auto az_0 = current_data[58];

      current_out[39] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[40] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[41] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyyyz
      const auto a0_0 = current_data[14];
      const auto ax_0 = current_data[51];
      const auto ay_0 = current_data[52];
      const auto az_0 = current_data[59];

      current_out[42] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[43] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[44] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxxzz
      const auto a0_0 = current_data[15];
      const auto ax_0 = current_data[53];
      const auto ay_0 = current_data[54];
      const auto az_0 = current_data[60];

      current_out[45] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[46] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[47] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxyzz
      const auto a0_0 = current_data[16];
      const auto ax_0 = current_data[54];
      const auto ay_0 = current_data[55];
      const auto az_0 = current_data[61];

      current_out[48] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[49] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[50] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyyzz
      const auto a0_0 = current_data[17];
      const auto ax_0 = current_data[55];
      const auto ay_0 = current_data[56];
      const auto az_0 = current_data[62];

      current_out[51] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[52] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[53] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyyzz
      const auto a0_0 = current_data[18];
      const auto ax_0 = current_data[56];
      const auto ay_0 = current_data[57];
      const auto az_0 = current_data[63];

      current_out[54] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[55] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[56] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyyzz
      const auto a0_0 = current_data[19];
      const auto ax_0 = current_data[57];
      const auto ay_0 = current_data[58];
      const auto az_0 = current_data[64];

      current_out[57] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[58] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[59] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyyzz
      const auto a0_0 = current_data[20];
      const auto ax_0 = current_data[58];
      const auto ay_0 = current_data[59];
      const auto az_0 = current_data[65];

      current_out[60] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[61] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[62] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxzzz
      const auto a0_0 = current_data[21];
      const auto ax_0 = current_data[60];
      const auto ay_0 = current_data[61];
      const auto az_0 = current_data[66];

      current_out[63] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[64] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[65] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyzzz
      const auto a0_0 = current_data[22];
      const auto ax_0 = current_data[61];
      const auto ay_0 = current_data[62];
      const auto az_0 = current_data[67];

      current_out[66] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[67] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[68] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyzzz
      const auto a0_0 = current_data[23];
      const auto ax_0 = current_data[62];
      const auto ay_0 = current_data[63];
      const auto az_0 = current_data[68];

      current_out[69] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[70] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[71] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyzzz
      const auto a0_0 = current_data[24];
      const auto ax_0 = current_data[63];
      const auto ay_0 = current_data[64];
      const auto az_0 = current_data[69];

      current_out[72] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[73] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[74] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyzzz
      const auto a0_0 = current_data[25];
      const auto ax_0 = current_data[64];
      const auto ay_0 = current_data[65];
      const auto az_0 = current_data[70];

      current_out[75] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[76] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[77] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxzzzz
      const auto a0_0 = current_data[26];
      const auto ax_0 = current_data[66];
      const auto ay_0 = current_data[67];
      const auto az_0 = current_data[71];

      current_out[78] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[79] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[80] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyzzzz
      const auto a0_0 = current_data[27];
      const auto ax_0 = current_data[67];
      const auto ay_0 = current_data[68];
      const auto az_0 = current_data[72];

      current_out[81] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[82] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[83] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyzzzz
      const auto a0_0 = current_data[28];
      const auto ax_0 = current_data[68];
      const auto ay_0 = current_data[69];
      const auto az_0 = current_data[73];

      current_out[84] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[85] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[86] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyzzzz
      const auto a0_0 = current_data[29];
      const auto ax_0 = current_data[69];
      const auto ay_0 = current_data[70];
      const auto az_0 = current_data[74];

      current_out[87] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[88] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[89] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxzzzzz
      const auto a0_0 = current_data[30];
      const auto ax_0 = current_data[71];
      const auto ay_0 = current_data[72];
      const auto az_0 = current_data[75];

      current_out[90] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[91] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[92] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyzzzzz
      const auto a0_0 = current_data[31];
      const auto ax_0 = current_data[72];
      const auto ay_0 = current_data[73];
      const auto az_0 = current_data[76];

      current_out[93] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[94] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[95] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyzzzzz
      const auto a0_0 = current_data[32];
      const auto ax_0 = current_data[73];
      const auto ay_0 = current_data[74];
      const auto az_0 = current_data[77];

      current_out[96] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[97] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[98] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xzzzzzz
      const auto a0_0 = current_data[33];
      const auto ax_0 = current_data[75];
      const auto ay_0 = current_data[76];
      const auto az_0 = current_data[78];

      current_out[99] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[100] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[101] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yzzzzzz
      const auto a0_0 = current_data[34];
      const auto ax_0 = current_data[76];
      const auto ay_0 = current_data[77];
      const auto az_0 = current_data[79];

      current_out[102] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[103] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[104] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: zzzzzzz
      const auto a0_0 = current_data[35];
      const auto ax_0 = current_data[78];
      const auto ay_0 = current_data[79];
      const auto az_0 = current_data[80];

      current_out[105] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[106] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[107] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}


void CHRRList::perform_HRR_80_71(const int nloop, const complex<double>* data_start, const array<double,3>& AB, complex<double>* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 81];
    auto current_out = &data_out[c * 108];
   {
     //current index a: xxxxxxx
      const auto a0_0 = current_data[0];
      const auto ax_0 = current_data[36];
      const auto ay_0 = current_data[37];
      const auto az_0 = current_data[45];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxxxy
      const auto a0_0 = current_data[1];
      const auto ax_0 = current_data[37];
      const auto ay_0 = current_data[38];
      const auto az_0 = current_data[46];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxxyy
      const auto a0_0 = current_data[2];
      const auto ax_0 = current_data[38];
      const auto ay_0 = current_data[39];
      const auto az_0 = current_data[47];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxyyy
      const auto a0_0 = current_data[3];
      const auto ax_0 = current_data[39];
      const auto ay_0 = current_data[40];
      const auto az_0 = current_data[48];

      current_out[9] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[10] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[11] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyyyy
      const auto a0_0 = current_data[4];
      const auto ax_0 = current_data[40];
      const auto ay_0 = current_data[41];
      const auto az_0 = current_data[49];

      current_out[12] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[13] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[14] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyyyy
      const auto a0_0 = current_data[5];
      const auto ax_0 = current_data[41];
      const auto ay_0 = current_data[42];
      const auto az_0 = current_data[50];

      current_out[15] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[16] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[17] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyyyy
      const auto a0_0 = current_data[6];
      const auto ax_0 = current_data[42];
      const auto ay_0 = current_data[43];
      const auto az_0 = current_data[51];

      current_out[18] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[19] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[20] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyyyy
      const auto a0_0 = current_data[7];
      const auto ax_0 = current_data[43];
      const auto ay_0 = current_data[44];
      const auto az_0 = current_data[52];

      current_out[21] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[22] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[23] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxxxz
      const auto a0_0 = current_data[8];
      const auto ax_0 = current_data[45];
      const auto ay_0 = current_data[46];
      const auto az_0 = current_data[53];

      current_out[24] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[25] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[26] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxxyz
      const auto a0_0 = current_data[9];
      const auto ax_0 = current_data[46];
      const auto ay_0 = current_data[47];
      const auto az_0 = current_data[54];

      current_out[27] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[28] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[29] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxyyz
      const auto a0_0 = current_data[10];
      const auto ax_0 = current_data[47];
      const auto ay_0 = current_data[48];
      const auto az_0 = current_data[55];

      current_out[30] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[31] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[32] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyyyz
      const auto a0_0 = current_data[11];
      const auto ax_0 = current_data[48];
      const auto ay_0 = current_data[49];
      const auto az_0 = current_data[56];

      current_out[33] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[34] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[35] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyyyz
      const auto a0_0 = current_data[12];
      const auto ax_0 = current_data[49];
      const auto ay_0 = current_data[50];
      const auto az_0 = current_data[57];

      current_out[36] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[37] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[38] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyyyz
      const auto a0_0 = current_data[13];
      const auto ax_0 = current_data[50];
      const auto ay_0 = current_data[51];
      const auto az_0 = current_data[58];

      current_out[39] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[40] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[41] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyyyz
      const auto a0_0 = current_data[14];
      const auto ax_0 = current_data[51];
      const auto ay_0 = current_data[52];
      const auto az_0 = current_data[59];

      current_out[42] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[43] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[44] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxxzz
      const auto a0_0 = current_data[15];
      const auto ax_0 = current_data[53];
      const auto ay_0 = current_data[54];
      const auto az_0 = current_data[60];

      current_out[45] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[46] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[47] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxyzz
      const auto a0_0 = current_data[16];
      const auto ax_0 = current_data[54];
      const auto ay_0 = current_data[55];
      const auto az_0 = current_data[61];

      current_out[48] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[49] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[50] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyyzz
      const auto a0_0 = current_data[17];
      const auto ax_0 = current_data[55];
      const auto ay_0 = current_data[56];
      const auto az_0 = current_data[62];

      current_out[51] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[52] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[53] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyyzz
      const auto a0_0 = current_data[18];
      const auto ax_0 = current_data[56];
      const auto ay_0 = current_data[57];
      const auto az_0 = current_data[63];

      current_out[54] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[55] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[56] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyyzz
      const auto a0_0 = current_data[19];
      const auto ax_0 = current_data[57];
      const auto ay_0 = current_data[58];
      const auto az_0 = current_data[64];

      current_out[57] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[58] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[59] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyyzz
      const auto a0_0 = current_data[20];
      const auto ax_0 = current_data[58];
      const auto ay_0 = current_data[59];
      const auto az_0 = current_data[65];

      current_out[60] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[61] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[62] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxxzzz
      const auto a0_0 = current_data[21];
      const auto ax_0 = current_data[60];
      const auto ay_0 = current_data[61];
      const auto az_0 = current_data[66];

      current_out[63] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[64] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[65] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxyzzz
      const auto a0_0 = current_data[22];
      const auto ax_0 = current_data[61];
      const auto ay_0 = current_data[62];
      const auto az_0 = current_data[67];

      current_out[66] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[67] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[68] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyyzzz
      const auto a0_0 = current_data[23];
      const auto ax_0 = current_data[62];
      const auto ay_0 = current_data[63];
      const auto az_0 = current_data[68];

      current_out[69] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[70] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[71] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyyzzz
      const auto a0_0 = current_data[24];
      const auto ax_0 = current_data[63];
      const auto ay_0 = current_data[64];
      const auto az_0 = current_data[69];

      current_out[72] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[73] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[74] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyyzzz
      const auto a0_0 = current_data[25];
      const auto ax_0 = current_data[64];
      const auto ay_0 = current_data[65];
      const auto az_0 = current_data[70];

      current_out[75] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[76] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[77] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxxzzzz
      const auto a0_0 = current_data[26];
      const auto ax_0 = current_data[66];
      const auto ay_0 = current_data[67];
      const auto az_0 = current_data[71];

      current_out[78] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[79] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[80] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxyzzzz
      const auto a0_0 = current_data[27];
      const auto ax_0 = current_data[67];
      const auto ay_0 = current_data[68];
      const auto az_0 = current_data[72];

      current_out[81] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[82] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[83] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyyzzzz
      const auto a0_0 = current_data[28];
      const auto ax_0 = current_data[68];
      const auto ay_0 = current_data[69];
      const auto az_0 = current_data[73];

      current_out[84] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[85] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[86] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyyzzzz
      const auto a0_0 = current_data[29];
      const auto ax_0 = current_data[69];
      const auto ay_0 = current_data[70];
      const auto az_0 = current_data[74];

      current_out[87] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[88] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[89] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxzzzzz
      const auto a0_0 = current_data[30];
      const auto ax_0 = current_data[71];
      const auto ay_0 = current_data[72];
      const auto az_0 = current_data[75];

      current_out[90] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[91] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[92] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyzzzzz
      const auto a0_0 = current_data[31];
      const auto ax_0 = current_data[72];
      const auto ay_0 = current_data[73];
      const auto az_0 = current_data[76];

      current_out[93] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[94] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[95] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyzzzzz
      const auto a0_0 = current_data[32];
      const auto ax_0 = current_data[73];
      const auto ay_0 = current_data[74];
      const auto az_0 = current_data[77];

      current_out[96] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[97] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[98] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xzzzzzz
      const auto a0_0 = current_data[33];
      const auto ax_0 = current_data[75];
      const auto ay_0 = current_data[76];
      const auto az_0 = current_data[78];

      current_out[99] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[100] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[101] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yzzzzzz
      const auto a0_0 = current_data[34];
      const auto ax_0 = current_data[76];
      const auto ay_0 = current_data[77];
      const auto az_0 = current_data[79];

      current_out[102] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[103] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[104] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: zzzzzzz
      const auto a0_0 = current_data[35];
      const auto ax_0 = current_data[78];
      const auto ay_0 = current_data[79];
      const auto az_0 = current_data[80];

      current_out[105] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[106] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[107] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}
#endif

