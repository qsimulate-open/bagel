//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _hrr_90_72.cc
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

void HRRList::perform_HRR_90_72(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 136];
    auto current_out = &data_out[c * 216];
   {
     //current index a: xxxxxxx
      const auto a0_0 = current_data[0];
      const auto ax_0 = current_data[36];
      const auto ay_0 = current_data[37];
      const auto az_0 = current_data[45];
      const auto axx_0 = current_data[81];
      const auto axy_0 = current_data[82];
      const auto ayy_0 = current_data[83];
      const auto axz_0 = current_data[91];
      const auto ayz_0 = current_data[92];
      const auto azz_0 = current_data[100];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[0] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[1] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[2] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[3] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[4] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[5] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxxxy
      const auto a0_0 = current_data[1];
      const auto ax_0 = current_data[37];
      const auto ay_0 = current_data[38];
      const auto az_0 = current_data[46];
      const auto axx_0 = current_data[82];
      const auto axy_0 = current_data[83];
      const auto ayy_0 = current_data[84];
      const auto axz_0 = current_data[92];
      const auto ayz_0 = current_data[93];
      const auto azz_0 = current_data[101];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[6] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[7] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[8] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[9] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[10] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[11] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxxyy
      const auto a0_0 = current_data[2];
      const auto ax_0 = current_data[38];
      const auto ay_0 = current_data[39];
      const auto az_0 = current_data[47];
      const auto axx_0 = current_data[83];
      const auto axy_0 = current_data[84];
      const auto ayy_0 = current_data[85];
      const auto axz_0 = current_data[93];
      const auto ayz_0 = current_data[94];
      const auto azz_0 = current_data[102];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[12] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[13] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[14] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[15] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[16] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[17] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxyyy
      const auto a0_0 = current_data[3];
      const auto ax_0 = current_data[39];
      const auto ay_0 = current_data[40];
      const auto az_0 = current_data[48];
      const auto axx_0 = current_data[84];
      const auto axy_0 = current_data[85];
      const auto ayy_0 = current_data[86];
      const auto axz_0 = current_data[94];
      const auto ayz_0 = current_data[95];
      const auto azz_0 = current_data[103];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[18] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[19] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[20] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[21] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[22] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[23] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxyyyy
      const auto a0_0 = current_data[4];
      const auto ax_0 = current_data[40];
      const auto ay_0 = current_data[41];
      const auto az_0 = current_data[49];
      const auto axx_0 = current_data[85];
      const auto axy_0 = current_data[86];
      const auto ayy_0 = current_data[87];
      const auto axz_0 = current_data[95];
      const auto ayz_0 = current_data[96];
      const auto azz_0 = current_data[104];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[24] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[25] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[26] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[27] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[28] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[29] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyyyyy
      const auto a0_0 = current_data[5];
      const auto ax_0 = current_data[41];
      const auto ay_0 = current_data[42];
      const auto az_0 = current_data[50];
      const auto axx_0 = current_data[86];
      const auto axy_0 = current_data[87];
      const auto ayy_0 = current_data[88];
      const auto axz_0 = current_data[96];
      const auto ayz_0 = current_data[97];
      const auto azz_0 = current_data[105];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[30] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[31] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[32] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[33] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[34] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[35] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyyyyy
      const auto a0_0 = current_data[6];
      const auto ax_0 = current_data[42];
      const auto ay_0 = current_data[43];
      const auto az_0 = current_data[51];
      const auto axx_0 = current_data[87];
      const auto axy_0 = current_data[88];
      const auto ayy_0 = current_data[89];
      const auto axz_0 = current_data[97];
      const auto ayz_0 = current_data[98];
      const auto azz_0 = current_data[106];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[36] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[37] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[38] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[39] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[40] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[41] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyyyyy
      const auto a0_0 = current_data[7];
      const auto ax_0 = current_data[43];
      const auto ay_0 = current_data[44];
      const auto az_0 = current_data[52];
      const auto axx_0 = current_data[88];
      const auto axy_0 = current_data[89];
      const auto ayy_0 = current_data[90];
      const auto axz_0 = current_data[98];
      const auto ayz_0 = current_data[99];
      const auto azz_0 = current_data[107];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[42] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[43] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[44] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[45] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[46] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[47] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxxxz
      const auto a0_0 = current_data[8];
      const auto ax_0 = current_data[45];
      const auto ay_0 = current_data[46];
      const auto az_0 = current_data[53];
      const auto axx_0 = current_data[91];
      const auto axy_0 = current_data[92];
      const auto ayy_0 = current_data[93];
      const auto axz_0 = current_data[100];
      const auto ayz_0 = current_data[101];
      const auto azz_0 = current_data[108];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[48] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[49] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[50] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[51] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[52] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[53] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxxyz
      const auto a0_0 = current_data[9];
      const auto ax_0 = current_data[46];
      const auto ay_0 = current_data[47];
      const auto az_0 = current_data[54];
      const auto axx_0 = current_data[92];
      const auto axy_0 = current_data[93];
      const auto ayy_0 = current_data[94];
      const auto axz_0 = current_data[101];
      const auto ayz_0 = current_data[102];
      const auto azz_0 = current_data[109];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[54] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[55] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[56] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[57] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[58] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[59] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxyyz
      const auto a0_0 = current_data[10];
      const auto ax_0 = current_data[47];
      const auto ay_0 = current_data[48];
      const auto az_0 = current_data[55];
      const auto axx_0 = current_data[93];
      const auto axy_0 = current_data[94];
      const auto ayy_0 = current_data[95];
      const auto axz_0 = current_data[102];
      const auto ayz_0 = current_data[103];
      const auto azz_0 = current_data[110];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[60] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[61] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[62] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[63] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[64] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[65] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxyyyz
      const auto a0_0 = current_data[11];
      const auto ax_0 = current_data[48];
      const auto ay_0 = current_data[49];
      const auto az_0 = current_data[56];
      const auto axx_0 = current_data[94];
      const auto axy_0 = current_data[95];
      const auto ayy_0 = current_data[96];
      const auto axz_0 = current_data[103];
      const auto ayz_0 = current_data[104];
      const auto azz_0 = current_data[111];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[66] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[67] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[68] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[69] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[70] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[71] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyyyyz
      const auto a0_0 = current_data[12];
      const auto ax_0 = current_data[49];
      const auto ay_0 = current_data[50];
      const auto az_0 = current_data[57];
      const auto axx_0 = current_data[95];
      const auto axy_0 = current_data[96];
      const auto ayy_0 = current_data[97];
      const auto axz_0 = current_data[104];
      const auto ayz_0 = current_data[105];
      const auto azz_0 = current_data[112];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[72] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[73] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[74] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[75] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[76] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[77] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyyyyz
      const auto a0_0 = current_data[13];
      const auto ax_0 = current_data[50];
      const auto ay_0 = current_data[51];
      const auto az_0 = current_data[58];
      const auto axx_0 = current_data[96];
      const auto axy_0 = current_data[97];
      const auto ayy_0 = current_data[98];
      const auto axz_0 = current_data[105];
      const auto ayz_0 = current_data[106];
      const auto azz_0 = current_data[113];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[78] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[79] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[80] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[81] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[82] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[83] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyyyyz
      const auto a0_0 = current_data[14];
      const auto ax_0 = current_data[51];
      const auto ay_0 = current_data[52];
      const auto az_0 = current_data[59];
      const auto axx_0 = current_data[97];
      const auto axy_0 = current_data[98];
      const auto ayy_0 = current_data[99];
      const auto axz_0 = current_data[106];
      const auto ayz_0 = current_data[107];
      const auto azz_0 = current_data[114];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[84] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[85] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[86] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[87] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[88] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[89] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxxzz
      const auto a0_0 = current_data[15];
      const auto ax_0 = current_data[53];
      const auto ay_0 = current_data[54];
      const auto az_0 = current_data[60];
      const auto axx_0 = current_data[100];
      const auto axy_0 = current_data[101];
      const auto ayy_0 = current_data[102];
      const auto axz_0 = current_data[108];
      const auto ayz_0 = current_data[109];
      const auto azz_0 = current_data[115];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[90] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[91] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[92] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[93] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[94] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[95] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxyzz
      const auto a0_0 = current_data[16];
      const auto ax_0 = current_data[54];
      const auto ay_0 = current_data[55];
      const auto az_0 = current_data[61];
      const auto axx_0 = current_data[101];
      const auto axy_0 = current_data[102];
      const auto ayy_0 = current_data[103];
      const auto axz_0 = current_data[109];
      const auto ayz_0 = current_data[110];
      const auto azz_0 = current_data[116];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[96] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[97] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[98] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[99] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[100] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[101] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxyyzz
      const auto a0_0 = current_data[17];
      const auto ax_0 = current_data[55];
      const auto ay_0 = current_data[56];
      const auto az_0 = current_data[62];
      const auto axx_0 = current_data[102];
      const auto axy_0 = current_data[103];
      const auto ayy_0 = current_data[104];
      const auto axz_0 = current_data[110];
      const auto ayz_0 = current_data[111];
      const auto azz_0 = current_data[117];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[102] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[103] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[104] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[105] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[106] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[107] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyyyzz
      const auto a0_0 = current_data[18];
      const auto ax_0 = current_data[56];
      const auto ay_0 = current_data[57];
      const auto az_0 = current_data[63];
      const auto axx_0 = current_data[103];
      const auto axy_0 = current_data[104];
      const auto ayy_0 = current_data[105];
      const auto axz_0 = current_data[111];
      const auto ayz_0 = current_data[112];
      const auto azz_0 = current_data[118];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[108] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[109] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[110] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[111] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[112] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[113] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyyyzz
      const auto a0_0 = current_data[19];
      const auto ax_0 = current_data[57];
      const auto ay_0 = current_data[58];
      const auto az_0 = current_data[64];
      const auto axx_0 = current_data[104];
      const auto axy_0 = current_data[105];
      const auto ayy_0 = current_data[106];
      const auto axz_0 = current_data[112];
      const auto ayz_0 = current_data[113];
      const auto azz_0 = current_data[119];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[114] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[115] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[116] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[117] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[118] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[119] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyyyzz
      const auto a0_0 = current_data[20];
      const auto ax_0 = current_data[58];
      const auto ay_0 = current_data[59];
      const auto az_0 = current_data[65];
      const auto axx_0 = current_data[105];
      const auto axy_0 = current_data[106];
      const auto ayy_0 = current_data[107];
      const auto axz_0 = current_data[113];
      const auto ayz_0 = current_data[114];
      const auto azz_0 = current_data[120];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[120] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[121] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[122] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[123] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[124] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[125] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxzzz
      const auto a0_0 = current_data[21];
      const auto ax_0 = current_data[60];
      const auto ay_0 = current_data[61];
      const auto az_0 = current_data[66];
      const auto axx_0 = current_data[108];
      const auto axy_0 = current_data[109];
      const auto ayy_0 = current_data[110];
      const auto axz_0 = current_data[115];
      const auto ayz_0 = current_data[116];
      const auto azz_0 = current_data[121];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[126] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[127] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[128] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[129] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[130] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[131] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxyzzz
      const auto a0_0 = current_data[22];
      const auto ax_0 = current_data[61];
      const auto ay_0 = current_data[62];
      const auto az_0 = current_data[67];
      const auto axx_0 = current_data[109];
      const auto axy_0 = current_data[110];
      const auto ayy_0 = current_data[111];
      const auto axz_0 = current_data[116];
      const auto ayz_0 = current_data[117];
      const auto azz_0 = current_data[122];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[132] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[133] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[134] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[135] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[136] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[137] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyyzzz
      const auto a0_0 = current_data[23];
      const auto ax_0 = current_data[62];
      const auto ay_0 = current_data[63];
      const auto az_0 = current_data[68];
      const auto axx_0 = current_data[110];
      const auto axy_0 = current_data[111];
      const auto ayy_0 = current_data[112];
      const auto axz_0 = current_data[117];
      const auto ayz_0 = current_data[118];
      const auto azz_0 = current_data[123];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[138] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[139] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[140] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[141] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[142] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[143] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyyzzz
      const auto a0_0 = current_data[24];
      const auto ax_0 = current_data[63];
      const auto ay_0 = current_data[64];
      const auto az_0 = current_data[69];
      const auto axx_0 = current_data[111];
      const auto axy_0 = current_data[112];
      const auto ayy_0 = current_data[113];
      const auto axz_0 = current_data[118];
      const auto ayz_0 = current_data[119];
      const auto azz_0 = current_data[124];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[144] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[145] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[146] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[147] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[148] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[149] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyyzzz
      const auto a0_0 = current_data[25];
      const auto ax_0 = current_data[64];
      const auto ay_0 = current_data[65];
      const auto az_0 = current_data[70];
      const auto axx_0 = current_data[112];
      const auto axy_0 = current_data[113];
      const auto ayy_0 = current_data[114];
      const auto axz_0 = current_data[119];
      const auto ayz_0 = current_data[120];
      const auto azz_0 = current_data[125];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[150] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[151] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[152] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[153] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[154] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[155] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxzzzz
      const auto a0_0 = current_data[26];
      const auto ax_0 = current_data[66];
      const auto ay_0 = current_data[67];
      const auto az_0 = current_data[71];
      const auto axx_0 = current_data[115];
      const auto axy_0 = current_data[116];
      const auto ayy_0 = current_data[117];
      const auto axz_0 = current_data[121];
      const auto ayz_0 = current_data[122];
      const auto azz_0 = current_data[126];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[156] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[157] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[158] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[159] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[160] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[161] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyzzzz
      const auto a0_0 = current_data[27];
      const auto ax_0 = current_data[67];
      const auto ay_0 = current_data[68];
      const auto az_0 = current_data[72];
      const auto axx_0 = current_data[116];
      const auto axy_0 = current_data[117];
      const auto ayy_0 = current_data[118];
      const auto axz_0 = current_data[122];
      const auto ayz_0 = current_data[123];
      const auto azz_0 = current_data[127];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[162] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[163] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[164] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[165] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[166] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[167] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyzzzz
      const auto a0_0 = current_data[28];
      const auto ax_0 = current_data[68];
      const auto ay_0 = current_data[69];
      const auto az_0 = current_data[73];
      const auto axx_0 = current_data[117];
      const auto axy_0 = current_data[118];
      const auto ayy_0 = current_data[119];
      const auto axz_0 = current_data[123];
      const auto ayz_0 = current_data[124];
      const auto azz_0 = current_data[128];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[168] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[169] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[170] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[171] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[172] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[173] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyzzzz
      const auto a0_0 = current_data[29];
      const auto ax_0 = current_data[69];
      const auto ay_0 = current_data[70];
      const auto az_0 = current_data[74];
      const auto axx_0 = current_data[118];
      const auto axy_0 = current_data[119];
      const auto ayy_0 = current_data[120];
      const auto axz_0 = current_data[124];
      const auto ayz_0 = current_data[125];
      const auto azz_0 = current_data[129];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[174] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[175] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[176] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[177] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[178] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[179] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxzzzzz
      const auto a0_0 = current_data[30];
      const auto ax_0 = current_data[71];
      const auto ay_0 = current_data[72];
      const auto az_0 = current_data[75];
      const auto axx_0 = current_data[121];
      const auto axy_0 = current_data[122];
      const auto ayy_0 = current_data[123];
      const auto axz_0 = current_data[126];
      const auto ayz_0 = current_data[127];
      const auto azz_0 = current_data[130];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[180] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[181] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[182] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[183] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[184] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[185] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyzzzzz
      const auto a0_0 = current_data[31];
      const auto ax_0 = current_data[72];
      const auto ay_0 = current_data[73];
      const auto az_0 = current_data[76];
      const auto axx_0 = current_data[122];
      const auto axy_0 = current_data[123];
      const auto ayy_0 = current_data[124];
      const auto axz_0 = current_data[127];
      const auto ayz_0 = current_data[128];
      const auto azz_0 = current_data[131];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[186] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[187] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[188] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[189] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[190] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[191] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyzzzzz
      const auto a0_0 = current_data[32];
      const auto ax_0 = current_data[73];
      const auto ay_0 = current_data[74];
      const auto az_0 = current_data[77];
      const auto axx_0 = current_data[123];
      const auto axy_0 = current_data[124];
      const auto ayy_0 = current_data[125];
      const auto axz_0 = current_data[128];
      const auto ayz_0 = current_data[129];
      const auto azz_0 = current_data[132];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[192] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[193] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[194] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[195] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[196] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[197] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xzzzzzz
      const auto a0_0 = current_data[33];
      const auto ax_0 = current_data[75];
      const auto ay_0 = current_data[76];
      const auto az_0 = current_data[78];
      const auto axx_0 = current_data[126];
      const auto axy_0 = current_data[127];
      const auto ayy_0 = current_data[128];
      const auto axz_0 = current_data[130];
      const auto ayz_0 = current_data[131];
      const auto azz_0 = current_data[133];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[198] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[199] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[200] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[201] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[202] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[203] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yzzzzzz
      const auto a0_0 = current_data[34];
      const auto ax_0 = current_data[76];
      const auto ay_0 = current_data[77];
      const auto az_0 = current_data[79];
      const auto axx_0 = current_data[127];
      const auto axy_0 = current_data[128];
      const auto ayy_0 = current_data[129];
      const auto axz_0 = current_data[131];
      const auto ayz_0 = current_data[132];
      const auto azz_0 = current_data[134];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[204] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[205] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[206] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[207] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[208] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[209] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: zzzzzzz
      const auto a0_0 = current_data[35];
      const auto ax_0 = current_data[78];
      const auto ay_0 = current_data[79];
      const auto az_0 = current_data[80];
      const auto axx_0 = current_data[130];
      const auto axy_0 = current_data[131];
      const auto ayy_0 = current_data[132];
      const auto axz_0 = current_data[133];
      const auto ayz_0 = current_data[134];
      const auto azz_0 = current_data[135];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[210] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[211] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[212] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[213] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[214] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[215] = az_z + AB[2] * a0_z; // a0_zz

    }
  }
}


void CHRRList::perform_HRR_90_72(const int nloop, const complex<double>* data_start, const array<double,3>& AB, complex<double>* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 136];
    auto current_out = &data_out[c * 216];
   {
     //current index a: xxxxxxx
      const auto a0_0 = current_data[0];
      const auto ax_0 = current_data[36];
      const auto ay_0 = current_data[37];
      const auto az_0 = current_data[45];
      const auto axx_0 = current_data[81];
      const auto axy_0 = current_data[82];
      const auto ayy_0 = current_data[83];
      const auto axz_0 = current_data[91];
      const auto ayz_0 = current_data[92];
      const auto azz_0 = current_data[100];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[0] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[1] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[2] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[3] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[4] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[5] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxxxy
      const auto a0_0 = current_data[1];
      const auto ax_0 = current_data[37];
      const auto ay_0 = current_data[38];
      const auto az_0 = current_data[46];
      const auto axx_0 = current_data[82];
      const auto axy_0 = current_data[83];
      const auto ayy_0 = current_data[84];
      const auto axz_0 = current_data[92];
      const auto ayz_0 = current_data[93];
      const auto azz_0 = current_data[101];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[6] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[7] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[8] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[9] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[10] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[11] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxxyy
      const auto a0_0 = current_data[2];
      const auto ax_0 = current_data[38];
      const auto ay_0 = current_data[39];
      const auto az_0 = current_data[47];
      const auto axx_0 = current_data[83];
      const auto axy_0 = current_data[84];
      const auto ayy_0 = current_data[85];
      const auto axz_0 = current_data[93];
      const auto ayz_0 = current_data[94];
      const auto azz_0 = current_data[102];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[12] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[13] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[14] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[15] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[16] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[17] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxyyy
      const auto a0_0 = current_data[3];
      const auto ax_0 = current_data[39];
      const auto ay_0 = current_data[40];
      const auto az_0 = current_data[48];
      const auto axx_0 = current_data[84];
      const auto axy_0 = current_data[85];
      const auto ayy_0 = current_data[86];
      const auto axz_0 = current_data[94];
      const auto ayz_0 = current_data[95];
      const auto azz_0 = current_data[103];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[18] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[19] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[20] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[21] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[22] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[23] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxyyyy
      const auto a0_0 = current_data[4];
      const auto ax_0 = current_data[40];
      const auto ay_0 = current_data[41];
      const auto az_0 = current_data[49];
      const auto axx_0 = current_data[85];
      const auto axy_0 = current_data[86];
      const auto ayy_0 = current_data[87];
      const auto axz_0 = current_data[95];
      const auto ayz_0 = current_data[96];
      const auto azz_0 = current_data[104];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[24] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[25] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[26] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[27] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[28] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[29] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyyyyy
      const auto a0_0 = current_data[5];
      const auto ax_0 = current_data[41];
      const auto ay_0 = current_data[42];
      const auto az_0 = current_data[50];
      const auto axx_0 = current_data[86];
      const auto axy_0 = current_data[87];
      const auto ayy_0 = current_data[88];
      const auto axz_0 = current_data[96];
      const auto ayz_0 = current_data[97];
      const auto azz_0 = current_data[105];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[30] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[31] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[32] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[33] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[34] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[35] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyyyyy
      const auto a0_0 = current_data[6];
      const auto ax_0 = current_data[42];
      const auto ay_0 = current_data[43];
      const auto az_0 = current_data[51];
      const auto axx_0 = current_data[87];
      const auto axy_0 = current_data[88];
      const auto ayy_0 = current_data[89];
      const auto axz_0 = current_data[97];
      const auto ayz_0 = current_data[98];
      const auto azz_0 = current_data[106];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[36] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[37] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[38] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[39] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[40] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[41] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyyyyy
      const auto a0_0 = current_data[7];
      const auto ax_0 = current_data[43];
      const auto ay_0 = current_data[44];
      const auto az_0 = current_data[52];
      const auto axx_0 = current_data[88];
      const auto axy_0 = current_data[89];
      const auto ayy_0 = current_data[90];
      const auto axz_0 = current_data[98];
      const auto ayz_0 = current_data[99];
      const auto azz_0 = current_data[107];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[42] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[43] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[44] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[45] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[46] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[47] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxxxz
      const auto a0_0 = current_data[8];
      const auto ax_0 = current_data[45];
      const auto ay_0 = current_data[46];
      const auto az_0 = current_data[53];
      const auto axx_0 = current_data[91];
      const auto axy_0 = current_data[92];
      const auto ayy_0 = current_data[93];
      const auto axz_0 = current_data[100];
      const auto ayz_0 = current_data[101];
      const auto azz_0 = current_data[108];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[48] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[49] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[50] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[51] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[52] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[53] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxxyz
      const auto a0_0 = current_data[9];
      const auto ax_0 = current_data[46];
      const auto ay_0 = current_data[47];
      const auto az_0 = current_data[54];
      const auto axx_0 = current_data[92];
      const auto axy_0 = current_data[93];
      const auto ayy_0 = current_data[94];
      const auto axz_0 = current_data[101];
      const auto ayz_0 = current_data[102];
      const auto azz_0 = current_data[109];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[54] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[55] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[56] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[57] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[58] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[59] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxyyz
      const auto a0_0 = current_data[10];
      const auto ax_0 = current_data[47];
      const auto ay_0 = current_data[48];
      const auto az_0 = current_data[55];
      const auto axx_0 = current_data[93];
      const auto axy_0 = current_data[94];
      const auto ayy_0 = current_data[95];
      const auto axz_0 = current_data[102];
      const auto ayz_0 = current_data[103];
      const auto azz_0 = current_data[110];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[60] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[61] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[62] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[63] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[64] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[65] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxyyyz
      const auto a0_0 = current_data[11];
      const auto ax_0 = current_data[48];
      const auto ay_0 = current_data[49];
      const auto az_0 = current_data[56];
      const auto axx_0 = current_data[94];
      const auto axy_0 = current_data[95];
      const auto ayy_0 = current_data[96];
      const auto axz_0 = current_data[103];
      const auto ayz_0 = current_data[104];
      const auto azz_0 = current_data[111];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[66] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[67] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[68] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[69] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[70] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[71] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyyyyz
      const auto a0_0 = current_data[12];
      const auto ax_0 = current_data[49];
      const auto ay_0 = current_data[50];
      const auto az_0 = current_data[57];
      const auto axx_0 = current_data[95];
      const auto axy_0 = current_data[96];
      const auto ayy_0 = current_data[97];
      const auto axz_0 = current_data[104];
      const auto ayz_0 = current_data[105];
      const auto azz_0 = current_data[112];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[72] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[73] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[74] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[75] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[76] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[77] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyyyyz
      const auto a0_0 = current_data[13];
      const auto ax_0 = current_data[50];
      const auto ay_0 = current_data[51];
      const auto az_0 = current_data[58];
      const auto axx_0 = current_data[96];
      const auto axy_0 = current_data[97];
      const auto ayy_0 = current_data[98];
      const auto axz_0 = current_data[105];
      const auto ayz_0 = current_data[106];
      const auto azz_0 = current_data[113];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[78] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[79] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[80] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[81] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[82] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[83] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyyyyz
      const auto a0_0 = current_data[14];
      const auto ax_0 = current_data[51];
      const auto ay_0 = current_data[52];
      const auto az_0 = current_data[59];
      const auto axx_0 = current_data[97];
      const auto axy_0 = current_data[98];
      const auto ayy_0 = current_data[99];
      const auto axz_0 = current_data[106];
      const auto ayz_0 = current_data[107];
      const auto azz_0 = current_data[114];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[84] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[85] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[86] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[87] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[88] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[89] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxxzz
      const auto a0_0 = current_data[15];
      const auto ax_0 = current_data[53];
      const auto ay_0 = current_data[54];
      const auto az_0 = current_data[60];
      const auto axx_0 = current_data[100];
      const auto axy_0 = current_data[101];
      const auto ayy_0 = current_data[102];
      const auto axz_0 = current_data[108];
      const auto ayz_0 = current_data[109];
      const auto azz_0 = current_data[115];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[90] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[91] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[92] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[93] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[94] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[95] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxyzz
      const auto a0_0 = current_data[16];
      const auto ax_0 = current_data[54];
      const auto ay_0 = current_data[55];
      const auto az_0 = current_data[61];
      const auto axx_0 = current_data[101];
      const auto axy_0 = current_data[102];
      const auto ayy_0 = current_data[103];
      const auto axz_0 = current_data[109];
      const auto ayz_0 = current_data[110];
      const auto azz_0 = current_data[116];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[96] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[97] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[98] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[99] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[100] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[101] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxyyzz
      const auto a0_0 = current_data[17];
      const auto ax_0 = current_data[55];
      const auto ay_0 = current_data[56];
      const auto az_0 = current_data[62];
      const auto axx_0 = current_data[102];
      const auto axy_0 = current_data[103];
      const auto ayy_0 = current_data[104];
      const auto axz_0 = current_data[110];
      const auto ayz_0 = current_data[111];
      const auto azz_0 = current_data[117];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[102] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[103] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[104] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[105] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[106] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[107] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyyyzz
      const auto a0_0 = current_data[18];
      const auto ax_0 = current_data[56];
      const auto ay_0 = current_data[57];
      const auto az_0 = current_data[63];
      const auto axx_0 = current_data[103];
      const auto axy_0 = current_data[104];
      const auto ayy_0 = current_data[105];
      const auto axz_0 = current_data[111];
      const auto ayz_0 = current_data[112];
      const auto azz_0 = current_data[118];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[108] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[109] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[110] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[111] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[112] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[113] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyyyzz
      const auto a0_0 = current_data[19];
      const auto ax_0 = current_data[57];
      const auto ay_0 = current_data[58];
      const auto az_0 = current_data[64];
      const auto axx_0 = current_data[104];
      const auto axy_0 = current_data[105];
      const auto ayy_0 = current_data[106];
      const auto axz_0 = current_data[112];
      const auto ayz_0 = current_data[113];
      const auto azz_0 = current_data[119];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[114] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[115] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[116] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[117] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[118] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[119] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyyyzz
      const auto a0_0 = current_data[20];
      const auto ax_0 = current_data[58];
      const auto ay_0 = current_data[59];
      const auto az_0 = current_data[65];
      const auto axx_0 = current_data[105];
      const auto axy_0 = current_data[106];
      const auto ayy_0 = current_data[107];
      const auto axz_0 = current_data[113];
      const auto ayz_0 = current_data[114];
      const auto azz_0 = current_data[120];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[120] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[121] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[122] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[123] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[124] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[125] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxzzz
      const auto a0_0 = current_data[21];
      const auto ax_0 = current_data[60];
      const auto ay_0 = current_data[61];
      const auto az_0 = current_data[66];
      const auto axx_0 = current_data[108];
      const auto axy_0 = current_data[109];
      const auto ayy_0 = current_data[110];
      const auto axz_0 = current_data[115];
      const auto ayz_0 = current_data[116];
      const auto azz_0 = current_data[121];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[126] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[127] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[128] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[129] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[130] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[131] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxyzzz
      const auto a0_0 = current_data[22];
      const auto ax_0 = current_data[61];
      const auto ay_0 = current_data[62];
      const auto az_0 = current_data[67];
      const auto axx_0 = current_data[109];
      const auto axy_0 = current_data[110];
      const auto ayy_0 = current_data[111];
      const auto axz_0 = current_data[116];
      const auto ayz_0 = current_data[117];
      const auto azz_0 = current_data[122];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[132] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[133] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[134] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[135] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[136] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[137] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyyzzz
      const auto a0_0 = current_data[23];
      const auto ax_0 = current_data[62];
      const auto ay_0 = current_data[63];
      const auto az_0 = current_data[68];
      const auto axx_0 = current_data[110];
      const auto axy_0 = current_data[111];
      const auto ayy_0 = current_data[112];
      const auto axz_0 = current_data[117];
      const auto ayz_0 = current_data[118];
      const auto azz_0 = current_data[123];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[138] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[139] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[140] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[141] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[142] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[143] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyyzzz
      const auto a0_0 = current_data[24];
      const auto ax_0 = current_data[63];
      const auto ay_0 = current_data[64];
      const auto az_0 = current_data[69];
      const auto axx_0 = current_data[111];
      const auto axy_0 = current_data[112];
      const auto ayy_0 = current_data[113];
      const auto axz_0 = current_data[118];
      const auto ayz_0 = current_data[119];
      const auto azz_0 = current_data[124];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[144] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[145] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[146] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[147] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[148] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[149] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyyzzz
      const auto a0_0 = current_data[25];
      const auto ax_0 = current_data[64];
      const auto ay_0 = current_data[65];
      const auto az_0 = current_data[70];
      const auto axx_0 = current_data[112];
      const auto axy_0 = current_data[113];
      const auto ayy_0 = current_data[114];
      const auto axz_0 = current_data[119];
      const auto ayz_0 = current_data[120];
      const auto azz_0 = current_data[125];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[150] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[151] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[152] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[153] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[154] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[155] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxzzzz
      const auto a0_0 = current_data[26];
      const auto ax_0 = current_data[66];
      const auto ay_0 = current_data[67];
      const auto az_0 = current_data[71];
      const auto axx_0 = current_data[115];
      const auto axy_0 = current_data[116];
      const auto ayy_0 = current_data[117];
      const auto axz_0 = current_data[121];
      const auto ayz_0 = current_data[122];
      const auto azz_0 = current_data[126];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[156] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[157] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[158] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[159] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[160] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[161] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyzzzz
      const auto a0_0 = current_data[27];
      const auto ax_0 = current_data[67];
      const auto ay_0 = current_data[68];
      const auto az_0 = current_data[72];
      const auto axx_0 = current_data[116];
      const auto axy_0 = current_data[117];
      const auto ayy_0 = current_data[118];
      const auto axz_0 = current_data[122];
      const auto ayz_0 = current_data[123];
      const auto azz_0 = current_data[127];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[162] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[163] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[164] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[165] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[166] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[167] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyzzzz
      const auto a0_0 = current_data[28];
      const auto ax_0 = current_data[68];
      const auto ay_0 = current_data[69];
      const auto az_0 = current_data[73];
      const auto axx_0 = current_data[117];
      const auto axy_0 = current_data[118];
      const auto ayy_0 = current_data[119];
      const auto axz_0 = current_data[123];
      const auto ayz_0 = current_data[124];
      const auto azz_0 = current_data[128];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[168] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[169] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[170] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[171] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[172] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[173] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyzzzz
      const auto a0_0 = current_data[29];
      const auto ax_0 = current_data[69];
      const auto ay_0 = current_data[70];
      const auto az_0 = current_data[74];
      const auto axx_0 = current_data[118];
      const auto axy_0 = current_data[119];
      const auto ayy_0 = current_data[120];
      const auto axz_0 = current_data[124];
      const auto ayz_0 = current_data[125];
      const auto azz_0 = current_data[129];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[174] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[175] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[176] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[177] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[178] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[179] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxzzzzz
      const auto a0_0 = current_data[30];
      const auto ax_0 = current_data[71];
      const auto ay_0 = current_data[72];
      const auto az_0 = current_data[75];
      const auto axx_0 = current_data[121];
      const auto axy_0 = current_data[122];
      const auto ayy_0 = current_data[123];
      const auto axz_0 = current_data[126];
      const auto ayz_0 = current_data[127];
      const auto azz_0 = current_data[130];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[180] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[181] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[182] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[183] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[184] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[185] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyzzzzz
      const auto a0_0 = current_data[31];
      const auto ax_0 = current_data[72];
      const auto ay_0 = current_data[73];
      const auto az_0 = current_data[76];
      const auto axx_0 = current_data[122];
      const auto axy_0 = current_data[123];
      const auto ayy_0 = current_data[124];
      const auto axz_0 = current_data[127];
      const auto ayz_0 = current_data[128];
      const auto azz_0 = current_data[131];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[186] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[187] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[188] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[189] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[190] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[191] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyzzzzz
      const auto a0_0 = current_data[32];
      const auto ax_0 = current_data[73];
      const auto ay_0 = current_data[74];
      const auto az_0 = current_data[77];
      const auto axx_0 = current_data[123];
      const auto axy_0 = current_data[124];
      const auto ayy_0 = current_data[125];
      const auto axz_0 = current_data[128];
      const auto ayz_0 = current_data[129];
      const auto azz_0 = current_data[132];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[192] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[193] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[194] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[195] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[196] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[197] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xzzzzzz
      const auto a0_0 = current_data[33];
      const auto ax_0 = current_data[75];
      const auto ay_0 = current_data[76];
      const auto az_0 = current_data[78];
      const auto axx_0 = current_data[126];
      const auto axy_0 = current_data[127];
      const auto ayy_0 = current_data[128];
      const auto axz_0 = current_data[130];
      const auto ayz_0 = current_data[131];
      const auto azz_0 = current_data[133];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[198] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[199] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[200] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[201] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[202] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[203] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yzzzzzz
      const auto a0_0 = current_data[34];
      const auto ax_0 = current_data[76];
      const auto ay_0 = current_data[77];
      const auto az_0 = current_data[79];
      const auto axx_0 = current_data[127];
      const auto axy_0 = current_data[128];
      const auto ayy_0 = current_data[129];
      const auto axz_0 = current_data[131];
      const auto ayz_0 = current_data[132];
      const auto azz_0 = current_data[134];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[204] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[205] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[206] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[207] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[208] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[209] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: zzzzzzz
      const auto a0_0 = current_data[35];
      const auto ax_0 = current_data[78];
      const auto ay_0 = current_data[79];
      const auto az_0 = current_data[80];
      const auto axx_0 = current_data[130];
      const auto axy_0 = current_data[131];
      const auto ayy_0 = current_data[132];
      const auto axz_0 = current_data[133];
      const auto ayz_0 = current_data[134];
      const auto azz_0 = current_data[135];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      current_out[210] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[211] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[212] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[213] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[214] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[215] = az_z + AB[2] * a0_z; // a0_zz

    }
  }
}
#endif

