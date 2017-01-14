//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _hrr_a0_73.cc
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

void HRRList::perform_HRR_a0_73(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 202];
    auto current_out = &data_out[c * 360];
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
      const auto axxx_0 = current_data[136];
      const auto axxy_0 = current_data[137];
      const auto axyy_0 = current_data[138];
      const auto ayyy_0 = current_data[139];
      const auto axxz_0 = current_data[147];
      const auto axyz_0 = current_data[148];
      const auto ayyz_0 = current_data[149];
      const auto axzz_0 = current_data[157];
      const auto ayzz_0 = current_data[158];
      const auto azzz_0 = current_data[166];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[0] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[1] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[2] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[3] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[4] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[5] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[6] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[7] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[8] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[9] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[137];
      const auto axxy_0 = current_data[138];
      const auto axyy_0 = current_data[139];
      const auto ayyy_0 = current_data[140];
      const auto axxz_0 = current_data[148];
      const auto axyz_0 = current_data[149];
      const auto ayyz_0 = current_data[150];
      const auto axzz_0 = current_data[158];
      const auto ayzz_0 = current_data[159];
      const auto azzz_0 = current_data[167];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[10] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[11] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[12] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[13] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[14] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[15] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[16] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[17] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[18] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[19] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[138];
      const auto axxy_0 = current_data[139];
      const auto axyy_0 = current_data[140];
      const auto ayyy_0 = current_data[141];
      const auto axxz_0 = current_data[149];
      const auto axyz_0 = current_data[150];
      const auto ayyz_0 = current_data[151];
      const auto axzz_0 = current_data[159];
      const auto ayzz_0 = current_data[160];
      const auto azzz_0 = current_data[168];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[20] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[21] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[22] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[23] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[24] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[25] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[26] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[27] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[28] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[29] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[139];
      const auto axxy_0 = current_data[140];
      const auto axyy_0 = current_data[141];
      const auto ayyy_0 = current_data[142];
      const auto axxz_0 = current_data[150];
      const auto axyz_0 = current_data[151];
      const auto ayyz_0 = current_data[152];
      const auto axzz_0 = current_data[160];
      const auto ayzz_0 = current_data[161];
      const auto azzz_0 = current_data[169];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[30] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[31] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[32] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[33] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[34] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[35] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[36] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[37] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[38] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[39] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[140];
      const auto axxy_0 = current_data[141];
      const auto axyy_0 = current_data[142];
      const auto ayyy_0 = current_data[143];
      const auto axxz_0 = current_data[151];
      const auto axyz_0 = current_data[152];
      const auto ayyz_0 = current_data[153];
      const auto axzz_0 = current_data[161];
      const auto ayzz_0 = current_data[162];
      const auto azzz_0 = current_data[170];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[40] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[41] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[42] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[43] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[44] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[45] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[46] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[47] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[48] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[49] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[141];
      const auto axxy_0 = current_data[142];
      const auto axyy_0 = current_data[143];
      const auto ayyy_0 = current_data[144];
      const auto axxz_0 = current_data[152];
      const auto axyz_0 = current_data[153];
      const auto ayyz_0 = current_data[154];
      const auto axzz_0 = current_data[162];
      const auto ayzz_0 = current_data[163];
      const auto azzz_0 = current_data[171];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[50] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[51] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[52] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[53] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[54] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[55] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[56] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[57] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[58] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[59] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[142];
      const auto axxy_0 = current_data[143];
      const auto axyy_0 = current_data[144];
      const auto ayyy_0 = current_data[145];
      const auto axxz_0 = current_data[153];
      const auto axyz_0 = current_data[154];
      const auto ayyz_0 = current_data[155];
      const auto axzz_0 = current_data[163];
      const auto ayzz_0 = current_data[164];
      const auto azzz_0 = current_data[172];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[60] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[61] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[62] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[63] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[64] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[65] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[66] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[67] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[68] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[69] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[143];
      const auto axxy_0 = current_data[144];
      const auto axyy_0 = current_data[145];
      const auto ayyy_0 = current_data[146];
      const auto axxz_0 = current_data[154];
      const auto axyz_0 = current_data[155];
      const auto ayyz_0 = current_data[156];
      const auto axzz_0 = current_data[164];
      const auto ayzz_0 = current_data[165];
      const auto azzz_0 = current_data[173];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[70] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[71] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[72] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[73] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[74] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[75] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[76] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[77] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[78] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[79] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[147];
      const auto axxy_0 = current_data[148];
      const auto axyy_0 = current_data[149];
      const auto ayyy_0 = current_data[150];
      const auto axxz_0 = current_data[157];
      const auto axyz_0 = current_data[158];
      const auto ayyz_0 = current_data[159];
      const auto axzz_0 = current_data[166];
      const auto ayzz_0 = current_data[167];
      const auto azzz_0 = current_data[174];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[80] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[81] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[82] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[83] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[84] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[85] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[86] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[87] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[88] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[89] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[148];
      const auto axxy_0 = current_data[149];
      const auto axyy_0 = current_data[150];
      const auto ayyy_0 = current_data[151];
      const auto axxz_0 = current_data[158];
      const auto axyz_0 = current_data[159];
      const auto ayyz_0 = current_data[160];
      const auto axzz_0 = current_data[167];
      const auto ayzz_0 = current_data[168];
      const auto azzz_0 = current_data[175];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[90] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[91] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[92] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[93] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[94] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[95] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[96] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[97] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[98] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[99] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[149];
      const auto axxy_0 = current_data[150];
      const auto axyy_0 = current_data[151];
      const auto ayyy_0 = current_data[152];
      const auto axxz_0 = current_data[159];
      const auto axyz_0 = current_data[160];
      const auto ayyz_0 = current_data[161];
      const auto axzz_0 = current_data[168];
      const auto ayzz_0 = current_data[169];
      const auto azzz_0 = current_data[176];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[100] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[101] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[102] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[103] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[104] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[105] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[106] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[107] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[108] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[109] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[150];
      const auto axxy_0 = current_data[151];
      const auto axyy_0 = current_data[152];
      const auto ayyy_0 = current_data[153];
      const auto axxz_0 = current_data[160];
      const auto axyz_0 = current_data[161];
      const auto ayyz_0 = current_data[162];
      const auto axzz_0 = current_data[169];
      const auto ayzz_0 = current_data[170];
      const auto azzz_0 = current_data[177];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[110] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[111] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[112] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[113] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[114] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[115] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[116] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[117] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[118] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[119] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[151];
      const auto axxy_0 = current_data[152];
      const auto axyy_0 = current_data[153];
      const auto ayyy_0 = current_data[154];
      const auto axxz_0 = current_data[161];
      const auto axyz_0 = current_data[162];
      const auto ayyz_0 = current_data[163];
      const auto axzz_0 = current_data[170];
      const auto ayzz_0 = current_data[171];
      const auto azzz_0 = current_data[178];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[120] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[121] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[122] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[123] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[124] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[125] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[126] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[127] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[128] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[129] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[152];
      const auto axxy_0 = current_data[153];
      const auto axyy_0 = current_data[154];
      const auto ayyy_0 = current_data[155];
      const auto axxz_0 = current_data[162];
      const auto axyz_0 = current_data[163];
      const auto ayyz_0 = current_data[164];
      const auto axzz_0 = current_data[171];
      const auto ayzz_0 = current_data[172];
      const auto azzz_0 = current_data[179];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[130] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[131] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[132] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[133] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[134] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[135] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[136] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[137] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[138] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[139] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[153];
      const auto axxy_0 = current_data[154];
      const auto axyy_0 = current_data[155];
      const auto ayyy_0 = current_data[156];
      const auto axxz_0 = current_data[163];
      const auto axyz_0 = current_data[164];
      const auto ayyz_0 = current_data[165];
      const auto axzz_0 = current_data[172];
      const auto ayzz_0 = current_data[173];
      const auto azzz_0 = current_data[180];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[140] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[141] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[142] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[143] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[144] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[145] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[146] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[147] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[148] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[149] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[157];
      const auto axxy_0 = current_data[158];
      const auto axyy_0 = current_data[159];
      const auto ayyy_0 = current_data[160];
      const auto axxz_0 = current_data[166];
      const auto axyz_0 = current_data[167];
      const auto ayyz_0 = current_data[168];
      const auto axzz_0 = current_data[174];
      const auto ayzz_0 = current_data[175];
      const auto azzz_0 = current_data[181];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[150] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[151] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[152] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[153] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[154] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[155] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[156] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[157] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[158] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[159] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[158];
      const auto axxy_0 = current_data[159];
      const auto axyy_0 = current_data[160];
      const auto ayyy_0 = current_data[161];
      const auto axxz_0 = current_data[167];
      const auto axyz_0 = current_data[168];
      const auto ayyz_0 = current_data[169];
      const auto axzz_0 = current_data[175];
      const auto ayzz_0 = current_data[176];
      const auto azzz_0 = current_data[182];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[160] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[161] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[162] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[163] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[164] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[165] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[166] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[167] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[168] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[169] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[159];
      const auto axxy_0 = current_data[160];
      const auto axyy_0 = current_data[161];
      const auto ayyy_0 = current_data[162];
      const auto axxz_0 = current_data[168];
      const auto axyz_0 = current_data[169];
      const auto ayyz_0 = current_data[170];
      const auto axzz_0 = current_data[176];
      const auto ayzz_0 = current_data[177];
      const auto azzz_0 = current_data[183];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[170] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[171] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[172] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[173] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[174] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[175] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[176] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[177] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[178] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[179] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[160];
      const auto axxy_0 = current_data[161];
      const auto axyy_0 = current_data[162];
      const auto ayyy_0 = current_data[163];
      const auto axxz_0 = current_data[169];
      const auto axyz_0 = current_data[170];
      const auto ayyz_0 = current_data[171];
      const auto axzz_0 = current_data[177];
      const auto ayzz_0 = current_data[178];
      const auto azzz_0 = current_data[184];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[180] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[181] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[182] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[183] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[184] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[185] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[186] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[187] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[188] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[189] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[161];
      const auto axxy_0 = current_data[162];
      const auto axyy_0 = current_data[163];
      const auto ayyy_0 = current_data[164];
      const auto axxz_0 = current_data[170];
      const auto axyz_0 = current_data[171];
      const auto ayyz_0 = current_data[172];
      const auto axzz_0 = current_data[178];
      const auto ayzz_0 = current_data[179];
      const auto azzz_0 = current_data[185];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[190] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[191] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[192] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[193] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[194] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[195] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[196] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[197] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[198] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[199] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[162];
      const auto axxy_0 = current_data[163];
      const auto axyy_0 = current_data[164];
      const auto ayyy_0 = current_data[165];
      const auto axxz_0 = current_data[171];
      const auto axyz_0 = current_data[172];
      const auto ayyz_0 = current_data[173];
      const auto axzz_0 = current_data[179];
      const auto ayzz_0 = current_data[180];
      const auto azzz_0 = current_data[186];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[200] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[201] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[202] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[203] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[204] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[205] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[206] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[207] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[208] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[209] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[166];
      const auto axxy_0 = current_data[167];
      const auto axyy_0 = current_data[168];
      const auto ayyy_0 = current_data[169];
      const auto axxz_0 = current_data[174];
      const auto axyz_0 = current_data[175];
      const auto ayyz_0 = current_data[176];
      const auto axzz_0 = current_data[181];
      const auto ayzz_0 = current_data[182];
      const auto azzz_0 = current_data[187];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[210] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[211] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[212] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[213] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[214] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[215] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[216] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[217] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[218] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[219] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[167];
      const auto axxy_0 = current_data[168];
      const auto axyy_0 = current_data[169];
      const auto ayyy_0 = current_data[170];
      const auto axxz_0 = current_data[175];
      const auto axyz_0 = current_data[176];
      const auto ayyz_0 = current_data[177];
      const auto axzz_0 = current_data[182];
      const auto ayzz_0 = current_data[183];
      const auto azzz_0 = current_data[188];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[220] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[221] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[222] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[223] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[224] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[225] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[226] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[227] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[228] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[229] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[168];
      const auto axxy_0 = current_data[169];
      const auto axyy_0 = current_data[170];
      const auto ayyy_0 = current_data[171];
      const auto axxz_0 = current_data[176];
      const auto axyz_0 = current_data[177];
      const auto ayyz_0 = current_data[178];
      const auto axzz_0 = current_data[183];
      const auto ayzz_0 = current_data[184];
      const auto azzz_0 = current_data[189];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[230] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[231] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[232] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[233] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[234] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[235] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[236] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[237] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[238] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[239] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[169];
      const auto axxy_0 = current_data[170];
      const auto axyy_0 = current_data[171];
      const auto ayyy_0 = current_data[172];
      const auto axxz_0 = current_data[177];
      const auto axyz_0 = current_data[178];
      const auto ayyz_0 = current_data[179];
      const auto axzz_0 = current_data[184];
      const auto ayzz_0 = current_data[185];
      const auto azzz_0 = current_data[190];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[240] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[241] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[242] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[243] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[244] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[245] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[246] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[247] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[248] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[249] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[170];
      const auto axxy_0 = current_data[171];
      const auto axyy_0 = current_data[172];
      const auto ayyy_0 = current_data[173];
      const auto axxz_0 = current_data[178];
      const auto axyz_0 = current_data[179];
      const auto ayyz_0 = current_data[180];
      const auto axzz_0 = current_data[185];
      const auto ayzz_0 = current_data[186];
      const auto azzz_0 = current_data[191];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[250] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[251] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[252] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[253] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[254] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[255] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[256] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[257] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[258] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[259] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[174];
      const auto axxy_0 = current_data[175];
      const auto axyy_0 = current_data[176];
      const auto ayyy_0 = current_data[177];
      const auto axxz_0 = current_data[181];
      const auto axyz_0 = current_data[182];
      const auto ayyz_0 = current_data[183];
      const auto axzz_0 = current_data[187];
      const auto ayzz_0 = current_data[188];
      const auto azzz_0 = current_data[192];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[260] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[261] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[262] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[263] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[264] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[265] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[266] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[267] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[268] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[269] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[175];
      const auto axxy_0 = current_data[176];
      const auto axyy_0 = current_data[177];
      const auto ayyy_0 = current_data[178];
      const auto axxz_0 = current_data[182];
      const auto axyz_0 = current_data[183];
      const auto ayyz_0 = current_data[184];
      const auto axzz_0 = current_data[188];
      const auto ayzz_0 = current_data[189];
      const auto azzz_0 = current_data[193];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[270] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[271] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[272] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[273] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[274] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[275] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[276] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[277] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[278] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[279] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[176];
      const auto axxy_0 = current_data[177];
      const auto axyy_0 = current_data[178];
      const auto ayyy_0 = current_data[179];
      const auto axxz_0 = current_data[183];
      const auto axyz_0 = current_data[184];
      const auto ayyz_0 = current_data[185];
      const auto axzz_0 = current_data[189];
      const auto ayzz_0 = current_data[190];
      const auto azzz_0 = current_data[194];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[280] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[281] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[282] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[283] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[284] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[285] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[286] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[287] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[288] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[289] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[177];
      const auto axxy_0 = current_data[178];
      const auto axyy_0 = current_data[179];
      const auto ayyy_0 = current_data[180];
      const auto axxz_0 = current_data[184];
      const auto axyz_0 = current_data[185];
      const auto ayyz_0 = current_data[186];
      const auto axzz_0 = current_data[190];
      const auto ayzz_0 = current_data[191];
      const auto azzz_0 = current_data[195];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[290] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[291] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[292] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[293] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[294] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[295] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[296] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[297] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[298] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[299] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[181];
      const auto axxy_0 = current_data[182];
      const auto axyy_0 = current_data[183];
      const auto ayyy_0 = current_data[184];
      const auto axxz_0 = current_data[187];
      const auto axyz_0 = current_data[188];
      const auto ayyz_0 = current_data[189];
      const auto axzz_0 = current_data[192];
      const auto ayzz_0 = current_data[193];
      const auto azzz_0 = current_data[196];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[300] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[301] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[302] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[303] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[304] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[305] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[306] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[307] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[308] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[309] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[182];
      const auto axxy_0 = current_data[183];
      const auto axyy_0 = current_data[184];
      const auto ayyy_0 = current_data[185];
      const auto axxz_0 = current_data[188];
      const auto axyz_0 = current_data[189];
      const auto ayyz_0 = current_data[190];
      const auto axzz_0 = current_data[193];
      const auto ayzz_0 = current_data[194];
      const auto azzz_0 = current_data[197];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[310] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[311] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[312] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[313] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[314] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[315] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[316] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[317] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[318] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[319] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[183];
      const auto axxy_0 = current_data[184];
      const auto axyy_0 = current_data[185];
      const auto ayyy_0 = current_data[186];
      const auto axxz_0 = current_data[189];
      const auto axyz_0 = current_data[190];
      const auto ayyz_0 = current_data[191];
      const auto axzz_0 = current_data[194];
      const auto ayzz_0 = current_data[195];
      const auto azzz_0 = current_data[198];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[320] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[321] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[322] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[323] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[324] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[325] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[326] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[327] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[328] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[329] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[187];
      const auto axxy_0 = current_data[188];
      const auto axyy_0 = current_data[189];
      const auto ayyy_0 = current_data[190];
      const auto axxz_0 = current_data[192];
      const auto axyz_0 = current_data[193];
      const auto ayyz_0 = current_data[194];
      const auto axzz_0 = current_data[196];
      const auto ayzz_0 = current_data[197];
      const auto azzz_0 = current_data[199];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[330] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[331] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[332] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[333] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[334] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[335] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[336] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[337] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[338] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[339] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[188];
      const auto axxy_0 = current_data[189];
      const auto axyy_0 = current_data[190];
      const auto ayyy_0 = current_data[191];
      const auto axxz_0 = current_data[193];
      const auto axyz_0 = current_data[194];
      const auto ayyz_0 = current_data[195];
      const auto axzz_0 = current_data[197];
      const auto ayzz_0 = current_data[198];
      const auto azzz_0 = current_data[200];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[340] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[341] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[342] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[343] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[344] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[345] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[346] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[347] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[348] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[349] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[192];
      const auto axxy_0 = current_data[193];
      const auto axyy_0 = current_data[194];
      const auto ayyy_0 = current_data[195];
      const auto axxz_0 = current_data[196];
      const auto axyz_0 = current_data[197];
      const auto ayyz_0 = current_data[198];
      const auto axzz_0 = current_data[199];
      const auto ayzz_0 = current_data[200];
      const auto azzz_0 = current_data[201];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[350] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[351] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[352] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[353] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[354] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[355] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[356] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[357] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[358] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[359] = az_zz + AB[2] * a0_zz; // a0_zzz

    }
  }
}


void CHRRList::perform_HRR_a0_73(const int nloop, const complex<double>* data_start, const array<double,3>& AB, complex<double>* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 202];
    auto current_out = &data_out[c * 360];
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
      const auto axxx_0 = current_data[136];
      const auto axxy_0 = current_data[137];
      const auto axyy_0 = current_data[138];
      const auto ayyy_0 = current_data[139];
      const auto axxz_0 = current_data[147];
      const auto axyz_0 = current_data[148];
      const auto ayyz_0 = current_data[149];
      const auto axzz_0 = current_data[157];
      const auto ayzz_0 = current_data[158];
      const auto azzz_0 = current_data[166];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[0] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[1] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[2] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[3] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[4] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[5] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[6] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[7] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[8] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[9] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[137];
      const auto axxy_0 = current_data[138];
      const auto axyy_0 = current_data[139];
      const auto ayyy_0 = current_data[140];
      const auto axxz_0 = current_data[148];
      const auto axyz_0 = current_data[149];
      const auto ayyz_0 = current_data[150];
      const auto axzz_0 = current_data[158];
      const auto ayzz_0 = current_data[159];
      const auto azzz_0 = current_data[167];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[10] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[11] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[12] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[13] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[14] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[15] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[16] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[17] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[18] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[19] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[138];
      const auto axxy_0 = current_data[139];
      const auto axyy_0 = current_data[140];
      const auto ayyy_0 = current_data[141];
      const auto axxz_0 = current_data[149];
      const auto axyz_0 = current_data[150];
      const auto ayyz_0 = current_data[151];
      const auto axzz_0 = current_data[159];
      const auto ayzz_0 = current_data[160];
      const auto azzz_0 = current_data[168];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[20] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[21] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[22] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[23] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[24] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[25] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[26] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[27] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[28] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[29] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[139];
      const auto axxy_0 = current_data[140];
      const auto axyy_0 = current_data[141];
      const auto ayyy_0 = current_data[142];
      const auto axxz_0 = current_data[150];
      const auto axyz_0 = current_data[151];
      const auto ayyz_0 = current_data[152];
      const auto axzz_0 = current_data[160];
      const auto ayzz_0 = current_data[161];
      const auto azzz_0 = current_data[169];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[30] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[31] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[32] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[33] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[34] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[35] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[36] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[37] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[38] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[39] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[140];
      const auto axxy_0 = current_data[141];
      const auto axyy_0 = current_data[142];
      const auto ayyy_0 = current_data[143];
      const auto axxz_0 = current_data[151];
      const auto axyz_0 = current_data[152];
      const auto ayyz_0 = current_data[153];
      const auto axzz_0 = current_data[161];
      const auto ayzz_0 = current_data[162];
      const auto azzz_0 = current_data[170];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[40] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[41] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[42] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[43] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[44] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[45] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[46] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[47] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[48] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[49] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[141];
      const auto axxy_0 = current_data[142];
      const auto axyy_0 = current_data[143];
      const auto ayyy_0 = current_data[144];
      const auto axxz_0 = current_data[152];
      const auto axyz_0 = current_data[153];
      const auto ayyz_0 = current_data[154];
      const auto axzz_0 = current_data[162];
      const auto ayzz_0 = current_data[163];
      const auto azzz_0 = current_data[171];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[50] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[51] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[52] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[53] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[54] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[55] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[56] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[57] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[58] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[59] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[142];
      const auto axxy_0 = current_data[143];
      const auto axyy_0 = current_data[144];
      const auto ayyy_0 = current_data[145];
      const auto axxz_0 = current_data[153];
      const auto axyz_0 = current_data[154];
      const auto ayyz_0 = current_data[155];
      const auto axzz_0 = current_data[163];
      const auto ayzz_0 = current_data[164];
      const auto azzz_0 = current_data[172];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[60] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[61] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[62] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[63] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[64] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[65] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[66] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[67] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[68] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[69] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[143];
      const auto axxy_0 = current_data[144];
      const auto axyy_0 = current_data[145];
      const auto ayyy_0 = current_data[146];
      const auto axxz_0 = current_data[154];
      const auto axyz_0 = current_data[155];
      const auto ayyz_0 = current_data[156];
      const auto axzz_0 = current_data[164];
      const auto ayzz_0 = current_data[165];
      const auto azzz_0 = current_data[173];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[70] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[71] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[72] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[73] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[74] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[75] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[76] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[77] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[78] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[79] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[147];
      const auto axxy_0 = current_data[148];
      const auto axyy_0 = current_data[149];
      const auto ayyy_0 = current_data[150];
      const auto axxz_0 = current_data[157];
      const auto axyz_0 = current_data[158];
      const auto ayyz_0 = current_data[159];
      const auto axzz_0 = current_data[166];
      const auto ayzz_0 = current_data[167];
      const auto azzz_0 = current_data[174];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[80] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[81] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[82] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[83] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[84] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[85] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[86] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[87] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[88] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[89] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[148];
      const auto axxy_0 = current_data[149];
      const auto axyy_0 = current_data[150];
      const auto ayyy_0 = current_data[151];
      const auto axxz_0 = current_data[158];
      const auto axyz_0 = current_data[159];
      const auto ayyz_0 = current_data[160];
      const auto axzz_0 = current_data[167];
      const auto ayzz_0 = current_data[168];
      const auto azzz_0 = current_data[175];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[90] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[91] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[92] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[93] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[94] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[95] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[96] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[97] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[98] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[99] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[149];
      const auto axxy_0 = current_data[150];
      const auto axyy_0 = current_data[151];
      const auto ayyy_0 = current_data[152];
      const auto axxz_0 = current_data[159];
      const auto axyz_0 = current_data[160];
      const auto ayyz_0 = current_data[161];
      const auto axzz_0 = current_data[168];
      const auto ayzz_0 = current_data[169];
      const auto azzz_0 = current_data[176];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[100] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[101] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[102] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[103] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[104] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[105] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[106] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[107] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[108] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[109] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[150];
      const auto axxy_0 = current_data[151];
      const auto axyy_0 = current_data[152];
      const auto ayyy_0 = current_data[153];
      const auto axxz_0 = current_data[160];
      const auto axyz_0 = current_data[161];
      const auto ayyz_0 = current_data[162];
      const auto axzz_0 = current_data[169];
      const auto ayzz_0 = current_data[170];
      const auto azzz_0 = current_data[177];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[110] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[111] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[112] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[113] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[114] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[115] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[116] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[117] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[118] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[119] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[151];
      const auto axxy_0 = current_data[152];
      const auto axyy_0 = current_data[153];
      const auto ayyy_0 = current_data[154];
      const auto axxz_0 = current_data[161];
      const auto axyz_0 = current_data[162];
      const auto ayyz_0 = current_data[163];
      const auto axzz_0 = current_data[170];
      const auto ayzz_0 = current_data[171];
      const auto azzz_0 = current_data[178];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[120] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[121] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[122] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[123] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[124] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[125] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[126] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[127] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[128] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[129] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[152];
      const auto axxy_0 = current_data[153];
      const auto axyy_0 = current_data[154];
      const auto ayyy_0 = current_data[155];
      const auto axxz_0 = current_data[162];
      const auto axyz_0 = current_data[163];
      const auto ayyz_0 = current_data[164];
      const auto axzz_0 = current_data[171];
      const auto ayzz_0 = current_data[172];
      const auto azzz_0 = current_data[179];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[130] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[131] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[132] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[133] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[134] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[135] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[136] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[137] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[138] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[139] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[153];
      const auto axxy_0 = current_data[154];
      const auto axyy_0 = current_data[155];
      const auto ayyy_0 = current_data[156];
      const auto axxz_0 = current_data[163];
      const auto axyz_0 = current_data[164];
      const auto ayyz_0 = current_data[165];
      const auto axzz_0 = current_data[172];
      const auto ayzz_0 = current_data[173];
      const auto azzz_0 = current_data[180];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[140] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[141] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[142] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[143] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[144] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[145] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[146] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[147] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[148] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[149] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[157];
      const auto axxy_0 = current_data[158];
      const auto axyy_0 = current_data[159];
      const auto ayyy_0 = current_data[160];
      const auto axxz_0 = current_data[166];
      const auto axyz_0 = current_data[167];
      const auto ayyz_0 = current_data[168];
      const auto axzz_0 = current_data[174];
      const auto ayzz_0 = current_data[175];
      const auto azzz_0 = current_data[181];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[150] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[151] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[152] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[153] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[154] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[155] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[156] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[157] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[158] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[159] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[158];
      const auto axxy_0 = current_data[159];
      const auto axyy_0 = current_data[160];
      const auto ayyy_0 = current_data[161];
      const auto axxz_0 = current_data[167];
      const auto axyz_0 = current_data[168];
      const auto ayyz_0 = current_data[169];
      const auto axzz_0 = current_data[175];
      const auto ayzz_0 = current_data[176];
      const auto azzz_0 = current_data[182];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[160] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[161] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[162] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[163] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[164] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[165] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[166] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[167] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[168] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[169] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[159];
      const auto axxy_0 = current_data[160];
      const auto axyy_0 = current_data[161];
      const auto ayyy_0 = current_data[162];
      const auto axxz_0 = current_data[168];
      const auto axyz_0 = current_data[169];
      const auto ayyz_0 = current_data[170];
      const auto axzz_0 = current_data[176];
      const auto ayzz_0 = current_data[177];
      const auto azzz_0 = current_data[183];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[170] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[171] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[172] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[173] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[174] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[175] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[176] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[177] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[178] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[179] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[160];
      const auto axxy_0 = current_data[161];
      const auto axyy_0 = current_data[162];
      const auto ayyy_0 = current_data[163];
      const auto axxz_0 = current_data[169];
      const auto axyz_0 = current_data[170];
      const auto ayyz_0 = current_data[171];
      const auto axzz_0 = current_data[177];
      const auto ayzz_0 = current_data[178];
      const auto azzz_0 = current_data[184];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[180] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[181] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[182] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[183] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[184] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[185] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[186] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[187] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[188] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[189] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[161];
      const auto axxy_0 = current_data[162];
      const auto axyy_0 = current_data[163];
      const auto ayyy_0 = current_data[164];
      const auto axxz_0 = current_data[170];
      const auto axyz_0 = current_data[171];
      const auto ayyz_0 = current_data[172];
      const auto axzz_0 = current_data[178];
      const auto ayzz_0 = current_data[179];
      const auto azzz_0 = current_data[185];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[190] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[191] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[192] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[193] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[194] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[195] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[196] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[197] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[198] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[199] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[162];
      const auto axxy_0 = current_data[163];
      const auto axyy_0 = current_data[164];
      const auto ayyy_0 = current_data[165];
      const auto axxz_0 = current_data[171];
      const auto axyz_0 = current_data[172];
      const auto ayyz_0 = current_data[173];
      const auto axzz_0 = current_data[179];
      const auto ayzz_0 = current_data[180];
      const auto azzz_0 = current_data[186];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[200] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[201] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[202] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[203] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[204] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[205] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[206] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[207] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[208] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[209] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[166];
      const auto axxy_0 = current_data[167];
      const auto axyy_0 = current_data[168];
      const auto ayyy_0 = current_data[169];
      const auto axxz_0 = current_data[174];
      const auto axyz_0 = current_data[175];
      const auto ayyz_0 = current_data[176];
      const auto axzz_0 = current_data[181];
      const auto ayzz_0 = current_data[182];
      const auto azzz_0 = current_data[187];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[210] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[211] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[212] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[213] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[214] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[215] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[216] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[217] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[218] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[219] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[167];
      const auto axxy_0 = current_data[168];
      const auto axyy_0 = current_data[169];
      const auto ayyy_0 = current_data[170];
      const auto axxz_0 = current_data[175];
      const auto axyz_0 = current_data[176];
      const auto ayyz_0 = current_data[177];
      const auto axzz_0 = current_data[182];
      const auto ayzz_0 = current_data[183];
      const auto azzz_0 = current_data[188];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[220] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[221] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[222] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[223] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[224] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[225] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[226] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[227] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[228] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[229] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[168];
      const auto axxy_0 = current_data[169];
      const auto axyy_0 = current_data[170];
      const auto ayyy_0 = current_data[171];
      const auto axxz_0 = current_data[176];
      const auto axyz_0 = current_data[177];
      const auto ayyz_0 = current_data[178];
      const auto axzz_0 = current_data[183];
      const auto ayzz_0 = current_data[184];
      const auto azzz_0 = current_data[189];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[230] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[231] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[232] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[233] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[234] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[235] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[236] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[237] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[238] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[239] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[169];
      const auto axxy_0 = current_data[170];
      const auto axyy_0 = current_data[171];
      const auto ayyy_0 = current_data[172];
      const auto axxz_0 = current_data[177];
      const auto axyz_0 = current_data[178];
      const auto ayyz_0 = current_data[179];
      const auto axzz_0 = current_data[184];
      const auto ayzz_0 = current_data[185];
      const auto azzz_0 = current_data[190];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[240] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[241] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[242] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[243] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[244] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[245] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[246] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[247] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[248] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[249] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[170];
      const auto axxy_0 = current_data[171];
      const auto axyy_0 = current_data[172];
      const auto ayyy_0 = current_data[173];
      const auto axxz_0 = current_data[178];
      const auto axyz_0 = current_data[179];
      const auto ayyz_0 = current_data[180];
      const auto axzz_0 = current_data[185];
      const auto ayzz_0 = current_data[186];
      const auto azzz_0 = current_data[191];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[250] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[251] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[252] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[253] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[254] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[255] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[256] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[257] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[258] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[259] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[174];
      const auto axxy_0 = current_data[175];
      const auto axyy_0 = current_data[176];
      const auto ayyy_0 = current_data[177];
      const auto axxz_0 = current_data[181];
      const auto axyz_0 = current_data[182];
      const auto ayyz_0 = current_data[183];
      const auto axzz_0 = current_data[187];
      const auto ayzz_0 = current_data[188];
      const auto azzz_0 = current_data[192];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[260] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[261] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[262] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[263] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[264] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[265] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[266] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[267] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[268] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[269] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[175];
      const auto axxy_0 = current_data[176];
      const auto axyy_0 = current_data[177];
      const auto ayyy_0 = current_data[178];
      const auto axxz_0 = current_data[182];
      const auto axyz_0 = current_data[183];
      const auto ayyz_0 = current_data[184];
      const auto axzz_0 = current_data[188];
      const auto ayzz_0 = current_data[189];
      const auto azzz_0 = current_data[193];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[270] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[271] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[272] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[273] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[274] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[275] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[276] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[277] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[278] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[279] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[176];
      const auto axxy_0 = current_data[177];
      const auto axyy_0 = current_data[178];
      const auto ayyy_0 = current_data[179];
      const auto axxz_0 = current_data[183];
      const auto axyz_0 = current_data[184];
      const auto ayyz_0 = current_data[185];
      const auto axzz_0 = current_data[189];
      const auto ayzz_0 = current_data[190];
      const auto azzz_0 = current_data[194];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[280] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[281] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[282] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[283] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[284] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[285] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[286] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[287] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[288] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[289] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[177];
      const auto axxy_0 = current_data[178];
      const auto axyy_0 = current_data[179];
      const auto ayyy_0 = current_data[180];
      const auto axxz_0 = current_data[184];
      const auto axyz_0 = current_data[185];
      const auto ayyz_0 = current_data[186];
      const auto axzz_0 = current_data[190];
      const auto ayzz_0 = current_data[191];
      const auto azzz_0 = current_data[195];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[290] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[291] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[292] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[293] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[294] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[295] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[296] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[297] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[298] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[299] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[181];
      const auto axxy_0 = current_data[182];
      const auto axyy_0 = current_data[183];
      const auto ayyy_0 = current_data[184];
      const auto axxz_0 = current_data[187];
      const auto axyz_0 = current_data[188];
      const auto ayyz_0 = current_data[189];
      const auto axzz_0 = current_data[192];
      const auto ayzz_0 = current_data[193];
      const auto azzz_0 = current_data[196];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[300] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[301] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[302] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[303] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[304] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[305] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[306] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[307] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[308] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[309] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[182];
      const auto axxy_0 = current_data[183];
      const auto axyy_0 = current_data[184];
      const auto ayyy_0 = current_data[185];
      const auto axxz_0 = current_data[188];
      const auto axyz_0 = current_data[189];
      const auto ayyz_0 = current_data[190];
      const auto axzz_0 = current_data[193];
      const auto ayzz_0 = current_data[194];
      const auto azzz_0 = current_data[197];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[310] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[311] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[312] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[313] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[314] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[315] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[316] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[317] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[318] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[319] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[183];
      const auto axxy_0 = current_data[184];
      const auto axyy_0 = current_data[185];
      const auto ayyy_0 = current_data[186];
      const auto axxz_0 = current_data[189];
      const auto axyz_0 = current_data[190];
      const auto ayyz_0 = current_data[191];
      const auto axzz_0 = current_data[194];
      const auto ayzz_0 = current_data[195];
      const auto azzz_0 = current_data[198];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[320] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[321] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[322] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[323] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[324] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[325] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[326] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[327] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[328] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[329] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[187];
      const auto axxy_0 = current_data[188];
      const auto axyy_0 = current_data[189];
      const auto ayyy_0 = current_data[190];
      const auto axxz_0 = current_data[192];
      const auto axyz_0 = current_data[193];
      const auto ayyz_0 = current_data[194];
      const auto axzz_0 = current_data[196];
      const auto ayzz_0 = current_data[197];
      const auto azzz_0 = current_data[199];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[330] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[331] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[332] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[333] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[334] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[335] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[336] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[337] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[338] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[339] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[188];
      const auto axxy_0 = current_data[189];
      const auto axyy_0 = current_data[190];
      const auto ayyy_0 = current_data[191];
      const auto axxz_0 = current_data[193];
      const auto axyz_0 = current_data[194];
      const auto ayyz_0 = current_data[195];
      const auto axzz_0 = current_data[197];
      const auto ayzz_0 = current_data[198];
      const auto azzz_0 = current_data[200];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[340] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[341] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[342] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[343] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[344] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[345] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[346] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[347] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[348] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[349] = az_zz + AB[2] * a0_zz; // a0_zzz

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
      const auto axxx_0 = current_data[192];
      const auto axxy_0 = current_data[193];
      const auto axyy_0 = current_data[194];
      const auto ayyy_0 = current_data[195];
      const auto axxz_0 = current_data[196];
      const auto axyz_0 = current_data[197];
      const auto ayyz_0 = current_data[198];
      const auto axzz_0 = current_data[199];
      const auto ayzz_0 = current_data[200];
      const auto azzz_0 = current_data[201];

      const auto a0_x = ax_0 + AB[0] * a0_0;
      const auto a0_y = ay_0 + AB[1] * a0_0;
      const auto a0_z = az_0 + AB[2] * a0_0;

      const auto ax_x = axx_0 + AB[0] * ax_0;
      const auto ax_y = axy_0 + AB[1] * ax_0;
      const auto ay_y = ayy_0 + AB[1] * ay_0;
      const auto ax_z = axz_0 + AB[2] * ax_0;
      const auto ay_z = ayz_0 + AB[2] * ay_0;
      const auto az_z = azz_0 + AB[2] * az_0;

      const auto axx_x = axxx_0 + AB[0] * axx_0;
      const auto axx_y = axxy_0 + AB[1] * axx_0;
      const auto axy_y = axyy_0 + AB[1] * axy_0;
      const auto ayy_y = ayyy_0 + AB[1] * ayy_0;
      const auto axx_z = axxz_0 + AB[2] * axx_0;
      const auto axy_z = axyz_0 + AB[2] * axy_0;
      const auto ayy_z = ayyz_0 + AB[2] * ayy_0;
      const auto axz_z = axzz_0 + AB[2] * axz_0;
      const auto ayz_z = ayzz_0 + AB[2] * ayz_0;
      const auto azz_z = azzz_0 + AB[2] * azz_0;

      const auto a0_xx = ax_x + AB[0] * a0_x;
      const auto a0_xy = ax_y + AB[0] * a0_y;
      const auto a0_yy = ay_y + AB[1] * a0_y;
      const auto a0_xz = ax_z + AB[0] * a0_z;
      const auto a0_yz = ay_z + AB[1] * a0_z;
      const auto a0_zz = az_z + AB[2] * a0_z;

      const auto ax_xx = axx_x + AB[0] * ax_x;
      const auto ax_xy = axx_y + AB[0] * ax_y;
      const auto ax_yy = axy_y + AB[1] * ax_y;
      const auto ay_yy = ayy_y + AB[1] * ay_y;
      const auto ax_xz = axx_z + AB[0] * ax_z;
      const auto ax_yz = axy_z + AB[1] * ax_z;
      const auto ay_yz = ayy_z + AB[1] * ay_z;
      const auto ax_zz = axz_z + AB[2] * ax_z;
      const auto ay_zz = ayz_z + AB[2] * ay_z;
      const auto az_zz = azz_z + AB[2] * az_z;

      current_out[350] = ax_xx + AB[0] * a0_xx; // a0_xxx
      current_out[351] = ax_xy + AB[0] * a0_xy; // a0_xxy
      current_out[352] = ax_yy + AB[0] * a0_yy; // a0_xyy
      current_out[353] = ay_yy + AB[1] * a0_yy; // a0_yyy
      current_out[354] = ax_xz + AB[0] * a0_xz; // a0_xxz
      current_out[355] = ax_yz + AB[0] * a0_yz; // a0_xyz
      current_out[356] = ay_yz + AB[1] * a0_yz; // a0_yyz
      current_out[357] = ax_zz + AB[0] * a0_zz; // a0_xzz
      current_out[358] = ay_zz + AB[1] * a0_zz; // a0_yzz
      current_out[359] = az_zz + AB[2] * a0_zz; // a0_zzz

    }
  }
}
#endif

