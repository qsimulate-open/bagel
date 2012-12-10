//
// BAGEL - Parallel electron correlation program.
// Filename: _hrr_80_44.cc
// Copyright (C) 2009 Toru Shiozaki
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

#include <src/rysint/hrrlist.h>
#include <array>
#include <algorithm>

using namespace std;
using namespace bagel;

void HRRList::perform_HRR_80_44(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    const double* current_data = &data_start[c * 145];
    double* current_out = &data_out[c * 225];
   {
     //current index a: xxxx
      const double a0_0 = current_data[0];
      const double ax_0 = current_data[15];
      const double ay_0 = current_data[16];
      const double az_0 = current_data[21];
      const double axx_0 = current_data[36];
      const double axy_0 = current_data[37];
      const double ayy_0 = current_data[38];
      const double axz_0 = current_data[43];
      const double ayz_0 = current_data[44];
      const double azz_0 = current_data[49];
      const double axxx_0 = current_data[64];
      const double axxy_0 = current_data[65];
      const double axyy_0 = current_data[66];
      const double ayyy_0 = current_data[67];
      const double axxz_0 = current_data[72];
      const double axyz_0 = current_data[73];
      const double ayyz_0 = current_data[74];
      const double axzz_0 = current_data[79];
      const double ayzz_0 = current_data[80];
      const double azzz_0 = current_data[85];
      const double axxxx_0 = current_data[100];
      const double axxxy_0 = current_data[101];
      const double axxyy_0 = current_data[102];
      const double axyyy_0 = current_data[103];
      const double ayyyy_0 = current_data[104];
      const double axxxz_0 = current_data[109];
      const double axxyz_0 = current_data[110];
      const double axyyz_0 = current_data[111];
      const double ayyyz_0 = current_data[112];
      const double axxzz_0 = current_data[117];
      const double axyzz_0 = current_data[118];
      const double ayyzz_0 = current_data[119];
      const double axzzz_0 = current_data[124];
      const double ayzzz_0 = current_data[125];
      const double azzzz_0 = current_data[130];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[0] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[1] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[2] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[3] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[4] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[5] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[6] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[7] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[8] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[9] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[10] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[11] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[12] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[13] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[14] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: xxxy
      const double a0_0 = current_data[1];
      const double ax_0 = current_data[16];
      const double ay_0 = current_data[17];
      const double az_0 = current_data[22];
      const double axx_0 = current_data[37];
      const double axy_0 = current_data[38];
      const double ayy_0 = current_data[39];
      const double axz_0 = current_data[44];
      const double ayz_0 = current_data[45];
      const double azz_0 = current_data[50];
      const double axxx_0 = current_data[65];
      const double axxy_0 = current_data[66];
      const double axyy_0 = current_data[67];
      const double ayyy_0 = current_data[68];
      const double axxz_0 = current_data[73];
      const double axyz_0 = current_data[74];
      const double ayyz_0 = current_data[75];
      const double axzz_0 = current_data[80];
      const double ayzz_0 = current_data[81];
      const double azzz_0 = current_data[86];
      const double axxxx_0 = current_data[101];
      const double axxxy_0 = current_data[102];
      const double axxyy_0 = current_data[103];
      const double axyyy_0 = current_data[104];
      const double ayyyy_0 = current_data[105];
      const double axxxz_0 = current_data[110];
      const double axxyz_0 = current_data[111];
      const double axyyz_0 = current_data[112];
      const double ayyyz_0 = current_data[113];
      const double axxzz_0 = current_data[118];
      const double axyzz_0 = current_data[119];
      const double ayyzz_0 = current_data[120];
      const double axzzz_0 = current_data[125];
      const double ayzzz_0 = current_data[126];
      const double azzzz_0 = current_data[131];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[15] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[16] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[17] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[18] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[19] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[20] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[21] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[22] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[23] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[24] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[25] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[26] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[27] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[28] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[29] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: xxyy
      const double a0_0 = current_data[2];
      const double ax_0 = current_data[17];
      const double ay_0 = current_data[18];
      const double az_0 = current_data[23];
      const double axx_0 = current_data[38];
      const double axy_0 = current_data[39];
      const double ayy_0 = current_data[40];
      const double axz_0 = current_data[45];
      const double ayz_0 = current_data[46];
      const double azz_0 = current_data[51];
      const double axxx_0 = current_data[66];
      const double axxy_0 = current_data[67];
      const double axyy_0 = current_data[68];
      const double ayyy_0 = current_data[69];
      const double axxz_0 = current_data[74];
      const double axyz_0 = current_data[75];
      const double ayyz_0 = current_data[76];
      const double axzz_0 = current_data[81];
      const double ayzz_0 = current_data[82];
      const double azzz_0 = current_data[87];
      const double axxxx_0 = current_data[102];
      const double axxxy_0 = current_data[103];
      const double axxyy_0 = current_data[104];
      const double axyyy_0 = current_data[105];
      const double ayyyy_0 = current_data[106];
      const double axxxz_0 = current_data[111];
      const double axxyz_0 = current_data[112];
      const double axyyz_0 = current_data[113];
      const double ayyyz_0 = current_data[114];
      const double axxzz_0 = current_data[119];
      const double axyzz_0 = current_data[120];
      const double ayyzz_0 = current_data[121];
      const double axzzz_0 = current_data[126];
      const double ayzzz_0 = current_data[127];
      const double azzzz_0 = current_data[132];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[30] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[31] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[32] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[33] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[34] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[35] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[36] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[37] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[38] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[39] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[40] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[41] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[42] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[43] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[44] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: xyyy
      const double a0_0 = current_data[3];
      const double ax_0 = current_data[18];
      const double ay_0 = current_data[19];
      const double az_0 = current_data[24];
      const double axx_0 = current_data[39];
      const double axy_0 = current_data[40];
      const double ayy_0 = current_data[41];
      const double axz_0 = current_data[46];
      const double ayz_0 = current_data[47];
      const double azz_0 = current_data[52];
      const double axxx_0 = current_data[67];
      const double axxy_0 = current_data[68];
      const double axyy_0 = current_data[69];
      const double ayyy_0 = current_data[70];
      const double axxz_0 = current_data[75];
      const double axyz_0 = current_data[76];
      const double ayyz_0 = current_data[77];
      const double axzz_0 = current_data[82];
      const double ayzz_0 = current_data[83];
      const double azzz_0 = current_data[88];
      const double axxxx_0 = current_data[103];
      const double axxxy_0 = current_data[104];
      const double axxyy_0 = current_data[105];
      const double axyyy_0 = current_data[106];
      const double ayyyy_0 = current_data[107];
      const double axxxz_0 = current_data[112];
      const double axxyz_0 = current_data[113];
      const double axyyz_0 = current_data[114];
      const double ayyyz_0 = current_data[115];
      const double axxzz_0 = current_data[120];
      const double axyzz_0 = current_data[121];
      const double ayyzz_0 = current_data[122];
      const double axzzz_0 = current_data[127];
      const double ayzzz_0 = current_data[128];
      const double azzzz_0 = current_data[133];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[45] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[46] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[47] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[48] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[49] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[50] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[51] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[52] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[53] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[54] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[55] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[56] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[57] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[58] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[59] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: yyyy
      const double a0_0 = current_data[4];
      const double ax_0 = current_data[19];
      const double ay_0 = current_data[20];
      const double az_0 = current_data[25];
      const double axx_0 = current_data[40];
      const double axy_0 = current_data[41];
      const double ayy_0 = current_data[42];
      const double axz_0 = current_data[47];
      const double ayz_0 = current_data[48];
      const double azz_0 = current_data[53];
      const double axxx_0 = current_data[68];
      const double axxy_0 = current_data[69];
      const double axyy_0 = current_data[70];
      const double ayyy_0 = current_data[71];
      const double axxz_0 = current_data[76];
      const double axyz_0 = current_data[77];
      const double ayyz_0 = current_data[78];
      const double axzz_0 = current_data[83];
      const double ayzz_0 = current_data[84];
      const double azzz_0 = current_data[89];
      const double axxxx_0 = current_data[104];
      const double axxxy_0 = current_data[105];
      const double axxyy_0 = current_data[106];
      const double axyyy_0 = current_data[107];
      const double ayyyy_0 = current_data[108];
      const double axxxz_0 = current_data[113];
      const double axxyz_0 = current_data[114];
      const double axyyz_0 = current_data[115];
      const double ayyyz_0 = current_data[116];
      const double axxzz_0 = current_data[121];
      const double axyzz_0 = current_data[122];
      const double ayyzz_0 = current_data[123];
      const double axzzz_0 = current_data[128];
      const double ayzzz_0 = current_data[129];
      const double azzzz_0 = current_data[134];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[60] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[61] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[62] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[63] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[64] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[65] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[66] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[67] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[68] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[69] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[70] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[71] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[72] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[73] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[74] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: xxxz
      const double a0_0 = current_data[5];
      const double ax_0 = current_data[21];
      const double ay_0 = current_data[22];
      const double az_0 = current_data[26];
      const double axx_0 = current_data[43];
      const double axy_0 = current_data[44];
      const double ayy_0 = current_data[45];
      const double axz_0 = current_data[49];
      const double ayz_0 = current_data[50];
      const double azz_0 = current_data[54];
      const double axxx_0 = current_data[72];
      const double axxy_0 = current_data[73];
      const double axyy_0 = current_data[74];
      const double ayyy_0 = current_data[75];
      const double axxz_0 = current_data[79];
      const double axyz_0 = current_data[80];
      const double ayyz_0 = current_data[81];
      const double axzz_0 = current_data[85];
      const double ayzz_0 = current_data[86];
      const double azzz_0 = current_data[90];
      const double axxxx_0 = current_data[109];
      const double axxxy_0 = current_data[110];
      const double axxyy_0 = current_data[111];
      const double axyyy_0 = current_data[112];
      const double ayyyy_0 = current_data[113];
      const double axxxz_0 = current_data[117];
      const double axxyz_0 = current_data[118];
      const double axyyz_0 = current_data[119];
      const double ayyyz_0 = current_data[120];
      const double axxzz_0 = current_data[124];
      const double axyzz_0 = current_data[125];
      const double ayyzz_0 = current_data[126];
      const double axzzz_0 = current_data[130];
      const double ayzzz_0 = current_data[131];
      const double azzzz_0 = current_data[135];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[75] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[76] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[77] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[78] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[79] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[80] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[81] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[82] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[83] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[84] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[85] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[86] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[87] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[88] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[89] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: xxyz
      const double a0_0 = current_data[6];
      const double ax_0 = current_data[22];
      const double ay_0 = current_data[23];
      const double az_0 = current_data[27];
      const double axx_0 = current_data[44];
      const double axy_0 = current_data[45];
      const double ayy_0 = current_data[46];
      const double axz_0 = current_data[50];
      const double ayz_0 = current_data[51];
      const double azz_0 = current_data[55];
      const double axxx_0 = current_data[73];
      const double axxy_0 = current_data[74];
      const double axyy_0 = current_data[75];
      const double ayyy_0 = current_data[76];
      const double axxz_0 = current_data[80];
      const double axyz_0 = current_data[81];
      const double ayyz_0 = current_data[82];
      const double axzz_0 = current_data[86];
      const double ayzz_0 = current_data[87];
      const double azzz_0 = current_data[91];
      const double axxxx_0 = current_data[110];
      const double axxxy_0 = current_data[111];
      const double axxyy_0 = current_data[112];
      const double axyyy_0 = current_data[113];
      const double ayyyy_0 = current_data[114];
      const double axxxz_0 = current_data[118];
      const double axxyz_0 = current_data[119];
      const double axyyz_0 = current_data[120];
      const double ayyyz_0 = current_data[121];
      const double axxzz_0 = current_data[125];
      const double axyzz_0 = current_data[126];
      const double ayyzz_0 = current_data[127];
      const double axzzz_0 = current_data[131];
      const double ayzzz_0 = current_data[132];
      const double azzzz_0 = current_data[136];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[90] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[91] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[92] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[93] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[94] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[95] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[96] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[97] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[98] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[99] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[100] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[101] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[102] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[103] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[104] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: xyyz
      const double a0_0 = current_data[7];
      const double ax_0 = current_data[23];
      const double ay_0 = current_data[24];
      const double az_0 = current_data[28];
      const double axx_0 = current_data[45];
      const double axy_0 = current_data[46];
      const double ayy_0 = current_data[47];
      const double axz_0 = current_data[51];
      const double ayz_0 = current_data[52];
      const double azz_0 = current_data[56];
      const double axxx_0 = current_data[74];
      const double axxy_0 = current_data[75];
      const double axyy_0 = current_data[76];
      const double ayyy_0 = current_data[77];
      const double axxz_0 = current_data[81];
      const double axyz_0 = current_data[82];
      const double ayyz_0 = current_data[83];
      const double axzz_0 = current_data[87];
      const double ayzz_0 = current_data[88];
      const double azzz_0 = current_data[92];
      const double axxxx_0 = current_data[111];
      const double axxxy_0 = current_data[112];
      const double axxyy_0 = current_data[113];
      const double axyyy_0 = current_data[114];
      const double ayyyy_0 = current_data[115];
      const double axxxz_0 = current_data[119];
      const double axxyz_0 = current_data[120];
      const double axyyz_0 = current_data[121];
      const double ayyyz_0 = current_data[122];
      const double axxzz_0 = current_data[126];
      const double axyzz_0 = current_data[127];
      const double ayyzz_0 = current_data[128];
      const double axzzz_0 = current_data[132];
      const double ayzzz_0 = current_data[133];
      const double azzzz_0 = current_data[137];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[105] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[106] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[107] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[108] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[109] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[110] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[111] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[112] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[113] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[114] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[115] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[116] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[117] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[118] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[119] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: yyyz
      const double a0_0 = current_data[8];
      const double ax_0 = current_data[24];
      const double ay_0 = current_data[25];
      const double az_0 = current_data[29];
      const double axx_0 = current_data[46];
      const double axy_0 = current_data[47];
      const double ayy_0 = current_data[48];
      const double axz_0 = current_data[52];
      const double ayz_0 = current_data[53];
      const double azz_0 = current_data[57];
      const double axxx_0 = current_data[75];
      const double axxy_0 = current_data[76];
      const double axyy_0 = current_data[77];
      const double ayyy_0 = current_data[78];
      const double axxz_0 = current_data[82];
      const double axyz_0 = current_data[83];
      const double ayyz_0 = current_data[84];
      const double axzz_0 = current_data[88];
      const double ayzz_0 = current_data[89];
      const double azzz_0 = current_data[93];
      const double axxxx_0 = current_data[112];
      const double axxxy_0 = current_data[113];
      const double axxyy_0 = current_data[114];
      const double axyyy_0 = current_data[115];
      const double ayyyy_0 = current_data[116];
      const double axxxz_0 = current_data[120];
      const double axxyz_0 = current_data[121];
      const double axyyz_0 = current_data[122];
      const double ayyyz_0 = current_data[123];
      const double axxzz_0 = current_data[127];
      const double axyzz_0 = current_data[128];
      const double ayyzz_0 = current_data[129];
      const double axzzz_0 = current_data[133];
      const double ayzzz_0 = current_data[134];
      const double azzzz_0 = current_data[138];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[120] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[121] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[122] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[123] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[124] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[125] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[126] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[127] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[128] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[129] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[130] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[131] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[132] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[133] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[134] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: xxzz
      const double a0_0 = current_data[9];
      const double ax_0 = current_data[26];
      const double ay_0 = current_data[27];
      const double az_0 = current_data[30];
      const double axx_0 = current_data[49];
      const double axy_0 = current_data[50];
      const double ayy_0 = current_data[51];
      const double axz_0 = current_data[54];
      const double ayz_0 = current_data[55];
      const double azz_0 = current_data[58];
      const double axxx_0 = current_data[79];
      const double axxy_0 = current_data[80];
      const double axyy_0 = current_data[81];
      const double ayyy_0 = current_data[82];
      const double axxz_0 = current_data[85];
      const double axyz_0 = current_data[86];
      const double ayyz_0 = current_data[87];
      const double axzz_0 = current_data[90];
      const double ayzz_0 = current_data[91];
      const double azzz_0 = current_data[94];
      const double axxxx_0 = current_data[117];
      const double axxxy_0 = current_data[118];
      const double axxyy_0 = current_data[119];
      const double axyyy_0 = current_data[120];
      const double ayyyy_0 = current_data[121];
      const double axxxz_0 = current_data[124];
      const double axxyz_0 = current_data[125];
      const double axyyz_0 = current_data[126];
      const double ayyyz_0 = current_data[127];
      const double axxzz_0 = current_data[130];
      const double axyzz_0 = current_data[131];
      const double ayyzz_0 = current_data[132];
      const double axzzz_0 = current_data[135];
      const double ayzzz_0 = current_data[136];
      const double azzzz_0 = current_data[139];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[135] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[136] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[137] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[138] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[139] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[140] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[141] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[142] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[143] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[144] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[145] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[146] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[147] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[148] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[149] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: xyzz
      const double a0_0 = current_data[10];
      const double ax_0 = current_data[27];
      const double ay_0 = current_data[28];
      const double az_0 = current_data[31];
      const double axx_0 = current_data[50];
      const double axy_0 = current_data[51];
      const double ayy_0 = current_data[52];
      const double axz_0 = current_data[55];
      const double ayz_0 = current_data[56];
      const double azz_0 = current_data[59];
      const double axxx_0 = current_data[80];
      const double axxy_0 = current_data[81];
      const double axyy_0 = current_data[82];
      const double ayyy_0 = current_data[83];
      const double axxz_0 = current_data[86];
      const double axyz_0 = current_data[87];
      const double ayyz_0 = current_data[88];
      const double axzz_0 = current_data[91];
      const double ayzz_0 = current_data[92];
      const double azzz_0 = current_data[95];
      const double axxxx_0 = current_data[118];
      const double axxxy_0 = current_data[119];
      const double axxyy_0 = current_data[120];
      const double axyyy_0 = current_data[121];
      const double ayyyy_0 = current_data[122];
      const double axxxz_0 = current_data[125];
      const double axxyz_0 = current_data[126];
      const double axyyz_0 = current_data[127];
      const double ayyyz_0 = current_data[128];
      const double axxzz_0 = current_data[131];
      const double axyzz_0 = current_data[132];
      const double ayyzz_0 = current_data[133];
      const double axzzz_0 = current_data[136];
      const double ayzzz_0 = current_data[137];
      const double azzzz_0 = current_data[140];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[150] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[151] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[152] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[153] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[154] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[155] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[156] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[157] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[158] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[159] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[160] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[161] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[162] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[163] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[164] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: yyzz
      const double a0_0 = current_data[11];
      const double ax_0 = current_data[28];
      const double ay_0 = current_data[29];
      const double az_0 = current_data[32];
      const double axx_0 = current_data[51];
      const double axy_0 = current_data[52];
      const double ayy_0 = current_data[53];
      const double axz_0 = current_data[56];
      const double ayz_0 = current_data[57];
      const double azz_0 = current_data[60];
      const double axxx_0 = current_data[81];
      const double axxy_0 = current_data[82];
      const double axyy_0 = current_data[83];
      const double ayyy_0 = current_data[84];
      const double axxz_0 = current_data[87];
      const double axyz_0 = current_data[88];
      const double ayyz_0 = current_data[89];
      const double axzz_0 = current_data[92];
      const double ayzz_0 = current_data[93];
      const double azzz_0 = current_data[96];
      const double axxxx_0 = current_data[119];
      const double axxxy_0 = current_data[120];
      const double axxyy_0 = current_data[121];
      const double axyyy_0 = current_data[122];
      const double ayyyy_0 = current_data[123];
      const double axxxz_0 = current_data[126];
      const double axxyz_0 = current_data[127];
      const double axyyz_0 = current_data[128];
      const double ayyyz_0 = current_data[129];
      const double axxzz_0 = current_data[132];
      const double axyzz_0 = current_data[133];
      const double ayyzz_0 = current_data[134];
      const double axzzz_0 = current_data[137];
      const double ayzzz_0 = current_data[138];
      const double azzzz_0 = current_data[141];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[165] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[166] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[167] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[168] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[169] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[170] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[171] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[172] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[173] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[174] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[175] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[176] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[177] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[178] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[179] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: xzzz
      const double a0_0 = current_data[12];
      const double ax_0 = current_data[30];
      const double ay_0 = current_data[31];
      const double az_0 = current_data[33];
      const double axx_0 = current_data[54];
      const double axy_0 = current_data[55];
      const double ayy_0 = current_data[56];
      const double axz_0 = current_data[58];
      const double ayz_0 = current_data[59];
      const double azz_0 = current_data[61];
      const double axxx_0 = current_data[85];
      const double axxy_0 = current_data[86];
      const double axyy_0 = current_data[87];
      const double ayyy_0 = current_data[88];
      const double axxz_0 = current_data[90];
      const double axyz_0 = current_data[91];
      const double ayyz_0 = current_data[92];
      const double axzz_0 = current_data[94];
      const double ayzz_0 = current_data[95];
      const double azzz_0 = current_data[97];
      const double axxxx_0 = current_data[124];
      const double axxxy_0 = current_data[125];
      const double axxyy_0 = current_data[126];
      const double axyyy_0 = current_data[127];
      const double ayyyy_0 = current_data[128];
      const double axxxz_0 = current_data[130];
      const double axxyz_0 = current_data[131];
      const double axyyz_0 = current_data[132];
      const double ayyyz_0 = current_data[133];
      const double axxzz_0 = current_data[135];
      const double axyzz_0 = current_data[136];
      const double ayyzz_0 = current_data[137];
      const double axzzz_0 = current_data[139];
      const double ayzzz_0 = current_data[140];
      const double azzzz_0 = current_data[142];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[180] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[181] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[182] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[183] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[184] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[185] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[186] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[187] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[188] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[189] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[190] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[191] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[192] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[193] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[194] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: yzzz
      const double a0_0 = current_data[13];
      const double ax_0 = current_data[31];
      const double ay_0 = current_data[32];
      const double az_0 = current_data[34];
      const double axx_0 = current_data[55];
      const double axy_0 = current_data[56];
      const double ayy_0 = current_data[57];
      const double axz_0 = current_data[59];
      const double ayz_0 = current_data[60];
      const double azz_0 = current_data[62];
      const double axxx_0 = current_data[86];
      const double axxy_0 = current_data[87];
      const double axyy_0 = current_data[88];
      const double ayyy_0 = current_data[89];
      const double axxz_0 = current_data[91];
      const double axyz_0 = current_data[92];
      const double ayyz_0 = current_data[93];
      const double axzz_0 = current_data[95];
      const double ayzz_0 = current_data[96];
      const double azzz_0 = current_data[98];
      const double axxxx_0 = current_data[125];
      const double axxxy_0 = current_data[126];
      const double axxyy_0 = current_data[127];
      const double axyyy_0 = current_data[128];
      const double ayyyy_0 = current_data[129];
      const double axxxz_0 = current_data[131];
      const double axxyz_0 = current_data[132];
      const double axyyz_0 = current_data[133];
      const double ayyyz_0 = current_data[134];
      const double axxzz_0 = current_data[136];
      const double axyzz_0 = current_data[137];
      const double ayyzz_0 = current_data[138];
      const double axzzz_0 = current_data[140];
      const double ayzzz_0 = current_data[141];
      const double azzzz_0 = current_data[143];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[195] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[196] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[197] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[198] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[199] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[200] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[201] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[202] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[203] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[204] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[205] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[206] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[207] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[208] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[209] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
   {
     //current index a: zzzz
      const double a0_0 = current_data[14];
      const double ax_0 = current_data[33];
      const double ay_0 = current_data[34];
      const double az_0 = current_data[35];
      const double axx_0 = current_data[58];
      const double axy_0 = current_data[59];
      const double ayy_0 = current_data[60];
      const double axz_0 = current_data[61];
      const double ayz_0 = current_data[62];
      const double azz_0 = current_data[63];
      const double axxx_0 = current_data[90];
      const double axxy_0 = current_data[91];
      const double axyy_0 = current_data[92];
      const double ayyy_0 = current_data[93];
      const double axxz_0 = current_data[94];
      const double axyz_0 = current_data[95];
      const double ayyz_0 = current_data[96];
      const double axzz_0 = current_data[97];
      const double ayzz_0 = current_data[98];
      const double azzz_0 = current_data[99];
      const double axxxx_0 = current_data[130];
      const double axxxy_0 = current_data[131];
      const double axxyy_0 = current_data[132];
      const double axyyy_0 = current_data[133];
      const double ayyyy_0 = current_data[134];
      const double axxxz_0 = current_data[135];
      const double axxyz_0 = current_data[136];
      const double axyyz_0 = current_data[137];
      const double ayyyz_0 = current_data[138];
      const double axxzz_0 = current_data[139];
      const double axyzz_0 = current_data[140];
      const double ayyzz_0 = current_data[141];
      const double axzzz_0 = current_data[142];
      const double ayzzz_0 = current_data[143];
      const double azzzz_0 = current_data[144];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      const double axx_x = axxx_0 + AB[0] * axx_0;
      const double axx_y = axxy_0 + AB[1] * axx_0;
      const double axy_y = axyy_0 + AB[1] * axy_0;
      const double ayy_y = ayyy_0 + AB[1] * ayy_0;
      const double axx_z = axxz_0 + AB[2] * axx_0;
      const double axy_z = axyz_0 + AB[2] * axy_0;
      const double ayy_z = ayyz_0 + AB[2] * ayy_0;
      const double axz_z = axzz_0 + AB[2] * axz_0;
      const double ayz_z = ayzz_0 + AB[2] * ayz_0;
      const double azz_z = azzz_0 + AB[2] * azz_0;

      const double axxx_x = axxxx_0 + AB[0] * axxx_0;
      const double axxx_y = axxxy_0 + AB[1] * axxx_0;
      const double axxy_y = axxyy_0 + AB[1] * axxy_0;
      const double axyy_y = axyyy_0 + AB[1] * axyy_0;
      const double ayyy_y = ayyyy_0 + AB[1] * ayyy_0;
      const double axxx_z = axxxz_0 + AB[2] * axxx_0;
      const double axxy_z = axxyz_0 + AB[2] * axxy_0;
      const double axyy_z = axyyz_0 + AB[2] * axyy_0;
      const double ayyy_z = ayyyz_0 + AB[2] * ayyy_0;
      const double axxz_z = axxzz_0 + AB[2] * axxz_0;
      const double axyz_z = axyzz_0 + AB[2] * axyz_0;
      const double ayyz_z = ayyzz_0 + AB[2] * ayyz_0;
      const double axzz_z = axzzz_0 + AB[2] * axzz_0;
      const double ayzz_z = ayzzz_0 + AB[2] * ayzz_0;
      const double azzz_z = azzzz_0 + AB[2] * azzz_0;

      const double a0_xx = ax_x + AB[0] * a0_x;
      const double a0_xy = ax_y + AB[0] * a0_y;
      const double a0_yy = ay_y + AB[1] * a0_y;
      const double a0_xz = ax_z + AB[0] * a0_z;
      const double a0_yz = ay_z + AB[1] * a0_z;
      const double a0_zz = az_z + AB[2] * a0_z;

      const double ax_xx = axx_x + AB[0] * ax_x;
      const double ax_xy = axx_y + AB[0] * ax_y;
      const double ax_yy = axy_y + AB[1] * ax_y;
      const double ay_yy = ayy_y + AB[1] * ay_y;
      const double ax_xz = axx_z + AB[0] * ax_z;
      const double ax_yz = axy_z + AB[1] * ax_z;
      const double ay_yz = ayy_z + AB[1] * ay_z;
      const double ax_zz = axz_z + AB[2] * ax_z;
      const double ay_zz = ayz_z + AB[2] * ay_z;
      const double az_zz = azz_z + AB[2] * az_z;

      const double axx_xx = axxx_x + AB[0] * axx_x;
      const double axx_xy = axxx_y + AB[0] * axx_y;
      const double axx_yy = axxy_y + AB[1] * axx_y;
      const double axy_yy = axyy_y + AB[1] * axy_y;
      const double ayy_yy = ayyy_y + AB[1] * ayy_y;
      const double axx_xz = axxx_z + AB[0] * axx_z;
      const double axx_yz = axxy_z + AB[1] * axx_z;
      const double axy_yz = axyy_z + AB[1] * axy_z;
      const double ayy_yz = ayyy_z + AB[1] * ayy_z;
      const double axx_zz = axxz_z + AB[2] * axx_z;
      const double axy_zz = axyz_z + AB[2] * axy_z;
      const double ayy_zz = ayyz_z + AB[2] * ayy_z;
      const double axz_zz = axzz_z + AB[2] * axz_z;
      const double ayz_zz = ayzz_z + AB[2] * ayz_z;
      const double azz_zz = azzz_z + AB[2] * azz_z;

      const double a0_xxx = ax_xx + AB[0] * a0_xx;
      const double a0_xxy = ax_xy + AB[0] * a0_xy;
      const double a0_xyy = ax_yy + AB[0] * a0_yy;
      const double a0_yyy = ay_yy + AB[1] * a0_yy;
      const double a0_xxz = ax_xz + AB[0] * a0_xz;
      const double a0_xyz = ax_yz + AB[0] * a0_yz;
      const double a0_yyz = ay_yz + AB[1] * a0_yz;
      const double a0_xzz = ax_zz + AB[0] * a0_zz;
      const double a0_yzz = ay_zz + AB[1] * a0_zz;
      const double a0_zzz = az_zz + AB[2] * a0_zz;

      const double ax_xxx = axx_xx + AB[0] * ax_xx;
      const double ax_xxy = axx_xy + AB[0] * ax_xy;
      const double ax_xyy = axx_yy + AB[0] * ax_yy;
      const double ax_yyy = axy_yy + AB[1] * ax_yy;
      const double ay_yyy = ayy_yy + AB[1] * ay_yy;
      const double ax_xxz = axx_xz + AB[0] * ax_xz;
      const double ax_xyz = axx_yz + AB[0] * ax_yz;
      const double ax_yyz = axy_yz + AB[1] * ax_yz;
      const double ay_yyz = ayy_yz + AB[1] * ay_yz;
      const double ax_xzz = axx_zz + AB[0] * ax_zz;
      const double ax_yzz = axy_zz + AB[1] * ax_zz;
      const double ay_yzz = ayy_zz + AB[1] * ay_zz;
      const double ax_zzz = axz_zz + AB[2] * ax_zz;
      const double ay_zzz = ayz_zz + AB[2] * ay_zz;
      const double az_zzz = azz_zz + AB[2] * az_zz;

      current_out[210] = ax_xxx + AB[0] * a0_xxx; // a0_xxxx
      current_out[211] = ax_xxy + AB[0] * a0_xxy; // a0_xxxy
      current_out[212] = ax_xyy + AB[0] * a0_xyy; // a0_xxyy
      current_out[213] = ax_yyy + AB[0] * a0_yyy; // a0_xyyy
      current_out[214] = ay_yyy + AB[1] * a0_yyy; // a0_yyyy
      current_out[215] = ax_xxz + AB[0] * a0_xxz; // a0_xxxz
      current_out[216] = ax_xyz + AB[0] * a0_xyz; // a0_xxyz
      current_out[217] = ax_yyz + AB[0] * a0_yyz; // a0_xyyz
      current_out[218] = ay_yyz + AB[1] * a0_yyz; // a0_yyyz
      current_out[219] = ax_xzz + AB[0] * a0_xzz; // a0_xxzz
      current_out[220] = ax_yzz + AB[0] * a0_yzz; // a0_xyzz
      current_out[221] = ay_yzz + AB[1] * a0_yzz; // a0_yyzz
      current_out[222] = ax_zzz + AB[0] * a0_zzz; // a0_xzzz
      current_out[223] = ay_zzz + AB[1] * a0_zzz; // a0_yzzz
      current_out[224] = az_zzz + AB[2] * a0_zzz; // a0_zzzz

    }
  }
}

