//
// Newint - Parallel electron correlation program.
// Filename: _hrr_70_43.cc
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

void HRRList::perform_HRR_70_43(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    const double* current_data = &data_start[c * 100];
    double* current_out = &data_out[c * 150];
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
  }
}

