//
// Newint - Parallel electron correlation program.
// Filename: _hrr_70_52.cc
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

void HRRList::perform_HRR_70_52(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    const double* current_data = &data_start[c * 85];
    double* current_out = &data_out[c * 126];
   {
     //current index a: xxxxx
      const double a0_0 = current_data[0];
      const double ax_0 = current_data[21];
      const double ay_0 = current_data[22];
      const double az_0 = current_data[28];
      const double axx_0 = current_data[49];
      const double axy_0 = current_data[50];
      const double ayy_0 = current_data[51];
      const double axz_0 = current_data[57];
      const double ayz_0 = current_data[58];
      const double azz_0 = current_data[64];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[0] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[1] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[2] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[3] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[4] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[5] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxy
      const double a0_0 = current_data[1];
      const double ax_0 = current_data[22];
      const double ay_0 = current_data[23];
      const double az_0 = current_data[29];
      const double axx_0 = current_data[50];
      const double axy_0 = current_data[51];
      const double ayy_0 = current_data[52];
      const double axz_0 = current_data[58];
      const double ayz_0 = current_data[59];
      const double azz_0 = current_data[65];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[6] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[7] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[8] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[9] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[10] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[11] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxyy
      const double a0_0 = current_data[2];
      const double ax_0 = current_data[23];
      const double ay_0 = current_data[24];
      const double az_0 = current_data[30];
      const double axx_0 = current_data[51];
      const double axy_0 = current_data[52];
      const double ayy_0 = current_data[53];
      const double axz_0 = current_data[59];
      const double ayz_0 = current_data[60];
      const double azz_0 = current_data[66];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[12] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[13] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[14] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[15] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[16] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[17] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyyy
      const double a0_0 = current_data[3];
      const double ax_0 = current_data[24];
      const double ay_0 = current_data[25];
      const double az_0 = current_data[31];
      const double axx_0 = current_data[52];
      const double axy_0 = current_data[53];
      const double ayy_0 = current_data[54];
      const double axz_0 = current_data[60];
      const double ayz_0 = current_data[61];
      const double azz_0 = current_data[67];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[18] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[19] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[20] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[21] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[22] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[23] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyyy
      const double a0_0 = current_data[4];
      const double ax_0 = current_data[25];
      const double ay_0 = current_data[26];
      const double az_0 = current_data[32];
      const double axx_0 = current_data[53];
      const double axy_0 = current_data[54];
      const double ayy_0 = current_data[55];
      const double axz_0 = current_data[61];
      const double ayz_0 = current_data[62];
      const double azz_0 = current_data[68];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[24] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[25] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[26] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[27] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[28] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[29] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyyy
      const double a0_0 = current_data[5];
      const double ax_0 = current_data[26];
      const double ay_0 = current_data[27];
      const double az_0 = current_data[33];
      const double axx_0 = current_data[54];
      const double axy_0 = current_data[55];
      const double ayy_0 = current_data[56];
      const double axz_0 = current_data[62];
      const double ayz_0 = current_data[63];
      const double azz_0 = current_data[69];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[30] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[31] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[32] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[33] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[34] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[35] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxz
      const double a0_0 = current_data[6];
      const double ax_0 = current_data[28];
      const double ay_0 = current_data[29];
      const double az_0 = current_data[34];
      const double axx_0 = current_data[57];
      const double axy_0 = current_data[58];
      const double ayy_0 = current_data[59];
      const double axz_0 = current_data[64];
      const double ayz_0 = current_data[65];
      const double azz_0 = current_data[70];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[36] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[37] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[38] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[39] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[40] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[41] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxyz
      const double a0_0 = current_data[7];
      const double ax_0 = current_data[29];
      const double ay_0 = current_data[30];
      const double az_0 = current_data[35];
      const double axx_0 = current_data[58];
      const double axy_0 = current_data[59];
      const double ayy_0 = current_data[60];
      const double axz_0 = current_data[65];
      const double ayz_0 = current_data[66];
      const double azz_0 = current_data[71];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[42] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[43] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[44] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[45] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[46] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[47] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyyz
      const double a0_0 = current_data[8];
      const double ax_0 = current_data[30];
      const double ay_0 = current_data[31];
      const double az_0 = current_data[36];
      const double axx_0 = current_data[59];
      const double axy_0 = current_data[60];
      const double ayy_0 = current_data[61];
      const double axz_0 = current_data[66];
      const double ayz_0 = current_data[67];
      const double azz_0 = current_data[72];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[48] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[49] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[50] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[51] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[52] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[53] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyyz
      const double a0_0 = current_data[9];
      const double ax_0 = current_data[31];
      const double ay_0 = current_data[32];
      const double az_0 = current_data[37];
      const double axx_0 = current_data[60];
      const double axy_0 = current_data[61];
      const double ayy_0 = current_data[62];
      const double axz_0 = current_data[67];
      const double ayz_0 = current_data[68];
      const double azz_0 = current_data[73];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[54] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[55] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[56] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[57] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[58] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[59] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyyz
      const double a0_0 = current_data[10];
      const double ax_0 = current_data[32];
      const double ay_0 = current_data[33];
      const double az_0 = current_data[38];
      const double axx_0 = current_data[61];
      const double axy_0 = current_data[62];
      const double ayy_0 = current_data[63];
      const double axz_0 = current_data[68];
      const double ayz_0 = current_data[69];
      const double azz_0 = current_data[74];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[60] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[61] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[62] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[63] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[64] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[65] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxzz
      const double a0_0 = current_data[11];
      const double ax_0 = current_data[34];
      const double ay_0 = current_data[35];
      const double az_0 = current_data[39];
      const double axx_0 = current_data[64];
      const double axy_0 = current_data[65];
      const double ayy_0 = current_data[66];
      const double axz_0 = current_data[70];
      const double ayz_0 = current_data[71];
      const double azz_0 = current_data[75];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[66] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[67] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[68] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[69] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[70] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[71] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyzz
      const double a0_0 = current_data[12];
      const double ax_0 = current_data[35];
      const double ay_0 = current_data[36];
      const double az_0 = current_data[40];
      const double axx_0 = current_data[65];
      const double axy_0 = current_data[66];
      const double ayy_0 = current_data[67];
      const double axz_0 = current_data[71];
      const double ayz_0 = current_data[72];
      const double azz_0 = current_data[76];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[72] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[73] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[74] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[75] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[76] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[77] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyyzz
      const double a0_0 = current_data[13];
      const double ax_0 = current_data[36];
      const double ay_0 = current_data[37];
      const double az_0 = current_data[41];
      const double axx_0 = current_data[66];
      const double axy_0 = current_data[67];
      const double ayy_0 = current_data[68];
      const double axz_0 = current_data[72];
      const double ayz_0 = current_data[73];
      const double azz_0 = current_data[77];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[78] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[79] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[80] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[81] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[82] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[83] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyyzz
      const double a0_0 = current_data[14];
      const double ax_0 = current_data[37];
      const double ay_0 = current_data[38];
      const double az_0 = current_data[42];
      const double axx_0 = current_data[67];
      const double axy_0 = current_data[68];
      const double ayy_0 = current_data[69];
      const double axz_0 = current_data[73];
      const double ayz_0 = current_data[74];
      const double azz_0 = current_data[78];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[84] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[85] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[86] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[87] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[88] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[89] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxzzz
      const double a0_0 = current_data[15];
      const double ax_0 = current_data[39];
      const double ay_0 = current_data[40];
      const double az_0 = current_data[43];
      const double axx_0 = current_data[70];
      const double axy_0 = current_data[71];
      const double ayy_0 = current_data[72];
      const double axz_0 = current_data[75];
      const double ayz_0 = current_data[76];
      const double azz_0 = current_data[79];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[90] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[91] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[92] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[93] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[94] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[95] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyzzz
      const double a0_0 = current_data[16];
      const double ax_0 = current_data[40];
      const double ay_0 = current_data[41];
      const double az_0 = current_data[44];
      const double axx_0 = current_data[71];
      const double axy_0 = current_data[72];
      const double ayy_0 = current_data[73];
      const double axz_0 = current_data[76];
      const double ayz_0 = current_data[77];
      const double azz_0 = current_data[80];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[96] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[97] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[98] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[99] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[100] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[101] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyzzz
      const double a0_0 = current_data[17];
      const double ax_0 = current_data[41];
      const double ay_0 = current_data[42];
      const double az_0 = current_data[45];
      const double axx_0 = current_data[72];
      const double axy_0 = current_data[73];
      const double ayy_0 = current_data[74];
      const double axz_0 = current_data[77];
      const double ayz_0 = current_data[78];
      const double azz_0 = current_data[81];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[102] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[103] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[104] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[105] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[106] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[107] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xzzzz
      const double a0_0 = current_data[18];
      const double ax_0 = current_data[43];
      const double ay_0 = current_data[44];
      const double az_0 = current_data[46];
      const double axx_0 = current_data[75];
      const double axy_0 = current_data[76];
      const double ayy_0 = current_data[77];
      const double axz_0 = current_data[79];
      const double ayz_0 = current_data[80];
      const double azz_0 = current_data[82];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[108] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[109] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[110] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[111] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[112] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[113] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yzzzz
      const double a0_0 = current_data[19];
      const double ax_0 = current_data[44];
      const double ay_0 = current_data[45];
      const double az_0 = current_data[47];
      const double axx_0 = current_data[76];
      const double axy_0 = current_data[77];
      const double ayy_0 = current_data[78];
      const double axz_0 = current_data[80];
      const double ayz_0 = current_data[81];
      const double azz_0 = current_data[83];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[114] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[115] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[116] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[117] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[118] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[119] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: zzzzz
      const double a0_0 = current_data[20];
      const double ax_0 = current_data[46];
      const double ay_0 = current_data[47];
      const double az_0 = current_data[48];
      const double axx_0 = current_data[79];
      const double axy_0 = current_data[80];
      const double ayy_0 = current_data[81];
      const double axz_0 = current_data[82];
      const double ayz_0 = current_data[83];
      const double azz_0 = current_data[84];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[120] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[121] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[122] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[123] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[124] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[125] = az_z + AB[2] * a0_z; // a0_zz

    }
  }
}

