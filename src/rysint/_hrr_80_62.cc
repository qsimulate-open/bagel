//
// BAGEL - Parallel electron correlation program.
// Filename: _hrr_80_62.cc
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

void HRRList::perform_HRR_80_62(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    const double* current_data = &data_start[c * 109];
    double* current_out = &data_out[c * 168];
   {
     //current index a: xxxxxx
      const double a0_0 = current_data[0];
      const double ax_0 = current_data[28];
      const double ay_0 = current_data[29];
      const double az_0 = current_data[36];
      const double axx_0 = current_data[64];
      const double axy_0 = current_data[65];
      const double ayy_0 = current_data[66];
      const double axz_0 = current_data[73];
      const double ayz_0 = current_data[74];
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

      current_out[0] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[1] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[2] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[3] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[4] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[5] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxxy
      const double a0_0 = current_data[1];
      const double ax_0 = current_data[29];
      const double ay_0 = current_data[30];
      const double az_0 = current_data[37];
      const double axx_0 = current_data[65];
      const double axy_0 = current_data[66];
      const double ayy_0 = current_data[67];
      const double axz_0 = current_data[74];
      const double ayz_0 = current_data[75];
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

      current_out[6] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[7] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[8] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[9] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[10] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[11] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxxyy
      const double a0_0 = current_data[2];
      const double ax_0 = current_data[30];
      const double ay_0 = current_data[31];
      const double az_0 = current_data[38];
      const double axx_0 = current_data[66];
      const double axy_0 = current_data[67];
      const double ayy_0 = current_data[68];
      const double axz_0 = current_data[75];
      const double ayz_0 = current_data[76];
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

      current_out[12] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[13] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[14] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[15] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[16] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[17] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxxyyy
      const double a0_0 = current_data[3];
      const double ax_0 = current_data[31];
      const double ay_0 = current_data[32];
      const double az_0 = current_data[39];
      const double axx_0 = current_data[67];
      const double axy_0 = current_data[68];
      const double ayy_0 = current_data[69];
      const double axz_0 = current_data[76];
      const double ayz_0 = current_data[77];
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

      current_out[18] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[19] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[20] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[21] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[22] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[23] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxyyyy
      const double a0_0 = current_data[4];
      const double ax_0 = current_data[32];
      const double ay_0 = current_data[33];
      const double az_0 = current_data[40];
      const double axx_0 = current_data[68];
      const double axy_0 = current_data[69];
      const double ayy_0 = current_data[70];
      const double axz_0 = current_data[77];
      const double ayz_0 = current_data[78];
      const double azz_0 = current_data[85];

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
     //current index a: xyyyyy
      const double a0_0 = current_data[5];
      const double ax_0 = current_data[33];
      const double ay_0 = current_data[34];
      const double az_0 = current_data[41];
      const double axx_0 = current_data[69];
      const double axy_0 = current_data[70];
      const double ayy_0 = current_data[71];
      const double axz_0 = current_data[78];
      const double ayz_0 = current_data[79];
      const double azz_0 = current_data[86];

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
     //current index a: yyyyyy
      const double a0_0 = current_data[6];
      const double ax_0 = current_data[34];
      const double ay_0 = current_data[35];
      const double az_0 = current_data[42];
      const double axx_0 = current_data[70];
      const double axy_0 = current_data[71];
      const double ayy_0 = current_data[72];
      const double axz_0 = current_data[79];
      const double ayz_0 = current_data[80];
      const double azz_0 = current_data[87];

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
     //current index a: xxxxxz
      const double a0_0 = current_data[7];
      const double ax_0 = current_data[36];
      const double ay_0 = current_data[37];
      const double az_0 = current_data[43];
      const double axx_0 = current_data[73];
      const double axy_0 = current_data[74];
      const double ayy_0 = current_data[75];
      const double axz_0 = current_data[81];
      const double ayz_0 = current_data[82];
      const double azz_0 = current_data[88];

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
     //current index a: xxxxyz
      const double a0_0 = current_data[8];
      const double ax_0 = current_data[37];
      const double ay_0 = current_data[38];
      const double az_0 = current_data[44];
      const double axx_0 = current_data[74];
      const double axy_0 = current_data[75];
      const double ayy_0 = current_data[76];
      const double axz_0 = current_data[82];
      const double ayz_0 = current_data[83];
      const double azz_0 = current_data[89];

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
     //current index a: xxxyyz
      const double a0_0 = current_data[9];
      const double ax_0 = current_data[38];
      const double ay_0 = current_data[39];
      const double az_0 = current_data[45];
      const double axx_0 = current_data[75];
      const double axy_0 = current_data[76];
      const double ayy_0 = current_data[77];
      const double axz_0 = current_data[83];
      const double ayz_0 = current_data[84];
      const double azz_0 = current_data[90];

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
     //current index a: xxyyyz
      const double a0_0 = current_data[10];
      const double ax_0 = current_data[39];
      const double ay_0 = current_data[40];
      const double az_0 = current_data[46];
      const double axx_0 = current_data[76];
      const double axy_0 = current_data[77];
      const double ayy_0 = current_data[78];
      const double axz_0 = current_data[84];
      const double ayz_0 = current_data[85];
      const double azz_0 = current_data[91];

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
     //current index a: xyyyyz
      const double a0_0 = current_data[11];
      const double ax_0 = current_data[40];
      const double ay_0 = current_data[41];
      const double az_0 = current_data[47];
      const double axx_0 = current_data[77];
      const double axy_0 = current_data[78];
      const double ayy_0 = current_data[79];
      const double axz_0 = current_data[85];
      const double ayz_0 = current_data[86];
      const double azz_0 = current_data[92];

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
     //current index a: yyyyyz
      const double a0_0 = current_data[12];
      const double ax_0 = current_data[41];
      const double ay_0 = current_data[42];
      const double az_0 = current_data[48];
      const double axx_0 = current_data[78];
      const double axy_0 = current_data[79];
      const double ayy_0 = current_data[80];
      const double axz_0 = current_data[86];
      const double ayz_0 = current_data[87];
      const double azz_0 = current_data[93];

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
     //current index a: xxxxzz
      const double a0_0 = current_data[13];
      const double ax_0 = current_data[43];
      const double ay_0 = current_data[44];
      const double az_0 = current_data[49];
      const double axx_0 = current_data[81];
      const double axy_0 = current_data[82];
      const double ayy_0 = current_data[83];
      const double axz_0 = current_data[88];
      const double ayz_0 = current_data[89];
      const double azz_0 = current_data[94];

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
     //current index a: xxxyzz
      const double a0_0 = current_data[14];
      const double ax_0 = current_data[44];
      const double ay_0 = current_data[45];
      const double az_0 = current_data[50];
      const double axx_0 = current_data[82];
      const double axy_0 = current_data[83];
      const double ayy_0 = current_data[84];
      const double axz_0 = current_data[89];
      const double ayz_0 = current_data[90];
      const double azz_0 = current_data[95];

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
     //current index a: xxyyzz
      const double a0_0 = current_data[15];
      const double ax_0 = current_data[45];
      const double ay_0 = current_data[46];
      const double az_0 = current_data[51];
      const double axx_0 = current_data[83];
      const double axy_0 = current_data[84];
      const double ayy_0 = current_data[85];
      const double axz_0 = current_data[90];
      const double ayz_0 = current_data[91];
      const double azz_0 = current_data[96];

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
     //current index a: xyyyzz
      const double a0_0 = current_data[16];
      const double ax_0 = current_data[46];
      const double ay_0 = current_data[47];
      const double az_0 = current_data[52];
      const double axx_0 = current_data[84];
      const double axy_0 = current_data[85];
      const double ayy_0 = current_data[86];
      const double axz_0 = current_data[91];
      const double ayz_0 = current_data[92];
      const double azz_0 = current_data[97];

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
     //current index a: yyyyzz
      const double a0_0 = current_data[17];
      const double ax_0 = current_data[47];
      const double ay_0 = current_data[48];
      const double az_0 = current_data[53];
      const double axx_0 = current_data[85];
      const double axy_0 = current_data[86];
      const double ayy_0 = current_data[87];
      const double axz_0 = current_data[92];
      const double ayz_0 = current_data[93];
      const double azz_0 = current_data[98];

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
     //current index a: xxxzzz
      const double a0_0 = current_data[18];
      const double ax_0 = current_data[49];
      const double ay_0 = current_data[50];
      const double az_0 = current_data[54];
      const double axx_0 = current_data[88];
      const double axy_0 = current_data[89];
      const double ayy_0 = current_data[90];
      const double axz_0 = current_data[94];
      const double ayz_0 = current_data[95];
      const double azz_0 = current_data[99];

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
     //current index a: xxyzzz
      const double a0_0 = current_data[19];
      const double ax_0 = current_data[50];
      const double ay_0 = current_data[51];
      const double az_0 = current_data[55];
      const double axx_0 = current_data[89];
      const double axy_0 = current_data[90];
      const double ayy_0 = current_data[91];
      const double axz_0 = current_data[95];
      const double ayz_0 = current_data[96];
      const double azz_0 = current_data[100];

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
     //current index a: xyyzzz
      const double a0_0 = current_data[20];
      const double ax_0 = current_data[51];
      const double ay_0 = current_data[52];
      const double az_0 = current_data[56];
      const double axx_0 = current_data[90];
      const double axy_0 = current_data[91];
      const double ayy_0 = current_data[92];
      const double axz_0 = current_data[96];
      const double ayz_0 = current_data[97];
      const double azz_0 = current_data[101];

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
   {
     //current index a: yyyzzz
      const double a0_0 = current_data[21];
      const double ax_0 = current_data[52];
      const double ay_0 = current_data[53];
      const double az_0 = current_data[57];
      const double axx_0 = current_data[91];
      const double axy_0 = current_data[92];
      const double ayy_0 = current_data[93];
      const double axz_0 = current_data[97];
      const double ayz_0 = current_data[98];
      const double azz_0 = current_data[102];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[126] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[127] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[128] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[129] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[130] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[131] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xxzzzz
      const double a0_0 = current_data[22];
      const double ax_0 = current_data[54];
      const double ay_0 = current_data[55];
      const double az_0 = current_data[58];
      const double axx_0 = current_data[94];
      const double axy_0 = current_data[95];
      const double ayy_0 = current_data[96];
      const double axz_0 = current_data[99];
      const double ayz_0 = current_data[100];
      const double azz_0 = current_data[103];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[132] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[133] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[134] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[135] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[136] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[137] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xyzzzz
      const double a0_0 = current_data[23];
      const double ax_0 = current_data[55];
      const double ay_0 = current_data[56];
      const double az_0 = current_data[59];
      const double axx_0 = current_data[95];
      const double axy_0 = current_data[96];
      const double ayy_0 = current_data[97];
      const double axz_0 = current_data[100];
      const double ayz_0 = current_data[101];
      const double azz_0 = current_data[104];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[138] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[139] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[140] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[141] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[142] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[143] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yyzzzz
      const double a0_0 = current_data[24];
      const double ax_0 = current_data[56];
      const double ay_0 = current_data[57];
      const double az_0 = current_data[60];
      const double axx_0 = current_data[96];
      const double axy_0 = current_data[97];
      const double ayy_0 = current_data[98];
      const double axz_0 = current_data[101];
      const double ayz_0 = current_data[102];
      const double azz_0 = current_data[105];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[144] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[145] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[146] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[147] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[148] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[149] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: xzzzzz
      const double a0_0 = current_data[25];
      const double ax_0 = current_data[58];
      const double ay_0 = current_data[59];
      const double az_0 = current_data[61];
      const double axx_0 = current_data[99];
      const double axy_0 = current_data[100];
      const double ayy_0 = current_data[101];
      const double axz_0 = current_data[103];
      const double ayz_0 = current_data[104];
      const double azz_0 = current_data[106];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[150] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[151] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[152] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[153] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[154] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[155] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: yzzzzz
      const double a0_0 = current_data[26];
      const double ax_0 = current_data[59];
      const double ay_0 = current_data[60];
      const double az_0 = current_data[62];
      const double axx_0 = current_data[100];
      const double axy_0 = current_data[101];
      const double ayy_0 = current_data[102];
      const double axz_0 = current_data[104];
      const double ayz_0 = current_data[105];
      const double azz_0 = current_data[107];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[156] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[157] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[158] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[159] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[160] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[161] = az_z + AB[2] * a0_z; // a0_zz

    }
   {
     //current index a: zzzzzz
      const double a0_0 = current_data[27];
      const double ax_0 = current_data[61];
      const double ay_0 = current_data[62];
      const double az_0 = current_data[63];
      const double axx_0 = current_data[103];
      const double axy_0 = current_data[104];
      const double ayy_0 = current_data[105];
      const double axz_0 = current_data[106];
      const double ayz_0 = current_data[107];
      const double azz_0 = current_data[108];

      const double a0_x = ax_0 + AB[0] * a0_0;
      const double a0_y = ay_0 + AB[1] * a0_0;
      const double a0_z = az_0 + AB[2] * a0_0;

      const double ax_x = axx_0 + AB[0] * ax_0;
      const double ax_y = axy_0 + AB[1] * ax_0;
      const double ay_y = ayy_0 + AB[1] * ay_0;
      const double ax_z = axz_0 + AB[2] * ax_0;
      const double ay_z = ayz_0 + AB[2] * ay_0;
      const double az_z = azz_0 + AB[2] * az_0;

      current_out[162] = ax_x + AB[0] * a0_x; // a0_xx
      current_out[163] = ax_y + AB[0] * a0_y; // a0_xy
      current_out[164] = ay_y + AB[1] * a0_y; // a0_yy
      current_out[165] = ax_z + AB[0] * a0_z; // a0_xz
      current_out[166] = ay_z + AB[1] * a0_z; // a0_yz
      current_out[167] = az_z + AB[2] * a0_z; // a0_zz

    }
  }
}

