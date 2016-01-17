//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _hrr_40_22.cc
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

void HRRList::perform_HRR_40_22(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 31];
    auto current_out = &data_out[c * 36];
   {
     //current index a: xx
      const auto a0_0 = current_data[0];
      const auto ax_0 = current_data[6];
      const auto ay_0 = current_data[7];
      const auto az_0 = current_data[10];
      const auto axx_0 = current_data[16];
      const auto axy_0 = current_data[17];
      const auto ayy_0 = current_data[18];
      const auto axz_0 = current_data[21];
      const auto ayz_0 = current_data[22];
      const auto azz_0 = current_data[25];

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
     //current index a: xy
      const auto a0_0 = current_data[1];
      const auto ax_0 = current_data[7];
      const auto ay_0 = current_data[8];
      const auto az_0 = current_data[11];
      const auto axx_0 = current_data[17];
      const auto axy_0 = current_data[18];
      const auto ayy_0 = current_data[19];
      const auto axz_0 = current_data[22];
      const auto ayz_0 = current_data[23];
      const auto azz_0 = current_data[26];

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
     //current index a: yy
      const auto a0_0 = current_data[2];
      const auto ax_0 = current_data[8];
      const auto ay_0 = current_data[9];
      const auto az_0 = current_data[12];
      const auto axx_0 = current_data[18];
      const auto axy_0 = current_data[19];
      const auto ayy_0 = current_data[20];
      const auto axz_0 = current_data[23];
      const auto ayz_0 = current_data[24];
      const auto azz_0 = current_data[27];

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
     //current index a: xz
      const auto a0_0 = current_data[3];
      const auto ax_0 = current_data[10];
      const auto ay_0 = current_data[11];
      const auto az_0 = current_data[13];
      const auto axx_0 = current_data[21];
      const auto axy_0 = current_data[22];
      const auto ayy_0 = current_data[23];
      const auto axz_0 = current_data[25];
      const auto ayz_0 = current_data[26];
      const auto azz_0 = current_data[28];

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
     //current index a: yz
      const auto a0_0 = current_data[4];
      const auto ax_0 = current_data[11];
      const auto ay_0 = current_data[12];
      const auto az_0 = current_data[14];
      const auto axx_0 = current_data[22];
      const auto axy_0 = current_data[23];
      const auto ayy_0 = current_data[24];
      const auto axz_0 = current_data[26];
      const auto ayz_0 = current_data[27];
      const auto azz_0 = current_data[29];

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
     //current index a: zz
      const auto a0_0 = current_data[5];
      const auto ax_0 = current_data[13];
      const auto ay_0 = current_data[14];
      const auto az_0 = current_data[15];
      const auto axx_0 = current_data[25];
      const auto axy_0 = current_data[26];
      const auto ayy_0 = current_data[27];
      const auto axz_0 = current_data[28];
      const auto ayz_0 = current_data[29];
      const auto azz_0 = current_data[30];

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
  }
}


void CHRRList::perform_HRR_40_22(const int nloop, const complex<double>* data_start, const array<double,3>& AB, complex<double>* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 31];
    auto current_out = &data_out[c * 36];
   {
     //current index a: xx
      const auto a0_0 = current_data[0];
      const auto ax_0 = current_data[6];
      const auto ay_0 = current_data[7];
      const auto az_0 = current_data[10];
      const auto axx_0 = current_data[16];
      const auto axy_0 = current_data[17];
      const auto ayy_0 = current_data[18];
      const auto axz_0 = current_data[21];
      const auto ayz_0 = current_data[22];
      const auto azz_0 = current_data[25];

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
     //current index a: xy
      const auto a0_0 = current_data[1];
      const auto ax_0 = current_data[7];
      const auto ay_0 = current_data[8];
      const auto az_0 = current_data[11];
      const auto axx_0 = current_data[17];
      const auto axy_0 = current_data[18];
      const auto ayy_0 = current_data[19];
      const auto axz_0 = current_data[22];
      const auto ayz_0 = current_data[23];
      const auto azz_0 = current_data[26];

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
     //current index a: yy
      const auto a0_0 = current_data[2];
      const auto ax_0 = current_data[8];
      const auto ay_0 = current_data[9];
      const auto az_0 = current_data[12];
      const auto axx_0 = current_data[18];
      const auto axy_0 = current_data[19];
      const auto ayy_0 = current_data[20];
      const auto axz_0 = current_data[23];
      const auto ayz_0 = current_data[24];
      const auto azz_0 = current_data[27];

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
     //current index a: xz
      const auto a0_0 = current_data[3];
      const auto ax_0 = current_data[10];
      const auto ay_0 = current_data[11];
      const auto az_0 = current_data[13];
      const auto axx_0 = current_data[21];
      const auto axy_0 = current_data[22];
      const auto ayy_0 = current_data[23];
      const auto axz_0 = current_data[25];
      const auto ayz_0 = current_data[26];
      const auto azz_0 = current_data[28];

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
     //current index a: yz
      const auto a0_0 = current_data[4];
      const auto ax_0 = current_data[11];
      const auto ay_0 = current_data[12];
      const auto az_0 = current_data[14];
      const auto axx_0 = current_data[22];
      const auto axy_0 = current_data[23];
      const auto ayy_0 = current_data[24];
      const auto axz_0 = current_data[26];
      const auto ayz_0 = current_data[27];
      const auto azz_0 = current_data[29];

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
     //current index a: zz
      const auto a0_0 = current_data[5];
      const auto ax_0 = current_data[13];
      const auto ay_0 = current_data[14];
      const auto az_0 = current_data[15];
      const auto axx_0 = current_data[25];
      const auto axy_0 = current_data[26];
      const auto ayy_0 = current_data[27];
      const auto axz_0 = current_data[28];
      const auto ayz_0 = current_data[29];
      const auto azz_0 = current_data[30];

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
  }
}

