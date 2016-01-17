//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _hrr_40_31.cc
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

void HRRList::perform_HRR_40_31(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 25];
    auto current_out = &data_out[c * 30];
   {
     //current index a: xxx
      const auto a0_0 = current_data[0];
      const auto ax_0 = current_data[10];
      const auto ay_0 = current_data[11];
      const auto az_0 = current_data[15];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxy
      const auto a0_0 = current_data[1];
      const auto ax_0 = current_data[11];
      const auto ay_0 = current_data[12];
      const auto az_0 = current_data[16];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyy
      const auto a0_0 = current_data[2];
      const auto ax_0 = current_data[12];
      const auto ay_0 = current_data[13];
      const auto az_0 = current_data[17];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyy
      const auto a0_0 = current_data[3];
      const auto ax_0 = current_data[13];
      const auto ay_0 = current_data[14];
      const auto az_0 = current_data[18];

      current_out[9] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[10] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[11] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxz
      const auto a0_0 = current_data[4];
      const auto ax_0 = current_data[15];
      const auto ay_0 = current_data[16];
      const auto az_0 = current_data[19];

      current_out[12] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[13] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[14] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyz
      const auto a0_0 = current_data[5];
      const auto ax_0 = current_data[16];
      const auto ay_0 = current_data[17];
      const auto az_0 = current_data[20];

      current_out[15] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[16] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[17] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyz
      const auto a0_0 = current_data[6];
      const auto ax_0 = current_data[17];
      const auto ay_0 = current_data[18];
      const auto az_0 = current_data[21];

      current_out[18] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[19] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[20] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xzz
      const auto a0_0 = current_data[7];
      const auto ax_0 = current_data[19];
      const auto ay_0 = current_data[20];
      const auto az_0 = current_data[22];

      current_out[21] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[22] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[23] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yzz
      const auto a0_0 = current_data[8];
      const auto ax_0 = current_data[20];
      const auto ay_0 = current_data[21];
      const auto az_0 = current_data[23];

      current_out[24] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[25] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[26] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: zzz
      const auto a0_0 = current_data[9];
      const auto ax_0 = current_data[22];
      const auto ay_0 = current_data[23];
      const auto az_0 = current_data[24];

      current_out[27] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[28] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[29] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}


void CHRRList::perform_HRR_40_31(const int nloop, const complex<double>* data_start, const array<double,3>& AB, complex<double>* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 25];
    auto current_out = &data_out[c * 30];
   {
     //current index a: xxx
      const auto a0_0 = current_data[0];
      const auto ax_0 = current_data[10];
      const auto ay_0 = current_data[11];
      const auto az_0 = current_data[15];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxy
      const auto a0_0 = current_data[1];
      const auto ax_0 = current_data[11];
      const auto ay_0 = current_data[12];
      const auto az_0 = current_data[16];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyy
      const auto a0_0 = current_data[2];
      const auto ax_0 = current_data[12];
      const auto ay_0 = current_data[13];
      const auto az_0 = current_data[17];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyy
      const auto a0_0 = current_data[3];
      const auto ax_0 = current_data[13];
      const auto ay_0 = current_data[14];
      const auto az_0 = current_data[18];

      current_out[9] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[10] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[11] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xxz
      const auto a0_0 = current_data[4];
      const auto ax_0 = current_data[15];
      const auto ay_0 = current_data[16];
      const auto az_0 = current_data[19];

      current_out[12] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[13] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[14] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xyz
      const auto a0_0 = current_data[5];
      const auto ax_0 = current_data[16];
      const auto ay_0 = current_data[17];
      const auto az_0 = current_data[20];

      current_out[15] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[16] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[17] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yyz
      const auto a0_0 = current_data[6];
      const auto ax_0 = current_data[17];
      const auto ay_0 = current_data[18];
      const auto az_0 = current_data[21];

      current_out[18] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[19] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[20] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xzz
      const auto a0_0 = current_data[7];
      const auto ax_0 = current_data[19];
      const auto ay_0 = current_data[20];
      const auto az_0 = current_data[22];

      current_out[21] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[22] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[23] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yzz
      const auto a0_0 = current_data[8];
      const auto ax_0 = current_data[20];
      const auto ay_0 = current_data[21];
      const auto az_0 = current_data[23];

      current_out[24] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[25] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[26] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: zzz
      const auto a0_0 = current_data[9];
      const auto ax_0 = current_data[22];
      const auto ay_0 = current_data[23];
      const auto az_0 = current_data[24];

      current_out[27] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[28] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[29] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}

