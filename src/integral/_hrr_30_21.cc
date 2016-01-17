//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _hrr_30_21.cc
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

void HRRList::perform_HRR_30_21(const int nloop, const double* data_start, const array<double,3>& AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 16];
    auto current_out = &data_out[c * 18];
   {
     //current index a: xx
      const auto a0_0 = current_data[0];
      const auto ax_0 = current_data[6];
      const auto ay_0 = current_data[7];
      const auto az_0 = current_data[10];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xy
      const auto a0_0 = current_data[1];
      const auto ax_0 = current_data[7];
      const auto ay_0 = current_data[8];
      const auto az_0 = current_data[11];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yy
      const auto a0_0 = current_data[2];
      const auto ax_0 = current_data[8];
      const auto ay_0 = current_data[9];
      const auto az_0 = current_data[12];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xz
      const auto a0_0 = current_data[3];
      const auto ax_0 = current_data[10];
      const auto ay_0 = current_data[11];
      const auto az_0 = current_data[13];

      current_out[9] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[10] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[11] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yz
      const auto a0_0 = current_data[4];
      const auto ax_0 = current_data[11];
      const auto ay_0 = current_data[12];
      const auto az_0 = current_data[14];

      current_out[12] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[13] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[14] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: zz
      const auto a0_0 = current_data[5];
      const auto ax_0 = current_data[13];
      const auto ay_0 = current_data[14];
      const auto az_0 = current_data[15];

      current_out[15] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[16] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[17] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}


void CHRRList::perform_HRR_30_21(const int nloop, const complex<double>* data_start, const array<double,3>& AB, complex<double>* data_out) {
  for (int c = 0; c != nloop; ++c) {
    auto current_data = &data_start[c * 16];
    auto current_out = &data_out[c * 18];
   {
     //current index a: xx
      const auto a0_0 = current_data[0];
      const auto ax_0 = current_data[6];
      const auto ay_0 = current_data[7];
      const auto az_0 = current_data[10];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xy
      const auto a0_0 = current_data[1];
      const auto ax_0 = current_data[7];
      const auto ay_0 = current_data[8];
      const auto az_0 = current_data[11];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yy
      const auto a0_0 = current_data[2];
      const auto ax_0 = current_data[8];
      const auto ay_0 = current_data[9];
      const auto az_0 = current_data[12];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xz
      const auto a0_0 = current_data[3];
      const auto ax_0 = current_data[10];
      const auto ay_0 = current_data[11];
      const auto az_0 = current_data[13];

      current_out[9] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[10] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[11] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yz
      const auto a0_0 = current_data[4];
      const auto ax_0 = current_data[11];
      const auto ay_0 = current_data[12];
      const auto az_0 = current_data[14];

      current_out[12] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[13] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[14] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: zz
      const auto a0_0 = current_data[5];
      const auto ax_0 = current_data[13];
      const auto ay_0 = current_data[14];
      const auto az_0 = current_data[15];

      current_out[15] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[16] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[17] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}

