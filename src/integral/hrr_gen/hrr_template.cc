//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: hrr_template.cc
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


#include "../rys/eribatch.h"
#include <algorithm>

using namespace std;

void ERIBatch::perform_HRR_40_22(const int nloop, const double* data_start, const int size) {
  // (ab): (40) -> (22)

  // loop for different contracted integral

  int offset = 0;
  for (int c = 0; c != nloop; ++c, offset += size) {
    const double* current_data = &data_start[offset];

    // loop of a
    for (int jz = 0; jz <= 2; ++jz) {
      const int rjz = 5 * jz;
      for (int jy = 0; jy <= 2 - jz; ++jy) {
        const int jx = 2 - jy - jz;
        if (jx >= 0) {
          const int place = jx + 5 * (jy + rjz);
          const int ja_0 = amapping_[place];
          const int jax_0 = amapping_[place + 1];
          const int jay_0 = amapping_[place + 5];
          const int jaz_0 = amapping_[place + 5 * 5];
          const int jaxx_0 = amapping_[place + 2];
          const int jaxy_0 = amapping_[place + 1 + 5];
          const int jaxz_0 = amapping_[place + 1 + 5 * 5 ];
          const int jayy_0 = amapping_[place + 5 * 2];
          const int jayz_0 = amapping_[place + 5 + 5 * 5];
          const int jazz_0 = amapping_[place + 5 * 5 * 2];

          const double a_0 = current_data[ja_0];
          const double ax_0 = current_data[jax_0];
          const double ay_0 = current_data[jay_0];
          const double az_0 = current_data[jaz_0];
          const double axx_0 = current_data[jaxx_0];
          const double axy_0 = current_data[jaxy_0];
          const double axz_0 = current_data[jaxz_0];
          const double ayy_0 = current_data[jayy_0];
          const double ayz_0 = current_data[jayz_0];
          const double azz_0 = current_data[jazz_0];

          const double a_x  = ax_0  + AB_[0] * a_0;
          const double a_y  = ay_0  + AB_[1] * a_0;
          const double a_z  = az_0  + AB_[2] * a_0;

          const double ax_x = axx_0 + AB_[0] * ax_0;
          const double ay_x = axy_0 + AB_[0] * ay_0;
          const double az_x = axz_0 + AB_[0] * az_0;
          const double ay_y = ayy_0 + AB_[1] * ay_0;
          const double az_y = ayz_0 + AB_[1] * az_0;
          const double az_z = azz_0 + AB_[2] * az_0;

          const double a_xx = ax_x  + AB_[0] * a_x;
          const double a_xy = ay_x  + AB_[1] * a_x;
          const double a_xz = az_x  + AB_[2] * a_x;
          const double a_yy = ay_y  + AB_[1] * a_y;
          const double a_yz = az_y  + AB_[2] * a_y;
          const double a_zz = az_z  + AB_[2] * a_z;

        }
      }
    }
  }
}

