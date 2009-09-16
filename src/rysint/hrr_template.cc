//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include <src/rysint/eribatch.h>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace std;

void ERIBatch::perform_HRR(const int nloop, const double* data0_start, const double* mapping, const double* AB, double* data0_out) {
  // (ab): (40) -> (22)

  // loop for different contracted integral

  int offset = 0;
  for (int c = 0; c != nloop; ++c) {
    const double* current_data = &data0_start[c * 31];
    double* current_out = &data0_out[c * 36];

    // loop of a
    for (int jz = 0; jz <= 2; ++jz) {
      const int rjz = 5 * jz;
      for (int jy = 0; jy <= 2 - jz; ++jy) {
        const int jx = 2 - jy - jz;
        if (jx >= 0) {
          current_out += 6;

          const int place = jx + 5 * (jy + rjz);
          const int ja0_0 = mapping[place];
          const int jax_0 = mapping[place + 1];
          const int jay_0 = mapping[place + 5];
          const int jaz_0 = mapping[place + 5 * 5];
          const int jaxx_0 = mapping[place + 2];
          const int jaxy_0 = mapping[place + 1 + 5];
          const int jaxz_0 = mapping[place + 1 + 5 * 5 ];
          const int jayy_0 = mapping[place + 5 * 2];
          const int jayz_0 = mapping[place + 5 + 5 * 5];
          const int jazz_0 = mapping[place + 5 * 5 * 2];

          const double a0_0 = current_data[ja0_0];
          const double ax_0 = current_data[jax_0];
          const double ay_0 = current_data[jay_0];
          const double az_0 = current_data[jaz_0];
          const double axx_0 = current_data[jaxx_0];
          const double axy_0 = current_data[jaxy_0];
          const double axz_0 = current_data[jaxz_0];
          const double ayy_0 = current_data[jayy_0];
          const double ayz_0 = current_data[jayz_0];
          const double azz_0 = current_data[jazz_0];

          const double a0_x  = ax_0 + AB[0] * a0_0;
          const double a0_y  = ay_0 + AB[1] * a0_0;
          const double a0_z  = az_0 + AB[2] * a0_0;

          const double ax_x = axx_0 + AB[0] * ax_0;
          const double ay_x = axy_0 + AB[0] * ay_0;
          const double az_x = axz_0 + AB[0] * az_0;
          const double ay_y = ayy_0 + AB[1] * ay_0;
          const double az_y = ayz_0 + AB[1] * az_0;
          const double az_z = azz_0 + AB[2] * az_0;

          const double a0_xx = ax_x + AB[0] * a0_x;
          const double a0_xy = ay_x + AB[1] * a0_x;
          const double a0_xz = az_x + AB[2] * a0_x;
          const double a0_yy = ay_y + AB[1] * a0_y;
          const double a0_yz = az_y + AB[2] * a0_y;
          const double a0_zz = az_z + AB[2] * a0_z;
 
// *********
if(jy == 1 && jx == 1) cout  <<  fixed << setprecision(10) << a0_xy << endl;

          current_out[0] = a0_xx;
          current_out[1] = a0_xy;
          current_out[2] = a0_xz;
          current_out[3] = a0_yy;
          current_out[4] = a0_yz;
          current_out[5] = a0_zz;

        }
      }
    }
  }
}

