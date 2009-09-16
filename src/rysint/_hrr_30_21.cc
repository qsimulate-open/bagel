//
// Author: Toru Shiozaki
// Machine Generated Code
//

#include "hrrlist.h"
#include <algorithm>

using namespace std;

void HRRList::perform_HRR_30_21(const int nloop, const double* data_start, const double* AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    const double* current_data = &data_start[c * 16];
    double* current_out = &data_out[c * 18];
   {
     //current index a: xx
      const double a0_0 = current_data[0];
      const double ax_0 = current_data[6];
      const double ay_0 = current_data[7];
      const double az_0 = current_data[10];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xy
      const double a0_0 = current_data[1];
      const double ax_0 = current_data[7];
      const double ay_0 = current_data[8];
      const double az_0 = current_data[11];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yy
      const double a0_0 = current_data[2];
      const double ax_0 = current_data[8];
      const double ay_0 = current_data[9];
      const double az_0 = current_data[12];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: xz
      const double a0_0 = current_data[3];
      const double ax_0 = current_data[10];
      const double ay_0 = current_data[11];
      const double az_0 = current_data[13];

      current_out[9] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[10] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[11] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: yz
      const double a0_0 = current_data[4];
      const double ax_0 = current_data[11];
      const double ay_0 = current_data[12];
      const double az_0 = current_data[14];

      current_out[12] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[13] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[14] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: zz
      const double a0_0 = current_data[5];
      const double ax_0 = current_data[13];
      const double ay_0 = current_data[14];
      const double az_0 = current_data[15];

      current_out[15] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[16] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[17] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}

