//
// Author: Toru Shiozaki
// Machine Generated Code
//

#include "hrrlist.h"
#include <algorithm>

using namespace std;

void HRRList::perform_HRR_20_11(const int nloop, const double* data_start, const double* AB, double* data_out) {
  for (int c = 0; c != nloop; ++c) {
    const double* current_data = &data_start[c * 9];
    double* current_out = &data_out[c * 9];
   {
     //current index a: x
      const double a0_0 = current_data[0];
      const double ax_0 = current_data[3];
      const double ay_0 = current_data[4];
      const double az_0 = current_data[6];

      current_out[0] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[1] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[2] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: y
      const double a0_0 = current_data[1];
      const double ax_0 = current_data[4];
      const double ay_0 = current_data[5];
      const double az_0 = current_data[7];

      current_out[3] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[4] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[5] = az_0 + AB[2] * a0_0; // a0_z

    }
   {
     //current index a: z
      const double a0_0 = current_data[2];
      const double ax_0 = current_data[6];
      const double ay_0 = current_data[7];
      const double az_0 = current_data[8];

      current_out[6] = ax_0 + AB[0] * a0_0; // a0_x
      current_out[7] = ay_0 + AB[1] * a0_0; // a0_y
      current_out[8] = az_0 + AB[2] * a0_0; // a0_z

    }
  }
}

