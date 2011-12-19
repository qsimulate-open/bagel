//
// Author: Toru Shiozaki
// Machine Generated Code
//

#include "hrrlist.h"
#include <algorithm>

using namespace std;

// modified by hand

void HRRList::perform_HRR_20_11(const int nloop, const double* data_start, const double* AB, double* data_out) {
  const double* current_data = data_start;
  double* current_out  = data_out;
  for (int c = 0; c != nloop; ++c, current_data += 9, current_out += 9) {
    current_out[0] = current_data[3] + AB[0] * current_data[0]; // a0_x
    current_out[1] = current_data[4] + AB[1] * current_data[0]; // a0_y
    current_out[2] = current_data[6] + AB[2] * current_data[0]; // a0_z
    current_out[3] = current_data[4] + AB[0] * current_data[1]; // a0_x
    current_out[4] = current_data[5] + AB[1] * current_data[1]; // a0_y
    current_out[5] = current_data[7] + AB[2] * current_data[1]; // a0_z
    current_out[6] = current_data[6] + AB[0] * current_data[2]; // a0_x
    current_out[7] = current_data[7] + AB[1] * current_data[2]; // a0_y
    current_out[8] = current_data[8] + AB[2] * current_data[2]; // a0_z
  }
}

