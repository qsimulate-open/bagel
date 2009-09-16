//
// Author: Toru Shiozaki
// Machine Generated Code in NewInt
//

#include "vrrlist.h"

// returns double array of length 15
void VRRList::_vrr_4000(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;

  data_[3] = C00_[0];
  data_[4] = C00_[1];
  data_[5] = C00_[2];

  double B10_current[3];
  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];

  data_[6] = C00_[0] * data_[3] + B10_current[0];
  data_[7] = C00_[1] * data_[4] + B10_current[1];
  data_[8] = C00_[2] * data_[5] + B10_current[2];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];

  data_[9] = C00_[0] * data_[6] + B10_current[0] * data_[3];
  data_[10] = C00_[1] * data_[7] + B10_current[1] * data_[4];
  data_[11] = C00_[2] * data_[8] + B10_current[2] * data_[5];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];

  data_[12] = C00_[0] * data_[9] + B10_current[0] * data_[6];
  data_[13] = C00_[1] * data_[10] + B10_current[1] * data_[7];
  data_[14] = C00_[2] * data_[11] + B10_current[2] * data_[8];
}

