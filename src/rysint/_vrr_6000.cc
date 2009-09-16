//
// Author: Toru Shiozaki
// Machine Generated Code in NewInt
//

#include "vrrlist.h"

// returns double array of length 28
void VRRList::_vrr_6000(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;
  data_[3] = 1.0;

  data_[4] = C00_[0];
  data_[5] = C00_[1];
  data_[6] = C00_[2];
  data_[7] = C00_[3];

  double B10_current[4];
  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];

  data_[8] = C00_[0] * data_[4] + B10_current[0];
  data_[9] = C00_[1] * data_[5] + B10_current[1];
  data_[10] = C00_[2] * data_[6] + B10_current[2];
  data_[11] = C00_[3] * data_[7] + B10_current[3];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];

  data_[12] = C00_[0] * data_[8] + B10_current[0] * data_[4];
  data_[13] = C00_[1] * data_[9] + B10_current[1] * data_[5];
  data_[14] = C00_[2] * data_[10] + B10_current[2] * data_[6];
  data_[15] = C00_[3] * data_[11] + B10_current[3] * data_[7];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];

  data_[16] = C00_[0] * data_[12] + B10_current[0] * data_[8];
  data_[17] = C00_[1] * data_[13] + B10_current[1] * data_[9];
  data_[18] = C00_[2] * data_[14] + B10_current[2] * data_[10];
  data_[19] = C00_[3] * data_[15] + B10_current[3] * data_[11];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];

  data_[20] = C00_[0] * data_[16] + B10_current[0] * data_[12];
  data_[21] = C00_[1] * data_[17] + B10_current[1] * data_[13];
  data_[22] = C00_[2] * data_[18] + B10_current[2] * data_[14];
  data_[23] = C00_[3] * data_[19] + B10_current[3] * data_[15];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];

  data_[24] = C00_[0] * data_[20] + B10_current[0] * data_[16];
  data_[25] = C00_[1] * data_[21] + B10_current[1] * data_[17];
  data_[26] = C00_[2] * data_[22] + B10_current[2] * data_[18];
  data_[27] = C00_[3] * data_[23] + B10_current[3] * data_[19];
}

