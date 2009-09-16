//
// Author: Toru Shiozaki
// Machine Generated Code in NewInt
//

#include "vrrlist.h"

// returns double array of length 36
void VRRList::_vrr_2030(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;

  data_[3] = C00_[0];
  data_[4] = C00_[1];
  data_[5] = C00_[2];

  data_[6] = C00_[0] * data_[3] + B10_[0];
  data_[7] = C00_[1] * data_[4] + B10_[1];
  data_[8] = C00_[2] * data_[5] + B10_[2];

  data_[9] = D00_[0];
  data_[10] = D00_[1];
  data_[11] = D00_[2];

  double cB00_current[3];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];

  data_[12] = C00_[0] * data_[9] + cB00_current[0];
  data_[13] = C00_[1] * data_[10] + cB00_current[1];
  data_[14] = C00_[2] * data_[11] + cB00_current[2];

  data_[15] = C00_[0] * data_[12] + B10_[0] * data_[9] + cB00_current[0] * data_[3];
  data_[16] = C00_[1] * data_[13] + B10_[1] * data_[10] + cB00_current[1] * data_[4];
  data_[17] = C00_[2] * data_[14] + B10_[2] * data_[11] + cB00_current[2] * data_[5];

  double B01_current[3];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];

  data_[18] = D00_[0] * data_[9] + B01_current[0];
  data_[19] = D00_[1] * data_[10] + B01_current[1];
  data_[20] = D00_[2] * data_[11] + B01_current[2];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];

  data_[21] = C00_[0] * data_[18] + cB00_current[0] * data_[9];
  data_[22] = C00_[1] * data_[19] + cB00_current[1] * data_[10];
  data_[23] = C00_[2] * data_[20] + cB00_current[2] * data_[11];

  data_[24] = C00_[0] * data_[21] + B10_[0] * data_[18] + cB00_current[0] * data_[12];
  data_[25] = C00_[1] * data_[22] + B10_[1] * data_[19] + cB00_current[1] * data_[13];
  data_[26] = C00_[2] * data_[23] + B10_[2] * data_[20] + cB00_current[2] * data_[14];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];

  data_[27] = D00_[0] * data_[18] + B01_current[0] * data_[9];
  data_[28] = D00_[1] * data_[19] + B01_current[1] * data_[10];
  data_[29] = D00_[2] * data_[20] + B01_current[2] * data_[11];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];

  data_[30] = C00_[0] * data_[27] + cB00_current[0] * data_[18];
  data_[31] = C00_[1] * data_[28] + cB00_current[1] * data_[19];
  data_[32] = C00_[2] * data_[29] + cB00_current[2] * data_[20];

  data_[33] = C00_[0] * data_[30] + B10_[0] * data_[27] + cB00_current[0] * data_[21];
  data_[34] = C00_[1] * data_[31] + B10_[1] * data_[28] + cB00_current[1] * data_[22];
  data_[35] = C00_[2] * data_[32] + B10_[2] * data_[29] + cB00_current[2] * data_[23];
}

