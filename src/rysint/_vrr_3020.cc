//
// Author: Toru Shiozaki
// Machine Generated Code in NewInt
//

#include "vrrlist.h"

// returns double array of length 36
void VRRList::_vrr_3020(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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

  data_[12] = D00_[0];
  data_[13] = D00_[1];
  data_[14] = D00_[2];

  double cB00_current[3];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];

  data_[15] = C00_[0] * data_[12] + cB00_current[0];
  data_[16] = C00_[1] * data_[13] + cB00_current[1];
  data_[17] = C00_[2] * data_[14] + cB00_current[2];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];

  data_[18] = C00_[0] * data_[15] + B10_current[0] * data_[12] + cB00_current[0] * data_[3];
  data_[19] = C00_[1] * data_[16] + B10_current[1] * data_[13] + cB00_current[1] * data_[4];
  data_[20] = C00_[2] * data_[17] + B10_current[2] * data_[14] + cB00_current[2] * data_[5];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];

  data_[21] = C00_[0] * data_[18] + B10_current[0] * data_[15] + cB00_current[0] * data_[6];
  data_[22] = C00_[1] * data_[19] + B10_current[1] * data_[16] + cB00_current[1] * data_[7];
  data_[23] = C00_[2] * data_[20] + B10_current[2] * data_[17] + cB00_current[2] * data_[8];


  data_[24] = D00_[0] * data_[12] + B01_[0];
  data_[25] = D00_[1] * data_[13] + B01_[1];
  data_[26] = D00_[2] * data_[14] + B01_[2];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];

  data_[27] = C00_[0] * data_[24] + cB00_current[0] * data_[12];
  data_[28] = C00_[1] * data_[25] + cB00_current[1] * data_[13];
  data_[29] = C00_[2] * data_[26] + cB00_current[2] * data_[14];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];

  data_[30] = C00_[0] * data_[27] + B10_current[0] * data_[24] + cB00_current[0] * data_[15];
  data_[31] = C00_[1] * data_[28] + B10_current[1] * data_[25] + cB00_current[1] * data_[16];
  data_[32] = C00_[2] * data_[29] + B10_current[2] * data_[26] + cB00_current[2] * data_[17];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];

  data_[33] = C00_[0] * data_[30] + B10_current[0] * data_[27] + cB00_current[0] * data_[18];
  data_[34] = C00_[1] * data_[31] + B10_current[1] * data_[28] + cB00_current[1] * data_[19];
  data_[35] = C00_[2] * data_[32] + B10_current[2] * data_[29] + cB00_current[2] * data_[20];
}

