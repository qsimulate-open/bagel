//
// Author: Toru Shiozaki
// Machine Generated Code
//

#include "svrrlist.h"

// returns double array of length 45
void SVRRList::_svrr_0080(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;
  data_[3] = 1.0;
  data_[4] = 1.0;

  data_[5] = D00_[0];
  data_[6] = D00_[1];
  data_[7] = D00_[2];
  data_[8] = D00_[3];
  data_[9] = D00_[4];

  double B01_current[5];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];
  B01_current[4] = B01_[4];

  data_[10] = D00_[0] * data_[5] + B01_current[0];
  data_[11] = D00_[1] * data_[6] + B01_current[1];
  data_[12] = D00_[2] * data_[7] + B01_current[2];
  data_[13] = D00_[3] * data_[8] + B01_current[3];
  data_[14] = D00_[4] * data_[9] + B01_current[4];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];

  data_[15] = D00_[0] * data_[10] + B01_current[0] * data_[5];
  data_[16] = D00_[1] * data_[11] + B01_current[1] * data_[6];
  data_[17] = D00_[2] * data_[12] + B01_current[2] * data_[7];
  data_[18] = D00_[3] * data_[13] + B01_current[3] * data_[8];
  data_[19] = D00_[4] * data_[14] + B01_current[4] * data_[9];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];

  data_[20] = D00_[0] * data_[15] + B01_current[0] * data_[10];
  data_[21] = D00_[1] * data_[16] + B01_current[1] * data_[11];
  data_[22] = D00_[2] * data_[17] + B01_current[2] * data_[12];
  data_[23] = D00_[3] * data_[18] + B01_current[3] * data_[13];
  data_[24] = D00_[4] * data_[19] + B01_current[4] * data_[14];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];

  data_[25] = D00_[0] * data_[20] + B01_current[0] * data_[15];
  data_[26] = D00_[1] * data_[21] + B01_current[1] * data_[16];
  data_[27] = D00_[2] * data_[22] + B01_current[2] * data_[17];
  data_[28] = D00_[3] * data_[23] + B01_current[3] * data_[18];
  data_[29] = D00_[4] * data_[24] + B01_current[4] * data_[19];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];

  data_[30] = D00_[0] * data_[25] + B01_current[0] * data_[20];
  data_[31] = D00_[1] * data_[26] + B01_current[1] * data_[21];
  data_[32] = D00_[2] * data_[27] + B01_current[2] * data_[22];
  data_[33] = D00_[3] * data_[28] + B01_current[3] * data_[23];
  data_[34] = D00_[4] * data_[29] + B01_current[4] * data_[24];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];

  data_[35] = D00_[0] * data_[30] + B01_current[0] * data_[25];
  data_[36] = D00_[1] * data_[31] + B01_current[1] * data_[26];
  data_[37] = D00_[2] * data_[32] + B01_current[2] * data_[27];
  data_[38] = D00_[3] * data_[33] + B01_current[3] * data_[28];
  data_[39] = D00_[4] * data_[34] + B01_current[4] * data_[29];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];

  data_[40] = D00_[0] * data_[35] + B01_current[0] * data_[30];
  data_[41] = D00_[1] * data_[36] + B01_current[1] * data_[31];
  data_[42] = D00_[2] * data_[37] + B01_current[2] * data_[32];
  data_[43] = D00_[3] * data_[38] + B01_current[3] * data_[33];
  data_[44] = D00_[4] * data_[39] + B01_current[4] * data_[34];
}

