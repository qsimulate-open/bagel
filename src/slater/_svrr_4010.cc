//
// Author: Toru Shiozaki
// Machine Generated Code
//

#include <src/slater/svrrlist.h>

// returns double array of length 40
void SVRRList::_svrr_4010(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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

  data_[20] = D00_[0];
  data_[21] = D00_[1];
  data_[22] = D00_[2];
  data_[23] = D00_[3];

  data_[24] = C00_[0] * data_[20] + B00_[0];
  data_[25] = C00_[1] * data_[21] + B00_[1];
  data_[26] = C00_[2] * data_[22] + B00_[2];
  data_[27] = C00_[3] * data_[23] + B00_[3];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];

  data_[28] = C00_[0] * data_[24] + B10_current[0] * data_[20] + B00_[0] * data_[4];
  data_[29] = C00_[1] * data_[25] + B10_current[1] * data_[21] + B00_[1] * data_[5];
  data_[30] = C00_[2] * data_[26] + B10_current[2] * data_[22] + B00_[2] * data_[6];
  data_[31] = C00_[3] * data_[27] + B10_current[3] * data_[23] + B00_[3] * data_[7];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];

  data_[32] = C00_[0] * data_[28] + B10_current[0] * data_[24] + B00_[0] * data_[8];
  data_[33] = C00_[1] * data_[29] + B10_current[1] * data_[25] + B00_[1] * data_[9];
  data_[34] = C00_[2] * data_[30] + B10_current[2] * data_[26] + B00_[2] * data_[10];
  data_[35] = C00_[3] * data_[31] + B10_current[3] * data_[27] + B00_[3] * data_[11];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];

  data_[36] = C00_[0] * data_[32] + B10_current[0] * data_[28] + B00_[0] * data_[12];
  data_[37] = C00_[1] * data_[33] + B10_current[1] * data_[29] + B00_[1] * data_[13];
  data_[38] = C00_[2] * data_[34] + B10_current[2] * data_[30] + B00_[2] * data_[14];
  data_[39] = C00_[3] * data_[35] + B10_current[3] * data_[31] + B00_[3] * data_[15];
}

