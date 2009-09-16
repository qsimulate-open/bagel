//
// Author: Toru Shiozaki
// Machine Generated Code
//

#include "svrrlist.h"

// returns double array of length 64
void SVRRList::_svrr_3030(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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

  data_[16] = D00_[0];
  data_[17] = D00_[1];
  data_[18] = D00_[2];
  data_[19] = D00_[3];

  double cB00_current[4];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];

  data_[20] = C00_[0] * data_[16] + cB00_current[0];
  data_[21] = C00_[1] * data_[17] + cB00_current[1];
  data_[22] = C00_[2] * data_[18] + cB00_current[2];
  data_[23] = C00_[3] * data_[19] + cB00_current[3];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];

  data_[24] = C00_[0] * data_[20] + B10_current[0] * data_[16] + cB00_current[0] * data_[4];
  data_[25] = C00_[1] * data_[21] + B10_current[1] * data_[17] + cB00_current[1] * data_[5];
  data_[26] = C00_[2] * data_[22] + B10_current[2] * data_[18] + cB00_current[2] * data_[6];
  data_[27] = C00_[3] * data_[23] + B10_current[3] * data_[19] + cB00_current[3] * data_[7];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];

  data_[28] = C00_[0] * data_[24] + B10_current[0] * data_[20] + cB00_current[0] * data_[8];
  data_[29] = C00_[1] * data_[25] + B10_current[1] * data_[21] + cB00_current[1] * data_[9];
  data_[30] = C00_[2] * data_[26] + B10_current[2] * data_[22] + cB00_current[2] * data_[10];
  data_[31] = C00_[3] * data_[27] + B10_current[3] * data_[23] + cB00_current[3] * data_[11];

  double B01_current[4];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];

  data_[32] = D00_[0] * data_[16] + B01_current[0];
  data_[33] = D00_[1] * data_[17] + B01_current[1];
  data_[34] = D00_[2] * data_[18] + B01_current[2];
  data_[35] = D00_[3] * data_[19] + B01_current[3];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];

  data_[36] = C00_[0] * data_[32] + cB00_current[0] * data_[16];
  data_[37] = C00_[1] * data_[33] + cB00_current[1] * data_[17];
  data_[38] = C00_[2] * data_[34] + cB00_current[2] * data_[18];
  data_[39] = C00_[3] * data_[35] + cB00_current[3] * data_[19];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];

  data_[40] = C00_[0] * data_[36] + B10_current[0] * data_[32] + cB00_current[0] * data_[20];
  data_[41] = C00_[1] * data_[37] + B10_current[1] * data_[33] + cB00_current[1] * data_[21];
  data_[42] = C00_[2] * data_[38] + B10_current[2] * data_[34] + cB00_current[2] * data_[22];
  data_[43] = C00_[3] * data_[39] + B10_current[3] * data_[35] + cB00_current[3] * data_[23];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];

  data_[44] = C00_[0] * data_[40] + B10_current[0] * data_[36] + cB00_current[0] * data_[24];
  data_[45] = C00_[1] * data_[41] + B10_current[1] * data_[37] + cB00_current[1] * data_[25];
  data_[46] = C00_[2] * data_[42] + B10_current[2] * data_[38] + cB00_current[2] * data_[26];
  data_[47] = C00_[3] * data_[43] + B10_current[3] * data_[39] + cB00_current[3] * data_[27];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];

  data_[48] = D00_[0] * data_[32] + B01_current[0] * data_[16];
  data_[49] = D00_[1] * data_[33] + B01_current[1] * data_[17];
  data_[50] = D00_[2] * data_[34] + B01_current[2] * data_[18];
  data_[51] = D00_[3] * data_[35] + B01_current[3] * data_[19];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];

  data_[52] = C00_[0] * data_[48] + cB00_current[0] * data_[32];
  data_[53] = C00_[1] * data_[49] + cB00_current[1] * data_[33];
  data_[54] = C00_[2] * data_[50] + cB00_current[2] * data_[34];
  data_[55] = C00_[3] * data_[51] + cB00_current[3] * data_[35];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];

  data_[56] = C00_[0] * data_[52] + B10_current[0] * data_[48] + cB00_current[0] * data_[36];
  data_[57] = C00_[1] * data_[53] + B10_current[1] * data_[49] + cB00_current[1] * data_[37];
  data_[58] = C00_[2] * data_[54] + B10_current[2] * data_[50] + cB00_current[2] * data_[38];
  data_[59] = C00_[3] * data_[55] + B10_current[3] * data_[51] + cB00_current[3] * data_[39];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];

  data_[60] = C00_[0] * data_[56] + B10_current[0] * data_[52] + cB00_current[0] * data_[40];
  data_[61] = C00_[1] * data_[57] + B10_current[1] * data_[53] + cB00_current[1] * data_[41];
  data_[62] = C00_[2] * data_[58] + B10_current[2] * data_[54] + cB00_current[2] * data_[42];
  data_[63] = C00_[3] * data_[59] + B10_current[3] * data_[55] + cB00_current[3] * data_[43];
}

