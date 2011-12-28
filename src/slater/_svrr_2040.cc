//
// Author: Toru Shiozaki
// Machine Generated Code
//

#include <src/slater/svrrlist.h>

// returns double array of length 60
void SVRRList::_svrr_2040(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;
  data_[3] = 1.0;

  data_[4] = C00_[0];
  data_[5] = C00_[1];
  data_[6] = C00_[2];
  data_[7] = C00_[3];

  data_[8] = C00_[0] * data_[4] + B10_[0];
  data_[9] = C00_[1] * data_[5] + B10_[1];
  data_[10] = C00_[2] * data_[6] + B10_[2];
  data_[11] = C00_[3] * data_[7] + B10_[3];

  data_[12] = D00_[0];
  data_[13] = D00_[1];
  data_[14] = D00_[2];
  data_[15] = D00_[3];

  double cB00_current[4];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];

  data_[16] = C00_[0] * data_[12] + cB00_current[0];
  data_[17] = C00_[1] * data_[13] + cB00_current[1];
  data_[18] = C00_[2] * data_[14] + cB00_current[2];
  data_[19] = C00_[3] * data_[15] + cB00_current[3];

  data_[20] = C00_[0] * data_[16] + B10_[0] * data_[12] + cB00_current[0] * data_[4];
  data_[21] = C00_[1] * data_[17] + B10_[1] * data_[13] + cB00_current[1] * data_[5];
  data_[22] = C00_[2] * data_[18] + B10_[2] * data_[14] + cB00_current[2] * data_[6];
  data_[23] = C00_[3] * data_[19] + B10_[3] * data_[15] + cB00_current[3] * data_[7];

  double B01_current[4];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];

  data_[24] = D00_[0] * data_[12] + B01_current[0];
  data_[25] = D00_[1] * data_[13] + B01_current[1];
  data_[26] = D00_[2] * data_[14] + B01_current[2];
  data_[27] = D00_[3] * data_[15] + B01_current[3];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];

  data_[28] = C00_[0] * data_[24] + cB00_current[0] * data_[12];
  data_[29] = C00_[1] * data_[25] + cB00_current[1] * data_[13];
  data_[30] = C00_[2] * data_[26] + cB00_current[2] * data_[14];
  data_[31] = C00_[3] * data_[27] + cB00_current[3] * data_[15];

  data_[32] = C00_[0] * data_[28] + B10_[0] * data_[24] + cB00_current[0] * data_[16];
  data_[33] = C00_[1] * data_[29] + B10_[1] * data_[25] + cB00_current[1] * data_[17];
  data_[34] = C00_[2] * data_[30] + B10_[2] * data_[26] + cB00_current[2] * data_[18];
  data_[35] = C00_[3] * data_[31] + B10_[3] * data_[27] + cB00_current[3] * data_[19];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];

  data_[36] = D00_[0] * data_[24] + B01_current[0] * data_[12];
  data_[37] = D00_[1] * data_[25] + B01_current[1] * data_[13];
  data_[38] = D00_[2] * data_[26] + B01_current[2] * data_[14];
  data_[39] = D00_[3] * data_[27] + B01_current[3] * data_[15];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];

  data_[40] = C00_[0] * data_[36] + cB00_current[0] * data_[24];
  data_[41] = C00_[1] * data_[37] + cB00_current[1] * data_[25];
  data_[42] = C00_[2] * data_[38] + cB00_current[2] * data_[26];
  data_[43] = C00_[3] * data_[39] + cB00_current[3] * data_[27];

  data_[44] = C00_[0] * data_[40] + B10_[0] * data_[36] + cB00_current[0] * data_[28];
  data_[45] = C00_[1] * data_[41] + B10_[1] * data_[37] + cB00_current[1] * data_[29];
  data_[46] = C00_[2] * data_[42] + B10_[2] * data_[38] + cB00_current[2] * data_[30];
  data_[47] = C00_[3] * data_[43] + B10_[3] * data_[39] + cB00_current[3] * data_[31];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];

  data_[48] = D00_[0] * data_[36] + B01_current[0] * data_[24];
  data_[49] = D00_[1] * data_[37] + B01_current[1] * data_[25];
  data_[50] = D00_[2] * data_[38] + B01_current[2] * data_[26];
  data_[51] = D00_[3] * data_[39] + B01_current[3] * data_[27];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];

  data_[52] = C00_[0] * data_[48] + cB00_current[0] * data_[36];
  data_[53] = C00_[1] * data_[49] + cB00_current[1] * data_[37];
  data_[54] = C00_[2] * data_[50] + cB00_current[2] * data_[38];
  data_[55] = C00_[3] * data_[51] + cB00_current[3] * data_[39];

  data_[56] = C00_[0] * data_[52] + B10_[0] * data_[48] + cB00_current[0] * data_[40];
  data_[57] = C00_[1] * data_[53] + B10_[1] * data_[49] + cB00_current[1] * data_[41];
  data_[58] = C00_[2] * data_[54] + B10_[2] * data_[50] + cB00_current[2] * data_[42];
  data_[59] = C00_[3] * data_[55] + B10_[3] * data_[51] + cB00_current[3] * data_[43];
}

