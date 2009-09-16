//
// Author: Toru Shiozaki
// Machine Generated Code in NewInt
//

#include "vrrlist.h"

// returns double array of length 105
void VRRList::_vrr_6020(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;
  data_[3] = 1.0;
  data_[4] = 1.0;

  data_[5] = C00_[0];
  data_[6] = C00_[1];
  data_[7] = C00_[2];
  data_[8] = C00_[3];
  data_[9] = C00_[4];

  double B10_current[5];
  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];

  data_[10] = C00_[0] * data_[5] + B10_current[0];
  data_[11] = C00_[1] * data_[6] + B10_current[1];
  data_[12] = C00_[2] * data_[7] + B10_current[2];
  data_[13] = C00_[3] * data_[8] + B10_current[3];
  data_[14] = C00_[4] * data_[9] + B10_current[4];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];

  data_[15] = C00_[0] * data_[10] + B10_current[0] * data_[5];
  data_[16] = C00_[1] * data_[11] + B10_current[1] * data_[6];
  data_[17] = C00_[2] * data_[12] + B10_current[2] * data_[7];
  data_[18] = C00_[3] * data_[13] + B10_current[3] * data_[8];
  data_[19] = C00_[4] * data_[14] + B10_current[4] * data_[9];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];

  data_[20] = C00_[0] * data_[15] + B10_current[0] * data_[10];
  data_[21] = C00_[1] * data_[16] + B10_current[1] * data_[11];
  data_[22] = C00_[2] * data_[17] + B10_current[2] * data_[12];
  data_[23] = C00_[3] * data_[18] + B10_current[3] * data_[13];
  data_[24] = C00_[4] * data_[19] + B10_current[4] * data_[14];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];

  data_[25] = C00_[0] * data_[20] + B10_current[0] * data_[15];
  data_[26] = C00_[1] * data_[21] + B10_current[1] * data_[16];
  data_[27] = C00_[2] * data_[22] + B10_current[2] * data_[17];
  data_[28] = C00_[3] * data_[23] + B10_current[3] * data_[18];
  data_[29] = C00_[4] * data_[24] + B10_current[4] * data_[19];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];

  data_[30] = C00_[0] * data_[25] + B10_current[0] * data_[20];
  data_[31] = C00_[1] * data_[26] + B10_current[1] * data_[21];
  data_[32] = C00_[2] * data_[27] + B10_current[2] * data_[22];
  data_[33] = C00_[3] * data_[28] + B10_current[3] * data_[23];
  data_[34] = C00_[4] * data_[29] + B10_current[4] * data_[24];

  data_[35] = D00_[0];
  data_[36] = D00_[1];
  data_[37] = D00_[2];
  data_[38] = D00_[3];
  data_[39] = D00_[4];

  double cB00_current[5];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];
  cB00_current[4] = B00_[4];

  data_[40] = C00_[0] * data_[35] + cB00_current[0];
  data_[41] = C00_[1] * data_[36] + cB00_current[1];
  data_[42] = C00_[2] * data_[37] + cB00_current[2];
  data_[43] = C00_[3] * data_[38] + cB00_current[3];
  data_[44] = C00_[4] * data_[39] + cB00_current[4];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];

  data_[45] = C00_[0] * data_[40] + B10_current[0] * data_[35] + cB00_current[0] * data_[5];
  data_[46] = C00_[1] * data_[41] + B10_current[1] * data_[36] + cB00_current[1] * data_[6];
  data_[47] = C00_[2] * data_[42] + B10_current[2] * data_[37] + cB00_current[2] * data_[7];
  data_[48] = C00_[3] * data_[43] + B10_current[3] * data_[38] + cB00_current[3] * data_[8];
  data_[49] = C00_[4] * data_[44] + B10_current[4] * data_[39] + cB00_current[4] * data_[9];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];

  data_[50] = C00_[0] * data_[45] + B10_current[0] * data_[40] + cB00_current[0] * data_[10];
  data_[51] = C00_[1] * data_[46] + B10_current[1] * data_[41] + cB00_current[1] * data_[11];
  data_[52] = C00_[2] * data_[47] + B10_current[2] * data_[42] + cB00_current[2] * data_[12];
  data_[53] = C00_[3] * data_[48] + B10_current[3] * data_[43] + cB00_current[3] * data_[13];
  data_[54] = C00_[4] * data_[49] + B10_current[4] * data_[44] + cB00_current[4] * data_[14];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];

  data_[55] = C00_[0] * data_[50] + B10_current[0] * data_[45] + cB00_current[0] * data_[15];
  data_[56] = C00_[1] * data_[51] + B10_current[1] * data_[46] + cB00_current[1] * data_[16];
  data_[57] = C00_[2] * data_[52] + B10_current[2] * data_[47] + cB00_current[2] * data_[17];
  data_[58] = C00_[3] * data_[53] + B10_current[3] * data_[48] + cB00_current[3] * data_[18];
  data_[59] = C00_[4] * data_[54] + B10_current[4] * data_[49] + cB00_current[4] * data_[19];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];

  data_[60] = C00_[0] * data_[55] + B10_current[0] * data_[50] + cB00_current[0] * data_[20];
  data_[61] = C00_[1] * data_[56] + B10_current[1] * data_[51] + cB00_current[1] * data_[21];
  data_[62] = C00_[2] * data_[57] + B10_current[2] * data_[52] + cB00_current[2] * data_[22];
  data_[63] = C00_[3] * data_[58] + B10_current[3] * data_[53] + cB00_current[3] * data_[23];
  data_[64] = C00_[4] * data_[59] + B10_current[4] * data_[54] + cB00_current[4] * data_[24];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];

  data_[65] = C00_[0] * data_[60] + B10_current[0] * data_[55] + cB00_current[0] * data_[25];
  data_[66] = C00_[1] * data_[61] + B10_current[1] * data_[56] + cB00_current[1] * data_[26];
  data_[67] = C00_[2] * data_[62] + B10_current[2] * data_[57] + cB00_current[2] * data_[27];
  data_[68] = C00_[3] * data_[63] + B10_current[3] * data_[58] + cB00_current[3] * data_[28];
  data_[69] = C00_[4] * data_[64] + B10_current[4] * data_[59] + cB00_current[4] * data_[29];


  data_[70] = D00_[0] * data_[35] + B01_[0];
  data_[71] = D00_[1] * data_[36] + B01_[1];
  data_[72] = D00_[2] * data_[37] + B01_[2];
  data_[73] = D00_[3] * data_[38] + B01_[3];
  data_[74] = D00_[4] * data_[39] + B01_[4];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];

  data_[75] = C00_[0] * data_[70] + cB00_current[0] * data_[35];
  data_[76] = C00_[1] * data_[71] + cB00_current[1] * data_[36];
  data_[77] = C00_[2] * data_[72] + cB00_current[2] * data_[37];
  data_[78] = C00_[3] * data_[73] + cB00_current[3] * data_[38];
  data_[79] = C00_[4] * data_[74] + cB00_current[4] * data_[39];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];

  data_[80] = C00_[0] * data_[75] + B10_current[0] * data_[70] + cB00_current[0] * data_[40];
  data_[81] = C00_[1] * data_[76] + B10_current[1] * data_[71] + cB00_current[1] * data_[41];
  data_[82] = C00_[2] * data_[77] + B10_current[2] * data_[72] + cB00_current[2] * data_[42];
  data_[83] = C00_[3] * data_[78] + B10_current[3] * data_[73] + cB00_current[3] * data_[43];
  data_[84] = C00_[4] * data_[79] + B10_current[4] * data_[74] + cB00_current[4] * data_[44];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];

  data_[85] = C00_[0] * data_[80] + B10_current[0] * data_[75] + cB00_current[0] * data_[45];
  data_[86] = C00_[1] * data_[81] + B10_current[1] * data_[76] + cB00_current[1] * data_[46];
  data_[87] = C00_[2] * data_[82] + B10_current[2] * data_[77] + cB00_current[2] * data_[47];
  data_[88] = C00_[3] * data_[83] + B10_current[3] * data_[78] + cB00_current[3] * data_[48];
  data_[89] = C00_[4] * data_[84] + B10_current[4] * data_[79] + cB00_current[4] * data_[49];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];

  data_[90] = C00_[0] * data_[85] + B10_current[0] * data_[80] + cB00_current[0] * data_[50];
  data_[91] = C00_[1] * data_[86] + B10_current[1] * data_[81] + cB00_current[1] * data_[51];
  data_[92] = C00_[2] * data_[87] + B10_current[2] * data_[82] + cB00_current[2] * data_[52];
  data_[93] = C00_[3] * data_[88] + B10_current[3] * data_[83] + cB00_current[3] * data_[53];
  data_[94] = C00_[4] * data_[89] + B10_current[4] * data_[84] + cB00_current[4] * data_[54];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];

  data_[95] = C00_[0] * data_[90] + B10_current[0] * data_[85] + cB00_current[0] * data_[55];
  data_[96] = C00_[1] * data_[91] + B10_current[1] * data_[86] + cB00_current[1] * data_[56];
  data_[97] = C00_[2] * data_[92] + B10_current[2] * data_[87] + cB00_current[2] * data_[57];
  data_[98] = C00_[3] * data_[93] + B10_current[3] * data_[88] + cB00_current[3] * data_[58];
  data_[99] = C00_[4] * data_[94] + B10_current[4] * data_[89] + cB00_current[4] * data_[59];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];

  data_[100] = C00_[0] * data_[95] + B10_current[0] * data_[90] + cB00_current[0] * data_[60];
  data_[101] = C00_[1] * data_[96] + B10_current[1] * data_[91] + cB00_current[1] * data_[61];
  data_[102] = C00_[2] * data_[97] + B10_current[2] * data_[92] + cB00_current[2] * data_[62];
  data_[103] = C00_[3] * data_[98] + B10_current[3] * data_[93] + cB00_current[3] * data_[63];
  data_[104] = C00_[4] * data_[99] + B10_current[4] * data_[94] + cB00_current[4] * data_[64];
}

