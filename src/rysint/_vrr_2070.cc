//
// Author: Toru Shiozaki
// Machine Generated Code in NewInt
//

#include "vrrlist.h"

// returns double array of length 120
void VRRList::_vrr_2070(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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

  data_[10] = C00_[0] * data_[5] + B10_[0];
  data_[11] = C00_[1] * data_[6] + B10_[1];
  data_[12] = C00_[2] * data_[7] + B10_[2];
  data_[13] = C00_[3] * data_[8] + B10_[3];
  data_[14] = C00_[4] * data_[9] + B10_[4];

  data_[15] = D00_[0];
  data_[16] = D00_[1];
  data_[17] = D00_[2];
  data_[18] = D00_[3];
  data_[19] = D00_[4];

  double cB00_current[5];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];
  cB00_current[4] = B00_[4];

  data_[20] = C00_[0] * data_[15] + cB00_current[0];
  data_[21] = C00_[1] * data_[16] + cB00_current[1];
  data_[22] = C00_[2] * data_[17] + cB00_current[2];
  data_[23] = C00_[3] * data_[18] + cB00_current[3];
  data_[24] = C00_[4] * data_[19] + cB00_current[4];

  data_[25] = C00_[0] * data_[20] + B10_[0] * data_[15] + cB00_current[0] * data_[5];
  data_[26] = C00_[1] * data_[21] + B10_[1] * data_[16] + cB00_current[1] * data_[6];
  data_[27] = C00_[2] * data_[22] + B10_[2] * data_[17] + cB00_current[2] * data_[7];
  data_[28] = C00_[3] * data_[23] + B10_[3] * data_[18] + cB00_current[3] * data_[8];
  data_[29] = C00_[4] * data_[24] + B10_[4] * data_[19] + cB00_current[4] * data_[9];

  double B01_current[5];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];
  B01_current[4] = B01_[4];

  data_[30] = D00_[0] * data_[15] + B01_current[0];
  data_[31] = D00_[1] * data_[16] + B01_current[1];
  data_[32] = D00_[2] * data_[17] + B01_current[2];
  data_[33] = D00_[3] * data_[18] + B01_current[3];
  data_[34] = D00_[4] * data_[19] + B01_current[4];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];

  data_[35] = C00_[0] * data_[30] + cB00_current[0] * data_[15];
  data_[36] = C00_[1] * data_[31] + cB00_current[1] * data_[16];
  data_[37] = C00_[2] * data_[32] + cB00_current[2] * data_[17];
  data_[38] = C00_[3] * data_[33] + cB00_current[3] * data_[18];
  data_[39] = C00_[4] * data_[34] + cB00_current[4] * data_[19];

  data_[40] = C00_[0] * data_[35] + B10_[0] * data_[30] + cB00_current[0] * data_[20];
  data_[41] = C00_[1] * data_[36] + B10_[1] * data_[31] + cB00_current[1] * data_[21];
  data_[42] = C00_[2] * data_[37] + B10_[2] * data_[32] + cB00_current[2] * data_[22];
  data_[43] = C00_[3] * data_[38] + B10_[3] * data_[33] + cB00_current[3] * data_[23];
  data_[44] = C00_[4] * data_[39] + B10_[4] * data_[34] + cB00_current[4] * data_[24];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];

  data_[45] = D00_[0] * data_[30] + B01_current[0] * data_[15];
  data_[46] = D00_[1] * data_[31] + B01_current[1] * data_[16];
  data_[47] = D00_[2] * data_[32] + B01_current[2] * data_[17];
  data_[48] = D00_[3] * data_[33] + B01_current[3] * data_[18];
  data_[49] = D00_[4] * data_[34] + B01_current[4] * data_[19];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];

  data_[50] = C00_[0] * data_[45] + cB00_current[0] * data_[30];
  data_[51] = C00_[1] * data_[46] + cB00_current[1] * data_[31];
  data_[52] = C00_[2] * data_[47] + cB00_current[2] * data_[32];
  data_[53] = C00_[3] * data_[48] + cB00_current[3] * data_[33];
  data_[54] = C00_[4] * data_[49] + cB00_current[4] * data_[34];

  data_[55] = C00_[0] * data_[50] + B10_[0] * data_[45] + cB00_current[0] * data_[35];
  data_[56] = C00_[1] * data_[51] + B10_[1] * data_[46] + cB00_current[1] * data_[36];
  data_[57] = C00_[2] * data_[52] + B10_[2] * data_[47] + cB00_current[2] * data_[37];
  data_[58] = C00_[3] * data_[53] + B10_[3] * data_[48] + cB00_current[3] * data_[38];
  data_[59] = C00_[4] * data_[54] + B10_[4] * data_[49] + cB00_current[4] * data_[39];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];

  data_[60] = D00_[0] * data_[45] + B01_current[0] * data_[30];
  data_[61] = D00_[1] * data_[46] + B01_current[1] * data_[31];
  data_[62] = D00_[2] * data_[47] + B01_current[2] * data_[32];
  data_[63] = D00_[3] * data_[48] + B01_current[3] * data_[33];
  data_[64] = D00_[4] * data_[49] + B01_current[4] * data_[34];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];

  data_[65] = C00_[0] * data_[60] + cB00_current[0] * data_[45];
  data_[66] = C00_[1] * data_[61] + cB00_current[1] * data_[46];
  data_[67] = C00_[2] * data_[62] + cB00_current[2] * data_[47];
  data_[68] = C00_[3] * data_[63] + cB00_current[3] * data_[48];
  data_[69] = C00_[4] * data_[64] + cB00_current[4] * data_[49];

  data_[70] = C00_[0] * data_[65] + B10_[0] * data_[60] + cB00_current[0] * data_[50];
  data_[71] = C00_[1] * data_[66] + B10_[1] * data_[61] + cB00_current[1] * data_[51];
  data_[72] = C00_[2] * data_[67] + B10_[2] * data_[62] + cB00_current[2] * data_[52];
  data_[73] = C00_[3] * data_[68] + B10_[3] * data_[63] + cB00_current[3] * data_[53];
  data_[74] = C00_[4] * data_[69] + B10_[4] * data_[64] + cB00_current[4] * data_[54];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];

  data_[75] = D00_[0] * data_[60] + B01_current[0] * data_[45];
  data_[76] = D00_[1] * data_[61] + B01_current[1] * data_[46];
  data_[77] = D00_[2] * data_[62] + B01_current[2] * data_[47];
  data_[78] = D00_[3] * data_[63] + B01_current[3] * data_[48];
  data_[79] = D00_[4] * data_[64] + B01_current[4] * data_[49];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];

  data_[80] = C00_[0] * data_[75] + cB00_current[0] * data_[60];
  data_[81] = C00_[1] * data_[76] + cB00_current[1] * data_[61];
  data_[82] = C00_[2] * data_[77] + cB00_current[2] * data_[62];
  data_[83] = C00_[3] * data_[78] + cB00_current[3] * data_[63];
  data_[84] = C00_[4] * data_[79] + cB00_current[4] * data_[64];

  data_[85] = C00_[0] * data_[80] + B10_[0] * data_[75] + cB00_current[0] * data_[65];
  data_[86] = C00_[1] * data_[81] + B10_[1] * data_[76] + cB00_current[1] * data_[66];
  data_[87] = C00_[2] * data_[82] + B10_[2] * data_[77] + cB00_current[2] * data_[67];
  data_[88] = C00_[3] * data_[83] + B10_[3] * data_[78] + cB00_current[3] * data_[68];
  data_[89] = C00_[4] * data_[84] + B10_[4] * data_[79] + cB00_current[4] * data_[69];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];

  data_[90] = D00_[0] * data_[75] + B01_current[0] * data_[60];
  data_[91] = D00_[1] * data_[76] + B01_current[1] * data_[61];
  data_[92] = D00_[2] * data_[77] + B01_current[2] * data_[62];
  data_[93] = D00_[3] * data_[78] + B01_current[3] * data_[63];
  data_[94] = D00_[4] * data_[79] + B01_current[4] * data_[64];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];

  data_[95] = C00_[0] * data_[90] + cB00_current[0] * data_[75];
  data_[96] = C00_[1] * data_[91] + cB00_current[1] * data_[76];
  data_[97] = C00_[2] * data_[92] + cB00_current[2] * data_[77];
  data_[98] = C00_[3] * data_[93] + cB00_current[3] * data_[78];
  data_[99] = C00_[4] * data_[94] + cB00_current[4] * data_[79];

  data_[100] = C00_[0] * data_[95] + B10_[0] * data_[90] + cB00_current[0] * data_[80];
  data_[101] = C00_[1] * data_[96] + B10_[1] * data_[91] + cB00_current[1] * data_[81];
  data_[102] = C00_[2] * data_[97] + B10_[2] * data_[92] + cB00_current[2] * data_[82];
  data_[103] = C00_[3] * data_[98] + B10_[3] * data_[93] + cB00_current[3] * data_[83];
  data_[104] = C00_[4] * data_[99] + B10_[4] * data_[94] + cB00_current[4] * data_[84];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];

  data_[105] = D00_[0] * data_[90] + B01_current[0] * data_[75];
  data_[106] = D00_[1] * data_[91] + B01_current[1] * data_[76];
  data_[107] = D00_[2] * data_[92] + B01_current[2] * data_[77];
  data_[108] = D00_[3] * data_[93] + B01_current[3] * data_[78];
  data_[109] = D00_[4] * data_[94] + B01_current[4] * data_[79];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];

  data_[110] = C00_[0] * data_[105] + cB00_current[0] * data_[90];
  data_[111] = C00_[1] * data_[106] + cB00_current[1] * data_[91];
  data_[112] = C00_[2] * data_[107] + cB00_current[2] * data_[92];
  data_[113] = C00_[3] * data_[108] + cB00_current[3] * data_[93];
  data_[114] = C00_[4] * data_[109] + cB00_current[4] * data_[94];

  data_[115] = C00_[0] * data_[110] + B10_[0] * data_[105] + cB00_current[0] * data_[95];
  data_[116] = C00_[1] * data_[111] + B10_[1] * data_[106] + cB00_current[1] * data_[96];
  data_[117] = C00_[2] * data_[112] + B10_[2] * data_[107] + cB00_current[2] * data_[97];
  data_[118] = C00_[3] * data_[113] + B10_[3] * data_[108] + cB00_current[3] * data_[98];
  data_[119] = C00_[4] * data_[114] + B10_[4] * data_[109] + cB00_current[4] * data_[99];
}

