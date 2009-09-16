//
// Author: Toru Shiozaki
// Machine Generated Code in NewInt
//

#include "vrrlist.h"

// returns double array of length 132
void VRRList::_vrr_10a0(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;
  data_[3] = 1.0;
  data_[4] = 1.0;
  data_[5] = 1.0;

  data_[6] = C00_[0];
  data_[7] = C00_[1];
  data_[8] = C00_[2];
  data_[9] = C00_[3];
  data_[10] = C00_[4];
  data_[11] = C00_[5];

  data_[12] = D00_[0];
  data_[13] = D00_[1];
  data_[14] = D00_[2];
  data_[15] = D00_[3];
  data_[16] = D00_[4];
  data_[17] = D00_[5];

  double cB00_current[6];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];
  cB00_current[4] = B00_[4];
  cB00_current[5] = B00_[5];

  data_[18] = C00_[0] * data_[12] + cB00_current[0];
  data_[19] = C00_[1] * data_[13] + cB00_current[1];
  data_[20] = C00_[2] * data_[14] + cB00_current[2];
  data_[21] = C00_[3] * data_[15] + cB00_current[3];
  data_[22] = C00_[4] * data_[16] + cB00_current[4];
  data_[23] = C00_[5] * data_[17] + cB00_current[5];

  double B01_current[6];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];
  B01_current[4] = B01_[4];
  B01_current[5] = B01_[5];

  data_[24] = D00_[0] * data_[12] + B01_current[0];
  data_[25] = D00_[1] * data_[13] + B01_current[1];
  data_[26] = D00_[2] * data_[14] + B01_current[2];
  data_[27] = D00_[3] * data_[15] + B01_current[3];
  data_[28] = D00_[4] * data_[16] + B01_current[4];
  data_[29] = D00_[5] * data_[17] + B01_current[5];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];

  data_[30] = C00_[0] * data_[24] + cB00_current[0] * data_[12];
  data_[31] = C00_[1] * data_[25] + cB00_current[1] * data_[13];
  data_[32] = C00_[2] * data_[26] + cB00_current[2] * data_[14];
  data_[33] = C00_[3] * data_[27] + cB00_current[3] * data_[15];
  data_[34] = C00_[4] * data_[28] + cB00_current[4] * data_[16];
  data_[35] = C00_[5] * data_[29] + cB00_current[5] * data_[17];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];

  data_[36] = D00_[0] * data_[24] + B01_current[0] * data_[12];
  data_[37] = D00_[1] * data_[25] + B01_current[1] * data_[13];
  data_[38] = D00_[2] * data_[26] + B01_current[2] * data_[14];
  data_[39] = D00_[3] * data_[27] + B01_current[3] * data_[15];
  data_[40] = D00_[4] * data_[28] + B01_current[4] * data_[16];
  data_[41] = D00_[5] * data_[29] + B01_current[5] * data_[17];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];

  data_[42] = C00_[0] * data_[36] + cB00_current[0] * data_[24];
  data_[43] = C00_[1] * data_[37] + cB00_current[1] * data_[25];
  data_[44] = C00_[2] * data_[38] + cB00_current[2] * data_[26];
  data_[45] = C00_[3] * data_[39] + cB00_current[3] * data_[27];
  data_[46] = C00_[4] * data_[40] + cB00_current[4] * data_[28];
  data_[47] = C00_[5] * data_[41] + cB00_current[5] * data_[29];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];

  data_[48] = D00_[0] * data_[36] + B01_current[0] * data_[24];
  data_[49] = D00_[1] * data_[37] + B01_current[1] * data_[25];
  data_[50] = D00_[2] * data_[38] + B01_current[2] * data_[26];
  data_[51] = D00_[3] * data_[39] + B01_current[3] * data_[27];
  data_[52] = D00_[4] * data_[40] + B01_current[4] * data_[28];
  data_[53] = D00_[5] * data_[41] + B01_current[5] * data_[29];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];

  data_[54] = C00_[0] * data_[48] + cB00_current[0] * data_[36];
  data_[55] = C00_[1] * data_[49] + cB00_current[1] * data_[37];
  data_[56] = C00_[2] * data_[50] + cB00_current[2] * data_[38];
  data_[57] = C00_[3] * data_[51] + cB00_current[3] * data_[39];
  data_[58] = C00_[4] * data_[52] + cB00_current[4] * data_[40];
  data_[59] = C00_[5] * data_[53] + cB00_current[5] * data_[41];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];

  data_[60] = D00_[0] * data_[48] + B01_current[0] * data_[36];
  data_[61] = D00_[1] * data_[49] + B01_current[1] * data_[37];
  data_[62] = D00_[2] * data_[50] + B01_current[2] * data_[38];
  data_[63] = D00_[3] * data_[51] + B01_current[3] * data_[39];
  data_[64] = D00_[4] * data_[52] + B01_current[4] * data_[40];
  data_[65] = D00_[5] * data_[53] + B01_current[5] * data_[41];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];

  data_[66] = C00_[0] * data_[60] + cB00_current[0] * data_[48];
  data_[67] = C00_[1] * data_[61] + cB00_current[1] * data_[49];
  data_[68] = C00_[2] * data_[62] + cB00_current[2] * data_[50];
  data_[69] = C00_[3] * data_[63] + cB00_current[3] * data_[51];
  data_[70] = C00_[4] * data_[64] + cB00_current[4] * data_[52];
  data_[71] = C00_[5] * data_[65] + cB00_current[5] * data_[53];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];

  data_[72] = D00_[0] * data_[60] + B01_current[0] * data_[48];
  data_[73] = D00_[1] * data_[61] + B01_current[1] * data_[49];
  data_[74] = D00_[2] * data_[62] + B01_current[2] * data_[50];
  data_[75] = D00_[3] * data_[63] + B01_current[3] * data_[51];
  data_[76] = D00_[4] * data_[64] + B01_current[4] * data_[52];
  data_[77] = D00_[5] * data_[65] + B01_current[5] * data_[53];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];

  data_[78] = C00_[0] * data_[72] + cB00_current[0] * data_[60];
  data_[79] = C00_[1] * data_[73] + cB00_current[1] * data_[61];
  data_[80] = C00_[2] * data_[74] + cB00_current[2] * data_[62];
  data_[81] = C00_[3] * data_[75] + cB00_current[3] * data_[63];
  data_[82] = C00_[4] * data_[76] + cB00_current[4] * data_[64];
  data_[83] = C00_[5] * data_[77] + cB00_current[5] * data_[65];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];

  data_[84] = D00_[0] * data_[72] + B01_current[0] * data_[60];
  data_[85] = D00_[1] * data_[73] + B01_current[1] * data_[61];
  data_[86] = D00_[2] * data_[74] + B01_current[2] * data_[62];
  data_[87] = D00_[3] * data_[75] + B01_current[3] * data_[63];
  data_[88] = D00_[4] * data_[76] + B01_current[4] * data_[64];
  data_[89] = D00_[5] * data_[77] + B01_current[5] * data_[65];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];

  data_[90] = C00_[0] * data_[84] + cB00_current[0] * data_[72];
  data_[91] = C00_[1] * data_[85] + cB00_current[1] * data_[73];
  data_[92] = C00_[2] * data_[86] + cB00_current[2] * data_[74];
  data_[93] = C00_[3] * data_[87] + cB00_current[3] * data_[75];
  data_[94] = C00_[4] * data_[88] + cB00_current[4] * data_[76];
  data_[95] = C00_[5] * data_[89] + cB00_current[5] * data_[77];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];

  data_[96] = D00_[0] * data_[84] + B01_current[0] * data_[72];
  data_[97] = D00_[1] * data_[85] + B01_current[1] * data_[73];
  data_[98] = D00_[2] * data_[86] + B01_current[2] * data_[74];
  data_[99] = D00_[3] * data_[87] + B01_current[3] * data_[75];
  data_[100] = D00_[4] * data_[88] + B01_current[4] * data_[76];
  data_[101] = D00_[5] * data_[89] + B01_current[5] * data_[77];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];

  data_[102] = C00_[0] * data_[96] + cB00_current[0] * data_[84];
  data_[103] = C00_[1] * data_[97] + cB00_current[1] * data_[85];
  data_[104] = C00_[2] * data_[98] + cB00_current[2] * data_[86];
  data_[105] = C00_[3] * data_[99] + cB00_current[3] * data_[87];
  data_[106] = C00_[4] * data_[100] + cB00_current[4] * data_[88];
  data_[107] = C00_[5] * data_[101] + cB00_current[5] * data_[89];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];

  data_[108] = D00_[0] * data_[96] + B01_current[0] * data_[84];
  data_[109] = D00_[1] * data_[97] + B01_current[1] * data_[85];
  data_[110] = D00_[2] * data_[98] + B01_current[2] * data_[86];
  data_[111] = D00_[3] * data_[99] + B01_current[3] * data_[87];
  data_[112] = D00_[4] * data_[100] + B01_current[4] * data_[88];
  data_[113] = D00_[5] * data_[101] + B01_current[5] * data_[89];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];

  data_[114] = C00_[0] * data_[108] + cB00_current[0] * data_[96];
  data_[115] = C00_[1] * data_[109] + cB00_current[1] * data_[97];
  data_[116] = C00_[2] * data_[110] + cB00_current[2] * data_[98];
  data_[117] = C00_[3] * data_[111] + cB00_current[3] * data_[99];
  data_[118] = C00_[4] * data_[112] + cB00_current[4] * data_[100];
  data_[119] = C00_[5] * data_[113] + cB00_current[5] * data_[101];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];

  data_[120] = D00_[0] * data_[108] + B01_current[0] * data_[96];
  data_[121] = D00_[1] * data_[109] + B01_current[1] * data_[97];
  data_[122] = D00_[2] * data_[110] + B01_current[2] * data_[98];
  data_[123] = D00_[3] * data_[111] + B01_current[3] * data_[99];
  data_[124] = D00_[4] * data_[112] + B01_current[4] * data_[100];
  data_[125] = D00_[5] * data_[113] + B01_current[5] * data_[101];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];

  data_[126] = C00_[0] * data_[120] + cB00_current[0] * data_[108];
  data_[127] = C00_[1] * data_[121] + cB00_current[1] * data_[109];
  data_[128] = C00_[2] * data_[122] + cB00_current[2] * data_[110];
  data_[129] = C00_[3] * data_[123] + cB00_current[3] * data_[111];
  data_[130] = C00_[4] * data_[124] + cB00_current[4] * data_[112];
  data_[131] = C00_[5] * data_[125] + cB00_current[5] * data_[113];
}

