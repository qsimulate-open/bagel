//
// Author: Toru Shiozaki
// Machine Generated Code
//

#include "svrrlist.h"

// returns double array of length 208
void SVRRList::_svrr_10c0(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;
  data_[3] = 1.0;
  data_[4] = 1.0;
  data_[5] = 1.0;
  data_[6] = 1.0;
  data_[7] = 1.0;

  data_[8] = C00_[0];
  data_[9] = C00_[1];
  data_[10] = C00_[2];
  data_[11] = C00_[3];
  data_[12] = C00_[4];
  data_[13] = C00_[5];
  data_[14] = C00_[6];
  data_[15] = C00_[7];

  data_[16] = D00_[0];
  data_[17] = D00_[1];
  data_[18] = D00_[2];
  data_[19] = D00_[3];
  data_[20] = D00_[4];
  data_[21] = D00_[5];
  data_[22] = D00_[6];
  data_[23] = D00_[7];

  double cB00_current[8];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];
  cB00_current[4] = B00_[4];
  cB00_current[5] = B00_[5];
  cB00_current[6] = B00_[6];
  cB00_current[7] = B00_[7];

  data_[24] = C00_[0] * data_[16] + cB00_current[0];
  data_[25] = C00_[1] * data_[17] + cB00_current[1];
  data_[26] = C00_[2] * data_[18] + cB00_current[2];
  data_[27] = C00_[3] * data_[19] + cB00_current[3];
  data_[28] = C00_[4] * data_[20] + cB00_current[4];
  data_[29] = C00_[5] * data_[21] + cB00_current[5];
  data_[30] = C00_[6] * data_[22] + cB00_current[6];
  data_[31] = C00_[7] * data_[23] + cB00_current[7];

  double B01_current[8];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];
  B01_current[4] = B01_[4];
  B01_current[5] = B01_[5];
  B01_current[6] = B01_[6];
  B01_current[7] = B01_[7];

  data_[32] = D00_[0] * data_[16] + B01_current[0];
  data_[33] = D00_[1] * data_[17] + B01_current[1];
  data_[34] = D00_[2] * data_[18] + B01_current[2];
  data_[35] = D00_[3] * data_[19] + B01_current[3];
  data_[36] = D00_[4] * data_[20] + B01_current[4];
  data_[37] = D00_[5] * data_[21] + B01_current[5];
  data_[38] = D00_[6] * data_[22] + B01_current[6];
  data_[39] = D00_[7] * data_[23] + B01_current[7];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[40] = C00_[0] * data_[32] + cB00_current[0] * data_[16];
  data_[41] = C00_[1] * data_[33] + cB00_current[1] * data_[17];
  data_[42] = C00_[2] * data_[34] + cB00_current[2] * data_[18];
  data_[43] = C00_[3] * data_[35] + cB00_current[3] * data_[19];
  data_[44] = C00_[4] * data_[36] + cB00_current[4] * data_[20];
  data_[45] = C00_[5] * data_[37] + cB00_current[5] * data_[21];
  data_[46] = C00_[6] * data_[38] + cB00_current[6] * data_[22];
  data_[47] = C00_[7] * data_[39] + cB00_current[7] * data_[23];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[48] = D00_[0] * data_[32] + B01_current[0] * data_[16];
  data_[49] = D00_[1] * data_[33] + B01_current[1] * data_[17];
  data_[50] = D00_[2] * data_[34] + B01_current[2] * data_[18];
  data_[51] = D00_[3] * data_[35] + B01_current[3] * data_[19];
  data_[52] = D00_[4] * data_[36] + B01_current[4] * data_[20];
  data_[53] = D00_[5] * data_[37] + B01_current[5] * data_[21];
  data_[54] = D00_[6] * data_[38] + B01_current[6] * data_[22];
  data_[55] = D00_[7] * data_[39] + B01_current[7] * data_[23];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[56] = C00_[0] * data_[48] + cB00_current[0] * data_[32];
  data_[57] = C00_[1] * data_[49] + cB00_current[1] * data_[33];
  data_[58] = C00_[2] * data_[50] + cB00_current[2] * data_[34];
  data_[59] = C00_[3] * data_[51] + cB00_current[3] * data_[35];
  data_[60] = C00_[4] * data_[52] + cB00_current[4] * data_[36];
  data_[61] = C00_[5] * data_[53] + cB00_current[5] * data_[37];
  data_[62] = C00_[6] * data_[54] + cB00_current[6] * data_[38];
  data_[63] = C00_[7] * data_[55] + cB00_current[7] * data_[39];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[64] = D00_[0] * data_[48] + B01_current[0] * data_[32];
  data_[65] = D00_[1] * data_[49] + B01_current[1] * data_[33];
  data_[66] = D00_[2] * data_[50] + B01_current[2] * data_[34];
  data_[67] = D00_[3] * data_[51] + B01_current[3] * data_[35];
  data_[68] = D00_[4] * data_[52] + B01_current[4] * data_[36];
  data_[69] = D00_[5] * data_[53] + B01_current[5] * data_[37];
  data_[70] = D00_[6] * data_[54] + B01_current[6] * data_[38];
  data_[71] = D00_[7] * data_[55] + B01_current[7] * data_[39];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[72] = C00_[0] * data_[64] + cB00_current[0] * data_[48];
  data_[73] = C00_[1] * data_[65] + cB00_current[1] * data_[49];
  data_[74] = C00_[2] * data_[66] + cB00_current[2] * data_[50];
  data_[75] = C00_[3] * data_[67] + cB00_current[3] * data_[51];
  data_[76] = C00_[4] * data_[68] + cB00_current[4] * data_[52];
  data_[77] = C00_[5] * data_[69] + cB00_current[5] * data_[53];
  data_[78] = C00_[6] * data_[70] + cB00_current[6] * data_[54];
  data_[79] = C00_[7] * data_[71] + cB00_current[7] * data_[55];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[80] = D00_[0] * data_[64] + B01_current[0] * data_[48];
  data_[81] = D00_[1] * data_[65] + B01_current[1] * data_[49];
  data_[82] = D00_[2] * data_[66] + B01_current[2] * data_[50];
  data_[83] = D00_[3] * data_[67] + B01_current[3] * data_[51];
  data_[84] = D00_[4] * data_[68] + B01_current[4] * data_[52];
  data_[85] = D00_[5] * data_[69] + B01_current[5] * data_[53];
  data_[86] = D00_[6] * data_[70] + B01_current[6] * data_[54];
  data_[87] = D00_[7] * data_[71] + B01_current[7] * data_[55];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[88] = C00_[0] * data_[80] + cB00_current[0] * data_[64];
  data_[89] = C00_[1] * data_[81] + cB00_current[1] * data_[65];
  data_[90] = C00_[2] * data_[82] + cB00_current[2] * data_[66];
  data_[91] = C00_[3] * data_[83] + cB00_current[3] * data_[67];
  data_[92] = C00_[4] * data_[84] + cB00_current[4] * data_[68];
  data_[93] = C00_[5] * data_[85] + cB00_current[5] * data_[69];
  data_[94] = C00_[6] * data_[86] + cB00_current[6] * data_[70];
  data_[95] = C00_[7] * data_[87] + cB00_current[7] * data_[71];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[96] = D00_[0] * data_[80] + B01_current[0] * data_[64];
  data_[97] = D00_[1] * data_[81] + B01_current[1] * data_[65];
  data_[98] = D00_[2] * data_[82] + B01_current[2] * data_[66];
  data_[99] = D00_[3] * data_[83] + B01_current[3] * data_[67];
  data_[100] = D00_[4] * data_[84] + B01_current[4] * data_[68];
  data_[101] = D00_[5] * data_[85] + B01_current[5] * data_[69];
  data_[102] = D00_[6] * data_[86] + B01_current[6] * data_[70];
  data_[103] = D00_[7] * data_[87] + B01_current[7] * data_[71];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[104] = C00_[0] * data_[96] + cB00_current[0] * data_[80];
  data_[105] = C00_[1] * data_[97] + cB00_current[1] * data_[81];
  data_[106] = C00_[2] * data_[98] + cB00_current[2] * data_[82];
  data_[107] = C00_[3] * data_[99] + cB00_current[3] * data_[83];
  data_[108] = C00_[4] * data_[100] + cB00_current[4] * data_[84];
  data_[109] = C00_[5] * data_[101] + cB00_current[5] * data_[85];
  data_[110] = C00_[6] * data_[102] + cB00_current[6] * data_[86];
  data_[111] = C00_[7] * data_[103] + cB00_current[7] * data_[87];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[112] = D00_[0] * data_[96] + B01_current[0] * data_[80];
  data_[113] = D00_[1] * data_[97] + B01_current[1] * data_[81];
  data_[114] = D00_[2] * data_[98] + B01_current[2] * data_[82];
  data_[115] = D00_[3] * data_[99] + B01_current[3] * data_[83];
  data_[116] = D00_[4] * data_[100] + B01_current[4] * data_[84];
  data_[117] = D00_[5] * data_[101] + B01_current[5] * data_[85];
  data_[118] = D00_[6] * data_[102] + B01_current[6] * data_[86];
  data_[119] = D00_[7] * data_[103] + B01_current[7] * data_[87];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[120] = C00_[0] * data_[112] + cB00_current[0] * data_[96];
  data_[121] = C00_[1] * data_[113] + cB00_current[1] * data_[97];
  data_[122] = C00_[2] * data_[114] + cB00_current[2] * data_[98];
  data_[123] = C00_[3] * data_[115] + cB00_current[3] * data_[99];
  data_[124] = C00_[4] * data_[116] + cB00_current[4] * data_[100];
  data_[125] = C00_[5] * data_[117] + cB00_current[5] * data_[101];
  data_[126] = C00_[6] * data_[118] + cB00_current[6] * data_[102];
  data_[127] = C00_[7] * data_[119] + cB00_current[7] * data_[103];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[128] = D00_[0] * data_[112] + B01_current[0] * data_[96];
  data_[129] = D00_[1] * data_[113] + B01_current[1] * data_[97];
  data_[130] = D00_[2] * data_[114] + B01_current[2] * data_[98];
  data_[131] = D00_[3] * data_[115] + B01_current[3] * data_[99];
  data_[132] = D00_[4] * data_[116] + B01_current[4] * data_[100];
  data_[133] = D00_[5] * data_[117] + B01_current[5] * data_[101];
  data_[134] = D00_[6] * data_[118] + B01_current[6] * data_[102];
  data_[135] = D00_[7] * data_[119] + B01_current[7] * data_[103];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[136] = C00_[0] * data_[128] + cB00_current[0] * data_[112];
  data_[137] = C00_[1] * data_[129] + cB00_current[1] * data_[113];
  data_[138] = C00_[2] * data_[130] + cB00_current[2] * data_[114];
  data_[139] = C00_[3] * data_[131] + cB00_current[3] * data_[115];
  data_[140] = C00_[4] * data_[132] + cB00_current[4] * data_[116];
  data_[141] = C00_[5] * data_[133] + cB00_current[5] * data_[117];
  data_[142] = C00_[6] * data_[134] + cB00_current[6] * data_[118];
  data_[143] = C00_[7] * data_[135] + cB00_current[7] * data_[119];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[144] = D00_[0] * data_[128] + B01_current[0] * data_[112];
  data_[145] = D00_[1] * data_[129] + B01_current[1] * data_[113];
  data_[146] = D00_[2] * data_[130] + B01_current[2] * data_[114];
  data_[147] = D00_[3] * data_[131] + B01_current[3] * data_[115];
  data_[148] = D00_[4] * data_[132] + B01_current[4] * data_[116];
  data_[149] = D00_[5] * data_[133] + B01_current[5] * data_[117];
  data_[150] = D00_[6] * data_[134] + B01_current[6] * data_[118];
  data_[151] = D00_[7] * data_[135] + B01_current[7] * data_[119];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[152] = C00_[0] * data_[144] + cB00_current[0] * data_[128];
  data_[153] = C00_[1] * data_[145] + cB00_current[1] * data_[129];
  data_[154] = C00_[2] * data_[146] + cB00_current[2] * data_[130];
  data_[155] = C00_[3] * data_[147] + cB00_current[3] * data_[131];
  data_[156] = C00_[4] * data_[148] + cB00_current[4] * data_[132];
  data_[157] = C00_[5] * data_[149] + cB00_current[5] * data_[133];
  data_[158] = C00_[6] * data_[150] + cB00_current[6] * data_[134];
  data_[159] = C00_[7] * data_[151] + cB00_current[7] * data_[135];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[160] = D00_[0] * data_[144] + B01_current[0] * data_[128];
  data_[161] = D00_[1] * data_[145] + B01_current[1] * data_[129];
  data_[162] = D00_[2] * data_[146] + B01_current[2] * data_[130];
  data_[163] = D00_[3] * data_[147] + B01_current[3] * data_[131];
  data_[164] = D00_[4] * data_[148] + B01_current[4] * data_[132];
  data_[165] = D00_[5] * data_[149] + B01_current[5] * data_[133];
  data_[166] = D00_[6] * data_[150] + B01_current[6] * data_[134];
  data_[167] = D00_[7] * data_[151] + B01_current[7] * data_[135];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[168] = C00_[0] * data_[160] + cB00_current[0] * data_[144];
  data_[169] = C00_[1] * data_[161] + cB00_current[1] * data_[145];
  data_[170] = C00_[2] * data_[162] + cB00_current[2] * data_[146];
  data_[171] = C00_[3] * data_[163] + cB00_current[3] * data_[147];
  data_[172] = C00_[4] * data_[164] + cB00_current[4] * data_[148];
  data_[173] = C00_[5] * data_[165] + cB00_current[5] * data_[149];
  data_[174] = C00_[6] * data_[166] + cB00_current[6] * data_[150];
  data_[175] = C00_[7] * data_[167] + cB00_current[7] * data_[151];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[176] = D00_[0] * data_[160] + B01_current[0] * data_[144];
  data_[177] = D00_[1] * data_[161] + B01_current[1] * data_[145];
  data_[178] = D00_[2] * data_[162] + B01_current[2] * data_[146];
  data_[179] = D00_[3] * data_[163] + B01_current[3] * data_[147];
  data_[180] = D00_[4] * data_[164] + B01_current[4] * data_[148];
  data_[181] = D00_[5] * data_[165] + B01_current[5] * data_[149];
  data_[182] = D00_[6] * data_[166] + B01_current[6] * data_[150];
  data_[183] = D00_[7] * data_[167] + B01_current[7] * data_[151];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[184] = C00_[0] * data_[176] + cB00_current[0] * data_[160];
  data_[185] = C00_[1] * data_[177] + cB00_current[1] * data_[161];
  data_[186] = C00_[2] * data_[178] + cB00_current[2] * data_[162];
  data_[187] = C00_[3] * data_[179] + cB00_current[3] * data_[163];
  data_[188] = C00_[4] * data_[180] + cB00_current[4] * data_[164];
  data_[189] = C00_[5] * data_[181] + cB00_current[5] * data_[165];
  data_[190] = C00_[6] * data_[182] + cB00_current[6] * data_[166];
  data_[191] = C00_[7] * data_[183] + cB00_current[7] * data_[167];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[192] = D00_[0] * data_[176] + B01_current[0] * data_[160];
  data_[193] = D00_[1] * data_[177] + B01_current[1] * data_[161];
  data_[194] = D00_[2] * data_[178] + B01_current[2] * data_[162];
  data_[195] = D00_[3] * data_[179] + B01_current[3] * data_[163];
  data_[196] = D00_[4] * data_[180] + B01_current[4] * data_[164];
  data_[197] = D00_[5] * data_[181] + B01_current[5] * data_[165];
  data_[198] = D00_[6] * data_[182] + B01_current[6] * data_[166];
  data_[199] = D00_[7] * data_[183] + B01_current[7] * data_[167];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[200] = C00_[0] * data_[192] + cB00_current[0] * data_[176];
  data_[201] = C00_[1] * data_[193] + cB00_current[1] * data_[177];
  data_[202] = C00_[2] * data_[194] + cB00_current[2] * data_[178];
  data_[203] = C00_[3] * data_[195] + cB00_current[3] * data_[179];
  data_[204] = C00_[4] * data_[196] + cB00_current[4] * data_[180];
  data_[205] = C00_[5] * data_[197] + cB00_current[5] * data_[181];
  data_[206] = C00_[6] * data_[198] + cB00_current[6] * data_[182];
  data_[207] = C00_[7] * data_[199] + cB00_current[7] * data_[183];
}

