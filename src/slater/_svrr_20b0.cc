//
// Newint - Parallel electron correlation program.
// Filename: _svrr_20b0.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/slater/svrrlist.h>

// returns double array of length 288
void SVRRList::_svrr_20b0(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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

  data_[16] = C00_[0] * data_[8] + B10_[0];
  data_[17] = C00_[1] * data_[9] + B10_[1];
  data_[18] = C00_[2] * data_[10] + B10_[2];
  data_[19] = C00_[3] * data_[11] + B10_[3];
  data_[20] = C00_[4] * data_[12] + B10_[4];
  data_[21] = C00_[5] * data_[13] + B10_[5];
  data_[22] = C00_[6] * data_[14] + B10_[6];
  data_[23] = C00_[7] * data_[15] + B10_[7];

  data_[24] = D00_[0];
  data_[25] = D00_[1];
  data_[26] = D00_[2];
  data_[27] = D00_[3];
  data_[28] = D00_[4];
  data_[29] = D00_[5];
  data_[30] = D00_[6];
  data_[31] = D00_[7];

  double cB00_current[8];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];
  cB00_current[4] = B00_[4];
  cB00_current[5] = B00_[5];
  cB00_current[6] = B00_[6];
  cB00_current[7] = B00_[7];

  data_[32] = C00_[0] * data_[24] + cB00_current[0];
  data_[33] = C00_[1] * data_[25] + cB00_current[1];
  data_[34] = C00_[2] * data_[26] + cB00_current[2];
  data_[35] = C00_[3] * data_[27] + cB00_current[3];
  data_[36] = C00_[4] * data_[28] + cB00_current[4];
  data_[37] = C00_[5] * data_[29] + cB00_current[5];
  data_[38] = C00_[6] * data_[30] + cB00_current[6];
  data_[39] = C00_[7] * data_[31] + cB00_current[7];

  data_[40] = C00_[0] * data_[32] + B10_[0] * data_[24] + cB00_current[0] * data_[8];
  data_[41] = C00_[1] * data_[33] + B10_[1] * data_[25] + cB00_current[1] * data_[9];
  data_[42] = C00_[2] * data_[34] + B10_[2] * data_[26] + cB00_current[2] * data_[10];
  data_[43] = C00_[3] * data_[35] + B10_[3] * data_[27] + cB00_current[3] * data_[11];
  data_[44] = C00_[4] * data_[36] + B10_[4] * data_[28] + cB00_current[4] * data_[12];
  data_[45] = C00_[5] * data_[37] + B10_[5] * data_[29] + cB00_current[5] * data_[13];
  data_[46] = C00_[6] * data_[38] + B10_[6] * data_[30] + cB00_current[6] * data_[14];
  data_[47] = C00_[7] * data_[39] + B10_[7] * data_[31] + cB00_current[7] * data_[15];

  double B01_current[8];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];
  B01_current[4] = B01_[4];
  B01_current[5] = B01_[5];
  B01_current[6] = B01_[6];
  B01_current[7] = B01_[7];

  data_[48] = D00_[0] * data_[24] + B01_current[0];
  data_[49] = D00_[1] * data_[25] + B01_current[1];
  data_[50] = D00_[2] * data_[26] + B01_current[2];
  data_[51] = D00_[3] * data_[27] + B01_current[3];
  data_[52] = D00_[4] * data_[28] + B01_current[4];
  data_[53] = D00_[5] * data_[29] + B01_current[5];
  data_[54] = D00_[6] * data_[30] + B01_current[6];
  data_[55] = D00_[7] * data_[31] + B01_current[7];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[56] = C00_[0] * data_[48] + cB00_current[0] * data_[24];
  data_[57] = C00_[1] * data_[49] + cB00_current[1] * data_[25];
  data_[58] = C00_[2] * data_[50] + cB00_current[2] * data_[26];
  data_[59] = C00_[3] * data_[51] + cB00_current[3] * data_[27];
  data_[60] = C00_[4] * data_[52] + cB00_current[4] * data_[28];
  data_[61] = C00_[5] * data_[53] + cB00_current[5] * data_[29];
  data_[62] = C00_[6] * data_[54] + cB00_current[6] * data_[30];
  data_[63] = C00_[7] * data_[55] + cB00_current[7] * data_[31];

  data_[64] = C00_[0] * data_[56] + B10_[0] * data_[48] + cB00_current[0] * data_[32];
  data_[65] = C00_[1] * data_[57] + B10_[1] * data_[49] + cB00_current[1] * data_[33];
  data_[66] = C00_[2] * data_[58] + B10_[2] * data_[50] + cB00_current[2] * data_[34];
  data_[67] = C00_[3] * data_[59] + B10_[3] * data_[51] + cB00_current[3] * data_[35];
  data_[68] = C00_[4] * data_[60] + B10_[4] * data_[52] + cB00_current[4] * data_[36];
  data_[69] = C00_[5] * data_[61] + B10_[5] * data_[53] + cB00_current[5] * data_[37];
  data_[70] = C00_[6] * data_[62] + B10_[6] * data_[54] + cB00_current[6] * data_[38];
  data_[71] = C00_[7] * data_[63] + B10_[7] * data_[55] + cB00_current[7] * data_[39];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[72] = D00_[0] * data_[48] + B01_current[0] * data_[24];
  data_[73] = D00_[1] * data_[49] + B01_current[1] * data_[25];
  data_[74] = D00_[2] * data_[50] + B01_current[2] * data_[26];
  data_[75] = D00_[3] * data_[51] + B01_current[3] * data_[27];
  data_[76] = D00_[4] * data_[52] + B01_current[4] * data_[28];
  data_[77] = D00_[5] * data_[53] + B01_current[5] * data_[29];
  data_[78] = D00_[6] * data_[54] + B01_current[6] * data_[30];
  data_[79] = D00_[7] * data_[55] + B01_current[7] * data_[31];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[80] = C00_[0] * data_[72] + cB00_current[0] * data_[48];
  data_[81] = C00_[1] * data_[73] + cB00_current[1] * data_[49];
  data_[82] = C00_[2] * data_[74] + cB00_current[2] * data_[50];
  data_[83] = C00_[3] * data_[75] + cB00_current[3] * data_[51];
  data_[84] = C00_[4] * data_[76] + cB00_current[4] * data_[52];
  data_[85] = C00_[5] * data_[77] + cB00_current[5] * data_[53];
  data_[86] = C00_[6] * data_[78] + cB00_current[6] * data_[54];
  data_[87] = C00_[7] * data_[79] + cB00_current[7] * data_[55];

  data_[88] = C00_[0] * data_[80] + B10_[0] * data_[72] + cB00_current[0] * data_[56];
  data_[89] = C00_[1] * data_[81] + B10_[1] * data_[73] + cB00_current[1] * data_[57];
  data_[90] = C00_[2] * data_[82] + B10_[2] * data_[74] + cB00_current[2] * data_[58];
  data_[91] = C00_[3] * data_[83] + B10_[3] * data_[75] + cB00_current[3] * data_[59];
  data_[92] = C00_[4] * data_[84] + B10_[4] * data_[76] + cB00_current[4] * data_[60];
  data_[93] = C00_[5] * data_[85] + B10_[5] * data_[77] + cB00_current[5] * data_[61];
  data_[94] = C00_[6] * data_[86] + B10_[6] * data_[78] + cB00_current[6] * data_[62];
  data_[95] = C00_[7] * data_[87] + B10_[7] * data_[79] + cB00_current[7] * data_[63];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[96] = D00_[0] * data_[72] + B01_current[0] * data_[48];
  data_[97] = D00_[1] * data_[73] + B01_current[1] * data_[49];
  data_[98] = D00_[2] * data_[74] + B01_current[2] * data_[50];
  data_[99] = D00_[3] * data_[75] + B01_current[3] * data_[51];
  data_[100] = D00_[4] * data_[76] + B01_current[4] * data_[52];
  data_[101] = D00_[5] * data_[77] + B01_current[5] * data_[53];
  data_[102] = D00_[6] * data_[78] + B01_current[6] * data_[54];
  data_[103] = D00_[7] * data_[79] + B01_current[7] * data_[55];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[104] = C00_[0] * data_[96] + cB00_current[0] * data_[72];
  data_[105] = C00_[1] * data_[97] + cB00_current[1] * data_[73];
  data_[106] = C00_[2] * data_[98] + cB00_current[2] * data_[74];
  data_[107] = C00_[3] * data_[99] + cB00_current[3] * data_[75];
  data_[108] = C00_[4] * data_[100] + cB00_current[4] * data_[76];
  data_[109] = C00_[5] * data_[101] + cB00_current[5] * data_[77];
  data_[110] = C00_[6] * data_[102] + cB00_current[6] * data_[78];
  data_[111] = C00_[7] * data_[103] + cB00_current[7] * data_[79];

  data_[112] = C00_[0] * data_[104] + B10_[0] * data_[96] + cB00_current[0] * data_[80];
  data_[113] = C00_[1] * data_[105] + B10_[1] * data_[97] + cB00_current[1] * data_[81];
  data_[114] = C00_[2] * data_[106] + B10_[2] * data_[98] + cB00_current[2] * data_[82];
  data_[115] = C00_[3] * data_[107] + B10_[3] * data_[99] + cB00_current[3] * data_[83];
  data_[116] = C00_[4] * data_[108] + B10_[4] * data_[100] + cB00_current[4] * data_[84];
  data_[117] = C00_[5] * data_[109] + B10_[5] * data_[101] + cB00_current[5] * data_[85];
  data_[118] = C00_[6] * data_[110] + B10_[6] * data_[102] + cB00_current[6] * data_[86];
  data_[119] = C00_[7] * data_[111] + B10_[7] * data_[103] + cB00_current[7] * data_[87];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[120] = D00_[0] * data_[96] + B01_current[0] * data_[72];
  data_[121] = D00_[1] * data_[97] + B01_current[1] * data_[73];
  data_[122] = D00_[2] * data_[98] + B01_current[2] * data_[74];
  data_[123] = D00_[3] * data_[99] + B01_current[3] * data_[75];
  data_[124] = D00_[4] * data_[100] + B01_current[4] * data_[76];
  data_[125] = D00_[5] * data_[101] + B01_current[5] * data_[77];
  data_[126] = D00_[6] * data_[102] + B01_current[6] * data_[78];
  data_[127] = D00_[7] * data_[103] + B01_current[7] * data_[79];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[128] = C00_[0] * data_[120] + cB00_current[0] * data_[96];
  data_[129] = C00_[1] * data_[121] + cB00_current[1] * data_[97];
  data_[130] = C00_[2] * data_[122] + cB00_current[2] * data_[98];
  data_[131] = C00_[3] * data_[123] + cB00_current[3] * data_[99];
  data_[132] = C00_[4] * data_[124] + cB00_current[4] * data_[100];
  data_[133] = C00_[5] * data_[125] + cB00_current[5] * data_[101];
  data_[134] = C00_[6] * data_[126] + cB00_current[6] * data_[102];
  data_[135] = C00_[7] * data_[127] + cB00_current[7] * data_[103];

  data_[136] = C00_[0] * data_[128] + B10_[0] * data_[120] + cB00_current[0] * data_[104];
  data_[137] = C00_[1] * data_[129] + B10_[1] * data_[121] + cB00_current[1] * data_[105];
  data_[138] = C00_[2] * data_[130] + B10_[2] * data_[122] + cB00_current[2] * data_[106];
  data_[139] = C00_[3] * data_[131] + B10_[3] * data_[123] + cB00_current[3] * data_[107];
  data_[140] = C00_[4] * data_[132] + B10_[4] * data_[124] + cB00_current[4] * data_[108];
  data_[141] = C00_[5] * data_[133] + B10_[5] * data_[125] + cB00_current[5] * data_[109];
  data_[142] = C00_[6] * data_[134] + B10_[6] * data_[126] + cB00_current[6] * data_[110];
  data_[143] = C00_[7] * data_[135] + B10_[7] * data_[127] + cB00_current[7] * data_[111];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[144] = D00_[0] * data_[120] + B01_current[0] * data_[96];
  data_[145] = D00_[1] * data_[121] + B01_current[1] * data_[97];
  data_[146] = D00_[2] * data_[122] + B01_current[2] * data_[98];
  data_[147] = D00_[3] * data_[123] + B01_current[3] * data_[99];
  data_[148] = D00_[4] * data_[124] + B01_current[4] * data_[100];
  data_[149] = D00_[5] * data_[125] + B01_current[5] * data_[101];
  data_[150] = D00_[6] * data_[126] + B01_current[6] * data_[102];
  data_[151] = D00_[7] * data_[127] + B01_current[7] * data_[103];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[152] = C00_[0] * data_[144] + cB00_current[0] * data_[120];
  data_[153] = C00_[1] * data_[145] + cB00_current[1] * data_[121];
  data_[154] = C00_[2] * data_[146] + cB00_current[2] * data_[122];
  data_[155] = C00_[3] * data_[147] + cB00_current[3] * data_[123];
  data_[156] = C00_[4] * data_[148] + cB00_current[4] * data_[124];
  data_[157] = C00_[5] * data_[149] + cB00_current[5] * data_[125];
  data_[158] = C00_[6] * data_[150] + cB00_current[6] * data_[126];
  data_[159] = C00_[7] * data_[151] + cB00_current[7] * data_[127];

  data_[160] = C00_[0] * data_[152] + B10_[0] * data_[144] + cB00_current[0] * data_[128];
  data_[161] = C00_[1] * data_[153] + B10_[1] * data_[145] + cB00_current[1] * data_[129];
  data_[162] = C00_[2] * data_[154] + B10_[2] * data_[146] + cB00_current[2] * data_[130];
  data_[163] = C00_[3] * data_[155] + B10_[3] * data_[147] + cB00_current[3] * data_[131];
  data_[164] = C00_[4] * data_[156] + B10_[4] * data_[148] + cB00_current[4] * data_[132];
  data_[165] = C00_[5] * data_[157] + B10_[5] * data_[149] + cB00_current[5] * data_[133];
  data_[166] = C00_[6] * data_[158] + B10_[6] * data_[150] + cB00_current[6] * data_[134];
  data_[167] = C00_[7] * data_[159] + B10_[7] * data_[151] + cB00_current[7] * data_[135];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[168] = D00_[0] * data_[144] + B01_current[0] * data_[120];
  data_[169] = D00_[1] * data_[145] + B01_current[1] * data_[121];
  data_[170] = D00_[2] * data_[146] + B01_current[2] * data_[122];
  data_[171] = D00_[3] * data_[147] + B01_current[3] * data_[123];
  data_[172] = D00_[4] * data_[148] + B01_current[4] * data_[124];
  data_[173] = D00_[5] * data_[149] + B01_current[5] * data_[125];
  data_[174] = D00_[6] * data_[150] + B01_current[6] * data_[126];
  data_[175] = D00_[7] * data_[151] + B01_current[7] * data_[127];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[176] = C00_[0] * data_[168] + cB00_current[0] * data_[144];
  data_[177] = C00_[1] * data_[169] + cB00_current[1] * data_[145];
  data_[178] = C00_[2] * data_[170] + cB00_current[2] * data_[146];
  data_[179] = C00_[3] * data_[171] + cB00_current[3] * data_[147];
  data_[180] = C00_[4] * data_[172] + cB00_current[4] * data_[148];
  data_[181] = C00_[5] * data_[173] + cB00_current[5] * data_[149];
  data_[182] = C00_[6] * data_[174] + cB00_current[6] * data_[150];
  data_[183] = C00_[7] * data_[175] + cB00_current[7] * data_[151];

  data_[184] = C00_[0] * data_[176] + B10_[0] * data_[168] + cB00_current[0] * data_[152];
  data_[185] = C00_[1] * data_[177] + B10_[1] * data_[169] + cB00_current[1] * data_[153];
  data_[186] = C00_[2] * data_[178] + B10_[2] * data_[170] + cB00_current[2] * data_[154];
  data_[187] = C00_[3] * data_[179] + B10_[3] * data_[171] + cB00_current[3] * data_[155];
  data_[188] = C00_[4] * data_[180] + B10_[4] * data_[172] + cB00_current[4] * data_[156];
  data_[189] = C00_[5] * data_[181] + B10_[5] * data_[173] + cB00_current[5] * data_[157];
  data_[190] = C00_[6] * data_[182] + B10_[6] * data_[174] + cB00_current[6] * data_[158];
  data_[191] = C00_[7] * data_[183] + B10_[7] * data_[175] + cB00_current[7] * data_[159];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[192] = D00_[0] * data_[168] + B01_current[0] * data_[144];
  data_[193] = D00_[1] * data_[169] + B01_current[1] * data_[145];
  data_[194] = D00_[2] * data_[170] + B01_current[2] * data_[146];
  data_[195] = D00_[3] * data_[171] + B01_current[3] * data_[147];
  data_[196] = D00_[4] * data_[172] + B01_current[4] * data_[148];
  data_[197] = D00_[5] * data_[173] + B01_current[5] * data_[149];
  data_[198] = D00_[6] * data_[174] + B01_current[6] * data_[150];
  data_[199] = D00_[7] * data_[175] + B01_current[7] * data_[151];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[200] = C00_[0] * data_[192] + cB00_current[0] * data_[168];
  data_[201] = C00_[1] * data_[193] + cB00_current[1] * data_[169];
  data_[202] = C00_[2] * data_[194] + cB00_current[2] * data_[170];
  data_[203] = C00_[3] * data_[195] + cB00_current[3] * data_[171];
  data_[204] = C00_[4] * data_[196] + cB00_current[4] * data_[172];
  data_[205] = C00_[5] * data_[197] + cB00_current[5] * data_[173];
  data_[206] = C00_[6] * data_[198] + cB00_current[6] * data_[174];
  data_[207] = C00_[7] * data_[199] + cB00_current[7] * data_[175];

  data_[208] = C00_[0] * data_[200] + B10_[0] * data_[192] + cB00_current[0] * data_[176];
  data_[209] = C00_[1] * data_[201] + B10_[1] * data_[193] + cB00_current[1] * data_[177];
  data_[210] = C00_[2] * data_[202] + B10_[2] * data_[194] + cB00_current[2] * data_[178];
  data_[211] = C00_[3] * data_[203] + B10_[3] * data_[195] + cB00_current[3] * data_[179];
  data_[212] = C00_[4] * data_[204] + B10_[4] * data_[196] + cB00_current[4] * data_[180];
  data_[213] = C00_[5] * data_[205] + B10_[5] * data_[197] + cB00_current[5] * data_[181];
  data_[214] = C00_[6] * data_[206] + B10_[6] * data_[198] + cB00_current[6] * data_[182];
  data_[215] = C00_[7] * data_[207] + B10_[7] * data_[199] + cB00_current[7] * data_[183];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[216] = D00_[0] * data_[192] + B01_current[0] * data_[168];
  data_[217] = D00_[1] * data_[193] + B01_current[1] * data_[169];
  data_[218] = D00_[2] * data_[194] + B01_current[2] * data_[170];
  data_[219] = D00_[3] * data_[195] + B01_current[3] * data_[171];
  data_[220] = D00_[4] * data_[196] + B01_current[4] * data_[172];
  data_[221] = D00_[5] * data_[197] + B01_current[5] * data_[173];
  data_[222] = D00_[6] * data_[198] + B01_current[6] * data_[174];
  data_[223] = D00_[7] * data_[199] + B01_current[7] * data_[175];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[224] = C00_[0] * data_[216] + cB00_current[0] * data_[192];
  data_[225] = C00_[1] * data_[217] + cB00_current[1] * data_[193];
  data_[226] = C00_[2] * data_[218] + cB00_current[2] * data_[194];
  data_[227] = C00_[3] * data_[219] + cB00_current[3] * data_[195];
  data_[228] = C00_[4] * data_[220] + cB00_current[4] * data_[196];
  data_[229] = C00_[5] * data_[221] + cB00_current[5] * data_[197];
  data_[230] = C00_[6] * data_[222] + cB00_current[6] * data_[198];
  data_[231] = C00_[7] * data_[223] + cB00_current[7] * data_[199];

  data_[232] = C00_[0] * data_[224] + B10_[0] * data_[216] + cB00_current[0] * data_[200];
  data_[233] = C00_[1] * data_[225] + B10_[1] * data_[217] + cB00_current[1] * data_[201];
  data_[234] = C00_[2] * data_[226] + B10_[2] * data_[218] + cB00_current[2] * data_[202];
  data_[235] = C00_[3] * data_[227] + B10_[3] * data_[219] + cB00_current[3] * data_[203];
  data_[236] = C00_[4] * data_[228] + B10_[4] * data_[220] + cB00_current[4] * data_[204];
  data_[237] = C00_[5] * data_[229] + B10_[5] * data_[221] + cB00_current[5] * data_[205];
  data_[238] = C00_[6] * data_[230] + B10_[6] * data_[222] + cB00_current[6] * data_[206];
  data_[239] = C00_[7] * data_[231] + B10_[7] * data_[223] + cB00_current[7] * data_[207];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[240] = D00_[0] * data_[216] + B01_current[0] * data_[192];
  data_[241] = D00_[1] * data_[217] + B01_current[1] * data_[193];
  data_[242] = D00_[2] * data_[218] + B01_current[2] * data_[194];
  data_[243] = D00_[3] * data_[219] + B01_current[3] * data_[195];
  data_[244] = D00_[4] * data_[220] + B01_current[4] * data_[196];
  data_[245] = D00_[5] * data_[221] + B01_current[5] * data_[197];
  data_[246] = D00_[6] * data_[222] + B01_current[6] * data_[198];
  data_[247] = D00_[7] * data_[223] + B01_current[7] * data_[199];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[248] = C00_[0] * data_[240] + cB00_current[0] * data_[216];
  data_[249] = C00_[1] * data_[241] + cB00_current[1] * data_[217];
  data_[250] = C00_[2] * data_[242] + cB00_current[2] * data_[218];
  data_[251] = C00_[3] * data_[243] + cB00_current[3] * data_[219];
  data_[252] = C00_[4] * data_[244] + cB00_current[4] * data_[220];
  data_[253] = C00_[5] * data_[245] + cB00_current[5] * data_[221];
  data_[254] = C00_[6] * data_[246] + cB00_current[6] * data_[222];
  data_[255] = C00_[7] * data_[247] + cB00_current[7] * data_[223];

  data_[256] = C00_[0] * data_[248] + B10_[0] * data_[240] + cB00_current[0] * data_[224];
  data_[257] = C00_[1] * data_[249] + B10_[1] * data_[241] + cB00_current[1] * data_[225];
  data_[258] = C00_[2] * data_[250] + B10_[2] * data_[242] + cB00_current[2] * data_[226];
  data_[259] = C00_[3] * data_[251] + B10_[3] * data_[243] + cB00_current[3] * data_[227];
  data_[260] = C00_[4] * data_[252] + B10_[4] * data_[244] + cB00_current[4] * data_[228];
  data_[261] = C00_[5] * data_[253] + B10_[5] * data_[245] + cB00_current[5] * data_[229];
  data_[262] = C00_[6] * data_[254] + B10_[6] * data_[246] + cB00_current[6] * data_[230];
  data_[263] = C00_[7] * data_[255] + B10_[7] * data_[247] + cB00_current[7] * data_[231];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[264] = D00_[0] * data_[240] + B01_current[0] * data_[216];
  data_[265] = D00_[1] * data_[241] + B01_current[1] * data_[217];
  data_[266] = D00_[2] * data_[242] + B01_current[2] * data_[218];
  data_[267] = D00_[3] * data_[243] + B01_current[3] * data_[219];
  data_[268] = D00_[4] * data_[244] + B01_current[4] * data_[220];
  data_[269] = D00_[5] * data_[245] + B01_current[5] * data_[221];
  data_[270] = D00_[6] * data_[246] + B01_current[6] * data_[222];
  data_[271] = D00_[7] * data_[247] + B01_current[7] * data_[223];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[272] = C00_[0] * data_[264] + cB00_current[0] * data_[240];
  data_[273] = C00_[1] * data_[265] + cB00_current[1] * data_[241];
  data_[274] = C00_[2] * data_[266] + cB00_current[2] * data_[242];
  data_[275] = C00_[3] * data_[267] + cB00_current[3] * data_[243];
  data_[276] = C00_[4] * data_[268] + cB00_current[4] * data_[244];
  data_[277] = C00_[5] * data_[269] + cB00_current[5] * data_[245];
  data_[278] = C00_[6] * data_[270] + cB00_current[6] * data_[246];
  data_[279] = C00_[7] * data_[271] + cB00_current[7] * data_[247];

  data_[280] = C00_[0] * data_[272] + B10_[0] * data_[264] + cB00_current[0] * data_[248];
  data_[281] = C00_[1] * data_[273] + B10_[1] * data_[265] + cB00_current[1] * data_[249];
  data_[282] = C00_[2] * data_[274] + B10_[2] * data_[266] + cB00_current[2] * data_[250];
  data_[283] = C00_[3] * data_[275] + B10_[3] * data_[267] + cB00_current[3] * data_[251];
  data_[284] = C00_[4] * data_[276] + B10_[4] * data_[268] + cB00_current[4] * data_[252];
  data_[285] = C00_[5] * data_[277] + B10_[5] * data_[269] + cB00_current[5] * data_[253];
  data_[286] = C00_[6] * data_[278] + B10_[6] * data_[270] + cB00_current[6] * data_[254];
  data_[287] = C00_[7] * data_[279] + B10_[7] * data_[271] + cB00_current[7] * data_[255];
}

