//
// Newint - Parallel electron correlation program.
// Filename: _svrr_8050.cc
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

// returns double array of length 432
void SVRRList::_svrr_8050(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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

  double B10_current[8];
  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[16] = C00_[0] * data_[8] + B10_current[0];
  data_[17] = C00_[1] * data_[9] + B10_current[1];
  data_[18] = C00_[2] * data_[10] + B10_current[2];
  data_[19] = C00_[3] * data_[11] + B10_current[3];
  data_[20] = C00_[4] * data_[12] + B10_current[4];
  data_[21] = C00_[5] * data_[13] + B10_current[5];
  data_[22] = C00_[6] * data_[14] + B10_current[6];
  data_[23] = C00_[7] * data_[15] + B10_current[7];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[24] = C00_[0] * data_[16] + B10_current[0] * data_[8];
  data_[25] = C00_[1] * data_[17] + B10_current[1] * data_[9];
  data_[26] = C00_[2] * data_[18] + B10_current[2] * data_[10];
  data_[27] = C00_[3] * data_[19] + B10_current[3] * data_[11];
  data_[28] = C00_[4] * data_[20] + B10_current[4] * data_[12];
  data_[29] = C00_[5] * data_[21] + B10_current[5] * data_[13];
  data_[30] = C00_[6] * data_[22] + B10_current[6] * data_[14];
  data_[31] = C00_[7] * data_[23] + B10_current[7] * data_[15];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[32] = C00_[0] * data_[24] + B10_current[0] * data_[16];
  data_[33] = C00_[1] * data_[25] + B10_current[1] * data_[17];
  data_[34] = C00_[2] * data_[26] + B10_current[2] * data_[18];
  data_[35] = C00_[3] * data_[27] + B10_current[3] * data_[19];
  data_[36] = C00_[4] * data_[28] + B10_current[4] * data_[20];
  data_[37] = C00_[5] * data_[29] + B10_current[5] * data_[21];
  data_[38] = C00_[6] * data_[30] + B10_current[6] * data_[22];
  data_[39] = C00_[7] * data_[31] + B10_current[7] * data_[23];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[40] = C00_[0] * data_[32] + B10_current[0] * data_[24];
  data_[41] = C00_[1] * data_[33] + B10_current[1] * data_[25];
  data_[42] = C00_[2] * data_[34] + B10_current[2] * data_[26];
  data_[43] = C00_[3] * data_[35] + B10_current[3] * data_[27];
  data_[44] = C00_[4] * data_[36] + B10_current[4] * data_[28];
  data_[45] = C00_[5] * data_[37] + B10_current[5] * data_[29];
  data_[46] = C00_[6] * data_[38] + B10_current[6] * data_[30];
  data_[47] = C00_[7] * data_[39] + B10_current[7] * data_[31];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[48] = C00_[0] * data_[40] + B10_current[0] * data_[32];
  data_[49] = C00_[1] * data_[41] + B10_current[1] * data_[33];
  data_[50] = C00_[2] * data_[42] + B10_current[2] * data_[34];
  data_[51] = C00_[3] * data_[43] + B10_current[3] * data_[35];
  data_[52] = C00_[4] * data_[44] + B10_current[4] * data_[36];
  data_[53] = C00_[5] * data_[45] + B10_current[5] * data_[37];
  data_[54] = C00_[6] * data_[46] + B10_current[6] * data_[38];
  data_[55] = C00_[7] * data_[47] + B10_current[7] * data_[39];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[56] = C00_[0] * data_[48] + B10_current[0] * data_[40];
  data_[57] = C00_[1] * data_[49] + B10_current[1] * data_[41];
  data_[58] = C00_[2] * data_[50] + B10_current[2] * data_[42];
  data_[59] = C00_[3] * data_[51] + B10_current[3] * data_[43];
  data_[60] = C00_[4] * data_[52] + B10_current[4] * data_[44];
  data_[61] = C00_[5] * data_[53] + B10_current[5] * data_[45];
  data_[62] = C00_[6] * data_[54] + B10_current[6] * data_[46];
  data_[63] = C00_[7] * data_[55] + B10_current[7] * data_[47];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[64] = C00_[0] * data_[56] + B10_current[0] * data_[48];
  data_[65] = C00_[1] * data_[57] + B10_current[1] * data_[49];
  data_[66] = C00_[2] * data_[58] + B10_current[2] * data_[50];
  data_[67] = C00_[3] * data_[59] + B10_current[3] * data_[51];
  data_[68] = C00_[4] * data_[60] + B10_current[4] * data_[52];
  data_[69] = C00_[5] * data_[61] + B10_current[5] * data_[53];
  data_[70] = C00_[6] * data_[62] + B10_current[6] * data_[54];
  data_[71] = C00_[7] * data_[63] + B10_current[7] * data_[55];

  data_[72] = D00_[0];
  data_[73] = D00_[1];
  data_[74] = D00_[2];
  data_[75] = D00_[3];
  data_[76] = D00_[4];
  data_[77] = D00_[5];
  data_[78] = D00_[6];
  data_[79] = D00_[7];

  double cB00_current[8];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];
  cB00_current[4] = B00_[4];
  cB00_current[5] = B00_[5];
  cB00_current[6] = B00_[6];
  cB00_current[7] = B00_[7];

  data_[80] = C00_[0] * data_[72] + cB00_current[0];
  data_[81] = C00_[1] * data_[73] + cB00_current[1];
  data_[82] = C00_[2] * data_[74] + cB00_current[2];
  data_[83] = C00_[3] * data_[75] + cB00_current[3];
  data_[84] = C00_[4] * data_[76] + cB00_current[4];
  data_[85] = C00_[5] * data_[77] + cB00_current[5];
  data_[86] = C00_[6] * data_[78] + cB00_current[6];
  data_[87] = C00_[7] * data_[79] + cB00_current[7];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[88] = C00_[0] * data_[80] + B10_current[0] * data_[72] + cB00_current[0] * data_[8];
  data_[89] = C00_[1] * data_[81] + B10_current[1] * data_[73] + cB00_current[1] * data_[9];
  data_[90] = C00_[2] * data_[82] + B10_current[2] * data_[74] + cB00_current[2] * data_[10];
  data_[91] = C00_[3] * data_[83] + B10_current[3] * data_[75] + cB00_current[3] * data_[11];
  data_[92] = C00_[4] * data_[84] + B10_current[4] * data_[76] + cB00_current[4] * data_[12];
  data_[93] = C00_[5] * data_[85] + B10_current[5] * data_[77] + cB00_current[5] * data_[13];
  data_[94] = C00_[6] * data_[86] + B10_current[6] * data_[78] + cB00_current[6] * data_[14];
  data_[95] = C00_[7] * data_[87] + B10_current[7] * data_[79] + cB00_current[7] * data_[15];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[96] = C00_[0] * data_[88] + B10_current[0] * data_[80] + cB00_current[0] * data_[16];
  data_[97] = C00_[1] * data_[89] + B10_current[1] * data_[81] + cB00_current[1] * data_[17];
  data_[98] = C00_[2] * data_[90] + B10_current[2] * data_[82] + cB00_current[2] * data_[18];
  data_[99] = C00_[3] * data_[91] + B10_current[3] * data_[83] + cB00_current[3] * data_[19];
  data_[100] = C00_[4] * data_[92] + B10_current[4] * data_[84] + cB00_current[4] * data_[20];
  data_[101] = C00_[5] * data_[93] + B10_current[5] * data_[85] + cB00_current[5] * data_[21];
  data_[102] = C00_[6] * data_[94] + B10_current[6] * data_[86] + cB00_current[6] * data_[22];
  data_[103] = C00_[7] * data_[95] + B10_current[7] * data_[87] + cB00_current[7] * data_[23];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[104] = C00_[0] * data_[96] + B10_current[0] * data_[88] + cB00_current[0] * data_[24];
  data_[105] = C00_[1] * data_[97] + B10_current[1] * data_[89] + cB00_current[1] * data_[25];
  data_[106] = C00_[2] * data_[98] + B10_current[2] * data_[90] + cB00_current[2] * data_[26];
  data_[107] = C00_[3] * data_[99] + B10_current[3] * data_[91] + cB00_current[3] * data_[27];
  data_[108] = C00_[4] * data_[100] + B10_current[4] * data_[92] + cB00_current[4] * data_[28];
  data_[109] = C00_[5] * data_[101] + B10_current[5] * data_[93] + cB00_current[5] * data_[29];
  data_[110] = C00_[6] * data_[102] + B10_current[6] * data_[94] + cB00_current[6] * data_[30];
  data_[111] = C00_[7] * data_[103] + B10_current[7] * data_[95] + cB00_current[7] * data_[31];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[112] = C00_[0] * data_[104] + B10_current[0] * data_[96] + cB00_current[0] * data_[32];
  data_[113] = C00_[1] * data_[105] + B10_current[1] * data_[97] + cB00_current[1] * data_[33];
  data_[114] = C00_[2] * data_[106] + B10_current[2] * data_[98] + cB00_current[2] * data_[34];
  data_[115] = C00_[3] * data_[107] + B10_current[3] * data_[99] + cB00_current[3] * data_[35];
  data_[116] = C00_[4] * data_[108] + B10_current[4] * data_[100] + cB00_current[4] * data_[36];
  data_[117] = C00_[5] * data_[109] + B10_current[5] * data_[101] + cB00_current[5] * data_[37];
  data_[118] = C00_[6] * data_[110] + B10_current[6] * data_[102] + cB00_current[6] * data_[38];
  data_[119] = C00_[7] * data_[111] + B10_current[7] * data_[103] + cB00_current[7] * data_[39];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[120] = C00_[0] * data_[112] + B10_current[0] * data_[104] + cB00_current[0] * data_[40];
  data_[121] = C00_[1] * data_[113] + B10_current[1] * data_[105] + cB00_current[1] * data_[41];
  data_[122] = C00_[2] * data_[114] + B10_current[2] * data_[106] + cB00_current[2] * data_[42];
  data_[123] = C00_[3] * data_[115] + B10_current[3] * data_[107] + cB00_current[3] * data_[43];
  data_[124] = C00_[4] * data_[116] + B10_current[4] * data_[108] + cB00_current[4] * data_[44];
  data_[125] = C00_[5] * data_[117] + B10_current[5] * data_[109] + cB00_current[5] * data_[45];
  data_[126] = C00_[6] * data_[118] + B10_current[6] * data_[110] + cB00_current[6] * data_[46];
  data_[127] = C00_[7] * data_[119] + B10_current[7] * data_[111] + cB00_current[7] * data_[47];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[128] = C00_[0] * data_[120] + B10_current[0] * data_[112] + cB00_current[0] * data_[48];
  data_[129] = C00_[1] * data_[121] + B10_current[1] * data_[113] + cB00_current[1] * data_[49];
  data_[130] = C00_[2] * data_[122] + B10_current[2] * data_[114] + cB00_current[2] * data_[50];
  data_[131] = C00_[3] * data_[123] + B10_current[3] * data_[115] + cB00_current[3] * data_[51];
  data_[132] = C00_[4] * data_[124] + B10_current[4] * data_[116] + cB00_current[4] * data_[52];
  data_[133] = C00_[5] * data_[125] + B10_current[5] * data_[117] + cB00_current[5] * data_[53];
  data_[134] = C00_[6] * data_[126] + B10_current[6] * data_[118] + cB00_current[6] * data_[54];
  data_[135] = C00_[7] * data_[127] + B10_current[7] * data_[119] + cB00_current[7] * data_[55];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[136] = C00_[0] * data_[128] + B10_current[0] * data_[120] + cB00_current[0] * data_[56];
  data_[137] = C00_[1] * data_[129] + B10_current[1] * data_[121] + cB00_current[1] * data_[57];
  data_[138] = C00_[2] * data_[130] + B10_current[2] * data_[122] + cB00_current[2] * data_[58];
  data_[139] = C00_[3] * data_[131] + B10_current[3] * data_[123] + cB00_current[3] * data_[59];
  data_[140] = C00_[4] * data_[132] + B10_current[4] * data_[124] + cB00_current[4] * data_[60];
  data_[141] = C00_[5] * data_[133] + B10_current[5] * data_[125] + cB00_current[5] * data_[61];
  data_[142] = C00_[6] * data_[134] + B10_current[6] * data_[126] + cB00_current[6] * data_[62];
  data_[143] = C00_[7] * data_[135] + B10_current[7] * data_[127] + cB00_current[7] * data_[63];

  double B01_current[8];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];
  B01_current[4] = B01_[4];
  B01_current[5] = B01_[5];
  B01_current[6] = B01_[6];
  B01_current[7] = B01_[7];

  data_[144] = D00_[0] * data_[72] + B01_current[0];
  data_[145] = D00_[1] * data_[73] + B01_current[1];
  data_[146] = D00_[2] * data_[74] + B01_current[2];
  data_[147] = D00_[3] * data_[75] + B01_current[3];
  data_[148] = D00_[4] * data_[76] + B01_current[4];
  data_[149] = D00_[5] * data_[77] + B01_current[5];
  data_[150] = D00_[6] * data_[78] + B01_current[6];
  data_[151] = D00_[7] * data_[79] + B01_current[7];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[152] = C00_[0] * data_[144] + cB00_current[0] * data_[72];
  data_[153] = C00_[1] * data_[145] + cB00_current[1] * data_[73];
  data_[154] = C00_[2] * data_[146] + cB00_current[2] * data_[74];
  data_[155] = C00_[3] * data_[147] + cB00_current[3] * data_[75];
  data_[156] = C00_[4] * data_[148] + cB00_current[4] * data_[76];
  data_[157] = C00_[5] * data_[149] + cB00_current[5] * data_[77];
  data_[158] = C00_[6] * data_[150] + cB00_current[6] * data_[78];
  data_[159] = C00_[7] * data_[151] + cB00_current[7] * data_[79];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[160] = C00_[0] * data_[152] + B10_current[0] * data_[144] + cB00_current[0] * data_[80];
  data_[161] = C00_[1] * data_[153] + B10_current[1] * data_[145] + cB00_current[1] * data_[81];
  data_[162] = C00_[2] * data_[154] + B10_current[2] * data_[146] + cB00_current[2] * data_[82];
  data_[163] = C00_[3] * data_[155] + B10_current[3] * data_[147] + cB00_current[3] * data_[83];
  data_[164] = C00_[4] * data_[156] + B10_current[4] * data_[148] + cB00_current[4] * data_[84];
  data_[165] = C00_[5] * data_[157] + B10_current[5] * data_[149] + cB00_current[5] * data_[85];
  data_[166] = C00_[6] * data_[158] + B10_current[6] * data_[150] + cB00_current[6] * data_[86];
  data_[167] = C00_[7] * data_[159] + B10_current[7] * data_[151] + cB00_current[7] * data_[87];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[168] = C00_[0] * data_[160] + B10_current[0] * data_[152] + cB00_current[0] * data_[88];
  data_[169] = C00_[1] * data_[161] + B10_current[1] * data_[153] + cB00_current[1] * data_[89];
  data_[170] = C00_[2] * data_[162] + B10_current[2] * data_[154] + cB00_current[2] * data_[90];
  data_[171] = C00_[3] * data_[163] + B10_current[3] * data_[155] + cB00_current[3] * data_[91];
  data_[172] = C00_[4] * data_[164] + B10_current[4] * data_[156] + cB00_current[4] * data_[92];
  data_[173] = C00_[5] * data_[165] + B10_current[5] * data_[157] + cB00_current[5] * data_[93];
  data_[174] = C00_[6] * data_[166] + B10_current[6] * data_[158] + cB00_current[6] * data_[94];
  data_[175] = C00_[7] * data_[167] + B10_current[7] * data_[159] + cB00_current[7] * data_[95];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[176] = C00_[0] * data_[168] + B10_current[0] * data_[160] + cB00_current[0] * data_[96];
  data_[177] = C00_[1] * data_[169] + B10_current[1] * data_[161] + cB00_current[1] * data_[97];
  data_[178] = C00_[2] * data_[170] + B10_current[2] * data_[162] + cB00_current[2] * data_[98];
  data_[179] = C00_[3] * data_[171] + B10_current[3] * data_[163] + cB00_current[3] * data_[99];
  data_[180] = C00_[4] * data_[172] + B10_current[4] * data_[164] + cB00_current[4] * data_[100];
  data_[181] = C00_[5] * data_[173] + B10_current[5] * data_[165] + cB00_current[5] * data_[101];
  data_[182] = C00_[6] * data_[174] + B10_current[6] * data_[166] + cB00_current[6] * data_[102];
  data_[183] = C00_[7] * data_[175] + B10_current[7] * data_[167] + cB00_current[7] * data_[103];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[184] = C00_[0] * data_[176] + B10_current[0] * data_[168] + cB00_current[0] * data_[104];
  data_[185] = C00_[1] * data_[177] + B10_current[1] * data_[169] + cB00_current[1] * data_[105];
  data_[186] = C00_[2] * data_[178] + B10_current[2] * data_[170] + cB00_current[2] * data_[106];
  data_[187] = C00_[3] * data_[179] + B10_current[3] * data_[171] + cB00_current[3] * data_[107];
  data_[188] = C00_[4] * data_[180] + B10_current[4] * data_[172] + cB00_current[4] * data_[108];
  data_[189] = C00_[5] * data_[181] + B10_current[5] * data_[173] + cB00_current[5] * data_[109];
  data_[190] = C00_[6] * data_[182] + B10_current[6] * data_[174] + cB00_current[6] * data_[110];
  data_[191] = C00_[7] * data_[183] + B10_current[7] * data_[175] + cB00_current[7] * data_[111];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[192] = C00_[0] * data_[184] + B10_current[0] * data_[176] + cB00_current[0] * data_[112];
  data_[193] = C00_[1] * data_[185] + B10_current[1] * data_[177] + cB00_current[1] * data_[113];
  data_[194] = C00_[2] * data_[186] + B10_current[2] * data_[178] + cB00_current[2] * data_[114];
  data_[195] = C00_[3] * data_[187] + B10_current[3] * data_[179] + cB00_current[3] * data_[115];
  data_[196] = C00_[4] * data_[188] + B10_current[4] * data_[180] + cB00_current[4] * data_[116];
  data_[197] = C00_[5] * data_[189] + B10_current[5] * data_[181] + cB00_current[5] * data_[117];
  data_[198] = C00_[6] * data_[190] + B10_current[6] * data_[182] + cB00_current[6] * data_[118];
  data_[199] = C00_[7] * data_[191] + B10_current[7] * data_[183] + cB00_current[7] * data_[119];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[200] = C00_[0] * data_[192] + B10_current[0] * data_[184] + cB00_current[0] * data_[120];
  data_[201] = C00_[1] * data_[193] + B10_current[1] * data_[185] + cB00_current[1] * data_[121];
  data_[202] = C00_[2] * data_[194] + B10_current[2] * data_[186] + cB00_current[2] * data_[122];
  data_[203] = C00_[3] * data_[195] + B10_current[3] * data_[187] + cB00_current[3] * data_[123];
  data_[204] = C00_[4] * data_[196] + B10_current[4] * data_[188] + cB00_current[4] * data_[124];
  data_[205] = C00_[5] * data_[197] + B10_current[5] * data_[189] + cB00_current[5] * data_[125];
  data_[206] = C00_[6] * data_[198] + B10_current[6] * data_[190] + cB00_current[6] * data_[126];
  data_[207] = C00_[7] * data_[199] + B10_current[7] * data_[191] + cB00_current[7] * data_[127];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[208] = C00_[0] * data_[200] + B10_current[0] * data_[192] + cB00_current[0] * data_[128];
  data_[209] = C00_[1] * data_[201] + B10_current[1] * data_[193] + cB00_current[1] * data_[129];
  data_[210] = C00_[2] * data_[202] + B10_current[2] * data_[194] + cB00_current[2] * data_[130];
  data_[211] = C00_[3] * data_[203] + B10_current[3] * data_[195] + cB00_current[3] * data_[131];
  data_[212] = C00_[4] * data_[204] + B10_current[4] * data_[196] + cB00_current[4] * data_[132];
  data_[213] = C00_[5] * data_[205] + B10_current[5] * data_[197] + cB00_current[5] * data_[133];
  data_[214] = C00_[6] * data_[206] + B10_current[6] * data_[198] + cB00_current[6] * data_[134];
  data_[215] = C00_[7] * data_[207] + B10_current[7] * data_[199] + cB00_current[7] * data_[135];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[216] = D00_[0] * data_[144] + B01_current[0] * data_[72];
  data_[217] = D00_[1] * data_[145] + B01_current[1] * data_[73];
  data_[218] = D00_[2] * data_[146] + B01_current[2] * data_[74];
  data_[219] = D00_[3] * data_[147] + B01_current[3] * data_[75];
  data_[220] = D00_[4] * data_[148] + B01_current[4] * data_[76];
  data_[221] = D00_[5] * data_[149] + B01_current[5] * data_[77];
  data_[222] = D00_[6] * data_[150] + B01_current[6] * data_[78];
  data_[223] = D00_[7] * data_[151] + B01_current[7] * data_[79];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[224] = C00_[0] * data_[216] + cB00_current[0] * data_[144];
  data_[225] = C00_[1] * data_[217] + cB00_current[1] * data_[145];
  data_[226] = C00_[2] * data_[218] + cB00_current[2] * data_[146];
  data_[227] = C00_[3] * data_[219] + cB00_current[3] * data_[147];
  data_[228] = C00_[4] * data_[220] + cB00_current[4] * data_[148];
  data_[229] = C00_[5] * data_[221] + cB00_current[5] * data_[149];
  data_[230] = C00_[6] * data_[222] + cB00_current[6] * data_[150];
  data_[231] = C00_[7] * data_[223] + cB00_current[7] * data_[151];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[232] = C00_[0] * data_[224] + B10_current[0] * data_[216] + cB00_current[0] * data_[152];
  data_[233] = C00_[1] * data_[225] + B10_current[1] * data_[217] + cB00_current[1] * data_[153];
  data_[234] = C00_[2] * data_[226] + B10_current[2] * data_[218] + cB00_current[2] * data_[154];
  data_[235] = C00_[3] * data_[227] + B10_current[3] * data_[219] + cB00_current[3] * data_[155];
  data_[236] = C00_[4] * data_[228] + B10_current[4] * data_[220] + cB00_current[4] * data_[156];
  data_[237] = C00_[5] * data_[229] + B10_current[5] * data_[221] + cB00_current[5] * data_[157];
  data_[238] = C00_[6] * data_[230] + B10_current[6] * data_[222] + cB00_current[6] * data_[158];
  data_[239] = C00_[7] * data_[231] + B10_current[7] * data_[223] + cB00_current[7] * data_[159];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[240] = C00_[0] * data_[232] + B10_current[0] * data_[224] + cB00_current[0] * data_[160];
  data_[241] = C00_[1] * data_[233] + B10_current[1] * data_[225] + cB00_current[1] * data_[161];
  data_[242] = C00_[2] * data_[234] + B10_current[2] * data_[226] + cB00_current[2] * data_[162];
  data_[243] = C00_[3] * data_[235] + B10_current[3] * data_[227] + cB00_current[3] * data_[163];
  data_[244] = C00_[4] * data_[236] + B10_current[4] * data_[228] + cB00_current[4] * data_[164];
  data_[245] = C00_[5] * data_[237] + B10_current[5] * data_[229] + cB00_current[5] * data_[165];
  data_[246] = C00_[6] * data_[238] + B10_current[6] * data_[230] + cB00_current[6] * data_[166];
  data_[247] = C00_[7] * data_[239] + B10_current[7] * data_[231] + cB00_current[7] * data_[167];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[248] = C00_[0] * data_[240] + B10_current[0] * data_[232] + cB00_current[0] * data_[168];
  data_[249] = C00_[1] * data_[241] + B10_current[1] * data_[233] + cB00_current[1] * data_[169];
  data_[250] = C00_[2] * data_[242] + B10_current[2] * data_[234] + cB00_current[2] * data_[170];
  data_[251] = C00_[3] * data_[243] + B10_current[3] * data_[235] + cB00_current[3] * data_[171];
  data_[252] = C00_[4] * data_[244] + B10_current[4] * data_[236] + cB00_current[4] * data_[172];
  data_[253] = C00_[5] * data_[245] + B10_current[5] * data_[237] + cB00_current[5] * data_[173];
  data_[254] = C00_[6] * data_[246] + B10_current[6] * data_[238] + cB00_current[6] * data_[174];
  data_[255] = C00_[7] * data_[247] + B10_current[7] * data_[239] + cB00_current[7] * data_[175];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[256] = C00_[0] * data_[248] + B10_current[0] * data_[240] + cB00_current[0] * data_[176];
  data_[257] = C00_[1] * data_[249] + B10_current[1] * data_[241] + cB00_current[1] * data_[177];
  data_[258] = C00_[2] * data_[250] + B10_current[2] * data_[242] + cB00_current[2] * data_[178];
  data_[259] = C00_[3] * data_[251] + B10_current[3] * data_[243] + cB00_current[3] * data_[179];
  data_[260] = C00_[4] * data_[252] + B10_current[4] * data_[244] + cB00_current[4] * data_[180];
  data_[261] = C00_[5] * data_[253] + B10_current[5] * data_[245] + cB00_current[5] * data_[181];
  data_[262] = C00_[6] * data_[254] + B10_current[6] * data_[246] + cB00_current[6] * data_[182];
  data_[263] = C00_[7] * data_[255] + B10_current[7] * data_[247] + cB00_current[7] * data_[183];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[264] = C00_[0] * data_[256] + B10_current[0] * data_[248] + cB00_current[0] * data_[184];
  data_[265] = C00_[1] * data_[257] + B10_current[1] * data_[249] + cB00_current[1] * data_[185];
  data_[266] = C00_[2] * data_[258] + B10_current[2] * data_[250] + cB00_current[2] * data_[186];
  data_[267] = C00_[3] * data_[259] + B10_current[3] * data_[251] + cB00_current[3] * data_[187];
  data_[268] = C00_[4] * data_[260] + B10_current[4] * data_[252] + cB00_current[4] * data_[188];
  data_[269] = C00_[5] * data_[261] + B10_current[5] * data_[253] + cB00_current[5] * data_[189];
  data_[270] = C00_[6] * data_[262] + B10_current[6] * data_[254] + cB00_current[6] * data_[190];
  data_[271] = C00_[7] * data_[263] + B10_current[7] * data_[255] + cB00_current[7] * data_[191];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[272] = C00_[0] * data_[264] + B10_current[0] * data_[256] + cB00_current[0] * data_[192];
  data_[273] = C00_[1] * data_[265] + B10_current[1] * data_[257] + cB00_current[1] * data_[193];
  data_[274] = C00_[2] * data_[266] + B10_current[2] * data_[258] + cB00_current[2] * data_[194];
  data_[275] = C00_[3] * data_[267] + B10_current[3] * data_[259] + cB00_current[3] * data_[195];
  data_[276] = C00_[4] * data_[268] + B10_current[4] * data_[260] + cB00_current[4] * data_[196];
  data_[277] = C00_[5] * data_[269] + B10_current[5] * data_[261] + cB00_current[5] * data_[197];
  data_[278] = C00_[6] * data_[270] + B10_current[6] * data_[262] + cB00_current[6] * data_[198];
  data_[279] = C00_[7] * data_[271] + B10_current[7] * data_[263] + cB00_current[7] * data_[199];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[280] = C00_[0] * data_[272] + B10_current[0] * data_[264] + cB00_current[0] * data_[200];
  data_[281] = C00_[1] * data_[273] + B10_current[1] * data_[265] + cB00_current[1] * data_[201];
  data_[282] = C00_[2] * data_[274] + B10_current[2] * data_[266] + cB00_current[2] * data_[202];
  data_[283] = C00_[3] * data_[275] + B10_current[3] * data_[267] + cB00_current[3] * data_[203];
  data_[284] = C00_[4] * data_[276] + B10_current[4] * data_[268] + cB00_current[4] * data_[204];
  data_[285] = C00_[5] * data_[277] + B10_current[5] * data_[269] + cB00_current[5] * data_[205];
  data_[286] = C00_[6] * data_[278] + B10_current[6] * data_[270] + cB00_current[6] * data_[206];
  data_[287] = C00_[7] * data_[279] + B10_current[7] * data_[271] + cB00_current[7] * data_[207];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[288] = D00_[0] * data_[216] + B01_current[0] * data_[144];
  data_[289] = D00_[1] * data_[217] + B01_current[1] * data_[145];
  data_[290] = D00_[2] * data_[218] + B01_current[2] * data_[146];
  data_[291] = D00_[3] * data_[219] + B01_current[3] * data_[147];
  data_[292] = D00_[4] * data_[220] + B01_current[4] * data_[148];
  data_[293] = D00_[5] * data_[221] + B01_current[5] * data_[149];
  data_[294] = D00_[6] * data_[222] + B01_current[6] * data_[150];
  data_[295] = D00_[7] * data_[223] + B01_current[7] * data_[151];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[296] = C00_[0] * data_[288] + cB00_current[0] * data_[216];
  data_[297] = C00_[1] * data_[289] + cB00_current[1] * data_[217];
  data_[298] = C00_[2] * data_[290] + cB00_current[2] * data_[218];
  data_[299] = C00_[3] * data_[291] + cB00_current[3] * data_[219];
  data_[300] = C00_[4] * data_[292] + cB00_current[4] * data_[220];
  data_[301] = C00_[5] * data_[293] + cB00_current[5] * data_[221];
  data_[302] = C00_[6] * data_[294] + cB00_current[6] * data_[222];
  data_[303] = C00_[7] * data_[295] + cB00_current[7] * data_[223];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[304] = C00_[0] * data_[296] + B10_current[0] * data_[288] + cB00_current[0] * data_[224];
  data_[305] = C00_[1] * data_[297] + B10_current[1] * data_[289] + cB00_current[1] * data_[225];
  data_[306] = C00_[2] * data_[298] + B10_current[2] * data_[290] + cB00_current[2] * data_[226];
  data_[307] = C00_[3] * data_[299] + B10_current[3] * data_[291] + cB00_current[3] * data_[227];
  data_[308] = C00_[4] * data_[300] + B10_current[4] * data_[292] + cB00_current[4] * data_[228];
  data_[309] = C00_[5] * data_[301] + B10_current[5] * data_[293] + cB00_current[5] * data_[229];
  data_[310] = C00_[6] * data_[302] + B10_current[6] * data_[294] + cB00_current[6] * data_[230];
  data_[311] = C00_[7] * data_[303] + B10_current[7] * data_[295] + cB00_current[7] * data_[231];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[312] = C00_[0] * data_[304] + B10_current[0] * data_[296] + cB00_current[0] * data_[232];
  data_[313] = C00_[1] * data_[305] + B10_current[1] * data_[297] + cB00_current[1] * data_[233];
  data_[314] = C00_[2] * data_[306] + B10_current[2] * data_[298] + cB00_current[2] * data_[234];
  data_[315] = C00_[3] * data_[307] + B10_current[3] * data_[299] + cB00_current[3] * data_[235];
  data_[316] = C00_[4] * data_[308] + B10_current[4] * data_[300] + cB00_current[4] * data_[236];
  data_[317] = C00_[5] * data_[309] + B10_current[5] * data_[301] + cB00_current[5] * data_[237];
  data_[318] = C00_[6] * data_[310] + B10_current[6] * data_[302] + cB00_current[6] * data_[238];
  data_[319] = C00_[7] * data_[311] + B10_current[7] * data_[303] + cB00_current[7] * data_[239];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[320] = C00_[0] * data_[312] + B10_current[0] * data_[304] + cB00_current[0] * data_[240];
  data_[321] = C00_[1] * data_[313] + B10_current[1] * data_[305] + cB00_current[1] * data_[241];
  data_[322] = C00_[2] * data_[314] + B10_current[2] * data_[306] + cB00_current[2] * data_[242];
  data_[323] = C00_[3] * data_[315] + B10_current[3] * data_[307] + cB00_current[3] * data_[243];
  data_[324] = C00_[4] * data_[316] + B10_current[4] * data_[308] + cB00_current[4] * data_[244];
  data_[325] = C00_[5] * data_[317] + B10_current[5] * data_[309] + cB00_current[5] * data_[245];
  data_[326] = C00_[6] * data_[318] + B10_current[6] * data_[310] + cB00_current[6] * data_[246];
  data_[327] = C00_[7] * data_[319] + B10_current[7] * data_[311] + cB00_current[7] * data_[247];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[328] = C00_[0] * data_[320] + B10_current[0] * data_[312] + cB00_current[0] * data_[248];
  data_[329] = C00_[1] * data_[321] + B10_current[1] * data_[313] + cB00_current[1] * data_[249];
  data_[330] = C00_[2] * data_[322] + B10_current[2] * data_[314] + cB00_current[2] * data_[250];
  data_[331] = C00_[3] * data_[323] + B10_current[3] * data_[315] + cB00_current[3] * data_[251];
  data_[332] = C00_[4] * data_[324] + B10_current[4] * data_[316] + cB00_current[4] * data_[252];
  data_[333] = C00_[5] * data_[325] + B10_current[5] * data_[317] + cB00_current[5] * data_[253];
  data_[334] = C00_[6] * data_[326] + B10_current[6] * data_[318] + cB00_current[6] * data_[254];
  data_[335] = C00_[7] * data_[327] + B10_current[7] * data_[319] + cB00_current[7] * data_[255];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[336] = C00_[0] * data_[328] + B10_current[0] * data_[320] + cB00_current[0] * data_[256];
  data_[337] = C00_[1] * data_[329] + B10_current[1] * data_[321] + cB00_current[1] * data_[257];
  data_[338] = C00_[2] * data_[330] + B10_current[2] * data_[322] + cB00_current[2] * data_[258];
  data_[339] = C00_[3] * data_[331] + B10_current[3] * data_[323] + cB00_current[3] * data_[259];
  data_[340] = C00_[4] * data_[332] + B10_current[4] * data_[324] + cB00_current[4] * data_[260];
  data_[341] = C00_[5] * data_[333] + B10_current[5] * data_[325] + cB00_current[5] * data_[261];
  data_[342] = C00_[6] * data_[334] + B10_current[6] * data_[326] + cB00_current[6] * data_[262];
  data_[343] = C00_[7] * data_[335] + B10_current[7] * data_[327] + cB00_current[7] * data_[263];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[344] = C00_[0] * data_[336] + B10_current[0] * data_[328] + cB00_current[0] * data_[264];
  data_[345] = C00_[1] * data_[337] + B10_current[1] * data_[329] + cB00_current[1] * data_[265];
  data_[346] = C00_[2] * data_[338] + B10_current[2] * data_[330] + cB00_current[2] * data_[266];
  data_[347] = C00_[3] * data_[339] + B10_current[3] * data_[331] + cB00_current[3] * data_[267];
  data_[348] = C00_[4] * data_[340] + B10_current[4] * data_[332] + cB00_current[4] * data_[268];
  data_[349] = C00_[5] * data_[341] + B10_current[5] * data_[333] + cB00_current[5] * data_[269];
  data_[350] = C00_[6] * data_[342] + B10_current[6] * data_[334] + cB00_current[6] * data_[270];
  data_[351] = C00_[7] * data_[343] + B10_current[7] * data_[335] + cB00_current[7] * data_[271];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[352] = C00_[0] * data_[344] + B10_current[0] * data_[336] + cB00_current[0] * data_[272];
  data_[353] = C00_[1] * data_[345] + B10_current[1] * data_[337] + cB00_current[1] * data_[273];
  data_[354] = C00_[2] * data_[346] + B10_current[2] * data_[338] + cB00_current[2] * data_[274];
  data_[355] = C00_[3] * data_[347] + B10_current[3] * data_[339] + cB00_current[3] * data_[275];
  data_[356] = C00_[4] * data_[348] + B10_current[4] * data_[340] + cB00_current[4] * data_[276];
  data_[357] = C00_[5] * data_[349] + B10_current[5] * data_[341] + cB00_current[5] * data_[277];
  data_[358] = C00_[6] * data_[350] + B10_current[6] * data_[342] + cB00_current[6] * data_[278];
  data_[359] = C00_[7] * data_[351] + B10_current[7] * data_[343] + cB00_current[7] * data_[279];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[360] = D00_[0] * data_[288] + B01_current[0] * data_[216];
  data_[361] = D00_[1] * data_[289] + B01_current[1] * data_[217];
  data_[362] = D00_[2] * data_[290] + B01_current[2] * data_[218];
  data_[363] = D00_[3] * data_[291] + B01_current[3] * data_[219];
  data_[364] = D00_[4] * data_[292] + B01_current[4] * data_[220];
  data_[365] = D00_[5] * data_[293] + B01_current[5] * data_[221];
  data_[366] = D00_[6] * data_[294] + B01_current[6] * data_[222];
  data_[367] = D00_[7] * data_[295] + B01_current[7] * data_[223];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[368] = C00_[0] * data_[360] + cB00_current[0] * data_[288];
  data_[369] = C00_[1] * data_[361] + cB00_current[1] * data_[289];
  data_[370] = C00_[2] * data_[362] + cB00_current[2] * data_[290];
  data_[371] = C00_[3] * data_[363] + cB00_current[3] * data_[291];
  data_[372] = C00_[4] * data_[364] + cB00_current[4] * data_[292];
  data_[373] = C00_[5] * data_[365] + cB00_current[5] * data_[293];
  data_[374] = C00_[6] * data_[366] + cB00_current[6] * data_[294];
  data_[375] = C00_[7] * data_[367] + cB00_current[7] * data_[295];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[376] = C00_[0] * data_[368] + B10_current[0] * data_[360] + cB00_current[0] * data_[296];
  data_[377] = C00_[1] * data_[369] + B10_current[1] * data_[361] + cB00_current[1] * data_[297];
  data_[378] = C00_[2] * data_[370] + B10_current[2] * data_[362] + cB00_current[2] * data_[298];
  data_[379] = C00_[3] * data_[371] + B10_current[3] * data_[363] + cB00_current[3] * data_[299];
  data_[380] = C00_[4] * data_[372] + B10_current[4] * data_[364] + cB00_current[4] * data_[300];
  data_[381] = C00_[5] * data_[373] + B10_current[5] * data_[365] + cB00_current[5] * data_[301];
  data_[382] = C00_[6] * data_[374] + B10_current[6] * data_[366] + cB00_current[6] * data_[302];
  data_[383] = C00_[7] * data_[375] + B10_current[7] * data_[367] + cB00_current[7] * data_[303];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[384] = C00_[0] * data_[376] + B10_current[0] * data_[368] + cB00_current[0] * data_[304];
  data_[385] = C00_[1] * data_[377] + B10_current[1] * data_[369] + cB00_current[1] * data_[305];
  data_[386] = C00_[2] * data_[378] + B10_current[2] * data_[370] + cB00_current[2] * data_[306];
  data_[387] = C00_[3] * data_[379] + B10_current[3] * data_[371] + cB00_current[3] * data_[307];
  data_[388] = C00_[4] * data_[380] + B10_current[4] * data_[372] + cB00_current[4] * data_[308];
  data_[389] = C00_[5] * data_[381] + B10_current[5] * data_[373] + cB00_current[5] * data_[309];
  data_[390] = C00_[6] * data_[382] + B10_current[6] * data_[374] + cB00_current[6] * data_[310];
  data_[391] = C00_[7] * data_[383] + B10_current[7] * data_[375] + cB00_current[7] * data_[311];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[392] = C00_[0] * data_[384] + B10_current[0] * data_[376] + cB00_current[0] * data_[312];
  data_[393] = C00_[1] * data_[385] + B10_current[1] * data_[377] + cB00_current[1] * data_[313];
  data_[394] = C00_[2] * data_[386] + B10_current[2] * data_[378] + cB00_current[2] * data_[314];
  data_[395] = C00_[3] * data_[387] + B10_current[3] * data_[379] + cB00_current[3] * data_[315];
  data_[396] = C00_[4] * data_[388] + B10_current[4] * data_[380] + cB00_current[4] * data_[316];
  data_[397] = C00_[5] * data_[389] + B10_current[5] * data_[381] + cB00_current[5] * data_[317];
  data_[398] = C00_[6] * data_[390] + B10_current[6] * data_[382] + cB00_current[6] * data_[318];
  data_[399] = C00_[7] * data_[391] + B10_current[7] * data_[383] + cB00_current[7] * data_[319];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[400] = C00_[0] * data_[392] + B10_current[0] * data_[384] + cB00_current[0] * data_[320];
  data_[401] = C00_[1] * data_[393] + B10_current[1] * data_[385] + cB00_current[1] * data_[321];
  data_[402] = C00_[2] * data_[394] + B10_current[2] * data_[386] + cB00_current[2] * data_[322];
  data_[403] = C00_[3] * data_[395] + B10_current[3] * data_[387] + cB00_current[3] * data_[323];
  data_[404] = C00_[4] * data_[396] + B10_current[4] * data_[388] + cB00_current[4] * data_[324];
  data_[405] = C00_[5] * data_[397] + B10_current[5] * data_[389] + cB00_current[5] * data_[325];
  data_[406] = C00_[6] * data_[398] + B10_current[6] * data_[390] + cB00_current[6] * data_[326];
  data_[407] = C00_[7] * data_[399] + B10_current[7] * data_[391] + cB00_current[7] * data_[327];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[408] = C00_[0] * data_[400] + B10_current[0] * data_[392] + cB00_current[0] * data_[328];
  data_[409] = C00_[1] * data_[401] + B10_current[1] * data_[393] + cB00_current[1] * data_[329];
  data_[410] = C00_[2] * data_[402] + B10_current[2] * data_[394] + cB00_current[2] * data_[330];
  data_[411] = C00_[3] * data_[403] + B10_current[3] * data_[395] + cB00_current[3] * data_[331];
  data_[412] = C00_[4] * data_[404] + B10_current[4] * data_[396] + cB00_current[4] * data_[332];
  data_[413] = C00_[5] * data_[405] + B10_current[5] * data_[397] + cB00_current[5] * data_[333];
  data_[414] = C00_[6] * data_[406] + B10_current[6] * data_[398] + cB00_current[6] * data_[334];
  data_[415] = C00_[7] * data_[407] + B10_current[7] * data_[399] + cB00_current[7] * data_[335];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[416] = C00_[0] * data_[408] + B10_current[0] * data_[400] + cB00_current[0] * data_[336];
  data_[417] = C00_[1] * data_[409] + B10_current[1] * data_[401] + cB00_current[1] * data_[337];
  data_[418] = C00_[2] * data_[410] + B10_current[2] * data_[402] + cB00_current[2] * data_[338];
  data_[419] = C00_[3] * data_[411] + B10_current[3] * data_[403] + cB00_current[3] * data_[339];
  data_[420] = C00_[4] * data_[412] + B10_current[4] * data_[404] + cB00_current[4] * data_[340];
  data_[421] = C00_[5] * data_[413] + B10_current[5] * data_[405] + cB00_current[5] * data_[341];
  data_[422] = C00_[6] * data_[414] + B10_current[6] * data_[406] + cB00_current[6] * data_[342];
  data_[423] = C00_[7] * data_[415] + B10_current[7] * data_[407] + cB00_current[7] * data_[343];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[424] = C00_[0] * data_[416] + B10_current[0] * data_[408] + cB00_current[0] * data_[344];
  data_[425] = C00_[1] * data_[417] + B10_current[1] * data_[409] + cB00_current[1] * data_[345];
  data_[426] = C00_[2] * data_[418] + B10_current[2] * data_[410] + cB00_current[2] * data_[346];
  data_[427] = C00_[3] * data_[419] + B10_current[3] * data_[411] + cB00_current[3] * data_[347];
  data_[428] = C00_[4] * data_[420] + B10_current[4] * data_[412] + cB00_current[4] * data_[348];
  data_[429] = C00_[5] * data_[421] + B10_current[5] * data_[413] + cB00_current[5] * data_[349];
  data_[430] = C00_[6] * data_[422] + B10_current[6] * data_[414] + cB00_current[6] * data_[350];
  data_[431] = C00_[7] * data_[423] + B10_current[7] * data_[415] + cB00_current[7] * data_[351];
}

