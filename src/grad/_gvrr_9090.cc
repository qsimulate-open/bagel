//
// Newint - Parallel electron correlation program.
// Filename: _gvrr_9090.cc
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


#include "gvrrlist.h"

// returns double array of length 900
void GVRRList::_gvrr_9090(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;
  data_[3] = 1.0;
  data_[4] = 1.0;
  data_[5] = 1.0;
  data_[6] = 1.0;
  data_[7] = 1.0;
  data_[8] = 1.0;

  data_[9] = C00_[0];
  data_[10] = C00_[1];
  data_[11] = C00_[2];
  data_[12] = C00_[3];
  data_[13] = C00_[4];
  data_[14] = C00_[5];
  data_[15] = C00_[6];
  data_[16] = C00_[7];
  data_[17] = C00_[8];

  double B10_current[9];
  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];
  B10_current[8] = B10_[8];

  data_[18] = C00_[0] * data_[9] + B10_current[0];
  data_[19] = C00_[1] * data_[10] + B10_current[1];
  data_[20] = C00_[2] * data_[11] + B10_current[2];
  data_[21] = C00_[3] * data_[12] + B10_current[3];
  data_[22] = C00_[4] * data_[13] + B10_current[4];
  data_[23] = C00_[5] * data_[14] + B10_current[5];
  data_[24] = C00_[6] * data_[15] + B10_current[6];
  data_[25] = C00_[7] * data_[16] + B10_current[7];
  data_[26] = C00_[8] * data_[17] + B10_current[8];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[27] = C00_[0] * data_[18] + B10_current[0] * data_[9];
  data_[28] = C00_[1] * data_[19] + B10_current[1] * data_[10];
  data_[29] = C00_[2] * data_[20] + B10_current[2] * data_[11];
  data_[30] = C00_[3] * data_[21] + B10_current[3] * data_[12];
  data_[31] = C00_[4] * data_[22] + B10_current[4] * data_[13];
  data_[32] = C00_[5] * data_[23] + B10_current[5] * data_[14];
  data_[33] = C00_[6] * data_[24] + B10_current[6] * data_[15];
  data_[34] = C00_[7] * data_[25] + B10_current[7] * data_[16];
  data_[35] = C00_[8] * data_[26] + B10_current[8] * data_[17];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[36] = C00_[0] * data_[27] + B10_current[0] * data_[18];
  data_[37] = C00_[1] * data_[28] + B10_current[1] * data_[19];
  data_[38] = C00_[2] * data_[29] + B10_current[2] * data_[20];
  data_[39] = C00_[3] * data_[30] + B10_current[3] * data_[21];
  data_[40] = C00_[4] * data_[31] + B10_current[4] * data_[22];
  data_[41] = C00_[5] * data_[32] + B10_current[5] * data_[23];
  data_[42] = C00_[6] * data_[33] + B10_current[6] * data_[24];
  data_[43] = C00_[7] * data_[34] + B10_current[7] * data_[25];
  data_[44] = C00_[8] * data_[35] + B10_current[8] * data_[26];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[45] = C00_[0] * data_[36] + B10_current[0] * data_[27];
  data_[46] = C00_[1] * data_[37] + B10_current[1] * data_[28];
  data_[47] = C00_[2] * data_[38] + B10_current[2] * data_[29];
  data_[48] = C00_[3] * data_[39] + B10_current[3] * data_[30];
  data_[49] = C00_[4] * data_[40] + B10_current[4] * data_[31];
  data_[50] = C00_[5] * data_[41] + B10_current[5] * data_[32];
  data_[51] = C00_[6] * data_[42] + B10_current[6] * data_[33];
  data_[52] = C00_[7] * data_[43] + B10_current[7] * data_[34];
  data_[53] = C00_[8] * data_[44] + B10_current[8] * data_[35];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[54] = C00_[0] * data_[45] + B10_current[0] * data_[36];
  data_[55] = C00_[1] * data_[46] + B10_current[1] * data_[37];
  data_[56] = C00_[2] * data_[47] + B10_current[2] * data_[38];
  data_[57] = C00_[3] * data_[48] + B10_current[3] * data_[39];
  data_[58] = C00_[4] * data_[49] + B10_current[4] * data_[40];
  data_[59] = C00_[5] * data_[50] + B10_current[5] * data_[41];
  data_[60] = C00_[6] * data_[51] + B10_current[6] * data_[42];
  data_[61] = C00_[7] * data_[52] + B10_current[7] * data_[43];
  data_[62] = C00_[8] * data_[53] + B10_current[8] * data_[44];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[63] = C00_[0] * data_[54] + B10_current[0] * data_[45];
  data_[64] = C00_[1] * data_[55] + B10_current[1] * data_[46];
  data_[65] = C00_[2] * data_[56] + B10_current[2] * data_[47];
  data_[66] = C00_[3] * data_[57] + B10_current[3] * data_[48];
  data_[67] = C00_[4] * data_[58] + B10_current[4] * data_[49];
  data_[68] = C00_[5] * data_[59] + B10_current[5] * data_[50];
  data_[69] = C00_[6] * data_[60] + B10_current[6] * data_[51];
  data_[70] = C00_[7] * data_[61] + B10_current[7] * data_[52];
  data_[71] = C00_[8] * data_[62] + B10_current[8] * data_[53];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[72] = C00_[0] * data_[63] + B10_current[0] * data_[54];
  data_[73] = C00_[1] * data_[64] + B10_current[1] * data_[55];
  data_[74] = C00_[2] * data_[65] + B10_current[2] * data_[56];
  data_[75] = C00_[3] * data_[66] + B10_current[3] * data_[57];
  data_[76] = C00_[4] * data_[67] + B10_current[4] * data_[58];
  data_[77] = C00_[5] * data_[68] + B10_current[5] * data_[59];
  data_[78] = C00_[6] * data_[69] + B10_current[6] * data_[60];
  data_[79] = C00_[7] * data_[70] + B10_current[7] * data_[61];
  data_[80] = C00_[8] * data_[71] + B10_current[8] * data_[62];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[81] = C00_[0] * data_[72] + B10_current[0] * data_[63];
  data_[82] = C00_[1] * data_[73] + B10_current[1] * data_[64];
  data_[83] = C00_[2] * data_[74] + B10_current[2] * data_[65];
  data_[84] = C00_[3] * data_[75] + B10_current[3] * data_[66];
  data_[85] = C00_[4] * data_[76] + B10_current[4] * data_[67];
  data_[86] = C00_[5] * data_[77] + B10_current[5] * data_[68];
  data_[87] = C00_[6] * data_[78] + B10_current[6] * data_[69];
  data_[88] = C00_[7] * data_[79] + B10_current[7] * data_[70];
  data_[89] = C00_[8] * data_[80] + B10_current[8] * data_[71];

  data_[90] = D00_[0];
  data_[91] = D00_[1];
  data_[92] = D00_[2];
  data_[93] = D00_[3];
  data_[94] = D00_[4];
  data_[95] = D00_[5];
  data_[96] = D00_[6];
  data_[97] = D00_[7];
  data_[98] = D00_[8];

  double cB00_current[9];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];
  cB00_current[4] = B00_[4];
  cB00_current[5] = B00_[5];
  cB00_current[6] = B00_[6];
  cB00_current[7] = B00_[7];
  cB00_current[8] = B00_[8];

  data_[99] = C00_[0] * data_[90] + cB00_current[0];
  data_[100] = C00_[1] * data_[91] + cB00_current[1];
  data_[101] = C00_[2] * data_[92] + cB00_current[2];
  data_[102] = C00_[3] * data_[93] + cB00_current[3];
  data_[103] = C00_[4] * data_[94] + cB00_current[4];
  data_[104] = C00_[5] * data_[95] + cB00_current[5];
  data_[105] = C00_[6] * data_[96] + cB00_current[6];
  data_[106] = C00_[7] * data_[97] + cB00_current[7];
  data_[107] = C00_[8] * data_[98] + cB00_current[8];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];
  B10_current[8] = B10_[8];

  data_[108] = C00_[0] * data_[99] + B10_current[0] * data_[90] + cB00_current[0] * data_[9];
  data_[109] = C00_[1] * data_[100] + B10_current[1] * data_[91] + cB00_current[1] * data_[10];
  data_[110] = C00_[2] * data_[101] + B10_current[2] * data_[92] + cB00_current[2] * data_[11];
  data_[111] = C00_[3] * data_[102] + B10_current[3] * data_[93] + cB00_current[3] * data_[12];
  data_[112] = C00_[4] * data_[103] + B10_current[4] * data_[94] + cB00_current[4] * data_[13];
  data_[113] = C00_[5] * data_[104] + B10_current[5] * data_[95] + cB00_current[5] * data_[14];
  data_[114] = C00_[6] * data_[105] + B10_current[6] * data_[96] + cB00_current[6] * data_[15];
  data_[115] = C00_[7] * data_[106] + B10_current[7] * data_[97] + cB00_current[7] * data_[16];
  data_[116] = C00_[8] * data_[107] + B10_current[8] * data_[98] + cB00_current[8] * data_[17];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[117] = C00_[0] * data_[108] + B10_current[0] * data_[99] + cB00_current[0] * data_[18];
  data_[118] = C00_[1] * data_[109] + B10_current[1] * data_[100] + cB00_current[1] * data_[19];
  data_[119] = C00_[2] * data_[110] + B10_current[2] * data_[101] + cB00_current[2] * data_[20];
  data_[120] = C00_[3] * data_[111] + B10_current[3] * data_[102] + cB00_current[3] * data_[21];
  data_[121] = C00_[4] * data_[112] + B10_current[4] * data_[103] + cB00_current[4] * data_[22];
  data_[122] = C00_[5] * data_[113] + B10_current[5] * data_[104] + cB00_current[5] * data_[23];
  data_[123] = C00_[6] * data_[114] + B10_current[6] * data_[105] + cB00_current[6] * data_[24];
  data_[124] = C00_[7] * data_[115] + B10_current[7] * data_[106] + cB00_current[7] * data_[25];
  data_[125] = C00_[8] * data_[116] + B10_current[8] * data_[107] + cB00_current[8] * data_[26];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[126] = C00_[0] * data_[117] + B10_current[0] * data_[108] + cB00_current[0] * data_[27];
  data_[127] = C00_[1] * data_[118] + B10_current[1] * data_[109] + cB00_current[1] * data_[28];
  data_[128] = C00_[2] * data_[119] + B10_current[2] * data_[110] + cB00_current[2] * data_[29];
  data_[129] = C00_[3] * data_[120] + B10_current[3] * data_[111] + cB00_current[3] * data_[30];
  data_[130] = C00_[4] * data_[121] + B10_current[4] * data_[112] + cB00_current[4] * data_[31];
  data_[131] = C00_[5] * data_[122] + B10_current[5] * data_[113] + cB00_current[5] * data_[32];
  data_[132] = C00_[6] * data_[123] + B10_current[6] * data_[114] + cB00_current[6] * data_[33];
  data_[133] = C00_[7] * data_[124] + B10_current[7] * data_[115] + cB00_current[7] * data_[34];
  data_[134] = C00_[8] * data_[125] + B10_current[8] * data_[116] + cB00_current[8] * data_[35];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[135] = C00_[0] * data_[126] + B10_current[0] * data_[117] + cB00_current[0] * data_[36];
  data_[136] = C00_[1] * data_[127] + B10_current[1] * data_[118] + cB00_current[1] * data_[37];
  data_[137] = C00_[2] * data_[128] + B10_current[2] * data_[119] + cB00_current[2] * data_[38];
  data_[138] = C00_[3] * data_[129] + B10_current[3] * data_[120] + cB00_current[3] * data_[39];
  data_[139] = C00_[4] * data_[130] + B10_current[4] * data_[121] + cB00_current[4] * data_[40];
  data_[140] = C00_[5] * data_[131] + B10_current[5] * data_[122] + cB00_current[5] * data_[41];
  data_[141] = C00_[6] * data_[132] + B10_current[6] * data_[123] + cB00_current[6] * data_[42];
  data_[142] = C00_[7] * data_[133] + B10_current[7] * data_[124] + cB00_current[7] * data_[43];
  data_[143] = C00_[8] * data_[134] + B10_current[8] * data_[125] + cB00_current[8] * data_[44];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[144] = C00_[0] * data_[135] + B10_current[0] * data_[126] + cB00_current[0] * data_[45];
  data_[145] = C00_[1] * data_[136] + B10_current[1] * data_[127] + cB00_current[1] * data_[46];
  data_[146] = C00_[2] * data_[137] + B10_current[2] * data_[128] + cB00_current[2] * data_[47];
  data_[147] = C00_[3] * data_[138] + B10_current[3] * data_[129] + cB00_current[3] * data_[48];
  data_[148] = C00_[4] * data_[139] + B10_current[4] * data_[130] + cB00_current[4] * data_[49];
  data_[149] = C00_[5] * data_[140] + B10_current[5] * data_[131] + cB00_current[5] * data_[50];
  data_[150] = C00_[6] * data_[141] + B10_current[6] * data_[132] + cB00_current[6] * data_[51];
  data_[151] = C00_[7] * data_[142] + B10_current[7] * data_[133] + cB00_current[7] * data_[52];
  data_[152] = C00_[8] * data_[143] + B10_current[8] * data_[134] + cB00_current[8] * data_[53];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[153] = C00_[0] * data_[144] + B10_current[0] * data_[135] + cB00_current[0] * data_[54];
  data_[154] = C00_[1] * data_[145] + B10_current[1] * data_[136] + cB00_current[1] * data_[55];
  data_[155] = C00_[2] * data_[146] + B10_current[2] * data_[137] + cB00_current[2] * data_[56];
  data_[156] = C00_[3] * data_[147] + B10_current[3] * data_[138] + cB00_current[3] * data_[57];
  data_[157] = C00_[4] * data_[148] + B10_current[4] * data_[139] + cB00_current[4] * data_[58];
  data_[158] = C00_[5] * data_[149] + B10_current[5] * data_[140] + cB00_current[5] * data_[59];
  data_[159] = C00_[6] * data_[150] + B10_current[6] * data_[141] + cB00_current[6] * data_[60];
  data_[160] = C00_[7] * data_[151] + B10_current[7] * data_[142] + cB00_current[7] * data_[61];
  data_[161] = C00_[8] * data_[152] + B10_current[8] * data_[143] + cB00_current[8] * data_[62];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[162] = C00_[0] * data_[153] + B10_current[0] * data_[144] + cB00_current[0] * data_[63];
  data_[163] = C00_[1] * data_[154] + B10_current[1] * data_[145] + cB00_current[1] * data_[64];
  data_[164] = C00_[2] * data_[155] + B10_current[2] * data_[146] + cB00_current[2] * data_[65];
  data_[165] = C00_[3] * data_[156] + B10_current[3] * data_[147] + cB00_current[3] * data_[66];
  data_[166] = C00_[4] * data_[157] + B10_current[4] * data_[148] + cB00_current[4] * data_[67];
  data_[167] = C00_[5] * data_[158] + B10_current[5] * data_[149] + cB00_current[5] * data_[68];
  data_[168] = C00_[6] * data_[159] + B10_current[6] * data_[150] + cB00_current[6] * data_[69];
  data_[169] = C00_[7] * data_[160] + B10_current[7] * data_[151] + cB00_current[7] * data_[70];
  data_[170] = C00_[8] * data_[161] + B10_current[8] * data_[152] + cB00_current[8] * data_[71];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[171] = C00_[0] * data_[162] + B10_current[0] * data_[153] + cB00_current[0] * data_[72];
  data_[172] = C00_[1] * data_[163] + B10_current[1] * data_[154] + cB00_current[1] * data_[73];
  data_[173] = C00_[2] * data_[164] + B10_current[2] * data_[155] + cB00_current[2] * data_[74];
  data_[174] = C00_[3] * data_[165] + B10_current[3] * data_[156] + cB00_current[3] * data_[75];
  data_[175] = C00_[4] * data_[166] + B10_current[4] * data_[157] + cB00_current[4] * data_[76];
  data_[176] = C00_[5] * data_[167] + B10_current[5] * data_[158] + cB00_current[5] * data_[77];
  data_[177] = C00_[6] * data_[168] + B10_current[6] * data_[159] + cB00_current[6] * data_[78];
  data_[178] = C00_[7] * data_[169] + B10_current[7] * data_[160] + cB00_current[7] * data_[79];
  data_[179] = C00_[8] * data_[170] + B10_current[8] * data_[161] + cB00_current[8] * data_[80];

  double B01_current[9];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];
  B01_current[4] = B01_[4];
  B01_current[5] = B01_[5];
  B01_current[6] = B01_[6];
  B01_current[7] = B01_[7];
  B01_current[8] = B01_[8];

  data_[180] = D00_[0] * data_[90] + B01_current[0];
  data_[181] = D00_[1] * data_[91] + B01_current[1];
  data_[182] = D00_[2] * data_[92] + B01_current[2];
  data_[183] = D00_[3] * data_[93] + B01_current[3];
  data_[184] = D00_[4] * data_[94] + B01_current[4];
  data_[185] = D00_[5] * data_[95] + B01_current[5];
  data_[186] = D00_[6] * data_[96] + B01_current[6];
  data_[187] = D00_[7] * data_[97] + B01_current[7];
  data_[188] = D00_[8] * data_[98] + B01_current[8];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];
  cB00_current[8] += B00_[8];

  data_[189] = C00_[0] * data_[180] + cB00_current[0] * data_[90];
  data_[190] = C00_[1] * data_[181] + cB00_current[1] * data_[91];
  data_[191] = C00_[2] * data_[182] + cB00_current[2] * data_[92];
  data_[192] = C00_[3] * data_[183] + cB00_current[3] * data_[93];
  data_[193] = C00_[4] * data_[184] + cB00_current[4] * data_[94];
  data_[194] = C00_[5] * data_[185] + cB00_current[5] * data_[95];
  data_[195] = C00_[6] * data_[186] + cB00_current[6] * data_[96];
  data_[196] = C00_[7] * data_[187] + cB00_current[7] * data_[97];
  data_[197] = C00_[8] * data_[188] + cB00_current[8] * data_[98];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];
  B10_current[8] = B10_[8];

  data_[198] = C00_[0] * data_[189] + B10_current[0] * data_[180] + cB00_current[0] * data_[99];
  data_[199] = C00_[1] * data_[190] + B10_current[1] * data_[181] + cB00_current[1] * data_[100];
  data_[200] = C00_[2] * data_[191] + B10_current[2] * data_[182] + cB00_current[2] * data_[101];
  data_[201] = C00_[3] * data_[192] + B10_current[3] * data_[183] + cB00_current[3] * data_[102];
  data_[202] = C00_[4] * data_[193] + B10_current[4] * data_[184] + cB00_current[4] * data_[103];
  data_[203] = C00_[5] * data_[194] + B10_current[5] * data_[185] + cB00_current[5] * data_[104];
  data_[204] = C00_[6] * data_[195] + B10_current[6] * data_[186] + cB00_current[6] * data_[105];
  data_[205] = C00_[7] * data_[196] + B10_current[7] * data_[187] + cB00_current[7] * data_[106];
  data_[206] = C00_[8] * data_[197] + B10_current[8] * data_[188] + cB00_current[8] * data_[107];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[207] = C00_[0] * data_[198] + B10_current[0] * data_[189] + cB00_current[0] * data_[108];
  data_[208] = C00_[1] * data_[199] + B10_current[1] * data_[190] + cB00_current[1] * data_[109];
  data_[209] = C00_[2] * data_[200] + B10_current[2] * data_[191] + cB00_current[2] * data_[110];
  data_[210] = C00_[3] * data_[201] + B10_current[3] * data_[192] + cB00_current[3] * data_[111];
  data_[211] = C00_[4] * data_[202] + B10_current[4] * data_[193] + cB00_current[4] * data_[112];
  data_[212] = C00_[5] * data_[203] + B10_current[5] * data_[194] + cB00_current[5] * data_[113];
  data_[213] = C00_[6] * data_[204] + B10_current[6] * data_[195] + cB00_current[6] * data_[114];
  data_[214] = C00_[7] * data_[205] + B10_current[7] * data_[196] + cB00_current[7] * data_[115];
  data_[215] = C00_[8] * data_[206] + B10_current[8] * data_[197] + cB00_current[8] * data_[116];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[216] = C00_[0] * data_[207] + B10_current[0] * data_[198] + cB00_current[0] * data_[117];
  data_[217] = C00_[1] * data_[208] + B10_current[1] * data_[199] + cB00_current[1] * data_[118];
  data_[218] = C00_[2] * data_[209] + B10_current[2] * data_[200] + cB00_current[2] * data_[119];
  data_[219] = C00_[3] * data_[210] + B10_current[3] * data_[201] + cB00_current[3] * data_[120];
  data_[220] = C00_[4] * data_[211] + B10_current[4] * data_[202] + cB00_current[4] * data_[121];
  data_[221] = C00_[5] * data_[212] + B10_current[5] * data_[203] + cB00_current[5] * data_[122];
  data_[222] = C00_[6] * data_[213] + B10_current[6] * data_[204] + cB00_current[6] * data_[123];
  data_[223] = C00_[7] * data_[214] + B10_current[7] * data_[205] + cB00_current[7] * data_[124];
  data_[224] = C00_[8] * data_[215] + B10_current[8] * data_[206] + cB00_current[8] * data_[125];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[225] = C00_[0] * data_[216] + B10_current[0] * data_[207] + cB00_current[0] * data_[126];
  data_[226] = C00_[1] * data_[217] + B10_current[1] * data_[208] + cB00_current[1] * data_[127];
  data_[227] = C00_[2] * data_[218] + B10_current[2] * data_[209] + cB00_current[2] * data_[128];
  data_[228] = C00_[3] * data_[219] + B10_current[3] * data_[210] + cB00_current[3] * data_[129];
  data_[229] = C00_[4] * data_[220] + B10_current[4] * data_[211] + cB00_current[4] * data_[130];
  data_[230] = C00_[5] * data_[221] + B10_current[5] * data_[212] + cB00_current[5] * data_[131];
  data_[231] = C00_[6] * data_[222] + B10_current[6] * data_[213] + cB00_current[6] * data_[132];
  data_[232] = C00_[7] * data_[223] + B10_current[7] * data_[214] + cB00_current[7] * data_[133];
  data_[233] = C00_[8] * data_[224] + B10_current[8] * data_[215] + cB00_current[8] * data_[134];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[234] = C00_[0] * data_[225] + B10_current[0] * data_[216] + cB00_current[0] * data_[135];
  data_[235] = C00_[1] * data_[226] + B10_current[1] * data_[217] + cB00_current[1] * data_[136];
  data_[236] = C00_[2] * data_[227] + B10_current[2] * data_[218] + cB00_current[2] * data_[137];
  data_[237] = C00_[3] * data_[228] + B10_current[3] * data_[219] + cB00_current[3] * data_[138];
  data_[238] = C00_[4] * data_[229] + B10_current[4] * data_[220] + cB00_current[4] * data_[139];
  data_[239] = C00_[5] * data_[230] + B10_current[5] * data_[221] + cB00_current[5] * data_[140];
  data_[240] = C00_[6] * data_[231] + B10_current[6] * data_[222] + cB00_current[6] * data_[141];
  data_[241] = C00_[7] * data_[232] + B10_current[7] * data_[223] + cB00_current[7] * data_[142];
  data_[242] = C00_[8] * data_[233] + B10_current[8] * data_[224] + cB00_current[8] * data_[143];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[243] = C00_[0] * data_[234] + B10_current[0] * data_[225] + cB00_current[0] * data_[144];
  data_[244] = C00_[1] * data_[235] + B10_current[1] * data_[226] + cB00_current[1] * data_[145];
  data_[245] = C00_[2] * data_[236] + B10_current[2] * data_[227] + cB00_current[2] * data_[146];
  data_[246] = C00_[3] * data_[237] + B10_current[3] * data_[228] + cB00_current[3] * data_[147];
  data_[247] = C00_[4] * data_[238] + B10_current[4] * data_[229] + cB00_current[4] * data_[148];
  data_[248] = C00_[5] * data_[239] + B10_current[5] * data_[230] + cB00_current[5] * data_[149];
  data_[249] = C00_[6] * data_[240] + B10_current[6] * data_[231] + cB00_current[6] * data_[150];
  data_[250] = C00_[7] * data_[241] + B10_current[7] * data_[232] + cB00_current[7] * data_[151];
  data_[251] = C00_[8] * data_[242] + B10_current[8] * data_[233] + cB00_current[8] * data_[152];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[252] = C00_[0] * data_[243] + B10_current[0] * data_[234] + cB00_current[0] * data_[153];
  data_[253] = C00_[1] * data_[244] + B10_current[1] * data_[235] + cB00_current[1] * data_[154];
  data_[254] = C00_[2] * data_[245] + B10_current[2] * data_[236] + cB00_current[2] * data_[155];
  data_[255] = C00_[3] * data_[246] + B10_current[3] * data_[237] + cB00_current[3] * data_[156];
  data_[256] = C00_[4] * data_[247] + B10_current[4] * data_[238] + cB00_current[4] * data_[157];
  data_[257] = C00_[5] * data_[248] + B10_current[5] * data_[239] + cB00_current[5] * data_[158];
  data_[258] = C00_[6] * data_[249] + B10_current[6] * data_[240] + cB00_current[6] * data_[159];
  data_[259] = C00_[7] * data_[250] + B10_current[7] * data_[241] + cB00_current[7] * data_[160];
  data_[260] = C00_[8] * data_[251] + B10_current[8] * data_[242] + cB00_current[8] * data_[161];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[261] = C00_[0] * data_[252] + B10_current[0] * data_[243] + cB00_current[0] * data_[162];
  data_[262] = C00_[1] * data_[253] + B10_current[1] * data_[244] + cB00_current[1] * data_[163];
  data_[263] = C00_[2] * data_[254] + B10_current[2] * data_[245] + cB00_current[2] * data_[164];
  data_[264] = C00_[3] * data_[255] + B10_current[3] * data_[246] + cB00_current[3] * data_[165];
  data_[265] = C00_[4] * data_[256] + B10_current[4] * data_[247] + cB00_current[4] * data_[166];
  data_[266] = C00_[5] * data_[257] + B10_current[5] * data_[248] + cB00_current[5] * data_[167];
  data_[267] = C00_[6] * data_[258] + B10_current[6] * data_[249] + cB00_current[6] * data_[168];
  data_[268] = C00_[7] * data_[259] + B10_current[7] * data_[250] + cB00_current[7] * data_[169];
  data_[269] = C00_[8] * data_[260] + B10_current[8] * data_[251] + cB00_current[8] * data_[170];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];
  B01_current[8] += B01_[8];

  data_[270] = D00_[0] * data_[180] + B01_current[0] * data_[90];
  data_[271] = D00_[1] * data_[181] + B01_current[1] * data_[91];
  data_[272] = D00_[2] * data_[182] + B01_current[2] * data_[92];
  data_[273] = D00_[3] * data_[183] + B01_current[3] * data_[93];
  data_[274] = D00_[4] * data_[184] + B01_current[4] * data_[94];
  data_[275] = D00_[5] * data_[185] + B01_current[5] * data_[95];
  data_[276] = D00_[6] * data_[186] + B01_current[6] * data_[96];
  data_[277] = D00_[7] * data_[187] + B01_current[7] * data_[97];
  data_[278] = D00_[8] * data_[188] + B01_current[8] * data_[98];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];
  cB00_current[8] += B00_[8];

  data_[279] = C00_[0] * data_[270] + cB00_current[0] * data_[180];
  data_[280] = C00_[1] * data_[271] + cB00_current[1] * data_[181];
  data_[281] = C00_[2] * data_[272] + cB00_current[2] * data_[182];
  data_[282] = C00_[3] * data_[273] + cB00_current[3] * data_[183];
  data_[283] = C00_[4] * data_[274] + cB00_current[4] * data_[184];
  data_[284] = C00_[5] * data_[275] + cB00_current[5] * data_[185];
  data_[285] = C00_[6] * data_[276] + cB00_current[6] * data_[186];
  data_[286] = C00_[7] * data_[277] + cB00_current[7] * data_[187];
  data_[287] = C00_[8] * data_[278] + cB00_current[8] * data_[188];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];
  B10_current[8] = B10_[8];

  data_[288] = C00_[0] * data_[279] + B10_current[0] * data_[270] + cB00_current[0] * data_[189];
  data_[289] = C00_[1] * data_[280] + B10_current[1] * data_[271] + cB00_current[1] * data_[190];
  data_[290] = C00_[2] * data_[281] + B10_current[2] * data_[272] + cB00_current[2] * data_[191];
  data_[291] = C00_[3] * data_[282] + B10_current[3] * data_[273] + cB00_current[3] * data_[192];
  data_[292] = C00_[4] * data_[283] + B10_current[4] * data_[274] + cB00_current[4] * data_[193];
  data_[293] = C00_[5] * data_[284] + B10_current[5] * data_[275] + cB00_current[5] * data_[194];
  data_[294] = C00_[6] * data_[285] + B10_current[6] * data_[276] + cB00_current[6] * data_[195];
  data_[295] = C00_[7] * data_[286] + B10_current[7] * data_[277] + cB00_current[7] * data_[196];
  data_[296] = C00_[8] * data_[287] + B10_current[8] * data_[278] + cB00_current[8] * data_[197];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[297] = C00_[0] * data_[288] + B10_current[0] * data_[279] + cB00_current[0] * data_[198];
  data_[298] = C00_[1] * data_[289] + B10_current[1] * data_[280] + cB00_current[1] * data_[199];
  data_[299] = C00_[2] * data_[290] + B10_current[2] * data_[281] + cB00_current[2] * data_[200];
  data_[300] = C00_[3] * data_[291] + B10_current[3] * data_[282] + cB00_current[3] * data_[201];
  data_[301] = C00_[4] * data_[292] + B10_current[4] * data_[283] + cB00_current[4] * data_[202];
  data_[302] = C00_[5] * data_[293] + B10_current[5] * data_[284] + cB00_current[5] * data_[203];
  data_[303] = C00_[6] * data_[294] + B10_current[6] * data_[285] + cB00_current[6] * data_[204];
  data_[304] = C00_[7] * data_[295] + B10_current[7] * data_[286] + cB00_current[7] * data_[205];
  data_[305] = C00_[8] * data_[296] + B10_current[8] * data_[287] + cB00_current[8] * data_[206];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[306] = C00_[0] * data_[297] + B10_current[0] * data_[288] + cB00_current[0] * data_[207];
  data_[307] = C00_[1] * data_[298] + B10_current[1] * data_[289] + cB00_current[1] * data_[208];
  data_[308] = C00_[2] * data_[299] + B10_current[2] * data_[290] + cB00_current[2] * data_[209];
  data_[309] = C00_[3] * data_[300] + B10_current[3] * data_[291] + cB00_current[3] * data_[210];
  data_[310] = C00_[4] * data_[301] + B10_current[4] * data_[292] + cB00_current[4] * data_[211];
  data_[311] = C00_[5] * data_[302] + B10_current[5] * data_[293] + cB00_current[5] * data_[212];
  data_[312] = C00_[6] * data_[303] + B10_current[6] * data_[294] + cB00_current[6] * data_[213];
  data_[313] = C00_[7] * data_[304] + B10_current[7] * data_[295] + cB00_current[7] * data_[214];
  data_[314] = C00_[8] * data_[305] + B10_current[8] * data_[296] + cB00_current[8] * data_[215];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[315] = C00_[0] * data_[306] + B10_current[0] * data_[297] + cB00_current[0] * data_[216];
  data_[316] = C00_[1] * data_[307] + B10_current[1] * data_[298] + cB00_current[1] * data_[217];
  data_[317] = C00_[2] * data_[308] + B10_current[2] * data_[299] + cB00_current[2] * data_[218];
  data_[318] = C00_[3] * data_[309] + B10_current[3] * data_[300] + cB00_current[3] * data_[219];
  data_[319] = C00_[4] * data_[310] + B10_current[4] * data_[301] + cB00_current[4] * data_[220];
  data_[320] = C00_[5] * data_[311] + B10_current[5] * data_[302] + cB00_current[5] * data_[221];
  data_[321] = C00_[6] * data_[312] + B10_current[6] * data_[303] + cB00_current[6] * data_[222];
  data_[322] = C00_[7] * data_[313] + B10_current[7] * data_[304] + cB00_current[7] * data_[223];
  data_[323] = C00_[8] * data_[314] + B10_current[8] * data_[305] + cB00_current[8] * data_[224];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[324] = C00_[0] * data_[315] + B10_current[0] * data_[306] + cB00_current[0] * data_[225];
  data_[325] = C00_[1] * data_[316] + B10_current[1] * data_[307] + cB00_current[1] * data_[226];
  data_[326] = C00_[2] * data_[317] + B10_current[2] * data_[308] + cB00_current[2] * data_[227];
  data_[327] = C00_[3] * data_[318] + B10_current[3] * data_[309] + cB00_current[3] * data_[228];
  data_[328] = C00_[4] * data_[319] + B10_current[4] * data_[310] + cB00_current[4] * data_[229];
  data_[329] = C00_[5] * data_[320] + B10_current[5] * data_[311] + cB00_current[5] * data_[230];
  data_[330] = C00_[6] * data_[321] + B10_current[6] * data_[312] + cB00_current[6] * data_[231];
  data_[331] = C00_[7] * data_[322] + B10_current[7] * data_[313] + cB00_current[7] * data_[232];
  data_[332] = C00_[8] * data_[323] + B10_current[8] * data_[314] + cB00_current[8] * data_[233];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[333] = C00_[0] * data_[324] + B10_current[0] * data_[315] + cB00_current[0] * data_[234];
  data_[334] = C00_[1] * data_[325] + B10_current[1] * data_[316] + cB00_current[1] * data_[235];
  data_[335] = C00_[2] * data_[326] + B10_current[2] * data_[317] + cB00_current[2] * data_[236];
  data_[336] = C00_[3] * data_[327] + B10_current[3] * data_[318] + cB00_current[3] * data_[237];
  data_[337] = C00_[4] * data_[328] + B10_current[4] * data_[319] + cB00_current[4] * data_[238];
  data_[338] = C00_[5] * data_[329] + B10_current[5] * data_[320] + cB00_current[5] * data_[239];
  data_[339] = C00_[6] * data_[330] + B10_current[6] * data_[321] + cB00_current[6] * data_[240];
  data_[340] = C00_[7] * data_[331] + B10_current[7] * data_[322] + cB00_current[7] * data_[241];
  data_[341] = C00_[8] * data_[332] + B10_current[8] * data_[323] + cB00_current[8] * data_[242];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[342] = C00_[0] * data_[333] + B10_current[0] * data_[324] + cB00_current[0] * data_[243];
  data_[343] = C00_[1] * data_[334] + B10_current[1] * data_[325] + cB00_current[1] * data_[244];
  data_[344] = C00_[2] * data_[335] + B10_current[2] * data_[326] + cB00_current[2] * data_[245];
  data_[345] = C00_[3] * data_[336] + B10_current[3] * data_[327] + cB00_current[3] * data_[246];
  data_[346] = C00_[4] * data_[337] + B10_current[4] * data_[328] + cB00_current[4] * data_[247];
  data_[347] = C00_[5] * data_[338] + B10_current[5] * data_[329] + cB00_current[5] * data_[248];
  data_[348] = C00_[6] * data_[339] + B10_current[6] * data_[330] + cB00_current[6] * data_[249];
  data_[349] = C00_[7] * data_[340] + B10_current[7] * data_[331] + cB00_current[7] * data_[250];
  data_[350] = C00_[8] * data_[341] + B10_current[8] * data_[332] + cB00_current[8] * data_[251];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[351] = C00_[0] * data_[342] + B10_current[0] * data_[333] + cB00_current[0] * data_[252];
  data_[352] = C00_[1] * data_[343] + B10_current[1] * data_[334] + cB00_current[1] * data_[253];
  data_[353] = C00_[2] * data_[344] + B10_current[2] * data_[335] + cB00_current[2] * data_[254];
  data_[354] = C00_[3] * data_[345] + B10_current[3] * data_[336] + cB00_current[3] * data_[255];
  data_[355] = C00_[4] * data_[346] + B10_current[4] * data_[337] + cB00_current[4] * data_[256];
  data_[356] = C00_[5] * data_[347] + B10_current[5] * data_[338] + cB00_current[5] * data_[257];
  data_[357] = C00_[6] * data_[348] + B10_current[6] * data_[339] + cB00_current[6] * data_[258];
  data_[358] = C00_[7] * data_[349] + B10_current[7] * data_[340] + cB00_current[7] * data_[259];
  data_[359] = C00_[8] * data_[350] + B10_current[8] * data_[341] + cB00_current[8] * data_[260];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];
  B01_current[8] += B01_[8];

  data_[360] = D00_[0] * data_[270] + B01_current[0] * data_[180];
  data_[361] = D00_[1] * data_[271] + B01_current[1] * data_[181];
  data_[362] = D00_[2] * data_[272] + B01_current[2] * data_[182];
  data_[363] = D00_[3] * data_[273] + B01_current[3] * data_[183];
  data_[364] = D00_[4] * data_[274] + B01_current[4] * data_[184];
  data_[365] = D00_[5] * data_[275] + B01_current[5] * data_[185];
  data_[366] = D00_[6] * data_[276] + B01_current[6] * data_[186];
  data_[367] = D00_[7] * data_[277] + B01_current[7] * data_[187];
  data_[368] = D00_[8] * data_[278] + B01_current[8] * data_[188];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];
  cB00_current[8] += B00_[8];

  data_[369] = C00_[0] * data_[360] + cB00_current[0] * data_[270];
  data_[370] = C00_[1] * data_[361] + cB00_current[1] * data_[271];
  data_[371] = C00_[2] * data_[362] + cB00_current[2] * data_[272];
  data_[372] = C00_[3] * data_[363] + cB00_current[3] * data_[273];
  data_[373] = C00_[4] * data_[364] + cB00_current[4] * data_[274];
  data_[374] = C00_[5] * data_[365] + cB00_current[5] * data_[275];
  data_[375] = C00_[6] * data_[366] + cB00_current[6] * data_[276];
  data_[376] = C00_[7] * data_[367] + cB00_current[7] * data_[277];
  data_[377] = C00_[8] * data_[368] + cB00_current[8] * data_[278];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];
  B10_current[8] = B10_[8];

  data_[378] = C00_[0] * data_[369] + B10_current[0] * data_[360] + cB00_current[0] * data_[279];
  data_[379] = C00_[1] * data_[370] + B10_current[1] * data_[361] + cB00_current[1] * data_[280];
  data_[380] = C00_[2] * data_[371] + B10_current[2] * data_[362] + cB00_current[2] * data_[281];
  data_[381] = C00_[3] * data_[372] + B10_current[3] * data_[363] + cB00_current[3] * data_[282];
  data_[382] = C00_[4] * data_[373] + B10_current[4] * data_[364] + cB00_current[4] * data_[283];
  data_[383] = C00_[5] * data_[374] + B10_current[5] * data_[365] + cB00_current[5] * data_[284];
  data_[384] = C00_[6] * data_[375] + B10_current[6] * data_[366] + cB00_current[6] * data_[285];
  data_[385] = C00_[7] * data_[376] + B10_current[7] * data_[367] + cB00_current[7] * data_[286];
  data_[386] = C00_[8] * data_[377] + B10_current[8] * data_[368] + cB00_current[8] * data_[287];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[387] = C00_[0] * data_[378] + B10_current[0] * data_[369] + cB00_current[0] * data_[288];
  data_[388] = C00_[1] * data_[379] + B10_current[1] * data_[370] + cB00_current[1] * data_[289];
  data_[389] = C00_[2] * data_[380] + B10_current[2] * data_[371] + cB00_current[2] * data_[290];
  data_[390] = C00_[3] * data_[381] + B10_current[3] * data_[372] + cB00_current[3] * data_[291];
  data_[391] = C00_[4] * data_[382] + B10_current[4] * data_[373] + cB00_current[4] * data_[292];
  data_[392] = C00_[5] * data_[383] + B10_current[5] * data_[374] + cB00_current[5] * data_[293];
  data_[393] = C00_[6] * data_[384] + B10_current[6] * data_[375] + cB00_current[6] * data_[294];
  data_[394] = C00_[7] * data_[385] + B10_current[7] * data_[376] + cB00_current[7] * data_[295];
  data_[395] = C00_[8] * data_[386] + B10_current[8] * data_[377] + cB00_current[8] * data_[296];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[396] = C00_[0] * data_[387] + B10_current[0] * data_[378] + cB00_current[0] * data_[297];
  data_[397] = C00_[1] * data_[388] + B10_current[1] * data_[379] + cB00_current[1] * data_[298];
  data_[398] = C00_[2] * data_[389] + B10_current[2] * data_[380] + cB00_current[2] * data_[299];
  data_[399] = C00_[3] * data_[390] + B10_current[3] * data_[381] + cB00_current[3] * data_[300];
  data_[400] = C00_[4] * data_[391] + B10_current[4] * data_[382] + cB00_current[4] * data_[301];
  data_[401] = C00_[5] * data_[392] + B10_current[5] * data_[383] + cB00_current[5] * data_[302];
  data_[402] = C00_[6] * data_[393] + B10_current[6] * data_[384] + cB00_current[6] * data_[303];
  data_[403] = C00_[7] * data_[394] + B10_current[7] * data_[385] + cB00_current[7] * data_[304];
  data_[404] = C00_[8] * data_[395] + B10_current[8] * data_[386] + cB00_current[8] * data_[305];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[405] = C00_[0] * data_[396] + B10_current[0] * data_[387] + cB00_current[0] * data_[306];
  data_[406] = C00_[1] * data_[397] + B10_current[1] * data_[388] + cB00_current[1] * data_[307];
  data_[407] = C00_[2] * data_[398] + B10_current[2] * data_[389] + cB00_current[2] * data_[308];
  data_[408] = C00_[3] * data_[399] + B10_current[3] * data_[390] + cB00_current[3] * data_[309];
  data_[409] = C00_[4] * data_[400] + B10_current[4] * data_[391] + cB00_current[4] * data_[310];
  data_[410] = C00_[5] * data_[401] + B10_current[5] * data_[392] + cB00_current[5] * data_[311];
  data_[411] = C00_[6] * data_[402] + B10_current[6] * data_[393] + cB00_current[6] * data_[312];
  data_[412] = C00_[7] * data_[403] + B10_current[7] * data_[394] + cB00_current[7] * data_[313];
  data_[413] = C00_[8] * data_[404] + B10_current[8] * data_[395] + cB00_current[8] * data_[314];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[414] = C00_[0] * data_[405] + B10_current[0] * data_[396] + cB00_current[0] * data_[315];
  data_[415] = C00_[1] * data_[406] + B10_current[1] * data_[397] + cB00_current[1] * data_[316];
  data_[416] = C00_[2] * data_[407] + B10_current[2] * data_[398] + cB00_current[2] * data_[317];
  data_[417] = C00_[3] * data_[408] + B10_current[3] * data_[399] + cB00_current[3] * data_[318];
  data_[418] = C00_[4] * data_[409] + B10_current[4] * data_[400] + cB00_current[4] * data_[319];
  data_[419] = C00_[5] * data_[410] + B10_current[5] * data_[401] + cB00_current[5] * data_[320];
  data_[420] = C00_[6] * data_[411] + B10_current[6] * data_[402] + cB00_current[6] * data_[321];
  data_[421] = C00_[7] * data_[412] + B10_current[7] * data_[403] + cB00_current[7] * data_[322];
  data_[422] = C00_[8] * data_[413] + B10_current[8] * data_[404] + cB00_current[8] * data_[323];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[423] = C00_[0] * data_[414] + B10_current[0] * data_[405] + cB00_current[0] * data_[324];
  data_[424] = C00_[1] * data_[415] + B10_current[1] * data_[406] + cB00_current[1] * data_[325];
  data_[425] = C00_[2] * data_[416] + B10_current[2] * data_[407] + cB00_current[2] * data_[326];
  data_[426] = C00_[3] * data_[417] + B10_current[3] * data_[408] + cB00_current[3] * data_[327];
  data_[427] = C00_[4] * data_[418] + B10_current[4] * data_[409] + cB00_current[4] * data_[328];
  data_[428] = C00_[5] * data_[419] + B10_current[5] * data_[410] + cB00_current[5] * data_[329];
  data_[429] = C00_[6] * data_[420] + B10_current[6] * data_[411] + cB00_current[6] * data_[330];
  data_[430] = C00_[7] * data_[421] + B10_current[7] * data_[412] + cB00_current[7] * data_[331];
  data_[431] = C00_[8] * data_[422] + B10_current[8] * data_[413] + cB00_current[8] * data_[332];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[432] = C00_[0] * data_[423] + B10_current[0] * data_[414] + cB00_current[0] * data_[333];
  data_[433] = C00_[1] * data_[424] + B10_current[1] * data_[415] + cB00_current[1] * data_[334];
  data_[434] = C00_[2] * data_[425] + B10_current[2] * data_[416] + cB00_current[2] * data_[335];
  data_[435] = C00_[3] * data_[426] + B10_current[3] * data_[417] + cB00_current[3] * data_[336];
  data_[436] = C00_[4] * data_[427] + B10_current[4] * data_[418] + cB00_current[4] * data_[337];
  data_[437] = C00_[5] * data_[428] + B10_current[5] * data_[419] + cB00_current[5] * data_[338];
  data_[438] = C00_[6] * data_[429] + B10_current[6] * data_[420] + cB00_current[6] * data_[339];
  data_[439] = C00_[7] * data_[430] + B10_current[7] * data_[421] + cB00_current[7] * data_[340];
  data_[440] = C00_[8] * data_[431] + B10_current[8] * data_[422] + cB00_current[8] * data_[341];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[441] = C00_[0] * data_[432] + B10_current[0] * data_[423] + cB00_current[0] * data_[342];
  data_[442] = C00_[1] * data_[433] + B10_current[1] * data_[424] + cB00_current[1] * data_[343];
  data_[443] = C00_[2] * data_[434] + B10_current[2] * data_[425] + cB00_current[2] * data_[344];
  data_[444] = C00_[3] * data_[435] + B10_current[3] * data_[426] + cB00_current[3] * data_[345];
  data_[445] = C00_[4] * data_[436] + B10_current[4] * data_[427] + cB00_current[4] * data_[346];
  data_[446] = C00_[5] * data_[437] + B10_current[5] * data_[428] + cB00_current[5] * data_[347];
  data_[447] = C00_[6] * data_[438] + B10_current[6] * data_[429] + cB00_current[6] * data_[348];
  data_[448] = C00_[7] * data_[439] + B10_current[7] * data_[430] + cB00_current[7] * data_[349];
  data_[449] = C00_[8] * data_[440] + B10_current[8] * data_[431] + cB00_current[8] * data_[350];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];
  B01_current[8] += B01_[8];

  data_[450] = D00_[0] * data_[360] + B01_current[0] * data_[270];
  data_[451] = D00_[1] * data_[361] + B01_current[1] * data_[271];
  data_[452] = D00_[2] * data_[362] + B01_current[2] * data_[272];
  data_[453] = D00_[3] * data_[363] + B01_current[3] * data_[273];
  data_[454] = D00_[4] * data_[364] + B01_current[4] * data_[274];
  data_[455] = D00_[5] * data_[365] + B01_current[5] * data_[275];
  data_[456] = D00_[6] * data_[366] + B01_current[6] * data_[276];
  data_[457] = D00_[7] * data_[367] + B01_current[7] * data_[277];
  data_[458] = D00_[8] * data_[368] + B01_current[8] * data_[278];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];
  cB00_current[8] += B00_[8];

  data_[459] = C00_[0] * data_[450] + cB00_current[0] * data_[360];
  data_[460] = C00_[1] * data_[451] + cB00_current[1] * data_[361];
  data_[461] = C00_[2] * data_[452] + cB00_current[2] * data_[362];
  data_[462] = C00_[3] * data_[453] + cB00_current[3] * data_[363];
  data_[463] = C00_[4] * data_[454] + cB00_current[4] * data_[364];
  data_[464] = C00_[5] * data_[455] + cB00_current[5] * data_[365];
  data_[465] = C00_[6] * data_[456] + cB00_current[6] * data_[366];
  data_[466] = C00_[7] * data_[457] + cB00_current[7] * data_[367];
  data_[467] = C00_[8] * data_[458] + cB00_current[8] * data_[368];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];
  B10_current[8] = B10_[8];

  data_[468] = C00_[0] * data_[459] + B10_current[0] * data_[450] + cB00_current[0] * data_[369];
  data_[469] = C00_[1] * data_[460] + B10_current[1] * data_[451] + cB00_current[1] * data_[370];
  data_[470] = C00_[2] * data_[461] + B10_current[2] * data_[452] + cB00_current[2] * data_[371];
  data_[471] = C00_[3] * data_[462] + B10_current[3] * data_[453] + cB00_current[3] * data_[372];
  data_[472] = C00_[4] * data_[463] + B10_current[4] * data_[454] + cB00_current[4] * data_[373];
  data_[473] = C00_[5] * data_[464] + B10_current[5] * data_[455] + cB00_current[5] * data_[374];
  data_[474] = C00_[6] * data_[465] + B10_current[6] * data_[456] + cB00_current[6] * data_[375];
  data_[475] = C00_[7] * data_[466] + B10_current[7] * data_[457] + cB00_current[7] * data_[376];
  data_[476] = C00_[8] * data_[467] + B10_current[8] * data_[458] + cB00_current[8] * data_[377];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[477] = C00_[0] * data_[468] + B10_current[0] * data_[459] + cB00_current[0] * data_[378];
  data_[478] = C00_[1] * data_[469] + B10_current[1] * data_[460] + cB00_current[1] * data_[379];
  data_[479] = C00_[2] * data_[470] + B10_current[2] * data_[461] + cB00_current[2] * data_[380];
  data_[480] = C00_[3] * data_[471] + B10_current[3] * data_[462] + cB00_current[3] * data_[381];
  data_[481] = C00_[4] * data_[472] + B10_current[4] * data_[463] + cB00_current[4] * data_[382];
  data_[482] = C00_[5] * data_[473] + B10_current[5] * data_[464] + cB00_current[5] * data_[383];
  data_[483] = C00_[6] * data_[474] + B10_current[6] * data_[465] + cB00_current[6] * data_[384];
  data_[484] = C00_[7] * data_[475] + B10_current[7] * data_[466] + cB00_current[7] * data_[385];
  data_[485] = C00_[8] * data_[476] + B10_current[8] * data_[467] + cB00_current[8] * data_[386];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[486] = C00_[0] * data_[477] + B10_current[0] * data_[468] + cB00_current[0] * data_[387];
  data_[487] = C00_[1] * data_[478] + B10_current[1] * data_[469] + cB00_current[1] * data_[388];
  data_[488] = C00_[2] * data_[479] + B10_current[2] * data_[470] + cB00_current[2] * data_[389];
  data_[489] = C00_[3] * data_[480] + B10_current[3] * data_[471] + cB00_current[3] * data_[390];
  data_[490] = C00_[4] * data_[481] + B10_current[4] * data_[472] + cB00_current[4] * data_[391];
  data_[491] = C00_[5] * data_[482] + B10_current[5] * data_[473] + cB00_current[5] * data_[392];
  data_[492] = C00_[6] * data_[483] + B10_current[6] * data_[474] + cB00_current[6] * data_[393];
  data_[493] = C00_[7] * data_[484] + B10_current[7] * data_[475] + cB00_current[7] * data_[394];
  data_[494] = C00_[8] * data_[485] + B10_current[8] * data_[476] + cB00_current[8] * data_[395];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[495] = C00_[0] * data_[486] + B10_current[0] * data_[477] + cB00_current[0] * data_[396];
  data_[496] = C00_[1] * data_[487] + B10_current[1] * data_[478] + cB00_current[1] * data_[397];
  data_[497] = C00_[2] * data_[488] + B10_current[2] * data_[479] + cB00_current[2] * data_[398];
  data_[498] = C00_[3] * data_[489] + B10_current[3] * data_[480] + cB00_current[3] * data_[399];
  data_[499] = C00_[4] * data_[490] + B10_current[4] * data_[481] + cB00_current[4] * data_[400];
  data_[500] = C00_[5] * data_[491] + B10_current[5] * data_[482] + cB00_current[5] * data_[401];
  data_[501] = C00_[6] * data_[492] + B10_current[6] * data_[483] + cB00_current[6] * data_[402];
  data_[502] = C00_[7] * data_[493] + B10_current[7] * data_[484] + cB00_current[7] * data_[403];
  data_[503] = C00_[8] * data_[494] + B10_current[8] * data_[485] + cB00_current[8] * data_[404];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[504] = C00_[0] * data_[495] + B10_current[0] * data_[486] + cB00_current[0] * data_[405];
  data_[505] = C00_[1] * data_[496] + B10_current[1] * data_[487] + cB00_current[1] * data_[406];
  data_[506] = C00_[2] * data_[497] + B10_current[2] * data_[488] + cB00_current[2] * data_[407];
  data_[507] = C00_[3] * data_[498] + B10_current[3] * data_[489] + cB00_current[3] * data_[408];
  data_[508] = C00_[4] * data_[499] + B10_current[4] * data_[490] + cB00_current[4] * data_[409];
  data_[509] = C00_[5] * data_[500] + B10_current[5] * data_[491] + cB00_current[5] * data_[410];
  data_[510] = C00_[6] * data_[501] + B10_current[6] * data_[492] + cB00_current[6] * data_[411];
  data_[511] = C00_[7] * data_[502] + B10_current[7] * data_[493] + cB00_current[7] * data_[412];
  data_[512] = C00_[8] * data_[503] + B10_current[8] * data_[494] + cB00_current[8] * data_[413];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[513] = C00_[0] * data_[504] + B10_current[0] * data_[495] + cB00_current[0] * data_[414];
  data_[514] = C00_[1] * data_[505] + B10_current[1] * data_[496] + cB00_current[1] * data_[415];
  data_[515] = C00_[2] * data_[506] + B10_current[2] * data_[497] + cB00_current[2] * data_[416];
  data_[516] = C00_[3] * data_[507] + B10_current[3] * data_[498] + cB00_current[3] * data_[417];
  data_[517] = C00_[4] * data_[508] + B10_current[4] * data_[499] + cB00_current[4] * data_[418];
  data_[518] = C00_[5] * data_[509] + B10_current[5] * data_[500] + cB00_current[5] * data_[419];
  data_[519] = C00_[6] * data_[510] + B10_current[6] * data_[501] + cB00_current[6] * data_[420];
  data_[520] = C00_[7] * data_[511] + B10_current[7] * data_[502] + cB00_current[7] * data_[421];
  data_[521] = C00_[8] * data_[512] + B10_current[8] * data_[503] + cB00_current[8] * data_[422];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[522] = C00_[0] * data_[513] + B10_current[0] * data_[504] + cB00_current[0] * data_[423];
  data_[523] = C00_[1] * data_[514] + B10_current[1] * data_[505] + cB00_current[1] * data_[424];
  data_[524] = C00_[2] * data_[515] + B10_current[2] * data_[506] + cB00_current[2] * data_[425];
  data_[525] = C00_[3] * data_[516] + B10_current[3] * data_[507] + cB00_current[3] * data_[426];
  data_[526] = C00_[4] * data_[517] + B10_current[4] * data_[508] + cB00_current[4] * data_[427];
  data_[527] = C00_[5] * data_[518] + B10_current[5] * data_[509] + cB00_current[5] * data_[428];
  data_[528] = C00_[6] * data_[519] + B10_current[6] * data_[510] + cB00_current[6] * data_[429];
  data_[529] = C00_[7] * data_[520] + B10_current[7] * data_[511] + cB00_current[7] * data_[430];
  data_[530] = C00_[8] * data_[521] + B10_current[8] * data_[512] + cB00_current[8] * data_[431];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[531] = C00_[0] * data_[522] + B10_current[0] * data_[513] + cB00_current[0] * data_[432];
  data_[532] = C00_[1] * data_[523] + B10_current[1] * data_[514] + cB00_current[1] * data_[433];
  data_[533] = C00_[2] * data_[524] + B10_current[2] * data_[515] + cB00_current[2] * data_[434];
  data_[534] = C00_[3] * data_[525] + B10_current[3] * data_[516] + cB00_current[3] * data_[435];
  data_[535] = C00_[4] * data_[526] + B10_current[4] * data_[517] + cB00_current[4] * data_[436];
  data_[536] = C00_[5] * data_[527] + B10_current[5] * data_[518] + cB00_current[5] * data_[437];
  data_[537] = C00_[6] * data_[528] + B10_current[6] * data_[519] + cB00_current[6] * data_[438];
  data_[538] = C00_[7] * data_[529] + B10_current[7] * data_[520] + cB00_current[7] * data_[439];
  data_[539] = C00_[8] * data_[530] + B10_current[8] * data_[521] + cB00_current[8] * data_[440];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];
  B01_current[8] += B01_[8];

  data_[540] = D00_[0] * data_[450] + B01_current[0] * data_[360];
  data_[541] = D00_[1] * data_[451] + B01_current[1] * data_[361];
  data_[542] = D00_[2] * data_[452] + B01_current[2] * data_[362];
  data_[543] = D00_[3] * data_[453] + B01_current[3] * data_[363];
  data_[544] = D00_[4] * data_[454] + B01_current[4] * data_[364];
  data_[545] = D00_[5] * data_[455] + B01_current[5] * data_[365];
  data_[546] = D00_[6] * data_[456] + B01_current[6] * data_[366];
  data_[547] = D00_[7] * data_[457] + B01_current[7] * data_[367];
  data_[548] = D00_[8] * data_[458] + B01_current[8] * data_[368];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];
  cB00_current[8] += B00_[8];

  data_[549] = C00_[0] * data_[540] + cB00_current[0] * data_[450];
  data_[550] = C00_[1] * data_[541] + cB00_current[1] * data_[451];
  data_[551] = C00_[2] * data_[542] + cB00_current[2] * data_[452];
  data_[552] = C00_[3] * data_[543] + cB00_current[3] * data_[453];
  data_[553] = C00_[4] * data_[544] + cB00_current[4] * data_[454];
  data_[554] = C00_[5] * data_[545] + cB00_current[5] * data_[455];
  data_[555] = C00_[6] * data_[546] + cB00_current[6] * data_[456];
  data_[556] = C00_[7] * data_[547] + cB00_current[7] * data_[457];
  data_[557] = C00_[8] * data_[548] + cB00_current[8] * data_[458];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];
  B10_current[8] = B10_[8];

  data_[558] = C00_[0] * data_[549] + B10_current[0] * data_[540] + cB00_current[0] * data_[459];
  data_[559] = C00_[1] * data_[550] + B10_current[1] * data_[541] + cB00_current[1] * data_[460];
  data_[560] = C00_[2] * data_[551] + B10_current[2] * data_[542] + cB00_current[2] * data_[461];
  data_[561] = C00_[3] * data_[552] + B10_current[3] * data_[543] + cB00_current[3] * data_[462];
  data_[562] = C00_[4] * data_[553] + B10_current[4] * data_[544] + cB00_current[4] * data_[463];
  data_[563] = C00_[5] * data_[554] + B10_current[5] * data_[545] + cB00_current[5] * data_[464];
  data_[564] = C00_[6] * data_[555] + B10_current[6] * data_[546] + cB00_current[6] * data_[465];
  data_[565] = C00_[7] * data_[556] + B10_current[7] * data_[547] + cB00_current[7] * data_[466];
  data_[566] = C00_[8] * data_[557] + B10_current[8] * data_[548] + cB00_current[8] * data_[467];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[567] = C00_[0] * data_[558] + B10_current[0] * data_[549] + cB00_current[0] * data_[468];
  data_[568] = C00_[1] * data_[559] + B10_current[1] * data_[550] + cB00_current[1] * data_[469];
  data_[569] = C00_[2] * data_[560] + B10_current[2] * data_[551] + cB00_current[2] * data_[470];
  data_[570] = C00_[3] * data_[561] + B10_current[3] * data_[552] + cB00_current[3] * data_[471];
  data_[571] = C00_[4] * data_[562] + B10_current[4] * data_[553] + cB00_current[4] * data_[472];
  data_[572] = C00_[5] * data_[563] + B10_current[5] * data_[554] + cB00_current[5] * data_[473];
  data_[573] = C00_[6] * data_[564] + B10_current[6] * data_[555] + cB00_current[6] * data_[474];
  data_[574] = C00_[7] * data_[565] + B10_current[7] * data_[556] + cB00_current[7] * data_[475];
  data_[575] = C00_[8] * data_[566] + B10_current[8] * data_[557] + cB00_current[8] * data_[476];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[576] = C00_[0] * data_[567] + B10_current[0] * data_[558] + cB00_current[0] * data_[477];
  data_[577] = C00_[1] * data_[568] + B10_current[1] * data_[559] + cB00_current[1] * data_[478];
  data_[578] = C00_[2] * data_[569] + B10_current[2] * data_[560] + cB00_current[2] * data_[479];
  data_[579] = C00_[3] * data_[570] + B10_current[3] * data_[561] + cB00_current[3] * data_[480];
  data_[580] = C00_[4] * data_[571] + B10_current[4] * data_[562] + cB00_current[4] * data_[481];
  data_[581] = C00_[5] * data_[572] + B10_current[5] * data_[563] + cB00_current[5] * data_[482];
  data_[582] = C00_[6] * data_[573] + B10_current[6] * data_[564] + cB00_current[6] * data_[483];
  data_[583] = C00_[7] * data_[574] + B10_current[7] * data_[565] + cB00_current[7] * data_[484];
  data_[584] = C00_[8] * data_[575] + B10_current[8] * data_[566] + cB00_current[8] * data_[485];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[585] = C00_[0] * data_[576] + B10_current[0] * data_[567] + cB00_current[0] * data_[486];
  data_[586] = C00_[1] * data_[577] + B10_current[1] * data_[568] + cB00_current[1] * data_[487];
  data_[587] = C00_[2] * data_[578] + B10_current[2] * data_[569] + cB00_current[2] * data_[488];
  data_[588] = C00_[3] * data_[579] + B10_current[3] * data_[570] + cB00_current[3] * data_[489];
  data_[589] = C00_[4] * data_[580] + B10_current[4] * data_[571] + cB00_current[4] * data_[490];
  data_[590] = C00_[5] * data_[581] + B10_current[5] * data_[572] + cB00_current[5] * data_[491];
  data_[591] = C00_[6] * data_[582] + B10_current[6] * data_[573] + cB00_current[6] * data_[492];
  data_[592] = C00_[7] * data_[583] + B10_current[7] * data_[574] + cB00_current[7] * data_[493];
  data_[593] = C00_[8] * data_[584] + B10_current[8] * data_[575] + cB00_current[8] * data_[494];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[594] = C00_[0] * data_[585] + B10_current[0] * data_[576] + cB00_current[0] * data_[495];
  data_[595] = C00_[1] * data_[586] + B10_current[1] * data_[577] + cB00_current[1] * data_[496];
  data_[596] = C00_[2] * data_[587] + B10_current[2] * data_[578] + cB00_current[2] * data_[497];
  data_[597] = C00_[3] * data_[588] + B10_current[3] * data_[579] + cB00_current[3] * data_[498];
  data_[598] = C00_[4] * data_[589] + B10_current[4] * data_[580] + cB00_current[4] * data_[499];
  data_[599] = C00_[5] * data_[590] + B10_current[5] * data_[581] + cB00_current[5] * data_[500];
  data_[600] = C00_[6] * data_[591] + B10_current[6] * data_[582] + cB00_current[6] * data_[501];
  data_[601] = C00_[7] * data_[592] + B10_current[7] * data_[583] + cB00_current[7] * data_[502];
  data_[602] = C00_[8] * data_[593] + B10_current[8] * data_[584] + cB00_current[8] * data_[503];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[603] = C00_[0] * data_[594] + B10_current[0] * data_[585] + cB00_current[0] * data_[504];
  data_[604] = C00_[1] * data_[595] + B10_current[1] * data_[586] + cB00_current[1] * data_[505];
  data_[605] = C00_[2] * data_[596] + B10_current[2] * data_[587] + cB00_current[2] * data_[506];
  data_[606] = C00_[3] * data_[597] + B10_current[3] * data_[588] + cB00_current[3] * data_[507];
  data_[607] = C00_[4] * data_[598] + B10_current[4] * data_[589] + cB00_current[4] * data_[508];
  data_[608] = C00_[5] * data_[599] + B10_current[5] * data_[590] + cB00_current[5] * data_[509];
  data_[609] = C00_[6] * data_[600] + B10_current[6] * data_[591] + cB00_current[6] * data_[510];
  data_[610] = C00_[7] * data_[601] + B10_current[7] * data_[592] + cB00_current[7] * data_[511];
  data_[611] = C00_[8] * data_[602] + B10_current[8] * data_[593] + cB00_current[8] * data_[512];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[612] = C00_[0] * data_[603] + B10_current[0] * data_[594] + cB00_current[0] * data_[513];
  data_[613] = C00_[1] * data_[604] + B10_current[1] * data_[595] + cB00_current[1] * data_[514];
  data_[614] = C00_[2] * data_[605] + B10_current[2] * data_[596] + cB00_current[2] * data_[515];
  data_[615] = C00_[3] * data_[606] + B10_current[3] * data_[597] + cB00_current[3] * data_[516];
  data_[616] = C00_[4] * data_[607] + B10_current[4] * data_[598] + cB00_current[4] * data_[517];
  data_[617] = C00_[5] * data_[608] + B10_current[5] * data_[599] + cB00_current[5] * data_[518];
  data_[618] = C00_[6] * data_[609] + B10_current[6] * data_[600] + cB00_current[6] * data_[519];
  data_[619] = C00_[7] * data_[610] + B10_current[7] * data_[601] + cB00_current[7] * data_[520];
  data_[620] = C00_[8] * data_[611] + B10_current[8] * data_[602] + cB00_current[8] * data_[521];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[621] = C00_[0] * data_[612] + B10_current[0] * data_[603] + cB00_current[0] * data_[522];
  data_[622] = C00_[1] * data_[613] + B10_current[1] * data_[604] + cB00_current[1] * data_[523];
  data_[623] = C00_[2] * data_[614] + B10_current[2] * data_[605] + cB00_current[2] * data_[524];
  data_[624] = C00_[3] * data_[615] + B10_current[3] * data_[606] + cB00_current[3] * data_[525];
  data_[625] = C00_[4] * data_[616] + B10_current[4] * data_[607] + cB00_current[4] * data_[526];
  data_[626] = C00_[5] * data_[617] + B10_current[5] * data_[608] + cB00_current[5] * data_[527];
  data_[627] = C00_[6] * data_[618] + B10_current[6] * data_[609] + cB00_current[6] * data_[528];
  data_[628] = C00_[7] * data_[619] + B10_current[7] * data_[610] + cB00_current[7] * data_[529];
  data_[629] = C00_[8] * data_[620] + B10_current[8] * data_[611] + cB00_current[8] * data_[530];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];
  B01_current[8] += B01_[8];

  data_[630] = D00_[0] * data_[540] + B01_current[0] * data_[450];
  data_[631] = D00_[1] * data_[541] + B01_current[1] * data_[451];
  data_[632] = D00_[2] * data_[542] + B01_current[2] * data_[452];
  data_[633] = D00_[3] * data_[543] + B01_current[3] * data_[453];
  data_[634] = D00_[4] * data_[544] + B01_current[4] * data_[454];
  data_[635] = D00_[5] * data_[545] + B01_current[5] * data_[455];
  data_[636] = D00_[6] * data_[546] + B01_current[6] * data_[456];
  data_[637] = D00_[7] * data_[547] + B01_current[7] * data_[457];
  data_[638] = D00_[8] * data_[548] + B01_current[8] * data_[458];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];
  cB00_current[8] += B00_[8];

  data_[639] = C00_[0] * data_[630] + cB00_current[0] * data_[540];
  data_[640] = C00_[1] * data_[631] + cB00_current[1] * data_[541];
  data_[641] = C00_[2] * data_[632] + cB00_current[2] * data_[542];
  data_[642] = C00_[3] * data_[633] + cB00_current[3] * data_[543];
  data_[643] = C00_[4] * data_[634] + cB00_current[4] * data_[544];
  data_[644] = C00_[5] * data_[635] + cB00_current[5] * data_[545];
  data_[645] = C00_[6] * data_[636] + cB00_current[6] * data_[546];
  data_[646] = C00_[7] * data_[637] + cB00_current[7] * data_[547];
  data_[647] = C00_[8] * data_[638] + cB00_current[8] * data_[548];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];
  B10_current[8] = B10_[8];

  data_[648] = C00_[0] * data_[639] + B10_current[0] * data_[630] + cB00_current[0] * data_[549];
  data_[649] = C00_[1] * data_[640] + B10_current[1] * data_[631] + cB00_current[1] * data_[550];
  data_[650] = C00_[2] * data_[641] + B10_current[2] * data_[632] + cB00_current[2] * data_[551];
  data_[651] = C00_[3] * data_[642] + B10_current[3] * data_[633] + cB00_current[3] * data_[552];
  data_[652] = C00_[4] * data_[643] + B10_current[4] * data_[634] + cB00_current[4] * data_[553];
  data_[653] = C00_[5] * data_[644] + B10_current[5] * data_[635] + cB00_current[5] * data_[554];
  data_[654] = C00_[6] * data_[645] + B10_current[6] * data_[636] + cB00_current[6] * data_[555];
  data_[655] = C00_[7] * data_[646] + B10_current[7] * data_[637] + cB00_current[7] * data_[556];
  data_[656] = C00_[8] * data_[647] + B10_current[8] * data_[638] + cB00_current[8] * data_[557];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[657] = C00_[0] * data_[648] + B10_current[0] * data_[639] + cB00_current[0] * data_[558];
  data_[658] = C00_[1] * data_[649] + B10_current[1] * data_[640] + cB00_current[1] * data_[559];
  data_[659] = C00_[2] * data_[650] + B10_current[2] * data_[641] + cB00_current[2] * data_[560];
  data_[660] = C00_[3] * data_[651] + B10_current[3] * data_[642] + cB00_current[3] * data_[561];
  data_[661] = C00_[4] * data_[652] + B10_current[4] * data_[643] + cB00_current[4] * data_[562];
  data_[662] = C00_[5] * data_[653] + B10_current[5] * data_[644] + cB00_current[5] * data_[563];
  data_[663] = C00_[6] * data_[654] + B10_current[6] * data_[645] + cB00_current[6] * data_[564];
  data_[664] = C00_[7] * data_[655] + B10_current[7] * data_[646] + cB00_current[7] * data_[565];
  data_[665] = C00_[8] * data_[656] + B10_current[8] * data_[647] + cB00_current[8] * data_[566];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[666] = C00_[0] * data_[657] + B10_current[0] * data_[648] + cB00_current[0] * data_[567];
  data_[667] = C00_[1] * data_[658] + B10_current[1] * data_[649] + cB00_current[1] * data_[568];
  data_[668] = C00_[2] * data_[659] + B10_current[2] * data_[650] + cB00_current[2] * data_[569];
  data_[669] = C00_[3] * data_[660] + B10_current[3] * data_[651] + cB00_current[3] * data_[570];
  data_[670] = C00_[4] * data_[661] + B10_current[4] * data_[652] + cB00_current[4] * data_[571];
  data_[671] = C00_[5] * data_[662] + B10_current[5] * data_[653] + cB00_current[5] * data_[572];
  data_[672] = C00_[6] * data_[663] + B10_current[6] * data_[654] + cB00_current[6] * data_[573];
  data_[673] = C00_[7] * data_[664] + B10_current[7] * data_[655] + cB00_current[7] * data_[574];
  data_[674] = C00_[8] * data_[665] + B10_current[8] * data_[656] + cB00_current[8] * data_[575];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[675] = C00_[0] * data_[666] + B10_current[0] * data_[657] + cB00_current[0] * data_[576];
  data_[676] = C00_[1] * data_[667] + B10_current[1] * data_[658] + cB00_current[1] * data_[577];
  data_[677] = C00_[2] * data_[668] + B10_current[2] * data_[659] + cB00_current[2] * data_[578];
  data_[678] = C00_[3] * data_[669] + B10_current[3] * data_[660] + cB00_current[3] * data_[579];
  data_[679] = C00_[4] * data_[670] + B10_current[4] * data_[661] + cB00_current[4] * data_[580];
  data_[680] = C00_[5] * data_[671] + B10_current[5] * data_[662] + cB00_current[5] * data_[581];
  data_[681] = C00_[6] * data_[672] + B10_current[6] * data_[663] + cB00_current[6] * data_[582];
  data_[682] = C00_[7] * data_[673] + B10_current[7] * data_[664] + cB00_current[7] * data_[583];
  data_[683] = C00_[8] * data_[674] + B10_current[8] * data_[665] + cB00_current[8] * data_[584];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[684] = C00_[0] * data_[675] + B10_current[0] * data_[666] + cB00_current[0] * data_[585];
  data_[685] = C00_[1] * data_[676] + B10_current[1] * data_[667] + cB00_current[1] * data_[586];
  data_[686] = C00_[2] * data_[677] + B10_current[2] * data_[668] + cB00_current[2] * data_[587];
  data_[687] = C00_[3] * data_[678] + B10_current[3] * data_[669] + cB00_current[3] * data_[588];
  data_[688] = C00_[4] * data_[679] + B10_current[4] * data_[670] + cB00_current[4] * data_[589];
  data_[689] = C00_[5] * data_[680] + B10_current[5] * data_[671] + cB00_current[5] * data_[590];
  data_[690] = C00_[6] * data_[681] + B10_current[6] * data_[672] + cB00_current[6] * data_[591];
  data_[691] = C00_[7] * data_[682] + B10_current[7] * data_[673] + cB00_current[7] * data_[592];
  data_[692] = C00_[8] * data_[683] + B10_current[8] * data_[674] + cB00_current[8] * data_[593];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[693] = C00_[0] * data_[684] + B10_current[0] * data_[675] + cB00_current[0] * data_[594];
  data_[694] = C00_[1] * data_[685] + B10_current[1] * data_[676] + cB00_current[1] * data_[595];
  data_[695] = C00_[2] * data_[686] + B10_current[2] * data_[677] + cB00_current[2] * data_[596];
  data_[696] = C00_[3] * data_[687] + B10_current[3] * data_[678] + cB00_current[3] * data_[597];
  data_[697] = C00_[4] * data_[688] + B10_current[4] * data_[679] + cB00_current[4] * data_[598];
  data_[698] = C00_[5] * data_[689] + B10_current[5] * data_[680] + cB00_current[5] * data_[599];
  data_[699] = C00_[6] * data_[690] + B10_current[6] * data_[681] + cB00_current[6] * data_[600];
  data_[700] = C00_[7] * data_[691] + B10_current[7] * data_[682] + cB00_current[7] * data_[601];
  data_[701] = C00_[8] * data_[692] + B10_current[8] * data_[683] + cB00_current[8] * data_[602];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[702] = C00_[0] * data_[693] + B10_current[0] * data_[684] + cB00_current[0] * data_[603];
  data_[703] = C00_[1] * data_[694] + B10_current[1] * data_[685] + cB00_current[1] * data_[604];
  data_[704] = C00_[2] * data_[695] + B10_current[2] * data_[686] + cB00_current[2] * data_[605];
  data_[705] = C00_[3] * data_[696] + B10_current[3] * data_[687] + cB00_current[3] * data_[606];
  data_[706] = C00_[4] * data_[697] + B10_current[4] * data_[688] + cB00_current[4] * data_[607];
  data_[707] = C00_[5] * data_[698] + B10_current[5] * data_[689] + cB00_current[5] * data_[608];
  data_[708] = C00_[6] * data_[699] + B10_current[6] * data_[690] + cB00_current[6] * data_[609];
  data_[709] = C00_[7] * data_[700] + B10_current[7] * data_[691] + cB00_current[7] * data_[610];
  data_[710] = C00_[8] * data_[701] + B10_current[8] * data_[692] + cB00_current[8] * data_[611];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[711] = C00_[0] * data_[702] + B10_current[0] * data_[693] + cB00_current[0] * data_[612];
  data_[712] = C00_[1] * data_[703] + B10_current[1] * data_[694] + cB00_current[1] * data_[613];
  data_[713] = C00_[2] * data_[704] + B10_current[2] * data_[695] + cB00_current[2] * data_[614];
  data_[714] = C00_[3] * data_[705] + B10_current[3] * data_[696] + cB00_current[3] * data_[615];
  data_[715] = C00_[4] * data_[706] + B10_current[4] * data_[697] + cB00_current[4] * data_[616];
  data_[716] = C00_[5] * data_[707] + B10_current[5] * data_[698] + cB00_current[5] * data_[617];
  data_[717] = C00_[6] * data_[708] + B10_current[6] * data_[699] + cB00_current[6] * data_[618];
  data_[718] = C00_[7] * data_[709] + B10_current[7] * data_[700] + cB00_current[7] * data_[619];
  data_[719] = C00_[8] * data_[710] + B10_current[8] * data_[701] + cB00_current[8] * data_[620];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];
  B01_current[8] += B01_[8];

  data_[720] = D00_[0] * data_[630] + B01_current[0] * data_[540];
  data_[721] = D00_[1] * data_[631] + B01_current[1] * data_[541];
  data_[722] = D00_[2] * data_[632] + B01_current[2] * data_[542];
  data_[723] = D00_[3] * data_[633] + B01_current[3] * data_[543];
  data_[724] = D00_[4] * data_[634] + B01_current[4] * data_[544];
  data_[725] = D00_[5] * data_[635] + B01_current[5] * data_[545];
  data_[726] = D00_[6] * data_[636] + B01_current[6] * data_[546];
  data_[727] = D00_[7] * data_[637] + B01_current[7] * data_[547];
  data_[728] = D00_[8] * data_[638] + B01_current[8] * data_[548];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];
  cB00_current[8] += B00_[8];

  data_[729] = C00_[0] * data_[720] + cB00_current[0] * data_[630];
  data_[730] = C00_[1] * data_[721] + cB00_current[1] * data_[631];
  data_[731] = C00_[2] * data_[722] + cB00_current[2] * data_[632];
  data_[732] = C00_[3] * data_[723] + cB00_current[3] * data_[633];
  data_[733] = C00_[4] * data_[724] + cB00_current[4] * data_[634];
  data_[734] = C00_[5] * data_[725] + cB00_current[5] * data_[635];
  data_[735] = C00_[6] * data_[726] + cB00_current[6] * data_[636];
  data_[736] = C00_[7] * data_[727] + cB00_current[7] * data_[637];
  data_[737] = C00_[8] * data_[728] + cB00_current[8] * data_[638];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];
  B10_current[8] = B10_[8];

  data_[738] = C00_[0] * data_[729] + B10_current[0] * data_[720] + cB00_current[0] * data_[639];
  data_[739] = C00_[1] * data_[730] + B10_current[1] * data_[721] + cB00_current[1] * data_[640];
  data_[740] = C00_[2] * data_[731] + B10_current[2] * data_[722] + cB00_current[2] * data_[641];
  data_[741] = C00_[3] * data_[732] + B10_current[3] * data_[723] + cB00_current[3] * data_[642];
  data_[742] = C00_[4] * data_[733] + B10_current[4] * data_[724] + cB00_current[4] * data_[643];
  data_[743] = C00_[5] * data_[734] + B10_current[5] * data_[725] + cB00_current[5] * data_[644];
  data_[744] = C00_[6] * data_[735] + B10_current[6] * data_[726] + cB00_current[6] * data_[645];
  data_[745] = C00_[7] * data_[736] + B10_current[7] * data_[727] + cB00_current[7] * data_[646];
  data_[746] = C00_[8] * data_[737] + B10_current[8] * data_[728] + cB00_current[8] * data_[647];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[747] = C00_[0] * data_[738] + B10_current[0] * data_[729] + cB00_current[0] * data_[648];
  data_[748] = C00_[1] * data_[739] + B10_current[1] * data_[730] + cB00_current[1] * data_[649];
  data_[749] = C00_[2] * data_[740] + B10_current[2] * data_[731] + cB00_current[2] * data_[650];
  data_[750] = C00_[3] * data_[741] + B10_current[3] * data_[732] + cB00_current[3] * data_[651];
  data_[751] = C00_[4] * data_[742] + B10_current[4] * data_[733] + cB00_current[4] * data_[652];
  data_[752] = C00_[5] * data_[743] + B10_current[5] * data_[734] + cB00_current[5] * data_[653];
  data_[753] = C00_[6] * data_[744] + B10_current[6] * data_[735] + cB00_current[6] * data_[654];
  data_[754] = C00_[7] * data_[745] + B10_current[7] * data_[736] + cB00_current[7] * data_[655];
  data_[755] = C00_[8] * data_[746] + B10_current[8] * data_[737] + cB00_current[8] * data_[656];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[756] = C00_[0] * data_[747] + B10_current[0] * data_[738] + cB00_current[0] * data_[657];
  data_[757] = C00_[1] * data_[748] + B10_current[1] * data_[739] + cB00_current[1] * data_[658];
  data_[758] = C00_[2] * data_[749] + B10_current[2] * data_[740] + cB00_current[2] * data_[659];
  data_[759] = C00_[3] * data_[750] + B10_current[3] * data_[741] + cB00_current[3] * data_[660];
  data_[760] = C00_[4] * data_[751] + B10_current[4] * data_[742] + cB00_current[4] * data_[661];
  data_[761] = C00_[5] * data_[752] + B10_current[5] * data_[743] + cB00_current[5] * data_[662];
  data_[762] = C00_[6] * data_[753] + B10_current[6] * data_[744] + cB00_current[6] * data_[663];
  data_[763] = C00_[7] * data_[754] + B10_current[7] * data_[745] + cB00_current[7] * data_[664];
  data_[764] = C00_[8] * data_[755] + B10_current[8] * data_[746] + cB00_current[8] * data_[665];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[765] = C00_[0] * data_[756] + B10_current[0] * data_[747] + cB00_current[0] * data_[666];
  data_[766] = C00_[1] * data_[757] + B10_current[1] * data_[748] + cB00_current[1] * data_[667];
  data_[767] = C00_[2] * data_[758] + B10_current[2] * data_[749] + cB00_current[2] * data_[668];
  data_[768] = C00_[3] * data_[759] + B10_current[3] * data_[750] + cB00_current[3] * data_[669];
  data_[769] = C00_[4] * data_[760] + B10_current[4] * data_[751] + cB00_current[4] * data_[670];
  data_[770] = C00_[5] * data_[761] + B10_current[5] * data_[752] + cB00_current[5] * data_[671];
  data_[771] = C00_[6] * data_[762] + B10_current[6] * data_[753] + cB00_current[6] * data_[672];
  data_[772] = C00_[7] * data_[763] + B10_current[7] * data_[754] + cB00_current[7] * data_[673];
  data_[773] = C00_[8] * data_[764] + B10_current[8] * data_[755] + cB00_current[8] * data_[674];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[774] = C00_[0] * data_[765] + B10_current[0] * data_[756] + cB00_current[0] * data_[675];
  data_[775] = C00_[1] * data_[766] + B10_current[1] * data_[757] + cB00_current[1] * data_[676];
  data_[776] = C00_[2] * data_[767] + B10_current[2] * data_[758] + cB00_current[2] * data_[677];
  data_[777] = C00_[3] * data_[768] + B10_current[3] * data_[759] + cB00_current[3] * data_[678];
  data_[778] = C00_[4] * data_[769] + B10_current[4] * data_[760] + cB00_current[4] * data_[679];
  data_[779] = C00_[5] * data_[770] + B10_current[5] * data_[761] + cB00_current[5] * data_[680];
  data_[780] = C00_[6] * data_[771] + B10_current[6] * data_[762] + cB00_current[6] * data_[681];
  data_[781] = C00_[7] * data_[772] + B10_current[7] * data_[763] + cB00_current[7] * data_[682];
  data_[782] = C00_[8] * data_[773] + B10_current[8] * data_[764] + cB00_current[8] * data_[683];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[783] = C00_[0] * data_[774] + B10_current[0] * data_[765] + cB00_current[0] * data_[684];
  data_[784] = C00_[1] * data_[775] + B10_current[1] * data_[766] + cB00_current[1] * data_[685];
  data_[785] = C00_[2] * data_[776] + B10_current[2] * data_[767] + cB00_current[2] * data_[686];
  data_[786] = C00_[3] * data_[777] + B10_current[3] * data_[768] + cB00_current[3] * data_[687];
  data_[787] = C00_[4] * data_[778] + B10_current[4] * data_[769] + cB00_current[4] * data_[688];
  data_[788] = C00_[5] * data_[779] + B10_current[5] * data_[770] + cB00_current[5] * data_[689];
  data_[789] = C00_[6] * data_[780] + B10_current[6] * data_[771] + cB00_current[6] * data_[690];
  data_[790] = C00_[7] * data_[781] + B10_current[7] * data_[772] + cB00_current[7] * data_[691];
  data_[791] = C00_[8] * data_[782] + B10_current[8] * data_[773] + cB00_current[8] * data_[692];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[792] = C00_[0] * data_[783] + B10_current[0] * data_[774] + cB00_current[0] * data_[693];
  data_[793] = C00_[1] * data_[784] + B10_current[1] * data_[775] + cB00_current[1] * data_[694];
  data_[794] = C00_[2] * data_[785] + B10_current[2] * data_[776] + cB00_current[2] * data_[695];
  data_[795] = C00_[3] * data_[786] + B10_current[3] * data_[777] + cB00_current[3] * data_[696];
  data_[796] = C00_[4] * data_[787] + B10_current[4] * data_[778] + cB00_current[4] * data_[697];
  data_[797] = C00_[5] * data_[788] + B10_current[5] * data_[779] + cB00_current[5] * data_[698];
  data_[798] = C00_[6] * data_[789] + B10_current[6] * data_[780] + cB00_current[6] * data_[699];
  data_[799] = C00_[7] * data_[790] + B10_current[7] * data_[781] + cB00_current[7] * data_[700];
  data_[800] = C00_[8] * data_[791] + B10_current[8] * data_[782] + cB00_current[8] * data_[701];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[801] = C00_[0] * data_[792] + B10_current[0] * data_[783] + cB00_current[0] * data_[702];
  data_[802] = C00_[1] * data_[793] + B10_current[1] * data_[784] + cB00_current[1] * data_[703];
  data_[803] = C00_[2] * data_[794] + B10_current[2] * data_[785] + cB00_current[2] * data_[704];
  data_[804] = C00_[3] * data_[795] + B10_current[3] * data_[786] + cB00_current[3] * data_[705];
  data_[805] = C00_[4] * data_[796] + B10_current[4] * data_[787] + cB00_current[4] * data_[706];
  data_[806] = C00_[5] * data_[797] + B10_current[5] * data_[788] + cB00_current[5] * data_[707];
  data_[807] = C00_[6] * data_[798] + B10_current[6] * data_[789] + cB00_current[6] * data_[708];
  data_[808] = C00_[7] * data_[799] + B10_current[7] * data_[790] + cB00_current[7] * data_[709];
  data_[809] = C00_[8] * data_[800] + B10_current[8] * data_[791] + cB00_current[8] * data_[710];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];
  B01_current[8] += B01_[8];

  data_[810] = D00_[0] * data_[720] + B01_current[0] * data_[630];
  data_[811] = D00_[1] * data_[721] + B01_current[1] * data_[631];
  data_[812] = D00_[2] * data_[722] + B01_current[2] * data_[632];
  data_[813] = D00_[3] * data_[723] + B01_current[3] * data_[633];
  data_[814] = D00_[4] * data_[724] + B01_current[4] * data_[634];
  data_[815] = D00_[5] * data_[725] + B01_current[5] * data_[635];
  data_[816] = D00_[6] * data_[726] + B01_current[6] * data_[636];
  data_[817] = D00_[7] * data_[727] + B01_current[7] * data_[637];
  data_[818] = D00_[8] * data_[728] + B01_current[8] * data_[638];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];
  cB00_current[8] += B00_[8];

  data_[819] = C00_[0] * data_[810] + cB00_current[0] * data_[720];
  data_[820] = C00_[1] * data_[811] + cB00_current[1] * data_[721];
  data_[821] = C00_[2] * data_[812] + cB00_current[2] * data_[722];
  data_[822] = C00_[3] * data_[813] + cB00_current[3] * data_[723];
  data_[823] = C00_[4] * data_[814] + cB00_current[4] * data_[724];
  data_[824] = C00_[5] * data_[815] + cB00_current[5] * data_[725];
  data_[825] = C00_[6] * data_[816] + cB00_current[6] * data_[726];
  data_[826] = C00_[7] * data_[817] + cB00_current[7] * data_[727];
  data_[827] = C00_[8] * data_[818] + cB00_current[8] * data_[728];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];
  B10_current[8] = B10_[8];

  data_[828] = C00_[0] * data_[819] + B10_current[0] * data_[810] + cB00_current[0] * data_[729];
  data_[829] = C00_[1] * data_[820] + B10_current[1] * data_[811] + cB00_current[1] * data_[730];
  data_[830] = C00_[2] * data_[821] + B10_current[2] * data_[812] + cB00_current[2] * data_[731];
  data_[831] = C00_[3] * data_[822] + B10_current[3] * data_[813] + cB00_current[3] * data_[732];
  data_[832] = C00_[4] * data_[823] + B10_current[4] * data_[814] + cB00_current[4] * data_[733];
  data_[833] = C00_[5] * data_[824] + B10_current[5] * data_[815] + cB00_current[5] * data_[734];
  data_[834] = C00_[6] * data_[825] + B10_current[6] * data_[816] + cB00_current[6] * data_[735];
  data_[835] = C00_[7] * data_[826] + B10_current[7] * data_[817] + cB00_current[7] * data_[736];
  data_[836] = C00_[8] * data_[827] + B10_current[8] * data_[818] + cB00_current[8] * data_[737];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[837] = C00_[0] * data_[828] + B10_current[0] * data_[819] + cB00_current[0] * data_[738];
  data_[838] = C00_[1] * data_[829] + B10_current[1] * data_[820] + cB00_current[1] * data_[739];
  data_[839] = C00_[2] * data_[830] + B10_current[2] * data_[821] + cB00_current[2] * data_[740];
  data_[840] = C00_[3] * data_[831] + B10_current[3] * data_[822] + cB00_current[3] * data_[741];
  data_[841] = C00_[4] * data_[832] + B10_current[4] * data_[823] + cB00_current[4] * data_[742];
  data_[842] = C00_[5] * data_[833] + B10_current[5] * data_[824] + cB00_current[5] * data_[743];
  data_[843] = C00_[6] * data_[834] + B10_current[6] * data_[825] + cB00_current[6] * data_[744];
  data_[844] = C00_[7] * data_[835] + B10_current[7] * data_[826] + cB00_current[7] * data_[745];
  data_[845] = C00_[8] * data_[836] + B10_current[8] * data_[827] + cB00_current[8] * data_[746];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[846] = C00_[0] * data_[837] + B10_current[0] * data_[828] + cB00_current[0] * data_[747];
  data_[847] = C00_[1] * data_[838] + B10_current[1] * data_[829] + cB00_current[1] * data_[748];
  data_[848] = C00_[2] * data_[839] + B10_current[2] * data_[830] + cB00_current[2] * data_[749];
  data_[849] = C00_[3] * data_[840] + B10_current[3] * data_[831] + cB00_current[3] * data_[750];
  data_[850] = C00_[4] * data_[841] + B10_current[4] * data_[832] + cB00_current[4] * data_[751];
  data_[851] = C00_[5] * data_[842] + B10_current[5] * data_[833] + cB00_current[5] * data_[752];
  data_[852] = C00_[6] * data_[843] + B10_current[6] * data_[834] + cB00_current[6] * data_[753];
  data_[853] = C00_[7] * data_[844] + B10_current[7] * data_[835] + cB00_current[7] * data_[754];
  data_[854] = C00_[8] * data_[845] + B10_current[8] * data_[836] + cB00_current[8] * data_[755];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[855] = C00_[0] * data_[846] + B10_current[0] * data_[837] + cB00_current[0] * data_[756];
  data_[856] = C00_[1] * data_[847] + B10_current[1] * data_[838] + cB00_current[1] * data_[757];
  data_[857] = C00_[2] * data_[848] + B10_current[2] * data_[839] + cB00_current[2] * data_[758];
  data_[858] = C00_[3] * data_[849] + B10_current[3] * data_[840] + cB00_current[3] * data_[759];
  data_[859] = C00_[4] * data_[850] + B10_current[4] * data_[841] + cB00_current[4] * data_[760];
  data_[860] = C00_[5] * data_[851] + B10_current[5] * data_[842] + cB00_current[5] * data_[761];
  data_[861] = C00_[6] * data_[852] + B10_current[6] * data_[843] + cB00_current[6] * data_[762];
  data_[862] = C00_[7] * data_[853] + B10_current[7] * data_[844] + cB00_current[7] * data_[763];
  data_[863] = C00_[8] * data_[854] + B10_current[8] * data_[845] + cB00_current[8] * data_[764];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[864] = C00_[0] * data_[855] + B10_current[0] * data_[846] + cB00_current[0] * data_[765];
  data_[865] = C00_[1] * data_[856] + B10_current[1] * data_[847] + cB00_current[1] * data_[766];
  data_[866] = C00_[2] * data_[857] + B10_current[2] * data_[848] + cB00_current[2] * data_[767];
  data_[867] = C00_[3] * data_[858] + B10_current[3] * data_[849] + cB00_current[3] * data_[768];
  data_[868] = C00_[4] * data_[859] + B10_current[4] * data_[850] + cB00_current[4] * data_[769];
  data_[869] = C00_[5] * data_[860] + B10_current[5] * data_[851] + cB00_current[5] * data_[770];
  data_[870] = C00_[6] * data_[861] + B10_current[6] * data_[852] + cB00_current[6] * data_[771];
  data_[871] = C00_[7] * data_[862] + B10_current[7] * data_[853] + cB00_current[7] * data_[772];
  data_[872] = C00_[8] * data_[863] + B10_current[8] * data_[854] + cB00_current[8] * data_[773];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[873] = C00_[0] * data_[864] + B10_current[0] * data_[855] + cB00_current[0] * data_[774];
  data_[874] = C00_[1] * data_[865] + B10_current[1] * data_[856] + cB00_current[1] * data_[775];
  data_[875] = C00_[2] * data_[866] + B10_current[2] * data_[857] + cB00_current[2] * data_[776];
  data_[876] = C00_[3] * data_[867] + B10_current[3] * data_[858] + cB00_current[3] * data_[777];
  data_[877] = C00_[4] * data_[868] + B10_current[4] * data_[859] + cB00_current[4] * data_[778];
  data_[878] = C00_[5] * data_[869] + B10_current[5] * data_[860] + cB00_current[5] * data_[779];
  data_[879] = C00_[6] * data_[870] + B10_current[6] * data_[861] + cB00_current[6] * data_[780];
  data_[880] = C00_[7] * data_[871] + B10_current[7] * data_[862] + cB00_current[7] * data_[781];
  data_[881] = C00_[8] * data_[872] + B10_current[8] * data_[863] + cB00_current[8] * data_[782];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[882] = C00_[0] * data_[873] + B10_current[0] * data_[864] + cB00_current[0] * data_[783];
  data_[883] = C00_[1] * data_[874] + B10_current[1] * data_[865] + cB00_current[1] * data_[784];
  data_[884] = C00_[2] * data_[875] + B10_current[2] * data_[866] + cB00_current[2] * data_[785];
  data_[885] = C00_[3] * data_[876] + B10_current[3] * data_[867] + cB00_current[3] * data_[786];
  data_[886] = C00_[4] * data_[877] + B10_current[4] * data_[868] + cB00_current[4] * data_[787];
  data_[887] = C00_[5] * data_[878] + B10_current[5] * data_[869] + cB00_current[5] * data_[788];
  data_[888] = C00_[6] * data_[879] + B10_current[6] * data_[870] + cB00_current[6] * data_[789];
  data_[889] = C00_[7] * data_[880] + B10_current[7] * data_[871] + cB00_current[7] * data_[790];
  data_[890] = C00_[8] * data_[881] + B10_current[8] * data_[872] + cB00_current[8] * data_[791];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];
  B10_current[8] += B10_[8];

  data_[891] = C00_[0] * data_[882] + B10_current[0] * data_[873] + cB00_current[0] * data_[792];
  data_[892] = C00_[1] * data_[883] + B10_current[1] * data_[874] + cB00_current[1] * data_[793];
  data_[893] = C00_[2] * data_[884] + B10_current[2] * data_[875] + cB00_current[2] * data_[794];
  data_[894] = C00_[3] * data_[885] + B10_current[3] * data_[876] + cB00_current[3] * data_[795];
  data_[895] = C00_[4] * data_[886] + B10_current[4] * data_[877] + cB00_current[4] * data_[796];
  data_[896] = C00_[5] * data_[887] + B10_current[5] * data_[878] + cB00_current[5] * data_[797];
  data_[897] = C00_[6] * data_[888] + B10_current[6] * data_[879] + cB00_current[6] * data_[798];
  data_[898] = C00_[7] * data_[889] + B10_current[7] * data_[880] + cB00_current[7] * data_[799];
  data_[899] = C00_[8] * data_[890] + B10_current[8] * data_[881] + cB00_current[8] * data_[800];
}

