//
// Newint - Parallel electron correlation program.
// Filename: _gvrr_10d0.cc
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

// returns double array of length 196
void GVRRList::_gvrr_10d0(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;
  data_[3] = 1.0;
  data_[4] = 1.0;
  data_[5] = 1.0;
  data_[6] = 1.0;

  data_[7] = C00_[0];
  data_[8] = C00_[1];
  data_[9] = C00_[2];
  data_[10] = C00_[3];
  data_[11] = C00_[4];
  data_[12] = C00_[5];
  data_[13] = C00_[6];

  data_[14] = D00_[0];
  data_[15] = D00_[1];
  data_[16] = D00_[2];
  data_[17] = D00_[3];
  data_[18] = D00_[4];
  data_[19] = D00_[5];
  data_[20] = D00_[6];

  double cB00_current[7];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];
  cB00_current[4] = B00_[4];
  cB00_current[5] = B00_[5];
  cB00_current[6] = B00_[6];

  data_[21] = C00_[0] * data_[14] + cB00_current[0];
  data_[22] = C00_[1] * data_[15] + cB00_current[1];
  data_[23] = C00_[2] * data_[16] + cB00_current[2];
  data_[24] = C00_[3] * data_[17] + cB00_current[3];
  data_[25] = C00_[4] * data_[18] + cB00_current[4];
  data_[26] = C00_[5] * data_[19] + cB00_current[5];
  data_[27] = C00_[6] * data_[20] + cB00_current[6];

  double B01_current[7];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];
  B01_current[4] = B01_[4];
  B01_current[5] = B01_[5];
  B01_current[6] = B01_[6];

  data_[28] = D00_[0] * data_[14] + B01_current[0];
  data_[29] = D00_[1] * data_[15] + B01_current[1];
  data_[30] = D00_[2] * data_[16] + B01_current[2];
  data_[31] = D00_[3] * data_[17] + B01_current[3];
  data_[32] = D00_[4] * data_[18] + B01_current[4];
  data_[33] = D00_[5] * data_[19] + B01_current[5];
  data_[34] = D00_[6] * data_[20] + B01_current[6];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[35] = C00_[0] * data_[28] + cB00_current[0] * data_[14];
  data_[36] = C00_[1] * data_[29] + cB00_current[1] * data_[15];
  data_[37] = C00_[2] * data_[30] + cB00_current[2] * data_[16];
  data_[38] = C00_[3] * data_[31] + cB00_current[3] * data_[17];
  data_[39] = C00_[4] * data_[32] + cB00_current[4] * data_[18];
  data_[40] = C00_[5] * data_[33] + cB00_current[5] * data_[19];
  data_[41] = C00_[6] * data_[34] + cB00_current[6] * data_[20];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[42] = D00_[0] * data_[28] + B01_current[0] * data_[14];
  data_[43] = D00_[1] * data_[29] + B01_current[1] * data_[15];
  data_[44] = D00_[2] * data_[30] + B01_current[2] * data_[16];
  data_[45] = D00_[3] * data_[31] + B01_current[3] * data_[17];
  data_[46] = D00_[4] * data_[32] + B01_current[4] * data_[18];
  data_[47] = D00_[5] * data_[33] + B01_current[5] * data_[19];
  data_[48] = D00_[6] * data_[34] + B01_current[6] * data_[20];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[49] = C00_[0] * data_[42] + cB00_current[0] * data_[28];
  data_[50] = C00_[1] * data_[43] + cB00_current[1] * data_[29];
  data_[51] = C00_[2] * data_[44] + cB00_current[2] * data_[30];
  data_[52] = C00_[3] * data_[45] + cB00_current[3] * data_[31];
  data_[53] = C00_[4] * data_[46] + cB00_current[4] * data_[32];
  data_[54] = C00_[5] * data_[47] + cB00_current[5] * data_[33];
  data_[55] = C00_[6] * data_[48] + cB00_current[6] * data_[34];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[56] = D00_[0] * data_[42] + B01_current[0] * data_[28];
  data_[57] = D00_[1] * data_[43] + B01_current[1] * data_[29];
  data_[58] = D00_[2] * data_[44] + B01_current[2] * data_[30];
  data_[59] = D00_[3] * data_[45] + B01_current[3] * data_[31];
  data_[60] = D00_[4] * data_[46] + B01_current[4] * data_[32];
  data_[61] = D00_[5] * data_[47] + B01_current[5] * data_[33];
  data_[62] = D00_[6] * data_[48] + B01_current[6] * data_[34];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[63] = C00_[0] * data_[56] + cB00_current[0] * data_[42];
  data_[64] = C00_[1] * data_[57] + cB00_current[1] * data_[43];
  data_[65] = C00_[2] * data_[58] + cB00_current[2] * data_[44];
  data_[66] = C00_[3] * data_[59] + cB00_current[3] * data_[45];
  data_[67] = C00_[4] * data_[60] + cB00_current[4] * data_[46];
  data_[68] = C00_[5] * data_[61] + cB00_current[5] * data_[47];
  data_[69] = C00_[6] * data_[62] + cB00_current[6] * data_[48];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[70] = D00_[0] * data_[56] + B01_current[0] * data_[42];
  data_[71] = D00_[1] * data_[57] + B01_current[1] * data_[43];
  data_[72] = D00_[2] * data_[58] + B01_current[2] * data_[44];
  data_[73] = D00_[3] * data_[59] + B01_current[3] * data_[45];
  data_[74] = D00_[4] * data_[60] + B01_current[4] * data_[46];
  data_[75] = D00_[5] * data_[61] + B01_current[5] * data_[47];
  data_[76] = D00_[6] * data_[62] + B01_current[6] * data_[48];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[77] = C00_[0] * data_[70] + cB00_current[0] * data_[56];
  data_[78] = C00_[1] * data_[71] + cB00_current[1] * data_[57];
  data_[79] = C00_[2] * data_[72] + cB00_current[2] * data_[58];
  data_[80] = C00_[3] * data_[73] + cB00_current[3] * data_[59];
  data_[81] = C00_[4] * data_[74] + cB00_current[4] * data_[60];
  data_[82] = C00_[5] * data_[75] + cB00_current[5] * data_[61];
  data_[83] = C00_[6] * data_[76] + cB00_current[6] * data_[62];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[84] = D00_[0] * data_[70] + B01_current[0] * data_[56];
  data_[85] = D00_[1] * data_[71] + B01_current[1] * data_[57];
  data_[86] = D00_[2] * data_[72] + B01_current[2] * data_[58];
  data_[87] = D00_[3] * data_[73] + B01_current[3] * data_[59];
  data_[88] = D00_[4] * data_[74] + B01_current[4] * data_[60];
  data_[89] = D00_[5] * data_[75] + B01_current[5] * data_[61];
  data_[90] = D00_[6] * data_[76] + B01_current[6] * data_[62];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[91] = C00_[0] * data_[84] + cB00_current[0] * data_[70];
  data_[92] = C00_[1] * data_[85] + cB00_current[1] * data_[71];
  data_[93] = C00_[2] * data_[86] + cB00_current[2] * data_[72];
  data_[94] = C00_[3] * data_[87] + cB00_current[3] * data_[73];
  data_[95] = C00_[4] * data_[88] + cB00_current[4] * data_[74];
  data_[96] = C00_[5] * data_[89] + cB00_current[5] * data_[75];
  data_[97] = C00_[6] * data_[90] + cB00_current[6] * data_[76];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[98] = D00_[0] * data_[84] + B01_current[0] * data_[70];
  data_[99] = D00_[1] * data_[85] + B01_current[1] * data_[71];
  data_[100] = D00_[2] * data_[86] + B01_current[2] * data_[72];
  data_[101] = D00_[3] * data_[87] + B01_current[3] * data_[73];
  data_[102] = D00_[4] * data_[88] + B01_current[4] * data_[74];
  data_[103] = D00_[5] * data_[89] + B01_current[5] * data_[75];
  data_[104] = D00_[6] * data_[90] + B01_current[6] * data_[76];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[105] = C00_[0] * data_[98] + cB00_current[0] * data_[84];
  data_[106] = C00_[1] * data_[99] + cB00_current[1] * data_[85];
  data_[107] = C00_[2] * data_[100] + cB00_current[2] * data_[86];
  data_[108] = C00_[3] * data_[101] + cB00_current[3] * data_[87];
  data_[109] = C00_[4] * data_[102] + cB00_current[4] * data_[88];
  data_[110] = C00_[5] * data_[103] + cB00_current[5] * data_[89];
  data_[111] = C00_[6] * data_[104] + cB00_current[6] * data_[90];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[112] = D00_[0] * data_[98] + B01_current[0] * data_[84];
  data_[113] = D00_[1] * data_[99] + B01_current[1] * data_[85];
  data_[114] = D00_[2] * data_[100] + B01_current[2] * data_[86];
  data_[115] = D00_[3] * data_[101] + B01_current[3] * data_[87];
  data_[116] = D00_[4] * data_[102] + B01_current[4] * data_[88];
  data_[117] = D00_[5] * data_[103] + B01_current[5] * data_[89];
  data_[118] = D00_[6] * data_[104] + B01_current[6] * data_[90];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[119] = C00_[0] * data_[112] + cB00_current[0] * data_[98];
  data_[120] = C00_[1] * data_[113] + cB00_current[1] * data_[99];
  data_[121] = C00_[2] * data_[114] + cB00_current[2] * data_[100];
  data_[122] = C00_[3] * data_[115] + cB00_current[3] * data_[101];
  data_[123] = C00_[4] * data_[116] + cB00_current[4] * data_[102];
  data_[124] = C00_[5] * data_[117] + cB00_current[5] * data_[103];
  data_[125] = C00_[6] * data_[118] + cB00_current[6] * data_[104];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[126] = D00_[0] * data_[112] + B01_current[0] * data_[98];
  data_[127] = D00_[1] * data_[113] + B01_current[1] * data_[99];
  data_[128] = D00_[2] * data_[114] + B01_current[2] * data_[100];
  data_[129] = D00_[3] * data_[115] + B01_current[3] * data_[101];
  data_[130] = D00_[4] * data_[116] + B01_current[4] * data_[102];
  data_[131] = D00_[5] * data_[117] + B01_current[5] * data_[103];
  data_[132] = D00_[6] * data_[118] + B01_current[6] * data_[104];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[133] = C00_[0] * data_[126] + cB00_current[0] * data_[112];
  data_[134] = C00_[1] * data_[127] + cB00_current[1] * data_[113];
  data_[135] = C00_[2] * data_[128] + cB00_current[2] * data_[114];
  data_[136] = C00_[3] * data_[129] + cB00_current[3] * data_[115];
  data_[137] = C00_[4] * data_[130] + cB00_current[4] * data_[116];
  data_[138] = C00_[5] * data_[131] + cB00_current[5] * data_[117];
  data_[139] = C00_[6] * data_[132] + cB00_current[6] * data_[118];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[140] = D00_[0] * data_[126] + B01_current[0] * data_[112];
  data_[141] = D00_[1] * data_[127] + B01_current[1] * data_[113];
  data_[142] = D00_[2] * data_[128] + B01_current[2] * data_[114];
  data_[143] = D00_[3] * data_[129] + B01_current[3] * data_[115];
  data_[144] = D00_[4] * data_[130] + B01_current[4] * data_[116];
  data_[145] = D00_[5] * data_[131] + B01_current[5] * data_[117];
  data_[146] = D00_[6] * data_[132] + B01_current[6] * data_[118];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[147] = C00_[0] * data_[140] + cB00_current[0] * data_[126];
  data_[148] = C00_[1] * data_[141] + cB00_current[1] * data_[127];
  data_[149] = C00_[2] * data_[142] + cB00_current[2] * data_[128];
  data_[150] = C00_[3] * data_[143] + cB00_current[3] * data_[129];
  data_[151] = C00_[4] * data_[144] + cB00_current[4] * data_[130];
  data_[152] = C00_[5] * data_[145] + cB00_current[5] * data_[131];
  data_[153] = C00_[6] * data_[146] + cB00_current[6] * data_[132];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[154] = D00_[0] * data_[140] + B01_current[0] * data_[126];
  data_[155] = D00_[1] * data_[141] + B01_current[1] * data_[127];
  data_[156] = D00_[2] * data_[142] + B01_current[2] * data_[128];
  data_[157] = D00_[3] * data_[143] + B01_current[3] * data_[129];
  data_[158] = D00_[4] * data_[144] + B01_current[4] * data_[130];
  data_[159] = D00_[5] * data_[145] + B01_current[5] * data_[131];
  data_[160] = D00_[6] * data_[146] + B01_current[6] * data_[132];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[161] = C00_[0] * data_[154] + cB00_current[0] * data_[140];
  data_[162] = C00_[1] * data_[155] + cB00_current[1] * data_[141];
  data_[163] = C00_[2] * data_[156] + cB00_current[2] * data_[142];
  data_[164] = C00_[3] * data_[157] + cB00_current[3] * data_[143];
  data_[165] = C00_[4] * data_[158] + cB00_current[4] * data_[144];
  data_[166] = C00_[5] * data_[159] + cB00_current[5] * data_[145];
  data_[167] = C00_[6] * data_[160] + cB00_current[6] * data_[146];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[168] = D00_[0] * data_[154] + B01_current[0] * data_[140];
  data_[169] = D00_[1] * data_[155] + B01_current[1] * data_[141];
  data_[170] = D00_[2] * data_[156] + B01_current[2] * data_[142];
  data_[171] = D00_[3] * data_[157] + B01_current[3] * data_[143];
  data_[172] = D00_[4] * data_[158] + B01_current[4] * data_[144];
  data_[173] = D00_[5] * data_[159] + B01_current[5] * data_[145];
  data_[174] = D00_[6] * data_[160] + B01_current[6] * data_[146];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[175] = C00_[0] * data_[168] + cB00_current[0] * data_[154];
  data_[176] = C00_[1] * data_[169] + cB00_current[1] * data_[155];
  data_[177] = C00_[2] * data_[170] + cB00_current[2] * data_[156];
  data_[178] = C00_[3] * data_[171] + cB00_current[3] * data_[157];
  data_[179] = C00_[4] * data_[172] + cB00_current[4] * data_[158];
  data_[180] = C00_[5] * data_[173] + cB00_current[5] * data_[159];
  data_[181] = C00_[6] * data_[174] + cB00_current[6] * data_[160];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[182] = D00_[0] * data_[168] + B01_current[0] * data_[154];
  data_[183] = D00_[1] * data_[169] + B01_current[1] * data_[155];
  data_[184] = D00_[2] * data_[170] + B01_current[2] * data_[156];
  data_[185] = D00_[3] * data_[171] + B01_current[3] * data_[157];
  data_[186] = D00_[4] * data_[172] + B01_current[4] * data_[158];
  data_[187] = D00_[5] * data_[173] + B01_current[5] * data_[159];
  data_[188] = D00_[6] * data_[174] + B01_current[6] * data_[160];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[189] = C00_[0] * data_[182] + cB00_current[0] * data_[168];
  data_[190] = C00_[1] * data_[183] + cB00_current[1] * data_[169];
  data_[191] = C00_[2] * data_[184] + cB00_current[2] * data_[170];
  data_[192] = C00_[3] * data_[185] + cB00_current[3] * data_[171];
  data_[193] = C00_[4] * data_[186] + cB00_current[4] * data_[172];
  data_[194] = C00_[5] * data_[187] + cB00_current[5] * data_[173];
  data_[195] = C00_[6] * data_[188] + cB00_current[6] * data_[174];
}

