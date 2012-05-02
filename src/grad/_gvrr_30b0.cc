//
// Newint - Parallel electron correlation program.
// Filename: _gvrr_30b0.cc
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

// returns double array of length 336
void GVRRList::_gvrr_30b0(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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

  double B10_current[7];
  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];

  data_[14] = C00_[0] * data_[7] + B10_current[0];
  data_[15] = C00_[1] * data_[8] + B10_current[1];
  data_[16] = C00_[2] * data_[9] + B10_current[2];
  data_[17] = C00_[3] * data_[10] + B10_current[3];
  data_[18] = C00_[4] * data_[11] + B10_current[4];
  data_[19] = C00_[5] * data_[12] + B10_current[5];
  data_[20] = C00_[6] * data_[13] + B10_current[6];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];

  data_[21] = C00_[0] * data_[14] + B10_current[0] * data_[7];
  data_[22] = C00_[1] * data_[15] + B10_current[1] * data_[8];
  data_[23] = C00_[2] * data_[16] + B10_current[2] * data_[9];
  data_[24] = C00_[3] * data_[17] + B10_current[3] * data_[10];
  data_[25] = C00_[4] * data_[18] + B10_current[4] * data_[11];
  data_[26] = C00_[5] * data_[19] + B10_current[5] * data_[12];
  data_[27] = C00_[6] * data_[20] + B10_current[6] * data_[13];

  data_[28] = D00_[0];
  data_[29] = D00_[1];
  data_[30] = D00_[2];
  data_[31] = D00_[3];
  data_[32] = D00_[4];
  data_[33] = D00_[5];
  data_[34] = D00_[6];

  double cB00_current[7];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];
  cB00_current[4] = B00_[4];
  cB00_current[5] = B00_[5];
  cB00_current[6] = B00_[6];

  data_[35] = C00_[0] * data_[28] + cB00_current[0];
  data_[36] = C00_[1] * data_[29] + cB00_current[1];
  data_[37] = C00_[2] * data_[30] + cB00_current[2];
  data_[38] = C00_[3] * data_[31] + cB00_current[3];
  data_[39] = C00_[4] * data_[32] + cB00_current[4];
  data_[40] = C00_[5] * data_[33] + cB00_current[5];
  data_[41] = C00_[6] * data_[34] + cB00_current[6];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];

  data_[42] = C00_[0] * data_[35] + B10_current[0] * data_[28] + cB00_current[0] * data_[7];
  data_[43] = C00_[1] * data_[36] + B10_current[1] * data_[29] + cB00_current[1] * data_[8];
  data_[44] = C00_[2] * data_[37] + B10_current[2] * data_[30] + cB00_current[2] * data_[9];
  data_[45] = C00_[3] * data_[38] + B10_current[3] * data_[31] + cB00_current[3] * data_[10];
  data_[46] = C00_[4] * data_[39] + B10_current[4] * data_[32] + cB00_current[4] * data_[11];
  data_[47] = C00_[5] * data_[40] + B10_current[5] * data_[33] + cB00_current[5] * data_[12];
  data_[48] = C00_[6] * data_[41] + B10_current[6] * data_[34] + cB00_current[6] * data_[13];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];

  data_[49] = C00_[0] * data_[42] + B10_current[0] * data_[35] + cB00_current[0] * data_[14];
  data_[50] = C00_[1] * data_[43] + B10_current[1] * data_[36] + cB00_current[1] * data_[15];
  data_[51] = C00_[2] * data_[44] + B10_current[2] * data_[37] + cB00_current[2] * data_[16];
  data_[52] = C00_[3] * data_[45] + B10_current[3] * data_[38] + cB00_current[3] * data_[17];
  data_[53] = C00_[4] * data_[46] + B10_current[4] * data_[39] + cB00_current[4] * data_[18];
  data_[54] = C00_[5] * data_[47] + B10_current[5] * data_[40] + cB00_current[5] * data_[19];
  data_[55] = C00_[6] * data_[48] + B10_current[6] * data_[41] + cB00_current[6] * data_[20];

  double B01_current[7];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];
  B01_current[4] = B01_[4];
  B01_current[5] = B01_[5];
  B01_current[6] = B01_[6];

  data_[56] = D00_[0] * data_[28] + B01_current[0];
  data_[57] = D00_[1] * data_[29] + B01_current[1];
  data_[58] = D00_[2] * data_[30] + B01_current[2];
  data_[59] = D00_[3] * data_[31] + B01_current[3];
  data_[60] = D00_[4] * data_[32] + B01_current[4];
  data_[61] = D00_[5] * data_[33] + B01_current[5];
  data_[62] = D00_[6] * data_[34] + B01_current[6];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[63] = C00_[0] * data_[56] + cB00_current[0] * data_[28];
  data_[64] = C00_[1] * data_[57] + cB00_current[1] * data_[29];
  data_[65] = C00_[2] * data_[58] + cB00_current[2] * data_[30];
  data_[66] = C00_[3] * data_[59] + cB00_current[3] * data_[31];
  data_[67] = C00_[4] * data_[60] + cB00_current[4] * data_[32];
  data_[68] = C00_[5] * data_[61] + cB00_current[5] * data_[33];
  data_[69] = C00_[6] * data_[62] + cB00_current[6] * data_[34];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];

  data_[70] = C00_[0] * data_[63] + B10_current[0] * data_[56] + cB00_current[0] * data_[35];
  data_[71] = C00_[1] * data_[64] + B10_current[1] * data_[57] + cB00_current[1] * data_[36];
  data_[72] = C00_[2] * data_[65] + B10_current[2] * data_[58] + cB00_current[2] * data_[37];
  data_[73] = C00_[3] * data_[66] + B10_current[3] * data_[59] + cB00_current[3] * data_[38];
  data_[74] = C00_[4] * data_[67] + B10_current[4] * data_[60] + cB00_current[4] * data_[39];
  data_[75] = C00_[5] * data_[68] + B10_current[5] * data_[61] + cB00_current[5] * data_[40];
  data_[76] = C00_[6] * data_[69] + B10_current[6] * data_[62] + cB00_current[6] * data_[41];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];

  data_[77] = C00_[0] * data_[70] + B10_current[0] * data_[63] + cB00_current[0] * data_[42];
  data_[78] = C00_[1] * data_[71] + B10_current[1] * data_[64] + cB00_current[1] * data_[43];
  data_[79] = C00_[2] * data_[72] + B10_current[2] * data_[65] + cB00_current[2] * data_[44];
  data_[80] = C00_[3] * data_[73] + B10_current[3] * data_[66] + cB00_current[3] * data_[45];
  data_[81] = C00_[4] * data_[74] + B10_current[4] * data_[67] + cB00_current[4] * data_[46];
  data_[82] = C00_[5] * data_[75] + B10_current[5] * data_[68] + cB00_current[5] * data_[47];
  data_[83] = C00_[6] * data_[76] + B10_current[6] * data_[69] + cB00_current[6] * data_[48];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[84] = D00_[0] * data_[56] + B01_current[0] * data_[28];
  data_[85] = D00_[1] * data_[57] + B01_current[1] * data_[29];
  data_[86] = D00_[2] * data_[58] + B01_current[2] * data_[30];
  data_[87] = D00_[3] * data_[59] + B01_current[3] * data_[31];
  data_[88] = D00_[4] * data_[60] + B01_current[4] * data_[32];
  data_[89] = D00_[5] * data_[61] + B01_current[5] * data_[33];
  data_[90] = D00_[6] * data_[62] + B01_current[6] * data_[34];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[91] = C00_[0] * data_[84] + cB00_current[0] * data_[56];
  data_[92] = C00_[1] * data_[85] + cB00_current[1] * data_[57];
  data_[93] = C00_[2] * data_[86] + cB00_current[2] * data_[58];
  data_[94] = C00_[3] * data_[87] + cB00_current[3] * data_[59];
  data_[95] = C00_[4] * data_[88] + cB00_current[4] * data_[60];
  data_[96] = C00_[5] * data_[89] + cB00_current[5] * data_[61];
  data_[97] = C00_[6] * data_[90] + cB00_current[6] * data_[62];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];

  data_[98] = C00_[0] * data_[91] + B10_current[0] * data_[84] + cB00_current[0] * data_[63];
  data_[99] = C00_[1] * data_[92] + B10_current[1] * data_[85] + cB00_current[1] * data_[64];
  data_[100] = C00_[2] * data_[93] + B10_current[2] * data_[86] + cB00_current[2] * data_[65];
  data_[101] = C00_[3] * data_[94] + B10_current[3] * data_[87] + cB00_current[3] * data_[66];
  data_[102] = C00_[4] * data_[95] + B10_current[4] * data_[88] + cB00_current[4] * data_[67];
  data_[103] = C00_[5] * data_[96] + B10_current[5] * data_[89] + cB00_current[5] * data_[68];
  data_[104] = C00_[6] * data_[97] + B10_current[6] * data_[90] + cB00_current[6] * data_[69];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];

  data_[105] = C00_[0] * data_[98] + B10_current[0] * data_[91] + cB00_current[0] * data_[70];
  data_[106] = C00_[1] * data_[99] + B10_current[1] * data_[92] + cB00_current[1] * data_[71];
  data_[107] = C00_[2] * data_[100] + B10_current[2] * data_[93] + cB00_current[2] * data_[72];
  data_[108] = C00_[3] * data_[101] + B10_current[3] * data_[94] + cB00_current[3] * data_[73];
  data_[109] = C00_[4] * data_[102] + B10_current[4] * data_[95] + cB00_current[4] * data_[74];
  data_[110] = C00_[5] * data_[103] + B10_current[5] * data_[96] + cB00_current[5] * data_[75];
  data_[111] = C00_[6] * data_[104] + B10_current[6] * data_[97] + cB00_current[6] * data_[76];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[112] = D00_[0] * data_[84] + B01_current[0] * data_[56];
  data_[113] = D00_[1] * data_[85] + B01_current[1] * data_[57];
  data_[114] = D00_[2] * data_[86] + B01_current[2] * data_[58];
  data_[115] = D00_[3] * data_[87] + B01_current[3] * data_[59];
  data_[116] = D00_[4] * data_[88] + B01_current[4] * data_[60];
  data_[117] = D00_[5] * data_[89] + B01_current[5] * data_[61];
  data_[118] = D00_[6] * data_[90] + B01_current[6] * data_[62];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[119] = C00_[0] * data_[112] + cB00_current[0] * data_[84];
  data_[120] = C00_[1] * data_[113] + cB00_current[1] * data_[85];
  data_[121] = C00_[2] * data_[114] + cB00_current[2] * data_[86];
  data_[122] = C00_[3] * data_[115] + cB00_current[3] * data_[87];
  data_[123] = C00_[4] * data_[116] + cB00_current[4] * data_[88];
  data_[124] = C00_[5] * data_[117] + cB00_current[5] * data_[89];
  data_[125] = C00_[6] * data_[118] + cB00_current[6] * data_[90];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];

  data_[126] = C00_[0] * data_[119] + B10_current[0] * data_[112] + cB00_current[0] * data_[91];
  data_[127] = C00_[1] * data_[120] + B10_current[1] * data_[113] + cB00_current[1] * data_[92];
  data_[128] = C00_[2] * data_[121] + B10_current[2] * data_[114] + cB00_current[2] * data_[93];
  data_[129] = C00_[3] * data_[122] + B10_current[3] * data_[115] + cB00_current[3] * data_[94];
  data_[130] = C00_[4] * data_[123] + B10_current[4] * data_[116] + cB00_current[4] * data_[95];
  data_[131] = C00_[5] * data_[124] + B10_current[5] * data_[117] + cB00_current[5] * data_[96];
  data_[132] = C00_[6] * data_[125] + B10_current[6] * data_[118] + cB00_current[6] * data_[97];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];

  data_[133] = C00_[0] * data_[126] + B10_current[0] * data_[119] + cB00_current[0] * data_[98];
  data_[134] = C00_[1] * data_[127] + B10_current[1] * data_[120] + cB00_current[1] * data_[99];
  data_[135] = C00_[2] * data_[128] + B10_current[2] * data_[121] + cB00_current[2] * data_[100];
  data_[136] = C00_[3] * data_[129] + B10_current[3] * data_[122] + cB00_current[3] * data_[101];
  data_[137] = C00_[4] * data_[130] + B10_current[4] * data_[123] + cB00_current[4] * data_[102];
  data_[138] = C00_[5] * data_[131] + B10_current[5] * data_[124] + cB00_current[5] * data_[103];
  data_[139] = C00_[6] * data_[132] + B10_current[6] * data_[125] + cB00_current[6] * data_[104];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[140] = D00_[0] * data_[112] + B01_current[0] * data_[84];
  data_[141] = D00_[1] * data_[113] + B01_current[1] * data_[85];
  data_[142] = D00_[2] * data_[114] + B01_current[2] * data_[86];
  data_[143] = D00_[3] * data_[115] + B01_current[3] * data_[87];
  data_[144] = D00_[4] * data_[116] + B01_current[4] * data_[88];
  data_[145] = D00_[5] * data_[117] + B01_current[5] * data_[89];
  data_[146] = D00_[6] * data_[118] + B01_current[6] * data_[90];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[147] = C00_[0] * data_[140] + cB00_current[0] * data_[112];
  data_[148] = C00_[1] * data_[141] + cB00_current[1] * data_[113];
  data_[149] = C00_[2] * data_[142] + cB00_current[2] * data_[114];
  data_[150] = C00_[3] * data_[143] + cB00_current[3] * data_[115];
  data_[151] = C00_[4] * data_[144] + cB00_current[4] * data_[116];
  data_[152] = C00_[5] * data_[145] + cB00_current[5] * data_[117];
  data_[153] = C00_[6] * data_[146] + cB00_current[6] * data_[118];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];

  data_[154] = C00_[0] * data_[147] + B10_current[0] * data_[140] + cB00_current[0] * data_[119];
  data_[155] = C00_[1] * data_[148] + B10_current[1] * data_[141] + cB00_current[1] * data_[120];
  data_[156] = C00_[2] * data_[149] + B10_current[2] * data_[142] + cB00_current[2] * data_[121];
  data_[157] = C00_[3] * data_[150] + B10_current[3] * data_[143] + cB00_current[3] * data_[122];
  data_[158] = C00_[4] * data_[151] + B10_current[4] * data_[144] + cB00_current[4] * data_[123];
  data_[159] = C00_[5] * data_[152] + B10_current[5] * data_[145] + cB00_current[5] * data_[124];
  data_[160] = C00_[6] * data_[153] + B10_current[6] * data_[146] + cB00_current[6] * data_[125];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];

  data_[161] = C00_[0] * data_[154] + B10_current[0] * data_[147] + cB00_current[0] * data_[126];
  data_[162] = C00_[1] * data_[155] + B10_current[1] * data_[148] + cB00_current[1] * data_[127];
  data_[163] = C00_[2] * data_[156] + B10_current[2] * data_[149] + cB00_current[2] * data_[128];
  data_[164] = C00_[3] * data_[157] + B10_current[3] * data_[150] + cB00_current[3] * data_[129];
  data_[165] = C00_[4] * data_[158] + B10_current[4] * data_[151] + cB00_current[4] * data_[130];
  data_[166] = C00_[5] * data_[159] + B10_current[5] * data_[152] + cB00_current[5] * data_[131];
  data_[167] = C00_[6] * data_[160] + B10_current[6] * data_[153] + cB00_current[6] * data_[132];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[168] = D00_[0] * data_[140] + B01_current[0] * data_[112];
  data_[169] = D00_[1] * data_[141] + B01_current[1] * data_[113];
  data_[170] = D00_[2] * data_[142] + B01_current[2] * data_[114];
  data_[171] = D00_[3] * data_[143] + B01_current[3] * data_[115];
  data_[172] = D00_[4] * data_[144] + B01_current[4] * data_[116];
  data_[173] = D00_[5] * data_[145] + B01_current[5] * data_[117];
  data_[174] = D00_[6] * data_[146] + B01_current[6] * data_[118];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[175] = C00_[0] * data_[168] + cB00_current[0] * data_[140];
  data_[176] = C00_[1] * data_[169] + cB00_current[1] * data_[141];
  data_[177] = C00_[2] * data_[170] + cB00_current[2] * data_[142];
  data_[178] = C00_[3] * data_[171] + cB00_current[3] * data_[143];
  data_[179] = C00_[4] * data_[172] + cB00_current[4] * data_[144];
  data_[180] = C00_[5] * data_[173] + cB00_current[5] * data_[145];
  data_[181] = C00_[6] * data_[174] + cB00_current[6] * data_[146];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];

  data_[182] = C00_[0] * data_[175] + B10_current[0] * data_[168] + cB00_current[0] * data_[147];
  data_[183] = C00_[1] * data_[176] + B10_current[1] * data_[169] + cB00_current[1] * data_[148];
  data_[184] = C00_[2] * data_[177] + B10_current[2] * data_[170] + cB00_current[2] * data_[149];
  data_[185] = C00_[3] * data_[178] + B10_current[3] * data_[171] + cB00_current[3] * data_[150];
  data_[186] = C00_[4] * data_[179] + B10_current[4] * data_[172] + cB00_current[4] * data_[151];
  data_[187] = C00_[5] * data_[180] + B10_current[5] * data_[173] + cB00_current[5] * data_[152];
  data_[188] = C00_[6] * data_[181] + B10_current[6] * data_[174] + cB00_current[6] * data_[153];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];

  data_[189] = C00_[0] * data_[182] + B10_current[0] * data_[175] + cB00_current[0] * data_[154];
  data_[190] = C00_[1] * data_[183] + B10_current[1] * data_[176] + cB00_current[1] * data_[155];
  data_[191] = C00_[2] * data_[184] + B10_current[2] * data_[177] + cB00_current[2] * data_[156];
  data_[192] = C00_[3] * data_[185] + B10_current[3] * data_[178] + cB00_current[3] * data_[157];
  data_[193] = C00_[4] * data_[186] + B10_current[4] * data_[179] + cB00_current[4] * data_[158];
  data_[194] = C00_[5] * data_[187] + B10_current[5] * data_[180] + cB00_current[5] * data_[159];
  data_[195] = C00_[6] * data_[188] + B10_current[6] * data_[181] + cB00_current[6] * data_[160];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[196] = D00_[0] * data_[168] + B01_current[0] * data_[140];
  data_[197] = D00_[1] * data_[169] + B01_current[1] * data_[141];
  data_[198] = D00_[2] * data_[170] + B01_current[2] * data_[142];
  data_[199] = D00_[3] * data_[171] + B01_current[3] * data_[143];
  data_[200] = D00_[4] * data_[172] + B01_current[4] * data_[144];
  data_[201] = D00_[5] * data_[173] + B01_current[5] * data_[145];
  data_[202] = D00_[6] * data_[174] + B01_current[6] * data_[146];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[203] = C00_[0] * data_[196] + cB00_current[0] * data_[168];
  data_[204] = C00_[1] * data_[197] + cB00_current[1] * data_[169];
  data_[205] = C00_[2] * data_[198] + cB00_current[2] * data_[170];
  data_[206] = C00_[3] * data_[199] + cB00_current[3] * data_[171];
  data_[207] = C00_[4] * data_[200] + cB00_current[4] * data_[172];
  data_[208] = C00_[5] * data_[201] + cB00_current[5] * data_[173];
  data_[209] = C00_[6] * data_[202] + cB00_current[6] * data_[174];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];

  data_[210] = C00_[0] * data_[203] + B10_current[0] * data_[196] + cB00_current[0] * data_[175];
  data_[211] = C00_[1] * data_[204] + B10_current[1] * data_[197] + cB00_current[1] * data_[176];
  data_[212] = C00_[2] * data_[205] + B10_current[2] * data_[198] + cB00_current[2] * data_[177];
  data_[213] = C00_[3] * data_[206] + B10_current[3] * data_[199] + cB00_current[3] * data_[178];
  data_[214] = C00_[4] * data_[207] + B10_current[4] * data_[200] + cB00_current[4] * data_[179];
  data_[215] = C00_[5] * data_[208] + B10_current[5] * data_[201] + cB00_current[5] * data_[180];
  data_[216] = C00_[6] * data_[209] + B10_current[6] * data_[202] + cB00_current[6] * data_[181];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];

  data_[217] = C00_[0] * data_[210] + B10_current[0] * data_[203] + cB00_current[0] * data_[182];
  data_[218] = C00_[1] * data_[211] + B10_current[1] * data_[204] + cB00_current[1] * data_[183];
  data_[219] = C00_[2] * data_[212] + B10_current[2] * data_[205] + cB00_current[2] * data_[184];
  data_[220] = C00_[3] * data_[213] + B10_current[3] * data_[206] + cB00_current[3] * data_[185];
  data_[221] = C00_[4] * data_[214] + B10_current[4] * data_[207] + cB00_current[4] * data_[186];
  data_[222] = C00_[5] * data_[215] + B10_current[5] * data_[208] + cB00_current[5] * data_[187];
  data_[223] = C00_[6] * data_[216] + B10_current[6] * data_[209] + cB00_current[6] * data_[188];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[224] = D00_[0] * data_[196] + B01_current[0] * data_[168];
  data_[225] = D00_[1] * data_[197] + B01_current[1] * data_[169];
  data_[226] = D00_[2] * data_[198] + B01_current[2] * data_[170];
  data_[227] = D00_[3] * data_[199] + B01_current[3] * data_[171];
  data_[228] = D00_[4] * data_[200] + B01_current[4] * data_[172];
  data_[229] = D00_[5] * data_[201] + B01_current[5] * data_[173];
  data_[230] = D00_[6] * data_[202] + B01_current[6] * data_[174];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[231] = C00_[0] * data_[224] + cB00_current[0] * data_[196];
  data_[232] = C00_[1] * data_[225] + cB00_current[1] * data_[197];
  data_[233] = C00_[2] * data_[226] + cB00_current[2] * data_[198];
  data_[234] = C00_[3] * data_[227] + cB00_current[3] * data_[199];
  data_[235] = C00_[4] * data_[228] + cB00_current[4] * data_[200];
  data_[236] = C00_[5] * data_[229] + cB00_current[5] * data_[201];
  data_[237] = C00_[6] * data_[230] + cB00_current[6] * data_[202];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];

  data_[238] = C00_[0] * data_[231] + B10_current[0] * data_[224] + cB00_current[0] * data_[203];
  data_[239] = C00_[1] * data_[232] + B10_current[1] * data_[225] + cB00_current[1] * data_[204];
  data_[240] = C00_[2] * data_[233] + B10_current[2] * data_[226] + cB00_current[2] * data_[205];
  data_[241] = C00_[3] * data_[234] + B10_current[3] * data_[227] + cB00_current[3] * data_[206];
  data_[242] = C00_[4] * data_[235] + B10_current[4] * data_[228] + cB00_current[4] * data_[207];
  data_[243] = C00_[5] * data_[236] + B10_current[5] * data_[229] + cB00_current[5] * data_[208];
  data_[244] = C00_[6] * data_[237] + B10_current[6] * data_[230] + cB00_current[6] * data_[209];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];

  data_[245] = C00_[0] * data_[238] + B10_current[0] * data_[231] + cB00_current[0] * data_[210];
  data_[246] = C00_[1] * data_[239] + B10_current[1] * data_[232] + cB00_current[1] * data_[211];
  data_[247] = C00_[2] * data_[240] + B10_current[2] * data_[233] + cB00_current[2] * data_[212];
  data_[248] = C00_[3] * data_[241] + B10_current[3] * data_[234] + cB00_current[3] * data_[213];
  data_[249] = C00_[4] * data_[242] + B10_current[4] * data_[235] + cB00_current[4] * data_[214];
  data_[250] = C00_[5] * data_[243] + B10_current[5] * data_[236] + cB00_current[5] * data_[215];
  data_[251] = C00_[6] * data_[244] + B10_current[6] * data_[237] + cB00_current[6] * data_[216];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[252] = D00_[0] * data_[224] + B01_current[0] * data_[196];
  data_[253] = D00_[1] * data_[225] + B01_current[1] * data_[197];
  data_[254] = D00_[2] * data_[226] + B01_current[2] * data_[198];
  data_[255] = D00_[3] * data_[227] + B01_current[3] * data_[199];
  data_[256] = D00_[4] * data_[228] + B01_current[4] * data_[200];
  data_[257] = D00_[5] * data_[229] + B01_current[5] * data_[201];
  data_[258] = D00_[6] * data_[230] + B01_current[6] * data_[202];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[259] = C00_[0] * data_[252] + cB00_current[0] * data_[224];
  data_[260] = C00_[1] * data_[253] + cB00_current[1] * data_[225];
  data_[261] = C00_[2] * data_[254] + cB00_current[2] * data_[226];
  data_[262] = C00_[3] * data_[255] + cB00_current[3] * data_[227];
  data_[263] = C00_[4] * data_[256] + cB00_current[4] * data_[228];
  data_[264] = C00_[5] * data_[257] + cB00_current[5] * data_[229];
  data_[265] = C00_[6] * data_[258] + cB00_current[6] * data_[230];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];

  data_[266] = C00_[0] * data_[259] + B10_current[0] * data_[252] + cB00_current[0] * data_[231];
  data_[267] = C00_[1] * data_[260] + B10_current[1] * data_[253] + cB00_current[1] * data_[232];
  data_[268] = C00_[2] * data_[261] + B10_current[2] * data_[254] + cB00_current[2] * data_[233];
  data_[269] = C00_[3] * data_[262] + B10_current[3] * data_[255] + cB00_current[3] * data_[234];
  data_[270] = C00_[4] * data_[263] + B10_current[4] * data_[256] + cB00_current[4] * data_[235];
  data_[271] = C00_[5] * data_[264] + B10_current[5] * data_[257] + cB00_current[5] * data_[236];
  data_[272] = C00_[6] * data_[265] + B10_current[6] * data_[258] + cB00_current[6] * data_[237];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];

  data_[273] = C00_[0] * data_[266] + B10_current[0] * data_[259] + cB00_current[0] * data_[238];
  data_[274] = C00_[1] * data_[267] + B10_current[1] * data_[260] + cB00_current[1] * data_[239];
  data_[275] = C00_[2] * data_[268] + B10_current[2] * data_[261] + cB00_current[2] * data_[240];
  data_[276] = C00_[3] * data_[269] + B10_current[3] * data_[262] + cB00_current[3] * data_[241];
  data_[277] = C00_[4] * data_[270] + B10_current[4] * data_[263] + cB00_current[4] * data_[242];
  data_[278] = C00_[5] * data_[271] + B10_current[5] * data_[264] + cB00_current[5] * data_[243];
  data_[279] = C00_[6] * data_[272] + B10_current[6] * data_[265] + cB00_current[6] * data_[244];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[280] = D00_[0] * data_[252] + B01_current[0] * data_[224];
  data_[281] = D00_[1] * data_[253] + B01_current[1] * data_[225];
  data_[282] = D00_[2] * data_[254] + B01_current[2] * data_[226];
  data_[283] = D00_[3] * data_[255] + B01_current[3] * data_[227];
  data_[284] = D00_[4] * data_[256] + B01_current[4] * data_[228];
  data_[285] = D00_[5] * data_[257] + B01_current[5] * data_[229];
  data_[286] = D00_[6] * data_[258] + B01_current[6] * data_[230];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[287] = C00_[0] * data_[280] + cB00_current[0] * data_[252];
  data_[288] = C00_[1] * data_[281] + cB00_current[1] * data_[253];
  data_[289] = C00_[2] * data_[282] + cB00_current[2] * data_[254];
  data_[290] = C00_[3] * data_[283] + cB00_current[3] * data_[255];
  data_[291] = C00_[4] * data_[284] + cB00_current[4] * data_[256];
  data_[292] = C00_[5] * data_[285] + cB00_current[5] * data_[257];
  data_[293] = C00_[6] * data_[286] + cB00_current[6] * data_[258];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];

  data_[294] = C00_[0] * data_[287] + B10_current[0] * data_[280] + cB00_current[0] * data_[259];
  data_[295] = C00_[1] * data_[288] + B10_current[1] * data_[281] + cB00_current[1] * data_[260];
  data_[296] = C00_[2] * data_[289] + B10_current[2] * data_[282] + cB00_current[2] * data_[261];
  data_[297] = C00_[3] * data_[290] + B10_current[3] * data_[283] + cB00_current[3] * data_[262];
  data_[298] = C00_[4] * data_[291] + B10_current[4] * data_[284] + cB00_current[4] * data_[263];
  data_[299] = C00_[5] * data_[292] + B10_current[5] * data_[285] + cB00_current[5] * data_[264];
  data_[300] = C00_[6] * data_[293] + B10_current[6] * data_[286] + cB00_current[6] * data_[265];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];

  data_[301] = C00_[0] * data_[294] + B10_current[0] * data_[287] + cB00_current[0] * data_[266];
  data_[302] = C00_[1] * data_[295] + B10_current[1] * data_[288] + cB00_current[1] * data_[267];
  data_[303] = C00_[2] * data_[296] + B10_current[2] * data_[289] + cB00_current[2] * data_[268];
  data_[304] = C00_[3] * data_[297] + B10_current[3] * data_[290] + cB00_current[3] * data_[269];
  data_[305] = C00_[4] * data_[298] + B10_current[4] * data_[291] + cB00_current[4] * data_[270];
  data_[306] = C00_[5] * data_[299] + B10_current[5] * data_[292] + cB00_current[5] * data_[271];
  data_[307] = C00_[6] * data_[300] + B10_current[6] * data_[293] + cB00_current[6] * data_[272];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];

  data_[308] = D00_[0] * data_[280] + B01_current[0] * data_[252];
  data_[309] = D00_[1] * data_[281] + B01_current[1] * data_[253];
  data_[310] = D00_[2] * data_[282] + B01_current[2] * data_[254];
  data_[311] = D00_[3] * data_[283] + B01_current[3] * data_[255];
  data_[312] = D00_[4] * data_[284] + B01_current[4] * data_[256];
  data_[313] = D00_[5] * data_[285] + B01_current[5] * data_[257];
  data_[314] = D00_[6] * data_[286] + B01_current[6] * data_[258];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];

  data_[315] = C00_[0] * data_[308] + cB00_current[0] * data_[280];
  data_[316] = C00_[1] * data_[309] + cB00_current[1] * data_[281];
  data_[317] = C00_[2] * data_[310] + cB00_current[2] * data_[282];
  data_[318] = C00_[3] * data_[311] + cB00_current[3] * data_[283];
  data_[319] = C00_[4] * data_[312] + cB00_current[4] * data_[284];
  data_[320] = C00_[5] * data_[313] + cB00_current[5] * data_[285];
  data_[321] = C00_[6] * data_[314] + cB00_current[6] * data_[286];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];

  data_[322] = C00_[0] * data_[315] + B10_current[0] * data_[308] + cB00_current[0] * data_[287];
  data_[323] = C00_[1] * data_[316] + B10_current[1] * data_[309] + cB00_current[1] * data_[288];
  data_[324] = C00_[2] * data_[317] + B10_current[2] * data_[310] + cB00_current[2] * data_[289];
  data_[325] = C00_[3] * data_[318] + B10_current[3] * data_[311] + cB00_current[3] * data_[290];
  data_[326] = C00_[4] * data_[319] + B10_current[4] * data_[312] + cB00_current[4] * data_[291];
  data_[327] = C00_[5] * data_[320] + B10_current[5] * data_[313] + cB00_current[5] * data_[292];
  data_[328] = C00_[6] * data_[321] + B10_current[6] * data_[314] + cB00_current[6] * data_[293];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];

  data_[329] = C00_[0] * data_[322] + B10_current[0] * data_[315] + cB00_current[0] * data_[294];
  data_[330] = C00_[1] * data_[323] + B10_current[1] * data_[316] + cB00_current[1] * data_[295];
  data_[331] = C00_[2] * data_[324] + B10_current[2] * data_[317] + cB00_current[2] * data_[296];
  data_[332] = C00_[3] * data_[325] + B10_current[3] * data_[318] + cB00_current[3] * data_[297];
  data_[333] = C00_[4] * data_[326] + B10_current[4] * data_[319] + cB00_current[4] * data_[298];
  data_[334] = C00_[5] * data_[327] + B10_current[5] * data_[320] + cB00_current[5] * data_[299];
  data_[335] = C00_[6] * data_[328] + B10_current[6] * data_[321] + cB00_current[6] * data_[300];
}

