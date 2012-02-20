//
// Newint - Parallel electron correlation program.
// Filename: _vrr_7070.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


#include "vrrlist.h"

// returns double array of length 512
void VRRList::_vrr_7070(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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

  data_[64] = D00_[0];
  data_[65] = D00_[1];
  data_[66] = D00_[2];
  data_[67] = D00_[3];
  data_[68] = D00_[4];
  data_[69] = D00_[5];
  data_[70] = D00_[6];
  data_[71] = D00_[7];

  double cB00_current[8];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];
  cB00_current[4] = B00_[4];
  cB00_current[5] = B00_[5];
  cB00_current[6] = B00_[6];
  cB00_current[7] = B00_[7];

  data_[72] = C00_[0] * data_[64] + cB00_current[0];
  data_[73] = C00_[1] * data_[65] + cB00_current[1];
  data_[74] = C00_[2] * data_[66] + cB00_current[2];
  data_[75] = C00_[3] * data_[67] + cB00_current[3];
  data_[76] = C00_[4] * data_[68] + cB00_current[4];
  data_[77] = C00_[5] * data_[69] + cB00_current[5];
  data_[78] = C00_[6] * data_[70] + cB00_current[6];
  data_[79] = C00_[7] * data_[71] + cB00_current[7];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[80] = C00_[0] * data_[72] + B10_current[0] * data_[64] + cB00_current[0] * data_[8];
  data_[81] = C00_[1] * data_[73] + B10_current[1] * data_[65] + cB00_current[1] * data_[9];
  data_[82] = C00_[2] * data_[74] + B10_current[2] * data_[66] + cB00_current[2] * data_[10];
  data_[83] = C00_[3] * data_[75] + B10_current[3] * data_[67] + cB00_current[3] * data_[11];
  data_[84] = C00_[4] * data_[76] + B10_current[4] * data_[68] + cB00_current[4] * data_[12];
  data_[85] = C00_[5] * data_[77] + B10_current[5] * data_[69] + cB00_current[5] * data_[13];
  data_[86] = C00_[6] * data_[78] + B10_current[6] * data_[70] + cB00_current[6] * data_[14];
  data_[87] = C00_[7] * data_[79] + B10_current[7] * data_[71] + cB00_current[7] * data_[15];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[88] = C00_[0] * data_[80] + B10_current[0] * data_[72] + cB00_current[0] * data_[16];
  data_[89] = C00_[1] * data_[81] + B10_current[1] * data_[73] + cB00_current[1] * data_[17];
  data_[90] = C00_[2] * data_[82] + B10_current[2] * data_[74] + cB00_current[2] * data_[18];
  data_[91] = C00_[3] * data_[83] + B10_current[3] * data_[75] + cB00_current[3] * data_[19];
  data_[92] = C00_[4] * data_[84] + B10_current[4] * data_[76] + cB00_current[4] * data_[20];
  data_[93] = C00_[5] * data_[85] + B10_current[5] * data_[77] + cB00_current[5] * data_[21];
  data_[94] = C00_[6] * data_[86] + B10_current[6] * data_[78] + cB00_current[6] * data_[22];
  data_[95] = C00_[7] * data_[87] + B10_current[7] * data_[79] + cB00_current[7] * data_[23];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[96] = C00_[0] * data_[88] + B10_current[0] * data_[80] + cB00_current[0] * data_[24];
  data_[97] = C00_[1] * data_[89] + B10_current[1] * data_[81] + cB00_current[1] * data_[25];
  data_[98] = C00_[2] * data_[90] + B10_current[2] * data_[82] + cB00_current[2] * data_[26];
  data_[99] = C00_[3] * data_[91] + B10_current[3] * data_[83] + cB00_current[3] * data_[27];
  data_[100] = C00_[4] * data_[92] + B10_current[4] * data_[84] + cB00_current[4] * data_[28];
  data_[101] = C00_[5] * data_[93] + B10_current[5] * data_[85] + cB00_current[5] * data_[29];
  data_[102] = C00_[6] * data_[94] + B10_current[6] * data_[86] + cB00_current[6] * data_[30];
  data_[103] = C00_[7] * data_[95] + B10_current[7] * data_[87] + cB00_current[7] * data_[31];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[104] = C00_[0] * data_[96] + B10_current[0] * data_[88] + cB00_current[0] * data_[32];
  data_[105] = C00_[1] * data_[97] + B10_current[1] * data_[89] + cB00_current[1] * data_[33];
  data_[106] = C00_[2] * data_[98] + B10_current[2] * data_[90] + cB00_current[2] * data_[34];
  data_[107] = C00_[3] * data_[99] + B10_current[3] * data_[91] + cB00_current[3] * data_[35];
  data_[108] = C00_[4] * data_[100] + B10_current[4] * data_[92] + cB00_current[4] * data_[36];
  data_[109] = C00_[5] * data_[101] + B10_current[5] * data_[93] + cB00_current[5] * data_[37];
  data_[110] = C00_[6] * data_[102] + B10_current[6] * data_[94] + cB00_current[6] * data_[38];
  data_[111] = C00_[7] * data_[103] + B10_current[7] * data_[95] + cB00_current[7] * data_[39];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[112] = C00_[0] * data_[104] + B10_current[0] * data_[96] + cB00_current[0] * data_[40];
  data_[113] = C00_[1] * data_[105] + B10_current[1] * data_[97] + cB00_current[1] * data_[41];
  data_[114] = C00_[2] * data_[106] + B10_current[2] * data_[98] + cB00_current[2] * data_[42];
  data_[115] = C00_[3] * data_[107] + B10_current[3] * data_[99] + cB00_current[3] * data_[43];
  data_[116] = C00_[4] * data_[108] + B10_current[4] * data_[100] + cB00_current[4] * data_[44];
  data_[117] = C00_[5] * data_[109] + B10_current[5] * data_[101] + cB00_current[5] * data_[45];
  data_[118] = C00_[6] * data_[110] + B10_current[6] * data_[102] + cB00_current[6] * data_[46];
  data_[119] = C00_[7] * data_[111] + B10_current[7] * data_[103] + cB00_current[7] * data_[47];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[120] = C00_[0] * data_[112] + B10_current[0] * data_[104] + cB00_current[0] * data_[48];
  data_[121] = C00_[1] * data_[113] + B10_current[1] * data_[105] + cB00_current[1] * data_[49];
  data_[122] = C00_[2] * data_[114] + B10_current[2] * data_[106] + cB00_current[2] * data_[50];
  data_[123] = C00_[3] * data_[115] + B10_current[3] * data_[107] + cB00_current[3] * data_[51];
  data_[124] = C00_[4] * data_[116] + B10_current[4] * data_[108] + cB00_current[4] * data_[52];
  data_[125] = C00_[5] * data_[117] + B10_current[5] * data_[109] + cB00_current[5] * data_[53];
  data_[126] = C00_[6] * data_[118] + B10_current[6] * data_[110] + cB00_current[6] * data_[54];
  data_[127] = C00_[7] * data_[119] + B10_current[7] * data_[111] + cB00_current[7] * data_[55];

  double B01_current[8];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];
  B01_current[4] = B01_[4];
  B01_current[5] = B01_[5];
  B01_current[6] = B01_[6];
  B01_current[7] = B01_[7];

  data_[128] = D00_[0] * data_[64] + B01_current[0];
  data_[129] = D00_[1] * data_[65] + B01_current[1];
  data_[130] = D00_[2] * data_[66] + B01_current[2];
  data_[131] = D00_[3] * data_[67] + B01_current[3];
  data_[132] = D00_[4] * data_[68] + B01_current[4];
  data_[133] = D00_[5] * data_[69] + B01_current[5];
  data_[134] = D00_[6] * data_[70] + B01_current[6];
  data_[135] = D00_[7] * data_[71] + B01_current[7];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[136] = C00_[0] * data_[128] + cB00_current[0] * data_[64];
  data_[137] = C00_[1] * data_[129] + cB00_current[1] * data_[65];
  data_[138] = C00_[2] * data_[130] + cB00_current[2] * data_[66];
  data_[139] = C00_[3] * data_[131] + cB00_current[3] * data_[67];
  data_[140] = C00_[4] * data_[132] + cB00_current[4] * data_[68];
  data_[141] = C00_[5] * data_[133] + cB00_current[5] * data_[69];
  data_[142] = C00_[6] * data_[134] + cB00_current[6] * data_[70];
  data_[143] = C00_[7] * data_[135] + cB00_current[7] * data_[71];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[144] = C00_[0] * data_[136] + B10_current[0] * data_[128] + cB00_current[0] * data_[72];
  data_[145] = C00_[1] * data_[137] + B10_current[1] * data_[129] + cB00_current[1] * data_[73];
  data_[146] = C00_[2] * data_[138] + B10_current[2] * data_[130] + cB00_current[2] * data_[74];
  data_[147] = C00_[3] * data_[139] + B10_current[3] * data_[131] + cB00_current[3] * data_[75];
  data_[148] = C00_[4] * data_[140] + B10_current[4] * data_[132] + cB00_current[4] * data_[76];
  data_[149] = C00_[5] * data_[141] + B10_current[5] * data_[133] + cB00_current[5] * data_[77];
  data_[150] = C00_[6] * data_[142] + B10_current[6] * data_[134] + cB00_current[6] * data_[78];
  data_[151] = C00_[7] * data_[143] + B10_current[7] * data_[135] + cB00_current[7] * data_[79];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[152] = C00_[0] * data_[144] + B10_current[0] * data_[136] + cB00_current[0] * data_[80];
  data_[153] = C00_[1] * data_[145] + B10_current[1] * data_[137] + cB00_current[1] * data_[81];
  data_[154] = C00_[2] * data_[146] + B10_current[2] * data_[138] + cB00_current[2] * data_[82];
  data_[155] = C00_[3] * data_[147] + B10_current[3] * data_[139] + cB00_current[3] * data_[83];
  data_[156] = C00_[4] * data_[148] + B10_current[4] * data_[140] + cB00_current[4] * data_[84];
  data_[157] = C00_[5] * data_[149] + B10_current[5] * data_[141] + cB00_current[5] * data_[85];
  data_[158] = C00_[6] * data_[150] + B10_current[6] * data_[142] + cB00_current[6] * data_[86];
  data_[159] = C00_[7] * data_[151] + B10_current[7] * data_[143] + cB00_current[7] * data_[87];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[160] = C00_[0] * data_[152] + B10_current[0] * data_[144] + cB00_current[0] * data_[88];
  data_[161] = C00_[1] * data_[153] + B10_current[1] * data_[145] + cB00_current[1] * data_[89];
  data_[162] = C00_[2] * data_[154] + B10_current[2] * data_[146] + cB00_current[2] * data_[90];
  data_[163] = C00_[3] * data_[155] + B10_current[3] * data_[147] + cB00_current[3] * data_[91];
  data_[164] = C00_[4] * data_[156] + B10_current[4] * data_[148] + cB00_current[4] * data_[92];
  data_[165] = C00_[5] * data_[157] + B10_current[5] * data_[149] + cB00_current[5] * data_[93];
  data_[166] = C00_[6] * data_[158] + B10_current[6] * data_[150] + cB00_current[6] * data_[94];
  data_[167] = C00_[7] * data_[159] + B10_current[7] * data_[151] + cB00_current[7] * data_[95];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[168] = C00_[0] * data_[160] + B10_current[0] * data_[152] + cB00_current[0] * data_[96];
  data_[169] = C00_[1] * data_[161] + B10_current[1] * data_[153] + cB00_current[1] * data_[97];
  data_[170] = C00_[2] * data_[162] + B10_current[2] * data_[154] + cB00_current[2] * data_[98];
  data_[171] = C00_[3] * data_[163] + B10_current[3] * data_[155] + cB00_current[3] * data_[99];
  data_[172] = C00_[4] * data_[164] + B10_current[4] * data_[156] + cB00_current[4] * data_[100];
  data_[173] = C00_[5] * data_[165] + B10_current[5] * data_[157] + cB00_current[5] * data_[101];
  data_[174] = C00_[6] * data_[166] + B10_current[6] * data_[158] + cB00_current[6] * data_[102];
  data_[175] = C00_[7] * data_[167] + B10_current[7] * data_[159] + cB00_current[7] * data_[103];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[176] = C00_[0] * data_[168] + B10_current[0] * data_[160] + cB00_current[0] * data_[104];
  data_[177] = C00_[1] * data_[169] + B10_current[1] * data_[161] + cB00_current[1] * data_[105];
  data_[178] = C00_[2] * data_[170] + B10_current[2] * data_[162] + cB00_current[2] * data_[106];
  data_[179] = C00_[3] * data_[171] + B10_current[3] * data_[163] + cB00_current[3] * data_[107];
  data_[180] = C00_[4] * data_[172] + B10_current[4] * data_[164] + cB00_current[4] * data_[108];
  data_[181] = C00_[5] * data_[173] + B10_current[5] * data_[165] + cB00_current[5] * data_[109];
  data_[182] = C00_[6] * data_[174] + B10_current[6] * data_[166] + cB00_current[6] * data_[110];
  data_[183] = C00_[7] * data_[175] + B10_current[7] * data_[167] + cB00_current[7] * data_[111];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[184] = C00_[0] * data_[176] + B10_current[0] * data_[168] + cB00_current[0] * data_[112];
  data_[185] = C00_[1] * data_[177] + B10_current[1] * data_[169] + cB00_current[1] * data_[113];
  data_[186] = C00_[2] * data_[178] + B10_current[2] * data_[170] + cB00_current[2] * data_[114];
  data_[187] = C00_[3] * data_[179] + B10_current[3] * data_[171] + cB00_current[3] * data_[115];
  data_[188] = C00_[4] * data_[180] + B10_current[4] * data_[172] + cB00_current[4] * data_[116];
  data_[189] = C00_[5] * data_[181] + B10_current[5] * data_[173] + cB00_current[5] * data_[117];
  data_[190] = C00_[6] * data_[182] + B10_current[6] * data_[174] + cB00_current[6] * data_[118];
  data_[191] = C00_[7] * data_[183] + B10_current[7] * data_[175] + cB00_current[7] * data_[119];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[192] = D00_[0] * data_[128] + B01_current[0] * data_[64];
  data_[193] = D00_[1] * data_[129] + B01_current[1] * data_[65];
  data_[194] = D00_[2] * data_[130] + B01_current[2] * data_[66];
  data_[195] = D00_[3] * data_[131] + B01_current[3] * data_[67];
  data_[196] = D00_[4] * data_[132] + B01_current[4] * data_[68];
  data_[197] = D00_[5] * data_[133] + B01_current[5] * data_[69];
  data_[198] = D00_[6] * data_[134] + B01_current[6] * data_[70];
  data_[199] = D00_[7] * data_[135] + B01_current[7] * data_[71];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[200] = C00_[0] * data_[192] + cB00_current[0] * data_[128];
  data_[201] = C00_[1] * data_[193] + cB00_current[1] * data_[129];
  data_[202] = C00_[2] * data_[194] + cB00_current[2] * data_[130];
  data_[203] = C00_[3] * data_[195] + cB00_current[3] * data_[131];
  data_[204] = C00_[4] * data_[196] + cB00_current[4] * data_[132];
  data_[205] = C00_[5] * data_[197] + cB00_current[5] * data_[133];
  data_[206] = C00_[6] * data_[198] + cB00_current[6] * data_[134];
  data_[207] = C00_[7] * data_[199] + cB00_current[7] * data_[135];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[208] = C00_[0] * data_[200] + B10_current[0] * data_[192] + cB00_current[0] * data_[136];
  data_[209] = C00_[1] * data_[201] + B10_current[1] * data_[193] + cB00_current[1] * data_[137];
  data_[210] = C00_[2] * data_[202] + B10_current[2] * data_[194] + cB00_current[2] * data_[138];
  data_[211] = C00_[3] * data_[203] + B10_current[3] * data_[195] + cB00_current[3] * data_[139];
  data_[212] = C00_[4] * data_[204] + B10_current[4] * data_[196] + cB00_current[4] * data_[140];
  data_[213] = C00_[5] * data_[205] + B10_current[5] * data_[197] + cB00_current[5] * data_[141];
  data_[214] = C00_[6] * data_[206] + B10_current[6] * data_[198] + cB00_current[6] * data_[142];
  data_[215] = C00_[7] * data_[207] + B10_current[7] * data_[199] + cB00_current[7] * data_[143];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[216] = C00_[0] * data_[208] + B10_current[0] * data_[200] + cB00_current[0] * data_[144];
  data_[217] = C00_[1] * data_[209] + B10_current[1] * data_[201] + cB00_current[1] * data_[145];
  data_[218] = C00_[2] * data_[210] + B10_current[2] * data_[202] + cB00_current[2] * data_[146];
  data_[219] = C00_[3] * data_[211] + B10_current[3] * data_[203] + cB00_current[3] * data_[147];
  data_[220] = C00_[4] * data_[212] + B10_current[4] * data_[204] + cB00_current[4] * data_[148];
  data_[221] = C00_[5] * data_[213] + B10_current[5] * data_[205] + cB00_current[5] * data_[149];
  data_[222] = C00_[6] * data_[214] + B10_current[6] * data_[206] + cB00_current[6] * data_[150];
  data_[223] = C00_[7] * data_[215] + B10_current[7] * data_[207] + cB00_current[7] * data_[151];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[224] = C00_[0] * data_[216] + B10_current[0] * data_[208] + cB00_current[0] * data_[152];
  data_[225] = C00_[1] * data_[217] + B10_current[1] * data_[209] + cB00_current[1] * data_[153];
  data_[226] = C00_[2] * data_[218] + B10_current[2] * data_[210] + cB00_current[2] * data_[154];
  data_[227] = C00_[3] * data_[219] + B10_current[3] * data_[211] + cB00_current[3] * data_[155];
  data_[228] = C00_[4] * data_[220] + B10_current[4] * data_[212] + cB00_current[4] * data_[156];
  data_[229] = C00_[5] * data_[221] + B10_current[5] * data_[213] + cB00_current[5] * data_[157];
  data_[230] = C00_[6] * data_[222] + B10_current[6] * data_[214] + cB00_current[6] * data_[158];
  data_[231] = C00_[7] * data_[223] + B10_current[7] * data_[215] + cB00_current[7] * data_[159];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[232] = C00_[0] * data_[224] + B10_current[0] * data_[216] + cB00_current[0] * data_[160];
  data_[233] = C00_[1] * data_[225] + B10_current[1] * data_[217] + cB00_current[1] * data_[161];
  data_[234] = C00_[2] * data_[226] + B10_current[2] * data_[218] + cB00_current[2] * data_[162];
  data_[235] = C00_[3] * data_[227] + B10_current[3] * data_[219] + cB00_current[3] * data_[163];
  data_[236] = C00_[4] * data_[228] + B10_current[4] * data_[220] + cB00_current[4] * data_[164];
  data_[237] = C00_[5] * data_[229] + B10_current[5] * data_[221] + cB00_current[5] * data_[165];
  data_[238] = C00_[6] * data_[230] + B10_current[6] * data_[222] + cB00_current[6] * data_[166];
  data_[239] = C00_[7] * data_[231] + B10_current[7] * data_[223] + cB00_current[7] * data_[167];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[240] = C00_[0] * data_[232] + B10_current[0] * data_[224] + cB00_current[0] * data_[168];
  data_[241] = C00_[1] * data_[233] + B10_current[1] * data_[225] + cB00_current[1] * data_[169];
  data_[242] = C00_[2] * data_[234] + B10_current[2] * data_[226] + cB00_current[2] * data_[170];
  data_[243] = C00_[3] * data_[235] + B10_current[3] * data_[227] + cB00_current[3] * data_[171];
  data_[244] = C00_[4] * data_[236] + B10_current[4] * data_[228] + cB00_current[4] * data_[172];
  data_[245] = C00_[5] * data_[237] + B10_current[5] * data_[229] + cB00_current[5] * data_[173];
  data_[246] = C00_[6] * data_[238] + B10_current[6] * data_[230] + cB00_current[6] * data_[174];
  data_[247] = C00_[7] * data_[239] + B10_current[7] * data_[231] + cB00_current[7] * data_[175];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[248] = C00_[0] * data_[240] + B10_current[0] * data_[232] + cB00_current[0] * data_[176];
  data_[249] = C00_[1] * data_[241] + B10_current[1] * data_[233] + cB00_current[1] * data_[177];
  data_[250] = C00_[2] * data_[242] + B10_current[2] * data_[234] + cB00_current[2] * data_[178];
  data_[251] = C00_[3] * data_[243] + B10_current[3] * data_[235] + cB00_current[3] * data_[179];
  data_[252] = C00_[4] * data_[244] + B10_current[4] * data_[236] + cB00_current[4] * data_[180];
  data_[253] = C00_[5] * data_[245] + B10_current[5] * data_[237] + cB00_current[5] * data_[181];
  data_[254] = C00_[6] * data_[246] + B10_current[6] * data_[238] + cB00_current[6] * data_[182];
  data_[255] = C00_[7] * data_[247] + B10_current[7] * data_[239] + cB00_current[7] * data_[183];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[256] = D00_[0] * data_[192] + B01_current[0] * data_[128];
  data_[257] = D00_[1] * data_[193] + B01_current[1] * data_[129];
  data_[258] = D00_[2] * data_[194] + B01_current[2] * data_[130];
  data_[259] = D00_[3] * data_[195] + B01_current[3] * data_[131];
  data_[260] = D00_[4] * data_[196] + B01_current[4] * data_[132];
  data_[261] = D00_[5] * data_[197] + B01_current[5] * data_[133];
  data_[262] = D00_[6] * data_[198] + B01_current[6] * data_[134];
  data_[263] = D00_[7] * data_[199] + B01_current[7] * data_[135];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[264] = C00_[0] * data_[256] + cB00_current[0] * data_[192];
  data_[265] = C00_[1] * data_[257] + cB00_current[1] * data_[193];
  data_[266] = C00_[2] * data_[258] + cB00_current[2] * data_[194];
  data_[267] = C00_[3] * data_[259] + cB00_current[3] * data_[195];
  data_[268] = C00_[4] * data_[260] + cB00_current[4] * data_[196];
  data_[269] = C00_[5] * data_[261] + cB00_current[5] * data_[197];
  data_[270] = C00_[6] * data_[262] + cB00_current[6] * data_[198];
  data_[271] = C00_[7] * data_[263] + cB00_current[7] * data_[199];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[272] = C00_[0] * data_[264] + B10_current[0] * data_[256] + cB00_current[0] * data_[200];
  data_[273] = C00_[1] * data_[265] + B10_current[1] * data_[257] + cB00_current[1] * data_[201];
  data_[274] = C00_[2] * data_[266] + B10_current[2] * data_[258] + cB00_current[2] * data_[202];
  data_[275] = C00_[3] * data_[267] + B10_current[3] * data_[259] + cB00_current[3] * data_[203];
  data_[276] = C00_[4] * data_[268] + B10_current[4] * data_[260] + cB00_current[4] * data_[204];
  data_[277] = C00_[5] * data_[269] + B10_current[5] * data_[261] + cB00_current[5] * data_[205];
  data_[278] = C00_[6] * data_[270] + B10_current[6] * data_[262] + cB00_current[6] * data_[206];
  data_[279] = C00_[7] * data_[271] + B10_current[7] * data_[263] + cB00_current[7] * data_[207];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[280] = C00_[0] * data_[272] + B10_current[0] * data_[264] + cB00_current[0] * data_[208];
  data_[281] = C00_[1] * data_[273] + B10_current[1] * data_[265] + cB00_current[1] * data_[209];
  data_[282] = C00_[2] * data_[274] + B10_current[2] * data_[266] + cB00_current[2] * data_[210];
  data_[283] = C00_[3] * data_[275] + B10_current[3] * data_[267] + cB00_current[3] * data_[211];
  data_[284] = C00_[4] * data_[276] + B10_current[4] * data_[268] + cB00_current[4] * data_[212];
  data_[285] = C00_[5] * data_[277] + B10_current[5] * data_[269] + cB00_current[5] * data_[213];
  data_[286] = C00_[6] * data_[278] + B10_current[6] * data_[270] + cB00_current[6] * data_[214];
  data_[287] = C00_[7] * data_[279] + B10_current[7] * data_[271] + cB00_current[7] * data_[215];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[288] = C00_[0] * data_[280] + B10_current[0] * data_[272] + cB00_current[0] * data_[216];
  data_[289] = C00_[1] * data_[281] + B10_current[1] * data_[273] + cB00_current[1] * data_[217];
  data_[290] = C00_[2] * data_[282] + B10_current[2] * data_[274] + cB00_current[2] * data_[218];
  data_[291] = C00_[3] * data_[283] + B10_current[3] * data_[275] + cB00_current[3] * data_[219];
  data_[292] = C00_[4] * data_[284] + B10_current[4] * data_[276] + cB00_current[4] * data_[220];
  data_[293] = C00_[5] * data_[285] + B10_current[5] * data_[277] + cB00_current[5] * data_[221];
  data_[294] = C00_[6] * data_[286] + B10_current[6] * data_[278] + cB00_current[6] * data_[222];
  data_[295] = C00_[7] * data_[287] + B10_current[7] * data_[279] + cB00_current[7] * data_[223];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[296] = C00_[0] * data_[288] + B10_current[0] * data_[280] + cB00_current[0] * data_[224];
  data_[297] = C00_[1] * data_[289] + B10_current[1] * data_[281] + cB00_current[1] * data_[225];
  data_[298] = C00_[2] * data_[290] + B10_current[2] * data_[282] + cB00_current[2] * data_[226];
  data_[299] = C00_[3] * data_[291] + B10_current[3] * data_[283] + cB00_current[3] * data_[227];
  data_[300] = C00_[4] * data_[292] + B10_current[4] * data_[284] + cB00_current[4] * data_[228];
  data_[301] = C00_[5] * data_[293] + B10_current[5] * data_[285] + cB00_current[5] * data_[229];
  data_[302] = C00_[6] * data_[294] + B10_current[6] * data_[286] + cB00_current[6] * data_[230];
  data_[303] = C00_[7] * data_[295] + B10_current[7] * data_[287] + cB00_current[7] * data_[231];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[304] = C00_[0] * data_[296] + B10_current[0] * data_[288] + cB00_current[0] * data_[232];
  data_[305] = C00_[1] * data_[297] + B10_current[1] * data_[289] + cB00_current[1] * data_[233];
  data_[306] = C00_[2] * data_[298] + B10_current[2] * data_[290] + cB00_current[2] * data_[234];
  data_[307] = C00_[3] * data_[299] + B10_current[3] * data_[291] + cB00_current[3] * data_[235];
  data_[308] = C00_[4] * data_[300] + B10_current[4] * data_[292] + cB00_current[4] * data_[236];
  data_[309] = C00_[5] * data_[301] + B10_current[5] * data_[293] + cB00_current[5] * data_[237];
  data_[310] = C00_[6] * data_[302] + B10_current[6] * data_[294] + cB00_current[6] * data_[238];
  data_[311] = C00_[7] * data_[303] + B10_current[7] * data_[295] + cB00_current[7] * data_[239];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[312] = C00_[0] * data_[304] + B10_current[0] * data_[296] + cB00_current[0] * data_[240];
  data_[313] = C00_[1] * data_[305] + B10_current[1] * data_[297] + cB00_current[1] * data_[241];
  data_[314] = C00_[2] * data_[306] + B10_current[2] * data_[298] + cB00_current[2] * data_[242];
  data_[315] = C00_[3] * data_[307] + B10_current[3] * data_[299] + cB00_current[3] * data_[243];
  data_[316] = C00_[4] * data_[308] + B10_current[4] * data_[300] + cB00_current[4] * data_[244];
  data_[317] = C00_[5] * data_[309] + B10_current[5] * data_[301] + cB00_current[5] * data_[245];
  data_[318] = C00_[6] * data_[310] + B10_current[6] * data_[302] + cB00_current[6] * data_[246];
  data_[319] = C00_[7] * data_[311] + B10_current[7] * data_[303] + cB00_current[7] * data_[247];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[320] = D00_[0] * data_[256] + B01_current[0] * data_[192];
  data_[321] = D00_[1] * data_[257] + B01_current[1] * data_[193];
  data_[322] = D00_[2] * data_[258] + B01_current[2] * data_[194];
  data_[323] = D00_[3] * data_[259] + B01_current[3] * data_[195];
  data_[324] = D00_[4] * data_[260] + B01_current[4] * data_[196];
  data_[325] = D00_[5] * data_[261] + B01_current[5] * data_[197];
  data_[326] = D00_[6] * data_[262] + B01_current[6] * data_[198];
  data_[327] = D00_[7] * data_[263] + B01_current[7] * data_[199];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[328] = C00_[0] * data_[320] + cB00_current[0] * data_[256];
  data_[329] = C00_[1] * data_[321] + cB00_current[1] * data_[257];
  data_[330] = C00_[2] * data_[322] + cB00_current[2] * data_[258];
  data_[331] = C00_[3] * data_[323] + cB00_current[3] * data_[259];
  data_[332] = C00_[4] * data_[324] + cB00_current[4] * data_[260];
  data_[333] = C00_[5] * data_[325] + cB00_current[5] * data_[261];
  data_[334] = C00_[6] * data_[326] + cB00_current[6] * data_[262];
  data_[335] = C00_[7] * data_[327] + cB00_current[7] * data_[263];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[336] = C00_[0] * data_[328] + B10_current[0] * data_[320] + cB00_current[0] * data_[264];
  data_[337] = C00_[1] * data_[329] + B10_current[1] * data_[321] + cB00_current[1] * data_[265];
  data_[338] = C00_[2] * data_[330] + B10_current[2] * data_[322] + cB00_current[2] * data_[266];
  data_[339] = C00_[3] * data_[331] + B10_current[3] * data_[323] + cB00_current[3] * data_[267];
  data_[340] = C00_[4] * data_[332] + B10_current[4] * data_[324] + cB00_current[4] * data_[268];
  data_[341] = C00_[5] * data_[333] + B10_current[5] * data_[325] + cB00_current[5] * data_[269];
  data_[342] = C00_[6] * data_[334] + B10_current[6] * data_[326] + cB00_current[6] * data_[270];
  data_[343] = C00_[7] * data_[335] + B10_current[7] * data_[327] + cB00_current[7] * data_[271];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[344] = C00_[0] * data_[336] + B10_current[0] * data_[328] + cB00_current[0] * data_[272];
  data_[345] = C00_[1] * data_[337] + B10_current[1] * data_[329] + cB00_current[1] * data_[273];
  data_[346] = C00_[2] * data_[338] + B10_current[2] * data_[330] + cB00_current[2] * data_[274];
  data_[347] = C00_[3] * data_[339] + B10_current[3] * data_[331] + cB00_current[3] * data_[275];
  data_[348] = C00_[4] * data_[340] + B10_current[4] * data_[332] + cB00_current[4] * data_[276];
  data_[349] = C00_[5] * data_[341] + B10_current[5] * data_[333] + cB00_current[5] * data_[277];
  data_[350] = C00_[6] * data_[342] + B10_current[6] * data_[334] + cB00_current[6] * data_[278];
  data_[351] = C00_[7] * data_[343] + B10_current[7] * data_[335] + cB00_current[7] * data_[279];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[352] = C00_[0] * data_[344] + B10_current[0] * data_[336] + cB00_current[0] * data_[280];
  data_[353] = C00_[1] * data_[345] + B10_current[1] * data_[337] + cB00_current[1] * data_[281];
  data_[354] = C00_[2] * data_[346] + B10_current[2] * data_[338] + cB00_current[2] * data_[282];
  data_[355] = C00_[3] * data_[347] + B10_current[3] * data_[339] + cB00_current[3] * data_[283];
  data_[356] = C00_[4] * data_[348] + B10_current[4] * data_[340] + cB00_current[4] * data_[284];
  data_[357] = C00_[5] * data_[349] + B10_current[5] * data_[341] + cB00_current[5] * data_[285];
  data_[358] = C00_[6] * data_[350] + B10_current[6] * data_[342] + cB00_current[6] * data_[286];
  data_[359] = C00_[7] * data_[351] + B10_current[7] * data_[343] + cB00_current[7] * data_[287];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[360] = C00_[0] * data_[352] + B10_current[0] * data_[344] + cB00_current[0] * data_[288];
  data_[361] = C00_[1] * data_[353] + B10_current[1] * data_[345] + cB00_current[1] * data_[289];
  data_[362] = C00_[2] * data_[354] + B10_current[2] * data_[346] + cB00_current[2] * data_[290];
  data_[363] = C00_[3] * data_[355] + B10_current[3] * data_[347] + cB00_current[3] * data_[291];
  data_[364] = C00_[4] * data_[356] + B10_current[4] * data_[348] + cB00_current[4] * data_[292];
  data_[365] = C00_[5] * data_[357] + B10_current[5] * data_[349] + cB00_current[5] * data_[293];
  data_[366] = C00_[6] * data_[358] + B10_current[6] * data_[350] + cB00_current[6] * data_[294];
  data_[367] = C00_[7] * data_[359] + B10_current[7] * data_[351] + cB00_current[7] * data_[295];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[368] = C00_[0] * data_[360] + B10_current[0] * data_[352] + cB00_current[0] * data_[296];
  data_[369] = C00_[1] * data_[361] + B10_current[1] * data_[353] + cB00_current[1] * data_[297];
  data_[370] = C00_[2] * data_[362] + B10_current[2] * data_[354] + cB00_current[2] * data_[298];
  data_[371] = C00_[3] * data_[363] + B10_current[3] * data_[355] + cB00_current[3] * data_[299];
  data_[372] = C00_[4] * data_[364] + B10_current[4] * data_[356] + cB00_current[4] * data_[300];
  data_[373] = C00_[5] * data_[365] + B10_current[5] * data_[357] + cB00_current[5] * data_[301];
  data_[374] = C00_[6] * data_[366] + B10_current[6] * data_[358] + cB00_current[6] * data_[302];
  data_[375] = C00_[7] * data_[367] + B10_current[7] * data_[359] + cB00_current[7] * data_[303];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[376] = C00_[0] * data_[368] + B10_current[0] * data_[360] + cB00_current[0] * data_[304];
  data_[377] = C00_[1] * data_[369] + B10_current[1] * data_[361] + cB00_current[1] * data_[305];
  data_[378] = C00_[2] * data_[370] + B10_current[2] * data_[362] + cB00_current[2] * data_[306];
  data_[379] = C00_[3] * data_[371] + B10_current[3] * data_[363] + cB00_current[3] * data_[307];
  data_[380] = C00_[4] * data_[372] + B10_current[4] * data_[364] + cB00_current[4] * data_[308];
  data_[381] = C00_[5] * data_[373] + B10_current[5] * data_[365] + cB00_current[5] * data_[309];
  data_[382] = C00_[6] * data_[374] + B10_current[6] * data_[366] + cB00_current[6] * data_[310];
  data_[383] = C00_[7] * data_[375] + B10_current[7] * data_[367] + cB00_current[7] * data_[311];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[384] = D00_[0] * data_[320] + B01_current[0] * data_[256];
  data_[385] = D00_[1] * data_[321] + B01_current[1] * data_[257];
  data_[386] = D00_[2] * data_[322] + B01_current[2] * data_[258];
  data_[387] = D00_[3] * data_[323] + B01_current[3] * data_[259];
  data_[388] = D00_[4] * data_[324] + B01_current[4] * data_[260];
  data_[389] = D00_[5] * data_[325] + B01_current[5] * data_[261];
  data_[390] = D00_[6] * data_[326] + B01_current[6] * data_[262];
  data_[391] = D00_[7] * data_[327] + B01_current[7] * data_[263];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[392] = C00_[0] * data_[384] + cB00_current[0] * data_[320];
  data_[393] = C00_[1] * data_[385] + cB00_current[1] * data_[321];
  data_[394] = C00_[2] * data_[386] + cB00_current[2] * data_[322];
  data_[395] = C00_[3] * data_[387] + cB00_current[3] * data_[323];
  data_[396] = C00_[4] * data_[388] + cB00_current[4] * data_[324];
  data_[397] = C00_[5] * data_[389] + cB00_current[5] * data_[325];
  data_[398] = C00_[6] * data_[390] + cB00_current[6] * data_[326];
  data_[399] = C00_[7] * data_[391] + cB00_current[7] * data_[327];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[400] = C00_[0] * data_[392] + B10_current[0] * data_[384] + cB00_current[0] * data_[328];
  data_[401] = C00_[1] * data_[393] + B10_current[1] * data_[385] + cB00_current[1] * data_[329];
  data_[402] = C00_[2] * data_[394] + B10_current[2] * data_[386] + cB00_current[2] * data_[330];
  data_[403] = C00_[3] * data_[395] + B10_current[3] * data_[387] + cB00_current[3] * data_[331];
  data_[404] = C00_[4] * data_[396] + B10_current[4] * data_[388] + cB00_current[4] * data_[332];
  data_[405] = C00_[5] * data_[397] + B10_current[5] * data_[389] + cB00_current[5] * data_[333];
  data_[406] = C00_[6] * data_[398] + B10_current[6] * data_[390] + cB00_current[6] * data_[334];
  data_[407] = C00_[7] * data_[399] + B10_current[7] * data_[391] + cB00_current[7] * data_[335];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[408] = C00_[0] * data_[400] + B10_current[0] * data_[392] + cB00_current[0] * data_[336];
  data_[409] = C00_[1] * data_[401] + B10_current[1] * data_[393] + cB00_current[1] * data_[337];
  data_[410] = C00_[2] * data_[402] + B10_current[2] * data_[394] + cB00_current[2] * data_[338];
  data_[411] = C00_[3] * data_[403] + B10_current[3] * data_[395] + cB00_current[3] * data_[339];
  data_[412] = C00_[4] * data_[404] + B10_current[4] * data_[396] + cB00_current[4] * data_[340];
  data_[413] = C00_[5] * data_[405] + B10_current[5] * data_[397] + cB00_current[5] * data_[341];
  data_[414] = C00_[6] * data_[406] + B10_current[6] * data_[398] + cB00_current[6] * data_[342];
  data_[415] = C00_[7] * data_[407] + B10_current[7] * data_[399] + cB00_current[7] * data_[343];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[416] = C00_[0] * data_[408] + B10_current[0] * data_[400] + cB00_current[0] * data_[344];
  data_[417] = C00_[1] * data_[409] + B10_current[1] * data_[401] + cB00_current[1] * data_[345];
  data_[418] = C00_[2] * data_[410] + B10_current[2] * data_[402] + cB00_current[2] * data_[346];
  data_[419] = C00_[3] * data_[411] + B10_current[3] * data_[403] + cB00_current[3] * data_[347];
  data_[420] = C00_[4] * data_[412] + B10_current[4] * data_[404] + cB00_current[4] * data_[348];
  data_[421] = C00_[5] * data_[413] + B10_current[5] * data_[405] + cB00_current[5] * data_[349];
  data_[422] = C00_[6] * data_[414] + B10_current[6] * data_[406] + cB00_current[6] * data_[350];
  data_[423] = C00_[7] * data_[415] + B10_current[7] * data_[407] + cB00_current[7] * data_[351];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[424] = C00_[0] * data_[416] + B10_current[0] * data_[408] + cB00_current[0] * data_[352];
  data_[425] = C00_[1] * data_[417] + B10_current[1] * data_[409] + cB00_current[1] * data_[353];
  data_[426] = C00_[2] * data_[418] + B10_current[2] * data_[410] + cB00_current[2] * data_[354];
  data_[427] = C00_[3] * data_[419] + B10_current[3] * data_[411] + cB00_current[3] * data_[355];
  data_[428] = C00_[4] * data_[420] + B10_current[4] * data_[412] + cB00_current[4] * data_[356];
  data_[429] = C00_[5] * data_[421] + B10_current[5] * data_[413] + cB00_current[5] * data_[357];
  data_[430] = C00_[6] * data_[422] + B10_current[6] * data_[414] + cB00_current[6] * data_[358];
  data_[431] = C00_[7] * data_[423] + B10_current[7] * data_[415] + cB00_current[7] * data_[359];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[432] = C00_[0] * data_[424] + B10_current[0] * data_[416] + cB00_current[0] * data_[360];
  data_[433] = C00_[1] * data_[425] + B10_current[1] * data_[417] + cB00_current[1] * data_[361];
  data_[434] = C00_[2] * data_[426] + B10_current[2] * data_[418] + cB00_current[2] * data_[362];
  data_[435] = C00_[3] * data_[427] + B10_current[3] * data_[419] + cB00_current[3] * data_[363];
  data_[436] = C00_[4] * data_[428] + B10_current[4] * data_[420] + cB00_current[4] * data_[364];
  data_[437] = C00_[5] * data_[429] + B10_current[5] * data_[421] + cB00_current[5] * data_[365];
  data_[438] = C00_[6] * data_[430] + B10_current[6] * data_[422] + cB00_current[6] * data_[366];
  data_[439] = C00_[7] * data_[431] + B10_current[7] * data_[423] + cB00_current[7] * data_[367];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[440] = C00_[0] * data_[432] + B10_current[0] * data_[424] + cB00_current[0] * data_[368];
  data_[441] = C00_[1] * data_[433] + B10_current[1] * data_[425] + cB00_current[1] * data_[369];
  data_[442] = C00_[2] * data_[434] + B10_current[2] * data_[426] + cB00_current[2] * data_[370];
  data_[443] = C00_[3] * data_[435] + B10_current[3] * data_[427] + cB00_current[3] * data_[371];
  data_[444] = C00_[4] * data_[436] + B10_current[4] * data_[428] + cB00_current[4] * data_[372];
  data_[445] = C00_[5] * data_[437] + B10_current[5] * data_[429] + cB00_current[5] * data_[373];
  data_[446] = C00_[6] * data_[438] + B10_current[6] * data_[430] + cB00_current[6] * data_[374];
  data_[447] = C00_[7] * data_[439] + B10_current[7] * data_[431] + cB00_current[7] * data_[375];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];
  B01_current[6] += B01_[6];
  B01_current[7] += B01_[7];

  data_[448] = D00_[0] * data_[384] + B01_current[0] * data_[320];
  data_[449] = D00_[1] * data_[385] + B01_current[1] * data_[321];
  data_[450] = D00_[2] * data_[386] + B01_current[2] * data_[322];
  data_[451] = D00_[3] * data_[387] + B01_current[3] * data_[323];
  data_[452] = D00_[4] * data_[388] + B01_current[4] * data_[324];
  data_[453] = D00_[5] * data_[389] + B01_current[5] * data_[325];
  data_[454] = D00_[6] * data_[390] + B01_current[6] * data_[326];
  data_[455] = D00_[7] * data_[391] + B01_current[7] * data_[327];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];
  cB00_current[6] += B00_[6];
  cB00_current[7] += B00_[7];

  data_[456] = C00_[0] * data_[448] + cB00_current[0] * data_[384];
  data_[457] = C00_[1] * data_[449] + cB00_current[1] * data_[385];
  data_[458] = C00_[2] * data_[450] + cB00_current[2] * data_[386];
  data_[459] = C00_[3] * data_[451] + cB00_current[3] * data_[387];
  data_[460] = C00_[4] * data_[452] + cB00_current[4] * data_[388];
  data_[461] = C00_[5] * data_[453] + cB00_current[5] * data_[389];
  data_[462] = C00_[6] * data_[454] + cB00_current[6] * data_[390];
  data_[463] = C00_[7] * data_[455] + cB00_current[7] * data_[391];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];
  B10_current[6] = B10_[6];
  B10_current[7] = B10_[7];

  data_[464] = C00_[0] * data_[456] + B10_current[0] * data_[448] + cB00_current[0] * data_[392];
  data_[465] = C00_[1] * data_[457] + B10_current[1] * data_[449] + cB00_current[1] * data_[393];
  data_[466] = C00_[2] * data_[458] + B10_current[2] * data_[450] + cB00_current[2] * data_[394];
  data_[467] = C00_[3] * data_[459] + B10_current[3] * data_[451] + cB00_current[3] * data_[395];
  data_[468] = C00_[4] * data_[460] + B10_current[4] * data_[452] + cB00_current[4] * data_[396];
  data_[469] = C00_[5] * data_[461] + B10_current[5] * data_[453] + cB00_current[5] * data_[397];
  data_[470] = C00_[6] * data_[462] + B10_current[6] * data_[454] + cB00_current[6] * data_[398];
  data_[471] = C00_[7] * data_[463] + B10_current[7] * data_[455] + cB00_current[7] * data_[399];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[472] = C00_[0] * data_[464] + B10_current[0] * data_[456] + cB00_current[0] * data_[400];
  data_[473] = C00_[1] * data_[465] + B10_current[1] * data_[457] + cB00_current[1] * data_[401];
  data_[474] = C00_[2] * data_[466] + B10_current[2] * data_[458] + cB00_current[2] * data_[402];
  data_[475] = C00_[3] * data_[467] + B10_current[3] * data_[459] + cB00_current[3] * data_[403];
  data_[476] = C00_[4] * data_[468] + B10_current[4] * data_[460] + cB00_current[4] * data_[404];
  data_[477] = C00_[5] * data_[469] + B10_current[5] * data_[461] + cB00_current[5] * data_[405];
  data_[478] = C00_[6] * data_[470] + B10_current[6] * data_[462] + cB00_current[6] * data_[406];
  data_[479] = C00_[7] * data_[471] + B10_current[7] * data_[463] + cB00_current[7] * data_[407];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[480] = C00_[0] * data_[472] + B10_current[0] * data_[464] + cB00_current[0] * data_[408];
  data_[481] = C00_[1] * data_[473] + B10_current[1] * data_[465] + cB00_current[1] * data_[409];
  data_[482] = C00_[2] * data_[474] + B10_current[2] * data_[466] + cB00_current[2] * data_[410];
  data_[483] = C00_[3] * data_[475] + B10_current[3] * data_[467] + cB00_current[3] * data_[411];
  data_[484] = C00_[4] * data_[476] + B10_current[4] * data_[468] + cB00_current[4] * data_[412];
  data_[485] = C00_[5] * data_[477] + B10_current[5] * data_[469] + cB00_current[5] * data_[413];
  data_[486] = C00_[6] * data_[478] + B10_current[6] * data_[470] + cB00_current[6] * data_[414];
  data_[487] = C00_[7] * data_[479] + B10_current[7] * data_[471] + cB00_current[7] * data_[415];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[488] = C00_[0] * data_[480] + B10_current[0] * data_[472] + cB00_current[0] * data_[416];
  data_[489] = C00_[1] * data_[481] + B10_current[1] * data_[473] + cB00_current[1] * data_[417];
  data_[490] = C00_[2] * data_[482] + B10_current[2] * data_[474] + cB00_current[2] * data_[418];
  data_[491] = C00_[3] * data_[483] + B10_current[3] * data_[475] + cB00_current[3] * data_[419];
  data_[492] = C00_[4] * data_[484] + B10_current[4] * data_[476] + cB00_current[4] * data_[420];
  data_[493] = C00_[5] * data_[485] + B10_current[5] * data_[477] + cB00_current[5] * data_[421];
  data_[494] = C00_[6] * data_[486] + B10_current[6] * data_[478] + cB00_current[6] * data_[422];
  data_[495] = C00_[7] * data_[487] + B10_current[7] * data_[479] + cB00_current[7] * data_[423];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[496] = C00_[0] * data_[488] + B10_current[0] * data_[480] + cB00_current[0] * data_[424];
  data_[497] = C00_[1] * data_[489] + B10_current[1] * data_[481] + cB00_current[1] * data_[425];
  data_[498] = C00_[2] * data_[490] + B10_current[2] * data_[482] + cB00_current[2] * data_[426];
  data_[499] = C00_[3] * data_[491] + B10_current[3] * data_[483] + cB00_current[3] * data_[427];
  data_[500] = C00_[4] * data_[492] + B10_current[4] * data_[484] + cB00_current[4] * data_[428];
  data_[501] = C00_[5] * data_[493] + B10_current[5] * data_[485] + cB00_current[5] * data_[429];
  data_[502] = C00_[6] * data_[494] + B10_current[6] * data_[486] + cB00_current[6] * data_[430];
  data_[503] = C00_[7] * data_[495] + B10_current[7] * data_[487] + cB00_current[7] * data_[431];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];
  B10_current[6] += B10_[6];
  B10_current[7] += B10_[7];

  data_[504] = C00_[0] * data_[496] + B10_current[0] * data_[488] + cB00_current[0] * data_[432];
  data_[505] = C00_[1] * data_[497] + B10_current[1] * data_[489] + cB00_current[1] * data_[433];
  data_[506] = C00_[2] * data_[498] + B10_current[2] * data_[490] + cB00_current[2] * data_[434];
  data_[507] = C00_[3] * data_[499] + B10_current[3] * data_[491] + cB00_current[3] * data_[435];
  data_[508] = C00_[4] * data_[500] + B10_current[4] * data_[492] + cB00_current[4] * data_[436];
  data_[509] = C00_[5] * data_[501] + B10_current[5] * data_[493] + cB00_current[5] * data_[437];
  data_[510] = C00_[6] * data_[502] + B10_current[6] * data_[494] + cB00_current[6] * data_[438];
  data_[511] = C00_[7] * data_[503] + B10_current[7] * data_[495] + cB00_current[7] * data_[439];
}

