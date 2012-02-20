//
// Newint - Parallel electron correlation program.
// Filename: _svrr_6040.cc
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


#include <src/slater/svrrlist.h>

// returns double array of length 210
void SVRRList::_svrr_6040(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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

  double B10_current[6];
  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];

  data_[12] = C00_[0] * data_[6] + B10_current[0];
  data_[13] = C00_[1] * data_[7] + B10_current[1];
  data_[14] = C00_[2] * data_[8] + B10_current[2];
  data_[15] = C00_[3] * data_[9] + B10_current[3];
  data_[16] = C00_[4] * data_[10] + B10_current[4];
  data_[17] = C00_[5] * data_[11] + B10_current[5];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[18] = C00_[0] * data_[12] + B10_current[0] * data_[6];
  data_[19] = C00_[1] * data_[13] + B10_current[1] * data_[7];
  data_[20] = C00_[2] * data_[14] + B10_current[2] * data_[8];
  data_[21] = C00_[3] * data_[15] + B10_current[3] * data_[9];
  data_[22] = C00_[4] * data_[16] + B10_current[4] * data_[10];
  data_[23] = C00_[5] * data_[17] + B10_current[5] * data_[11];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[24] = C00_[0] * data_[18] + B10_current[0] * data_[12];
  data_[25] = C00_[1] * data_[19] + B10_current[1] * data_[13];
  data_[26] = C00_[2] * data_[20] + B10_current[2] * data_[14];
  data_[27] = C00_[3] * data_[21] + B10_current[3] * data_[15];
  data_[28] = C00_[4] * data_[22] + B10_current[4] * data_[16];
  data_[29] = C00_[5] * data_[23] + B10_current[5] * data_[17];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[30] = C00_[0] * data_[24] + B10_current[0] * data_[18];
  data_[31] = C00_[1] * data_[25] + B10_current[1] * data_[19];
  data_[32] = C00_[2] * data_[26] + B10_current[2] * data_[20];
  data_[33] = C00_[3] * data_[27] + B10_current[3] * data_[21];
  data_[34] = C00_[4] * data_[28] + B10_current[4] * data_[22];
  data_[35] = C00_[5] * data_[29] + B10_current[5] * data_[23];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[36] = C00_[0] * data_[30] + B10_current[0] * data_[24];
  data_[37] = C00_[1] * data_[31] + B10_current[1] * data_[25];
  data_[38] = C00_[2] * data_[32] + B10_current[2] * data_[26];
  data_[39] = C00_[3] * data_[33] + B10_current[3] * data_[27];
  data_[40] = C00_[4] * data_[34] + B10_current[4] * data_[28];
  data_[41] = C00_[5] * data_[35] + B10_current[5] * data_[29];

  data_[42] = D00_[0];
  data_[43] = D00_[1];
  data_[44] = D00_[2];
  data_[45] = D00_[3];
  data_[46] = D00_[4];
  data_[47] = D00_[5];

  double cB00_current[6];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];
  cB00_current[3] = B00_[3];
  cB00_current[4] = B00_[4];
  cB00_current[5] = B00_[5];

  data_[48] = C00_[0] * data_[42] + cB00_current[0];
  data_[49] = C00_[1] * data_[43] + cB00_current[1];
  data_[50] = C00_[2] * data_[44] + cB00_current[2];
  data_[51] = C00_[3] * data_[45] + cB00_current[3];
  data_[52] = C00_[4] * data_[46] + cB00_current[4];
  data_[53] = C00_[5] * data_[47] + cB00_current[5];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];

  data_[54] = C00_[0] * data_[48] + B10_current[0] * data_[42] + cB00_current[0] * data_[6];
  data_[55] = C00_[1] * data_[49] + B10_current[1] * data_[43] + cB00_current[1] * data_[7];
  data_[56] = C00_[2] * data_[50] + B10_current[2] * data_[44] + cB00_current[2] * data_[8];
  data_[57] = C00_[3] * data_[51] + B10_current[3] * data_[45] + cB00_current[3] * data_[9];
  data_[58] = C00_[4] * data_[52] + B10_current[4] * data_[46] + cB00_current[4] * data_[10];
  data_[59] = C00_[5] * data_[53] + B10_current[5] * data_[47] + cB00_current[5] * data_[11];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[60] = C00_[0] * data_[54] + B10_current[0] * data_[48] + cB00_current[0] * data_[12];
  data_[61] = C00_[1] * data_[55] + B10_current[1] * data_[49] + cB00_current[1] * data_[13];
  data_[62] = C00_[2] * data_[56] + B10_current[2] * data_[50] + cB00_current[2] * data_[14];
  data_[63] = C00_[3] * data_[57] + B10_current[3] * data_[51] + cB00_current[3] * data_[15];
  data_[64] = C00_[4] * data_[58] + B10_current[4] * data_[52] + cB00_current[4] * data_[16];
  data_[65] = C00_[5] * data_[59] + B10_current[5] * data_[53] + cB00_current[5] * data_[17];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[66] = C00_[0] * data_[60] + B10_current[0] * data_[54] + cB00_current[0] * data_[18];
  data_[67] = C00_[1] * data_[61] + B10_current[1] * data_[55] + cB00_current[1] * data_[19];
  data_[68] = C00_[2] * data_[62] + B10_current[2] * data_[56] + cB00_current[2] * data_[20];
  data_[69] = C00_[3] * data_[63] + B10_current[3] * data_[57] + cB00_current[3] * data_[21];
  data_[70] = C00_[4] * data_[64] + B10_current[4] * data_[58] + cB00_current[4] * data_[22];
  data_[71] = C00_[5] * data_[65] + B10_current[5] * data_[59] + cB00_current[5] * data_[23];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[72] = C00_[0] * data_[66] + B10_current[0] * data_[60] + cB00_current[0] * data_[24];
  data_[73] = C00_[1] * data_[67] + B10_current[1] * data_[61] + cB00_current[1] * data_[25];
  data_[74] = C00_[2] * data_[68] + B10_current[2] * data_[62] + cB00_current[2] * data_[26];
  data_[75] = C00_[3] * data_[69] + B10_current[3] * data_[63] + cB00_current[3] * data_[27];
  data_[76] = C00_[4] * data_[70] + B10_current[4] * data_[64] + cB00_current[4] * data_[28];
  data_[77] = C00_[5] * data_[71] + B10_current[5] * data_[65] + cB00_current[5] * data_[29];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[78] = C00_[0] * data_[72] + B10_current[0] * data_[66] + cB00_current[0] * data_[30];
  data_[79] = C00_[1] * data_[73] + B10_current[1] * data_[67] + cB00_current[1] * data_[31];
  data_[80] = C00_[2] * data_[74] + B10_current[2] * data_[68] + cB00_current[2] * data_[32];
  data_[81] = C00_[3] * data_[75] + B10_current[3] * data_[69] + cB00_current[3] * data_[33];
  data_[82] = C00_[4] * data_[76] + B10_current[4] * data_[70] + cB00_current[4] * data_[34];
  data_[83] = C00_[5] * data_[77] + B10_current[5] * data_[71] + cB00_current[5] * data_[35];

  double B01_current[6];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];
  B01_current[3] = B01_[3];
  B01_current[4] = B01_[4];
  B01_current[5] = B01_[5];

  data_[84] = D00_[0] * data_[42] + B01_current[0];
  data_[85] = D00_[1] * data_[43] + B01_current[1];
  data_[86] = D00_[2] * data_[44] + B01_current[2];
  data_[87] = D00_[3] * data_[45] + B01_current[3];
  data_[88] = D00_[4] * data_[46] + B01_current[4];
  data_[89] = D00_[5] * data_[47] + B01_current[5];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];

  data_[90] = C00_[0] * data_[84] + cB00_current[0] * data_[42];
  data_[91] = C00_[1] * data_[85] + cB00_current[1] * data_[43];
  data_[92] = C00_[2] * data_[86] + cB00_current[2] * data_[44];
  data_[93] = C00_[3] * data_[87] + cB00_current[3] * data_[45];
  data_[94] = C00_[4] * data_[88] + cB00_current[4] * data_[46];
  data_[95] = C00_[5] * data_[89] + cB00_current[5] * data_[47];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];

  data_[96] = C00_[0] * data_[90] + B10_current[0] * data_[84] + cB00_current[0] * data_[48];
  data_[97] = C00_[1] * data_[91] + B10_current[1] * data_[85] + cB00_current[1] * data_[49];
  data_[98] = C00_[2] * data_[92] + B10_current[2] * data_[86] + cB00_current[2] * data_[50];
  data_[99] = C00_[3] * data_[93] + B10_current[3] * data_[87] + cB00_current[3] * data_[51];
  data_[100] = C00_[4] * data_[94] + B10_current[4] * data_[88] + cB00_current[4] * data_[52];
  data_[101] = C00_[5] * data_[95] + B10_current[5] * data_[89] + cB00_current[5] * data_[53];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[102] = C00_[0] * data_[96] + B10_current[0] * data_[90] + cB00_current[0] * data_[54];
  data_[103] = C00_[1] * data_[97] + B10_current[1] * data_[91] + cB00_current[1] * data_[55];
  data_[104] = C00_[2] * data_[98] + B10_current[2] * data_[92] + cB00_current[2] * data_[56];
  data_[105] = C00_[3] * data_[99] + B10_current[3] * data_[93] + cB00_current[3] * data_[57];
  data_[106] = C00_[4] * data_[100] + B10_current[4] * data_[94] + cB00_current[4] * data_[58];
  data_[107] = C00_[5] * data_[101] + B10_current[5] * data_[95] + cB00_current[5] * data_[59];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[108] = C00_[0] * data_[102] + B10_current[0] * data_[96] + cB00_current[0] * data_[60];
  data_[109] = C00_[1] * data_[103] + B10_current[1] * data_[97] + cB00_current[1] * data_[61];
  data_[110] = C00_[2] * data_[104] + B10_current[2] * data_[98] + cB00_current[2] * data_[62];
  data_[111] = C00_[3] * data_[105] + B10_current[3] * data_[99] + cB00_current[3] * data_[63];
  data_[112] = C00_[4] * data_[106] + B10_current[4] * data_[100] + cB00_current[4] * data_[64];
  data_[113] = C00_[5] * data_[107] + B10_current[5] * data_[101] + cB00_current[5] * data_[65];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[114] = C00_[0] * data_[108] + B10_current[0] * data_[102] + cB00_current[0] * data_[66];
  data_[115] = C00_[1] * data_[109] + B10_current[1] * data_[103] + cB00_current[1] * data_[67];
  data_[116] = C00_[2] * data_[110] + B10_current[2] * data_[104] + cB00_current[2] * data_[68];
  data_[117] = C00_[3] * data_[111] + B10_current[3] * data_[105] + cB00_current[3] * data_[69];
  data_[118] = C00_[4] * data_[112] + B10_current[4] * data_[106] + cB00_current[4] * data_[70];
  data_[119] = C00_[5] * data_[113] + B10_current[5] * data_[107] + cB00_current[5] * data_[71];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[120] = C00_[0] * data_[114] + B10_current[0] * data_[108] + cB00_current[0] * data_[72];
  data_[121] = C00_[1] * data_[115] + B10_current[1] * data_[109] + cB00_current[1] * data_[73];
  data_[122] = C00_[2] * data_[116] + B10_current[2] * data_[110] + cB00_current[2] * data_[74];
  data_[123] = C00_[3] * data_[117] + B10_current[3] * data_[111] + cB00_current[3] * data_[75];
  data_[124] = C00_[4] * data_[118] + B10_current[4] * data_[112] + cB00_current[4] * data_[76];
  data_[125] = C00_[5] * data_[119] + B10_current[5] * data_[113] + cB00_current[5] * data_[77];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];

  data_[126] = D00_[0] * data_[84] + B01_current[0] * data_[42];
  data_[127] = D00_[1] * data_[85] + B01_current[1] * data_[43];
  data_[128] = D00_[2] * data_[86] + B01_current[2] * data_[44];
  data_[129] = D00_[3] * data_[87] + B01_current[3] * data_[45];
  data_[130] = D00_[4] * data_[88] + B01_current[4] * data_[46];
  data_[131] = D00_[5] * data_[89] + B01_current[5] * data_[47];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];

  data_[132] = C00_[0] * data_[126] + cB00_current[0] * data_[84];
  data_[133] = C00_[1] * data_[127] + cB00_current[1] * data_[85];
  data_[134] = C00_[2] * data_[128] + cB00_current[2] * data_[86];
  data_[135] = C00_[3] * data_[129] + cB00_current[3] * data_[87];
  data_[136] = C00_[4] * data_[130] + cB00_current[4] * data_[88];
  data_[137] = C00_[5] * data_[131] + cB00_current[5] * data_[89];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];

  data_[138] = C00_[0] * data_[132] + B10_current[0] * data_[126] + cB00_current[0] * data_[90];
  data_[139] = C00_[1] * data_[133] + B10_current[1] * data_[127] + cB00_current[1] * data_[91];
  data_[140] = C00_[2] * data_[134] + B10_current[2] * data_[128] + cB00_current[2] * data_[92];
  data_[141] = C00_[3] * data_[135] + B10_current[3] * data_[129] + cB00_current[3] * data_[93];
  data_[142] = C00_[4] * data_[136] + B10_current[4] * data_[130] + cB00_current[4] * data_[94];
  data_[143] = C00_[5] * data_[137] + B10_current[5] * data_[131] + cB00_current[5] * data_[95];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[144] = C00_[0] * data_[138] + B10_current[0] * data_[132] + cB00_current[0] * data_[96];
  data_[145] = C00_[1] * data_[139] + B10_current[1] * data_[133] + cB00_current[1] * data_[97];
  data_[146] = C00_[2] * data_[140] + B10_current[2] * data_[134] + cB00_current[2] * data_[98];
  data_[147] = C00_[3] * data_[141] + B10_current[3] * data_[135] + cB00_current[3] * data_[99];
  data_[148] = C00_[4] * data_[142] + B10_current[4] * data_[136] + cB00_current[4] * data_[100];
  data_[149] = C00_[5] * data_[143] + B10_current[5] * data_[137] + cB00_current[5] * data_[101];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[150] = C00_[0] * data_[144] + B10_current[0] * data_[138] + cB00_current[0] * data_[102];
  data_[151] = C00_[1] * data_[145] + B10_current[1] * data_[139] + cB00_current[1] * data_[103];
  data_[152] = C00_[2] * data_[146] + B10_current[2] * data_[140] + cB00_current[2] * data_[104];
  data_[153] = C00_[3] * data_[147] + B10_current[3] * data_[141] + cB00_current[3] * data_[105];
  data_[154] = C00_[4] * data_[148] + B10_current[4] * data_[142] + cB00_current[4] * data_[106];
  data_[155] = C00_[5] * data_[149] + B10_current[5] * data_[143] + cB00_current[5] * data_[107];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[156] = C00_[0] * data_[150] + B10_current[0] * data_[144] + cB00_current[0] * data_[108];
  data_[157] = C00_[1] * data_[151] + B10_current[1] * data_[145] + cB00_current[1] * data_[109];
  data_[158] = C00_[2] * data_[152] + B10_current[2] * data_[146] + cB00_current[2] * data_[110];
  data_[159] = C00_[3] * data_[153] + B10_current[3] * data_[147] + cB00_current[3] * data_[111];
  data_[160] = C00_[4] * data_[154] + B10_current[4] * data_[148] + cB00_current[4] * data_[112];
  data_[161] = C00_[5] * data_[155] + B10_current[5] * data_[149] + cB00_current[5] * data_[113];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[162] = C00_[0] * data_[156] + B10_current[0] * data_[150] + cB00_current[0] * data_[114];
  data_[163] = C00_[1] * data_[157] + B10_current[1] * data_[151] + cB00_current[1] * data_[115];
  data_[164] = C00_[2] * data_[158] + B10_current[2] * data_[152] + cB00_current[2] * data_[116];
  data_[165] = C00_[3] * data_[159] + B10_current[3] * data_[153] + cB00_current[3] * data_[117];
  data_[166] = C00_[4] * data_[160] + B10_current[4] * data_[154] + cB00_current[4] * data_[118];
  data_[167] = C00_[5] * data_[161] + B10_current[5] * data_[155] + cB00_current[5] * data_[119];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];
  B01_current[3] += B01_[3];
  B01_current[4] += B01_[4];
  B01_current[5] += B01_[5];

  data_[168] = D00_[0] * data_[126] + B01_current[0] * data_[84];
  data_[169] = D00_[1] * data_[127] + B01_current[1] * data_[85];
  data_[170] = D00_[2] * data_[128] + B01_current[2] * data_[86];
  data_[171] = D00_[3] * data_[129] + B01_current[3] * data_[87];
  data_[172] = D00_[4] * data_[130] + B01_current[4] * data_[88];
  data_[173] = D00_[5] * data_[131] + B01_current[5] * data_[89];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];
  cB00_current[3] += B00_[3];
  cB00_current[4] += B00_[4];
  cB00_current[5] += B00_[5];

  data_[174] = C00_[0] * data_[168] + cB00_current[0] * data_[126];
  data_[175] = C00_[1] * data_[169] + cB00_current[1] * data_[127];
  data_[176] = C00_[2] * data_[170] + cB00_current[2] * data_[128];
  data_[177] = C00_[3] * data_[171] + cB00_current[3] * data_[129];
  data_[178] = C00_[4] * data_[172] + cB00_current[4] * data_[130];
  data_[179] = C00_[5] * data_[173] + cB00_current[5] * data_[131];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];
  B10_current[3] = B10_[3];
  B10_current[4] = B10_[4];
  B10_current[5] = B10_[5];

  data_[180] = C00_[0] * data_[174] + B10_current[0] * data_[168] + cB00_current[0] * data_[132];
  data_[181] = C00_[1] * data_[175] + B10_current[1] * data_[169] + cB00_current[1] * data_[133];
  data_[182] = C00_[2] * data_[176] + B10_current[2] * data_[170] + cB00_current[2] * data_[134];
  data_[183] = C00_[3] * data_[177] + B10_current[3] * data_[171] + cB00_current[3] * data_[135];
  data_[184] = C00_[4] * data_[178] + B10_current[4] * data_[172] + cB00_current[4] * data_[136];
  data_[185] = C00_[5] * data_[179] + B10_current[5] * data_[173] + cB00_current[5] * data_[137];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[186] = C00_[0] * data_[180] + B10_current[0] * data_[174] + cB00_current[0] * data_[138];
  data_[187] = C00_[1] * data_[181] + B10_current[1] * data_[175] + cB00_current[1] * data_[139];
  data_[188] = C00_[2] * data_[182] + B10_current[2] * data_[176] + cB00_current[2] * data_[140];
  data_[189] = C00_[3] * data_[183] + B10_current[3] * data_[177] + cB00_current[3] * data_[141];
  data_[190] = C00_[4] * data_[184] + B10_current[4] * data_[178] + cB00_current[4] * data_[142];
  data_[191] = C00_[5] * data_[185] + B10_current[5] * data_[179] + cB00_current[5] * data_[143];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[192] = C00_[0] * data_[186] + B10_current[0] * data_[180] + cB00_current[0] * data_[144];
  data_[193] = C00_[1] * data_[187] + B10_current[1] * data_[181] + cB00_current[1] * data_[145];
  data_[194] = C00_[2] * data_[188] + B10_current[2] * data_[182] + cB00_current[2] * data_[146];
  data_[195] = C00_[3] * data_[189] + B10_current[3] * data_[183] + cB00_current[3] * data_[147];
  data_[196] = C00_[4] * data_[190] + B10_current[4] * data_[184] + cB00_current[4] * data_[148];
  data_[197] = C00_[5] * data_[191] + B10_current[5] * data_[185] + cB00_current[5] * data_[149];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[198] = C00_[0] * data_[192] + B10_current[0] * data_[186] + cB00_current[0] * data_[150];
  data_[199] = C00_[1] * data_[193] + B10_current[1] * data_[187] + cB00_current[1] * data_[151];
  data_[200] = C00_[2] * data_[194] + B10_current[2] * data_[188] + cB00_current[2] * data_[152];
  data_[201] = C00_[3] * data_[195] + B10_current[3] * data_[189] + cB00_current[3] * data_[153];
  data_[202] = C00_[4] * data_[196] + B10_current[4] * data_[190] + cB00_current[4] * data_[154];
  data_[203] = C00_[5] * data_[197] + B10_current[5] * data_[191] + cB00_current[5] * data_[155];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];
  B10_current[3] += B10_[3];
  B10_current[4] += B10_[4];
  B10_current[5] += B10_[5];

  data_[204] = C00_[0] * data_[198] + B10_current[0] * data_[192] + cB00_current[0] * data_[156];
  data_[205] = C00_[1] * data_[199] + B10_current[1] * data_[193] + cB00_current[1] * data_[157];
  data_[206] = C00_[2] * data_[200] + B10_current[2] * data_[194] + cB00_current[2] * data_[158];
  data_[207] = C00_[3] * data_[201] + B10_current[3] * data_[195] + cB00_current[3] * data_[159];
  data_[208] = C00_[4] * data_[202] + B10_current[4] * data_[196] + cB00_current[4] * data_[160];
  data_[209] = C00_[5] * data_[203] + B10_current[5] * data_[197] + cB00_current[5] * data_[161];
}

