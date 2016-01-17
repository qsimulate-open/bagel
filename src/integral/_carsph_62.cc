//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_62.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/integral/carsphlist.h>
#include <algorithm>

using namespace std;
using namespace bagel;


void CarSphList::carsph_62(const int nloop, const double* source, double* target) {
  const double c25 = 51.553976568253198;
  const double c17 = 40.301597362883768;
  const double c49 = 37.649701194033398;
  const double c38 = 34.369317712168801;
  const double c30 = 29.764702249476645;
  const double c21 = 25.776988284126599;
  const double c66 = 25.099800796022265;
  const double c11 = 23.268138086232856;
  const double c54 = 21.737065119284157;
  const double c16 = 20.150798681441884;
  const double c40 = 19.843134832984429;
  const double c86 = 19.48557158514987;
  const double c45 = 18.824850597016699;
  const double c3 = 17.451103564674643;
  const double c36 = 17.1846588560844;
  const double c31 = 14.882351124738323;
  const double c68 = 14.491376746189438;
  const double c51 = 14.118637947762524;
  const double c14 = 13.433865787627923;
  const double c85 = 12.99038105676658;
  const double c48 = 12.549900398011133;
  const double c9 = 11.634069043116428;
  const double c74 = 11.456439237389599;
  const double c94 = 11.25;
  const double c55 = 10.868532559642079;
  const double c6 = 10.075399340720942;
  const double c41 = 9.9215674164922145;
  const double c81 = 9.742785792574935;
  const double c1 = 8.7255517823373214;
  const double c24 = 8.5923294280422002;
  const double c58 = 8.1513994197315593;
  const double c71 = 7.9372539331937721;
  const double c90 = 7.5;
  const double c52 = 7.245688373094719;
  const double c47 = 7.0593189738812621;
  const double c10 = 6.9804414258698566;
  const double c15 = 6.7169328938139614;
  const double c79 = 6.49519052838329;
  const double c44 = 6.2749501990055663;
  const double c19 = 5.817034521558214;
  const double c75 = 5.7282196186947996;
  const double c92 = 5.625;
  const double c27 = 5.1553976568253201;
  const double c7 = 5.0376996703604711;
  const double c28 = 4.9607837082461073;
  const double c80 = 4.8713928962874675;
  const double c50 = 4.7062126492541747;
  const double c72 = 4.5825756949558398;
  const double c20 = 4.2961647140211001;
  const double c59 = 4.0756997098657797;
  const double c12 = 4.0301597362883772;
  const double c70 = 3.9686269665968861;
  const double c91 = 3.75;
  const double c53 = 3.6228441865473595;
  const double c8 = 3.4902207129349283;
  const double c39 = 3.4369317712168801;
  const double c67 = 3.1374750995027831;
  const double c34 = 2.9764702249476644;
  const double c76 = 2.8641098093473998;
  const double c93 = 2.8125;
  const double c56 = 2.7171331399105196;
  const double c23 = 2.5776988284126601;
  const double c29 = 2.4803918541230536;
  const double c46 = 2.3531063246270874;
  const double c18 = 2.3268138086232857;
  const double c73 = 2.2912878474779199;
  const double c13 = 2.0150798681441886;
  const double c42 = 1.984313483298443;
  const double c69 = 1.8114220932736798;
  const double c84 = 1.7320508075688772;
  const double c37 = 1.71846588560844;
  const double c88 = 1.6237976320958225;
  const double c62 = 1.5687375497513916;
  const double c35 = 1.4882351124738322;
  const double c77 = 1.4320549046736999;
  const double c57 = 1.3585665699552598;
  const double c2 = 1.1634069043116428;
  const double c43 = 0.99215674164922152;
  const double c97 = 0.9375;
  const double c65 = 0.90571104663683988;
  const double c78 = 0.8660254037844386;
  const double c26 = 0.85923294280422002;
  const double c83 = 0.81189881604791125;
  const double c61 = 0.78436877487569578;
  const double c4 = 0.67169328938139616;
  const double c0 = 0.58170345215582142;
  const double c87 = 0.54126587736527421;
  const double c89 = 0.5;
  const double c32 = 0.49607837082461076;
  const double c98 = 0.46875;
  const double c63 = 0.45285552331841994;
  const double c22 = 0.42961647140211001;
  const double c60 = 0.39218438743784789;
  const double c5 = 0.33584664469069808;
  const double c95 = 0.3125;
  const double c82 = 0.2706329386826371;
  const double c33 = 0.24803918541230538;
  const double c64 = 0.22642776165920997;
  const double c96 = 0.15625;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 65, source += 168) {
    target[0] =  c0 * source[0] - c0 * source[2] - c1 * source[12]
                  + c1 * source[14] + c1 * source[24] - c1 * source[26]
                  - c0 * source[36] + c0 * source[38];
    target[1] =  c2 * source[1] - c3 * source[13] + c3 * source[25]
                  - c2 * source[37];
    target[2] =  c2 * source[3] - c3 * source[15] + c3 * source[27]
                  - c2 * source[39];
    target[3] =  c2 * source[4] - c3 * source[16] + c3 * source[28]
                  - c2 * source[40];
    target[4] =  c4 * source[5] - c5 * source[0] - c5 * source[2]
                  - c6 * source[17] + c7 * source[12] + c7 * source[14]
                  + c6 * source[29] - c7 * source[24] - c7 * source[26]
                  - c4 * source[41] + c5 * source[36] + c5 * source[38];
    target[5] =  c8 * source[6] - c8 * source[8] - c9 * source[18]
                  + c9 * source[20] + c8 * source[30] - c8 * source[32];
    target[6] =  c10 * source[7] - c11 * source[19] + c10 * source[31];
    target[7] =  c10 * source[9] - c11 * source[21] + c10 * source[33];
    target[8] =  c10 * source[10] - c11 * source[22] + c10 * source[34];
    target[9] =  c12 * source[11] - c13 * source[6] - c13 * source[8]
                  - c14 * source[23] + c15 * source[18] + c15 * source[20]
                  + c12 * source[35] - c13 * source[30] - c13 * source[32];
    target[10] =  c13 * source[42] - c13 * source[44] - c16 * source[54]
                  + c16 * source[56] + c6 * source[66] - c6 * source[68];
    target[11] =  c12 * source[43] - c17 * source[55] + c16 * source[67];
    target[12] =  c12 * source[45] - c17 * source[57] + c16 * source[69];
    target[13] =  c12 * source[46] - c17 * source[58] + c16 * source[70];
    target[14] =  c18 * source[47] - c2 * source[42] - c2 * source[44]
                  - c11 * source[59] + c9 * source[54] + c9 * source[56]
                  + c9 * source[71] - c19 * source[66] - c19 * source[68];
    target[15] =  c6 * source[48] - c6 * source[50] - c16 * source[60]
                  + c16 * source[62] + c13 * source[72] - c13 * source[74];
    target[16] =  c16 * source[49] - c17 * source[61] + c12 * source[73];
    target[17] =  c16 * source[51] - c17 * source[63] + c12 * source[75];
    target[18] =  c16 * source[52] - c17 * source[64] + c12 * source[76];
    target[19] =  c9 * source[53] - c19 * source[48] - c19 * source[50]
                  - c11 * source[65] + c9 * source[60] + c9 * source[62]
                  + c18 * source[77] - c2 * source[72] - c2 * source[74];
    target[20] =  c20 * source[78] - c20 * source[80] - c21 * source[90]
                  + c21 * source[92] + c20 * source[102] - c20 * source[104]
                  - c22 * source[0] + c22 * source[2] + c23 * source[12]
                  - c23 * source[14] - c22 * source[24] + c22 * source[26]
                  - c22 * source[12] + c22 * source[14] + c23 * source[24]
                  - c23 * source[26] - c22 * source[36] + c22 * source[38];
    target[21] =  c24 * source[79] - c25 * source[91] + c24 * source[103]
                  - c26 * source[1] + c27 * source[13] - c26 * source[25]
                  - c26 * source[13] + c27 * source[25] - c26 * source[37];
    target[22] =  c24 * source[81] - c25 * source[93] + c24 * source[105]
                  - c26 * source[3] + c27 * source[15] - c26 * source[27]
                  - c26 * source[15] + c27 * source[27] - c26 * source[39];
    target[23] =  c24 * source[82] - c25 * source[94] + c24 * source[106]
                  - c26 * source[4] + c27 * source[16] - c26 * source[28]
                  - c26 * source[16] + c27 * source[28] - c26 * source[40];
    target[24] =  c28 * source[83] - c29 * source[78] - c29 * source[80]
                  - c30 * source[95] + c31 * source[90] + c31 * source[92]
                  + c28 * source[107] - c29 * source[102] - c29 * source[104]
                  - c32 * source[5] + c33 * source[0] + c33 * source[2]
                  + c34 * source[17] - c35 * source[12] - c35 * source[14]
                  - c32 * source[29] + c33 * source[24] + c33 * source[26]
                  - c32 * source[17] + c33 * source[12] + c33 * source[14]
                  + c34 * source[29] - c35 * source[24] - c35 * source[26]
                  - c32 * source[41] + c33 * source[36] + c33 * source[38];
    target[25] =  c36 * source[84] - c36 * source[86] - c36 * source[96]
                  + c36 * source[98] - c37 * source[6] + c37 * source[8]
                  + c37 * source[18] - c37 * source[20] - c37 * source[18]
                  + c37 * source[20] + c37 * source[30] - c37 * source[32];
    target[26] =  c38 * source[85] - c38 * source[97] - c39 * source[7]
                  + c39 * source[19] - c39 * source[19] + c39 * source[31];
    target[27] =  c38 * source[87] - c38 * source[99] - c39 * source[9]
                  + c39 * source[21] - c39 * source[21] + c39 * source[33];
    target[28] =  c38 * source[88] - c38 * source[100] - c39 * source[10]
                  + c39 * source[22] - c39 * source[22] + c39 * source[34];
    target[29] =  c40 * source[89] - c41 * source[84] - c41 * source[86]
                  - c40 * source[101] + c41 * source[96] + c41 * source[98]
                  - c42 * source[11] + c43 * source[6] + c43 * source[8]
                  + c42 * source[23] - c43 * source[18] - c43 * source[20]
                  - c42 * source[23] + c43 * source[18] + c43 * source[20]
                  + c42 * source[35] - c43 * source[30] - c43 * source[32];
    target[30] =  c44 * source[108] - c44 * source[110] - c45 * source[120]
                  + c45 * source[122] - c46 * source[42] + c46 * source[44]
                  + c47 * source[54] - c47 * source[56] - c46 * source[54]
                  + c46 * source[56] + c47 * source[66] - c47 * source[68];
    target[31] =  c48 * source[109] - c49 * source[121] - c50 * source[43]
                  + c51 * source[55] - c50 * source[55] + c51 * source[67];
    target[32] =  c48 * source[111] - c49 * source[123] - c50 * source[45]
                  + c51 * source[57] - c50 * source[57] + c51 * source[69];
    target[33] =  c48 * source[112] - c49 * source[124] - c50 * source[46]
                  + c51 * source[58] - c50 * source[58] + c51 * source[70];
    target[34] =  c52 * source[113] - c53 * source[108] - c53 * source[110]
                  - c54 * source[125] + c55 * source[120] + c55 * source[122]
                  - c56 * source[47] + c57 * source[42] + c57 * source[44]
                  + c58 * source[59] - c59 * source[54] - c59 * source[56]
                  - c56 * source[59] + c57 * source[54] + c57 * source[56]
                  + c58 * source[71] - c59 * source[66] - c59 * source[68];
    target[35] =  c45 * source[114] - c45 * source[116] - c44 * source[126]
                  + c44 * source[128] - c47 * source[48] + c47 * source[50]
                  + c46 * source[60] - c46 * source[62] - c47 * source[60]
                  + c47 * source[62] + c46 * source[72] - c46 * source[74];
    target[36] =  c49 * source[115] - c48 * source[127] - c51 * source[49]
                  + c50 * source[61] - c51 * source[61] + c50 * source[73];
    target[37] =  c49 * source[117] - c48 * source[129] - c51 * source[51]
                  + c50 * source[63] - c51 * source[63] + c50 * source[75];
    target[38] =  c49 * source[118] - c48 * source[130] - c51 * source[52]
                  + c50 * source[64] - c51 * source[64] + c50 * source[76];
    target[39] =  c54 * source[119] - c55 * source[114] - c55 * source[116]
                  - c52 * source[131] + c53 * source[126] + c53 * source[128]
                  - c58 * source[53] + c59 * source[48] + c59 * source[50]
                  + c56 * source[65] - c57 * source[60] - c57 * source[62]
                  - c58 * source[65] + c59 * source[60] + c59 * source[62]
                  + c56 * source[77] - c57 * source[72] - c57 * source[74];
    target[40] =  c44 * source[132] - c44 * source[134] - c44 * source[144]
                  + c44 * source[146] - c44 * source[78] + c44 * source[80]
                  + c44 * source[90] - c44 * source[92] - c44 * source[90]
                  + c44 * source[92] + c44 * source[102] - c44 * source[104]
                  + c60 * source[0] - c60 * source[2] - c60 * source[12]
                  + c60 * source[14] + c61 * source[12] - c61 * source[14]
                  - c61 * source[24] + c61 * source[26] + c60 * source[24]
                  - c60 * source[26] - c60 * source[36] + c60 * source[38];
    target[41] =  c48 * source[133] - c48 * source[145] - c48 * source[79]
                  + c48 * source[91] - c48 * source[91] + c48 * source[103]
                  + c61 * source[1] - c61 * source[13] + c62 * source[13]
                  - c62 * source[25] + c61 * source[25] - c61 * source[37];
    target[42] =  c48 * source[135] - c48 * source[147] - c48 * source[81]
                  + c48 * source[93] - c48 * source[93] + c48 * source[105]
                  + c61 * source[3] - c61 * source[15] + c62 * source[15]
                  - c62 * source[27] + c61 * source[27] - c61 * source[39];
    target[43] =  c48 * source[136] - c48 * source[148] - c48 * source[82]
                  + c48 * source[94] - c48 * source[94] + c48 * source[106]
                  + c61 * source[4] - c61 * source[16] + c62 * source[16]
                  - c62 * source[28] + c61 * source[28] - c61 * source[40];
    target[44] =  c52 * source[137] - c53 * source[132] - c53 * source[134]
                  - c52 * source[149] + c53 * source[144] + c53 * source[146]
                  - c52 * source[83] + c53 * source[78] + c53 * source[80]
                  + c52 * source[95] - c53 * source[90] - c53 * source[92]
                  - c52 * source[95] + c53 * source[90] + c53 * source[92]
                  + c52 * source[107] - c53 * source[102] - c53 * source[104]
                  + c63 * source[5] - c64 * source[0] - c64 * source[2]
                  - c63 * source[17] + c64 * source[12] + c64 * source[14]
                  + c65 * source[17] - c63 * source[12] - c63 * source[14]
                  - c65 * source[29] + c63 * source[24] + c63 * source[26]
                  + c63 * source[29] - c64 * source[24] - c64 * source[26]
                  - c63 * source[41] + c64 * source[36] + c64 * source[38];
    target[45] =  c48 * source[138] - c48 * source[140] - c48 * source[84]
                  + c48 * source[86] - c48 * source[96] + c48 * source[98]
                  + c61 * source[6] - c61 * source[8] + c62 * source[18]
                  - c62 * source[20] + c61 * source[30] - c61 * source[32];
    target[46] =  c66 * source[139] - c66 * source[85] - c66 * source[97]
                  + c62 * source[7] + c67 * source[19] + c62 * source[31];
    target[47] =  c66 * source[141] - c66 * source[87] - c66 * source[99]
                  + c62 * source[9] + c67 * source[21] + c62 * source[33];
    target[48] =  c66 * source[142] - c66 * source[88] - c66 * source[100]
                  + c62 * source[10] + c67 * source[22] + c62 * source[34];
    target[49] =  c68 * source[143] - c52 * source[138] - c52 * source[140]
                  - c68 * source[89] + c52 * source[84] + c52 * source[86]
                  - c68 * source[101] + c52 * source[96] + c52 * source[98]
                  + c65 * source[11] - c63 * source[6] - c63 * source[8]
                  + c69 * source[23] - c65 * source[18] - c65 * source[20]
                  + c65 * source[35] - c63 * source[30] - c63 * source[32];
    target[50] =  c70 * source[150] - c70 * source[152] - c41 * source[108]
                  + c41 * source[110] - c41 * source[120] + c41 * source[122]
                  + c29 * source[42] - c29 * source[44] + c28 * source[54]
                  - c28 * source[56] + c29 * source[66] - c29 * source[68];
    target[51] =  c71 * source[151] - c40 * source[109] - c40 * source[121]
                  + c28 * source[43] + c41 * source[55] + c28 * source[67];
    target[52] =  c71 * source[153] - c40 * source[111] - c40 * source[123]
                  + c28 * source[45] + c41 * source[57] + c28 * source[69];
    target[53] =  c71 * source[154] - c40 * source[112] - c40 * source[124]
                  + c28 * source[46] + c41 * source[58] + c28 * source[70];
    target[54] =  c72 * source[155] - c73 * source[150] - c73 * source[152]
                  - c74 * source[113] + c75 * source[108] + c75 * source[110]
                  - c74 * source[125] + c75 * source[120] + c75 * source[122]
                  + c76 * source[47] - c77 * source[42] - c77 * source[44]
                  + c75 * source[59] - c76 * source[54] - c76 * source[56]
                  + c76 * source[71] - c77 * source[66] - c77 * source[68];
    target[55] =  c70 * source[156] - c70 * source[158] - c41 * source[114]
                  + c41 * source[116] - c41 * source[126] + c41 * source[128]
                  + c29 * source[48] - c29 * source[50] + c28 * source[60]
                  - c28 * source[62] + c29 * source[72] - c29 * source[74];
    target[56] =  c71 * source[157] - c40 * source[115] - c40 * source[127]
                  + c28 * source[49] + c41 * source[61] + c28 * source[73];
    target[57] =  c71 * source[159] - c40 * source[117] - c40 * source[129]
                  + c28 * source[51] + c41 * source[63] + c28 * source[75];
    target[58] =  c71 * source[160] - c40 * source[118] - c40 * source[130]
                  + c28 * source[52] + c41 * source[64] + c28 * source[76];
    target[59] =  c72 * source[161] - c73 * source[156] - c73 * source[158]
                  - c74 * source[119] + c75 * source[114] + c75 * source[116]
                  - c74 * source[131] + c75 * source[126] + c75 * source[128]
                  + c76 * source[53] - c77 * source[48] - c77 * source[50]
                  + c75 * source[65] - c76 * source[60] - c76 * source[62]
                  + c76 * source[77] - c77 * source[72] - c77 * source[74];
    target[60] =  c78 * source[162] - c78 * source[164] - c79 * source[132]
                  + c79 * source[134] - c79 * source[144] + c79 * source[146]
                  + c80 * source[78] - c80 * source[80] + c81 * source[90]
                  - c81 * source[92] + c80 * source[102] - c80 * source[104]
                  - c82 * source[0] + c82 * source[2] - c83 * source[12]
                  + c83 * source[14] - c83 * source[24] + c83 * source[26]
                  - c82 * source[36] + c82 * source[38];
    target[61] =  c84 * source[163] - c85 * source[133] - c85 * source[145]
                  + c81 * source[79] + c86 * source[91] + c81 * source[103]
                  - c87 * source[1] - c88 * source[13] - c88 * source[25]
                  - c87 * source[37];
    target[62] =  c84 * source[165] - c85 * source[135] - c85 * source[147]
                  + c81 * source[81] + c86 * source[93] + c81 * source[105]
                  - c87 * source[3] - c88 * source[15] - c88 * source[27]
                  - c87 * source[39];
    target[63] =  c84 * source[166] - c85 * source[136] - c85 * source[148]
                  + c81 * source[82] + c86 * source[94] + c81 * source[106]
                  - c87 * source[4] - c88 * source[16] - c88 * source[28]
                  - c87 * source[40];
    target[64] =  source[167] - c89 * source[162] - c89 * source[164]
                  - c90 * source[137] + c91 * source[132] + c91 * source[134]
                  - c90 * source[149] + c91 * source[144] + c91 * source[146]
                  + c92 * source[83] - c93 * source[78] - c93 * source[80]
                  + c94 * source[95] - c92 * source[90] - c92 * source[92]
                  + c92 * source[107] - c93 * source[102] - c93 * source[104]
                  - c95 * source[5] + c96 * source[0] + c96 * source[2]
                  - c97 * source[17] + c98 * source[12] + c98 * source[14]
                  - c97 * source[29] + c98 * source[24] + c98 * source[26]
                  - c95 * source[41] + c96 * source[36] + c96 * source[38];
  }
}

void CCarSphList::carsph_62(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c25 = 51.553976568253198;
  const double c17 = 40.301597362883768;
  const double c49 = 37.649701194033398;
  const double c38 = 34.369317712168801;
  const double c30 = 29.764702249476645;
  const double c21 = 25.776988284126599;
  const double c66 = 25.099800796022265;
  const double c11 = 23.268138086232856;
  const double c54 = 21.737065119284157;
  const double c16 = 20.150798681441884;
  const double c40 = 19.843134832984429;
  const double c86 = 19.48557158514987;
  const double c45 = 18.824850597016699;
  const double c3 = 17.451103564674643;
  const double c36 = 17.1846588560844;
  const double c31 = 14.882351124738323;
  const double c68 = 14.491376746189438;
  const double c51 = 14.118637947762524;
  const double c14 = 13.433865787627923;
  const double c85 = 12.99038105676658;
  const double c48 = 12.549900398011133;
  const double c9 = 11.634069043116428;
  const double c74 = 11.456439237389599;
  const double c94 = 11.25;
  const double c55 = 10.868532559642079;
  const double c6 = 10.075399340720942;
  const double c41 = 9.9215674164922145;
  const double c81 = 9.742785792574935;
  const double c1 = 8.7255517823373214;
  const double c24 = 8.5923294280422002;
  const double c58 = 8.1513994197315593;
  const double c71 = 7.9372539331937721;
  const double c90 = 7.5;
  const double c52 = 7.245688373094719;
  const double c47 = 7.0593189738812621;
  const double c10 = 6.9804414258698566;
  const double c15 = 6.7169328938139614;
  const double c79 = 6.49519052838329;
  const double c44 = 6.2749501990055663;
  const double c19 = 5.817034521558214;
  const double c75 = 5.7282196186947996;
  const double c92 = 5.625;
  const double c27 = 5.1553976568253201;
  const double c7 = 5.0376996703604711;
  const double c28 = 4.9607837082461073;
  const double c80 = 4.8713928962874675;
  const double c50 = 4.7062126492541747;
  const double c72 = 4.5825756949558398;
  const double c20 = 4.2961647140211001;
  const double c59 = 4.0756997098657797;
  const double c12 = 4.0301597362883772;
  const double c70 = 3.9686269665968861;
  const double c91 = 3.75;
  const double c53 = 3.6228441865473595;
  const double c8 = 3.4902207129349283;
  const double c39 = 3.4369317712168801;
  const double c67 = 3.1374750995027831;
  const double c34 = 2.9764702249476644;
  const double c76 = 2.8641098093473998;
  const double c93 = 2.8125;
  const double c56 = 2.7171331399105196;
  const double c23 = 2.5776988284126601;
  const double c29 = 2.4803918541230536;
  const double c46 = 2.3531063246270874;
  const double c18 = 2.3268138086232857;
  const double c73 = 2.2912878474779199;
  const double c13 = 2.0150798681441886;
  const double c42 = 1.984313483298443;
  const double c69 = 1.8114220932736798;
  const double c84 = 1.7320508075688772;
  const double c37 = 1.71846588560844;
  const double c88 = 1.6237976320958225;
  const double c62 = 1.5687375497513916;
  const double c35 = 1.4882351124738322;
  const double c77 = 1.4320549046736999;
  const double c57 = 1.3585665699552598;
  const double c2 = 1.1634069043116428;
  const double c43 = 0.99215674164922152;
  const double c97 = 0.9375;
  const double c65 = 0.90571104663683988;
  const double c78 = 0.8660254037844386;
  const double c26 = 0.85923294280422002;
  const double c83 = 0.81189881604791125;
  const double c61 = 0.78436877487569578;
  const double c4 = 0.67169328938139616;
  const double c0 = 0.58170345215582142;
  const double c87 = 0.54126587736527421;
  const double c89 = 0.5;
  const double c32 = 0.49607837082461076;
  const double c98 = 0.46875;
  const double c63 = 0.45285552331841994;
  const double c22 = 0.42961647140211001;
  const double c60 = 0.39218438743784789;
  const double c5 = 0.33584664469069808;
  const double c95 = 0.3125;
  const double c82 = 0.2706329386826371;
  const double c33 = 0.24803918541230538;
  const double c64 = 0.22642776165920997;
  const double c96 = 0.15625;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 65, source += 168) {
    target[0] =  c0 * source[0] - c0 * source[2] - c1 * source[12]
                  + c1 * source[14] + c1 * source[24] - c1 * source[26]
                  - c0 * source[36] + c0 * source[38];
    target[1] =  c2 * source[1] - c3 * source[13] + c3 * source[25]
                  - c2 * source[37];
    target[2] =  c2 * source[3] - c3 * source[15] + c3 * source[27]
                  - c2 * source[39];
    target[3] =  c2 * source[4] - c3 * source[16] + c3 * source[28]
                  - c2 * source[40];
    target[4] =  c4 * source[5] - c5 * source[0] - c5 * source[2]
                  - c6 * source[17] + c7 * source[12] + c7 * source[14]
                  + c6 * source[29] - c7 * source[24] - c7 * source[26]
                  - c4 * source[41] + c5 * source[36] + c5 * source[38];
    target[5] =  c8 * source[6] - c8 * source[8] - c9 * source[18]
                  + c9 * source[20] + c8 * source[30] - c8 * source[32];
    target[6] =  c10 * source[7] - c11 * source[19] + c10 * source[31];
    target[7] =  c10 * source[9] - c11 * source[21] + c10 * source[33];
    target[8] =  c10 * source[10] - c11 * source[22] + c10 * source[34];
    target[9] =  c12 * source[11] - c13 * source[6] - c13 * source[8]
                  - c14 * source[23] + c15 * source[18] + c15 * source[20]
                  + c12 * source[35] - c13 * source[30] - c13 * source[32];
    target[10] =  c13 * source[42] - c13 * source[44] - c16 * source[54]
                  + c16 * source[56] + c6 * source[66] - c6 * source[68];
    target[11] =  c12 * source[43] - c17 * source[55] + c16 * source[67];
    target[12] =  c12 * source[45] - c17 * source[57] + c16 * source[69];
    target[13] =  c12 * source[46] - c17 * source[58] + c16 * source[70];
    target[14] =  c18 * source[47] - c2 * source[42] - c2 * source[44]
                  - c11 * source[59] + c9 * source[54] + c9 * source[56]
                  + c9 * source[71] - c19 * source[66] - c19 * source[68];
    target[15] =  c6 * source[48] - c6 * source[50] - c16 * source[60]
                  + c16 * source[62] + c13 * source[72] - c13 * source[74];
    target[16] =  c16 * source[49] - c17 * source[61] + c12 * source[73];
    target[17] =  c16 * source[51] - c17 * source[63] + c12 * source[75];
    target[18] =  c16 * source[52] - c17 * source[64] + c12 * source[76];
    target[19] =  c9 * source[53] - c19 * source[48] - c19 * source[50]
                  - c11 * source[65] + c9 * source[60] + c9 * source[62]
                  + c18 * source[77] - c2 * source[72] - c2 * source[74];
    target[20] =  c20 * source[78] - c20 * source[80] - c21 * source[90]
                  + c21 * source[92] + c20 * source[102] - c20 * source[104]
                  - c22 * source[0] + c22 * source[2] + c23 * source[12]
                  - c23 * source[14] - c22 * source[24] + c22 * source[26]
                  - c22 * source[12] + c22 * source[14] + c23 * source[24]
                  - c23 * source[26] - c22 * source[36] + c22 * source[38];
    target[21] =  c24 * source[79] - c25 * source[91] + c24 * source[103]
                  - c26 * source[1] + c27 * source[13] - c26 * source[25]
                  - c26 * source[13] + c27 * source[25] - c26 * source[37];
    target[22] =  c24 * source[81] - c25 * source[93] + c24 * source[105]
                  - c26 * source[3] + c27 * source[15] - c26 * source[27]
                  - c26 * source[15] + c27 * source[27] - c26 * source[39];
    target[23] =  c24 * source[82] - c25 * source[94] + c24 * source[106]
                  - c26 * source[4] + c27 * source[16] - c26 * source[28]
                  - c26 * source[16] + c27 * source[28] - c26 * source[40];
    target[24] =  c28 * source[83] - c29 * source[78] - c29 * source[80]
                  - c30 * source[95] + c31 * source[90] + c31 * source[92]
                  + c28 * source[107] - c29 * source[102] - c29 * source[104]
                  - c32 * source[5] + c33 * source[0] + c33 * source[2]
                  + c34 * source[17] - c35 * source[12] - c35 * source[14]
                  - c32 * source[29] + c33 * source[24] + c33 * source[26]
                  - c32 * source[17] + c33 * source[12] + c33 * source[14]
                  + c34 * source[29] - c35 * source[24] - c35 * source[26]
                  - c32 * source[41] + c33 * source[36] + c33 * source[38];
    target[25] =  c36 * source[84] - c36 * source[86] - c36 * source[96]
                  + c36 * source[98] - c37 * source[6] + c37 * source[8]
                  + c37 * source[18] - c37 * source[20] - c37 * source[18]
                  + c37 * source[20] + c37 * source[30] - c37 * source[32];
    target[26] =  c38 * source[85] - c38 * source[97] - c39 * source[7]
                  + c39 * source[19] - c39 * source[19] + c39 * source[31];
    target[27] =  c38 * source[87] - c38 * source[99] - c39 * source[9]
                  + c39 * source[21] - c39 * source[21] + c39 * source[33];
    target[28] =  c38 * source[88] - c38 * source[100] - c39 * source[10]
                  + c39 * source[22] - c39 * source[22] + c39 * source[34];
    target[29] =  c40 * source[89] - c41 * source[84] - c41 * source[86]
                  - c40 * source[101] + c41 * source[96] + c41 * source[98]
                  - c42 * source[11] + c43 * source[6] + c43 * source[8]
                  + c42 * source[23] - c43 * source[18] - c43 * source[20]
                  - c42 * source[23] + c43 * source[18] + c43 * source[20]
                  + c42 * source[35] - c43 * source[30] - c43 * source[32];
    target[30] =  c44 * source[108] - c44 * source[110] - c45 * source[120]
                  + c45 * source[122] - c46 * source[42] + c46 * source[44]
                  + c47 * source[54] - c47 * source[56] - c46 * source[54]
                  + c46 * source[56] + c47 * source[66] - c47 * source[68];
    target[31] =  c48 * source[109] - c49 * source[121] - c50 * source[43]
                  + c51 * source[55] - c50 * source[55] + c51 * source[67];
    target[32] =  c48 * source[111] - c49 * source[123] - c50 * source[45]
                  + c51 * source[57] - c50 * source[57] + c51 * source[69];
    target[33] =  c48 * source[112] - c49 * source[124] - c50 * source[46]
                  + c51 * source[58] - c50 * source[58] + c51 * source[70];
    target[34] =  c52 * source[113] - c53 * source[108] - c53 * source[110]
                  - c54 * source[125] + c55 * source[120] + c55 * source[122]
                  - c56 * source[47] + c57 * source[42] + c57 * source[44]
                  + c58 * source[59] - c59 * source[54] - c59 * source[56]
                  - c56 * source[59] + c57 * source[54] + c57 * source[56]
                  + c58 * source[71] - c59 * source[66] - c59 * source[68];
    target[35] =  c45 * source[114] - c45 * source[116] - c44 * source[126]
                  + c44 * source[128] - c47 * source[48] + c47 * source[50]
                  + c46 * source[60] - c46 * source[62] - c47 * source[60]
                  + c47 * source[62] + c46 * source[72] - c46 * source[74];
    target[36] =  c49 * source[115] - c48 * source[127] - c51 * source[49]
                  + c50 * source[61] - c51 * source[61] + c50 * source[73];
    target[37] =  c49 * source[117] - c48 * source[129] - c51 * source[51]
                  + c50 * source[63] - c51 * source[63] + c50 * source[75];
    target[38] =  c49 * source[118] - c48 * source[130] - c51 * source[52]
                  + c50 * source[64] - c51 * source[64] + c50 * source[76];
    target[39] =  c54 * source[119] - c55 * source[114] - c55 * source[116]
                  - c52 * source[131] + c53 * source[126] + c53 * source[128]
                  - c58 * source[53] + c59 * source[48] + c59 * source[50]
                  + c56 * source[65] - c57 * source[60] - c57 * source[62]
                  - c58 * source[65] + c59 * source[60] + c59 * source[62]
                  + c56 * source[77] - c57 * source[72] - c57 * source[74];
    target[40] =  c44 * source[132] - c44 * source[134] - c44 * source[144]
                  + c44 * source[146] - c44 * source[78] + c44 * source[80]
                  + c44 * source[90] - c44 * source[92] - c44 * source[90]
                  + c44 * source[92] + c44 * source[102] - c44 * source[104]
                  + c60 * source[0] - c60 * source[2] - c60 * source[12]
                  + c60 * source[14] + c61 * source[12] - c61 * source[14]
                  - c61 * source[24] + c61 * source[26] + c60 * source[24]
                  - c60 * source[26] - c60 * source[36] + c60 * source[38];
    target[41] =  c48 * source[133] - c48 * source[145] - c48 * source[79]
                  + c48 * source[91] - c48 * source[91] + c48 * source[103]
                  + c61 * source[1] - c61 * source[13] + c62 * source[13]
                  - c62 * source[25] + c61 * source[25] - c61 * source[37];
    target[42] =  c48 * source[135] - c48 * source[147] - c48 * source[81]
                  + c48 * source[93] - c48 * source[93] + c48 * source[105]
                  + c61 * source[3] - c61 * source[15] + c62 * source[15]
                  - c62 * source[27] + c61 * source[27] - c61 * source[39];
    target[43] =  c48 * source[136] - c48 * source[148] - c48 * source[82]
                  + c48 * source[94] - c48 * source[94] + c48 * source[106]
                  + c61 * source[4] - c61 * source[16] + c62 * source[16]
                  - c62 * source[28] + c61 * source[28] - c61 * source[40];
    target[44] =  c52 * source[137] - c53 * source[132] - c53 * source[134]
                  - c52 * source[149] + c53 * source[144] + c53 * source[146]
                  - c52 * source[83] + c53 * source[78] + c53 * source[80]
                  + c52 * source[95] - c53 * source[90] - c53 * source[92]
                  - c52 * source[95] + c53 * source[90] + c53 * source[92]
                  + c52 * source[107] - c53 * source[102] - c53 * source[104]
                  + c63 * source[5] - c64 * source[0] - c64 * source[2]
                  - c63 * source[17] + c64 * source[12] + c64 * source[14]
                  + c65 * source[17] - c63 * source[12] - c63 * source[14]
                  - c65 * source[29] + c63 * source[24] + c63 * source[26]
                  + c63 * source[29] - c64 * source[24] - c64 * source[26]
                  - c63 * source[41] + c64 * source[36] + c64 * source[38];
    target[45] =  c48 * source[138] - c48 * source[140] - c48 * source[84]
                  + c48 * source[86] - c48 * source[96] + c48 * source[98]
                  + c61 * source[6] - c61 * source[8] + c62 * source[18]
                  - c62 * source[20] + c61 * source[30] - c61 * source[32];
    target[46] =  c66 * source[139] - c66 * source[85] - c66 * source[97]
                  + c62 * source[7] + c67 * source[19] + c62 * source[31];
    target[47] =  c66 * source[141] - c66 * source[87] - c66 * source[99]
                  + c62 * source[9] + c67 * source[21] + c62 * source[33];
    target[48] =  c66 * source[142] - c66 * source[88] - c66 * source[100]
                  + c62 * source[10] + c67 * source[22] + c62 * source[34];
    target[49] =  c68 * source[143] - c52 * source[138] - c52 * source[140]
                  - c68 * source[89] + c52 * source[84] + c52 * source[86]
                  - c68 * source[101] + c52 * source[96] + c52 * source[98]
                  + c65 * source[11] - c63 * source[6] - c63 * source[8]
                  + c69 * source[23] - c65 * source[18] - c65 * source[20]
                  + c65 * source[35] - c63 * source[30] - c63 * source[32];
    target[50] =  c70 * source[150] - c70 * source[152] - c41 * source[108]
                  + c41 * source[110] - c41 * source[120] + c41 * source[122]
                  + c29 * source[42] - c29 * source[44] + c28 * source[54]
                  - c28 * source[56] + c29 * source[66] - c29 * source[68];
    target[51] =  c71 * source[151] - c40 * source[109] - c40 * source[121]
                  + c28 * source[43] + c41 * source[55] + c28 * source[67];
    target[52] =  c71 * source[153] - c40 * source[111] - c40 * source[123]
                  + c28 * source[45] + c41 * source[57] + c28 * source[69];
    target[53] =  c71 * source[154] - c40 * source[112] - c40 * source[124]
                  + c28 * source[46] + c41 * source[58] + c28 * source[70];
    target[54] =  c72 * source[155] - c73 * source[150] - c73 * source[152]
                  - c74 * source[113] + c75 * source[108] + c75 * source[110]
                  - c74 * source[125] + c75 * source[120] + c75 * source[122]
                  + c76 * source[47] - c77 * source[42] - c77 * source[44]
                  + c75 * source[59] - c76 * source[54] - c76 * source[56]
                  + c76 * source[71] - c77 * source[66] - c77 * source[68];
    target[55] =  c70 * source[156] - c70 * source[158] - c41 * source[114]
                  + c41 * source[116] - c41 * source[126] + c41 * source[128]
                  + c29 * source[48] - c29 * source[50] + c28 * source[60]
                  - c28 * source[62] + c29 * source[72] - c29 * source[74];
    target[56] =  c71 * source[157] - c40 * source[115] - c40 * source[127]
                  + c28 * source[49] + c41 * source[61] + c28 * source[73];
    target[57] =  c71 * source[159] - c40 * source[117] - c40 * source[129]
                  + c28 * source[51] + c41 * source[63] + c28 * source[75];
    target[58] =  c71 * source[160] - c40 * source[118] - c40 * source[130]
                  + c28 * source[52] + c41 * source[64] + c28 * source[76];
    target[59] =  c72 * source[161] - c73 * source[156] - c73 * source[158]
                  - c74 * source[119] + c75 * source[114] + c75 * source[116]
                  - c74 * source[131] + c75 * source[126] + c75 * source[128]
                  + c76 * source[53] - c77 * source[48] - c77 * source[50]
                  + c75 * source[65] - c76 * source[60] - c76 * source[62]
                  + c76 * source[77] - c77 * source[72] - c77 * source[74];
    target[60] =  c78 * source[162] - c78 * source[164] - c79 * source[132]
                  + c79 * source[134] - c79 * source[144] + c79 * source[146]
                  + c80 * source[78] - c80 * source[80] + c81 * source[90]
                  - c81 * source[92] + c80 * source[102] - c80 * source[104]
                  - c82 * source[0] + c82 * source[2] - c83 * source[12]
                  + c83 * source[14] - c83 * source[24] + c83 * source[26]
                  - c82 * source[36] + c82 * source[38];
    target[61] =  c84 * source[163] - c85 * source[133] - c85 * source[145]
                  + c81 * source[79] + c86 * source[91] + c81 * source[103]
                  - c87 * source[1] - c88 * source[13] - c88 * source[25]
                  - c87 * source[37];
    target[62] =  c84 * source[165] - c85 * source[135] - c85 * source[147]
                  + c81 * source[81] + c86 * source[93] + c81 * source[105]
                  - c87 * source[3] - c88 * source[15] - c88 * source[27]
                  - c87 * source[39];
    target[63] =  c84 * source[166] - c85 * source[136] - c85 * source[148]
                  + c81 * source[82] + c86 * source[94] + c81 * source[106]
                  - c87 * source[4] - c88 * source[16] - c88 * source[28]
                  - c87 * source[40];
    target[64] =  source[167] - c89 * source[162] - c89 * source[164]
                  - c90 * source[137] + c91 * source[132] + c91 * source[134]
                  - c90 * source[149] + c91 * source[144] + c91 * source[146]
                  + c92 * source[83] - c93 * source[78] - c93 * source[80]
                  + c94 * source[95] - c92 * source[90] - c92 * source[92]
                  + c92 * source[107] - c93 * source[102] - c93 * source[104]
                  - c95 * source[5] + c96 * source[0] + c96 * source[2]
                  - c97 * source[17] + c98 * source[12] + c98 * source[14]
                  - c97 * source[29] + c98 * source[24] + c98 * source[26]
                  - c95 * source[41] + c96 * source[36] + c96 * source[38];
  }
}

