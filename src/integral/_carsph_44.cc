//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_44.cc
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


void CarSphList::carsph_44(const int nloop, const double* source, double* target) {
  const double c71 = 45;
  const double c44 = 42.093645601206838;
  const double c40 = 39.375;
  const double c14 = 29.764702249476645;
  const double c8 = 27.84232950922031;
  const double c56 = 22.5;
  const double c74 = 21.213203435596427;
  const double c42 = 21.046822800603419;
  const double c79 = 20.124611797498108;
  const double c32 = 19.843134832984429;
  const double c2 = 19.6875;
  const double c50 = 18.824850597016699;
  const double c29 = 18.561553006146873;
  const double c75 = 15.90990257669732;
  const double c11 = 14.882351124738323;
  const double c17 = 14.031215200402281;
  const double c24 = 13.311179511974137;
  const double c4 = 13.125;
  const double c53 = 11.25;
  const double c59 = 10.606601717798213;
  const double c18 = 10.523411400301709;
  const double c64 = 10.062305898749054;
  const double c83 = 10;
  const double c30 = 9.9215674164922145;
  const double c86 = 9.4868329805051381;
  const double c34 = 9.354143466934854;
  const double c7 = 9.2807765030734366;
  const double c95 = 9;
  const double c37 = 8.8741196746494246;
  const double c27 = 8.75;
  const double c60 = 7.9549512883486599;
  const double c72 = 7.5;
  const double c89 = 7.1151247353788536;
  const double c35 = 7.0156076002011405;
  const double c78 = 6.7082039324993694;
  const double c45 = 6.6143782776614763;
  const double c47 = 6.2749501990055663;
  const double c28 = 6.1871843353822911;
  const double c84 = 5.625;
  const double c80 = 5.0311529493745271;
  const double c12 = 4.9607837082461073;
  const double c52 = 4.7062126492541747;
  const double c6 = 4.6403882515367183;
  const double c23 = 4.4370598373247123;
  const double c39 = 4.375;
  const double c57 = 3.75;
  const double c76 = 3.5355339059327378;
  const double c43 = 3.5078038001005702;
  const double c63 = 3.3541019662496847;
  const double c26 = 3.3277948779935342;
  const double c33 = 3.3071891388307382;
  const double c1 = 3.28125;
  const double c85 = 3.1622776601683795;
  const double c92 = 3;
  const double c36 = 2.9580398915498081;
  const double c77 = 2.6516504294495533;
  const double c66 = 2.5155764746872635;
  const double c9 = 2.4803918541230536;
  const double c88 = 2.3717082451262845;
  const double c51 = 2.3531063246270874;
  const double c15 = 2.3385358667337135;
  const double c97 = 2.25;
  const double c20 = 2.2185299186623562;
  const double c3 = 2.1875;
  const double c46 = 2.0916500663351889;
  const double c54 = 1.875;
  const double c91 = 1.7787811838447134;
  const double c61 = 1.7677669529663689;
  const double c16 = 1.7539019000502851;
  const double c68 = 1.6770509831248424;
  const double c25 = 1.6638974389967671;
  const double c31 = 1.6535945694153691;
  const double c49 = 1.5687375497513916;
  const double c5 = 1.5467960838455728;
  const double c62 = 1.3258252147247767;
  const double c65 = 1.2577882373436318;
  const double c73 = 1.25;
  const double c87 = 1.1858541225631423;
  const double c41 = 1.1692679333668567;
  const double c96 = 1.125;
  const double c81 = 1.1180339887498949;
  const double c38 = 1.1092649593311781;
  const double c90 = 0.8893905919223567;
  const double c82 = 0.83852549156242118;
  const double c13 = 0.82679728470768454;
  const double c48 = 0.78436877487569578;
  const double c94 = 0.75;
  const double c19 = 0.73950997288745202;
  const double c58 = 0.625;
  const double c100 = 0.5625;
  const double c67 = 0.55901699437494745;
  const double c22 = 0.55463247966558904;
  const double c0 = 0.546875;
  const double c70 = 0.41926274578121059;
  const double c10 = 0.41339864235384227;
  const double c93 = 0.375;
  const double c55 = 0.3125;
  const double c99 = 0.28125;
  const double c21 = 0.27731623983279452;
  const double c69 = 0.20963137289060529;
  const double c98 = 0.140625;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 81, source += 225) {
    target[0] =  c0 * source[0] - c1 * source[2] + c0 * source[4]
                  - c1 * source[30] + c2 * source[32] - c1 * source[34]
                  + c0 * source[60] - c1 * source[62] + c0 * source[64];
    target[1] =  c3 * source[1] - c3 * source[3] - c4 * source[31]
                  + c4 * source[33] + c3 * source[61] - c3 * source[63];
    target[2] =  c5 * source[5] - c6 * source[7] - c7 * source[35]
                  + c8 * source[37] + c5 * source[65] - c6 * source[67];
    target[3] =  c6 * source[6] - c5 * source[8] - c8 * source[36]
                  + c7 * source[38] + c6 * source[66] - c5 * source[68];
    target[4] =  c9 * source[9] - c9 * source[11] - c10 * source[0]
                  + c10 * source[2] - c10 * source[2] + c10 * source[4]
                  - c11 * source[39] + c11 * source[41] + c9 * source[30]
                  - c9 * source[32] + c9 * source[32] - c9 * source[34]
                  + c9 * source[69] - c9 * source[71] - c10 * source[60]
                  + c10 * source[62] - c10 * source[62] + c10 * source[64];
    target[5] =  c12 * source[10] - c13 * source[1] - c13 * source[3]
                  - c14 * source[40] + c12 * source[31] + c12 * source[33]
                  + c12 * source[70] - c13 * source[61] - c13 * source[63];
    target[6] =  c15 * source[12] - c16 * source[5] - c16 * source[7]
                  - c17 * source[42] + c18 * source[35] + c18 * source[37]
                  + c15 * source[72] - c16 * source[65] - c16 * source[67];
    target[7] =  c15 * source[13] - c16 * source[6] - c16 * source[8]
                  - c17 * source[43] + c18 * source[36] + c18 * source[38]
                  + c15 * source[73] - c16 * source[66] - c16 * source[68];
    target[8] =  c19 * source[14] - c20 * source[9] - c20 * source[11]
                  + c21 * source[0] + c22 * source[2] + c21 * source[4]
                  - c23 * source[44] + c24 * source[39] + c24 * source[41]
                  - c25 * source[30] - c26 * source[32] - c25 * source[34]
                  + c19 * source[74] - c20 * source[69] - c20 * source[71]
                  + c21 * source[60] + c22 * source[62] + c21 * source[64];
    target[9] =  c3 * source[15] - c4 * source[17] + c3 * source[19]
                  - c3 * source[45] + c4 * source[47] - c3 * source[49];
    target[10] =  c27 * source[16] - c27 * source[18] - c27 * source[46]
                  + c27 * source[48];
    target[11] =  c28 * source[20] - c29 * source[22] - c28 * source[50]
                  + c29 * source[52];
    target[12] =  c29 * source[21] - c28 * source[23] - c29 * source[51]
                  + c28 * source[53];
    target[13] =  c30 * source[24] - c30 * source[26] - c31 * source[15]
                  + c31 * source[17] - c31 * source[17] + c31 * source[19]
                  - c30 * source[54] + c30 * source[56] + c31 * source[45]
                  - c31 * source[47] + c31 * source[47] - c31 * source[49];
    target[14] =  c32 * source[25] - c33 * source[16] - c33 * source[18]
                  - c32 * source[55] + c33 * source[46] + c33 * source[48];
    target[15] =  c34 * source[27] - c35 * source[20] - c35 * source[22]
                  - c34 * source[57] + c35 * source[50] + c35 * source[52];
    target[16] =  c34 * source[28] - c35 * source[21] - c35 * source[23]
                  - c34 * source[58] + c35 * source[51] + c35 * source[53];
    target[17] =  c36 * source[29] - c37 * source[24] - c37 * source[26]
                  + c38 * source[15] + c20 * source[17] + c38 * source[19]
                  - c36 * source[59] + c37 * source[54] + c37 * source[56]
                  - c38 * source[45] - c20 * source[47] - c38 * source[49];
    target[18] =  c5 * source[75] - c7 * source[77] + c5 * source[79]
                  - c6 * source[105] + c8 * source[107] - c6 * source[109];
    target[19] =  c28 * source[76] - c28 * source[78] - c29 * source[106]
                  + c29 * source[108];
    target[20] =  c39 * source[80] - c4 * source[82] - c4 * source[110]
                  + c40 * source[112];
    target[21] =  c4 * source[81] - c39 * source[83] - c40 * source[111]
                  + c4 * source[113];
    target[22] =  c35 * source[84] - c35 * source[86] - c41 * source[75]
                  + c41 * source[77] - c41 * source[77] + c41 * source[79]
                  - c42 * source[114] + c42 * source[116] + c43 * source[105]
                  - c43 * source[107] + c43 * source[107] - c43 * source[109];
    target[23] =  c17 * source[85] - c15 * source[76] - c15 * source[78]
                  - c44 * source[115] + c35 * source[106] + c35 * source[108];
    target[24] =  c45 * source[87] - c12 * source[80] - c12 * source[82]
                  - c32 * source[117] + c11 * source[110] + c11 * source[112];
    target[25] =  c45 * source[88] - c12 * source[81] - c12 * source[83]
                  - c32 * source[118] + c11 * source[111] + c11 * source[113];
    target[26] =  c46 * source[89] - c47 * source[84] - c47 * source[86]
                  + c48 * source[75] + c49 * source[77] + c48 * source[79]
                  - c47 * source[119] + c50 * source[114] + c50 * source[116]
                  - c51 * source[105] - c52 * source[107] - c51 * source[109];
    target[27] =  c6 * source[90] - c8 * source[92] + c6 * source[94]
                  - c5 * source[120] + c7 * source[122] - c5 * source[124];
    target[28] =  c29 * source[91] - c29 * source[93] - c28 * source[121]
                  + c28 * source[123];
    target[29] =  c4 * source[95] - c40 * source[97] - c39 * source[125]
                  + c4 * source[127];
    target[30] =  c40 * source[96] - c4 * source[98] - c4 * source[126]
                  + c39 * source[128];
    target[31] =  c42 * source[99] - c42 * source[101] - c43 * source[90]
                  + c43 * source[92] - c43 * source[92] + c43 * source[94]
                  - c35 * source[129] + c35 * source[131] + c41 * source[120]
                  - c41 * source[122] + c41 * source[122] - c41 * source[124];
    target[32] =  c44 * source[100] - c35 * source[91] - c35 * source[93]
                  - c17 * source[130] + c15 * source[121] + c15 * source[123];
    target[33] =  c32 * source[102] - c11 * source[95] - c11 * source[97]
                  - c45 * source[132] + c12 * source[125] + c12 * source[127];
    target[34] =  c32 * source[103] - c11 * source[96] - c11 * source[98]
                  - c45 * source[133] + c12 * source[126] + c12 * source[128];
    target[35] =  c47 * source[104] - c50 * source[99] - c50 * source[101]
                  + c51 * source[90] + c52 * source[92] + c51 * source[94]
                  - c46 * source[134] + c47 * source[129] + c47 * source[131]
                  - c48 * source[120] - c49 * source[122] - c48 * source[124];
    target[36] =  c9 * source[135] - c11 * source[137] + c9 * source[139]
                  - c9 * source[165] + c11 * source[167] - c9 * source[169]
                  - c10 * source[0] + c9 * source[2] - c10 * source[4]
                  + c10 * source[30] - c9 * source[32] + c10 * source[34]
                  - c10 * source[30] + c9 * source[32] - c10 * source[34]
                  + c10 * source[60] - c9 * source[62] + c10 * source[64];
    target[37] =  c30 * source[136] - c30 * source[138] - c30 * source[166]
                  + c30 * source[168] - c31 * source[1] + c31 * source[3]
                  + c31 * source[31] - c31 * source[33] - c31 * source[31]
                  + c31 * source[33] + c31 * source[61] - c31 * source[63];
    target[38] =  c35 * source[140] - c42 * source[142] - c35 * source[170]
                  + c42 * source[172] - c41 * source[5] + c43 * source[7]
                  + c41 * source[35] - c43 * source[37] - c41 * source[35]
                  + c43 * source[37] + c41 * source[65] - c43 * source[67];
    target[39] =  c42 * source[141] - c35 * source[143] - c42 * source[171]
                  + c35 * source[173] - c43 * source[6] + c41 * source[8]
                  + c43 * source[36] - c41 * source[38] - c43 * source[36]
                  + c41 * source[38] + c43 * source[66] - c41 * source[68];
    target[40] =  c53 * source[144] - c53 * source[146] - c54 * source[135]
                  + c54 * source[137] - c54 * source[137] + c54 * source[139]
                  - c53 * source[174] + c53 * source[176] + c54 * source[165]
                  - c54 * source[167] + c54 * source[167] - c54 * source[169]
                  - c54 * source[9] + c54 * source[11] + c55 * source[0]
                  - c55 * source[2] + c55 * source[2] - c55 * source[4]
                  + c54 * source[39] - c54 * source[41] - c55 * source[30]
                  + c55 * source[32] - c55 * source[32] + c55 * source[34]
                  - c54 * source[39] + c54 * source[41] + c55 * source[30]
                  - c55 * source[32] + c55 * source[32] - c55 * source[34]
                  + c54 * source[69] - c54 * source[71] - c55 * source[60]
                  + c55 * source[62] - c55 * source[62] + c55 * source[64];
    target[41] =  c56 * source[145] - c57 * source[136] - c57 * source[138]
                  - c56 * source[175] + c57 * source[166] + c57 * source[168]
                  - c57 * source[10] + c58 * source[1] + c58 * source[3]
                  + c57 * source[40] - c58 * source[31] - c58 * source[33]
                  - c57 * source[40] + c58 * source[31] + c58 * source[33]
                  + c57 * source[70] - c58 * source[61] - c58 * source[63];
    target[42] =  c59 * source[147] - c60 * source[140] - c60 * source[142]
                  - c59 * source[177] + c60 * source[170] + c60 * source[172]
                  - c61 * source[12] + c62 * source[5] + c62 * source[7]
                  + c61 * source[42] - c62 * source[35] - c62 * source[37]
                  - c61 * source[42] + c62 * source[35] + c62 * source[37]
                  + c61 * source[72] - c62 * source[65] - c62 * source[67];
    target[43] =  c59 * source[148] - c60 * source[141] - c60 * source[143]
                  - c59 * source[178] + c60 * source[171] + c60 * source[173]
                  - c61 * source[13] + c62 * source[6] + c62 * source[8]
                  + c61 * source[43] - c62 * source[36] - c62 * source[38]
                  - c61 * source[43] + c62 * source[36] + c62 * source[38]
                  + c61 * source[73] - c62 * source[66] - c62 * source[68];
    target[44] =  c63 * source[149] - c64 * source[144] - c64 * source[146]
                  + c65 * source[135] + c66 * source[137] + c65 * source[139]
                  - c63 * source[179] + c64 * source[174] + c64 * source[176]
                  - c65 * source[165] - c66 * source[167] - c65 * source[169]
                  - c67 * source[14] + c68 * source[9] + c68 * source[11]
                  - c69 * source[0] - c70 * source[2] - c69 * source[4]
                  + c67 * source[44] - c68 * source[39] - c68 * source[41]
                  + c69 * source[30] + c70 * source[32] + c69 * source[34]
                  - c67 * source[44] + c68 * source[39] + c68 * source[41]
                  - c69 * source[30] - c70 * source[32] - c69 * source[34]
                  + c67 * source[74] - c68 * source[69] - c68 * source[71]
                  + c69 * source[60] + c70 * source[62] + c69 * source[64];
    target[45] =  c12 * source[150] - c14 * source[152] + c12 * source[154]
                  - c13 * source[15] + c12 * source[17] - c13 * source[19]
                  - c13 * source[45] + c12 * source[47] - c13 * source[49];
    target[46] =  c32 * source[151] - c32 * source[153] - c33 * source[16]
                  + c33 * source[18] - c33 * source[46] + c33 * source[48];
    target[47] =  c17 * source[155] - c44 * source[157] - c15 * source[20]
                  + c35 * source[22] - c15 * source[50] + c35 * source[52];
    target[48] =  c44 * source[156] - c17 * source[158] - c35 * source[21]
                  + c15 * source[23] - c35 * source[51] + c15 * source[53];
    target[49] =  c56 * source[159] - c56 * source[161] - c57 * source[150]
                  + c57 * source[152] - c57 * source[152] + c57 * source[154]
                  - c57 * source[24] + c57 * source[26] + c58 * source[15]
                  - c58 * source[17] + c58 * source[17] - c58 * source[19]
                  - c57 * source[54] + c57 * source[56] + c58 * source[45]
                  - c58 * source[47] + c58 * source[47] - c58 * source[49];
    target[50] =  c71 * source[160] - c72 * source[151] - c72 * source[153]
                  - c72 * source[25] + c73 * source[16] + c73 * source[18]
                  - c72 * source[55] + c73 * source[46] + c73 * source[48];
    target[51] =  c74 * source[162] - c75 * source[155] - c75 * source[157]
                  - c76 * source[27] + c77 * source[20] + c77 * source[22]
                  - c76 * source[57] + c77 * source[50] + c77 * source[52];
    target[52] =  c74 * source[163] - c75 * source[156] - c75 * source[158]
                  - c76 * source[28] + c77 * source[21] + c77 * source[23]
                  - c76 * source[58] + c77 * source[51] + c77 * source[53];
    target[53] =  c78 * source[164] - c79 * source[159] - c79 * source[161]
                  + c66 * source[150] + c80 * source[152] + c66 * source[154]
                  - c81 * source[29] + c63 * source[24] + c63 * source[26]
                  - c70 * source[15] - c82 * source[17] - c70 * source[19]
                  - c81 * source[59] + c63 * source[54] + c63 * source[56]
                  - c70 * source[45] - c82 * source[47] - c70 * source[49];
    target[54] =  c15 * source[180] - c17 * source[182] + c15 * source[184]
                  - c16 * source[75] + c18 * source[77] - c16 * source[79]
                  - c16 * source[105] + c18 * source[107] - c16 * source[109];
    target[55] =  c34 * source[181] - c34 * source[183] - c35 * source[76]
                  + c35 * source[78] - c35 * source[106] + c35 * source[108];
    target[56] =  c45 * source[185] - c32 * source[187] - c12 * source[80]
                  + c11 * source[82] - c12 * source[110] + c11 * source[112];
    target[57] =  c32 * source[186] - c45 * source[188] - c11 * source[81]
                  + c12 * source[83] - c11 * source[111] + c12 * source[113];
    target[58] =  c59 * source[189] - c59 * source[191] - c61 * source[180]
                  + c61 * source[182] - c61 * source[182] + c61 * source[184]
                  - c60 * source[84] + c60 * source[86] + c62 * source[75]
                  - c62 * source[77] + c62 * source[77] - c62 * source[79]
                  - c60 * source[114] + c60 * source[116] + c62 * source[105]
                  - c62 * source[107] + c62 * source[107] - c62 * source[109];
    target[59] =  c74 * source[190] - c76 * source[181] - c76 * source[183]
                  - c75 * source[85] + c77 * source[76] + c77 * source[78]
                  - c75 * source[115] + c77 * source[106] + c77 * source[108];
    target[60] =  c83 * source[192] - c72 * source[185] - c72 * source[187]
                  - c72 * source[87] + c84 * source[80] + c84 * source[82]
                  - c72 * source[117] + c84 * source[110] + c84 * source[112];
    target[61] =  c83 * source[193] - c72 * source[186] - c72 * source[188]
                  - c72 * source[88] + c84 * source[81] + c84 * source[83]
                  - c72 * source[118] + c84 * source[111] + c84 * source[113];
    target[62] =  c85 * source[194] - c86 * source[189] - c86 * source[191]
                  + c87 * source[180] + c88 * source[182] + c87 * source[184]
                  - c88 * source[89] + c89 * source[84] + c89 * source[86]
                  - c90 * source[75] - c91 * source[77] - c90 * source[79]
                  - c88 * source[119] + c89 * source[114] + c89 * source[116]
                  - c90 * source[105] - c91 * source[107] - c90 * source[109];
    target[63] =  c15 * source[195] - c17 * source[197] + c15 * source[199]
                  - c16 * source[90] + c18 * source[92] - c16 * source[94]
                  - c16 * source[120] + c18 * source[122] - c16 * source[124];
    target[64] =  c34 * source[196] - c34 * source[198] - c35 * source[91]
                  + c35 * source[93] - c35 * source[121] + c35 * source[123];
    target[65] =  c45 * source[200] - c32 * source[202] - c12 * source[95]
                  + c11 * source[97] - c12 * source[125] + c11 * source[127];
    target[66] =  c32 * source[201] - c45 * source[203] - c11 * source[96]
                  + c12 * source[98] - c11 * source[126] + c12 * source[128];
    target[67] =  c59 * source[204] - c59 * source[206] - c61 * source[195]
                  + c61 * source[197] - c61 * source[197] + c61 * source[199]
                  - c60 * source[99] + c60 * source[101] + c62 * source[90]
                  - c62 * source[92] + c62 * source[92] - c62 * source[94]
                  - c60 * source[129] + c60 * source[131] + c62 * source[120]
                  - c62 * source[122] + c62 * source[122] - c62 * source[124];
    target[68] =  c74 * source[205] - c76 * source[196] - c76 * source[198]
                  - c75 * source[100] + c77 * source[91] + c77 * source[93]
                  - c75 * source[130] + c77 * source[121] + c77 * source[123];
    target[69] =  c83 * source[207] - c72 * source[200] - c72 * source[202]
                  - c72 * source[102] + c84 * source[95] + c84 * source[97]
                  - c72 * source[132] + c84 * source[125] + c84 * source[127];
    target[70] =  c83 * source[208] - c72 * source[201] - c72 * source[203]
                  - c72 * source[103] + c84 * source[96] + c84 * source[98]
                  - c72 * source[133] + c84 * source[126] + c84 * source[128];
    target[71] =  c85 * source[209] - c86 * source[204] - c86 * source[206]
                  + c87 * source[195] + c88 * source[197] + c87 * source[199]
                  - c88 * source[104] + c89 * source[99] + c89 * source[101]
                  - c90 * source[90] - c91 * source[92] - c90 * source[94]
                  - c88 * source[134] + c89 * source[129] + c89 * source[131]
                  - c90 * source[120] - c91 * source[122] - c90 * source[124];
    target[72] =  c19 * source[210] - c23 * source[212] + c19 * source[214]
                  - c20 * source[135] + c24 * source[137] - c20 * source[139]
                  - c20 * source[165] + c24 * source[167] - c20 * source[169]
                  + c21 * source[0] - c25 * source[2] + c21 * source[4]
                  + c22 * source[30] - c26 * source[32] + c22 * source[34]
                  + c21 * source[60] - c25 * source[62] + c21 * source[64];
    target[73] =  c36 * source[211] - c36 * source[213] - c37 * source[136]
                  + c37 * source[138] - c37 * source[166] + c37 * source[168]
                  + c38 * source[1] - c38 * source[3] + c20 * source[31]
                  - c20 * source[33] + c38 * source[61] - c38 * source[63];
    target[74] =  c46 * source[215] - c47 * source[217] - c47 * source[140]
                  + c50 * source[142] - c47 * source[170] + c50 * source[172]
                  + c48 * source[5] - c51 * source[7] + c49 * source[35]
                  - c52 * source[37] + c48 * source[65] - c51 * source[67];
    target[75] =  c47 * source[216] - c46 * source[218] - c50 * source[141]
                  + c47 * source[143] - c50 * source[171] + c47 * source[173]
                  + c51 * source[6] - c48 * source[8] + c52 * source[36]
                  - c49 * source[38] + c51 * source[66] - c48 * source[68];
    target[76] =  c63 * source[219] - c63 * source[221] - c67 * source[210]
                  + c67 * source[212] - c67 * source[212] + c67 * source[214]
                  - c64 * source[144] + c64 * source[146] + c68 * source[135]
                  - c68 * source[137] + c68 * source[137] - c68 * source[139]
                  - c64 * source[174] + c64 * source[176] + c68 * source[165]
                  - c68 * source[167] + c68 * source[167] - c68 * source[169]
                  + c65 * source[9] - c65 * source[11] - c69 * source[0]
                  + c69 * source[2] - c69 * source[2] + c69 * source[4]
                  + c66 * source[39] - c66 * source[41] - c70 * source[30]
                  + c70 * source[32] - c70 * source[32] + c70 * source[34]
                  + c65 * source[69] - c65 * source[71] - c69 * source[60]
                  + c69 * source[62] - c69 * source[62] + c69 * source[64];
    target[77] =  c78 * source[220] - c81 * source[211] - c81 * source[213]
                  - c79 * source[145] + c63 * source[136] + c63 * source[138]
                  - c79 * source[175] + c63 * source[166] + c63 * source[168]
                  + c66 * source[10] - c70 * source[1] - c70 * source[3]
                  + c80 * source[40] - c82 * source[31] - c82 * source[33]
                  + c66 * source[70] - c70 * source[61] - c70 * source[63];
    target[78] =  c85 * source[222] - c88 * source[215] - c88 * source[217]
                  - c86 * source[147] + c89 * source[140] + c89 * source[142]
                  - c86 * source[177] + c89 * source[170] + c89 * source[172]
                  + c87 * source[12] - c90 * source[5] - c90 * source[7]
                  + c88 * source[42] - c91 * source[35] - c91 * source[37]
                  + c87 * source[72] - c90 * source[65] - c90 * source[67];
    target[79] =  c85 * source[223] - c88 * source[216] - c88 * source[218]
                  - c86 * source[148] + c89 * source[141] + c89 * source[143]
                  - c86 * source[178] + c89 * source[171] + c89 * source[173]
                  + c87 * source[13] - c90 * source[6] - c90 * source[8]
                  + c88 * source[43] - c91 * source[36] - c91 * source[38]
                  + c87 * source[73] - c90 * source[66] - c90 * source[68];
    target[80] =  source[224] - c92 * source[219] - c92 * source[221]
                  + c93 * source[210] + c94 * source[212] + c93 * source[214]
                  - c92 * source[149] + c95 * source[144] + c95 * source[146]
                  - c96 * source[135] - c97 * source[137] - c96 * source[139]
                  - c92 * source[179] + c95 * source[174] + c95 * source[176]
                  - c96 * source[165] - c97 * source[167] - c96 * source[169]
                  + c93 * source[14] - c96 * source[9] - c96 * source[11]
                  + c98 * source[0] + c99 * source[2] + c98 * source[4]
                  + c94 * source[44] - c97 * source[39] - c97 * source[41]
                  + c99 * source[30] + c100 * source[32] + c99 * source[34]
                  + c93 * source[74] - c96 * source[69] - c96 * source[71]
                  + c98 * source[60] + c99 * source[62] + c98 * source[64];
  }
}

void CCarSphList::carsph_44(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c71 = 45;
  const double c44 = 42.093645601206838;
  const double c40 = 39.375;
  const double c14 = 29.764702249476645;
  const double c8 = 27.84232950922031;
  const double c56 = 22.5;
  const double c74 = 21.213203435596427;
  const double c42 = 21.046822800603419;
  const double c79 = 20.124611797498108;
  const double c32 = 19.843134832984429;
  const double c2 = 19.6875;
  const double c50 = 18.824850597016699;
  const double c29 = 18.561553006146873;
  const double c75 = 15.90990257669732;
  const double c11 = 14.882351124738323;
  const double c17 = 14.031215200402281;
  const double c24 = 13.311179511974137;
  const double c4 = 13.125;
  const double c53 = 11.25;
  const double c59 = 10.606601717798213;
  const double c18 = 10.523411400301709;
  const double c64 = 10.062305898749054;
  const double c83 = 10;
  const double c30 = 9.9215674164922145;
  const double c86 = 9.4868329805051381;
  const double c34 = 9.354143466934854;
  const double c7 = 9.2807765030734366;
  const double c95 = 9;
  const double c37 = 8.8741196746494246;
  const double c27 = 8.75;
  const double c60 = 7.9549512883486599;
  const double c72 = 7.5;
  const double c89 = 7.1151247353788536;
  const double c35 = 7.0156076002011405;
  const double c78 = 6.7082039324993694;
  const double c45 = 6.6143782776614763;
  const double c47 = 6.2749501990055663;
  const double c28 = 6.1871843353822911;
  const double c84 = 5.625;
  const double c80 = 5.0311529493745271;
  const double c12 = 4.9607837082461073;
  const double c52 = 4.7062126492541747;
  const double c6 = 4.6403882515367183;
  const double c23 = 4.4370598373247123;
  const double c39 = 4.375;
  const double c57 = 3.75;
  const double c76 = 3.5355339059327378;
  const double c43 = 3.5078038001005702;
  const double c63 = 3.3541019662496847;
  const double c26 = 3.3277948779935342;
  const double c33 = 3.3071891388307382;
  const double c1 = 3.28125;
  const double c85 = 3.1622776601683795;
  const double c92 = 3;
  const double c36 = 2.9580398915498081;
  const double c77 = 2.6516504294495533;
  const double c66 = 2.5155764746872635;
  const double c9 = 2.4803918541230536;
  const double c88 = 2.3717082451262845;
  const double c51 = 2.3531063246270874;
  const double c15 = 2.3385358667337135;
  const double c97 = 2.25;
  const double c20 = 2.2185299186623562;
  const double c3 = 2.1875;
  const double c46 = 2.0916500663351889;
  const double c54 = 1.875;
  const double c91 = 1.7787811838447134;
  const double c61 = 1.7677669529663689;
  const double c16 = 1.7539019000502851;
  const double c68 = 1.6770509831248424;
  const double c25 = 1.6638974389967671;
  const double c31 = 1.6535945694153691;
  const double c49 = 1.5687375497513916;
  const double c5 = 1.5467960838455728;
  const double c62 = 1.3258252147247767;
  const double c65 = 1.2577882373436318;
  const double c73 = 1.25;
  const double c87 = 1.1858541225631423;
  const double c41 = 1.1692679333668567;
  const double c96 = 1.125;
  const double c81 = 1.1180339887498949;
  const double c38 = 1.1092649593311781;
  const double c90 = 0.8893905919223567;
  const double c82 = 0.83852549156242118;
  const double c13 = 0.82679728470768454;
  const double c48 = 0.78436877487569578;
  const double c94 = 0.75;
  const double c19 = 0.73950997288745202;
  const double c58 = 0.625;
  const double c100 = 0.5625;
  const double c67 = 0.55901699437494745;
  const double c22 = 0.55463247966558904;
  const double c0 = 0.546875;
  const double c70 = 0.41926274578121059;
  const double c10 = 0.41339864235384227;
  const double c93 = 0.375;
  const double c55 = 0.3125;
  const double c99 = 0.28125;
  const double c21 = 0.27731623983279452;
  const double c69 = 0.20963137289060529;
  const double c98 = 0.140625;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 81, source += 225) {
    target[0] =  c0 * source[0] - c1 * source[2] + c0 * source[4]
                  - c1 * source[30] + c2 * source[32] - c1 * source[34]
                  + c0 * source[60] - c1 * source[62] + c0 * source[64];
    target[1] =  c3 * source[1] - c3 * source[3] - c4 * source[31]
                  + c4 * source[33] + c3 * source[61] - c3 * source[63];
    target[2] =  c5 * source[5] - c6 * source[7] - c7 * source[35]
                  + c8 * source[37] + c5 * source[65] - c6 * source[67];
    target[3] =  c6 * source[6] - c5 * source[8] - c8 * source[36]
                  + c7 * source[38] + c6 * source[66] - c5 * source[68];
    target[4] =  c9 * source[9] - c9 * source[11] - c10 * source[0]
                  + c10 * source[2] - c10 * source[2] + c10 * source[4]
                  - c11 * source[39] + c11 * source[41] + c9 * source[30]
                  - c9 * source[32] + c9 * source[32] - c9 * source[34]
                  + c9 * source[69] - c9 * source[71] - c10 * source[60]
                  + c10 * source[62] - c10 * source[62] + c10 * source[64];
    target[5] =  c12 * source[10] - c13 * source[1] - c13 * source[3]
                  - c14 * source[40] + c12 * source[31] + c12 * source[33]
                  + c12 * source[70] - c13 * source[61] - c13 * source[63];
    target[6] =  c15 * source[12] - c16 * source[5] - c16 * source[7]
                  - c17 * source[42] + c18 * source[35] + c18 * source[37]
                  + c15 * source[72] - c16 * source[65] - c16 * source[67];
    target[7] =  c15 * source[13] - c16 * source[6] - c16 * source[8]
                  - c17 * source[43] + c18 * source[36] + c18 * source[38]
                  + c15 * source[73] - c16 * source[66] - c16 * source[68];
    target[8] =  c19 * source[14] - c20 * source[9] - c20 * source[11]
                  + c21 * source[0] + c22 * source[2] + c21 * source[4]
                  - c23 * source[44] + c24 * source[39] + c24 * source[41]
                  - c25 * source[30] - c26 * source[32] - c25 * source[34]
                  + c19 * source[74] - c20 * source[69] - c20 * source[71]
                  + c21 * source[60] + c22 * source[62] + c21 * source[64];
    target[9] =  c3 * source[15] - c4 * source[17] + c3 * source[19]
                  - c3 * source[45] + c4 * source[47] - c3 * source[49];
    target[10] =  c27 * source[16] - c27 * source[18] - c27 * source[46]
                  + c27 * source[48];
    target[11] =  c28 * source[20] - c29 * source[22] - c28 * source[50]
                  + c29 * source[52];
    target[12] =  c29 * source[21] - c28 * source[23] - c29 * source[51]
                  + c28 * source[53];
    target[13] =  c30 * source[24] - c30 * source[26] - c31 * source[15]
                  + c31 * source[17] - c31 * source[17] + c31 * source[19]
                  - c30 * source[54] + c30 * source[56] + c31 * source[45]
                  - c31 * source[47] + c31 * source[47] - c31 * source[49];
    target[14] =  c32 * source[25] - c33 * source[16] - c33 * source[18]
                  - c32 * source[55] + c33 * source[46] + c33 * source[48];
    target[15] =  c34 * source[27] - c35 * source[20] - c35 * source[22]
                  - c34 * source[57] + c35 * source[50] + c35 * source[52];
    target[16] =  c34 * source[28] - c35 * source[21] - c35 * source[23]
                  - c34 * source[58] + c35 * source[51] + c35 * source[53];
    target[17] =  c36 * source[29] - c37 * source[24] - c37 * source[26]
                  + c38 * source[15] + c20 * source[17] + c38 * source[19]
                  - c36 * source[59] + c37 * source[54] + c37 * source[56]
                  - c38 * source[45] - c20 * source[47] - c38 * source[49];
    target[18] =  c5 * source[75] - c7 * source[77] + c5 * source[79]
                  - c6 * source[105] + c8 * source[107] - c6 * source[109];
    target[19] =  c28 * source[76] - c28 * source[78] - c29 * source[106]
                  + c29 * source[108];
    target[20] =  c39 * source[80] - c4 * source[82] - c4 * source[110]
                  + c40 * source[112];
    target[21] =  c4 * source[81] - c39 * source[83] - c40 * source[111]
                  + c4 * source[113];
    target[22] =  c35 * source[84] - c35 * source[86] - c41 * source[75]
                  + c41 * source[77] - c41 * source[77] + c41 * source[79]
                  - c42 * source[114] + c42 * source[116] + c43 * source[105]
                  - c43 * source[107] + c43 * source[107] - c43 * source[109];
    target[23] =  c17 * source[85] - c15 * source[76] - c15 * source[78]
                  - c44 * source[115] + c35 * source[106] + c35 * source[108];
    target[24] =  c45 * source[87] - c12 * source[80] - c12 * source[82]
                  - c32 * source[117] + c11 * source[110] + c11 * source[112];
    target[25] =  c45 * source[88] - c12 * source[81] - c12 * source[83]
                  - c32 * source[118] + c11 * source[111] + c11 * source[113];
    target[26] =  c46 * source[89] - c47 * source[84] - c47 * source[86]
                  + c48 * source[75] + c49 * source[77] + c48 * source[79]
                  - c47 * source[119] + c50 * source[114] + c50 * source[116]
                  - c51 * source[105] - c52 * source[107] - c51 * source[109];
    target[27] =  c6 * source[90] - c8 * source[92] + c6 * source[94]
                  - c5 * source[120] + c7 * source[122] - c5 * source[124];
    target[28] =  c29 * source[91] - c29 * source[93] - c28 * source[121]
                  + c28 * source[123];
    target[29] =  c4 * source[95] - c40 * source[97] - c39 * source[125]
                  + c4 * source[127];
    target[30] =  c40 * source[96] - c4 * source[98] - c4 * source[126]
                  + c39 * source[128];
    target[31] =  c42 * source[99] - c42 * source[101] - c43 * source[90]
                  + c43 * source[92] - c43 * source[92] + c43 * source[94]
                  - c35 * source[129] + c35 * source[131] + c41 * source[120]
                  - c41 * source[122] + c41 * source[122] - c41 * source[124];
    target[32] =  c44 * source[100] - c35 * source[91] - c35 * source[93]
                  - c17 * source[130] + c15 * source[121] + c15 * source[123];
    target[33] =  c32 * source[102] - c11 * source[95] - c11 * source[97]
                  - c45 * source[132] + c12 * source[125] + c12 * source[127];
    target[34] =  c32 * source[103] - c11 * source[96] - c11 * source[98]
                  - c45 * source[133] + c12 * source[126] + c12 * source[128];
    target[35] =  c47 * source[104] - c50 * source[99] - c50 * source[101]
                  + c51 * source[90] + c52 * source[92] + c51 * source[94]
                  - c46 * source[134] + c47 * source[129] + c47 * source[131]
                  - c48 * source[120] - c49 * source[122] - c48 * source[124];
    target[36] =  c9 * source[135] - c11 * source[137] + c9 * source[139]
                  - c9 * source[165] + c11 * source[167] - c9 * source[169]
                  - c10 * source[0] + c9 * source[2] - c10 * source[4]
                  + c10 * source[30] - c9 * source[32] + c10 * source[34]
                  - c10 * source[30] + c9 * source[32] - c10 * source[34]
                  + c10 * source[60] - c9 * source[62] + c10 * source[64];
    target[37] =  c30 * source[136] - c30 * source[138] - c30 * source[166]
                  + c30 * source[168] - c31 * source[1] + c31 * source[3]
                  + c31 * source[31] - c31 * source[33] - c31 * source[31]
                  + c31 * source[33] + c31 * source[61] - c31 * source[63];
    target[38] =  c35 * source[140] - c42 * source[142] - c35 * source[170]
                  + c42 * source[172] - c41 * source[5] + c43 * source[7]
                  + c41 * source[35] - c43 * source[37] - c41 * source[35]
                  + c43 * source[37] + c41 * source[65] - c43 * source[67];
    target[39] =  c42 * source[141] - c35 * source[143] - c42 * source[171]
                  + c35 * source[173] - c43 * source[6] + c41 * source[8]
                  + c43 * source[36] - c41 * source[38] - c43 * source[36]
                  + c41 * source[38] + c43 * source[66] - c41 * source[68];
    target[40] =  c53 * source[144] - c53 * source[146] - c54 * source[135]
                  + c54 * source[137] - c54 * source[137] + c54 * source[139]
                  - c53 * source[174] + c53 * source[176] + c54 * source[165]
                  - c54 * source[167] + c54 * source[167] - c54 * source[169]
                  - c54 * source[9] + c54 * source[11] + c55 * source[0]
                  - c55 * source[2] + c55 * source[2] - c55 * source[4]
                  + c54 * source[39] - c54 * source[41] - c55 * source[30]
                  + c55 * source[32] - c55 * source[32] + c55 * source[34]
                  - c54 * source[39] + c54 * source[41] + c55 * source[30]
                  - c55 * source[32] + c55 * source[32] - c55 * source[34]
                  + c54 * source[69] - c54 * source[71] - c55 * source[60]
                  + c55 * source[62] - c55 * source[62] + c55 * source[64];
    target[41] =  c56 * source[145] - c57 * source[136] - c57 * source[138]
                  - c56 * source[175] + c57 * source[166] + c57 * source[168]
                  - c57 * source[10] + c58 * source[1] + c58 * source[3]
                  + c57 * source[40] - c58 * source[31] - c58 * source[33]
                  - c57 * source[40] + c58 * source[31] + c58 * source[33]
                  + c57 * source[70] - c58 * source[61] - c58 * source[63];
    target[42] =  c59 * source[147] - c60 * source[140] - c60 * source[142]
                  - c59 * source[177] + c60 * source[170] + c60 * source[172]
                  - c61 * source[12] + c62 * source[5] + c62 * source[7]
                  + c61 * source[42] - c62 * source[35] - c62 * source[37]
                  - c61 * source[42] + c62 * source[35] + c62 * source[37]
                  + c61 * source[72] - c62 * source[65] - c62 * source[67];
    target[43] =  c59 * source[148] - c60 * source[141] - c60 * source[143]
                  - c59 * source[178] + c60 * source[171] + c60 * source[173]
                  - c61 * source[13] + c62 * source[6] + c62 * source[8]
                  + c61 * source[43] - c62 * source[36] - c62 * source[38]
                  - c61 * source[43] + c62 * source[36] + c62 * source[38]
                  + c61 * source[73] - c62 * source[66] - c62 * source[68];
    target[44] =  c63 * source[149] - c64 * source[144] - c64 * source[146]
                  + c65 * source[135] + c66 * source[137] + c65 * source[139]
                  - c63 * source[179] + c64 * source[174] + c64 * source[176]
                  - c65 * source[165] - c66 * source[167] - c65 * source[169]
                  - c67 * source[14] + c68 * source[9] + c68 * source[11]
                  - c69 * source[0] - c70 * source[2] - c69 * source[4]
                  + c67 * source[44] - c68 * source[39] - c68 * source[41]
                  + c69 * source[30] + c70 * source[32] + c69 * source[34]
                  - c67 * source[44] + c68 * source[39] + c68 * source[41]
                  - c69 * source[30] - c70 * source[32] - c69 * source[34]
                  + c67 * source[74] - c68 * source[69] - c68 * source[71]
                  + c69 * source[60] + c70 * source[62] + c69 * source[64];
    target[45] =  c12 * source[150] - c14 * source[152] + c12 * source[154]
                  - c13 * source[15] + c12 * source[17] - c13 * source[19]
                  - c13 * source[45] + c12 * source[47] - c13 * source[49];
    target[46] =  c32 * source[151] - c32 * source[153] - c33 * source[16]
                  + c33 * source[18] - c33 * source[46] + c33 * source[48];
    target[47] =  c17 * source[155] - c44 * source[157] - c15 * source[20]
                  + c35 * source[22] - c15 * source[50] + c35 * source[52];
    target[48] =  c44 * source[156] - c17 * source[158] - c35 * source[21]
                  + c15 * source[23] - c35 * source[51] + c15 * source[53];
    target[49] =  c56 * source[159] - c56 * source[161] - c57 * source[150]
                  + c57 * source[152] - c57 * source[152] + c57 * source[154]
                  - c57 * source[24] + c57 * source[26] + c58 * source[15]
                  - c58 * source[17] + c58 * source[17] - c58 * source[19]
                  - c57 * source[54] + c57 * source[56] + c58 * source[45]
                  - c58 * source[47] + c58 * source[47] - c58 * source[49];
    target[50] =  c71 * source[160] - c72 * source[151] - c72 * source[153]
                  - c72 * source[25] + c73 * source[16] + c73 * source[18]
                  - c72 * source[55] + c73 * source[46] + c73 * source[48];
    target[51] =  c74 * source[162] - c75 * source[155] - c75 * source[157]
                  - c76 * source[27] + c77 * source[20] + c77 * source[22]
                  - c76 * source[57] + c77 * source[50] + c77 * source[52];
    target[52] =  c74 * source[163] - c75 * source[156] - c75 * source[158]
                  - c76 * source[28] + c77 * source[21] + c77 * source[23]
                  - c76 * source[58] + c77 * source[51] + c77 * source[53];
    target[53] =  c78 * source[164] - c79 * source[159] - c79 * source[161]
                  + c66 * source[150] + c80 * source[152] + c66 * source[154]
                  - c81 * source[29] + c63 * source[24] + c63 * source[26]
                  - c70 * source[15] - c82 * source[17] - c70 * source[19]
                  - c81 * source[59] + c63 * source[54] + c63 * source[56]
                  - c70 * source[45] - c82 * source[47] - c70 * source[49];
    target[54] =  c15 * source[180] - c17 * source[182] + c15 * source[184]
                  - c16 * source[75] + c18 * source[77] - c16 * source[79]
                  - c16 * source[105] + c18 * source[107] - c16 * source[109];
    target[55] =  c34 * source[181] - c34 * source[183] - c35 * source[76]
                  + c35 * source[78] - c35 * source[106] + c35 * source[108];
    target[56] =  c45 * source[185] - c32 * source[187] - c12 * source[80]
                  + c11 * source[82] - c12 * source[110] + c11 * source[112];
    target[57] =  c32 * source[186] - c45 * source[188] - c11 * source[81]
                  + c12 * source[83] - c11 * source[111] + c12 * source[113];
    target[58] =  c59 * source[189] - c59 * source[191] - c61 * source[180]
                  + c61 * source[182] - c61 * source[182] + c61 * source[184]
                  - c60 * source[84] + c60 * source[86] + c62 * source[75]
                  - c62 * source[77] + c62 * source[77] - c62 * source[79]
                  - c60 * source[114] + c60 * source[116] + c62 * source[105]
                  - c62 * source[107] + c62 * source[107] - c62 * source[109];
    target[59] =  c74 * source[190] - c76 * source[181] - c76 * source[183]
                  - c75 * source[85] + c77 * source[76] + c77 * source[78]
                  - c75 * source[115] + c77 * source[106] + c77 * source[108];
    target[60] =  c83 * source[192] - c72 * source[185] - c72 * source[187]
                  - c72 * source[87] + c84 * source[80] + c84 * source[82]
                  - c72 * source[117] + c84 * source[110] + c84 * source[112];
    target[61] =  c83 * source[193] - c72 * source[186] - c72 * source[188]
                  - c72 * source[88] + c84 * source[81] + c84 * source[83]
                  - c72 * source[118] + c84 * source[111] + c84 * source[113];
    target[62] =  c85 * source[194] - c86 * source[189] - c86 * source[191]
                  + c87 * source[180] + c88 * source[182] + c87 * source[184]
                  - c88 * source[89] + c89 * source[84] + c89 * source[86]
                  - c90 * source[75] - c91 * source[77] - c90 * source[79]
                  - c88 * source[119] + c89 * source[114] + c89 * source[116]
                  - c90 * source[105] - c91 * source[107] - c90 * source[109];
    target[63] =  c15 * source[195] - c17 * source[197] + c15 * source[199]
                  - c16 * source[90] + c18 * source[92] - c16 * source[94]
                  - c16 * source[120] + c18 * source[122] - c16 * source[124];
    target[64] =  c34 * source[196] - c34 * source[198] - c35 * source[91]
                  + c35 * source[93] - c35 * source[121] + c35 * source[123];
    target[65] =  c45 * source[200] - c32 * source[202] - c12 * source[95]
                  + c11 * source[97] - c12 * source[125] + c11 * source[127];
    target[66] =  c32 * source[201] - c45 * source[203] - c11 * source[96]
                  + c12 * source[98] - c11 * source[126] + c12 * source[128];
    target[67] =  c59 * source[204] - c59 * source[206] - c61 * source[195]
                  + c61 * source[197] - c61 * source[197] + c61 * source[199]
                  - c60 * source[99] + c60 * source[101] + c62 * source[90]
                  - c62 * source[92] + c62 * source[92] - c62 * source[94]
                  - c60 * source[129] + c60 * source[131] + c62 * source[120]
                  - c62 * source[122] + c62 * source[122] - c62 * source[124];
    target[68] =  c74 * source[205] - c76 * source[196] - c76 * source[198]
                  - c75 * source[100] + c77 * source[91] + c77 * source[93]
                  - c75 * source[130] + c77 * source[121] + c77 * source[123];
    target[69] =  c83 * source[207] - c72 * source[200] - c72 * source[202]
                  - c72 * source[102] + c84 * source[95] + c84 * source[97]
                  - c72 * source[132] + c84 * source[125] + c84 * source[127];
    target[70] =  c83 * source[208] - c72 * source[201] - c72 * source[203]
                  - c72 * source[103] + c84 * source[96] + c84 * source[98]
                  - c72 * source[133] + c84 * source[126] + c84 * source[128];
    target[71] =  c85 * source[209] - c86 * source[204] - c86 * source[206]
                  + c87 * source[195] + c88 * source[197] + c87 * source[199]
                  - c88 * source[104] + c89 * source[99] + c89 * source[101]
                  - c90 * source[90] - c91 * source[92] - c90 * source[94]
                  - c88 * source[134] + c89 * source[129] + c89 * source[131]
                  - c90 * source[120] - c91 * source[122] - c90 * source[124];
    target[72] =  c19 * source[210] - c23 * source[212] + c19 * source[214]
                  - c20 * source[135] + c24 * source[137] - c20 * source[139]
                  - c20 * source[165] + c24 * source[167] - c20 * source[169]
                  + c21 * source[0] - c25 * source[2] + c21 * source[4]
                  + c22 * source[30] - c26 * source[32] + c22 * source[34]
                  + c21 * source[60] - c25 * source[62] + c21 * source[64];
    target[73] =  c36 * source[211] - c36 * source[213] - c37 * source[136]
                  + c37 * source[138] - c37 * source[166] + c37 * source[168]
                  + c38 * source[1] - c38 * source[3] + c20 * source[31]
                  - c20 * source[33] + c38 * source[61] - c38 * source[63];
    target[74] =  c46 * source[215] - c47 * source[217] - c47 * source[140]
                  + c50 * source[142] - c47 * source[170] + c50 * source[172]
                  + c48 * source[5] - c51 * source[7] + c49 * source[35]
                  - c52 * source[37] + c48 * source[65] - c51 * source[67];
    target[75] =  c47 * source[216] - c46 * source[218] - c50 * source[141]
                  + c47 * source[143] - c50 * source[171] + c47 * source[173]
                  + c51 * source[6] - c48 * source[8] + c52 * source[36]
                  - c49 * source[38] + c51 * source[66] - c48 * source[68];
    target[76] =  c63 * source[219] - c63 * source[221] - c67 * source[210]
                  + c67 * source[212] - c67 * source[212] + c67 * source[214]
                  - c64 * source[144] + c64 * source[146] + c68 * source[135]
                  - c68 * source[137] + c68 * source[137] - c68 * source[139]
                  - c64 * source[174] + c64 * source[176] + c68 * source[165]
                  - c68 * source[167] + c68 * source[167] - c68 * source[169]
                  + c65 * source[9] - c65 * source[11] - c69 * source[0]
                  + c69 * source[2] - c69 * source[2] + c69 * source[4]
                  + c66 * source[39] - c66 * source[41] - c70 * source[30]
                  + c70 * source[32] - c70 * source[32] + c70 * source[34]
                  + c65 * source[69] - c65 * source[71] - c69 * source[60]
                  + c69 * source[62] - c69 * source[62] + c69 * source[64];
    target[77] =  c78 * source[220] - c81 * source[211] - c81 * source[213]
                  - c79 * source[145] + c63 * source[136] + c63 * source[138]
                  - c79 * source[175] + c63 * source[166] + c63 * source[168]
                  + c66 * source[10] - c70 * source[1] - c70 * source[3]
                  + c80 * source[40] - c82 * source[31] - c82 * source[33]
                  + c66 * source[70] - c70 * source[61] - c70 * source[63];
    target[78] =  c85 * source[222] - c88 * source[215] - c88 * source[217]
                  - c86 * source[147] + c89 * source[140] + c89 * source[142]
                  - c86 * source[177] + c89 * source[170] + c89 * source[172]
                  + c87 * source[12] - c90 * source[5] - c90 * source[7]
                  + c88 * source[42] - c91 * source[35] - c91 * source[37]
                  + c87 * source[72] - c90 * source[65] - c90 * source[67];
    target[79] =  c85 * source[223] - c88 * source[216] - c88 * source[218]
                  - c86 * source[148] + c89 * source[141] + c89 * source[143]
                  - c86 * source[178] + c89 * source[171] + c89 * source[173]
                  + c87 * source[13] - c90 * source[6] - c90 * source[8]
                  + c88 * source[43] - c91 * source[36] - c91 * source[38]
                  + c87 * source[73] - c90 * source[66] - c90 * source[68];
    target[80] =  source[224] - c92 * source[219] - c92 * source[221]
                  + c93 * source[210] + c94 * source[212] + c93 * source[214]
                  - c92 * source[149] + c95 * source[144] + c95 * source[146]
                  - c96 * source[135] - c97 * source[137] - c96 * source[139]
                  - c92 * source[179] + c95 * source[174] + c95 * source[176]
                  - c96 * source[165] - c97 * source[167] - c96 * source[169]
                  + c93 * source[14] - c96 * source[9] - c96 * source[11]
                  + c98 * source[0] + c99 * source[2] + c98 * source[4]
                  + c94 * source[44] - c97 * source[39] - c97 * source[41]
                  + c99 * source[30] + c100 * source[32] + c99 * source[34]
                  + c93 * source[74] - c96 * source[69] - c96 * source[71]
                  + c98 * source[60] + c99 * source[62] + c98 * source[64];
  }
}

