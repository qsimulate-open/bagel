//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_53.cc
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


void CarSphList::carsph_53(const int nloop, const double* source, double* target) {
  const double c26 = 51.553976568253198;
  const double c49 = 48.605555238058955;
  const double c74 = 39.686269665968858;
  const double c35 = 34.369317712168801;
  const double c28 = 32.605597678926237;
  const double c24 = 31.57023420090513;
  const double c54 = 30.740852297878796;
  const double c40 = 29.764702249476645;
  const double c10 = 27.171331399105199;
  const double c25 = 25.776988284126599;
  const double c75 = 25.099800796022265;
  const double c45 = 24.302777619029477;
  const double c90 = 22.5;
  const double c36 = 21.737065119284157;
  const double c34 = 21.046822800603419;
  const double c33 = 19.966769267961205;
  const double c71 = 19.843134832984429;
  const double c117 = 19.364916731037084;
  const double c63 = 18.824850597016699;
  const double c13 = 17.1846588560844;
  const double c3 = 16.638974389967672;
  const double c48 = 16.201851746019649;
  const double c76 = 15.370426148939398;
  const double c89 = 15;
  const double c118 = 14.523687548277813;
  const double c94 = 14.230249470757707;
  const double c80 = 13.778379803155376;
  const double c7 = 13.5856656995526;
  const double c32 = 13.311179511974137;
  const double c62 = 12.549900398011133;
  const double c121 = 12.24744871391589;
  const double c69 = 12.151388809514739;
  const double c108 = 11.858541225631422;
  const double c86 = 11.25;
  const double c20 = 10.523411400301709;
  const double c52 = 10.246950765959598;
  const double c39 = 9.9215674164922145;
  const double c114 = 9.6824583655185421;
  const double c92 = 9.4868329805051381;
  const double c78 = 9.1855865354369186;
  const double c112 = 8.8939059192235668;
  const double c37 = 8.8741196746494246;
  const double c101 = 8.7142125289666872;
  const double c15 = 8.5923294280422002;
  const double c5 = 8.3194871949838358;
  const double c29 = 8.1513994197315593;
  const double c44 = 8.1009258730098246;
  const double c55 = 7.6852130744696989;
  const double c85 = 7.5;
  const double c116 = 7.2618437741389066;
  const double c19 = 7.0156076002011405;
  const double c8 = 6.7928328497762998;
  const double c61 = 6.2749501990055663;
  const double c51 = 6.0756944047573693;
  const double c100 = 5.8094750193111251;
  const double c125 = 5.625;
  const double c2 = 5.54632479665589;
  const double c27 = 5.4342662798210393;
  const double c22 = 5.2617057001508547;
  const double c73 = 5.123475382979799;
  const double c123 = 5;
  const double c70 = 4.9607837082461073;
  const double c79 = 4.5927932677184593;
  const double c110 = 4.4469529596117834;
  const double c14 = 4.2961647140211001;
  const double c60 = 4.1833001326703778;
  const double c68 = 4.0504629365049123;
  const double c107 = 3.9528470752104741;
  const double c99 = 3.872983346207417;
  const double c58 = 3.8426065372348495;
  const double c91 = 3.75;
  const double c43 = 3.7205877811845807;
  const double c115 = 3.6309218870694533;
  const double c95 = 3.5575623676894268;
  const double c21 = 3.5078038001005702;
  const double c31 = 3.3277948779935342;
  const double c38 = 3.3071891388307382;
  const double c72 = 3.1374750995027831;
  const double c77 = 3.0618621784789726;
  const double c47 = 3.0378472023786847;
  const double c111 = 2.9646353064078554;
  const double c124 = 2.8125;
  const double c4 = 2.773162398327945;
  const double c9 = 2.7171331399105196;
  const double c53 = 2.5617376914898995;
  const double c119 = 2.4494897427831779;
  const double c93 = 2.3717082451262845;
  const double c67 = 2.3531063246270874;
  const double c84 = 2.2963966338592297;
  const double c30 = 2.2185299186623562;
  const double c16 = 2.1480823570105501;
  const double c50 = 2.0252314682524561;
  const double c113 = 1.9364916731037085;
  const double c88 = 1.875;
  const double c23 = 1.7539019000502851;
  const double c11 = 1.71846588560844;
  const double c1 = 1.6638974389967671;
  const double c66 = 1.5687375497513916;
  const double c122 = 1.5;
  const double c109 = 1.4823176532039277;
  const double c105 = 1.4523687548277813;
  const double c6 = 1.3585665699552598;
  const double c56 = 1.2808688457449497;
  const double c42 = 1.2401959270615268;
  const double c96 = 1.1858541225631423;
  const double c82 = 1.1481983169296148;
  const double c18 = 1.052341140030171;
  const double c46 = 1.0126157341262281;
  const double c104 = 0.96824583655185426;
  const double c59 = 0.96065163430871237;
  const double c87 = 0.9375;
  const double c106 = 0.79056941504209488;
  const double c65 = 0.78436877487569578;
  const double c83 = 0.76546554461974314;
  const double c103 = 0.72618437741389064;
  const double c17 = 0.70156076002011403;
  const double c120 = 0.61237243569579447;
  const double c98 = 0.59292706128157113;
  const double c0 = 0.55463247966558904;
  const double c64 = 0.52291251658379723;
  const double c102 = 0.48412291827592713;
  const double c12 = 0.42961647140211001;
  const double c41 = 0.41339864235384227;
  const double c81 = 0.38273277230987157;
  const double c57 = 0.32021721143623744;
  const double c97 = 0.29646353064078557;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 77, source += 210) {
    target[0] =  c0 * source[0] - c1 * source[2] - c2 * source[20]
                  + c3 * source[22] + c4 * source[40] - c5 * source[42];
    target[1] =  c1 * source[1] - c0 * source[3] - c3 * source[21]
                  + c2 * source[23] + c5 * source[41] - c4 * source[43];
    target[2] =  c6 * source[4] - c6 * source[6] - c7 * source[24]
                  + c7 * source[26] + c8 * source[44] - c8 * source[46];
    target[3] =  c9 * source[5] - c10 * source[25] + c7 * source[45];
    target[4] =  c11 * source[7] - c12 * source[0] - c12 * source[2]
                  - c13 * source[27] + c14 * source[20] + c14 * source[22]
                  + c15 * source[47] - c16 * source[40] - c16 * source[42];
    target[5] =  c11 * source[8] - c12 * source[1] - c12 * source[3]
                  - c13 * source[28] + c14 * source[21] + c14 * source[23]
                  + c15 * source[48] - c16 * source[41] - c16 * source[43];
    target[6] =  c17 * source[9] - c18 * source[4] - c18 * source[6]
                  - c19 * source[29] + c20 * source[24] + c20 * source[26]
                  + c21 * source[49] - c22 * source[44] - c22 * source[46];
    target[7] =  c4 * source[10] - c5 * source[12] - c2 * source[30]
                  + c3 * source[32] + c0 * source[50] - c1 * source[52];
    target[8] =  c5 * source[11] - c4 * source[13] - c3 * source[31]
                  + c2 * source[33] + c1 * source[51] - c0 * source[53];
    target[9] =  c8 * source[14] - c8 * source[16] - c7 * source[34]
                  + c7 * source[36] + c6 * source[54] - c6 * source[56];
    target[10] =  c7 * source[15] - c10 * source[35] + c9 * source[55];
    target[11] =  c15 * source[17] - c16 * source[10] - c16 * source[12]
                  - c13 * source[37] + c14 * source[30] + c14 * source[32]
                  + c11 * source[57] - c12 * source[50] - c12 * source[52];
    target[12] =  c15 * source[18] - c16 * source[11] - c16 * source[13]
                  - c13 * source[38] + c14 * source[31] + c14 * source[33]
                  + c11 * source[58] - c12 * source[51] - c12 * source[53];
    target[13] =  c21 * source[19] - c22 * source[14] - c22 * source[16]
                  - c19 * source[39] + c20 * source[34] + c20 * source[36]
                  + c17 * source[59] - c18 * source[54] - c18 * source[56];
    target[14] =  c23 * source[60] - c22 * source[62] - c20 * source[80]
                  + c24 * source[82] + c23 * source[100] - c22 * source[102];
    target[15] =  c22 * source[61] - c23 * source[63] - c24 * source[81]
                  + c20 * source[83] + c22 * source[101] - c23 * source[103];
    target[16] =  c14 * source[64] - c14 * source[66] - c25 * source[84]
                  + c25 * source[86] + c14 * source[104] - c14 * source[106];
    target[17] =  c15 * source[65] - c26 * source[85] + c15 * source[105];
    target[18] =  c27 * source[67] - c6 * source[60] - c6 * source[62]
                  - c28 * source[87] + c29 * source[80] + c29 * source[82]
                  + c27 * source[107] - c6 * source[100] - c6 * source[102];
    target[19] =  c27 * source[68] - c6 * source[61] - c6 * source[63]
                  - c28 * source[88] + c29 * source[81] + c29 * source[83]
                  + c27 * source[108] - c6 * source[101] - c6 * source[103];
    target[20] =  c30 * source[69] - c31 * source[64] - c31 * source[66]
                  - c32 * source[89] + c33 * source[84] + c33 * source[86]
                  + c30 * source[109] - c31 * source[104] - c31 * source[106];
    target[21] =  c19 * source[70] - c34 * source[72] - c19 * source[90]
                  + c34 * source[92];
    target[22] =  c34 * source[71] - c19 * source[73] - c34 * source[91]
                  + c19 * source[93];
    target[23] =  c13 * source[74] - c13 * source[76] - c13 * source[94]
                  + c13 * source[96];
    target[24] =  c35 * source[75] - c35 * source[95];
    target[25] =  c36 * source[77] - c27 * source[70] - c27 * source[72]
                  - c36 * source[97] + c27 * source[90] + c27 * source[92];
    target[26] =  c36 * source[78] - c27 * source[71] - c27 * source[73]
                  - c36 * source[98] + c27 * source[91] + c27 * source[93];
    target[27] =  c37 * source[79] - c32 * source[74] - c32 * source[76]
                  - c37 * source[99] + c32 * source[94] + c32 * source[96];
    target[28] =  c38 * source[110] - c39 * source[112] - c39 * source[130]
                  + c40 * source[132] - c41 * source[0] + c42 * source[2]
                  + c42 * source[20] - c43 * source[22] - c41 * source[20]
                  + c42 * source[22] + c42 * source[40] - c43 * source[42];
    target[29] =  c39 * source[111] - c38 * source[113] - c40 * source[131]
                  + c39 * source[133] - c42 * source[1] + c41 * source[3]
                  + c43 * source[21] - c42 * source[23] - c42 * source[21]
                  + c41 * source[23] + c43 * source[41] - c42 * source[43];
    target[30] =  c44 * source[114] - c44 * source[116] - c45 * source[134]
                  + c45 * source[136] - c46 * source[4] + c46 * source[6]
                  + c47 * source[24] - c47 * source[26] - c46 * source[24]
                  + c46 * source[26] + c47 * source[44] - c47 * source[46];
    target[31] =  c48 * source[115] - c49 * source[135] - c50 * source[5]
                  + c51 * source[25] - c50 * source[25] + c51 * source[45];
    target[32] =  c52 * source[117] - c53 * source[110] - c53 * source[112]
                  - c54 * source[137] + c55 * source[130] + c55 * source[132]
                  - c56 * source[7] + c57 * source[0] + c57 * source[2]
                  + c58 * source[27] - c59 * source[20] - c59 * source[22]
                  - c56 * source[27] + c57 * source[20] + c57 * source[22]
                  + c58 * source[47] - c59 * source[40] - c59 * source[42];
    target[33] =  c52 * source[118] - c53 * source[111] - c53 * source[113]
                  - c54 * source[138] + c55 * source[131] + c55 * source[133]
                  - c56 * source[8] + c57 * source[1] + c57 * source[3]
                  + c58 * source[28] - c59 * source[21] - c59 * source[23]
                  - c56 * source[28] + c57 * source[21] + c57 * source[23]
                  + c58 * source[48] - c59 * source[41] - c59 * source[43];
    target[34] =  c60 * source[119] - c61 * source[114] - c61 * source[116]
                  - c62 * source[139] + c63 * source[134] + c63 * source[136]
                  - c64 * source[9] + c65 * source[4] + c65 * source[6]
                  + c66 * source[29] - c67 * source[24] - c67 * source[26]
                  - c64 * source[29] + c65 * source[24] + c65 * source[26]
                  + c66 * source[49] - c67 * source[44] - c67 * source[46];
    target[35] =  c39 * source[120] - c40 * source[122] - c38 * source[140]
                  + c39 * source[142] - c42 * source[10] + c43 * source[12]
                  + c41 * source[30] - c42 * source[32] - c42 * source[30]
                  + c43 * source[32] + c41 * source[50] - c42 * source[52];
    target[36] =  c40 * source[121] - c39 * source[123] - c39 * source[141]
                  + c38 * source[143] - c43 * source[11] + c42 * source[13]
                  + c42 * source[31] - c41 * source[33] - c43 * source[31]
                  + c42 * source[33] + c42 * source[51] - c41 * source[53];
    target[37] =  c45 * source[124] - c45 * source[126] - c44 * source[144]
                  + c44 * source[146] - c47 * source[14] + c47 * source[16]
                  + c46 * source[34] - c46 * source[36] - c47 * source[34]
                  + c47 * source[36] + c46 * source[54] - c46 * source[56];
    target[38] =  c49 * source[125] - c48 * source[145] - c51 * source[15]
                  + c50 * source[35] - c51 * source[35] + c50 * source[55];
    target[39] =  c54 * source[127] - c55 * source[120] - c55 * source[122]
                  - c52 * source[147] + c53 * source[140] + c53 * source[142]
                  - c58 * source[17] + c59 * source[10] + c59 * source[12]
                  + c56 * source[37] - c57 * source[30] - c57 * source[32]
                  - c58 * source[37] + c59 * source[30] + c59 * source[32]
                  + c56 * source[57] - c57 * source[50] - c57 * source[52];
    target[40] =  c54 * source[128] - c55 * source[121] - c55 * source[123]
                  - c52 * source[148] + c53 * source[141] + c53 * source[143]
                  - c58 * source[18] + c59 * source[11] + c59 * source[13]
                  + c56 * source[38] - c57 * source[31] - c57 * source[33]
                  - c58 * source[38] + c59 * source[31] + c59 * source[33]
                  + c56 * source[58] - c57 * source[51] - c57 * source[53];
    target[41] =  c62 * source[129] - c63 * source[124] - c63 * source[126]
                  - c60 * source[149] + c61 * source[144] + c61 * source[146]
                  - c66 * source[19] + c67 * source[14] + c67 * source[16]
                  + c64 * source[39] - c65 * source[34] - c65 * source[36]
                  - c66 * source[39] + c67 * source[34] + c67 * source[36]
                  + c64 * source[59] - c65 * source[54] - c65 * source[56];
    target[42] =  c68 * source[150] - c69 * source[152] - c68 * source[170]
                  + c69 * source[172] - c50 * source[60] + c51 * source[62]
                  + c50 * source[80] - c51 * source[82] - c50 * source[80]
                  + c51 * source[82] + c50 * source[100] - c51 * source[102];
    target[43] =  c69 * source[151] - c68 * source[153] - c69 * source[171]
                  + c68 * source[173] - c51 * source[61] + c50 * source[63]
                  + c51 * source[81] - c50 * source[83] - c51 * source[81]
                  + c50 * source[83] + c51 * source[101] - c50 * source[103];
    target[44] =  c39 * source[154] - c39 * source[156] - c39 * source[174]
                  + c39 * source[176] - c70 * source[64] + c70 * source[66]
                  + c70 * source[84] - c70 * source[86] - c70 * source[84]
                  + c70 * source[86] + c70 * source[104] - c70 * source[106];
    target[45] =  c71 * source[155] - c71 * source[175] - c39 * source[65]
                  + c39 * source[85] - c39 * source[85] + c39 * source[105];
    target[46] =  c62 * source[157] - c72 * source[150] - c72 * source[152]
                  - c62 * source[177] + c72 * source[170] + c72 * source[172]
                  - c61 * source[67] + c66 * source[60] + c66 * source[62]
                  + c61 * source[87] - c66 * source[80] - c66 * source[82]
                  - c61 * source[87] + c66 * source[80] + c66 * source[82]
                  + c61 * source[107] - c66 * source[100] - c66 * source[102];
    target[47] =  c62 * source[158] - c72 * source[151] - c72 * source[153]
                  - c62 * source[178] + c72 * source[171] + c72 * source[173]
                  - c61 * source[68] + c66 * source[61] + c66 * source[63]
                  + c61 * source[88] - c66 * source[81] - c66 * source[83]
                  - c61 * source[88] + c66 * source[81] + c66 * source[83]
                  + c61 * source[108] - c66 * source[101] - c66 * source[103];
    target[48] =  c73 * source[159] - c55 * source[154] - c55 * source[156]
                  - c73 * source[179] + c55 * source[174] + c55 * source[176]
                  - c53 * source[69] + c58 * source[64] + c58 * source[66]
                  + c53 * source[89] - c58 * source[84] - c58 * source[86]
                  - c53 * source[89] + c58 * source[84] + c58 * source[86]
                  + c53 * source[109] - c58 * source[104] - c58 * source[106];
    target[49] =  c44 * source[160] - c45 * source[162] - c68 * source[70]
                  + c69 * source[72] - c68 * source[90] + c69 * source[92];
    target[50] =  c45 * source[161] - c44 * source[163] - c69 * source[71]
                  + c68 * source[73] - c69 * source[91] + c68 * source[93];
    target[51] =  c71 * source[164] - c71 * source[166] - c39 * source[74]
                  + c39 * source[76] - c39 * source[94] + c39 * source[96];
    target[52] =  c74 * source[165] - c71 * source[75] - c71 * source[95];
    target[53] =  c75 * source[167] - c61 * source[160] - c61 * source[162]
                  - c62 * source[77] + c72 * source[70] + c72 * source[72]
                  - c62 * source[97] + c72 * source[90] + c72 * source[92];
    target[54] =  c75 * source[168] - c61 * source[161] - c61 * source[163]
                  - c62 * source[78] + c72 * source[71] + c72 * source[73]
                  - c62 * source[98] + c72 * source[91] + c72 * source[93];
    target[55] =  c52 * source[169] - c76 * source[164] - c76 * source[166]
                  - c73 * source[79] + c55 * source[74] + c55 * source[76]
                  - c73 * source[99] + c55 * source[94] + c55 * source[96];
    target[56] =  c77 * source[180] - c78 * source[182] - c79 * source[110]
                  + c80 * source[112] - c79 * source[130] + c80 * source[132]
                  + c81 * source[0] - c82 * source[2] + c83 * source[20]
                  - c84 * source[22] + c81 * source[40] - c82 * source[42];
    target[57] =  c78 * source[181] - c77 * source[183] - c80 * source[111]
                  + c79 * source[113] - c80 * source[131] + c79 * source[133]
                  + c82 * source[1] - c81 * source[3] + c84 * source[21]
                  - c83 * source[23] + c82 * source[41] - c81 * source[43];
    target[58] =  c85 * source[184] - c85 * source[186] - c86 * source[114]
                  + c86 * source[116] - c86 * source[134] + c86 * source[136]
                  + c87 * source[4] - c87 * source[6] + c88 * source[24]
                  - c88 * source[26] + c87 * source[44] - c87 * source[46];
    target[59] =  c89 * source[185] - c90 * source[115] - c90 * source[135]
                  + c88 * source[5] + c91 * source[25] + c88 * source[45];
    target[60] =  c92 * source[187] - c93 * source[180] - c93 * source[182]
                  - c94 * source[117] + c95 * source[110] + c95 * source[112]
                  - c94 * source[137] + c95 * source[130] + c95 * source[132]
                  + c96 * source[7] - c97 * source[0] - c97 * source[2]
                  + c93 * source[27] - c98 * source[20] - c98 * source[22]
                  + c96 * source[47] - c97 * source[40] - c97 * source[42];
    target[61] =  c92 * source[188] - c93 * source[181] - c93 * source[183]
                  - c94 * source[118] + c95 * source[111] + c95 * source[113]
                  - c94 * source[138] + c95 * source[131] + c95 * source[133]
                  + c96 * source[8] - c97 * source[1] - c97 * source[3]
                  + c93 * source[28] - c98 * source[21] - c98 * source[23]
                  + c96 * source[48] - c97 * source[41] - c97 * source[43];
    target[62] =  c99 * source[189] - c100 * source[184] - c100 * source[186]
                  - c100 * source[119] + c101 * source[114] + c101 * source[116]
                  - c100 * source[139] + c101 * source[134] + c101 * source[136]
                  + c102 * source[9] - c103 * source[4] - c103 * source[6]
                  + c104 * source[29] - c105 * source[24] - c105 * source[26]
                  + c102 * source[49] - c103 * source[44] - c103 * source[46];
    target[63] =  c77 * source[190] - c78 * source[192] - c79 * source[120]
                  + c80 * source[122] - c79 * source[140] + c80 * source[142]
                  + c81 * source[10] - c82 * source[12] + c83 * source[30]
                  - c84 * source[32] + c81 * source[50] - c82 * source[52];
    target[64] =  c78 * source[191] - c77 * source[193] - c80 * source[121]
                  + c79 * source[123] - c80 * source[141] + c79 * source[143]
                  + c82 * source[11] - c81 * source[13] + c84 * source[31]
                  - c83 * source[33] + c82 * source[51] - c81 * source[53];
    target[65] =  c85 * source[194] - c85 * source[196] - c86 * source[124]
                  + c86 * source[126] - c86 * source[144] + c86 * source[146]
                  + c87 * source[14] - c87 * source[16] + c88 * source[34]
                  - c88 * source[36] + c87 * source[54] - c87 * source[56];
    target[66] =  c89 * source[195] - c90 * source[125] - c90 * source[145]
                  + c88 * source[15] + c91 * source[35] + c88 * source[55];
    target[67] =  c92 * source[197] - c93 * source[190] - c93 * source[192]
                  - c94 * source[127] + c95 * source[120] + c95 * source[122]
                  - c94 * source[147] + c95 * source[140] + c95 * source[142]
                  + c96 * source[17] - c97 * source[10] - c97 * source[12]
                  + c93 * source[37] - c98 * source[30] - c98 * source[32]
                  + c96 * source[57] - c97 * source[50] - c97 * source[52];
    target[68] =  c92 * source[198] - c93 * source[191] - c93 * source[193]
                  - c94 * source[128] + c95 * source[121] + c95 * source[123]
                  - c94 * source[148] + c95 * source[141] + c95 * source[143]
                  + c96 * source[18] - c97 * source[11] - c97 * source[13]
                  + c93 * source[38] - c98 * source[31] - c98 * source[33]
                  + c96 * source[58] - c97 * source[51] - c97 * source[53];
    target[69] =  c99 * source[199] - c100 * source[194] - c100 * source[196]
                  - c100 * source[129] + c101 * source[124] + c101 * source[126]
                  - c100 * source[149] + c101 * source[144] + c101 * source[146]
                  + c102 * source[19] - c103 * source[14] - c103 * source[16]
                  + c104 * source[39] - c105 * source[34] - c105 * source[36]
                  + c102 * source[59] - c103 * source[54] - c103 * source[56];
    target[70] =  c106 * source[200] - c93 * source[202] - c107 * source[150]
                  + c108 * source[152] - c107 * source[170] + c108 * source[172]
                  + c109 * source[60] - c110 * source[62] + c111 * source[80]
                  - c112 * source[82] + c109 * source[100] - c110 * source[102];
    target[71] =  c93 * source[201] - c106 * source[203] - c108 * source[151]
                  + c107 * source[153] - c108 * source[171] + c107 * source[173]
                  + c110 * source[61] - c109 * source[63] + c112 * source[81]
                  - c111 * source[83] + c110 * source[101] - c109 * source[103];
    target[72] =  c113 * source[204] - c113 * source[206] - c114 * source[154]
                  + c114 * source[156] - c114 * source[174] + c114 * source[176]
                  + c115 * source[64] - c115 * source[66] + c116 * source[84]
                  - c116 * source[86] + c115 * source[104] - c115 * source[106];
    target[73] =  c99 * source[205] - c117 * source[155] - c117 * source[175]
                  + c116 * source[65] + c118 * source[85] + c116 * source[105];
    target[74] =  c119 * source[207] - c120 * source[200] - c120 * source[202]
                  - c121 * source[157] + c77 * source[150] + c77 * source[152]
                  - c121 * source[177] + c77 * source[170] + c77 * source[172]
                  + c79 * source[67] - c82 * source[60] - c82 * source[62]
                  + c78 * source[87] - c84 * source[80] - c84 * source[82]
                  + c79 * source[107] - c82 * source[100] - c82 * source[102];
    target[75] =  c119 * source[208] - c120 * source[201] - c120 * source[203]
                  - c121 * source[158] + c77 * source[151] + c77 * source[153]
                  - c121 * source[178] + c77 * source[171] + c77 * source[173]
                  + c79 * source[68] - c82 * source[61] - c82 * source[63]
                  + c78 * source[88] - c84 * source[81] - c84 * source[83]
                  + c79 * source[108] - c82 * source[101] - c82 * source[103];
    target[76] =  source[209] - c122 * source[204] - c122 * source[206]
                  - c123 * source[159] + c85 * source[154] + c85 * source[156]
                  - c123 * source[179] + c85 * source[174] + c85 * source[176]
                  + c88 * source[69] - c124 * source[64] - c124 * source[66]
                  + c91 * source[89] - c125 * source[84] - c125 * source[86]
                  + c88 * source[109] - c124 * source[104] - c124 * source[106];
  }
}

void CCarSphList::carsph_53(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c26 = 51.553976568253198;
  const double c49 = 48.605555238058955;
  const double c74 = 39.686269665968858;
  const double c35 = 34.369317712168801;
  const double c28 = 32.605597678926237;
  const double c24 = 31.57023420090513;
  const double c54 = 30.740852297878796;
  const double c40 = 29.764702249476645;
  const double c10 = 27.171331399105199;
  const double c25 = 25.776988284126599;
  const double c75 = 25.099800796022265;
  const double c45 = 24.302777619029477;
  const double c90 = 22.5;
  const double c36 = 21.737065119284157;
  const double c34 = 21.046822800603419;
  const double c33 = 19.966769267961205;
  const double c71 = 19.843134832984429;
  const double c117 = 19.364916731037084;
  const double c63 = 18.824850597016699;
  const double c13 = 17.1846588560844;
  const double c3 = 16.638974389967672;
  const double c48 = 16.201851746019649;
  const double c76 = 15.370426148939398;
  const double c89 = 15;
  const double c118 = 14.523687548277813;
  const double c94 = 14.230249470757707;
  const double c80 = 13.778379803155376;
  const double c7 = 13.5856656995526;
  const double c32 = 13.311179511974137;
  const double c62 = 12.549900398011133;
  const double c121 = 12.24744871391589;
  const double c69 = 12.151388809514739;
  const double c108 = 11.858541225631422;
  const double c86 = 11.25;
  const double c20 = 10.523411400301709;
  const double c52 = 10.246950765959598;
  const double c39 = 9.9215674164922145;
  const double c114 = 9.6824583655185421;
  const double c92 = 9.4868329805051381;
  const double c78 = 9.1855865354369186;
  const double c112 = 8.8939059192235668;
  const double c37 = 8.8741196746494246;
  const double c101 = 8.7142125289666872;
  const double c15 = 8.5923294280422002;
  const double c5 = 8.3194871949838358;
  const double c29 = 8.1513994197315593;
  const double c44 = 8.1009258730098246;
  const double c55 = 7.6852130744696989;
  const double c85 = 7.5;
  const double c116 = 7.2618437741389066;
  const double c19 = 7.0156076002011405;
  const double c8 = 6.7928328497762998;
  const double c61 = 6.2749501990055663;
  const double c51 = 6.0756944047573693;
  const double c100 = 5.8094750193111251;
  const double c125 = 5.625;
  const double c2 = 5.54632479665589;
  const double c27 = 5.4342662798210393;
  const double c22 = 5.2617057001508547;
  const double c73 = 5.123475382979799;
  const double c123 = 5;
  const double c70 = 4.9607837082461073;
  const double c79 = 4.5927932677184593;
  const double c110 = 4.4469529596117834;
  const double c14 = 4.2961647140211001;
  const double c60 = 4.1833001326703778;
  const double c68 = 4.0504629365049123;
  const double c107 = 3.9528470752104741;
  const double c99 = 3.872983346207417;
  const double c58 = 3.8426065372348495;
  const double c91 = 3.75;
  const double c43 = 3.7205877811845807;
  const double c115 = 3.6309218870694533;
  const double c95 = 3.5575623676894268;
  const double c21 = 3.5078038001005702;
  const double c31 = 3.3277948779935342;
  const double c38 = 3.3071891388307382;
  const double c72 = 3.1374750995027831;
  const double c77 = 3.0618621784789726;
  const double c47 = 3.0378472023786847;
  const double c111 = 2.9646353064078554;
  const double c124 = 2.8125;
  const double c4 = 2.773162398327945;
  const double c9 = 2.7171331399105196;
  const double c53 = 2.5617376914898995;
  const double c119 = 2.4494897427831779;
  const double c93 = 2.3717082451262845;
  const double c67 = 2.3531063246270874;
  const double c84 = 2.2963966338592297;
  const double c30 = 2.2185299186623562;
  const double c16 = 2.1480823570105501;
  const double c50 = 2.0252314682524561;
  const double c113 = 1.9364916731037085;
  const double c88 = 1.875;
  const double c23 = 1.7539019000502851;
  const double c11 = 1.71846588560844;
  const double c1 = 1.6638974389967671;
  const double c66 = 1.5687375497513916;
  const double c122 = 1.5;
  const double c109 = 1.4823176532039277;
  const double c105 = 1.4523687548277813;
  const double c6 = 1.3585665699552598;
  const double c56 = 1.2808688457449497;
  const double c42 = 1.2401959270615268;
  const double c96 = 1.1858541225631423;
  const double c82 = 1.1481983169296148;
  const double c18 = 1.052341140030171;
  const double c46 = 1.0126157341262281;
  const double c104 = 0.96824583655185426;
  const double c59 = 0.96065163430871237;
  const double c87 = 0.9375;
  const double c106 = 0.79056941504209488;
  const double c65 = 0.78436877487569578;
  const double c83 = 0.76546554461974314;
  const double c103 = 0.72618437741389064;
  const double c17 = 0.70156076002011403;
  const double c120 = 0.61237243569579447;
  const double c98 = 0.59292706128157113;
  const double c0 = 0.55463247966558904;
  const double c64 = 0.52291251658379723;
  const double c102 = 0.48412291827592713;
  const double c12 = 0.42961647140211001;
  const double c41 = 0.41339864235384227;
  const double c81 = 0.38273277230987157;
  const double c57 = 0.32021721143623744;
  const double c97 = 0.29646353064078557;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 77, source += 210) {
    target[0] =  c0 * source[0] - c1 * source[2] - c2 * source[20]
                  + c3 * source[22] + c4 * source[40] - c5 * source[42];
    target[1] =  c1 * source[1] - c0 * source[3] - c3 * source[21]
                  + c2 * source[23] + c5 * source[41] - c4 * source[43];
    target[2] =  c6 * source[4] - c6 * source[6] - c7 * source[24]
                  + c7 * source[26] + c8 * source[44] - c8 * source[46];
    target[3] =  c9 * source[5] - c10 * source[25] + c7 * source[45];
    target[4] =  c11 * source[7] - c12 * source[0] - c12 * source[2]
                  - c13 * source[27] + c14 * source[20] + c14 * source[22]
                  + c15 * source[47] - c16 * source[40] - c16 * source[42];
    target[5] =  c11 * source[8] - c12 * source[1] - c12 * source[3]
                  - c13 * source[28] + c14 * source[21] + c14 * source[23]
                  + c15 * source[48] - c16 * source[41] - c16 * source[43];
    target[6] =  c17 * source[9] - c18 * source[4] - c18 * source[6]
                  - c19 * source[29] + c20 * source[24] + c20 * source[26]
                  + c21 * source[49] - c22 * source[44] - c22 * source[46];
    target[7] =  c4 * source[10] - c5 * source[12] - c2 * source[30]
                  + c3 * source[32] + c0 * source[50] - c1 * source[52];
    target[8] =  c5 * source[11] - c4 * source[13] - c3 * source[31]
                  + c2 * source[33] + c1 * source[51] - c0 * source[53];
    target[9] =  c8 * source[14] - c8 * source[16] - c7 * source[34]
                  + c7 * source[36] + c6 * source[54] - c6 * source[56];
    target[10] =  c7 * source[15] - c10 * source[35] + c9 * source[55];
    target[11] =  c15 * source[17] - c16 * source[10] - c16 * source[12]
                  - c13 * source[37] + c14 * source[30] + c14 * source[32]
                  + c11 * source[57] - c12 * source[50] - c12 * source[52];
    target[12] =  c15 * source[18] - c16 * source[11] - c16 * source[13]
                  - c13 * source[38] + c14 * source[31] + c14 * source[33]
                  + c11 * source[58] - c12 * source[51] - c12 * source[53];
    target[13] =  c21 * source[19] - c22 * source[14] - c22 * source[16]
                  - c19 * source[39] + c20 * source[34] + c20 * source[36]
                  + c17 * source[59] - c18 * source[54] - c18 * source[56];
    target[14] =  c23 * source[60] - c22 * source[62] - c20 * source[80]
                  + c24 * source[82] + c23 * source[100] - c22 * source[102];
    target[15] =  c22 * source[61] - c23 * source[63] - c24 * source[81]
                  + c20 * source[83] + c22 * source[101] - c23 * source[103];
    target[16] =  c14 * source[64] - c14 * source[66] - c25 * source[84]
                  + c25 * source[86] + c14 * source[104] - c14 * source[106];
    target[17] =  c15 * source[65] - c26 * source[85] + c15 * source[105];
    target[18] =  c27 * source[67] - c6 * source[60] - c6 * source[62]
                  - c28 * source[87] + c29 * source[80] + c29 * source[82]
                  + c27 * source[107] - c6 * source[100] - c6 * source[102];
    target[19] =  c27 * source[68] - c6 * source[61] - c6 * source[63]
                  - c28 * source[88] + c29 * source[81] + c29 * source[83]
                  + c27 * source[108] - c6 * source[101] - c6 * source[103];
    target[20] =  c30 * source[69] - c31 * source[64] - c31 * source[66]
                  - c32 * source[89] + c33 * source[84] + c33 * source[86]
                  + c30 * source[109] - c31 * source[104] - c31 * source[106];
    target[21] =  c19 * source[70] - c34 * source[72] - c19 * source[90]
                  + c34 * source[92];
    target[22] =  c34 * source[71] - c19 * source[73] - c34 * source[91]
                  + c19 * source[93];
    target[23] =  c13 * source[74] - c13 * source[76] - c13 * source[94]
                  + c13 * source[96];
    target[24] =  c35 * source[75] - c35 * source[95];
    target[25] =  c36 * source[77] - c27 * source[70] - c27 * source[72]
                  - c36 * source[97] + c27 * source[90] + c27 * source[92];
    target[26] =  c36 * source[78] - c27 * source[71] - c27 * source[73]
                  - c36 * source[98] + c27 * source[91] + c27 * source[93];
    target[27] =  c37 * source[79] - c32 * source[74] - c32 * source[76]
                  - c37 * source[99] + c32 * source[94] + c32 * source[96];
    target[28] =  c38 * source[110] - c39 * source[112] - c39 * source[130]
                  + c40 * source[132] - c41 * source[0] + c42 * source[2]
                  + c42 * source[20] - c43 * source[22] - c41 * source[20]
                  + c42 * source[22] + c42 * source[40] - c43 * source[42];
    target[29] =  c39 * source[111] - c38 * source[113] - c40 * source[131]
                  + c39 * source[133] - c42 * source[1] + c41 * source[3]
                  + c43 * source[21] - c42 * source[23] - c42 * source[21]
                  + c41 * source[23] + c43 * source[41] - c42 * source[43];
    target[30] =  c44 * source[114] - c44 * source[116] - c45 * source[134]
                  + c45 * source[136] - c46 * source[4] + c46 * source[6]
                  + c47 * source[24] - c47 * source[26] - c46 * source[24]
                  + c46 * source[26] + c47 * source[44] - c47 * source[46];
    target[31] =  c48 * source[115] - c49 * source[135] - c50 * source[5]
                  + c51 * source[25] - c50 * source[25] + c51 * source[45];
    target[32] =  c52 * source[117] - c53 * source[110] - c53 * source[112]
                  - c54 * source[137] + c55 * source[130] + c55 * source[132]
                  - c56 * source[7] + c57 * source[0] + c57 * source[2]
                  + c58 * source[27] - c59 * source[20] - c59 * source[22]
                  - c56 * source[27] + c57 * source[20] + c57 * source[22]
                  + c58 * source[47] - c59 * source[40] - c59 * source[42];
    target[33] =  c52 * source[118] - c53 * source[111] - c53 * source[113]
                  - c54 * source[138] + c55 * source[131] + c55 * source[133]
                  - c56 * source[8] + c57 * source[1] + c57 * source[3]
                  + c58 * source[28] - c59 * source[21] - c59 * source[23]
                  - c56 * source[28] + c57 * source[21] + c57 * source[23]
                  + c58 * source[48] - c59 * source[41] - c59 * source[43];
    target[34] =  c60 * source[119] - c61 * source[114] - c61 * source[116]
                  - c62 * source[139] + c63 * source[134] + c63 * source[136]
                  - c64 * source[9] + c65 * source[4] + c65 * source[6]
                  + c66 * source[29] - c67 * source[24] - c67 * source[26]
                  - c64 * source[29] + c65 * source[24] + c65 * source[26]
                  + c66 * source[49] - c67 * source[44] - c67 * source[46];
    target[35] =  c39 * source[120] - c40 * source[122] - c38 * source[140]
                  + c39 * source[142] - c42 * source[10] + c43 * source[12]
                  + c41 * source[30] - c42 * source[32] - c42 * source[30]
                  + c43 * source[32] + c41 * source[50] - c42 * source[52];
    target[36] =  c40 * source[121] - c39 * source[123] - c39 * source[141]
                  + c38 * source[143] - c43 * source[11] + c42 * source[13]
                  + c42 * source[31] - c41 * source[33] - c43 * source[31]
                  + c42 * source[33] + c42 * source[51] - c41 * source[53];
    target[37] =  c45 * source[124] - c45 * source[126] - c44 * source[144]
                  + c44 * source[146] - c47 * source[14] + c47 * source[16]
                  + c46 * source[34] - c46 * source[36] - c47 * source[34]
                  + c47 * source[36] + c46 * source[54] - c46 * source[56];
    target[38] =  c49 * source[125] - c48 * source[145] - c51 * source[15]
                  + c50 * source[35] - c51 * source[35] + c50 * source[55];
    target[39] =  c54 * source[127] - c55 * source[120] - c55 * source[122]
                  - c52 * source[147] + c53 * source[140] + c53 * source[142]
                  - c58 * source[17] + c59 * source[10] + c59 * source[12]
                  + c56 * source[37] - c57 * source[30] - c57 * source[32]
                  - c58 * source[37] + c59 * source[30] + c59 * source[32]
                  + c56 * source[57] - c57 * source[50] - c57 * source[52];
    target[40] =  c54 * source[128] - c55 * source[121] - c55 * source[123]
                  - c52 * source[148] + c53 * source[141] + c53 * source[143]
                  - c58 * source[18] + c59 * source[11] + c59 * source[13]
                  + c56 * source[38] - c57 * source[31] - c57 * source[33]
                  - c58 * source[38] + c59 * source[31] + c59 * source[33]
                  + c56 * source[58] - c57 * source[51] - c57 * source[53];
    target[41] =  c62 * source[129] - c63 * source[124] - c63 * source[126]
                  - c60 * source[149] + c61 * source[144] + c61 * source[146]
                  - c66 * source[19] + c67 * source[14] + c67 * source[16]
                  + c64 * source[39] - c65 * source[34] - c65 * source[36]
                  - c66 * source[39] + c67 * source[34] + c67 * source[36]
                  + c64 * source[59] - c65 * source[54] - c65 * source[56];
    target[42] =  c68 * source[150] - c69 * source[152] - c68 * source[170]
                  + c69 * source[172] - c50 * source[60] + c51 * source[62]
                  + c50 * source[80] - c51 * source[82] - c50 * source[80]
                  + c51 * source[82] + c50 * source[100] - c51 * source[102];
    target[43] =  c69 * source[151] - c68 * source[153] - c69 * source[171]
                  + c68 * source[173] - c51 * source[61] + c50 * source[63]
                  + c51 * source[81] - c50 * source[83] - c51 * source[81]
                  + c50 * source[83] + c51 * source[101] - c50 * source[103];
    target[44] =  c39 * source[154] - c39 * source[156] - c39 * source[174]
                  + c39 * source[176] - c70 * source[64] + c70 * source[66]
                  + c70 * source[84] - c70 * source[86] - c70 * source[84]
                  + c70 * source[86] + c70 * source[104] - c70 * source[106];
    target[45] =  c71 * source[155] - c71 * source[175] - c39 * source[65]
                  + c39 * source[85] - c39 * source[85] + c39 * source[105];
    target[46] =  c62 * source[157] - c72 * source[150] - c72 * source[152]
                  - c62 * source[177] + c72 * source[170] + c72 * source[172]
                  - c61 * source[67] + c66 * source[60] + c66 * source[62]
                  + c61 * source[87] - c66 * source[80] - c66 * source[82]
                  - c61 * source[87] + c66 * source[80] + c66 * source[82]
                  + c61 * source[107] - c66 * source[100] - c66 * source[102];
    target[47] =  c62 * source[158] - c72 * source[151] - c72 * source[153]
                  - c62 * source[178] + c72 * source[171] + c72 * source[173]
                  - c61 * source[68] + c66 * source[61] + c66 * source[63]
                  + c61 * source[88] - c66 * source[81] - c66 * source[83]
                  - c61 * source[88] + c66 * source[81] + c66 * source[83]
                  + c61 * source[108] - c66 * source[101] - c66 * source[103];
    target[48] =  c73 * source[159] - c55 * source[154] - c55 * source[156]
                  - c73 * source[179] + c55 * source[174] + c55 * source[176]
                  - c53 * source[69] + c58 * source[64] + c58 * source[66]
                  + c53 * source[89] - c58 * source[84] - c58 * source[86]
                  - c53 * source[89] + c58 * source[84] + c58 * source[86]
                  + c53 * source[109] - c58 * source[104] - c58 * source[106];
    target[49] =  c44 * source[160] - c45 * source[162] - c68 * source[70]
                  + c69 * source[72] - c68 * source[90] + c69 * source[92];
    target[50] =  c45 * source[161] - c44 * source[163] - c69 * source[71]
                  + c68 * source[73] - c69 * source[91] + c68 * source[93];
    target[51] =  c71 * source[164] - c71 * source[166] - c39 * source[74]
                  + c39 * source[76] - c39 * source[94] + c39 * source[96];
    target[52] =  c74 * source[165] - c71 * source[75] - c71 * source[95];
    target[53] =  c75 * source[167] - c61 * source[160] - c61 * source[162]
                  - c62 * source[77] + c72 * source[70] + c72 * source[72]
                  - c62 * source[97] + c72 * source[90] + c72 * source[92];
    target[54] =  c75 * source[168] - c61 * source[161] - c61 * source[163]
                  - c62 * source[78] + c72 * source[71] + c72 * source[73]
                  - c62 * source[98] + c72 * source[91] + c72 * source[93];
    target[55] =  c52 * source[169] - c76 * source[164] - c76 * source[166]
                  - c73 * source[79] + c55 * source[74] + c55 * source[76]
                  - c73 * source[99] + c55 * source[94] + c55 * source[96];
    target[56] =  c77 * source[180] - c78 * source[182] - c79 * source[110]
                  + c80 * source[112] - c79 * source[130] + c80 * source[132]
                  + c81 * source[0] - c82 * source[2] + c83 * source[20]
                  - c84 * source[22] + c81 * source[40] - c82 * source[42];
    target[57] =  c78 * source[181] - c77 * source[183] - c80 * source[111]
                  + c79 * source[113] - c80 * source[131] + c79 * source[133]
                  + c82 * source[1] - c81 * source[3] + c84 * source[21]
                  - c83 * source[23] + c82 * source[41] - c81 * source[43];
    target[58] =  c85 * source[184] - c85 * source[186] - c86 * source[114]
                  + c86 * source[116] - c86 * source[134] + c86 * source[136]
                  + c87 * source[4] - c87 * source[6] + c88 * source[24]
                  - c88 * source[26] + c87 * source[44] - c87 * source[46];
    target[59] =  c89 * source[185] - c90 * source[115] - c90 * source[135]
                  + c88 * source[5] + c91 * source[25] + c88 * source[45];
    target[60] =  c92 * source[187] - c93 * source[180] - c93 * source[182]
                  - c94 * source[117] + c95 * source[110] + c95 * source[112]
                  - c94 * source[137] + c95 * source[130] + c95 * source[132]
                  + c96 * source[7] - c97 * source[0] - c97 * source[2]
                  + c93 * source[27] - c98 * source[20] - c98 * source[22]
                  + c96 * source[47] - c97 * source[40] - c97 * source[42];
    target[61] =  c92 * source[188] - c93 * source[181] - c93 * source[183]
                  - c94 * source[118] + c95 * source[111] + c95 * source[113]
                  - c94 * source[138] + c95 * source[131] + c95 * source[133]
                  + c96 * source[8] - c97 * source[1] - c97 * source[3]
                  + c93 * source[28] - c98 * source[21] - c98 * source[23]
                  + c96 * source[48] - c97 * source[41] - c97 * source[43];
    target[62] =  c99 * source[189] - c100 * source[184] - c100 * source[186]
                  - c100 * source[119] + c101 * source[114] + c101 * source[116]
                  - c100 * source[139] + c101 * source[134] + c101 * source[136]
                  + c102 * source[9] - c103 * source[4] - c103 * source[6]
                  + c104 * source[29] - c105 * source[24] - c105 * source[26]
                  + c102 * source[49] - c103 * source[44] - c103 * source[46];
    target[63] =  c77 * source[190] - c78 * source[192] - c79 * source[120]
                  + c80 * source[122] - c79 * source[140] + c80 * source[142]
                  + c81 * source[10] - c82 * source[12] + c83 * source[30]
                  - c84 * source[32] + c81 * source[50] - c82 * source[52];
    target[64] =  c78 * source[191] - c77 * source[193] - c80 * source[121]
                  + c79 * source[123] - c80 * source[141] + c79 * source[143]
                  + c82 * source[11] - c81 * source[13] + c84 * source[31]
                  - c83 * source[33] + c82 * source[51] - c81 * source[53];
    target[65] =  c85 * source[194] - c85 * source[196] - c86 * source[124]
                  + c86 * source[126] - c86 * source[144] + c86 * source[146]
                  + c87 * source[14] - c87 * source[16] + c88 * source[34]
                  - c88 * source[36] + c87 * source[54] - c87 * source[56];
    target[66] =  c89 * source[195] - c90 * source[125] - c90 * source[145]
                  + c88 * source[15] + c91 * source[35] + c88 * source[55];
    target[67] =  c92 * source[197] - c93 * source[190] - c93 * source[192]
                  - c94 * source[127] + c95 * source[120] + c95 * source[122]
                  - c94 * source[147] + c95 * source[140] + c95 * source[142]
                  + c96 * source[17] - c97 * source[10] - c97 * source[12]
                  + c93 * source[37] - c98 * source[30] - c98 * source[32]
                  + c96 * source[57] - c97 * source[50] - c97 * source[52];
    target[68] =  c92 * source[198] - c93 * source[191] - c93 * source[193]
                  - c94 * source[128] + c95 * source[121] + c95 * source[123]
                  - c94 * source[148] + c95 * source[141] + c95 * source[143]
                  + c96 * source[18] - c97 * source[11] - c97 * source[13]
                  + c93 * source[38] - c98 * source[31] - c98 * source[33]
                  + c96 * source[58] - c97 * source[51] - c97 * source[53];
    target[69] =  c99 * source[199] - c100 * source[194] - c100 * source[196]
                  - c100 * source[129] + c101 * source[124] + c101 * source[126]
                  - c100 * source[149] + c101 * source[144] + c101 * source[146]
                  + c102 * source[19] - c103 * source[14] - c103 * source[16]
                  + c104 * source[39] - c105 * source[34] - c105 * source[36]
                  + c102 * source[59] - c103 * source[54] - c103 * source[56];
    target[70] =  c106 * source[200] - c93 * source[202] - c107 * source[150]
                  + c108 * source[152] - c107 * source[170] + c108 * source[172]
                  + c109 * source[60] - c110 * source[62] + c111 * source[80]
                  - c112 * source[82] + c109 * source[100] - c110 * source[102];
    target[71] =  c93 * source[201] - c106 * source[203] - c108 * source[151]
                  + c107 * source[153] - c108 * source[171] + c107 * source[173]
                  + c110 * source[61] - c109 * source[63] + c112 * source[81]
                  - c111 * source[83] + c110 * source[101] - c109 * source[103];
    target[72] =  c113 * source[204] - c113 * source[206] - c114 * source[154]
                  + c114 * source[156] - c114 * source[174] + c114 * source[176]
                  + c115 * source[64] - c115 * source[66] + c116 * source[84]
                  - c116 * source[86] + c115 * source[104] - c115 * source[106];
    target[73] =  c99 * source[205] - c117 * source[155] - c117 * source[175]
                  + c116 * source[65] + c118 * source[85] + c116 * source[105];
    target[74] =  c119 * source[207] - c120 * source[200] - c120 * source[202]
                  - c121 * source[157] + c77 * source[150] + c77 * source[152]
                  - c121 * source[177] + c77 * source[170] + c77 * source[172]
                  + c79 * source[67] - c82 * source[60] - c82 * source[62]
                  + c78 * source[87] - c84 * source[80] - c84 * source[82]
                  + c79 * source[107] - c82 * source[100] - c82 * source[102];
    target[75] =  c119 * source[208] - c120 * source[201] - c120 * source[203]
                  - c121 * source[158] + c77 * source[151] + c77 * source[153]
                  - c121 * source[178] + c77 * source[171] + c77 * source[173]
                  + c79 * source[68] - c82 * source[61] - c82 * source[63]
                  + c78 * source[88] - c84 * source[81] - c84 * source[83]
                  + c79 * source[108] - c82 * source[101] - c82 * source[103];
    target[76] =  source[209] - c122 * source[204] - c122 * source[206]
                  - c123 * source[159] + c85 * source[154] + c85 * source[156]
                  - c123 * source[179] + c85 * source[174] + c85 * source[176]
                  + c88 * source[69] - c124 * source[64] - c124 * source[66]
                  + c91 * source[89] - c125 * source[84] - c125 * source[86]
                  + c88 * source[109] - c124 * source[104] - c124 * source[106];
  }
}

