//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_72.cc
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

#ifdef COMPILE_J_ORB
#include <src/integral/carsphlist.h>
#include <algorithm>

using namespace std;
using namespace bagel;


void CarSphList::carsph_72(const int nloop, const double* source, double* target) {
  const double c39 = 98.718349358161376;
  const double c27 = 83.894464954489095;
  const double c61 = 65.812232905440922;
  const double c19 = 62.920848715866825;
  const double c75 = 59.529404498953291;
  const double c44 = 56.995065575889988;
  const double c109 = 56.124860801609124;
  const double c33 = 49.359174679080688;
  const double c30 = 48.436491924993909;
  const double c135 = 45.46633369868303;
  const double c77 = 44.647053374214963;
  const double c25 = 41.947232477244548;
  const double c6 = 39.238043063301461;
  const double c63 = 37.996710383926661;
  const double c22 = 36.327368943745434;
  const double c82 = 34.369317712168801;
  const double c108 = 33.674916480965472;
  const double c60 = 32.906116452720461;
  const double c112 = 32.403703492039298;
  const double c17 = 31.460424357933412;
  const double c67 = 29.764702249476645;
  const double c55 = 29.615504807448414;
  const double c45 = 28.497532787944994;
  const double c99 = 28.062430400804562;
  const double c142 = 26.25;
  const double c86 = 25.776988284126599;
  const double c26 = 25.168339486346731;
  const double c34 = 24.679587339540344;
  const double c31 = 24.218245962496955;
  const double c5 = 23.542825837980878;
  const double c130 = 22.733166849341515;
  const double c12 = 22.654094725071229;
  const double c69 = 22.323526687107481;
  const double c110 = 21.046822800603419;
  const double c74 = 19.843134832984429;
  const double c62 = 19.743669871632274;
  const double c2 = 19.619021531650731;
  const double c111 = 19.442222095223581;
  const double c64 = 18.998355191963331;
  const double c134 = 18.186533479473212;
  const double c23 = 18.163684471872717;
  const double c83 = 17.1846588560844;
  const double c58 = 17.098519672766997;
  const double c98 = 16.837458240482736;
  const double c54 = 16.453058226360231;
  const double c103 = 16.201851746019649;
  const double c76 = 14.882351124738323;
  const double c53 = 14.807752403724207;
  const double c28 = 14.530947577498171;
  const double c46 = 14.248766393972497;
  const double c95 = 14.031215200402281;
  const double c10 = 13.592456835042736;
  const double c140 = 13.125;
  const double c87 = 12.888494142063299;
  const double c24 = 12.584169743173366;
  const double c113 = 12.151388809514739;
  const double c1 = 11.771412918990439;
  const double c80 = 11.456439237389599;
  const double c65 = 11.399013115177997;
  const double c129 = 11.366583424670758;
  const double c13 = 11.327047362535614;
  const double c100 = 10.523411400301709;
  const double c138 = 10.5;
  const double c66 = 9.9215674164922145;
  const double c38 = 9.8718349358161372;
  const double c101 = 9.7211110476117906;
  const double c56 = 9.4991775959816653;
  const double c117 = 9.1651513899116797;
  const double c128 = 9.093266739736606;
  const double c84 = 8.5923294280422002;
  const double c59 = 8.5492598363834986;
  const double c94 = 8.4187291202413679;
  const double c41 = 8.2265291131801153;
  const double c104 = 8.1009258730098246;
  const double c7 = 7.8476086126602924;
  const double c68 = 7.4411755623691613;
  const double c29 = 7.2654737887490857;
  const double c11 = 6.7962284175213679;
  const double c141 = 6.5625;
  const double c107 = 6.0756944047573693;
  const double c81 = 5.7282196186947996;
  const double c42 = 5.6995065575889985;
  const double c132 = 5.6832917123353788;
  const double c120 = 5.2915026221291814;
  const double c97 = 5.2617057001508547;
  const double c139 = 5.25;
  const double c122 = 4.9607837082461073;
  const double c32 = 4.9359174679080686;
  const double c102 = 4.8605555238058953;
  const double c49 = 4.7495887979908327;
  const double c114 = 4.5825756949558398;
  const double c14 = 4.5308189450142455;
  const double c79 = 4.4647053374214964;
  const double c85 = 4.2961647140211001;
  const double c18 = 4.1947232477244549;
  const double c36 = 4.1132645565900576;
  const double c3 = 3.9238043063301462;
  const double c136 = 3.7888611415569189;
  const double c145 = 3.28125;
  const double c105 = 3.0378472023786847;
  const double c43 = 2.8497532787944992;
  const double c121 = 2.6457513110645907;
  const double c96 = 2.6308528500754274;
  const double c93 = 2.5776988284126601;
  const double c52 = 2.4679587339540343;
  const double c20 = 2.4218245962496954;
  const double c50 = 2.3747943989954163;
  const double c15 = 2.2654094725071228;
  const double c73 = 2.2323526687107482;
  const double c143 = 2.1875;
  const double c119 = 2.1480823570105501;
  const double c16 = 2.0973616238622275;
  const double c37 = 2.0566322782950288;
  const double c131 = 1.8944305707784594;
  const double c133 = 1.7320508075688772;
  const double c106 = 1.5189236011893423;
  const double c78 = 1.4882351124738322;
  const double c57 = 1.4248766393972496;
  const double c90 = 1.28884941420633;
  const double c125 = 1.2401959270615268;
  const double c21 = 1.2109122981248477;
  const double c51 = 1.1873971994977082;
  const double c4 = 1.1210869446657561;
  const double c71 = 1.1161763343553741;
  const double c144 = 1.09375;
  const double c116 = 1.074041178505275;
  const double c127 = 0.8660254037844386;
  const double c92 = 0.85923294280422002;
  const double c40 = 0.82265291131801144;
  const double c72 = 0.74411755623691611;
  const double c118 = 0.71602745233684995;
  const double c8 = 0.64725984928774938;
  const double c91 = 0.64442470710316502;
  const double c126 = 0.62009796353076341;
  const double c0 = 0.56054347233287805;
  const double c137 = 0.5;
  const double c47 = 0.47495887979908324;
  const double c88 = 0.42961647140211001;
  const double c123 = 0.41339864235384227;
  const double c35 = 0.41132645565900572;
  const double c70 = 0.37205877811845806;
  const double c115 = 0.35801372616842497;
  const double c9 = 0.32362992464387469;
  const double c48 = 0.23747943989954162;
  const double c89 = 0.21480823570105501;
  const double c124 = 0.20669932117692114;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 75, source += 216) {
    target[0] =  c0 * source[0] - c0 * source[2] - c1 * source[12]
                  + c1 * source[14] + c2 * source[24] - c2 * source[26]
                  - c3 * source[36] + c3 * source[38];
    target[1] =  c4 * source[1] - c5 * source[13] + c6 * source[25]
                  - c7 * source[37];
    target[2] =  c4 * source[3] - c5 * source[15] + c6 * source[27]
                  - c7 * source[39];
    target[3] =  c4 * source[4] - c5 * source[16] + c6 * source[28]
                  - c7 * source[40];
    target[4] =  c8 * source[5] - c9 * source[0] - c9 * source[2]
                  - c10 * source[17] + c11 * source[12] + c11 * source[14]
                  + c12 * source[29] - c13 * source[24] - c13 * source[26]
                  - c14 * source[41] + c15 * source[36] + c15 * source[38];
    target[5] =  c3 * source[6] - c3 * source[8] - c2 * source[18]
                  + c2 * source[20] + c1 * source[30] - c1 * source[32]
                  - c0 * source[42] + c0 * source[44];
    target[6] =  c7 * source[7] - c6 * source[19] + c5 * source[31]
                  - c4 * source[43];
    target[7] =  c7 * source[9] - c6 * source[21] + c5 * source[33]
                  - c4 * source[45];
    target[8] =  c7 * source[10] - c6 * source[22] + c5 * source[34]
                  - c4 * source[46];
    target[9] =  c14 * source[11] - c15 * source[6] - c15 * source[8]
                  - c12 * source[23] + c13 * source[18] + c13 * source[20]
                  + c10 * source[35] - c11 * source[30] - c11 * source[32]
                  - c8 * source[47] + c9 * source[42] + c9 * source[44];
    target[10] =  c16 * source[48] - c16 * source[50] - c17 * source[60]
                  + c17 * source[62] + c17 * source[72] - c17 * source[74]
                  - c16 * source[84] + c16 * source[86];
    target[11] =  c18 * source[49] - c19 * source[61] + c19 * source[73]
                  - c18 * source[85];
    target[12] =  c18 * source[51] - c19 * source[63] + c19 * source[75]
                  - c18 * source[87];
    target[13] =  c18 * source[52] - c19 * source[64] + c19 * source[76]
                  - c18 * source[88];
    target[14] =  c20 * source[53] - c21 * source[48] - c21 * source[50]
                  - c22 * source[65] + c23 * source[60] + c23 * source[62]
                  + c22 * source[77] - c23 * source[72] - c23 * source[74]
                  - c20 * source[89] + c21 * source[84] + c21 * source[86];
    target[15] =  c24 * source[54] - c24 * source[56] - c25 * source[66]
                  + c25 * source[68] + c24 * source[78] - c24 * source[80];
    target[16] =  c26 * source[55] - c27 * source[67] + c26 * source[79];
    target[17] =  c26 * source[57] - c27 * source[69] + c26 * source[81];
    target[18] =  c26 * source[58] - c27 * source[70] + c26 * source[82];
    target[19] =  c28 * source[59] - c29 * source[54] - c29 * source[56]
                  - c30 * source[71] + c31 * source[66] + c31 * source[68]
                  + c28 * source[83] - c29 * source[78] - c29 * source[80];
    target[20] =  c32 * source[90] - c32 * source[92] - c33 * source[102]
                  + c33 * source[104] + c34 * source[114] - c34 * source[116]
                  - c35 * source[0] + c35 * source[2] + c36 * source[12]
                  - c36 * source[14] - c37 * source[24] + c37 * source[26]
                  - c35 * source[12] + c35 * source[14] + c36 * source[24]
                  - c36 * source[26] - c37 * source[36] + c37 * source[38];
    target[21] =  c38 * source[91] - c39 * source[103] + c33 * source[115]
                  - c40 * source[1] + c41 * source[13] - c36 * source[25]
                  - c40 * source[13] + c41 * source[25] - c36 * source[37];
    target[22] =  c38 * source[93] - c39 * source[105] + c33 * source[117]
                  - c40 * source[3] + c41 * source[15] - c36 * source[27]
                  - c40 * source[15] + c41 * source[27] - c36 * source[39];
    target[23] =  c38 * source[94] - c39 * source[106] + c33 * source[118]
                  - c40 * source[4] + c41 * source[16] - c36 * source[28]
                  - c40 * source[16] + c41 * source[28] - c36 * source[40];
    target[24] =  c42 * source[95] - c43 * source[90] - c43 * source[92]
                  - c44 * source[107] + c45 * source[102] + c45 * source[104]
                  + c45 * source[119] - c46 * source[114] - c46 * source[116]
                  - c47 * source[5] + c48 * source[0] + c48 * source[2]
                  + c49 * source[17] - c50 * source[12] - c50 * source[14]
                  - c50 * source[29] + c51 * source[24] + c51 * source[26]
                  - c47 * source[17] + c48 * source[12] + c48 * source[14]
                  + c49 * source[29] - c50 * source[24] - c50 * source[26]
                  - c50 * source[41] + c51 * source[36] + c51 * source[38];
    target[25] =  c34 * source[96] - c34 * source[98] - c33 * source[108]
                  + c33 * source[110] + c32 * source[120] - c32 * source[122]
                  - c37 * source[6] + c37 * source[8] + c36 * source[18]
                  - c36 * source[20] - c35 * source[30] + c35 * source[32]
                  - c37 * source[18] + c37 * source[20] + c36 * source[30]
                  - c36 * source[32] - c35 * source[42] + c35 * source[44];
    target[26] =  c33 * source[97] - c39 * source[109] + c38 * source[121]
                  - c36 * source[7] + c41 * source[19] - c40 * source[31]
                  - c36 * source[19] + c41 * source[31] - c40 * source[43];
    target[27] =  c33 * source[99] - c39 * source[111] + c38 * source[123]
                  - c36 * source[9] + c41 * source[21] - c40 * source[33]
                  - c36 * source[21] + c41 * source[33] - c40 * source[45];
    target[28] =  c33 * source[100] - c39 * source[112] + c38 * source[124]
                  - c36 * source[10] + c41 * source[22] - c40 * source[34]
                  - c36 * source[22] + c41 * source[34] - c40 * source[46];
    target[29] =  c45 * source[101] - c46 * source[96] - c46 * source[98]
                  - c44 * source[113] + c45 * source[108] + c45 * source[110]
                  + c42 * source[125] - c43 * source[120] - c43 * source[122]
                  - c50 * source[11] + c51 * source[6] + c51 * source[8]
                  + c49 * source[23] - c50 * source[18] - c50 * source[20]
                  - c47 * source[35] + c48 * source[30] + c48 * source[32]
                  - c50 * source[23] + c51 * source[18] + c51 * source[20]
                  + c49 * source[35] - c50 * source[30] - c50 * source[32]
                  - c47 * source[47] + c48 * source[42] + c48 * source[44];
    target[30] =  c41 * source[126] - c41 * source[128] - c33 * source[138]
                  + c33 * source[140] + c41 * source[150] - c41 * source[152]
                  - c52 * source[48] + c52 * source[50] + c53 * source[60]
                  - c53 * source[62] - c52 * source[72] + c52 * source[74]
                  - c52 * source[60] + c52 * source[62] + c53 * source[72]
                  - c53 * source[74] - c52 * source[84] + c52 * source[86];
    target[31] =  c54 * source[127] - c39 * source[139] + c54 * source[151]
                  - c32 * source[49] + c55 * source[61] - c32 * source[73]
                  - c32 * source[61] + c55 * source[73] - c32 * source[85];
    target[32] =  c54 * source[129] - c39 * source[141] + c54 * source[153]
                  - c32 * source[51] + c55 * source[63] - c32 * source[75]
                  - c32 * source[63] + c55 * source[75] - c32 * source[87];
    target[33] =  c54 * source[130] - c39 * source[142] + c54 * source[154]
                  - c32 * source[52] + c55 * source[64] - c32 * source[76]
                  - c32 * source[64] + c55 * source[76] - c32 * source[88];
    target[34] =  c56 * source[131] - c49 * source[126] - c49 * source[128]
                  - c44 * source[143] + c45 * source[138] + c45 * source[140]
                  + c56 * source[155] - c49 * source[150] - c49 * source[152]
                  - c43 * source[53] + c57 * source[48] + c57 * source[50]
                  + c58 * source[65] - c59 * source[60] - c59 * source[62]
                  - c43 * source[77] + c57 * source[72] + c57 * source[74]
                  - c43 * source[65] + c57 * source[60] + c57 * source[62]
                  + c58 * source[77] - c59 * source[72] - c59 * source[74]
                  - c43 * source[89] + c57 * source[84] + c57 * source[86];
    target[35] =  c60 * source[132] - c60 * source[134] - c60 * source[144]
                  + c60 * source[146] - c38 * source[54] + c38 * source[56]
                  + c38 * source[66] - c38 * source[68] - c38 * source[66]
                  + c38 * source[68] + c38 * source[78] - c38 * source[80];
    target[36] =  c61 * source[133] - c61 * source[145] - c62 * source[55]
                  + c62 * source[67] - c62 * source[67] + c62 * source[79];
    target[37] =  c61 * source[135] - c61 * source[147] - c62 * source[57]
                  + c62 * source[69] - c62 * source[69] + c62 * source[81];
    target[38] =  c61 * source[136] - c61 * source[148] - c62 * source[58]
                  + c62 * source[70] - c62 * source[70] + c62 * source[82];
    target[39] =  c63 * source[137] - c64 * source[132] - c64 * source[134]
                  - c63 * source[149] + c64 * source[144] + c64 * source[146]
                  - c65 * source[59] + c42 * source[54] + c42 * source[56]
                  + c65 * source[71] - c42 * source[66] - c42 * source[68]
                  - c65 * source[71] + c42 * source[66] + c42 * source[68]
                  + c65 * source[83] - c42 * source[78] - c42 * source[80];
    target[40] =  c66 * source[156] - c66 * source[158] - c67 * source[168]
                  + c67 * source[170] - c68 * source[90] + c68 * source[92]
                  + c69 * source[102] - c69 * source[104] - c68 * source[102]
                  + c68 * source[104] + c69 * source[114] - c69 * source[116]
                  + c70 * source[0] - c70 * source[2] - c71 * source[12]
                  + c71 * source[14] + c72 * source[12] - c72 * source[14]
                  - c73 * source[24] + c73 * source[26] + c70 * source[24]
                  - c70 * source[26] - c71 * source[36] + c71 * source[38];
    target[41] =  c74 * source[157] - c75 * source[169] - c76 * source[91]
                  + c77 * source[103] - c76 * source[103] + c77 * source[115]
                  + c72 * source[1] - c73 * source[13] + c78 * source[13]
                  - c79 * source[25] + c72 * source[25] - c73 * source[37];
    target[42] =  c74 * source[159] - c75 * source[171] - c76 * source[93]
                  + c77 * source[105] - c76 * source[105] + c77 * source[117]
                  + c72 * source[3] - c73 * source[15] + c78 * source[15]
                  - c79 * source[27] + c72 * source[27] - c73 * source[39];
    target[43] =  c74 * source[160] - c75 * source[172] - c76 * source[94]
                  + c77 * source[106] - c76 * source[106] + c77 * source[118]
                  + c72 * source[4] - c73 * source[16] + c78 * source[16]
                  - c79 * source[28] + c72 * source[28] - c73 * source[40];
    target[44] =  c80 * source[161] - c81 * source[156] - c81 * source[158]
                  - c82 * source[173] + c83 * source[168] + c83 * source[170]
                  - c84 * source[95] + c85 * source[90] + c85 * source[92]
                  + c86 * source[107] - c87 * source[102] - c87 * source[104]
                  - c84 * source[107] + c85 * source[102] + c85 * source[104]
                  + c86 * source[119] - c87 * source[114] - c87 * source[116]
                  + c88 * source[5] - c89 * source[0] - c89 * source[2]
                  - c90 * source[17] + c91 * source[12] + c91 * source[14]
                  + c92 * source[17] - c88 * source[12] - c88 * source[14]
                  - c93 * source[29] + c90 * source[24] + c90 * source[26]
                  + c88 * source[29] - c89 * source[24] - c89 * source[26]
                  - c90 * source[41] + c91 * source[36] + c91 * source[38];
    target[45] =  c67 * source[162] - c67 * source[164] - c66 * source[174]
                  + c66 * source[176] - c69 * source[96] + c69 * source[98]
                  + c68 * source[108] - c68 * source[110] - c69 * source[108]
                  + c69 * source[110] + c68 * source[120] - c68 * source[122]
                  + c71 * source[6] - c71 * source[8] - c70 * source[18]
                  + c70 * source[20] + c73 * source[18] - c73 * source[20]
                  - c72 * source[30] + c72 * source[32] + c71 * source[30]
                  - c71 * source[32] - c70 * source[42] + c70 * source[44];
    target[46] =  c75 * source[163] - c74 * source[175] - c77 * source[97]
                  + c76 * source[109] - c77 * source[109] + c76 * source[121]
                  + c73 * source[7] - c72 * source[19] + c79 * source[19]
                  - c78 * source[31] + c73 * source[31] - c72 * source[43];
    target[47] =  c75 * source[165] - c74 * source[177] - c77 * source[99]
                  + c76 * source[111] - c77 * source[111] + c76 * source[123]
                  + c73 * source[9] - c72 * source[21] + c79 * source[21]
                  - c78 * source[33] + c73 * source[33] - c72 * source[45];
    target[48] =  c75 * source[166] - c74 * source[178] - c77 * source[100]
                  + c76 * source[112] - c77 * source[112] + c76 * source[124]
                  + c73 * source[10] - c72 * source[22] + c79 * source[22]
                  - c78 * source[34] + c73 * source[34] - c72 * source[46];
    target[49] =  c82 * source[167] - c83 * source[162] - c83 * source[164]
                  - c80 * source[179] + c81 * source[174] + c81 * source[176]
                  - c86 * source[101] + c87 * source[96] + c87 * source[98]
                  + c84 * source[113] - c85 * source[108] - c85 * source[110]
                  - c86 * source[113] + c87 * source[108] + c87 * source[110]
                  + c84 * source[125] - c85 * source[120] - c85 * source[122]
                  + c90 * source[11] - c91 * source[6] - c91 * source[8]
                  - c88 * source[23] + c89 * source[18] + c89 * source[20]
                  + c93 * source[23] - c90 * source[18] - c90 * source[20]
                  - c92 * source[35] + c88 * source[30] + c88 * source[32]
                  + c90 * source[35] - c91 * source[30] - c91 * source[32]
                  - c88 * source[47] + c89 * source[42] + c89 * source[44];
    target[50] =  c94 * source[180] - c94 * source[182] - c94 * source[192]
                  + c94 * source[194] - c95 * source[126] + c95 * source[128]
                  + c95 * source[138] - c95 * source[140] - c95 * source[138]
                  + c95 * source[140] + c95 * source[150] - c95 * source[152]
                  + c96 * source[48] - c96 * source[50] - c96 * source[60]
                  + c96 * source[62] + c97 * source[60] - c97 * source[62]
                  - c97 * source[72] + c97 * source[74] + c96 * source[72]
                  - c96 * source[74] - c96 * source[84] + c96 * source[86];
    target[51] =  c98 * source[181] - c98 * source[193] - c99 * source[127]
                  + c99 * source[139] - c99 * source[139] + c99 * source[151]
                  + c97 * source[49] - c97 * source[61] + c100 * source[61]
                  - c100 * source[73] + c97 * source[73] - c97 * source[85];
    target[52] =  c98 * source[183] - c98 * source[195] - c99 * source[129]
                  + c99 * source[141] - c99 * source[141] + c99 * source[153]
                  + c97 * source[51] - c97 * source[63] + c100 * source[63]
                  - c100 * source[75] + c97 * source[75] - c97 * source[87];
    target[53] =  c98 * source[184] - c98 * source[196] - c99 * source[130]
                  + c99 * source[142] - c99 * source[142] + c99 * source[154]
                  + c97 * source[52] - c97 * source[64] + c100 * source[64]
                  - c100 * source[76] + c97 * source[76] - c97 * source[88];
    target[54] =  c101 * source[185] - c102 * source[180] - c102 * source[182]
                  - c101 * source[197] + c102 * source[192] + c102 * source[194]
                  - c103 * source[131] + c104 * source[126] + c104 * source[128]
                  + c103 * source[143] - c104 * source[138] - c104 * source[140]
                  - c103 * source[143] + c104 * source[138] + c104 * source[140]
                  + c103 * source[155] - c104 * source[150] - c104 * source[152]
                  + c105 * source[53] - c106 * source[48] - c106 * source[50]
                  - c105 * source[65] + c106 * source[60] + c106 * source[62]
                  + c107 * source[65] - c105 * source[60] - c105 * source[62]
                  - c107 * source[77] + c105 * source[72] + c105 * source[74]
                  + c105 * source[77] - c106 * source[72] - c106 * source[74]
                  - c105 * source[89] + c106 * source[84] + c106 * source[86];
    target[55] =  c98 * source[186] - c98 * source[188] - c99 * source[132]
                  + c99 * source[134] - c99 * source[144] + c99 * source[146]
                  + c97 * source[54] - c97 * source[56] + c100 * source[66]
                  - c100 * source[68] + c97 * source[78] - c97 * source[80];
    target[56] =  c108 * source[187] - c109 * source[133] - c109 * source[145]
                  + c100 * source[55] + c110 * source[67] + c100 * source[79];
    target[57] =  c108 * source[189] - c109 * source[135] - c109 * source[147]
                  + c100 * source[57] + c110 * source[69] + c100 * source[81];
    target[58] =  c108 * source[190] - c109 * source[136] - c109 * source[148]
                  + c100 * source[58] + c110 * source[70] + c100 * source[82];
    target[59] =  c111 * source[191] - c101 * source[186] - c101 * source[188]
                  - c112 * source[137] + c103 * source[132] + c103 * source[134]
                  - c112 * source[149] + c103 * source[144] + c103 * source[146]
                  + c107 * source[59] - c105 * source[54] - c105 * source[56]
                  + c113 * source[71] - c107 * source[66] - c107 * source[68]
                  + c107 * source[83] - c105 * source[78] - c105 * source[80];
    target[60] =  c114 * source[198] - c114 * source[200] - c83 * source[156]
                  + c83 * source[158] - c83 * source[168] + c83 * source[170]
                  + c84 * source[90] - c84 * source[92] + c83 * source[102]
                  - c83 * source[104] + c84 * source[114] - c84 * source[116]
                  - c115 * source[0] + c115 * source[2] - c116 * source[12]
                  + c116 * source[14] - c116 * source[24] + c116 * source[26]
                  - c115 * source[36] + c115 * source[38];
    target[61] =  c117 * source[199] - c82 * source[157] - c82 * source[169]
                  + c83 * source[91] + c82 * source[103] + c83 * source[115]
                  - c118 * source[1] - c119 * source[13] - c119 * source[25]
                  - c118 * source[37];
    target[62] =  c117 * source[201] - c82 * source[159] - c82 * source[171]
                  + c83 * source[93] + c82 * source[105] + c83 * source[117]
                  - c118 * source[3] - c119 * source[15] - c119 * source[27]
                  - c118 * source[39];
    target[63] =  c117 * source[202] - c82 * source[160] - c82 * source[172]
                  + c83 * source[94] + c82 * source[106] + c83 * source[118]
                  - c118 * source[4] - c119 * source[16] - c119 * source[28]
                  - c118 * source[40];
    target[64] =  c120 * source[203] - c121 * source[198] - c121 * source[200]
                  - c74 * source[161] + c66 * source[156] + c66 * source[158]
                  - c74 * source[173] + c66 * source[168] + c66 * source[170]
                  + c66 * source[95] - c122 * source[90] - c122 * source[92]
                  + c74 * source[107] - c66 * source[102] - c66 * source[104]
                  + c66 * source[119] - c122 * source[114] - c122 * source[116]
                  - c123 * source[5] + c124 * source[0] + c124 * source[2]
                  - c125 * source[17] + c126 * source[12] + c126 * source[14]
                  - c125 * source[29] + c126 * source[24] + c126 * source[26]
                  - c123 * source[41] + c124 * source[36] + c124 * source[38];
    target[65] =  c114 * source[204] - c114 * source[206] - c83 * source[162]
                  + c83 * source[164] - c83 * source[174] + c83 * source[176]
                  + c84 * source[96] - c84 * source[98] + c83 * source[108]
                  - c83 * source[110] + c84 * source[120] - c84 * source[122]
                  - c115 * source[6] + c115 * source[8] - c116 * source[18]
                  + c116 * source[20] - c116 * source[30] + c116 * source[32]
                  - c115 * source[42] + c115 * source[44];
    target[66] =  c117 * source[205] - c82 * source[163] - c82 * source[175]
                  + c83 * source[97] + c82 * source[109] + c83 * source[121]
                  - c118 * source[7] - c119 * source[19] - c119 * source[31]
                  - c118 * source[43];
    target[67] =  c117 * source[207] - c82 * source[165] - c82 * source[177]
                  + c83 * source[99] + c82 * source[111] + c83 * source[123]
                  - c118 * source[9] - c119 * source[21] - c119 * source[33]
                  - c118 * source[45];
    target[68] =  c117 * source[208] - c82 * source[166] - c82 * source[178]
                  + c83 * source[100] + c82 * source[112] + c83 * source[124]
                  - c118 * source[10] - c119 * source[22] - c119 * source[34]
                  - c118 * source[46];
    target[69] =  c120 * source[209] - c121 * source[204] - c121 * source[206]
                  - c74 * source[167] + c66 * source[162] + c66 * source[164]
                  - c74 * source[179] + c66 * source[174] + c66 * source[176]
                  + c66 * source[101] - c122 * source[96] - c122 * source[98]
                  + c74 * source[113] - c66 * source[108] - c66 * source[110]
                  + c66 * source[125] - c122 * source[120] - c122 * source[122]
                  - c123 * source[11] + c124 * source[6] + c124 * source[8]
                  - c125 * source[23] + c126 * source[18] + c126 * source[20]
                  - c125 * source[35] + c126 * source[30] + c126 * source[32]
                  - c123 * source[47] + c124 * source[42] + c124 * source[44];
    target[70] =  c127 * source[210] - c127 * source[212] - c128 * source[180]
                  + c128 * source[182] - c128 * source[192] + c128 * source[194]
                  + c129 * source[126] - c129 * source[128] + c130 * source[138]
                  - c130 * source[140] + c129 * source[150] - c129 * source[152]
                  - c131 * source[48] + c131 * source[50] - c132 * source[60]
                  + c132 * source[62] - c132 * source[72] + c132 * source[74]
                  - c131 * source[84] + c131 * source[86];
    target[71] =  c133 * source[211] - c134 * source[181] - c134 * source[193]
                  + c130 * source[127] + c135 * source[139] + c130 * source[151]
                  - c136 * source[49] - c129 * source[61] - c129 * source[73]
                  - c136 * source[85];
    target[72] =  c133 * source[213] - c134 * source[183] - c134 * source[195]
                  + c130 * source[129] + c135 * source[141] + c130 * source[153]
                  - c136 * source[51] - c129 * source[63] - c129 * source[75]
                  - c136 * source[87];
    target[73] =  c133 * source[214] - c134 * source[184] - c134 * source[196]
                  + c130 * source[130] + c135 * source[142] + c130 * source[154]
                  - c136 * source[52] - c129 * source[64] - c129 * source[76]
                  - c136 * source[88];
    target[74] =  source[215] - c137 * source[210] - c137 * source[212]
                  - c138 * source[185] + c139 * source[180] + c139 * source[182]
                  - c138 * source[197] + c139 * source[192] + c139 * source[194]
                  + c140 * source[131] - c141 * source[126] - c141 * source[128]
                  + c142 * source[143] - c140 * source[138] - c140 * source[140]
                  + c140 * source[155] - c141 * source[150] - c141 * source[152]
                  - c143 * source[53] + c144 * source[48] + c144 * source[50]
                  - c141 * source[65] + c145 * source[60] + c145 * source[62]
                  - c141 * source[77] + c145 * source[72] + c145 * source[74]
                  - c143 * source[89] + c144 * source[84] + c144 * source[86];
  }
}

void CCarSphList::carsph_72(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c39 = 98.718349358161376;
  const double c27 = 83.894464954489095;
  const double c61 = 65.812232905440922;
  const double c19 = 62.920848715866825;
  const double c75 = 59.529404498953291;
  const double c44 = 56.995065575889988;
  const double c109 = 56.124860801609124;
  const double c33 = 49.359174679080688;
  const double c30 = 48.436491924993909;
  const double c135 = 45.46633369868303;
  const double c77 = 44.647053374214963;
  const double c25 = 41.947232477244548;
  const double c6 = 39.238043063301461;
  const double c63 = 37.996710383926661;
  const double c22 = 36.327368943745434;
  const double c82 = 34.369317712168801;
  const double c108 = 33.674916480965472;
  const double c60 = 32.906116452720461;
  const double c112 = 32.403703492039298;
  const double c17 = 31.460424357933412;
  const double c67 = 29.764702249476645;
  const double c55 = 29.615504807448414;
  const double c45 = 28.497532787944994;
  const double c99 = 28.062430400804562;
  const double c142 = 26.25;
  const double c86 = 25.776988284126599;
  const double c26 = 25.168339486346731;
  const double c34 = 24.679587339540344;
  const double c31 = 24.218245962496955;
  const double c5 = 23.542825837980878;
  const double c130 = 22.733166849341515;
  const double c12 = 22.654094725071229;
  const double c69 = 22.323526687107481;
  const double c110 = 21.046822800603419;
  const double c74 = 19.843134832984429;
  const double c62 = 19.743669871632274;
  const double c2 = 19.619021531650731;
  const double c111 = 19.442222095223581;
  const double c64 = 18.998355191963331;
  const double c134 = 18.186533479473212;
  const double c23 = 18.163684471872717;
  const double c83 = 17.1846588560844;
  const double c58 = 17.098519672766997;
  const double c98 = 16.837458240482736;
  const double c54 = 16.453058226360231;
  const double c103 = 16.201851746019649;
  const double c76 = 14.882351124738323;
  const double c53 = 14.807752403724207;
  const double c28 = 14.530947577498171;
  const double c46 = 14.248766393972497;
  const double c95 = 14.031215200402281;
  const double c10 = 13.592456835042736;
  const double c140 = 13.125;
  const double c87 = 12.888494142063299;
  const double c24 = 12.584169743173366;
  const double c113 = 12.151388809514739;
  const double c1 = 11.771412918990439;
  const double c80 = 11.456439237389599;
  const double c65 = 11.399013115177997;
  const double c129 = 11.366583424670758;
  const double c13 = 11.327047362535614;
  const double c100 = 10.523411400301709;
  const double c138 = 10.5;
  const double c66 = 9.9215674164922145;
  const double c38 = 9.8718349358161372;
  const double c101 = 9.7211110476117906;
  const double c56 = 9.4991775959816653;
  const double c117 = 9.1651513899116797;
  const double c128 = 9.093266739736606;
  const double c84 = 8.5923294280422002;
  const double c59 = 8.5492598363834986;
  const double c94 = 8.4187291202413679;
  const double c41 = 8.2265291131801153;
  const double c104 = 8.1009258730098246;
  const double c7 = 7.8476086126602924;
  const double c68 = 7.4411755623691613;
  const double c29 = 7.2654737887490857;
  const double c11 = 6.7962284175213679;
  const double c141 = 6.5625;
  const double c107 = 6.0756944047573693;
  const double c81 = 5.7282196186947996;
  const double c42 = 5.6995065575889985;
  const double c132 = 5.6832917123353788;
  const double c120 = 5.2915026221291814;
  const double c97 = 5.2617057001508547;
  const double c139 = 5.25;
  const double c122 = 4.9607837082461073;
  const double c32 = 4.9359174679080686;
  const double c102 = 4.8605555238058953;
  const double c49 = 4.7495887979908327;
  const double c114 = 4.5825756949558398;
  const double c14 = 4.5308189450142455;
  const double c79 = 4.4647053374214964;
  const double c85 = 4.2961647140211001;
  const double c18 = 4.1947232477244549;
  const double c36 = 4.1132645565900576;
  const double c3 = 3.9238043063301462;
  const double c136 = 3.7888611415569189;
  const double c145 = 3.28125;
  const double c105 = 3.0378472023786847;
  const double c43 = 2.8497532787944992;
  const double c121 = 2.6457513110645907;
  const double c96 = 2.6308528500754274;
  const double c93 = 2.5776988284126601;
  const double c52 = 2.4679587339540343;
  const double c20 = 2.4218245962496954;
  const double c50 = 2.3747943989954163;
  const double c15 = 2.2654094725071228;
  const double c73 = 2.2323526687107482;
  const double c143 = 2.1875;
  const double c119 = 2.1480823570105501;
  const double c16 = 2.0973616238622275;
  const double c37 = 2.0566322782950288;
  const double c131 = 1.8944305707784594;
  const double c133 = 1.7320508075688772;
  const double c106 = 1.5189236011893423;
  const double c78 = 1.4882351124738322;
  const double c57 = 1.4248766393972496;
  const double c90 = 1.28884941420633;
  const double c125 = 1.2401959270615268;
  const double c21 = 1.2109122981248477;
  const double c51 = 1.1873971994977082;
  const double c4 = 1.1210869446657561;
  const double c71 = 1.1161763343553741;
  const double c144 = 1.09375;
  const double c116 = 1.074041178505275;
  const double c127 = 0.8660254037844386;
  const double c92 = 0.85923294280422002;
  const double c40 = 0.82265291131801144;
  const double c72 = 0.74411755623691611;
  const double c118 = 0.71602745233684995;
  const double c8 = 0.64725984928774938;
  const double c91 = 0.64442470710316502;
  const double c126 = 0.62009796353076341;
  const double c0 = 0.56054347233287805;
  const double c137 = 0.5;
  const double c47 = 0.47495887979908324;
  const double c88 = 0.42961647140211001;
  const double c123 = 0.41339864235384227;
  const double c35 = 0.41132645565900572;
  const double c70 = 0.37205877811845806;
  const double c115 = 0.35801372616842497;
  const double c9 = 0.32362992464387469;
  const double c48 = 0.23747943989954162;
  const double c89 = 0.21480823570105501;
  const double c124 = 0.20669932117692114;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 75, source += 216) {
    target[0] =  c0 * source[0] - c0 * source[2] - c1 * source[12]
                  + c1 * source[14] + c2 * source[24] - c2 * source[26]
                  - c3 * source[36] + c3 * source[38];
    target[1] =  c4 * source[1] - c5 * source[13] + c6 * source[25]
                  - c7 * source[37];
    target[2] =  c4 * source[3] - c5 * source[15] + c6 * source[27]
                  - c7 * source[39];
    target[3] =  c4 * source[4] - c5 * source[16] + c6 * source[28]
                  - c7 * source[40];
    target[4] =  c8 * source[5] - c9 * source[0] - c9 * source[2]
                  - c10 * source[17] + c11 * source[12] + c11 * source[14]
                  + c12 * source[29] - c13 * source[24] - c13 * source[26]
                  - c14 * source[41] + c15 * source[36] + c15 * source[38];
    target[5] =  c3 * source[6] - c3 * source[8] - c2 * source[18]
                  + c2 * source[20] + c1 * source[30] - c1 * source[32]
                  - c0 * source[42] + c0 * source[44];
    target[6] =  c7 * source[7] - c6 * source[19] + c5 * source[31]
                  - c4 * source[43];
    target[7] =  c7 * source[9] - c6 * source[21] + c5 * source[33]
                  - c4 * source[45];
    target[8] =  c7 * source[10] - c6 * source[22] + c5 * source[34]
                  - c4 * source[46];
    target[9] =  c14 * source[11] - c15 * source[6] - c15 * source[8]
                  - c12 * source[23] + c13 * source[18] + c13 * source[20]
                  + c10 * source[35] - c11 * source[30] - c11 * source[32]
                  - c8 * source[47] + c9 * source[42] + c9 * source[44];
    target[10] =  c16 * source[48] - c16 * source[50] - c17 * source[60]
                  + c17 * source[62] + c17 * source[72] - c17 * source[74]
                  - c16 * source[84] + c16 * source[86];
    target[11] =  c18 * source[49] - c19 * source[61] + c19 * source[73]
                  - c18 * source[85];
    target[12] =  c18 * source[51] - c19 * source[63] + c19 * source[75]
                  - c18 * source[87];
    target[13] =  c18 * source[52] - c19 * source[64] + c19 * source[76]
                  - c18 * source[88];
    target[14] =  c20 * source[53] - c21 * source[48] - c21 * source[50]
                  - c22 * source[65] + c23 * source[60] + c23 * source[62]
                  + c22 * source[77] - c23 * source[72] - c23 * source[74]
                  - c20 * source[89] + c21 * source[84] + c21 * source[86];
    target[15] =  c24 * source[54] - c24 * source[56] - c25 * source[66]
                  + c25 * source[68] + c24 * source[78] - c24 * source[80];
    target[16] =  c26 * source[55] - c27 * source[67] + c26 * source[79];
    target[17] =  c26 * source[57] - c27 * source[69] + c26 * source[81];
    target[18] =  c26 * source[58] - c27 * source[70] + c26 * source[82];
    target[19] =  c28 * source[59] - c29 * source[54] - c29 * source[56]
                  - c30 * source[71] + c31 * source[66] + c31 * source[68]
                  + c28 * source[83] - c29 * source[78] - c29 * source[80];
    target[20] =  c32 * source[90] - c32 * source[92] - c33 * source[102]
                  + c33 * source[104] + c34 * source[114] - c34 * source[116]
                  - c35 * source[0] + c35 * source[2] + c36 * source[12]
                  - c36 * source[14] - c37 * source[24] + c37 * source[26]
                  - c35 * source[12] + c35 * source[14] + c36 * source[24]
                  - c36 * source[26] - c37 * source[36] + c37 * source[38];
    target[21] =  c38 * source[91] - c39 * source[103] + c33 * source[115]
                  - c40 * source[1] + c41 * source[13] - c36 * source[25]
                  - c40 * source[13] + c41 * source[25] - c36 * source[37];
    target[22] =  c38 * source[93] - c39 * source[105] + c33 * source[117]
                  - c40 * source[3] + c41 * source[15] - c36 * source[27]
                  - c40 * source[15] + c41 * source[27] - c36 * source[39];
    target[23] =  c38 * source[94] - c39 * source[106] + c33 * source[118]
                  - c40 * source[4] + c41 * source[16] - c36 * source[28]
                  - c40 * source[16] + c41 * source[28] - c36 * source[40];
    target[24] =  c42 * source[95] - c43 * source[90] - c43 * source[92]
                  - c44 * source[107] + c45 * source[102] + c45 * source[104]
                  + c45 * source[119] - c46 * source[114] - c46 * source[116]
                  - c47 * source[5] + c48 * source[0] + c48 * source[2]
                  + c49 * source[17] - c50 * source[12] - c50 * source[14]
                  - c50 * source[29] + c51 * source[24] + c51 * source[26]
                  - c47 * source[17] + c48 * source[12] + c48 * source[14]
                  + c49 * source[29] - c50 * source[24] - c50 * source[26]
                  - c50 * source[41] + c51 * source[36] + c51 * source[38];
    target[25] =  c34 * source[96] - c34 * source[98] - c33 * source[108]
                  + c33 * source[110] + c32 * source[120] - c32 * source[122]
                  - c37 * source[6] + c37 * source[8] + c36 * source[18]
                  - c36 * source[20] - c35 * source[30] + c35 * source[32]
                  - c37 * source[18] + c37 * source[20] + c36 * source[30]
                  - c36 * source[32] - c35 * source[42] + c35 * source[44];
    target[26] =  c33 * source[97] - c39 * source[109] + c38 * source[121]
                  - c36 * source[7] + c41 * source[19] - c40 * source[31]
                  - c36 * source[19] + c41 * source[31] - c40 * source[43];
    target[27] =  c33 * source[99] - c39 * source[111] + c38 * source[123]
                  - c36 * source[9] + c41 * source[21] - c40 * source[33]
                  - c36 * source[21] + c41 * source[33] - c40 * source[45];
    target[28] =  c33 * source[100] - c39 * source[112] + c38 * source[124]
                  - c36 * source[10] + c41 * source[22] - c40 * source[34]
                  - c36 * source[22] + c41 * source[34] - c40 * source[46];
    target[29] =  c45 * source[101] - c46 * source[96] - c46 * source[98]
                  - c44 * source[113] + c45 * source[108] + c45 * source[110]
                  + c42 * source[125] - c43 * source[120] - c43 * source[122]
                  - c50 * source[11] + c51 * source[6] + c51 * source[8]
                  + c49 * source[23] - c50 * source[18] - c50 * source[20]
                  - c47 * source[35] + c48 * source[30] + c48 * source[32]
                  - c50 * source[23] + c51 * source[18] + c51 * source[20]
                  + c49 * source[35] - c50 * source[30] - c50 * source[32]
                  - c47 * source[47] + c48 * source[42] + c48 * source[44];
    target[30] =  c41 * source[126] - c41 * source[128] - c33 * source[138]
                  + c33 * source[140] + c41 * source[150] - c41 * source[152]
                  - c52 * source[48] + c52 * source[50] + c53 * source[60]
                  - c53 * source[62] - c52 * source[72] + c52 * source[74]
                  - c52 * source[60] + c52 * source[62] + c53 * source[72]
                  - c53 * source[74] - c52 * source[84] + c52 * source[86];
    target[31] =  c54 * source[127] - c39 * source[139] + c54 * source[151]
                  - c32 * source[49] + c55 * source[61] - c32 * source[73]
                  - c32 * source[61] + c55 * source[73] - c32 * source[85];
    target[32] =  c54 * source[129] - c39 * source[141] + c54 * source[153]
                  - c32 * source[51] + c55 * source[63] - c32 * source[75]
                  - c32 * source[63] + c55 * source[75] - c32 * source[87];
    target[33] =  c54 * source[130] - c39 * source[142] + c54 * source[154]
                  - c32 * source[52] + c55 * source[64] - c32 * source[76]
                  - c32 * source[64] + c55 * source[76] - c32 * source[88];
    target[34] =  c56 * source[131] - c49 * source[126] - c49 * source[128]
                  - c44 * source[143] + c45 * source[138] + c45 * source[140]
                  + c56 * source[155] - c49 * source[150] - c49 * source[152]
                  - c43 * source[53] + c57 * source[48] + c57 * source[50]
                  + c58 * source[65] - c59 * source[60] - c59 * source[62]
                  - c43 * source[77] + c57 * source[72] + c57 * source[74]
                  - c43 * source[65] + c57 * source[60] + c57 * source[62]
                  + c58 * source[77] - c59 * source[72] - c59 * source[74]
                  - c43 * source[89] + c57 * source[84] + c57 * source[86];
    target[35] =  c60 * source[132] - c60 * source[134] - c60 * source[144]
                  + c60 * source[146] - c38 * source[54] + c38 * source[56]
                  + c38 * source[66] - c38 * source[68] - c38 * source[66]
                  + c38 * source[68] + c38 * source[78] - c38 * source[80];
    target[36] =  c61 * source[133] - c61 * source[145] - c62 * source[55]
                  + c62 * source[67] - c62 * source[67] + c62 * source[79];
    target[37] =  c61 * source[135] - c61 * source[147] - c62 * source[57]
                  + c62 * source[69] - c62 * source[69] + c62 * source[81];
    target[38] =  c61 * source[136] - c61 * source[148] - c62 * source[58]
                  + c62 * source[70] - c62 * source[70] + c62 * source[82];
    target[39] =  c63 * source[137] - c64 * source[132] - c64 * source[134]
                  - c63 * source[149] + c64 * source[144] + c64 * source[146]
                  - c65 * source[59] + c42 * source[54] + c42 * source[56]
                  + c65 * source[71] - c42 * source[66] - c42 * source[68]
                  - c65 * source[71] + c42 * source[66] + c42 * source[68]
                  + c65 * source[83] - c42 * source[78] - c42 * source[80];
    target[40] =  c66 * source[156] - c66 * source[158] - c67 * source[168]
                  + c67 * source[170] - c68 * source[90] + c68 * source[92]
                  + c69 * source[102] - c69 * source[104] - c68 * source[102]
                  + c68 * source[104] + c69 * source[114] - c69 * source[116]
                  + c70 * source[0] - c70 * source[2] - c71 * source[12]
                  + c71 * source[14] + c72 * source[12] - c72 * source[14]
                  - c73 * source[24] + c73 * source[26] + c70 * source[24]
                  - c70 * source[26] - c71 * source[36] + c71 * source[38];
    target[41] =  c74 * source[157] - c75 * source[169] - c76 * source[91]
                  + c77 * source[103] - c76 * source[103] + c77 * source[115]
                  + c72 * source[1] - c73 * source[13] + c78 * source[13]
                  - c79 * source[25] + c72 * source[25] - c73 * source[37];
    target[42] =  c74 * source[159] - c75 * source[171] - c76 * source[93]
                  + c77 * source[105] - c76 * source[105] + c77 * source[117]
                  + c72 * source[3] - c73 * source[15] + c78 * source[15]
                  - c79 * source[27] + c72 * source[27] - c73 * source[39];
    target[43] =  c74 * source[160] - c75 * source[172] - c76 * source[94]
                  + c77 * source[106] - c76 * source[106] + c77 * source[118]
                  + c72 * source[4] - c73 * source[16] + c78 * source[16]
                  - c79 * source[28] + c72 * source[28] - c73 * source[40];
    target[44] =  c80 * source[161] - c81 * source[156] - c81 * source[158]
                  - c82 * source[173] + c83 * source[168] + c83 * source[170]
                  - c84 * source[95] + c85 * source[90] + c85 * source[92]
                  + c86 * source[107] - c87 * source[102] - c87 * source[104]
                  - c84 * source[107] + c85 * source[102] + c85 * source[104]
                  + c86 * source[119] - c87 * source[114] - c87 * source[116]
                  + c88 * source[5] - c89 * source[0] - c89 * source[2]
                  - c90 * source[17] + c91 * source[12] + c91 * source[14]
                  + c92 * source[17] - c88 * source[12] - c88 * source[14]
                  - c93 * source[29] + c90 * source[24] + c90 * source[26]
                  + c88 * source[29] - c89 * source[24] - c89 * source[26]
                  - c90 * source[41] + c91 * source[36] + c91 * source[38];
    target[45] =  c67 * source[162] - c67 * source[164] - c66 * source[174]
                  + c66 * source[176] - c69 * source[96] + c69 * source[98]
                  + c68 * source[108] - c68 * source[110] - c69 * source[108]
                  + c69 * source[110] + c68 * source[120] - c68 * source[122]
                  + c71 * source[6] - c71 * source[8] - c70 * source[18]
                  + c70 * source[20] + c73 * source[18] - c73 * source[20]
                  - c72 * source[30] + c72 * source[32] + c71 * source[30]
                  - c71 * source[32] - c70 * source[42] + c70 * source[44];
    target[46] =  c75 * source[163] - c74 * source[175] - c77 * source[97]
                  + c76 * source[109] - c77 * source[109] + c76 * source[121]
                  + c73 * source[7] - c72 * source[19] + c79 * source[19]
                  - c78 * source[31] + c73 * source[31] - c72 * source[43];
    target[47] =  c75 * source[165] - c74 * source[177] - c77 * source[99]
                  + c76 * source[111] - c77 * source[111] + c76 * source[123]
                  + c73 * source[9] - c72 * source[21] + c79 * source[21]
                  - c78 * source[33] + c73 * source[33] - c72 * source[45];
    target[48] =  c75 * source[166] - c74 * source[178] - c77 * source[100]
                  + c76 * source[112] - c77 * source[112] + c76 * source[124]
                  + c73 * source[10] - c72 * source[22] + c79 * source[22]
                  - c78 * source[34] + c73 * source[34] - c72 * source[46];
    target[49] =  c82 * source[167] - c83 * source[162] - c83 * source[164]
                  - c80 * source[179] + c81 * source[174] + c81 * source[176]
                  - c86 * source[101] + c87 * source[96] + c87 * source[98]
                  + c84 * source[113] - c85 * source[108] - c85 * source[110]
                  - c86 * source[113] + c87 * source[108] + c87 * source[110]
                  + c84 * source[125] - c85 * source[120] - c85 * source[122]
                  + c90 * source[11] - c91 * source[6] - c91 * source[8]
                  - c88 * source[23] + c89 * source[18] + c89 * source[20]
                  + c93 * source[23] - c90 * source[18] - c90 * source[20]
                  - c92 * source[35] + c88 * source[30] + c88 * source[32]
                  + c90 * source[35] - c91 * source[30] - c91 * source[32]
                  - c88 * source[47] + c89 * source[42] + c89 * source[44];
    target[50] =  c94 * source[180] - c94 * source[182] - c94 * source[192]
                  + c94 * source[194] - c95 * source[126] + c95 * source[128]
                  + c95 * source[138] - c95 * source[140] - c95 * source[138]
                  + c95 * source[140] + c95 * source[150] - c95 * source[152]
                  + c96 * source[48] - c96 * source[50] - c96 * source[60]
                  + c96 * source[62] + c97 * source[60] - c97 * source[62]
                  - c97 * source[72] + c97 * source[74] + c96 * source[72]
                  - c96 * source[74] - c96 * source[84] + c96 * source[86];
    target[51] =  c98 * source[181] - c98 * source[193] - c99 * source[127]
                  + c99 * source[139] - c99 * source[139] + c99 * source[151]
                  + c97 * source[49] - c97 * source[61] + c100 * source[61]
                  - c100 * source[73] + c97 * source[73] - c97 * source[85];
    target[52] =  c98 * source[183] - c98 * source[195] - c99 * source[129]
                  + c99 * source[141] - c99 * source[141] + c99 * source[153]
                  + c97 * source[51] - c97 * source[63] + c100 * source[63]
                  - c100 * source[75] + c97 * source[75] - c97 * source[87];
    target[53] =  c98 * source[184] - c98 * source[196] - c99 * source[130]
                  + c99 * source[142] - c99 * source[142] + c99 * source[154]
                  + c97 * source[52] - c97 * source[64] + c100 * source[64]
                  - c100 * source[76] + c97 * source[76] - c97 * source[88];
    target[54] =  c101 * source[185] - c102 * source[180] - c102 * source[182]
                  - c101 * source[197] + c102 * source[192] + c102 * source[194]
                  - c103 * source[131] + c104 * source[126] + c104 * source[128]
                  + c103 * source[143] - c104 * source[138] - c104 * source[140]
                  - c103 * source[143] + c104 * source[138] + c104 * source[140]
                  + c103 * source[155] - c104 * source[150] - c104 * source[152]
                  + c105 * source[53] - c106 * source[48] - c106 * source[50]
                  - c105 * source[65] + c106 * source[60] + c106 * source[62]
                  + c107 * source[65] - c105 * source[60] - c105 * source[62]
                  - c107 * source[77] + c105 * source[72] + c105 * source[74]
                  + c105 * source[77] - c106 * source[72] - c106 * source[74]
                  - c105 * source[89] + c106 * source[84] + c106 * source[86];
    target[55] =  c98 * source[186] - c98 * source[188] - c99 * source[132]
                  + c99 * source[134] - c99 * source[144] + c99 * source[146]
                  + c97 * source[54] - c97 * source[56] + c100 * source[66]
                  - c100 * source[68] + c97 * source[78] - c97 * source[80];
    target[56] =  c108 * source[187] - c109 * source[133] - c109 * source[145]
                  + c100 * source[55] + c110 * source[67] + c100 * source[79];
    target[57] =  c108 * source[189] - c109 * source[135] - c109 * source[147]
                  + c100 * source[57] + c110 * source[69] + c100 * source[81];
    target[58] =  c108 * source[190] - c109 * source[136] - c109 * source[148]
                  + c100 * source[58] + c110 * source[70] + c100 * source[82];
    target[59] =  c111 * source[191] - c101 * source[186] - c101 * source[188]
                  - c112 * source[137] + c103 * source[132] + c103 * source[134]
                  - c112 * source[149] + c103 * source[144] + c103 * source[146]
                  + c107 * source[59] - c105 * source[54] - c105 * source[56]
                  + c113 * source[71] - c107 * source[66] - c107 * source[68]
                  + c107 * source[83] - c105 * source[78] - c105 * source[80];
    target[60] =  c114 * source[198] - c114 * source[200] - c83 * source[156]
                  + c83 * source[158] - c83 * source[168] + c83 * source[170]
                  + c84 * source[90] - c84 * source[92] + c83 * source[102]
                  - c83 * source[104] + c84 * source[114] - c84 * source[116]
                  - c115 * source[0] + c115 * source[2] - c116 * source[12]
                  + c116 * source[14] - c116 * source[24] + c116 * source[26]
                  - c115 * source[36] + c115 * source[38];
    target[61] =  c117 * source[199] - c82 * source[157] - c82 * source[169]
                  + c83 * source[91] + c82 * source[103] + c83 * source[115]
                  - c118 * source[1] - c119 * source[13] - c119 * source[25]
                  - c118 * source[37];
    target[62] =  c117 * source[201] - c82 * source[159] - c82 * source[171]
                  + c83 * source[93] + c82 * source[105] + c83 * source[117]
                  - c118 * source[3] - c119 * source[15] - c119 * source[27]
                  - c118 * source[39];
    target[63] =  c117 * source[202] - c82 * source[160] - c82 * source[172]
                  + c83 * source[94] + c82 * source[106] + c83 * source[118]
                  - c118 * source[4] - c119 * source[16] - c119 * source[28]
                  - c118 * source[40];
    target[64] =  c120 * source[203] - c121 * source[198] - c121 * source[200]
                  - c74 * source[161] + c66 * source[156] + c66 * source[158]
                  - c74 * source[173] + c66 * source[168] + c66 * source[170]
                  + c66 * source[95] - c122 * source[90] - c122 * source[92]
                  + c74 * source[107] - c66 * source[102] - c66 * source[104]
                  + c66 * source[119] - c122 * source[114] - c122 * source[116]
                  - c123 * source[5] + c124 * source[0] + c124 * source[2]
                  - c125 * source[17] + c126 * source[12] + c126 * source[14]
                  - c125 * source[29] + c126 * source[24] + c126 * source[26]
                  - c123 * source[41] + c124 * source[36] + c124 * source[38];
    target[65] =  c114 * source[204] - c114 * source[206] - c83 * source[162]
                  + c83 * source[164] - c83 * source[174] + c83 * source[176]
                  + c84 * source[96] - c84 * source[98] + c83 * source[108]
                  - c83 * source[110] + c84 * source[120] - c84 * source[122]
                  - c115 * source[6] + c115 * source[8] - c116 * source[18]
                  + c116 * source[20] - c116 * source[30] + c116 * source[32]
                  - c115 * source[42] + c115 * source[44];
    target[66] =  c117 * source[205] - c82 * source[163] - c82 * source[175]
                  + c83 * source[97] + c82 * source[109] + c83 * source[121]
                  - c118 * source[7] - c119 * source[19] - c119 * source[31]
                  - c118 * source[43];
    target[67] =  c117 * source[207] - c82 * source[165] - c82 * source[177]
                  + c83 * source[99] + c82 * source[111] + c83 * source[123]
                  - c118 * source[9] - c119 * source[21] - c119 * source[33]
                  - c118 * source[45];
    target[68] =  c117 * source[208] - c82 * source[166] - c82 * source[178]
                  + c83 * source[100] + c82 * source[112] + c83 * source[124]
                  - c118 * source[10] - c119 * source[22] - c119 * source[34]
                  - c118 * source[46];
    target[69] =  c120 * source[209] - c121 * source[204] - c121 * source[206]
                  - c74 * source[167] + c66 * source[162] + c66 * source[164]
                  - c74 * source[179] + c66 * source[174] + c66 * source[176]
                  + c66 * source[101] - c122 * source[96] - c122 * source[98]
                  + c74 * source[113] - c66 * source[108] - c66 * source[110]
                  + c66 * source[125] - c122 * source[120] - c122 * source[122]
                  - c123 * source[11] + c124 * source[6] + c124 * source[8]
                  - c125 * source[23] + c126 * source[18] + c126 * source[20]
                  - c125 * source[35] + c126 * source[30] + c126 * source[32]
                  - c123 * source[47] + c124 * source[42] + c124 * source[44];
    target[70] =  c127 * source[210] - c127 * source[212] - c128 * source[180]
                  + c128 * source[182] - c128 * source[192] + c128 * source[194]
                  + c129 * source[126] - c129 * source[128] + c130 * source[138]
                  - c130 * source[140] + c129 * source[150] - c129 * source[152]
                  - c131 * source[48] + c131 * source[50] - c132 * source[60]
                  + c132 * source[62] - c132 * source[72] + c132 * source[74]
                  - c131 * source[84] + c131 * source[86];
    target[71] =  c133 * source[211] - c134 * source[181] - c134 * source[193]
                  + c130 * source[127] + c135 * source[139] + c130 * source[151]
                  - c136 * source[49] - c129 * source[61] - c129 * source[73]
                  - c136 * source[85];
    target[72] =  c133 * source[213] - c134 * source[183] - c134 * source[195]
                  + c130 * source[129] + c135 * source[141] + c130 * source[153]
                  - c136 * source[51] - c129 * source[63] - c129 * source[75]
                  - c136 * source[87];
    target[73] =  c133 * source[214] - c134 * source[184] - c134 * source[196]
                  + c130 * source[130] + c135 * source[142] + c130 * source[154]
                  - c136 * source[52] - c129 * source[64] - c129 * source[76]
                  - c136 * source[88];
    target[74] =  source[215] - c137 * source[210] - c137 * source[212]
                  - c138 * source[185] + c139 * source[180] + c139 * source[182]
                  - c138 * source[197] + c139 * source[192] + c139 * source[194]
                  + c140 * source[131] - c141 * source[126] - c141 * source[128]
                  + c142 * source[143] - c140 * source[138] - c140 * source[140]
                  + c140 * source[155] - c141 * source[150] - c141 * source[152]
                  - c143 * source[53] + c144 * source[48] + c144 * source[50]
                  - c141 * source[65] + c145 * source[60] + c145 * source[62]
                  - c141 * source[77] + c145 * source[72] + c145 * source[74]
                  - c143 * source[89] + c144 * source[84] + c144 * source[86];
  }
}

#endif
