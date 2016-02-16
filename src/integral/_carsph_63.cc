//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_63.cc
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


void CarSphList::carsph_63(const int nloop, const double* source, double* target) {
  const double c68 = 115.27819611704548;
  const double c42 = 90.117111305234374;
  const double c110 = 84.187291202413675;
  const double c93 = 76.852130744696993;
  const double c73 = 72.908332857088425;
  const double c58 = 70.593189738812626;
  const double c64 = 57.639098058522741;
  const double c45 = 56.995065575889988;
  const double c147 = 56.124860801609124;
  const double c35 = 55.185234947855392;
  const double c115 = 53.244718047896548;
  const double c23 = 52.029138470668528;
  const double c101 = 51.553976568253198;
  const double c95 = 48.605555238058955;
  const double c88 = 47.062126492541751;
  const double c39 = 45.058555652617187;
  const double c82 = 44.647053374214963;
  const double c163 = 44.37059837324712;
  const double c189 = 43.571062644833439;
  const double c106 = 42.093645601206838;
  const double c7 = 39.021853853001396;
  const double c91 = 38.426065372348496;
  const double c149 = 35.496478698597699;
  const double c52 = 34.902207129349286;
  const double c144 = 34.369317712168801;
  const double c26 = 32.906116452720461;
  const double c124 = 32.605597678926237;
  const double c19 = 31.861210252437054;
  const double c112 = 31.57023420090513;
  const double c81 = 29.764702249476645;
  const double c188 = 29.047375096555626;
  const double c47 = 28.497532787944994;
  const double c109 = 28.062430400804562;
  const double c37 = 27.592617473927696;
  const double c198 = 27.556759606310752;
  const double c155 = 27.171331399105199;
  const double c177 = 26.6817177576707;
  const double c21 = 26.014569235334264;
  const double c10 = 24.679587339540344;
  const double c3 = 23.895907689327789;
  const double c57 = 23.531063246270875;
  const double c51 = 23.268138086232856;
  const double c40 = 22.529277826308594;
  const double c160 = 22.18529918662356;
  const double c184 = 21.78553132241672;
  const double c123 = 21.737065119284157;
  const double c31 = 20.150798681441884;
  const double c119 = 19.966769267961205;
  const double c97 = 19.843134832984429;
  const double c5 = 19.510926926500698;
  const double c104 = 19.332741213094948;
  const double c67 = 19.213032686174248;
  const double c34 = 18.395078315951796;
  const double c194 = 18.371173070873837;
  const double c74 = 18.227083214272106;
  const double c173 = 17.787811838447134;
  const double c113 = 17.748239349298849;
  const double c54 = 17.451103564674643;
  const double c100 = 17.1846588560844;
  const double c209 = 16.875;
  const double c108 = 15.785117100452565;
  const double c87 = 15.687375497513916;
  const double c22 = 15.608741541200558;
  const double c15 = 15.113099011081413;
  const double c182 = 14.523687548277813;
  const double c151 = 14.491376746189438;
  const double c46 = 14.248766393972497;
  const double c105 = 14.031215200402281;
  const double c196 = 13.778379803155376;
  const double c159 = 13.5856656995526;
  const double c30 = 13.433865787627923;
  const double c175 = 13.34085887883535;
  const double c116 = 13.311179511974137;
  const double c128 = 12.227099129597338;
  const double c71 = 12.151388809514739;
  const double c56 = 11.765531623135438;
  const double c53 = 11.634069043116428;
  const double c70 = 11.527819611704547;
  const double c143 = 11.456439237389599;
  const double c206 = 11.25;
  const double c164 = 11.224972160321824;
  const double c162 = 11.09264959331178;
  const double c183 = 10.89276566120836;
  const double c122 = 10.868532559642079;
  const double c18 = 10.620403417479018;
  const double c111 = 10.523411400301709;
  const double c14 = 10.075399340720942;
  const double c24 = 9.8718349358161372;
  const double c63 = 9.6065163430871241;
  const double c17 = 9.5583630757311155;
  const double c36 = 9.197539157975898;
  const double c154 = 9.0571104663683997;
  const double c41 = 9.0117111305234374;
  const double c176 = 8.8939059192235668;
  const double c150 = 8.8741196746494246;
  const double c169 = 8.5923294280422002;
  const double c208 = 8.4375;
  const double c27 = 8.2265291131801153;
  const double c127 = 8.1513994197315593;
  const double c2 = 7.9653025631092635;
  const double c20 = 7.8043707706002792;
  const double c94 = 7.6852130744696989;
  const double c205 = 7.5;
  const double c80 = 7.4411755623691613;
  const double c77 = 7.2908332857088425;
  const double c121 = 7.245688373094719;
  const double c48 = 7.1243831969862486;
  const double c62 = 7.0593189738812621;
  const double c148 = 7.0156076002011405;
  const double c199 = 6.8891899015776881;
  const double c167 = 6.8738635424337602;
  const double c157 = 6.7928328497762998;
  const double c117 = 6.6555897559870685;
  const double c103 = 6.4442470710316497;
  const double c11 = 6.169896834885086;
  const double c29 = 6.0452396044325658;
  const double c172 = 5.9292706128157109;
  const double c66 = 5.7639098058522737;
  const double c99 = 5.7282196186947996;
  const double c43 = 5.6995065575889985;
  const double c207 = 5.625;
  const double c161 = 5.54632479665589;
  const double c33 = 5.5185234947855388;
  const double c107 = 5.2617057001508547;
  const double c120 = 4.9916923169903011;
  const double c79 = 4.9607837082461073;
  const double c96 = 4.8605555238058953;
  const double c90 = 4.7062126492541747;
  const double c195 = 4.5927932677184593;
  const double c166 = 4.5825756949558398;
  const double c158 = 4.5285552331841998;
  const double c38 = 4.5058555652617187;
  const double c86 = 4.4647053374214964;
  const double c174 = 4.4469529596117834;
  const double c114 = 4.4370598373247123;
  const double c146 = 4.2961647140211001;
  const double c126 = 4.0756997098657797;
  const double c28 = 4.0301597362883772;
  const double c55 = 3.9218438743784789;
  const double c187 = 3.872983346207417;
  const double c92 = 3.8426065372348495;
  const double c191 = 3.6309218870694533;
  const double c153 = 3.6228441865473595;
  const double c134 = 3.5078038001005702;
  const double c50 = 3.4902207129349283;
  const double c197 = 3.444594950788844;
  const double c16 = 3.1861210252437053;
  const double c72 = 3.0378472023786847;
  const double c85 = 2.9764702249476644;
  const double c168 = 2.8641098093473998;
  const double c165 = 2.8062430400804561;
  const double c125 = 2.7171331399105196;
  const double c6 = 2.6014569235334264;
  const double c25 = 2.4679587339540343;
  const double c192 = 2.4494897427831779;
  const double c171 = 2.3717082451262845;
  const double c61 = 2.3531063246270874;
  const double c49 = 2.3268138086232857;
  const double c202 = 2.2963966338592297;
  const double c156 = 2.2642776165920999;
  const double c180 = 2.2234764798058917;
  const double c137 = 2.2185299186623562;
  const double c102 = 2.1480823570105501;
  const double c98 = 1.984313483298443;
  const double c181 = 1.9364916731037085;
  const double c69 = 1.9213032686174247;
  const double c32 = 1.8395078315951796;
  const double c78 = 1.8227083214272106;
  const double c186 = 1.8154609435347266;
  const double c152 = 1.8114220932736798;
  const double c133 = 1.7539019000502851;
  const double c118 = 1.6638974389967671;
  const double c8 = 1.6453058226360229;
  const double c1 = 1.5930605126218527;
  const double c89 = 1.5687375497513916;
  const double c204 = 1.5;
  const double c145 = 1.4320549046736999;
  const double c44 = 1.4248766393972496;
  const double c213 = 1.40625;
  const double c142 = 1.3585665699552598;
  const double c4 = 1.3007284617667132;
  const double c75 = 1.2151388809514738;
  const double c190 = 1.2103072956898178;
  const double c60 = 1.1765531623135437;
  const double c135 = 1.1092649593311781;
  const double c130 = 1.074041178505275;
  const double c13 = 1.0075399340720943;
  const double c65 = 0.96065163430871237;
  const double c212 = 0.9375;
  const double c141 = 0.90571104663683988;
  const double c132 = 0.87695095002514256;
  const double c170 = 0.79056941504209488;
  const double c200 = 0.76546554461974314;
  const double c84 = 0.74411755623691611;
  const double c179 = 0.74115882660196386;
  const double c131 = 0.71602745233684995;
  const double c140 = 0.67928328497762991;
  const double c12 = 0.67169328938139616;
  const double c193 = 0.61237243569579447;
  const double c185 = 0.60515364784490888;
  const double c203 = 0.57409915846480741;
  const double c138 = 0.55463247966558904;
  const double c0 = 0.53102017087395093;
  const double c83 = 0.49607837082461076;
  const double c211 = 0.46875;
  const double c139 = 0.45285552331841994;
  const double c9 = 0.41132645565900572;
  const double c59 = 0.39218438743784789;
  const double c129 = 0.35801372616842497;
  const double c210 = 0.3125;
  const double c76 = 0.30378472023786846;
  const double c136 = 0.27731623983279452;
  const double c178 = 0.24705294220065463;
  const double c201 = 0.19136638615493579;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 91, source += 280) {
    target[0] =  c0 * source[0] - c1 * source[2] - c2 * source[20]
                  + c3 * source[22] + c2 * source[40] - c3 * source[42]
                  - c0 * source[60] + c1 * source[62];
    target[1] =  c1 * source[1] - c0 * source[3] - c3 * source[21]
                  + c2 * source[23] + c3 * source[41] - c2 * source[43]
                  - c1 * source[61] + c0 * source[63];
    target[2] =  c4 * source[4] - c4 * source[6] - c5 * source[24]
                  + c5 * source[26] + c5 * source[44] - c5 * source[46]
                  - c4 * source[64] + c4 * source[66];
    target[3] =  c6 * source[5] - c7 * source[25] + c7 * source[45]
                  - c6 * source[65];
    target[4] =  c8 * source[7] - c9 * source[0] - c9 * source[2]
                  - c10 * source[27] + c11 * source[20] + c11 * source[22]
                  + c10 * source[47] - c11 * source[40] - c11 * source[42]
                  - c8 * source[67] + c9 * source[60] + c9 * source[62];
    target[5] =  c8 * source[8] - c9 * source[1] - c9 * source[3]
                  - c10 * source[28] + c11 * source[21] + c11 * source[23]
                  + c10 * source[48] - c11 * source[41] - c11 * source[43]
                  - c8 * source[68] + c9 * source[61] + c9 * source[63];
    target[6] =  c12 * source[9] - c13 * source[4] - c13 * source[6]
                  - c14 * source[29] + c15 * source[24] + c15 * source[26]
                  + c14 * source[49] - c15 * source[44] - c15 * source[46]
                  - c12 * source[69] + c13 * source[64] + c13 * source[66];
    target[7] =  c16 * source[10] - c17 * source[12] - c18 * source[30]
                  + c19 * source[32] + c16 * source[50] - c17 * source[52];
    target[8] =  c17 * source[11] - c16 * source[13] - c19 * source[31]
                  + c18 * source[33] + c17 * source[51] - c16 * source[53];
    target[9] =  c20 * source[14] - c20 * source[16] - c21 * source[34]
                  + c21 * source[36] + c20 * source[54] - c20 * source[56];
    target[10] =  c22 * source[15] - c23 * source[35] + c22 * source[55];
    target[11] =  c24 * source[17] - c25 * source[10] - c25 * source[12]
                  - c26 * source[37] + c27 * source[30] + c27 * source[32]
                  + c24 * source[57] - c25 * source[50] - c25 * source[52];
    target[12] =  c24 * source[18] - c25 * source[11] - c25 * source[13]
                  - c26 * source[38] + c27 * source[31] + c27 * source[33]
                  + c24 * source[58] - c25 * source[51] - c25 * source[53];
    target[13] =  c28 * source[19] - c29 * source[14] - c29 * source[16]
                  - c30 * source[39] + c31 * source[34] + c31 * source[36]
                  + c28 * source[59] - c29 * source[54] - c29 * source[56];
    target[14] =  c32 * source[70] - c33 * source[72] - c34 * source[90]
                  + c35 * source[92] + c36 * source[110] - c37 * source[112];
    target[15] =  c33 * source[71] - c32 * source[73] - c35 * source[91]
                  + c34 * source[93] + c37 * source[111] - c36 * source[113];
    target[16] =  c38 * source[74] - c38 * source[76] - c39 * source[94]
                  + c39 * source[96] + c40 * source[114] - c40 * source[116];
    target[17] =  c41 * source[75] - c42 * source[95] + c39 * source[115];
    target[18] =  c43 * source[77] - c44 * source[70] - c44 * source[72]
                  - c45 * source[97] + c46 * source[90] + c46 * source[92]
                  + c47 * source[117] - c48 * source[110] - c48 * source[112];
    target[19] =  c43 * source[78] - c44 * source[71] - c44 * source[73]
                  - c45 * source[98] + c46 * source[91] + c46 * source[93]
                  + c47 * source[118] - c48 * source[111] - c48 * source[113];
    target[20] =  c49 * source[79] - c50 * source[74] - c50 * source[76]
                  - c51 * source[99] + c52 * source[94] + c52 * source[96]
                  + c53 * source[119] - c54 * source[114] - c54 * source[116];
    target[21] =  c36 * source[80] - c37 * source[82] - c34 * source[100]
                  + c35 * source[102] + c32 * source[120] - c33 * source[122];
    target[22] =  c37 * source[81] - c36 * source[83] - c35 * source[101]
                  + c34 * source[103] + c33 * source[121] - c32 * source[123];
    target[23] =  c40 * source[84] - c40 * source[86] - c39 * source[104]
                  + c39 * source[106] + c38 * source[124] - c38 * source[126];
    target[24] =  c39 * source[85] - c42 * source[105] + c41 * source[125];
    target[25] =  c47 * source[87] - c48 * source[80] - c48 * source[82]
                  - c45 * source[107] + c46 * source[100] + c46 * source[102]
                  + c43 * source[127] - c44 * source[120] - c44 * source[122];
    target[26] =  c47 * source[88] - c48 * source[81] - c48 * source[83]
                  - c45 * source[108] + c46 * source[101] + c46 * source[103]
                  + c43 * source[128] - c44 * source[121] - c44 * source[123];
    target[27] =  c53 * source[89] - c54 * source[84] - c54 * source[86]
                  - c51 * source[109] + c52 * source[104] + c52 * source[106]
                  + c49 * source[129] - c50 * source[124] - c50 * source[126];
    target[28] =  c55 * source[130] - c56 * source[132] - c57 * source[150]
                  + c58 * source[152] + c55 * source[170] - c56 * source[172]
                  - c59 * source[0] + c60 * source[2] + c61 * source[20]
                  - c62 * source[22] - c59 * source[40] + c60 * source[42]
                  - c59 * source[20] + c60 * source[22] + c61 * source[40]
                  - c62 * source[42] - c59 * source[60] + c60 * source[62];
    target[29] =  c56 * source[131] - c55 * source[133] - c58 * source[151]
                  + c57 * source[153] + c56 * source[171] - c55 * source[173]
                  - c60 * source[1] + c59 * source[3] + c62 * source[21]
                  - c61 * source[23] - c60 * source[41] + c59 * source[43]
                  - c60 * source[21] + c59 * source[23] + c62 * source[41]
                  - c61 * source[43] - c60 * source[61] + c59 * source[63];
    target[30] =  c63 * source[134] - c63 * source[136] - c64 * source[154]
                  + c64 * source[156] + c63 * source[174] - c63 * source[176]
                  - c65 * source[4] + c65 * source[6] + c66 * source[24]
                  - c66 * source[26] - c65 * source[44] + c65 * source[46]
                  - c65 * source[24] + c65 * source[26] + c66 * source[44]
                  - c66 * source[46] - c65 * source[64] + c65 * source[66];
    target[31] =  c67 * source[135] - c68 * source[155] + c67 * source[175]
                  - c69 * source[5] + c70 * source[25] - c69 * source[45]
                  - c69 * source[25] + c70 * source[45] - c69 * source[65];
    target[32] =  c71 * source[137] - c72 * source[130] - c72 * source[132]
                  - c73 * source[157] + c74 * source[150] + c74 * source[152]
                  + c71 * source[177] - c72 * source[170] - c72 * source[172]
                  - c75 * source[7] + c76 * source[0] + c76 * source[2]
                  + c77 * source[27] - c78 * source[20] - c78 * source[22]
                  - c75 * source[47] + c76 * source[40] + c76 * source[42]
                  - c75 * source[27] + c76 * source[20] + c76 * source[22]
                  + c77 * source[47] - c78 * source[40] - c78 * source[42]
                  - c75 * source[67] + c76 * source[60] + c76 * source[62];
    target[33] =  c71 * source[138] - c72 * source[131] - c72 * source[133]
                  - c73 * source[158] + c74 * source[151] + c74 * source[153]
                  + c71 * source[178] - c72 * source[171] - c72 * source[173]
                  - c75 * source[8] + c76 * source[1] + c76 * source[3]
                  + c77 * source[28] - c78 * source[21] - c78 * source[23]
                  - c75 * source[48] + c76 * source[41] + c76 * source[43]
                  - c75 * source[28] + c76 * source[21] + c76 * source[23]
                  + c77 * source[48] - c78 * source[41] - c78 * source[43]
                  - c75 * source[68] + c76 * source[61] + c76 * source[63];
    target[34] =  c79 * source[139] - c80 * source[134] - c80 * source[136]
                  - c81 * source[159] + c82 * source[154] + c82 * source[156]
                  + c79 * source[179] - c80 * source[174] - c80 * source[176]
                  - c83 * source[9] + c84 * source[4] + c84 * source[6]
                  + c85 * source[29] - c86 * source[24] - c86 * source[26]
                  - c83 * source[49] + c84 * source[44] + c84 * source[46]
                  - c83 * source[29] + c84 * source[24] + c84 * source[26]
                  + c85 * source[49] - c86 * source[44] - c86 * source[46]
                  - c83 * source[69] + c84 * source[64] + c84 * source[66];
    target[35] =  c87 * source[140] - c88 * source[142] - c87 * source[160]
                  + c88 * source[162] - c89 * source[10] + c90 * source[12]
                  + c89 * source[30] - c90 * source[32] - c89 * source[30]
                  + c90 * source[32] + c89 * source[50] - c90 * source[52];
    target[36] =  c88 * source[141] - c87 * source[143] - c88 * source[161]
                  + c87 * source[163] - c90 * source[11] + c89 * source[13]
                  + c90 * source[31] - c89 * source[33] - c90 * source[31]
                  + c89 * source[33] + c90 * source[51] - c89 * source[53];
    target[37] =  c91 * source[144] - c91 * source[146] - c91 * source[164]
                  + c91 * source[166] - c92 * source[14] + c92 * source[16]
                  + c92 * source[34] - c92 * source[36] - c92 * source[34]
                  + c92 * source[36] + c92 * source[54] - c92 * source[56];
    target[38] =  c93 * source[145] - c93 * source[165] - c94 * source[15]
                  + c94 * source[35] - c94 * source[35] + c94 * source[55];
    target[39] =  c95 * source[147] - c71 * source[140] - c71 * source[142]
                  - c95 * source[167] + c71 * source[160] + c71 * source[162]
                  - c96 * source[17] + c75 * source[10] + c75 * source[12]
                  + c96 * source[37] - c75 * source[30] - c75 * source[32]
                  - c96 * source[37] + c75 * source[30] + c75 * source[32]
                  + c96 * source[57] - c75 * source[50] - c75 * source[52];
    target[40] =  c95 * source[148] - c71 * source[141] - c71 * source[143]
                  - c95 * source[168] + c71 * source[161] + c71 * source[163]
                  - c96 * source[18] + c75 * source[11] + c75 * source[13]
                  + c96 * source[38] - c75 * source[31] - c75 * source[33]
                  - c96 * source[38] + c75 * source[31] + c75 * source[33]
                  + c96 * source[58] - c75 * source[51] - c75 * source[53];
    target[41] =  c97 * source[149] - c81 * source[144] - c81 * source[146]
                  - c97 * source[169] + c81 * source[164] + c81 * source[166]
                  - c98 * source[19] + c85 * source[14] + c85 * source[16]
                  + c98 * source[39] - c85 * source[34] - c85 * source[36]
                  - c98 * source[39] + c85 * source[34] + c85 * source[36]
                  + c98 * source[59] - c85 * source[54] - c85 * source[56];
    target[42] =  c99 * source[180] - c100 * source[182] - c100 * source[200]
                  + c101 * source[202] - c102 * source[70] + c103 * source[72]
                  + c103 * source[90] - c104 * source[92] - c102 * source[90]
                  + c103 * source[92] + c103 * source[110] - c104 * source[112];
    target[43] =  c100 * source[181] - c99 * source[183] - c101 * source[201]
                  + c100 * source[203] - c103 * source[71] + c102 * source[73]
                  + c104 * source[91] - c103 * source[93] - c103 * source[91]
                  + c102 * source[93] + c104 * source[111] - c103 * source[113];
    target[44] =  c105 * source[184] - c105 * source[186] - c106 * source[204]
                  + c106 * source[206] - c107 * source[74] + c107 * source[76]
                  + c108 * source[94] - c108 * source[96] - c107 * source[94]
                  + c107 * source[96] + c108 * source[114] - c108 * source[116];
    target[45] =  c109 * source[185] - c110 * source[205] - c111 * source[75]
                  + c112 * source[95] - c111 * source[95] + c112 * source[115];
    target[46] =  c113 * source[187] - c114 * source[180] - c114 * source[182]
                  - c115 * source[207] + c116 * source[200] + c116 * source[202]
                  - c117 * source[77] + c118 * source[70] + c118 * source[72]
                  + c119 * source[97] - c120 * source[90] - c120 * source[92]
                  - c117 * source[97] + c118 * source[90] + c118 * source[92]
                  + c119 * source[117] - c120 * source[110] - c120 * source[112];
    target[47] =  c113 * source[188] - c114 * source[181] - c114 * source[183]
                  - c115 * source[208] + c116 * source[201] + c116 * source[203]
                  - c117 * source[78] + c118 * source[71] + c118 * source[73]
                  + c119 * source[98] - c120 * source[91] - c120 * source[93]
                  - c117 * source[98] + c118 * source[91] + c118 * source[93]
                  + c119 * source[118] - c120 * source[111] - c120 * source[113];
    target[48] =  c121 * source[189] - c122 * source[184] - c122 * source[186]
                  - c123 * source[209] + c124 * source[204] + c124 * source[206]
                  - c125 * source[79] + c126 * source[74] + c126 * source[76]
                  + c127 * source[99] - c128 * source[94] - c128 * source[96]
                  - c125 * source[99] + c126 * source[94] + c126 * source[96]
                  + c127 * source[119] - c128 * source[114] - c128 * source[116];
    target[49] =  c100 * source[190] - c101 * source[192] - c99 * source[210]
                  + c100 * source[212] - c103 * source[80] + c104 * source[82]
                  + c102 * source[100] - c103 * source[102] - c103 * source[100]
                  + c104 * source[102] + c102 * source[120] - c103 * source[122];
    target[50] =  c101 * source[191] - c100 * source[193] - c100 * source[211]
                  + c99 * source[213] - c104 * source[81] + c103 * source[83]
                  + c103 * source[101] - c102 * source[103] - c104 * source[101]
                  + c103 * source[103] + c103 * source[121] - c102 * source[123];
    target[51] =  c106 * source[194] - c106 * source[196] - c105 * source[214]
                  + c105 * source[216] - c108 * source[84] + c108 * source[86]
                  + c107 * source[104] - c107 * source[106] - c108 * source[104]
                  + c108 * source[106] + c107 * source[124] - c107 * source[126];
    target[52] =  c110 * source[195] - c109 * source[215] - c112 * source[85]
                  + c111 * source[105] - c112 * source[105] + c111 * source[125];
    target[53] =  c115 * source[197] - c116 * source[190] - c116 * source[192]
                  - c113 * source[217] + c114 * source[210] + c114 * source[212]
                  - c119 * source[87] + c120 * source[80] + c120 * source[82]
                  + c117 * source[107] - c118 * source[100] - c118 * source[102]
                  - c119 * source[107] + c120 * source[100] + c120 * source[102]
                  + c117 * source[127] - c118 * source[120] - c118 * source[122];
    target[54] =  c115 * source[198] - c116 * source[191] - c116 * source[193]
                  - c113 * source[218] + c114 * source[211] + c114 * source[213]
                  - c119 * source[88] + c120 * source[81] + c120 * source[83]
                  + c117 * source[108] - c118 * source[101] - c118 * source[103]
                  - c119 * source[108] + c120 * source[101] + c120 * source[103]
                  + c117 * source[128] - c118 * source[121] - c118 * source[123];
    target[55] =  c123 * source[199] - c124 * source[194] - c124 * source[196]
                  - c121 * source[219] + c122 * source[214] + c122 * source[216]
                  - c127 * source[89] + c128 * source[84] + c128 * source[86]
                  + c125 * source[109] - c126 * source[104] - c126 * source[106]
                  - c127 * source[109] + c128 * source[104] + c128 * source[106]
                  + c125 * source[129] - c126 * source[124] - c126 * source[126];
    target[56] =  c99 * source[220] - c100 * source[222] - c99 * source[240]
                  + c100 * source[242] - c99 * source[130] + c100 * source[132]
                  + c99 * source[150] - c100 * source[152] - c99 * source[150]
                  + c100 * source[152] + c99 * source[170] - c100 * source[172]
                  + c129 * source[0] - c130 * source[2] - c129 * source[20]
                  + c130 * source[22] + c131 * source[20] - c102 * source[22]
                  - c131 * source[40] + c102 * source[42] + c129 * source[40]
                  - c130 * source[42] - c129 * source[60] + c130 * source[62];
    target[57] =  c100 * source[221] - c99 * source[223] - c100 * source[241]
                  + c99 * source[243] - c100 * source[131] + c99 * source[133]
                  + c100 * source[151] - c99 * source[153] - c100 * source[151]
                  + c99 * source[153] + c100 * source[171] - c99 * source[173]
                  + c130 * source[1] - c129 * source[3] - c130 * source[21]
                  + c129 * source[23] + c102 * source[21] - c131 * source[23]
                  - c102 * source[41] + c131 * source[43] + c130 * source[41]
                  - c129 * source[43] - c130 * source[61] + c129 * source[63];
    target[58] =  c105 * source[224] - c105 * source[226] - c105 * source[244]
                  + c105 * source[246] - c105 * source[134] + c105 * source[136]
                  + c105 * source[154] - c105 * source[156] - c105 * source[154]
                  + c105 * source[156] + c105 * source[174] - c105 * source[176]
                  + c132 * source[4] - c132 * source[6] - c132 * source[24]
                  + c132 * source[26] + c133 * source[24] - c133 * source[26]
                  - c133 * source[44] + c133 * source[46] + c132 * source[44]
                  - c132 * source[46] - c132 * source[64] + c132 * source[66];
    target[59] =  c109 * source[225] - c109 * source[245] - c109 * source[135]
                  + c109 * source[155] - c109 * source[155] + c109 * source[175]
                  + c133 * source[5] - c133 * source[25] + c134 * source[25]
                  - c134 * source[45] + c133 * source[45] - c133 * source[65];
    target[60] =  c113 * source[227] - c114 * source[220] - c114 * source[222]
                  - c113 * source[247] + c114 * source[240] + c114 * source[242]
                  - c113 * source[137] + c114 * source[130] + c114 * source[132]
                  + c113 * source[157] - c114 * source[150] - c114 * source[152]
                  - c113 * source[157] + c114 * source[150] + c114 * source[152]
                  + c113 * source[177] - c114 * source[170] - c114 * source[172]
                  + c135 * source[7] - c136 * source[0] - c136 * source[2]
                  - c135 * source[27] + c136 * source[20] + c136 * source[22]
                  + c137 * source[27] - c138 * source[20] - c138 * source[22]
                  - c137 * source[47] + c138 * source[40] + c138 * source[42]
                  + c135 * source[47] - c136 * source[40] - c136 * source[42]
                  - c135 * source[67] + c136 * source[60] + c136 * source[62];
    target[61] =  c113 * source[228] - c114 * source[221] - c114 * source[223]
                  - c113 * source[248] + c114 * source[241] + c114 * source[243]
                  - c113 * source[138] + c114 * source[131] + c114 * source[133]
                  + c113 * source[158] - c114 * source[151] - c114 * source[153]
                  - c113 * source[158] + c114 * source[151] + c114 * source[153]
                  + c113 * source[178] - c114 * source[171] - c114 * source[173]
                  + c135 * source[8] - c136 * source[1] - c136 * source[3]
                  - c135 * source[28] + c136 * source[21] + c136 * source[23]
                  + c137 * source[28] - c138 * source[21] - c138 * source[23]
                  - c137 * source[48] + c138 * source[41] + c138 * source[43]
                  + c135 * source[48] - c136 * source[41] - c136 * source[43]
                  - c135 * source[68] + c136 * source[61] + c136 * source[63];
    target[62] =  c121 * source[229] - c122 * source[224] - c122 * source[226]
                  - c121 * source[249] + c122 * source[244] + c122 * source[246]
                  - c121 * source[139] + c122 * source[134] + c122 * source[136]
                  + c121 * source[159] - c122 * source[154] - c122 * source[156]
                  - c121 * source[159] + c122 * source[154] + c122 * source[156]
                  + c121 * source[179] - c122 * source[174] - c122 * source[176]
                  + c139 * source[9] - c140 * source[4] - c140 * source[6]
                  - c139 * source[29] + c140 * source[24] + c140 * source[26]
                  + c141 * source[29] - c142 * source[24] - c142 * source[26]
                  - c141 * source[49] + c142 * source[44] + c142 * source[46]
                  + c139 * source[49] - c140 * source[44] - c140 * source[46]
                  - c139 * source[69] + c140 * source[64] + c140 * source[66];
    target[63] =  c143 * source[230] - c144 * source[232] - c143 * source[140]
                  + c144 * source[142] - c143 * source[160] + c144 * source[162]
                  + c131 * source[10] - c102 * source[12] + c145 * source[30]
                  - c146 * source[32] + c131 * source[50] - c102 * source[52];
    target[64] =  c144 * source[231] - c143 * source[233] - c144 * source[141]
                  + c143 * source[143] - c144 * source[161] + c143 * source[163]
                  + c102 * source[11] - c131 * source[13] + c146 * source[31]
                  - c145 * source[33] + c102 * source[51] - c131 * source[53];
    target[65] =  c109 * source[234] - c109 * source[236] - c109 * source[144]
                  + c109 * source[146] - c109 * source[164] + c109 * source[166]
                  + c133 * source[14] - c133 * source[16] + c134 * source[34]
                  - c134 * source[36] + c133 * source[54] - c133 * source[56];
    target[66] =  c147 * source[235] - c147 * source[145] - c147 * source[165]
                  + c134 * source[15] + c148 * source[35] + c134 * source[55];
    target[67] =  c149 * source[237] - c150 * source[230] - c150 * source[232]
                  - c149 * source[147] + c150 * source[140] + c150 * source[142]
                  - c149 * source[167] + c150 * source[160] + c150 * source[162]
                  + c137 * source[17] - c138 * source[10] - c138 * source[12]
                  + c114 * source[37] - c135 * source[30] - c135 * source[32]
                  + c137 * source[57] - c138 * source[50] - c138 * source[52];
    target[68] =  c149 * source[238] - c150 * source[231] - c150 * source[233]
                  - c149 * source[148] + c150 * source[141] + c150 * source[143]
                  - c149 * source[168] + c150 * source[161] + c150 * source[163]
                  + c137 * source[18] - c138 * source[11] - c138 * source[13]
                  + c114 * source[38] - c135 * source[31] - c135 * source[33]
                  + c137 * source[58] - c138 * source[51] - c138 * source[53];
    target[69] =  c151 * source[239] - c123 * source[234] - c123 * source[236]
                  - c151 * source[149] + c123 * source[144] + c123 * source[146]
                  - c151 * source[169] + c123 * source[164] + c123 * source[166]
                  + c141 * source[19] - c142 * source[14] - c142 * source[16]
                  + c152 * source[39] - c125 * source[34] - c125 * source[36]
                  + c141 * source[59] - c142 * source[54] - c142 * source[56];
    target[70] =  c153 * source[250] - c122 * source[252] - c154 * source[180]
                  + c155 * source[182] - c154 * source[200] + c155 * source[202]
                  + c156 * source[70] - c157 * source[72] + c158 * source[90]
                  - c159 * source[92] + c156 * source[110] - c157 * source[112];
    target[71] =  c122 * source[251] - c153 * source[253] - c155 * source[181]
                  + c154 * source[183] - c155 * source[201] + c154 * source[203]
                  + c157 * source[71] - c156 * source[73] + c159 * source[91]
                  - c158 * source[93] + c157 * source[111] - c156 * source[113];
    target[72] =  c150 * source[254] - c150 * source[256] - c160 * source[184]
                  + c160 * source[186] - c160 * source[204] + c160 * source[206]
                  + c161 * source[74] - c161 * source[76] + c162 * source[94]
                  - c162 * source[96] + c161 * source[114] - c161 * source[116];
    target[73] =  c113 * source[255] - c163 * source[185] - c163 * source[205]
                  + c162 * source[75] + c160 * source[95] + c162 * source[115];
    target[74] =  c164 * source[257] - c165 * source[250] - c165 * source[252]
                  - c109 * source[187] + c148 * source[180] + c148 * source[182]
                  - c109 * source[207] + c148 * source[200] + c148 * source[202]
                  + c148 * source[77] - c133 * source[70] - c133 * source[72]
                  + c105 * source[97] - c134 * source[90] - c134 * source[92]
                  + c148 * source[117] - c133 * source[110] - c133 * source[112];
    target[75] =  c164 * source[258] - c165 * source[251] - c165 * source[253]
                  - c109 * source[188] + c148 * source[181] + c148 * source[183]
                  - c109 * source[208] + c148 * source[201] + c148 * source[203]
                  + c148 * source[78] - c133 * source[71] - c133 * source[73]
                  + c105 * source[98] - c134 * source[91] - c134 * source[93]
                  + c148 * source[118] - c133 * source[111] - c133 * source[113];
    target[76] =  c166 * source[259] - c167 * source[254] - c167 * source[256]
                  - c143 * source[189] + c100 * source[184] + c100 * source[186]
                  - c143 * source[209] + c100 * source[204] + c100 * source[206]
                  + c168 * source[79] - c146 * source[74] - c146 * source[76]
                  + c99 * source[99] - c169 * source[94] - c169 * source[96]
                  + c168 * source[119] - c146 * source[114] - c146 * source[116];
    target[77] =  c153 * source[260] - c122 * source[262] - c154 * source[190]
                  + c155 * source[192] - c154 * source[210] + c155 * source[212]
                  + c156 * source[80] - c157 * source[82] + c158 * source[100]
                  - c159 * source[102] + c156 * source[120] - c157 * source[122];
    target[78] =  c122 * source[261] - c153 * source[263] - c155 * source[191]
                  + c154 * source[193] - c155 * source[211] + c154 * source[213]
                  + c157 * source[81] - c156 * source[83] + c159 * source[101]
                  - c158 * source[103] + c157 * source[121] - c156 * source[123];
    target[79] =  c150 * source[264] - c150 * source[266] - c160 * source[194]
                  + c160 * source[196] - c160 * source[214] + c160 * source[216]
                  + c161 * source[84] - c161 * source[86] + c162 * source[104]
                  - c162 * source[106] + c161 * source[124] - c161 * source[126];
    target[80] =  c113 * source[265] - c163 * source[195] - c163 * source[215]
                  + c162 * source[85] + c160 * source[105] + c162 * source[125];
    target[81] =  c164 * source[267] - c165 * source[260] - c165 * source[262]
                  - c109 * source[197] + c148 * source[190] + c148 * source[192]
                  - c109 * source[217] + c148 * source[210] + c148 * source[212]
                  + c148 * source[87] - c133 * source[80] - c133 * source[82]
                  + c105 * source[107] - c134 * source[100] - c134 * source[102]
                  + c148 * source[127] - c133 * source[120] - c133 * source[122];
    target[82] =  c164 * source[268] - c165 * source[261] - c165 * source[263]
                  - c109 * source[198] + c148 * source[191] + c148 * source[193]
                  - c109 * source[218] + c148 * source[211] + c148 * source[213]
                  + c148 * source[88] - c133 * source[81] - c133 * source[83]
                  + c105 * source[108] - c134 * source[101] - c134 * source[103]
                  + c148 * source[128] - c133 * source[121] - c133 * source[123];
    target[83] =  c166 * source[269] - c167 * source[264] - c167 * source[266]
                  - c143 * source[199] + c100 * source[194] + c100 * source[196]
                  - c143 * source[219] + c100 * source[214] + c100 * source[216]
                  + c168 * source[89] - c146 * source[84] - c146 * source[86]
                  + c99 * source[109] - c169 * source[104] - c169 * source[106]
                  + c168 * source[129] - c146 * source[124] - c146 * source[126];
    target[84] =  c170 * source[270] - c171 * source[272] - c172 * source[220]
                  + c173 * source[222] - c172 * source[240] + c173 * source[242]
                  + c174 * source[130] - c175 * source[132] + c176 * source[150]
                  - c177 * source[152] + c174 * source[170] - c175 * source[172]
                  - c178 * source[0] + c179 * source[2] - c179 * source[20]
                  + c180 * source[22] - c179 * source[40] + c180 * source[42]
                  - c178 * source[60] + c179 * source[62];
    target[85] =  c171 * source[271] - c170 * source[273] - c173 * source[221]
                  + c172 * source[223] - c173 * source[241] + c172 * source[243]
                  + c175 * source[131] - c174 * source[133] + c177 * source[151]
                  - c176 * source[153] + c175 * source[171] - c174 * source[173]
                  - c179 * source[1] + c178 * source[3] - c180 * source[21]
                  + c179 * source[23] - c180 * source[41] + c179 * source[43]
                  - c179 * source[61] + c178 * source[63];
    target[86] =  c181 * source[274] - c181 * source[276] - c182 * source[224]
                  + c182 * source[226] - c182 * source[244] + c182 * source[246]
                  + c183 * source[134] - c183 * source[136] + c184 * source[154]
                  - c184 * source[156] + c183 * source[174] - c183 * source[176]
                  - c185 * source[4] + c185 * source[6] - c186 * source[24]
                  + c186 * source[26] - c186 * source[44] + c186 * source[46]
                  - c185 * source[64] + c185 * source[66];
    target[87] =  c187 * source[275] - c188 * source[225] - c188 * source[245]
                  + c184 * source[135] + c189 * source[155] + c184 * source[175]
                  - c190 * source[5] - c191 * source[25] - c191 * source[45]
                  - c190 * source[65];
    target[88] =  c192 * source[277] - c193 * source[270] - c193 * source[272]
                  - c194 * source[227] + c195 * source[220] + c195 * source[222]
                  - c194 * source[247] + c195 * source[240] + c195 * source[242]
                  + c196 * source[137] - c197 * source[130] - c197 * source[132]
                  + c198 * source[157] - c199 * source[150] - c199 * source[152]
                  + c196 * source[177] - c197 * source[170] - c197 * source[172]
                  - c200 * source[7] + c201 * source[0] + c201 * source[2]
                  - c202 * source[27] + c203 * source[20] + c203 * source[22]
                  - c202 * source[47] + c203 * source[40] + c203 * source[42]
                  - c200 * source[67] + c201 * source[60] + c201 * source[62];
    target[89] =  c192 * source[278] - c193 * source[271] - c193 * source[273]
                  - c194 * source[228] + c195 * source[221] + c195 * source[223]
                  - c194 * source[248] + c195 * source[241] + c195 * source[243]
                  + c196 * source[138] - c197 * source[131] - c197 * source[133]
                  + c198 * source[158] - c199 * source[151] - c199 * source[153]
                  + c196 * source[178] - c197 * source[171] - c197 * source[173]
                  - c200 * source[8] + c201 * source[1] + c201 * source[3]
                  - c202 * source[28] + c203 * source[21] + c203 * source[23]
                  - c202 * source[48] + c203 * source[41] + c203 * source[43]
                  - c200 * source[68] + c201 * source[61] + c201 * source[63];
    target[90] =  source[279] - c204 * source[274] - c204 * source[276]
                  - c205 * source[229] + c206 * source[224] + c206 * source[226]
                  - c205 * source[249] + c206 * source[244] + c206 * source[246]
                  + c207 * source[139] - c208 * source[134] - c208 * source[136]
                  + c206 * source[159] - c209 * source[154] - c209 * source[156]
                  + c207 * source[179] - c208 * source[174] - c208 * source[176]
                  - c210 * source[9] + c211 * source[4] + c211 * source[6]
                  - c212 * source[29] + c213 * source[24] + c213 * source[26]
                  - c212 * source[49] + c213 * source[44] + c213 * source[46]
                  - c210 * source[69] + c211 * source[64] + c211 * source[66];
  }
}

void CCarSphList::carsph_63(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c68 = 115.27819611704548;
  const double c42 = 90.117111305234374;
  const double c110 = 84.187291202413675;
  const double c93 = 76.852130744696993;
  const double c73 = 72.908332857088425;
  const double c58 = 70.593189738812626;
  const double c64 = 57.639098058522741;
  const double c45 = 56.995065575889988;
  const double c147 = 56.124860801609124;
  const double c35 = 55.185234947855392;
  const double c115 = 53.244718047896548;
  const double c23 = 52.029138470668528;
  const double c101 = 51.553976568253198;
  const double c95 = 48.605555238058955;
  const double c88 = 47.062126492541751;
  const double c39 = 45.058555652617187;
  const double c82 = 44.647053374214963;
  const double c163 = 44.37059837324712;
  const double c189 = 43.571062644833439;
  const double c106 = 42.093645601206838;
  const double c7 = 39.021853853001396;
  const double c91 = 38.426065372348496;
  const double c149 = 35.496478698597699;
  const double c52 = 34.902207129349286;
  const double c144 = 34.369317712168801;
  const double c26 = 32.906116452720461;
  const double c124 = 32.605597678926237;
  const double c19 = 31.861210252437054;
  const double c112 = 31.57023420090513;
  const double c81 = 29.764702249476645;
  const double c188 = 29.047375096555626;
  const double c47 = 28.497532787944994;
  const double c109 = 28.062430400804562;
  const double c37 = 27.592617473927696;
  const double c198 = 27.556759606310752;
  const double c155 = 27.171331399105199;
  const double c177 = 26.6817177576707;
  const double c21 = 26.014569235334264;
  const double c10 = 24.679587339540344;
  const double c3 = 23.895907689327789;
  const double c57 = 23.531063246270875;
  const double c51 = 23.268138086232856;
  const double c40 = 22.529277826308594;
  const double c160 = 22.18529918662356;
  const double c184 = 21.78553132241672;
  const double c123 = 21.737065119284157;
  const double c31 = 20.150798681441884;
  const double c119 = 19.966769267961205;
  const double c97 = 19.843134832984429;
  const double c5 = 19.510926926500698;
  const double c104 = 19.332741213094948;
  const double c67 = 19.213032686174248;
  const double c34 = 18.395078315951796;
  const double c194 = 18.371173070873837;
  const double c74 = 18.227083214272106;
  const double c173 = 17.787811838447134;
  const double c113 = 17.748239349298849;
  const double c54 = 17.451103564674643;
  const double c100 = 17.1846588560844;
  const double c209 = 16.875;
  const double c108 = 15.785117100452565;
  const double c87 = 15.687375497513916;
  const double c22 = 15.608741541200558;
  const double c15 = 15.113099011081413;
  const double c182 = 14.523687548277813;
  const double c151 = 14.491376746189438;
  const double c46 = 14.248766393972497;
  const double c105 = 14.031215200402281;
  const double c196 = 13.778379803155376;
  const double c159 = 13.5856656995526;
  const double c30 = 13.433865787627923;
  const double c175 = 13.34085887883535;
  const double c116 = 13.311179511974137;
  const double c128 = 12.227099129597338;
  const double c71 = 12.151388809514739;
  const double c56 = 11.765531623135438;
  const double c53 = 11.634069043116428;
  const double c70 = 11.527819611704547;
  const double c143 = 11.456439237389599;
  const double c206 = 11.25;
  const double c164 = 11.224972160321824;
  const double c162 = 11.09264959331178;
  const double c183 = 10.89276566120836;
  const double c122 = 10.868532559642079;
  const double c18 = 10.620403417479018;
  const double c111 = 10.523411400301709;
  const double c14 = 10.075399340720942;
  const double c24 = 9.8718349358161372;
  const double c63 = 9.6065163430871241;
  const double c17 = 9.5583630757311155;
  const double c36 = 9.197539157975898;
  const double c154 = 9.0571104663683997;
  const double c41 = 9.0117111305234374;
  const double c176 = 8.8939059192235668;
  const double c150 = 8.8741196746494246;
  const double c169 = 8.5923294280422002;
  const double c208 = 8.4375;
  const double c27 = 8.2265291131801153;
  const double c127 = 8.1513994197315593;
  const double c2 = 7.9653025631092635;
  const double c20 = 7.8043707706002792;
  const double c94 = 7.6852130744696989;
  const double c205 = 7.5;
  const double c80 = 7.4411755623691613;
  const double c77 = 7.2908332857088425;
  const double c121 = 7.245688373094719;
  const double c48 = 7.1243831969862486;
  const double c62 = 7.0593189738812621;
  const double c148 = 7.0156076002011405;
  const double c199 = 6.8891899015776881;
  const double c167 = 6.8738635424337602;
  const double c157 = 6.7928328497762998;
  const double c117 = 6.6555897559870685;
  const double c103 = 6.4442470710316497;
  const double c11 = 6.169896834885086;
  const double c29 = 6.0452396044325658;
  const double c172 = 5.9292706128157109;
  const double c66 = 5.7639098058522737;
  const double c99 = 5.7282196186947996;
  const double c43 = 5.6995065575889985;
  const double c207 = 5.625;
  const double c161 = 5.54632479665589;
  const double c33 = 5.5185234947855388;
  const double c107 = 5.2617057001508547;
  const double c120 = 4.9916923169903011;
  const double c79 = 4.9607837082461073;
  const double c96 = 4.8605555238058953;
  const double c90 = 4.7062126492541747;
  const double c195 = 4.5927932677184593;
  const double c166 = 4.5825756949558398;
  const double c158 = 4.5285552331841998;
  const double c38 = 4.5058555652617187;
  const double c86 = 4.4647053374214964;
  const double c174 = 4.4469529596117834;
  const double c114 = 4.4370598373247123;
  const double c146 = 4.2961647140211001;
  const double c126 = 4.0756997098657797;
  const double c28 = 4.0301597362883772;
  const double c55 = 3.9218438743784789;
  const double c187 = 3.872983346207417;
  const double c92 = 3.8426065372348495;
  const double c191 = 3.6309218870694533;
  const double c153 = 3.6228441865473595;
  const double c134 = 3.5078038001005702;
  const double c50 = 3.4902207129349283;
  const double c197 = 3.444594950788844;
  const double c16 = 3.1861210252437053;
  const double c72 = 3.0378472023786847;
  const double c85 = 2.9764702249476644;
  const double c168 = 2.8641098093473998;
  const double c165 = 2.8062430400804561;
  const double c125 = 2.7171331399105196;
  const double c6 = 2.6014569235334264;
  const double c25 = 2.4679587339540343;
  const double c192 = 2.4494897427831779;
  const double c171 = 2.3717082451262845;
  const double c61 = 2.3531063246270874;
  const double c49 = 2.3268138086232857;
  const double c202 = 2.2963966338592297;
  const double c156 = 2.2642776165920999;
  const double c180 = 2.2234764798058917;
  const double c137 = 2.2185299186623562;
  const double c102 = 2.1480823570105501;
  const double c98 = 1.984313483298443;
  const double c181 = 1.9364916731037085;
  const double c69 = 1.9213032686174247;
  const double c32 = 1.8395078315951796;
  const double c78 = 1.8227083214272106;
  const double c186 = 1.8154609435347266;
  const double c152 = 1.8114220932736798;
  const double c133 = 1.7539019000502851;
  const double c118 = 1.6638974389967671;
  const double c8 = 1.6453058226360229;
  const double c1 = 1.5930605126218527;
  const double c89 = 1.5687375497513916;
  const double c204 = 1.5;
  const double c145 = 1.4320549046736999;
  const double c44 = 1.4248766393972496;
  const double c213 = 1.40625;
  const double c142 = 1.3585665699552598;
  const double c4 = 1.3007284617667132;
  const double c75 = 1.2151388809514738;
  const double c190 = 1.2103072956898178;
  const double c60 = 1.1765531623135437;
  const double c135 = 1.1092649593311781;
  const double c130 = 1.074041178505275;
  const double c13 = 1.0075399340720943;
  const double c65 = 0.96065163430871237;
  const double c212 = 0.9375;
  const double c141 = 0.90571104663683988;
  const double c132 = 0.87695095002514256;
  const double c170 = 0.79056941504209488;
  const double c200 = 0.76546554461974314;
  const double c84 = 0.74411755623691611;
  const double c179 = 0.74115882660196386;
  const double c131 = 0.71602745233684995;
  const double c140 = 0.67928328497762991;
  const double c12 = 0.67169328938139616;
  const double c193 = 0.61237243569579447;
  const double c185 = 0.60515364784490888;
  const double c203 = 0.57409915846480741;
  const double c138 = 0.55463247966558904;
  const double c0 = 0.53102017087395093;
  const double c83 = 0.49607837082461076;
  const double c211 = 0.46875;
  const double c139 = 0.45285552331841994;
  const double c9 = 0.41132645565900572;
  const double c59 = 0.39218438743784789;
  const double c129 = 0.35801372616842497;
  const double c210 = 0.3125;
  const double c76 = 0.30378472023786846;
  const double c136 = 0.27731623983279452;
  const double c178 = 0.24705294220065463;
  const double c201 = 0.19136638615493579;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 91, source += 280) {
    target[0] =  c0 * source[0] - c1 * source[2] - c2 * source[20]
                  + c3 * source[22] + c2 * source[40] - c3 * source[42]
                  - c0 * source[60] + c1 * source[62];
    target[1] =  c1 * source[1] - c0 * source[3] - c3 * source[21]
                  + c2 * source[23] + c3 * source[41] - c2 * source[43]
                  - c1 * source[61] + c0 * source[63];
    target[2] =  c4 * source[4] - c4 * source[6] - c5 * source[24]
                  + c5 * source[26] + c5 * source[44] - c5 * source[46]
                  - c4 * source[64] + c4 * source[66];
    target[3] =  c6 * source[5] - c7 * source[25] + c7 * source[45]
                  - c6 * source[65];
    target[4] =  c8 * source[7] - c9 * source[0] - c9 * source[2]
                  - c10 * source[27] + c11 * source[20] + c11 * source[22]
                  + c10 * source[47] - c11 * source[40] - c11 * source[42]
                  - c8 * source[67] + c9 * source[60] + c9 * source[62];
    target[5] =  c8 * source[8] - c9 * source[1] - c9 * source[3]
                  - c10 * source[28] + c11 * source[21] + c11 * source[23]
                  + c10 * source[48] - c11 * source[41] - c11 * source[43]
                  - c8 * source[68] + c9 * source[61] + c9 * source[63];
    target[6] =  c12 * source[9] - c13 * source[4] - c13 * source[6]
                  - c14 * source[29] + c15 * source[24] + c15 * source[26]
                  + c14 * source[49] - c15 * source[44] - c15 * source[46]
                  - c12 * source[69] + c13 * source[64] + c13 * source[66];
    target[7] =  c16 * source[10] - c17 * source[12] - c18 * source[30]
                  + c19 * source[32] + c16 * source[50] - c17 * source[52];
    target[8] =  c17 * source[11] - c16 * source[13] - c19 * source[31]
                  + c18 * source[33] + c17 * source[51] - c16 * source[53];
    target[9] =  c20 * source[14] - c20 * source[16] - c21 * source[34]
                  + c21 * source[36] + c20 * source[54] - c20 * source[56];
    target[10] =  c22 * source[15] - c23 * source[35] + c22 * source[55];
    target[11] =  c24 * source[17] - c25 * source[10] - c25 * source[12]
                  - c26 * source[37] + c27 * source[30] + c27 * source[32]
                  + c24 * source[57] - c25 * source[50] - c25 * source[52];
    target[12] =  c24 * source[18] - c25 * source[11] - c25 * source[13]
                  - c26 * source[38] + c27 * source[31] + c27 * source[33]
                  + c24 * source[58] - c25 * source[51] - c25 * source[53];
    target[13] =  c28 * source[19] - c29 * source[14] - c29 * source[16]
                  - c30 * source[39] + c31 * source[34] + c31 * source[36]
                  + c28 * source[59] - c29 * source[54] - c29 * source[56];
    target[14] =  c32 * source[70] - c33 * source[72] - c34 * source[90]
                  + c35 * source[92] + c36 * source[110] - c37 * source[112];
    target[15] =  c33 * source[71] - c32 * source[73] - c35 * source[91]
                  + c34 * source[93] + c37 * source[111] - c36 * source[113];
    target[16] =  c38 * source[74] - c38 * source[76] - c39 * source[94]
                  + c39 * source[96] + c40 * source[114] - c40 * source[116];
    target[17] =  c41 * source[75] - c42 * source[95] + c39 * source[115];
    target[18] =  c43 * source[77] - c44 * source[70] - c44 * source[72]
                  - c45 * source[97] + c46 * source[90] + c46 * source[92]
                  + c47 * source[117] - c48 * source[110] - c48 * source[112];
    target[19] =  c43 * source[78] - c44 * source[71] - c44 * source[73]
                  - c45 * source[98] + c46 * source[91] + c46 * source[93]
                  + c47 * source[118] - c48 * source[111] - c48 * source[113];
    target[20] =  c49 * source[79] - c50 * source[74] - c50 * source[76]
                  - c51 * source[99] + c52 * source[94] + c52 * source[96]
                  + c53 * source[119] - c54 * source[114] - c54 * source[116];
    target[21] =  c36 * source[80] - c37 * source[82] - c34 * source[100]
                  + c35 * source[102] + c32 * source[120] - c33 * source[122];
    target[22] =  c37 * source[81] - c36 * source[83] - c35 * source[101]
                  + c34 * source[103] + c33 * source[121] - c32 * source[123];
    target[23] =  c40 * source[84] - c40 * source[86] - c39 * source[104]
                  + c39 * source[106] + c38 * source[124] - c38 * source[126];
    target[24] =  c39 * source[85] - c42 * source[105] + c41 * source[125];
    target[25] =  c47 * source[87] - c48 * source[80] - c48 * source[82]
                  - c45 * source[107] + c46 * source[100] + c46 * source[102]
                  + c43 * source[127] - c44 * source[120] - c44 * source[122];
    target[26] =  c47 * source[88] - c48 * source[81] - c48 * source[83]
                  - c45 * source[108] + c46 * source[101] + c46 * source[103]
                  + c43 * source[128] - c44 * source[121] - c44 * source[123];
    target[27] =  c53 * source[89] - c54 * source[84] - c54 * source[86]
                  - c51 * source[109] + c52 * source[104] + c52 * source[106]
                  + c49 * source[129] - c50 * source[124] - c50 * source[126];
    target[28] =  c55 * source[130] - c56 * source[132] - c57 * source[150]
                  + c58 * source[152] + c55 * source[170] - c56 * source[172]
                  - c59 * source[0] + c60 * source[2] + c61 * source[20]
                  - c62 * source[22] - c59 * source[40] + c60 * source[42]
                  - c59 * source[20] + c60 * source[22] + c61 * source[40]
                  - c62 * source[42] - c59 * source[60] + c60 * source[62];
    target[29] =  c56 * source[131] - c55 * source[133] - c58 * source[151]
                  + c57 * source[153] + c56 * source[171] - c55 * source[173]
                  - c60 * source[1] + c59 * source[3] + c62 * source[21]
                  - c61 * source[23] - c60 * source[41] + c59 * source[43]
                  - c60 * source[21] + c59 * source[23] + c62 * source[41]
                  - c61 * source[43] - c60 * source[61] + c59 * source[63];
    target[30] =  c63 * source[134] - c63 * source[136] - c64 * source[154]
                  + c64 * source[156] + c63 * source[174] - c63 * source[176]
                  - c65 * source[4] + c65 * source[6] + c66 * source[24]
                  - c66 * source[26] - c65 * source[44] + c65 * source[46]
                  - c65 * source[24] + c65 * source[26] + c66 * source[44]
                  - c66 * source[46] - c65 * source[64] + c65 * source[66];
    target[31] =  c67 * source[135] - c68 * source[155] + c67 * source[175]
                  - c69 * source[5] + c70 * source[25] - c69 * source[45]
                  - c69 * source[25] + c70 * source[45] - c69 * source[65];
    target[32] =  c71 * source[137] - c72 * source[130] - c72 * source[132]
                  - c73 * source[157] + c74 * source[150] + c74 * source[152]
                  + c71 * source[177] - c72 * source[170] - c72 * source[172]
                  - c75 * source[7] + c76 * source[0] + c76 * source[2]
                  + c77 * source[27] - c78 * source[20] - c78 * source[22]
                  - c75 * source[47] + c76 * source[40] + c76 * source[42]
                  - c75 * source[27] + c76 * source[20] + c76 * source[22]
                  + c77 * source[47] - c78 * source[40] - c78 * source[42]
                  - c75 * source[67] + c76 * source[60] + c76 * source[62];
    target[33] =  c71 * source[138] - c72 * source[131] - c72 * source[133]
                  - c73 * source[158] + c74 * source[151] + c74 * source[153]
                  + c71 * source[178] - c72 * source[171] - c72 * source[173]
                  - c75 * source[8] + c76 * source[1] + c76 * source[3]
                  + c77 * source[28] - c78 * source[21] - c78 * source[23]
                  - c75 * source[48] + c76 * source[41] + c76 * source[43]
                  - c75 * source[28] + c76 * source[21] + c76 * source[23]
                  + c77 * source[48] - c78 * source[41] - c78 * source[43]
                  - c75 * source[68] + c76 * source[61] + c76 * source[63];
    target[34] =  c79 * source[139] - c80 * source[134] - c80 * source[136]
                  - c81 * source[159] + c82 * source[154] + c82 * source[156]
                  + c79 * source[179] - c80 * source[174] - c80 * source[176]
                  - c83 * source[9] + c84 * source[4] + c84 * source[6]
                  + c85 * source[29] - c86 * source[24] - c86 * source[26]
                  - c83 * source[49] + c84 * source[44] + c84 * source[46]
                  - c83 * source[29] + c84 * source[24] + c84 * source[26]
                  + c85 * source[49] - c86 * source[44] - c86 * source[46]
                  - c83 * source[69] + c84 * source[64] + c84 * source[66];
    target[35] =  c87 * source[140] - c88 * source[142] - c87 * source[160]
                  + c88 * source[162] - c89 * source[10] + c90 * source[12]
                  + c89 * source[30] - c90 * source[32] - c89 * source[30]
                  + c90 * source[32] + c89 * source[50] - c90 * source[52];
    target[36] =  c88 * source[141] - c87 * source[143] - c88 * source[161]
                  + c87 * source[163] - c90 * source[11] + c89 * source[13]
                  + c90 * source[31] - c89 * source[33] - c90 * source[31]
                  + c89 * source[33] + c90 * source[51] - c89 * source[53];
    target[37] =  c91 * source[144] - c91 * source[146] - c91 * source[164]
                  + c91 * source[166] - c92 * source[14] + c92 * source[16]
                  + c92 * source[34] - c92 * source[36] - c92 * source[34]
                  + c92 * source[36] + c92 * source[54] - c92 * source[56];
    target[38] =  c93 * source[145] - c93 * source[165] - c94 * source[15]
                  + c94 * source[35] - c94 * source[35] + c94 * source[55];
    target[39] =  c95 * source[147] - c71 * source[140] - c71 * source[142]
                  - c95 * source[167] + c71 * source[160] + c71 * source[162]
                  - c96 * source[17] + c75 * source[10] + c75 * source[12]
                  + c96 * source[37] - c75 * source[30] - c75 * source[32]
                  - c96 * source[37] + c75 * source[30] + c75 * source[32]
                  + c96 * source[57] - c75 * source[50] - c75 * source[52];
    target[40] =  c95 * source[148] - c71 * source[141] - c71 * source[143]
                  - c95 * source[168] + c71 * source[161] + c71 * source[163]
                  - c96 * source[18] + c75 * source[11] + c75 * source[13]
                  + c96 * source[38] - c75 * source[31] - c75 * source[33]
                  - c96 * source[38] + c75 * source[31] + c75 * source[33]
                  + c96 * source[58] - c75 * source[51] - c75 * source[53];
    target[41] =  c97 * source[149] - c81 * source[144] - c81 * source[146]
                  - c97 * source[169] + c81 * source[164] + c81 * source[166]
                  - c98 * source[19] + c85 * source[14] + c85 * source[16]
                  + c98 * source[39] - c85 * source[34] - c85 * source[36]
                  - c98 * source[39] + c85 * source[34] + c85 * source[36]
                  + c98 * source[59] - c85 * source[54] - c85 * source[56];
    target[42] =  c99 * source[180] - c100 * source[182] - c100 * source[200]
                  + c101 * source[202] - c102 * source[70] + c103 * source[72]
                  + c103 * source[90] - c104 * source[92] - c102 * source[90]
                  + c103 * source[92] + c103 * source[110] - c104 * source[112];
    target[43] =  c100 * source[181] - c99 * source[183] - c101 * source[201]
                  + c100 * source[203] - c103 * source[71] + c102 * source[73]
                  + c104 * source[91] - c103 * source[93] - c103 * source[91]
                  + c102 * source[93] + c104 * source[111] - c103 * source[113];
    target[44] =  c105 * source[184] - c105 * source[186] - c106 * source[204]
                  + c106 * source[206] - c107 * source[74] + c107 * source[76]
                  + c108 * source[94] - c108 * source[96] - c107 * source[94]
                  + c107 * source[96] + c108 * source[114] - c108 * source[116];
    target[45] =  c109 * source[185] - c110 * source[205] - c111 * source[75]
                  + c112 * source[95] - c111 * source[95] + c112 * source[115];
    target[46] =  c113 * source[187] - c114 * source[180] - c114 * source[182]
                  - c115 * source[207] + c116 * source[200] + c116 * source[202]
                  - c117 * source[77] + c118 * source[70] + c118 * source[72]
                  + c119 * source[97] - c120 * source[90] - c120 * source[92]
                  - c117 * source[97] + c118 * source[90] + c118 * source[92]
                  + c119 * source[117] - c120 * source[110] - c120 * source[112];
    target[47] =  c113 * source[188] - c114 * source[181] - c114 * source[183]
                  - c115 * source[208] + c116 * source[201] + c116 * source[203]
                  - c117 * source[78] + c118 * source[71] + c118 * source[73]
                  + c119 * source[98] - c120 * source[91] - c120 * source[93]
                  - c117 * source[98] + c118 * source[91] + c118 * source[93]
                  + c119 * source[118] - c120 * source[111] - c120 * source[113];
    target[48] =  c121 * source[189] - c122 * source[184] - c122 * source[186]
                  - c123 * source[209] + c124 * source[204] + c124 * source[206]
                  - c125 * source[79] + c126 * source[74] + c126 * source[76]
                  + c127 * source[99] - c128 * source[94] - c128 * source[96]
                  - c125 * source[99] + c126 * source[94] + c126 * source[96]
                  + c127 * source[119] - c128 * source[114] - c128 * source[116];
    target[49] =  c100 * source[190] - c101 * source[192] - c99 * source[210]
                  + c100 * source[212] - c103 * source[80] + c104 * source[82]
                  + c102 * source[100] - c103 * source[102] - c103 * source[100]
                  + c104 * source[102] + c102 * source[120] - c103 * source[122];
    target[50] =  c101 * source[191] - c100 * source[193] - c100 * source[211]
                  + c99 * source[213] - c104 * source[81] + c103 * source[83]
                  + c103 * source[101] - c102 * source[103] - c104 * source[101]
                  + c103 * source[103] + c103 * source[121] - c102 * source[123];
    target[51] =  c106 * source[194] - c106 * source[196] - c105 * source[214]
                  + c105 * source[216] - c108 * source[84] + c108 * source[86]
                  + c107 * source[104] - c107 * source[106] - c108 * source[104]
                  + c108 * source[106] + c107 * source[124] - c107 * source[126];
    target[52] =  c110 * source[195] - c109 * source[215] - c112 * source[85]
                  + c111 * source[105] - c112 * source[105] + c111 * source[125];
    target[53] =  c115 * source[197] - c116 * source[190] - c116 * source[192]
                  - c113 * source[217] + c114 * source[210] + c114 * source[212]
                  - c119 * source[87] + c120 * source[80] + c120 * source[82]
                  + c117 * source[107] - c118 * source[100] - c118 * source[102]
                  - c119 * source[107] + c120 * source[100] + c120 * source[102]
                  + c117 * source[127] - c118 * source[120] - c118 * source[122];
    target[54] =  c115 * source[198] - c116 * source[191] - c116 * source[193]
                  - c113 * source[218] + c114 * source[211] + c114 * source[213]
                  - c119 * source[88] + c120 * source[81] + c120 * source[83]
                  + c117 * source[108] - c118 * source[101] - c118 * source[103]
                  - c119 * source[108] + c120 * source[101] + c120 * source[103]
                  + c117 * source[128] - c118 * source[121] - c118 * source[123];
    target[55] =  c123 * source[199] - c124 * source[194] - c124 * source[196]
                  - c121 * source[219] + c122 * source[214] + c122 * source[216]
                  - c127 * source[89] + c128 * source[84] + c128 * source[86]
                  + c125 * source[109] - c126 * source[104] - c126 * source[106]
                  - c127 * source[109] + c128 * source[104] + c128 * source[106]
                  + c125 * source[129] - c126 * source[124] - c126 * source[126];
    target[56] =  c99 * source[220] - c100 * source[222] - c99 * source[240]
                  + c100 * source[242] - c99 * source[130] + c100 * source[132]
                  + c99 * source[150] - c100 * source[152] - c99 * source[150]
                  + c100 * source[152] + c99 * source[170] - c100 * source[172]
                  + c129 * source[0] - c130 * source[2] - c129 * source[20]
                  + c130 * source[22] + c131 * source[20] - c102 * source[22]
                  - c131 * source[40] + c102 * source[42] + c129 * source[40]
                  - c130 * source[42] - c129 * source[60] + c130 * source[62];
    target[57] =  c100 * source[221] - c99 * source[223] - c100 * source[241]
                  + c99 * source[243] - c100 * source[131] + c99 * source[133]
                  + c100 * source[151] - c99 * source[153] - c100 * source[151]
                  + c99 * source[153] + c100 * source[171] - c99 * source[173]
                  + c130 * source[1] - c129 * source[3] - c130 * source[21]
                  + c129 * source[23] + c102 * source[21] - c131 * source[23]
                  - c102 * source[41] + c131 * source[43] + c130 * source[41]
                  - c129 * source[43] - c130 * source[61] + c129 * source[63];
    target[58] =  c105 * source[224] - c105 * source[226] - c105 * source[244]
                  + c105 * source[246] - c105 * source[134] + c105 * source[136]
                  + c105 * source[154] - c105 * source[156] - c105 * source[154]
                  + c105 * source[156] + c105 * source[174] - c105 * source[176]
                  + c132 * source[4] - c132 * source[6] - c132 * source[24]
                  + c132 * source[26] + c133 * source[24] - c133 * source[26]
                  - c133 * source[44] + c133 * source[46] + c132 * source[44]
                  - c132 * source[46] - c132 * source[64] + c132 * source[66];
    target[59] =  c109 * source[225] - c109 * source[245] - c109 * source[135]
                  + c109 * source[155] - c109 * source[155] + c109 * source[175]
                  + c133 * source[5] - c133 * source[25] + c134 * source[25]
                  - c134 * source[45] + c133 * source[45] - c133 * source[65];
    target[60] =  c113 * source[227] - c114 * source[220] - c114 * source[222]
                  - c113 * source[247] + c114 * source[240] + c114 * source[242]
                  - c113 * source[137] + c114 * source[130] + c114 * source[132]
                  + c113 * source[157] - c114 * source[150] - c114 * source[152]
                  - c113 * source[157] + c114 * source[150] + c114 * source[152]
                  + c113 * source[177] - c114 * source[170] - c114 * source[172]
                  + c135 * source[7] - c136 * source[0] - c136 * source[2]
                  - c135 * source[27] + c136 * source[20] + c136 * source[22]
                  + c137 * source[27] - c138 * source[20] - c138 * source[22]
                  - c137 * source[47] + c138 * source[40] + c138 * source[42]
                  + c135 * source[47] - c136 * source[40] - c136 * source[42]
                  - c135 * source[67] + c136 * source[60] + c136 * source[62];
    target[61] =  c113 * source[228] - c114 * source[221] - c114 * source[223]
                  - c113 * source[248] + c114 * source[241] + c114 * source[243]
                  - c113 * source[138] + c114 * source[131] + c114 * source[133]
                  + c113 * source[158] - c114 * source[151] - c114 * source[153]
                  - c113 * source[158] + c114 * source[151] + c114 * source[153]
                  + c113 * source[178] - c114 * source[171] - c114 * source[173]
                  + c135 * source[8] - c136 * source[1] - c136 * source[3]
                  - c135 * source[28] + c136 * source[21] + c136 * source[23]
                  + c137 * source[28] - c138 * source[21] - c138 * source[23]
                  - c137 * source[48] + c138 * source[41] + c138 * source[43]
                  + c135 * source[48] - c136 * source[41] - c136 * source[43]
                  - c135 * source[68] + c136 * source[61] + c136 * source[63];
    target[62] =  c121 * source[229] - c122 * source[224] - c122 * source[226]
                  - c121 * source[249] + c122 * source[244] + c122 * source[246]
                  - c121 * source[139] + c122 * source[134] + c122 * source[136]
                  + c121 * source[159] - c122 * source[154] - c122 * source[156]
                  - c121 * source[159] + c122 * source[154] + c122 * source[156]
                  + c121 * source[179] - c122 * source[174] - c122 * source[176]
                  + c139 * source[9] - c140 * source[4] - c140 * source[6]
                  - c139 * source[29] + c140 * source[24] + c140 * source[26]
                  + c141 * source[29] - c142 * source[24] - c142 * source[26]
                  - c141 * source[49] + c142 * source[44] + c142 * source[46]
                  + c139 * source[49] - c140 * source[44] - c140 * source[46]
                  - c139 * source[69] + c140 * source[64] + c140 * source[66];
    target[63] =  c143 * source[230] - c144 * source[232] - c143 * source[140]
                  + c144 * source[142] - c143 * source[160] + c144 * source[162]
                  + c131 * source[10] - c102 * source[12] + c145 * source[30]
                  - c146 * source[32] + c131 * source[50] - c102 * source[52];
    target[64] =  c144 * source[231] - c143 * source[233] - c144 * source[141]
                  + c143 * source[143] - c144 * source[161] + c143 * source[163]
                  + c102 * source[11] - c131 * source[13] + c146 * source[31]
                  - c145 * source[33] + c102 * source[51] - c131 * source[53];
    target[65] =  c109 * source[234] - c109 * source[236] - c109 * source[144]
                  + c109 * source[146] - c109 * source[164] + c109 * source[166]
                  + c133 * source[14] - c133 * source[16] + c134 * source[34]
                  - c134 * source[36] + c133 * source[54] - c133 * source[56];
    target[66] =  c147 * source[235] - c147 * source[145] - c147 * source[165]
                  + c134 * source[15] + c148 * source[35] + c134 * source[55];
    target[67] =  c149 * source[237] - c150 * source[230] - c150 * source[232]
                  - c149 * source[147] + c150 * source[140] + c150 * source[142]
                  - c149 * source[167] + c150 * source[160] + c150 * source[162]
                  + c137 * source[17] - c138 * source[10] - c138 * source[12]
                  + c114 * source[37] - c135 * source[30] - c135 * source[32]
                  + c137 * source[57] - c138 * source[50] - c138 * source[52];
    target[68] =  c149 * source[238] - c150 * source[231] - c150 * source[233]
                  - c149 * source[148] + c150 * source[141] + c150 * source[143]
                  - c149 * source[168] + c150 * source[161] + c150 * source[163]
                  + c137 * source[18] - c138 * source[11] - c138 * source[13]
                  + c114 * source[38] - c135 * source[31] - c135 * source[33]
                  + c137 * source[58] - c138 * source[51] - c138 * source[53];
    target[69] =  c151 * source[239] - c123 * source[234] - c123 * source[236]
                  - c151 * source[149] + c123 * source[144] + c123 * source[146]
                  - c151 * source[169] + c123 * source[164] + c123 * source[166]
                  + c141 * source[19] - c142 * source[14] - c142 * source[16]
                  + c152 * source[39] - c125 * source[34] - c125 * source[36]
                  + c141 * source[59] - c142 * source[54] - c142 * source[56];
    target[70] =  c153 * source[250] - c122 * source[252] - c154 * source[180]
                  + c155 * source[182] - c154 * source[200] + c155 * source[202]
                  + c156 * source[70] - c157 * source[72] + c158 * source[90]
                  - c159 * source[92] + c156 * source[110] - c157 * source[112];
    target[71] =  c122 * source[251] - c153 * source[253] - c155 * source[181]
                  + c154 * source[183] - c155 * source[201] + c154 * source[203]
                  + c157 * source[71] - c156 * source[73] + c159 * source[91]
                  - c158 * source[93] + c157 * source[111] - c156 * source[113];
    target[72] =  c150 * source[254] - c150 * source[256] - c160 * source[184]
                  + c160 * source[186] - c160 * source[204] + c160 * source[206]
                  + c161 * source[74] - c161 * source[76] + c162 * source[94]
                  - c162 * source[96] + c161 * source[114] - c161 * source[116];
    target[73] =  c113 * source[255] - c163 * source[185] - c163 * source[205]
                  + c162 * source[75] + c160 * source[95] + c162 * source[115];
    target[74] =  c164 * source[257] - c165 * source[250] - c165 * source[252]
                  - c109 * source[187] + c148 * source[180] + c148 * source[182]
                  - c109 * source[207] + c148 * source[200] + c148 * source[202]
                  + c148 * source[77] - c133 * source[70] - c133 * source[72]
                  + c105 * source[97] - c134 * source[90] - c134 * source[92]
                  + c148 * source[117] - c133 * source[110] - c133 * source[112];
    target[75] =  c164 * source[258] - c165 * source[251] - c165 * source[253]
                  - c109 * source[188] + c148 * source[181] + c148 * source[183]
                  - c109 * source[208] + c148 * source[201] + c148 * source[203]
                  + c148 * source[78] - c133 * source[71] - c133 * source[73]
                  + c105 * source[98] - c134 * source[91] - c134 * source[93]
                  + c148 * source[118] - c133 * source[111] - c133 * source[113];
    target[76] =  c166 * source[259] - c167 * source[254] - c167 * source[256]
                  - c143 * source[189] + c100 * source[184] + c100 * source[186]
                  - c143 * source[209] + c100 * source[204] + c100 * source[206]
                  + c168 * source[79] - c146 * source[74] - c146 * source[76]
                  + c99 * source[99] - c169 * source[94] - c169 * source[96]
                  + c168 * source[119] - c146 * source[114] - c146 * source[116];
    target[77] =  c153 * source[260] - c122 * source[262] - c154 * source[190]
                  + c155 * source[192] - c154 * source[210] + c155 * source[212]
                  + c156 * source[80] - c157 * source[82] + c158 * source[100]
                  - c159 * source[102] + c156 * source[120] - c157 * source[122];
    target[78] =  c122 * source[261] - c153 * source[263] - c155 * source[191]
                  + c154 * source[193] - c155 * source[211] + c154 * source[213]
                  + c157 * source[81] - c156 * source[83] + c159 * source[101]
                  - c158 * source[103] + c157 * source[121] - c156 * source[123];
    target[79] =  c150 * source[264] - c150 * source[266] - c160 * source[194]
                  + c160 * source[196] - c160 * source[214] + c160 * source[216]
                  + c161 * source[84] - c161 * source[86] + c162 * source[104]
                  - c162 * source[106] + c161 * source[124] - c161 * source[126];
    target[80] =  c113 * source[265] - c163 * source[195] - c163 * source[215]
                  + c162 * source[85] + c160 * source[105] + c162 * source[125];
    target[81] =  c164 * source[267] - c165 * source[260] - c165 * source[262]
                  - c109 * source[197] + c148 * source[190] + c148 * source[192]
                  - c109 * source[217] + c148 * source[210] + c148 * source[212]
                  + c148 * source[87] - c133 * source[80] - c133 * source[82]
                  + c105 * source[107] - c134 * source[100] - c134 * source[102]
                  + c148 * source[127] - c133 * source[120] - c133 * source[122];
    target[82] =  c164 * source[268] - c165 * source[261] - c165 * source[263]
                  - c109 * source[198] + c148 * source[191] + c148 * source[193]
                  - c109 * source[218] + c148 * source[211] + c148 * source[213]
                  + c148 * source[88] - c133 * source[81] - c133 * source[83]
                  + c105 * source[108] - c134 * source[101] - c134 * source[103]
                  + c148 * source[128] - c133 * source[121] - c133 * source[123];
    target[83] =  c166 * source[269] - c167 * source[264] - c167 * source[266]
                  - c143 * source[199] + c100 * source[194] + c100 * source[196]
                  - c143 * source[219] + c100 * source[214] + c100 * source[216]
                  + c168 * source[89] - c146 * source[84] - c146 * source[86]
                  + c99 * source[109] - c169 * source[104] - c169 * source[106]
                  + c168 * source[129] - c146 * source[124] - c146 * source[126];
    target[84] =  c170 * source[270] - c171 * source[272] - c172 * source[220]
                  + c173 * source[222] - c172 * source[240] + c173 * source[242]
                  + c174 * source[130] - c175 * source[132] + c176 * source[150]
                  - c177 * source[152] + c174 * source[170] - c175 * source[172]
                  - c178 * source[0] + c179 * source[2] - c179 * source[20]
                  + c180 * source[22] - c179 * source[40] + c180 * source[42]
                  - c178 * source[60] + c179 * source[62];
    target[85] =  c171 * source[271] - c170 * source[273] - c173 * source[221]
                  + c172 * source[223] - c173 * source[241] + c172 * source[243]
                  + c175 * source[131] - c174 * source[133] + c177 * source[151]
                  - c176 * source[153] + c175 * source[171] - c174 * source[173]
                  - c179 * source[1] + c178 * source[3] - c180 * source[21]
                  + c179 * source[23] - c180 * source[41] + c179 * source[43]
                  - c179 * source[61] + c178 * source[63];
    target[86] =  c181 * source[274] - c181 * source[276] - c182 * source[224]
                  + c182 * source[226] - c182 * source[244] + c182 * source[246]
                  + c183 * source[134] - c183 * source[136] + c184 * source[154]
                  - c184 * source[156] + c183 * source[174] - c183 * source[176]
                  - c185 * source[4] + c185 * source[6] - c186 * source[24]
                  + c186 * source[26] - c186 * source[44] + c186 * source[46]
                  - c185 * source[64] + c185 * source[66];
    target[87] =  c187 * source[275] - c188 * source[225] - c188 * source[245]
                  + c184 * source[135] + c189 * source[155] + c184 * source[175]
                  - c190 * source[5] - c191 * source[25] - c191 * source[45]
                  - c190 * source[65];
    target[88] =  c192 * source[277] - c193 * source[270] - c193 * source[272]
                  - c194 * source[227] + c195 * source[220] + c195 * source[222]
                  - c194 * source[247] + c195 * source[240] + c195 * source[242]
                  + c196 * source[137] - c197 * source[130] - c197 * source[132]
                  + c198 * source[157] - c199 * source[150] - c199 * source[152]
                  + c196 * source[177] - c197 * source[170] - c197 * source[172]
                  - c200 * source[7] + c201 * source[0] + c201 * source[2]
                  - c202 * source[27] + c203 * source[20] + c203 * source[22]
                  - c202 * source[47] + c203 * source[40] + c203 * source[42]
                  - c200 * source[67] + c201 * source[60] + c201 * source[62];
    target[89] =  c192 * source[278] - c193 * source[271] - c193 * source[273]
                  - c194 * source[228] + c195 * source[221] + c195 * source[223]
                  - c194 * source[248] + c195 * source[241] + c195 * source[243]
                  + c196 * source[138] - c197 * source[131] - c197 * source[133]
                  + c198 * source[158] - c199 * source[151] - c199 * source[153]
                  + c196 * source[178] - c197 * source[171] - c197 * source[173]
                  - c200 * source[8] + c201 * source[1] + c201 * source[3]
                  - c202 * source[28] + c203 * source[21] + c203 * source[23]
                  - c202 * source[48] + c203 * source[41] + c203 * source[43]
                  - c200 * source[68] + c201 * source[61] + c201 * source[63];
    target[90] =  source[279] - c204 * source[274] - c204 * source[276]
                  - c205 * source[229] + c206 * source[224] + c206 * source[226]
                  - c205 * source[249] + c206 * source[244] + c206 * source[246]
                  + c207 * source[139] - c208 * source[134] - c208 * source[136]
                  + c206 * source[159] - c209 * source[154] - c209 * source[156]
                  + c207 * source[179] - c208 * source[174] - c208 * source[176]
                  - c210 * source[9] + c211 * source[4] + c211 * source[6]
                  - c212 * source[29] + c213 * source[24] + c213 * source[26]
                  - c212 * source[49] + c213 * source[44] + c213 * source[46]
                  - c210 * source[69] + c211 * source[64] + c211 * source[66];
  }
}

