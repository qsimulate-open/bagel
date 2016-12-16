//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_54.cc
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


void CarSphList::carsph_54(const int nloop, const double* source, double* target) {
  const double c56 = 89.294106748429925;
  const double c95 = 84.187291202413675;
  const double c50 = 83.526988527660933;
  const double c86 = 78.75;
  const double c142 = 68.738635424337602;
  const double c141 = 64.299105748058423;
  const double c70 = 59.529404498953291;
  const double c44 = 59.0625;
  const double c67 = 55.68465901844062;
  const double c23 = 47.062126492541751;
  const double c138 = 45.46633369868303;
  const double c53 = 44.647053374214963;
  const double c12 = 44.022588307027107;
  const double c57 = 42.093645601206838;
  const double c62 = 39.933538535922409;
  const double c98 = 39.686269665968858;
  const double c46 = 39.375;
  const double c166 = 38.97114317029974;
  const double c105 = 37.649701194033398;
  const double c83 = 37.123106012293746;
  const double c153 = 36.454166428544212;
  const double c125 = 34.369317712168801;
  const double c212 = 33.541019662496844;
  const double c144 = 32.403703492039298;
  const double c118 = 32.149552874029212;
  const double c58 = 31.57023420090513;
  const double c201 = 31.374750995027831;
  const double c3 = 31.128670717282485;
  const double c147 = 30.740852297878796;
  const double c139 = 30.310889132455351;
  const double c68 = 29.764702249476645;
  const double c72 = 28.062430400804562;
  const double c49 = 27.84232950922031;
  const double c74 = 26.622359023948274;
  const double c65 = 26.25;
  const double c164 = 25.98076211353316;
  const double c149 = 25.776988284126599;
  const double c214 = 25.155764746872634;
  const double c145 = 24.302777619029477;
  const double c17 = 23.531063246270875;
  const double c112 = 22.733166849341515;
  const double c27 = 22.18529918662356;
  const double c14 = 22.011294153513553;
  const double c140 = 21.433035249352809;
  const double c36 = 21.046822800603419;
  const double c7 = 20.75244714485499;
  const double c159 = 19.48557158514987;
  const double c66 = 18.561553006146873;
  const double c171 = 18.371173070873837;
  const double c182 = 17.428425057933374;
  const double c121 = 17.1846588560844;
  const double c204 = 16.770509831248422;
  const double c28 = 16.638974389967672;
  const double c127 = 16.201851746019649;
  const double c120 = 16.074776437014606;
  const double c218 = 15.811388300841896;
  const double c5 = 15.564335358641243;
  const double c132 = 15.370426148939398;
  const double c115 = 15.155444566227676;
  const double c227 = 15;
  const double c54 = 14.882351124738323;
  const double c196 = 14.79019945774904;
  const double c11 = 14.674196102342369;
  const double c89 = 14.031215200402281;
  const double c48 = 13.921164754610155;
  const double c172 = 13.778379803155376;
  const double c61 = 13.311179511974137;
  const double c97 = 13.228756555322953;
  const double c157 = 12.99038105676658;
  const double c208 = 12.577882373436317;
  const double c102 = 12.549900398011133;
  const double c82 = 12.374368670764582;
  const double c169 = 12.24744871391589;
  const double c128 = 12.151388809514739;
  const double c219 = 11.858541225631422;
  const double c19 = 11.765531623135438;
  const double c178 = 11.61895003862225;
  const double c143 = 11.456439237389599;
  const double c114 = 11.366583424670758;
  const double c233 = 11.25;
  const double c29 = 11.09264959331178;
  const double c117 = 10.716517624676404;
  const double c40 = 10.523411400301709;
  const double c200 = 10.458250331675945;
  const double c8 = 10.376223572427495;
  const double c146 = 10.246950765959598;
  const double c64 = 9.9833846339806023;
  const double c71 = 9.9215674164922145;
  const double c43 = 9.84375;
  const double c106 = 9.4124252985083494;
  const double c77 = 9.2807765030734366;
  const double c170 = 9.1855865354369186;
  const double c222 = 8.8939059192235668;
  const double c73 = 8.8741196746494246;
  const double c85 = 8.75;
  const double c123 = 8.5923294280422002;
  const double c30 = 8.3194871949838358;
  const double c129 = 8.1009258730098246;
  const double c24 = 7.8436877487569578;
  const double c136 = 7.6852130744696989;
  const double c116 = 7.5777222831138378;
  const double c51 = 7.4411755623691613;
  const double c13 = 7.3370980511711847;
  const double c35 = 7.0156076002011405;
  const double c81 = 6.9605823773050775;
  const double c210 = 6.7082039324993694;
  const double c59 = 6.6555897559870685;
  const double c45 = 6.5625;
  const double c167 = 6.49519052838329;
  const double c206 = 6.2889411867181586;
  const double c199 = 6.2749501990055663;
  const double c130 = 6.0756944047573693;
  const double c220 = 5.9292706128157109;
  const double c181 = 5.8094750193111251;
  const double c126 = 5.7282196186947996;
  const double c230 = 5.625;
  const double c213 = 5.5901699437494745;
  const double c197 = 5.54632479665589;
  const double c119 = 5.3582588123382022;
  const double c38 = 5.2617057001508547;
  const double c2 = 5.1881117862137476;
  const double c131 = 5.123475382979799;
  const double c226 = 5;
  const double c63 = 4.9916923169903011;
  const double c69 = 4.9607837082461073;
  const double c21 = 4.7062126492541747;
  const double c94 = 4.677071733467427;
  const double c47 = 4.6403882515367183;
  const double c221 = 4.4469529596117834;
  const double c191 = 4.4370598373247123;
  const double c10 = 4.4022588307027108;
  const double c184 = 4.3571062644833436;
  const double c165 = 4.3301270189221936;
  const double c148 = 4.2961647140211001;
  const double c215 = 4.1926274578121054;
  const double c101 = 4.1833001326703778;
  const double c18 = 3.9218438743784789;
  const double c177 = 3.872983346207417;
  const double c134 = 3.8426065372348495;
  const double c111 = 3.7888611415569189;
  const double c229 = 3.75;
  const double c100 = 3.7205877811845807;
  const double c192 = 3.69754986443726;
  const double c39 = 3.5078038001005702;
  const double c202 = 3.3541019662496847;
  const double c75 = 3.3277948779935342;
  const double c88 = 3.28125;
  const double c160 = 3.247595264191645;
  const double c216 = 3.1622776601683795;
  const double c104 = 3.1374750995027831;
  const double c1 = 3.1128670717282483;
  const double c76 = 3.0935921676911455;
  const double c175 = 3.0618621784789726;
  const double c155 = 3.0378472023786847;
  const double c223 = 3;
  const double c195 = 2.9580398915498081;
  const double c180 = 2.9047375096555625;
  const double c122 = 2.8641098093473998;
  const double c234 = 2.8125;
  const double c205 = 2.7950849718747373;
  const double c194 = 2.773162398327945;
  const double c37 = 2.6308528500754274;
  const double c4 = 2.5940558931068738;
  const double c135 = 2.5617376914898995;
  const double c55 = 2.4803918541230536;
  const double c217 = 2.3717082451262845;
  const double c15 = 2.3531063246270874;
  const double c90 = 2.3385358667337135;
  const double c79 = 2.3201941257683592;
  const double c176 = 2.2963966338592297;
  const double c25 = 2.2185299186623562;
  const double c183 = 2.1785531322416718;
  const double c158 = 2.1650635094610968;
  const double c151 = 2.1480823570105501;
  const double c32 = 2.104682280060342;
  const double c209 = 2.0963137289060527;
  const double c198 = 2.0916500663351889;
  const double c6 = 2.0752447144854989;
  const double c156 = 2.0252314682524561;
  const double c20 = 1.9609219371892395;
  const double c133 = 1.9213032686174247;
  const double c113 = 1.8944305707784594;
  const double c228 = 1.875;
  const double c91 = 1.7539019000502851;
  const double c26 = 1.6638974389967671;
  const double c99 = 1.6535945694153691;
  const double c42 = 1.640625;
  const double c161 = 1.6237976320958225;
  const double c103 = 1.5687375497513916;
  const double c84 = 1.5467960838455728;
  const double c173 = 1.5309310892394863;
  const double c9 = 1.4674196102342369;
  const double c179 = 1.4523687548277813;
  const double c124 = 1.4320549046736999;
  const double c232 = 1.40625;
  const double c193 = 1.3865811991639725;
  const double c41 = 1.3154264250377137;
  const double c52 = 1.2401959270615268;
  const double c110 = 1.1765531623135437;
  const double c80 = 1.1600970628841796;
  const double c174 = 1.1481983169296148;
  const double c211 = 1.1180339887498949;
  const double c87 = 1.09375;
  const double c168 = 1.0825317547305484;
  const double c207 = 1.0481568644530264;
  const double c154 = 1.0126157341262281;
  const double c188 = 0.96824583655185426;
  const double c137 = 0.96065163430871237;
  const double c93 = 0.87695095002514256;
  const double c60 = 0.83194871949838356;
  const double c22 = 0.78436877487569578;
  const double c225 = 0.75;
  const double c190 = 0.73950997288745202;
  const double c189 = 0.72618437741389064;
  const double c152 = 0.71602745233684995;
  const double c231 = 0.703125;
  const double c31 = 0.70156076002011403;
  const double c109 = 0.58827658115677184;
  const double c96 = 0.58463396668342837;
  const double c203 = 0.55901699437494745;
  const double c163 = 0.54126587736527421;
  const double c34 = 0.52617057001508549;
  const double c107 = 0.52291251658379723;
  const double c0 = 0.51881117862137471;
  const double c185 = 0.48412291827592713;
  const double c16 = 0.39218438743784789;
  const double c78 = 0.38669902096139319;
  const double c224 = 0.375;
  const double c187 = 0.36309218870694532;
  const double c150 = 0.35801372616842497;
  const double c92 = 0.29231698334171419;
  const double c162 = 0.2706329386826371;
  const double c33 = 0.26308528500754275;
  const double c108 = 0.19609219371892395;
  const double c186 = 0.18154609435347266;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 99, source += 315) {
    target[0] =  c0 * source[0] - c1 * source[2] + c0 * source[4]
                  - c2 * source[30] + c3 * source[32] - c2 * source[34]
                  + c4 * source[60] - c5 * source[62] + c4 * source[64];
    target[1] =  c6 * source[1] - c6 * source[3] - c7 * source[31]
                  + c7 * source[33] + c8 * source[61] - c8 * source[63];
    target[2] =  c9 * source[5] - c10 * source[7] - c11 * source[35]
                  + c12 * source[37] + c13 * source[65] - c14 * source[67];
    target[3] =  c10 * source[6] - c9 * source[8] - c12 * source[36]
                  + c11 * source[38] + c14 * source[66] - c13 * source[68];
    target[4] =  c15 * source[9] - c15 * source[11] - c16 * source[0]
                  + c16 * source[2] - c16 * source[2] + c16 * source[4]
                  - c17 * source[39] + c17 * source[41] + c18 * source[30]
                  - c18 * source[32] + c18 * source[32] - c18 * source[34]
                  + c19 * source[69] - c19 * source[71] - c20 * source[60]
                  + c20 * source[62] - c20 * source[62] + c20 * source[64];
    target[5] =  c21 * source[10] - c22 * source[1] - c22 * source[3]
                  - c23 * source[40] + c24 * source[31] + c24 * source[33]
                  + c17 * source[70] - c18 * source[61] - c18 * source[63];
    target[6] =  c25 * source[12] - c26 * source[5] - c26 * source[7]
                  - c27 * source[42] + c28 * source[35] + c28 * source[37]
                  + c29 * source[72] - c30 * source[65] - c30 * source[67];
    target[7] =  c25 * source[13] - c26 * source[6] - c26 * source[8]
                  - c27 * source[43] + c28 * source[36] + c28 * source[38]
                  + c29 * source[73] - c30 * source[66] - c30 * source[68];
    target[8] =  c31 * source[14] - c32 * source[9] - c32 * source[11]
                  + c33 * source[0] + c34 * source[2] + c33 * source[4]
                  - c35 * source[44] + c36 * source[39] + c36 * source[41]
                  - c37 * source[30] - c38 * source[32] - c37 * source[34]
                  + c39 * source[74] - c40 * source[69] - c40 * source[71]
                  + c41 * source[60] + c37 * source[62] + c41 * source[64];
    target[9] =  c4 * source[15] - c5 * source[17] + c4 * source[19]
                  - c2 * source[45] + c3 * source[47] - c2 * source[49]
                  + c0 * source[75] - c1 * source[77] + c0 * source[79];
    target[10] =  c8 * source[16] - c8 * source[18] - c7 * source[46]
                  + c7 * source[48] + c6 * source[76] - c6 * source[78];
    target[11] =  c13 * source[20] - c14 * source[22] - c11 * source[50]
                  + c12 * source[52] + c9 * source[80] - c10 * source[82];
    target[12] =  c14 * source[21] - c13 * source[23] - c12 * source[51]
                  + c11 * source[53] + c10 * source[81] - c9 * source[83];
    target[13] =  c19 * source[24] - c19 * source[26] - c20 * source[15]
                  + c20 * source[17] - c20 * source[17] + c20 * source[19]
                  - c17 * source[54] + c17 * source[56] + c18 * source[45]
                  - c18 * source[47] + c18 * source[47] - c18 * source[49]
                  + c15 * source[84] - c15 * source[86] - c16 * source[75]
                  + c16 * source[77] - c16 * source[77] + c16 * source[79];
    target[14] =  c17 * source[25] - c18 * source[16] - c18 * source[18]
                  - c23 * source[55] + c24 * source[46] + c24 * source[48]
                  + c21 * source[85] - c22 * source[76] - c22 * source[78];
    target[15] =  c29 * source[27] - c30 * source[20] - c30 * source[22]
                  - c27 * source[57] + c28 * source[50] + c28 * source[52]
                  + c25 * source[87] - c26 * source[80] - c26 * source[82];
    target[16] =  c29 * source[28] - c30 * source[21] - c30 * source[23]
                  - c27 * source[58] + c28 * source[51] + c28 * source[53]
                  + c25 * source[88] - c26 * source[81] - c26 * source[83];
    target[17] =  c39 * source[29] - c40 * source[24] - c40 * source[26]
                  + c41 * source[15] + c37 * source[17] + c41 * source[19]
                  - c35 * source[59] + c36 * source[54] + c36 * source[56]
                  - c37 * source[45] - c38 * source[47] - c37 * source[49]
                  + c31 * source[89] - c32 * source[84] - c32 * source[86]
                  + c33 * source[75] + c34 * source[77] + c33 * source[79];
    target[18] =  c42 * source[90] - c43 * source[92] + c42 * source[94]
                  - c43 * source[120] + c44 * source[122] - c43 * source[124]
                  + c42 * source[150] - c43 * source[152] + c42 * source[154];
    target[19] =  c45 * source[91] - c45 * source[93] - c46 * source[121]
                  + c46 * source[123] + c45 * source[151] - c45 * source[153];
    target[20] =  c47 * source[95] - c48 * source[97] - c49 * source[125]
                  + c50 * source[127] + c47 * source[155] - c48 * source[157];
    target[21] =  c48 * source[96] - c47 * source[98] - c50 * source[126]
                  + c49 * source[128] + c48 * source[156] - c47 * source[158];
    target[22] =  c51 * source[99] - c51 * source[101] - c52 * source[90]
                  + c52 * source[92] - c52 * source[92] + c52 * source[94]
                  - c53 * source[129] + c53 * source[131] + c51 * source[120]
                  - c51 * source[122] + c51 * source[122] - c51 * source[124]
                  + c51 * source[159] - c51 * source[161] - c52 * source[150]
                  + c52 * source[152] - c52 * source[152] + c52 * source[154];
    target[23] =  c54 * source[100] - c55 * source[91] - c55 * source[93]
                  - c56 * source[130] + c54 * source[121] + c54 * source[123]
                  + c54 * source[160] - c55 * source[151] - c55 * source[153];
    target[24] =  c35 * source[102] - c38 * source[95] - c38 * source[97]
                  - c57 * source[132] + c58 * source[125] + c58 * source[127]
                  + c35 * source[162] - c38 * source[155] - c38 * source[157];
    target[25] =  c35 * source[103] - c38 * source[96] - c38 * source[98]
                  - c57 * source[133] + c58 * source[126] + c58 * source[128]
                  + c35 * source[163] - c38 * source[156] - c38 * source[158];
    target[26] =  c25 * source[104] - c59 * source[99] - c59 * source[101]
                  + c60 * source[90] + c26 * source[92] + c60 * source[94]
                  - c61 * source[134] + c62 * source[129] + c62 * source[131]
                  - c63 * source[120] - c64 * source[122] - c63 * source[124]
                  + c25 * source[164] - c59 * source[159] - c59 * source[161]
                  + c60 * source[150] + c26 * source[152] + c60 * source[154];
    target[27] =  c45 * source[105] - c46 * source[107] + c45 * source[109]
                  - c45 * source[135] + c46 * source[137] - c45 * source[139];
    target[28] =  c65 * source[106] - c65 * source[108] - c65 * source[136]
                  + c65 * source[138];
    target[29] =  c66 * source[110] - c67 * source[112] - c66 * source[140]
                  + c67 * source[142];
    target[30] =  c67 * source[111] - c66 * source[113] - c67 * source[141]
                  + c66 * source[143];
    target[31] =  c68 * source[114] - c68 * source[116] - c69 * source[105]
                  + c69 * source[107] - c69 * source[107] + c69 * source[109]
                  - c68 * source[144] + c68 * source[146] + c69 * source[135]
                  - c69 * source[137] + c69 * source[137] - c69 * source[139];
    target[32] =  c70 * source[115] - c71 * source[106] - c71 * source[108]
                  - c70 * source[145] + c71 * source[136] + c71 * source[138];
    target[33] =  c72 * source[117] - c36 * source[110] - c36 * source[112]
                  - c72 * source[147] + c36 * source[140] + c36 * source[142];
    target[34] =  c72 * source[118] - c36 * source[111] - c36 * source[113]
                  - c72 * source[148] + c36 * source[141] + c36 * source[143];
    target[35] =  c73 * source[119] - c74 * source[114] - c74 * source[116]
                  + c75 * source[105] + c59 * source[107] + c75 * source[109]
                  - c73 * source[149] + c74 * source[144] + c74 * source[146]
                  - c75 * source[135] - c59 * source[137] - c75 * source[139];
    target[36] =  c76 * source[165] - c66 * source[167] + c76 * source[169]
                  - c77 * source[195] + c67 * source[197] - c77 * source[199]
                  - c78 * source[0] + c79 * source[2] - c78 * source[4]
                  + c80 * source[30] - c81 * source[32] + c80 * source[34]
                  - c78 * source[30] + c79 * source[32] - c78 * source[34]
                  + c80 * source[60] - c81 * source[62] + c80 * source[64];
    target[37] =  c82 * source[166] - c82 * source[168] - c83 * source[196]
                  + c83 * source[198] - c84 * source[1] + c84 * source[3]
                  + c47 * source[31] - c47 * source[33] - c84 * source[31]
                  + c84 * source[33] + c47 * source[61] - c47 * source[63];
    target[38] =  c85 * source[170] - c65 * source[172] - c65 * source[200]
                  + c86 * source[202] - c87 * source[5] + c88 * source[7]
                  + c88 * source[35] - c43 * source[37] - c87 * source[35]
                  + c88 * source[37] + c88 * source[65] - c43 * source[67];
    target[39] =  c65 * source[171] - c85 * source[173] - c86 * source[201]
                  + c65 * source[203] - c88 * source[6] + c87 * source[8]
                  + c43 * source[36] - c88 * source[38] - c88 * source[36]
                  + c87 * source[38] + c43 * source[66] - c88 * source[68];
    target[40] =  c89 * source[174] - c89 * source[176] - c90 * source[165]
                  + c90 * source[167] - c90 * source[167] + c90 * source[169]
                  - c57 * source[204] + c57 * source[206] + c35 * source[195]
                  - c35 * source[197] + c35 * source[197] - c35 * source[199]
                  - c91 * source[9] + c91 * source[11] + c92 * source[0]
                  - c92 * source[2] + c92 * source[2] - c92 * source[4]
                  + c38 * source[39] - c38 * source[41] - c93 * source[30]
                  + c93 * source[32] - c93 * source[32] + c93 * source[34]
                  - c91 * source[39] + c91 * source[41] + c92 * source[30]
                  - c92 * source[32] + c92 * source[32] - c92 * source[34]
                  + c38 * source[69] - c38 * source[71] - c93 * source[60]
                  + c93 * source[62] - c93 * source[62] + c93 * source[64];
    target[41] =  c72 * source[175] - c94 * source[166] - c94 * source[168]
                  - c95 * source[205] + c89 * source[196] + c89 * source[198]
                  - c39 * source[10] + c96 * source[1] + c96 * source[3]
                  + c40 * source[40] - c91 * source[31] - c91 * source[33]
                  - c39 * source[40] + c96 * source[31] + c96 * source[33]
                  + c40 * source[70] - c91 * source[61] - c91 * source[63];
    target[42] =  c97 * source[177] - c71 * source[170] - c71 * source[172]
                  - c98 * source[207] + c68 * source[200] + c68 * source[202]
                  - c99 * source[12] + c52 * source[5] + c52 * source[7]
                  + c69 * source[42] - c100 * source[35] - c100 * source[37]
                  - c99 * source[42] + c52 * source[35] + c52 * source[37]
                  + c69 * source[72] - c100 * source[65] - c100 * source[67];
    target[43] =  c97 * source[178] - c71 * source[171] - c71 * source[173]
                  - c98 * source[208] + c68 * source[201] + c68 * source[203]
                  - c99 * source[13] + c52 * source[6] + c52 * source[8]
                  + c69 * source[43] - c100 * source[36] - c100 * source[38]
                  - c99 * source[43] + c52 * source[36] + c52 * source[38]
                  + c69 * source[73] - c100 * source[66] - c100 * source[68];
    target[44] =  c101 * source[179] - c102 * source[174] - c102 * source[176]
                  + c103 * source[165] + c104 * source[167] + c103 * source[169]
                  - c102 * source[209] + c105 * source[204] + c105 * source[206]
                  - c21 * source[195] - c106 * source[197] - c21 * source[199]
                  - c107 * source[14] + c103 * source[9] + c103 * source[11]
                  - c108 * source[0] - c16 * source[2] - c108 * source[4]
                  + c103 * source[44] - c21 * source[39] - c21 * source[41]
                  + c109 * source[30] + c110 * source[32] + c109 * source[34]
                  - c107 * source[44] + c103 * source[39] + c103 * source[41]
                  - c108 * source[30] - c16 * source[32] - c108 * source[34]
                  + c103 * source[74] - c21 * source[69] - c21 * source[71]
                  + c109 * source[60] + c110 * source[62] + c109 * source[64];
    target[45] =  c77 * source[180] - c67 * source[182] + c77 * source[184]
                  - c76 * source[210] + c66 * source[212] - c76 * source[214]
                  - c80 * source[15] + c81 * source[17] - c80 * source[19]
                  + c78 * source[45] - c79 * source[47] + c78 * source[49]
                  - c80 * source[45] + c81 * source[47] - c80 * source[49]
                  + c78 * source[75] - c79 * source[77] + c78 * source[79];
    target[46] =  c83 * source[181] - c83 * source[183] - c82 * source[211]
                  + c82 * source[213] - c47 * source[16] + c47 * source[18]
                  + c84 * source[46] - c84 * source[48] - c47 * source[46]
                  + c47 * source[48] + c84 * source[76] - c84 * source[78];
    target[47] =  c65 * source[185] - c86 * source[187] - c85 * source[215]
                  + c65 * source[217] - c88 * source[20] + c43 * source[22]
                  + c87 * source[50] - c88 * source[52] - c88 * source[50]
                  + c43 * source[52] + c87 * source[80] - c88 * source[82];
    target[48] =  c86 * source[186] - c65 * source[188] - c65 * source[216]
                  + c85 * source[218] - c43 * source[21] + c88 * source[23]
                  + c88 * source[51] - c87 * source[53] - c43 * source[51]
                  + c88 * source[53] + c88 * source[81] - c87 * source[83];
    target[49] =  c57 * source[189] - c57 * source[191] - c35 * source[180]
                  + c35 * source[182] - c35 * source[182] + c35 * source[184]
                  - c89 * source[219] + c89 * source[221] + c90 * source[210]
                  - c90 * source[212] + c90 * source[212] - c90 * source[214]
                  - c38 * source[24] + c38 * source[26] + c93 * source[15]
                  - c93 * source[17] + c93 * source[17] - c93 * source[19]
                  + c91 * source[54] - c91 * source[56] - c92 * source[45]
                  + c92 * source[47] - c92 * source[47] + c92 * source[49]
                  - c38 * source[54] + c38 * source[56] + c93 * source[45]
                  - c93 * source[47] + c93 * source[47] - c93 * source[49]
                  + c91 * source[84] - c91 * source[86] - c92 * source[75]
                  + c92 * source[77] - c92 * source[77] + c92 * source[79];
    target[50] =  c95 * source[190] - c89 * source[181] - c89 * source[183]
                  - c72 * source[220] + c94 * source[211] + c94 * source[213]
                  - c40 * source[25] + c91 * source[16] + c91 * source[18]
                  + c39 * source[55] - c96 * source[46] - c96 * source[48]
                  - c40 * source[55] + c91 * source[46] + c91 * source[48]
                  + c39 * source[85] - c96 * source[76] - c96 * source[78];
    target[51] =  c98 * source[192] - c68 * source[185] - c68 * source[187]
                  - c97 * source[222] + c71 * source[215] + c71 * source[217]
                  - c69 * source[27] + c100 * source[20] + c100 * source[22]
                  + c99 * source[57] - c52 * source[50] - c52 * source[52]
                  - c69 * source[57] + c100 * source[50] + c100 * source[52]
                  + c99 * source[87] - c52 * source[80] - c52 * source[82];
    target[52] =  c98 * source[193] - c68 * source[186] - c68 * source[188]
                  - c97 * source[223] + c71 * source[216] + c71 * source[218]
                  - c69 * source[28] + c100 * source[21] + c100 * source[23]
                  + c99 * source[58] - c52 * source[51] - c52 * source[53]
                  - c69 * source[58] + c100 * source[51] + c100 * source[53]
                  + c99 * source[88] - c52 * source[81] - c52 * source[83];
    target[53] =  c102 * source[194] - c105 * source[189] - c105 * source[191]
                  + c21 * source[180] + c106 * source[182] + c21 * source[184]
                  - c101 * source[224] + c102 * source[219] + c102 * source[221]
                  - c103 * source[210] - c104 * source[212] - c103 * source[214]
                  - c103 * source[29] + c21 * source[24] + c21 * source[26]
                  - c109 * source[15] - c110 * source[17] - c109 * source[19]
                  + c107 * source[59] - c103 * source[54] - c103 * source[56]
                  + c108 * source[45] + c16 * source[47] + c108 * source[49]
                  - c103 * source[59] + c21 * source[54] + c21 * source[56]
                  - c109 * source[45] - c110 * source[47] - c109 * source[49]
                  + c107 * source[89] - c103 * source[84] - c103 * source[86]
                  + c108 * source[75] + c16 * source[77] + c108 * source[79];
    target[54] =  c111 * source[225] - c112 * source[227] + c111 * source[229]
                  - c111 * source[255] + c112 * source[257] - c111 * source[259]
                  - c113 * source[90] + c114 * source[92] - c113 * source[94]
                  + c113 * source[120] - c114 * source[122] + c113 * source[124]
                  - c113 * source[120] + c114 * source[122] - c113 * source[124]
                  + c113 * source[150] - c114 * source[152] + c113 * source[154];
    target[55] =  c115 * source[226] - c115 * source[228] - c115 * source[256]
                  + c115 * source[258] - c116 * source[91] + c116 * source[93]
                  + c116 * source[121] - c116 * source[123] - c116 * source[121]
                  + c116 * source[123] + c116 * source[151] - c116 * source[153];
    target[56] =  c117 * source[230] - c118 * source[232] - c117 * source[260]
                  + c118 * source[262] - c119 * source[95] + c120 * source[97]
                  + c119 * source[125] - c120 * source[127] - c119 * source[125]
                  + c120 * source[127] + c119 * source[155] - c120 * source[157];
    target[57] =  c118 * source[231] - c117 * source[233] - c118 * source[261]
                  + c117 * source[263] - c120 * source[96] + c119 * source[98]
                  + c120 * source[126] - c119 * source[128] - c120 * source[126]
                  + c119 * source[128] + c120 * source[156] - c119 * source[158];
    target[58] =  c121 * source[234] - c121 * source[236] - c122 * source[225]
                  + c122 * source[227] - c122 * source[227] + c122 * source[229]
                  - c121 * source[264] + c121 * source[266] + c122 * source[255]
                  - c122 * source[257] + c122 * source[257] - c122 * source[259]
                  - c123 * source[99] + c123 * source[101] + c124 * source[90]
                  - c124 * source[92] + c124 * source[92] - c124 * source[94]
                  + c123 * source[129] - c123 * source[131] - c124 * source[120]
                  + c124 * source[122] - c124 * source[122] + c124 * source[124]
                  - c123 * source[129] + c123 * source[131] + c124 * source[120]
                  - c124 * source[122] + c124 * source[122] - c124 * source[124]
                  + c123 * source[159] - c123 * source[161] - c124 * source[150]
                  + c124 * source[152] - c124 * source[152] + c124 * source[154];
    target[59] =  c125 * source[235] - c126 * source[226] - c126 * source[228]
                  - c125 * source[265] + c126 * source[256] + c126 * source[258]
                  - c121 * source[100] + c122 * source[91] + c122 * source[93]
                  + c121 * source[130] - c122 * source[121] - c122 * source[123]
                  - c121 * source[130] + c122 * source[121] + c122 * source[123]
                  + c121 * source[160] - c122 * source[151] - c122 * source[153];
    target[60] =  c127 * source[237] - c128 * source[230] - c128 * source[232]
                  - c127 * source[267] + c128 * source[260] + c128 * source[262]
                  - c129 * source[102] + c130 * source[95] + c130 * source[97]
                  + c129 * source[132] - c130 * source[125] - c130 * source[127]
                  - c129 * source[132] + c130 * source[125] + c130 * source[127]
                  + c129 * source[162] - c130 * source[155] - c130 * source[157];
    target[61] =  c127 * source[238] - c128 * source[231] - c128 * source[233]
                  - c127 * source[268] + c128 * source[261] + c128 * source[263]
                  - c129 * source[103] + c130 * source[96] + c130 * source[98]
                  + c129 * source[133] - c130 * source[126] - c130 * source[128]
                  - c129 * source[133] + c130 * source[126] + c130 * source[128]
                  + c129 * source[163] - c130 * source[156] - c130 * source[158];
    target[62] =  c131 * source[239] - c132 * source[234] - c132 * source[236]
                  + c133 * source[225] + c134 * source[227] + c133 * source[229]
                  - c131 * source[269] + c132 * source[264] + c132 * source[266]
                  - c133 * source[255] - c134 * source[257] - c133 * source[259]
                  - c135 * source[104] + c136 * source[99] + c136 * source[101]
                  - c137 * source[90] - c133 * source[92] - c137 * source[94]
                  + c135 * source[134] - c136 * source[129] - c136 * source[131]
                  + c137 * source[120] + c133 * source[122] + c137 * source[124]
                  - c135 * source[134] + c136 * source[129] + c136 * source[131]
                  - c137 * source[120] - c133 * source[122] - c137 * source[124]
                  + c135 * source[164] - c136 * source[159] - c136 * source[161]
                  + c137 * source[150] + c133 * source[152] + c137 * source[154];
    target[63] =  c116 * source[240] - c138 * source[242] + c116 * source[244]
                  - c111 * source[105] + c112 * source[107] - c111 * source[109]
                  - c111 * source[135] + c112 * source[137] - c111 * source[139];
    target[64] =  c139 * source[241] - c139 * source[243] - c115 * source[106]
                  + c115 * source[108] - c115 * source[136] + c115 * source[138];
    target[65] =  c140 * source[245] - c141 * source[247] - c117 * source[110]
                  + c118 * source[112] - c117 * source[140] + c118 * source[142];
    target[66] =  c141 * source[246] - c140 * source[248] - c118 * source[111]
                  + c117 * source[113] - c118 * source[141] + c117 * source[143];
    target[67] =  c125 * source[249] - c125 * source[251] - c126 * source[240]
                  + c126 * source[242] - c126 * source[242] + c126 * source[244]
                  - c121 * source[114] + c121 * source[116] + c122 * source[105]
                  - c122 * source[107] + c122 * source[107] - c122 * source[109]
                  - c121 * source[144] + c121 * source[146] + c122 * source[135]
                  - c122 * source[137] + c122 * source[137] - c122 * source[139];
    target[68] =  c142 * source[250] - c143 * source[241] - c143 * source[243]
                  - c125 * source[115] + c126 * source[106] + c126 * source[108]
                  - c125 * source[145] + c126 * source[136] + c126 * source[138];
    target[69] =  c144 * source[252] - c145 * source[245] - c145 * source[247]
                  - c127 * source[117] + c128 * source[110] + c128 * source[112]
                  - c127 * source[147] + c128 * source[140] + c128 * source[142];
    target[70] =  c144 * source[253] - c145 * source[246] - c145 * source[248]
                  - c127 * source[118] + c128 * source[111] + c128 * source[113]
                  - c127 * source[148] + c128 * source[141] + c128 * source[143];
    target[71] =  c146 * source[254] - c147 * source[249] - c147 * source[251]
                  + c134 * source[240] + c136 * source[242] + c134 * source[244]
                  - c131 * source[119] + c132 * source[114] + c132 * source[116]
                  - c133 * source[105] - c134 * source[107] - c133 * source[109]
                  - c131 * source[149] + c132 * source[144] + c132 * source[146]
                  - c133 * source[135] - c134 * source[137] - c133 * source[139];
    target[72] =  c122 * source[270] - c121 * source[272] + c122 * source[274]
                  - c148 * source[165] + c149 * source[167] - c148 * source[169]
                  - c148 * source[195] + c149 * source[197] - c148 * source[199]
                  + c150 * source[0] - c151 * source[2] + c150 * source[4]
                  + c152 * source[30] - c148 * source[32] + c152 * source[34]
                  + c150 * source[60] - c151 * source[62] + c150 * source[64];
    target[73] =  c143 * source[271] - c143 * source[273] - c121 * source[166]
                  + c121 * source[168] - c121 * source[196] + c121 * source[198]
                  + c124 * source[1] - c124 * source[3] + c122 * source[31]
                  - c122 * source[33] + c124 * source[61] - c124 * source[63];
    target[74] =  c129 * source[275] - c145 * source[277] - c128 * source[170]
                  + c153 * source[172] - c128 * source[200] + c153 * source[202]
                  + c154 * source[5] - c155 * source[7] + c156 * source[35]
                  - c130 * source[37] + c154 * source[65] - c155 * source[67];
    target[75] =  c145 * source[276] - c129 * source[278] - c153 * source[171]
                  + c128 * source[173] - c153 * source[201] + c128 * source[203]
                  + c155 * source[6] - c154 * source[8] + c130 * source[36]
                  - c156 * source[38] + c155 * source[66] - c154 * source[68];
    target[76] =  c157 * source[279] - c157 * source[281] - c158 * source[270]
                  + c158 * source[272] - c158 * source[272] + c158 * source[274]
                  - c159 * source[174] + c159 * source[176] + c160 * source[165]
                  - c160 * source[167] + c160 * source[167] - c160 * source[169]
                  - c159 * source[204] + c159 * source[206] + c160 * source[195]
                  - c160 * source[197] + c160 * source[197] - c160 * source[199]
                  + c161 * source[9] - c161 * source[11] - c162 * source[0]
                  + c162 * source[2] - c162 * source[2] + c162 * source[4]
                  + c160 * source[39] - c160 * source[41] - c163 * source[30]
                  + c163 * source[32] - c163 * source[32] + c163 * source[34]
                  + c161 * source[69] - c161 * source[71] - c162 * source[60]
                  + c162 * source[62] - c162 * source[62] + c162 * source[64];
    target[77] =  c164 * source[280] - c165 * source[271] - c165 * source[273]
                  - c166 * source[175] + c167 * source[166] + c167 * source[168]
                  - c166 * source[205] + c167 * source[196] + c167 * source[198]
                  + c160 * source[10] - c163 * source[1] - c163 * source[3]
                  + c167 * source[40] - c168 * source[31] - c168 * source[33]
                  + c160 * source[70] - c163 * source[61] - c163 * source[63];
    target[78] =  c169 * source[282] - c170 * source[275] - c170 * source[277]
                  - c171 * source[177] + c172 * source[170] + c172 * source[172]
                  - c171 * source[207] + c172 * source[200] + c172 * source[202]
                  + c173 * source[12] - c174 * source[5] - c174 * source[7]
                  + c175 * source[42] - c176 * source[35] - c176 * source[37]
                  + c173 * source[72] - c174 * source[65] - c174 * source[67];
    target[79] =  c169 * source[283] - c170 * source[276] - c170 * source[278]
                  - c171 * source[178] + c172 * source[171] + c172 * source[173]
                  - c171 * source[208] + c172 * source[201] + c172 * source[203]
                  + c173 * source[13] - c174 * source[6] - c174 * source[8]
                  + c175 * source[43] - c176 * source[36] - c176 * source[38]
                  + c173 * source[73] - c174 * source[66] - c174 * source[68];
    target[80] =  c177 * source[284] - c178 * source[279] - c178 * source[281]
                  + c179 * source[270] + c180 * source[272] + c179 * source[274]
                  - c181 * source[179] + c182 * source[174] + c182 * source[176]
                  - c183 * source[165] - c184 * source[167] - c183 * source[169]
                  - c181 * source[209] + c182 * source[204] + c182 * source[206]
                  - c183 * source[195] - c184 * source[197] - c183 * source[199]
                  + c185 * source[14] - c179 * source[9] - c179 * source[11]
                  + c186 * source[0] + c187 * source[2] + c186 * source[4]
                  + c188 * source[44] - c180 * source[39] - c180 * source[41]
                  + c187 * source[30] + c189 * source[32] + c187 * source[34]
                  + c185 * source[74] - c179 * source[69] - c179 * source[71]
                  + c186 * source[60] + c187 * source[62] + c186 * source[64];
    target[81] =  c122 * source[285] - c121 * source[287] + c122 * source[289]
                  - c148 * source[180] + c149 * source[182] - c148 * source[184]
                  - c148 * source[210] + c149 * source[212] - c148 * source[214]
                  + c150 * source[15] - c151 * source[17] + c150 * source[19]
                  + c152 * source[45] - c148 * source[47] + c152 * source[49]
                  + c150 * source[75] - c151 * source[77] + c150 * source[79];
    target[82] =  c143 * source[286] - c143 * source[288] - c121 * source[181]
                  + c121 * source[183] - c121 * source[211] + c121 * source[213]
                  + c124 * source[16] - c124 * source[18] + c122 * source[46]
                  - c122 * source[48] + c124 * source[76] - c124 * source[78];
    target[83] =  c129 * source[290] - c145 * source[292] - c128 * source[185]
                  + c153 * source[187] - c128 * source[215] + c153 * source[217]
                  + c154 * source[20] - c155 * source[22] + c156 * source[50]
                  - c130 * source[52] + c154 * source[80] - c155 * source[82];
    target[84] =  c145 * source[291] - c129 * source[293] - c153 * source[186]
                  + c128 * source[188] - c153 * source[216] + c128 * source[218]
                  + c155 * source[21] - c154 * source[23] + c130 * source[51]
                  - c156 * source[53] + c155 * source[81] - c154 * source[83];
    target[85] =  c157 * source[294] - c157 * source[296] - c158 * source[285]
                  + c158 * source[287] - c158 * source[287] + c158 * source[289]
                  - c159 * source[189] + c159 * source[191] + c160 * source[180]
                  - c160 * source[182] + c160 * source[182] - c160 * source[184]
                  - c159 * source[219] + c159 * source[221] + c160 * source[210]
                  - c160 * source[212] + c160 * source[212] - c160 * source[214]
                  + c161 * source[24] - c161 * source[26] - c162 * source[15]
                  + c162 * source[17] - c162 * source[17] + c162 * source[19]
                  + c160 * source[54] - c160 * source[56] - c163 * source[45]
                  + c163 * source[47] - c163 * source[47] + c163 * source[49]
                  + c161 * source[84] - c161 * source[86] - c162 * source[75]
                  + c162 * source[77] - c162 * source[77] + c162 * source[79];
    target[86] =  c164 * source[295] - c165 * source[286] - c165 * source[288]
                  - c166 * source[190] + c167 * source[181] + c167 * source[183]
                  - c166 * source[220] + c167 * source[211] + c167 * source[213]
                  + c160 * source[25] - c163 * source[16] - c163 * source[18]
                  + c167 * source[55] - c168 * source[46] - c168 * source[48]
                  + c160 * source[85] - c163 * source[76] - c163 * source[78];
    target[87] =  c169 * source[297] - c170 * source[290] - c170 * source[292]
                  - c171 * source[192] + c172 * source[185] + c172 * source[187]
                  - c171 * source[222] + c172 * source[215] + c172 * source[217]
                  + c173 * source[27] - c174 * source[20] - c174 * source[22]
                  + c175 * source[57] - c176 * source[50] - c176 * source[52]
                  + c173 * source[87] - c174 * source[80] - c174 * source[82];
    target[88] =  c169 * source[298] - c170 * source[291] - c170 * source[293]
                  - c171 * source[193] + c172 * source[186] + c172 * source[188]
                  - c171 * source[223] + c172 * source[216] + c172 * source[218]
                  + c173 * source[28] - c174 * source[21] - c174 * source[23]
                  + c175 * source[58] - c176 * source[51] - c176 * source[53]
                  + c173 * source[88] - c174 * source[81] - c174 * source[83];
    target[89] =  c177 * source[299] - c178 * source[294] - c178 * source[296]
                  + c179 * source[285] + c180 * source[287] + c179 * source[289]
                  - c181 * source[194] + c182 * source[189] + c182 * source[191]
                  - c183 * source[180] - c184 * source[182] - c183 * source[184]
                  - c181 * source[224] + c182 * source[219] + c182 * source[221]
                  - c183 * source[210] - c184 * source[212] - c183 * source[214]
                  + c185 * source[29] - c179 * source[24] - c179 * source[26]
                  + c186 * source[15] + c187 * source[17] + c186 * source[19]
                  + c188 * source[59] - c180 * source[54] - c180 * source[56]
                  + c187 * source[45] + c189 * source[47] + c187 * source[49]
                  + c185 * source[89] - c179 * source[84] - c179 * source[86]
                  + c186 * source[75] + c187 * source[77] + c186 * source[79];
    target[90] =  c190 * source[300] - c191 * source[302] + c190 * source[304]
                  - c192 * source[225] + c27 * source[227] - c192 * source[229]
                  - c192 * source[255] + c27 * source[257] - c192 * source[259]
                  + c193 * source[90] - c30 * source[92] + c193 * source[94]
                  + c194 * source[120] - c28 * source[122] + c194 * source[124]
                  + c193 * source[150] - c30 * source[152] + c193 * source[154];
    target[91] =  c195 * source[301] - c195 * source[303] - c196 * source[226]
                  + c196 * source[228] - c196 * source[256] + c196 * source[258]
                  + c197 * source[91] - c197 * source[93] + c29 * source[121]
                  - c29 * source[123] + c197 * source[151] - c197 * source[153];
    target[92] =  c198 * source[305] - c199 * source[307] - c200 * source[230]
                  + c201 * source[232] - c200 * source[260] + c201 * source[262]
                  + c18 * source[95] - c19 * source[97] + c24 * source[125]
                  - c17 * source[127] + c18 * source[155] - c19 * source[157];
    target[93] =  c199 * source[306] - c198 * source[308] - c201 * source[231]
                  + c200 * source[233] - c201 * source[261] + c200 * source[263]
                  + c19 * source[96] - c18 * source[98] + c17 * source[126]
                  - c24 * source[128] + c19 * source[156] - c18 * source[158];
    target[94] =  c202 * source[309] - c202 * source[311] - c203 * source[300]
                  + c203 * source[302] - c203 * source[302] + c203 * source[304]
                  - c204 * source[234] + c204 * source[236] + c205 * source[225]
                  - c205 * source[227] + c205 * source[227] - c205 * source[229]
                  - c204 * source[264] + c204 * source[266] + c205 * source[255]
                  - c205 * source[257] + c205 * source[257] - c205 * source[259]
                  + c206 * source[99] - c206 * source[101] - c207 * source[90]
                  + c207 * source[92] - c207 * source[92] + c207 * source[94]
                  + c208 * source[129] - c208 * source[131] - c209 * source[120]
                  + c209 * source[122] - c209 * source[122] + c209 * source[124]
                  + c206 * source[159] - c206 * source[161] - c207 * source[150]
                  + c207 * source[152] - c207 * source[152] + c207 * source[154];
    target[95] =  c210 * source[310] - c211 * source[301] - c211 * source[303]
                  - c212 * source[235] + c213 * source[226] + c213 * source[228]
                  - c212 * source[265] + c213 * source[256] + c213 * source[258]
                  + c208 * source[100] - c209 * source[91] - c209 * source[93]
                  + c214 * source[130] - c215 * source[121] - c215 * source[123]
                  + c208 * source[160] - c209 * source[151] - c209 * source[153];
    target[96] =  c216 * source[312] - c217 * source[305] - c217 * source[307]
                  - c218 * source[237] + c219 * source[230] + c219 * source[232]
                  - c218 * source[267] + c219 * source[260] + c219 * source[262]
                  + c220 * source[102] - c221 * source[95] - c221 * source[97]
                  + c219 * source[132] - c222 * source[125] - c222 * source[127]
                  + c220 * source[162] - c221 * source[155] - c221 * source[157];
    target[97] =  c216 * source[313] - c217 * source[306] - c217 * source[308]
                  - c218 * source[238] + c219 * source[231] + c219 * source[233]
                  - c218 * source[268] + c219 * source[261] + c219 * source[263]
                  + c220 * source[103] - c221 * source[96] - c221 * source[98]
                  + c219 * source[133] - c222 * source[126] - c222 * source[128]
                  + c220 * source[163] - c221 * source[156] - c221 * source[158];
    target[98] =  source[314] - c223 * source[309] - c223 * source[311]
                  + c224 * source[300] + c225 * source[302] + c224 * source[304]
                  - c226 * source[239] + c227 * source[234] + c227 * source[236]
                  - c228 * source[225] - c229 * source[227] - c228 * source[229]
                  - c226 * source[269] + c227 * source[264] + c227 * source[266]
                  - c228 * source[255] - c229 * source[257] - c228 * source[259]
                  + c228 * source[104] - c230 * source[99] - c230 * source[101]
                  + c231 * source[90] + c232 * source[92] + c231 * source[94]
                  + c229 * source[134] - c233 * source[129] - c233 * source[131]
                  + c232 * source[120] + c234 * source[122] + c232 * source[124]
                  + c228 * source[164] - c230 * source[159] - c230 * source[161]
                  + c231 * source[150] + c232 * source[152] + c231 * source[154];
  }
}

void CCarSphList::carsph_54(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c56 = 89.294106748429925;
  const double c95 = 84.187291202413675;
  const double c50 = 83.526988527660933;
  const double c86 = 78.75;
  const double c142 = 68.738635424337602;
  const double c141 = 64.299105748058423;
  const double c70 = 59.529404498953291;
  const double c44 = 59.0625;
  const double c67 = 55.68465901844062;
  const double c23 = 47.062126492541751;
  const double c138 = 45.46633369868303;
  const double c53 = 44.647053374214963;
  const double c12 = 44.022588307027107;
  const double c57 = 42.093645601206838;
  const double c62 = 39.933538535922409;
  const double c98 = 39.686269665968858;
  const double c46 = 39.375;
  const double c166 = 38.97114317029974;
  const double c105 = 37.649701194033398;
  const double c83 = 37.123106012293746;
  const double c153 = 36.454166428544212;
  const double c125 = 34.369317712168801;
  const double c212 = 33.541019662496844;
  const double c144 = 32.403703492039298;
  const double c118 = 32.149552874029212;
  const double c58 = 31.57023420090513;
  const double c201 = 31.374750995027831;
  const double c3 = 31.128670717282485;
  const double c147 = 30.740852297878796;
  const double c139 = 30.310889132455351;
  const double c68 = 29.764702249476645;
  const double c72 = 28.062430400804562;
  const double c49 = 27.84232950922031;
  const double c74 = 26.622359023948274;
  const double c65 = 26.25;
  const double c164 = 25.98076211353316;
  const double c149 = 25.776988284126599;
  const double c214 = 25.155764746872634;
  const double c145 = 24.302777619029477;
  const double c17 = 23.531063246270875;
  const double c112 = 22.733166849341515;
  const double c27 = 22.18529918662356;
  const double c14 = 22.011294153513553;
  const double c140 = 21.433035249352809;
  const double c36 = 21.046822800603419;
  const double c7 = 20.75244714485499;
  const double c159 = 19.48557158514987;
  const double c66 = 18.561553006146873;
  const double c171 = 18.371173070873837;
  const double c182 = 17.428425057933374;
  const double c121 = 17.1846588560844;
  const double c204 = 16.770509831248422;
  const double c28 = 16.638974389967672;
  const double c127 = 16.201851746019649;
  const double c120 = 16.074776437014606;
  const double c218 = 15.811388300841896;
  const double c5 = 15.564335358641243;
  const double c132 = 15.370426148939398;
  const double c115 = 15.155444566227676;
  const double c227 = 15;
  const double c54 = 14.882351124738323;
  const double c196 = 14.79019945774904;
  const double c11 = 14.674196102342369;
  const double c89 = 14.031215200402281;
  const double c48 = 13.921164754610155;
  const double c172 = 13.778379803155376;
  const double c61 = 13.311179511974137;
  const double c97 = 13.228756555322953;
  const double c157 = 12.99038105676658;
  const double c208 = 12.577882373436317;
  const double c102 = 12.549900398011133;
  const double c82 = 12.374368670764582;
  const double c169 = 12.24744871391589;
  const double c128 = 12.151388809514739;
  const double c219 = 11.858541225631422;
  const double c19 = 11.765531623135438;
  const double c178 = 11.61895003862225;
  const double c143 = 11.456439237389599;
  const double c114 = 11.366583424670758;
  const double c233 = 11.25;
  const double c29 = 11.09264959331178;
  const double c117 = 10.716517624676404;
  const double c40 = 10.523411400301709;
  const double c200 = 10.458250331675945;
  const double c8 = 10.376223572427495;
  const double c146 = 10.246950765959598;
  const double c64 = 9.9833846339806023;
  const double c71 = 9.9215674164922145;
  const double c43 = 9.84375;
  const double c106 = 9.4124252985083494;
  const double c77 = 9.2807765030734366;
  const double c170 = 9.1855865354369186;
  const double c222 = 8.8939059192235668;
  const double c73 = 8.8741196746494246;
  const double c85 = 8.75;
  const double c123 = 8.5923294280422002;
  const double c30 = 8.3194871949838358;
  const double c129 = 8.1009258730098246;
  const double c24 = 7.8436877487569578;
  const double c136 = 7.6852130744696989;
  const double c116 = 7.5777222831138378;
  const double c51 = 7.4411755623691613;
  const double c13 = 7.3370980511711847;
  const double c35 = 7.0156076002011405;
  const double c81 = 6.9605823773050775;
  const double c210 = 6.7082039324993694;
  const double c59 = 6.6555897559870685;
  const double c45 = 6.5625;
  const double c167 = 6.49519052838329;
  const double c206 = 6.2889411867181586;
  const double c199 = 6.2749501990055663;
  const double c130 = 6.0756944047573693;
  const double c220 = 5.9292706128157109;
  const double c181 = 5.8094750193111251;
  const double c126 = 5.7282196186947996;
  const double c230 = 5.625;
  const double c213 = 5.5901699437494745;
  const double c197 = 5.54632479665589;
  const double c119 = 5.3582588123382022;
  const double c38 = 5.2617057001508547;
  const double c2 = 5.1881117862137476;
  const double c131 = 5.123475382979799;
  const double c226 = 5;
  const double c63 = 4.9916923169903011;
  const double c69 = 4.9607837082461073;
  const double c21 = 4.7062126492541747;
  const double c94 = 4.677071733467427;
  const double c47 = 4.6403882515367183;
  const double c221 = 4.4469529596117834;
  const double c191 = 4.4370598373247123;
  const double c10 = 4.4022588307027108;
  const double c184 = 4.3571062644833436;
  const double c165 = 4.3301270189221936;
  const double c148 = 4.2961647140211001;
  const double c215 = 4.1926274578121054;
  const double c101 = 4.1833001326703778;
  const double c18 = 3.9218438743784789;
  const double c177 = 3.872983346207417;
  const double c134 = 3.8426065372348495;
  const double c111 = 3.7888611415569189;
  const double c229 = 3.75;
  const double c100 = 3.7205877811845807;
  const double c192 = 3.69754986443726;
  const double c39 = 3.5078038001005702;
  const double c202 = 3.3541019662496847;
  const double c75 = 3.3277948779935342;
  const double c88 = 3.28125;
  const double c160 = 3.247595264191645;
  const double c216 = 3.1622776601683795;
  const double c104 = 3.1374750995027831;
  const double c1 = 3.1128670717282483;
  const double c76 = 3.0935921676911455;
  const double c175 = 3.0618621784789726;
  const double c155 = 3.0378472023786847;
  const double c223 = 3;
  const double c195 = 2.9580398915498081;
  const double c180 = 2.9047375096555625;
  const double c122 = 2.8641098093473998;
  const double c234 = 2.8125;
  const double c205 = 2.7950849718747373;
  const double c194 = 2.773162398327945;
  const double c37 = 2.6308528500754274;
  const double c4 = 2.5940558931068738;
  const double c135 = 2.5617376914898995;
  const double c55 = 2.4803918541230536;
  const double c217 = 2.3717082451262845;
  const double c15 = 2.3531063246270874;
  const double c90 = 2.3385358667337135;
  const double c79 = 2.3201941257683592;
  const double c176 = 2.2963966338592297;
  const double c25 = 2.2185299186623562;
  const double c183 = 2.1785531322416718;
  const double c158 = 2.1650635094610968;
  const double c151 = 2.1480823570105501;
  const double c32 = 2.104682280060342;
  const double c209 = 2.0963137289060527;
  const double c198 = 2.0916500663351889;
  const double c6 = 2.0752447144854989;
  const double c156 = 2.0252314682524561;
  const double c20 = 1.9609219371892395;
  const double c133 = 1.9213032686174247;
  const double c113 = 1.8944305707784594;
  const double c228 = 1.875;
  const double c91 = 1.7539019000502851;
  const double c26 = 1.6638974389967671;
  const double c99 = 1.6535945694153691;
  const double c42 = 1.640625;
  const double c161 = 1.6237976320958225;
  const double c103 = 1.5687375497513916;
  const double c84 = 1.5467960838455728;
  const double c173 = 1.5309310892394863;
  const double c9 = 1.4674196102342369;
  const double c179 = 1.4523687548277813;
  const double c124 = 1.4320549046736999;
  const double c232 = 1.40625;
  const double c193 = 1.3865811991639725;
  const double c41 = 1.3154264250377137;
  const double c52 = 1.2401959270615268;
  const double c110 = 1.1765531623135437;
  const double c80 = 1.1600970628841796;
  const double c174 = 1.1481983169296148;
  const double c211 = 1.1180339887498949;
  const double c87 = 1.09375;
  const double c168 = 1.0825317547305484;
  const double c207 = 1.0481568644530264;
  const double c154 = 1.0126157341262281;
  const double c188 = 0.96824583655185426;
  const double c137 = 0.96065163430871237;
  const double c93 = 0.87695095002514256;
  const double c60 = 0.83194871949838356;
  const double c22 = 0.78436877487569578;
  const double c225 = 0.75;
  const double c190 = 0.73950997288745202;
  const double c189 = 0.72618437741389064;
  const double c152 = 0.71602745233684995;
  const double c231 = 0.703125;
  const double c31 = 0.70156076002011403;
  const double c109 = 0.58827658115677184;
  const double c96 = 0.58463396668342837;
  const double c203 = 0.55901699437494745;
  const double c163 = 0.54126587736527421;
  const double c34 = 0.52617057001508549;
  const double c107 = 0.52291251658379723;
  const double c0 = 0.51881117862137471;
  const double c185 = 0.48412291827592713;
  const double c16 = 0.39218438743784789;
  const double c78 = 0.38669902096139319;
  const double c224 = 0.375;
  const double c187 = 0.36309218870694532;
  const double c150 = 0.35801372616842497;
  const double c92 = 0.29231698334171419;
  const double c162 = 0.2706329386826371;
  const double c33 = 0.26308528500754275;
  const double c108 = 0.19609219371892395;
  const double c186 = 0.18154609435347266;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 99, source += 315) {
    target[0] =  c0 * source[0] - c1 * source[2] + c0 * source[4]
                  - c2 * source[30] + c3 * source[32] - c2 * source[34]
                  + c4 * source[60] - c5 * source[62] + c4 * source[64];
    target[1] =  c6 * source[1] - c6 * source[3] - c7 * source[31]
                  + c7 * source[33] + c8 * source[61] - c8 * source[63];
    target[2] =  c9 * source[5] - c10 * source[7] - c11 * source[35]
                  + c12 * source[37] + c13 * source[65] - c14 * source[67];
    target[3] =  c10 * source[6] - c9 * source[8] - c12 * source[36]
                  + c11 * source[38] + c14 * source[66] - c13 * source[68];
    target[4] =  c15 * source[9] - c15 * source[11] - c16 * source[0]
                  + c16 * source[2] - c16 * source[2] + c16 * source[4]
                  - c17 * source[39] + c17 * source[41] + c18 * source[30]
                  - c18 * source[32] + c18 * source[32] - c18 * source[34]
                  + c19 * source[69] - c19 * source[71] - c20 * source[60]
                  + c20 * source[62] - c20 * source[62] + c20 * source[64];
    target[5] =  c21 * source[10] - c22 * source[1] - c22 * source[3]
                  - c23 * source[40] + c24 * source[31] + c24 * source[33]
                  + c17 * source[70] - c18 * source[61] - c18 * source[63];
    target[6] =  c25 * source[12] - c26 * source[5] - c26 * source[7]
                  - c27 * source[42] + c28 * source[35] + c28 * source[37]
                  + c29 * source[72] - c30 * source[65] - c30 * source[67];
    target[7] =  c25 * source[13] - c26 * source[6] - c26 * source[8]
                  - c27 * source[43] + c28 * source[36] + c28 * source[38]
                  + c29 * source[73] - c30 * source[66] - c30 * source[68];
    target[8] =  c31 * source[14] - c32 * source[9] - c32 * source[11]
                  + c33 * source[0] + c34 * source[2] + c33 * source[4]
                  - c35 * source[44] + c36 * source[39] + c36 * source[41]
                  - c37 * source[30] - c38 * source[32] - c37 * source[34]
                  + c39 * source[74] - c40 * source[69] - c40 * source[71]
                  + c41 * source[60] + c37 * source[62] + c41 * source[64];
    target[9] =  c4 * source[15] - c5 * source[17] + c4 * source[19]
                  - c2 * source[45] + c3 * source[47] - c2 * source[49]
                  + c0 * source[75] - c1 * source[77] + c0 * source[79];
    target[10] =  c8 * source[16] - c8 * source[18] - c7 * source[46]
                  + c7 * source[48] + c6 * source[76] - c6 * source[78];
    target[11] =  c13 * source[20] - c14 * source[22] - c11 * source[50]
                  + c12 * source[52] + c9 * source[80] - c10 * source[82];
    target[12] =  c14 * source[21] - c13 * source[23] - c12 * source[51]
                  + c11 * source[53] + c10 * source[81] - c9 * source[83];
    target[13] =  c19 * source[24] - c19 * source[26] - c20 * source[15]
                  + c20 * source[17] - c20 * source[17] + c20 * source[19]
                  - c17 * source[54] + c17 * source[56] + c18 * source[45]
                  - c18 * source[47] + c18 * source[47] - c18 * source[49]
                  + c15 * source[84] - c15 * source[86] - c16 * source[75]
                  + c16 * source[77] - c16 * source[77] + c16 * source[79];
    target[14] =  c17 * source[25] - c18 * source[16] - c18 * source[18]
                  - c23 * source[55] + c24 * source[46] + c24 * source[48]
                  + c21 * source[85] - c22 * source[76] - c22 * source[78];
    target[15] =  c29 * source[27] - c30 * source[20] - c30 * source[22]
                  - c27 * source[57] + c28 * source[50] + c28 * source[52]
                  + c25 * source[87] - c26 * source[80] - c26 * source[82];
    target[16] =  c29 * source[28] - c30 * source[21] - c30 * source[23]
                  - c27 * source[58] + c28 * source[51] + c28 * source[53]
                  + c25 * source[88] - c26 * source[81] - c26 * source[83];
    target[17] =  c39 * source[29] - c40 * source[24] - c40 * source[26]
                  + c41 * source[15] + c37 * source[17] + c41 * source[19]
                  - c35 * source[59] + c36 * source[54] + c36 * source[56]
                  - c37 * source[45] - c38 * source[47] - c37 * source[49]
                  + c31 * source[89] - c32 * source[84] - c32 * source[86]
                  + c33 * source[75] + c34 * source[77] + c33 * source[79];
    target[18] =  c42 * source[90] - c43 * source[92] + c42 * source[94]
                  - c43 * source[120] + c44 * source[122] - c43 * source[124]
                  + c42 * source[150] - c43 * source[152] + c42 * source[154];
    target[19] =  c45 * source[91] - c45 * source[93] - c46 * source[121]
                  + c46 * source[123] + c45 * source[151] - c45 * source[153];
    target[20] =  c47 * source[95] - c48 * source[97] - c49 * source[125]
                  + c50 * source[127] + c47 * source[155] - c48 * source[157];
    target[21] =  c48 * source[96] - c47 * source[98] - c50 * source[126]
                  + c49 * source[128] + c48 * source[156] - c47 * source[158];
    target[22] =  c51 * source[99] - c51 * source[101] - c52 * source[90]
                  + c52 * source[92] - c52 * source[92] + c52 * source[94]
                  - c53 * source[129] + c53 * source[131] + c51 * source[120]
                  - c51 * source[122] + c51 * source[122] - c51 * source[124]
                  + c51 * source[159] - c51 * source[161] - c52 * source[150]
                  + c52 * source[152] - c52 * source[152] + c52 * source[154];
    target[23] =  c54 * source[100] - c55 * source[91] - c55 * source[93]
                  - c56 * source[130] + c54 * source[121] + c54 * source[123]
                  + c54 * source[160] - c55 * source[151] - c55 * source[153];
    target[24] =  c35 * source[102] - c38 * source[95] - c38 * source[97]
                  - c57 * source[132] + c58 * source[125] + c58 * source[127]
                  + c35 * source[162] - c38 * source[155] - c38 * source[157];
    target[25] =  c35 * source[103] - c38 * source[96] - c38 * source[98]
                  - c57 * source[133] + c58 * source[126] + c58 * source[128]
                  + c35 * source[163] - c38 * source[156] - c38 * source[158];
    target[26] =  c25 * source[104] - c59 * source[99] - c59 * source[101]
                  + c60 * source[90] + c26 * source[92] + c60 * source[94]
                  - c61 * source[134] + c62 * source[129] + c62 * source[131]
                  - c63 * source[120] - c64 * source[122] - c63 * source[124]
                  + c25 * source[164] - c59 * source[159] - c59 * source[161]
                  + c60 * source[150] + c26 * source[152] + c60 * source[154];
    target[27] =  c45 * source[105] - c46 * source[107] + c45 * source[109]
                  - c45 * source[135] + c46 * source[137] - c45 * source[139];
    target[28] =  c65 * source[106] - c65 * source[108] - c65 * source[136]
                  + c65 * source[138];
    target[29] =  c66 * source[110] - c67 * source[112] - c66 * source[140]
                  + c67 * source[142];
    target[30] =  c67 * source[111] - c66 * source[113] - c67 * source[141]
                  + c66 * source[143];
    target[31] =  c68 * source[114] - c68 * source[116] - c69 * source[105]
                  + c69 * source[107] - c69 * source[107] + c69 * source[109]
                  - c68 * source[144] + c68 * source[146] + c69 * source[135]
                  - c69 * source[137] + c69 * source[137] - c69 * source[139];
    target[32] =  c70 * source[115] - c71 * source[106] - c71 * source[108]
                  - c70 * source[145] + c71 * source[136] + c71 * source[138];
    target[33] =  c72 * source[117] - c36 * source[110] - c36 * source[112]
                  - c72 * source[147] + c36 * source[140] + c36 * source[142];
    target[34] =  c72 * source[118] - c36 * source[111] - c36 * source[113]
                  - c72 * source[148] + c36 * source[141] + c36 * source[143];
    target[35] =  c73 * source[119] - c74 * source[114] - c74 * source[116]
                  + c75 * source[105] + c59 * source[107] + c75 * source[109]
                  - c73 * source[149] + c74 * source[144] + c74 * source[146]
                  - c75 * source[135] - c59 * source[137] - c75 * source[139];
    target[36] =  c76 * source[165] - c66 * source[167] + c76 * source[169]
                  - c77 * source[195] + c67 * source[197] - c77 * source[199]
                  - c78 * source[0] + c79 * source[2] - c78 * source[4]
                  + c80 * source[30] - c81 * source[32] + c80 * source[34]
                  - c78 * source[30] + c79 * source[32] - c78 * source[34]
                  + c80 * source[60] - c81 * source[62] + c80 * source[64];
    target[37] =  c82 * source[166] - c82 * source[168] - c83 * source[196]
                  + c83 * source[198] - c84 * source[1] + c84 * source[3]
                  + c47 * source[31] - c47 * source[33] - c84 * source[31]
                  + c84 * source[33] + c47 * source[61] - c47 * source[63];
    target[38] =  c85 * source[170] - c65 * source[172] - c65 * source[200]
                  + c86 * source[202] - c87 * source[5] + c88 * source[7]
                  + c88 * source[35] - c43 * source[37] - c87 * source[35]
                  + c88 * source[37] + c88 * source[65] - c43 * source[67];
    target[39] =  c65 * source[171] - c85 * source[173] - c86 * source[201]
                  + c65 * source[203] - c88 * source[6] + c87 * source[8]
                  + c43 * source[36] - c88 * source[38] - c88 * source[36]
                  + c87 * source[38] + c43 * source[66] - c88 * source[68];
    target[40] =  c89 * source[174] - c89 * source[176] - c90 * source[165]
                  + c90 * source[167] - c90 * source[167] + c90 * source[169]
                  - c57 * source[204] + c57 * source[206] + c35 * source[195]
                  - c35 * source[197] + c35 * source[197] - c35 * source[199]
                  - c91 * source[9] + c91 * source[11] + c92 * source[0]
                  - c92 * source[2] + c92 * source[2] - c92 * source[4]
                  + c38 * source[39] - c38 * source[41] - c93 * source[30]
                  + c93 * source[32] - c93 * source[32] + c93 * source[34]
                  - c91 * source[39] + c91 * source[41] + c92 * source[30]
                  - c92 * source[32] + c92 * source[32] - c92 * source[34]
                  + c38 * source[69] - c38 * source[71] - c93 * source[60]
                  + c93 * source[62] - c93 * source[62] + c93 * source[64];
    target[41] =  c72 * source[175] - c94 * source[166] - c94 * source[168]
                  - c95 * source[205] + c89 * source[196] + c89 * source[198]
                  - c39 * source[10] + c96 * source[1] + c96 * source[3]
                  + c40 * source[40] - c91 * source[31] - c91 * source[33]
                  - c39 * source[40] + c96 * source[31] + c96 * source[33]
                  + c40 * source[70] - c91 * source[61] - c91 * source[63];
    target[42] =  c97 * source[177] - c71 * source[170] - c71 * source[172]
                  - c98 * source[207] + c68 * source[200] + c68 * source[202]
                  - c99 * source[12] + c52 * source[5] + c52 * source[7]
                  + c69 * source[42] - c100 * source[35] - c100 * source[37]
                  - c99 * source[42] + c52 * source[35] + c52 * source[37]
                  + c69 * source[72] - c100 * source[65] - c100 * source[67];
    target[43] =  c97 * source[178] - c71 * source[171] - c71 * source[173]
                  - c98 * source[208] + c68 * source[201] + c68 * source[203]
                  - c99 * source[13] + c52 * source[6] + c52 * source[8]
                  + c69 * source[43] - c100 * source[36] - c100 * source[38]
                  - c99 * source[43] + c52 * source[36] + c52 * source[38]
                  + c69 * source[73] - c100 * source[66] - c100 * source[68];
    target[44] =  c101 * source[179] - c102 * source[174] - c102 * source[176]
                  + c103 * source[165] + c104 * source[167] + c103 * source[169]
                  - c102 * source[209] + c105 * source[204] + c105 * source[206]
                  - c21 * source[195] - c106 * source[197] - c21 * source[199]
                  - c107 * source[14] + c103 * source[9] + c103 * source[11]
                  - c108 * source[0] - c16 * source[2] - c108 * source[4]
                  + c103 * source[44] - c21 * source[39] - c21 * source[41]
                  + c109 * source[30] + c110 * source[32] + c109 * source[34]
                  - c107 * source[44] + c103 * source[39] + c103 * source[41]
                  - c108 * source[30] - c16 * source[32] - c108 * source[34]
                  + c103 * source[74] - c21 * source[69] - c21 * source[71]
                  + c109 * source[60] + c110 * source[62] + c109 * source[64];
    target[45] =  c77 * source[180] - c67 * source[182] + c77 * source[184]
                  - c76 * source[210] + c66 * source[212] - c76 * source[214]
                  - c80 * source[15] + c81 * source[17] - c80 * source[19]
                  + c78 * source[45] - c79 * source[47] + c78 * source[49]
                  - c80 * source[45] + c81 * source[47] - c80 * source[49]
                  + c78 * source[75] - c79 * source[77] + c78 * source[79];
    target[46] =  c83 * source[181] - c83 * source[183] - c82 * source[211]
                  + c82 * source[213] - c47 * source[16] + c47 * source[18]
                  + c84 * source[46] - c84 * source[48] - c47 * source[46]
                  + c47 * source[48] + c84 * source[76] - c84 * source[78];
    target[47] =  c65 * source[185] - c86 * source[187] - c85 * source[215]
                  + c65 * source[217] - c88 * source[20] + c43 * source[22]
                  + c87 * source[50] - c88 * source[52] - c88 * source[50]
                  + c43 * source[52] + c87 * source[80] - c88 * source[82];
    target[48] =  c86 * source[186] - c65 * source[188] - c65 * source[216]
                  + c85 * source[218] - c43 * source[21] + c88 * source[23]
                  + c88 * source[51] - c87 * source[53] - c43 * source[51]
                  + c88 * source[53] + c88 * source[81] - c87 * source[83];
    target[49] =  c57 * source[189] - c57 * source[191] - c35 * source[180]
                  + c35 * source[182] - c35 * source[182] + c35 * source[184]
                  - c89 * source[219] + c89 * source[221] + c90 * source[210]
                  - c90 * source[212] + c90 * source[212] - c90 * source[214]
                  - c38 * source[24] + c38 * source[26] + c93 * source[15]
                  - c93 * source[17] + c93 * source[17] - c93 * source[19]
                  + c91 * source[54] - c91 * source[56] - c92 * source[45]
                  + c92 * source[47] - c92 * source[47] + c92 * source[49]
                  - c38 * source[54] + c38 * source[56] + c93 * source[45]
                  - c93 * source[47] + c93 * source[47] - c93 * source[49]
                  + c91 * source[84] - c91 * source[86] - c92 * source[75]
                  + c92 * source[77] - c92 * source[77] + c92 * source[79];
    target[50] =  c95 * source[190] - c89 * source[181] - c89 * source[183]
                  - c72 * source[220] + c94 * source[211] + c94 * source[213]
                  - c40 * source[25] + c91 * source[16] + c91 * source[18]
                  + c39 * source[55] - c96 * source[46] - c96 * source[48]
                  - c40 * source[55] + c91 * source[46] + c91 * source[48]
                  + c39 * source[85] - c96 * source[76] - c96 * source[78];
    target[51] =  c98 * source[192] - c68 * source[185] - c68 * source[187]
                  - c97 * source[222] + c71 * source[215] + c71 * source[217]
                  - c69 * source[27] + c100 * source[20] + c100 * source[22]
                  + c99 * source[57] - c52 * source[50] - c52 * source[52]
                  - c69 * source[57] + c100 * source[50] + c100 * source[52]
                  + c99 * source[87] - c52 * source[80] - c52 * source[82];
    target[52] =  c98 * source[193] - c68 * source[186] - c68 * source[188]
                  - c97 * source[223] + c71 * source[216] + c71 * source[218]
                  - c69 * source[28] + c100 * source[21] + c100 * source[23]
                  + c99 * source[58] - c52 * source[51] - c52 * source[53]
                  - c69 * source[58] + c100 * source[51] + c100 * source[53]
                  + c99 * source[88] - c52 * source[81] - c52 * source[83];
    target[53] =  c102 * source[194] - c105 * source[189] - c105 * source[191]
                  + c21 * source[180] + c106 * source[182] + c21 * source[184]
                  - c101 * source[224] + c102 * source[219] + c102 * source[221]
                  - c103 * source[210] - c104 * source[212] - c103 * source[214]
                  - c103 * source[29] + c21 * source[24] + c21 * source[26]
                  - c109 * source[15] - c110 * source[17] - c109 * source[19]
                  + c107 * source[59] - c103 * source[54] - c103 * source[56]
                  + c108 * source[45] + c16 * source[47] + c108 * source[49]
                  - c103 * source[59] + c21 * source[54] + c21 * source[56]
                  - c109 * source[45] - c110 * source[47] - c109 * source[49]
                  + c107 * source[89] - c103 * source[84] - c103 * source[86]
                  + c108 * source[75] + c16 * source[77] + c108 * source[79];
    target[54] =  c111 * source[225] - c112 * source[227] + c111 * source[229]
                  - c111 * source[255] + c112 * source[257] - c111 * source[259]
                  - c113 * source[90] + c114 * source[92] - c113 * source[94]
                  + c113 * source[120] - c114 * source[122] + c113 * source[124]
                  - c113 * source[120] + c114 * source[122] - c113 * source[124]
                  + c113 * source[150] - c114 * source[152] + c113 * source[154];
    target[55] =  c115 * source[226] - c115 * source[228] - c115 * source[256]
                  + c115 * source[258] - c116 * source[91] + c116 * source[93]
                  + c116 * source[121] - c116 * source[123] - c116 * source[121]
                  + c116 * source[123] + c116 * source[151] - c116 * source[153];
    target[56] =  c117 * source[230] - c118 * source[232] - c117 * source[260]
                  + c118 * source[262] - c119 * source[95] + c120 * source[97]
                  + c119 * source[125] - c120 * source[127] - c119 * source[125]
                  + c120 * source[127] + c119 * source[155] - c120 * source[157];
    target[57] =  c118 * source[231] - c117 * source[233] - c118 * source[261]
                  + c117 * source[263] - c120 * source[96] + c119 * source[98]
                  + c120 * source[126] - c119 * source[128] - c120 * source[126]
                  + c119 * source[128] + c120 * source[156] - c119 * source[158];
    target[58] =  c121 * source[234] - c121 * source[236] - c122 * source[225]
                  + c122 * source[227] - c122 * source[227] + c122 * source[229]
                  - c121 * source[264] + c121 * source[266] + c122 * source[255]
                  - c122 * source[257] + c122 * source[257] - c122 * source[259]
                  - c123 * source[99] + c123 * source[101] + c124 * source[90]
                  - c124 * source[92] + c124 * source[92] - c124 * source[94]
                  + c123 * source[129] - c123 * source[131] - c124 * source[120]
                  + c124 * source[122] - c124 * source[122] + c124 * source[124]
                  - c123 * source[129] + c123 * source[131] + c124 * source[120]
                  - c124 * source[122] + c124 * source[122] - c124 * source[124]
                  + c123 * source[159] - c123 * source[161] - c124 * source[150]
                  + c124 * source[152] - c124 * source[152] + c124 * source[154];
    target[59] =  c125 * source[235] - c126 * source[226] - c126 * source[228]
                  - c125 * source[265] + c126 * source[256] + c126 * source[258]
                  - c121 * source[100] + c122 * source[91] + c122 * source[93]
                  + c121 * source[130] - c122 * source[121] - c122 * source[123]
                  - c121 * source[130] + c122 * source[121] + c122 * source[123]
                  + c121 * source[160] - c122 * source[151] - c122 * source[153];
    target[60] =  c127 * source[237] - c128 * source[230] - c128 * source[232]
                  - c127 * source[267] + c128 * source[260] + c128 * source[262]
                  - c129 * source[102] + c130 * source[95] + c130 * source[97]
                  + c129 * source[132] - c130 * source[125] - c130 * source[127]
                  - c129 * source[132] + c130 * source[125] + c130 * source[127]
                  + c129 * source[162] - c130 * source[155] - c130 * source[157];
    target[61] =  c127 * source[238] - c128 * source[231] - c128 * source[233]
                  - c127 * source[268] + c128 * source[261] + c128 * source[263]
                  - c129 * source[103] + c130 * source[96] + c130 * source[98]
                  + c129 * source[133] - c130 * source[126] - c130 * source[128]
                  - c129 * source[133] + c130 * source[126] + c130 * source[128]
                  + c129 * source[163] - c130 * source[156] - c130 * source[158];
    target[62] =  c131 * source[239] - c132 * source[234] - c132 * source[236]
                  + c133 * source[225] + c134 * source[227] + c133 * source[229]
                  - c131 * source[269] + c132 * source[264] + c132 * source[266]
                  - c133 * source[255] - c134 * source[257] - c133 * source[259]
                  - c135 * source[104] + c136 * source[99] + c136 * source[101]
                  - c137 * source[90] - c133 * source[92] - c137 * source[94]
                  + c135 * source[134] - c136 * source[129] - c136 * source[131]
                  + c137 * source[120] + c133 * source[122] + c137 * source[124]
                  - c135 * source[134] + c136 * source[129] + c136 * source[131]
                  - c137 * source[120] - c133 * source[122] - c137 * source[124]
                  + c135 * source[164] - c136 * source[159] - c136 * source[161]
                  + c137 * source[150] + c133 * source[152] + c137 * source[154];
    target[63] =  c116 * source[240] - c138 * source[242] + c116 * source[244]
                  - c111 * source[105] + c112 * source[107] - c111 * source[109]
                  - c111 * source[135] + c112 * source[137] - c111 * source[139];
    target[64] =  c139 * source[241] - c139 * source[243] - c115 * source[106]
                  + c115 * source[108] - c115 * source[136] + c115 * source[138];
    target[65] =  c140 * source[245] - c141 * source[247] - c117 * source[110]
                  + c118 * source[112] - c117 * source[140] + c118 * source[142];
    target[66] =  c141 * source[246] - c140 * source[248] - c118 * source[111]
                  + c117 * source[113] - c118 * source[141] + c117 * source[143];
    target[67] =  c125 * source[249] - c125 * source[251] - c126 * source[240]
                  + c126 * source[242] - c126 * source[242] + c126 * source[244]
                  - c121 * source[114] + c121 * source[116] + c122 * source[105]
                  - c122 * source[107] + c122 * source[107] - c122 * source[109]
                  - c121 * source[144] + c121 * source[146] + c122 * source[135]
                  - c122 * source[137] + c122 * source[137] - c122 * source[139];
    target[68] =  c142 * source[250] - c143 * source[241] - c143 * source[243]
                  - c125 * source[115] + c126 * source[106] + c126 * source[108]
                  - c125 * source[145] + c126 * source[136] + c126 * source[138];
    target[69] =  c144 * source[252] - c145 * source[245] - c145 * source[247]
                  - c127 * source[117] + c128 * source[110] + c128 * source[112]
                  - c127 * source[147] + c128 * source[140] + c128 * source[142];
    target[70] =  c144 * source[253] - c145 * source[246] - c145 * source[248]
                  - c127 * source[118] + c128 * source[111] + c128 * source[113]
                  - c127 * source[148] + c128 * source[141] + c128 * source[143];
    target[71] =  c146 * source[254] - c147 * source[249] - c147 * source[251]
                  + c134 * source[240] + c136 * source[242] + c134 * source[244]
                  - c131 * source[119] + c132 * source[114] + c132 * source[116]
                  - c133 * source[105] - c134 * source[107] - c133 * source[109]
                  - c131 * source[149] + c132 * source[144] + c132 * source[146]
                  - c133 * source[135] - c134 * source[137] - c133 * source[139];
    target[72] =  c122 * source[270] - c121 * source[272] + c122 * source[274]
                  - c148 * source[165] + c149 * source[167] - c148 * source[169]
                  - c148 * source[195] + c149 * source[197] - c148 * source[199]
                  + c150 * source[0] - c151 * source[2] + c150 * source[4]
                  + c152 * source[30] - c148 * source[32] + c152 * source[34]
                  + c150 * source[60] - c151 * source[62] + c150 * source[64];
    target[73] =  c143 * source[271] - c143 * source[273] - c121 * source[166]
                  + c121 * source[168] - c121 * source[196] + c121 * source[198]
                  + c124 * source[1] - c124 * source[3] + c122 * source[31]
                  - c122 * source[33] + c124 * source[61] - c124 * source[63];
    target[74] =  c129 * source[275] - c145 * source[277] - c128 * source[170]
                  + c153 * source[172] - c128 * source[200] + c153 * source[202]
                  + c154 * source[5] - c155 * source[7] + c156 * source[35]
                  - c130 * source[37] + c154 * source[65] - c155 * source[67];
    target[75] =  c145 * source[276] - c129 * source[278] - c153 * source[171]
                  + c128 * source[173] - c153 * source[201] + c128 * source[203]
                  + c155 * source[6] - c154 * source[8] + c130 * source[36]
                  - c156 * source[38] + c155 * source[66] - c154 * source[68];
    target[76] =  c157 * source[279] - c157 * source[281] - c158 * source[270]
                  + c158 * source[272] - c158 * source[272] + c158 * source[274]
                  - c159 * source[174] + c159 * source[176] + c160 * source[165]
                  - c160 * source[167] + c160 * source[167] - c160 * source[169]
                  - c159 * source[204] + c159 * source[206] + c160 * source[195]
                  - c160 * source[197] + c160 * source[197] - c160 * source[199]
                  + c161 * source[9] - c161 * source[11] - c162 * source[0]
                  + c162 * source[2] - c162 * source[2] + c162 * source[4]
                  + c160 * source[39] - c160 * source[41] - c163 * source[30]
                  + c163 * source[32] - c163 * source[32] + c163 * source[34]
                  + c161 * source[69] - c161 * source[71] - c162 * source[60]
                  + c162 * source[62] - c162 * source[62] + c162 * source[64];
    target[77] =  c164 * source[280] - c165 * source[271] - c165 * source[273]
                  - c166 * source[175] + c167 * source[166] + c167 * source[168]
                  - c166 * source[205] + c167 * source[196] + c167 * source[198]
                  + c160 * source[10] - c163 * source[1] - c163 * source[3]
                  + c167 * source[40] - c168 * source[31] - c168 * source[33]
                  + c160 * source[70] - c163 * source[61] - c163 * source[63];
    target[78] =  c169 * source[282] - c170 * source[275] - c170 * source[277]
                  - c171 * source[177] + c172 * source[170] + c172 * source[172]
                  - c171 * source[207] + c172 * source[200] + c172 * source[202]
                  + c173 * source[12] - c174 * source[5] - c174 * source[7]
                  + c175 * source[42] - c176 * source[35] - c176 * source[37]
                  + c173 * source[72] - c174 * source[65] - c174 * source[67];
    target[79] =  c169 * source[283] - c170 * source[276] - c170 * source[278]
                  - c171 * source[178] + c172 * source[171] + c172 * source[173]
                  - c171 * source[208] + c172 * source[201] + c172 * source[203]
                  + c173 * source[13] - c174 * source[6] - c174 * source[8]
                  + c175 * source[43] - c176 * source[36] - c176 * source[38]
                  + c173 * source[73] - c174 * source[66] - c174 * source[68];
    target[80] =  c177 * source[284] - c178 * source[279] - c178 * source[281]
                  + c179 * source[270] + c180 * source[272] + c179 * source[274]
                  - c181 * source[179] + c182 * source[174] + c182 * source[176]
                  - c183 * source[165] - c184 * source[167] - c183 * source[169]
                  - c181 * source[209] + c182 * source[204] + c182 * source[206]
                  - c183 * source[195] - c184 * source[197] - c183 * source[199]
                  + c185 * source[14] - c179 * source[9] - c179 * source[11]
                  + c186 * source[0] + c187 * source[2] + c186 * source[4]
                  + c188 * source[44] - c180 * source[39] - c180 * source[41]
                  + c187 * source[30] + c189 * source[32] + c187 * source[34]
                  + c185 * source[74] - c179 * source[69] - c179 * source[71]
                  + c186 * source[60] + c187 * source[62] + c186 * source[64];
    target[81] =  c122 * source[285] - c121 * source[287] + c122 * source[289]
                  - c148 * source[180] + c149 * source[182] - c148 * source[184]
                  - c148 * source[210] + c149 * source[212] - c148 * source[214]
                  + c150 * source[15] - c151 * source[17] + c150 * source[19]
                  + c152 * source[45] - c148 * source[47] + c152 * source[49]
                  + c150 * source[75] - c151 * source[77] + c150 * source[79];
    target[82] =  c143 * source[286] - c143 * source[288] - c121 * source[181]
                  + c121 * source[183] - c121 * source[211] + c121 * source[213]
                  + c124 * source[16] - c124 * source[18] + c122 * source[46]
                  - c122 * source[48] + c124 * source[76] - c124 * source[78];
    target[83] =  c129 * source[290] - c145 * source[292] - c128 * source[185]
                  + c153 * source[187] - c128 * source[215] + c153 * source[217]
                  + c154 * source[20] - c155 * source[22] + c156 * source[50]
                  - c130 * source[52] + c154 * source[80] - c155 * source[82];
    target[84] =  c145 * source[291] - c129 * source[293] - c153 * source[186]
                  + c128 * source[188] - c153 * source[216] + c128 * source[218]
                  + c155 * source[21] - c154 * source[23] + c130 * source[51]
                  - c156 * source[53] + c155 * source[81] - c154 * source[83];
    target[85] =  c157 * source[294] - c157 * source[296] - c158 * source[285]
                  + c158 * source[287] - c158 * source[287] + c158 * source[289]
                  - c159 * source[189] + c159 * source[191] + c160 * source[180]
                  - c160 * source[182] + c160 * source[182] - c160 * source[184]
                  - c159 * source[219] + c159 * source[221] + c160 * source[210]
                  - c160 * source[212] + c160 * source[212] - c160 * source[214]
                  + c161 * source[24] - c161 * source[26] - c162 * source[15]
                  + c162 * source[17] - c162 * source[17] + c162 * source[19]
                  + c160 * source[54] - c160 * source[56] - c163 * source[45]
                  + c163 * source[47] - c163 * source[47] + c163 * source[49]
                  + c161 * source[84] - c161 * source[86] - c162 * source[75]
                  + c162 * source[77] - c162 * source[77] + c162 * source[79];
    target[86] =  c164 * source[295] - c165 * source[286] - c165 * source[288]
                  - c166 * source[190] + c167 * source[181] + c167 * source[183]
                  - c166 * source[220] + c167 * source[211] + c167 * source[213]
                  + c160 * source[25] - c163 * source[16] - c163 * source[18]
                  + c167 * source[55] - c168 * source[46] - c168 * source[48]
                  + c160 * source[85] - c163 * source[76] - c163 * source[78];
    target[87] =  c169 * source[297] - c170 * source[290] - c170 * source[292]
                  - c171 * source[192] + c172 * source[185] + c172 * source[187]
                  - c171 * source[222] + c172 * source[215] + c172 * source[217]
                  + c173 * source[27] - c174 * source[20] - c174 * source[22]
                  + c175 * source[57] - c176 * source[50] - c176 * source[52]
                  + c173 * source[87] - c174 * source[80] - c174 * source[82];
    target[88] =  c169 * source[298] - c170 * source[291] - c170 * source[293]
                  - c171 * source[193] + c172 * source[186] + c172 * source[188]
                  - c171 * source[223] + c172 * source[216] + c172 * source[218]
                  + c173 * source[28] - c174 * source[21] - c174 * source[23]
                  + c175 * source[58] - c176 * source[51] - c176 * source[53]
                  + c173 * source[88] - c174 * source[81] - c174 * source[83];
    target[89] =  c177 * source[299] - c178 * source[294] - c178 * source[296]
                  + c179 * source[285] + c180 * source[287] + c179 * source[289]
                  - c181 * source[194] + c182 * source[189] + c182 * source[191]
                  - c183 * source[180] - c184 * source[182] - c183 * source[184]
                  - c181 * source[224] + c182 * source[219] + c182 * source[221]
                  - c183 * source[210] - c184 * source[212] - c183 * source[214]
                  + c185 * source[29] - c179 * source[24] - c179 * source[26]
                  + c186 * source[15] + c187 * source[17] + c186 * source[19]
                  + c188 * source[59] - c180 * source[54] - c180 * source[56]
                  + c187 * source[45] + c189 * source[47] + c187 * source[49]
                  + c185 * source[89] - c179 * source[84] - c179 * source[86]
                  + c186 * source[75] + c187 * source[77] + c186 * source[79];
    target[90] =  c190 * source[300] - c191 * source[302] + c190 * source[304]
                  - c192 * source[225] + c27 * source[227] - c192 * source[229]
                  - c192 * source[255] + c27 * source[257] - c192 * source[259]
                  + c193 * source[90] - c30 * source[92] + c193 * source[94]
                  + c194 * source[120] - c28 * source[122] + c194 * source[124]
                  + c193 * source[150] - c30 * source[152] + c193 * source[154];
    target[91] =  c195 * source[301] - c195 * source[303] - c196 * source[226]
                  + c196 * source[228] - c196 * source[256] + c196 * source[258]
                  + c197 * source[91] - c197 * source[93] + c29 * source[121]
                  - c29 * source[123] + c197 * source[151] - c197 * source[153];
    target[92] =  c198 * source[305] - c199 * source[307] - c200 * source[230]
                  + c201 * source[232] - c200 * source[260] + c201 * source[262]
                  + c18 * source[95] - c19 * source[97] + c24 * source[125]
                  - c17 * source[127] + c18 * source[155] - c19 * source[157];
    target[93] =  c199 * source[306] - c198 * source[308] - c201 * source[231]
                  + c200 * source[233] - c201 * source[261] + c200 * source[263]
                  + c19 * source[96] - c18 * source[98] + c17 * source[126]
                  - c24 * source[128] + c19 * source[156] - c18 * source[158];
    target[94] =  c202 * source[309] - c202 * source[311] - c203 * source[300]
                  + c203 * source[302] - c203 * source[302] + c203 * source[304]
                  - c204 * source[234] + c204 * source[236] + c205 * source[225]
                  - c205 * source[227] + c205 * source[227] - c205 * source[229]
                  - c204 * source[264] + c204 * source[266] + c205 * source[255]
                  - c205 * source[257] + c205 * source[257] - c205 * source[259]
                  + c206 * source[99] - c206 * source[101] - c207 * source[90]
                  + c207 * source[92] - c207 * source[92] + c207 * source[94]
                  + c208 * source[129] - c208 * source[131] - c209 * source[120]
                  + c209 * source[122] - c209 * source[122] + c209 * source[124]
                  + c206 * source[159] - c206 * source[161] - c207 * source[150]
                  + c207 * source[152] - c207 * source[152] + c207 * source[154];
    target[95] =  c210 * source[310] - c211 * source[301] - c211 * source[303]
                  - c212 * source[235] + c213 * source[226] + c213 * source[228]
                  - c212 * source[265] + c213 * source[256] + c213 * source[258]
                  + c208 * source[100] - c209 * source[91] - c209 * source[93]
                  + c214 * source[130] - c215 * source[121] - c215 * source[123]
                  + c208 * source[160] - c209 * source[151] - c209 * source[153];
    target[96] =  c216 * source[312] - c217 * source[305] - c217 * source[307]
                  - c218 * source[237] + c219 * source[230] + c219 * source[232]
                  - c218 * source[267] + c219 * source[260] + c219 * source[262]
                  + c220 * source[102] - c221 * source[95] - c221 * source[97]
                  + c219 * source[132] - c222 * source[125] - c222 * source[127]
                  + c220 * source[162] - c221 * source[155] - c221 * source[157];
    target[97] =  c216 * source[313] - c217 * source[306] - c217 * source[308]
                  - c218 * source[238] + c219 * source[231] + c219 * source[233]
                  - c218 * source[268] + c219 * source[261] + c219 * source[263]
                  + c220 * source[103] - c221 * source[96] - c221 * source[98]
                  + c219 * source[133] - c222 * source[126] - c222 * source[128]
                  + c220 * source[163] - c221 * source[156] - c221 * source[158];
    target[98] =  source[314] - c223 * source[309] - c223 * source[311]
                  + c224 * source[300] + c225 * source[302] + c224 * source[304]
                  - c226 * source[239] + c227 * source[234] + c227 * source[236]
                  - c228 * source[225] - c229 * source[227] - c228 * source[229]
                  - c226 * source[269] + c227 * source[264] + c227 * source[266]
                  - c228 * source[255] - c229 * source[257] - c228 * source[259]
                  + c228 * source[104] - c230 * source[99] - c230 * source[101]
                  + c231 * source[90] + c232 * source[92] + c231 * source[94]
                  + c229 * source[134] - c233 * source[129] - c233 * source[131]
                  + c232 * source[120] + c234 * source[122] + c232 * source[124]
                  + c228 * source[164] - c230 * source[159] - c230 * source[161]
                  + c231 * source[150] + c232 * source[152] + c231 * source[154];
  }
}

