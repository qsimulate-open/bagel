//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_43.cc
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


void CarSphList::carsph_43(const int nloop, const double* source, double* target) {
  const double c56 = 25.98076211353316;
  const double c28 = 24.302777619029477;
  const double c7 = 17.1846588560844;
  const double c58 = 16.431676725154983;
  const double c54 = 15.90990257669732;
  const double c31 = 15.370426148939398;
  const double c24 = 14.882351124738323;
  const double c43 = 12.99038105676658;
  const double c72 = 12.24744871391589;
  const double c26 = 12.151388809514739;
  const double c90 = 11.61895003862225;
  const double c19 = 11.456439237389599;
  const double c10 = 10.868532559642079;
  const double c3 = 10.523411400301709;
  const double c63 = 10.062305898749054;
  const double c36 = 9.4124252985083494;
  const double c73 = 9.1855865354369186;
  const double c5 = 8.5923294280422002;
  const double c45 = 8.2158383625774913;
  const double c27 = 8.1009258730098246;
  const double c38 = 7.9549512883486599;
  const double c74 = 7.745966692414834;
  const double c67 = 7.5;
  const double c94 = 7.3484692283495345;
  const double c20 = 7.245688373094719;
  const double c83 = 7.1151247353788536;
  const double c17 = 7.0156076002011405;
  const double c62 = 6.7082039324993694;
  const double c15 = 6.6555897559870685;
  const double c41 = 6.49519052838329;
  const double c35 = 6.2749501990055663;
  const double c70 = 6.1237243569579451;
  const double c76 = 5.8094750193111251;
  const double c18 = 5.7282196186947996;
  const double c69 = 5.625;
  const double c53 = 5.3033008588991066;
  const double c29 = 5.123475382979799;
  const double c50 = 5.0311529493745271;
  const double c23 = 4.9607837082461073;
  const double c79 = 4.7434164902525691;
  const double c71 = 4.5927932677184593;
  const double c101 = 4.5;
  const double c14 = 4.4370598373247123;
  const double c57 = 4.3301270189221936;
  const double c59 = 4.1079191812887457;
  const double c25 = 4.0504629365049123;
  const double c89 = 3.872983346207417;
  const double c32 = 3.8426065372348495;
  const double c81 = 3.5575623676894268;
  const double c2 = 3.5078038001005702;
  const double c49 = 3.3541019662496847;
  const double c78 = 3.1622776601683795;
  const double c34 = 3.1374750995027831;
  const double c100 = 3;
  const double c21 = 2.9580398915498081;
  const double c91 = 2.9047375096555625;
  const double c6 = 2.8641098093473998;
  const double c60 = 2.7386127875258306;
  const double c11 = 2.7171331399105196;
  const double c37 = 2.6516504294495533;
  const double c66 = 2.5;
  const double c92 = 2.4494897427831779;
  const double c80 = 2.3717082451262845;
  const double c16 = 2.3385358667337135;
  const double c44 = 2.1650635094610968;
  const double c33 = 2.0916500663351889;
  const double c46 = 2.0539595906443728;
  const double c75 = 1.9364916731037085;
  const double c68 = 1.875;
  const double c95 = 1.8371173070873836;
  const double c8 = 1.8114220932736798;
  const double c87 = 1.7787811838447134;
  const double c1 = 1.7539019000502851;
  const double c65 = 1.6770509831248424;
  const double c22 = 1.6535945694153691;
  const double c99 = 1.5;
  const double c77 = 1.4523687548277813;
  const double c4 = 1.4320549046736999;
  const double c47 = 1.3693063937629153;
  const double c40 = 1.3258252147247767;
  const double c30 = 1.2808688457449497;
  const double c105 = 1.125;
  const double c64 = 1.1180339887498949;
  const double c13 = 1.1092649593311781;
  const double c42 = 1.0825317547305484;
  const double c96 = 0.91855865354369182;
  const double c85 = 0.8893905919223567;
  const double c55 = 0.88388347648318444;
  const double c52 = 0.83852549156242118;
  const double c82 = 0.79056941504209488;
  const double c104 = 0.75;
  const double c12 = 0.73950997288745202;
  const double c88 = 0.72618437741389064;
  const double c61 = 0.68465319688145765;
  const double c93 = 0.61237243569579447;
  const double c86 = 0.59292706128157113;
  const double c0 = 0.58463396668342837;
  const double c103 = 0.5625;
  const double c51 = 0.55901699437494745;
  const double c98 = 0.45927932677184591;
  const double c9 = 0.45285552331841994;
  const double c39 = 0.44194173824159222;
  const double c102 = 0.375;
  const double c48 = 0.34232659844072882;
  const double c84 = 0.29646353064078557;
  const double c97 = 0.22963966338592295;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 63, source += 150) {
    target[0] =  c0 * source[0] - c1 * source[2] - c2 * source[20]
                  + c3 * source[22] + c0 * source[40] - c1 * source[42];
    target[1] =  c1 * source[1] - c0 * source[3] - c3 * source[21]
                  + c2 * source[23] + c1 * source[41] - c0 * source[43];
    target[2] =  c4 * source[4] - c4 * source[6] - c5 * source[24]
                  + c5 * source[26] + c4 * source[44] - c4 * source[46];
    target[3] =  c6 * source[5] - c7 * source[25] + c6 * source[45];
    target[4] =  c8 * source[7] - c9 * source[0] - c9 * source[2]
                  - c10 * source[27] + c11 * source[20] + c11 * source[22]
                  + c8 * source[47] - c9 * source[40] - c9 * source[42];
    target[5] =  c8 * source[8] - c9 * source[1] - c9 * source[3]
                  - c10 * source[28] + c11 * source[21] + c11 * source[23]
                  + c8 * source[48] - c9 * source[41] - c9 * source[43];
    target[6] =  c12 * source[9] - c13 * source[4] - c13 * source[6]
                  - c14 * source[29] + c15 * source[24] + c15 * source[26]
                  + c12 * source[49] - c13 * source[44] - c13 * source[46];
    target[7] =  c16 * source[10] - c17 * source[12] - c16 * source[30]
                  + c17 * source[32];
    target[8] =  c17 * source[11] - c16 * source[13] - c17 * source[31]
                  + c16 * source[33];
    target[9] =  c18 * source[14] - c18 * source[16] - c18 * source[34]
                  + c18 * source[36];
    target[10] =  c19 * source[15] - c19 * source[35];
    target[11] =  c20 * source[17] - c8 * source[10] - c8 * source[12]
                  - c20 * source[37] + c8 * source[30] + c8 * source[32];
    target[12] =  c20 * source[18] - c8 * source[11] - c8 * source[13]
                  - c20 * source[38] + c8 * source[31] + c8 * source[33];
    target[13] =  c21 * source[19] - c14 * source[14] - c14 * source[16]
                  - c21 * source[39] + c14 * source[34] + c14 * source[36];
    target[14] =  c22 * source[50] - c23 * source[52] - c23 * source[70]
                  + c24 * source[72];
    target[15] =  c23 * source[51] - c22 * source[53] - c24 * source[71]
                  + c23 * source[73];
    target[16] =  c25 * source[54] - c25 * source[56] - c26 * source[74]
                  + c26 * source[76];
    target[17] =  c27 * source[55] - c28 * source[75];
    target[18] =  c29 * source[57] - c30 * source[50] - c30 * source[52]
                  - c31 * source[77] + c32 * source[70] + c32 * source[72];
    target[19] =  c29 * source[58] - c30 * source[51] - c30 * source[53]
                  - c31 * source[78] + c32 * source[71] + c32 * source[73];
    target[20] =  c33 * source[59] - c34 * source[54] - c34 * source[56]
                  - c35 * source[79] + c36 * source[74] + c36 * source[76];
    target[21] =  c23 * source[60] - c24 * source[62] - c22 * source[80]
                  + c23 * source[82];
    target[22] =  c24 * source[61] - c23 * source[63] - c23 * source[81]
                  + c22 * source[83];
    target[23] =  c26 * source[64] - c26 * source[66] - c25 * source[84]
                  + c25 * source[86];
    target[24] =  c28 * source[65] - c27 * source[85];
    target[25] =  c31 * source[67] - c32 * source[60] - c32 * source[62]
                  - c29 * source[87] + c30 * source[80] + c30 * source[82];
    target[26] =  c31 * source[68] - c32 * source[61] - c32 * source[63]
                  - c29 * source[88] + c30 * source[81] + c30 * source[83];
    target[27] =  c35 * source[69] - c36 * source[64] - c36 * source[66]
                  - c33 * source[89] + c34 * source[84] + c34 * source[86];
    target[28] =  c37 * source[90] - c38 * source[92] - c37 * source[110]
                  + c38 * source[112] - c39 * source[0] + c40 * source[2]
                  + c39 * source[20] - c40 * source[22] - c39 * source[20]
                  + c40 * source[22] + c39 * source[40] - c40 * source[42];
    target[29] =  c38 * source[91] - c37 * source[93] - c38 * source[111]
                  + c37 * source[113] - c40 * source[1] + c39 * source[3]
                  + c40 * source[21] - c39 * source[23] - c40 * source[21]
                  + c39 * source[23] + c40 * source[41] - c39 * source[43];
    target[30] =  c41 * source[94] - c41 * source[96] - c41 * source[114]
                  + c41 * source[116] - c42 * source[4] + c42 * source[6]
                  + c42 * source[24] - c42 * source[26] - c42 * source[24]
                  + c42 * source[26] + c42 * source[44] - c42 * source[46];
    target[31] =  c43 * source[95] - c43 * source[115] - c44 * source[5]
                  + c44 * source[25] - c44 * source[25] + c44 * source[45];
    target[32] =  c45 * source[97] - c46 * source[90] - c46 * source[92]
                  - c45 * source[117] + c46 * source[110] + c46 * source[112]
                  - c47 * source[7] + c48 * source[0] + c48 * source[2]
                  + c47 * source[27] - c48 * source[20] - c48 * source[22]
                  - c47 * source[27] + c48 * source[20] + c48 * source[22]
                  + c47 * source[47] - c48 * source[40] - c48 * source[42];
    target[33] =  c45 * source[98] - c46 * source[91] - c46 * source[93]
                  - c45 * source[118] + c46 * source[111] + c46 * source[113]
                  - c47 * source[8] + c48 * source[1] + c48 * source[3]
                  + c47 * source[28] - c48 * source[21] - c48 * source[23]
                  - c47 * source[28] + c48 * source[21] + c48 * source[23]
                  + c47 * source[48] - c48 * source[41] - c48 * source[43];
    target[34] =  c49 * source[99] - c50 * source[94] - c50 * source[96]
                  - c49 * source[119] + c50 * source[114] + c50 * source[116]
                  - c51 * source[9] + c52 * source[4] + c52 * source[6]
                  + c51 * source[29] - c52 * source[24] - c52 * source[26]
                  - c51 * source[29] + c52 * source[24] + c52 * source[26]
                  + c51 * source[49] - c52 * source[44] - c52 * source[46];
    target[35] =  c53 * source[100] - c54 * source[102] - c55 * source[10]
                  + c37 * source[12] - c55 * source[30] + c37 * source[32];
    target[36] =  c54 * source[101] - c53 * source[103] - c37 * source[11]
                  + c55 * source[13] - c37 * source[31] + c55 * source[33];
    target[37] =  c43 * source[104] - c43 * source[106] - c44 * source[14]
                  + c44 * source[16] - c44 * source[34] + c44 * source[36];
    target[38] =  c56 * source[105] - c57 * source[15] - c57 * source[35];
    target[39] =  c58 * source[107] - c59 * source[100] - c59 * source[102]
                  - c60 * source[17] + c61 * source[10] + c61 * source[12]
                  - c60 * source[37] + c61 * source[30] + c61 * source[32];
    target[40] =  c58 * source[108] - c59 * source[101] - c59 * source[103]
                  - c60 * source[18] + c61 * source[11] + c61 * source[13]
                  - c60 * source[38] + c61 * source[31] + c61 * source[33];
    target[41] =  c62 * source[109] - c63 * source[104] - c63 * source[106]
                  - c64 * source[19] + c65 * source[14] + c65 * source[16]
                  - c64 * source[39] + c65 * source[34] + c65 * source[36];
    target[42] =  c66 * source[120] - c67 * source[122] - c68 * source[50]
                  + c69 * source[52] - c68 * source[70] + c69 * source[72];
    target[43] =  c67 * source[121] - c66 * source[123] - c69 * source[51]
                  + c68 * source[53] - c69 * source[71] + c68 * source[73];
    target[44] =  c70 * source[124] - c70 * source[126] - c71 * source[54]
                  + c71 * source[56] - c71 * source[74] + c71 * source[76];
    target[45] =  c72 * source[125] - c73 * source[55] - c73 * source[75];
    target[46] =  c74 * source[127] - c75 * source[120] - c75 * source[122]
                  - c76 * source[57] + c77 * source[50] + c77 * source[52]
                  - c76 * source[77] + c77 * source[70] + c77 * source[72];
    target[47] =  c74 * source[128] - c75 * source[121] - c75 * source[123]
                  - c76 * source[58] + c77 * source[51] + c77 * source[53]
                  - c76 * source[78] + c77 * source[71] + c77 * source[73];
    target[48] =  c78 * source[129] - c79 * source[124] - c79 * source[126]
                  - c80 * source[59] + c81 * source[54] + c81 * source[56]
                  - c80 * source[79] + c81 * source[74] + c81 * source[76];
    target[49] =  c66 * source[130] - c67 * source[132] - c68 * source[60]
                  + c69 * source[62] - c68 * source[80] + c69 * source[82];
    target[50] =  c67 * source[131] - c66 * source[133] - c69 * source[61]
                  + c68 * source[63] - c69 * source[81] + c68 * source[83];
    target[51] =  c70 * source[134] - c70 * source[136] - c71 * source[64]
                  + c71 * source[66] - c71 * source[84] + c71 * source[86];
    target[52] =  c72 * source[135] - c73 * source[65] - c73 * source[85];
    target[53] =  c74 * source[137] - c75 * source[130] - c75 * source[132]
                  - c76 * source[67] + c77 * source[60] + c77 * source[62]
                  - c76 * source[87] + c77 * source[80] + c77 * source[82];
    target[54] =  c74 * source[138] - c75 * source[131] - c75 * source[133]
                  - c76 * source[68] + c77 * source[61] + c77 * source[63]
                  - c76 * source[88] + c77 * source[81] + c77 * source[83];
    target[55] =  c78 * source[139] - c79 * source[134] - c79 * source[136]
                  - c80 * source[69] + c81 * source[64] + c81 * source[66]
                  - c80 * source[89] + c81 * source[84] + c81 * source[86];
    target[56] =  c82 * source[140] - c80 * source[142] - c80 * source[90]
                  + c83 * source[92] - c80 * source[110] + c83 * source[112]
                  + c84 * source[0] - c85 * source[2] + c86 * source[20]
                  - c87 * source[22] + c84 * source[40] - c85 * source[42];
    target[57] =  c80 * source[141] - c82 * source[143] - c83 * source[91]
                  + c80 * source[93] - c83 * source[111] + c80 * source[113]
                  + c85 * source[1] - c84 * source[3] + c87 * source[21]
                  - c86 * source[23] + c85 * source[41] - c84 * source[43];
    target[58] =  c75 * source[144] - c75 * source[146] - c76 * source[94]
                  + c76 * source[96] - c76 * source[114] + c76 * source[116]
                  + c88 * source[4] - c88 * source[6] + c77 * source[24]
                  - c77 * source[26] + c88 * source[44] - c88 * source[46];
    target[59] =  c89 * source[145] - c90 * source[95] - c90 * source[115]
                  + c77 * source[5] + c91 * source[25] + c77 * source[45];
    target[60] =  c92 * source[147] - c93 * source[140] - c93 * source[142]
                  - c94 * source[97] + c95 * source[90] + c95 * source[92]
                  - c94 * source[117] + c95 * source[110] + c95 * source[112]
                  + c96 * source[7] - c97 * source[0] - c97 * source[2]
                  + c95 * source[27] - c98 * source[20] - c98 * source[22]
                  + c96 * source[47] - c97 * source[40] - c97 * source[42];
    target[61] =  c92 * source[148] - c93 * source[141] - c93 * source[143]
                  - c94 * source[98] + c95 * source[91] + c95 * source[93]
                  - c94 * source[118] + c95 * source[111] + c95 * source[113]
                  + c96 * source[8] - c97 * source[1] - c97 * source[3]
                  + c95 * source[28] - c98 * source[21] - c98 * source[23]
                  + c96 * source[48] - c97 * source[41] - c97 * source[43];
    target[62] =  source[149] - c99 * source[144] - c99 * source[146]
                  - c100 * source[99] + c101 * source[94] + c101 * source[96]
                  - c100 * source[119] + c101 * source[114] + c101 * source[116]
                  + c102 * source[9] - c103 * source[4] - c103 * source[6]
                  + c104 * source[29] - c105 * source[24] - c105 * source[26]
                  + c102 * source[49] - c103 * source[44] - c103 * source[46];
  }
}

void CCarSphList::carsph_43(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c56 = 25.98076211353316;
  const double c28 = 24.302777619029477;
  const double c7 = 17.1846588560844;
  const double c58 = 16.431676725154983;
  const double c54 = 15.90990257669732;
  const double c31 = 15.370426148939398;
  const double c24 = 14.882351124738323;
  const double c43 = 12.99038105676658;
  const double c72 = 12.24744871391589;
  const double c26 = 12.151388809514739;
  const double c90 = 11.61895003862225;
  const double c19 = 11.456439237389599;
  const double c10 = 10.868532559642079;
  const double c3 = 10.523411400301709;
  const double c63 = 10.062305898749054;
  const double c36 = 9.4124252985083494;
  const double c73 = 9.1855865354369186;
  const double c5 = 8.5923294280422002;
  const double c45 = 8.2158383625774913;
  const double c27 = 8.1009258730098246;
  const double c38 = 7.9549512883486599;
  const double c74 = 7.745966692414834;
  const double c67 = 7.5;
  const double c94 = 7.3484692283495345;
  const double c20 = 7.245688373094719;
  const double c83 = 7.1151247353788536;
  const double c17 = 7.0156076002011405;
  const double c62 = 6.7082039324993694;
  const double c15 = 6.6555897559870685;
  const double c41 = 6.49519052838329;
  const double c35 = 6.2749501990055663;
  const double c70 = 6.1237243569579451;
  const double c76 = 5.8094750193111251;
  const double c18 = 5.7282196186947996;
  const double c69 = 5.625;
  const double c53 = 5.3033008588991066;
  const double c29 = 5.123475382979799;
  const double c50 = 5.0311529493745271;
  const double c23 = 4.9607837082461073;
  const double c79 = 4.7434164902525691;
  const double c71 = 4.5927932677184593;
  const double c101 = 4.5;
  const double c14 = 4.4370598373247123;
  const double c57 = 4.3301270189221936;
  const double c59 = 4.1079191812887457;
  const double c25 = 4.0504629365049123;
  const double c89 = 3.872983346207417;
  const double c32 = 3.8426065372348495;
  const double c81 = 3.5575623676894268;
  const double c2 = 3.5078038001005702;
  const double c49 = 3.3541019662496847;
  const double c78 = 3.1622776601683795;
  const double c34 = 3.1374750995027831;
  const double c100 = 3;
  const double c21 = 2.9580398915498081;
  const double c91 = 2.9047375096555625;
  const double c6 = 2.8641098093473998;
  const double c60 = 2.7386127875258306;
  const double c11 = 2.7171331399105196;
  const double c37 = 2.6516504294495533;
  const double c66 = 2.5;
  const double c92 = 2.4494897427831779;
  const double c80 = 2.3717082451262845;
  const double c16 = 2.3385358667337135;
  const double c44 = 2.1650635094610968;
  const double c33 = 2.0916500663351889;
  const double c46 = 2.0539595906443728;
  const double c75 = 1.9364916731037085;
  const double c68 = 1.875;
  const double c95 = 1.8371173070873836;
  const double c8 = 1.8114220932736798;
  const double c87 = 1.7787811838447134;
  const double c1 = 1.7539019000502851;
  const double c65 = 1.6770509831248424;
  const double c22 = 1.6535945694153691;
  const double c99 = 1.5;
  const double c77 = 1.4523687548277813;
  const double c4 = 1.4320549046736999;
  const double c47 = 1.3693063937629153;
  const double c40 = 1.3258252147247767;
  const double c30 = 1.2808688457449497;
  const double c105 = 1.125;
  const double c64 = 1.1180339887498949;
  const double c13 = 1.1092649593311781;
  const double c42 = 1.0825317547305484;
  const double c96 = 0.91855865354369182;
  const double c85 = 0.8893905919223567;
  const double c55 = 0.88388347648318444;
  const double c52 = 0.83852549156242118;
  const double c82 = 0.79056941504209488;
  const double c104 = 0.75;
  const double c12 = 0.73950997288745202;
  const double c88 = 0.72618437741389064;
  const double c61 = 0.68465319688145765;
  const double c93 = 0.61237243569579447;
  const double c86 = 0.59292706128157113;
  const double c0 = 0.58463396668342837;
  const double c103 = 0.5625;
  const double c51 = 0.55901699437494745;
  const double c98 = 0.45927932677184591;
  const double c9 = 0.45285552331841994;
  const double c39 = 0.44194173824159222;
  const double c102 = 0.375;
  const double c48 = 0.34232659844072882;
  const double c84 = 0.29646353064078557;
  const double c97 = 0.22963966338592295;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 63, source += 150) {
    target[0] =  c0 * source[0] - c1 * source[2] - c2 * source[20]
                  + c3 * source[22] + c0 * source[40] - c1 * source[42];
    target[1] =  c1 * source[1] - c0 * source[3] - c3 * source[21]
                  + c2 * source[23] + c1 * source[41] - c0 * source[43];
    target[2] =  c4 * source[4] - c4 * source[6] - c5 * source[24]
                  + c5 * source[26] + c4 * source[44] - c4 * source[46];
    target[3] =  c6 * source[5] - c7 * source[25] + c6 * source[45];
    target[4] =  c8 * source[7] - c9 * source[0] - c9 * source[2]
                  - c10 * source[27] + c11 * source[20] + c11 * source[22]
                  + c8 * source[47] - c9 * source[40] - c9 * source[42];
    target[5] =  c8 * source[8] - c9 * source[1] - c9 * source[3]
                  - c10 * source[28] + c11 * source[21] + c11 * source[23]
                  + c8 * source[48] - c9 * source[41] - c9 * source[43];
    target[6] =  c12 * source[9] - c13 * source[4] - c13 * source[6]
                  - c14 * source[29] + c15 * source[24] + c15 * source[26]
                  + c12 * source[49] - c13 * source[44] - c13 * source[46];
    target[7] =  c16 * source[10] - c17 * source[12] - c16 * source[30]
                  + c17 * source[32];
    target[8] =  c17 * source[11] - c16 * source[13] - c17 * source[31]
                  + c16 * source[33];
    target[9] =  c18 * source[14] - c18 * source[16] - c18 * source[34]
                  + c18 * source[36];
    target[10] =  c19 * source[15] - c19 * source[35];
    target[11] =  c20 * source[17] - c8 * source[10] - c8 * source[12]
                  - c20 * source[37] + c8 * source[30] + c8 * source[32];
    target[12] =  c20 * source[18] - c8 * source[11] - c8 * source[13]
                  - c20 * source[38] + c8 * source[31] + c8 * source[33];
    target[13] =  c21 * source[19] - c14 * source[14] - c14 * source[16]
                  - c21 * source[39] + c14 * source[34] + c14 * source[36];
    target[14] =  c22 * source[50] - c23 * source[52] - c23 * source[70]
                  + c24 * source[72];
    target[15] =  c23 * source[51] - c22 * source[53] - c24 * source[71]
                  + c23 * source[73];
    target[16] =  c25 * source[54] - c25 * source[56] - c26 * source[74]
                  + c26 * source[76];
    target[17] =  c27 * source[55] - c28 * source[75];
    target[18] =  c29 * source[57] - c30 * source[50] - c30 * source[52]
                  - c31 * source[77] + c32 * source[70] + c32 * source[72];
    target[19] =  c29 * source[58] - c30 * source[51] - c30 * source[53]
                  - c31 * source[78] + c32 * source[71] + c32 * source[73];
    target[20] =  c33 * source[59] - c34 * source[54] - c34 * source[56]
                  - c35 * source[79] + c36 * source[74] + c36 * source[76];
    target[21] =  c23 * source[60] - c24 * source[62] - c22 * source[80]
                  + c23 * source[82];
    target[22] =  c24 * source[61] - c23 * source[63] - c23 * source[81]
                  + c22 * source[83];
    target[23] =  c26 * source[64] - c26 * source[66] - c25 * source[84]
                  + c25 * source[86];
    target[24] =  c28 * source[65] - c27 * source[85];
    target[25] =  c31 * source[67] - c32 * source[60] - c32 * source[62]
                  - c29 * source[87] + c30 * source[80] + c30 * source[82];
    target[26] =  c31 * source[68] - c32 * source[61] - c32 * source[63]
                  - c29 * source[88] + c30 * source[81] + c30 * source[83];
    target[27] =  c35 * source[69] - c36 * source[64] - c36 * source[66]
                  - c33 * source[89] + c34 * source[84] + c34 * source[86];
    target[28] =  c37 * source[90] - c38 * source[92] - c37 * source[110]
                  + c38 * source[112] - c39 * source[0] + c40 * source[2]
                  + c39 * source[20] - c40 * source[22] - c39 * source[20]
                  + c40 * source[22] + c39 * source[40] - c40 * source[42];
    target[29] =  c38 * source[91] - c37 * source[93] - c38 * source[111]
                  + c37 * source[113] - c40 * source[1] + c39 * source[3]
                  + c40 * source[21] - c39 * source[23] - c40 * source[21]
                  + c39 * source[23] + c40 * source[41] - c39 * source[43];
    target[30] =  c41 * source[94] - c41 * source[96] - c41 * source[114]
                  + c41 * source[116] - c42 * source[4] + c42 * source[6]
                  + c42 * source[24] - c42 * source[26] - c42 * source[24]
                  + c42 * source[26] + c42 * source[44] - c42 * source[46];
    target[31] =  c43 * source[95] - c43 * source[115] - c44 * source[5]
                  + c44 * source[25] - c44 * source[25] + c44 * source[45];
    target[32] =  c45 * source[97] - c46 * source[90] - c46 * source[92]
                  - c45 * source[117] + c46 * source[110] + c46 * source[112]
                  - c47 * source[7] + c48 * source[0] + c48 * source[2]
                  + c47 * source[27] - c48 * source[20] - c48 * source[22]
                  - c47 * source[27] + c48 * source[20] + c48 * source[22]
                  + c47 * source[47] - c48 * source[40] - c48 * source[42];
    target[33] =  c45 * source[98] - c46 * source[91] - c46 * source[93]
                  - c45 * source[118] + c46 * source[111] + c46 * source[113]
                  - c47 * source[8] + c48 * source[1] + c48 * source[3]
                  + c47 * source[28] - c48 * source[21] - c48 * source[23]
                  - c47 * source[28] + c48 * source[21] + c48 * source[23]
                  + c47 * source[48] - c48 * source[41] - c48 * source[43];
    target[34] =  c49 * source[99] - c50 * source[94] - c50 * source[96]
                  - c49 * source[119] + c50 * source[114] + c50 * source[116]
                  - c51 * source[9] + c52 * source[4] + c52 * source[6]
                  + c51 * source[29] - c52 * source[24] - c52 * source[26]
                  - c51 * source[29] + c52 * source[24] + c52 * source[26]
                  + c51 * source[49] - c52 * source[44] - c52 * source[46];
    target[35] =  c53 * source[100] - c54 * source[102] - c55 * source[10]
                  + c37 * source[12] - c55 * source[30] + c37 * source[32];
    target[36] =  c54 * source[101] - c53 * source[103] - c37 * source[11]
                  + c55 * source[13] - c37 * source[31] + c55 * source[33];
    target[37] =  c43 * source[104] - c43 * source[106] - c44 * source[14]
                  + c44 * source[16] - c44 * source[34] + c44 * source[36];
    target[38] =  c56 * source[105] - c57 * source[15] - c57 * source[35];
    target[39] =  c58 * source[107] - c59 * source[100] - c59 * source[102]
                  - c60 * source[17] + c61 * source[10] + c61 * source[12]
                  - c60 * source[37] + c61 * source[30] + c61 * source[32];
    target[40] =  c58 * source[108] - c59 * source[101] - c59 * source[103]
                  - c60 * source[18] + c61 * source[11] + c61 * source[13]
                  - c60 * source[38] + c61 * source[31] + c61 * source[33];
    target[41] =  c62 * source[109] - c63 * source[104] - c63 * source[106]
                  - c64 * source[19] + c65 * source[14] + c65 * source[16]
                  - c64 * source[39] + c65 * source[34] + c65 * source[36];
    target[42] =  c66 * source[120] - c67 * source[122] - c68 * source[50]
                  + c69 * source[52] - c68 * source[70] + c69 * source[72];
    target[43] =  c67 * source[121] - c66 * source[123] - c69 * source[51]
                  + c68 * source[53] - c69 * source[71] + c68 * source[73];
    target[44] =  c70 * source[124] - c70 * source[126] - c71 * source[54]
                  + c71 * source[56] - c71 * source[74] + c71 * source[76];
    target[45] =  c72 * source[125] - c73 * source[55] - c73 * source[75];
    target[46] =  c74 * source[127] - c75 * source[120] - c75 * source[122]
                  - c76 * source[57] + c77 * source[50] + c77 * source[52]
                  - c76 * source[77] + c77 * source[70] + c77 * source[72];
    target[47] =  c74 * source[128] - c75 * source[121] - c75 * source[123]
                  - c76 * source[58] + c77 * source[51] + c77 * source[53]
                  - c76 * source[78] + c77 * source[71] + c77 * source[73];
    target[48] =  c78 * source[129] - c79 * source[124] - c79 * source[126]
                  - c80 * source[59] + c81 * source[54] + c81 * source[56]
                  - c80 * source[79] + c81 * source[74] + c81 * source[76];
    target[49] =  c66 * source[130] - c67 * source[132] - c68 * source[60]
                  + c69 * source[62] - c68 * source[80] + c69 * source[82];
    target[50] =  c67 * source[131] - c66 * source[133] - c69 * source[61]
                  + c68 * source[63] - c69 * source[81] + c68 * source[83];
    target[51] =  c70 * source[134] - c70 * source[136] - c71 * source[64]
                  + c71 * source[66] - c71 * source[84] + c71 * source[86];
    target[52] =  c72 * source[135] - c73 * source[65] - c73 * source[85];
    target[53] =  c74 * source[137] - c75 * source[130] - c75 * source[132]
                  - c76 * source[67] + c77 * source[60] + c77 * source[62]
                  - c76 * source[87] + c77 * source[80] + c77 * source[82];
    target[54] =  c74 * source[138] - c75 * source[131] - c75 * source[133]
                  - c76 * source[68] + c77 * source[61] + c77 * source[63]
                  - c76 * source[88] + c77 * source[81] + c77 * source[83];
    target[55] =  c78 * source[139] - c79 * source[134] - c79 * source[136]
                  - c80 * source[69] + c81 * source[64] + c81 * source[66]
                  - c80 * source[89] + c81 * source[84] + c81 * source[86];
    target[56] =  c82 * source[140] - c80 * source[142] - c80 * source[90]
                  + c83 * source[92] - c80 * source[110] + c83 * source[112]
                  + c84 * source[0] - c85 * source[2] + c86 * source[20]
                  - c87 * source[22] + c84 * source[40] - c85 * source[42];
    target[57] =  c80 * source[141] - c82 * source[143] - c83 * source[91]
                  + c80 * source[93] - c83 * source[111] + c80 * source[113]
                  + c85 * source[1] - c84 * source[3] + c87 * source[21]
                  - c86 * source[23] + c85 * source[41] - c84 * source[43];
    target[58] =  c75 * source[144] - c75 * source[146] - c76 * source[94]
                  + c76 * source[96] - c76 * source[114] + c76 * source[116]
                  + c88 * source[4] - c88 * source[6] + c77 * source[24]
                  - c77 * source[26] + c88 * source[44] - c88 * source[46];
    target[59] =  c89 * source[145] - c90 * source[95] - c90 * source[115]
                  + c77 * source[5] + c91 * source[25] + c77 * source[45];
    target[60] =  c92 * source[147] - c93 * source[140] - c93 * source[142]
                  - c94 * source[97] + c95 * source[90] + c95 * source[92]
                  - c94 * source[117] + c95 * source[110] + c95 * source[112]
                  + c96 * source[7] - c97 * source[0] - c97 * source[2]
                  + c95 * source[27] - c98 * source[20] - c98 * source[22]
                  + c96 * source[47] - c97 * source[40] - c97 * source[42];
    target[61] =  c92 * source[148] - c93 * source[141] - c93 * source[143]
                  - c94 * source[98] + c95 * source[91] + c95 * source[93]
                  - c94 * source[118] + c95 * source[111] + c95 * source[113]
                  + c96 * source[8] - c97 * source[1] - c97 * source[3]
                  + c95 * source[28] - c98 * source[21] - c98 * source[23]
                  + c96 * source[48] - c97 * source[41] - c97 * source[43];
    target[62] =  source[149] - c99 * source[144] - c99 * source[146]
                  - c100 * source[99] + c101 * source[94] + c101 * source[96]
                  - c100 * source[119] + c101 * source[114] + c101 * source[116]
                  + c102 * source[9] - c103 * source[4] - c103 * source[6]
                  + c104 * source[29] - c105 * source[24] - c105 * source[26]
                  + c102 * source[49] - c103 * source[44] - c103 * source[46];
  }
}

