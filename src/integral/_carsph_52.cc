//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_52.cc
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


void CarSphList::carsph_52(const int nloop, const double* source, double* target) {
  const double c13 = 23.055639223409095;
  const double c27 = 21.737065119284157;
  const double c41 = 17.748239349298849;
  const double c19 = 15.370426148939398;
  const double c16 = 13.311179511974137;
  const double c32 = 12.549900398011133;
  const double c4 = 12.151388809514739;
  const double c11 = 11.527819611704547;
  const double c23 = 10.868532559642079;
  const double c42 = 10.246950765959598;
  const double c48 = 10.062305898749054;
  const double c20 = 8.8741196746494246;
  const double c62 = 8.6602540378443873;
  const double c18 = 7.6852130744696989;
  const double c26 = 7.245688373094719;
  const double c7 = 7.0156076002011405;
  const double c47 = 6.7082039324993694;
  const double c17 = 6.6555897559870685;
  const double c63 = 6.49519052838329;
  const double c33 = 6.2749501990055663;
  const double c1 = 6.0756944047573693;
  const double c52 = 5.8094750193111251;
  const double c38 = 5.123475382979799;
  const double c44 = 5.0311529493745271;
  const double c65 = 5;
  const double c21 = 4.4370598373247123;
  const double c58 = 4.3301270189221936;
  const double c30 = 4.1833001326703778;
  const double c50 = 3.872983346207417;
  const double c12 = 3.8426065372348495;
  const double c69 = 3.75;
  const double c22 = 3.6228441865473595;
  const double c8 = 3.5078038001005702;
  const double c43 = 3.3541019662496847;
  const double c60 = 3.247595264191645;
  const double c2 = 3.0378472023786847;
  const double c53 = 2.9047375096555625;
  const double c29 = 2.7171331399105196;
  const double c39 = 2.5617376914898995;
  const double c66 = 2.5;
  const double c14 = 2.2185299186623562;
  const double c31 = 2.0916500663351889;
  const double c51 = 1.9364916731037085;
  const double c10 = 1.9213032686174247;
  const double c67 = 1.875;
  const double c9 = 1.7539019000502851;
  const double c61 = 1.7320508075688772;
  const double c49 = 1.6770509831248424;
  const double c59 = 1.6237976320958225;
  const double c36 = 1.5687375497513916;
  const double c25 = 1.3585665699552598;
  const double c40 = 1.2808688457449497;
  const double c3 = 1.2151388809514738;
  const double c15 = 1.1092649593311781;
  const double c56 = 0.96824583655185426;
  const double c68 = 0.9375;
  const double c28 = 0.90571104663683988;
  const double c57 = 0.8660254037844386;
  const double c46 = 0.83852549156242118;
  const double c37 = 0.78436877487569578;
  const double c5 = 0.70156076002011403;
  const double c0 = 0.60756944047573691;
  const double c34 = 0.52291251658379723;
  const double c64 = 0.5;
  const double c54 = 0.48412291827592713;
  const double c24 = 0.45285552331841994;
  const double c45 = 0.41926274578121059;
  const double c6 = 0.35078038001005701;
  const double c35 = 0.26145625829189861;
  const double c55 = 0.24206145913796356;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 55, source += 126) {
    target[0] =  c0 * source[0] - c0 * source[2] - c1 * source[12]
                  + c1 * source[14] + c2 * source[24] - c2 * source[26];
    target[1] =  c3 * source[1] - c4 * source[13] + c1 * source[25];
    target[2] =  c3 * source[3] - c4 * source[15] + c1 * source[27];
    target[3] =  c3 * source[4] - c4 * source[16] + c1 * source[28];
    target[4] =  c5 * source[5] - c6 * source[0] - c6 * source[2]
                  - c7 * source[17] + c8 * source[12] + c8 * source[14]
                  + c8 * source[29] - c9 * source[24] - c9 * source[26];
    target[5] =  c2 * source[6] - c2 * source[8] - c1 * source[18]
                  + c1 * source[20] + c0 * source[30] - c0 * source[32];
    target[6] =  c1 * source[7] - c4 * source[19] + c3 * source[31];
    target[7] =  c1 * source[9] - c4 * source[21] + c3 * source[33];
    target[8] =  c1 * source[10] - c4 * source[22] + c3 * source[34];
    target[9] =  c8 * source[11] - c9 * source[6] - c9 * source[8]
                  - c7 * source[23] + c8 * source[18] + c8 * source[20]
                  + c5 * source[35] - c6 * source[30] - c6 * source[32];
    target[10] =  c10 * source[36] - c10 * source[38] - c11 * source[48]
                  + c11 * source[50] + c10 * source[60] - c10 * source[62];
    target[11] =  c12 * source[37] - c13 * source[49] + c12 * source[61];
    target[12] =  c12 * source[39] - c13 * source[51] + c12 * source[63];
    target[13] =  c12 * source[40] - c13 * source[52] + c12 * source[64];
    target[14] =  c14 * source[41] - c15 * source[36] - c15 * source[38]
                  - c16 * source[53] + c17 * source[48] + c17 * source[50]
                  + c14 * source[65] - c15 * source[60] - c15 * source[62];
    target[15] =  c18 * source[42] - c18 * source[44] - c18 * source[54]
                  + c18 * source[56];
    target[16] =  c19 * source[43] - c19 * source[55];
    target[17] =  c19 * source[45] - c19 * source[57];
    target[18] =  c19 * source[46] - c19 * source[58];
    target[19] =  c20 * source[47] - c21 * source[42] - c21 * source[44]
                  - c20 * source[59] + c21 * source[54] + c21 * source[56];
    target[20] =  c22 * source[66] - c22 * source[68] - c23 * source[78]
                  + c23 * source[80] - c24 * source[0] + c24 * source[2]
                  + c25 * source[12] - c25 * source[14] - c24 * source[12]
                  + c24 * source[14] + c25 * source[24] - c25 * source[26];
    target[21] =  c26 * source[67] - c27 * source[79] - c28 * source[1]
                  + c29 * source[13] - c28 * source[13] + c29 * source[25];
    target[22] =  c26 * source[69] - c27 * source[81] - c28 * source[3]
                  + c29 * source[15] - c28 * source[15] + c29 * source[27];
    target[23] =  c26 * source[70] - c27 * source[82] - c28 * source[4]
                  + c29 * source[16] - c28 * source[16] + c29 * source[28];
    target[24] =  c30 * source[71] - c31 * source[66] - c31 * source[68]
                  - c32 * source[83] + c33 * source[78] + c33 * source[80]
                  - c34 * source[5] + c35 * source[0] + c35 * source[2]
                  + c36 * source[17] - c37 * source[12] - c37 * source[14]
                  - c34 * source[17] + c35 * source[12] + c35 * source[14]
                  + c36 * source[29] - c37 * source[24] - c37 * source[26];
    target[25] =  c23 * source[72] - c23 * source[74] - c22 * source[84]
                  + c22 * source[86] - c25 * source[6] + c25 * source[8]
                  + c24 * source[18] - c24 * source[20] - c25 * source[18]
                  + c25 * source[20] + c24 * source[30] - c24 * source[32];
    target[26] =  c27 * source[73] - c26 * source[85] - c29 * source[7]
                  + c28 * source[19] - c29 * source[19] + c28 * source[31];
    target[27] =  c27 * source[75] - c26 * source[87] - c29 * source[9]
                  + c28 * source[21] - c29 * source[21] + c28 * source[33];
    target[28] =  c27 * source[76] - c26 * source[88] - c29 * source[10]
                  + c28 * source[22] - c29 * source[22] + c28 * source[34];
    target[29] =  c32 * source[77] - c33 * source[72] - c33 * source[74]
                  - c30 * source[89] + c31 * source[84] + c31 * source[86]
                  - c36 * source[11] + c37 * source[6] + c37 * source[8]
                  + c34 * source[23] - c35 * source[18] - c35 * source[20]
                  - c36 * source[23] + c37 * source[18] + c37 * source[20]
                  + c34 * source[35] - c35 * source[30] - c35 * source[32];
    target[30] =  c21 * source[90] - c21 * source[92] - c21 * source[102]
                  + c21 * source[104] - c14 * source[36] + c14 * source[38]
                  + c14 * source[48] - c14 * source[50] - c14 * source[48]
                  + c14 * source[50] + c14 * source[60] - c14 * source[62];
    target[31] =  c20 * source[91] - c20 * source[103] - c21 * source[37]
                  + c21 * source[49] - c21 * source[49] + c21 * source[61];
    target[32] =  c20 * source[93] - c20 * source[105] - c21 * source[39]
                  + c21 * source[51] - c21 * source[51] + c21 * source[63];
    target[33] =  c20 * source[94] - c20 * source[106] - c21 * source[40]
                  + c21 * source[52] - c21 * source[52] + c21 * source[64];
    target[34] =  c38 * source[95] - c39 * source[90] - c39 * source[92]
                  - c38 * source[107] + c39 * source[102] + c39 * source[104]
                  - c39 * source[41] + c40 * source[36] + c40 * source[38]
                  + c39 * source[53] - c40 * source[48] - c40 * source[50]
                  - c39 * source[53] + c40 * source[48] + c40 * source[50]
                  + c39 * source[65] - c40 * source[60] - c40 * source[62];
    target[35] =  c20 * source[96] - c20 * source[98] - c21 * source[42]
                  + c21 * source[44] - c21 * source[54] + c21 * source[56];
    target[36] =  c41 * source[97] - c20 * source[43] - c20 * source[55];
    target[37] =  c41 * source[99] - c20 * source[45] - c20 * source[57];
    target[38] =  c41 * source[100] - c20 * source[46] - c20 * source[58];
    target[39] =  c42 * source[101] - c38 * source[96] - c38 * source[98]
                  - c38 * source[47] + c39 * source[42] + c39 * source[44]
                  - c38 * source[59] + c39 * source[54] + c39 * source[56];
    target[40] =  c43 * source[108] - c43 * source[110] - c44 * source[66]
                  + c44 * source[68] - c44 * source[78] + c44 * source[80]
                  + c45 * source[0] - c45 * source[2] + c46 * source[12]
                  - c46 * source[14] + c45 * source[24] - c45 * source[26];
    target[41] =  c47 * source[109] - c48 * source[67] - c48 * source[79]
                  + c46 * source[1] + c49 * source[13] + c46 * source[25];
    target[42] =  c47 * source[111] - c48 * source[69] - c48 * source[81]
                  + c46 * source[3] + c49 * source[15] + c46 * source[27];
    target[43] =  c47 * source[112] - c48 * source[70] - c48 * source[82]
                  + c46 * source[4] + c49 * source[16] + c46 * source[28];
    target[44] =  c50 * source[113] - c51 * source[108] - c51 * source[110]
                  - c52 * source[71] + c53 * source[66] + c53 * source[68]
                  - c52 * source[83] + c53 * source[78] + c53 * source[80]
                  + c54 * source[5] - c55 * source[0] - c55 * source[2]
                  + c56 * source[17] - c54 * source[12] - c54 * source[14]
                  + c54 * source[29] - c55 * source[24] - c55 * source[26];
    target[45] =  c43 * source[114] - c43 * source[116] - c44 * source[72]
                  + c44 * source[74] - c44 * source[84] + c44 * source[86]
                  + c45 * source[6] - c45 * source[8] + c46 * source[18]
                  - c46 * source[20] + c45 * source[30] - c45 * source[32];
    target[46] =  c47 * source[115] - c48 * source[73] - c48 * source[85]
                  + c46 * source[7] + c49 * source[19] + c46 * source[31];
    target[47] =  c47 * source[117] - c48 * source[75] - c48 * source[87]
                  + c46 * source[9] + c49 * source[21] + c46 * source[33];
    target[48] =  c47 * source[118] - c48 * source[76] - c48 * source[88]
                  + c46 * source[10] + c49 * source[22] + c46 * source[34];
    target[49] =  c50 * source[119] - c51 * source[114] - c51 * source[116]
                  - c52 * source[77] + c53 * source[72] + c53 * source[74]
                  - c52 * source[89] + c53 * source[84] + c53 * source[86]
                  + c54 * source[11] - c55 * source[6] - c55 * source[8]
                  + c56 * source[23] - c54 * source[18] - c54 * source[20]
                  + c54 * source[35] - c55 * source[30] - c55 * source[32];
    target[50] =  c57 * source[120] - c57 * source[122] - c58 * source[90]
                  + c58 * source[92] - c58 * source[102] + c58 * source[104]
                  + c59 * source[36] - c59 * source[38] + c60 * source[48]
                  - c60 * source[50] + c59 * source[60] - c59 * source[62];
    target[51] =  c61 * source[121] - c62 * source[91] - c62 * source[103]
                  + c60 * source[37] + c63 * source[49] + c60 * source[61];
    target[52] =  c61 * source[123] - c62 * source[93] - c62 * source[105]
                  + c60 * source[39] + c63 * source[51] + c60 * source[63];
    target[53] =  c61 * source[124] - c62 * source[94] - c62 * source[106]
                  + c60 * source[40] + c63 * source[52] + c60 * source[64];
    target[54] =  source[125] - c64 * source[120] - c64 * source[122]
                  - c65 * source[95] + c66 * source[90] + c66 * source[92]
                  - c65 * source[107] + c66 * source[102] + c66 * source[104]
                  + c67 * source[41] - c68 * source[36] - c68 * source[38]
                  + c69 * source[53] - c67 * source[48] - c67 * source[50]
                  + c67 * source[65] - c68 * source[60] - c68 * source[62];
  }
}

void CCarSphList::carsph_52(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c13 = 23.055639223409095;
  const double c27 = 21.737065119284157;
  const double c41 = 17.748239349298849;
  const double c19 = 15.370426148939398;
  const double c16 = 13.311179511974137;
  const double c32 = 12.549900398011133;
  const double c4 = 12.151388809514739;
  const double c11 = 11.527819611704547;
  const double c23 = 10.868532559642079;
  const double c42 = 10.246950765959598;
  const double c48 = 10.062305898749054;
  const double c20 = 8.8741196746494246;
  const double c62 = 8.6602540378443873;
  const double c18 = 7.6852130744696989;
  const double c26 = 7.245688373094719;
  const double c7 = 7.0156076002011405;
  const double c47 = 6.7082039324993694;
  const double c17 = 6.6555897559870685;
  const double c63 = 6.49519052838329;
  const double c33 = 6.2749501990055663;
  const double c1 = 6.0756944047573693;
  const double c52 = 5.8094750193111251;
  const double c38 = 5.123475382979799;
  const double c44 = 5.0311529493745271;
  const double c65 = 5;
  const double c21 = 4.4370598373247123;
  const double c58 = 4.3301270189221936;
  const double c30 = 4.1833001326703778;
  const double c50 = 3.872983346207417;
  const double c12 = 3.8426065372348495;
  const double c69 = 3.75;
  const double c22 = 3.6228441865473595;
  const double c8 = 3.5078038001005702;
  const double c43 = 3.3541019662496847;
  const double c60 = 3.247595264191645;
  const double c2 = 3.0378472023786847;
  const double c53 = 2.9047375096555625;
  const double c29 = 2.7171331399105196;
  const double c39 = 2.5617376914898995;
  const double c66 = 2.5;
  const double c14 = 2.2185299186623562;
  const double c31 = 2.0916500663351889;
  const double c51 = 1.9364916731037085;
  const double c10 = 1.9213032686174247;
  const double c67 = 1.875;
  const double c9 = 1.7539019000502851;
  const double c61 = 1.7320508075688772;
  const double c49 = 1.6770509831248424;
  const double c59 = 1.6237976320958225;
  const double c36 = 1.5687375497513916;
  const double c25 = 1.3585665699552598;
  const double c40 = 1.2808688457449497;
  const double c3 = 1.2151388809514738;
  const double c15 = 1.1092649593311781;
  const double c56 = 0.96824583655185426;
  const double c68 = 0.9375;
  const double c28 = 0.90571104663683988;
  const double c57 = 0.8660254037844386;
  const double c46 = 0.83852549156242118;
  const double c37 = 0.78436877487569578;
  const double c5 = 0.70156076002011403;
  const double c0 = 0.60756944047573691;
  const double c34 = 0.52291251658379723;
  const double c64 = 0.5;
  const double c54 = 0.48412291827592713;
  const double c24 = 0.45285552331841994;
  const double c45 = 0.41926274578121059;
  const double c6 = 0.35078038001005701;
  const double c35 = 0.26145625829189861;
  const double c55 = 0.24206145913796356;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 55, source += 126) {
    target[0] =  c0 * source[0] - c0 * source[2] - c1 * source[12]
                  + c1 * source[14] + c2 * source[24] - c2 * source[26];
    target[1] =  c3 * source[1] - c4 * source[13] + c1 * source[25];
    target[2] =  c3 * source[3] - c4 * source[15] + c1 * source[27];
    target[3] =  c3 * source[4] - c4 * source[16] + c1 * source[28];
    target[4] =  c5 * source[5] - c6 * source[0] - c6 * source[2]
                  - c7 * source[17] + c8 * source[12] + c8 * source[14]
                  + c8 * source[29] - c9 * source[24] - c9 * source[26];
    target[5] =  c2 * source[6] - c2 * source[8] - c1 * source[18]
                  + c1 * source[20] + c0 * source[30] - c0 * source[32];
    target[6] =  c1 * source[7] - c4 * source[19] + c3 * source[31];
    target[7] =  c1 * source[9] - c4 * source[21] + c3 * source[33];
    target[8] =  c1 * source[10] - c4 * source[22] + c3 * source[34];
    target[9] =  c8 * source[11] - c9 * source[6] - c9 * source[8]
                  - c7 * source[23] + c8 * source[18] + c8 * source[20]
                  + c5 * source[35] - c6 * source[30] - c6 * source[32];
    target[10] =  c10 * source[36] - c10 * source[38] - c11 * source[48]
                  + c11 * source[50] + c10 * source[60] - c10 * source[62];
    target[11] =  c12 * source[37] - c13 * source[49] + c12 * source[61];
    target[12] =  c12 * source[39] - c13 * source[51] + c12 * source[63];
    target[13] =  c12 * source[40] - c13 * source[52] + c12 * source[64];
    target[14] =  c14 * source[41] - c15 * source[36] - c15 * source[38]
                  - c16 * source[53] + c17 * source[48] + c17 * source[50]
                  + c14 * source[65] - c15 * source[60] - c15 * source[62];
    target[15] =  c18 * source[42] - c18 * source[44] - c18 * source[54]
                  + c18 * source[56];
    target[16] =  c19 * source[43] - c19 * source[55];
    target[17] =  c19 * source[45] - c19 * source[57];
    target[18] =  c19 * source[46] - c19 * source[58];
    target[19] =  c20 * source[47] - c21 * source[42] - c21 * source[44]
                  - c20 * source[59] + c21 * source[54] + c21 * source[56];
    target[20] =  c22 * source[66] - c22 * source[68] - c23 * source[78]
                  + c23 * source[80] - c24 * source[0] + c24 * source[2]
                  + c25 * source[12] - c25 * source[14] - c24 * source[12]
                  + c24 * source[14] + c25 * source[24] - c25 * source[26];
    target[21] =  c26 * source[67] - c27 * source[79] - c28 * source[1]
                  + c29 * source[13] - c28 * source[13] + c29 * source[25];
    target[22] =  c26 * source[69] - c27 * source[81] - c28 * source[3]
                  + c29 * source[15] - c28 * source[15] + c29 * source[27];
    target[23] =  c26 * source[70] - c27 * source[82] - c28 * source[4]
                  + c29 * source[16] - c28 * source[16] + c29 * source[28];
    target[24] =  c30 * source[71] - c31 * source[66] - c31 * source[68]
                  - c32 * source[83] + c33 * source[78] + c33 * source[80]
                  - c34 * source[5] + c35 * source[0] + c35 * source[2]
                  + c36 * source[17] - c37 * source[12] - c37 * source[14]
                  - c34 * source[17] + c35 * source[12] + c35 * source[14]
                  + c36 * source[29] - c37 * source[24] - c37 * source[26];
    target[25] =  c23 * source[72] - c23 * source[74] - c22 * source[84]
                  + c22 * source[86] - c25 * source[6] + c25 * source[8]
                  + c24 * source[18] - c24 * source[20] - c25 * source[18]
                  + c25 * source[20] + c24 * source[30] - c24 * source[32];
    target[26] =  c27 * source[73] - c26 * source[85] - c29 * source[7]
                  + c28 * source[19] - c29 * source[19] + c28 * source[31];
    target[27] =  c27 * source[75] - c26 * source[87] - c29 * source[9]
                  + c28 * source[21] - c29 * source[21] + c28 * source[33];
    target[28] =  c27 * source[76] - c26 * source[88] - c29 * source[10]
                  + c28 * source[22] - c29 * source[22] + c28 * source[34];
    target[29] =  c32 * source[77] - c33 * source[72] - c33 * source[74]
                  - c30 * source[89] + c31 * source[84] + c31 * source[86]
                  - c36 * source[11] + c37 * source[6] + c37 * source[8]
                  + c34 * source[23] - c35 * source[18] - c35 * source[20]
                  - c36 * source[23] + c37 * source[18] + c37 * source[20]
                  + c34 * source[35] - c35 * source[30] - c35 * source[32];
    target[30] =  c21 * source[90] - c21 * source[92] - c21 * source[102]
                  + c21 * source[104] - c14 * source[36] + c14 * source[38]
                  + c14 * source[48] - c14 * source[50] - c14 * source[48]
                  + c14 * source[50] + c14 * source[60] - c14 * source[62];
    target[31] =  c20 * source[91] - c20 * source[103] - c21 * source[37]
                  + c21 * source[49] - c21 * source[49] + c21 * source[61];
    target[32] =  c20 * source[93] - c20 * source[105] - c21 * source[39]
                  + c21 * source[51] - c21 * source[51] + c21 * source[63];
    target[33] =  c20 * source[94] - c20 * source[106] - c21 * source[40]
                  + c21 * source[52] - c21 * source[52] + c21 * source[64];
    target[34] =  c38 * source[95] - c39 * source[90] - c39 * source[92]
                  - c38 * source[107] + c39 * source[102] + c39 * source[104]
                  - c39 * source[41] + c40 * source[36] + c40 * source[38]
                  + c39 * source[53] - c40 * source[48] - c40 * source[50]
                  - c39 * source[53] + c40 * source[48] + c40 * source[50]
                  + c39 * source[65] - c40 * source[60] - c40 * source[62];
    target[35] =  c20 * source[96] - c20 * source[98] - c21 * source[42]
                  + c21 * source[44] - c21 * source[54] + c21 * source[56];
    target[36] =  c41 * source[97] - c20 * source[43] - c20 * source[55];
    target[37] =  c41 * source[99] - c20 * source[45] - c20 * source[57];
    target[38] =  c41 * source[100] - c20 * source[46] - c20 * source[58];
    target[39] =  c42 * source[101] - c38 * source[96] - c38 * source[98]
                  - c38 * source[47] + c39 * source[42] + c39 * source[44]
                  - c38 * source[59] + c39 * source[54] + c39 * source[56];
    target[40] =  c43 * source[108] - c43 * source[110] - c44 * source[66]
                  + c44 * source[68] - c44 * source[78] + c44 * source[80]
                  + c45 * source[0] - c45 * source[2] + c46 * source[12]
                  - c46 * source[14] + c45 * source[24] - c45 * source[26];
    target[41] =  c47 * source[109] - c48 * source[67] - c48 * source[79]
                  + c46 * source[1] + c49 * source[13] + c46 * source[25];
    target[42] =  c47 * source[111] - c48 * source[69] - c48 * source[81]
                  + c46 * source[3] + c49 * source[15] + c46 * source[27];
    target[43] =  c47 * source[112] - c48 * source[70] - c48 * source[82]
                  + c46 * source[4] + c49 * source[16] + c46 * source[28];
    target[44] =  c50 * source[113] - c51 * source[108] - c51 * source[110]
                  - c52 * source[71] + c53 * source[66] + c53 * source[68]
                  - c52 * source[83] + c53 * source[78] + c53 * source[80]
                  + c54 * source[5] - c55 * source[0] - c55 * source[2]
                  + c56 * source[17] - c54 * source[12] - c54 * source[14]
                  + c54 * source[29] - c55 * source[24] - c55 * source[26];
    target[45] =  c43 * source[114] - c43 * source[116] - c44 * source[72]
                  + c44 * source[74] - c44 * source[84] + c44 * source[86]
                  + c45 * source[6] - c45 * source[8] + c46 * source[18]
                  - c46 * source[20] + c45 * source[30] - c45 * source[32];
    target[46] =  c47 * source[115] - c48 * source[73] - c48 * source[85]
                  + c46 * source[7] + c49 * source[19] + c46 * source[31];
    target[47] =  c47 * source[117] - c48 * source[75] - c48 * source[87]
                  + c46 * source[9] + c49 * source[21] + c46 * source[33];
    target[48] =  c47 * source[118] - c48 * source[76] - c48 * source[88]
                  + c46 * source[10] + c49 * source[22] + c46 * source[34];
    target[49] =  c50 * source[119] - c51 * source[114] - c51 * source[116]
                  - c52 * source[77] + c53 * source[72] + c53 * source[74]
                  - c52 * source[89] + c53 * source[84] + c53 * source[86]
                  + c54 * source[11] - c55 * source[6] - c55 * source[8]
                  + c56 * source[23] - c54 * source[18] - c54 * source[20]
                  + c54 * source[35] - c55 * source[30] - c55 * source[32];
    target[50] =  c57 * source[120] - c57 * source[122] - c58 * source[90]
                  + c58 * source[92] - c58 * source[102] + c58 * source[104]
                  + c59 * source[36] - c59 * source[38] + c60 * source[48]
                  - c60 * source[50] + c59 * source[60] - c59 * source[62];
    target[51] =  c61 * source[121] - c62 * source[91] - c62 * source[103]
                  + c60 * source[37] + c63 * source[49] + c60 * source[61];
    target[52] =  c61 * source[123] - c62 * source[93] - c62 * source[105]
                  + c60 * source[39] + c63 * source[51] + c60 * source[63];
    target[53] =  c61 * source[124] - c62 * source[94] - c62 * source[106]
                  + c60 * source[40] + c63 * source[52] + c60 * source[64];
    target[54] =  source[125] - c64 * source[120] - c64 * source[122]
                  - c65 * source[95] + c66 * source[90] + c66 * source[92]
                  - c65 * source[107] + c66 * source[102] + c66 * source[104]
                  + c67 * source[41] - c68 * source[36] - c68 * source[38]
                  + c69 * source[53] - c67 * source[48] - c67 * source[50]
                  + c67 * source[65] - c68 * source[60] - c68 * source[62];
  }
}

