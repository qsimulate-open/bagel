//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_42.cc
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


void CarSphList::carsph_42(const int nloop, const double* source, double* target) {
  const double c28 = 11.61895003862225;
  const double c15 = 10.868532559642079;
  const double c3 = 7.6852130744696989;
  const double c30 = 6.7082039324993694;
  const double c18 = 6.2749501990055663;
  const double c22 = 5.8094750193111251;
  const double c34 = 5.4772255750516612;
  const double c13 = 5.4342662798210393;
  const double c45 = 5.196152422706632;
  const double c9 = 5.123475382979799;
  const double c6 = 4.4370598373247123;
  const double c35 = 4.1079191812887457;
  const double c1 = 3.8426065372348495;
  const double c14 = 3.6228441865473595;
  const double c24 = 3.3541019662496847;
  const double c36 = 3.1622776601683795;
  const double c19 = 3.1374750995027831;
  const double c48 = 3;
  const double c10 = 2.9580398915498081;
  const double c20 = 2.9047375096555625;
  const double c32 = 2.7386127875258306;
  const double c41 = 2.598076211353316;
  const double c8 = 2.5617376914898995;
  const double c38 = 2.3717082451262845;
  const double c7 = 2.2185299186623562;
  const double c16 = 2.0916500663351889;
  const double c33 = 2.0539595906443728;
  const double c29 = 1.9364916731037085;
  const double c12 = 1.8114220932736798;
  const double c44 = 1.7320508075688772;
  const double c25 = 1.6770509831248424;
  const double c37 = 1.5811388300841898;
  const double c49 = 1.5;
  const double c11 = 1.479019945774904;
  const double c46 = 1.299038105676658;
  const double c2 = 1.2808688457449497;
  const double c39 = 1.1858541225631423;
  const double c31 = 1.1180339887498949;
  const double c17 = 1.0458250331675945;
  const double c23 = 0.96824583655185426;
  const double c40 = 0.8660254037844386;
  const double c52 = 0.75;
  const double c4 = 0.73950997288745202;
  const double c43 = 0.649519052838329;
  const double c0 = 0.64043442287247487;
  const double c26 = 0.55901699437494745;
  const double c47 = 0.5;
  const double c21 = 0.48412291827592713;
  const double c50 = 0.375;
  const double c5 = 0.36975498644372601;
  const double c42 = 0.3247595264191645;
  const double c27 = 0.27950849718747373;
  const double c51 = 0.1875;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 45, source += 90) {
    target[0] =  c0 * source[0] - c0 * source[2] - c1 * source[12]
                  + c1 * source[14] + c0 * source[24] - c0 * source[26];
    target[1] =  c2 * source[1] - c3 * source[13] + c2 * source[25];
    target[2] =  c2 * source[3] - c3 * source[15] + c2 * source[27];
    target[3] =  c2 * source[4] - c3 * source[16] + c2 * source[28];
    target[4] =  c4 * source[5] - c5 * source[0] - c5 * source[2]
                  - c6 * source[17] + c7 * source[12] + c7 * source[14]
                  + c4 * source[29] - c5 * source[24] - c5 * source[26];
    target[5] =  c8 * source[6] - c8 * source[8] - c8 * source[18]
                  + c8 * source[20];
    target[6] =  c9 * source[7] - c9 * source[19];
    target[7] =  c9 * source[9] - c9 * source[21];
    target[8] =  c9 * source[10] - c9 * source[22];
    target[9] =  c10 * source[11] - c11 * source[6] - c11 * source[8]
                  - c10 * source[23] + c11 * source[18] + c11 * source[20];
    target[10] =  c12 * source[30] - c12 * source[32] - c13 * source[42]
                  + c13 * source[44];
    target[11] =  c14 * source[31] - c15 * source[43];
    target[12] =  c14 * source[33] - c15 * source[45];
    target[13] =  c14 * source[34] - c15 * source[46];
    target[14] =  c16 * source[35] - c17 * source[30] - c17 * source[32]
                  - c18 * source[47] + c19 * source[42] + c19 * source[44];
    target[15] =  c13 * source[36] - c13 * source[38] - c12 * source[48]
                  + c12 * source[50];
    target[16] =  c15 * source[37] - c14 * source[49];
    target[17] =  c15 * source[39] - c14 * source[51];
    target[18] =  c15 * source[40] - c14 * source[52];
    target[19] =  c18 * source[41] - c19 * source[36] - c19 * source[38]
                  - c16 * source[53] + c17 * source[48] + c17 * source[50];
    target[20] =  c20 * source[54] - c20 * source[56] - c20 * source[66]
                  + c20 * source[68] - c21 * source[0] + c21 * source[2]
                  + c21 * source[12] - c21 * source[14] - c21 * source[12]
                  + c21 * source[14] + c21 * source[24] - c21 * source[26];
    target[21] =  c22 * source[55] - c22 * source[67] - c23 * source[1]
                  + c23 * source[13] - c23 * source[13] + c23 * source[25];
    target[22] =  c22 * source[57] - c22 * source[69] - c23 * source[3]
                  + c23 * source[15] - c23 * source[15] + c23 * source[27];
    target[23] =  c22 * source[58] - c22 * source[70] - c23 * source[4]
                  + c23 * source[16] - c23 * source[16] + c23 * source[28];
    target[24] =  c24 * source[59] - c25 * source[54] - c25 * source[56]
                  - c24 * source[71] + c25 * source[66] + c25 * source[68]
                  - c26 * source[5] + c27 * source[0] + c27 * source[2]
                  + c26 * source[17] - c27 * source[12] - c27 * source[14]
                  - c26 * source[17] + c27 * source[12] + c27 * source[14]
                  + c26 * source[29] - c27 * source[24] - c27 * source[26];
    target[25] =  c22 * source[60] - c22 * source[62] - c23 * source[6]
                  + c23 * source[8] - c23 * source[18] + c23 * source[20];
    target[26] =  c28 * source[61] - c29 * source[7] - c29 * source[19];
    target[27] =  c28 * source[63] - c29 * source[9] - c29 * source[21];
    target[28] =  c28 * source[64] - c29 * source[10] - c29 * source[22];
    target[29] =  c30 * source[65] - c24 * source[60] - c24 * source[62]
                  - c31 * source[11] + c26 * source[6] + c26 * source[8]
                  - c31 * source[23] + c26 * source[18] + c26 * source[20];
    target[30] =  c32 * source[72] - c32 * source[74] - c33 * source[30]
                  + c33 * source[32] - c33 * source[42] + c33 * source[44];
    target[31] =  c34 * source[73] - c35 * source[31] - c35 * source[43];
    target[32] =  c34 * source[75] - c35 * source[33] - c35 * source[45];
    target[33] =  c34 * source[76] - c35 * source[34] - c35 * source[46];
    target[34] =  c36 * source[77] - c37 * source[72] - c37 * source[74]
                  - c38 * source[35] + c39 * source[30] + c39 * source[32]
                  - c38 * source[47] + c39 * source[42] + c39 * source[44];
    target[35] =  c32 * source[78] - c32 * source[80] - c33 * source[36]
                  + c33 * source[38] - c33 * source[48] + c33 * source[50];
    target[36] =  c34 * source[79] - c35 * source[37] - c35 * source[49];
    target[37] =  c34 * source[81] - c35 * source[39] - c35 * source[51];
    target[38] =  c34 * source[82] - c35 * source[40] - c35 * source[52];
    target[39] =  c36 * source[83] - c37 * source[78] - c37 * source[80]
                  - c38 * source[41] + c39 * source[36] + c39 * source[38]
                  - c38 * source[53] + c39 * source[48] + c39 * source[50];
    target[40] =  c40 * source[84] - c40 * source[86] - c41 * source[54]
                  + c41 * source[56] - c41 * source[66] + c41 * source[68]
                  + c42 * source[0] - c42 * source[2] + c43 * source[12]
                  - c43 * source[14] + c42 * source[24] - c42 * source[26];
    target[41] =  c44 * source[85] - c45 * source[55] - c45 * source[67]
                  + c43 * source[1] + c46 * source[13] + c43 * source[25];
    target[42] =  c44 * source[87] - c45 * source[57] - c45 * source[69]
                  + c43 * source[3] + c46 * source[15] + c43 * source[27];
    target[43] =  c44 * source[88] - c45 * source[58] - c45 * source[70]
                  + c43 * source[4] + c46 * source[16] + c43 * source[28];
    target[44] =  source[89] - c47 * source[84] - c47 * source[86]
                  - c48 * source[59] + c49 * source[54] + c49 * source[56]
                  - c48 * source[71] + c49 * source[66] + c49 * source[68]
                  + c50 * source[5] - c51 * source[0] - c51 * source[2]
                  + c52 * source[17] - c50 * source[12] - c50 * source[14]
                  + c50 * source[29] - c51 * source[24] - c51 * source[26];
  }
}

void CCarSphList::carsph_42(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c28 = 11.61895003862225;
  const double c15 = 10.868532559642079;
  const double c3 = 7.6852130744696989;
  const double c30 = 6.7082039324993694;
  const double c18 = 6.2749501990055663;
  const double c22 = 5.8094750193111251;
  const double c34 = 5.4772255750516612;
  const double c13 = 5.4342662798210393;
  const double c45 = 5.196152422706632;
  const double c9 = 5.123475382979799;
  const double c6 = 4.4370598373247123;
  const double c35 = 4.1079191812887457;
  const double c1 = 3.8426065372348495;
  const double c14 = 3.6228441865473595;
  const double c24 = 3.3541019662496847;
  const double c36 = 3.1622776601683795;
  const double c19 = 3.1374750995027831;
  const double c48 = 3;
  const double c10 = 2.9580398915498081;
  const double c20 = 2.9047375096555625;
  const double c32 = 2.7386127875258306;
  const double c41 = 2.598076211353316;
  const double c8 = 2.5617376914898995;
  const double c38 = 2.3717082451262845;
  const double c7 = 2.2185299186623562;
  const double c16 = 2.0916500663351889;
  const double c33 = 2.0539595906443728;
  const double c29 = 1.9364916731037085;
  const double c12 = 1.8114220932736798;
  const double c44 = 1.7320508075688772;
  const double c25 = 1.6770509831248424;
  const double c37 = 1.5811388300841898;
  const double c49 = 1.5;
  const double c11 = 1.479019945774904;
  const double c46 = 1.299038105676658;
  const double c2 = 1.2808688457449497;
  const double c39 = 1.1858541225631423;
  const double c31 = 1.1180339887498949;
  const double c17 = 1.0458250331675945;
  const double c23 = 0.96824583655185426;
  const double c40 = 0.8660254037844386;
  const double c52 = 0.75;
  const double c4 = 0.73950997288745202;
  const double c43 = 0.649519052838329;
  const double c0 = 0.64043442287247487;
  const double c26 = 0.55901699437494745;
  const double c47 = 0.5;
  const double c21 = 0.48412291827592713;
  const double c50 = 0.375;
  const double c5 = 0.36975498644372601;
  const double c42 = 0.3247595264191645;
  const double c27 = 0.27950849718747373;
  const double c51 = 0.1875;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 45, source += 90) {
    target[0] =  c0 * source[0] - c0 * source[2] - c1 * source[12]
                  + c1 * source[14] + c0 * source[24] - c0 * source[26];
    target[1] =  c2 * source[1] - c3 * source[13] + c2 * source[25];
    target[2] =  c2 * source[3] - c3 * source[15] + c2 * source[27];
    target[3] =  c2 * source[4] - c3 * source[16] + c2 * source[28];
    target[4] =  c4 * source[5] - c5 * source[0] - c5 * source[2]
                  - c6 * source[17] + c7 * source[12] + c7 * source[14]
                  + c4 * source[29] - c5 * source[24] - c5 * source[26];
    target[5] =  c8 * source[6] - c8 * source[8] - c8 * source[18]
                  + c8 * source[20];
    target[6] =  c9 * source[7] - c9 * source[19];
    target[7] =  c9 * source[9] - c9 * source[21];
    target[8] =  c9 * source[10] - c9 * source[22];
    target[9] =  c10 * source[11] - c11 * source[6] - c11 * source[8]
                  - c10 * source[23] + c11 * source[18] + c11 * source[20];
    target[10] =  c12 * source[30] - c12 * source[32] - c13 * source[42]
                  + c13 * source[44];
    target[11] =  c14 * source[31] - c15 * source[43];
    target[12] =  c14 * source[33] - c15 * source[45];
    target[13] =  c14 * source[34] - c15 * source[46];
    target[14] =  c16 * source[35] - c17 * source[30] - c17 * source[32]
                  - c18 * source[47] + c19 * source[42] + c19 * source[44];
    target[15] =  c13 * source[36] - c13 * source[38] - c12 * source[48]
                  + c12 * source[50];
    target[16] =  c15 * source[37] - c14 * source[49];
    target[17] =  c15 * source[39] - c14 * source[51];
    target[18] =  c15 * source[40] - c14 * source[52];
    target[19] =  c18 * source[41] - c19 * source[36] - c19 * source[38]
                  - c16 * source[53] + c17 * source[48] + c17 * source[50];
    target[20] =  c20 * source[54] - c20 * source[56] - c20 * source[66]
                  + c20 * source[68] - c21 * source[0] + c21 * source[2]
                  + c21 * source[12] - c21 * source[14] - c21 * source[12]
                  + c21 * source[14] + c21 * source[24] - c21 * source[26];
    target[21] =  c22 * source[55] - c22 * source[67] - c23 * source[1]
                  + c23 * source[13] - c23 * source[13] + c23 * source[25];
    target[22] =  c22 * source[57] - c22 * source[69] - c23 * source[3]
                  + c23 * source[15] - c23 * source[15] + c23 * source[27];
    target[23] =  c22 * source[58] - c22 * source[70] - c23 * source[4]
                  + c23 * source[16] - c23 * source[16] + c23 * source[28];
    target[24] =  c24 * source[59] - c25 * source[54] - c25 * source[56]
                  - c24 * source[71] + c25 * source[66] + c25 * source[68]
                  - c26 * source[5] + c27 * source[0] + c27 * source[2]
                  + c26 * source[17] - c27 * source[12] - c27 * source[14]
                  - c26 * source[17] + c27 * source[12] + c27 * source[14]
                  + c26 * source[29] - c27 * source[24] - c27 * source[26];
    target[25] =  c22 * source[60] - c22 * source[62] - c23 * source[6]
                  + c23 * source[8] - c23 * source[18] + c23 * source[20];
    target[26] =  c28 * source[61] - c29 * source[7] - c29 * source[19];
    target[27] =  c28 * source[63] - c29 * source[9] - c29 * source[21];
    target[28] =  c28 * source[64] - c29 * source[10] - c29 * source[22];
    target[29] =  c30 * source[65] - c24 * source[60] - c24 * source[62]
                  - c31 * source[11] + c26 * source[6] + c26 * source[8]
                  - c31 * source[23] + c26 * source[18] + c26 * source[20];
    target[30] =  c32 * source[72] - c32 * source[74] - c33 * source[30]
                  + c33 * source[32] - c33 * source[42] + c33 * source[44];
    target[31] =  c34 * source[73] - c35 * source[31] - c35 * source[43];
    target[32] =  c34 * source[75] - c35 * source[33] - c35 * source[45];
    target[33] =  c34 * source[76] - c35 * source[34] - c35 * source[46];
    target[34] =  c36 * source[77] - c37 * source[72] - c37 * source[74]
                  - c38 * source[35] + c39 * source[30] + c39 * source[32]
                  - c38 * source[47] + c39 * source[42] + c39 * source[44];
    target[35] =  c32 * source[78] - c32 * source[80] - c33 * source[36]
                  + c33 * source[38] - c33 * source[48] + c33 * source[50];
    target[36] =  c34 * source[79] - c35 * source[37] - c35 * source[49];
    target[37] =  c34 * source[81] - c35 * source[39] - c35 * source[51];
    target[38] =  c34 * source[82] - c35 * source[40] - c35 * source[52];
    target[39] =  c36 * source[83] - c37 * source[78] - c37 * source[80]
                  - c38 * source[41] + c39 * source[36] + c39 * source[38]
                  - c38 * source[53] + c39 * source[48] + c39 * source[50];
    target[40] =  c40 * source[84] - c40 * source[86] - c41 * source[54]
                  + c41 * source[56] - c41 * source[66] + c41 * source[68]
                  + c42 * source[0] - c42 * source[2] + c43 * source[12]
                  - c43 * source[14] + c42 * source[24] - c42 * source[26];
    target[41] =  c44 * source[85] - c45 * source[55] - c45 * source[67]
                  + c43 * source[1] + c46 * source[13] + c43 * source[25];
    target[42] =  c44 * source[87] - c45 * source[57] - c45 * source[69]
                  + c43 * source[3] + c46 * source[15] + c43 * source[27];
    target[43] =  c44 * source[88] - c45 * source[58] - c45 * source[70]
                  + c43 * source[4] + c46 * source[16] + c43 * source[28];
    target[44] =  source[89] - c47 * source[84] - c47 * source[86]
                  - c48 * source[59] + c49 * source[54] + c49 * source[56]
                  - c48 * source[71] + c49 * source[66] + c49 * source[68]
                  + c50 * source[5] - c51 * source[0] - c51 * source[2]
                  + c52 * source[17] - c50 * source[12] - c50 * source[14]
                  + c50 * source[29] - c51 * source[24] - c51 * source[26];
  }
}

