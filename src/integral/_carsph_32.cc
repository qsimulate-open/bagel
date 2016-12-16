//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_32.cc
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


void CarSphList::carsph_32(const int nloop, const double* source, double* target) {
  const double c12 = 6.7082039324993694;
  const double c16 = 4.2426406871192848;
  const double c3 = 4.1079191812887457;
  const double c13 = 3.872983346207417;
  const double c9 = 3.3541019662496847;
  const double c25 = 2.598076211353316;
  const double c18 = 2.4494897427831779;
  const double c6 = 2.3717082451262845;
  const double c14 = 2.1213203435596424;
  const double c1 = 2.0539595906443728;
  const double c10 = 1.9364916731037085;
  const double c24 = 1.7320508075688772;
  const double c8 = 1.6770509831248424;
  const double c27 = 1.5;
  const double c2 = 1.3693063937629153;
  const double c23 = 1.299038105676658;
  const double c19 = 1.2247448713915889;
  const double c7 = 1.1858541225631423;
  const double c17 = 1.0606601717798212;
  const double c11 = 0.96824583655185426;
  const double c22 = 0.8660254037844386;
  const double c4 = 0.79056941504209488;
  const double c28 = 0.75;
  const double c0 = 0.68465319688145765;
  const double c20 = 0.61237243569579447;
  const double c15 = 0.5303300858899106;
  const double c26 = 0.5;
  const double c5 = 0.39528470752104744;
  const double c21 = 0.30618621784789724;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 35, source += 60) {
    target[0] =  c0 * source[0] - c0 * source[2] - c1 * source[12]
                  + c1 * source[14];
    target[1] =  c2 * source[1] - c3 * source[13];
    target[2] =  c2 * source[3] - c3 * source[15];
    target[3] =  c2 * source[4] - c3 * source[16];
    target[4] =  c4 * source[5] - c5 * source[0] - c5 * source[2]
                  - c6 * source[17] + c7 * source[12] + c7 * source[14];
    target[5] =  c1 * source[6] - c1 * source[8] - c0 * source[18]
                  + c0 * source[20];
    target[6] =  c3 * source[7] - c2 * source[19];
    target[7] =  c3 * source[9] - c2 * source[21];
    target[8] =  c3 * source[10] - c2 * source[22];
    target[9] =  c6 * source[11] - c7 * source[6] - c7 * source[8]
                  - c4 * source[23] + c5 * source[18] + c5 * source[20];
    target[10] =  c8 * source[24] - c8 * source[26] - c8 * source[36]
                  + c8 * source[38];
    target[11] =  c9 * source[25] - c9 * source[37];
    target[12] =  c9 * source[27] - c9 * source[39];
    target[13] =  c9 * source[28] - c9 * source[40];
    target[14] =  c10 * source[29] - c11 * source[24] - c11 * source[26]
                  - c10 * source[41] + c11 * source[36] + c11 * source[38];
    target[15] =  c9 * source[30] - c9 * source[32];
    target[16] =  c12 * source[31];
    target[17] =  c12 * source[33];
    target[18] =  c12 * source[34];
    target[19] =  c13 * source[35] - c10 * source[30] - c10 * source[32];
    target[20] =  c14 * source[42] - c14 * source[44] - c15 * source[0]
                  + c15 * source[2] - c15 * source[12] + c15 * source[14];
    target[21] =  c16 * source[43] - c17 * source[1] - c17 * source[13];
    target[22] =  c16 * source[45] - c17 * source[3] - c17 * source[15];
    target[23] =  c16 * source[46] - c17 * source[4] - c17 * source[16];
    target[24] =  c18 * source[47] - c19 * source[42] - c19 * source[44]
                  - c20 * source[5] + c21 * source[0] + c21 * source[2]
                  - c20 * source[17] + c21 * source[12] + c21 * source[14];
    target[25] =  c14 * source[48] - c14 * source[50] - c15 * source[6]
                  + c15 * source[8] - c15 * source[18] + c15 * source[20];
    target[26] =  c16 * source[49] - c17 * source[7] - c17 * source[19];
    target[27] =  c16 * source[51] - c17 * source[9] - c17 * source[21];
    target[28] =  c16 * source[52] - c17 * source[10] - c17 * source[22];
    target[29] =  c18 * source[53] - c19 * source[48] - c19 * source[50]
                  - c20 * source[11] + c21 * source[6] + c21 * source[8]
                  - c20 * source[23] + c21 * source[18] + c21 * source[20];
    target[30] =  c22 * source[54] - c22 * source[56] - c23 * source[24]
                  + c23 * source[26] - c23 * source[36] + c23 * source[38];
    target[31] =  c24 * source[55] - c25 * source[25] - c25 * source[37];
    target[32] =  c24 * source[57] - c25 * source[27] - c25 * source[39];
    target[33] =  c24 * source[58] - c25 * source[28] - c25 * source[40];
    target[34] =  source[59] - c26 * source[54] - c26 * source[56]
                  - c27 * source[29] + c28 * source[24] + c28 * source[26]
                  - c27 * source[41] + c28 * source[36] + c28 * source[38];
  }
}

void CCarSphList::carsph_32(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c12 = 6.7082039324993694;
  const double c16 = 4.2426406871192848;
  const double c3 = 4.1079191812887457;
  const double c13 = 3.872983346207417;
  const double c9 = 3.3541019662496847;
  const double c25 = 2.598076211353316;
  const double c18 = 2.4494897427831779;
  const double c6 = 2.3717082451262845;
  const double c14 = 2.1213203435596424;
  const double c1 = 2.0539595906443728;
  const double c10 = 1.9364916731037085;
  const double c24 = 1.7320508075688772;
  const double c8 = 1.6770509831248424;
  const double c27 = 1.5;
  const double c2 = 1.3693063937629153;
  const double c23 = 1.299038105676658;
  const double c19 = 1.2247448713915889;
  const double c7 = 1.1858541225631423;
  const double c17 = 1.0606601717798212;
  const double c11 = 0.96824583655185426;
  const double c22 = 0.8660254037844386;
  const double c4 = 0.79056941504209488;
  const double c28 = 0.75;
  const double c0 = 0.68465319688145765;
  const double c20 = 0.61237243569579447;
  const double c15 = 0.5303300858899106;
  const double c26 = 0.5;
  const double c5 = 0.39528470752104744;
  const double c21 = 0.30618621784789724;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 35, source += 60) {
    target[0] =  c0 * source[0] - c0 * source[2] - c1 * source[12]
                  + c1 * source[14];
    target[1] =  c2 * source[1] - c3 * source[13];
    target[2] =  c2 * source[3] - c3 * source[15];
    target[3] =  c2 * source[4] - c3 * source[16];
    target[4] =  c4 * source[5] - c5 * source[0] - c5 * source[2]
                  - c6 * source[17] + c7 * source[12] + c7 * source[14];
    target[5] =  c1 * source[6] - c1 * source[8] - c0 * source[18]
                  + c0 * source[20];
    target[6] =  c3 * source[7] - c2 * source[19];
    target[7] =  c3 * source[9] - c2 * source[21];
    target[8] =  c3 * source[10] - c2 * source[22];
    target[9] =  c6 * source[11] - c7 * source[6] - c7 * source[8]
                  - c4 * source[23] + c5 * source[18] + c5 * source[20];
    target[10] =  c8 * source[24] - c8 * source[26] - c8 * source[36]
                  + c8 * source[38];
    target[11] =  c9 * source[25] - c9 * source[37];
    target[12] =  c9 * source[27] - c9 * source[39];
    target[13] =  c9 * source[28] - c9 * source[40];
    target[14] =  c10 * source[29] - c11 * source[24] - c11 * source[26]
                  - c10 * source[41] + c11 * source[36] + c11 * source[38];
    target[15] =  c9 * source[30] - c9 * source[32];
    target[16] =  c12 * source[31];
    target[17] =  c12 * source[33];
    target[18] =  c12 * source[34];
    target[19] =  c13 * source[35] - c10 * source[30] - c10 * source[32];
    target[20] =  c14 * source[42] - c14 * source[44] - c15 * source[0]
                  + c15 * source[2] - c15 * source[12] + c15 * source[14];
    target[21] =  c16 * source[43] - c17 * source[1] - c17 * source[13];
    target[22] =  c16 * source[45] - c17 * source[3] - c17 * source[15];
    target[23] =  c16 * source[46] - c17 * source[4] - c17 * source[16];
    target[24] =  c18 * source[47] - c19 * source[42] - c19 * source[44]
                  - c20 * source[5] + c21 * source[0] + c21 * source[2]
                  - c20 * source[17] + c21 * source[12] + c21 * source[14];
    target[25] =  c14 * source[48] - c14 * source[50] - c15 * source[6]
                  + c15 * source[8] - c15 * source[18] + c15 * source[20];
    target[26] =  c16 * source[49] - c17 * source[7] - c17 * source[19];
    target[27] =  c16 * source[51] - c17 * source[9] - c17 * source[21];
    target[28] =  c16 * source[52] - c17 * source[10] - c17 * source[22];
    target[29] =  c18 * source[53] - c19 * source[48] - c19 * source[50]
                  - c20 * source[11] + c21 * source[6] + c21 * source[8]
                  - c20 * source[23] + c21 * source[18] + c21 * source[20];
    target[30] =  c22 * source[54] - c22 * source[56] - c23 * source[24]
                  + c23 * source[26] - c23 * source[36] + c23 * source[38];
    target[31] =  c24 * source[55] - c25 * source[25] - c25 * source[37];
    target[32] =  c24 * source[57] - c25 * source[27] - c25 * source[39];
    target[33] =  c24 * source[58] - c25 * source[28] - c25 * source[40];
    target[34] =  source[59] - c26 * source[54] - c26 * source[56]
                  - c27 * source[29] + c28 * source[24] + c28 * source[26]
                  - c27 * source[41] + c28 * source[36] + c28 * source[38];
  }
}

