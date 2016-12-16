//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_61.cc
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


void CarSphList::carsph_61(const int nloop, const double* source, double* target) {
  const double c8 = 29.764702249476645;
  const double c5 = 23.268138086232856;
  const double c14 = 21.737065119284157;
  const double c11 = 19.843134832984429;
  const double c19 = 14.491376746189438;
  const double c3 = 13.433865787627923;
  const double c6 = 11.634069043116428;
  const double c22 = 11.456439237389599;
  const double c27 = 11.25;
  const double c1 = 10.075399340720942;
  const double c16 = 8.1513994197315593;
  const double c25 = 7.5;
  const double c13 = 7.245688373094719;
  const double c24 = 5.7282196186947996;
  const double c26 = 5.625;
  const double c7 = 4.9607837082461073;
  const double c21 = 4.5825756949558398;
  const double c2 = 4.0301597362883772;
  const double c10 = 2.9764702249476644;
  const double c23 = 2.8641098093473998;
  const double c15 = 2.7171331399105196;
  const double c4 = 2.3268138086232857;
  const double c12 = 1.984313483298443;
  const double c20 = 1.8114220932736798;
  const double c29 = 0.9375;
  const double c18 = 0.90571104663683988;
  const double c0 = 0.67169328938139616;
  const double c9 = 0.49607837082461076;
  const double c17 = 0.45285552331841994;
  const double c28 = 0.3125;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 39, source += 84) {
    target[0] =  c0 * source[0] - c1 * source[6] + c1 * source[12]
                  - c0 * source[18];
    target[1] =  c0 * source[1] - c1 * source[7] + c1 * source[13]
                  - c0 * source[19];
    target[2] =  c0 * source[2] - c1 * source[8] + c1 * source[14]
                  - c0 * source[20];
    target[3] =  c2 * source[3] - c3 * source[9] + c2 * source[15];
    target[4] =  c2 * source[4] - c3 * source[10] + c2 * source[16];
    target[5] =  c2 * source[5] - c3 * source[11] + c2 * source[17];
    target[6] =  c4 * source[21] - c5 * source[27] + c6 * source[33];
    target[7] =  c4 * source[22] - c5 * source[28] + c6 * source[34];
    target[8] =  c4 * source[23] - c5 * source[29] + c6 * source[35];
    target[9] =  c6 * source[24] - c5 * source[30] + c4 * source[36];
    target[10] =  c6 * source[25] - c5 * source[31] + c4 * source[37];
    target[11] =  c6 * source[26] - c5 * source[32] + c4 * source[38];
    target[12] =  c7 * source[39] - c8 * source[45] + c7 * source[51]
                  - c9 * source[0] + c10 * source[6] - c9 * source[12]
                  - c9 * source[6] + c10 * source[12] - c9 * source[18];
    target[13] =  c7 * source[40] - c8 * source[46] + c7 * source[52]
                  - c9 * source[1] + c10 * source[7] - c9 * source[13]
                  - c9 * source[7] + c10 * source[13] - c9 * source[19];
    target[14] =  c7 * source[41] - c8 * source[47] + c7 * source[53]
                  - c9 * source[2] + c10 * source[8] - c9 * source[14]
                  - c9 * source[8] + c10 * source[14] - c9 * source[20];
    target[15] =  c11 * source[42] - c11 * source[48] - c12 * source[3]
                  + c12 * source[9] - c12 * source[9] + c12 * source[15];
    target[16] =  c11 * source[43] - c11 * source[49] - c12 * source[4]
                  + c12 * source[10] - c12 * source[10] + c12 * source[16];
    target[17] =  c11 * source[44] - c11 * source[50] - c12 * source[5]
                  + c12 * source[11] - c12 * source[11] + c12 * source[17];
    target[18] =  c13 * source[54] - c14 * source[60] - c15 * source[21]
                  + c16 * source[27] - c15 * source[27] + c16 * source[33];
    target[19] =  c13 * source[55] - c14 * source[61] - c15 * source[22]
                  + c16 * source[28] - c15 * source[28] + c16 * source[34];
    target[20] =  c13 * source[56] - c14 * source[62] - c15 * source[23]
                  + c16 * source[29] - c15 * source[29] + c16 * source[35];
    target[21] =  c14 * source[57] - c13 * source[63] - c16 * source[24]
                  + c15 * source[30] - c16 * source[30] + c15 * source[36];
    target[22] =  c14 * source[58] - c13 * source[64] - c16 * source[25]
                  + c15 * source[31] - c16 * source[31] + c15 * source[37];
    target[23] =  c14 * source[59] - c13 * source[65] - c16 * source[26]
                  + c15 * source[32] - c16 * source[32] + c15 * source[38];
    target[24] =  c13 * source[66] - c13 * source[72] - c13 * source[39]
                  + c13 * source[45] - c13 * source[45] + c13 * source[51]
                  + c17 * source[0] - c17 * source[6] + c18 * source[6]
                  - c18 * source[12] + c17 * source[12] - c17 * source[18];
    target[25] =  c13 * source[67] - c13 * source[73] - c13 * source[40]
                  + c13 * source[46] - c13 * source[46] + c13 * source[52]
                  + c17 * source[1] - c17 * source[7] + c18 * source[7]
                  - c18 * source[13] + c17 * source[13] - c17 * source[19];
    target[26] =  c13 * source[68] - c13 * source[74] - c13 * source[41]
                  + c13 * source[47] - c13 * source[47] + c13 * source[53]
                  + c17 * source[2] - c17 * source[8] + c18 * source[8]
                  - c18 * source[14] + c17 * source[14] - c17 * source[20];
    target[27] =  c19 * source[69] - c19 * source[42] - c19 * source[48]
                  + c18 * source[3] + c20 * source[9] + c18 * source[15];
    target[28] =  c19 * source[70] - c19 * source[43] - c19 * source[49]
                  + c18 * source[4] + c20 * source[10] + c18 * source[16];
    target[29] =  c19 * source[71] - c19 * source[44] - c19 * source[50]
                  + c18 * source[5] + c20 * source[11] + c18 * source[17];
    target[30] =  c21 * source[75] - c22 * source[54] - c22 * source[60]
                  + c23 * source[21] + c24 * source[27] + c23 * source[33];
    target[31] =  c21 * source[76] - c22 * source[55] - c22 * source[61]
                  + c23 * source[22] + c24 * source[28] + c23 * source[34];
    target[32] =  c21 * source[77] - c22 * source[56] - c22 * source[62]
                  + c23 * source[23] + c24 * source[29] + c23 * source[35];
    target[33] =  c21 * source[78] - c22 * source[57] - c22 * source[63]
                  + c23 * source[24] + c24 * source[30] + c23 * source[36];
    target[34] =  c21 * source[79] - c22 * source[58] - c22 * source[64]
                  + c23 * source[25] + c24 * source[31] + c23 * source[37];
    target[35] =  c21 * source[80] - c22 * source[59] - c22 * source[65]
                  + c23 * source[26] + c24 * source[32] + c23 * source[38];
    target[36] =  source[81] - c25 * source[66] - c25 * source[72]
                  + c26 * source[39] + c27 * source[45] + c26 * source[51]
                  - c28 * source[0] - c29 * source[6] - c29 * source[12]
                  - c28 * source[18];
    target[37] =  source[82] - c25 * source[67] - c25 * source[73]
                  + c26 * source[40] + c27 * source[46] + c26 * source[52]
                  - c28 * source[1] - c29 * source[7] - c29 * source[13]
                  - c28 * source[19];
    target[38] =  source[83] - c25 * source[68] - c25 * source[74]
                  + c26 * source[41] + c27 * source[47] + c26 * source[53]
                  - c28 * source[2] - c29 * source[8] - c29 * source[14]
                  - c28 * source[20];
  }
}

void CCarSphList::carsph_61(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c8 = 29.764702249476645;
  const double c5 = 23.268138086232856;
  const double c14 = 21.737065119284157;
  const double c11 = 19.843134832984429;
  const double c19 = 14.491376746189438;
  const double c3 = 13.433865787627923;
  const double c6 = 11.634069043116428;
  const double c22 = 11.456439237389599;
  const double c27 = 11.25;
  const double c1 = 10.075399340720942;
  const double c16 = 8.1513994197315593;
  const double c25 = 7.5;
  const double c13 = 7.245688373094719;
  const double c24 = 5.7282196186947996;
  const double c26 = 5.625;
  const double c7 = 4.9607837082461073;
  const double c21 = 4.5825756949558398;
  const double c2 = 4.0301597362883772;
  const double c10 = 2.9764702249476644;
  const double c23 = 2.8641098093473998;
  const double c15 = 2.7171331399105196;
  const double c4 = 2.3268138086232857;
  const double c12 = 1.984313483298443;
  const double c20 = 1.8114220932736798;
  const double c29 = 0.9375;
  const double c18 = 0.90571104663683988;
  const double c0 = 0.67169328938139616;
  const double c9 = 0.49607837082461076;
  const double c17 = 0.45285552331841994;
  const double c28 = 0.3125;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 39, source += 84) {
    target[0] =  c0 * source[0] - c1 * source[6] + c1 * source[12]
                  - c0 * source[18];
    target[1] =  c0 * source[1] - c1 * source[7] + c1 * source[13]
                  - c0 * source[19];
    target[2] =  c0 * source[2] - c1 * source[8] + c1 * source[14]
                  - c0 * source[20];
    target[3] =  c2 * source[3] - c3 * source[9] + c2 * source[15];
    target[4] =  c2 * source[4] - c3 * source[10] + c2 * source[16];
    target[5] =  c2 * source[5] - c3 * source[11] + c2 * source[17];
    target[6] =  c4 * source[21] - c5 * source[27] + c6 * source[33];
    target[7] =  c4 * source[22] - c5 * source[28] + c6 * source[34];
    target[8] =  c4 * source[23] - c5 * source[29] + c6 * source[35];
    target[9] =  c6 * source[24] - c5 * source[30] + c4 * source[36];
    target[10] =  c6 * source[25] - c5 * source[31] + c4 * source[37];
    target[11] =  c6 * source[26] - c5 * source[32] + c4 * source[38];
    target[12] =  c7 * source[39] - c8 * source[45] + c7 * source[51]
                  - c9 * source[0] + c10 * source[6] - c9 * source[12]
                  - c9 * source[6] + c10 * source[12] - c9 * source[18];
    target[13] =  c7 * source[40] - c8 * source[46] + c7 * source[52]
                  - c9 * source[1] + c10 * source[7] - c9 * source[13]
                  - c9 * source[7] + c10 * source[13] - c9 * source[19];
    target[14] =  c7 * source[41] - c8 * source[47] + c7 * source[53]
                  - c9 * source[2] + c10 * source[8] - c9 * source[14]
                  - c9 * source[8] + c10 * source[14] - c9 * source[20];
    target[15] =  c11 * source[42] - c11 * source[48] - c12 * source[3]
                  + c12 * source[9] - c12 * source[9] + c12 * source[15];
    target[16] =  c11 * source[43] - c11 * source[49] - c12 * source[4]
                  + c12 * source[10] - c12 * source[10] + c12 * source[16];
    target[17] =  c11 * source[44] - c11 * source[50] - c12 * source[5]
                  + c12 * source[11] - c12 * source[11] + c12 * source[17];
    target[18] =  c13 * source[54] - c14 * source[60] - c15 * source[21]
                  + c16 * source[27] - c15 * source[27] + c16 * source[33];
    target[19] =  c13 * source[55] - c14 * source[61] - c15 * source[22]
                  + c16 * source[28] - c15 * source[28] + c16 * source[34];
    target[20] =  c13 * source[56] - c14 * source[62] - c15 * source[23]
                  + c16 * source[29] - c15 * source[29] + c16 * source[35];
    target[21] =  c14 * source[57] - c13 * source[63] - c16 * source[24]
                  + c15 * source[30] - c16 * source[30] + c15 * source[36];
    target[22] =  c14 * source[58] - c13 * source[64] - c16 * source[25]
                  + c15 * source[31] - c16 * source[31] + c15 * source[37];
    target[23] =  c14 * source[59] - c13 * source[65] - c16 * source[26]
                  + c15 * source[32] - c16 * source[32] + c15 * source[38];
    target[24] =  c13 * source[66] - c13 * source[72] - c13 * source[39]
                  + c13 * source[45] - c13 * source[45] + c13 * source[51]
                  + c17 * source[0] - c17 * source[6] + c18 * source[6]
                  - c18 * source[12] + c17 * source[12] - c17 * source[18];
    target[25] =  c13 * source[67] - c13 * source[73] - c13 * source[40]
                  + c13 * source[46] - c13 * source[46] + c13 * source[52]
                  + c17 * source[1] - c17 * source[7] + c18 * source[7]
                  - c18 * source[13] + c17 * source[13] - c17 * source[19];
    target[26] =  c13 * source[68] - c13 * source[74] - c13 * source[41]
                  + c13 * source[47] - c13 * source[47] + c13 * source[53]
                  + c17 * source[2] - c17 * source[8] + c18 * source[8]
                  - c18 * source[14] + c17 * source[14] - c17 * source[20];
    target[27] =  c19 * source[69] - c19 * source[42] - c19 * source[48]
                  + c18 * source[3] + c20 * source[9] + c18 * source[15];
    target[28] =  c19 * source[70] - c19 * source[43] - c19 * source[49]
                  + c18 * source[4] + c20 * source[10] + c18 * source[16];
    target[29] =  c19 * source[71] - c19 * source[44] - c19 * source[50]
                  + c18 * source[5] + c20 * source[11] + c18 * source[17];
    target[30] =  c21 * source[75] - c22 * source[54] - c22 * source[60]
                  + c23 * source[21] + c24 * source[27] + c23 * source[33];
    target[31] =  c21 * source[76] - c22 * source[55] - c22 * source[61]
                  + c23 * source[22] + c24 * source[28] + c23 * source[34];
    target[32] =  c21 * source[77] - c22 * source[56] - c22 * source[62]
                  + c23 * source[23] + c24 * source[29] + c23 * source[35];
    target[33] =  c21 * source[78] - c22 * source[57] - c22 * source[63]
                  + c23 * source[24] + c24 * source[30] + c23 * source[36];
    target[34] =  c21 * source[79] - c22 * source[58] - c22 * source[64]
                  + c23 * source[25] + c24 * source[31] + c23 * source[37];
    target[35] =  c21 * source[80] - c22 * source[59] - c22 * source[65]
                  + c23 * source[26] + c24 * source[32] + c23 * source[38];
    target[36] =  source[81] - c25 * source[66] - c25 * source[72]
                  + c26 * source[39] + c27 * source[45] + c26 * source[51]
                  - c28 * source[0] - c29 * source[6] - c29 * source[12]
                  - c28 * source[18];
    target[37] =  source[82] - c25 * source[67] - c25 * source[73]
                  + c26 * source[40] + c27 * source[46] + c26 * source[52]
                  - c28 * source[1] - c29 * source[7] - c29 * source[13]
                  - c28 * source[19];
    target[38] =  source[83] - c25 * source[68] - c25 * source[74]
                  + c26 * source[41] + c27 * source[47] + c26 * source[53]
                  - c28 * source[2] - c29 * source[8] - c29 * source[14]
                  - c28 * source[20];
  }
}

