//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_33.cc
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


void CarSphList::carsph_33(const int nloop, const double* source, double* target) {
  const double c19 = 15;
  const double c20 = 9.4868329805051381;
  const double c6 = 9.1855865354369186;
  const double c16 = 7.5;
  const double c22 = 6;
  const double c9 = 5.8094750193111251;
  const double c2 = 5.625;
  const double c17 = 4.7434164902525691;
  const double c4 = 4.5927932677184593;
  const double c21 = 3.872983346207417;
  const double c15 = 3.75;
  const double c26 = 3.6742346141747673;
  const double c14 = 3.5575623676894268;
  const double c5 = 3.0618621784789726;
  const double c18 = 2.9047375096555625;
  const double c25 = 2.4494897427831779;
  const double c13 = 2.3717082451262845;
  const double c29 = 2.25;
  const double c7 = 1.9364916731037085;
  const double c1 = 1.875;
  const double c3 = 1.5309310892394863;
  const double c23 = 1.5;
  const double c10 = 1.4523687548277813;
  const double c12 = 1.1858541225631423;
  const double c28 = 0.91855865354369182;
  const double c11 = 0.79056941504209488;
  const double c0 = 0.625;
  const double c27 = 0.61237243569579447;
  const double c8 = 0.48412291827592713;
  const double c24 = 0.375;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 49, source += 100) {
    target[0] =  c0 * source[0] - c1 * source[2] - c1 * source[20]
                  + c2 * source[22];
    target[1] =  c1 * source[1] - c0 * source[3] - c2 * source[21]
                  + c1 * source[23];
    target[2] =  c3 * source[4] - c3 * source[6] - c4 * source[24]
                  + c4 * source[26];
    target[3] =  c5 * source[5] - c6 * source[25];
    target[4] =  c7 * source[7] - c8 * source[0] - c8 * source[2]
                  - c9 * source[27] + c10 * source[20] + c10 * source[22];
    target[5] =  c7 * source[8] - c8 * source[1] - c8 * source[3]
                  - c9 * source[28] + c10 * source[21] + c10 * source[23];
    target[6] =  c11 * source[9] - c12 * source[4] - c12 * source[6]
                  - c13 * source[29] + c14 * source[24] + c14 * source[26];
    target[7] =  c1 * source[10] - c2 * source[12] - c0 * source[30]
                  + c1 * source[32];
    target[8] =  c2 * source[11] - c1 * source[13] - c1 * source[31]
                  + c0 * source[33];
    target[9] =  c4 * source[14] - c4 * source[16] - c3 * source[34]
                  + c3 * source[36];
    target[10] =  c6 * source[15] - c5 * source[35];
    target[11] =  c9 * source[17] - c10 * source[10] - c10 * source[12]
                  - c7 * source[37] + c8 * source[30] + c8 * source[32];
    target[12] =  c9 * source[18] - c10 * source[11] - c10 * source[13]
                  - c7 * source[38] + c8 * source[31] + c8 * source[33];
    target[13] =  c13 * source[19] - c14 * source[14] - c14 * source[16]
                  - c11 * source[39] + c12 * source[34] + c12 * source[36];
    target[14] =  c3 * source[40] - c4 * source[42] - c3 * source[60]
                  + c4 * source[62];
    target[15] =  c4 * source[41] - c3 * source[43] - c4 * source[61]
                  + c3 * source[63];
    target[16] =  c15 * source[44] - c15 * source[46] - c15 * source[64]
                  + c15 * source[66];
    target[17] =  c16 * source[45] - c16 * source[65];
    target[18] =  c17 * source[47] - c12 * source[40] - c12 * source[42]
                  - c17 * source[67] + c12 * source[60] + c12 * source[62];
    target[19] =  c17 * source[48] - c12 * source[41] - c12 * source[43]
                  - c17 * source[68] + c12 * source[61] + c12 * source[63];
    target[20] =  c7 * source[49] - c18 * source[44] - c18 * source[46]
                  - c7 * source[69] + c18 * source[64] + c18 * source[66];
    target[21] =  c5 * source[50] - c6 * source[52];
    target[22] =  c6 * source[51] - c5 * source[53];
    target[23] =  c16 * source[54] - c16 * source[56];
    target[24] =  c19 * source[55];
    target[25] =  c20 * source[57] - c13 * source[50] - c13 * source[52];
    target[26] =  c20 * source[58] - c13 * source[51] - c13 * source[53];
    target[27] =  c21 * source[59] - c9 * source[54] - c9 * source[56];
    target[28] =  c7 * source[70] - c9 * source[72] - c8 * source[0]
                  + c10 * source[2] - c8 * source[20] + c10 * source[22];
    target[29] =  c9 * source[71] - c7 * source[73] - c10 * source[1]
                  + c8 * source[3] - c10 * source[21] + c8 * source[23];
    target[30] =  c17 * source[74] - c17 * source[76] - c12 * source[4]
                  + c12 * source[6] - c12 * source[24] + c12 * source[26];
    target[31] =  c20 * source[75] - c13 * source[5] - c13 * source[25];
    target[32] =  c22 * source[77] - c23 * source[70] - c23 * source[72]
                  - c23 * source[7] + c24 * source[0] + c24 * source[2]
                  - c23 * source[27] + c24 * source[20] + c24 * source[22];
    target[33] =  c22 * source[78] - c23 * source[71] - c23 * source[73]
                  - c23 * source[8] + c24 * source[1] + c24 * source[3]
                  - c23 * source[28] + c24 * source[21] + c24 * source[23];
    target[34] =  c25 * source[79] - c26 * source[74] - c26 * source[76]
                  - c27 * source[9] + c28 * source[4] + c28 * source[6]
                  - c27 * source[29] + c28 * source[24] + c28 * source[26];
    target[35] =  c7 * source[80] - c9 * source[82] - c8 * source[10]
                  + c10 * source[12] - c8 * source[30] + c10 * source[32];
    target[36] =  c9 * source[81] - c7 * source[83] - c10 * source[11]
                  + c8 * source[13] - c10 * source[31] + c8 * source[33];
    target[37] =  c17 * source[84] - c17 * source[86] - c12 * source[14]
                  + c12 * source[16] - c12 * source[34] + c12 * source[36];
    target[38] =  c20 * source[85] - c13 * source[15] - c13 * source[35];
    target[39] =  c22 * source[87] - c23 * source[80] - c23 * source[82]
                  - c23 * source[17] + c24 * source[10] + c24 * source[12]
                  - c23 * source[37] + c24 * source[30] + c24 * source[32];
    target[40] =  c22 * source[88] - c23 * source[81] - c23 * source[83]
                  - c23 * source[18] + c24 * source[11] + c24 * source[13]
                  - c23 * source[38] + c24 * source[31] + c24 * source[33];
    target[41] =  c25 * source[89] - c26 * source[84] - c26 * source[86]
                  - c27 * source[19] + c28 * source[14] + c28 * source[16]
                  - c27 * source[39] + c28 * source[34] + c28 * source[36];
    target[42] =  c11 * source[90] - c13 * source[92] - c12 * source[40]
                  + c14 * source[42] - c12 * source[60] + c14 * source[62];
    target[43] =  c13 * source[91] - c11 * source[93] - c14 * source[41]
                  + c12 * source[43] - c14 * source[61] + c12 * source[63];
    target[44] =  c7 * source[94] - c7 * source[96] - c18 * source[44]
                  + c18 * source[46] - c18 * source[64] + c18 * source[66];
    target[45] =  c21 * source[95] - c9 * source[45] - c9 * source[65];
    target[46] =  c25 * source[97] - c27 * source[90] - c27 * source[92]
                  - c26 * source[47] + c28 * source[40] + c28 * source[42]
                  - c26 * source[67] + c28 * source[60] + c28 * source[62];
    target[47] =  c25 * source[98] - c27 * source[91] - c27 * source[93]
                  - c26 * source[48] + c28 * source[41] + c28 * source[43]
                  - c26 * source[68] + c28 * source[61] + c28 * source[63];
    target[48] =  source[99] - c23 * source[94] - c23 * source[96]
                  - c23 * source[49] + c29 * source[44] + c29 * source[46]
                  - c23 * source[69] + c29 * source[64] + c29 * source[66];
  }
}

void CCarSphList::carsph_33(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c19 = 15;
  const double c20 = 9.4868329805051381;
  const double c6 = 9.1855865354369186;
  const double c16 = 7.5;
  const double c22 = 6;
  const double c9 = 5.8094750193111251;
  const double c2 = 5.625;
  const double c17 = 4.7434164902525691;
  const double c4 = 4.5927932677184593;
  const double c21 = 3.872983346207417;
  const double c15 = 3.75;
  const double c26 = 3.6742346141747673;
  const double c14 = 3.5575623676894268;
  const double c5 = 3.0618621784789726;
  const double c18 = 2.9047375096555625;
  const double c25 = 2.4494897427831779;
  const double c13 = 2.3717082451262845;
  const double c29 = 2.25;
  const double c7 = 1.9364916731037085;
  const double c1 = 1.875;
  const double c3 = 1.5309310892394863;
  const double c23 = 1.5;
  const double c10 = 1.4523687548277813;
  const double c12 = 1.1858541225631423;
  const double c28 = 0.91855865354369182;
  const double c11 = 0.79056941504209488;
  const double c0 = 0.625;
  const double c27 = 0.61237243569579447;
  const double c8 = 0.48412291827592713;
  const double c24 = 0.375;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 49, source += 100) {
    target[0] =  c0 * source[0] - c1 * source[2] - c1 * source[20]
                  + c2 * source[22];
    target[1] =  c1 * source[1] - c0 * source[3] - c2 * source[21]
                  + c1 * source[23];
    target[2] =  c3 * source[4] - c3 * source[6] - c4 * source[24]
                  + c4 * source[26];
    target[3] =  c5 * source[5] - c6 * source[25];
    target[4] =  c7 * source[7] - c8 * source[0] - c8 * source[2]
                  - c9 * source[27] + c10 * source[20] + c10 * source[22];
    target[5] =  c7 * source[8] - c8 * source[1] - c8 * source[3]
                  - c9 * source[28] + c10 * source[21] + c10 * source[23];
    target[6] =  c11 * source[9] - c12 * source[4] - c12 * source[6]
                  - c13 * source[29] + c14 * source[24] + c14 * source[26];
    target[7] =  c1 * source[10] - c2 * source[12] - c0 * source[30]
                  + c1 * source[32];
    target[8] =  c2 * source[11] - c1 * source[13] - c1 * source[31]
                  + c0 * source[33];
    target[9] =  c4 * source[14] - c4 * source[16] - c3 * source[34]
                  + c3 * source[36];
    target[10] =  c6 * source[15] - c5 * source[35];
    target[11] =  c9 * source[17] - c10 * source[10] - c10 * source[12]
                  - c7 * source[37] + c8 * source[30] + c8 * source[32];
    target[12] =  c9 * source[18] - c10 * source[11] - c10 * source[13]
                  - c7 * source[38] + c8 * source[31] + c8 * source[33];
    target[13] =  c13 * source[19] - c14 * source[14] - c14 * source[16]
                  - c11 * source[39] + c12 * source[34] + c12 * source[36];
    target[14] =  c3 * source[40] - c4 * source[42] - c3 * source[60]
                  + c4 * source[62];
    target[15] =  c4 * source[41] - c3 * source[43] - c4 * source[61]
                  + c3 * source[63];
    target[16] =  c15 * source[44] - c15 * source[46] - c15 * source[64]
                  + c15 * source[66];
    target[17] =  c16 * source[45] - c16 * source[65];
    target[18] =  c17 * source[47] - c12 * source[40] - c12 * source[42]
                  - c17 * source[67] + c12 * source[60] + c12 * source[62];
    target[19] =  c17 * source[48] - c12 * source[41] - c12 * source[43]
                  - c17 * source[68] + c12 * source[61] + c12 * source[63];
    target[20] =  c7 * source[49] - c18 * source[44] - c18 * source[46]
                  - c7 * source[69] + c18 * source[64] + c18 * source[66];
    target[21] =  c5 * source[50] - c6 * source[52];
    target[22] =  c6 * source[51] - c5 * source[53];
    target[23] =  c16 * source[54] - c16 * source[56];
    target[24] =  c19 * source[55];
    target[25] =  c20 * source[57] - c13 * source[50] - c13 * source[52];
    target[26] =  c20 * source[58] - c13 * source[51] - c13 * source[53];
    target[27] =  c21 * source[59] - c9 * source[54] - c9 * source[56];
    target[28] =  c7 * source[70] - c9 * source[72] - c8 * source[0]
                  + c10 * source[2] - c8 * source[20] + c10 * source[22];
    target[29] =  c9 * source[71] - c7 * source[73] - c10 * source[1]
                  + c8 * source[3] - c10 * source[21] + c8 * source[23];
    target[30] =  c17 * source[74] - c17 * source[76] - c12 * source[4]
                  + c12 * source[6] - c12 * source[24] + c12 * source[26];
    target[31] =  c20 * source[75] - c13 * source[5] - c13 * source[25];
    target[32] =  c22 * source[77] - c23 * source[70] - c23 * source[72]
                  - c23 * source[7] + c24 * source[0] + c24 * source[2]
                  - c23 * source[27] + c24 * source[20] + c24 * source[22];
    target[33] =  c22 * source[78] - c23 * source[71] - c23 * source[73]
                  - c23 * source[8] + c24 * source[1] + c24 * source[3]
                  - c23 * source[28] + c24 * source[21] + c24 * source[23];
    target[34] =  c25 * source[79] - c26 * source[74] - c26 * source[76]
                  - c27 * source[9] + c28 * source[4] + c28 * source[6]
                  - c27 * source[29] + c28 * source[24] + c28 * source[26];
    target[35] =  c7 * source[80] - c9 * source[82] - c8 * source[10]
                  + c10 * source[12] - c8 * source[30] + c10 * source[32];
    target[36] =  c9 * source[81] - c7 * source[83] - c10 * source[11]
                  + c8 * source[13] - c10 * source[31] + c8 * source[33];
    target[37] =  c17 * source[84] - c17 * source[86] - c12 * source[14]
                  + c12 * source[16] - c12 * source[34] + c12 * source[36];
    target[38] =  c20 * source[85] - c13 * source[15] - c13 * source[35];
    target[39] =  c22 * source[87] - c23 * source[80] - c23 * source[82]
                  - c23 * source[17] + c24 * source[10] + c24 * source[12]
                  - c23 * source[37] + c24 * source[30] + c24 * source[32];
    target[40] =  c22 * source[88] - c23 * source[81] - c23 * source[83]
                  - c23 * source[18] + c24 * source[11] + c24 * source[13]
                  - c23 * source[38] + c24 * source[31] + c24 * source[33];
    target[41] =  c25 * source[89] - c26 * source[84] - c26 * source[86]
                  - c27 * source[19] + c28 * source[14] + c28 * source[16]
                  - c27 * source[39] + c28 * source[34] + c28 * source[36];
    target[42] =  c11 * source[90] - c13 * source[92] - c12 * source[40]
                  + c14 * source[42] - c12 * source[60] + c14 * source[62];
    target[43] =  c13 * source[91] - c11 * source[93] - c14 * source[41]
                  + c12 * source[43] - c14 * source[61] + c12 * source[63];
    target[44] =  c7 * source[94] - c7 * source[96] - c18 * source[44]
                  + c18 * source[46] - c18 * source[64] + c18 * source[66];
    target[45] =  c21 * source[95] - c9 * source[45] - c9 * source[65];
    target[46] =  c25 * source[97] - c27 * source[90] - c27 * source[92]
                  - c26 * source[47] + c28 * source[40] + c28 * source[42]
                  - c26 * source[67] + c28 * source[60] + c28 * source[62];
    target[47] =  c25 * source[98] - c27 * source[91] - c27 * source[93]
                  - c26 * source[48] + c28 * source[41] + c28 * source[43]
                  - c26 * source[68] + c28 * source[61] + c28 * source[63];
    target[48] =  source[99] - c23 * source[94] - c23 * source[96]
                  - c23 * source[49] + c29 * source[44] + c29 * source[46]
                  - c23 * source[69] + c29 * source[64] + c29 * source[66];
  }
}

