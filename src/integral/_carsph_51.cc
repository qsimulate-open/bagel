//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_51.cc
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


void CarSphList::carsph_51(const int nloop, const double* source, double* target) {
  const double c4 = 13.311179511974137;
  const double c7 = 12.549900398011133;
  const double c12 = 10.246950765959598;
  const double c5 = 8.8741196746494246;
  const double c1 = 7.0156076002011405;
  const double c14 = 5.8094750193111251;
  const double c10 = 5.123475382979799;
  const double c17 = 5;
  const double c6 = 4.1833001326703778;
  const double c13 = 3.872983346207417;
  const double c19 = 3.75;
  const double c2 = 3.5078038001005702;
  const double c11 = 2.5617376914898995;
  const double c3 = 2.2185299186623562;
  const double c18 = 1.875;
  const double c9 = 1.5687375497513916;
  const double c16 = 0.96824583655185426;
  const double c0 = 0.70156076002011403;
  const double c8 = 0.52291251658379723;
  const double c15 = 0.48412291827592713;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 33, source += 63) {
    target[0] =  c0 * source[0] - c1 * source[6] + c2 * source[12];
    target[1] =  c0 * source[1] - c1 * source[7] + c2 * source[13];
    target[2] =  c0 * source[2] - c1 * source[8] + c2 * source[14];
    target[3] =  c2 * source[3] - c1 * source[9] + c0 * source[15];
    target[4] =  c2 * source[4] - c1 * source[10] + c0 * source[16];
    target[5] =  c2 * source[5] - c1 * source[11] + c0 * source[17];
    target[6] =  c3 * source[18] - c4 * source[24] + c3 * source[30];
    target[7] =  c3 * source[19] - c4 * source[25] + c3 * source[31];
    target[8] =  c3 * source[20] - c4 * source[26] + c3 * source[32];
    target[9] =  c5 * source[21] - c5 * source[27];
    target[10] =  c5 * source[22] - c5 * source[28];
    target[11] =  c5 * source[23] - c5 * source[29];
    target[12] =  c6 * source[33] - c7 * source[39] - c8 * source[0]
                  + c9 * source[6] - c8 * source[6] + c9 * source[12];
    target[13] =  c6 * source[34] - c7 * source[40] - c8 * source[1]
                  + c9 * source[7] - c8 * source[7] + c9 * source[13];
    target[14] =  c6 * source[35] - c7 * source[41] - c8 * source[2]
                  + c9 * source[8] - c8 * source[8] + c9 * source[14];
    target[15] =  c7 * source[36] - c6 * source[42] - c9 * source[3]
                  + c8 * source[9] - c9 * source[9] + c8 * source[15];
    target[16] =  c7 * source[37] - c6 * source[43] - c9 * source[4]
                  + c8 * source[10] - c9 * source[10] + c8 * source[16];
    target[17] =  c7 * source[38] - c6 * source[44] - c9 * source[5]
                  + c8 * source[11] - c9 * source[11] + c8 * source[17];
    target[18] =  c10 * source[45] - c10 * source[51] - c11 * source[18]
                  + c11 * source[24] - c11 * source[24] + c11 * source[30];
    target[19] =  c10 * source[46] - c10 * source[52] - c11 * source[19]
                  + c11 * source[25] - c11 * source[25] + c11 * source[31];
    target[20] =  c10 * source[47] - c10 * source[53] - c11 * source[20]
                  + c11 * source[26] - c11 * source[26] + c11 * source[32];
    target[21] =  c12 * source[48] - c10 * source[21] - c10 * source[27];
    target[22] =  c12 * source[49] - c10 * source[22] - c10 * source[28];
    target[23] =  c12 * source[50] - c10 * source[23] - c10 * source[29];
    target[24] =  c13 * source[54] - c14 * source[33] - c14 * source[39]
                  + c15 * source[0] + c16 * source[6] + c15 * source[12];
    target[25] =  c13 * source[55] - c14 * source[34] - c14 * source[40]
                  + c15 * source[1] + c16 * source[7] + c15 * source[13];
    target[26] =  c13 * source[56] - c14 * source[35] - c14 * source[41]
                  + c15 * source[2] + c16 * source[8] + c15 * source[14];
    target[27] =  c13 * source[57] - c14 * source[36] - c14 * source[42]
                  + c15 * source[3] + c16 * source[9] + c15 * source[15];
    target[28] =  c13 * source[58] - c14 * source[37] - c14 * source[43]
                  + c15 * source[4] + c16 * source[10] + c15 * source[16];
    target[29] =  c13 * source[59] - c14 * source[38] - c14 * source[44]
                  + c15 * source[5] + c16 * source[11] + c15 * source[17];
    target[30] =  source[60] - c17 * source[45] - c17 * source[51]
                  + c18 * source[18] + c19 * source[24] + c18 * source[30];
    target[31] =  source[61] - c17 * source[46] - c17 * source[52]
                  + c18 * source[19] + c19 * source[25] + c18 * source[31];
    target[32] =  source[62] - c17 * source[47] - c17 * source[53]
                  + c18 * source[20] + c19 * source[26] + c18 * source[32];
  }
}

void CCarSphList::carsph_51(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c4 = 13.311179511974137;
  const double c7 = 12.549900398011133;
  const double c12 = 10.246950765959598;
  const double c5 = 8.8741196746494246;
  const double c1 = 7.0156076002011405;
  const double c14 = 5.8094750193111251;
  const double c10 = 5.123475382979799;
  const double c17 = 5;
  const double c6 = 4.1833001326703778;
  const double c13 = 3.872983346207417;
  const double c19 = 3.75;
  const double c2 = 3.5078038001005702;
  const double c11 = 2.5617376914898995;
  const double c3 = 2.2185299186623562;
  const double c18 = 1.875;
  const double c9 = 1.5687375497513916;
  const double c16 = 0.96824583655185426;
  const double c0 = 0.70156076002011403;
  const double c8 = 0.52291251658379723;
  const double c15 = 0.48412291827592713;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 33, source += 63) {
    target[0] =  c0 * source[0] - c1 * source[6] + c2 * source[12];
    target[1] =  c0 * source[1] - c1 * source[7] + c2 * source[13];
    target[2] =  c0 * source[2] - c1 * source[8] + c2 * source[14];
    target[3] =  c2 * source[3] - c1 * source[9] + c0 * source[15];
    target[4] =  c2 * source[4] - c1 * source[10] + c0 * source[16];
    target[5] =  c2 * source[5] - c1 * source[11] + c0 * source[17];
    target[6] =  c3 * source[18] - c4 * source[24] + c3 * source[30];
    target[7] =  c3 * source[19] - c4 * source[25] + c3 * source[31];
    target[8] =  c3 * source[20] - c4 * source[26] + c3 * source[32];
    target[9] =  c5 * source[21] - c5 * source[27];
    target[10] =  c5 * source[22] - c5 * source[28];
    target[11] =  c5 * source[23] - c5 * source[29];
    target[12] =  c6 * source[33] - c7 * source[39] - c8 * source[0]
                  + c9 * source[6] - c8 * source[6] + c9 * source[12];
    target[13] =  c6 * source[34] - c7 * source[40] - c8 * source[1]
                  + c9 * source[7] - c8 * source[7] + c9 * source[13];
    target[14] =  c6 * source[35] - c7 * source[41] - c8 * source[2]
                  + c9 * source[8] - c8 * source[8] + c9 * source[14];
    target[15] =  c7 * source[36] - c6 * source[42] - c9 * source[3]
                  + c8 * source[9] - c9 * source[9] + c8 * source[15];
    target[16] =  c7 * source[37] - c6 * source[43] - c9 * source[4]
                  + c8 * source[10] - c9 * source[10] + c8 * source[16];
    target[17] =  c7 * source[38] - c6 * source[44] - c9 * source[5]
                  + c8 * source[11] - c9 * source[11] + c8 * source[17];
    target[18] =  c10 * source[45] - c10 * source[51] - c11 * source[18]
                  + c11 * source[24] - c11 * source[24] + c11 * source[30];
    target[19] =  c10 * source[46] - c10 * source[52] - c11 * source[19]
                  + c11 * source[25] - c11 * source[25] + c11 * source[31];
    target[20] =  c10 * source[47] - c10 * source[53] - c11 * source[20]
                  + c11 * source[26] - c11 * source[26] + c11 * source[32];
    target[21] =  c12 * source[48] - c10 * source[21] - c10 * source[27];
    target[22] =  c12 * source[49] - c10 * source[22] - c10 * source[28];
    target[23] =  c12 * source[50] - c10 * source[23] - c10 * source[29];
    target[24] =  c13 * source[54] - c14 * source[33] - c14 * source[39]
                  + c15 * source[0] + c16 * source[6] + c15 * source[12];
    target[25] =  c13 * source[55] - c14 * source[34] - c14 * source[40]
                  + c15 * source[1] + c16 * source[7] + c15 * source[13];
    target[26] =  c13 * source[56] - c14 * source[35] - c14 * source[41]
                  + c15 * source[2] + c16 * source[8] + c15 * source[14];
    target[27] =  c13 * source[57] - c14 * source[36] - c14 * source[42]
                  + c15 * source[3] + c16 * source[9] + c15 * source[15];
    target[28] =  c13 * source[58] - c14 * source[37] - c14 * source[43]
                  + c15 * source[4] + c16 * source[10] + c15 * source[16];
    target[29] =  c13 * source[59] - c14 * source[38] - c14 * source[44]
                  + c15 * source[5] + c16 * source[11] + c15 * source[17];
    target[30] =  source[60] - c17 * source[45] - c17 * source[51]
                  + c18 * source[18] + c19 * source[24] + c18 * source[30];
    target[31] =  source[61] - c17 * source[46] - c17 * source[52]
                  + c18 * source[19] + c19 * source[25] + c18 * source[31];
    target[32] =  source[62] - c17 * source[47] - c17 * source[53]
                  + c18 * source[20] + c19 * source[26] + c18 * source[32];
  }
}

