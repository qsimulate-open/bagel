//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_41.cc
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


void CarSphList::carsph_41(const int nloop, const double* source, double* target) {
  const double c7 = 6.7082039324993694;
  const double c4 = 6.2749501990055663;
  const double c1 = 4.4370598373247123;
  const double c5 = 3.3541019662496847;
  const double c9 = 3.1622776601683795;
  const double c11 = 3;
  const double c2 = 2.9580398915498081;
  const double c10 = 2.3717082451262845;
  const double c3 = 2.0916500663351889;
  const double c8 = 1.1180339887498949;
  const double c13 = 0.75;
  const double c0 = 0.73950997288745202;
  const double c6 = 0.55901699437494745;
  const double c12 = 0.375;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 27, source += 45) {
    target[0] =  c0 * source[0] - c1 * source[6] + c0 * source[12];
    target[1] =  c0 * source[1] - c1 * source[7] + c0 * source[13];
    target[2] =  c0 * source[2] - c1 * source[8] + c0 * source[14];
    target[3] =  c2 * source[3] - c2 * source[9];
    target[4] =  c2 * source[4] - c2 * source[10];
    target[5] =  c2 * source[5] - c2 * source[11];
    target[6] =  c3 * source[15] - c4 * source[21];
    target[7] =  c3 * source[16] - c4 * source[22];
    target[8] =  c3 * source[17] - c4 * source[23];
    target[9] =  c4 * source[18] - c3 * source[24];
    target[10] =  c4 * source[19] - c3 * source[25];
    target[11] =  c4 * source[20] - c3 * source[26];
    target[12] =  c5 * source[27] - c5 * source[33] - c6 * source[0]
                  + c6 * source[6] - c6 * source[6] + c6 * source[12];
    target[13] =  c5 * source[28] - c5 * source[34] - c6 * source[1]
                  + c6 * source[7] - c6 * source[7] + c6 * source[13];
    target[14] =  c5 * source[29] - c5 * source[35] - c6 * source[2]
                  + c6 * source[8] - c6 * source[8] + c6 * source[14];
    target[15] =  c7 * source[30] - c8 * source[3] - c8 * source[9];
    target[16] =  c7 * source[31] - c8 * source[4] - c8 * source[10];
    target[17] =  c7 * source[32] - c8 * source[5] - c8 * source[11];
    target[18] =  c9 * source[36] - c10 * source[15] - c10 * source[21];
    target[19] =  c9 * source[37] - c10 * source[16] - c10 * source[22];
    target[20] =  c9 * source[38] - c10 * source[17] - c10 * source[23];
    target[21] =  c9 * source[39] - c10 * source[18] - c10 * source[24];
    target[22] =  c9 * source[40] - c10 * source[19] - c10 * source[25];
    target[23] =  c9 * source[41] - c10 * source[20] - c10 * source[26];
    target[24] =  source[42] - c11 * source[27] - c11 * source[33]
                  + c12 * source[0] + c13 * source[6] + c12 * source[12];
    target[25] =  source[43] - c11 * source[28] - c11 * source[34]
                  + c12 * source[1] + c13 * source[7] + c12 * source[13];
    target[26] =  source[44] - c11 * source[29] - c11 * source[35]
                  + c12 * source[2] + c13 * source[8] + c12 * source[14];
  }
}

void CCarSphList::carsph_41(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c7 = 6.7082039324993694;
  const double c4 = 6.2749501990055663;
  const double c1 = 4.4370598373247123;
  const double c5 = 3.3541019662496847;
  const double c9 = 3.1622776601683795;
  const double c11 = 3;
  const double c2 = 2.9580398915498081;
  const double c10 = 2.3717082451262845;
  const double c3 = 2.0916500663351889;
  const double c8 = 1.1180339887498949;
  const double c13 = 0.75;
  const double c0 = 0.73950997288745202;
  const double c6 = 0.55901699437494745;
  const double c12 = 0.375;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 27, source += 45) {
    target[0] =  c0 * source[0] - c1 * source[6] + c0 * source[12];
    target[1] =  c0 * source[1] - c1 * source[7] + c0 * source[13];
    target[2] =  c0 * source[2] - c1 * source[8] + c0 * source[14];
    target[3] =  c2 * source[3] - c2 * source[9];
    target[4] =  c2 * source[4] - c2 * source[10];
    target[5] =  c2 * source[5] - c2 * source[11];
    target[6] =  c3 * source[15] - c4 * source[21];
    target[7] =  c3 * source[16] - c4 * source[22];
    target[8] =  c3 * source[17] - c4 * source[23];
    target[9] =  c4 * source[18] - c3 * source[24];
    target[10] =  c4 * source[19] - c3 * source[25];
    target[11] =  c4 * source[20] - c3 * source[26];
    target[12] =  c5 * source[27] - c5 * source[33] - c6 * source[0]
                  + c6 * source[6] - c6 * source[6] + c6 * source[12];
    target[13] =  c5 * source[28] - c5 * source[34] - c6 * source[1]
                  + c6 * source[7] - c6 * source[7] + c6 * source[13];
    target[14] =  c5 * source[29] - c5 * source[35] - c6 * source[2]
                  + c6 * source[8] - c6 * source[8] + c6 * source[14];
    target[15] =  c7 * source[30] - c8 * source[3] - c8 * source[9];
    target[16] =  c7 * source[31] - c8 * source[4] - c8 * source[10];
    target[17] =  c7 * source[32] - c8 * source[5] - c8 * source[11];
    target[18] =  c9 * source[36] - c10 * source[15] - c10 * source[21];
    target[19] =  c9 * source[37] - c10 * source[16] - c10 * source[22];
    target[20] =  c9 * source[38] - c10 * source[17] - c10 * source[23];
    target[21] =  c9 * source[39] - c10 * source[18] - c10 * source[24];
    target[22] =  c9 * source[40] - c10 * source[19] - c10 * source[25];
    target[23] =  c9 * source[41] - c10 * source[20] - c10 * source[26];
    target[24] =  source[42] - c11 * source[27] - c11 * source[33]
                  + c12 * source[0] + c13 * source[6] + c12 * source[12];
    target[25] =  source[43] - c11 * source[28] - c11 * source[34]
                  + c12 * source[1] + c13 * source[7] + c12 * source[13];
    target[26] =  source[44] - c11 * source[29] - c11 * source[35]
                  + c12 * source[2] + c13 * source[8] + c12 * source[14];
  }
}

