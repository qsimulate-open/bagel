//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_60.cc
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


void CarSphList::carsph_60(const int nloop, const double* source, double* target) {
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
  for (int iloop = 0; iloop != nloop; ++iloop, target += 13, source += 28) {
    target[0] =  c0 * source[0] - c1 * source[2] + c1 * source[4]
                  - c0 * source[6];
    target[1] =  c2 * source[1] - c3 * source[3] + c2 * source[5];
    target[2] =  c4 * source[7] - c5 * source[9] + c6 * source[11];
    target[3] =  c6 * source[8] - c5 * source[10] + c4 * source[12];
    target[4] =  c7 * source[13] - c8 * source[15] + c7 * source[17]
                  - c9 * source[0] + c10 * source[2] - c9 * source[4]
                  - c9 * source[2] + c10 * source[4] - c9 * source[6];
    target[5] =  c11 * source[14] - c11 * source[16] - c12 * source[1]
                  + c12 * source[3] - c12 * source[3] + c12 * source[5];
    target[6] =  c13 * source[18] - c14 * source[20] - c15 * source[7]
                  + c16 * source[9] - c15 * source[9] + c16 * source[11];
    target[7] =  c14 * source[19] - c13 * source[21] - c16 * source[8]
                  + c15 * source[10] - c16 * source[10] + c15 * source[12];
    target[8] =  c13 * source[22] - c13 * source[24] - c13 * source[13]
                  + c13 * source[15] - c13 * source[15] + c13 * source[17]
                  + c17 * source[0] - c17 * source[2] + c18 * source[2]
                  - c18 * source[4] + c17 * source[4] - c17 * source[6];
    target[9] =  c19 * source[23] - c19 * source[14] - c19 * source[16]
                  + c18 * source[1] + c20 * source[3] + c18 * source[5];
    target[10] =  c21 * source[25] - c22 * source[18] - c22 * source[20]
                  + c23 * source[7] + c24 * source[9] + c23 * source[11];
    target[11] =  c21 * source[26] - c22 * source[19] - c22 * source[21]
                  + c23 * source[8] + c24 * source[10] + c23 * source[12];
    target[12] =  source[27] - c25 * source[22] - c25 * source[24]
                  + c26 * source[13] + c27 * source[15] + c26 * source[17]
                  - c28 * source[0] - c29 * source[2] - c29 * source[4]
                  - c28 * source[6];
  }
}

void CCarSphList::carsph_60(const int nloop, const complex<double>* source, complex<double>* target) {
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
  for (int iloop = 0; iloop != nloop; ++iloop, target += 13, source += 28) {
    target[0] =  c0 * source[0] - c1 * source[2] + c1 * source[4]
                  - c0 * source[6];
    target[1] =  c2 * source[1] - c3 * source[3] + c2 * source[5];
    target[2] =  c4 * source[7] - c5 * source[9] + c6 * source[11];
    target[3] =  c6 * source[8] - c5 * source[10] + c4 * source[12];
    target[4] =  c7 * source[13] - c8 * source[15] + c7 * source[17]
                  - c9 * source[0] + c10 * source[2] - c9 * source[4]
                  - c9 * source[2] + c10 * source[4] - c9 * source[6];
    target[5] =  c11 * source[14] - c11 * source[16] - c12 * source[1]
                  + c12 * source[3] - c12 * source[3] + c12 * source[5];
    target[6] =  c13 * source[18] - c14 * source[20] - c15 * source[7]
                  + c16 * source[9] - c15 * source[9] + c16 * source[11];
    target[7] =  c14 * source[19] - c13 * source[21] - c16 * source[8]
                  + c15 * source[10] - c16 * source[10] + c15 * source[12];
    target[8] =  c13 * source[22] - c13 * source[24] - c13 * source[13]
                  + c13 * source[15] - c13 * source[15] + c13 * source[17]
                  + c17 * source[0] - c17 * source[2] + c18 * source[2]
                  - c18 * source[4] + c17 * source[4] - c17 * source[6];
    target[9] =  c19 * source[23] - c19 * source[14] - c19 * source[16]
                  + c18 * source[1] + c20 * source[3] + c18 * source[5];
    target[10] =  c21 * source[25] - c22 * source[18] - c22 * source[20]
                  + c23 * source[7] + c24 * source[9] + c23 * source[11];
    target[11] =  c21 * source[26] - c22 * source[19] - c22 * source[21]
                  + c23 * source[8] + c24 * source[10] + c23 * source[12];
    target[12] =  source[27] - c25 * source[22] - c25 * source[24]
                  + c26 * source[13] + c27 * source[15] + c26 * source[17]
                  - c28 * source[0] - c29 * source[2] - c29 * source[4]
                  - c28 * source[6];
  }
}

