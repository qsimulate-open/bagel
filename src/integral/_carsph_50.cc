//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_50.cc
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


void CarSphList::carsph_50(const int nloop, const double* source, double* target) {
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
  for (int iloop = 0; iloop != nloop; ++iloop, target += 11, source += 21) {
    target[0] =  c0 * source[0] - c1 * source[2] + c2 * source[4];
    target[1] =  c2 * source[1] - c1 * source[3] + c0 * source[5];
    target[2] =  c3 * source[6] - c4 * source[8] + c3 * source[10];
    target[3] =  c5 * source[7] - c5 * source[9];
    target[4] =  c6 * source[11] - c7 * source[13] - c8 * source[0]
                  + c9 * source[2] - c8 * source[2] + c9 * source[4];
    target[5] =  c7 * source[12] - c6 * source[14] - c9 * source[1]
                  + c8 * source[3] - c9 * source[3] + c8 * source[5];
    target[6] =  c10 * source[15] - c10 * source[17] - c11 * source[6]
                  + c11 * source[8] - c11 * source[8] + c11 * source[10];
    target[7] =  c12 * source[16] - c10 * source[7] - c10 * source[9];
    target[8] =  c13 * source[18] - c14 * source[11] - c14 * source[13]
                  + c15 * source[0] + c16 * source[2] + c15 * source[4];
    target[9] =  c13 * source[19] - c14 * source[12] - c14 * source[14]
                  + c15 * source[1] + c16 * source[3] + c15 * source[5];
    target[10] =  source[20] - c17 * source[15] - c17 * source[17]
                  + c18 * source[6] + c19 * source[8] + c18 * source[10];
  }
}

void CCarSphList::carsph_50(const int nloop, const complex<double>* source, complex<double>* target) {
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
  for (int iloop = 0; iloop != nloop; ++iloop, target += 11, source += 21) {
    target[0] =  c0 * source[0] - c1 * source[2] + c2 * source[4];
    target[1] =  c2 * source[1] - c1 * source[3] + c0 * source[5];
    target[2] =  c3 * source[6] - c4 * source[8] + c3 * source[10];
    target[3] =  c5 * source[7] - c5 * source[9];
    target[4] =  c6 * source[11] - c7 * source[13] - c8 * source[0]
                  + c9 * source[2] - c8 * source[2] + c9 * source[4];
    target[5] =  c7 * source[12] - c6 * source[14] - c9 * source[1]
                  + c8 * source[3] - c9 * source[3] + c8 * source[5];
    target[6] =  c10 * source[15] - c10 * source[17] - c11 * source[6]
                  + c11 * source[8] - c11 * source[8] + c11 * source[10];
    target[7] =  c12 * source[16] - c10 * source[7] - c10 * source[9];
    target[8] =  c13 * source[18] - c14 * source[11] - c14 * source[13]
                  + c15 * source[0] + c16 * source[2] + c15 * source[4];
    target[9] =  c13 * source[19] - c14 * source[12] - c14 * source[14]
                  + c15 * source[1] + c16 * source[3] + c15 * source[5];
    target[10] =  source[20] - c17 * source[15] - c17 * source[17]
                  + c18 * source[6] + c19 * source[8] + c18 * source[10];
  }
}

