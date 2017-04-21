//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_70.cc
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

#ifdef COMPILE_J_ORB
#include <src/integral/carsphlist.h>
#include <algorithm>

using namespace std;
using namespace bagel;


void CarSphList::carsph_70(const int nloop, const double* source, double* target) {
  const double c9 = 56.995065575889988;
  const double c7 = 48.436491924993909;
  const double c17 = 37.996710383926661;
  const double c5 = 36.327368943745434;
  const double c20 = 34.369317712168801;
  const double c32 = 32.403703492039298;
  const double c10 = 28.497532787944994;
  const double c41 = 26.25;
  const double c22 = 25.776988284126599;
  const double c2 = 22.654094725071229;
  const double c35 = 19.843134832984429;
  const double c31 = 19.442222095223581;
  const double c16 = 17.098519672766997;
  const double c28 = 16.201851746019649;
  const double c6 = 14.530947577498171;
  const double c1 = 13.592456835042736;
  const double c40 = 13.125;
  const double c33 = 12.151388809514739;
  const double c19 = 11.456439237389599;
  const double c18 = 11.399013115177997;
  const double c39 = 10.5;
  const double c36 = 9.9215674164922145;
  const double c27 = 9.7211110476117906;
  const double c14 = 9.4991775959816653;
  const double c21 = 8.5923294280422002;
  const double c43 = 6.5625;
  const double c30 = 6.0756944047573693;
  const double c8 = 5.6995065575889985;
  const double c34 = 5.2915026221291814;
  const double c12 = 4.7495887979908327;
  const double c3 = 4.5308189450142455;
  const double c29 = 3.0378472023786847;
  const double c15 = 2.8497532787944992;
  const double c26 = 2.5776988284126601;
  const double c4 = 2.4218245962496954;
  const double c13 = 2.3747943989954163;
  const double c42 = 2.1875;
  const double c24 = 1.28884941420633;
  const double c38 = 1.2401959270615268;
  const double c25 = 0.85923294280422002;
  const double c0 = 0.64725984928774938;
  const double c11 = 0.47495887979908324;
  const double c23 = 0.42961647140211001;
  const double c37 = 0.41339864235384227;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 15, source += 36) {
    target[0] =  c0 * source[0] - c1 * source[2] + c2 * source[4]
                  - c3 * source[6];
    target[1] =  c3 * source[1] - c2 * source[3] + c1 * source[5]
                  - c0 * source[7];
    target[2] =  c4 * source[8] - c5 * source[10] + c5 * source[12]
                  - c4 * source[14];
    target[3] =  c6 * source[9] - c7 * source[11] + c6 * source[13];
    target[4] =  c8 * source[15] - c9 * source[17] + c10 * source[19]
                  - c11 * source[0] + c12 * source[2] - c13 * source[4]
                  - c11 * source[2] + c12 * source[4] - c13 * source[6];
    target[5] =  c10 * source[16] - c9 * source[18] + c8 * source[20]
                  - c13 * source[1] + c12 * source[3] - c11 * source[5]
                  - c13 * source[3] + c12 * source[5] - c11 * source[7];
    target[6] =  c14 * source[21] - c9 * source[23] + c14 * source[25]
                  - c15 * source[8] + c16 * source[10] - c15 * source[12]
                  - c15 * source[10] + c16 * source[12] - c15 * source[14];
    target[7] =  c17 * source[22] - c17 * source[24] - c18 * source[9]
                  + c18 * source[11] - c18 * source[11] + c18 * source[13];
    target[8] =  c19 * source[26] - c20 * source[28] - c21 * source[15]
                  + c22 * source[17] - c21 * source[17] + c22 * source[19]
                  + c23 * source[0] - c24 * source[2] + c25 * source[2]
                  - c26 * source[4] + c23 * source[4] - c24 * source[6];
    target[9] =  c20 * source[27] - c19 * source[29] - c22 * source[16]
                  + c21 * source[18] - c22 * source[18] + c21 * source[20]
                  + c24 * source[1] - c23 * source[3] + c26 * source[3]
                  - c25 * source[5] + c24 * source[5] - c23 * source[7];
    target[10] =  c27 * source[30] - c27 * source[32] - c28 * source[21]
                  + c28 * source[23] - c28 * source[23] + c28 * source[25]
                  + c29 * source[8] - c29 * source[10] + c30 * source[10]
                  - c30 * source[12] + c29 * source[12] - c29 * source[14];
    target[11] =  c31 * source[31] - c32 * source[22] - c32 * source[24]
                  + c30 * source[9] + c33 * source[11] + c30 * source[13];
    target[12] =  c34 * source[33] - c35 * source[26] - c35 * source[28]
                  + c36 * source[15] + c35 * source[17] + c36 * source[19]
                  - c37 * source[0] - c38 * source[2] - c38 * source[4]
                  - c37 * source[6];
    target[13] =  c34 * source[34] - c35 * source[27] - c35 * source[29]
                  + c36 * source[16] + c35 * source[18] + c36 * source[20]
                  - c37 * source[1] - c38 * source[3] - c38 * source[5]
                  - c37 * source[7];
    target[14] =  source[35] - c39 * source[30] - c39 * source[32]
                  + c40 * source[21] + c41 * source[23] + c40 * source[25]
                  - c42 * source[8] - c43 * source[10] - c43 * source[12]
                  - c42 * source[14];
  }
}

void CCarSphList::carsph_70(const int nloop, const complex<double>* source, complex<double>* target) {
  const double c9 = 56.995065575889988;
  const double c7 = 48.436491924993909;
  const double c17 = 37.996710383926661;
  const double c5 = 36.327368943745434;
  const double c20 = 34.369317712168801;
  const double c32 = 32.403703492039298;
  const double c10 = 28.497532787944994;
  const double c41 = 26.25;
  const double c22 = 25.776988284126599;
  const double c2 = 22.654094725071229;
  const double c35 = 19.843134832984429;
  const double c31 = 19.442222095223581;
  const double c16 = 17.098519672766997;
  const double c28 = 16.201851746019649;
  const double c6 = 14.530947577498171;
  const double c1 = 13.592456835042736;
  const double c40 = 13.125;
  const double c33 = 12.151388809514739;
  const double c19 = 11.456439237389599;
  const double c18 = 11.399013115177997;
  const double c39 = 10.5;
  const double c36 = 9.9215674164922145;
  const double c27 = 9.7211110476117906;
  const double c14 = 9.4991775959816653;
  const double c21 = 8.5923294280422002;
  const double c43 = 6.5625;
  const double c30 = 6.0756944047573693;
  const double c8 = 5.6995065575889985;
  const double c34 = 5.2915026221291814;
  const double c12 = 4.7495887979908327;
  const double c3 = 4.5308189450142455;
  const double c29 = 3.0378472023786847;
  const double c15 = 2.8497532787944992;
  const double c26 = 2.5776988284126601;
  const double c4 = 2.4218245962496954;
  const double c13 = 2.3747943989954163;
  const double c42 = 2.1875;
  const double c24 = 1.28884941420633;
  const double c38 = 1.2401959270615268;
  const double c25 = 0.85923294280422002;
  const double c0 = 0.64725984928774938;
  const double c11 = 0.47495887979908324;
  const double c23 = 0.42961647140211001;
  const double c37 = 0.41339864235384227;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 15, source += 36) {
    target[0] =  c0 * source[0] - c1 * source[2] + c2 * source[4]
                  - c3 * source[6];
    target[1] =  c3 * source[1] - c2 * source[3] + c1 * source[5]
                  - c0 * source[7];
    target[2] =  c4 * source[8] - c5 * source[10] + c5 * source[12]
                  - c4 * source[14];
    target[3] =  c6 * source[9] - c7 * source[11] + c6 * source[13];
    target[4] =  c8 * source[15] - c9 * source[17] + c10 * source[19]
                  - c11 * source[0] + c12 * source[2] - c13 * source[4]
                  - c11 * source[2] + c12 * source[4] - c13 * source[6];
    target[5] =  c10 * source[16] - c9 * source[18] + c8 * source[20]
                  - c13 * source[1] + c12 * source[3] - c11 * source[5]
                  - c13 * source[3] + c12 * source[5] - c11 * source[7];
    target[6] =  c14 * source[21] - c9 * source[23] + c14 * source[25]
                  - c15 * source[8] + c16 * source[10] - c15 * source[12]
                  - c15 * source[10] + c16 * source[12] - c15 * source[14];
    target[7] =  c17 * source[22] - c17 * source[24] - c18 * source[9]
                  + c18 * source[11] - c18 * source[11] + c18 * source[13];
    target[8] =  c19 * source[26] - c20 * source[28] - c21 * source[15]
                  + c22 * source[17] - c21 * source[17] + c22 * source[19]
                  + c23 * source[0] - c24 * source[2] + c25 * source[2]
                  - c26 * source[4] + c23 * source[4] - c24 * source[6];
    target[9] =  c20 * source[27] - c19 * source[29] - c22 * source[16]
                  + c21 * source[18] - c22 * source[18] + c21 * source[20]
                  + c24 * source[1] - c23 * source[3] + c26 * source[3]
                  - c25 * source[5] + c24 * source[5] - c23 * source[7];
    target[10] =  c27 * source[30] - c27 * source[32] - c28 * source[21]
                  + c28 * source[23] - c28 * source[23] + c28 * source[25]
                  + c29 * source[8] - c29 * source[10] + c30 * source[10]
                  - c30 * source[12] + c29 * source[12] - c29 * source[14];
    target[11] =  c31 * source[31] - c32 * source[22] - c32 * source[24]
                  + c30 * source[9] + c33 * source[11] + c30 * source[13];
    target[12] =  c34 * source[33] - c35 * source[26] - c35 * source[28]
                  + c36 * source[15] + c35 * source[17] + c36 * source[19]
                  - c37 * source[0] - c38 * source[2] - c38 * source[4]
                  - c37 * source[6];
    target[13] =  c34 * source[34] - c35 * source[27] - c35 * source[29]
                  + c36 * source[16] + c35 * source[18] + c36 * source[20]
                  - c37 * source[1] - c38 * source[3] - c38 * source[5]
                  - c37 * source[7];
    target[14] =  source[35] - c39 * source[30] - c39 * source[32]
                  + c40 * source[21] + c41 * source[23] + c40 * source[25]
                  - c42 * source[8] - c43 * source[10] - c43 * source[12]
                  - c42 * source[14];
  }
}

#endif
