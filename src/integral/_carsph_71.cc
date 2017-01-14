//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _carsph_71.cc
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


void CarSphList::carsph_71(const int nloop, const double* source, double* target) {
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
  for (int iloop = 0; iloop != nloop; ++iloop, target += 45, source += 108) {
    target[0] =  c0 * source[0] - c1 * source[6] + c2 * source[12]
                  - c3 * source[18];
    target[1] =  c0 * source[1] - c1 * source[7] + c2 * source[13]
                  - c3 * source[19];
    target[2] =  c0 * source[2] - c1 * source[8] + c2 * source[14]
                  - c3 * source[20];
    target[3] =  c3 * source[3] - c2 * source[9] + c1 * source[15]
                  - c0 * source[21];
    target[4] =  c3 * source[4] - c2 * source[10] + c1 * source[16]
                  - c0 * source[22];
    target[5] =  c3 * source[5] - c2 * source[11] + c1 * source[17]
                  - c0 * source[23];
    target[6] =  c4 * source[24] - c5 * source[30] + c5 * source[36]
                  - c4 * source[42];
    target[7] =  c4 * source[25] - c5 * source[31] + c5 * source[37]
                  - c4 * source[43];
    target[8] =  c4 * source[26] - c5 * source[32] + c5 * source[38]
                  - c4 * source[44];
    target[9] =  c6 * source[27] - c7 * source[33] + c6 * source[39];
    target[10] =  c6 * source[28] - c7 * source[34] + c6 * source[40];
    target[11] =  c6 * source[29] - c7 * source[35] + c6 * source[41];
    target[12] =  c8 * source[45] - c9 * source[51] + c10 * source[57]
                  - c11 * source[0] + c12 * source[6] - c13 * source[12]
                  - c11 * source[6] + c12 * source[12] - c13 * source[18];
    target[13] =  c8 * source[46] - c9 * source[52] + c10 * source[58]
                  - c11 * source[1] + c12 * source[7] - c13 * source[13]
                  - c11 * source[7] + c12 * source[13] - c13 * source[19];
    target[14] =  c8 * source[47] - c9 * source[53] + c10 * source[59]
                  - c11 * source[2] + c12 * source[8] - c13 * source[14]
                  - c11 * source[8] + c12 * source[14] - c13 * source[20];
    target[15] =  c10 * source[48] - c9 * source[54] + c8 * source[60]
                  - c13 * source[3] + c12 * source[9] - c11 * source[15]
                  - c13 * source[9] + c12 * source[15] - c11 * source[21];
    target[16] =  c10 * source[49] - c9 * source[55] + c8 * source[61]
                  - c13 * source[4] + c12 * source[10] - c11 * source[16]
                  - c13 * source[10] + c12 * source[16] - c11 * source[22];
    target[17] =  c10 * source[50] - c9 * source[56] + c8 * source[62]
                  - c13 * source[5] + c12 * source[11] - c11 * source[17]
                  - c13 * source[11] + c12 * source[17] - c11 * source[23];
    target[18] =  c14 * source[63] - c9 * source[69] + c14 * source[75]
                  - c15 * source[24] + c16 * source[30] - c15 * source[36]
                  - c15 * source[30] + c16 * source[36] - c15 * source[42];
    target[19] =  c14 * source[64] - c9 * source[70] + c14 * source[76]
                  - c15 * source[25] + c16 * source[31] - c15 * source[37]
                  - c15 * source[31] + c16 * source[37] - c15 * source[43];
    target[20] =  c14 * source[65] - c9 * source[71] + c14 * source[77]
                  - c15 * source[26] + c16 * source[32] - c15 * source[38]
                  - c15 * source[32] + c16 * source[38] - c15 * source[44];
    target[21] =  c17 * source[66] - c17 * source[72] - c18 * source[27]
                  + c18 * source[33] - c18 * source[33] + c18 * source[39];
    target[22] =  c17 * source[67] - c17 * source[73] - c18 * source[28]
                  + c18 * source[34] - c18 * source[34] + c18 * source[40];
    target[23] =  c17 * source[68] - c17 * source[74] - c18 * source[29]
                  + c18 * source[35] - c18 * source[35] + c18 * source[41];
    target[24] =  c19 * source[78] - c20 * source[84] - c21 * source[45]
                  + c22 * source[51] - c21 * source[51] + c22 * source[57]
                  + c23 * source[0] - c24 * source[6] + c25 * source[6]
                  - c26 * source[12] + c23 * source[12] - c24 * source[18];
    target[25] =  c19 * source[79] - c20 * source[85] - c21 * source[46]
                  + c22 * source[52] - c21 * source[52] + c22 * source[58]
                  + c23 * source[1] - c24 * source[7] + c25 * source[7]
                  - c26 * source[13] + c23 * source[13] - c24 * source[19];
    target[26] =  c19 * source[80] - c20 * source[86] - c21 * source[47]
                  + c22 * source[53] - c21 * source[53] + c22 * source[59]
                  + c23 * source[2] - c24 * source[8] + c25 * source[8]
                  - c26 * source[14] + c23 * source[14] - c24 * source[20];
    target[27] =  c20 * source[81] - c19 * source[87] - c22 * source[48]
                  + c21 * source[54] - c22 * source[54] + c21 * source[60]
                  + c24 * source[3] - c23 * source[9] + c26 * source[9]
                  - c25 * source[15] + c24 * source[15] - c23 * source[21];
    target[28] =  c20 * source[82] - c19 * source[88] - c22 * source[49]
                  + c21 * source[55] - c22 * source[55] + c21 * source[61]
                  + c24 * source[4] - c23 * source[10] + c26 * source[10]
                  - c25 * source[16] + c24 * source[16] - c23 * source[22];
    target[29] =  c20 * source[83] - c19 * source[89] - c22 * source[50]
                  + c21 * source[56] - c22 * source[56] + c21 * source[62]
                  + c24 * source[5] - c23 * source[11] + c26 * source[11]
                  - c25 * source[17] + c24 * source[17] - c23 * source[23];
    target[30] =  c27 * source[90] - c27 * source[96] - c28 * source[63]
                  + c28 * source[69] - c28 * source[69] + c28 * source[75]
                  + c29 * source[24] - c29 * source[30] + c30 * source[30]
                  - c30 * source[36] + c29 * source[36] - c29 * source[42];
    target[31] =  c27 * source[91] - c27 * source[97] - c28 * source[64]
                  + c28 * source[70] - c28 * source[70] + c28 * source[76]
                  + c29 * source[25] - c29 * source[31] + c30 * source[31]
                  - c30 * source[37] + c29 * source[37] - c29 * source[43];
    target[32] =  c27 * source[92] - c27 * source[98] - c28 * source[65]
                  + c28 * source[71] - c28 * source[71] + c28 * source[77]
                  + c29 * source[26] - c29 * source[32] + c30 * source[32]
                  - c30 * source[38] + c29 * source[38] - c29 * source[44];
    target[33] =  c31 * source[93] - c32 * source[66] - c32 * source[72]
                  + c30 * source[27] + c33 * source[33] + c30 * source[39];
    target[34] =  c31 * source[94] - c32 * source[67] - c32 * source[73]
                  + c30 * source[28] + c33 * source[34] + c30 * source[40];
    target[35] =  c31 * source[95] - c32 * source[68] - c32 * source[74]
                  + c30 * source[29] + c33 * source[35] + c30 * source[41];
    target[36] =  c34 * source[99] - c35 * source[78] - c35 * source[84]
                  + c36 * source[45] + c35 * source[51] + c36 * source[57]
                  - c37 * source[0] - c38 * source[6] - c38 * source[12]
                  - c37 * source[18];
    target[37] =  c34 * source[100] - c35 * source[79] - c35 * source[85]
                  + c36 * source[46] + c35 * source[52] + c36 * source[58]
                  - c37 * source[1] - c38 * source[7] - c38 * source[13]
                  - c37 * source[19];
    target[38] =  c34 * source[101] - c35 * source[80] - c35 * source[86]
                  + c36 * source[47] + c35 * source[53] + c36 * source[59]
                  - c37 * source[2] - c38 * source[8] - c38 * source[14]
                  - c37 * source[20];
    target[39] =  c34 * source[102] - c35 * source[81] - c35 * source[87]
                  + c36 * source[48] + c35 * source[54] + c36 * source[60]
                  - c37 * source[3] - c38 * source[9] - c38 * source[15]
                  - c37 * source[21];
    target[40] =  c34 * source[103] - c35 * source[82] - c35 * source[88]
                  + c36 * source[49] + c35 * source[55] + c36 * source[61]
                  - c37 * source[4] - c38 * source[10] - c38 * source[16]
                  - c37 * source[22];
    target[41] =  c34 * source[104] - c35 * source[83] - c35 * source[89]
                  + c36 * source[50] + c35 * source[56] + c36 * source[62]
                  - c37 * source[5] - c38 * source[11] - c38 * source[17]
                  - c37 * source[23];
    target[42] =  source[105] - c39 * source[90] - c39 * source[96]
                  + c40 * source[63] + c41 * source[69] + c40 * source[75]
                  - c42 * source[24] - c43 * source[30] - c43 * source[36]
                  - c42 * source[42];
    target[43] =  source[106] - c39 * source[91] - c39 * source[97]
                  + c40 * source[64] + c41 * source[70] + c40 * source[76]
                  - c42 * source[25] - c43 * source[31] - c43 * source[37]
                  - c42 * source[43];
    target[44] =  source[107] - c39 * source[92] - c39 * source[98]
                  + c40 * source[65] + c41 * source[71] + c40 * source[77]
                  - c42 * source[26] - c43 * source[32] - c43 * source[38]
                  - c42 * source[44];
  }
}

void CCarSphList::carsph_71(const int nloop, const complex<double>* source, complex<double>* target) {
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
  for (int iloop = 0; iloop != nloop; ++iloop, target += 45, source += 108) {
    target[0] =  c0 * source[0] - c1 * source[6] + c2 * source[12]
                  - c3 * source[18];
    target[1] =  c0 * source[1] - c1 * source[7] + c2 * source[13]
                  - c3 * source[19];
    target[2] =  c0 * source[2] - c1 * source[8] + c2 * source[14]
                  - c3 * source[20];
    target[3] =  c3 * source[3] - c2 * source[9] + c1 * source[15]
                  - c0 * source[21];
    target[4] =  c3 * source[4] - c2 * source[10] + c1 * source[16]
                  - c0 * source[22];
    target[5] =  c3 * source[5] - c2 * source[11] + c1 * source[17]
                  - c0 * source[23];
    target[6] =  c4 * source[24] - c5 * source[30] + c5 * source[36]
                  - c4 * source[42];
    target[7] =  c4 * source[25] - c5 * source[31] + c5 * source[37]
                  - c4 * source[43];
    target[8] =  c4 * source[26] - c5 * source[32] + c5 * source[38]
                  - c4 * source[44];
    target[9] =  c6 * source[27] - c7 * source[33] + c6 * source[39];
    target[10] =  c6 * source[28] - c7 * source[34] + c6 * source[40];
    target[11] =  c6 * source[29] - c7 * source[35] + c6 * source[41];
    target[12] =  c8 * source[45] - c9 * source[51] + c10 * source[57]
                  - c11 * source[0] + c12 * source[6] - c13 * source[12]
                  - c11 * source[6] + c12 * source[12] - c13 * source[18];
    target[13] =  c8 * source[46] - c9 * source[52] + c10 * source[58]
                  - c11 * source[1] + c12 * source[7] - c13 * source[13]
                  - c11 * source[7] + c12 * source[13] - c13 * source[19];
    target[14] =  c8 * source[47] - c9 * source[53] + c10 * source[59]
                  - c11 * source[2] + c12 * source[8] - c13 * source[14]
                  - c11 * source[8] + c12 * source[14] - c13 * source[20];
    target[15] =  c10 * source[48] - c9 * source[54] + c8 * source[60]
                  - c13 * source[3] + c12 * source[9] - c11 * source[15]
                  - c13 * source[9] + c12 * source[15] - c11 * source[21];
    target[16] =  c10 * source[49] - c9 * source[55] + c8 * source[61]
                  - c13 * source[4] + c12 * source[10] - c11 * source[16]
                  - c13 * source[10] + c12 * source[16] - c11 * source[22];
    target[17] =  c10 * source[50] - c9 * source[56] + c8 * source[62]
                  - c13 * source[5] + c12 * source[11] - c11 * source[17]
                  - c13 * source[11] + c12 * source[17] - c11 * source[23];
    target[18] =  c14 * source[63] - c9 * source[69] + c14 * source[75]
                  - c15 * source[24] + c16 * source[30] - c15 * source[36]
                  - c15 * source[30] + c16 * source[36] - c15 * source[42];
    target[19] =  c14 * source[64] - c9 * source[70] + c14 * source[76]
                  - c15 * source[25] + c16 * source[31] - c15 * source[37]
                  - c15 * source[31] + c16 * source[37] - c15 * source[43];
    target[20] =  c14 * source[65] - c9 * source[71] + c14 * source[77]
                  - c15 * source[26] + c16 * source[32] - c15 * source[38]
                  - c15 * source[32] + c16 * source[38] - c15 * source[44];
    target[21] =  c17 * source[66] - c17 * source[72] - c18 * source[27]
                  + c18 * source[33] - c18 * source[33] + c18 * source[39];
    target[22] =  c17 * source[67] - c17 * source[73] - c18 * source[28]
                  + c18 * source[34] - c18 * source[34] + c18 * source[40];
    target[23] =  c17 * source[68] - c17 * source[74] - c18 * source[29]
                  + c18 * source[35] - c18 * source[35] + c18 * source[41];
    target[24] =  c19 * source[78] - c20 * source[84] - c21 * source[45]
                  + c22 * source[51] - c21 * source[51] + c22 * source[57]
                  + c23 * source[0] - c24 * source[6] + c25 * source[6]
                  - c26 * source[12] + c23 * source[12] - c24 * source[18];
    target[25] =  c19 * source[79] - c20 * source[85] - c21 * source[46]
                  + c22 * source[52] - c21 * source[52] + c22 * source[58]
                  + c23 * source[1] - c24 * source[7] + c25 * source[7]
                  - c26 * source[13] + c23 * source[13] - c24 * source[19];
    target[26] =  c19 * source[80] - c20 * source[86] - c21 * source[47]
                  + c22 * source[53] - c21 * source[53] + c22 * source[59]
                  + c23 * source[2] - c24 * source[8] + c25 * source[8]
                  - c26 * source[14] + c23 * source[14] - c24 * source[20];
    target[27] =  c20 * source[81] - c19 * source[87] - c22 * source[48]
                  + c21 * source[54] - c22 * source[54] + c21 * source[60]
                  + c24 * source[3] - c23 * source[9] + c26 * source[9]
                  - c25 * source[15] + c24 * source[15] - c23 * source[21];
    target[28] =  c20 * source[82] - c19 * source[88] - c22 * source[49]
                  + c21 * source[55] - c22 * source[55] + c21 * source[61]
                  + c24 * source[4] - c23 * source[10] + c26 * source[10]
                  - c25 * source[16] + c24 * source[16] - c23 * source[22];
    target[29] =  c20 * source[83] - c19 * source[89] - c22 * source[50]
                  + c21 * source[56] - c22 * source[56] + c21 * source[62]
                  + c24 * source[5] - c23 * source[11] + c26 * source[11]
                  - c25 * source[17] + c24 * source[17] - c23 * source[23];
    target[30] =  c27 * source[90] - c27 * source[96] - c28 * source[63]
                  + c28 * source[69] - c28 * source[69] + c28 * source[75]
                  + c29 * source[24] - c29 * source[30] + c30 * source[30]
                  - c30 * source[36] + c29 * source[36] - c29 * source[42];
    target[31] =  c27 * source[91] - c27 * source[97] - c28 * source[64]
                  + c28 * source[70] - c28 * source[70] + c28 * source[76]
                  + c29 * source[25] - c29 * source[31] + c30 * source[31]
                  - c30 * source[37] + c29 * source[37] - c29 * source[43];
    target[32] =  c27 * source[92] - c27 * source[98] - c28 * source[65]
                  + c28 * source[71] - c28 * source[71] + c28 * source[77]
                  + c29 * source[26] - c29 * source[32] + c30 * source[32]
                  - c30 * source[38] + c29 * source[38] - c29 * source[44];
    target[33] =  c31 * source[93] - c32 * source[66] - c32 * source[72]
                  + c30 * source[27] + c33 * source[33] + c30 * source[39];
    target[34] =  c31 * source[94] - c32 * source[67] - c32 * source[73]
                  + c30 * source[28] + c33 * source[34] + c30 * source[40];
    target[35] =  c31 * source[95] - c32 * source[68] - c32 * source[74]
                  + c30 * source[29] + c33 * source[35] + c30 * source[41];
    target[36] =  c34 * source[99] - c35 * source[78] - c35 * source[84]
                  + c36 * source[45] + c35 * source[51] + c36 * source[57]
                  - c37 * source[0] - c38 * source[6] - c38 * source[12]
                  - c37 * source[18];
    target[37] =  c34 * source[100] - c35 * source[79] - c35 * source[85]
                  + c36 * source[46] + c35 * source[52] + c36 * source[58]
                  - c37 * source[1] - c38 * source[7] - c38 * source[13]
                  - c37 * source[19];
    target[38] =  c34 * source[101] - c35 * source[80] - c35 * source[86]
                  + c36 * source[47] + c35 * source[53] + c36 * source[59]
                  - c37 * source[2] - c38 * source[8] - c38 * source[14]
                  - c37 * source[20];
    target[39] =  c34 * source[102] - c35 * source[81] - c35 * source[87]
                  + c36 * source[48] + c35 * source[54] + c36 * source[60]
                  - c37 * source[3] - c38 * source[9] - c38 * source[15]
                  - c37 * source[21];
    target[40] =  c34 * source[103] - c35 * source[82] - c35 * source[88]
                  + c36 * source[49] + c35 * source[55] + c36 * source[61]
                  - c37 * source[4] - c38 * source[10] - c38 * source[16]
                  - c37 * source[22];
    target[41] =  c34 * source[104] - c35 * source[83] - c35 * source[89]
                  + c36 * source[50] + c35 * source[56] + c36 * source[62]
                  - c37 * source[5] - c38 * source[11] - c38 * source[17]
                  - c37 * source[23];
    target[42] =  source[105] - c39 * source[90] - c39 * source[96]
                  + c40 * source[63] + c41 * source[69] + c40 * source[75]
                  - c42 * source[24] - c43 * source[30] - c43 * source[36]
                  - c42 * source[42];
    target[43] =  source[106] - c39 * source[91] - c39 * source[97]
                  + c40 * source[64] + c41 * source[70] + c40 * source[76]
                  - c42 * source[25] - c43 * source[31] - c43 * source[37]
                  - c42 * source[43];
    target[44] =  source[107] - c39 * source[92] - c39 * source[98]
                  + c40 * source[65] + c41 * source[71] + c40 * source[77]
                  - c42 * source[26] - c43 * source[32] - c43 * source[38]
                  - c42 * source[44];
  }
}

#endif
