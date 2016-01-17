//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: vrr_template.cc
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


#include <src/integral/rys/eribatch.h>

using namespace std;
using namespace bagel;

void ERIBatch::perform_VRR3() {

// (4, 0, 0, 0) and (0, 0, 0, 4) hasn't been rewritten (using the generic code);
// (3, 0, 1, 0) and (1, 0, 3, 0) hasn't been rewritten yet
// (2, 1, 1, 0) and (1, 0, 2, 1) hasn't been rewritten yet

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);

  const unsigned int ang0 = basisinfo_[0]->angular_number();
  const unsigned int ang1 = basisinfo_[1]->angular_number();
  const unsigned int ang2 = basisinfo_[2]->angular_number();
  const unsigned int ang3 = basisinfo_[3]->angular_number();
  const unsigned int ang = (((((ang0 << 5) + ang1) << 5) + ang2) << 5) + ang3; // offsets are 2^5 = 32

  if (ang == 65600) { // (2, 0, 2, 0)
    double x[27];
    double y[27];
    double z[27];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;

      const double oxp2 = 0.5 / xp_[ii];
      const double opq = 1.0 / (xp_[ii] + xq_[ii]);
      const double xqopq = xq_[ii] * opq;
      const double c00i0x = P_[ii3]     - ax;
      const double c00i0y = P_[ii3 + 1] - ay;
      const double c00i0z = P_[ii3 + 2] - az;
      const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
      const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
      const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
      const double b10i0 = xqopq * oxp2;
      const double b10_0 = oxp2 - b10i0 * roots_[offset];
      const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
      const double b10_2 = oxp2 - b10i0 * roots_[offset + 2];
      const double oxq2 = 0.5 / xq_[ii];
      const double xpopq = xp_[ii] * opq;
      const double d00i0x = Q_[ii3] - cx;
      const double d00i0y = Q_[ii3 + 1] - cy;
      const double d00i0z = Q_[ii3 + 2] - cz;
      const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
      const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
      const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
      const double b00int = 0.5 * opq;
      const double b00_0 = b00int * roots_[offset];
      const double b00_1 = b00int * roots_[offset + 1];
      const double b00_2 = b00int * roots_[offset + 2];
      const double b01i0 = xpopq * oxq2;
      const double b01_0 = oxq2 - b01i0 * roots_[offset];
      const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
      const double b01_2 = oxq2 - b01i0 * roots_[offset + 2];

      x[3]  = c00i0x - c00i1x * roots_[offset];
      x[4]  = c00i0x - c00i1x * roots_[offset + 1];
      x[5]  = c00i0x - c00i1x * roots_[offset + 2];
      x[6]  = x[3] * x[3] + b10_0;
      x[7]  = x[4] * x[4] + b10_1;
      x[8]  = x[5] * x[5] + b10_2;
      x[9]  = d00i0x + d00i1x * roots_[offset];
      x[10] = d00i0x + d00i1x * roots_[offset + 1];
      x[11] = d00i0x + d00i1x * roots_[offset + 2];
      x[12] = x[3] * x[9]  + b00_0;
      x[13] = x[4] * x[10] + b00_1;
      x[14] = x[5] * x[11] + b00_2;
      x[15] = x[3] * (x[12] + b00_0) + b10_0 * x[9];
      x[16] = x[4] * (x[13] + b00_1) + b10_1 * x[10];
      x[17] = x[5] * (x[14] + b00_2) + b10_2 * x[11];
      x[18] = x[9]  * x[9]  + b01_0;
      x[19] = x[10] * x[10] + b01_1;
      x[20] = x[11] * x[11] + b01_2;
      x[21] = x[3] * x[18] + 2.0 * b00_0 * x[9];
      x[22] = x[4] * x[19] + 2.0 * b00_1 * x[10];
      x[23] = x[5] * x[20] + 2.0 * b00_2 * x[11];
      x[24] = x[3] * x[21] + 2.0 * b00_0 * x[12] + b10_0 * x[18];
      x[25] = x[4] * x[22] + 2.0 * b00_1 * x[13] + b10_1 * x[19];
      x[26] = x[5] * x[23] + 2.0 * b00_2 * x[14] + b10_2 * x[20];

      y[3]  = c00i0y - c00i1y * roots_[offset];
      y[4]  = c00i0y - c00i1y * roots_[offset + 1];
      y[5]  = c00i0y - c00i1y * roots_[offset + 2];
      y[6]  = y[3] * y[3] + b10_0;
      y[7]  = y[4] * y[4] + b10_1;
      y[8]  = y[5] * y[5] + b10_2;
      y[9]  = d00i0y + d00i1y * roots_[offset];
      y[10] = d00i0y + d00i1y * roots_[offset + 1];
      y[11] = d00i0y + d00i1y * roots_[offset + 2];
      y[12] = y[3] * y[9]  + b00_0;
      y[13] = y[4] * y[10] + b00_1;
      y[14] = y[5] * y[11] + b00_2;
      y[15] = y[3] * (y[12] + b00_0) + b10_0 * y[9];
      y[16] = y[4] * (y[13] + b00_1) + b10_1 * y[10];
      y[17] = y[5] * (y[14] + b00_2) + b10_2 * y[11];
      y[18] = y[9]  * y[9]  + b01_0;
      y[19] = y[10] * y[10] + b01_1;
      y[20] = y[11] * y[11] + b01_2;
      y[21] = y[3] * y[18] + 2.0 * b00_0 * y[9];
      y[22] = y[4] * y[19] + 2.0 * b00_1 * y[10];
      y[23] = y[5] * y[20] + 2.0 * b00_2 * y[11];
      y[24] = y[3] * y[21] + 2.0 * b00_0 * y[12] + b10_0 * y[18];
      y[25] = y[4] * y[22] + 2.0 * b00_1 * y[13] + b10_1 * y[19];
      y[26] = y[5] * y[23] + 2.0 * b00_2 * y[14] + b10_2 * y[20];

      z[3]  = c00i0z - c00i1z * roots_[offset];
      z[4]  = c00i0z - c00i1z * roots_[offset + 1];
      z[5]  = c00i0z - c00i1z * roots_[offset + 2];
      z[6]  = z[3] * z[3] + b10_0;
      z[7]  = z[4] * z[4] + b10_1;
      z[8]  = z[5] * z[5] + b10_2;
      z[9]  = d00i0z + d00i1z * roots_[offset];
      z[10] = d00i0z + d00i1z * roots_[offset + 1];
      z[11] = d00i0z + d00i1z * roots_[offset + 2];
      z[12] = z[3] * z[9]  + b00_0;
      z[13] = z[4] * z[10] + b00_1;
      z[14] = z[5] * z[11] + b00_2;
      z[15] = z[3] * (z[12] + b00_0) + b10_0 * z[9];
      z[16] = z[4] * (z[13] + b00_1) + b10_1 * z[10];
      z[17] = z[5] * (z[14] + b00_2) + b10_2 * z[11];
      z[18] = z[9]  * z[9]  + b01_0;
      z[19] = z[10] * z[10] + b01_1;
      z[20] = z[11] * z[11] + b01_2;
      z[21] = z[3] * z[18] + 2.0 * b00_0 * z[9];
      z[22] = z[4] * z[19] + 2.0 * b00_1 * z[10];
      z[23] = z[5] * z[20] + 2.0 * b00_2 * z[11];
      z[24] = z[3] * z[21] + 2.0 * b00_0 * z[12] + b10_0 * z[18];
      z[25] = z[4] * z[22] + 2.0 * b00_1 * z[13] + b10_1 * z[19];
      z[26] = z[5] * z[23] + 2.0 * b00_2 * z[14] + b10_2 * z[20];

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double scale2 = coeff_[ii] * weights_[offset + 2];

      const double x3s = x[3] * scale0;
      const double x4s = x[4] * scale1;
      const double x5s = x[5] * scale2;
      const double x6s = x[6] * scale0;
      const double x7s = x[7] * scale1;
      const double x8s = x[8] * scale2;
      const double x9s = x[9] * scale0;
      const double x10s = x[10] * scale1;
      const double x11s = x[11] * scale2;
      const double y3s = y[3] * scale0;
      const double y4s = y[4] * scale1;
      const double y5s = y[5] * scale2;
      const double y6s = y[6] * scale0;
      const double y7s = y[7] * scale1;
      const double y8s = y[8] * scale2;
      const double y9s = y[9] * scale0;
      const double y10s = y[10] * scale1;
      const double y11s = y[11] * scale2;
      const double z3s = z[3] * scale0;
      const double z4s = z[4] * scale1;
      const double z5s = z[5] * scale2;
      const double z6s = z[6] * scale0;
      const double z7s = z[7] * scale1;
      const double z8s = z[8] * scale2;
      const double z9s = z[9] * scale0;
      const double z10s = z[10] * scale1;
      const double z11s = z[11] * scale2;

      current_data[0]  = x[24] * scale0               + x[25] * scale1                + x[26] * scale2; // xx|xx
      current_data[1]  = x[21] * y3s                  + x[22] * y4s                   + x[23] * y5s; // xx|xy
      current_data[2]  = x[18] * y6s                  + x[19] * y7s                   + x[20] * y8s; // xx|yy
      current_data[3]  = x[21] * z3s                  + x[22] * z4s                   + x[23] * z5s; // xx|xz
      current_data[4]  = x[18] * y3s * z[3]           + x[19] * y4s * z[4]            + x[20] * y5s * z[5];// xx|yz
      current_data[5]  = x[18] * z6s                  + x[19] * z7s                   + x[20] * z8s; // xx|zz

      current_data[6]  = x[15] * y9s                  + x[16] * y10s                  + x[17] * y11s; // xy|xx
      current_data[7]  = x[12] * y[12] * scale0       + x[13] * y[13] * scale1        + x[14] * y[14] * scale2; // xy|xy
      current_data[8]  = x9s * y[15]                  + x10s * y[16]                  + x11s * y[17]; // xy|yy
      current_data[9]  = x[12] * y[9] * z3s           + x[13] * y[10] * z4s           + x[14] * y[11] * z5s; // xy|xz
      current_data[10] = x[9] * y[12] * z3s           + x[10] * y[13] * z4s           + x[11] * y[14] * z5s; // xy|yz
      current_data[11] = x[9] * y[9] * z6s            + x[10] * y[10] * z7s           + x[11] * y[11] * z8s; // xy|zz

      current_data[12] = y[18] * x6s                  + y[19] * x7s                   + y[20] * x8s; // yy|xx
      current_data[13] = y[21] * x3s                  + y[22] * x4s                   + y[23] * x5s; // yy|xy
      current_data[14] = y[24] * scale0               + y[25] * scale1                + y[26] * scale2; // yy|yy
      current_data[15] = y[18] * x[3] * z3s           + y[19] * x[4] * z4s            + y[20] * x[5] * z5s;// yy|xz
      current_data[16] = y[21] * z3s                  + y[22] * z4s                   + y[23] * z5s; // yy|yz
      current_data[17] = y[18] * z6s                  + y[19] * z7s                   + y[20] * z8s; // yy|zz

      current_data[18] = x[15] * z9s                  + x[16] * z10s                  + x[17] * z11s; // xz|xx
      current_data[19] = x[12] * z[9] * y3s           + x[13] * z[10] * y4s           + x[14] * z[11] * y5s; // xz|xy
      current_data[20] = x[9] * z[9] * y6s            + x[10] * z[10] * y7s           + x[11] * z[11] * y8s; // xz|yy
      current_data[21] = x[12] * z[12] * scale0       + x[13] * z[13] * scale1        + x[14] * z[14] * scale2; // xz|xz
      current_data[22] = x[9] * z[12] * y3s           + x[10] * z[13] * y4s           + x[11] * z[14] * y5s; // xz|yz
      current_data[23] = x9s * z[15]                  + x10s * z[16]                  + x11s * z[17]; // xz|zz

      current_data[24] = z[9] * y[9] * x6s            + z[10] * y[10] * x7s           + z[11] * y[11] * x8s; // yz|xx
      current_data[25] = z[9] * y[12] * x3s           + z[10] * y[13] * x4s           + z[11] * y[14] * x5s; // yz|xy
      current_data[26] = z9s * y[15]                  + z10s * y[16]                  + z11s * y[17]; // yz|yy
      current_data[27] = z[12] * y[9] * x3s           + z[13] * y[10] * x4s           + z[14] * y[11] * x5s; // yz|xz
      current_data[28] = z[12] * y[12] * scale0       + z[13] * y[13] * scale1        + z[14] * y[14] * scale2; // yz|yz
      current_data[29] = y9s * z[15]                  + y10s * z[16]                  + y11s * z[17]; // yz|zz

      current_data[30] = z[18] * x6s                  + z[19] * x7s                   + z[20] * x8s; // zz|xx
      current_data[31] = z[18] * y3s * x[3]           + z[19] * y4s * x[4]            + z[20] * y5s * x[5];// zz|xy
      current_data[32] = z[18] * y6s                  + z[19] * y7s                   + z[20] * y8s; // zz|yy
      current_data[33] = z[21] * x3s                  + z[22] * x4s                   + z[23] * x5s; // zz|xz
      current_data[34] = z[21] * y3s                  + z[22] * y4s                   + z[23] * y5s; // zz|yz
      current_data[35] = z[24] * scale0               + z[25] * scale1                + z[26] * scale2; // xx|xx

    }
  } else if (ang == 65569) { // (2, 0, 1, 1)
    double x[27];
    double y[27];
    double z[27];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;

      const double oxp2 = 0.5 / xp_[ii];
      const double opq = 1.0 / (xp_[ii] + xq_[ii]);
      const double xqopq = xq_[ii] * opq;
      const double c00i0x = P_[ii3]     - ax;
      const double c00i0y = P_[ii3 + 1] - ay;
      const double c00i0z = P_[ii3 + 2] - az;
      const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
      const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
      const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
      const double b10i0 = xqopq * oxp2;
      const double b10_0 = oxp2 - b10i0 * roots_[offset];
      const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
      const double b10_2 = oxp2 - b10i0 * roots_[offset + 2];
      const double oxq2 = 0.5 / xq_[ii];
      const double xpopq = xp_[ii] * opq;
      const double d00i0x = Q_[ii3] - cx;
      const double d00i0y = Q_[ii3 + 1] - cy;
      const double d00i0z = Q_[ii3 + 2] - cz;
      const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
      const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
      const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
      const double b00int = 0.5 * opq;
      const double b00_0 = b00int * roots_[offset];
      const double b00_1 = b00int * roots_[offset + 1];
      const double b00_2 = b00int * roots_[offset + 2];
      const double b01i0 = xpopq * oxq2;
      const double b01_0 = oxq2 - b01i0 * roots_[offset];
      const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
      const double b01_2 = oxq2 - b01i0 * roots_[offset + 2];

      x[3]  = c00i0x - c00i1x * roots_[offset];
      x[4]  = c00i0x - c00i1x * roots_[offset + 1];
      x[5]  = c00i0x - c00i1x * roots_[offset + 2];
      x[6]  = x[3] * x[3] + b10_0;
      x[7]  = x[4] * x[4] + b10_1;
      x[8]  = x[5] * x[5] + b10_2;
      x[9]  = d00i0x + d00i1x * roots_[offset];
      x[10] = d00i0x + d00i1x * roots_[offset + 1];
      x[11] = d00i0x + d00i1x * roots_[offset + 2];
      x[12] = x[3] * x[9]  + b00_0;
      x[13] = x[4] * x[10] + b00_1;
      x[14] = x[5] * x[11] + b00_2;
      x[15] = x[3] * (x[12] + b00_0) + b10_0 * x[9];
      x[16] = x[4] * (x[13] + b00_1) + b10_1 * x[10];
      x[17] = x[5] * (x[14] + b00_2) + b10_2 * x[11];
      x[18] = x[9]  * x[9]  + b01_0;
      x[19] = x[10] * x[10] + b01_1;
      x[20] = x[11] * x[11] + b01_2;
      x[21] = x[3] * x[18] + 2.0 * b00_0 * x[9];
      x[22] = x[4] * x[19] + 2.0 * b00_1 * x[10];
      x[23] = x[5] * x[20] + 2.0 * b00_2 * x[11];
      x[24] = x[3] * x[21] + 2.0 * b00_0 * x[12] + b10_0 * x[18];
      x[25] = x[4] * x[22] + 2.0 * b00_1 * x[13] + b10_1 * x[19];
      x[26] = x[5] * x[23] + 2.0 * b00_2 * x[14] + b10_2 * x[20];

      y[3]  = c00i0y - c00i1y * roots_[offset];
      y[4]  = c00i0y - c00i1y * roots_[offset + 1];
      y[5]  = c00i0y - c00i1y * roots_[offset + 2];
      y[6]  = y[3] * y[3] + b10_0;
      y[7]  = y[4] * y[4] + b10_1;
      y[8]  = y[5] * y[5] + b10_2;
      y[9]  = d00i0y + d00i1y * roots_[offset];
      y[10] = d00i0y + d00i1y * roots_[offset + 1];
      y[11] = d00i0y + d00i1y * roots_[offset + 2];
      y[12] = y[3] * y[9]  + b00_0;
      y[13] = y[4] * y[10] + b00_1;
      y[14] = y[5] * y[11] + b00_2;
      y[15] = y[3] * (y[12] + b00_0) + b10_0 * y[9];
      y[16] = y[4] * (y[13] + b00_1) + b10_1 * y[10];
      y[17] = y[5] * (y[14] + b00_2) + b10_2 * y[11];
      y[18] = y[9]  * y[9]  + b01_0;
      y[19] = y[10] * y[10] + b01_1;
      y[20] = y[11] * y[11] + b01_2;
      y[21] = y[3] * y[18] + 2.0 * b00_0 * y[9];
      y[22] = y[4] * y[19] + 2.0 * b00_1 * y[10];
      y[23] = y[5] * y[20] + 2.0 * b00_2 * y[11];
      y[24] = y[3] * y[21] + 2.0 * b00_0 * y[12] + b10_0 * y[18];
      y[25] = y[4] * y[22] + 2.0 * b00_1 * y[13] + b10_1 * y[19];
      y[26] = y[5] * y[23] + 2.0 * b00_2 * y[14] + b10_2 * y[20];

      z[3]  = c00i0z - c00i1z * roots_[offset];
      z[4]  = c00i0z - c00i1z * roots_[offset + 1];
      z[5]  = c00i0z - c00i1z * roots_[offset + 2];
      z[6]  = z[3] * z[3] + b10_0;
      z[7]  = z[4] * z[4] + b10_1;
      z[8]  = z[5] * z[5] + b10_2;
      z[9]  = d00i0z + d00i1z * roots_[offset];
      z[10] = d00i0z + d00i1z * roots_[offset + 1];
      z[11] = d00i0z + d00i1z * roots_[offset + 2];
      z[12] = z[3] * z[9]  + b00_0;
      z[13] = z[4] * z[10] + b00_1;
      z[14] = z[5] * z[11] + b00_2;
      z[15] = z[3] * (z[12] + b00_0) + b10_0 * z[9];
      z[16] = z[4] * (z[13] + b00_1) + b10_1 * z[10];
      z[17] = z[5] * (z[14] + b00_2) + b10_2 * z[11];
      z[18] = z[9]  * z[9]  + b01_0;
      z[19] = z[10] * z[10] + b01_1;
      z[20] = z[11] * z[11] + b01_2;
      z[21] = z[3] * z[18] + 2.0 * b00_0 * z[9];
      z[22] = z[4] * z[19] + 2.0 * b00_1 * z[10];
      z[23] = z[5] * z[20] + 2.0 * b00_2 * z[11];
      z[24] = z[3] * z[21] + 2.0 * b00_0 * z[12] + b10_0 * z[18];
      z[25] = z[4] * z[22] + 2.0 * b00_1 * z[13] + b10_1 * z[19];
      z[26] = z[5] * z[23] + 2.0 * b00_2 * z[14] + b10_2 * z[20];

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double scale2 = coeff_[ii] * weights_[offset + 2];

      const double x3s = x[3] * scale0;
      const double x4s = x[4] * scale1;
      const double x5s = x[5] * scale2;
      const double x6s = x[6] * scale0;
      const double x7s = x[7] * scale1;
      const double x8s = x[8] * scale2;
      const double x9s = x[9] * scale0;
      const double x10s = x[10] * scale1;
      const double x11s = x[11] * scale2;
      const double y3s = y[3] * scale0;
      const double y4s = y[4] * scale1;
      const double y5s = y[5] * scale2;
      const double y6s = y[6] * scale0;
      const double y7s = y[7] * scale1;
      const double y8s = y[8] * scale2;
      const double y9s = y[9] * scale0;
      const double y10s = y[10] * scale1;
      const double y11s = y[11] * scale2;
      const double z3s = z[3] * scale0;
      const double z4s = z[4] * scale1;
      const double z5s = z[5] * scale2;
      const double z6s = z[6] * scale0;
      const double z7s = z[7] * scale1;
      const double z8s = z[8] * scale2;
      const double z9s = z[9] * scale0;
      const double z10s = z[10] * scale1;
      const double z11s = z[11] * scale2;

      current_data[0]  = x[15] * scale0               + x[16] * scale1                + x[17] * scale2; // x|xx
      current_data[1]  = x[12] * y3s                  + x[13] * y4s                   + x[14] * y5s; // x|xy
      current_data[2]  = x[9] * y6s                   + x[10] * y7s                   + x[11] * y8s; // x|yy
      current_data[3]  = x[12] * z3s                  + x[13] * z4s                   + x[14] * z5s; // x|xz
      current_data[4]  = x[9] * y3s * z[3]            + x[10] * y4s * z[4]            + x[11] * y5s * z[5]; // x|yz
      current_data[5]  = x[9] * z6s                   + x[10] * z7s                   + x[11] * z8s; // x|zz

      current_data[6]  = y[9] * x6s                   + y[10] * x7s                   + y[11] * x8s; // y|xx
      current_data[7]  = y[12] * x3s                  + y[13] * x4s                   + y[14] * x5s; // y|xy
      current_data[8]  = y[15] * scale0               + y[16] * scale1                + y[17] * scale2; // y|yy
      current_data[9]  = y[9] * x3s * z[3]            + y[10] * x4s * z[4]            + y[11] * x5s * z[5]; // y|xz
      current_data[10] = y[12] * z3s                  + y[13] * z4s                   + y[14] * z5s; // y|yz
      current_data[11] = y[9] * z6s                   + y[10] * z7s                   + y[11] * z8s; // y|zz

      current_data[12] = z[9] * x6s                   + z[10] * x7s                   + z[11] * x8s; // z|xx
      current_data[13] = z[9] * x3s * y[3]            + z[10] * x4s * y[4]            + z[11] * x5s * y[5]; // z|xy
      current_data[14] = z[9] * y6s                   + z[10] * y7s                   + z[11] * y8s; // z|yy
      current_data[15] = z[12] * x3s                  + z[13] * x4s                   + z[14] * x5s; // z|xz
      current_data[16] = z[12] * y3s                  + z[13] * y4s                   + z[14] * y5s; // z|yz
      current_data[17] = z[15] * scale0               + z[16] * scale1                + z[17] * scale2; // z|zz


      current_data[18] = x[24] * scale0               + x[25] * scale1                + x[26] * scale2; // xx|xx
      current_data[19] = x[21] * y3s                  + x[22] * y4s                   + x[23] * y5s; // xx|xy
      current_data[20] = x[18] * y6s                  + x[19] * y7s                   + x[20] * y8s; // xx|yy
      current_data[21] = x[21] * z3s                  + x[22] * z4s                   + x[23] * z5s; // xx|xz
      current_data[22] = x[18] * y3s * z[3]           + x[19] * y4s * z[4]            + x[20] * y5s * z[5];// xx|yz
      current_data[23] = x[18] * z6s                  + x[19] * z7s                   + x[20] * z8s; // xx|zz

      current_data[24] = x[15] * y9s                  + x[16] * y10s                  + x[17] * y11s; // xy|xx
      current_data[25] = x[12] * y[12] * scale0       + x[13] * y[13] * scale1        + x[14] * y[14] * scale2; // xy|xy
      current_data[26] = x9s * y[15]                  + x10s * y[16]                  + x11s * y[17]; // xy|yy
      current_data[27] = x[12] * y[9] * z3s           + x[13] * y[10] * z4s           + x[14] * y[11] * z5s; // xy|xz
      current_data[28] = x[9] * y[12] * z3s           + x[10] * y[13] * z4s           + x[11] * y[14] * z5s; // xy|yz
      current_data[29] = x[9] * y[9] * z6s            + x[10] * y[10] * z7s           + x[11] * y[11] * z8s; // xy|zz

      current_data[30] = y[18] * x6s                  + y[19] * x7s                   + y[20] * x8s; // yy|xx
      current_data[31] = y[21] * x3s                  + y[22] * x4s                   + y[23] * x5s; // yy|xy
      current_data[32] = y[24] * scale0               + y[25] * scale1                + y[26] * scale2; // yy|yy
      current_data[33] = y[18] * x[3] * z3s           + y[19] * x[4] * z4s            + y[20] * x[5] * z5s;// yy|xz
      current_data[34] = y[21] * z3s                  + y[22] * z4s                   + y[23] * z5s; // yy|yz
      current_data[35] = y[18] * z6s                  + y[19] * z7s                   + y[20] * z8s; // yy|zz

      current_data[36] = x[15] * z9s                  + x[16] * z10s                  + x[17] * z11s; // xz|xx
      current_data[37] = x[12] * z[9] * y3s           + x[13] * z[10] * y4s           + x[14] * z[11] * y5s; // xz|xy
      current_data[38] = x[9] * z[9] * y6s            + x[10] * z[10] * y7s           + x[11] * z[11] * y8s; // xz|yy
      current_data[39] = x[12] * z[12] * scale0       + x[13] * z[13] * scale1        + x[14] * z[14] * scale2; // xz|xz
      current_data[40] = x[9] * z[12] * y3s           + x[10] * z[13] * y4s           + x[11] * z[14] * y5s; // xz|yz
      current_data[41] = x9s * z[15]                  + x10s * z[16]                  + x11s * z[17]; // xz|zz

      current_data[42] = z[9] * y[9] * x6s            + z[10] * y[10] * x7s           + z[11] * y[11] * x8s; // yz|xx
      current_data[43] = z[9] * y[12] * x3s           + z[10] * y[13] * x4s           + z[11] * y[14] * x5s; // yz|xy
      current_data[44] = z9s * y[15]                  + z10s * y[16]                  + z11s * y[17]; // yz|yy
      current_data[45] = z[12] * y[9] * x3s           + z[13] * y[10] * x4s           + z[14] * y[11] * x5s; // yz|xz
      current_data[46] = z[12] * y[12] * scale0       + z[13] * y[13] * scale1        + z[14] * y[14] * scale2; // yz|yz
      current_data[47] = y9s * z[15]                  + y10s * z[16]                  + y11s * z[17]; // yz|zz

      current_data[48] = z[18] * x6s                  + z[19] * x7s                   + z[20] * x8s; // zz|xx
      current_data[49] = z[18] * y3s * x[3]           + z[19] * y4s * x[4]            + z[20] * y5s * x[5];// zz|xy
      current_data[50] = z[18] * y6s                  + z[19] * y7s                   + z[20] * y8s; // zz|yy
      current_data[51] = z[21] * x3s                  + z[22] * x4s                   + z[23] * x5s; // zz|xz
      current_data[52] = z[21] * y3s                  + z[22] * y4s                   + z[23] * y5s; // zz|yz
      current_data[53] = z[24] * scale0               + z[25] * scale1                + z[26] * scale2; // xx|xx

    }
  } else if (ang == 33856) { // (1, 1, 2, 0)
    double x[27];
    double y[27];
    double z[27];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;

      const double oxp2 = 0.5 / xp_[ii];
      const double opq = 1.0 / (xp_[ii] + xq_[ii]);
      const double xqopq = xq_[ii] * opq;
      const double c00i0x = P_[ii3]     - ax;
      const double c00i0y = P_[ii3 + 1] - ay;
      const double c00i0z = P_[ii3 + 2] - az;
      const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
      const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
      const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
      const double b10i0 = xqopq * oxp2;
      const double b10_0 = oxp2 - b10i0 * roots_[offset];
      const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
      const double b10_2 = oxp2 - b10i0 * roots_[offset + 2];
      const double oxq2 = 0.5 / xq_[ii];
      const double xpopq = xp_[ii] * opq;
      const double d00i0x = Q_[ii3] - cx;
      const double d00i0y = Q_[ii3 + 1] - cy;
      const double d00i0z = Q_[ii3 + 2] - cz;
      const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
      const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
      const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
      const double b00int = 0.5 * opq;
      const double b00_0 = b00int * roots_[offset];
      const double b00_1 = b00int * roots_[offset + 1];
      const double b00_2 = b00int * roots_[offset + 2];
      const double b01i0 = xpopq * oxq2;
      const double b01_0 = oxq2 - b01i0 * roots_[offset];
      const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
      const double b01_2 = oxq2 - b01i0 * roots_[offset + 2];

      x[3]  = c00i0x - c00i1x * roots_[offset];
      x[4]  = c00i0x - c00i1x * roots_[offset + 1];
      x[5]  = c00i0x - c00i1x * roots_[offset + 2];
      x[6]  = x[3] * x[3] + b10_0;
      x[7]  = x[4] * x[4] + b10_1;
      x[8]  = x[5] * x[5] + b10_2;
      x[9]  = d00i0x + d00i1x * roots_[offset];
      x[10] = d00i0x + d00i1x * roots_[offset + 1];
      x[11] = d00i0x + d00i1x * roots_[offset + 2];
      x[12] = x[3] * x[9]  + b00_0;
      x[13] = x[4] * x[10] + b00_1;
      x[14] = x[5] * x[11] + b00_2;
      x[15] = x[3] * (x[12] + b00_0) + b10_0 * x[9];
      x[16] = x[4] * (x[13] + b00_1) + b10_1 * x[10];
      x[17] = x[5] * (x[14] + b00_2) + b10_2 * x[11];
      x[18] = x[9]  * x[9]  + b01_0;
      x[19] = x[10] * x[10] + b01_1;
      x[20] = x[11] * x[11] + b01_2;
      x[21] = x[3] * x[18] + 2.0 * b00_0 * x[9];
      x[22] = x[4] * x[19] + 2.0 * b00_1 * x[10];
      x[23] = x[5] * x[20] + 2.0 * b00_2 * x[11];
      x[24] = x[3] * x[21] + 2.0 * b00_0 * x[12] + b10_0 * x[18];
      x[25] = x[4] * x[22] + 2.0 * b00_1 * x[13] + b10_1 * x[19];
      x[26] = x[5] * x[23] + 2.0 * b00_2 * x[14] + b10_2 * x[20];

      y[3]  = c00i0y - c00i1y * roots_[offset];
      y[4]  = c00i0y - c00i1y * roots_[offset + 1];
      y[5]  = c00i0y - c00i1y * roots_[offset + 2];
      y[6]  = y[3] * y[3] + b10_0;
      y[7]  = y[4] * y[4] + b10_1;
      y[8]  = y[5] * y[5] + b10_2;
      y[9]  = d00i0y + d00i1y * roots_[offset];
      y[10] = d00i0y + d00i1y * roots_[offset + 1];
      y[11] = d00i0y + d00i1y * roots_[offset + 2];
      y[12] = y[3] * y[9]  + b00_0;
      y[13] = y[4] * y[10] + b00_1;
      y[14] = y[5] * y[11] + b00_2;
      y[15] = y[3] * (y[12] + b00_0) + b10_0 * y[9];
      y[16] = y[4] * (y[13] + b00_1) + b10_1 * y[10];
      y[17] = y[5] * (y[14] + b00_2) + b10_2 * y[11];
      y[18] = y[9]  * y[9]  + b01_0;
      y[19] = y[10] * y[10] + b01_1;
      y[20] = y[11] * y[11] + b01_2;
      y[21] = y[3] * y[18] + 2.0 * b00_0 * y[9];
      y[22] = y[4] * y[19] + 2.0 * b00_1 * y[10];
      y[23] = y[5] * y[20] + 2.0 * b00_2 * y[11];
      y[24] = y[3] * y[21] + 2.0 * b00_0 * y[12] + b10_0 * y[18];
      y[25] = y[4] * y[22] + 2.0 * b00_1 * y[13] + b10_1 * y[19];
      y[26] = y[5] * y[23] + 2.0 * b00_2 * y[14] + b10_2 * y[20];

      z[3]  = c00i0z - c00i1z * roots_[offset];
      z[4]  = c00i0z - c00i1z * roots_[offset + 1];
      z[5]  = c00i0z - c00i1z * roots_[offset + 2];
      z[6]  = z[3] * z[3] + b10_0;
      z[7]  = z[4] * z[4] + b10_1;
      z[8]  = z[5] * z[5] + b10_2;
      z[9]  = d00i0z + d00i1z * roots_[offset];
      z[10] = d00i0z + d00i1z * roots_[offset + 1];
      z[11] = d00i0z + d00i1z * roots_[offset + 2];
      z[12] = z[3] * z[9]  + b00_0;
      z[13] = z[4] * z[10] + b00_1;
      z[14] = z[5] * z[11] + b00_2;
      z[15] = z[3] * (z[12] + b00_0) + b10_0 * z[9];
      z[16] = z[4] * (z[13] + b00_1) + b10_1 * z[10];
      z[17] = z[5] * (z[14] + b00_2) + b10_2 * z[11];
      z[18] = z[9]  * z[9]  + b01_0;
      z[19] = z[10] * z[10] + b01_1;
      z[20] = z[11] * z[11] + b01_2;
      z[21] = z[3] * z[18] + 2.0 * b00_0 * z[9];
      z[22] = z[4] * z[19] + 2.0 * b00_1 * z[10];
      z[23] = z[5] * z[20] + 2.0 * b00_2 * z[11];
      z[24] = z[3] * z[21] + 2.0 * b00_0 * z[12] + b10_0 * z[18];
      z[25] = z[4] * z[22] + 2.0 * b00_1 * z[13] + b10_1 * z[19];
      z[26] = z[5] * z[23] + 2.0 * b00_2 * z[14] + b10_2 * z[20];

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double scale2 = coeff_[ii] * weights_[offset + 2];

      const double x3s = x[3] * scale0;
      const double x4s = x[4] * scale1;
      const double x5s = x[5] * scale2;
      const double x6s = x[6] * scale0;
      const double x7s = x[7] * scale1;
      const double x8s = x[8] * scale2;
      const double x9s = x[9] * scale0;
      const double x10s = x[10] * scale1;
      const double x11s = x[11] * scale2;
      const double y3s = y[3] * scale0;
      const double y4s = y[4] * scale1;
      const double y5s = y[5] * scale2;
      const double y6s = y[6] * scale0;
      const double y7s = y[7] * scale1;
      const double y8s = y[8] * scale2;
      const double y9s = y[9] * scale0;
      const double y10s = y[10] * scale1;
      const double y11s = y[11] * scale2;
      const double z3s = z[3] * scale0;
      const double z4s = z[4] * scale1;
      const double z5s = z[5] * scale2;
      const double z6s = z[6] * scale0;
      const double z7s = z[7] * scale1;
      const double z8s = z[8] * scale2;
      const double z9s = z[9] * scale0;
      const double z10s = z[10] * scale1;
      const double z11s = z[11] * scale2;

      current_data[ 0] = x[21] * scale0               + x[22] * scale1                + x[23] * scale2; // xx|x
      current_data[ 1] = x[18] * y3s                  + x[19] * y4s                   + x[20] * y5s; // xx|y
      current_data[ 2] = x[18] * z3s                  + x[19] * z4s                   + x[20] * z5s; // xx|z
      current_data[ 3] = x[24] * scale0               + x[25] * scale1                + x[26] * scale2; // xx|xx
      current_data[ 4] = x[21] * y3s                  + x[22] * y4s                   + x[23] * y5s; // xx|xy
      current_data[ 5] = x[18] * y6s                  + x[19] * y7s                   + x[20] * y8s; // xx|yy
      current_data[ 6] = x[21] * z3s                  + x[22] * z4s                   + x[23] * z5s; // xx|xz
      current_data[ 7] = x[18] * y3s * z[3]           + x[19] * y4s * z[4]            + x[20] * y5s * z[5];// xx|yz
      current_data[ 8] = x[18] * z6s                  + x[19] * z7s                   + x[20] * z8s; // xx|zz

      const double x9y9 = x[9] * y[9];
      const double x10y10 = x[10] * y[10];
      const double x11y11 = x[11] * y[11];
      current_data[ 9] = x[12] * y9s                  + x[13] * y10s                  + x[14] * y11s; // xy|x
      current_data[10] = y[12] * x9s                  + y[13] * x10s                  + y[14] * x11s; // xy|y
      current_data[11] = x9y9 * z3s                   + x10y10 * z4s                  + x11y11 * z5s; // xy|z
      current_data[12] = x[15] * y9s                  + x[16] * y10s                  + x[17] * y11s; // xy|xx
      current_data[13] = x[12] * y[12] * scale0       + x[13] * y[13] * scale1        + x[14] * y[14] * scale2; // xy|xy
      current_data[14] = x9s * y[15]                  + x10s * y[16]                  + x11s * y[17]; // xy|yy
      current_data[15] = x[12] * y[9] * z3s           + x[13] * y[10] * z4s           + x[14] * y[11] * z5s; // xy|xz
      current_data[16] = x[9] * y[12] * z3s           + x[10] * y[13] * z4s           + x[11] * y[14] * z5s; // xy|yz
      current_data[17] = x9y9 * z6s                   + x10y10 * z7s                  + x11y11 * z8s; // xy|zz

      current_data[18] = y[18] * x3s                  + y[19] * x4s                   + y[20] * x5s; // yy|x
      current_data[19] = y[21] * scale0               + y[22] * scale1                + y[23] * scale2; // yy|y
      current_data[20] = y[18] * z3s                  + y[19] * z4s                   + y[20] * z5s; // yy|z
      current_data[21] = y[18] * x6s                  + y[19] * x7s                   + y[20] * x8s; // yy|xx
      current_data[22] = y[21] * x3s                  + y[22] * x4s                   + y[23] * x5s; // yy|xy
      current_data[23] = y[24] * scale0               + y[25] * scale1                + y[26] * scale2; // yy|yy
      current_data[24] = y[18] * x[3] * z3s           + y[19] * x[4] * z4s            + y[20] * x[5] * z5s;// yy|xz
      current_data[25] = y[21] * z3s                  + y[22] * z4s                   + y[23] * z5s; // yy|yz
      current_data[26] = y[18] * z6s                  + y[19] * z7s                   + y[20] * z8s; // yy|zz

      const double x9z9 = x[9] * z[9];
      const double x10z10 = x[10] * z[10];
      const double x11z11 = x[11] * z[11];
      current_data[27] = x[12] * z9s                  + x[13] * z10s                  + x[14] * z11s; // xz|x
      current_data[28] = x9z9 * y3s                   + x10z10 * y4s                  + x11z11 * y5s; // xz|y
      current_data[29] = z[12] * x9s                  + z[13] * x10s                  + z[14] * x11s; // xz|z
      current_data[30] = x[15] * z9s                  + x[16] * z10s                  + x[17] * z11s; // xz|xx
      current_data[31] = x[12] * z[9] * y3s           + x[13] * z[10] * y4s           + x[14] * z[11] * y5s; // xz|xy
      current_data[32] = x9z9 * y6s                   + x10z10 * y7s                  + x11z11 * y8s; // xz|yy
      current_data[33] = x[12] * z[12] * scale0       + x[13] * z[13] * scale1        + x[14] * z[14] * scale2; // xz|xz
      current_data[34] = x[9] * z[12] * y3s           + x[10] * z[13] * y4s           + x[11] * z[14] * y5s; // xz|yz
      current_data[35] = x9s * z[15]                  + x10s * z[16]                  + x11s * z[17]; // xz|zz

      const double y9z9 = y[9] * z[9];
      const double y10z10 = y[10] * z[10];
      const double y11z11 = y[11] * z[11];
      current_data[36] = y9z9 * x3s                   + y10z10 * x4s                  + y11z11 * x5s; // yz|x
      current_data[37] = y[12] * z9s                  + y[13] * z10s                  + y[14] * z11s; // yz|y
      current_data[38] = z[12] * y9s                  + z[13] * y10s                  + z[14] * y11s; // yz|z
      current_data[39] = y9z9 * x6s                   + y10z10 * x7s                  + y11z11 * x8s; // yz|xx
      current_data[40] = z[9] * y[12] * x3s           + z[10] * y[13] * x4s           + z[11] * y[14] * x5s; // yz|xy
      current_data[41] = z9s * y[15]                  + z10s * y[16]                  + z11s * y[17]; // yz|yy
      current_data[42] = z[12] * y[9] * x3s           + z[13] * y[10] * x4s           + z[14] * y[11] * x5s; // yz|xz
      current_data[43] = z[12] * y[12] * scale0       + z[13] * y[13] * scale1        + z[14] * y[14] * scale2; // yz|yz
      current_data[44] = y9s * z[15]                  + y10s * z[16]                  + y11s * z[17]; // yz|zz

      current_data[45] = z[18] * x3s                  + z[19] * x4s                   + z[20] * x5s; // zz|x
      current_data[46] = z[18] * y3s                  + z[19] * y4s                   + z[20] * y5s; // zz|y
      current_data[47] = z[21] * scale0               + z[22] * scale1                + z[23] * scale2; // zz|z
      current_data[48] = z[18] * x6s                  + z[19] * x7s                   + z[20] * x8s; // zz|xx
      current_data[49] = z[18] * y3s * x[3]           + z[19] * y4s * x[4]            + z[20] * y5s * x[5];// zz|xy
      current_data[50] = z[18] * y6s                  + z[19] * y7s                   + z[20] * y8s; // zz|yy
      current_data[51] = z[21] * x3s                  + z[22] * x4s                   + z[23] * x5s; // zz|xz
      current_data[52] = z[21] * y3s                  + z[22] * y4s                   + z[23] * y5s; // zz|yz
      current_data[53] = z[24] * scale0               + z[25] * scale1                + z[26] * scale2; // xx|xx

    }
  } else if (ang == 33825) { // (1, 1, 1, 1)
    double x[27];
    double y[27];
    double z[27];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;

      const double oxp2 = 0.5 / xp_[ii];
      const double opq = 1.0 / (xp_[ii] + xq_[ii]);
      const double xqopq = xq_[ii] * opq;
      const double c00i0x = P_[ii3]     - ax;
      const double c00i0y = P_[ii3 + 1] - ay;
      const double c00i0z = P_[ii3 + 2] - az;
      const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
      const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
      const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
      const double b10i0 = xqopq * oxp2;
      const double b10_0 = oxp2 - b10i0 * roots_[offset];
      const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
      const double b10_2 = oxp2 - b10i0 * roots_[offset + 2];
      const double oxq2 = 0.5 / xq_[ii];
      const double xpopq = xp_[ii] * opq;
      const double d00i0x = Q_[ii3] - cx;
      const double d00i0y = Q_[ii3 + 1] - cy;
      const double d00i0z = Q_[ii3 + 2] - cz;
      const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
      const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
      const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
      const double b00int = 0.5 * opq;
      const double b00_0 = b00int * roots_[offset];
      const double b00_1 = b00int * roots_[offset + 1];
      const double b00_2 = b00int * roots_[offset + 2];
      const double b01i0 = xpopq * oxq2;
      const double b01_0 = oxq2 - b01i0 * roots_[offset];
      const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
      const double b01_2 = oxq2 - b01i0 * roots_[offset + 2];

      x[3]  = c00i0x - c00i1x * roots_[offset];
      x[4]  = c00i0x - c00i1x * roots_[offset + 1];
      x[5]  = c00i0x - c00i1x * roots_[offset + 2];
      x[6]  = x[3] * x[3] + b10_0;
      x[7]  = x[4] * x[4] + b10_1;
      x[8]  = x[5] * x[5] + b10_2;
      x[9]  = d00i0x + d00i1x * roots_[offset];
      x[10] = d00i0x + d00i1x * roots_[offset + 1];
      x[11] = d00i0x + d00i1x * roots_[offset + 2];
      x[12] = x[3] * x[9]  + b00_0;
      x[13] = x[4] * x[10] + b00_1;
      x[14] = x[5] * x[11] + b00_2;
      x[15] = x[3] * (x[12] + b00_0) + b10_0 * x[9];
      x[16] = x[4] * (x[13] + b00_1) + b10_1 * x[10];
      x[17] = x[5] * (x[14] + b00_2) + b10_2 * x[11];
      x[18] = x[9]  * x[9]  + b01_0;
      x[19] = x[10] * x[10] + b01_1;
      x[20] = x[11] * x[11] + b01_2;
      x[21] = x[3] * x[18] + 2.0 * b00_0 * x[9];
      x[22] = x[4] * x[19] + 2.0 * b00_1 * x[10];
      x[23] = x[5] * x[20] + 2.0 * b00_2 * x[11];
      x[24] = x[3] * x[21] + 2.0 * b00_0 * x[12] + b10_0 * x[18];
      x[25] = x[4] * x[22] + 2.0 * b00_1 * x[13] + b10_1 * x[19];
      x[26] = x[5] * x[23] + 2.0 * b00_2 * x[14] + b10_2 * x[20];

      y[3]  = c00i0y - c00i1y * roots_[offset];
      y[4]  = c00i0y - c00i1y * roots_[offset + 1];
      y[5]  = c00i0y - c00i1y * roots_[offset + 2];
      y[6]  = y[3] * y[3] + b10_0;
      y[7]  = y[4] * y[4] + b10_1;
      y[8]  = y[5] * y[5] + b10_2;
      y[9]  = d00i0y + d00i1y * roots_[offset];
      y[10] = d00i0y + d00i1y * roots_[offset + 1];
      y[11] = d00i0y + d00i1y * roots_[offset + 2];
      y[12] = y[3] * y[9]  + b00_0;
      y[13] = y[4] * y[10] + b00_1;
      y[14] = y[5] * y[11] + b00_2;
      y[15] = y[3] * (y[12] + b00_0) + b10_0 * y[9];
      y[16] = y[4] * (y[13] + b00_1) + b10_1 * y[10];
      y[17] = y[5] * (y[14] + b00_2) + b10_2 * y[11];
      y[18] = y[9]  * y[9]  + b01_0;
      y[19] = y[10] * y[10] + b01_1;
      y[20] = y[11] * y[11] + b01_2;
      y[21] = y[3] * y[18] + 2.0 * b00_0 * y[9];
      y[22] = y[4] * y[19] + 2.0 * b00_1 * y[10];
      y[23] = y[5] * y[20] + 2.0 * b00_2 * y[11];
      y[24] = y[3] * y[21] + 2.0 * b00_0 * y[12] + b10_0 * y[18];
      y[25] = y[4] * y[22] + 2.0 * b00_1 * y[13] + b10_1 * y[19];
      y[26] = y[5] * y[23] + 2.0 * b00_2 * y[14] + b10_2 * y[20];

      z[3]  = c00i0z - c00i1z * roots_[offset];
      z[4]  = c00i0z - c00i1z * roots_[offset + 1];
      z[5]  = c00i0z - c00i1z * roots_[offset + 2];
      z[6]  = z[3] * z[3] + b10_0;
      z[7]  = z[4] * z[4] + b10_1;
      z[8]  = z[5] * z[5] + b10_2;
      z[9]  = d00i0z + d00i1z * roots_[offset];
      z[10] = d00i0z + d00i1z * roots_[offset + 1];
      z[11] = d00i0z + d00i1z * roots_[offset + 2];
      z[12] = z[3] * z[9]  + b00_0;
      z[13] = z[4] * z[10] + b00_1;
      z[14] = z[5] * z[11] + b00_2;
      z[15] = z[3] * (z[12] + b00_0) + b10_0 * z[9];
      z[16] = z[4] * (z[13] + b00_1) + b10_1 * z[10];
      z[17] = z[5] * (z[14] + b00_2) + b10_2 * z[11];
      z[18] = z[9]  * z[9]  + b01_0;
      z[19] = z[10] * z[10] + b01_1;
      z[20] = z[11] * z[11] + b01_2;
      z[21] = z[3] * z[18] + 2.0 * b00_0 * z[9];
      z[22] = z[4] * z[19] + 2.0 * b00_1 * z[10];
      z[23] = z[5] * z[20] + 2.0 * b00_2 * z[11];
      z[24] = z[3] * z[21] + 2.0 * b00_0 * z[12] + b10_0 * z[18];
      z[25] = z[4] * z[22] + 2.0 * b00_1 * z[13] + b10_1 * z[19];
      z[26] = z[5] * z[23] + 2.0 * b00_2 * z[14] + b10_2 * z[20];

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double scale2 = coeff_[ii] * weights_[offset + 2];

      const double x3s = x[3] * scale0;
      const double x4s = x[4] * scale1;
      const double x5s = x[5] * scale2;
      const double x6s = x[6] * scale0;
      const double x7s = x[7] * scale1;
      const double x8s = x[8] * scale2;
      const double x9s = x[9] * scale0;
      const double x10s = x[10] * scale1;
      const double x11s = x[11] * scale2;
      const double y3s = y[3] * scale0;
      const double y4s = y[4] * scale1;
      const double y5s = y[5] * scale2;
      const double y6s = y[6] * scale0;
      const double y7s = y[7] * scale1;
      const double y8s = y[8] * scale2;
      const double y9s = y[9] * scale0;
      const double y10s = y[10] * scale1;
      const double y11s = y[11] * scale2;
      const double y12s = y[12] * scale0;
      const double y13s = y[13] * scale1;
      const double y14s = y[14] * scale2;
      const double z3s = z[3] * scale0;
      const double z4s = z[4] * scale1;
      const double z5s = z[5] * scale2;
      const double z6s = z[6] * scale0;
      const double z7s = z[7] * scale1;
      const double z8s = z[8] * scale2;
      const double z9s = z[9] * scale0;
      const double z10s = z[10] * scale1;
      const double z11s = z[11] * scale2;
      const double z12s = z[12] * scale0;
      const double z13s = z[13] * scale1;
      const double z14s = z[14] * scale2;

      current_data[0]  = x[12] * scale0               + x[13] * scale1                + x[14] * scale2; // x|x
      current_data[1]  = x[9] * y3s                   + x[10] * y4s                   + x[11] * y5s; // x|y
      current_data[2]  = x[9] * z3s                   + x[10] * z4s                   + x[11] * z5s; // x|z
      current_data[3]  = x[15] * scale0               + x[16] * scale1                + x[17] * scale2; // x|xx
      current_data[4]  = x[12] * y3s                  + x[13] * y4s                   + x[14] * y5s; // x|xy
      current_data[5]  = x[9] * y6s                   + x[10] * y7s                   + x[11] * y8s; // x|yy
      current_data[6]  = x[12] * z3s                  + x[13] * z4s                   + x[14] * z5s; // x|xz
      current_data[7]  = x[9] * y3s * z[3]            + x[10] * y4s * z[4]            + x[11] * y5s * z[5]; // x|yz
      current_data[8]  = x[9] * z6s                   + x[10] * z7s                   + x[11] * z8s; // x|zz

      current_data[9]  = y[9] * x3s                   + y[10] * x4s                   + y[11] * x5s; // y|x
      current_data[10] = y12s                         + y13s                          + y14s; // y|y
      current_data[11] = y[9] * z3s                   + y[10] * z4s                   + y[11] * z5s; // y|z
      current_data[12] = y[9] * x6s                   + y[10] * x7s                   + y[11] * x8s; // y|xx
      current_data[13] = y[12] * x3s                  + y[13] * x4s                   + y[14] * x5s; // y|xy
      current_data[14] = y[15] * scale0               + y[16] * scale1                + y[17] * scale2; // y|yy
      current_data[15] = y[9] * x3s * z[3]            + y[10] * x4s * z[4]            + y[11] * x5s * z[5]; // y|xz
      current_data[16] = y[12] * z3s                  + y[13] * z4s                   + y[14] * z5s; // y|yz
      current_data[17] = y[9] * z6s                   + y[10] * z7s                   + y[11] * z8s; // y|zz

      current_data[18] = z[9] * x3s                   + z[10] * x4s                   + z[11] * x5s; // z|x
      current_data[19] = z[9] * y3s                   + z[10] * y4s                   + z[11] * y5s; // z|y
      current_data[20] = z12s                         + z13s                          + z14s; // z|z
      current_data[21] = z[9] * x6s                   + z[10] * x7s                   + z[11] * x8s; // z|xx
      current_data[22] = z[9] * x3s * y[3]            + z[10] * x4s * y[4]            + z[11] * x5s * y[5]; // z|xy
      current_data[23] = z[9] * y6s                   + z[10] * y7s                   + z[11] * y8s; // z|yy
      current_data[24] = z[12] * x3s                  + z[13] * x4s                   + z[14] * x5s; // z|xz
      current_data[25] = z[12] * y3s                  + z[13] * y4s                   + z[14] * y5s; // z|yz
      current_data[26] = z[15] * scale0               + z[16] * scale1                + z[17] * scale2; // z|zz

      current_data[27] = x[21] * scale0               + x[22] * scale1                + x[23] * scale2; // xx|x
      current_data[28] = x[18] * y3s                  + x[19] * y4s                   + x[20] * y5s; // xx|y
      current_data[29] = x[18] * z3s                  + x[19] * z4s                   + x[20] * z5s; // xx|z
      current_data[30] = x[24] * scale0               + x[25] * scale1                + x[26] * scale2; // xx|xx
      current_data[31] = x[21] * y3s                  + x[22] * y4s                   + x[23] * y5s; // xx|xy
      current_data[32] = x[18] * y6s                  + x[19] * y7s                   + x[20] * y8s; // xx|yy
      current_data[33] = x[21] * z3s                  + x[22] * z4s                   + x[23] * z5s; // xx|xz
      current_data[34] = x[18] * y3s * z[3]           + x[19] * y4s * z[4]            + x[20] * y5s * z[5];// xx|yz
      current_data[35] = x[18] * z6s                  + x[19] * z7s                   + x[20] * z8s; // xx|zz

      const double x9y9 = x[9] * y[9];
      const double x10y10 = x[10] * y[10];
      const double x11y11 = x[11] * y[11];
      current_data[36] = x[12] * y9s                  + x[13] * y10s                  + x[14] * y11s; // xy|x
      current_data[37] = y[12] * x9s                  + y[13] * x10s                  + y[14] * x11s; // xy|y
      current_data[38] = x9y9 * z3s                   + x10y10 * z4s                  + x11y11 * z5s; // xy|z
      current_data[39] = x[15] * y9s                  + x[16] * y10s                  + x[17] * y11s; // xy|xx
      current_data[40] = x[12] * y12s                 + x[13] * y13s                  + x[14] * y14s; // xy|xy
      current_data[41] = x9s * y[15]                  + x10s * y[16]                  + x11s * y[17]; // xy|yy
      current_data[42] = x[12] * y[9] * z3s           + x[13] * y[10] * z4s           + x[14] * y[11] * z5s; // xy|xz
      current_data[43] = x[9] * y[12] * z3s           + x[10] * y[13] * z4s           + x[11] * y[14] * z5s; // xy|yz
      current_data[44] = x9y9 * z6s                   + x10y10 * z7s                  + x11y11 * z8s; // xy|zz

      current_data[45] = y[18] * x3s                  + y[19] * x4s                   + y[20] * x5s; // yy|x
      current_data[46] = y[21] * scale0               + y[22] * scale1                + y[23] * scale2; // yy|y
      current_data[47] = y[18] * z3s                  + y[19] * z4s                   + y[20] * z5s; // yy|z
      current_data[48] = y[18] * x6s                  + y[19] * x7s                   + y[20] * x8s; // yy|xx
      current_data[49] = y[21] * x3s                  + y[22] * x4s                   + y[23] * x5s; // yy|xy
      current_data[50] = y[24] * scale0               + y[25] * scale1                + y[26] * scale2; // yy|yy
      current_data[51] = y[18] * x[3] * z3s           + y[19] * x[4] * z4s            + y[20] * x[5] * z5s;// yy|xz
      current_data[52] = y[21] * z3s                  + y[22] * z4s                   + y[23] * z5s; // yy|yz
      current_data[53] = y[18] * z6s                  + y[19] * z7s                   + y[20] * z8s; // yy|zz

      const double x9z9 = x[9] * z[9];
      const double x10z10 = x[10] * z[10];
      const double x11z11 = x[11] * z[11];
      current_data[54] = x[12] * z9s                  + x[13] * z10s                  + x[14] * z11s; // xz|x
      current_data[55] = x9z9 * y3s                   + x10z10 * y4s                  + x11z11 * y5s; // xz|y
      current_data[56] = z[12] * x9s                  + z[13] * x10s                  + z[14] * x11s; // xz|z
      current_data[57] = x[15] * z9s                  + x[16] * z10s                  + x[17] * z11s; // xz|xx
      current_data[58] = x[12] * z[9] * y3s           + x[13] * z[10] * y4s           + x[14] * z[11] * y5s; // xz|xy
      current_data[59] = x9z9 * y6s                   + x10z10 * y7s                  + x11z11 * y8s; // xz|yy
      current_data[60] = x[12] * z12s                 + x[13] * z13s                  + x[14] * z14s; // xz|xz
      current_data[61] = x[9] * z[12] * y3s           + x[10] * z[13] * y4s           + x[11] * z[14] * y5s; // xz|yz
      current_data[62] = x9s * z[15]                  + x10s * z[16]                  + x11s * z[17]; // xz|zz

      const double y9z9 = y[9] * z[9];
      const double y10z10 = y[10] * z[10];
      const double y11z11 = y[11] * z[11];
      current_data[63] = y9z9 * x3s                   + y10z10 * x4s                  + y11z11 * x5s; // yz|x
      current_data[64] = y[12] * z9s                  + y[13] * z10s                  + y[14] * z11s; // yz|y
      current_data[65] = z[12] * y9s                  + z[13] * y10s                  + z[14] * y11s; // yz|z
      current_data[66] = y9z9 * x6s                   + y10z10 * x7s                  + y11z11 * x8s; // yz|xx
      current_data[67] = z[9] * y[12] * x3s           + z[10] * y[13] * x4s           + z[11] * y[14] * x5s; // yz|xy
      current_data[68] = z9s * y[15]                  + z10s * y[16]                  + z11s * y[17]; // yz|yy
      current_data[69] = z[12] * y[9] * x3s           + z[13] * y[10] * x4s           + z[14] * y[11] * x5s; // yz|xz
      current_data[70] = z12s * y[12]                 + z13s * y[13]                  + z14s * y[14]; // yz|yz
      current_data[71] = y9s * z[15]                  + y10s * z[16]                  + y11s * z[17]; // yz|zz

      current_data[72] = z[18] * x3s                  + z[19] * x4s                   + z[20] * x5s; // zz|x
      current_data[73] = z[18] * y3s                  + z[19] * y4s                   + z[20] * y5s; // zz|y
      current_data[74] = z[21] * scale0               + z[22] * scale1                + z[23] * scale2; // zz|z
      current_data[75] = z[18] * x6s                  + z[19] * x7s                   + z[20] * x8s; // zz|xx
      current_data[76] = z[18] * y3s * x[3]           + z[19] * y4s * x[4]            + z[20] * y5s * x[5];// zz|xy
      current_data[77] = z[18] * y6s                  + z[19] * y7s                   + z[20] * y8s; // zz|yy
      current_data[78] = z[21] * x3s                  + z[22] * x4s                   + z[23] * x5s; // zz|xz
      current_data[79] = z[21] * y3s                  + z[22] * y4s                   + z[23] * y5s; // zz|yz
      current_data[80] = z[24] * scale0               + z[25] * scale1                + z[26] * scale2; // xx|xx

    }
  } else if (ang == 67584) { // (2, 2, 0, 0)
    double x[15];
    double y[15];
    double z[15];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;

      const double oxp2 = 0.5 / xp_[ii];
      const double opq = 1.0 / (xp_[ii] + xq_[ii]);
      const double xqopq = xq_[ii] * opq;
      const double c00i0x = P_[ii3]     - ax;
      const double c00i0y = P_[ii3 + 1] - ay;
      const double c00i0z = P_[ii3 + 2] - az;
      const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
      const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
      const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
      const double b10i0 = xqopq * oxp2;
      const double b10_0 = oxp2 - b10i0 * roots_[offset];
      const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
      const double b10_2 = oxp2 - b10i0 * roots_[offset + 2];

      x[3]  = c00i0x - c00i1x * roots_[offset];
      x[4]  = c00i0x - c00i1x * roots_[offset + 1];
      x[5]  = c00i0x - c00i1x * roots_[offset + 2];
      x[6]  = x[3] * x[3] + b10_0;
      x[7]  = x[4] * x[4] + b10_1;
      x[8]  = x[5] * x[5] + b10_2;
      x[9]  = x[3] * (x[6] + 2.0 * b10_0);
      x[10] = x[4] * (x[7] + 2.0 * b10_1);
      x[11] = x[5] * (x[8] + 2.0 * b10_2);
      x[12] = x[3] * x[9] + 3.0 * b10_0 * x[6];
      x[13] = x[4] * x[10] + 3.0 * b10_1 * x[7];
      x[14] = x[5] * x[11] + 3.0 * b10_2 * x[8];

      y[3]  = c00i0y - c00i1y * roots_[offset];
      y[4]  = c00i0y - c00i1y * roots_[offset + 1];
      y[5]  = c00i0y - c00i1y * roots_[offset + 2];
      y[6]  = y[3] * y[3] + b10_0;
      y[7]  = y[4] * y[4] + b10_1;
      y[8]  = y[5] * y[5] + b10_2;
      y[9]  = y[3] * (y[6] + 2.0 * b10_0);
      y[10] = y[4] * (y[7] + 2.0 * b10_1);
      y[11] = y[5] * (y[8] + 2.0 * b10_2);
      y[12] = y[3] * y[9]  + 3.0 * b10_0 * y[6];
      y[13] = y[4] * y[10] + 3.0 * b10_1 * y[7];
      y[14] = y[5] * y[11] + 3.0 * b10_2 * y[8];

      z[3]  = c00i0z - c00i1z * roots_[offset];
      z[4]  = c00i0z - c00i1z * roots_[offset + 1];
      z[5]  = c00i0z - c00i1z * roots_[offset + 2];
      z[6]  = z[3] * z[3] + b10_0;
      z[7]  = z[4] * z[4] + b10_1;
      z[8]  = z[5] * z[5] + b10_2;
      z[9]  = z[3] * (z[6] + 2.0 * b10_0);
      z[10] = z[4] * (z[7] + 2.0 * b10_1);
      z[11] = z[5] * (z[8] + 2.0 * b10_2);
      z[12] = z[3] * z[9]  + 3.0 * b10_0 * z[6];
      z[13] = z[4] * z[10] + 3.0 * b10_1 * z[7];
      z[14] = z[5] * z[11] + 3.0 * b10_2 * z[8];

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double scale2 = coeff_[ii] * weights_[offset + 2];

      const double x3s = x[3] * scale0;
      const double x4s = x[4] * scale1;
      const double x5s = x[5] * scale2;
      const double y3s = y[3] * scale0;
      const double y4s = y[4] * scale1;
      const double y5s = y[5] * scale2;
      const double y6s = y[6] * scale0;
      const double y7s = y[7] * scale1;
      const double y8s = y[8] * scale2;
      const double z3s = z[3] * scale0;
      const double z4s = z[4] * scale1;
      const double z5s = z[5] * scale2;
      const double z6s = z[6] * scale0;
      const double z7s = z[7] * scale1;
      const double z8s = z[8] * scale2;

      const double x3y3 = x[3] * y[3];
      const double x4y4 = x[4] * y[4];
      const double x5y5 = x[5] * y[5];
      const double x3y6 = x[3] * y[6];
      const double x4y7 = x[4] * y[7];
      const double x5y8 = x[5] * y[8];
      const double x6y3 = y[3] * x[6];
      const double x7y4 = y[4] * x[7];
      const double x8y5 = y[5] * x[8];

      current_data[0]  = x[6] * scale0        + x[7] * scale1        + x[8] * scale2; // xx
      current_data[1]  = x[3] * y3s           + x[4] * y4s           + x[5] * y5s; // xy
      current_data[2]  = y6s                  + y7s                  + y8s; // yy
      current_data[3]  = x[3] * z3s           + x[4] * z4s           + x[5] * z5s; // xz
      current_data[4]  = y[3] * z3s           + y[4] * z4s           + y[5] * z5s; // yz
      current_data[5]  = z6s                  + z7s                  + z8s; // zz

      current_data[6] = x[9] * scale0               + x[10] * scale1              + x[11] * scale2; // xxx
      current_data[7] = x6y3 * scale0               + x7y4 * scale1               + x8y5 * scale2; // xxy
      current_data[8] = x3y6 * scale0               + x4y7 * scale1               + x5y8 * scale2; // xyy
      current_data[9] = y[9] * scale0               + y[10] * scale1              + y[11] * scale2; // yyy
      current_data[10] = x[6] * z3s                  + x[7] * z4s                  + x[8] * z5s; // xxz
      current_data[11] = x3y3 * z3s                  + x4y4 * z4s                  + x5y5 * z5s; // xyz
      current_data[12] = y[6] * z3s                  + y[7] * z4s                  + y[8] * z5s; // yyz
      current_data[13] = x[3] * z6s                  + x[4] * z7s                  + x[5] * z8s; // xzz
      current_data[14] = y[3] * z6s                  + y[4] * z7s                  + y[5] * z8s; // yzz
      current_data[15] = z[9] * scale0               + z[10] * scale1              + z[11] * scale2; // zzz

      current_data[16] = x[12] * scale0              + x[13] * scale1              + x[14] * scale2; // xxxx
      current_data[17] = x[9] * y3s                  + x[10] * y4s                 + x[11] * y5s; // xxxy
      current_data[18] = x[6] * y6s                  + x[7] * y7s                  + x[8] * y8s; // xxyy
      current_data[19] = y[9] * x3s                  + y[10] * x4s                 + y[11] * x5s; // xyyy
      current_data[20] = y[12] * scale0              + y[13] * scale1              + y[14] * scale2; // yyyy
      current_data[21] = x[9] * z3s                  + x[10] * z4s                 + x[11] * z5s; // xxxz
      current_data[22] = x6y3 * z3s                  + x7y4 * z4s                  + x8y5 * z5s; // xxyz
      current_data[23] = x3y6 * z3s                  + x4y7 * z4s                  + x5y8 * z5s; // xyyz
      current_data[24] = y[9] * z3s                  + y[10] * z4s                 + y[11] * z5s; // yyyz
      current_data[25] = x[6] * z6s                  + x[7] * z7s                  + x[8] * z8s; // xxzz
      current_data[26] = x3y3 * z6s                  + x4y4 * z7s                  + x5y5 * z8s; // xyzz
      current_data[27] = y[6] * z6s                  + y[7] * z7s                  + y[8] * z8s; // yyzz
      current_data[28] = z[9] * x3s                  + z[10] * x4s                 + z[11] * x5s; // xzzz
      current_data[29] = z[9] * y3s                  + z[10] * y4s                 + z[11] * y5s; // yzzz
      current_data[30] = z[12] * scale0              + z[13] * scale1              + z[14] * scale2; // yyyy

    }
  } else if (ang == 99328) { // (3, 1, 0, 0)
    double x[15];
    double y[15];
    double z[15];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;

      const double oxp2 = 0.5 / xp_[ii];
      const double opq = 1.0 / (xp_[ii] + xq_[ii]);
      const double xqopq = xq_[ii] * opq;
      const double c00i0x = P_[ii3]     - ax;
      const double c00i0y = P_[ii3 + 1] - ay;
      const double c00i0z = P_[ii3 + 2] - az;
      const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
      const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
      const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
      const double b10i0 = xqopq * oxp2;
      const double b10_0 = oxp2 - b10i0 * roots_[offset];
      const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
      const double b10_2 = oxp2 - b10i0 * roots_[offset + 2];

      x[3]  = c00i0x - c00i1x * roots_[offset];
      x[4]  = c00i0x - c00i1x * roots_[offset + 1];
      x[5]  = c00i0x - c00i1x * roots_[offset + 2];
      x[6]  = x[3] * x[3] + b10_0;
      x[7]  = x[4] * x[4] + b10_1;
      x[8]  = x[5] * x[5] + b10_2;
      x[9]  = x[3] * (x[6] + 2.0 * b10_0);
      x[10] = x[4] * (x[7] + 2.0 * b10_1);
      x[11] = x[5] * (x[8] + 2.0 * b10_2);
      x[12] = x[3] * x[9] + 3.0 * b10_0 * x[6];
      x[13] = x[4] * x[10] + 3.0 * b10_1 * x[7];
      x[14] = x[5] * x[11] + 3.0 * b10_2 * x[8];

      y[3]  = c00i0y - c00i1y * roots_[offset];
      y[4]  = c00i0y - c00i1y * roots_[offset + 1];
      y[5]  = c00i0y - c00i1y * roots_[offset + 2];
      y[6]  = y[3] * y[3] + b10_0;
      y[7]  = y[4] * y[4] + b10_1;
      y[8]  = y[5] * y[5] + b10_2;
      y[9]  = y[3] * (y[6] + 2.0 * b10_0);
      y[10] = y[4] * (y[7] + 2.0 * b10_1);
      y[11] = y[5] * (y[8] + 2.0 * b10_2);
      y[12] = y[3] * y[9]  + 3.0 * b10_0 * y[6];
      y[13] = y[4] * y[10] + 3.0 * b10_1 * y[7];
      y[14] = y[5] * y[11] + 3.0 * b10_2 * y[8];

      z[3]  = c00i0z - c00i1z * roots_[offset];
      z[4]  = c00i0z - c00i1z * roots_[offset + 1];
      z[5]  = c00i0z - c00i1z * roots_[offset + 2];
      z[6]  = z[3] * z[3] + b10_0;
      z[7]  = z[4] * z[4] + b10_1;
      z[8]  = z[5] * z[5] + b10_2;
      z[9]  = z[3] * (z[6] + 2.0 * b10_0);
      z[10] = z[4] * (z[7] + 2.0 * b10_1);
      z[11] = z[5] * (z[8] + 2.0 * b10_2);
      z[12] = z[3] * z[9]  + 3.0 * b10_0 * z[6];
      z[13] = z[4] * z[10] + 3.0 * b10_1 * z[7];
      z[14] = z[5] * z[11] + 3.0 * b10_2 * z[8];

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double scale2 = coeff_[ii] * weights_[offset + 2];

      const double x3s = x[3] * scale0;
      const double x4s = x[4] * scale1;
      const double x5s = x[5] * scale2;
      const double y3s = y[3] * scale0;
      const double y4s = y[4] * scale1;
      const double y5s = y[5] * scale2;
      const double y6s = y[6] * scale0;
      const double y7s = y[7] * scale1;
      const double y8s = y[8] * scale2;
      const double z3s = z[3] * scale0;
      const double z4s = z[4] * scale1;
      const double z5s = z[5] * scale2;
      const double z6s = z[6] * scale0;
      const double z7s = z[7] * scale1;
      const double z8s = z[8] * scale2;

      const double x3y3 = x[3] * y[3];
      const double x4y4 = x[4] * y[4];
      const double x5y5 = x[5] * y[5];
      const double x3y6 = x[3] * y[6];
      const double x4y7 = x[4] * y[7];
      const double x5y8 = x[5] * y[8];
      const double x6y3 = y[3] * x[6];
      const double x7y4 = y[4] * x[7];
      const double x8y5 = y[5] * x[8];

      current_data[0] = x[9] * scale0               + x[10] * scale1              + x[11] * scale2; // xxx
      current_data[1] = x6y3 * scale0               + x7y4 * scale1               + x8y5 * scale2; // xxy
      current_data[2] = x3y6 * scale0               + x4y7 * scale1               + x5y8 * scale2; // xyy
      current_data[3] = y[9] * scale0               + y[10] * scale1              + y[11] * scale2; // yyy
      current_data[4] = x[6] * z3s                  + x[7] * z4s                  + x[8] * z5s; // xxz
      current_data[5] = x3y3 * z3s                  + x4y4 * z4s                  + x5y5 * z5s; // xyz
      current_data[6] = y[6] * z3s                  + y[7] * z4s                  + y[8] * z5s; // yyz
      current_data[7] = x[3] * z6s                  + x[4] * z7s                  + x[5] * z8s; // xzz
      current_data[8] = y[3] * z6s                  + y[4] * z7s                  + y[5] * z8s; // yzz
      current_data[9] = z[9] * scale0               + z[10] * scale1              + z[11] * scale2; // zzz

      current_data[10] = x[12] * scale0              + x[13] * scale1              + x[14] * scale2; // xxxx
      current_data[11] = x[9] * y3s                  + x[10] * y4s                 + x[11] * y5s; // xxxy
      current_data[12] = x[6] * y6s                  + x[7] * y7s                  + x[8] * y8s; // xxyy
      current_data[13] = y[9] * x3s                  + y[10] * x4s                 + y[11] * x5s; // xyyy
      current_data[14] = y[12] * scale0              + y[13] * scale1              + y[14] * scale2; // yyyy
      current_data[15] = x[9] * z3s                  + x[10] * z4s                 + x[11] * z5s; // xxxz
      current_data[16] = x6y3 * z3s                  + x7y4 * z4s                  + x8y5 * z5s; // xxyz
      current_data[17] = x3y6 * z3s                  + x4y7 * z4s                  + x5y8 * z5s; // xyyz
      current_data[18] = y[9] * z3s                  + y[10] * z4s                 + y[11] * z5s; // yyyz
      current_data[19] = x[6] * z6s                  + x[7] * z7s                  + x[8] * z8s; // xxzz
      current_data[20] = x3y3 * z6s                  + x4y4 * z7s                  + x5y5 * z8s; // xyzz
      current_data[21] = y[6] * z6s                  + y[7] * z7s                  + y[8] * z8s; // yyzz
      current_data[22] = z[9] * x3s                  + z[10] * x4s                 + z[11] * x5s; // xzzz
      current_data[23] = z[9] * y3s                  + z[10] * y4s                 + z[11] * y5s; // yzzz
      current_data[24] = z[12] * scale0              + z[13] * scale1              + z[14] * scale2; // yyyy

    }
  } else if (ang == 97) { // (0, 0, 3, 1)
    double x[15];
    double y[15];
    double z[15];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;

      const double opq = 1.0 / (xp_[ii] + xq_[ii]);
      const double oxq2 = 0.5 / xq_[ii];
      const double xpopq = xp_[ii] * opq;
      const double d00i0x = Q_[ii3] - cx;
      const double d00i0y = Q_[ii3 + 1] - cy;
      const double d00i0z = Q_[ii3 + 2] - cz;
      const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
      const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
      const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
      const double b01i0 = xpopq * oxq2;
      const double b01_0 = oxq2 - b01i0 * roots_[offset];
      const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
      const double b01_2 = oxq2 - b01i0 * roots_[offset + 2];

      x[3]  = d00i0x + d00i1x * roots_[offset];
      x[4]  = d00i0x + d00i1x * roots_[offset + 1];
      x[5]  = d00i0x + d00i1x * roots_[offset + 2];
      x[6] = x[3] * x[3] + b01_0;
      x[7] = x[4] * x[4] + b01_1;
      x[8] = x[5] * x[5] + b01_2;
      x[9]  = x[3] * (x[6] + 2.0 * b01_0);
      x[10] = x[4] * (x[7] + 2.0 * b01_1);
      x[11] = x[5] * (x[8] + 2.0 * b01_2);
      x[12] = x[3] * x[9]  + 3.0 * b01_0 * x[6];
      x[13] = x[4] * x[10] + 3.0 * b01_1 * x[7];
      x[14] = x[5] * x[11] + 3.0 * b01_2 * x[8];

      y[3]  = d00i0y + d00i1y * roots_[offset];
      y[4]  = d00i0y + d00i1y * roots_[offset + 1];
      y[5]  = d00i0y + d00i1y * roots_[offset + 2];
      y[6] = y[3] * y[3] + b01_0;
      y[7] = y[4] * y[4] + b01_1;
      y[8] = y[5] * y[5] + b01_2;
      y[9]  = y[3] * (y[6] + 2.0 * b01_0);
      y[10] = y[4] * (y[7] + 2.0 * b01_1);
      y[11] = y[5] * (y[8] + 2.0 * b01_2);
      y[12] = y[3] * y[9]  + 3.0 * b01_0 * y[6];
      y[13] = y[4] * y[10] + 3.0 * b01_1 * y[7];
      y[14] = y[5] * y[11] + 3.0 * b01_2 * y[8];

      z[3]  = d00i0z + d00i1z * roots_[offset];
      z[4]  = d00i0z + d00i1z * roots_[offset + 1];
      z[5]  = d00i0z + d00i1z * roots_[offset + 2];
      z[6] = z[3] * z[3] + b01_0;
      z[7] = z[4] * z[4] + b01_1;
      z[8] = z[5] * z[5] + b01_2;
      z[9]  = z[3] * (z[6] + 2.0 * b01_0);
      z[10] = z[4] * (z[7] + 2.0 * b01_1);
      z[11] = z[5] * (z[8] + 2.0 * b01_2);
      z[12] = z[3] * z[9]  + 3.0 * b01_0 * z[6];
      z[13] = z[4] * z[10] + 3.0 * b01_1 * z[7];
      z[14] = z[5] * z[11] + 3.0 * b01_2 * z[8];

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double scale2 = coeff_[ii] * weights_[offset + 2];

      const double x3s = x[3] * scale0;
      const double x4s = x[4] * scale1;
      const double x5s = x[5] * scale2;
      const double y3s = y[3] * scale0;
      const double y4s = y[4] * scale1;
      const double y5s = y[5] * scale2;
      const double y6s = y[6] * scale0;
      const double y7s = y[7] * scale1;
      const double y8s = y[8] * scale2;
      const double z3s = z[3] * scale0;
      const double z4s = z[4] * scale1;
      const double z5s = z[5] * scale2;
      const double z6s = z[6] * scale0;
      const double z7s = z[7] * scale1;
      const double z8s = z[8] * scale2;

      const double x3y3 = x[3] * y[3];
      const double x4y4 = x[4] * y[4];
      const double x5y5 = x[5] * y[5];
      const double x3y6 = x[3] * y[6];
      const double x4y7 = x[4] * y[7];
      const double x5y8 = x[5] * y[8];
      const double x6y3 = y[3] * x[6];
      const double x7y4 = y[4] * x[7];
      const double x8y5 = y[5] * x[8];

      current_data[0] = x[9] * scale0               + x[10] * scale1              + x[11] * scale2; // xxx
      current_data[1] = x6y3 * scale0               + x7y4 * scale1               + x8y5 * scale2; // xxy
      current_data[2] = x3y6 * scale0               + x4y7 * scale1               + x5y8 * scale2; // xyy
      current_data[3] = y[9] * scale0               + y[10] * scale1              + y[11] * scale2; // yyy
      current_data[4] = x[6] * z3s                  + x[7] * z4s                  + x[8] * z5s; // xxz
      current_data[5] = x3y3 * z3s                  + x4y4 * z4s                  + x5y5 * z5s; // xyz
      current_data[6] = y[6] * z3s                  + y[7] * z4s                  + y[8] * z5s; // yyz
      current_data[7] = x[3] * z6s                  + x[4] * z7s                  + x[5] * z8s; // xzz
      current_data[8] = y[3] * z6s                  + y[4] * z7s                  + y[5] * z8s; // yzz
      current_data[9] = z[9] * scale0               + z[10] * scale1              + z[11] * scale2; // zzz

      current_data[10] = x[12] * scale0              + x[13] * scale1              + x[14] * scale2; // xxxx
      current_data[11] = x[9] * y3s                  + x[10] * y4s                 + x[11] * y5s; // xxxy
      current_data[12] = x[6] * y6s                  + x[7] * y7s                  + x[8] * y8s; // xxyy
      current_data[13] = y[9] * x3s                  + y[10] * x4s                 + y[11] * x5s; // xyyy
      current_data[14] = y[12] * scale0              + y[13] * scale1              + y[14] * scale2; // yyyy
      current_data[15] = x[9] * z3s                  + x[10] * z4s                 + x[11] * z5s; // xxxz
      current_data[16] = x6y3 * z3s                  + x7y4 * z4s                  + x8y5 * z5s; // xxyz
      current_data[17] = x3y6 * z3s                  + x4y7 * z4s                  + x5y8 * z5s; // xyyz
      current_data[18] = y[9] * z3s                  + y[10] * z4s                 + y[11] * z5s; // yyyz
      current_data[19] = x[6] * z6s                  + x[7] * z7s                  + x[8] * z8s; // xxzz
      current_data[20] = x3y3 * z6s                  + x4y4 * z7s                  + x5y5 * z8s; // xyzz
      current_data[21] = y[6] * z6s                  + y[7] * z7s                  + y[8] * z8s; // yyzz
      current_data[22] = z[9] * x3s                  + z[10] * x4s                 + z[11] * x5s; // xzzz
      current_data[23] = z[9] * y3s                  + z[10] * y4s                 + z[11] * y5s; // yzzz
      current_data[24] = z[12] * scale0              + z[13] * scale1              + z[14] * scale2; // yyyy

    }
  } else if (ang == 66) { // (0, 0, 2, 2)
    double x[15];
    double y[27];
    double z[27];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;

      const double opq = 1.0 / (xp_[ii] + xq_[ii]);
      const double oxq2 = 0.5 / xq_[ii];
      const double xpopq = xp_[ii] * opq;
      const double d00i0x = Q_[ii3] - cx;
      const double d00i0y = Q_[ii3 + 1] - cy;
      const double d00i0z = Q_[ii3 + 2] - cz;
      const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
      const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
      const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
      const double b01i0 = xpopq * oxq2;
      const double b01_0 = oxq2 - b01i0 * roots_[offset];
      const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
      const double b01_2 = oxq2 - b01i0 * roots_[offset + 2];

      x[3]  = d00i0x + d00i1x * roots_[offset];
      x[4]  = d00i0x + d00i1x * roots_[offset + 1];
      x[5]  = d00i0x + d00i1x * roots_[offset + 2];
      x[6] = x[3] * x[3] + b01_0;
      x[7] = x[4] * x[4] + b01_1;
      x[8] = x[5] * x[5] + b01_2;
      x[9]  = x[3] * (x[6] + 2.0 * b01_0);
      x[10] = x[4] * (x[7] + 2.0 * b01_1);
      x[11] = x[5] * (x[8] + 2.0 * b01_2);
      x[12] = x[3] * x[9]  + 3.0 * b01_0 * x[6];
      x[13] = x[4] * x[10] + 3.0 * b01_1 * x[7];
      x[14] = x[5] * x[11] + 3.0 * b01_2 * x[8];

      y[3]  = d00i0y + d00i1y * roots_[offset];
      y[4]  = d00i0y + d00i1y * roots_[offset + 1];
      y[5]  = d00i0y + d00i1y * roots_[offset + 2];
      y[6] = y[3] * y[3] + b01_0;
      y[7] = y[4] * y[4] + b01_1;
      y[8] = y[5] * y[5] + b01_2;
      y[9]  = y[3] * (y[6] + 2.0 * b01_0);
      y[10] = y[4] * (y[7] + 2.0 * b01_1);
      y[11] = y[5] * (y[8] + 2.0 * b01_2);
      y[12] = y[3] * y[9]  + 3.0 * b01_0 * y[6];
      y[13] = y[4] * y[10] + 3.0 * b01_1 * y[7];
      y[14] = y[5] * y[11] + 3.0 * b01_2 * y[8];

      z[3]  = d00i0z + d00i1z * roots_[offset];
      z[4]  = d00i0z + d00i1z * roots_[offset + 1];
      z[5]  = d00i0z + d00i1z * roots_[offset + 2];
      z[6] = z[3] * z[3] + b01_0;
      z[7] = z[4] * z[4] + b01_1;
      z[8] = z[5] * z[5] + b01_2;
      z[9]  = z[3] * (z[6] + 2.0 * b01_0);
      z[10] = z[4] * (z[7] + 2.0 * b01_1);
      z[11] = z[5] * (z[8] + 2.0 * b01_2);
      z[12] = z[3] * z[9]  + 3.0 * b01_0 * z[6];
      z[13] = z[4] * z[10] + 3.0 * b01_1 * z[7];
      z[14] = z[5] * z[11] + 3.0 * b01_2 * z[8];

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double scale2 = coeff_[ii] * weights_[offset + 2];

      const double x3s = x[3] * scale0;
      const double x4s = x[4] * scale1;
      const double x5s = x[5] * scale2;
      const double y3s = y[3] * scale0;
      const double y4s = y[4] * scale1;
      const double y5s = y[5] * scale2;
      const double y6s = y[6] * scale0;
      const double y7s = y[7] * scale1;
      const double y8s = y[8] * scale2;
      const double z3s = z[3] * scale0;
      const double z4s = z[4] * scale1;
      const double z5s = z[5] * scale2;
      const double z6s = z[6] * scale0;
      const double z7s = z[7] * scale1;
      const double z8s = z[8] * scale2;

      const double x3y3 = x[3] * y[3];
      const double x4y4 = x[4] * y[4];
      const double x5y5 = x[5] * y[5];
      const double x3y6 = x[3] * y[6];
      const double x4y7 = x[4] * y[7];
      const double x5y8 = x[5] * y[8];
      const double x6y3 = y[3] * x[6];
      const double x7y4 = y[4] * x[7];
      const double x8y5 = y[5] * x[8];

      current_data[0]  = x[6] * scale0        + x[7] * scale1        + x[8] * scale2; // xx
      current_data[1]  = x[3] * y3s           + x[4] * y4s           + x[5] * y5s; // xy
      current_data[2]  = y6s                  + y7s                  + y8s; // yy
      current_data[3]  = x[3] * z3s           + x[4] * z4s           + x[5] * z5s; // xz
      current_data[4]  = y[3] * z3s           + y[4] * z4s           + y[5] * z5s; // yz
      current_data[5]  = z6s                  + z7s                  + z8s; // zz

      current_data[6] = x[9] * scale0               + x[10] * scale1              + x[11] * scale2; // xxx
      current_data[7] = x6y3 * scale0               + x7y4 * scale1               + x8y5 * scale2; // xxy
      current_data[8] = x3y6 * scale0               + x4y7 * scale1               + x5y8 * scale2; // xyy
      current_data[9] = y[9] * scale0               + y[10] * scale1              + y[11] * scale2; // yyy
      current_data[10] = x[6] * z3s                  + x[7] * z4s                  + x[8] * z5s; // xxz
      current_data[11] = x3y3 * z3s                  + x4y4 * z4s                  + x5y5 * z5s; // xyz
      current_data[12] = y[6] * z3s                  + y[7] * z4s                  + y[8] * z5s; // yyz
      current_data[13] = x[3] * z6s                  + x[4] * z7s                  + x[5] * z8s; // xzz
      current_data[14] = y[3] * z6s                  + y[4] * z7s                  + y[5] * z8s; // yzz
      current_data[15] = z[9] * scale0               + z[10] * scale1              + z[11] * scale2; // zzz

      current_data[16] = x[12] * scale0              + x[13] * scale1              + x[14] * scale2; // xxxx
      current_data[17] = x[9] * y3s                  + x[10] * y4s                 + x[11] * y5s; // xxxy
      current_data[18] = x[6] * y6s                  + x[7] * y7s                  + x[8] * y8s; // xxyy
      current_data[19] = y[9] * x3s                  + y[10] * x4s                 + y[11] * x5s; // xyyy
      current_data[20] = y[12] * scale0              + y[13] * scale1              + y[14] * scale2; // yyyy
      current_data[21] = x[9] * z3s                  + x[10] * z4s                 + x[11] * z5s; // xxxz
      current_data[22] = x6y3 * z3s                  + x7y4 * z4s                  + x8y5 * z5s; // xxyz
      current_data[23] = x3y6 * z3s                  + x4y7 * z4s                  + x5y8 * z5s; // xyyz
      current_data[24] = y[9] * z3s                  + y[10] * z4s                 + y[11] * z5s; // yyyz
      current_data[25] = x[6] * z6s                  + x[7] * z7s                  + x[8] * z8s; // xxzz
      current_data[26] = x3y3 * z6s                  + x4y4 * z7s                  + x5y5 * z8s; // xyzz
      current_data[27] = y[6] * z6s                  + y[7] * z7s                  + y[8] * z8s; // yyzz
      current_data[28] = z[9] * x3s                  + z[10] * x4s                 + z[11] * x5s; // xzzz
      current_data[29] = z[9] * y3s                  + z[10] * y4s                 + z[11] * y5s; // yzzz
      current_data[30] = z[12] * scale0              + z[13] * scale1              + z[14] * scale2; // yyyy

    }

  } else {
// (4, 0, 0, 0) and (0, 0, 0, 4) hasn't been rewritten (using the generic code);
// (3, 0, 1, 0) and (1, 0, 3, 0) hasn't been rewritten yet
// (2, 1, 1, 0) and (1, 0, 2, 1) hasn't been rewritten yet

    perform_VRR();

  }
}


void ERIBatch::perform_VRR2() { // fully hand-written

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);

  const unsigned int ang0 = basisinfo_[0]->angular_number();
  const unsigned int ang1 = basisinfo_[1]->angular_number();
  const unsigned int ang2 = basisinfo_[2]->angular_number();
  const unsigned int ang3 = basisinfo_[3]->angular_number();
  const unsigned int ang = (((((ang0 << 5) + ang1) << 5) + ang2) << 5) + ang3; // offsets are 2^5 = 32

  if (ang == 65536) { // (2, 0, 0, 0)
    double x2, x3, x4, x5;
    double y2, y3, y4, y5;
    double z2, z3, z4, z5;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      {
        const double oxp2 = 0.5 / xp_[ii];
        const double xqopq = xq_[ii] / (xp_[ii] + xq_[ii]);
        const double c00i0x = P_[ii3]     - ax;
        const double c00i0y = P_[ii3 + 1] - ay;
        const double c00i0z = P_[ii3 + 2] - az;
        const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
        const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
        const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
        const double b10i0 = xqopq * oxp2;
        const double b10_0 = oxp2 - b10i0 * roots_[offset];
        const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
        x2 = c00i0x - c00i1x * roots_[offset];
        x3 = c00i0x - c00i1x * roots_[offset + 1];
        x4 = x2 * x2 + b10_0;
        x5 = x3 * x3 + b10_1;

        y2 = c00i0y - c00i1y * roots_[offset];
        y3 = c00i0y - c00i1y * roots_[offset + 1];
        y4 = y2 * y2 + b10_0;
        y5 = y3 * y3 + b10_1;

        z2 = c00i0z - c00i1z * roots_[offset];
        z3 = c00i0z - c00i1z * roots_[offset + 1];
        z4 = z2 * z2 + b10_0;
        z5 = z3 * z3 + b10_1;
      }

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double y2s0 = y2 * scale0;
      const double y3s1 = y3 * scale1;
      current_data[0] = scale0 * x4      + scale1 * x5      ; // xx
      current_data[1] = y2s0   * x2      + y3s1   * x3      ; // xy
      current_data[2] = scale0 * y4      + scale1 * y5      ; // yy
      current_data[3] = scale0 * x2 * z2 + scale1 * x3 * z3 ; // xz
      current_data[4] = y2s0   * z2      + y3s1   * z3      ; // yz
      current_data[5] = scale0 * z4      + scale1 * z5      ; // zz
    }
  } else if (ang == 64) { // (0, 0, 2, 0)
    double x2, x3, x4, x5;
    double y2, y3, y4, y5;
    double z2, z3, z4, z5;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      {
        const double oxq2 = 0.5 / xq_[ii];
        const double xpopq = xp_[ii] / (xp_[ii] + xq_[ii]);
        const double d00i0x = Q_[ii3] - cx;
        const double d00i0y = Q_[ii3 + 1] - cy;
        const double d00i0z = Q_[ii3 + 2] - cz;
        const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
        const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
        const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
        const double b01i0 = xpopq * oxq2;
        const double b01_0 = oxq2 - b01i0 * roots_[offset];
        const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
        x2 = d00i0x + d00i1x * roots_[offset];
        x3 = d00i0x + d00i1x * roots_[offset + 1];
        x4 = x2 * x2 + b01_0;
        x5 = x3 * x3 + b01_1;

        y2 = d00i0y + d00i1y * roots_[offset];
        y3 = d00i0y + d00i1y * roots_[offset + 1];
        y4 = y2 * y2 + b01_0;
        y5 = y3 * y3 + b01_1;

        z2 = d00i0z + d00i1z * roots_[offset];
        z3 = d00i0z + d00i1z * roots_[offset + 1];
        z4 = z2 * z2 + b01_0;
        z5 = z3 * z3 + b01_1;
      }
      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double y2s0 = y2 * scale0;
      const double y3s1 = y3 * scale1;
      const double z2s0 = z2 * scale0;
      const double z3s1 = z3 * scale1;
      current_data[0] = scale0 * x4      + scale1 * x5      ; // xx
      current_data[1] = y2s0   * x2      + y3s1   * x3      ; // xy
      current_data[2] = scale0 * y4      + scale1 * y5      ; // yy
      current_data[3] = z2s0   * x2      + z3s1   * x3      ; // xz
      current_data[4] = y2s0   * z2      + y3s1   * z3      ; // yz
      current_data[5] = scale0 * z4      + scale1 * z5      ; // zz
    }
  } else if (ang == 32800) { // (1, 0, 1, 0)
    double x2, x3, x4, x5, x6, x7;
    double y2, y3, y4, y5, y6, y7;
    double z2, z3, z4, z5, z6, z7;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      {
        const double opq = 1.0 / (xp_[ii] + xq_[ii]);
        const double xqopq = xq_[ii] * opq;
        const double c00i0x = P_[ii3]     - ax;
        const double c00i0y = P_[ii3 + 1] - ay;
        const double c00i0z = P_[ii3 + 2] - az;
        const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
        const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
        const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;

        const double xpopq = xp_[ii] * opq;
        const double d00i0x = Q_[ii3] - cx;
        const double d00i0y = Q_[ii3 + 1] - cy;
        const double d00i0z = Q_[ii3 + 2] - cz;
        const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
        const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
        const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;

        const double b00_0 = 0.5 * opq * roots_[offset];
        const double b00_1 = 0.5 * opq * roots_[offset + 1];

        x2 = c00i0x - c00i1x * roots_[offset];
        x3 = c00i0x - c00i1x * roots_[offset + 1];
        x4 = d00i0x + d00i1x * roots_[offset];
        x5 = d00i0x + d00i1x * roots_[offset + 1];
        x6  = x2 * x4 + b00_0;
        x7 = x3 * x5 + b00_1;

        y2 = c00i0y - c00i1y * roots_[offset];
        y3 = c00i0y - c00i1y * roots_[offset + 1];
        y4 = d00i0y + d00i1y * roots_[offset];
        y5 = d00i0y + d00i1y * roots_[offset + 1];
        y6 = y2 * y4 + b00_0;
        y7 = y3 * y5 + b00_1;

        z2 = c00i0z - c00i1z * roots_[offset];
        z3 = c00i0z - c00i1z * roots_[offset + 1];
        z4 = d00i0z + d00i1z * roots_[offset];
        z5 = d00i0z + d00i1z * roots_[offset + 1];
        z6  = z2 * z4 + b00_0;
        z7 = z3 * z5 + b00_1;
      }

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double x2s0 = x2 * scale0;
      const double x3s1 = x3 * scale1;
      const double x4s0 = x4 * scale0;
      const double x5s1 = x5 * scale1;
      current_data[0] = x6 *      scale0 + x7 *      scale1 ; // x|x
      current_data[1] = x4s0 * y2        + x5s1 * y3        ; // x|y
      current_data[2] = x4s0 * z2        + x5s1 * z3        ; // x|z
      current_data[3] = y4 * x2s0        + y5 * x3s1        ; // y|x
      current_data[4] = y6 *      scale0 + y7 *      scale1 ; // y|y
      current_data[5] = y4 * z2 * scale0 + y5 * z3 * scale1 ; // y|z
      current_data[6] = z4 * x2s0        + z5 * x3s1        ; // z|x
      current_data[7] = z4 * y2 * scale0 + z5 * y3 * scale1 ; // z|y
      current_data[8] = z6 *      scale0 + z7 *      scale1 ; // z|z
    }
  } else if (ang == 33792) { // (1, 1, 0, 0)
    double x2, x3, x4, x5;
    double y2, y3, y4, y5;
    double z2, z3, z4, z5;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];
      const int ii3 = 3 * ii;
      {
        const double oxp2 = 0.5 / xp_[ii];
        const double xqopq = xq_[ii] / (xp_[ii] + xq_[ii]);
        const double c00i0x = P_[ii3]     - ax;
        const double c00i0y = P_[ii3 + 1] - ay;
        const double c00i0z = P_[ii3 + 2] - az;
        const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
        const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
        const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
        const double b10i0 = xqopq * oxp2;
        const double b10_0 = oxp2 - b10i0 * roots_[offset];
        const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
        x2 = c00i0x - c00i1x * roots_[offset];
        x3 = c00i0x - c00i1x * roots_[offset + 1];
        x4 = x2 * x2 + b10_0;
        x5 = x3 * x3 + b10_1;

        y2 = c00i0y - c00i1y * roots_[offset];
        y3 = c00i0y - c00i1y * roots_[offset + 1];
        y4 = y2 * y2 + b10_0;
        y5 = y3 * y3 + b10_1;

        z2 = c00i0z - c00i1z * roots_[offset];
        z3 = c00i0z - c00i1z * roots_[offset + 1];
        z4 = z2 * z2 + b10_0;
        z5 = z3 * z3 + b10_1;
      }

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double y2s0 = y2 * scale0;
      const double y3s1 = y3 * scale1;
      const double z2s0 = z2 * scale0;
      const double z3s1 = z3 * scale1;
      current_data[0] = scale0 * x2      + scale1 * x3      ; // x
      current_data[1] = y2s0             + y3s1             ; // y
      current_data[2] = z2s0             + z3s1             ; // z
      current_data[3] = scale0 * x4      + scale1 * x5      ; // xx
      current_data[4] = y2s0   * x2      + y3s1   * x3      ; // xy
      current_data[5] = scale0 * y4      + scale1 * y5      ; // yy
      current_data[6] = scale0 * x2 * z2 + scale1 * x3 * z3 ; // xz
      current_data[7] = y2s0   * z2      + y3s1   * z3      ; // yz
      current_data[8] = scale0 * z4      + scale1 * z5      ; // zz
    }
  } else if (ang == 33) { // (0, 0, 1, 1)
    double x2, x3, x4, x5;
    double y2, y3, y4, y5;
    double z2, z3, z4, z5;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      {
        const double oxq2 = 0.5 / xq_[ii];
        const double xpopq = xp_[ii] / (xp_[ii] + xq_[ii]);
        const double d00i0x = Q_[ii3] - cx;
        const double d00i0y = Q_[ii3 + 1] - cy;
        const double d00i0z = Q_[ii3 + 2] - cz;
        const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
        const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
        const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
        const double b01i0 = xpopq * oxq2;
        const double b01_0 = oxq2 - b01i0 * roots_[offset];
        const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
        x2 = d00i0x + d00i1x * roots_[offset];
        x3 = d00i0x + d00i1x * roots_[offset + 1];
        x4 = x2 * x2 + b01_0;
        x5 = x3 * x3 + b01_1;

        y2 = d00i0y + d00i1y * roots_[offset];
        y3 = d00i0y + d00i1y * roots_[offset + 1];
        y4 = y2 * y2 + b01_0;
        y5 = y3 * y3 + b01_1;

        z2 = d00i0z + d00i1z * roots_[offset];
        z3 = d00i0z + d00i1z * roots_[offset + 1];
        z4 = z2 * z2 + b01_0;
        z5 = z3 * z3 + b01_1;
      }
      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double y2s0 = y2 * scale0;
      const double y3s1 = y3 * scale1;
      const double z2s0 = z2 * scale0;
      const double z3s1 = z3 * scale1;
      current_data[0] = scale0 * x2      + scale1 * x3      ; // x
      current_data[1] = y2s0             + y3s1             ; // y
      current_data[2] = z2s0             + z3s1             ; // z
      current_data[3] = scale0 * x4      + scale1 * x5      ; // xx
      current_data[4] = y2s0   * x2      + y3s1   * x3      ; // xy
      current_data[5] = scale0 * y4      + scale1 * y5      ; // yy
      current_data[6] = z2s0   * x2      + z3s1   * x3      ; // xz
      current_data[7] = y2s0   * z2      + y3s1   * z3      ; // yz
      current_data[8] = scale0 * z4      + scale1 * z5      ; // zz
    }
  } else if (ang == 65568) { // (2, 0, 1, 0)
    double x2, x3, x4, x5, x6, x7, x8, x9, x10, x11;
    double y2, y3, y4, y5, y6, y7, y8, y9, y10, y11;
    double z2, z3, z4, z5, z6, z7, z8, z9, z10, z11;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      {
        const double oxp2 = 0.5 / xp_[ii];
        const double opq = 1.0 / (xp_[ii] + xq_[ii]);
        const double xqopq = xq_[ii] * opq;
        const double c00i0x = P_[ii3]     - ax;
        const double c00i0y = P_[ii3 + 1] - ay;
        const double c00i0z = P_[ii3 + 2] - az;
        const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
        const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
        const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
        const double b10i0 = xqopq * oxp2;
        const double b10_0 = oxp2 - b10i0 * roots_[offset];
        const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
        const double xpopq = xp_[ii] * opq;
        const double d00i0x = Q_[ii3] - cx;
        const double d00i0y = Q_[ii3 + 1] - cy;
        const double d00i0z = Q_[ii3 + 2] - cz;
        const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
        const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
        const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
        const double b00int = 0.5 * opq;
        const double b00_0 = b00int * roots_[offset];
        const double b00_1 = b00int * roots_[offset + 1];
        x2 = c00i0x - c00i1x * roots_[offset];
        x3 = c00i0x - c00i1x * roots_[offset + 1];
        x4 = x2 * x2 + b10_0;
        x5 = x3 * x3 + b10_1;
        x6 = d00i0x + d00i1x * roots_[offset];
        x7 = d00i0x + d00i1x * roots_[offset + 1];
        x8 = x2 * x6 + b00_0;
        x9 = x3 * x7 + b00_1;
        x10 = x2 * (x8 + b00_0) + b10_0 * x6;
        x11 = x3 * (x9 + b00_1) + b10_1 * x7;

        y2 = c00i0y - c00i1y * roots_[offset];
        y3 = c00i0y - c00i1y * roots_[offset + 1];
        y4 = y2 * y2 + b10_0;
        y5 = y3 * y3 + b10_1;
        y6 = d00i0y + d00i1y * roots_[offset];
        y7 = d00i0y + d00i1y * roots_[offset + 1];
        y8 = y2 * y6 + b00_0;
        y9 = y3 * y7 + b00_1;
        y10 = y2 * (y8 + b00_0) + b10_0 * y6;
        y11 = y3 * (y9 + b00_1) + b10_1 * y7;

        z2 = c00i0z - c00i1z * roots_[offset];
        z3 = c00i0z - c00i1z * roots_[offset + 1];
        z4 = z2 * z2 + b10_0;
        z5 = z3 * z3 + b10_1;
        z6 = d00i0z + d00i1z * roots_[offset];
        z7 = d00i0z + d00i1z * roots_[offset + 1];
        z8 = z2 * z6 + b00_0;
        z9 = z3 * z7 + b00_1;
        z10 = z2 * (z8 + b00_0) + b10_0 * z6;
        z11 = z3 * (z9 + b00_1) + b10_1 * z7;

        const double s0 = coeff_[ii] * weights_[offset];
        const double s1 = coeff_[ii] * weights_[offset + 1];

        const double s0x6 = s0 * x6;
        const double s1x7 = s1 * x7;
        const double s0x8 = s0 * x8;
        const double s1x9 = s1 * x9;
        const double s0x10 = s0 * x10;
        const double s1x11 = s1 * x11;
        current_data[0] = s0x10            + s1x11            ; // x|xx
        current_data[1] = s0x8 * y2        + s1x9 * y3        ; // x|xy
        current_data[2] = s0x6 * y4        + s1x7 * y5        ; // x|yy
        current_data[3] = s0x8 * z2        + s1x9 * z3        ; // x|xz
        current_data[4] = s0x6 * y2 * z2   + s1x7 * y3 * z3   ; // x|yz
        current_data[5] = s0x6 * z4        + s1x7 * z5        ; // x|zz

        const double s0y6 = s0 * y6;
        const double s1y7 = s1 * y7;
        const double s0y8 = s0 * y8;
        const double s1y9 = s1 * y9;
        const double s0y10 = s0 * y10;
        const double s1y11 = s1 * y11;
        current_data[6] =  s0y6 * x4        + s1y7 * x5        ; // y|xx
        current_data[7] =  s0y8 * x2        + s1y9 * x3        ; // y|xy
        current_data[8] =  s0y10            + s1y11            ; // y|yy
        current_data[9] =  s0y6 * x2 * z2   + s1y7 * x3 * z3   ; // y|xz
        current_data[10] = s0y8 * z2        + s1y9 * z3        ; // y|yz
        current_data[11] = s0y6 * z4        + s1y7 * z5        ; // y|zz

        const double s0z6 = s0 * z6;
        const double s1z7 = s1 * z7;
        const double s0z8 = s0 * z8;
        const double s1z9 = s1 * z9;
        const double s0z10 = s0 * z10;
        const double s1z11 = s1 * z11;
        current_data[12] = s0z6 * x4        + s1z7 * x5        ; // z|xx
        current_data[13] = s0z6 * x2 * y2   + s1z7 * x3 * y3   ; // z|xy
        current_data[14] = s0z6 * y4        + s1z7 * y5        ; // z|yy
        current_data[15] = s0z8 * x2        + s1z9 * x3        ; // z|xz
        current_data[16] = s0z8 * y2        + s1z9 * y3        ; // z|yz
        current_data[17] = s0z10            + s1z11            ; // z|zz
      }
    }
  } else if (ang == 33824) { // (1, 1, 1, 0)
    double x[12];
    double y[12];
    double z[12];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      {
        const double oxp2 = 0.5 / xp_[ii];
        const double opq = 1.0 / (xp_[ii] + xq_[ii]);
        const double xqopq = xq_[ii] * opq;
        const double c00i0x = P_[ii3]     - ax;
        const double c00i0y = P_[ii3 + 1] - ay;
        const double c00i0z = P_[ii3 + 2] - az;
        const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
        const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
        const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
        const double b10i0 = xqopq * oxp2;
        const double b10_0 = oxp2 - b10i0 * roots_[offset];
        const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
        const double xpopq = xp_[ii] * opq;
        const double d00i0x = Q_[ii3] - cx;
        const double d00i0y = Q_[ii3 + 1] - cy;
        const double d00i0z = Q_[ii3 + 2] - cz;
        const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
        const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
        const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
        const double b00int = 0.5 * opq;
        const double b00_0 = b00int * roots_[offset];
        const double b00_1 = b00int * roots_[offset + 1];
        x[2] = c00i0x - c00i1x * roots_[offset];
        x[3] = c00i0x - c00i1x * roots_[offset + 1];
        x[4] = x[2] * x[2] + b10_0;
        x[5] = x[3] * x[3] + b10_1;
        x[6] = d00i0x + d00i1x * roots_[offset];
        x[7] = d00i0x + d00i1x * roots_[offset + 1];
        x[8] = x[2] * x[6] + b00_0;
        x[9] = x[3] * x[7] + b00_1;
        x[10] = x[2] * (x[8] + b00_0) + b10_0 * x[6];
        x[11] = x[3] * (x[9] + b00_1) + b10_1 * x[7];

        y[2] = c00i0y - c00i1y * roots_[offset];
        y[3] = c00i0y - c00i1y * roots_[offset + 1];
        y[4] = y[2] * y[2] + b10_0;
        y[5] = y[3] * y[3] + b10_1;
        y[6] = d00i0y + d00i1y * roots_[offset];
        y[7] = d00i0y + d00i1y * roots_[offset + 1];
        y[8] = y[2] * y[6] + b00_0;
        y[9] = y[3] * y[7] + b00_1;
        y[10] = y[2] * (y[8] + b00_0) + b10_0 * y[6];
        y[11] = y[3] * (y[9] + b00_1) + b10_1 * y[7];

        z[2] = c00i0z - c00i1z * roots_[offset];
        z[3] = c00i0z - c00i1z * roots_[offset + 1];
        z[4] = z[2] * z[2] + b10_0;
        z[5] = z[3] * z[3] + b10_1;
        z[6] = d00i0z + d00i1z * roots_[offset];
        z[7] = d00i0z + d00i1z * roots_[offset + 1];
        z[8] = z[2] * z[6] + b00_0;
        z[9] = z[3] * z[7] + b00_1;
        z[10] = z[2] * (z[8] + b00_0) + b10_0 * z[6];
        z[11] = z[3] * (z[9] + b00_1) + b10_1 * z[7];

        const double s0 = coeff_[ii] * weights_[offset];
        const double s1 = coeff_[ii] * weights_[offset + 1];
        {
        const double s0x6 = s0 * x[6];
        const double s1x7 = s1 * x[7];
        const double s0x8 = s0 * x[8];
        const double s1x9 = s1 * x[9];
        const double s0x10 = s0 * x[10];
        const double s1x11 = s1 * x[11];
        current_data[0] = s0x8               + s1x9               ; // x|x
        current_data[1] = s0x6 * y[2]        + s1x7 * y[3]        ; // x|y
        current_data[2] = s0x6 * z[2]        + s1x7 * z[3]        ; // x|z
        current_data[3] = s0x10              + s1x11              ; // x|xx
        current_data[4] = s0x8 * y[2]        + s1x9 * y[3]        ; // x|xy
        current_data[5] = s0x6 * y[4]        + s1x7 * y[5]        ; // x|yy
        current_data[6] = s0x8 * z[2]        + s1x9 * z[3]        ; // x|xz
        current_data[7] = s0x6 * y[2] * z[2] + s1x7 * y[3] * z[3] ; // x|yz
        current_data[8] = s0x6 * z[4]        + s1x7 * z[5]        ; // x|zz
        }{
        const double s0y6 = s0 * y[6];
        const double s1y7 = s1 * y[7];
        const double s0y8 = s0 * y[8];
        const double s1y9 = s1 * y[9];
        const double s0y10 = s0 * y[10];
        const double s1y11 = s1 * y[11];
        current_data[9] =  s0y6 * x[2]        + s1y7 * x[3]        ; // y|x
        current_data[10]=  s0y8               + s1y9               ; // y|y
        current_data[11]=  s0y6 * z[2]        + s1y7 * z[3]        ; // y|z
        current_data[12]=  s0y6 * x[4]        + s1y7 * x[5]        ; // y|xx
        current_data[13]=  s0y8 * x[2]        + s1y9 * x[3]        ; // y|xy
        current_data[14]=  s0y10              + s1y11              ; // y|yy
        current_data[15]=  s0y6 * x[2] * z[2] + s1y7 * x[3] * z[3] ; // y|xz
        current_data[16] = s0y8 * z[2]        + s1y9 * z[3]        ; // y|yz
        current_data[17] = s0y6 * z[4]        + s1y7 * z[5]        ; // y|zz
        }{
        const double s0z6 = s0 * z[6];
        const double s1z7 = s1 * z[7];
        const double s0z8 = s0 * z[8];
        const double s1z9 = s1 * z[9];
        const double s0z10 = s0 * z[10];
        const double s1z11 = s1 * z[11];
        current_data[18]=  s0z6 * x[2]        + s1z7 * x[3]        ; // z|x
        current_data[19]=  s0z6 * y[2]        + s1z7 * y[3]        ; // z|z
        current_data[20]=  s0z8               + s1z9               ; // z|y
        current_data[21] = s0z6 * x[4]        + s1z7 * x[5]        ; // z|xx
        current_data[22] = s0z6 * x[2] * y[2] + s1z7 * x[3] * y[3] ; // z|xy
        current_data[23] = s0z6 * y[4]        + s1z7 * y[5]        ; // z|yy
        current_data[24] = s0z8 * x[2]        + s1z9 * x[3]        ; // z|xz
        current_data[25] = s0z8 * y[2]        + s1z9 * y[3]        ; // z|yz
        current_data[26] = s0z10              + s1z11              ; // z|zz
        }
      }
    }
  } else if (ang == 32801) { // (1, 0, 1, 1)
    double x[12];
    double y[12];
    double z[12];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;

      const double opq = 1.0 / (xp_[ii] + xq_[ii]);
      const double xqopq = xq_[ii] * opq;
      const double c00i0x = P_[ii3]     - ax;
      const double c00i0y = P_[ii3 + 1] - ay;
      const double c00i0z = P_[ii3 + 2] - az;
      const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
      const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
      const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
      const double oxq2 = 0.5 / xq_[ii];
      const double xpopq = xp_[ii] * opq;
      const double d00i0x = Q_[ii3] - cx;
      const double d00i0y = Q_[ii3 + 1] - cy;
      const double d00i0z = Q_[ii3 + 2] - cz;
      const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
      const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
      const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
      const double b00int = 0.5 * opq;
      const double b00_0 = b00int * roots_[offset];
      const double b00_1 = b00int * roots_[offset + 1];
      const double b01i0 = xpopq * oxq2;
      const double b01_0 = oxq2 - b01i0 * roots_[offset];
      const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];

      x[2] = c00i0x - c00i1x * roots_[offset];
      x[3] = c00i0x - c00i1x * roots_[offset + 1];
      x[4] = d00i0x + d00i1x * roots_[offset];
      x[5] = d00i0x + d00i1x * roots_[offset + 1];
      x[6] = x[2] * x[4]   + b00_0;
      x[7] = x[3] * x[5]   + b00_1;
      x[8] = x[4] * x[4]   + b01_0;
      x[9] = x[5] * x[5]   + b01_1;
      x[10] = x[2] * x[8]  + 2.0 * b00_0 * x[4];
      x[11] = x[3] * x[9]  + 2.0 * b00_1 * x[5];

      y[2] = c00i0y - c00i1y * roots_[offset];
      y[3] = c00i0y - c00i1y * roots_[offset + 1];
      y[4] = d00i0y + d00i1y * roots_[offset];
      y[5] = d00i0y + d00i1y * roots_[offset + 1];
      y[6] = y[2] * y[4]   + b00_0;
      y[7] = y[3] * y[5]   + b00_1;
      y[8] = y[4] * y[4]   + b01_0;
      y[9] = y[5] * y[5]   + b01_1;
      y[10] = y[2] * y[8]  + 2.0 * b00_0 * y[4];
      y[11] = y[3] * y[9]  + 2.0 * b00_1 * y[5];

      z[2] = c00i0z - c00i1z * roots_[offset];
      z[3] = c00i0z - c00i1z * roots_[offset + 1];
      z[4] = d00i0z + d00i1z * roots_[offset];
      z[5] = d00i0z + d00i1z * roots_[offset + 1];
      z[6] = z[2] * z[4]   + b00_0;
      z[7] = z[3] * z[5]   + b00_1;
      z[8] = z[4] * z[4]   + b01_0;
      z[9] = z[5] * z[5]   + b01_1;
      z[10] = z[2] * z[8]  + 2.0 * b00_0 * z[4];
      z[11] = z[3] * z[9]  + 2.0 * b00_1 * z[5];

      const double s0 = coeff_[ii] * weights_[offset];
      const double s1 = coeff_[ii] * weights_[offset + 1];
      const double s0x2 = s0 * x[2];
      const double s1x3 = s1 * x[3];
      const double s0x4 = s0 * x[4];
      const double s1x5 = s1 * x[5];
      const double s0x6 = s0 * x[6];
      const double s1x7 = s1 * x[7];
      const double s0x8 = s0 * x[8];
      const double s1x9 = s1 * x[9];
      const double s0x10 = s0 * x[10];
      const double s1x11 = s1 * x[11];
      const double s0x2y4 = s0x2 * y[4];
      const double s1x3y5 = s1x3 * y[5];
      const double s0y6 = s0 * y[6];
      const double s1y7 = s1 * y[7];
      const double s0y4 = s0 * y[4];
      const double s1y5 = s1 * y[5];
      const double s0y2 = s0 * y[2];
      const double s1y3 = s1 * y[3];
      current_data[0] = s0x6              + s1x7              ; // x|x
      current_data[1] = s0x4 * y[2]       + s1x5 * y[3]       ; // x|y
      current_data[2] = s0x4 * z[2]       + s1x5 * z[3]       ; // x|z
      current_data[3] = s0x2y4            + s1x3y5            ; // y|x
      current_data[4] = s0y6              + s1y7              ; // y|y
      current_data[5] = s0y4 * z[2]       + s1y5 * z[3]       ; // y|z
      current_data[6] = s0x2 * z[4]       + s1x3 * z[5]       ; // z|x
      current_data[7] = s0y2 * z[4]       + s1y3 * z[5]       ; // z|y
      current_data[8] = s0 * z[6]         + s1 * z[7]         ; // z|z

      current_data[ 9] = s0x10              + s1x11             ; // xx|x
      current_data[10] = s0x8 * y[2]        + s1x9 * y[3]       ; // xx|y
      current_data[11] = s0x8 * z[2]        + s1x9 * z[3]       ; // xx|z
      current_data[12] = s0x6 * y[4]        + s1x7 * y[5]       ; // xy|x
      current_data[13] = s0x4 * y[6]        + s1x5 * y[7]       ; // xy|y
      current_data[14] = s0x4 * y[4] * z[2] + s1x5 * y[5] * z[3]; // xy|z
      current_data[15] = s0x2 * y[8]        + s1x3 * y[9]       ; // yy|x
      current_data[16] = s0 * y[10]         + s1 * y[11]        ; // yy|y
      current_data[17] = s0 * y[8] * z[2]   + s1 * y[9] * z[3]  ; // yy|z
      current_data[18] = s0x6 * z[4]        + s1x7 * z[5]       ; // xz|x
      current_data[19] = s0x4 * z[4] * y[2] + s1x5 * z[5] * y[3]; // xz|y
      current_data[20] = s0x4 * z[6]        + s1x5 * z[7]       ; // xz|z
      current_data[21] = s0x2y4 * z[4]      + s1x3y5 * z[5]     ; // yz|x
      current_data[22] = s0y6 * z[4]        + s1y7 * z[5]       ; // yz|y
      current_data[23] = s0y4 * z[6]        + s1y5 * z[7]       ; // yz|z
      current_data[24] = s0x2 * z[8]        + s1x3 * z[9]       ; // zz|x
      current_data[25] = s0y2 * z[8]        + s1y3 * z[9]       ; // zz|y
      current_data[26] = s0 * z[10]         + s1 * z[11]        ; // zz|z
    }
  } else if (ang == 32832) { // (1, 0, 2, 0)
    double x[12];
    double y[12];
    double z[12];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;

      const double opq = 1.0 / (xp_[ii] + xq_[ii]);
      const double xqopq = xq_[ii] * opq;
      const double c00i0x = P_[ii3]     - ax;
      const double c00i0y = P_[ii3 + 1] - ay;
      const double c00i0z = P_[ii3 + 2] - az;
      const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
      const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
      const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
      const double oxq2 = 0.5 / xq_[ii];
      const double xpopq = xp_[ii] * opq;
      const double d00i0x = Q_[ii3] - cx;
      const double d00i0y = Q_[ii3 + 1] - cy;
      const double d00i0z = Q_[ii3 + 2] - cz;
      const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
      const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
      const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
      const double b00int = 0.5 * opq;
      const double b00_0 = b00int * roots_[offset];
      const double b00_1 = b00int * roots_[offset + 1];
      const double b01i0 = xpopq * oxq2;
      const double b01_0 = oxq2 - b01i0 * roots_[offset];
      const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];

      x[2] = c00i0x - c00i1x * roots_[offset];
      x[3] = c00i0x - c00i1x * roots_[offset + 1];
      x[4] = d00i0x + d00i1x * roots_[offset];
      x[5] = d00i0x + d00i1x * roots_[offset + 1];
      x[6] = x[2] * x[4]   + b00_0;
      x[7] = x[3] * x[5]   + b00_1;
      x[8] = x[4] * x[4]   + b01_0;
      x[9] = x[5] * x[5]   + b01_1;
      x[10] = x[2] * x[8]  + 2.0 * b00_0 * x[4];
      x[11] = x[3] * x[9]  + 2.0 * b00_1 * x[5];

      y[2] = c00i0y - c00i1y * roots_[offset];
      y[3] = c00i0y - c00i1y * roots_[offset + 1];
      y[4] = d00i0y + d00i1y * roots_[offset];
      y[5] = d00i0y + d00i1y * roots_[offset + 1];
      y[6] = y[2] * y[4]   + b00_0;
      y[7] = y[3] * y[5]   + b00_1;
      y[8] = y[4] * y[4]   + b01_0;
      y[9] = y[5] * y[5]   + b01_1;
      y[10] = y[2] * y[8]  + 2.0 * b00_0 * y[4];
      y[11] = y[3] * y[9]  + 2.0 * b00_1 * y[5];

      z[2] = c00i0z - c00i1z * roots_[offset];
      z[3] = c00i0z - c00i1z * roots_[offset + 1];
      z[4] = d00i0z + d00i1z * roots_[offset];
      z[5] = d00i0z + d00i1z * roots_[offset + 1];
      z[6] = z[2] * z[4]   + b00_0;
      z[7] = z[3] * z[5]   + b00_1;
      z[8] = z[4] * z[4]   + b01_0;
      z[9] = z[5] * z[5]   + b01_1;
      z[10] = z[2] * z[8]  + 2.0 * b00_0 * z[4];
      z[11] = z[3] * z[9]  + 2.0 * b00_1 * z[5];

      const double s0 = coeff_[ii] * weights_[offset];
      const double s1 = coeff_[ii] * weights_[offset + 1];
      const double s0x2 = s0 * x[2];
      const double s1x3 = s1 * x[3];
      const double s0x4 = s0 * x[4];
      const double s1x5 = s1 * x[5];
      const double s0x6 = s0 * x[6];
      const double s1x7 = s1 * x[7];
      const double s0x8 = s0 * x[8];
      const double s1x9 = s1 * x[9];
      const double s0x10 = s0 * x[10];
      const double s1x11 = s1 * x[11];
      const double s0x2y4 = s0x2 * y[4];
      const double s1x3y5 = s1x3 * y[5];
      const double s0y6 = s0 * y[6];
      const double s1y7 = s1 * y[7];
      const double s0y4 = s0 * y[4];
      const double s1y5 = s1 * y[5];
      const double s0y2 = s0 * y[2];
      const double s1y3 = s1 * y[3];

      current_data[ 0] = s0x10              + s1x11             ; // xx|x
      current_data[ 1] = s0x8 * y[2]        + s1x9 * y[3]       ; // xx|y
      current_data[ 2] = s0x8 * z[2]        + s1x9 * z[3]       ; // xx|z
      current_data[ 3] = s0x6 * y[4]        + s1x7 * y[5]       ; // xy|x
      current_data[ 4] = s0x4 * y[6]        + s1x5 * y[7]       ; // xy|y
      current_data[ 5] = s0x4 * y[4] * z[2] + s1x5 * y[5] * z[3]; // xy|z
      current_data[ 6] = s0x2 * y[8]        + s1x3 * y[9]       ; // yy|x
      current_data[ 7] = s0 * y[10]         + s1 * y[11]        ; // yy|y
      current_data[ 8] = s0 * y[8] * z[2]   + s1 * y[9] * z[3]  ; // yy|z
      current_data[ 9] = s0x6 * z[4]        + s1x7 * z[5]       ; // xz|x
      current_data[10] = s0x4 * z[4] * y[2] + s1x5 * z[5] * y[3]; // xz|y
      current_data[11] = s0x4 * z[6]        + s1x5 * z[7]       ; // xz|z
      current_data[12] = s0x2y4 * z[4]      + s1x3y5 * z[5]     ; // yz|x
      current_data[13] = s0y6 * z[4]        + s1y7 * z[5]       ; // yz|y
      current_data[14] = s0y4 * z[6]        + s1y5 * z[7]       ; // yz|z
      current_data[15] = s0x2 * z[8]        + s1x3 * z[9]       ; // zz|x
      current_data[16] = s0y2 * z[8]        + s1y3 * z[9]       ; // zz|y
      current_data[17] = s0 * z[10]         + s1 * z[11]        ; // zz|z
    }
  } else if (ang == 98304) { // (3, 0, 0, 0)
    double x2, x3, x4, x5, x6, x7;
    double y2, y3, y4, y5, y6, y7;
    double z2, z3, z4, z5, z6, z7;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      {
        const double oxp2 = 0.5 / xp_[ii];
        const double xqopq = xq_[ii] / (xp_[ii] + xq_[ii]);
        const double c00i0x = P_[ii3]     - ax;
        const double c00i0y = P_[ii3 + 1] - ay;
        const double c00i0z = P_[ii3 + 2] - az;
        const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
        const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
        const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
        const double b10i0 = xqopq * oxp2;
        const double b10_0 = oxp2 - b10i0 * roots_[offset];
        const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
        x2 = c00i0x - c00i1x * roots_[offset];
        x3 = c00i0x - c00i1x * roots_[offset + 1];
        x4 = x2 * x2 + b10_0;
        x5 = x3 * x3 + b10_1;
        x6 = x2 * (x4 + 2.0 * b10_0);
        x7 = x3 * (x5 + 2.0 * b10_1);

        y2 = c00i0y - c00i1y * roots_[offset];
        y3 = c00i0y - c00i1y * roots_[offset + 1];
        y4 = y2 * y2 + b10_0;
        y5 = y3 * y3 + b10_1;
        y6 = y2 * (y4 + 2.0 * b10_0);
        y7 = y3 * (y5 + 2.0 * b10_1);

        z2 = c00i0z - c00i1z * roots_[offset];
        z3 = c00i0z - c00i1z * roots_[offset + 1];
        z4 = z2 * z2 + b10_0;
        z5 = z3 * z3 + b10_1;
        z6 = z2 * (z4 + 2.0 * b10_0);
        z7 = z3 * (z5 + 2.0 * b10_1);
      }

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double s0y2 = y2 * scale0;
      const double s1y3 = y3 * scale1;
      const double s0y4 = y4 * scale0;
      const double s1y5 = y5 * scale1;
      current_data[0] = x6      * scale0  + x7      * scale1  ; // xxx
      current_data[1] = x4 * s0y2         + x5 * s1y3         ; // xxy
      current_data[2] = x2 * s0y4         + x3 * s1y5         ; // xyy
      current_data[3] = y6      * scale0  + y7      * scale1  ; // yyy
      current_data[4] = x4 * z2 * scale0  + x5 * z3 * scale1  ; // xxz
      current_data[5] = x2 * s0y2 * z2    + x3 * s1y3 * z3    ; // xyz
      current_data[6] = s0y4 * z2         + s1y5 * z3         ; // yyz
      current_data[7] = x2 * z4 * scale0  + x3 * z5 * scale1  ; // xzz
      current_data[8] = s0y2 * z4         + s1y3 * z5         ; // yzz
      current_data[9] = z6      * scale0  + z7      * scale1  ; // zzz
    }
  } else if (ang == 66560) { // (2, 1, 0, 0)
    double x2, x3, x4, x5, x6, x7;
    double y2, y3, y4, y5, y6, y7;
    double z2, z3, z4, z5, z6, z7;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      {
        const double oxp2 = 0.5 / xp_[ii];
        const double xqopq = xq_[ii] / (xp_[ii] + xq_[ii]);
        const double c00i0x = P_[ii3]     - ax;
        const double c00i0y = P_[ii3 + 1] - ay;
        const double c00i0z = P_[ii3 + 2] - az;
        const double c00i1x = (P_[ii3]     - Q_[ii3])     * xqopq;
        const double c00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xqopq;
        const double c00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xqopq;
        const double b10i0 = xqopq * oxp2;
        const double b10_0 = oxp2 - b10i0 * roots_[offset];
        const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
        x2 = c00i0x - c00i1x * roots_[offset];
        x3 = c00i0x - c00i1x * roots_[offset + 1];
        x4 = x2 * x2 + b10_0;
        x5 = x3 * x3 + b10_1;
        x6 = x2 * (x4 + 2.0 * b10_0);
        x7 = x3 * (x5 + 2.0 * b10_1);

        y2 = c00i0y - c00i1y * roots_[offset];
        y3 = c00i0y - c00i1y * roots_[offset + 1];
        y4 = y2 * y2 + b10_0;
        y5 = y3 * y3 + b10_1;
        y6 = y2 * (y4 + 2.0 * b10_0);
        y7 = y3 * (y5 + 2.0 * b10_1);

        z2 = c00i0z - c00i1z * roots_[offset];
        z3 = c00i0z - c00i1z * roots_[offset + 1];
        z4 = z2 * z2 + b10_0;
        z5 = z3 * z3 + b10_1;
        z6 = z2 * (z4 + 2.0 * b10_0);
        z7 = z3 * (z5 + 2.0 * b10_1);
      }

      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double s0y2 = y2 * scale0;
      const double s1y3 = y3 * scale1;
      const double s0y4 = y4 * scale0;
      const double s1y5 = y5 * scale1;
      const double s0z2 = z2 * scale0;
      const double s1z3 = z3 * scale1;
      const double s0z4 = z4 * scale0;
      const double s1z5 = z5 * scale1;
      current_data[0] = x4 * scale0       + x5 * scale1       ; // xx
      current_data[1] = x2 * s0y2         + x3 * s1y3         ; // xy
      current_data[2] = s0y4              + s1y5              ; // yy
      current_data[3] = x2 * s0z2         + x3 * s1z3         ; // xz
      current_data[4] = s0y2 * z2         + s1y3 * z3         ; // yz
      current_data[5] = s0z4              + s1z5              ; // zz
      current_data[6] = x6      * scale0  + x7      * scale1  ; // xxx
      current_data[7] = x4 * s0y2         + x5 * s1y3         ; // xxy
      current_data[8] = x2 * s0y4         + x3 * s1y5         ; // xyy
      current_data[9] = y6      * scale0  + y7      * scale1  ; // yyy
      current_data[10] = x4 * s0z2         + x5 * s1z3         ; // xxz
      current_data[11] = x2 * s0y2 * z2    + x3 * s1y3 * z3    ; // xyz
      current_data[12] = s0y4 * z2         + s1y5 * z3         ; // yyz
      current_data[13] = x2 * s0z4         + x3 * s1z5         ; // xzz
      current_data[14] = s0y2 * z4         + s1y3 * z5         ; // yzz
      current_data[15] = z6      * scale0  + z7      * scale1  ; // zzz
    }
  } else if (ang == 96) { // (0, 0, 3, 0)
    double x2, x3, x4, x5, x6, x7;
    double y2, y3, y4, y5, y6, y7;
    double z2, z3, z4, z5, z6, z7;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      {
        const double oxq2 = 0.5 / xq_[ii];
        const double xpopq = xp_[ii] / (xp_[ii] + xq_[ii]);
        const double d00i0x = Q_[ii3] - cx;
        const double d00i0y = Q_[ii3 + 1] - cy;
        const double d00i0z = Q_[ii3 + 2] - cz;
        const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
        const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
        const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
        const double b01i0 = xpopq * oxq2;
        const double b01_0 = oxq2 - b01i0 * roots_[offset];
        const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
        x2 = d00i0x + d00i1x * roots_[offset];
        x3 = d00i0x + d00i1x * roots_[offset + 1];
        x4 = x2 * x2 + b01_0;
        x5 = x3 * x3 + b01_1;
        x6 = x2 * (x4 + 2.0 * b01_0);
        x7 = x3 * (x5 + 2.0 * b01_1);

        y2 = d00i0y + d00i1y * roots_[offset];
        y3 = d00i0y + d00i1y * roots_[offset + 1];
        y4 = y2 * y2 + b01_0;
        y5 = y3 * y3 + b01_1;
        y6 = y2 * (y4 + 2.0 * b01_0);
        y7 = y3 * (y5 + 2.0 * b01_1);

        z2 = d00i0z + d00i1z * roots_[offset];
        z3 = d00i0z + d00i1z * roots_[offset + 1];
        z4 = z2 * z2 + b01_0;
        z5 = z3 * z3 + b01_1;
        z6 = z2 * (z4 + 2.0 * b01_0);
        z7 = z3 * (z5 + 2.0 * b01_1);
      }
      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double s0y2 = y2 * scale0;
      const double s1y3 = y3 * scale1;
      const double s0y4 = y4 * scale0;
      const double s1y5 = y5 * scale1;
      current_data[0] = x6      * scale0  + x7      * scale1  ; // xxx
      current_data[1] = x4 * s0y2         + x5 * s1y3         ; // xxy
      current_data[2] = x2 * s0y4         + x3 * s1y5         ; // xyy
      current_data[3] = y6      * scale0  + y7      * scale1  ; // yyy
      current_data[4] = x4 * z2 * scale0  + x5 * z3 * scale1  ; // xxz
      current_data[5] = x2 * s0y2 * z2    + x3 * s1y3 * z3    ; // xyz
      current_data[6] = s0y4 * z2         + s1y5 * z3         ; // yyz
      current_data[7] = x2 * z4 * scale0  + x3 * z5 * scale1  ; // xzz
      current_data[8] = s0y2 * z4         + s1y3 * z5         ; // yzz
      current_data[9] = z6      * scale0  + z7      * scale1  ; // zzz
    }
  } else if (ang == 65) { // (0, 0, 2, 1)
    double x2, x3, x4, x5, x6, x7;
    double y2, y3, y4, y5, y6, y7;
    double z2, z3, z4, z5, z6, z7;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      {
        const double oxq2 = 0.5 / xq_[ii];
        const double xpopq = xp_[ii] / (xp_[ii] + xq_[ii]);
        const double d00i0x = Q_[ii3] - cx;
        const double d00i0y = Q_[ii3 + 1] - cy;
        const double d00i0z = Q_[ii3 + 2] - cz;
        const double d00i1x = (P_[ii3] - Q_[ii3]) * xpopq;
        const double d00i1y = (P_[ii3 + 1] - Q_[ii3 + 1]) * xpopq;
        const double d00i1z = (P_[ii3 + 2] - Q_[ii3 + 2]) * xpopq;
        const double b01i0 = xpopq * oxq2;
        const double b01_0 = oxq2 - b01i0 * roots_[offset];
        const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
        x2 = d00i0x + d00i1x * roots_[offset];
        x3 = d00i0x + d00i1x * roots_[offset + 1];
        x4 = x2 * x2 + b01_0;
        x5 = x3 * x3 + b01_1;
        x6 = x2 * (x4 + 2.0 * b01_0);
        x7 = x3 * (x5 + 2.0 * b01_1);

        y2 = d00i0y + d00i1y * roots_[offset];
        y3 = d00i0y + d00i1y * roots_[offset + 1];
        y4 = y2 * y2 + b01_0;
        y5 = y3 * y3 + b01_1;
        y6 = y2 * (y4 + 2.0 * b01_0);
        y7 = y3 * (y5 + 2.0 * b01_1);

        z2 = d00i0z + d00i1z * roots_[offset];
        z3 = d00i0z + d00i1z * roots_[offset + 1];
        z4 = z2 * z2 + b01_0;
        z5 = z3 * z3 + b01_1;
        z6 = z2 * (z4 + 2.0 * b01_0);
        z7 = z3 * (z5 + 2.0 * b01_1);
      }
      const double scale0 = coeff_[ii] * weights_[offset];
      const double scale1 = coeff_[ii] * weights_[offset + 1];
      const double s0y2 = y2 * scale0;
      const double s1y3 = y3 * scale1;
      const double s0y4 = y4 * scale0;
      const double s1y5 = y5 * scale1;
      const double s0z2 = z2 * scale0;
      const double s1z3 = z3 * scale1;
      const double s0z4 = z4 * scale0;
      const double s1z5 = z5 * scale1;
      current_data[0] = x4 * scale0       + x5 * scale1       ; // xx
      current_data[1] = x2 * s0y2         + x3 * s1y3         ; // xy
      current_data[2] = s0y4              + s1y5              ; // yy
      current_data[3] = x2 * s0z2         + x3 * s1z3         ; // xz
      current_data[4] = s0y2 * z2         + s1y3 * z3         ; // yz
      current_data[5] = s0z4              + s1z5              ; // zz
      current_data[6] = x6      * scale0  + x7      * scale1  ; // xxx
      current_data[7] = x4 * s0y2         + x5 * s1y3         ; // xxy
      current_data[8] = x2 * s0y4         + x3 * s1y5         ; // xyy
      current_data[9] = y6      * scale0  + y7      * scale1  ; // yyy
      current_data[10] = x4 * s0z2         + x5 * s1z3         ; // xxz
      current_data[11] = x2 * s0y2 * z2    + x3 * s1y3 * z3    ; // xyz
      current_data[12] = s0y4 * z2         + s1y5 * z3         ; // yyz
      current_data[13] = x2 * s0z4         + x3 * s1z5         ; // xzz
      current_data[14] = s0y2 * z4         + s1y3 * z5         ; // yzz
      current_data[15] = z6      * scale0  + z7      * scale1  ; // zzz
    }
  }
}


void ERIBatch::perform_VRR1() {
  const int ang0 = basisinfo_[0]->angular_number();
  const int ang2 = basisinfo_[2]->angular_number();
  if (!ang0 && !ang2) {
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      data_[ii] = coeff_[ii] * weights_[ii];
    }
  } else if (ang0) { // (1, 0, 0, 0)
    const double ax = basisinfo_[0]->position(0);
    const double ay = basisinfo_[0]->position(1);
    const double az = basisinfo_[0]->position(2);
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int ii3 = 3 * ii;

      const double cxp = xp_[ii];
      const double cxq = xq_[ii];
      const double one_pq = 1.0 / (cxp + cxq);
      const double xqopq = cxq * one_pq;
      double x1, y1, z1;
      {
        const double px = P_[ii3];
        const double py = P_[ii3 + 1];
        const double pz = P_[ii3 + 2];
        const double c00i0x = px - ax;
        const double c00i1x = (px - Q_[ii3]) * xqopq;
        const double c00i0y = py - ay;
        const double c00i1y = (py - Q_[ii3 + 1]) * xqopq;
        const double c00i0z = pz - az;
        const double c00i1z = (pz - Q_[ii3 + 2]) * xqopq;
        const double root = roots_[ii];
        x1 = c00i0x - c00i1x * root;
        y1 = c00i0y - c00i1y * root;
        z1 = c00i0z - c00i1z * root;
      }

      const double scale0 = coeff_[ii] * weights_[ii];
      data_[ii3] = scale0 * x1;
      data_[ii3 + 1] = scale0 * y1;
      data_[ii3 + 2] = scale0 * z1;
    }
  } else if (ang2) { // (0, 0, 1, 0)
    const double cx = basisinfo_[2]->position(0);
    const double cy = basisinfo_[2]->position(1);
    const double cz = basisinfo_[2]->position(2);
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];

      const int ii3 = 3 * ii;
      const double cxp = xp_[ii];
      const double cxq = xq_[ii];
      const double one_pq = 1.0 / (cxp + cxq);
      const double xpopq = cxp * one_pq;
      double x1, y1, z1;
      {
        const double qx = Q_[ii3];
        const double qy = Q_[ii3 + 1];
        const double qz = Q_[ii3 + 2];
        const double d00i0x = qx - cx;
        const double d00i1x = (P_[ii3] - qx) * xpopq;
        const double d00i0y = qy - cy;
        const double d00i1y = (P_[ii3 + 1] - qy) * xpopq;
        const double d00i0z = qz - cz;
        const double d00i1z = (P_[ii3 + 2] - qz) * xpopq;
        const double root = roots_[ii];
        x1 = d00i0x + d00i1x * root;
        y1 = d00i0y + d00i1y * root;
        z1 = d00i0z + d00i1z * root;
      }
      const double scale0 = coeff_[ii] * weights_[ii];
      data_[ii3] = scale0 * x1;
      data_[ii3 + 1] = scale0 * y1;
      data_[ii3 + 2] = scale0 * z1;
    }
  }

}

