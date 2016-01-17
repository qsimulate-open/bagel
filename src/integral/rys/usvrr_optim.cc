//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: usvrr_template.cc
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


#include <src/integral/rys/slaterbatch.h>

using namespace std;
using namespace bagel;

#ifdef HAVE_LIBSLATER

void SlaterBatch::perform_USVRR1() {
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    const double cw = weights_[ii] * roots_[ii];
    data2_[ii] = coeffy_[ii] * cw;
    data_[ii]  = coeff_[ii] * (weights_[ii] - cw);

  }
}


void SlaterBatch::perform_USVRR2() {

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
  if (ang == 32768) { // (1, 0, 0, 0)
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int ii3 = 3 * ii;
      const int offset = ii * rank_;

      const double cxp = xp_[ii];
      const double cxq = xq_[ii];
      const double one_pq = 1.0 / (cxp + cxq);
      const double xqopq = cxq * one_pq;
      double x0, y0, z0, x1, y1, z1;
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
        const double root0 = roots_[offset];
        const double root1 = roots_[offset + 1];
        x0 = c00i0x - c00i1x * root0;
        y0 = c00i0y - c00i1y * root0;
        z0 = c00i0z - c00i1z * root0;
        x1 = c00i0x - c00i1x * root1;
        y1 = c00i0y - c00i1y * root1;
        z1 = c00i0z - c00i1z * root1;
      }

      const double wr0 = weights_[offset] * roots_[offset];
      const double wr1 = weights_[offset + 1] * roots_[offset + 1];

      {
        const double scale0 = coeff_[ii] * (weights_[offset] - wr0);
        const double scale1 = coeff_[ii] * (weights_[offset + 1] -wr1);
        data_[ii3    ] = scale0 * x0 + scale1 * x1;
        data_[ii3 + 1] = scale0 * y0 + scale1 * y1;
        data_[ii3 + 2] = scale0 * z0 + scale1 * z1;
      }
      {
        const double scale0 = coeffy_[ii] * wr0;
        const double scale1 = coeffy_[ii] * wr1;
        data2_[ii3    ] = scale0 * x0 + scale1 * x1;
        data2_[ii3 + 1] = scale0 * y0 + scale1 * y1;
        data2_[ii3 + 2] = scale0 * z0 + scale1 * z1;
      }
    }
  } else if (ang == 32) { // (0, 0, 1, 0)
    const double cx = basisinfo_[2]->position(0);
    const double cy = basisinfo_[2]->position(1);
    const double cz = basisinfo_[2]->position(2);
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int ii3 = 3 * ii;
      const double cxp = xp_[ii];
      const double cxq = xq_[ii];
      const double one_pq = 1.0 / (cxp + cxq);
      const double xpopq = cxp * one_pq;
      double x0, y0, z0, x1, y1, z1;
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
        const double root0 = roots_[offset    ];
        const double root1 = roots_[offset + 1];
        x0 = d00i0x + d00i1x * root0;
        y0 = d00i0y + d00i1y * root0;
        z0 = d00i0z + d00i1z * root0;
        x1 = d00i0x + d00i1x * root1;
        y1 = d00i0y + d00i1y * root1;
        z1 = d00i0z + d00i1z * root1;
      }
      const double wr0 = weights_[offset] * roots_[offset];
      const double wr1 = weights_[offset + 1] * roots_[offset + 1];

      {
        const double scale0 = coeff_[ii] * (weights_[offset] - wr0);
        const double scale1 = coeff_[ii] * (weights_[offset + 1] - wr1);
        data_[ii3    ] = scale0 * x0 + scale1 * x1;
        data_[ii3 + 1] = scale0 * y0 + scale1 * y1;
        data_[ii3 + 2] = scale0 * z0 + scale1 * z1;
      }
      {
        const double scale0 = coeffy_[ii] * wr0;
        const double scale1 = coeffy_[ii] * wr1;
        data2_[ii3    ] = scale0 * x0 + scale1 * x1;
        data2_[ii3 + 1] = scale0 * y0 + scale1 * y1;
        data2_[ii3 + 2] = scale0 * z0 + scale1 * z1;
      }
    }
  } else if (ang == 65536) {
    double x2, x3, x4, x5;
    double y2, y3, y4, y5;
    double z2, z3, z4, z5;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

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

      const double wr0 = weights_[offset] * roots_[offset];
      const double wr1 = weights_[offset + 1] * roots_[offset + 1];

      {
        double* current_data = &data_[data_offset_ii];
        const double scale0 = coeff_[ii] * (weights_[offset] - wr0);
        const double scale1 = coeff_[ii] * (weights_[offset + 1] - wr1);
        const double y2s0 = y2 * scale0;
        const double y3s1 = y3 * scale1;
        current_data[0] = scale0 * x4      + scale1 * x5      ; // xx
        current_data[1] = y2s0   * x2      + y3s1   * x3      ; // xy
        current_data[2] = scale0 * y4      + scale1 * y5      ; // yy
        current_data[3] = scale0 * x2 * z2 + scale1 * x3 * z3 ; // xz
        current_data[4] = y2s0   * z2      + y3s1   * z3      ; // yz
        current_data[5] = scale0 * z4      + scale1 * z5      ; // zz
      }
      {
        double* current_data2 = &data2_[data_offset_ii];
        const double scale0 = coeffy_[ii] * wr0;
        const double scale1 = coeffy_[ii] * wr1;
        const double y2s0 = y2 * scale0;
        const double y3s1 = y3 * scale1;
        current_data2[0] = scale0 * x4      + scale1 * x5      ; // xx
        current_data2[1] = y2s0   * x2      + y3s1   * x3      ; // xy
        current_data2[2] = scale0 * y4      + scale1 * y5      ; // yy
        current_data2[3] = scale0 * x2 * z2 + scale1 * x3 * z3 ; // xz
        current_data2[4] = y2s0   * z2      + y3s1   * z3      ; // yz
        current_data2[5] = scale0 * z4      + scale1 * z5      ; // zz
      }
    }
  } else if (ang == 64) {
    double x2, x3, x4, x5;
    double y2, y3, y4, y5;
    double z2, z3, z4, z5;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;


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
      const double wr0 = weights_[offset] * roots_[offset];
      const double wr1 = weights_[offset + 1] * roots_[offset + 1];

      {
        double* current_data = &data_[data_offset_ii];
        const double scale0 = coeff_[ii] * (weights_[offset    ] - wr0);
        const double scale1 = coeff_[ii] * (weights_[offset + 1] - wr1);
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
      {
        double* current_data = &data2_[data_offset_ii];
        const double scale0 = coeffy_[ii] * wr0;
        const double scale1 = coeffy_[ii] * wr1;
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
    }
  } else if (ang == 33792) { // (1, 1, 0, 0)
    double x2, x3, x4, x5;
    double y2, y3, y4, y5;
    double z2, z3, z4, z5;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

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

      const double wr0 = weights_[offset] * roots_[offset];
      const double wr1 = weights_[offset + 1] * roots_[offset + 1];

      {
        double* current_data = &data_[data_offset_ii];
        const double scale0 = coeff_[ii] * (weights_[offset] - wr0);
        const double scale1 = coeff_[ii] * (weights_[offset + 1] - wr1);
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
      {
        double* current_data = &data2_[data_offset_ii];
        const double scale0 = coeffy_[ii] * wr0;
        const double scale1 = coeffy_[ii] * wr1;
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
    }
  } else if (ang == 33) { // (0, 0, 1, 1)
    double x2, x3, x4, x5;
    double y2, y3, y4, y5;
    double z2, z3, z4, z5;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

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

      const double wr0 = weights_[offset] * roots_[offset];
      const double wr1 = weights_[offset + 1] * roots_[offset + 1];

      {
        double* current_data = &data_[data_offset_ii];
        const double scale0 = coeff_[ii] * (weights_[offset] - wr0);
        const double scale1 = coeff_[ii] * (weights_[offset + 1] - wr1);
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
      {
        double* current_data = &data2_[data_offset_ii];
        const double scale0 = coeffy_[ii] * wr0;
        const double scale1 = coeffy_[ii] * wr1;
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
    }
  } else if (ang == 32800) { // (1, 0, 1, 0)
    double x2, x3, x4, x5, x6, x7;
    double y2, y3, y4, y5, y6, y7;
    double z2, z3, z4, z5, z6, z7;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

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

      const double wr0 = weights_[offset] * roots_[offset];
      const double wr1 = weights_[offset + 1] * roots_[offset + 1];

      {
        double* current_data = &data_[data_offset_ii];
        const double scale0 = coeff_[ii] * (weights_[offset] - wr0);
        const double scale1 = coeff_[ii] * (weights_[offset + 1] - wr1);
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
      {
        double* current_data = &data2_[data_offset_ii];
        const double scale0 = coeffy_[ii] * wr0;
        const double scale1 = coeffy_[ii] * wr1;
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
    }
  }
}

#endif
