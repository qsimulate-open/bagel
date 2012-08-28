//
// BAGEL - Parallel electron correlation program.
// Filename: svrr_template.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/slater/slaterbatch.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <cstring>

//#define LOCAL_DEBUG2
//#define LOCAL_DEBUG3

using namespace std;
using namespace bagel;

void SlaterBatch::perform_SVRR1() {
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    data_[ii] = coeff_[ii] * weights_[ii] * (1.0 - roots_[ii]); 
  }
}


void SlaterBatch::perform_SVRR2() {

  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);

  const unsigned int ang0 = basisinfo_[0]->angular_number();
  const unsigned int ang1 = basisinfo_[1]->angular_number();
  const unsigned int ang2 = basisinfo_[2]->angular_number();
  const unsigned int ang3 = basisinfo_[3]->angular_number();
  const unsigned int ang = (((((ang0 << 5) + ang1) << 5) + ang2) << 5) + ang3; // offsets are 2^5 = 32
#ifndef LOCAL_DEBUG2
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
        const double px = p_[ii3];
        const double py = p_[ii3 + 1];
        const double pz = p_[ii3 + 2];
        const double c00i0x = px - ax;
        const double c00i1x = (px - q_[ii3]) * xqopq; 
        const double c00i0y = py - ay;
        const double c00i1y = (py - q_[ii3 + 1]) * xqopq; 
        const double c00i0z = pz - az;
        const double c00i1z = (pz - q_[ii3 + 2]) * xqopq; 
        const double root0 = roots_[offset];
        const double root1 = roots_[offset + 1];
        x0 = c00i0x - c00i1x * root0; 
        y0 = c00i0y - c00i1y * root0;
        z0 = c00i0z - c00i1z * root0; 
        x1 = c00i0x - c00i1x * root1; 
        y1 = c00i0y - c00i1y * root1;
        z1 = c00i0z - c00i1z * root1; 
      }
      
      const double scale0 = coeff_[ii] * weights_[offset] * (1.0 - roots_[offset]);
      const double scale1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
      data_[ii3    ] = scale0 * x0 + scale1 * x1;
      data_[ii3 + 1] = scale0 * y0 + scale1 * y1;
      data_[ii3 + 2] = scale0 * z0 + scale1 * z1;
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
        const double qx = q_[ii3];
        const double qy = q_[ii3 + 1];
        const double qz = q_[ii3 + 2];
        const double d00i0x = qx - cx; 
        const double d00i1x = (p_[ii3] - qx) * xpopq; 
        const double d00i0y = qy - cy; 
        const double d00i1y = (p_[ii3 + 1] - qy) * xpopq; 
        const double d00i0z = qz - cz; 
        const double d00i1z = (p_[ii3 + 2] - qz) * xpopq; 
        const double root0 = roots_[offset    ];
        const double root1 = roots_[offset + 1];
        x0 = d00i0x + d00i1x * root0;
        y0 = d00i0y + d00i1y * root0;
        z0 = d00i0z + d00i1z * root0;
        x1 = d00i0x + d00i1x * root1;
        y1 = d00i0y + d00i1y * root1;
        z1 = d00i0z + d00i1z * root1;
      }
      const double scale0 = coeff_[ii] * weights_[offset] * (1.0 - roots_[offset]);
      const double scale1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
      data_[ii3    ] = scale0 * x0 + scale1 * x1;
      data_[ii3 + 1] = scale0 * y0 + scale1 * y1;
      data_[ii3 + 2] = scale0 * z0 + scale1 * z1;
    }
  } else if (ang == 65536) {
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
        const double c00i0x = p_[ii3]     - ax;
        const double c00i0y = p_[ii3 + 1] - ay;
        const double c00i0z = p_[ii3 + 2] - az;
        const double c00i1x = (p_[ii3]     - q_[ii3])     * xqopq; 
        const double c00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xqopq; 
        const double c00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xqopq; 
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

      const double scale0 = coeff_[ii] * weights_[offset] * (1.0 - roots_[offset]);
      const double scale1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
      const double y2s0 = y2 * scale0;
      const double y3s1 = y3 * scale1;
      current_data[0] = scale0 * x4      + scale1 * x5      ; // xx 
      current_data[1] = y2s0   * x2      + y3s1   * x3      ; // xy 
      current_data[2] = scale0 * y4      + scale1 * y5      ; // yy 
      current_data[3] = scale0 * x2 * z2 + scale1 * x3 * z3 ; // xz
      current_data[4] = y2s0   * z2      + y3s1   * z3      ; // yz 
      current_data[5] = scale0 * z4      + scale1 * z5      ; // zz 
    }
  } else if (ang == 64) {
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
        const double d00i0x = q_[ii3] - cx;
        const double d00i0y = q_[ii3 + 1] - cy;
        const double d00i0z = q_[ii3 + 2] - cz;
        const double d00i1x = (p_[ii3] - q_[ii3]) * xpopq; 
        const double d00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xpopq; 
        const double d00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xpopq; 
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
      const double scale0 = coeff_[ii] * weights_[offset    ] * (1.0 - roots_[offset]);
      const double scale1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
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
        const double c00i0x = p_[ii3]     - ax;
        const double c00i0y = p_[ii3 + 1] - ay;
        const double c00i0z = p_[ii3 + 2] - az;
        const double c00i1x = (p_[ii3]     - q_[ii3])     * xqopq; 
        const double c00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xqopq; 
        const double c00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xqopq; 
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

      const double scale0 = coeff_[ii] * weights_[offset] * (1.0 - roots_[offset]);
      const double scale1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
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
        const double d00i0x = q_[ii3] - cx;
        const double d00i0y = q_[ii3 + 1] - cy;
        const double d00i0z = q_[ii3 + 2] - cz;
        const double d00i1x = (p_[ii3] - q_[ii3]) * xpopq; 
        const double d00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xpopq; 
        const double d00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xpopq; 
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
      const double scale0 = coeff_[ii] * weights_[offset] * (1.0 - roots_[offset]);
      const double scale1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
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
        const double oxp2 = 0.5 / xp_[ii];
        const double opq = 1.0 / (xp_[ii] + xq_[ii]);
        const double xqopq = xq_[ii] * opq;
        const double c00i0x = p_[ii3]     - ax;
        const double c00i0y = p_[ii3 + 1] - ay;
        const double c00i0z = p_[ii3 + 2] - az;
        const double c00i1x = (p_[ii3]     - q_[ii3])     * xqopq; 
        const double c00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xqopq; 
        const double c00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xqopq; 
  
        const double oxq2 = 0.5 / xq_[ii];
        const double xpopq = xp_[ii] * opq;
        const double d00i0x = q_[ii3] - cx;
        const double d00i0y = q_[ii3 + 1] - cy;
        const double d00i0z = q_[ii3 + 2] - cz;
        const double d00i1x = (p_[ii3] - q_[ii3]) * xpopq; 
        const double d00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xpopq; 
        const double d00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xpopq; 

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
 
      const double scale0 = coeff_[ii] * weights_[offset] * (1.0 - roots_[offset]);
      const double scale1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
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
#else
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 2 * isize;
  const int svrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = new double[worksize];
  double* worky = new double[worksize];
  double* workz = new double[worksize];
  double iyiz[2];

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D cix(dparamx, &roots_[offset], 2, worksize, workx, vrr_->vrrfunc[svrr_index]);
    double womt[2];
    double wt[2];
    wt[0] = weights_[offset + 0] * roots_[offset + 0];
    wt[1] = weights_[offset + 1] * roots_[offset + 1];
    womt[0] = weights_[offset + 0]  - wt[0];
    womt[1] = weights_[offset + 1]  - wt[1];
    cix.scale_data(womt, coeff_[ii]);
 
    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], 2, worksize, worky, vrr_->vrrfunc[svrr_index]);
 
    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], 2, worksize, workz, vrr_->vrrfunc[svrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
            iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
              }
            }
          }
        }
      }
    }

  }

  delete[] workz;
  delete[] worky;
  delete[] workx;
#endif
}


void SlaterBatch::perform_SVRR3() {

#ifndef LOCAL_DEBUG3
  const int acsize = asize_ * csize_;
  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  const unsigned int ang0 = basisinfo_[0]->angular_number();
  const unsigned int ang1 = basisinfo_[1]->angular_number();
  const unsigned int ang2 = basisinfo_[2]->angular_number();
  const unsigned int ang3 = basisinfo_[3]->angular_number();
  const unsigned int ang = (((((ang0 << 5) + ang1) << 5) + ang2) << 5) + ang3; // offsets are 2^5 = 32
  if (ang == 66560) { // (2, 1, 0, 0)
    double x3, x4, x5, x6, x7, x8, x9, x10, x11;
    double y3, y4, y5, y6, y7, y8, y9, y10, y11;
    double z3, z4, z5, z6, z7, z8, z9, z10, z11;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;
   
      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      {
        const double oxp2 = 0.5 / xp_[ii];
        const double xqopq = xq_[ii] / (xp_[ii] + xq_[ii]);
        const double c00i0x = p_[ii3]     - ax;
        const double c00i0y = p_[ii3 + 1] - ay;
        const double c00i0z = p_[ii3 + 2] - az;
        const double c00i1x = (p_[ii3]     - q_[ii3])     * xqopq; 
        const double c00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xqopq; 
        const double c00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xqopq; 
        const double b10i0 = xqopq * oxp2; 
        const double b10_0 = oxp2 - b10i0 * roots_[offset];
        const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
        const double b10_2 = oxp2 - b10i0 * roots_[offset + 2];
        x3 = c00i0x - c00i1x * roots_[offset]; 
        x4 = c00i0x - c00i1x * roots_[offset + 1]; 
        x5 = c00i0x - c00i1x * roots_[offset + 2]; 
        x6 = x3 * x3 + b10_0;
        x7 = x4 * x4 + b10_1;
        x8 = x5 * x5 + b10_2;
        x9  = x3 * (x6 + 2.0 * b10_0);
        x10 = x4 * (x7 + 2.0 * b10_1);
        x11 = x5 * (x8 + 2.0 * b10_2);

        y3 = c00i0y - c00i1y * roots_[offset]; 
        y4 = c00i0y - c00i1y * roots_[offset + 1]; 
        y5 = c00i0y - c00i1y * roots_[offset + 2]; 
        y6 = y3 * y3 + b10_0;
        y7 = y4 * y4 + b10_1;
        y8 = y5 * y5 + b10_2;
        y9  = y3 * (y6 + 2.0 * b10_0);
        y10 = y4 * (y7 + 2.0 * b10_1);
        y11 = y5 * (y8 + 2.0 * b10_2);

        z3 = c00i0z - c00i1z * roots_[offset]; 
        z4 = c00i0z - c00i1z * roots_[offset + 1]; 
        z5 = c00i0z - c00i1z * roots_[offset + 2]; 
        z6 = z3 * z3 + b10_0;
        z7 = z4 * z4 + b10_1;
        z8 = z5 * z5+ b10_2;
        z9  = z3 * (z6 + 2.0 * b10_0);
        z10 = z4 * (z7 + 2.0 * b10_1);
        z11 = z5 * (z8 + 2.0 * b10_2);
      }

      const double scale0 = coeff_[ii] * weights_[offset    ] * (1.0 - roots_[offset]);
      const double scale1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
      const double scale2 = coeff_[ii] * weights_[offset + 2] * (1.0 - roots_[offset + 2]);
      const double s0z3 = z3 * scale0;
      const double s1z4 = z4 * scale1;
      const double s2z5 = z5 * scale2;
      const double s0z6 = z6 * scale0;
      const double s1z7 = z7 * scale1;
      const double s2z8 = z8 * scale2;
      const double s0y3 = y3 * scale0;
      const double s1y4 = y4 * scale1;
      const double s2y5 = y5 * scale2;
      const double s0y6 = y6 * scale0;
      const double s1y7 = y7 * scale1;
      const double s2y8 = y8 * scale2;
      current_data[0] = x6 * scale0       + x7 * scale1       + x8 * scale2; // xx 
      current_data[1] = x3 * s0y3         + x4 * s1y4         + x5 * s2y5; // xy 
      current_data[2] = s0y6              + s1y7              + s2y8; // yy 
      current_data[3] = x3 * s0z3         + x4 * s1z4         + x5 * s2z5; // xz 
      current_data[4] = s0y3 * z3         + s1y4 * z4         + s2y5 * z5; // yz 
      current_data[5] = s0z6              + s1z7              + s2z8; // zz 
      current_data[6] = x9      * scale0  + x10     * scale1  + x11 * scale2; // x,x,x
      current_data[7] = x6 * s0y3         + x7 * s1y4         + x8 * s2y5; // xxy
      current_data[8] = x3 * s0y6         + x4 * s1y7         + x5 * s2y8; // xyy
      current_data[9] = y9      * scale0  + y10      * scale1  + y11 * scale2; // yyy 
      current_data[10] = x6 * s0z3        + x7 * s1z4          + x8 * s2z5; // xxz
      current_data[11] = x3 * s0y3 * z3    + x4 * s1y4 * z4    + x5 * s2y5 * z5; // xyz
      current_data[12] = s0y6 * z3         + s1y7 * z4         + s2y8 * z5; // yyz
      current_data[13] = x3 * s0z6         + x4 * s1z7         + x5 * s2z8; // xzz
      current_data[14] = s0y3 * z6         + s1y4 * z7         + s2y5 * z8; // yzz
      current_data[15] = z9      * scale0  + z10      * scale1  + z11 * scale2; // zzz 
    }
  } else if (ang == 65) { // (0, 0, 2, 1)
    double x3, x4, x5, x6, x7, x8, x9, x10, x11;
    double y3, y4, y5, y6, y7, y8, y9, y10, y11;
    double z3, z4, z5, z6, z7, z8, z9, z10, z11;
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      {
        const double oxq2 = 0.5 / xq_[ii];
        const double xpopq = xp_[ii] / (xp_[ii] + xq_[ii]);
        const double d00i0x = q_[ii3] - cx;
        const double d00i0y = q_[ii3 + 1] - cy;
        const double d00i0z = q_[ii3 + 2] - cz;
        const double d00i1x = (p_[ii3] - q_[ii3]) * xpopq; 
        const double d00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xpopq; 
        const double d00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xpopq; 
        const double b01i0 = xpopq * oxq2; 
        const double b01_0 = oxq2 - b01i0 * roots_[offset];
        const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
        const double b01_2 = oxq2 - b01i0 * roots_[offset + 2];
        x3 = d00i0x + d00i1x * roots_[offset]; 
        x4 = d00i0x + d00i1x * roots_[offset + 1]; 
        x5 = d00i0x + d00i1x * roots_[offset + 2]; 
        x6 = x3 * x3 + b01_0;
        x7 = x4 * x4 + b01_1;
        x8 = x5 * x5 + b01_2;
        x9 = x3 * (x6 + 2.0 * b01_0);
        x10 = x4 * (x7 + 2.0 * b01_1);
        x11 = x5 * (x8 + 2.0 * b01_2);

        y3 = d00i0y + d00i1y * roots_[offset]; 
        y4 = d00i0y + d00i1y * roots_[offset + 1]; 
        y5 = d00i0y + d00i1y * roots_[offset + 2]; 
        y6 = y3 * y3 + b01_0;
        y7 = y4 * y4 + b01_1;
        y8 = y5 * y5 + b01_2;
        y9 = y3 * (y6 + 2.0 * b01_0);
        y10 = y4 * (y7 + 2.0 * b01_1);
        y11 = y5 * (y8 + 2.0 * b01_2);

        z3 = d00i0z + d00i1z * roots_[offset]; 
        z4 = d00i0z + d00i1z * roots_[offset + 1]; 
        z5 = d00i0z + d00i1z * roots_[offset + 2]; 
        z6 = z3 * z3 + b01_0;
        z7 = z4 * z4 + b01_1;
        z8 = z5 * z5 + b01_2;
        z9 = z3 * (z6 + 2.0 * b01_0);
        z10 = z4 * (z7 + 2.0 * b01_1);
        z11 = z5 * (z8 + 2.0 * b01_2);
      }
      const double scale0 = coeff_[ii] * weights_[offset    ] * (1.0 - roots_[offset]);
      const double scale1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
      const double scale2 = coeff_[ii] * weights_[offset + 2] * (1.0 - roots_[offset + 2]);
      const double s0y3 = y3 * scale0;
      const double s1y4 = y4 * scale1;
      const double s2y5 = y5 * scale2;
      const double s0y6 = y6 * scale0;
      const double s1y7 = y7 * scale1;
      const double s2y8 = y8 * scale2;
      const double s0z3 = z3 * scale0;
      const double s1z4 = z4 * scale1;
      const double s2z5 = z5 * scale2;
      const double s0z6 = z6 * scale0;
      const double s1z7 = z7 * scale1;
      const double s2z8 = z8 * scale2;
      current_data[0] = x6 * scale0       + x7 * scale1       + x8 * scale2; // xx 
      current_data[1] = x3 * s0y3         + x4 * s1y4         + x5 * s2y5; // xy 
      current_data[2] = s0y6              + s1y7              + s2y8; // yy 
      current_data[3] = x3 * s0z3         + x4 * s1z4         + x5 * s2z5; // xz 
      current_data[4] = s0y3 * z3         + s1y4 * z4         + s2y5 * z5; // yz 
      current_data[5] = s0z6              + s1z7              + s2z8; // zz 
      current_data[6] = x9      * scale0  + x10      * scale1 + x11 * scale2; // x,x,x
      current_data[7] = x6 * s0y3         + x7 * s1y4         + x8 * s2y5; // xxy
      current_data[8] = x3 * s0y6         + x4 * s1y7         + x5 * s2y8; // xyy
      current_data[9] = y9      * scale0  + y10      * scale1  + y11 * scale2; // yyy 
      current_data[10] = x6 * s0z3         + x7 * s1z4         + x8 * s2z5; // xxz
      current_data[11] = x3 * s0y3 * z3    + x4 * s1y4 * z4    + x5 * s2y5 * z5; // xyz
      current_data[12] = s0y6 * z3         + s1y7 * z4         + s2y8 * z5; // yyz
      current_data[13] = x3 * s0z6         + x4 * s1z7         + x5 * s2z8; // xzz
      current_data[14] = s0y3 * z6         + s1y4 * z7         + s2y5 * z8; // yzz
      current_data[15] = z9      * scale0  + z10      * scale1 + z11 * scale2; // zzz 
    }
  } else if (ang == 65568) { // (2, 0, 1, 0)
    double x3, x4, x6, x7, x9, x10, x12, x13, x15, x16;
    double y3, y4, y6, y7, y9, y10, y12, y13, y15, y16;
    double z3, z4, z6, z7, z9, z10, z12, z13, z15, z16;
    double x5, x8, x11, x14, x17;
    double y5, y8, y11, y14, y17;
    double z5, z8, z11, z14, z17;
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
        const double c00i0x = p_[ii3]     - ax;
        const double c00i0y = p_[ii3 + 1] - ay;
        const double c00i0z = p_[ii3 + 2] - az;
        const double c00i1x = (p_[ii3]     - q_[ii3])     * xqopq; 
        const double c00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xqopq; 
        const double c00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xqopq; 
        const double b10i0 = xqopq * oxp2; 
        const double b10_0 = oxp2 - b10i0 * roots_[offset];
        const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
        const double b10_2 = oxp2 - b10i0 * roots_[offset + 2];
        const double oxq2 = 0.5 / xq_[ii];
        const double xpopq = xp_[ii] * opq;
        const double d00i0x = q_[ii3] - cx;
        const double d00i0y = q_[ii3 + 1] - cy;
        const double d00i0z = q_[ii3 + 2] - cz;
        const double d00i1x = (p_[ii3] - q_[ii3]) * xpopq; 
        const double d00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xpopq; 
        const double d00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xpopq; 
        const double b00int = 0.5 * opq;
        const double b00_0 = b00int * roots_[offset];
        const double b00_1 = b00int * roots_[offset + 1];
        const double b00_2 = b00int * roots_[offset + 2];
        x3 = c00i0x - c00i1x * roots_[offset]; 
        x4 = c00i0x - c00i1x * roots_[offset + 1]; 
        x5 = c00i0x - c00i1x * roots_[offset + 2]; 
        x6 = x3 * x3 + b10_0;
        x7 = x4 * x4 + b10_1;
        x8 = x5 * x5 + b10_2;
        x9  = d00i0x + d00i1x * roots_[offset]; 
        x10 = d00i0x + d00i1x * roots_[offset + 1]; 
        x11 = d00i0x + d00i1x * roots_[offset + 2]; 
        x12 = x3 * x9  + b00_0;
        x13 = x4 * x10 + b00_1;
        x14 = x5 * x11 + b00_2;
        x15 = x3 * (x12 + b00_0) + b10_0 * x9;
        x16 = x4 * (x13 + b00_1) + b10_1 * x10;
        x17 = x5 * (x14 + b00_2) + b10_2 * x11;

        y3 = c00i0y - c00i1y * roots_[offset]; 
        y4 = c00i0y - c00i1y * roots_[offset + 1]; 
        y5 = c00i0y - c00i1y * roots_[offset + 2]; 
        y6 = y3 * y3 + b10_0;
        y7 = y4 * y4 + b10_1;
        y8 = y5 * y5 + b10_2;
        y9  = d00i0y + d00i1y * roots_[offset]; 
        y10 = d00i0y + d00i1y * roots_[offset + 1]; 
        y11 = d00i0y + d00i1y * roots_[offset + 2]; 
        y12 = y3 * y9  + b00_0;
        y13 = y4 * y10 + b00_1;
        y14 = y5 * y11 + b00_2;
        y15 = y3 * (y12 + b00_0) + b10_0 * y9;
        y16 = y4 * (y13 + b00_1) + b10_1 * y10;
        y17 = y5 * (y14 + b00_2) + b10_2 * y11;

        z3 = c00i0z - c00i1z * roots_[offset]; 
        z4 = c00i0z - c00i1z * roots_[offset + 1]; 
        z5 = c00i0z - c00i1z * roots_[offset + 2]; 
        z6 = z3 * z3 + b10_0;
        z7 = z4 * z4 + b10_1;
        z8 = z5 * z5 + b10_2;
        z9  = d00i0z + d00i1z * roots_[offset]; 
        z10 = d00i0z + d00i1z * roots_[offset + 1]; 
        z11 = d00i0z + d00i1z * roots_[offset + 2]; 
        z12 = z3 * z9  + b00_0;
        z13 = z4 * z10 + b00_1;
        z14 = z5 * z11 + b00_2;
        z15 = z3 * (z12 + b00_0) + b10_0 * z9;
        z16 = z4 * (z13 + b00_1) + b10_1 * z10;
        z17 = z5 * (z14 + b00_2) + b10_2 * z11;

        const double s0 = coeff_[ii] * weights_[offset    ] * (1.0 - roots_[offset]);
        const double s1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
        const double s2 = coeff_[ii] * weights_[offset + 2] * (1.0 - roots_[offset + 2]);

        const double s0x9 = s0 * x9;
        const double s1x10 = s1 * x10;
        const double s2x11 = s2 * x11;
        const double s0x12 = s0 * x12;
        const double s1x13 = s1 * x13;
        const double s2x14 = s2 * x14;
        const double s0x15 = s0 * x15;
        const double s1x16 = s1 * x16;
        const double s2x17 = s2 * x17;
        current_data[0] = s0x15            + s1x16             + s2x17; // x|xx 
        current_data[1] = s0x12 * y3       + s1x13 * y4        + s2x14 * y5; // x|xy 
        current_data[2] = s0x9 * y6        + s1x10 * y7        + s2x11 * y8; // x|yy 
        current_data[3] = s0x12 * z3       + s1x13 * z4        + s2x14 * z5; // x|xz 
        current_data[4] = s0x9 * y3 * z3   + s1x10 * y4 * z4   + s2x11 * y5 * z5; // x|yz 
        current_data[5] = s0x9 * z6        + s1x10 * z7        + s2x11 * z8; // x|zz 

        const double s0y9 = s0 * y9;
        const double s1y10 = s1 * y10;
        const double s2y11 = s2 * y11;
        const double s0y12 = s0 * y12;
        const double s1y13 = s1 * y13;
        const double s2y14 = s2 * y14;
        const double s0y15 = s0 * y15;
        const double s1y16 = s1 * y16;
        const double s2y17 = s2 * y17;
        current_data[6] =  s0y9 * x6        + s1y10 * x7        + s2y11 * x8; // y|xx
        current_data[7] =  s0y12 * x3       + s1y13 * x4        + s2y14 * x5; // y|xy 
        current_data[8] =  s0y15            + s1y16             + s2y17; // y|yy 
        current_data[9] =  s0y9 * x3 * z3   + s1y10 * x4 * z4   + s2y11 * x5 * z5; // y|xz 
        current_data[10] = s0y12 * z3       + s1y13 * z4        + s2y14 * z5; // y|yz 
        current_data[11] = s0y9 * z6        + s1y10 * z7        + s2y11 * z8; // y|zz 

        const double s0z9 = s0 * z9;
        const double s1z10 = s1 * z10;
        const double s2z11 = s2 * z11;
        const double s0z12 = s0 * z12;
        const double s1z13 = s1 * z13;
        const double s2z14 = s2 * z14;
        const double s0z15 = s0 * z15;
        const double s1z16 = s1 * z16;
        const double s2z17 = s2 * z17;
        current_data[12] = s0z9 * x6        + s1z10 * x7        + s2z11 * x8; // z|xx
        current_data[13] = s0z9 * x3 * y3   + s1z10 * x4 * y4   + s2z11 * x5 * y5; // z|xy 
        current_data[14] = s0z9 * y6        + s1z10 * y7        + s2z11 * y8; // z|yy 
        current_data[15] = s0z12 * x3       + s1z13 * x4        + s2z14 * x5; // z|xz 
        current_data[16] = s0z12 * y3       + s1z13 * y4        + s2z14 * y5; // z|yz 
        current_data[17] = s0z15            + s1z16             + s2z17; // z|zz 
      }
    } 
  } else if (ang == 33824) { // (1, 1, 1, 0)
    double x3, x4, x6, x7, x9, x10, x12, x13, x15, x16;
    double y3, y4, y6, y7, y9, y10, y12, y13, y15, y16;
    double z3, z4, z6, z7, z9, z10, z12, z13, z15, z16;
    double x5, x8, x11, x14, x17;
    double y5, y8, y11, y14, y17;
    double z5, z8, z11, z14, z17;
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
        const double c00i0x = p_[ii3]     - ax;
        const double c00i0y = p_[ii3 + 1] - ay;
        const double c00i0z = p_[ii3 + 2] - az;
        const double c00i1x = (p_[ii3]     - q_[ii3])     * xqopq; 
        const double c00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xqopq; 
        const double c00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xqopq; 
        const double b10i0 = xqopq * oxp2; 
        const double b10_0 = oxp2 - b10i0 * roots_[offset];
        const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
        const double b10_2 = oxp2 - b10i0 * roots_[offset + 2];
        const double oxq2 = 0.5 / xq_[ii];
        const double xpopq = xp_[ii] * opq;
        const double d00i0x = q_[ii3] - cx;
        const double d00i0y = q_[ii3 + 1] - cy;
        const double d00i0z = q_[ii3 + 2] - cz;
        const double d00i1x = (p_[ii3] - q_[ii3]) * xpopq; 
        const double d00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xpopq; 
        const double d00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xpopq; 
        const double b00int = 0.5 * opq;
        const double b00_0 = b00int * roots_[offset];
        const double b00_1 = b00int * roots_[offset + 1];
        const double b00_2 = b00int * roots_[offset + 2];
        x3 = c00i0x - c00i1x * roots_[offset]; 
        x4 = c00i0x - c00i1x * roots_[offset + 1]; 
        x5 = c00i0x - c00i1x * roots_[offset + 2]; 
        x6 = x3 * x3 + b10_0;
        x7 = x4 * x4 + b10_1;
        x8 = x5 * x5 + b10_2;
        x9  = d00i0x + d00i1x * roots_[offset]; 
        x10 = d00i0x + d00i1x * roots_[offset + 1]; 
        x11 = d00i0x + d00i1x * roots_[offset + 2]; 
        x12 = x3 * x9  + b00_0;
        x13 = x4 * x10 + b00_1;
        x14 = x5 * x11 + b00_2;
        x15 = x3 * (x12 + b00_0) + b10_0 * x9;
        x16 = x4 * (x13 + b00_1) + b10_1 * x10;
        x17 = x5 * (x14 + b00_2) + b10_2 * x11;

        y3 = c00i0y - c00i1y * roots_[offset]; 
        y4 = c00i0y - c00i1y * roots_[offset + 1]; 
        y5 = c00i0y - c00i1y * roots_[offset + 2]; 
        y6 = y3 * y3 + b10_0;
        y7 = y4 * y4 + b10_1;
        y8 = y5 * y5 + b10_2;
        y9  = d00i0y + d00i1y * roots_[offset]; 
        y10 = d00i0y + d00i1y * roots_[offset + 1]; 
        y11 = d00i0y + d00i1y * roots_[offset + 2]; 
        y12 = y3 * y9  + b00_0;
        y13 = y4 * y10 + b00_1;
        y14 = y5 * y11 + b00_2;
        y15 = y3 * (y12 + b00_0) + b10_0 * y9;
        y16 = y4 * (y13 + b00_1) + b10_1 * y10;
        y17 = y5 * (y14 + b00_2) + b10_2 * y11;

        z3 = c00i0z - c00i1z * roots_[offset]; 
        z4 = c00i0z - c00i1z * roots_[offset + 1]; 
        z5 = c00i0z - c00i1z * roots_[offset + 2]; 
        z6 = z3 * z3 + b10_0;
        z7 = z4 * z4 + b10_1;
        z8 = z5 * z5 + b10_2;
        z9  = d00i0z + d00i1z * roots_[offset]; 
        z10 = d00i0z + d00i1z * roots_[offset + 1]; 
        z11 = d00i0z + d00i1z * roots_[offset + 2]; 
        z12 = z3 * z9  + b00_0;
        z13 = z4 * z10 + b00_1;
        z14 = z5 * z11 + b00_2;
        z15 = z3 * (z12 + b00_0) + b10_0 * z9;
        z16 = z4 * (z13 + b00_1) + b10_1 * z10;
        z17 = z5 * (z14 + b00_2) + b10_2 * z11;

        const double s0 = coeff_[ii] * weights_[offset    ] * (1.0 - roots_[offset]);
        const double s1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
        const double s2 = coeff_[ii] * weights_[offset + 2] * (1.0 - roots_[offset + 2]);

        const double s0x9 = s0 * x9;
        const double s1x10 = s1 * x10;
        const double s2x11 = s2 * x11;
        const double s0x12 = s0 * x12;
        const double s1x13 = s1 * x13;
        const double s2x14 = s2 * x14;
        const double s0x15 = s0 * x15;
        const double s1x16 = s1 * x16;
        const double s2x17 = s2 * x17;
        current_data[0] = s0x12            + s1x13             + s2x14; // x|x
        current_data[1] = s0x9 * y3        + s1x10 * y4        + s2x11 * y5; // x|y
        current_data[2] = s0x9 * z3        + s1x10 * z4        + s2x11 * z5; // x|z
        current_data[3] = s0x15            + s1x16             + s2x17; // x|xx 
        current_data[4] = s0x12 * y3       + s1x13 * y4        + s2x14 * y5; // x|xy 
        current_data[5] = s0x9 * y6        + s1x10 * y7        + s2x11 * y8; // x|yy 
        current_data[6] = s0x12 * z3       + s1x13 * z4        + s2x14 * z5; // x|xz 
        current_data[7] = s0x9 * y3 * z3   + s1x10 * y4 * z4   + s2x11 * y5 * z5; // x|yz 
        current_data[8] = s0x9 * z6        + s1x10 * z7        + s2x11 * z8; // x|zz 

        const double s0y9 = s0 * y9;
        const double s1y10 = s1 * y10;
        const double s2y11 = s2 * y11;
        const double s0y12 = s0 * y12;
        const double s1y13 = s1 * y13;
        const double s2y14 = s2 * y14;
        const double s0y15 = s0 * y15;
        const double s1y16 = s1 * y16;
        const double s2y17 = s2 * y17;
        current_data[ 9] = s0y9 * x3        + s1y10 * x4        + s2y11 * x5; // y|x
        current_data[10] = s0y12            + s1y13             + s2y14; // y|y
        current_data[11] = s0y9 * z3        + s1y10 * z4        + s2y11 * z5; // x|z
        current_data[12] = s0y9 * x6        + s1y10 * x7        + s2y11 * x8; // y|xx
        current_data[13] = s0y12 * x3       + s1y13 * x4        + s2y14 * x5; // y|xy 
        current_data[14] = s0y15            + s1y16             + s2y17; // y|yy 
        current_data[15] = s0y9 * x3 * z3   + s1y10 * x4 * z4   + s2y11 * x5 * z5; // y|xz 
        current_data[16] = s0y12 * z3       + s1y13 * z4        + s2y14 * z5; // y|yz 
        current_data[17] = s0y9 * z6        + s1y10 * z7        + s2y11 * z8; // y|zz 

        const double s0z9 = s0 * z9;
        const double s1z10 = s1 * z10;
        const double s2z11 = s2 * z11;
        const double s0z12 = s0 * z12;
        const double s1z13 = s1 * z13;
        const double s2z14 = s2 * z14;
        const double s0z15 = s0 * z15;
        const double s1z16 = s1 * z16;
        const double s2z17 = s2 * z17;
        current_data[18] = s0z9 * x3        + s1z10 * x4        + s2z11 * x5; // z|x
        current_data[19] = s0z9 * y3        + s1z10 * y4        + s2z11 * y5; // z|y
        current_data[20] = s0z12            + s1z13             + s2z14; // z|z
        current_data[21] = s0z9 * x6        + s1z10 * x7        + s2z11 * x8; // z|xx
        current_data[22] = s0z9 * x3 * y3   + s1z10 * x4 * y4   + s2z11 * x5 * y5; // z|xy 
        current_data[23] = s0z9 * y6        + s1z10 * y7        + s2z11 * y8; // z|yy 
        current_data[24] = s0z12 * x3       + s1z13 * x4        + s2z14 * x5; // z|xz 
        current_data[25] = s0z12 * y3       + s1z13 * y4        + s2z14 * y5; // z|yz 
        current_data[26] = s0z15            + s1z16             + s2z17; // z|zz 
      }
    }
// most likely have a bug
#if 0
  } else if (ang == 32801) { // (1, 0, 1, 1)
    double x[18];
    double y[18];
    double z[18];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = data_ + data_offset_ii;

      const int ii3 = 3 * ii;

      const double oxp2 = 0.5 / xp_[ii];
      const double opq = 1.0 / (xp_[ii] + xq_[ii]);
      const double xqopq = xq_[ii] * opq;
      const double c00i0x = p_[ii3]     - ax;
      const double c00i0y = p_[ii3 + 1] - ay;
      const double c00i0z = p_[ii3 + 2] - az;
      const double c00i1x = (p_[ii3]     - q_[ii3])     * xqopq; 
      const double c00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xqopq; 
      const double c00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xqopq; 
      const double b10i0 = xqopq * oxp2; 
      const double b10_0 = oxp2 - b10i0 * roots_[offset];
      const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
      const double b10_2 = oxp2 - b10i0 * roots_[offset + 2];
      const double oxq2 = 0.5 / xq_[ii];
      const double xpopq = xp_[ii] * opq;
      const double d00i0x = q_[ii3] - cx;
      const double d00i0y = q_[ii3 + 1] - cy;
      const double d00i0z = q_[ii3 + 2] - cz;
      const double d00i1x = (p_[ii3] - q_[ii3]) * xpopq; 
      const double d00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xpopq; 
      const double d00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xpopq; 
      const double b00int = 0.5 * opq;
      const double b00_0 = b00int * roots_[offset];
      const double b00_1 = b00int * roots_[offset + 1];
      const double b00_2 = b00int * roots_[offset + 2];
      const double b01i0 = xpopq * oxq2; 
      const double b01_0 = oxq2 - b01i0 * roots_[offset];
      const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
      const double b01_2 = oxq2 - b01i0 * roots_[offset + 2];

      x[3] = c00i0x - c00i1x * roots_[offset]; 
      x[4] = c00i0x - c00i1x * roots_[offset + 1]; 
      x[5] = c00i0x - c00i1x * roots_[offset + 2]; 
      x[6] = d00i0x + d00i1x * roots_[offset]; 
      x[7] = d00i0x + d00i1x * roots_[offset + 1]; 
      x[8] = d00i0x + d00i1x * roots_[offset + 2]; 
      x[9] = x[3] * x[6]   + b00_0;
      x[10] = x[4] * x[7]   + b00_1;
      x[11] = x[5] * x[8]   + b00_2;
      x[12] = x[6] * x[6]   + b01_0; 
      x[13] = x[7] * x[7]   + b01_1; 
      x[14] = x[8] * x[8]   + b01_2; 
      x[15] = x[3] * x[12]  + 2.0 * b00_0 * x[6];
      x[16] = x[4] * x[13]  + 2.0 * b00_1 * x[7];
      x[17] = x[5] * x[14]  + 2.0 * b00_2 * x[8];

      y[3] = c00i0y - c00i1y * roots_[offset]; 
      y[4] = c00i0y - c00i1y * roots_[offset + 1]; 
      y[5] = c00i0y - c00i1y * roots_[offset + 2]; 
      y[6] = d00i0y + d00i1y * roots_[offset]; 
      y[7] = d00i0y + d00i1y * roots_[offset + 1]; 
      y[8] = d00i0y + d00i1y * roots_[offset + 2]; 
      y[9]  = y[3] * y[6]   + b00_0;
      y[10] = y[4] * y[7]   + b00_1;
      y[11] = y[5] * y[8]   + b00_2;
      y[12] = y[6] * y[6]   + b01_0; 
      y[13] = y[7] * y[7]   + b01_1; 
      y[14] = y[8] * y[8]   + b01_2; 
      y[15] = y[3] * y[12]  + 2.0 * b00_0 * y[6];
      y[16] = y[4] * y[13]  + 2.0 * b00_1 * y[7];
      y[17] = y[5] * y[14]  + 2.0 * b00_2 * y[8];

      z[3] = c00i0z - c00i1z * roots_[offset]; 
      z[4] = c00i0z - c00i1z * roots_[offset + 1]; 
      z[5] = c00i0z - c00i1z * roots_[offset + 2]; 
      z[6] = d00i0z + d00i1z * roots_[offset]; 
      z[7] = d00i0z + d00i1z * roots_[offset + 1]; 
      z[8] = d00i0z + d00i1z * roots_[offset + 2]; 
      z[9] = z[3] * z[6]   + b00_0;
      z[10] = z[4] * z[7]   + b00_1;
      z[11] = z[5] * z[8]   + b00_2;
      z[12] = z[6] * z[6]   + b01_0; 
      z[13] = z[7] * z[7]   + b01_1; 
      z[14] = z[8] * z[8]   + b01_2; 
      z[15] = z[3] * z[12]  + 2.0 * b00_0 * z[6];
      z[16] = z[4] * z[13]  + 2.0 * b00_1 * z[7];
      z[17] = z[5] * z[14]  + 2.0 * b00_2 * z[8];

      const double s0 = coeff_[ii] * weights_[offset    ] * (1.0 - roots_[offset]);
      const double s1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
      const double s2 = coeff_[ii] * weights_[offset + 2] * (1.0 - roots_[offset + 2]);
      const double s0x3 = s0 * x[3];
      const double s1x4 = s1 * x[4];
      const double s2x5 = s2 * x[5];
      const double s0x6 = s0 * x[6];
      const double s1x7 = s1 * x[7];
      const double s2x8 = s2 * x[8];
      const double s0x9 = s0 * x[9];
      const double s1x10 = s1 * x[10];
      const double s2x11 = s2 * x[11];
      const double s0x12 = s0 * x[12];
      const double s1x13 = s1 * x[13];
      const double s2x14 = s2 * x[14];
      const double s0x15 = s0 * x[15];
      const double s1x16 = s1 * x[16];
      const double s2x17 = s2 * x[17];
      const double s0x3y6 = s0x3 * y[6]; 
      const double s1x4y7 = s1x4 * y[7]; 
      const double s2x5y8 = s2x5 * y[8]; 
      const double s0y2 = s0 * y[3];
      const double s1y4 = s1 * y[4];
      const double s2y5 = s2 * y[5];
      const double s0y6 = s0 * y[6];
      const double s1y7 = s1 * y[7];
      const double s2y8 = s2 * y[8];
      const double s0y9  = s0 * y[9];
      const double s1y10 = s1 * y[10];
      const double s2y11 = s2 * y[11];

      current_data[0] = s0x9              + s1x10             + s2x11; // x|x
      current_data[1] = s0x6 * y[3]       + s1x7 * y[4]       + s2x8 * y[5]; // x|y
      current_data[2] = s0x6 * z[3]       + s1x7 * z[4]       + s2x8 * z[5]; // x|z
      current_data[3] = s0x3y6            + s1x4y7            + s2x5y8; // y|x
      current_data[4] = s0y9              + s1y10             + s2y11; // y|y
      current_data[5] = s0y6 * z[3]       + s1y7 * z[4]       + s2y8 * z[5]; // y|z 
      current_data[6] = s0x3 * z[6]       + s1x4 * z[7]       + s2x5 * z[8]; // z|x
      current_data[7] = s0y2 * z[6]       + s1y4 * z[7]       + s2y5 * z[8]; // z|y
      current_data[8] = s0 * z[9]         + s1 * z[10]         + s2 * z[11]; // z|z
   
      current_data[ 9] = s0x15              + s1x16               + s2x17; // xx|x
      current_data[10] = s0x12 * y[3]       + s1x13 * y[4]        + s2x14 * y[5]; // xx|y
      current_data[11] = s0x12 * z[3]       + s1x13 * z[4]        + s2x14 * z[5]; // xx|z
      current_data[12] = s0x9 * y[6]        + s1x10 * y[7]        + s2x11 * y[8]; // xy|x
      current_data[13] = s0x6 * y[9]        + s1x7 * y[10]        + s2x8 * y[11]; // xy|y
      current_data[14] = s0x6 * y[6] * z[3] + s1x7 * y[7] * z[4]  + s2x8 * y[8] * z[5]; // xy|z
      current_data[15] = s0x3 * y[12]       + s1x4 * y[13]        + s2x5 * y[14]; // yy|x
      current_data[16] = s0 * y[15]         + s1 * y[16]          + s2 * y[17]; // yy|y
      current_data[17] = s0 * y[12] * z[3]  + s1 * y[13] * z[4]   + s2 * y[14] * z[5]; // yy|z
      current_data[18] = s0x9 * z[6]        + s1x10 * z[7]        + s2x11 * z[8]; // xz|x
      current_data[19] = s0x6 * z[6] * y[3] + s1x7 * z[7] * y[4]  + s2x8 * z[8] * y[5]; // xz|y
      current_data[20] = s0x6 * z[9]        + s1x7 * z[10]         + s2x8 * z[11]; // xz|z
      current_data[21] = s0x3y6 * z[6]      + s1x4y7 * z[7]       + s2x5y8; // yz|x
      current_data[22] = s0y9 * z[6]        + s1y10 * z[7]        + s2y11 * z[8]; // yz|y
      current_data[23] = s0y6 * z[9]        + s1y7 * z[10]         + s2y8 * z[11]; // yz|z
      current_data[24] = s0x3 * z[12]        + s1x4 * z[13]         + s2x5 * z[14]; // zz|x
      current_data[25] = s0y2 * z[12]        + s1y4 * z[13]         + s2y5 * z[14]; // zz|y
      current_data[26] = s0 * z[15]         + s1 * z[16]          + s2 * z[17]; // zz|z
    } 
  } else if (ang == 32832) { // (1, 0, 2, 0)
    double x[18];
    double y[18];
    double z[18];
    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      const int offset = ii * rank_;
      const int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;

      const double oxp2 = 0.5 / xp_[ii];
      const double opq = 1.0 / (xp_[ii] + xq_[ii]);
      const double xqopq = xq_[ii] * opq;
      const double c00i0x = p_[ii3]     - ax;
      const double c00i0y = p_[ii3 + 1] - ay;
      const double c00i0z = p_[ii3 + 2] - az;
      const double c00i1x = (p_[ii3]     - q_[ii3])     * xqopq; 
      const double c00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xqopq; 
      const double c00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xqopq; 
      const double b10i0 = xqopq * oxp2; 
      const double b10_0 = oxp2 - b10i0 * roots_[offset];
      const double b10_1 = oxp2 - b10i0 * roots_[offset + 1];
      const double b10_2 = oxp2 - b10i0 * roots_[offset + 2];
      const double oxq2 = 0.5 / xq_[ii];
      const double xpopq = xp_[ii] * opq;
      const double d00i0x = q_[ii3] - cx;
      const double d00i0y = q_[ii3 + 1] - cy;
      const double d00i0z = q_[ii3 + 2] - cz;
      const double d00i1x = (p_[ii3] - q_[ii3]) * xpopq; 
      const double d00i1y = (p_[ii3 + 1] - q_[ii3 + 1]) * xpopq; 
      const double d00i1z = (p_[ii3 + 2] - q_[ii3 + 2]) * xpopq; 
      const double b00int = 0.5 * opq;
      const double b00_0 = b00int * roots_[offset];
      const double b00_1 = b00int * roots_[offset + 1];
      const double b00_2 = b00int * roots_[offset + 2];
      const double b01i0 = xpopq * oxq2; 
      const double b01_0 = oxq2 - b01i0 * roots_[offset];
      const double b01_1 = oxq2 - b01i0 * roots_[offset + 1];
      const double b01_2 = oxq2 - b01i0 * roots_[offset + 2];

      x[3] = c00i0x - c00i1x * roots_[offset]; 
      x[4] = c00i0x - c00i1x * roots_[offset + 1]; 
      x[5] = c00i0x - c00i1x * roots_[offset + 2]; 
      x[6] = d00i0x + d00i1x * roots_[offset]; 
      x[7] = d00i0x + d00i1x * roots_[offset + 1]; 
      x[8] = d00i0x + d00i1x * roots_[offset + 2]; 
      x[9] = x[3] * x[6]   + b00_0;
      x[10] = x[4] * x[7]   + b00_1;
      x[11] = x[5] * x[8]   + b00_2;
      x[12] = x[6] * x[6]   + b01_0; 
      x[13] = x[7] * x[7]   + b01_1; 
      x[14] = x[8] * x[8]   + b01_2; 
      x[15] = x[3] * x[12]  + 2.0 * b00_0 * x[6];
      x[16] = x[4] * x[13]  + 2.0 * b00_1 * x[7];
      x[17] = x[5] * x[14]  + 2.0 * b00_2 * x[8];

      y[3] = c00i0y - c00i1y * roots_[offset]; 
      y[4] = c00i0y - c00i1y * roots_[offset + 1]; 
      y[5] = c00i0y - c00i1y * roots_[offset + 2]; 
      y[6] = d00i0y + d00i1y * roots_[offset]; 
      y[7] = d00i0y + d00i1y * roots_[offset + 1]; 
      y[8] = d00i0y + d00i1y * roots_[offset + 2]; 
      y[9]  = y[3] * y[6]   + b00_0;
      y[10] = y[4] * y[7]   + b00_1;
      y[11] = y[5] * y[8]   + b00_2;
      y[12] = y[6] * y[6]   + b01_0; 
      y[13] = y[7] * y[7]   + b01_1; 
      y[14] = y[8] * y[8]   + b01_2; 
      y[15] = y[3] * y[12]  + 2.0 * b00_0 * y[6];
      y[16] = y[4] * y[13]  + 2.0 * b00_1 * y[7];
      y[17] = y[5] * y[14]  + 2.0 * b00_2 * y[8];

      z[3] = c00i0z - c00i1z * roots_[offset]; 
      z[4] = c00i0z - c00i1z * roots_[offset + 1]; 
      z[5] = c00i0z - c00i1z * roots_[offset + 2]; 
      z[6] = d00i0z + d00i1z * roots_[offset]; 
      z[7] = d00i0z + d00i1z * roots_[offset + 1]; 
      z[8] = d00i0z + d00i1z * roots_[offset + 2]; 
      z[9] = z[3] * z[6]   + b00_0;
      z[10] = z[4] * z[7]   + b00_1;
      z[11] = z[5] * z[8]   + b00_2;
      z[12] = z[6] * z[6]   + b01_0; 
      z[13] = z[7] * z[7]   + b01_1; 
      z[14] = z[8] * z[8]   + b01_2; 
      z[15] = z[3] * z[12]  + 2.0 * b00_0 * z[6];
      z[16] = z[4] * z[13]  + 2.0 * b00_1 * z[7];
      z[17] = z[5] * z[14]  + 2.0 * b00_2 * z[8];

      const double s0 = coeff_[ii] * weights_[offset    ] * (1.0 - roots_[offset]);
      const double s1 = coeff_[ii] * weights_[offset + 1] * (1.0 - roots_[offset + 1]);
      const double s2 = coeff_[ii] * weights_[offset + 2] * (1.0 - roots_[offset + 2]);
      const double s0x3 = s0 * x[3];
      const double s1x4 = s1 * x[4];
      const double s2x5 = s2 * x[5];
      const double s0x6 = s0 * x[6];
      const double s1x7 = s1 * x[7];
      const double s2x8 = s2 * x[8];
      const double s0x9 = s0 * x[9];
      const double s1x10 = s1 * x[10];
      const double s2x11 = s2 * x[11];
      const double s0x12 = s0 * x[12];
      const double s1x13 = s1 * x[13];
      const double s2x14 = s2 * x[14];
      const double s0x15 = s0 * x[15];
      const double s1x16 = s1 * x[16];
      const double s2x17 = s2 * x[17];
      const double s0x3y6 = s0x3 * y[6]; 
      const double s1x4y7 = s1x4 * y[7]; 
      const double s2x5y8 = s2x5 * y[8]; 
      const double s0y2 = s0 * y[3];
      const double s1y4 = s1 * y[4];
      const double s2y5 = s2 * y[5];
      const double s0y6 = s0 * y[6];
      const double s1y7 = s1 * y[7];
      const double s2y8 = s2 * y[8];
      const double s0y9  = s0 * y[9];
      const double s1y10 = s1 * y[10];
      const double s2y11 = s2 * y[11];

      current_data[ 0] = s0x15              + s1x16               + s2x17; // xx|x
      current_data[ 1] = s0x12 * y[3]       + s1x13 * y[4]        + s2x14 * y[5]; // xx|y
      current_data[ 2] = s0x12 * z[3]       + s1x13 * z[4]        + s2x14 * z[5]; // xx|z
      current_data[ 3] = s0x9 * y[6]        + s1x10 * y[7]        + s2x11 * y[8]; // xy|x
      current_data[ 4] = s0x6 * y[9]        + s1x7 * y[10]        + s2x8 * y[11]; // xy|y
      current_data[ 5] = s0x6 * y[6] * z[3] + s1x7 * y[7] * z[4]  + s2x8 * y[8] * z[5]; // xy|z
      current_data[ 6] = s0x3 * y[12]       + s1x4 * y[13]        + s2x5 * y[14]; // yy|x
      current_data[ 7] = s0 * y[15]         + s1 * y[16]          + s2 * y[17]; // yy|y
      current_data[ 8] = s0 * y[12] * z[3]  + s1 * y[13] * z[4]   + s2 * y[14] * z[5]; // yy|z
      current_data[ 9] = s0x9 * z[6]        + s1x10 * z[7]        + s2x11 * z[8]; // xz|x
      current_data[10] = s0x6 * z[6] * y[3] + s1x7 * z[7] * y[4]  + s2x8 * z[8] * y[5]; // xz|y
      current_data[11] = s0x6 * z[9]        + s1x7 * z[10]        + s2x8 * z[11]; // xz|z
      current_data[12] = s0x3y6 * z[6]      + s1x4y7 * z[7]       + s2x5y8; // yz|x
      current_data[13] = s0y9 * z[6]        + s1y10 * z[7]        + s2y11 * z[8]; // yz|y
      current_data[14] = s0y6 * z[9]        + s1y7 * z[10]        + s2y8 * z[11]; // yz|z
      current_data[15] = s0x3 * z[12]       + s1x4 * z[13]        + s2x5 * z[14]; // zz|x
      current_data[16] = s0y2 * z[12]       + s1y4 * z[13]        + s2y5 * z[14]; // zz|y
      current_data[17] = s0 * z[15]         + s1 * z[16]          + s2 * z[17]; // zz|z
    }
   
  } else {
#endif
#else
  {
    const int acsize = asize_ * csize_;
    const double ax = basisinfo_[0]->position(0);
    const double ay = basisinfo_[0]->position(1);
    const double az = basisinfo_[0]->position(2);
    const double bx = basisinfo_[1]->position(0);
    const double by = basisinfo_[1]->position(1);
    const double bz = basisinfo_[1]->position(2);
    const double cx = basisinfo_[2]->position(0);
    const double cy = basisinfo_[2]->position(1);
    const double cz = basisinfo_[2]->position(2);
    const double dx = basisinfo_[3]->position(0);
    const double dy = basisinfo_[3]->position(1);
    const double dz = basisinfo_[3]->position(2);
#endif
    const int isize = (amax_ + 1) * (cmax_ + 1);
    const int worksize = 3 * isize;
    const int svrr_index = amax_ * ANG_VRR_END + cmax_;

    double* workx = new double[worksize];
    double* worky = new double[worksize];
    double* workz = new double[worksize];
    double iyiz[3];

    for (int j = 0; j != screening_size_; ++j) {
      const int ii = screening_[j];
      int offset = ii * rank_;
      int data_offset_ii = ii * acsize;

      double* current_data = &data_[data_offset_ii];

      const int ii3 = 3 * ii;
      const double cxp = xp_[ii];
      const double cxq = xq_[ii];
      const double oxp2 = 0.5 / cxp;
      const double oxq2 = 0.5 / cxq;
      const double opq = 1.0 / (cxp + cxq);
      const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
      Int2D cix(dparamx, &roots_[offset], 3, worksize, workx, vrr_->vrrfunc[svrr_index]);
      double womt[3];
      double wt[3];
      wt[0] = weights_[offset + 0] * roots_[offset + 0];
      wt[1] = weights_[offset + 1] * roots_[offset + 1];
      wt[2] = weights_[offset + 2] * roots_[offset + 2];
      womt[0] = weights_[offset + 0]  - wt[0];
      womt[1] = weights_[offset + 1]  - wt[1];
      womt[2] = weights_[offset + 2]  - wt[2];
      cix.scale_data(womt, coeff_[ii]);

      const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
      Int2D ciy(dparamy, &roots_[offset], 3, worksize, worky, vrr_->vrrfunc[svrr_index]);

      const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
      Int2D ciz(dparamz, &roots_[offset], 3, worksize, workz, vrr_->vrrfunc[svrr_index]);

      for (int iz = 0; iz <= cmax_; ++iz) { 
        for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
          const int iyz = cmax1_ * (iy + cmax1_ * iz);
          for (int jz = 0; jz <= amax_; ++jz) { 
            const int offsetz = rank_ * (amax1_ * iz + jz);
            for (int jy = 0; jy <= amax_ - jz; ++jy) { 
              const int offsety = rank_ * (amax1_ * iy + jy);
              const int jyz = amax1_ * (jy + amax1_ * jz);
              iyiz[0] = worky[offsety + 0] * workz[offsetz + 0]; 
              iyiz[1] = worky[offsety + 1] * workz[offsetz + 1]; 
              iyiz[2] = worky[offsety + 2] * workz[offsetz + 2]; 
              for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
                const int iposition = cmapping_[ix + iyz];
                const int ipos_asize = iposition * asize_;
                for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                  const int offsetx = rank_ * (amax1_ * ix + jx);
                  const int jposition = amapping_[jx + jyz];
                  const int ijposition = jposition + ipos_asize;

                  current_data[ijposition] = iyiz[0] * workx[offsetx + 0]; 
                  current_data[ijposition] += iyiz[1] * workx[offsetx + 1]; 
                  current_data[ijposition] += iyiz[2] * workx[offsetx + 2]; 
                }
              }
            }
          }
        }
      }

    }

    delete[] workz;
    delete[] worky;
    delete[] workx;
  }
}


