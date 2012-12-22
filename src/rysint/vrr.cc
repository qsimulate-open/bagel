//
// BAGEL - Parallel electron correlation program.
// Filename: vrr.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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

#include <src/rysint/eribatch.h>
#include <src/rysint/int2d.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <src/parallel/resources.h>
#include <src/util/f77.h>

using namespace std;
using namespace bagel;


void ERIBatch::perform_VRR4() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 4 * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

  double* const workx = stack_->get(worksize);
  double* const worky = stack_->get(worksize);
  double* const workz = stack_->get(worksize);
  double iyiz[4];

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
  for (int j = 0; j != screening_size_; ++j) {
    int ii = screening_[j];
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
    Int2D cix(dparamx, &roots_[offset], 4, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], 4, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], 4, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) {
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            for (int i = 0; i != 4; ++i)
              iyiz[i] = worky[offsety + i] * workz[offsetz + i];
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;
                current_data[ijposition] = ddot_(4, iyiz, 1, workx+offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}


void ERIBatch::perform_VRR5() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 5 * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

  double* const workx = stack_->get(worksize);
  double* const worky = stack_->get(worksize);
  double* const workz = stack_->get(worksize);
  double iyiz[5];

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
  for (int j = 0; j != screening_size_; ++j) {
    int ii = screening_[j];
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
    Int2D cix(dparamx, &roots_[offset], 5, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], 5, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], 5, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) {
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            for (int i = 0; i != 5; ++i)
              iyiz[i] = worky[offsety + i] * workz[offsetz + i];
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;
                current_data[ijposition] = ddot_(5, iyiz, 1, workx+offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}


void ERIBatch::perform_VRR6() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 6 * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

  double* const workx = stack_->get(worksize);
  double* const worky = stack_->get(worksize);
  double* const workz = stack_->get(worksize);
  double iyiz[6];

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
  for (int j = 0; j != screening_size_; ++j) {
    int ii = screening_[j];
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
    Int2D cix(dparamx, &roots_[offset], 6, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], 6, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], 6, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) {
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            for (int i = 0; i != 6; ++i)
              iyiz[i] = worky[offsety + i] * workz[offsetz + i];
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;
                current_data[ijposition] = ddot_(6, iyiz, 1, workx+offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}


void ERIBatch::perform_VRR7() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 7 * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

  double* const workx = stack_->get(worksize);
  double* const worky = stack_->get(worksize);
  double* const workz = stack_->get(worksize);
  double iyiz[7];

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
  for (int j = 0; j != screening_size_; ++j) {
    int ii = screening_[j];
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
    Int2D cix(dparamx, &roots_[offset], 7, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], 7, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], 7, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) {
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            for (int i = 0; i != 7; ++i)
              iyiz[i] = worky[offsety + i] * workz[offsetz + i];
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;
                current_data[ijposition] = ddot_(7, iyiz, 1, workx+offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}


void ERIBatch::perform_VRR8() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 8 * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

  double* const workx = stack_->get(worksize);
  double* const worky = stack_->get(worksize);
  double* const workz = stack_->get(worksize);
  double iyiz[8];

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
  for (int j = 0; j != screening_size_; ++j) {
    int ii = screening_[j];
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
    Int2D cix(dparamx, &roots_[offset], 8, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], 8, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], 8, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) {
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            for (int i = 0; i != 8; ++i)
              iyiz[i] = worky[offsety + i] * workz[offsetz + i];
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;
                current_data[ijposition] = ddot_(8, iyiz, 1, workx+offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}


void ERIBatch::perform_VRR9() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 9 * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

  double* const workx = stack_->get(worksize);
  double* const worky = stack_->get(worksize);
  double* const workz = stack_->get(worksize);
  double iyiz[9];

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
  for (int j = 0; j != screening_size_; ++j) {
    int ii = screening_[j];
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
    Int2D cix(dparamx, &roots_[offset], 9, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], 9, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], 9, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) {
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            for (int i = 0; i != 9; ++i)
              iyiz[i] = worky[offsety + i] * workz[offsetz + i];
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;
                current_data[ijposition] = ddot_(9, iyiz, 1, workx+offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}


void ERIBatch::perform_VRR10() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 10 * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

  double* const workx = stack_->get(worksize);
  double* const worky = stack_->get(worksize);
  double* const workz = stack_->get(worksize);
  double iyiz[10];

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
  for (int j = 0; j != screening_size_; ++j) {
    int ii = screening_[j];
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
    Int2D cix(dparamx, &roots_[offset], 10, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], 10, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], 10, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) {
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            for (int i = 0; i != 10; ++i)
              iyiz[i] = worky[offsety + i] * workz[offsetz + i];
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;
                current_data[ijposition] = ddot_(10, iyiz, 1, workx+offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}


void ERIBatch::perform_VRR11() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 11 * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

  double* const workx = stack_->get(worksize);
  double* const worky = stack_->get(worksize);
  double* const workz = stack_->get(worksize);
  double iyiz[11];

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
  for (int j = 0; j != screening_size_; ++j) {
    int ii = screening_[j];
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
    Int2D cix(dparamx, &roots_[offset], 11, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], 11, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], 11, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) {
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            for (int i = 0; i != 11; ++i)
              iyiz[i] = worky[offsety + i] * workz[offsetz + i];
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;
                current_data[ijposition] = ddot_(11, iyiz, 1, workx+offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}


void ERIBatch::perform_VRR12() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 12 * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

  double* const workx = stack_->get(worksize);
  double* const worky = stack_->get(worksize);
  double* const workz = stack_->get(worksize);
  double iyiz[12];

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
  for (int j = 0; j != screening_size_; ++j) {
    int ii = screening_[j];
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
    Int2D cix(dparamx, &roots_[offset], 12, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], 12, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], 12, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) {
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            for (int i = 0; i != 12; ++i)
              iyiz[i] = worky[offsety + i] * workz[offsetz + i];
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;
                current_data[ijposition] = ddot_(12, iyiz, 1, workx+offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}


void ERIBatch::perform_VRR13() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = 13 * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

  double* const workx = stack_->get(worksize);
  double* const worky = stack_->get(worksize);
  double* const workz = stack_->get(worksize);
  double iyiz[13];

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
  for (int j = 0; j != screening_size_; ++j) {
    int ii = screening_[j];
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
    Int2D cix(dparamx, &roots_[offset], 13, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], 13, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], 13, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) {
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            for (int i = 0; i != 13; ++i)
              iyiz[i] = worky[offsety + i] * workz[offsetz + i];
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;
                current_data[ijposition] = ddot_(13, iyiz, 1, workx+offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}


