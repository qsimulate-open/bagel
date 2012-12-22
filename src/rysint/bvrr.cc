//
// BAGEL - Parallel electron correlation program.
// Filename: bvrr.cc
// Copyright (C) 2012 Toru Shiozaki
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


#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/rysint/breitbatch.h>
#include <src/rysint/int2d.h>

using namespace std;
using namespace bagel;


void BreitBatch::perform_VRR() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = rank_ * isize;
  // CAUTION we need up to amax1_, cmax1_
  const int vrr_index = amax1_ * ANG_VRR_END + cmax1_;

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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const worktx = stack_->get(worksize*3);
  double* const workty = worktx + worksize;
  double* const worktz = workty + worksize;
  double* const worksx = stack_->get(worksize*3);
  double* const worksy = worksx + worksize;
  double* const worksz = worksy + worksize;

  const int acsize = size_block_ / primsize_;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * rank_;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D cix(dparamx, &roots_[offset], rank_, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], rank_, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], rank_, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute \tidle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != rank_; ++i) {
          worktx[i+rank_*(ia+amax2*ic)] = pq[0]*workx[i+rank_*(ia+amax2*ic)] + (ia==0 ? 0.0 : oxp2*workx[i+rank_*(ia+amax2*ic)]) - (ic==0 ? 0.0 : oxq2*workx[i+rank_*(ia+amax2*(ic-1))]); 
          workty[i+rank_*(ia+amax2*ic)] = pq[1]*worky[i+rank_*(ia+amax2*ic)] + (ia==0 ? 0.0 : oxp2*worky[i+rank_*(ia+amax2*ic)]) - (ic==0 ? 0.0 : oxq2*worky[i+rank_*(ia+amax2*(ic-1))]); 
          worktz[i+rank_*(ia+amax2*ic)] = pq[2]*workz[i+rank_*(ia+amax2*ic)] + (ia==0 ? 0.0 : oxp2*workz[i+rank_*(ia+amax2*ic)]) - (ic==0 ? 0.0 : oxq2*workz[i+rank_*(ia+amax2*(ic-1))]); 
        }
    // then compute \tilde{\tilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != rank_; ++i) {
          worksx[i+rank_*(ia+amax2*ic)] = worktx[i+rank_*((ia+1)+amax2*ic)] - worktx[i+rank_*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+rank_*(ia+amax2*ic)]; 
          worksy[i+rank_*(ia+amax2*ic)] = workty[i+rank_*((ia+1)+amax2*ic)] - workty[i+rank_*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+rank_*(ia+amax2*ic)]; 
          worksz[i+rank_*(ia+amax2*ic)] = worktz[i+rank_*((ia+1)+amax2*ic)] - worktz[i+rank_*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+rank_*(ia+amax2*ic)]; 
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_; 
    double* const datayy = dataxy + size_block_; 
    double* const dataxz = datayy + size_block_; 
    double* const datayz = dataxz + size_block_; 
    double* const datazz = datayz + size_block_; 

    double* const iyiz_nn = stack_->get(rank_*6);
    double* const iyiz_tn = iyiz_nn + rank_;
    double* const iyiz_nt = iyiz_tn + rank_;
    double* const iyiz_tt = iyiz_nt + rank_;
    double* const iyiz_sn = iyiz_tt + rank_;
    double* const iyiz_ns = iyiz_sn + rank_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = rank_ * (amax2 * iz + jz);
            const int offsety = rank_ * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != rank_; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i];
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i];
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i];
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = rank_ * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(rank_, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(rank_, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(rank_, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(rank_, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(rank_, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(rank_, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

    stack_->release(rank_*6, iyiz_nn);
  }

  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}
