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

#include <src/rysint/breitbatch.h>
#include <src/rysint/int2d.h>
#include <src/parallel/resources.h>
#include <src/util/f77.h>

using namespace std;
using namespace bagel;


void BreitBatch::perform_VRR1() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 1 * isize;
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

  double* const iyiz_nn = stack_->get(1*6);
  double* const iyiz_tn = iyiz_nn + 1;
  double* const iyiz_nt = iyiz_tn + 1;
  double* const iyiz_tt = iyiz_nt + 1;
  double* const iyiz_sn = iyiz_tt + 1;
  double* const iyiz_ns = iyiz_sn + 1;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 1;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(1, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<1> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<1> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<1> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 1; ++i) {
          worktx[i+1*(ia+amax2*ic)] = pq[0]*workx[i+1*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+1*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+1*(ia+amax2*(ic-1))]);
          workty[i+1*(ia+amax2*ic)] = pq[1]*worky[i+1*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+1*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+1*(ia+amax2*(ic-1))]);
          worktz[i+1*(ia+amax2*ic)] = pq[2]*workz[i+1*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+1*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+1*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 1; ++i) {
          worksx[i+1*(ia+amax2*ic)] = worktx[i+1*((ia+1)+amax2*ic)] - worktx[i+1*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+1*(ia+amax2*ic)];
          worksy[i+1*(ia+amax2*ic)] = workty[i+1*((ia+1)+amax2*ic)] - workty[i+1*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+1*(ia+amax2*ic)];
          worksz[i+1*(ia+amax2*ic)] = worktz[i+1*((ia+1)+amax2*ic)] - worktz[i+1*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+1*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 1 * (amax2 * iz + jz);
            const int offsety = 1 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 1; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 1 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(1, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(1, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(1, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(1, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(1, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(1, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(1*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


void BreitBatch::perform_VRR2() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 2 * isize;
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

  double* const iyiz_nn = stack_->get(2*6);
  double* const iyiz_tn = iyiz_nn + 2;
  double* const iyiz_nt = iyiz_tn + 2;
  double* const iyiz_tt = iyiz_nt + 2;
  double* const iyiz_sn = iyiz_tt + 2;
  double* const iyiz_ns = iyiz_sn + 2;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 2;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(2, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<2> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<2> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<2> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 2; ++i) {
          worktx[i+2*(ia+amax2*ic)] = pq[0]*workx[i+2*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+2*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+2*(ia+amax2*(ic-1))]);
          workty[i+2*(ia+amax2*ic)] = pq[1]*worky[i+2*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+2*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+2*(ia+amax2*(ic-1))]);
          worktz[i+2*(ia+amax2*ic)] = pq[2]*workz[i+2*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+2*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+2*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 2; ++i) {
          worksx[i+2*(ia+amax2*ic)] = worktx[i+2*((ia+1)+amax2*ic)] - worktx[i+2*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+2*(ia+amax2*ic)];
          worksy[i+2*(ia+amax2*ic)] = workty[i+2*((ia+1)+amax2*ic)] - workty[i+2*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+2*(ia+amax2*ic)];
          worksz[i+2*(ia+amax2*ic)] = worktz[i+2*((ia+1)+amax2*ic)] - worktz[i+2*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+2*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 2 * (amax2 * iz + jz);
            const int offsety = 2 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 2; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 2 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(2, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(2, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(2, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(2, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(2, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(2, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(2*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


void BreitBatch::perform_VRR3() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 3 * isize;
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

  double* const iyiz_nn = stack_->get(3*6);
  double* const iyiz_tn = iyiz_nn + 3;
  double* const iyiz_nt = iyiz_tn + 3;
  double* const iyiz_tt = iyiz_nt + 3;
  double* const iyiz_sn = iyiz_tt + 3;
  double* const iyiz_ns = iyiz_sn + 3;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 3;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(3, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<3> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<3> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<3> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 3; ++i) {
          worktx[i+3*(ia+amax2*ic)] = pq[0]*workx[i+3*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+3*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+3*(ia+amax2*(ic-1))]);
          workty[i+3*(ia+amax2*ic)] = pq[1]*worky[i+3*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+3*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+3*(ia+amax2*(ic-1))]);
          worktz[i+3*(ia+amax2*ic)] = pq[2]*workz[i+3*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+3*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+3*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 3; ++i) {
          worksx[i+3*(ia+amax2*ic)] = worktx[i+3*((ia+1)+amax2*ic)] - worktx[i+3*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+3*(ia+amax2*ic)];
          worksy[i+3*(ia+amax2*ic)] = workty[i+3*((ia+1)+amax2*ic)] - workty[i+3*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+3*(ia+amax2*ic)];
          worksz[i+3*(ia+amax2*ic)] = worktz[i+3*((ia+1)+amax2*ic)] - worktz[i+3*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+3*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 3 * (amax2 * iz + jz);
            const int offsety = 3 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 3; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 3 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(3, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(3, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(3, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(3, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(3, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(3, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(3*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


void BreitBatch::perform_VRR4() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 4 * isize;
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

  double* const iyiz_nn = stack_->get(4*6);
  double* const iyiz_tn = iyiz_nn + 4;
  double* const iyiz_nt = iyiz_tn + 4;
  double* const iyiz_tt = iyiz_nt + 4;
  double* const iyiz_sn = iyiz_tt + 4;
  double* const iyiz_ns = iyiz_sn + 4;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 4;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(4, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<4> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<4> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<4> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 4; ++i) {
          worktx[i+4*(ia+amax2*ic)] = pq[0]*workx[i+4*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+4*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+4*(ia+amax2*(ic-1))]);
          workty[i+4*(ia+amax2*ic)] = pq[1]*worky[i+4*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+4*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+4*(ia+amax2*(ic-1))]);
          worktz[i+4*(ia+amax2*ic)] = pq[2]*workz[i+4*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+4*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+4*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 4; ++i) {
          worksx[i+4*(ia+amax2*ic)] = worktx[i+4*((ia+1)+amax2*ic)] - worktx[i+4*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+4*(ia+amax2*ic)];
          worksy[i+4*(ia+amax2*ic)] = workty[i+4*((ia+1)+amax2*ic)] - workty[i+4*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+4*(ia+amax2*ic)];
          worksz[i+4*(ia+amax2*ic)] = worktz[i+4*((ia+1)+amax2*ic)] - worktz[i+4*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+4*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 4 * (amax2 * iz + jz);
            const int offsety = 4 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 4; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 4 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(4, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(4, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(4, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(4, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(4, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(4, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(4*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


void BreitBatch::perform_VRR5() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 5 * isize;
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

  double* const iyiz_nn = stack_->get(5*6);
  double* const iyiz_tn = iyiz_nn + 5;
  double* const iyiz_nt = iyiz_tn + 5;
  double* const iyiz_tt = iyiz_nt + 5;
  double* const iyiz_sn = iyiz_tt + 5;
  double* const iyiz_ns = iyiz_sn + 5;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 5;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(5, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<5> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<5> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<5> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 5; ++i) {
          worktx[i+5*(ia+amax2*ic)] = pq[0]*workx[i+5*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+5*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+5*(ia+amax2*(ic-1))]);
          workty[i+5*(ia+amax2*ic)] = pq[1]*worky[i+5*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+5*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+5*(ia+amax2*(ic-1))]);
          worktz[i+5*(ia+amax2*ic)] = pq[2]*workz[i+5*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+5*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+5*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 5; ++i) {
          worksx[i+5*(ia+amax2*ic)] = worktx[i+5*((ia+1)+amax2*ic)] - worktx[i+5*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+5*(ia+amax2*ic)];
          worksy[i+5*(ia+amax2*ic)] = workty[i+5*((ia+1)+amax2*ic)] - workty[i+5*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+5*(ia+amax2*ic)];
          worksz[i+5*(ia+amax2*ic)] = worktz[i+5*((ia+1)+amax2*ic)] - worktz[i+5*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+5*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 5 * (amax2 * iz + jz);
            const int offsety = 5 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 5; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 5 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(5, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(5, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(5, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(5, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(5, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(5, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(5*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


void BreitBatch::perform_VRR6() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 6 * isize;
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

  double* const iyiz_nn = stack_->get(6*6);
  double* const iyiz_tn = iyiz_nn + 6;
  double* const iyiz_nt = iyiz_tn + 6;
  double* const iyiz_tt = iyiz_nt + 6;
  double* const iyiz_sn = iyiz_tt + 6;
  double* const iyiz_ns = iyiz_sn + 6;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 6;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(6, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<6> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<6> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<6> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 6; ++i) {
          worktx[i+6*(ia+amax2*ic)] = pq[0]*workx[i+6*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+6*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+6*(ia+amax2*(ic-1))]);
          workty[i+6*(ia+amax2*ic)] = pq[1]*worky[i+6*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+6*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+6*(ia+amax2*(ic-1))]);
          worktz[i+6*(ia+amax2*ic)] = pq[2]*workz[i+6*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+6*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+6*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 6; ++i) {
          worksx[i+6*(ia+amax2*ic)] = worktx[i+6*((ia+1)+amax2*ic)] - worktx[i+6*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+6*(ia+amax2*ic)];
          worksy[i+6*(ia+amax2*ic)] = workty[i+6*((ia+1)+amax2*ic)] - workty[i+6*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+6*(ia+amax2*ic)];
          worksz[i+6*(ia+amax2*ic)] = worktz[i+6*((ia+1)+amax2*ic)] - worktz[i+6*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+6*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 6 * (amax2 * iz + jz);
            const int offsety = 6 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 6; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 6 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(6, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(6, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(6, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(6, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(6, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(6, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(6*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


void BreitBatch::perform_VRR7() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 7 * isize;
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

  double* const iyiz_nn = stack_->get(7*6);
  double* const iyiz_tn = iyiz_nn + 7;
  double* const iyiz_nt = iyiz_tn + 7;
  double* const iyiz_tt = iyiz_nt + 7;
  double* const iyiz_sn = iyiz_tt + 7;
  double* const iyiz_ns = iyiz_sn + 7;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 7;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(7, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<7> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<7> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<7> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 7; ++i) {
          worktx[i+7*(ia+amax2*ic)] = pq[0]*workx[i+7*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+7*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+7*(ia+amax2*(ic-1))]);
          workty[i+7*(ia+amax2*ic)] = pq[1]*worky[i+7*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+7*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+7*(ia+amax2*(ic-1))]);
          worktz[i+7*(ia+amax2*ic)] = pq[2]*workz[i+7*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+7*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+7*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 7; ++i) {
          worksx[i+7*(ia+amax2*ic)] = worktx[i+7*((ia+1)+amax2*ic)] - worktx[i+7*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+7*(ia+amax2*ic)];
          worksy[i+7*(ia+amax2*ic)] = workty[i+7*((ia+1)+amax2*ic)] - workty[i+7*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+7*(ia+amax2*ic)];
          worksz[i+7*(ia+amax2*ic)] = worktz[i+7*((ia+1)+amax2*ic)] - worktz[i+7*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+7*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 7 * (amax2 * iz + jz);
            const int offsety = 7 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 7; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 7 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(7, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(7, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(7, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(7, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(7, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(7, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(7*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


void BreitBatch::perform_VRR8() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 8 * isize;
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

  double* const iyiz_nn = stack_->get(8*6);
  double* const iyiz_tn = iyiz_nn + 8;
  double* const iyiz_nt = iyiz_tn + 8;
  double* const iyiz_tt = iyiz_nt + 8;
  double* const iyiz_sn = iyiz_tt + 8;
  double* const iyiz_ns = iyiz_sn + 8;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 8;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(8, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<8> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<8> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<8> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 8; ++i) {
          worktx[i+8*(ia+amax2*ic)] = pq[0]*workx[i+8*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+8*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+8*(ia+amax2*(ic-1))]);
          workty[i+8*(ia+amax2*ic)] = pq[1]*worky[i+8*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+8*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+8*(ia+amax2*(ic-1))]);
          worktz[i+8*(ia+amax2*ic)] = pq[2]*workz[i+8*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+8*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+8*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 8; ++i) {
          worksx[i+8*(ia+amax2*ic)] = worktx[i+8*((ia+1)+amax2*ic)] - worktx[i+8*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+8*(ia+amax2*ic)];
          worksy[i+8*(ia+amax2*ic)] = workty[i+8*((ia+1)+amax2*ic)] - workty[i+8*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+8*(ia+amax2*ic)];
          worksz[i+8*(ia+amax2*ic)] = worktz[i+8*((ia+1)+amax2*ic)] - worktz[i+8*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+8*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 8 * (amax2 * iz + jz);
            const int offsety = 8 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 8; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 8 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(8, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(8, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(8, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(8, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(8, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(8, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(8*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


void BreitBatch::perform_VRR9() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 9 * isize;
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

  double* const iyiz_nn = stack_->get(9*6);
  double* const iyiz_tn = iyiz_nn + 9;
  double* const iyiz_nt = iyiz_tn + 9;
  double* const iyiz_tt = iyiz_nt + 9;
  double* const iyiz_sn = iyiz_tt + 9;
  double* const iyiz_ns = iyiz_sn + 9;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 9;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(9, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<9> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<9> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<9> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 9; ++i) {
          worktx[i+9*(ia+amax2*ic)] = pq[0]*workx[i+9*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+9*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+9*(ia+amax2*(ic-1))]);
          workty[i+9*(ia+amax2*ic)] = pq[1]*worky[i+9*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+9*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+9*(ia+amax2*(ic-1))]);
          worktz[i+9*(ia+amax2*ic)] = pq[2]*workz[i+9*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+9*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+9*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 9; ++i) {
          worksx[i+9*(ia+amax2*ic)] = worktx[i+9*((ia+1)+amax2*ic)] - worktx[i+9*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+9*(ia+amax2*ic)];
          worksy[i+9*(ia+amax2*ic)] = workty[i+9*((ia+1)+amax2*ic)] - workty[i+9*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+9*(ia+amax2*ic)];
          worksz[i+9*(ia+amax2*ic)] = worktz[i+9*((ia+1)+amax2*ic)] - worktz[i+9*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+9*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 9 * (amax2 * iz + jz);
            const int offsety = 9 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 9; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 9 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(9, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(9, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(9, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(9, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(9, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(9, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(9*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


void BreitBatch::perform_VRR10() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 10 * isize;
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

  double* const iyiz_nn = stack_->get(10*6);
  double* const iyiz_tn = iyiz_nn + 10;
  double* const iyiz_nt = iyiz_tn + 10;
  double* const iyiz_tt = iyiz_nt + 10;
  double* const iyiz_sn = iyiz_tt + 10;
  double* const iyiz_ns = iyiz_sn + 10;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 10;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(10, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<10> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<10> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<10> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 10; ++i) {
          worktx[i+10*(ia+amax2*ic)] = pq[0]*workx[i+10*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+10*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+10*(ia+amax2*(ic-1))]);
          workty[i+10*(ia+amax2*ic)] = pq[1]*worky[i+10*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+10*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+10*(ia+amax2*(ic-1))]);
          worktz[i+10*(ia+amax2*ic)] = pq[2]*workz[i+10*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+10*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+10*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 10; ++i) {
          worksx[i+10*(ia+amax2*ic)] = worktx[i+10*((ia+1)+amax2*ic)] - worktx[i+10*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+10*(ia+amax2*ic)];
          worksy[i+10*(ia+amax2*ic)] = workty[i+10*((ia+1)+amax2*ic)] - workty[i+10*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+10*(ia+amax2*ic)];
          worksz[i+10*(ia+amax2*ic)] = worktz[i+10*((ia+1)+amax2*ic)] - worktz[i+10*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+10*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 10 * (amax2 * iz + jz);
            const int offsety = 10 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 10; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 10 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(10, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(10, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(10, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(10, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(10, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(10, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(10*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


void BreitBatch::perform_VRR11() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 11 * isize;
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

  double* const iyiz_nn = stack_->get(11*6);
  double* const iyiz_tn = iyiz_nn + 11;
  double* const iyiz_nt = iyiz_tn + 11;
  double* const iyiz_tt = iyiz_nt + 11;
  double* const iyiz_sn = iyiz_tt + 11;
  double* const iyiz_ns = iyiz_sn + 11;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 11;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(11, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<11> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<11> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<11> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 11; ++i) {
          worktx[i+11*(ia+amax2*ic)] = pq[0]*workx[i+11*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+11*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+11*(ia+amax2*(ic-1))]);
          workty[i+11*(ia+amax2*ic)] = pq[1]*worky[i+11*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+11*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+11*(ia+amax2*(ic-1))]);
          worktz[i+11*(ia+amax2*ic)] = pq[2]*workz[i+11*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+11*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+11*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 11; ++i) {
          worksx[i+11*(ia+amax2*ic)] = worktx[i+11*((ia+1)+amax2*ic)] - worktx[i+11*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+11*(ia+amax2*ic)];
          worksy[i+11*(ia+amax2*ic)] = workty[i+11*((ia+1)+amax2*ic)] - workty[i+11*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+11*(ia+amax2*ic)];
          worksz[i+11*(ia+amax2*ic)] = worktz[i+11*((ia+1)+amax2*ic)] - worktz[i+11*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+11*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 11 * (amax2 * iz + jz);
            const int offsety = 11 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 11; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 11 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(11, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(11, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(11, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(11, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(11, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(11, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(11*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


void BreitBatch::perform_VRR12() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 12 * isize;
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

  double* const iyiz_nn = stack_->get(12*6);
  double* const iyiz_tn = iyiz_nn + 12;
  double* const iyiz_nt = iyiz_tn + 12;
  double* const iyiz_tt = iyiz_nt + 12;
  double* const iyiz_sn = iyiz_tt + 12;
  double* const iyiz_ns = iyiz_sn + 12;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 12;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(12, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<12> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<12> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<12> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 12; ++i) {
          worktx[i+12*(ia+amax2*ic)] = pq[0]*workx[i+12*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+12*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+12*(ia+amax2*(ic-1))]);
          workty[i+12*(ia+amax2*ic)] = pq[1]*worky[i+12*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+12*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+12*(ia+amax2*(ic-1))]);
          worktz[i+12*(ia+amax2*ic)] = pq[2]*workz[i+12*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+12*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+12*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 12; ++i) {
          worksx[i+12*(ia+amax2*ic)] = worktx[i+12*((ia+1)+amax2*ic)] - worktx[i+12*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+12*(ia+amax2*ic)];
          worksy[i+12*(ia+amax2*ic)] = workty[i+12*((ia+1)+amax2*ic)] - workty[i+12*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+12*(ia+amax2*ic)];
          worksz[i+12*(ia+amax2*ic)] = worktz[i+12*((ia+1)+amax2*ic)] - worktz[i+12*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+12*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 12 * (amax2 * iz + jz);
            const int offsety = 12 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 12; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 12 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(12, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(12, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(12, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(12, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(12, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(12, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(12*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


void BreitBatch::perform_VRR13() {
  const int amax2 = amax1_+1;
  const int cmax2 = cmax1_+1;
  const int isize = amax2 * cmax2;
  const int worksize = 13 * isize;
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

  double* const iyiz_nn = stack_->get(13*6);
  double* const iyiz_tn = iyiz_nn + 13;
  double* const iyiz_nt = iyiz_tn + 13;
  double* const iyiz_tt = iyiz_nt + 13;
  double* const iyiz_sn = iyiz_tt + 13;
  double* const iyiz_ns = iyiz_sn + 13;

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    const size_t offset = ii * 13;
    const size_t data_offset_ii = ii * acsize;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    dscal_(13, cxp*cxq*2.0*opq, weights_+offset, 1);

    const array<double, 11> dparamx = {{p_[ii3],   q_[ii3],   ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<13> cix(dparamx, roots_+offset, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(weights_+offset, coeff_[ii]);

    const array<double, 11> dparamy = {{p_[ii3+1], q_[ii3+1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<13> ciy(dparamy, roots_+offset, worksize, worky, vrr_->vrrfunc[vrr_index]);

    const array<double, 11> dparamz = {{p_[ii3+2], q_[ii3+2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<13> ciz(dparamz, roots_+offset, worksize, workz, vrr_->vrrfunc[vrr_index]);

    const double pq[3] = {p_[ii3]-q_[ii3], p_[ii3+1]-q_[ii3+1], p_[ii3+2]-q_[ii3+2]};

    // next compute 	idle{I}_x,y,z up to amax_, cmax_
    for (int ic = 0; ic <= cmax1_; ++ic)
      for (int ia = 0; ia <= amax1_; ++ia)
        for (int i = 0; i != 13; ++i) {
          worktx[i+13*(ia+amax2*ic)] = pq[0]*workx[i+13*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+13*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+13*(ia+amax2*(ic-1))]);
          workty[i+13*(ia+amax2*ic)] = pq[1]*worky[i+13*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+13*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+13*(ia+amax2*(ic-1))]);
          worktz[i+13*(ia+amax2*ic)] = pq[2]*workz[i+13*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+13*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+13*(ia+amax2*(ic-1))]);
        }
    // then compute 	ilde{	ilde{I}}_x,y,z up to amax_-1, cmax_-1
    for (int ic = 0; ic != cmax1_; ++ic)
      for (int ia = 0; ia != amax1_; ++ia)
        for (int i = 0; i != 13; ++i) {
          worksx[i+13*(ia+amax2*ic)] = worktx[i+13*((ia+1)+amax2*ic)] - worktx[i+13*(ia+amax2*(ic+1))] + (ax - cx)*worktx[i+13*(ia+amax2*ic)];
          worksy[i+13*(ia+amax2*ic)] = workty[i+13*((ia+1)+amax2*ic)] - workty[i+13*(ia+amax2*(ic+1))] + (ay - cy)*workty[i+13*(ia+amax2*ic)];
          worksz[i+13*(ia+amax2*ic)] = worktz[i+13*((ia+1)+amax2*ic)] - worktz[i+13*(ia+amax2*(ic+1))] + (az - cz)*worktz[i+13*(ia+amax2*ic)];
        }

    double* const dataxx = &data_[data_offset_ii];
    double* const dataxy = dataxx + size_block_;
    double* const datayy = dataxy + size_block_;
    double* const dataxz = datayy + size_block_;
    double* const datayz = dataxz + size_block_;
    double* const datazz = datayz + size_block_;


    // assemble up to amax_, cmax_
    for (int iz = 0; iz <= cmax_; ++iz) {
      for (int iy = 0; iy <= cmax_ - iz; ++iy) {
        for (int jz = 0; jz <= amax_; ++jz) {
          for (int jy = 0; jy <= amax_ - jz; ++jy) {
            const int offsetz = 13 * (amax2 * iz + jz);
            const int offsety = 13 * (amax2 * iy + jy);

            const int iyz = cmax1_ * (iy + cmax1_ * iz);
            const int jyz = amax1_ * (jy + amax1_ * jz);

            for (int i = 0; i != 13; ++i) {
              iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
              iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots_[offset+i]);
              iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
              iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
            }

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
                const int offsetx = 13 * (amax2 * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                dataxx[ijposition] = ddot_(13, iyiz_nn, 1, worksx+offsetx, 1);
                dataxy[ijposition] = ddot_(13, iyiz_tn, 1, worktx+offsetx, 1);
                dataxz[ijposition] = ddot_(13, iyiz_nt, 1, worktx+offsetx, 1);
                datayy[ijposition] = ddot_(13, iyiz_sn, 1, workx +offsetx, 1);
                datazz[ijposition] = ddot_(13, iyiz_ns, 1, workx +offsetx, 1);
                datayz[ijposition] = ddot_(13, iyiz_tt, 1, workx +offsetx, 1);
              }
            }
          }
        }
      }
    }

  }

  stack_->release(13*6, iyiz_nn);
  stack_->release(worksize*3, worksx);
  stack_->release(worksize*3, worktx);
  stack_->release(worksize*3, workx);

}


