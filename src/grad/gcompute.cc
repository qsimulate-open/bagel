//
// Newint - Parallel electron correlation program.
// Filename: gcompute.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/stackmem.h>
#include <src/grad/gradbatch.h>

using namespace std;

extern StackMem* stack;

void GradBatch::compute() {

  bkup_ = stack->get(size_alloc_);
  fill(data_, data_+data_size(), 0.0);

  // perform VRR
  // data_ will contain the intermediates: prim01{ prim23{ xyz{ } } } 
#if 0
  switch (rank_) {
    case 1: perform_VRR1(); break;
    case 2: perform_VRR2(); break;
    case 3: perform_VRR3(); break;
    case 4: perform_VRR4(); break;
    case 5: perform_VRR5(); break;
    case 6: perform_VRR6(); break;
    case 7: perform_VRR7(); break;
    case 8: perform_VRR8(); break;
    case 9: perform_VRR9(); break;  
    case 10: perform_VRR10(); break;  
    case 11: perform_VRR11(); break;  
    case 12: perform_VRR12(); break;  
    case 13: perform_VRR13(); break;  
    default: assert(false); break;
  }
#endif
  perform_VRR();

  stack->release(size_alloc_);
}


void GradBatch::perform_VRR() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = rank_ * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = stack->get(worksize*3);
  double* worky = workx + worksize; 
  double* workz = worky + worksize; 
  double iyiz[rank_];

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
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], rank_, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], rank_, worksize, worky, vrr_->vrrfunc[vrr_index]);
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], rank_, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz);
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * (amax1_ * iz + jz);
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * (amax1_ * iy + jy);
            const int jyz = amax1_ * (jy + amax1_ * jz);
            for (int i = 0; i != rank_; ++i)
              iyiz[i] = worky[offsety + i] * workz[offsetz + i]; 
            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_;
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int offsetx = rank_ * (amax1_ * ix + jx);
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                current_data[ijposition] = 0.0;
                for (int i = 0; i != rank_; ++i)
                  current_data[ijposition + i] += iyiz[i] * workx[offsetx + i]; 
              }
            }
          }
        }
      }
    }

  }

  stack->release(worksize*3);
}
