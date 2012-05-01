//
// Newint - Parallel electron correlation program.
// Filename: vrr_template.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <src/rysint/eribatch.h>
#include <iostream>
#include <iomanip>
#include <src/stackmem.h>

using namespace std;

extern StackMem* stack;

void ERIBatch::perform_VRR() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = rank_ * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_; 

  const int acsize = asize_ * csize_;

  double* workx = stack->get(worksize*3);
  double* worky = workx + worksize;
  double* workz = worky + worksize;
  double iyiz[RYS_MAX];

  // Perform VRR
  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data = &data_[data_offset_ii];

    /// workx, worky, workz would be Ix, Iy, and Iz data
    const int ii3 = 3 * ii;
    const double dparamx[8] = {p_[ii3], q_[ii3], basisinfo_[0]->position(0), basisinfo_[1]->position(0), 
                                                 basisinfo_[2]->position(0), basisinfo_[3]->position(0), xp_[ii], xq_[ii]};
    Int2D cix(dparamx, &roots_[offset], rank_, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);
  
    const double dparamy[8] = {p_[ii3 + 1], q_[ii3 + 1], basisinfo_[0]->position(1), basisinfo_[1]->position(1), 
                                                         basisinfo_[2]->position(1), basisinfo_[3]->position(1), xp_[ii], xq_[ii] };
    Int2D ciy(dparamy, &roots_[offset], rank_, worksize, worky, vrr_->vrrfunc[vrr_index]);
  
    const double dparamz[8] = {p_[ii3 + 2], q_[ii3 + 2], basisinfo_[0]->position(2), basisinfo_[1]->position(2), 
                                                         basisinfo_[2]->position(2), basisinfo_[3]->position(2), xp_[ii], xq_[ii] };
    Int2D ciz(dparamz, &roots_[offset], rank_, worksize, workz, vrr_->vrrfunc[vrr_index]);

    /// assembly process

    for (int iz = 0; iz <= cmax_; ++iz) { 
      for (int iy = 0; iy <= cmax_ - iz; ++iy) { 
        const int iyz = cmax1_ * (iy + cmax1_ * iz); 
        for (int jz = 0; jz <= amax_; ++jz) { 
          const int offsetz = rank_ * ((amax_ + 1) * iz + jz); 
          for (int jy = 0; jy <= amax_ - jz; ++jy) { 
            const int offsety = rank_ * ((amax_ + 1) * iy + jy); 
            const int jyz = amax1_ * (jy + amax1_ * jz); 
            for (int i = 0; i != rank_; ++i) iyiz[i] = worky[offsety + i] * workz[offsetz + i]; 

            for (int ix = max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) { 
              const int iposition = cmapping_[ix + iyz];
              const int ipos_asize = iposition * asize_; 
              for (int jx = max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) { 
                const int jposition = amapping_[jx + jyz];
                const int ijposition = jposition + ipos_asize;

                const int offsetx = rank_ * ((amax_ + 1) * ix + jx);
                current_data[ijposition] = 0.0;
                for (int i = 0; i != rank_; ++i) {
                  current_data[ijposition] += iyiz[i] * workx[offsetx + i]; 
                } 
              }
            }
          }
        }
      }
    }

  } // end of primsize loop

  stack->release(worksize*3);
}



