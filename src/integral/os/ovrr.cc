//
// BAGEL - Parallel electron correlation program.
// Filename: ovrr.cc
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


#include <iostream>
#include <iomanip>
#include <src/integral/os/overlapbatch.h>

using namespace std;
using namespace bagel;

// private functions
void OverlapBatch::perform_VRR(double* intermediate) {

  const int worksize = amax1_;
  double* workx = stack_->get(worksize);
  double* worky = stack_->get(worksize);
  double* workz = stack_->get(worksize);

  // Perform VRR
  for (int ii = 0; ii != prim0_ * prim1_; ++ii) {
    const int offset_ii = ii * asize_;
    double* current_data = &intermediate[offset_ii];

    /// Sx(0 : i + j, 0) etc will be made here
    workx[0] = coeffsx_[ii];
    worky[0] = coeffsy_[ii];
    workz[0] = coeffsz_[ii];
    if (ang0_ + ang1_ > 0) {
      workx[1] = (p_[ii * 3    ] - basisinfo_[0]->position(0)) * workx[0];
      worky[1] = (p_[ii * 3 + 1] - basisinfo_[0]->position(1)) * worky[0];
      workz[1] = (p_[ii * 3 + 2] - basisinfo_[0]->position(2)) * workz[0];
      for (int i = 2; i != amax1_; ++i) {
        workx[i] = (p_[ii * 3    ] - basisinfo_[0]->position(0)) * workx[i - 1] + 0.5 * (i - 1) / xp_[ii] * workx[i - 2];
        worky[i] = (p_[ii * 3 + 1] - basisinfo_[0]->position(1)) * worky[i - 1] + 0.5 * (i - 1) / xp_[ii] * worky[i - 2];
        workz[i] = (p_[ii * 3 + 2] - basisinfo_[0]->position(2)) * workz[i - 1] + 0.5 * (i - 1) / xp_[ii] * workz[i - 2];
      }
    }

    /// assembly process

    for (int iz = 0; iz <= amax_; ++iz) {
      for (int iy = 0; iy <= amax_ - iz; ++iy) {
        const double iyiz = workz[iz] * worky[iy];
        for (int ix = max(0, amin_ - iy - iz); ix <= amax_ - iy - iz; ++ix) {
          int pos = amapping_[ix + amax1_ * (iy + amax1_ * iz)];
          current_data[pos] = workx[ix] * iyiz;
        }
      }
    }

  } // end of primsize loop

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);
}

