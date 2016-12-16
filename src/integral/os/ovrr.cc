//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ovrr.cc
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
      workx[1] = (P_[ii * 3    ] - basisinfo_[0]->position(0)) * workx[0];
      worky[1] = (P_[ii * 3 + 1] - basisinfo_[0]->position(1)) * worky[0];
      workz[1] = (P_[ii * 3 + 2] - basisinfo_[0]->position(2)) * workz[0];
      for (int i = 2; i != amax1_; ++i) {
        workx[i] = (P_[ii * 3    ] - basisinfo_[0]->position(0)) * workx[i - 1] + 0.5 * (i - 1) / xp_[ii] * workx[i - 2];
        worky[i] = (P_[ii * 3 + 1] - basisinfo_[0]->position(1)) * worky[i - 1] + 0.5 * (i - 1) / xp_[ii] * worky[i - 2];
        workz[i] = (P_[ii * 3 + 2] - basisinfo_[0]->position(2)) * workz[i - 1] + 0.5 * (i - 1) / xp_[ii] * workz[i - 2];
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

