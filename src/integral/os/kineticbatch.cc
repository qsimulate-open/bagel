//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: kineticbatch.cc
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


#include <src/integral/carsphlist.h>
#include <src/integral/os/kineticbatch.h>

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;


void KineticBatch::compute() {

  double* const intermediate_p = stack_->get(prim0_ * prim1_ * asize_intermediate_);
  perform_VRR(intermediate_p);

  double* const intermediate_c = stack_->get(cont0_ * cont1_ * asize_intermediate_);
  perform_contraction(asize_intermediate_, intermediate_p, prim0_, prim1_, intermediate_c,
                      basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0_,
                      basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1_);

  const SortList sort(spherical_);
  if (spherical_) {
    double* const intermediate_i = stack_->get(cont0_ * cont1_ * asize_final_);
    const unsigned int carsph_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = cont0_ * cont1_;
    carsphlist.carsphfunc_call(carsph_index, nloops, intermediate_c, intermediate_i);

    const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort.sortfunc_call(sort_index, data_, intermediate_i, cont1_, cont0_, 1, swap01_);
    stack_->release(cont0_ * cont1_ * asize_final_, intermediate_i);
  } else {
    const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort.sortfunc_call(sort_index, data_, intermediate_c, cont1_, cont0_, 1, swap01_);
  }

  stack_->release(cont0_*cont1_*asize_intermediate_, intermediate_c);
  stack_->release(prim0_*prim1_*asize_intermediate_, intermediate_p);
}


void KineticBatch::perform_VRR(double* intermediate) {
  const int worksize = amax1_;
  double* worktx = stack_->get(worksize * worksize);
  double* workty = stack_->get(worksize * worksize);
  double* worktz = stack_->get(worksize * worksize);
  double* worksx = stack_->get(worksize * worksize);
  double* worksy = stack_->get(worksize * worksize);
  double* worksz = stack_->get(worksize * worksize);

  for (int ii = 0; ii != prim0_ * prim1_; ++ii) {
    // Perform VRR
    int offset_ii = ii * asize_intermediate_;
    const double cop = 1.0 / xp_[ii];
    const double ca = xa_[ii];
    const double cb = xb_[ii];
    const double tabop = 2.0 * ca * cb * cop;
    const double cxpa = P_[ii * 3    ] - basisinfo_[0]->position(0);
    const double cypa = P_[ii * 3 + 1] - basisinfo_[0]->position(1);
    const double czpa = P_[ii * 3 + 2] - basisinfo_[0]->position(2);
    double* current_data = &intermediate[offset_ii];
    worksx[0] = coeffsx_[ii];
    worksy[0] = coeffsy_[ii];
    worksz[0] = coeffsz_[ii];

    worktx[0] = coefftx_[ii];
    workty[0] = coeffty_[ii];
    worktz[0] = coefftz_[ii];
    if (ang0_ + ang1_ > 0) {
      worksx[1] = cxpa * worksx[0];
      worksy[1] = cypa * worksy[0];
      worksz[1] = czpa * worksz[0];

      worktx[1] = cxpa * worktx[0] + tabop * worksx[1];
      workty[1] = cypa * workty[0] + tabop * worksy[1];
      worktz[1] = czpa * worktz[0] + tabop * worksz[1];

      for (int i = 2; i != amax1_; ++i) {
        worksx[i] = cxpa * worksx[i-1] + 0.5 * (i-1) * cop * worksx[i-2];
        worksy[i] = cypa * worksy[i-1] + 0.5 * (i-1) * cop * worksy[i-2];
        worksz[i] = czpa * worksz[i-1] + 0.5 * (i-1) * cop * worksz[i-2];

        worktx[i] = cxpa * worktx[i-1] + 0.5 * (i-1) * cop * worktx[i-2] + tabop * worksx[i] - cb * cop * (i-1) * worksx[i-2];
        workty[i] = cypa * workty[i-1] + 0.5 * (i-1) * cop * workty[i-2] + tabop * worksy[i] - cb * cop * (i-1) * worksy[i-2];
        worktz[i] = czpa * worktz[i-1] + 0.5 * (i-1) * cop * worktz[i-2] + tabop * worksz[i] - cb * cop * (i-1) * worksz[i-2];
      }
    }

    // peform HRR to obtain S(1, j)
    if (ang1_ > 0) {
      // obtain S(0, 1) and T(0, 1)
      {
        const int j = 1;
        worksx[j * amax1_] = AB_[0] * worksx[(j-1) * amax1_] + worksx[(j-1) * amax1_ + 1];
        worksy[j * amax1_] = AB_[1] * worksy[(j-1) * amax1_] + worksy[(j-1) * amax1_ + 1];
        worksz[j * amax1_] = AB_[2] * worksz[(j-1) * amax1_] + worksz[(j-1) * amax1_ + 1];

        worktx[j * amax1_] = AB_[0] * worktx[(j-1) * amax1_] + worktx[(j-1) * amax1_ + 1] + tabop * (worksx[j * amax1_] - worksx[(j-1) * amax1_ + 1]);
        workty[j * amax1_] = AB_[1] * workty[(j-1) * amax1_] + workty[(j-1) * amax1_ + 1] + tabop * (worksy[j * amax1_] - worksy[(j-1) * amax1_ + 1]);
        worktz[j * amax1_] = AB_[2] * worktz[(j-1) * amax1_] + worktz[(j-1) * amax1_ + 1] + tabop * (worksz[j * amax1_] - worksz[(j-1) * amax1_ + 1]);

        for (int i = 1; i != amax1_ - j; ++i) {
          worksx[j * amax1_ + i] = AB_[0] * worksx[(j-1) * amax1_ + i] + worksx[(j-1) * amax1_ + i + 1];
          worksy[j * amax1_ + i] = AB_[1] * worksy[(j-1) * amax1_ + i] + worksy[(j-1) * amax1_ + i + 1];
          worksz[j * amax1_ + i] = AB_[2] * worksz[(j-1) * amax1_ + i] + worksz[(j-1) * amax1_ + i + 1];

          worktx[j * amax1_ + i] = AB_[0] * worktx[(j-1) * amax1_ + i] + worktx[(j-1) * amax1_ + i + 1]
                              + tabop * (worksx[j * amax1_ + i] - worksx[(j-1) * amax1_ + i + 1])
                              - cop * (- i * cb * worksx[(j-1) * amax1_ + i-1]);
          workty[j * amax1_ + i] = AB_[1] * workty[(j-1) * amax1_ + i] + workty[(j-1) * amax1_ + i + 1]
                              + tabop * (worksy[j * amax1_ + i] - worksy[(j-1) * amax1_ + i + 1])
                              - cop * (- i * cb * worksy[(j-1) * amax1_ + i-1]);
          worktz[j * amax1_ + i] = AB_[2] * worktz[(j-1) * amax1_ + i] + worktz[(j-1) * amax1_ + i + 1]
                              + tabop * (worksz[j * amax1_ + i] - worksz[(j-1) * amax1_ + i + 1])
                              - cop * (- i * cb * worksz[(j-1) * amax1_ + i-1]);
        }
      }
      for (int j = 2; j <= ang1_; ++j) {
        worksx[j * amax1_] = AB_[0] * worksx[(j-1) * amax1_] + worksx[(j-1) * amax1_ + 1];
        worksy[j * amax1_] = AB_[1] * worksy[(j-1) * amax1_] + worksy[(j-1) * amax1_ + 1];
        worksz[j * amax1_] = AB_[2] * worksz[(j-1) * amax1_] + worksz[(j-1) * amax1_ + 1];

        worktx[j * amax1_] = AB_[0] * worktx[(j-1) * amax1_] + worktx[(j-1) * amax1_ + 1] + tabop * (worksx[j * amax1_] - worksx[(j-1) * amax1_ + 1])
                              - cop * ((j-1) * ca * worksx[(j-2) * amax1_]);
        workty[j * amax1_] = AB_[1] * workty[(j-1) * amax1_] + workty[(j-1) * amax1_ + 1] + tabop * (worksy[j * amax1_] - worksy[(j-1) * amax1_ + 1])
                              - cop * ((j-1) * ca * worksy[(j-2) * amax1_]);
        worktz[j * amax1_] = AB_[2] * worktz[(j-1) * amax1_] + worktz[(j-1) * amax1_ + 1] + tabop * (worksz[j * amax1_] - worksz[(j-1) * amax1_ + 1])
                              - cop * ((j-1) * ca * worksz[(j-2) * amax1_]);
        for (int i = 1; i != amax1_ - j; ++i) {
          worksx[j * amax1_ + i] = AB_[0] * worksx[(j-1) * amax1_ + i] + worksx[(j-1) * amax1_ + i + 1];
          worksy[j * amax1_ + i] = AB_[1] * worksy[(j-1) * amax1_ + i] + worksy[(j-1) * amax1_ + i + 1];
          worksz[j * amax1_ + i] = AB_[2] * worksz[(j-1) * amax1_ + i] + worksz[(j-1) * amax1_ + i + 1];

          worktx[j * amax1_ + i] = AB_[0] * worktx[(j-1) * amax1_ + i] + worktx[(j-1) * amax1_ + i + 1]
                              + tabop * (worksx[j * amax1_ + i] - worksx[(j-1) * amax1_ + i + 1])
                              - cop * ((j-1) * ca * worksx[(j-2) * amax1_ + i] - i * cb * worksx[(j-1) * amax1_ + i-1]);
          workty[j * amax1_ + i] = AB_[1] * workty[(j-1) * amax1_ + i] + workty[(j-1) * amax1_ + i + 1]
                              + tabop * (worksy[j * amax1_ + i] - worksy[(j-1) * amax1_ + i + 1])
                              - cop * ((j-1) * ca * worksy[(j-2) * amax1_ + i] - i * cb * worksy[(j-1) * amax1_ + i-1]);
          worktz[j * amax1_ + i] = AB_[2] * worktz[(j-1) * amax1_ + i] + worktz[(j-1) * amax1_ + i + 1]
                              + tabop * (worksz[j * amax1_ + i] - worksz[(j-1) * amax1_ + i + 1])
                              - cop * ((j-1) * ca * worksz[(j-2) * amax1_ + i] - i * cb * worksz[(j-1) * amax1_ + i-1]);
        }
      }
    }

    // now we obtain the output
    int cnt = 0;
    for (int iz = 0; iz <= ang0_; ++iz) {
      for (int iy = 0; iy <= ang0_ - iz; ++iy) {
        const int ix = ang0_ - iy - iz;
        if (ix >= 0) {
          for (int jz = 0; jz <= ang1_; ++jz) {
            for (int jy = 0; jy <= ang1_ - jz; ++jy) {
              const int jx = ang1_ - jy - jz;
              if (jx >= 0) {
                current_data[cnt] =  worktx[ix + amax1_ * jx] * worksy[iy + amax1_ * jy] * worksz[iz + amax1_ * jz];
                current_data[cnt] += worksx[ix + amax1_ * jx] * workty[iy + amax1_ * jy] * worksz[iz + amax1_ * jz];
                current_data[cnt] += worksx[ix + amax1_ * jx] * worksy[iy + amax1_ * jy] * worktz[iz + amax1_ * jz];
                ++cnt;
              }
            }
          }
        }
      }
    }
    assert(cnt == asize_intermediate_);

  } // end of prim exponent loop

  stack_->release(worksize * worksize, worksz);
  stack_->release(worksize * worksize, worksy);
  stack_->release(worksize * worksize, worksx);
  stack_->release(worksize * worksize, worktz);
  stack_->release(worksize * worksize, workty);
  stack_->release(worksize * worksize, worktx);
}
