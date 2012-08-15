//
// Newint - Parallel electron correlation program.
// Filename: momentbatch.cc
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


#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/stackmem.h>
#include <src/osint/momentbatch.h>
#include <src/rysint/macros.h>
#include <src/rysint/carsphlist.h>

using namespace std;

extern StackMem* stack;

MomentBatch::MomentBatch(const vector<std::shared_ptr<const Shell> >& _basis) 
 : OSInt(_basis, -1) {

}  


MomentBatch::~MomentBatch() {
}


void MomentBatch::compute() {
  double* stack_save = stack->get(0);

  double* intermediate_p = stack->get(prim0_*prim1_*asize_intermediate_*3);
  perform_VRR(intermediate_p);

  for (int i = 0; i != 3; ++i) {
    double* cdata = data_ + i*size_block_;
    const double* csource = intermediate_p + i*prim0_*prim1_*asize_intermediate_;
    double* intermediate_c = stack->get(cont0_ * cont1_ * asize_intermediate_);
    fill(intermediate_c, intermediate_c + cont0_ * cont1_ * asize_intermediate_, 0.0);
    perform_contraction(asize_intermediate_, csource, prim0_, prim1_, intermediate_c, 
                        basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0_, 
                        basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1_);

    if (spherical_) {
      struct CarSphList carsphlist;
      double* intermediate_i = stack->get(cont0_ * cont1_ * asize_final_);
      const unsigned int carsph_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      const int nloops = cont0_ * cont1_;
      carsphlist.carsphfunc_call(carsph_index, nloops, intermediate_c, intermediate_i); 

      const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort_.sortfunc_call(sort_index, cdata, intermediate_i, cont1_, cont0_, 1, swap01_);
      stack->release(cont0_ * cont1_ * asize_final_);
    } else {
      const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort_.sortfunc_call(sort_index, cdata, intermediate_c, cont1_, cont0_, 1, swap01_);
    }

    stack->release(cont0_ * cont1_ * asize_intermediate_);
  }
  stack->release(prim0_*prim1_*asize_intermediate_*3);
  assert(stack->get(0) == stack_save);

  if (swap01_) dscal_(size_alloc_, -1.0, data_, 1); 

}


void MomentBatch::perform_VRR(double* intermediate) {
  const int worksize = amax1_;
  double* worktx = stack->get(worksize * worksize);
  double* workty = stack->get(worksize * worksize);
  double* worktz = stack->get(worksize * worksize);
  double* worksx = stack->get(worksize * worksize);
  double* worksy = stack->get(worksize * worksize);
  double* worksz = stack->get(worksize * worksize);

  for (int ii = 0; ii != prim0_ * prim1_; ++ii) {
    // Perform VRR
    int offset_ii = ii * asize_intermediate_;
    const double cop = 1.0 / xp_[ii];
    const double ca = xa_[ii];
    const double cb = xb_[ii];
    const double mbop = - cb * cop;
    const double pabop = (ca + cb) * cop;
    const double cxpa = p_[ii * 3    ] - basisinfo_[0]->position(0); 
    const double cypa = p_[ii * 3 + 1] - basisinfo_[0]->position(1); 
    const double czpa = p_[ii * 3 + 2] - basisinfo_[0]->position(2); 
    double* current_data = &intermediate[offset_ii];
    worksx[0] = coeffsx_[ii]; 
    worksy[0] = coeffsy_[ii];
    worksz[0] = coeffsz_[ii];

    worktx[0] = 2*ca*cxpa*worksx[0];
    workty[0] = 2*ca*cypa*worksy[0];
    worktz[0] = 2*ca*czpa*worksz[0];

    if (ang0_ + ang1_ > 0) {
      worksx[1] = cxpa * worksx[0];
      worksy[1] = cypa * worksy[0];
      worksz[1] = czpa * worksz[0];

      worktx[1] = cxpa * worktx[0] + mbop * worksx[0];
      workty[1] = cypa * workty[0] + mbop * worksy[0];
      worktz[1] = czpa * worktz[0] + mbop * worksz[0];

      for (int i = 2; i != amax1_; ++i) {
        worksx[i] = cxpa * worksx[i-1] + 0.5 * (i-1) * cop * worksx[i-2];
        worksy[i] = cypa * worksy[i-1] + 0.5 * (i-1) * cop * worksy[i-2];
        worksz[i] = czpa * worksz[i-1] + 0.5 * (i-1) * cop * worksz[i-2];

        worktx[i] = cxpa * worktx[i-1] + 0.5 * (i-1) * cop * worktx[i-2] + mbop * worksx[i-1];
        workty[i] = cypa * workty[i-1] + 0.5 * (i-1) * cop * workty[i-2] + mbop * worksy[i-1];
        worktz[i] = czpa * worktz[i-1] + 0.5 * (i-1) * cop * worktz[i-2] + mbop * worksz[i-1];
      }
    }

    // peform HRR to obtain S(1, j)
    if (ang1_ > 0) {
      for (int j = 1; j <= ang1_; ++j) { 
        for (int i = 0; i != amax1_ - j; ++i) {
          worksx[j * amax1_ + i] = AB_[0] * worksx[(j-1) * amax1_ + i] + worksx[(j-1) * amax1_ + i + 1];
          worksy[j * amax1_ + i] = AB_[1] * worksy[(j-1) * amax1_ + i] + worksy[(j-1) * amax1_ + i + 1];
          worksz[j * amax1_ + i] = AB_[2] * worksz[(j-1) * amax1_ + i] + worksz[(j-1) * amax1_ + i + 1];
 
          worktx[j * amax1_ + i] = AB_[0] * worktx[(j-1) * amax1_ + i] + worktx[(j-1) * amax1_ + i + 1] + pabop * worksx[(j-1) * amax1_ + i]; 
          workty[j * amax1_ + i] = AB_[1] * workty[(j-1) * amax1_ + i] + workty[(j-1) * amax1_ + i + 1] + pabop * worksy[(j-1) * amax1_ + i];
          worktz[j * amax1_ + i] = AB_[2] * worktz[(j-1) * amax1_ + i] + worktz[(j-1) * amax1_ + i + 1] + pabop * worksz[(j-1) * amax1_ + i];
        }
      }
    } 

    // now we obtain the output
    const int isize = (ang0_ + 1) * (ang0_ + 2) / 2;
    const int jsize = (ang1_ + 1) * (ang1_ + 2) / 2;
    assert(isize * jsize == asize_intermediate_);

    int cnt = 0;
    for (int iz = 0; iz <= ang0_; ++iz) { 
      for (int iy = 0; iy <= ang0_ - iz; ++iy) { 
        const int ix = ang0_ - iy - iz;
        if (ix >= 0) {
          for (int jz = 0; jz <= ang1_; ++jz) {
            for (int jy = 0; jy <= ang1_ - jz; ++jy) {
              const int jx = ang1_ - jy - jz;
              if (jx >= 0) {
                current_data[cnt              ] = worktx[ix + amax1_ * jx] * worksy[iy + amax1_ * jy] * worksz[iz + amax1_ * jz]; 
                current_data[cnt+size_block_  ] = worksx[ix + amax1_ * jx] * workty[iy + amax1_ * jy] * worksz[iz + amax1_ * jz]; 
                current_data[cnt+size_block_*2] = worksx[ix + amax1_ * jx] * worksy[iy + amax1_ * jy] * worktz[iz + amax1_ * jz]; 
                ++cnt;
              }
            }
          }
        }
      }
    }
    assert(cnt == asize_intermediate_);

  } // end of prim exponent loop

  stack->release(worksize * worksize * 6);
}
