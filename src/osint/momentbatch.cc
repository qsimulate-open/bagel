//
// BAGEL - Parallel electron correlation program.
// Filename: momentbatch.cc
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


#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/osint/momentbatch.h>
#include <src/rysint/macros.h>
#include <src/rysint/carsphlist.h>

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;

MomentBatch::MomentBatch(const array<std::shared_ptr<const Shell>,2>& _basis, shared_ptr<StackMem> stack)
 : OSInt(_basis, -1, stack) {

}


MomentBatch::~MomentBatch() {
}


void MomentBatch::compute() {

  double* const intermediate_p = stack_->get(prim0_*prim1_*asize_intermediate_*3);
  perform_VRR(intermediate_p);

  for (int i = 0; i != 3; ++i) {
    double* cdata = data_ + i*size_block_;
    const double* csource = intermediate_p + i*prim0_*prim1_*asize_intermediate_;
    double* const intermediate_c = stack_->get(cont0_ * cont1_ * asize_intermediate_);
    fill(intermediate_c, intermediate_c + cont0_ * cont1_ * asize_intermediate_, 0.0);
    perform_contraction(asize_intermediate_, csource, prim0_, prim1_, intermediate_c,
                        basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0_,
                        basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1_);

    if (basisinfo_[0]->spherical() && basisinfo_[1]->spherical()) {
      // transform both indices to spherical
      double* const intermediate_i = stack_->get(cont0_ * cont1_ * asize_final_);
      const unsigned int carsph_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      const int nloops = cont0_ * cont1_;
      carsphlist.carsphfunc_call(carsph_index, nloops, intermediate_c, intermediate_i);

      const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort_.sortfunc_call(sort_index, cdata, intermediate_i, cont1_, cont0_, 1, swap01_);
      stack_->release(cont0_ * cont1_ * asize_final_, intermediate_i);
    } else {
      const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort_.sortfunc_call(sort_index, cdata, intermediate_c, cont1_, cont0_, 1, swap01_);
    }

    stack_->release(cont0_ * cont1_ * asize_intermediate_, intermediate_c);
  }
  stack_->release(prim0_*prim1_*asize_intermediate_*3, intermediate_p);

  if (swap01_) dscal_(size_alloc_*3, -1.0, data_, 1);

}


void MomentBatch::perform_VRR(double* intermediate) {
  const int size_b = prim0_*prim1_*asize_intermediate_;

  const int dim = amax1_+1;
  const int worksize = dim;
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
    const double mbop = - cb * cop;
    const double pabop = (ca + cb) * cop;
    const double cxpa = p_[ii * 3    ] - basisinfo_[0]->position(0);
    const double cypa = p_[ii * 3 + 1] - basisinfo_[0]->position(1);
    const double czpa = p_[ii * 3 + 2] - basisinfo_[0]->position(2);
    double* current_data = &intermediate[offset_ii];
    worksx[0] = coeffsx_[ii];
    worksy[0] = coeffsy_[ii];
    worksz[0] = coeffsz_[ii];

    for (int i = 1; i != dim; ++i) {
      worksx[i] = cxpa * worksx[i-1] + 0.5 * (i-1) * cop * worksx[i-2];
      worksy[i] = cypa * worksy[i-1] + 0.5 * (i-1) * cop * worksy[i-2];
      worksz[i] = czpa * worksz[i-1] + 0.5 * (i-1) * cop * worksz[i-2];
    }

    // peform HRR to obtain S(1, j)
    for (int j = 1; j <= ang1_; ++j) {
      for (int i = 0; i != dim - j; ++i) {
        worksx[j * dim + i] = AB_[0] * worksx[(j-1) * dim + i] + worksx[(j-1) * dim + i + 1];
        worksy[j * dim + i] = AB_[1] * worksy[(j-1) * dim + i] + worksy[(j-1) * dim + i + 1];
        worksz[j * dim + i] = AB_[2] * worksz[(j-1) * dim + i] + worksz[(j-1) * dim + i + 1];
      }
    }

    for (int j = 0; j <= ang1_; ++j) {
      for (int i = 0; i <= ang0_; ++i) {
        worktx[j * dim + i] = 2.0 * ca * worksx[j * dim + i+1] - i * worksx[j * dim + i-1];
        workty[j * dim + i] = 2.0 * ca * worksy[j * dim + i+1] - i * worksy[j * dim + i-1];
        worktz[j * dim + i] = 2.0 * ca * worksz[j * dim + i+1] - i * worksz[j * dim + i-1];
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
                current_data[cnt         ] = worktx[ix + dim * jx] * worksy[iy + dim * jy] * worksz[iz + dim * jz];
                current_data[cnt+size_b  ] = worksx[ix + dim * jx] * workty[iy + dim * jy] * worksz[iz + dim * jz];
                current_data[cnt+size_b*2] = worksx[ix + dim * jx] * worksy[iy + dim * jy] * worktz[iz + dim * jz];
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
