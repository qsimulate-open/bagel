//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: momentumbatch.cc
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
#include <src/integral/os/momentumbatch.h>

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;


void MomentumBatch::compute() {

  const SortList sort(spherical_);

  double* const intermediate_p = stack_->get(size_block_*3);
  perform_VRR(intermediate_p);

  for (int i = 0; i != 3; ++i) {
    double* const cdata = data_ + i*size_block_;
    const double* const csource = intermediate_p + i*size_block_;
    double* const intermediate_c = stack_->get(cont0_ * cont1_ * asize_intermediate_);
    fill_n(intermediate_c, cont0_ * cont1_ * asize_intermediate_, 0.0);
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
      sort.sortfunc_call(sort_index, cdata, intermediate_i, cont1_, cont0_, 1, swap01_);
      stack_->release(cont0_ * cont1_ * asize_final_, intermediate_i);
    } else {
      const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort.sortfunc_call(sort_index, cdata, intermediate_c, cont1_, cont0_, 1, swap01_);
    }

    stack_->release(cont0_ * cont1_ * asize_intermediate_, intermediate_c);
  }
  stack_->release(size_block_*3, intermediate_p);

  if (swap01_) blas::scale_n(-1.0, data_, size_block_*3);

}


void MomentumBatch::perform_VRR(double* intermediate) {
  const int worksize = amax1_;
  double* workpx = stack_->get(worksize*worksize);
  double* workpy = stack_->get(worksize*worksize);
  double* workpz = stack_->get(worksize*worksize);
  double* worksx = stack_->get(worksize*worksize);
  double* worksy = stack_->get(worksize*worksize);
  double* worksz = stack_->get(worksize*worksize);

  for (int ii = 0; ii != prim0_*prim1_; ++ii) {
    // Perform VRR
    int offset_ii = ii*asize_intermediate_;
    const double cop = 1.0 / xp_[ii];
    const double ca = xa_[ii];
    const double cb = xb_[ii];
    const double bop = cb*cop;
    const double cxpa = P_[ii*3  ]-basisinfo_[0]->position(0);
    const double cypa = P_[ii*3+1]-basisinfo_[0]->position(1);
    const double czpa = P_[ii*3+2]-basisinfo_[0]->position(2);
    double* current_data = &intermediate[offset_ii];
    worksx[0] = coeffsx_[ii];
    worksy[0] = coeffsy_[ii];
    worksz[0] = coeffsz_[ii];

    workpx[0] = coeffsx_[ii]*2.0*ca*cxpa;
    workpy[0] = coeffsy_[ii]*2.0*ca*cypa;
    workpz[0] = coeffsz_[ii]*2.0*ca*czpa;

    if (ang0_+ang1_ > 0) {
      worksx[1] = cxpa*worksx[0];
      worksy[1] = cypa*worksy[0];
      worksz[1] = czpa*worksz[0];

      workpx[1] = cxpa*workpx[0]-bop*worksx[0];
      workpy[1] = cypa*workpy[0]-bop*worksy[0];
      workpz[1] = czpa*workpz[0]-bop*worksz[0];

      for (int i = 2; i != amax1_; ++i) {
        worksx[i] = cxpa*worksx[i-1]+0.5*(i-1)*cop*worksx[i-2];
        worksy[i] = cypa*worksy[i-1]+0.5*(i-1)*cop*worksy[i-2];
        worksz[i] = czpa*worksz[i-1]+0.5*(i-1)*cop*worksz[i-2];

        workpx[i] = cxpa*workpx[i-1]+0.5*(i-1)*cop*workpx[i-2]-bop*worksx[i-1];
        workpy[i] = cypa*workpy[i-1]+0.5*(i-1)*cop*workpy[i-2]-bop*worksy[i-1];
        workpz[i] = czpa*workpz[i-1]+0.5*(i-1)*cop*workpz[i-2]-bop*worksz[i-1];
      }
    }

    // peform HRR to obtain S(1, j)
    if (ang1_ > 0) {
      for (int j = 1; j <= ang1_; ++j) {
        for (int i = 0; i != amax1_-j; ++i) {
          worksx[j*amax1_+i] = AB_[0]*worksx[(j-1)*amax1_+i]+worksx[(j-1)*amax1_+i+1];
          worksy[j*amax1_+i] = AB_[1]*worksy[(j-1)*amax1_+i]+worksy[(j-1)*amax1_+i+1];
          worksz[j*amax1_+i] = AB_[2]*worksz[(j-1)*amax1_+i]+worksz[(j-1)*amax1_+i+1];

          workpx[j*amax1_+i] = AB_[0]*workpx[(j-1)*amax1_+i]+workpx[(j-1)*amax1_+i+1]+worksx[(j-1)*amax1_+i];
          workpy[j*amax1_+i] = AB_[1]*workpy[(j-1)*amax1_+i]+workpy[(j-1)*amax1_+i+1]+worksy[(j-1)*amax1_+i];
          workpz[j*amax1_+i] = AB_[2]*workpz[(j-1)*amax1_+i]+workpz[(j-1)*amax1_+i+1]+worksz[(j-1)*amax1_+i];
        }
      }
    }

    // now we obtain the output
    assert((ang0_+1)*(ang0_+2)*(ang1_+1)*(ang1_+2)/4 == asize_intermediate_);

    int cnt = 0;
    for (int iz = 0; iz <= ang0_; ++iz) {
      for (int iy = 0; iy <= ang0_-iz; ++iy) {
        const int ix = ang0_-iy-iz;
        if (ix >= 0) {
          for (int jz = 0; jz <= ang1_; ++jz) {
            for (int jy = 0; jy <= ang1_-jz; ++jy) {
              const int jx = ang1_-jy-jz;
              if (jx >= 0) {
                current_data[cnt              ] = workpx[ix+amax1_*jx]*worksy[iy+amax1_*jy]*worksz[iz+amax1_*jz];
                current_data[cnt+size_block_  ] = worksx[ix+amax1_*jx]*workpy[iy+amax1_*jy]*worksz[iz+amax1_*jz];
                current_data[cnt+size_block_*2] = worksx[ix+amax1_*jx]*worksy[iy+amax1_*jy]*workpz[iz+amax1_*jz];
                ++cnt;
              }
            }
          }
        }
      }
    }
    assert(cnt == asize_intermediate_);

  } // end of prim exponent loop

  stack_->release(worksize*worksize, worksz);
  stack_->release(worksize*worksize, worksy);
  stack_->release(worksize*worksize, worksx);
  stack_->release(worksize*worksize, workpz);
  stack_->release(worksize*worksize, workpy);
  stack_->release(worksize*worksize, workpx);
}
