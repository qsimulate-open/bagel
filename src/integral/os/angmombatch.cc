//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: angmombatch.cc
// Copyright (C) 2015 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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
#include <src/integral/os/angmombatch.h>

using namespace std;
using namespace bagel;

const static CarSphList carsphlist;


void AngMomBatch::compute() {

  const SortList sort_(spherical_);

  double* const intermediate_p = stack_->get<double>(size_block_*3);
  fill_n(intermediate_p, size_block_*3, 0.0);
  perform_VRR(intermediate_p);

  for (int i = 0; i != 3; ++i) {
    double* const cdata = data_ + i*size_block_;
    const double* const csource = intermediate_p + i*size_block_;
    double* const intermediate_c = stack_->get<double>(cont0_ * cont1_ * asize_intermediate_);
    fill_n(intermediate_c, cont0_ * cont1_ * asize_intermediate_, 0.0);
    perform_contraction(asize_intermediate_, csource, prim0_, prim1_, intermediate_c,
                        basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0_,
                        basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1_);

    if (basisinfo_[0]->spherical() && basisinfo_[1]->spherical()) {
      // transform both indices to spherical
      double* const intermediate_i = stack_->get<double>(cont0_ * cont1_ * asize_final_);
      fill_n(intermediate_i, cont0_ * cont1_ * asize_final_, 0.0);
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
  stack_->release(size_block_*3, intermediate_p);

}


void AngMomBatch::perform_VRR(double* intermediate) {
  const int amax3 = amax1_+2;
  const int worksize = amax3;
  double* worksx = stack_->get<double>(worksize * worksize);
  double* worksy = stack_->get<double>(worksize * worksize);
  double* worksz = stack_->get<double>(worksize * worksize);

  const double imag = swap01_ ? -1.0 : 1.0;
  const double Bx = basisinfo_[1]->position(0);
  const double By = basisinfo_[1]->position(1);
  const double Bz = basisinfo_[1]->position(2);

  for (int ii = 0; ii != prim0_ * prim1_; ++ii) {
    // Perform VRR
    int offset_ii = ii * asize_intermediate_;
    const double cop = 1.0 / xp_[ii];
    const double cb = xb_[ii];
    const double cxpa = P_[ii * 3    ] - basisinfo_[0]->position(0);
    const double cypa = P_[ii * 3 + 1] - basisinfo_[0]->position(1);
    const double czpa = P_[ii * 3 + 2] - basisinfo_[0]->position(2);
    double* current_data = &intermediate[offset_ii];

    // obtain S(0, 0)
    worksx[0] = coeffsx_[ii];
    worksy[0] = coeffsy_[ii];
    worksz[0] = coeffsz_[ii];

    // obtain S(1, 0)
    worksx[1] = cxpa * worksx[0];
    worksy[1] = cypa * worksy[0];
    worksz[1] = czpa * worksz[0];

    for (int i = 2; i != amax3; ++i) {
      // obtain S(i, 0)
      worksx[i] = cxpa * worksx[i-1] + 0.5 * (i-1) * cop * worksx[i-2];
      worksy[i] = cypa * worksy[i-1] + 0.5 * (i-1) * cop * worksy[i-2];
      worksz[i] = czpa * worksz[i-1] + 0.5 * (i-1) * cop * worksz[i-2];
    }

    // peform HRR
    {
      const int j = 1;
      // obtain S(0, 1)
      worksx[j * amax3] = AB_[0] * worksx[(j-1) * amax3] + worksx[(j-1) * amax3 + 1];
      worksy[j * amax3] = AB_[1] * worksy[(j-1) * amax3] + worksy[(j-1) * amax3 + 1];
      worksz[j * amax3] = AB_[2] * worksz[(j-1) * amax3] + worksz[(j-1) * amax3 + 1];

      for (int i = 1; i != amax3 - j; ++i) {
        // obtain S(i, 1)
        worksx[j * amax3 + i] = AB_[0] * worksx[(j-1) * amax3 + i] + worksx[(j-1) * amax3 + i + 1];
        worksy[j * amax3 + i] = AB_[1] * worksy[(j-1) * amax3 + i] + worksy[(j-1) * amax3 + i + 1];
        worksz[j * amax3 + i] = AB_[2] * worksz[(j-1) * amax3 + i] + worksz[(j-1) * amax3 + i + 1];
      }
    }
    for (int j = 2; j <= ang1_ + 2; ++j) {
      // obtain S(0, j)
      worksx[j * amax3] = AB_[0] * worksx[(j-1) * amax3] + worksx[(j-1) * amax3 + 1];
      worksy[j * amax3] = AB_[1] * worksy[(j-1) * amax3] + worksy[(j-1) * amax3 + 1];
      worksz[j * amax3] = AB_[2] * worksz[(j-1) * amax3] + worksz[(j-1) * amax3 + 1];

      for (int i = 1; i != amax3 - j; ++i) {
        // obtain S(i, j)
        worksx[j * amax3 + i] = AB_[0] * worksx[(j-1) * amax3 + i] + worksx[(j-1) * amax3 + i + 1];
        worksy[j * amax3 + i] = AB_[1] * worksy[(j-1) * amax3 + i] + worksy[(j-1) * amax3 + i + 1];
        worksz[j * amax3 + i] = AB_[2] * worksz[(j-1) * amax3 + i] + worksz[(j-1) * amax3 + i + 1];
      }
    }

    double Sx, Sy, Sz;    // operator = 1
    double ox, oy, oz;    // operator = x
    double dx, dy, dz;    // operator = d/dx

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

                const double jxd = jx;
                const double jyd = jy;
                const double jzd = jz;
                double minus1x = 0.0;
                double minus1y = 0.0;
                double minus1z = 0.0;
                if (jx > 0) minus1x = jxd * worksx[ix + amax3 * (jx-1)];
                if (jy > 0) minus1y = jyd * worksy[iy + amax3 * (jy-1)];
                if (jz > 0) minus1z = jzd * worksz[iz + amax3 * (jz-1)];

                Sx = worksx[ix + amax3 * jx];
                Sy = worksy[iy + amax3 * jy];
                Sz = worksz[iz + amax3 * jz];

                ox = worksx[ix + amax3 * (jx+1)] + Bx * worksx[ix + amax3 * jx];
                oy = worksy[iy + amax3 * (jy+1)] + By * worksy[iy + amax3 * jy];
                oz = worksz[iz + amax3 * (jz+1)] + Bz * worksz[iz + amax3 * jz];

                dx = minus1x - 2.0 * cb * worksx[ix+amax3 * (jx+1)];
                dy = minus1y - 2.0 * cb * worksy[iy+amax3 * (jy+1)];
                dz = minus1z - 2.0 * cb * worksz[iz+amax3 * (jz+1)];

                current_data[cnt              ] = -imag * mcoord_[2] * Sx * dy * Sz;
                current_data[cnt+size_block_  ] = -imag * mcoord_[0] * Sx * Sy * dz;
                current_data[cnt+size_block_*2] = -imag * mcoord_[1] * dx * Sy * Sz;

                current_data[cnt              ] += imag * mcoord_[1] * Sx * Sy * dz;
                current_data[cnt+size_block_  ] += imag * mcoord_[2] * dx * Sy * Sz;
                current_data[cnt+size_block_*2] += imag * mcoord_[0] * Sx * dy * Sz;

                current_data[cnt              ] -= imag * Sx * oy * dz;
                current_data[cnt+size_block_  ] -= imag * dx * Sy * oz;
                current_data[cnt+size_block_*2] -= imag * ox * dy * Sz;

                current_data[cnt              ] += imag * Sx * dy * oz;
                current_data[cnt+size_block_  ] += imag * ox * Sy * dz;
                current_data[cnt+size_block_*2] += imag * dx * oy * Sz;

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
}

