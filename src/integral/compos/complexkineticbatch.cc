//
// BAGEL - Parallel electron correlation program.
// Filename: complexkineticbatch.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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


#include <complex>
#include <src/integral/carsphlist.h>
#include <src/integral/compos/complexkineticbatch.h>

using namespace std;
using namespace bagel;

const static CCarSphList carsphlist;


void ComplexKineticBatch::compute() {

  const CSortList sort_ (spherical_);

  complex<double>* const intermediate_p = stack_->get<complex<double>>(prim0_ * prim1_ * asize_intermediate_);
  perform_VRR(intermediate_p);

  complex<double>* const intermediate_c = stack_->get<complex<double>>(cont0_ * cont1_ * asize_intermediate_);
  perform_contraction(asize_intermediate_, intermediate_p, prim0_, prim1_, intermediate_c,
                      basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0_,
                      basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1_);

  if (spherical_) {
    complex<double>* const intermediate_i = stack_->get<complex<double>>(cont0_ * cont1_ * asize_final_);
    const unsigned int carsph_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = cont0_ * cont1_;
    carsphlist.carsphfunc_call(carsph_index, nloops, intermediate_c, intermediate_i);

    const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort_.sortfunc_call(sort_index, data_, intermediate_i, cont1_, cont0_, 1, swap01_);
    stack_->release(cont0_ * cont1_ * asize_final_, intermediate_i);
  } else {
    const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort_.sortfunc_call(sort_index, data_, intermediate_c, cont1_, cont0_, 1, swap01_);
  }

  stack_->release(cont0_*cont1_*asize_intermediate_, intermediate_c);
  stack_->release(prim0_*prim1_*asize_intermediate_, intermediate_p);
}


void ComplexKineticBatch::perform_VRR(complex<double>* intermediate) {
  const int amax3 = amax1_+2;
  const int worksize = amax3;
  complex<double>* worksx = stack_->get<complex<double>>(worksize * worksize);
  complex<double>* worksy = stack_->get<complex<double>>(worksize * worksize);
  complex<double>* worksz = stack_->get<complex<double>>(worksize * worksize);

  const complex<double> imag (0.0, 1.0);
  const double Bx = basisinfo_[1]->position(0);
  const double By = basisinfo_[1]->position(1);
  const double Bz = basisinfo_[1]->position(2);
  const double A_Bx = basisinfo_[1]->vector_potential(0);
  const double A_By = basisinfo_[1]->vector_potential(1);
  const double A_Bz = basisinfo_[1]->vector_potential(2);
  const double field_x = basisinfo_[0]->magnetic_field(0);
  const double field_y = basisinfo_[0]->magnetic_field(1);
  const double field_z = basisinfo_[0]->magnetic_field(2);
  assert (basisinfo_[1]->magnetic_field(0)==field_x);
  assert (basisinfo_[1]->magnetic_field(1)==field_y);
  assert (basisinfo_[1]->magnetic_field(2)==field_z);

  for (int ii = 0; ii != prim0_ * prim1_; ++ii) {
    // Perform VRR
    int offset_ii = ii * asize_intermediate_;
    const double cop = 1.0 / xp_[ii];
    const double cb = xb_[ii];
    const complex<double> cxpa = P_[ii * 3    ] - basisinfo_[0]->position(0);
    const complex<double> cypa = P_[ii * 3 + 1] - basisinfo_[0]->position(1);
    const complex<double> czpa = P_[ii * 3 + 2] - basisinfo_[0]->position(2);
    complex<double>* current_data = &intermediate[offset_ii];
    worksx[0] = 0.5 * coeffsx_[ii];
    worksy[0] = coeffsy_[ii];
    worksz[0] = coeffsz_[ii];

    if (ang0_ + ang1_ + 2 > 0) {
      worksx[1] = cxpa * worksx[0];
      worksy[1] = cypa * worksy[0];
      worksz[1] = czpa * worksz[0];

      for (int i = 2; i != amax3; ++i) {
        worksx[i] = cxpa * worksx[i-1] + 0.5 * (i-1) * cop * worksx[i-2];
        worksy[i] = cypa * worksy[i-1] + 0.5 * (i-1) * cop * worksy[i-2];
        worksz[i] = czpa * worksz[i-1] + 0.5 * (i-1) * cop * worksz[i-2];
      }
    }

    // peform HRR to obtain S(1, j)
    if (ang1_ + 2 > 0) {
      // obtain S(0, 1) and T(0, 1)
      {
        const int j = 1;
        worksx[j * amax3] = AB_[0] * worksx[(j-1) * amax3] + worksx[(j-1) * amax3 + 1];
        worksy[j * amax3] = AB_[1] * worksy[(j-1) * amax3] + worksy[(j-1) * amax3 + 1];
        worksz[j * amax3] = AB_[2] * worksz[(j-1) * amax3] + worksz[(j-1) * amax3 + 1];

        for (int i = 1; i != amax3 - j; ++i) {
          worksx[j * amax3 + i] = AB_[0] * worksx[(j-1) * amax3 + i] + worksx[(j-1) * amax3 + i + 1];
          worksy[j * amax3 + i] = AB_[1] * worksy[(j-1) * amax3 + i] + worksy[(j-1) * amax3 + i + 1];
          worksz[j * amax3 + i] = AB_[2] * worksz[(j-1) * amax3 + i] + worksz[(j-1) * amax3 + i + 1];
        }
      }
      for (int j = 2; j <= ang1_ + 2; ++j) {
        worksx[j * amax3] = AB_[0] * worksx[(j-1) * amax3] + worksx[(j-1) * amax3 + 1];
        worksy[j * amax3] = AB_[1] * worksy[(j-1) * amax3] + worksy[(j-1) * amax3 + 1];
        worksz[j * amax3] = AB_[2] * worksz[(j-1) * amax3] + worksz[(j-1) * amax3 + 1];

        for (int i = 1; i != amax3 - j; ++i) {
          worksx[j * amax3 + i] = AB_[0] * worksx[(j-1) * amax3 + i] + worksx[(j-1) * amax3 + i + 1];
          worksy[j * amax3 + i] = AB_[1] * worksy[(j-1) * amax3 + i] + worksy[(j-1) * amax3 + i + 1];
          worksz[j * amax3 + i] = AB_[2] * worksz[(j-1) * amax3 + i] + worksz[(j-1) * amax3 + i + 1];
        }
      }
    }

    complex<double> Sx, Sy, Sz;    // operator = 1
    complex<double> ox, oy, oz;    // operator = x
    complex<double> tx, ty, tz;    // operator = x^2
    complex<double> dx, dy, dz;    // operator = d/dx
    complex<double> Tx, Ty, Tz;    // operator = d^2/dx^2

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
                complex<double> minus1x = 0.0;
                complex<double> minus1y = 0.0;
                complex<double> minus1z = 0.0;
                complex<double> minus2x = 0.0;
                complex<double> minus2y = 0.0;
                complex<double> minus2z = 0.0;
                if (jx > 0) minus1x = jxd * worksx[ix + amax3 * (jx-1)];
                if (jy > 0) minus1y = jyd * worksy[iy + amax3 * (jy-1)];
                if (jz > 0) minus1z = jzd * worksz[iz + amax3 * (jz-1)];
                if (jx > 1) minus2x = jxd * (jxd-1.0) * worksx[ix + amax3 * (jx-2)];
                if (jy > 1) minus2y = jyd * (jyd-1.0) * worksy[iy + amax3 * (jy-2)];
                if (jz > 1) minus2z = jzd * (jzd-1.0) * worksz[iz + amax3 * (jz-2)];

                Sx = worksx[ix + amax3 * jx];
                Sy = worksy[iy + amax3 * jy];
                Sz = worksz[iz + amax3 * jz];

                ox = worksx[ix + amax3 * (jx+1)] + Bx * worksx[ix + amax3 * jx];
                oy = worksy[iy + amax3 * (jy+1)] + By * worksy[iy + amax3 * jy];
                oz = worksz[iz + amax3 * (jz+1)] + Bz * worksz[iz + amax3 * jz];

                tx = worksx[ix + amax3 * (jx+2)] + 2.0 * Bx * worksx[ix + amax3 * (jx+1)] + Bx * Bx * worksx[ix + amax3 * jx];
                ty = worksy[iy + amax3 * (jy+2)] + 2.0 * By * worksy[iy + amax3 * (jy+1)] + By * By * worksy[iy + amax3 * jy];
                tz = worksz[iz + amax3 * (jz+2)] + 2.0 * Bz * worksz[iz + amax3 * (jz+1)] + Bz * Bz * worksz[iz + amax3 * jz];

                dx = minus1x - 2.0 * cb * worksx[ix+amax3 * (jx+1)] - imag * A_Bx * worksx[ix + amax3 * jx];
                dy = minus1y - 2.0 * cb * worksy[iy+amax3 * (jy+1)] - imag * A_By * worksy[iy + amax3 * jy];
                dz = minus1z - 2.0 * cb * worksz[iz+amax3 * (jz+1)] - imag * A_Bz * worksz[iz + amax3 * jz];

                Tx = minus2x - imag * 2.0 * A_Bx * minus1x - (2*cb*(2*jx+1) - A_Bx*A_Bx) * worksx[ix+amax3 * jx]
                     + imag * 4.0 * cb * A_Bx * worksx[ix + amax3 * jx] + 4 * cb * cb * worksx[ix + amax3 * (jx + 2)];
                Ty = minus2y - imag * 2.0 * A_By * minus1y - (2*cb*(2*jy+1) - A_By*A_By) * worksy[iy+amax3 * jy]
                     + imag * 4.0 * cb * A_By * worksy[iy + amax3 * jy] + 4 * cb * cb * worksy[iy + amax3 * (jy + 2)];
                Tz = minus2z - imag * 2.0 * A_Bz * minus1z - (2*cb*(2*jz+1) - A_Bz*A_Bz) * worksz[iz+amax3 * jz]
                     + imag * 4.0 * cb * A_Bz * worksz[iz + amax3 * jz] + 4 * cb * cb * worksz[iz + amax3 * (jz + 2)];

                current_data[cnt] = -Tx * Sy * Sz;
                current_data[cnt] -= Sx * Ty * Sz;
                current_data[cnt] -= Sx * Sy * Tz;
                current_data[cnt] += 0.25 * (field_y * field_y + field_z * field_z) * tx * Sy * Sz;
                current_data[cnt] += 0.25 * (field_z * field_z + field_x * field_x) * Sx * ty * Sz;
                current_data[cnt] += 0.25 * (field_x * field_x + field_y * field_y) * Sx * Sy * tz;
                current_data[cnt] -= 0.5 * field_x * field_y * ox * oy * Sz;
                current_data[cnt] -= 0.5 * field_y * field_z * Sx * oy * oz;
                current_data[cnt] -= 0.5 * field_z * field_x * ox * Sy * oz;
                current_data[cnt] += imag * field_z * dx * oy * Sz;
                current_data[cnt] -= imag * field_y * dx * Sy * oz;
                current_data[cnt] += imag * field_x * Sx * dy * oz;
                current_data[cnt] -= imag * field_z * ox * dy * Sz;
                current_data[cnt] += imag * field_y * ox * Sy * dz;
                current_data[cnt] -= imag * field_x * Sx * oy * dz;
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
