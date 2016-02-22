//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: complexkineticbatch.cc
// Copyright (C) 2014 Toru Shiozaki
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

  complex<double> imag (0.0, 1.0);
  if (swap01_) imag *= -1.0;
  const double Bx = basisinfo_[1]->position(0);
  const double By = basisinfo_[1]->position(1);
  const double Bz = basisinfo_[1]->position(2);
  const double A_Bx = basisinfo_[1]->vector_potential(0);
  const double A_By = basisinfo_[1]->vector_potential(1);
  const double A_Bz = basisinfo_[1]->vector_potential(2);

  for (int ii = 0; ii != prim0_ * prim1_; ++ii) {
    // Perform VRR
    int offset_ii = ii * asize_intermediate_;
    const double cop = 1.0 / xp_[ii];
    const double cb = xb_[ii];
    const complex<double> cxpa = P_[ii * 3    ] - basisinfo_[0]->position(0);
    const complex<double> cypa = P_[ii * 3 + 1] - basisinfo_[0]->position(1);
    const complex<double> czpa = P_[ii * 3 + 2] - basisinfo_[0]->position(2);
    complex<double>* current_data = &intermediate[offset_ii];

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

                Tx = minus2x - imag * 2.0 * A_Bx * minus1x - (2.0*cb*(2.0*jxd+1.0) + A_Bx*A_Bx) * worksx[ix+amax3 * jx]
                     + imag * 4.0 * cb * A_Bx * worksx[ix + amax3 * (jx+1)] + 4.0 * cb * cb * worksx[ix + amax3 * (jx + 2)];
                Ty = minus2y - imag * 2.0 * A_By * minus1y - (2.0*cb*(2.0*jyd+1.0) + A_By*A_By) * worksy[iy+amax3 * jy]
                     + imag * 4.0 * cb * A_By * worksy[iy + amax3 * (jy+1)] + 4.0 * cb * cb * worksy[iy + amax3 * (jy + 2)];
                Tz = minus2z - imag * 2.0 * A_Bz * minus1z - (2.0*cb*(2.0*jzd+1.0) + A_Bz*A_Bz) * worksz[iz+amax3 * jz]
                     + imag * 4.0 * cb * A_Bz * worksz[iz + amax3 * (jz+1)] + 4.0 * cb * cb * worksz[iz + amax3 * (jz + 2)];

                current_data[cnt] = -Tx * Sy * Sz;
                current_data[cnt] -= Sx * Ty * Sz;
                current_data[cnt] -= Sx * Sy * Tz;
                current_data[cnt] += 0.25 * (magnetic_field_[1] * magnetic_field_[1] + magnetic_field_[2] * magnetic_field_[2]) * tx * Sy * Sz;
                current_data[cnt] += 0.25 * (magnetic_field_[2] * magnetic_field_[2] + magnetic_field_[0] * magnetic_field_[0]) * Sx * ty * Sz;
                current_data[cnt] += 0.25 * (magnetic_field_[0] * magnetic_field_[0] + magnetic_field_[1] * magnetic_field_[1]) * Sx * Sy * tz;
                current_data[cnt] -= 0.5 * magnetic_field_[0] * magnetic_field_[1] * ox * oy * Sz;
                current_data[cnt] -= 0.5 * magnetic_field_[1] * magnetic_field_[2] * Sx * oy * oz;
                current_data[cnt] -= 0.5 * magnetic_field_[2] * magnetic_field_[0] * ox * Sy * oz;
                current_data[cnt] += imag * magnetic_field_[2] * dx * oy * Sz;
                current_data[cnt] -= imag * magnetic_field_[1] * dx * Sy * oz;
                current_data[cnt] += imag * magnetic_field_[0] * Sx * dy * oz;
                current_data[cnt] -= imag * magnetic_field_[2] * ox * dy * Sz;
                current_data[cnt] += imag * magnetic_field_[1] * ox * Sy * dz;
                current_data[cnt] -= imag * magnetic_field_[0] * Sx * oy * dz;
                current_data[cnt] *= 0.5;
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

complex<double> ComplexKineticBatch::get_P(const double coord1, const double coord2, const double exp1, const double exp2, const double one12,
                                                const int dim, const bool swap) {
  const double Areal = coord1*exp1;
  const double Breal = coord2*exp2;
  const double Aimag = basisinfo_[0]->vector_potential(dim);
  const double Bimag = basisinfo_[1]->vector_potential(dim);
  const double imag = 0.5 * (swap ? Bimag - Aimag : Aimag - Bimag );
  const complex<double> num (Areal + Breal, imag);
  return num * one12;
}

