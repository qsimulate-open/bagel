//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: point_complexmomentumbatch.cc
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
#include <src/util/timer.h>
#include <src/integral/compos/point_complexmomentumbatch.h>

using namespace std;
using namespace bagel;

const static CCarSphList carsphlist;

static const double pisqrtinv__ = 1.0 / sqrt(pi__);

void Point_ComplexMomentumBatch::compute() {


  const CSortList sort_(spherical_);

  complex<double>* const intermediate_p = stack_->get<complex<double>>(size_block_*3);
  fill_n(intermediate_p, size_block_*3, 0.0);
  perform_VRR(intermediate_p);

  for (int i = 0; i != 3; ++i) {
    complex<double>* const cdata = data_ + i*size_block_;
    const complex<double>* const csource = intermediate_p + i*size_block_;
    complex<double>* const intermediate_c = stack_->get<complex<double>>(cont0_ * cont1_ * asize_intermediate_);
    fill_n(intermediate_c, cont0_ * cont1_ * asize_intermediate_, 0.0);
    perform_contraction(asize_intermediate_, csource, prim0_, prim1_, intermediate_c,
                        basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0_,
                        basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1_);

    if (basisinfo_[0]->spherical() && basisinfo_[1]->spherical()) {
      // transform both indices to spherical
      complex<double>* const intermediate_i = stack_->get<complex<double>>(cont0_ * cont1_ * asize_final_);
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


void Point_ComplexMomentumBatch::perform_VRR(complex<double>* intermediate) {

  const int ang0s = ang0_+1;
  const int worksize = ang0s*(ang1_+2);
  complex<double>* worksx = stack_->get<complex<double>>(worksize);
  complex<double>* worksy = stack_->get<complex<double>>(worksize);
  complex<double>* worksz = stack_->get<complex<double>>(worksize);

  complex<double> imag (0.0, 1.0);
  if (swap01_) imag *= -1.0;
  const double Bx = basisinfo_[1]->position(0);
  const double By = basisinfo_[1]->position(1);
  const double Bz = basisinfo_[1]->position(2);
  const double A_Bx = basisinfo_[1]->vector_potential(0);
  const double A_By = basisinfo_[1]->vector_potential(1);
  const double A_Bz = basisinfo_[1]->vector_potential(2);
  const double field_x_half = 0.5 * magnetic_field_[0];
  const double field_y_half = 0.5 * magnetic_field_[1];
  const double field_z_half = 0.5 * magnetic_field_[2];

  const double pAx = location_[0] - basisinfo_[0]->position(0);
  const double pAy = location_[1] - basisinfo_[0]->position(1);
  const double pAz = location_[2] - basisinfo_[0]->position(2);
  const double pBx = location_[0] - basisinfo_[1]->position(0);
  const double pBy = location_[1] - basisinfo_[1]->position(1);
  const double pBz = location_[2] - basisinfo_[1]->position(2);

  for (int ii = 0; ii != prim0_*prim1_; ++ii) {
    // Perform VRR
    int offset_ii = ii*asize_intermediate_;
    const double cb = xb_[ii];
    complex<double>* current_data = &intermediate[offset_ii];

    const double cxp = xp_[ii];
    const complex<double> cxpl = location_[0] - P_[ii*3  ];
    const complex<double> cypl = location_[1] - P_[ii*3+1];
    const complex<double> czpl = location_[2] - P_[ii*3+2];

    // obtain S(0, 0)
    const double tmp = pisqrtinv__ * sqrt(cxp);
    worksx[0] = coeffsx_[ii] * tmp * exp(- cxp * (cxpl * cxpl));
    worksy[0] = coeffsy_[ii] * tmp * exp(- cxp * (cypl * cypl));
    worksz[0] = coeffsz_[ii] * tmp * exp(- cxp * (czpl * czpl));

    for (int i = 1; i <= ang0_; ++i) {
      // obtain S(i, 0)
      worksx[i] = pAx*worksx[i-1];
      worksy[i] = pAy*worksy[i-1];
      worksz[i] = pAz*worksz[i-1];
    }

    // peform HRR
    for (int j = 1; j <= ang1_ + 1; ++j) {
      for (int i = 0; i <= ang0_; ++i) {
        // obtain S(i, j)
        worksx[j*ang0s+i] = pBx*worksx[(j-1)*ang0s+i];
        worksy[j*ang0s+i] = pBy*worksy[(j-1)*ang0s+i];
        worksz[j*ang0s+i] = pBz*worksz[(j-1)*ang0s+i];
      }
    }


    complex<double> Sx, Sy, Sz;    // operator = 1
    complex<double> ox, oy, oz;    // operator = x
    complex<double> dx, dy, dz;    // operator = d/dx

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
                const double jxd = jx;
                const double jyd = jy;
                const double jzd = jz;
                complex<double> minus1x = 0.0;
                complex<double> minus1y = 0.0;
                complex<double> minus1z = 0.0;
                if (jx > 0) minus1x = jxd * worksx[ix+ang0s*(jx-1)];
                if (jy > 0) minus1y = jyd * worksy[iy+ang0s*(jy-1)];
                if (jz > 0) minus1z = jzd * worksz[iz+ang0s*(jz-1)];

                Sx = worksx[ix + ang0s * jx];
                Sy = worksy[iy + ang0s * jy];
                Sz = worksz[iz + ang0s * jz];

                ox = worksx[ix + ang0s * (jx+1)] + Bx * worksx[ix + ang0s * jx];
                oy = worksy[iy + ang0s * (jy+1)] + By * worksy[iy + ang0s * jy];
                oz = worksz[iz + ang0s * (jz+1)] + Bz * worksz[iz + ang0s * jz];

                dx = minus1x - 2.0 * cb * worksx[ix+ang0s * (jx+1)] - imag * A_Bx * worksx[ix + ang0s * jx];
                dy = minus1y - 2.0 * cb * worksy[iy+ang0s * (jy+1)] - imag * A_By * worksy[iy + ang0s * jy];
                dz = minus1z - 2.0 * cb * worksz[iz+ang0s * (jz+1)] - imag * A_Bz * worksz[iz + ang0s * jz];

                current_data[cnt              ] = -imag * dx * Sy * Sz;
                current_data[cnt+size_block_  ] = -imag * Sx * dy * Sz;
                current_data[cnt+size_block_*2] = -imag * Sx * Sy * dz;

                current_data[cnt              ] += field_y_half * Sx * Sy * oz;
                current_data[cnt+size_block_  ] += field_z_half * ox * Sy * Sz;
                current_data[cnt+size_block_*2] += field_x_half * Sx * oy * Sz;
                current_data[cnt              ] -= field_z_half * Sx * oy * Sz;
                current_data[cnt+size_block_  ] -= field_x_half * Sx * Sy * oz;
                current_data[cnt+size_block_*2] -= field_y_half * ox * Sy * Sz;

                ++cnt;
              }
            }
          }
        }
      }
    }
    assert(cnt == asize_intermediate_);

  } // end of prim exponent loop

  stack_->release(worksize, worksz);
  stack_->release(worksize, worksy);
  stack_->release(worksize, worksx);
}

complex<double> Point_ComplexMomentumBatch::get_P(const double coord1, const double coord2, const double exp1, const double exp2, const double one12,
                                                       const int dim, const bool swap) {
  const double Areal = coord1*exp1;
  const double Breal = coord2*exp2;
  const double Aimag = basisinfo_[0]->vector_potential(dim);
  const double Bimag = basisinfo_[1]->vector_potential(dim);
  const double imag = 0.5 * (swap ? Bimag - Aimag : Aimag - Bimag );
  const complex<double> num (Areal + Breal, imag);
  return num * one12;
}

