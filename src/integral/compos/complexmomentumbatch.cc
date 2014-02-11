//
// BAGEL - Parallel electron correlation program.
// Filename: complexmomentumbatch.cc
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
#include <src/integral/compos/complexmomentumbatch.h>

using namespace std;
using namespace bagel;

const static CCarSphList carsphlist;


void ComplexMomentumBatch::compute() {

  const CSortList sort_ (spherical_);

  complex<double>* const intermediate_p = stack_->get<complex<double>>(size_block_*3);
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

  if (swap01_) {for (int n=0; n!=size_block_*3; n++) data_[n] *= -1.0;}

}


void ComplexMomentumBatch::perform_VRR(complex<double>* intermediate) {
  const int amax2 = amax1_+1;
  const int worksize = amax2;
  complex<double>* worksx = stack_->get<complex<double>>(worksize*worksize);
  complex<double>* worksy = stack_->get<complex<double>>(worksize*worksize);
  complex<double>* worksz = stack_->get<complex<double>>(worksize*worksize);

  for (int ii = 0; ii != prim0_*prim1_; ++ii) {
    // Perform VRR
    int offset_ii = ii*asize_intermediate_;
    const double cop = 1.0 / xp_[ii];
    const complex<double> cxpa = P_[ii*3  ]-basisinfo_[0]->position(0);
    const complex<double> cypa = P_[ii*3+1]-basisinfo_[0]->position(1);
    const complex<double> czpa = P_[ii*3+2]-basisinfo_[0]->position(2);
    complex<double>* current_data = &intermediate[offset_ii];

    // obtain S(0, 0)
    worksx[0] = coeffsx_[ii];
    worksy[0] = coeffsy_[ii];
    worksz[0] = coeffsz_[ii];

    // obtain S(1, 0)
    worksx[1] = cxpa*worksx[0];
    worksy[1] = cypa*worksy[0];
    worksz[1] = czpa*worksz[0];

    for (int i = 2; i != amax2; ++i) {
      // obtain S(i, 0)
      worksx[i] = cxpa*worksx[i-1]+0.5*(i-1)*cop*worksx[i-2];
      worksy[i] = cypa*worksy[i-1]+0.5*(i-1)*cop*worksy[i-2];
      worksz[i] = czpa*worksz[i-1]+0.5*(i-1)*cop*worksz[i-2];
    }

    // peform HRR
    if (ang1_ > 0) {
      for (int j = 1; j <= ang1_ + 1; ++j) {
        for (int i = 0; i != amax2 - j; ++i) {
          // obtain S(i, j)
          worksx[j*amax2+i] = AB_[0]*worksx[(j-1)*amax2+i]+worksx[(j-1)*amax2+i+1];
          worksy[j*amax2+i] = AB_[1]*worksy[(j-1)*amax2+i]+worksy[(j-1)*amax2+i+1];
          worksz[j*amax2+i] = AB_[2]*worksz[(j-1)*amax2+i]+worksz[(j-1)*amax2+i+1];
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
                current_data[cnt              ] = 0.0; //workpx[ix+amax1_*jx]*worksy[iy+amax1_*jy]*worksz[iz+amax1_*jz];
                current_data[cnt+size_block_  ] = 0.0; //worksx[ix+amax1_*jx]*workpy[iy+amax1_*jy]*worksz[iz+amax1_*jz];
                current_data[cnt+size_block_*2] = 0.0; //worksx[ix+amax1_*jx]*worksy[iy+amax1_*jy]*workpz[iz+amax1_*jz];
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
}

// TODO For efficiency's sake, it's probably best to find a way to avoid repeatedly running basisinfo_[i]->vector_potential(j) each time we want to get a P or Q value
std::complex<double> ComplexMomentumBatch::get_P(const double coord1, const double coord2, const double exp1, const double exp2, const double one12,
                                                 const int dim, const bool swap) {
  const double Areal = coord1*exp1;
  const double Breal = coord2*exp2;
  const double Aimag = basisinfo_[0]->vector_potential(dim);
  const double Bimag = basisinfo_[1]->vector_potential(dim);
  double imag;
  if (swap) imag = 0.5*(Bimag - Aimag);
  else imag = 0.5*(Aimag - Bimag);
  const std::complex<double> num (Areal + Breal, imag);
  return num * one12;
}

