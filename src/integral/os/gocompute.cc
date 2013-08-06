//
// BAGEL - Parallel electron correlation program.
// Filename: gocompute.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#include <stddef.h>
#include <src/integral/carsphlist.h>
#include <src/integral/os/goverlapbatch.h>
#include <src/math/comb.h>

using namespace std;
using namespace bagel;

const static Comb comb;

void GOverlapBatch::compute() {

  fill_n(data_, size_alloc_, 0.0);

  const int a = basisinfo_[0]->angular_number();
  const int b = basisinfo_[1]->angular_number();
  const int a2 = a+2;
  const int b2 = b+2;
  assert(amax_ == a+b+1);
  const int acsize = (a+1)*(a+2)*(b+1)*(b+2)/4;

  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  fill_n(transx, (amax_+1)*a2*b2, 0.0);
  fill_n(transy, (amax_+1)*a2*b2, 0.0);
  fill_n(transz, (amax_+1)*a2*b2, 0.0);
  for (int ib = 0, k = 0; ib <= b+1; ++ib) {
    for (int ia = 0; ia <= a+1; ++ia, ++k) {
      if (ia == a+1 && ib == b+1) continue;
      for (int i = ia; i <= ia+ib; ++i) {
        transx[i + (amax_+1)*k] = comb.c(ib, ia+ib-i) * pow(AB_[0], ia+ib-i);
        transy[i + (amax_+1)*k] = comb.c(ib, ia+ib-i) * pow(AB_[1], ia+ib-i);
        transz[i + (amax_+1)*k] = comb.c(ib, ia+ib-i) * pow(AB_[2], ia+ib-i);
      }
    }
  }
  const int worksize = amax_+1;
  double* const workx = stack_->get(worksize);
  double* const worky = stack_->get(worksize);
  double* const workz = stack_->get(worksize);

  double* const bufx = stack_->get(a2*b2);
  double* const bufy = stack_->get(a2*b2);
  double* const bufz = stack_->get(a2*b2);
  double* const bufx_a = stack_->get(a2*b2);
  double* const bufx_b = stack_->get(a2*b2);
  double* const bufy_a = stack_->get(a2*b2);
  double* const bufy_b = stack_->get(a2*b2);
  double* const bufz_a = stack_->get(a2*b2);
  double* const bufz_b = stack_->get(a2*b2);

  // Perform VRR
  for (int ii = 0; ii != prim0_ * prim1_; ++ii) {

    /// Sx(0 : i+j+1, 0) will be made here
    workx[0] = coeffsx_[ii];
    worky[0] = coeffsy_[ii];
    workz[0] = coeffsz_[ii];
    workx[1] = (p_[ii * 3    ] - basisinfo_[0]->position(0)) * workx[0];
    worky[1] = (p_[ii * 3 + 1] - basisinfo_[0]->position(1)) * worky[0];
    workz[1] = (p_[ii * 3 + 2] - basisinfo_[0]->position(2)) * workz[0];
    for (int i = 2; i <= amax_; ++i) {
      workx[i] = (p_[ii * 3    ] - basisinfo_[0]->position(0)) * workx[i - 1] + 0.5 * (i - 1) / xp_[ii] * workx[i - 2];
      worky[i] = (p_[ii * 3 + 1] - basisinfo_[0]->position(1)) * worky[i - 1] + 0.5 * (i - 1) / xp_[ii] * worky[i - 2];
      workz[i] = (p_[ii * 3 + 2] - basisinfo_[0]->position(2)) * workz[i - 1] + 0.5 * (i - 1) / xp_[ii] * workz[i - 2];
    }
    // HRR is done in one shot
    dgemv_("T", amax_+1, a2*b2, 1.0, transx, amax_+1, workx, 1, 0.0, bufx, 1);
    dgemv_("T", amax_+1, a2*b2, 1.0, transy, amax_+1, worky, 1, 0.0, bufy, 1);
    dgemv_("T", amax_+1, a2*b2, 1.0, transz, amax_+1, workz, 1, 0.0, bufz, 1);

    const double alpha = xa_[ii];
    const double beta_ = xb_[ii];
    for (int ib = 0; ib <= b; ++ib) {
      for (int ia = 0; ia <= a; ++ia) {
        bufx_a[ia+a2*ib] = 2.0*alpha*bufx[ia+1+a2*(ib)] - ia*bufx[ia-1+a2*(ib)];
        bufx_b[ia+a2*ib] = 2.0*beta_*bufx[ia+a2*(ib+1)] - ib*bufx[ia+a2*(ib-1)];
        bufy_a[ia+a2*ib] = 2.0*alpha*bufy[ia+1+a2*(ib)] - ia*bufy[ia-1+a2*(ib)];
        bufy_b[ia+a2*ib] = 2.0*beta_*bufy[ia+a2*(ib+1)] - ib*bufy[ia+a2*(ib-1)];
        bufz_a[ia+a2*ib] = 2.0*alpha*bufz[ia+1+a2*(ib)] - ia*bufz[ia-1+a2*(ib)];
        bufz_b[ia+a2*ib] = 2.0*beta_*bufz[ia+a2*(ib+1)] - ib*bufz[ia+a2*(ib-1)];
      }
    }

    /// assembly process
    const int offset_ii = ii * acsize;
    double* current_data0 = data_ + offset_ii;
    double* current_data1 = data_ + offset_ii + size_block_;
    double* current_data2 = data_ + offset_ii + size_block_*2;
    double* current_data3 = data_ + offset_ii + size_block_*3;
    double* current_data4 = data_ + offset_ii + size_block_*4;
    double* current_data5 = data_ + offset_ii + size_block_*5;

    for (int iaz = 0; iaz <= a; ++iaz) {
      for (int iay = 0; iay <= a - iaz; ++iay) {
        const int iax = a - iaz - iay;
        for (int ibz = 0; ibz <= b; ++ibz) {
          for (int iby = 0; iby <= b - ibz; ++iby) {
            const int ibx = b - ibz - iby;

            *current_data0++ += bufx_a[iax+a2*ibx] * bufy  [iay+a2*iby] * bufz  [iaz+a2*ibz];
            *current_data1++ += bufx  [iax+a2*ibx] * bufy_a[iay+a2*iby] * bufz  [iaz+a2*ibz];
            *current_data2++ += bufx  [iax+a2*ibx] * bufy  [iay+a2*iby] * bufz_a[iaz+a2*ibz];
            *current_data3++ += bufx_b[iax+a2*ibx] * bufy  [iay+a2*iby] * bufz  [iaz+a2*ibz];
            *current_data4++ += bufx  [iax+a2*ibx] * bufy_b[iay+a2*iby] * bufz  [iaz+a2*ibz];
            *current_data5++ += bufx  [iax+a2*ibx] * bufy  [iay+a2*iby] * bufz_b[iaz+a2*ibz];
          }
        }
      }
    }

  } // end of primsize loop

  stack_->release(a2*b2, bufz_b);
  stack_->release(a2*b2, bufz_a);
  stack_->release(a2*b2, bufy_b);
  stack_->release(a2*b2, bufy_a);
  stack_->release(a2*b2, bufx_b);
  stack_->release(a2*b2, bufx_a);
  stack_->release(a2*b2, bufz);
  stack_->release(a2*b2, bufy);
  stack_->release(a2*b2, bufx);

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  const size_t acpsize = acsize*cont0_*cont1_;
  double* const bkup = stack_->get(acpsize);
  double* cdata = data_;
  for (int i = 0; i != 6; ++i, cdata += size_block_) {
    // first, contraction.
    const double* source = cdata;
    double* target = bkup;
    perform_contraction(acsize, source, prim0_, prim1_, target,
                        basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0_,
                        basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1_);

    if (spherical_) {
      struct CarSphList carsphlist;
      const unsigned int carsph_index = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      const int nloops = cont0_ * cont1_;
      source = bkup;
      target = cdata;
      carsphlist.carsphfunc_call(carsph_index, nloops, source, target);

      const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      source = cdata;
      target = bkup;
      sort_.sortfunc_call(sort_index, target, source, cont1_, cont0_, 1, swap01_);
      copy_n(bkup, acpsize, cdata);
    } else {
      const unsigned int sort_index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      source = bkup;
      target = cdata;
      sort_.sortfunc_call(sort_index, target, source, cont1_, cont0_, 1, swap01_);
    }
  }

  stack_->release(acpsize, bkup);

}

