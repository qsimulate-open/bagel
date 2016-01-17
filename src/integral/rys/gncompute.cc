//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gncompute.cc
// Copyright (C) 2012 Toru Shiozaki
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
#include <src/integral/sortlist.h>
#include <src/integral/rys/gnaibatch.h>
#include <src/util/math/comb.h>

using namespace std;
using namespace bagel;

const static Comb comb;
const static CarSphList carsphlist;

void GNAIBatch::compute() {

  fill_n(data_, size_alloc_, 0.0);

  // temp are for VRR
  const int worksize = rank_ * amax1_;
  double* workx = stack_->get(worksize);
  double* worky = stack_->get(worksize);
  double* workz = stack_->get(worksize);
  double r1x[RYS_MAX];
  double r1y[RYS_MAX];
  double r1z[RYS_MAX];
  double r2[RYS_MAX];


  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);

  // transformation matrix
  const int a = basisinfo_[0]->angular_number();
  const int b = basisinfo_[1]->angular_number();

  const int a2 = a+2;
  const int b2 = b+2;
  assert(a+b+1 == amax_);

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
        transx[i + (amax_+1)*k] = comb(ib, ia+ib-i) * pow(AB_[0], ia+ib-i);
        transy[i + (amax_+1)*k] = comb(ib, ia+ib-i) * pow(AB_[1], ia+ib-i);
        transz[i + (amax_+1)*k] = comb(ib, ia+ib-i) * pow(AB_[2], ia+ib-i);
      }
    }
  }
  double* const bufx = stack_->get(rank_*a2*b2);
  double* const bufy = stack_->get(rank_*a2*b2);
  double* const bufz = stack_->get(rank_*a2*b2);
  double* const bufx_a = stack_->get(rank_*a2*b2);
  double* const bufx_b = stack_->get(rank_*a2*b2);
  double* const bufx_c = stack_->get(rank_*a2*b2);
  double* const bufy_a = stack_->get(rank_*a2*b2);
  double* const bufy_b = stack_->get(rank_*a2*b2);
  double* const bufy_c = stack_->get(rank_*a2*b2);
  double* const bufz_a = stack_->get(rank_*a2*b2);
  double* const bufz_b = stack_->get(rank_*a2*b2);
  double* const bufz_c = stack_->get(rank_*a2*b2);

  const int acsize = (a+1)*(a+2)*(b+1)*(b+2)/4;
  assert(acsize*primsize_ == size_block_);

  const SortList sort(spherical1_);

  // perform VRR
  for (int xj = 0; xj != screening_size_; ++xj) {
    const int i = screening_[xj];
    const int iprim = i / natom_;
    const int iatom = i % natom_;
    const int offset_iprim = iprim * acsize;

    const double* croots = &roots_[i * rank_];
    const double* cweights = &weights_[i * rank_];
    double PC[3];
    PC[0] = P_[i*3  ] - mol_->atoms(iatom)->position(0);
    PC[1] = P_[i*3+1] - mol_->atoms(iatom)->position(1);
    PC[2] = P_[i*3+2] - mol_->atoms(iatom)->position(2);
    for (int r = 0; r != rank_; ++r) {
      r1x[r] = P_[i*3  ] - ax - PC[0] * croots[r];
      r1y[r] = P_[i*3+1] - ay - PC[1] * croots[r];
      r1z[r] = P_[i*3+2] - az - PC[2] * croots[r];
      r2[r] = (1.0 - croots[r]) * 0.5 / xp_[i];
      workx[r] = cweights[r] * coeff_[i];
      worky[r] = 1.0;
      workz[r] = 1.0;
    }

    if (amax_ > 0) {
      for (int r = 0; r != rank_; ++r) {
        workx[rank_ + r] = r1x[r] * workx[r];
        worky[rank_ + r] = r1y[r] * worky[r];
        workz[rank_ + r] = r1z[r] * workz[r];
      }
      for (int j = 2; j <= amax_; ++j) {
        const int offset = j * rank_;
        for (int r = 0; r != rank_; ++r) {
          workx[offset + r] = r1x[r] * workx[offset - rank_ + r] + (j - 1) * r2[r] * workx[offset - rank_ * 2 + r];
          worky[offset + r] = r1y[r] * worky[offset - rank_ + r] + (j - 1) * r2[r] * worky[offset - rank_ * 2 + r];
          workz[offset + r] = r1z[r] * workz[offset - rank_ + r] + (j - 1) * r2[r] * workz[offset - rank_ * 2 + r];
        }
      }
    } else {
      assert(false);
    }

    // HRR is done in one shot
    dgemm_("N", "N", rank_, b2*a2, amax_+1, 1.0, workx, rank_, transx, amax_+1, 0.0, bufx, rank_);
    dgemm_("N", "N", rank_, b2*a2, amax_+1, 1.0, worky, rank_, transy, amax_+1, 0.0, bufy, rank_);
    dgemm_("N", "N", rank_, b2*a2, amax_+1, 1.0, workz, rank_, transz, amax_+1, 0.0, bufz, rank_);

    unsigned int aatom, batom;
    tie(aatom, batom) = iatom_;
    const unsigned int catom = iatom;
    const double alpha = exponents_[iprim*2+0];
    const double beta_ = exponents_[iprim*2+1];

    double xpPC2[3];
    xpPC2[0] = 2.0 * xp_[i] * PC[0];
    xpPC2[1] = 2.0 * xp_[i] * PC[1];
    xpPC2[2] = 2.0 * xp_[i] * PC[2];
    for (int ib = 0; ib <= b; ++ib) {
      for (int ia = 0; ia <= a; ++ia) {
        for (int r = 0; r != rank_; ++r) {
          bufx_a[r+rank_*(ia+a2*ib)] = 2.0*alpha*bufx[r+rank_*(ia+1+a2*(ib))] - ia*bufx[r+rank_*(ia-1+a2*(ib))];
          bufx_b[r+rank_*(ia+a2*ib)] = 2.0*beta_*bufx[r+rank_*(ia+a2*(ib+1))] - ib*bufx[r+rank_*(ia+a2*(ib-1))];
          bufx_c[r+rank_*(ia+a2*ib)] = croots[r] * (xpPC2[0] * bufx[r+rank_*(ia+a2*ib)]
                                     + ia*bufx[r+rank_*(ia-1+a2*(ib))] + ib*bufx[r+rank_*(ia+a2*(ib-1))]);
          bufy_a[r+rank_*(ia+a2*ib)] = 2.0*alpha*bufy[r+rank_*(ia+1+a2*(ib))] - ia*bufy[r+rank_*(ia-1+a2*(ib))];
          bufy_b[r+rank_*(ia+a2*ib)] = 2.0*beta_*bufy[r+rank_*(ia+a2*(ib+1))] - ib*bufy[r+rank_*(ia+a2*(ib-1))];
          bufy_c[r+rank_*(ia+a2*ib)] = croots[r] * (xpPC2[1] * bufy[r+rank_*(ia+a2*ib)]
                                     + ia*bufy[r+rank_*(ia-1+a2*(ib))] + ib*bufy[r+rank_*(ia+a2*(ib-1))]);
          bufz_a[r+rank_*(ia+a2*ib)] = 2.0*alpha*bufz[r+rank_*(ia+1+a2*(ib))] - ia*bufz[r+rank_*(ia-1+a2*(ib))];
          bufz_b[r+rank_*(ia+a2*ib)] = 2.0*beta_*bufz[r+rank_*(ia+a2*(ib+1))] - ib*bufz[r+rank_*(ia+a2*(ib-1))];
          bufz_c[r+rank_*(ia+a2*ib)] = croots[r] * (xpPC2[2] * bufz[r+rank_*(ia+a2*ib)]
                                     + ia*bufz[r+rank_*(ia-1+a2*(ib))] + ib*bufz[r+rank_*(ia+a2*(ib-1))]);
        }
      }
    }

    // assembly step
    double* current_data0 = data_ + size_block_*(3*aatom+0) + offset_iprim;
    double* current_data1 = data_ + size_block_*(3*aatom+1) + offset_iprim;
    double* current_data2 = data_ + size_block_*(3*aatom+2) + offset_iprim;
    double* current_data3 = data_ + size_block_*(3*batom+0) + offset_iprim;
    double* current_data4 = data_ + size_block_*(3*batom+1) + offset_iprim;
    double* current_data5 = data_ + size_block_*(3*batom+2) + offset_iprim;
    double* current_data6 = data_ + size_block_*(3*catom+0) + offset_iprim;
    double* current_data7 = data_ + size_block_*(3*catom+1) + offset_iprim;
    double* current_data8 = data_ + size_block_*(3*catom+2) + offset_iprim;

    for (int iaz = 0; iaz <= a; ++iaz) {
      for (int iay = 0; iay <= a - iaz; ++iay) {
        const int iax = a - iaz - iay;

        for (int ibz = 0; ibz <= b; ++ibz) {
          for (int iby = 0; iby <= b - ibz; ++iby) {
            const int ibx = b - ibz - iby;
            for (int r = 0; r != rank_; ++r) {
              *current_data0 += bufx_a[r+rank_*(iax+a2*ibx)] * bufy  [r+rank_*(iay+a2*iby)] * bufz  [r+rank_*(iaz+a2*ibz)];
              *current_data1 += bufx  [r+rank_*(iax+a2*ibx)] * bufy_a[r+rank_*(iay+a2*iby)] * bufz  [r+rank_*(iaz+a2*ibz)];
              *current_data2 += bufx  [r+rank_*(iax+a2*ibx)] * bufy  [r+rank_*(iay+a2*iby)] * bufz_a[r+rank_*(iaz+a2*ibz)];
              *current_data3 += bufx_b[r+rank_*(iax+a2*ibx)] * bufy  [r+rank_*(iay+a2*iby)] * bufz  [r+rank_*(iaz+a2*ibz)];
              *current_data4 += bufx  [r+rank_*(iax+a2*ibx)] * bufy_b[r+rank_*(iay+a2*iby)] * bufz  [r+rank_*(iaz+a2*ibz)];
              *current_data5 += bufx  [r+rank_*(iax+a2*ibx)] * bufy  [r+rank_*(iay+a2*iby)] * bufz_b[r+rank_*(iaz+a2*ibz)];
              *current_data6 += bufx_c[r+rank_*(iax+a2*ibx)] * bufy  [r+rank_*(iay+a2*iby)] * bufz  [r+rank_*(iaz+a2*ibz)];
              *current_data7 += bufx  [r+rank_*(iax+a2*ibx)] * bufy_c[r+rank_*(iay+a2*iby)] * bufz  [r+rank_*(iaz+a2*ibz)];
              *current_data8 += bufx  [r+rank_*(iax+a2*ibx)] * bufy  [r+rank_*(iay+a2*iby)] * bufz_c[r+rank_*(iaz+a2*ibz)];
            }
            ++current_data0;
            ++current_data1;
            ++current_data2;
            ++current_data3;
            ++current_data4;
            ++current_data5;
            ++current_data6;
            ++current_data7;
            ++current_data8;
          }
        }
      }
    }
  }

  stack_->release(rank_*a2*b2, bufz_c);
  stack_->release(rank_*a2*b2, bufz_b);
  stack_->release(rank_*a2*b2, bufz_a);
  stack_->release(rank_*a2*b2, bufy_c);
  stack_->release(rank_*a2*b2, bufy_b);
  stack_->release(rank_*a2*b2, bufy_a);
  stack_->release(rank_*a2*b2, bufx_c);
  stack_->release(rank_*a2*b2, bufx_b);
  stack_->release(rank_*a2*b2, bufx_a);
  stack_->release(rank_*a2*b2, bufz);
  stack_->release(rank_*a2*b2, bufy);
  stack_->release(rank_*a2*b2, bufx);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize, workz);
  stack_->release(worksize, worky);
  stack_->release(worksize, workx);

  double* const buf = stack_->get(size_block_);
  bkup_ = buf;
  double* cdata = data_;
  for (int i = 0; i != natom_*3; ++i, cdata += size_block_) {
    double* target = bkup_;
    const double* source = cdata;
    // contract indices 01
    // data will be stored in bkup_: cont01{ xyz{ } }
    const int m = acsize;
    perform_contraction(m, source, prim0size_, prim1size_, target,
                        basisinfo_[0]->contractions(), basisinfo_[0]->contraction_ranges(), cont0size_,
                        basisinfo_[1]->contractions(), basisinfo_[1]->contraction_ranges(), cont1size_);

    // Cartesian to spherical 01 if necesarry
    // data will be stored in bkup_
    target = cdata;
    source = bkup_;
    if (spherical1_) {
      const int carsphindex = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      const int nloops = contsize_;
      carsphlist.carsphfunc_call(carsphindex, nloops, source, target);
    }

    // Sort cont01 and xyzab
    // data will be stored in data_: cont1b{ cont0a{ } }
    if (spherical1_) {
      target = bkup_;
      source = cdata;
      const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort.sortfunc_call(index, target, source, cont1size_, cont0size_, 1, swap01_);
      copy(bkup_, bkup_+size_block_, cdata);
    } else {
      target = cdata;
      source = bkup_;
      const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort.sortfunc_call(index, target, source, cont1size_, cont0size_, 1, swap01_);
    }
  }

  stack_->release(size_block_, buf);

}


