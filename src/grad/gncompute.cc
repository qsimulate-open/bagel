//
// Newint - Parallel electron correlation program.
// Filename: gncompute.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/stackmem.h>
#include <src/grad/gnaibatch.h>
#include <src/rysint/carsphlist.h>
#include <src/fci/comb.h>

using namespace std;

extern StackMem* stack;
static Comb comb;

void GNAIBatch::compute() {
  double* const stack_save = stack->get(0);

  bkup_ = stack->get(size_alloc_/12);

  const int worksize = rank_ * amax1_;
  
  double* workx = stack->get(worksize);
  double* worky = stack->get(worksize);
  double* workz = stack->get(worksize);

  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);

  fill(data_, data_ + size_alloc_, 0.0);

  double r1x[20];
  double r1y[20];
  double r1z[20];
  double r2[20];

  // transformation matrix
  const int a = basisinfo_[0]->angular_number();
  const int b = basisinfo_[1]->angular_number();
  assert(a+b+1 == amax_); 

  const int a2 = a+2;
  const int b2 = b+2;

  double* transx = stack->get((amax_+1)*a2*b2);
  double* transy = stack->get((amax_+1)*a2*b2);
  double* transz = stack->get((amax_+1)*a2*b2);
  fill(transx, transx+(amax_+1)*a2*b2, 0.0);
  fill(transy, transy+(amax_+1)*a2*b2, 0.0);
  fill(transz, transz+(amax_+1)*a2*b2, 0.0);
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
  double* bufx = stack->get(rank_*a2*b2);
  double* bufy = stack->get(rank_*a2*b2);
  double* bufz = stack->get(rank_*a2*b2);
  double* bufx_a = stack->get(rank_*a2*b2);
  double* bufx_b = stack->get(rank_*a2*b2);
  double* bufx_c = stack->get(rank_*a2*b2);
  double* bufy_a = stack->get(rank_*a2*b2);
  double* bufy_b = stack->get(rank_*a2*b2);
  double* bufy_c = stack->get(rank_*a2*b2);
  double* bufz_a = stack->get(rank_*a2*b2);
  double* bufz_b = stack->get(rank_*a2*b2);
  double* bufz_c = stack->get(rank_*a2*b2);

  const int acsize = size_alloc_ / 12 / primsize_;
  const size_t acpsize = acsize*primsize_;
  assert(acpsize == (a+1)*(a+2)*(b+1)*(b+2)/4*primsize_); 
  assert(primsize_ == prim0size_*prim1size_);

  double* const data0 = data_;
#if 0
  double* const data1 = data_ + acpsize;
  double* const data2 = data_ + acpsize*2;
  double* const data3 = data_ + acpsize*3;
  double* const data4 = data_ + acpsize*4;
  double* const data5 = data_ + acpsize*5;
  double* const data6 = data_ + acpsize*6;
  double* const data7 = data_ + acpsize*7;
  double* const data8 = data_ + acpsize*8;
#endif


  // perform VRR
  const int natom_unit = natom_ / (2 * L_ + 1);
  assert(natom_ % (2 * L_ + 1) == 0);
  for (int xj = 0; xj != screening_size_; ++xj) {
    const int i = screening_[xj];
    const int iprim = i / natom_;
    const int resid = i % natom_;
    const int cell  = resid / natom_unit - L_;
    const int iatom = resid % natom_unit;
    double disp[3];
    disp[0] = disp[1] = 0.0;
    disp[2] = A_ * cell;
    const int offset_iprim = iprim * acsize;

    if (cell != 0) throw logic_error("I haven't thought about periodic cases");

    const double* croots = &roots_[i * rank_]; 
    const double* cweights = &weights_[i * rank_]; 
    double PC[3];
    PC[0] = p_[i*3  ] - geom_->atoms(iatom)->position(0) - disp[0];
    PC[1] = p_[i*3+1] - geom_->atoms(iatom)->position(1) - disp[1];
    PC[2] = p_[i*3+2] - geom_->atoms(iatom)->position(2) - disp[2];
    for (int r = 0; r != rank_; ++r) {
      r1x[r] = p_[i*3  ] - ax - PC[0] * croots[r];
      r1y[r] = p_[i*3+1] - ay - PC[1] * croots[r];
      r1z[r] = p_[i*3+2] - az - PC[2] * croots[r];
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

    double* expo = exponents_.get() + i*2;
    double xpPC2[3];
    xpPC2[0] = 2.0 * xp_[i] * PC[0];
    xpPC2[1] = 2.0 * xp_[i] * PC[1];
    xpPC2[2] = 2.0 * xp_[i] * PC[2];
#if 0
    for (int ib = 0; ib <= b; ++ib) {
      const double ibb = ib;
      for (int ia = 0; ia <= a; ++ia) {
        const double iaa = ia;
        for (int r = 0; r != rank_; ++r) {
          bufx_a[r+rank_*(ia+a2*ib)] = 2.0*expo[0]*bufx[r+rank_*(ia+1+a2*(ib))] - ia*bufx[r+rank_*(ia-1+a2*(ib))];
          bufx_b[r+rank_*(ia+a2*ib)] = 2.0*expo[1]*bufx[r+rank_*(ia+a2*(ib+1))] - ia*bufx[r+rank_*(ia+a2*(ib-1))];
          bufx_c[r+rank_*(ia+a2*ib)] = croots[r] * (xpPC2[0] * bufx[r+rank_*(ia+a2*ib)]
                                     + ia*bufx[r+rank_*(ia-1+a2*(ib))] + ib*bufx[r+rank_*(ia+a2*(ib-1))]);
          bufy_a[r+rank_*(ia+a2*ib)] = 2.0*expo[0]*bufy[r+rank_*(ia+1+a2*(ib))] - ia*bufy[r+rank_*(ia-1+a2*(ib))];
          bufy_b[r+rank_*(ia+a2*ib)] = 2.0*expo[1]*bufy[r+rank_*(ia+a2*(ib+1))] - ia*bufy[r+rank_*(ia+a2*(ib-1))];
          bufy_c[r+rank_*(ia+a2*ib)] = croots[r] * (xpPC2[1] * bufy[r+rank_*(ia+a2*ib)]
                                     + ia*bufy[r+rank_*(ia-1+a2*(ib))] + ib*bufy[r+rank_*(ia+a2*(ib-1))]);
          bufz_a[r+rank_*(ia+a2*ib)] = 2.0*expo[0]*bufz[r+rank_*(ia+1+a2*(ib))] - ia*bufz[r+rank_*(ia-1+a2*(ib))];
          bufz_b[r+rank_*(ia+a2*ib)] = 2.0*expo[1]*bufz[r+rank_*(ia+a2*(ib+1))] - ia*bufz[r+rank_*(ia+a2*(ib-1))];
          bufz_c[r+rank_*(ia+a2*ib)] = croots[r] * (xpPC2[2] * bufz[r+rank_*(ia+a2*ib)]
                                     + ia*bufz[r+rank_*(ia-1+a2*(ib))] + ib*bufz[r+rank_*(ia+a2*(ib-1))]);
        }
      }
    }
#endif

    // assembly step
    double* current_data0 = data0 + offset_iprim;
#if 0
    double* current_data1 = data1 + offset_iprim;
    double* current_data2 = data2 + offset_iprim;
    double* current_data3 = data3 + offset_iprim;
    double* current_data4 = data4 + offset_iprim;
    double* current_data5 = data5 + offset_iprim;
    double* current_data6 = data6 + offset_iprim;
    double* current_data7 = data7 + offset_iprim;
    double* current_data8 = data8 + offset_iprim;
#endif

    for (int iaz = 0; iaz <= a; ++iaz) { 
      for (int iay = 0; iay <= a - iaz; ++iay) { 
        const int iax = a - iaz - iay; 

        for (int ibz = 0; ibz <= b; ++ibz) { 
          for (int iby = 0; iby <= b - ibz; ++iby) { 
            const int ibx = b - ibz - iby;
            for (int r = 0; r != rank_; ++r) {
#define NAIDEBUG
#ifdef  NAIDEBUG
              *current_data0 += bufx[r+rank_*(iax+a2*ibx)] * bufy[r+rank_*(iay+a2*iby)] * bufz[r+rank_*(iaz+a2*ibz)];
#else
//            *current_data0 += bufx_a[r+rank_*(iax+a2*ibx)] * bufy[r+rank_*(iay+a2*iby)] * bufz[r+rank_*(iaz+a2*ibz)];
#endif
#if 0
              *current_data1 += bufx[r+rank_*(iax+a2*ibx)] * bufy_a[r+rank_*(iay+a2*iby)] * bufz[r+rank_*(iaz+a2*ibz)];
              *current_data2 += bufx[r+rank_*(iax+a2*ibx)] * bufy[r+rank_*(iay+a2*iby)] * bufz_a[r+rank_*(iaz+a2*ibz)];
              *current_data3 += bufx_b[r+rank_*(iax+a2*ibx)] * bufy[r+rank_*(iay+a2*iby)] * bufz[r+rank_*(iaz+a2*ibz)];
              *current_data4 += bufx[r+rank_*(iax+a2*ibx)] * bufy_b[r+rank_*(iay+a2*iby)] * bufz[r+rank_*(iaz+a2*ibz)];
              *current_data5 += bufx[r+rank_*(iax+a2*ibx)] * bufy[r+rank_*(iay+a2*iby)] * bufz_b[r+rank_*(iaz+a2*ibz)];
              *current_data6 += bufx_c[r+rank_*(iax+a2*ibx)] * bufy[r+rank_*(iay+a2*iby)] * bufz[r+rank_*(iaz+a2*ibz)];
              *current_data7 += bufx[r+rank_*(iax+a2*ibx)] * bufy_c[r+rank_*(iay+a2*iby)] * bufz[r+rank_*(iaz+a2*ibz)];
              *current_data8 += bufx[r+rank_*(iax+a2*ibx)] * bufy[r+rank_*(iay+a2*iby)] * bufz_c[r+rank_*(iaz+a2*ibz)];
#endif
            }
            ++current_data0;
#if 0
            ++current_data1;
            ++current_data2;
            ++current_data3;
            ++current_data4;
            ++current_data5;
            ++current_data6;
            ++current_data7;
            ++current_data8;
#endif
          }
        }
      }
    }
  }

  double* cdata = data_;
  for (int i = 0; i != 9; ++i, cdata += acpsize) {
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
    if (spherical_) {
      struct CarSphList carsphlist;
      const int carsphindex = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
      const int nloops = contsize_;
assert(contsize_ == cont0size_ * cont1size_);
      carsphlist.carsphfunc_call(carsphindex, nloops, source, target); 
    }

    // Sort cont01 and xyzab
    // data will be stored in data_: cont1b{ cont0a{ } }
    if (spherical_) {
      target = bkup_; 
      source = cdata;
      const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort_->sortfunc_call(index, target, source, cont1size_, cont0size_, 1, swap01_);
      copy(bkup_, bkup_+acpsize, cdata);
    } else {
      target = cdata; 
      source = bkup_;
      const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
      sort_->sortfunc_call(index, target, source, cont1size_, cont0size_, 1, swap01_);
    }
  }

  stack->release((amax_+1)*a2*b2*3);
  stack->release(rank_*a2*b2*12);
  stack->release(size_alloc_/12 + worksize*3);

  if (stack->get(0) != stack_save) throw logic_error("memory is not completely deallocated in GNAIBatch::vrr()");
}


