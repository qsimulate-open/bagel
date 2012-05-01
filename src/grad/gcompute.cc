//
// Newint - Parallel electron correlation program.
// Filename: gcompute.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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

#include <cstring>
#include <src/stackmem.h>
#include <src/grad/gradbatch.h>
#include <src/fci/comb.h>
#include <src/util/f77.h>
#include <src/rysint/carsphlist.h>

using namespace std;

static const Comb comb;
extern StackMem* stack;

void GradBatch::compute() {

  bkup_ = stack->get(size_alloc_);
  bool swapped = false;

  const int ang0 = basisinfo_[0]->angular_number();
  const int ang1 = basisinfo_[1]->angular_number();
  const int ang2 = basisinfo_[2]->angular_number();
  const int ang3 = basisinfo_[3]->angular_number();
  const int absize_cart = (ang0+1) * (ang0+2) * (ang1+1) * (ang1+2) / 4;
  const int cdsize_cart = (ang2+1) * (ang2+2) * (ang3+1) * (ang3+2) / 4;
  const int absize_sph = spherical_ ? (2*ang0+1)*(2*ang1+1) : absize_cart;
  const int cdsize_sph = spherical_ ? (2*ang2+1)*(2*ang3+1) : cdsize_cart;

  assert(12* prim0size_ * prim1size_ * prim2size_ * prim3size_ * absize_cart * cdsize_cart == size_alloc_);

  // perform VRR
  // data_ will contain the intermediates: prim01{ prim23{ xyz{ } } } 
#if 0
  switch (rank_) {
    case 1: perform_VRR1(); break;
    case 2: perform_VRR2(); break;
    case 3: perform_VRR3(); break;
    case 4: perform_VRR4(); break;
    case 5: perform_VRR5(); break;
    case 6: perform_VRR6(); break;
    case 7: perform_VRR7(); break;
    case 8: perform_VRR8(); break;
    case 9: perform_VRR9(); break;  
    case 10: perform_VRR10(); break;  
    case 11: perform_VRR11(); break;  
    case 12: perform_VRR12(); break;  
    case 13: perform_VRR13(); break;  
    default: assert(false); break;
  }
#endif
  // CAUTION!
  // integrals in the 0(1(2(3(x2(x3(x0(x1))))))) order 
  perform_VRR();

  // contract indices 01 
  // data will be stored in bkup_: cont01{ prim23{ xyz{ } } }
  {
    const int m = prim2size_ * prim3size_ * absize_cart * cdsize_cart;
    perform_contraction_new_outer(m, data_, prim0size_, prim1size_, bkup_, 
      basisinfo_[0]->contractions(), basisinfo_[0]->contraction_upper(), basisinfo_[0]->contraction_lower(), cont0size_, 
      basisinfo_[1]->contractions(), basisinfo_[1]->contraction_upper(), basisinfo_[1]->contraction_lower(), cont1size_);
  }

  // contract indices 23 
  // data will be stored in data_: cont01{ cont23{ xyz{ } } }
  {
    const int n = cont0size_ * cont1size_;
    perform_contraction_new_inner(n, bkup_, prim2size_, prim3size_, data_, 
      basisinfo_[2]->contractions(), basisinfo_[2]->contraction_upper(), basisinfo_[2]->contraction_lower(), cont2size_, 
      basisinfo_[3]->contractions(), basisinfo_[3]->contraction_upper(), basisinfo_[3]->contraction_lower(), cont3size_);
  }

  // Cartesian to spherical 01 if necesarry
  // integrals in the 0(1(2(3(x2(x3(x0(x1))))))) order 
  struct CarSphList carsphlist;
  const bool need_sph01 = basisinfo_[0]->angular_number() > 1;
  if (spherical_ && need_sph01) {
    const int carsphindex = basisinfo_[0]->angular_number() * ANG_HRR_END + basisinfo_[1]->angular_number();
    const int nloops = contsize_ * cdsize_cart;
    carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_); 
  } else {
    swapped = (swapped ^ true); 
  }

  int a = (ang0 + 1) * (ang0 + 2) / 2;
  int b = (ang1 + 1) * (ang1 + 2) / 2;
  int c = (ang2 + 1) * (ang2 + 2) / 2;
  int d = (ang3 + 1) * (ang3 + 2) / 2;
  const int asph = 2 * ang0 + 1;
  const int bsph = 2 * ang1 + 1;
  const int csph = 2 * ang2 + 1;
  const int dsph = 2 * ang3 + 1;

  // the result will be 0(1(x0(x1(2(3(x2(x3)))))))
  if (basisinfo_[0]->angular_number() != 0) {
    const int m = spherical_ ? (asph * bsph) : (a * b);
    const int n = cont2size_ * cont3size_ * cdsize_cart;
    const int nloop = cont0size_ * cont1size_;
    int offset = 0;
    if (swapped) {
      for (int i = 0; i != nloop; ++i, offset += m * n)  
        mytranspose_(&data_[offset], &m, &n, &bkup_[offset]);
    } else {
      for (int i = 0; i != nloop; ++i, offset += m * n)  
        mytranspose_(&bkup_[offset], &m, &n, &data_[offset]);
    }
  } else {
    swapped = (swapped ^ true);
  }

  // Cartesian to spherical 23 if necesarry
  // data will be stored in bkup_
  // the result will be 0(1(x0(x1(2(3(x2(x3)))))))
  const bool need_sph23 = basisinfo_[2]->angular_number() > 1;
  if (spherical_ && need_sph23) {
    const int carsphindex = basisinfo_[2]->angular_number() * ANG_HRR_END + basisinfo_[3]->angular_number();
    const int nloops = contsize_ * asph * bsph;
    if (swapped)
      carsphlist.carsphfunc_call(carsphindex, nloops, bkup_, data_); 
    else
      carsphlist.carsphfunc_call(carsphindex, nloops, data_, bkup_); 
  } else {
    swapped = (swapped ^ true);
  }

  double *target_now = swapped ? bkup_ : data_;
  double *source_now = swapped ? data_ : bkup_;

  // Sort cont23 and xyzcd
  // data will be stored in data_: cont01{ xyzab{ cont3d{ cont2c{ } } } }
  if (basisinfo_[2]->angular_number() != 0) {
    const int nloop = a * b * cont0size_ * cont1size_;
    const unsigned int index = basisinfo_[3]->angular_number() * ANG_HRR_END + basisinfo_[2]->angular_number();
    sort_->sortfunc_call(index, target_now, source_now, cont3size_, cont2size_, nloop, swap23_);
  } else {
    swapped = (swapped ^ true);
  }

  target_now = swapped ? data_ : bkup_;
  source_now = swapped ? bkup_ : data_;
  // transpose batch
  // data will be stored in bkup_: cont3d{ cont2c{ cont01{ xyzab{ } } } } 
  if (!no_transpose_) {
    const int m = c * d * cont2size_ * cont3size_;
    const int n = a * b * cont0size_ * cont1size_; 
    mytranspose_(source_now, &m, &n, target_now);
  } else {
    swapped = (swapped ^ true);
  } 

  target_now = swapped ? bkup_ : data_;
  source_now = swapped ? data_ : bkup_;
  // Sort cont01 and xyzab
  // data will be stored in data_: cont3d{ cont2c{ cont1b{ cont0a{ } } } }
  if (basisinfo_[0]->angular_number() != 0) {
    const int nloop = c * d * cont2size_ * cont3size_;
    const unsigned int index = basisinfo_[1]->angular_number() * ANG_HRR_END + basisinfo_[0]->angular_number();
    sort_->sortfunc_call(index, target_now, source_now, cont1size_, cont0size_, nloop, swap01_);
  } else {
    swapped = (swapped ^ true);
  }

  if (swapped) ::memcpy(data_, bkup_, size_alloc_ * sizeof(double)); 

  stack->release(size_alloc_);
}


size_t GradBatch::m(int i, int a, int b, int c, int d) {
  const int la = basisinfo_[0]->angular_number()+2;
  const int lb = basisinfo_[1]->angular_number()+2;
  const int lc = basisinfo_[2]->angular_number()+2;
  const int ld = basisinfo_[3]->angular_number()+2;
  return i+rank_*(a+la*(b+lb*(c+lc*d)));
}

void GradBatch::perform_VRR() {
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = rank_ * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

  double* workx = stack->get(worksize*3);
  double* worky = workx + worksize; 
  double* workz = worky + worksize; 
  double iyiz[rank_];

  const double ax = basisinfo_[0]->position(0);
  const double ay = basisinfo_[0]->position(1);
  const double az = basisinfo_[0]->position(2);
  const double bx = basisinfo_[1]->position(0);
  const double by = basisinfo_[1]->position(1);
  const double bz = basisinfo_[1]->position(2);
  const double cx = basisinfo_[2]->position(0);
  const double cy = basisinfo_[2]->position(1);
  const double cz = basisinfo_[2]->position(2);
  const double dx = basisinfo_[3]->position(0);
  const double dy = basisinfo_[3]->position(1);
  const double dz = basisinfo_[3]->position(2);
  // transformation matrix
  const int a = basisinfo_[0]->angular_number();
  const int b = basisinfo_[1]->angular_number();
  const int c = basisinfo_[2]->angular_number();
  const int d = basisinfo_[3]->angular_number();
  assert(a+b+1 == amax_ && c+d+1 == cmax_); 

  double* trans  = stack->get((amax_+1)*(a+2)*(b+2));
  double* trans2 = stack->get((cmax_+1)*(c+2)*(d+2));
  fill(trans,  trans +(amax_+1)*(a+2)*(b+2), 0.0);
  fill(trans2, trans2+(cmax_+1)*(c+2)*(d+2), 0.0);
  double* tmp = trans;
  // for usual integrals
  for (int ib = 0, k = 0; ib <= b+1; ++ib) { 
    for (int ia = 0; ia <= a+1; ++ia, ++k) {
      if (ia == a+1 && ib == b+1) continue;
      for (int i = ia; i <= ia+ib; ++i) {
        trans[i + (amax_+1)*k] += comb.c(ib, ia+ib-i);
      }
    }
  }
  for (int id = 0, k = 0; id <= d+1; ++id) { 
    for (int ic = 0; ic <= c+1; ++ic, ++k) {
      if (ic == c+1 && id == d+1) continue;
      for (int i = ic; i <= ic+id; ++i) {
        trans2[i + (cmax_+1)*k] += comb.c(id, ic+id-i);
      }
    }
  }
  double* intermediate = stack->get(((b+2)*(a+2)-1)*(cmax_+1)*rank_);
  double* final_x  = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_xa = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_xb = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_xc = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_xd = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_y  = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_ya = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_yb = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_yc = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_yd = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_z  = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_za = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_zb = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_zc = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);
  double* final_zd = stack->get((b+2)*(a+2)*(c+2)*(d+2)*rank_);

  double* expo = exponents_.get();
  assert((size_alloc_/12)*12 == size_alloc_);
  const int acsize = size_alloc_ / 12 / primsize_;

  for (int j = 0; j != screening_size_; ++j, expo+=4) {
    const int ii = screening_[j];
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;

    double* current_data0 = &data_[data_offset_ii];
    double* current_data1 = &data_[data_offset_ii + acsize];
    double* current_data2 = &data_[data_offset_ii + acsize*2];
    double* current_data3 = &data_[data_offset_ii + acsize*3];
    double* current_data4 = &data_[data_offset_ii + acsize*4];
    double* current_data5 = &data_[data_offset_ii + acsize*5];
    double* current_data6 = &data_[data_offset_ii + acsize*6];
    double* current_data7 = &data_[data_offset_ii + acsize*7];
    double* current_data8 = &data_[data_offset_ii + acsize*8];
    double* current_data9 = &data_[data_offset_ii + acsize*9];
    double* current_data10 = &data_[data_offset_ii + acsize*10];
    double* current_data11 = &data_[data_offset_ii + acsize*11];

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const double dparamx[11] = {p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq};
    Int2D cix(dparamx, &roots_[offset], rank_, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1) (except a+1,b+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", rank_, (b+2)*(a+2), amax_+1, 1.0, workx+i*rank_*(amax_+1), rank_, trans, amax_+1, 0.0, intermediate+i*rank_*(b+2)*(a+2), rank_);
    dgemm_("N", "N", rank_*(b+2)*(a+2), (c+2)*(d+2), cmax_+1, 1.0, intermediate, rank_*(b+2)*(a+2), trans2, cmax_+1, 0.0, final_x, rank_*(b+2)*(a+2));

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != rank_; ++r) {
              final_xa[m(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m(r,ia+1,ib,ic,id)] - (ia == 0 ? 0.0 : ia*final_x[m(r,ia-1,ib,ic,id)]);
              final_xb[m(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m(r,ia,ib+1,ic,id)] - (ib == 0 ? 0.0 : ib*final_x[m(r,ia,ib-1,ic,id)]);
              final_xc[m(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m(r,ia,ib,ic+1,id)] - (ic == 0 ? 0.0 : ic*final_x[m(r,ia,ib,ic-1,id)]);
              final_xd[m(r,ia,ib,ic,id)] = 2.0*expo[3]*final_x[m(r,ia,ib,ic,id+1)] - (id == 0 ? 0.0 : id*final_x[m(r,ia,ib,ic,id-1)]);
            }
          }
        }
      }
    }
 
    const double dparamy[11] = {p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciy(dparamy, &roots_[offset], rank_, worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", rank_, (b+2)*(a+2), amax_+1, 1.0, worky+i*rank_*(amax_+1), rank_, trans, amax_+1, 0.0, intermediate+i*rank_*(b+2)*(a+2), rank_);
    dgemm_("N", "N", rank_*(b+2)*(a+2), (c+2)*(d+2), cmax_+1, 1.0, intermediate, rank_*(b+2)*(a+2), trans2, cmax_+1, 0.0, final_y, rank_*(b+2)*(a+2));

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != rank_; ++r) {
              final_ya[m(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m(r,ia+1,ib,ic,id)] - (ia == 0 ? 0.0 : ia*final_y[m(r,ia-1,ib,ic,id)]);
              final_yb[m(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m(r,ia,ib+1,ic,id)] - (ib == 0 ? 0.0 : ib*final_y[m(r,ia,ib-1,ic,id)]);
              final_yc[m(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m(r,ia,ib,ic+1,id)] - (ic == 0 ? 0.0 : ic*final_y[m(r,ia,ib,ic-1,id)]);
              final_yd[m(r,ia,ib,ic,id)] = 2.0*expo[3]*final_y[m(r,ia,ib,ic,id+1)] - (id == 0 ? 0.0 : id*final_y[m(r,ia,ib,ic,id-1)]);
            }
          }
        }
      }
    }
 
    const double dparamz[11] = {p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq};
    Int2D ciz(dparamz, &roots_[offset], rank_, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", rank_, (b+2)*(a+2), amax_+1, 1.0, workz+i*rank_*(amax_+1), rank_, trans, amax_+1, 0.0, intermediate+i*rank_*(b+2)*(a+2), rank_);
    dgemm_("N", "N", rank_*(b+2)*(a+2), (c+2)*(d+2), cmax_+1, 1.0, intermediate, rank_*(b+2)*(a+2), trans2, cmax_+1, 0.0, final_z, rank_*(b+2)*(a+2));

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != rank_; ++r) {
              final_za[m(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m(r,ia+1,ib,ic,id)] - (ia == 0 ? 0.0 : ia*final_z[m(r,ia-1,ib,ic,id)]);
              final_zb[m(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m(r,ia,ib+1,ic,id)] - (ib == 0 ? 0.0 : ib*final_z[m(r,ia,ib-1,ic,id)]);
              final_zc[m(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m(r,ia,ib,ic+1,id)] - (ic == 0 ? 0.0 : ic*final_z[m(r,ia,ib,ic-1,id)]);
              final_zd[m(r,ia,ib,ic,id)] = 2.0*expo[3]*final_z[m(r,ia,ib,ic,id+1)] - (id == 0 ? 0.0 : id*final_z[m(r,ia,ib,ic,id-1)]);
            }
          }
        }
      }
    }
 
#if 1
    // CAUTION!
    // integrals in the 0(1(2(3(x2(x3(x0(x1))))))) order 
    for (int icx = 0; icx <= c; ++icx) { 
    for (int icy = 0; icy <= c - icx; ++icy) { 
    const int icz = c - icx - icy; 
    for (int idx = 0; idx <= d; ++idx) { 
    for (int idy = 0; idy <= d - idx; ++idy) { 
    const int idz = d - idx - idy; 
    for (int iax = 0; iax <= a; ++iax) { 
    for (int iay = 0; iay <= a - iax; ++iay) { 
    const int iaz = a - iax - iay; 
    for (int ibx = 0; ibx <= b; ++ibx) { 
    for (int iby = 0; iby <= b - ibx; ++iby) { 
    const int ibz = b - ibx - iby; 

    *current_data0 = 0.0;
    *current_data1 = 0.0;
    *current_data2 = 0.0;
    *current_data3 = 0.0;
    *current_data4 = 0.0;
    *current_data5 = 0.0;
    *current_data6 = 0.0;
    *current_data7 = 0.0;
    *current_data8 = 0.0;
    *current_data9 = 0.0;
    *current_data10 = 0.0;
    *current_data11 = 0.0;
    for (int i = 0; i != rank_; ++i) {
      *current_data0 += final_xa[m(i, iax, ibx, icx, idx)] * final_y[m(i, iay, iby, icy, idy)] * final_z[m(i, iaz, ibz, icz, idz)];
      *current_data1 += final_x[m(i, iax, ibx, icx, idx)] * final_ya[m(i, iay, iby, icy, idy)] * final_z[m(i, iaz, ibz, icz, idz)];
      *current_data2 += final_x[m(i, iax, ibx, icx, idx)] * final_y[m(i, iay, iby, icy, idy)] * final_za[m(i, iaz, ibz, icz, idz)];
      *current_data3 += final_xb[m(i, iax, ibx, icx, idx)] * final_y[m(i, iay, iby, icy, idy)] * final_z[m(i, iaz, ibz, icz, idz)];
      *current_data4 += final_x[m(i, iax, ibx, icx, idx)] * final_yb[m(i, iay, iby, icy, idy)] * final_z[m(i, iaz, ibz, icz, idz)];
      *current_data5 += final_x[m(i, iax, ibx, icx, idx)] * final_y[m(i, iay, iby, icy, idy)] * final_zb[m(i, iaz, ibz, icz, idz)];
      *current_data6 += final_xc[m(i, iax, ibx, icx, idx)] * final_y[m(i, iay, iby, icy, idy)] * final_z[m(i, iaz, ibz, icz, idz)];
      *current_data7 += final_x[m(i, iax, ibx, icx, idx)] * final_yc[m(i, iay, iby, icy, idy)] * final_z[m(i, iaz, ibz, icz, idz)];
      *current_data8 += final_x[m(i, iax, ibx, icx, idx)] * final_y[m(i, iay, iby, icy, idy)] * final_zc[m(i, iaz, ibz, icz, idz)];
      *current_data9 += final_xd[m(i, iax, ibx, icx, idx)] * final_y[m(i, iay, iby, icy, idy)] * final_z[m(i, iaz, ibz, icz, idz)];
      *current_data10+= final_x[m(i, iax, ibx, icx, idx)] * final_yd[m(i, iay, iby, icy, idy)] * final_z[m(i, iaz, ibz, icz, idz)];
      *current_data11+= final_x[m(i, iax, ibx, icx, idx)] * final_y[m(i, iay, iby, icy, idy)] * final_zd[m(i, iaz, ibz, icz, idz)];
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
    ++current_data9;
    ++current_data10;
    ++current_data11;

    }} }}
    }} }}
#endif

  }

  stack->release((b+2)*(a+2)*(c+2)*(d+2)*rank_ * 15);
  stack->release(((b+2)*(a+2)-1)*(cmax_+1)*rank_);
  stack->release((cmax_+1)*(c+2)*(d+2));
  stack->release((amax_+1)*(b+2)*(a+2));
  stack->release(worksize*3);
}
