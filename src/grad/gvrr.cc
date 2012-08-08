//
// Newint - Parallel electron correlation program.
// Filename: gvrr.cc
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


#include <src/grad/gradbatch.h>
#include <src/rysint/int2d.h>
#include <src/stackmem.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <src/util/f77.h>
#include <src/util/comb.h>

using namespace std;

extern StackMem* stack;

static const Comb comb;

inline size_t GradBatch::m(int i, int a, int b, int c, int d) {
  const int la = basisinfo_[0]->angular_number()+2;
  const int lb = basisinfo_[1]->angular_number()+2;
  const int lc = basisinfo_[2]->angular_number()+2;
  const int ld = basisinfo_[3]->angular_number()+2;
  return i+rank_*(a+la*(b+lb*(c+lc*d)));
}

void GradBatch::perform_VRR() {
  double* const start = stack->get(0);
  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = rank_ * isize;
  const int vrr_index = amax_ * ANG_VRR_END + cmax_;

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

  const int a2 = a+2;
  const int b2 = b+2;
  const int c2 = c+2;
  const int d2 = d+2;

  double* const workx = stack->get(worksize*3);
  double* const worky = workx + worksize; 
  double* const workz = worky + worksize; 
  double* transx = stack->get((amax_+1)*a2*b2);
  double* transy = stack->get((amax_+1)*a2*b2);
  double* transz = stack->get((amax_+1)*a2*b2);
  double* trans2x = stack->get((cmax_+1)*c2*d2);
  double* trans2y = stack->get((cmax_+1)*c2*d2);
  double* trans2z = stack->get((cmax_+1)*c2*d2);
  fill(transx,  transx +(amax_+1)*a2*b2, 0.0);
  fill(transy,  transy +(amax_+1)*a2*b2, 0.0);
  fill(transz,  transz +(amax_+1)*a2*b2, 0.0);
  fill(trans2x, trans2x+(cmax_+1)*c2*d2, 0.0);
  fill(trans2y, trans2y+(cmax_+1)*c2*d2, 0.0);
  fill(trans2z, trans2z+(cmax_+1)*c2*d2, 0.0);
  // for usual integrals
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
  for (int id = 0, k = 0; id <= d+1; ++id) { 
    for (int ic = 0; ic <= c+1; ++ic, ++k) {
      if (ic == c+1 && id == d+1) continue;
      for (int i = ic; i <= ic+id; ++i) {
        trans2x[i + (cmax_+1)*k] = comb.c(id, ic+id-i) * pow(CD_[0], ic+id-i);
        trans2y[i + (cmax_+1)*k] = comb.c(id, ic+id-i) * pow(CD_[1], ic+id-i);
        trans2z[i + (cmax_+1)*k] = comb.c(id, ic+id-i) * pow(CD_[2], ic+id-i);
      }
    }
  }
  double* const intermediate = stack->get(b2*a2*(cmax_+1)*rank_);
  double* const final_x  = stack->get(b2*a2*c2*d2*rank_);
  double* const final_y  = stack->get(b2*a2*c2*d2*rank_);
  double* const final_z  = stack->get(b2*a2*c2*d2*rank_);
  double* const final_xa = stack->get(b2*a2*c2*d2*rank_);
  double* const final_xb = stack->get(b2*a2*c2*d2*rank_);
  double* const final_xc = stack->get(b2*a2*c2*d2*rank_);
  double* const final_xd = stack->get(b2*a2*c2*d2*rank_);
  double* const final_ya = stack->get(b2*a2*c2*d2*rank_);
  double* const final_yb = stack->get(b2*a2*c2*d2*rank_);
  double* const final_yc = stack->get(b2*a2*c2*d2*rank_);
  double* const final_yd = stack->get(b2*a2*c2*d2*rank_);
  double* const final_za = stack->get(b2*a2*c2*d2*rank_);
  double* const final_zb = stack->get(b2*a2*c2*d2*rank_);
  double* const final_zc = stack->get(b2*a2*c2*d2*rank_);
  double* const final_zd = stack->get(b2*a2*c2*d2*rank_);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];
    
    int offset = ii * rank_;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii;
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D cix(dparamx, &roots_[offset], rank_, worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", rank_, b2*a2, amax_+1, 1.0, workx+i*rank_*(amax_+1), rank_, transx, amax_+1, 0.0, intermediate+i*rank_*b2*a2, rank_);
    dgemm_("N", "N", rank_*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, rank_*b2*a2, trans2x, cmax_+1, 0.0, final_x, rank_*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciy(dparamy, &roots_[offset], rank_, worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", rank_, b2*a2, amax_+1, 1.0, worky+i*rank_*(amax_+1), rank_, transy, amax_+1, 0.0, intermediate+i*rank_*b2*a2, rank_);
    dgemm_("N", "N", rank_*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, rank_*b2*a2, trans2y, cmax_+1, 0.0, final_y, rank_*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D ciz(dparamz, &roots_[offset], rank_, worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", rank_, b2*a2, amax_+1, 1.0, workz+i*rank_*(amax_+1), rank_, transz, amax_+1, 0.0, intermediate+i*rank_*b2*a2, rank_);
    dgemm_("N", "N", rank_*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, rank_*b2*a2, trans2z, cmax_+1, 0.0, final_z, rank_*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != rank_; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m(r,ia+1,ib,ic,id)] - ia*final_x[m(r,ia-1,ib,ic,id)];
              final_xb[m(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m(r,ia,ib+1,ic,id)] - ib*final_x[m(r,ia,ib-1,ic,id)];
              final_xc[m(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m(r,ia,ib,ic+1,id)] - ic*final_x[m(r,ia,ib,ic-1,id)];
//            final_xd[m(r,ia,ib,ic,id)] = 2.0*expo[3]*final_x[m(r,ia,ib,ic,id+1)] - id*final_x[m(r,ia,ib,ic,id-1)];

              final_ya[m(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m(r,ia+1,ib,ic,id)] - ia*final_y[m(r,ia-1,ib,ic,id)];
              final_yb[m(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m(r,ia,ib+1,ic,id)] - ib*final_y[m(r,ia,ib-1,ic,id)];
              final_yc[m(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m(r,ia,ib,ic+1,id)] - ic*final_y[m(r,ia,ib,ic-1,id)];
//            final_yd[m(r,ia,ib,ic,id)] = 2.0*expo[3]*final_y[m(r,ia,ib,ic,id+1)] - id*final_y[m(r,ia,ib,ic,id-1)];

              final_za[m(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m(r,ia+1,ib,ic,id)] - ia*final_z[m(r,ia-1,ib,ic,id)];
              final_zb[m(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m(r,ia,ib+1,ic,id)] - ib*final_z[m(r,ia,ib-1,ic,id)];
              final_zc[m(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m(r,ia,ib,ic+1,id)] - ic*final_z[m(r,ia,ib,ic-1,id)];
//            final_zd[m(r,ia,ib,ic,id)] = 2.0*expo[3]*final_z[m(r,ia,ib,ic,id+1)] - id*final_z[m(r,ia,ib,ic,id-1)];
            }
          }
        }
      }
    }
 
    double* current_data0  = data_ + data_offset_ii;
    double* current_data1  = data_ + data_offset_ii + size_block_;
    double* current_data2  = data_ + data_offset_ii + size_block_* 2;
    double* current_data3  = data_ + data_offset_ii + size_block_* 3;
    double* current_data4  = data_ + data_offset_ii + size_block_* 4;
    double* current_data5  = data_ + data_offset_ii + size_block_* 5;
    double* current_data6  = data_ + data_offset_ii + size_block_* 6;
    double* current_data7  = data_ + data_offset_ii + size_block_* 7;
    double* current_data8  = data_ + data_offset_ii + size_block_* 8;
//  double* current_data9  = data_ + data_offset_ii + size_block_* 9;
//  double* current_data10 = data_ + data_offset_ii + size_block_*10;
//  double* current_data11 = data_ + data_offset_ii + size_block_*11;

    // CAUTION!
    // integrals in the 0(1(2(3(x2(x3(x0(x1))))))) order 
    for (int icz = 0; icz <= c; ++icz) { 
    for (int icy = 0; icy <= c - icz; ++icy) { 
    const int icx = c - icz - icy; 

      for (int idz = 0; idz <= d; ++idz) { 
      for (int idy = 0; idy <= d - idz; ++idy) { 
      const int idx = d - idz - idy; 

        for (int iaz = 0; iaz <= a; ++iaz) { 
        for (int iay = 0; iay <= a - iaz; ++iay) { 
        const int iax = a - iaz - iay; 

          for (int ibz = 0; ibz <= b; ++ibz) { 
          for (int iby = 0; iby <= b - ibz; ++iby) { 
          const int ibx = b - ibz - iby;

            for (int i = 0; i != rank_; ++i) {
              *current_data0  += final_xa[m(i, iax, ibx, icx, idx)] * final_y [m(i, iay, iby, icy, idy)] * final_z [m(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m(i, iax, ibx, icx, idx)] * final_ya[m(i, iay, iby, icy, idy)] * final_z [m(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m(i, iax, ibx, icx, idx)] * final_y [m(i, iay, iby, icy, idy)] * final_za[m(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m(i, iax, ibx, icx, idx)] * final_y [m(i, iay, iby, icy, idy)] * final_z [m(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m(i, iax, ibx, icx, idx)] * final_yb[m(i, iay, iby, icy, idy)] * final_z [m(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m(i, iax, ibx, icx, idx)] * final_y [m(i, iay, iby, icy, idy)] * final_zb[m(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m(i, iax, ibx, icx, idx)] * final_y [m(i, iay, iby, icy, idy)] * final_z [m(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m(i, iax, ibx, icx, idx)] * final_yc[m(i, iay, iby, icy, idy)] * final_z [m(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m(i, iax, ibx, icx, idx)] * final_y [m(i, iay, iby, icy, idy)] * final_zc[m(i, iaz, ibz, icz, idz)];
//            *current_data9  += final_xd[m(i, iax, ibx, icx, idx)] * final_y [m(i, iay, iby, icy, idy)] * final_z [m(i, iaz, ibz, icz, idz)];
//            *current_data10 += final_x [m(i, iax, ibx, icx, idx)] * final_yd[m(i, iay, iby, icy, idy)] * final_z [m(i, iaz, ibz, icz, idz)];
//            *current_data11 += final_x [m(i, iax, ibx, icx, idx)] * final_y [m(i, iay, iby, icy, idy)] * final_zd[m(i, iaz, ibz, icz, idz)];
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
//          ++current_data9;
//          ++current_data10;
//          ++current_data11;

          }}
        }}
      }}
    }}

  }

  stack->release(a2*b2*c2*d2*rank_ * 15);
  stack->release(a2*b2*(cmax_+1)*rank_);
  stack->release((cmax_+1)*c2*d2*3);
  stack->release((amax_+1)*b2*a2*3);
  stack->release(worksize*3);
  // checking memory leaks
  assert(start == stack->get(0));
}
