//
// BAGEL - Parallel electron correlation program.
// Filename: gvrr.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#include <src/grad/gradbatch.h>
#include <src/rysint/int2d.h>
#include <src/util/f77.h>
#include <src/util/comb.h>

using namespace std;
using namespace bagel;

static const Comb comb;

void GradBatch::perform_VRR1() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 1 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*1);
  double* const final_x  = stack_->get(b2*a2*c2*d2*1);
  double* const final_y  = stack_->get(b2*a2*c2*d2*1);
  double* const final_z  = stack_->get(b2*a2*c2*d2*1);
  double* const final_xa = stack_->get(b2*a2*c2*d2*1);
  double* const final_xb = stack_->get(b2*a2*c2*d2*1);
  double* const final_xc = stack_->get(b2*a2*c2*d2*1);
  double* const final_ya = stack_->get(b2*a2*c2*d2*1);
  double* const final_yb = stack_->get(b2*a2*c2*d2*1);
  double* const final_yc = stack_->get(b2*a2*c2*d2*1);
  double* const final_za = stack_->get(b2*a2*c2*d2*1);
  double* const final_zb = stack_->get(b2*a2*c2*d2*1);
  double* const final_zc = stack_->get(b2*a2*c2*d2*1);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 1;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<1> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 1, b2*a2, amax_+1, 1.0, workx+i*1*(amax_+1), 1, transx, amax_+1, 0.0, intermediate+i*1*b2*a2, 1);
    dgemm_("N", "N", 1*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 1*b2*a2, trans2x, cmax_+1, 0.0, final_x, 1*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<1> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 1, b2*a2, amax_+1, 1.0, worky+i*1*(amax_+1), 1, transy, amax_+1, 0.0, intermediate+i*1*b2*a2, 1);
    dgemm_("N", "N", 1*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 1*b2*a2, trans2y, cmax_+1, 0.0, final_y, 1*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<1> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 1, b2*a2, amax_+1, 1.0, workz+i*1*(amax_+1), 1, transz, amax_+1, 0.0, intermediate+i*1*b2*a2, 1);
    dgemm_("N", "N", 1*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 1*b2*a2, trans2z, cmax_+1, 0.0, final_z, 1*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 1; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<1>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<1>(r,ia+1,ib,ic,id)] - ia*final_x[m<1>(r,ia-1,ib,ic,id)];
              final_xb[m<1>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<1>(r,ia,ib+1,ic,id)] - ib*final_x[m<1>(r,ia,ib-1,ic,id)];
              final_xc[m<1>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<1>(r,ia,ib,ic+1,id)] - ic*final_x[m<1>(r,ia,ib,ic-1,id)];

              final_ya[m<1>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<1>(r,ia+1,ib,ic,id)] - ia*final_y[m<1>(r,ia-1,ib,ic,id)];
              final_yb[m<1>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<1>(r,ia,ib+1,ic,id)] - ib*final_y[m<1>(r,ia,ib-1,ic,id)];
              final_yc[m<1>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<1>(r,ia,ib,ic+1,id)] - ic*final_y[m<1>(r,ia,ib,ic-1,id)];

              final_za[m<1>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<1>(r,ia+1,ib,ic,id)] - ia*final_z[m<1>(r,ia-1,ib,ic,id)];
              final_zb[m<1>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<1>(r,ia,ib+1,ic,id)] - ib*final_z[m<1>(r,ia,ib-1,ic,id)];
              final_zc[m<1>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<1>(r,ia,ib,ic+1,id)] - ic*final_z[m<1>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 1; ++i) {
              *current_data0  += final_xa[m<1>(i, iax, ibx, icx, idx)] * final_y [m<1>(i, iay, iby, icy, idy)] * final_z [m<1>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<1>(i, iax, ibx, icx, idx)] * final_ya[m<1>(i, iay, iby, icy, idy)] * final_z [m<1>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<1>(i, iax, ibx, icx, idx)] * final_y [m<1>(i, iay, iby, icy, idy)] * final_za[m<1>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<1>(i, iax, ibx, icx, idx)] * final_y [m<1>(i, iay, iby, icy, idy)] * final_z [m<1>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<1>(i, iax, ibx, icx, idx)] * final_yb[m<1>(i, iay, iby, icy, idy)] * final_z [m<1>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<1>(i, iax, ibx, icx, idx)] * final_y [m<1>(i, iay, iby, icy, idy)] * final_zb[m<1>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<1>(i, iax, ibx, icx, idx)] * final_y [m<1>(i, iay, iby, icy, idy)] * final_z [m<1>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<1>(i, iax, ibx, icx, idx)] * final_yc[m<1>(i, iay, iby, icy, idy)] * final_z [m<1>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<1>(i, iax, ibx, icx, idx)] * final_y [m<1>(i, iay, iby, icy, idy)] * final_zc[m<1>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*1, final_zc);
  stack_->release(b2*a2*c2*d2*1, final_zb);
  stack_->release(b2*a2*c2*d2*1, final_za);
  stack_->release(b2*a2*c2*d2*1, final_yc);
  stack_->release(b2*a2*c2*d2*1, final_yb);
  stack_->release(b2*a2*c2*d2*1, final_ya);
  stack_->release(b2*a2*c2*d2*1, final_xc);
  stack_->release(b2*a2*c2*d2*1, final_xb);
  stack_->release(b2*a2*c2*d2*1, final_xa);
  stack_->release(b2*a2*c2*d2*1, final_z);
  stack_->release(b2*a2*c2*d2*1, final_y);
  stack_->release(b2*a2*c2*d2*1, final_x);

  stack_->release(b2*a2*(cmax_+1)*1, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

void GradBatch::perform_VRR2() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 2 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*2);
  double* const final_x  = stack_->get(b2*a2*c2*d2*2);
  double* const final_y  = stack_->get(b2*a2*c2*d2*2);
  double* const final_z  = stack_->get(b2*a2*c2*d2*2);
  double* const final_xa = stack_->get(b2*a2*c2*d2*2);
  double* const final_xb = stack_->get(b2*a2*c2*d2*2);
  double* const final_xc = stack_->get(b2*a2*c2*d2*2);
  double* const final_ya = stack_->get(b2*a2*c2*d2*2);
  double* const final_yb = stack_->get(b2*a2*c2*d2*2);
  double* const final_yc = stack_->get(b2*a2*c2*d2*2);
  double* const final_za = stack_->get(b2*a2*c2*d2*2);
  double* const final_zb = stack_->get(b2*a2*c2*d2*2);
  double* const final_zc = stack_->get(b2*a2*c2*d2*2);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 2;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<2> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 2, b2*a2, amax_+1, 1.0, workx+i*2*(amax_+1), 2, transx, amax_+1, 0.0, intermediate+i*2*b2*a2, 2);
    dgemm_("N", "N", 2*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 2*b2*a2, trans2x, cmax_+1, 0.0, final_x, 2*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<2> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 2, b2*a2, amax_+1, 1.0, worky+i*2*(amax_+1), 2, transy, amax_+1, 0.0, intermediate+i*2*b2*a2, 2);
    dgemm_("N", "N", 2*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 2*b2*a2, trans2y, cmax_+1, 0.0, final_y, 2*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<2> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 2, b2*a2, amax_+1, 1.0, workz+i*2*(amax_+1), 2, transz, amax_+1, 0.0, intermediate+i*2*b2*a2, 2);
    dgemm_("N", "N", 2*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 2*b2*a2, trans2z, cmax_+1, 0.0, final_z, 2*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 2; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<2>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<2>(r,ia+1,ib,ic,id)] - ia*final_x[m<2>(r,ia-1,ib,ic,id)];
              final_xb[m<2>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<2>(r,ia,ib+1,ic,id)] - ib*final_x[m<2>(r,ia,ib-1,ic,id)];
              final_xc[m<2>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<2>(r,ia,ib,ic+1,id)] - ic*final_x[m<2>(r,ia,ib,ic-1,id)];

              final_ya[m<2>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<2>(r,ia+1,ib,ic,id)] - ia*final_y[m<2>(r,ia-1,ib,ic,id)];
              final_yb[m<2>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<2>(r,ia,ib+1,ic,id)] - ib*final_y[m<2>(r,ia,ib-1,ic,id)];
              final_yc[m<2>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<2>(r,ia,ib,ic+1,id)] - ic*final_y[m<2>(r,ia,ib,ic-1,id)];

              final_za[m<2>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<2>(r,ia+1,ib,ic,id)] - ia*final_z[m<2>(r,ia-1,ib,ic,id)];
              final_zb[m<2>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<2>(r,ia,ib+1,ic,id)] - ib*final_z[m<2>(r,ia,ib-1,ic,id)];
              final_zc[m<2>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<2>(r,ia,ib,ic+1,id)] - ic*final_z[m<2>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 2; ++i) {
              *current_data0  += final_xa[m<2>(i, iax, ibx, icx, idx)] * final_y [m<2>(i, iay, iby, icy, idy)] * final_z [m<2>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<2>(i, iax, ibx, icx, idx)] * final_ya[m<2>(i, iay, iby, icy, idy)] * final_z [m<2>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<2>(i, iax, ibx, icx, idx)] * final_y [m<2>(i, iay, iby, icy, idy)] * final_za[m<2>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<2>(i, iax, ibx, icx, idx)] * final_y [m<2>(i, iay, iby, icy, idy)] * final_z [m<2>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<2>(i, iax, ibx, icx, idx)] * final_yb[m<2>(i, iay, iby, icy, idy)] * final_z [m<2>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<2>(i, iax, ibx, icx, idx)] * final_y [m<2>(i, iay, iby, icy, idy)] * final_zb[m<2>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<2>(i, iax, ibx, icx, idx)] * final_y [m<2>(i, iay, iby, icy, idy)] * final_z [m<2>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<2>(i, iax, ibx, icx, idx)] * final_yc[m<2>(i, iay, iby, icy, idy)] * final_z [m<2>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<2>(i, iax, ibx, icx, idx)] * final_y [m<2>(i, iay, iby, icy, idy)] * final_zc[m<2>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*2, final_zc);
  stack_->release(b2*a2*c2*d2*2, final_zb);
  stack_->release(b2*a2*c2*d2*2, final_za);
  stack_->release(b2*a2*c2*d2*2, final_yc);
  stack_->release(b2*a2*c2*d2*2, final_yb);
  stack_->release(b2*a2*c2*d2*2, final_ya);
  stack_->release(b2*a2*c2*d2*2, final_xc);
  stack_->release(b2*a2*c2*d2*2, final_xb);
  stack_->release(b2*a2*c2*d2*2, final_xa);
  stack_->release(b2*a2*c2*d2*2, final_z);
  stack_->release(b2*a2*c2*d2*2, final_y);
  stack_->release(b2*a2*c2*d2*2, final_x);

  stack_->release(b2*a2*(cmax_+1)*2, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

void GradBatch::perform_VRR3() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 3 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*3);
  double* const final_x  = stack_->get(b2*a2*c2*d2*3);
  double* const final_y  = stack_->get(b2*a2*c2*d2*3);
  double* const final_z  = stack_->get(b2*a2*c2*d2*3);
  double* const final_xa = stack_->get(b2*a2*c2*d2*3);
  double* const final_xb = stack_->get(b2*a2*c2*d2*3);
  double* const final_xc = stack_->get(b2*a2*c2*d2*3);
  double* const final_ya = stack_->get(b2*a2*c2*d2*3);
  double* const final_yb = stack_->get(b2*a2*c2*d2*3);
  double* const final_yc = stack_->get(b2*a2*c2*d2*3);
  double* const final_za = stack_->get(b2*a2*c2*d2*3);
  double* const final_zb = stack_->get(b2*a2*c2*d2*3);
  double* const final_zc = stack_->get(b2*a2*c2*d2*3);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 3;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<3> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 3, b2*a2, amax_+1, 1.0, workx+i*3*(amax_+1), 3, transx, amax_+1, 0.0, intermediate+i*3*b2*a2, 3);
    dgemm_("N", "N", 3*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 3*b2*a2, trans2x, cmax_+1, 0.0, final_x, 3*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<3> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 3, b2*a2, amax_+1, 1.0, worky+i*3*(amax_+1), 3, transy, amax_+1, 0.0, intermediate+i*3*b2*a2, 3);
    dgemm_("N", "N", 3*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 3*b2*a2, trans2y, cmax_+1, 0.0, final_y, 3*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<3> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 3, b2*a2, amax_+1, 1.0, workz+i*3*(amax_+1), 3, transz, amax_+1, 0.0, intermediate+i*3*b2*a2, 3);
    dgemm_("N", "N", 3*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 3*b2*a2, trans2z, cmax_+1, 0.0, final_z, 3*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 3; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<3>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<3>(r,ia+1,ib,ic,id)] - ia*final_x[m<3>(r,ia-1,ib,ic,id)];
              final_xb[m<3>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<3>(r,ia,ib+1,ic,id)] - ib*final_x[m<3>(r,ia,ib-1,ic,id)];
              final_xc[m<3>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<3>(r,ia,ib,ic+1,id)] - ic*final_x[m<3>(r,ia,ib,ic-1,id)];

              final_ya[m<3>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<3>(r,ia+1,ib,ic,id)] - ia*final_y[m<3>(r,ia-1,ib,ic,id)];
              final_yb[m<3>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<3>(r,ia,ib+1,ic,id)] - ib*final_y[m<3>(r,ia,ib-1,ic,id)];
              final_yc[m<3>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<3>(r,ia,ib,ic+1,id)] - ic*final_y[m<3>(r,ia,ib,ic-1,id)];

              final_za[m<3>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<3>(r,ia+1,ib,ic,id)] - ia*final_z[m<3>(r,ia-1,ib,ic,id)];
              final_zb[m<3>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<3>(r,ia,ib+1,ic,id)] - ib*final_z[m<3>(r,ia,ib-1,ic,id)];
              final_zc[m<3>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<3>(r,ia,ib,ic+1,id)] - ic*final_z[m<3>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 3; ++i) {
              *current_data0  += final_xa[m<3>(i, iax, ibx, icx, idx)] * final_y [m<3>(i, iay, iby, icy, idy)] * final_z [m<3>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<3>(i, iax, ibx, icx, idx)] * final_ya[m<3>(i, iay, iby, icy, idy)] * final_z [m<3>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<3>(i, iax, ibx, icx, idx)] * final_y [m<3>(i, iay, iby, icy, idy)] * final_za[m<3>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<3>(i, iax, ibx, icx, idx)] * final_y [m<3>(i, iay, iby, icy, idy)] * final_z [m<3>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<3>(i, iax, ibx, icx, idx)] * final_yb[m<3>(i, iay, iby, icy, idy)] * final_z [m<3>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<3>(i, iax, ibx, icx, idx)] * final_y [m<3>(i, iay, iby, icy, idy)] * final_zb[m<3>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<3>(i, iax, ibx, icx, idx)] * final_y [m<3>(i, iay, iby, icy, idy)] * final_z [m<3>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<3>(i, iax, ibx, icx, idx)] * final_yc[m<3>(i, iay, iby, icy, idy)] * final_z [m<3>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<3>(i, iax, ibx, icx, idx)] * final_y [m<3>(i, iay, iby, icy, idy)] * final_zc[m<3>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*3, final_zc);
  stack_->release(b2*a2*c2*d2*3, final_zb);
  stack_->release(b2*a2*c2*d2*3, final_za);
  stack_->release(b2*a2*c2*d2*3, final_yc);
  stack_->release(b2*a2*c2*d2*3, final_yb);
  stack_->release(b2*a2*c2*d2*3, final_ya);
  stack_->release(b2*a2*c2*d2*3, final_xc);
  stack_->release(b2*a2*c2*d2*3, final_xb);
  stack_->release(b2*a2*c2*d2*3, final_xa);
  stack_->release(b2*a2*c2*d2*3, final_z);
  stack_->release(b2*a2*c2*d2*3, final_y);
  stack_->release(b2*a2*c2*d2*3, final_x);

  stack_->release(b2*a2*(cmax_+1)*3, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

void GradBatch::perform_VRR4() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 4 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*4);
  double* const final_x  = stack_->get(b2*a2*c2*d2*4);
  double* const final_y  = stack_->get(b2*a2*c2*d2*4);
  double* const final_z  = stack_->get(b2*a2*c2*d2*4);
  double* const final_xa = stack_->get(b2*a2*c2*d2*4);
  double* const final_xb = stack_->get(b2*a2*c2*d2*4);
  double* const final_xc = stack_->get(b2*a2*c2*d2*4);
  double* const final_ya = stack_->get(b2*a2*c2*d2*4);
  double* const final_yb = stack_->get(b2*a2*c2*d2*4);
  double* const final_yc = stack_->get(b2*a2*c2*d2*4);
  double* const final_za = stack_->get(b2*a2*c2*d2*4);
  double* const final_zb = stack_->get(b2*a2*c2*d2*4);
  double* const final_zc = stack_->get(b2*a2*c2*d2*4);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 4;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<4> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 4, b2*a2, amax_+1, 1.0, workx+i*4*(amax_+1), 4, transx, amax_+1, 0.0, intermediate+i*4*b2*a2, 4);
    dgemm_("N", "N", 4*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 4*b2*a2, trans2x, cmax_+1, 0.0, final_x, 4*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<4> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 4, b2*a2, amax_+1, 1.0, worky+i*4*(amax_+1), 4, transy, amax_+1, 0.0, intermediate+i*4*b2*a2, 4);
    dgemm_("N", "N", 4*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 4*b2*a2, trans2y, cmax_+1, 0.0, final_y, 4*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<4> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 4, b2*a2, amax_+1, 1.0, workz+i*4*(amax_+1), 4, transz, amax_+1, 0.0, intermediate+i*4*b2*a2, 4);
    dgemm_("N", "N", 4*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 4*b2*a2, trans2z, cmax_+1, 0.0, final_z, 4*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 4; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<4>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<4>(r,ia+1,ib,ic,id)] - ia*final_x[m<4>(r,ia-1,ib,ic,id)];
              final_xb[m<4>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<4>(r,ia,ib+1,ic,id)] - ib*final_x[m<4>(r,ia,ib-1,ic,id)];
              final_xc[m<4>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<4>(r,ia,ib,ic+1,id)] - ic*final_x[m<4>(r,ia,ib,ic-1,id)];

              final_ya[m<4>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<4>(r,ia+1,ib,ic,id)] - ia*final_y[m<4>(r,ia-1,ib,ic,id)];
              final_yb[m<4>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<4>(r,ia,ib+1,ic,id)] - ib*final_y[m<4>(r,ia,ib-1,ic,id)];
              final_yc[m<4>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<4>(r,ia,ib,ic+1,id)] - ic*final_y[m<4>(r,ia,ib,ic-1,id)];

              final_za[m<4>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<4>(r,ia+1,ib,ic,id)] - ia*final_z[m<4>(r,ia-1,ib,ic,id)];
              final_zb[m<4>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<4>(r,ia,ib+1,ic,id)] - ib*final_z[m<4>(r,ia,ib-1,ic,id)];
              final_zc[m<4>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<4>(r,ia,ib,ic+1,id)] - ic*final_z[m<4>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 4; ++i) {
              *current_data0  += final_xa[m<4>(i, iax, ibx, icx, idx)] * final_y [m<4>(i, iay, iby, icy, idy)] * final_z [m<4>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<4>(i, iax, ibx, icx, idx)] * final_ya[m<4>(i, iay, iby, icy, idy)] * final_z [m<4>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<4>(i, iax, ibx, icx, idx)] * final_y [m<4>(i, iay, iby, icy, idy)] * final_za[m<4>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<4>(i, iax, ibx, icx, idx)] * final_y [m<4>(i, iay, iby, icy, idy)] * final_z [m<4>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<4>(i, iax, ibx, icx, idx)] * final_yb[m<4>(i, iay, iby, icy, idy)] * final_z [m<4>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<4>(i, iax, ibx, icx, idx)] * final_y [m<4>(i, iay, iby, icy, idy)] * final_zb[m<4>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<4>(i, iax, ibx, icx, idx)] * final_y [m<4>(i, iay, iby, icy, idy)] * final_z [m<4>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<4>(i, iax, ibx, icx, idx)] * final_yc[m<4>(i, iay, iby, icy, idy)] * final_z [m<4>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<4>(i, iax, ibx, icx, idx)] * final_y [m<4>(i, iay, iby, icy, idy)] * final_zc[m<4>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*4, final_zc);
  stack_->release(b2*a2*c2*d2*4, final_zb);
  stack_->release(b2*a2*c2*d2*4, final_za);
  stack_->release(b2*a2*c2*d2*4, final_yc);
  stack_->release(b2*a2*c2*d2*4, final_yb);
  stack_->release(b2*a2*c2*d2*4, final_ya);
  stack_->release(b2*a2*c2*d2*4, final_xc);
  stack_->release(b2*a2*c2*d2*4, final_xb);
  stack_->release(b2*a2*c2*d2*4, final_xa);
  stack_->release(b2*a2*c2*d2*4, final_z);
  stack_->release(b2*a2*c2*d2*4, final_y);
  stack_->release(b2*a2*c2*d2*4, final_x);

  stack_->release(b2*a2*(cmax_+1)*4, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

void GradBatch::perform_VRR5() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 5 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*5);
  double* const final_x  = stack_->get(b2*a2*c2*d2*5);
  double* const final_y  = stack_->get(b2*a2*c2*d2*5);
  double* const final_z  = stack_->get(b2*a2*c2*d2*5);
  double* const final_xa = stack_->get(b2*a2*c2*d2*5);
  double* const final_xb = stack_->get(b2*a2*c2*d2*5);
  double* const final_xc = stack_->get(b2*a2*c2*d2*5);
  double* const final_ya = stack_->get(b2*a2*c2*d2*5);
  double* const final_yb = stack_->get(b2*a2*c2*d2*5);
  double* const final_yc = stack_->get(b2*a2*c2*d2*5);
  double* const final_za = stack_->get(b2*a2*c2*d2*5);
  double* const final_zb = stack_->get(b2*a2*c2*d2*5);
  double* const final_zc = stack_->get(b2*a2*c2*d2*5);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 5;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<5> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 5, b2*a2, amax_+1, 1.0, workx+i*5*(amax_+1), 5, transx, amax_+1, 0.0, intermediate+i*5*b2*a2, 5);
    dgemm_("N", "N", 5*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 5*b2*a2, trans2x, cmax_+1, 0.0, final_x, 5*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<5> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 5, b2*a2, amax_+1, 1.0, worky+i*5*(amax_+1), 5, transy, amax_+1, 0.0, intermediate+i*5*b2*a2, 5);
    dgemm_("N", "N", 5*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 5*b2*a2, trans2y, cmax_+1, 0.0, final_y, 5*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<5> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 5, b2*a2, amax_+1, 1.0, workz+i*5*(amax_+1), 5, transz, amax_+1, 0.0, intermediate+i*5*b2*a2, 5);
    dgemm_("N", "N", 5*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 5*b2*a2, trans2z, cmax_+1, 0.0, final_z, 5*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 5; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<5>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<5>(r,ia+1,ib,ic,id)] - ia*final_x[m<5>(r,ia-1,ib,ic,id)];
              final_xb[m<5>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<5>(r,ia,ib+1,ic,id)] - ib*final_x[m<5>(r,ia,ib-1,ic,id)];
              final_xc[m<5>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<5>(r,ia,ib,ic+1,id)] - ic*final_x[m<5>(r,ia,ib,ic-1,id)];

              final_ya[m<5>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<5>(r,ia+1,ib,ic,id)] - ia*final_y[m<5>(r,ia-1,ib,ic,id)];
              final_yb[m<5>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<5>(r,ia,ib+1,ic,id)] - ib*final_y[m<5>(r,ia,ib-1,ic,id)];
              final_yc[m<5>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<5>(r,ia,ib,ic+1,id)] - ic*final_y[m<5>(r,ia,ib,ic-1,id)];

              final_za[m<5>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<5>(r,ia+1,ib,ic,id)] - ia*final_z[m<5>(r,ia-1,ib,ic,id)];
              final_zb[m<5>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<5>(r,ia,ib+1,ic,id)] - ib*final_z[m<5>(r,ia,ib-1,ic,id)];
              final_zc[m<5>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<5>(r,ia,ib,ic+1,id)] - ic*final_z[m<5>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 5; ++i) {
              *current_data0  += final_xa[m<5>(i, iax, ibx, icx, idx)] * final_y [m<5>(i, iay, iby, icy, idy)] * final_z [m<5>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<5>(i, iax, ibx, icx, idx)] * final_ya[m<5>(i, iay, iby, icy, idy)] * final_z [m<5>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<5>(i, iax, ibx, icx, idx)] * final_y [m<5>(i, iay, iby, icy, idy)] * final_za[m<5>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<5>(i, iax, ibx, icx, idx)] * final_y [m<5>(i, iay, iby, icy, idy)] * final_z [m<5>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<5>(i, iax, ibx, icx, idx)] * final_yb[m<5>(i, iay, iby, icy, idy)] * final_z [m<5>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<5>(i, iax, ibx, icx, idx)] * final_y [m<5>(i, iay, iby, icy, idy)] * final_zb[m<5>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<5>(i, iax, ibx, icx, idx)] * final_y [m<5>(i, iay, iby, icy, idy)] * final_z [m<5>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<5>(i, iax, ibx, icx, idx)] * final_yc[m<5>(i, iay, iby, icy, idy)] * final_z [m<5>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<5>(i, iax, ibx, icx, idx)] * final_y [m<5>(i, iay, iby, icy, idy)] * final_zc[m<5>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*5, final_zc);
  stack_->release(b2*a2*c2*d2*5, final_zb);
  stack_->release(b2*a2*c2*d2*5, final_za);
  stack_->release(b2*a2*c2*d2*5, final_yc);
  stack_->release(b2*a2*c2*d2*5, final_yb);
  stack_->release(b2*a2*c2*d2*5, final_ya);
  stack_->release(b2*a2*c2*d2*5, final_xc);
  stack_->release(b2*a2*c2*d2*5, final_xb);
  stack_->release(b2*a2*c2*d2*5, final_xa);
  stack_->release(b2*a2*c2*d2*5, final_z);
  stack_->release(b2*a2*c2*d2*5, final_y);
  stack_->release(b2*a2*c2*d2*5, final_x);

  stack_->release(b2*a2*(cmax_+1)*5, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

void GradBatch::perform_VRR6() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 6 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*6);
  double* const final_x  = stack_->get(b2*a2*c2*d2*6);
  double* const final_y  = stack_->get(b2*a2*c2*d2*6);
  double* const final_z  = stack_->get(b2*a2*c2*d2*6);
  double* const final_xa = stack_->get(b2*a2*c2*d2*6);
  double* const final_xb = stack_->get(b2*a2*c2*d2*6);
  double* const final_xc = stack_->get(b2*a2*c2*d2*6);
  double* const final_ya = stack_->get(b2*a2*c2*d2*6);
  double* const final_yb = stack_->get(b2*a2*c2*d2*6);
  double* const final_yc = stack_->get(b2*a2*c2*d2*6);
  double* const final_za = stack_->get(b2*a2*c2*d2*6);
  double* const final_zb = stack_->get(b2*a2*c2*d2*6);
  double* const final_zc = stack_->get(b2*a2*c2*d2*6);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 6;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<6> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 6, b2*a2, amax_+1, 1.0, workx+i*6*(amax_+1), 6, transx, amax_+1, 0.0, intermediate+i*6*b2*a2, 6);
    dgemm_("N", "N", 6*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 6*b2*a2, trans2x, cmax_+1, 0.0, final_x, 6*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<6> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 6, b2*a2, amax_+1, 1.0, worky+i*6*(amax_+1), 6, transy, amax_+1, 0.0, intermediate+i*6*b2*a2, 6);
    dgemm_("N", "N", 6*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 6*b2*a2, trans2y, cmax_+1, 0.0, final_y, 6*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<6> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 6, b2*a2, amax_+1, 1.0, workz+i*6*(amax_+1), 6, transz, amax_+1, 0.0, intermediate+i*6*b2*a2, 6);
    dgemm_("N", "N", 6*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 6*b2*a2, trans2z, cmax_+1, 0.0, final_z, 6*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 6; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<6>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<6>(r,ia+1,ib,ic,id)] - ia*final_x[m<6>(r,ia-1,ib,ic,id)];
              final_xb[m<6>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<6>(r,ia,ib+1,ic,id)] - ib*final_x[m<6>(r,ia,ib-1,ic,id)];
              final_xc[m<6>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<6>(r,ia,ib,ic+1,id)] - ic*final_x[m<6>(r,ia,ib,ic-1,id)];

              final_ya[m<6>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<6>(r,ia+1,ib,ic,id)] - ia*final_y[m<6>(r,ia-1,ib,ic,id)];
              final_yb[m<6>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<6>(r,ia,ib+1,ic,id)] - ib*final_y[m<6>(r,ia,ib-1,ic,id)];
              final_yc[m<6>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<6>(r,ia,ib,ic+1,id)] - ic*final_y[m<6>(r,ia,ib,ic-1,id)];

              final_za[m<6>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<6>(r,ia+1,ib,ic,id)] - ia*final_z[m<6>(r,ia-1,ib,ic,id)];
              final_zb[m<6>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<6>(r,ia,ib+1,ic,id)] - ib*final_z[m<6>(r,ia,ib-1,ic,id)];
              final_zc[m<6>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<6>(r,ia,ib,ic+1,id)] - ic*final_z[m<6>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 6; ++i) {
              *current_data0  += final_xa[m<6>(i, iax, ibx, icx, idx)] * final_y [m<6>(i, iay, iby, icy, idy)] * final_z [m<6>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<6>(i, iax, ibx, icx, idx)] * final_ya[m<6>(i, iay, iby, icy, idy)] * final_z [m<6>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<6>(i, iax, ibx, icx, idx)] * final_y [m<6>(i, iay, iby, icy, idy)] * final_za[m<6>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<6>(i, iax, ibx, icx, idx)] * final_y [m<6>(i, iay, iby, icy, idy)] * final_z [m<6>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<6>(i, iax, ibx, icx, idx)] * final_yb[m<6>(i, iay, iby, icy, idy)] * final_z [m<6>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<6>(i, iax, ibx, icx, idx)] * final_y [m<6>(i, iay, iby, icy, idy)] * final_zb[m<6>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<6>(i, iax, ibx, icx, idx)] * final_y [m<6>(i, iay, iby, icy, idy)] * final_z [m<6>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<6>(i, iax, ibx, icx, idx)] * final_yc[m<6>(i, iay, iby, icy, idy)] * final_z [m<6>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<6>(i, iax, ibx, icx, idx)] * final_y [m<6>(i, iay, iby, icy, idy)] * final_zc[m<6>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*6, final_zc);
  stack_->release(b2*a2*c2*d2*6, final_zb);
  stack_->release(b2*a2*c2*d2*6, final_za);
  stack_->release(b2*a2*c2*d2*6, final_yc);
  stack_->release(b2*a2*c2*d2*6, final_yb);
  stack_->release(b2*a2*c2*d2*6, final_ya);
  stack_->release(b2*a2*c2*d2*6, final_xc);
  stack_->release(b2*a2*c2*d2*6, final_xb);
  stack_->release(b2*a2*c2*d2*6, final_xa);
  stack_->release(b2*a2*c2*d2*6, final_z);
  stack_->release(b2*a2*c2*d2*6, final_y);
  stack_->release(b2*a2*c2*d2*6, final_x);

  stack_->release(b2*a2*(cmax_+1)*6, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

void GradBatch::perform_VRR7() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 7 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*7);
  double* const final_x  = stack_->get(b2*a2*c2*d2*7);
  double* const final_y  = stack_->get(b2*a2*c2*d2*7);
  double* const final_z  = stack_->get(b2*a2*c2*d2*7);
  double* const final_xa = stack_->get(b2*a2*c2*d2*7);
  double* const final_xb = stack_->get(b2*a2*c2*d2*7);
  double* const final_xc = stack_->get(b2*a2*c2*d2*7);
  double* const final_ya = stack_->get(b2*a2*c2*d2*7);
  double* const final_yb = stack_->get(b2*a2*c2*d2*7);
  double* const final_yc = stack_->get(b2*a2*c2*d2*7);
  double* const final_za = stack_->get(b2*a2*c2*d2*7);
  double* const final_zb = stack_->get(b2*a2*c2*d2*7);
  double* const final_zc = stack_->get(b2*a2*c2*d2*7);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 7;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<7> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 7, b2*a2, amax_+1, 1.0, workx+i*7*(amax_+1), 7, transx, amax_+1, 0.0, intermediate+i*7*b2*a2, 7);
    dgemm_("N", "N", 7*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 7*b2*a2, trans2x, cmax_+1, 0.0, final_x, 7*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<7> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 7, b2*a2, amax_+1, 1.0, worky+i*7*(amax_+1), 7, transy, amax_+1, 0.0, intermediate+i*7*b2*a2, 7);
    dgemm_("N", "N", 7*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 7*b2*a2, trans2y, cmax_+1, 0.0, final_y, 7*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<7> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 7, b2*a2, amax_+1, 1.0, workz+i*7*(amax_+1), 7, transz, amax_+1, 0.0, intermediate+i*7*b2*a2, 7);
    dgemm_("N", "N", 7*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 7*b2*a2, trans2z, cmax_+1, 0.0, final_z, 7*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 7; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<7>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<7>(r,ia+1,ib,ic,id)] - ia*final_x[m<7>(r,ia-1,ib,ic,id)];
              final_xb[m<7>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<7>(r,ia,ib+1,ic,id)] - ib*final_x[m<7>(r,ia,ib-1,ic,id)];
              final_xc[m<7>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<7>(r,ia,ib,ic+1,id)] - ic*final_x[m<7>(r,ia,ib,ic-1,id)];

              final_ya[m<7>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<7>(r,ia+1,ib,ic,id)] - ia*final_y[m<7>(r,ia-1,ib,ic,id)];
              final_yb[m<7>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<7>(r,ia,ib+1,ic,id)] - ib*final_y[m<7>(r,ia,ib-1,ic,id)];
              final_yc[m<7>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<7>(r,ia,ib,ic+1,id)] - ic*final_y[m<7>(r,ia,ib,ic-1,id)];

              final_za[m<7>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<7>(r,ia+1,ib,ic,id)] - ia*final_z[m<7>(r,ia-1,ib,ic,id)];
              final_zb[m<7>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<7>(r,ia,ib+1,ic,id)] - ib*final_z[m<7>(r,ia,ib-1,ic,id)];
              final_zc[m<7>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<7>(r,ia,ib,ic+1,id)] - ic*final_z[m<7>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 7; ++i) {
              *current_data0  += final_xa[m<7>(i, iax, ibx, icx, idx)] * final_y [m<7>(i, iay, iby, icy, idy)] * final_z [m<7>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<7>(i, iax, ibx, icx, idx)] * final_ya[m<7>(i, iay, iby, icy, idy)] * final_z [m<7>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<7>(i, iax, ibx, icx, idx)] * final_y [m<7>(i, iay, iby, icy, idy)] * final_za[m<7>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<7>(i, iax, ibx, icx, idx)] * final_y [m<7>(i, iay, iby, icy, idy)] * final_z [m<7>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<7>(i, iax, ibx, icx, idx)] * final_yb[m<7>(i, iay, iby, icy, idy)] * final_z [m<7>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<7>(i, iax, ibx, icx, idx)] * final_y [m<7>(i, iay, iby, icy, idy)] * final_zb[m<7>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<7>(i, iax, ibx, icx, idx)] * final_y [m<7>(i, iay, iby, icy, idy)] * final_z [m<7>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<7>(i, iax, ibx, icx, idx)] * final_yc[m<7>(i, iay, iby, icy, idy)] * final_z [m<7>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<7>(i, iax, ibx, icx, idx)] * final_y [m<7>(i, iay, iby, icy, idy)] * final_zc[m<7>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*7, final_zc);
  stack_->release(b2*a2*c2*d2*7, final_zb);
  stack_->release(b2*a2*c2*d2*7, final_za);
  stack_->release(b2*a2*c2*d2*7, final_yc);
  stack_->release(b2*a2*c2*d2*7, final_yb);
  stack_->release(b2*a2*c2*d2*7, final_ya);
  stack_->release(b2*a2*c2*d2*7, final_xc);
  stack_->release(b2*a2*c2*d2*7, final_xb);
  stack_->release(b2*a2*c2*d2*7, final_xa);
  stack_->release(b2*a2*c2*d2*7, final_z);
  stack_->release(b2*a2*c2*d2*7, final_y);
  stack_->release(b2*a2*c2*d2*7, final_x);

  stack_->release(b2*a2*(cmax_+1)*7, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

void GradBatch::perform_VRR8() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 8 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*8);
  double* const final_x  = stack_->get(b2*a2*c2*d2*8);
  double* const final_y  = stack_->get(b2*a2*c2*d2*8);
  double* const final_z  = stack_->get(b2*a2*c2*d2*8);
  double* const final_xa = stack_->get(b2*a2*c2*d2*8);
  double* const final_xb = stack_->get(b2*a2*c2*d2*8);
  double* const final_xc = stack_->get(b2*a2*c2*d2*8);
  double* const final_ya = stack_->get(b2*a2*c2*d2*8);
  double* const final_yb = stack_->get(b2*a2*c2*d2*8);
  double* const final_yc = stack_->get(b2*a2*c2*d2*8);
  double* const final_za = stack_->get(b2*a2*c2*d2*8);
  double* const final_zb = stack_->get(b2*a2*c2*d2*8);
  double* const final_zc = stack_->get(b2*a2*c2*d2*8);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 8;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<8> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 8, b2*a2, amax_+1, 1.0, workx+i*8*(amax_+1), 8, transx, amax_+1, 0.0, intermediate+i*8*b2*a2, 8);
    dgemm_("N", "N", 8*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 8*b2*a2, trans2x, cmax_+1, 0.0, final_x, 8*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<8> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 8, b2*a2, amax_+1, 1.0, worky+i*8*(amax_+1), 8, transy, amax_+1, 0.0, intermediate+i*8*b2*a2, 8);
    dgemm_("N", "N", 8*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 8*b2*a2, trans2y, cmax_+1, 0.0, final_y, 8*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<8> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 8, b2*a2, amax_+1, 1.0, workz+i*8*(amax_+1), 8, transz, amax_+1, 0.0, intermediate+i*8*b2*a2, 8);
    dgemm_("N", "N", 8*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 8*b2*a2, trans2z, cmax_+1, 0.0, final_z, 8*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 8; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<8>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<8>(r,ia+1,ib,ic,id)] - ia*final_x[m<8>(r,ia-1,ib,ic,id)];
              final_xb[m<8>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<8>(r,ia,ib+1,ic,id)] - ib*final_x[m<8>(r,ia,ib-1,ic,id)];
              final_xc[m<8>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<8>(r,ia,ib,ic+1,id)] - ic*final_x[m<8>(r,ia,ib,ic-1,id)];

              final_ya[m<8>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<8>(r,ia+1,ib,ic,id)] - ia*final_y[m<8>(r,ia-1,ib,ic,id)];
              final_yb[m<8>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<8>(r,ia,ib+1,ic,id)] - ib*final_y[m<8>(r,ia,ib-1,ic,id)];
              final_yc[m<8>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<8>(r,ia,ib,ic+1,id)] - ic*final_y[m<8>(r,ia,ib,ic-1,id)];

              final_za[m<8>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<8>(r,ia+1,ib,ic,id)] - ia*final_z[m<8>(r,ia-1,ib,ic,id)];
              final_zb[m<8>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<8>(r,ia,ib+1,ic,id)] - ib*final_z[m<8>(r,ia,ib-1,ic,id)];
              final_zc[m<8>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<8>(r,ia,ib,ic+1,id)] - ic*final_z[m<8>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 8; ++i) {
              *current_data0  += final_xa[m<8>(i, iax, ibx, icx, idx)] * final_y [m<8>(i, iay, iby, icy, idy)] * final_z [m<8>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<8>(i, iax, ibx, icx, idx)] * final_ya[m<8>(i, iay, iby, icy, idy)] * final_z [m<8>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<8>(i, iax, ibx, icx, idx)] * final_y [m<8>(i, iay, iby, icy, idy)] * final_za[m<8>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<8>(i, iax, ibx, icx, idx)] * final_y [m<8>(i, iay, iby, icy, idy)] * final_z [m<8>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<8>(i, iax, ibx, icx, idx)] * final_yb[m<8>(i, iay, iby, icy, idy)] * final_z [m<8>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<8>(i, iax, ibx, icx, idx)] * final_y [m<8>(i, iay, iby, icy, idy)] * final_zb[m<8>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<8>(i, iax, ibx, icx, idx)] * final_y [m<8>(i, iay, iby, icy, idy)] * final_z [m<8>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<8>(i, iax, ibx, icx, idx)] * final_yc[m<8>(i, iay, iby, icy, idy)] * final_z [m<8>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<8>(i, iax, ibx, icx, idx)] * final_y [m<8>(i, iay, iby, icy, idy)] * final_zc[m<8>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*8, final_zc);
  stack_->release(b2*a2*c2*d2*8, final_zb);
  stack_->release(b2*a2*c2*d2*8, final_za);
  stack_->release(b2*a2*c2*d2*8, final_yc);
  stack_->release(b2*a2*c2*d2*8, final_yb);
  stack_->release(b2*a2*c2*d2*8, final_ya);
  stack_->release(b2*a2*c2*d2*8, final_xc);
  stack_->release(b2*a2*c2*d2*8, final_xb);
  stack_->release(b2*a2*c2*d2*8, final_xa);
  stack_->release(b2*a2*c2*d2*8, final_z);
  stack_->release(b2*a2*c2*d2*8, final_y);
  stack_->release(b2*a2*c2*d2*8, final_x);

  stack_->release(b2*a2*(cmax_+1)*8, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

void GradBatch::perform_VRR9() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 9 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*9);
  double* const final_x  = stack_->get(b2*a2*c2*d2*9);
  double* const final_y  = stack_->get(b2*a2*c2*d2*9);
  double* const final_z  = stack_->get(b2*a2*c2*d2*9);
  double* const final_xa = stack_->get(b2*a2*c2*d2*9);
  double* const final_xb = stack_->get(b2*a2*c2*d2*9);
  double* const final_xc = stack_->get(b2*a2*c2*d2*9);
  double* const final_ya = stack_->get(b2*a2*c2*d2*9);
  double* const final_yb = stack_->get(b2*a2*c2*d2*9);
  double* const final_yc = stack_->get(b2*a2*c2*d2*9);
  double* const final_za = stack_->get(b2*a2*c2*d2*9);
  double* const final_zb = stack_->get(b2*a2*c2*d2*9);
  double* const final_zc = stack_->get(b2*a2*c2*d2*9);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 9;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<9> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 9, b2*a2, amax_+1, 1.0, workx+i*9*(amax_+1), 9, transx, amax_+1, 0.0, intermediate+i*9*b2*a2, 9);
    dgemm_("N", "N", 9*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 9*b2*a2, trans2x, cmax_+1, 0.0, final_x, 9*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<9> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 9, b2*a2, amax_+1, 1.0, worky+i*9*(amax_+1), 9, transy, amax_+1, 0.0, intermediate+i*9*b2*a2, 9);
    dgemm_("N", "N", 9*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 9*b2*a2, trans2y, cmax_+1, 0.0, final_y, 9*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<9> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 9, b2*a2, amax_+1, 1.0, workz+i*9*(amax_+1), 9, transz, amax_+1, 0.0, intermediate+i*9*b2*a2, 9);
    dgemm_("N", "N", 9*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 9*b2*a2, trans2z, cmax_+1, 0.0, final_z, 9*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 9; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<9>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<9>(r,ia+1,ib,ic,id)] - ia*final_x[m<9>(r,ia-1,ib,ic,id)];
              final_xb[m<9>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<9>(r,ia,ib+1,ic,id)] - ib*final_x[m<9>(r,ia,ib-1,ic,id)];
              final_xc[m<9>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<9>(r,ia,ib,ic+1,id)] - ic*final_x[m<9>(r,ia,ib,ic-1,id)];

              final_ya[m<9>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<9>(r,ia+1,ib,ic,id)] - ia*final_y[m<9>(r,ia-1,ib,ic,id)];
              final_yb[m<9>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<9>(r,ia,ib+1,ic,id)] - ib*final_y[m<9>(r,ia,ib-1,ic,id)];
              final_yc[m<9>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<9>(r,ia,ib,ic+1,id)] - ic*final_y[m<9>(r,ia,ib,ic-1,id)];

              final_za[m<9>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<9>(r,ia+1,ib,ic,id)] - ia*final_z[m<9>(r,ia-1,ib,ic,id)];
              final_zb[m<9>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<9>(r,ia,ib+1,ic,id)] - ib*final_z[m<9>(r,ia,ib-1,ic,id)];
              final_zc[m<9>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<9>(r,ia,ib,ic+1,id)] - ic*final_z[m<9>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 9; ++i) {
              *current_data0  += final_xa[m<9>(i, iax, ibx, icx, idx)] * final_y [m<9>(i, iay, iby, icy, idy)] * final_z [m<9>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<9>(i, iax, ibx, icx, idx)] * final_ya[m<9>(i, iay, iby, icy, idy)] * final_z [m<9>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<9>(i, iax, ibx, icx, idx)] * final_y [m<9>(i, iay, iby, icy, idy)] * final_za[m<9>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<9>(i, iax, ibx, icx, idx)] * final_y [m<9>(i, iay, iby, icy, idy)] * final_z [m<9>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<9>(i, iax, ibx, icx, idx)] * final_yb[m<9>(i, iay, iby, icy, idy)] * final_z [m<9>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<9>(i, iax, ibx, icx, idx)] * final_y [m<9>(i, iay, iby, icy, idy)] * final_zb[m<9>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<9>(i, iax, ibx, icx, idx)] * final_y [m<9>(i, iay, iby, icy, idy)] * final_z [m<9>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<9>(i, iax, ibx, icx, idx)] * final_yc[m<9>(i, iay, iby, icy, idy)] * final_z [m<9>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<9>(i, iax, ibx, icx, idx)] * final_y [m<9>(i, iay, iby, icy, idy)] * final_zc[m<9>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*9, final_zc);
  stack_->release(b2*a2*c2*d2*9, final_zb);
  stack_->release(b2*a2*c2*d2*9, final_za);
  stack_->release(b2*a2*c2*d2*9, final_yc);
  stack_->release(b2*a2*c2*d2*9, final_yb);
  stack_->release(b2*a2*c2*d2*9, final_ya);
  stack_->release(b2*a2*c2*d2*9, final_xc);
  stack_->release(b2*a2*c2*d2*9, final_xb);
  stack_->release(b2*a2*c2*d2*9, final_xa);
  stack_->release(b2*a2*c2*d2*9, final_z);
  stack_->release(b2*a2*c2*d2*9, final_y);
  stack_->release(b2*a2*c2*d2*9, final_x);

  stack_->release(b2*a2*(cmax_+1)*9, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

void GradBatch::perform_VRR10() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 10 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*10);
  double* const final_x  = stack_->get(b2*a2*c2*d2*10);
  double* const final_y  = stack_->get(b2*a2*c2*d2*10);
  double* const final_z  = stack_->get(b2*a2*c2*d2*10);
  double* const final_xa = stack_->get(b2*a2*c2*d2*10);
  double* const final_xb = stack_->get(b2*a2*c2*d2*10);
  double* const final_xc = stack_->get(b2*a2*c2*d2*10);
  double* const final_ya = stack_->get(b2*a2*c2*d2*10);
  double* const final_yb = stack_->get(b2*a2*c2*d2*10);
  double* const final_yc = stack_->get(b2*a2*c2*d2*10);
  double* const final_za = stack_->get(b2*a2*c2*d2*10);
  double* const final_zb = stack_->get(b2*a2*c2*d2*10);
  double* const final_zc = stack_->get(b2*a2*c2*d2*10);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 10;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<10> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 10, b2*a2, amax_+1, 1.0, workx+i*10*(amax_+1), 10, transx, amax_+1, 0.0, intermediate+i*10*b2*a2, 10);
    dgemm_("N", "N", 10*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 10*b2*a2, trans2x, cmax_+1, 0.0, final_x, 10*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<10> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 10, b2*a2, amax_+1, 1.0, worky+i*10*(amax_+1), 10, transy, amax_+1, 0.0, intermediate+i*10*b2*a2, 10);
    dgemm_("N", "N", 10*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 10*b2*a2, trans2y, cmax_+1, 0.0, final_y, 10*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<10> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 10, b2*a2, amax_+1, 1.0, workz+i*10*(amax_+1), 10, transz, amax_+1, 0.0, intermediate+i*10*b2*a2, 10);
    dgemm_("N", "N", 10*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 10*b2*a2, trans2z, cmax_+1, 0.0, final_z, 10*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 10; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<10>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<10>(r,ia+1,ib,ic,id)] - ia*final_x[m<10>(r,ia-1,ib,ic,id)];
              final_xb[m<10>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<10>(r,ia,ib+1,ic,id)] - ib*final_x[m<10>(r,ia,ib-1,ic,id)];
              final_xc[m<10>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<10>(r,ia,ib,ic+1,id)] - ic*final_x[m<10>(r,ia,ib,ic-1,id)];

              final_ya[m<10>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<10>(r,ia+1,ib,ic,id)] - ia*final_y[m<10>(r,ia-1,ib,ic,id)];
              final_yb[m<10>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<10>(r,ia,ib+1,ic,id)] - ib*final_y[m<10>(r,ia,ib-1,ic,id)];
              final_yc[m<10>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<10>(r,ia,ib,ic+1,id)] - ic*final_y[m<10>(r,ia,ib,ic-1,id)];

              final_za[m<10>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<10>(r,ia+1,ib,ic,id)] - ia*final_z[m<10>(r,ia-1,ib,ic,id)];
              final_zb[m<10>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<10>(r,ia,ib+1,ic,id)] - ib*final_z[m<10>(r,ia,ib-1,ic,id)];
              final_zc[m<10>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<10>(r,ia,ib,ic+1,id)] - ic*final_z[m<10>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 10; ++i) {
              *current_data0  += final_xa[m<10>(i, iax, ibx, icx, idx)] * final_y [m<10>(i, iay, iby, icy, idy)] * final_z [m<10>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<10>(i, iax, ibx, icx, idx)] * final_ya[m<10>(i, iay, iby, icy, idy)] * final_z [m<10>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<10>(i, iax, ibx, icx, idx)] * final_y [m<10>(i, iay, iby, icy, idy)] * final_za[m<10>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<10>(i, iax, ibx, icx, idx)] * final_y [m<10>(i, iay, iby, icy, idy)] * final_z [m<10>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<10>(i, iax, ibx, icx, idx)] * final_yb[m<10>(i, iay, iby, icy, idy)] * final_z [m<10>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<10>(i, iax, ibx, icx, idx)] * final_y [m<10>(i, iay, iby, icy, idy)] * final_zb[m<10>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<10>(i, iax, ibx, icx, idx)] * final_y [m<10>(i, iay, iby, icy, idy)] * final_z [m<10>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<10>(i, iax, ibx, icx, idx)] * final_yc[m<10>(i, iay, iby, icy, idy)] * final_z [m<10>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<10>(i, iax, ibx, icx, idx)] * final_y [m<10>(i, iay, iby, icy, idy)] * final_zc[m<10>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*10, final_zc);
  stack_->release(b2*a2*c2*d2*10, final_zb);
  stack_->release(b2*a2*c2*d2*10, final_za);
  stack_->release(b2*a2*c2*d2*10, final_yc);
  stack_->release(b2*a2*c2*d2*10, final_yb);
  stack_->release(b2*a2*c2*d2*10, final_ya);
  stack_->release(b2*a2*c2*d2*10, final_xc);
  stack_->release(b2*a2*c2*d2*10, final_xb);
  stack_->release(b2*a2*c2*d2*10, final_xa);
  stack_->release(b2*a2*c2*d2*10, final_z);
  stack_->release(b2*a2*c2*d2*10, final_y);
  stack_->release(b2*a2*c2*d2*10, final_x);

  stack_->release(b2*a2*(cmax_+1)*10, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

void GradBatch::perform_VRR11() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 11 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*11);
  double* const final_x  = stack_->get(b2*a2*c2*d2*11);
  double* const final_y  = stack_->get(b2*a2*c2*d2*11);
  double* const final_z  = stack_->get(b2*a2*c2*d2*11);
  double* const final_xa = stack_->get(b2*a2*c2*d2*11);
  double* const final_xb = stack_->get(b2*a2*c2*d2*11);
  double* const final_xc = stack_->get(b2*a2*c2*d2*11);
  double* const final_ya = stack_->get(b2*a2*c2*d2*11);
  double* const final_yb = stack_->get(b2*a2*c2*d2*11);
  double* const final_yc = stack_->get(b2*a2*c2*d2*11);
  double* const final_za = stack_->get(b2*a2*c2*d2*11);
  double* const final_zb = stack_->get(b2*a2*c2*d2*11);
  double* const final_zc = stack_->get(b2*a2*c2*d2*11);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 11;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<11> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 11, b2*a2, amax_+1, 1.0, workx+i*11*(amax_+1), 11, transx, amax_+1, 0.0, intermediate+i*11*b2*a2, 11);
    dgemm_("N", "N", 11*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 11*b2*a2, trans2x, cmax_+1, 0.0, final_x, 11*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<11> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 11, b2*a2, amax_+1, 1.0, worky+i*11*(amax_+1), 11, transy, amax_+1, 0.0, intermediate+i*11*b2*a2, 11);
    dgemm_("N", "N", 11*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 11*b2*a2, trans2y, cmax_+1, 0.0, final_y, 11*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<11> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 11, b2*a2, amax_+1, 1.0, workz+i*11*(amax_+1), 11, transz, amax_+1, 0.0, intermediate+i*11*b2*a2, 11);
    dgemm_("N", "N", 11*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 11*b2*a2, trans2z, cmax_+1, 0.0, final_z, 11*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 11; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<11>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<11>(r,ia+1,ib,ic,id)] - ia*final_x[m<11>(r,ia-1,ib,ic,id)];
              final_xb[m<11>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<11>(r,ia,ib+1,ic,id)] - ib*final_x[m<11>(r,ia,ib-1,ic,id)];
              final_xc[m<11>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<11>(r,ia,ib,ic+1,id)] - ic*final_x[m<11>(r,ia,ib,ic-1,id)];

              final_ya[m<11>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<11>(r,ia+1,ib,ic,id)] - ia*final_y[m<11>(r,ia-1,ib,ic,id)];
              final_yb[m<11>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<11>(r,ia,ib+1,ic,id)] - ib*final_y[m<11>(r,ia,ib-1,ic,id)];
              final_yc[m<11>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<11>(r,ia,ib,ic+1,id)] - ic*final_y[m<11>(r,ia,ib,ic-1,id)];

              final_za[m<11>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<11>(r,ia+1,ib,ic,id)] - ia*final_z[m<11>(r,ia-1,ib,ic,id)];
              final_zb[m<11>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<11>(r,ia,ib+1,ic,id)] - ib*final_z[m<11>(r,ia,ib-1,ic,id)];
              final_zc[m<11>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<11>(r,ia,ib,ic+1,id)] - ic*final_z[m<11>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 11; ++i) {
              *current_data0  += final_xa[m<11>(i, iax, ibx, icx, idx)] * final_y [m<11>(i, iay, iby, icy, idy)] * final_z [m<11>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<11>(i, iax, ibx, icx, idx)] * final_ya[m<11>(i, iay, iby, icy, idy)] * final_z [m<11>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<11>(i, iax, ibx, icx, idx)] * final_y [m<11>(i, iay, iby, icy, idy)] * final_za[m<11>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<11>(i, iax, ibx, icx, idx)] * final_y [m<11>(i, iay, iby, icy, idy)] * final_z [m<11>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<11>(i, iax, ibx, icx, idx)] * final_yb[m<11>(i, iay, iby, icy, idy)] * final_z [m<11>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<11>(i, iax, ibx, icx, idx)] * final_y [m<11>(i, iay, iby, icy, idy)] * final_zb[m<11>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<11>(i, iax, ibx, icx, idx)] * final_y [m<11>(i, iay, iby, icy, idy)] * final_z [m<11>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<11>(i, iax, ibx, icx, idx)] * final_yc[m<11>(i, iay, iby, icy, idy)] * final_z [m<11>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<11>(i, iax, ibx, icx, idx)] * final_y [m<11>(i, iay, iby, icy, idy)] * final_zc[m<11>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*11, final_zc);
  stack_->release(b2*a2*c2*d2*11, final_zb);
  stack_->release(b2*a2*c2*d2*11, final_za);
  stack_->release(b2*a2*c2*d2*11, final_yc);
  stack_->release(b2*a2*c2*d2*11, final_yb);
  stack_->release(b2*a2*c2*d2*11, final_ya);
  stack_->release(b2*a2*c2*d2*11, final_xc);
  stack_->release(b2*a2*c2*d2*11, final_xb);
  stack_->release(b2*a2*c2*d2*11, final_xa);
  stack_->release(b2*a2*c2*d2*11, final_z);
  stack_->release(b2*a2*c2*d2*11, final_y);
  stack_->release(b2*a2*c2*d2*11, final_x);

  stack_->release(b2*a2*(cmax_+1)*11, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

void GradBatch::perform_VRR12() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 12 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*12);
  double* const final_x  = stack_->get(b2*a2*c2*d2*12);
  double* const final_y  = stack_->get(b2*a2*c2*d2*12);
  double* const final_z  = stack_->get(b2*a2*c2*d2*12);
  double* const final_xa = stack_->get(b2*a2*c2*d2*12);
  double* const final_xb = stack_->get(b2*a2*c2*d2*12);
  double* const final_xc = stack_->get(b2*a2*c2*d2*12);
  double* const final_ya = stack_->get(b2*a2*c2*d2*12);
  double* const final_yb = stack_->get(b2*a2*c2*d2*12);
  double* const final_yc = stack_->get(b2*a2*c2*d2*12);
  double* const final_za = stack_->get(b2*a2*c2*d2*12);
  double* const final_zb = stack_->get(b2*a2*c2*d2*12);
  double* const final_zc = stack_->get(b2*a2*c2*d2*12);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 12;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<12> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 12, b2*a2, amax_+1, 1.0, workx+i*12*(amax_+1), 12, transx, amax_+1, 0.0, intermediate+i*12*b2*a2, 12);
    dgemm_("N", "N", 12*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 12*b2*a2, trans2x, cmax_+1, 0.0, final_x, 12*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<12> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 12, b2*a2, amax_+1, 1.0, worky+i*12*(amax_+1), 12, transy, amax_+1, 0.0, intermediate+i*12*b2*a2, 12);
    dgemm_("N", "N", 12*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 12*b2*a2, trans2y, cmax_+1, 0.0, final_y, 12*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<12> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 12, b2*a2, amax_+1, 1.0, workz+i*12*(amax_+1), 12, transz, amax_+1, 0.0, intermediate+i*12*b2*a2, 12);
    dgemm_("N", "N", 12*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 12*b2*a2, trans2z, cmax_+1, 0.0, final_z, 12*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 12; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<12>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<12>(r,ia+1,ib,ic,id)] - ia*final_x[m<12>(r,ia-1,ib,ic,id)];
              final_xb[m<12>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<12>(r,ia,ib+1,ic,id)] - ib*final_x[m<12>(r,ia,ib-1,ic,id)];
              final_xc[m<12>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<12>(r,ia,ib,ic+1,id)] - ic*final_x[m<12>(r,ia,ib,ic-1,id)];

              final_ya[m<12>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<12>(r,ia+1,ib,ic,id)] - ia*final_y[m<12>(r,ia-1,ib,ic,id)];
              final_yb[m<12>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<12>(r,ia,ib+1,ic,id)] - ib*final_y[m<12>(r,ia,ib-1,ic,id)];
              final_yc[m<12>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<12>(r,ia,ib,ic+1,id)] - ic*final_y[m<12>(r,ia,ib,ic-1,id)];

              final_za[m<12>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<12>(r,ia+1,ib,ic,id)] - ia*final_z[m<12>(r,ia-1,ib,ic,id)];
              final_zb[m<12>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<12>(r,ia,ib+1,ic,id)] - ib*final_z[m<12>(r,ia,ib-1,ic,id)];
              final_zc[m<12>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<12>(r,ia,ib,ic+1,id)] - ic*final_z[m<12>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 12; ++i) {
              *current_data0  += final_xa[m<12>(i, iax, ibx, icx, idx)] * final_y [m<12>(i, iay, iby, icy, idy)] * final_z [m<12>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<12>(i, iax, ibx, icx, idx)] * final_ya[m<12>(i, iay, iby, icy, idy)] * final_z [m<12>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<12>(i, iax, ibx, icx, idx)] * final_y [m<12>(i, iay, iby, icy, idy)] * final_za[m<12>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<12>(i, iax, ibx, icx, idx)] * final_y [m<12>(i, iay, iby, icy, idy)] * final_z [m<12>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<12>(i, iax, ibx, icx, idx)] * final_yb[m<12>(i, iay, iby, icy, idy)] * final_z [m<12>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<12>(i, iax, ibx, icx, idx)] * final_y [m<12>(i, iay, iby, icy, idy)] * final_zb[m<12>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<12>(i, iax, ibx, icx, idx)] * final_y [m<12>(i, iay, iby, icy, idy)] * final_z [m<12>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<12>(i, iax, ibx, icx, idx)] * final_yc[m<12>(i, iay, iby, icy, idy)] * final_z [m<12>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<12>(i, iax, ibx, icx, idx)] * final_y [m<12>(i, iay, iby, icy, idy)] * final_zc[m<12>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*12, final_zc);
  stack_->release(b2*a2*c2*d2*12, final_zb);
  stack_->release(b2*a2*c2*d2*12, final_za);
  stack_->release(b2*a2*c2*d2*12, final_yc);
  stack_->release(b2*a2*c2*d2*12, final_yb);
  stack_->release(b2*a2*c2*d2*12, final_ya);
  stack_->release(b2*a2*c2*d2*12, final_xc);
  stack_->release(b2*a2*c2*d2*12, final_xb);
  stack_->release(b2*a2*c2*d2*12, final_xa);
  stack_->release(b2*a2*c2*d2*12, final_z);
  stack_->release(b2*a2*c2*d2*12, final_y);
  stack_->release(b2*a2*c2*d2*12, final_x);

  stack_->release(b2*a2*(cmax_+1)*12, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

void GradBatch::perform_VRR13() {
  const int isize = (amax_ + 1) * (cmax_ + 1); 
  const int worksize = 13 * isize;
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

  double* const workx = stack_->get(worksize*3);
  double* const worky = workx + worksize;
  double* const workz = worky + worksize;
  double* const transx = stack_->get((amax_+1)*a2*b2);
  double* const transy = stack_->get((amax_+1)*a2*b2);
  double* const transz = stack_->get((amax_+1)*a2*b2);
  double* const trans2x = stack_->get((cmax_+1)*c2*d2);
  double* const trans2y = stack_->get((cmax_+1)*c2*d2);
  double* const trans2z = stack_->get((cmax_+1)*c2*d2);
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
  double* const intermediate = stack_->get(b2*a2*(cmax_+1)*13);
  double* const final_x  = stack_->get(b2*a2*c2*d2*13);
  double* const final_y  = stack_->get(b2*a2*c2*d2*13);
  double* const final_z  = stack_->get(b2*a2*c2*d2*13);
  double* const final_xa = stack_->get(b2*a2*c2*d2*13);
  double* const final_xb = stack_->get(b2*a2*c2*d2*13);
  double* const final_xc = stack_->get(b2*a2*c2*d2*13);
  double* const final_ya = stack_->get(b2*a2*c2*d2*13);
  double* const final_yb = stack_->get(b2*a2*c2*d2*13);
  double* const final_yc = stack_->get(b2*a2*c2*d2*13);
  double* const final_za = stack_->get(b2*a2*c2*d2*13);
  double* const final_zb = stack_->get(b2*a2*c2*d2*13);
  double* const final_zc = stack_->get(b2*a2*c2*d2*13);

  const int acsize = size_block_ / primsize_;
  assert(acsize == (a+1)*(b+1)*(c+1)*(d+1)*a2*b2*c2*d2/16 && size_block_*12 == size_alloc_);

  for (int j = 0; j != screening_size_; ++j) {
    const int ii = screening_[j];

    int offset = ii * 13;
    int data_offset_ii = ii * acsize;
    double* expo = exponents_.get() + ii*4;

    const int ii3 = 3 * ii; 
    const double cxp = xp_[ii];
    const double cxq = xq_[ii];
    const double oxp2 = 0.5 / cxp;
    const double oxq2 = 0.5 / cxq;
    const double opq = 1.0 / (cxp + cxq);
    const array<double, 11> dparamx = {{p_[ii3], q_[ii3], ax, bx, cx, dx, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<13> cix(dparamx, &roots_[offset], worksize, workx, vrr_->vrrfunc[vrr_index]);
    cix.scale_data(&weights_[offset], coeff_[ii]);

    // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 13, b2*a2, amax_+1, 1.0, workx+i*13*(amax_+1), 13, transx, amax_+1, 0.0, intermediate+i*13*b2*a2, 13);
    dgemm_("N", "N", 13*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 13*b2*a2, trans2x, cmax_+1, 0.0, final_x, 13*b2*a2);

    const array<double, 11> dparamy = {{p_[ii3 + 1], q_[ii3 + 1], ay, by, cy, dy, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<13> ciy(dparamy, &roots_[offset], worksize, worky, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 13, b2*a2, amax_+1, 1.0, worky+i*13*(amax_+1), 13, transy, amax_+1, 0.0, intermediate+i*13*b2*a2, 13);
    dgemm_("N", "N", 13*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 13*b2*a2, trans2y, cmax_+1, 0.0, final_y, 13*b2*a2);

    const array<double, 11> dparamz = {{p_[ii3 + 2], q_[ii3 + 2], az, bz, cz, dz, cxp, cxq, oxp2, oxq2, opq}};
    Int2D<13> ciz(dparamz, &roots_[offset], worksize, workz, vrr_->vrrfunc[vrr_index]);

    for (int i = 0; i <= cmax_; ++i)
      dgemm_("N", "N", 13, b2*a2, amax_+1, 1.0, workz+i*13*(amax_+1), 13, transz, amax_+1, 0.0, intermediate+i*13*b2*a2, 13);
    dgemm_("N", "N", 13*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, 13*b2*a2, trans2z, cmax_+1, 0.0, final_z, 13*b2*a2);

    for (int id = 0; id <= d; ++id) {
      for (int ic = 0; ic <= c; ++ic) {
        for (int ib = 0; ib <= b; ++ib) {
          for (int ia = 0; ia <= a; ++ia) {
            for (int r = 0; r != 13; ++r) {
                                                                                // v- this is a little dangerous, but perhaps the best
              final_xa[m<13>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_x[m<13>(r,ia+1,ib,ic,id)] - ia*final_x[m<13>(r,ia-1,ib,ic,id)];
              final_xb[m<13>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_x[m<13>(r,ia,ib+1,ic,id)] - ib*final_x[m<13>(r,ia,ib-1,ic,id)];
              final_xc[m<13>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_x[m<13>(r,ia,ib,ic+1,id)] - ic*final_x[m<13>(r,ia,ib,ic-1,id)];

              final_ya[m<13>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_y[m<13>(r,ia+1,ib,ic,id)] - ia*final_y[m<13>(r,ia-1,ib,ic,id)];
              final_yb[m<13>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_y[m<13>(r,ia,ib+1,ic,id)] - ib*final_y[m<13>(r,ia,ib-1,ic,id)];
              final_yc[m<13>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_y[m<13>(r,ia,ib,ic+1,id)] - ic*final_y[m<13>(r,ia,ib,ic-1,id)];

              final_za[m<13>(r,ia,ib,ic,id)] = 2.0*expo[0]*final_z[m<13>(r,ia+1,ib,ic,id)] - ia*final_z[m<13>(r,ia-1,ib,ic,id)];
              final_zb[m<13>(r,ia,ib,ic,id)] = 2.0*expo[1]*final_z[m<13>(r,ia,ib+1,ic,id)] - ib*final_z[m<13>(r,ia,ib-1,ic,id)];
              final_zc[m<13>(r,ia,ib,ic,id)] = 2.0*expo[2]*final_z[m<13>(r,ia,ib,ic+1,id)] - ic*final_z[m<13>(r,ia,ib,ic-1,id)];
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

            for (int i = 0; i != 13; ++i) {
              *current_data0  += final_xa[m<13>(i, iax, ibx, icx, idx)] * final_y [m<13>(i, iay, iby, icy, idy)] * final_z [m<13>(i, iaz, ibz, icz, idz)];
              *current_data1  += final_x [m<13>(i, iax, ibx, icx, idx)] * final_ya[m<13>(i, iay, iby, icy, idy)] * final_z [m<13>(i, iaz, ibz, icz, idz)];
              *current_data2  += final_x [m<13>(i, iax, ibx, icx, idx)] * final_y [m<13>(i, iay, iby, icy, idy)] * final_za[m<13>(i, iaz, ibz, icz, idz)];
              *current_data3  += final_xb[m<13>(i, iax, ibx, icx, idx)] * final_y [m<13>(i, iay, iby, icy, idy)] * final_z [m<13>(i, iaz, ibz, icz, idz)];
              *current_data4  += final_x [m<13>(i, iax, ibx, icx, idx)] * final_yb[m<13>(i, iay, iby, icy, idy)] * final_z [m<13>(i, iaz, ibz, icz, idz)];
              *current_data5  += final_x [m<13>(i, iax, ibx, icx, idx)] * final_y [m<13>(i, iay, iby, icy, idy)] * final_zb[m<13>(i, iaz, ibz, icz, idz)];
              *current_data6  += final_xc[m<13>(i, iax, ibx, icx, idx)] * final_y [m<13>(i, iay, iby, icy, idy)] * final_z [m<13>(i, iaz, ibz, icz, idz)];
              *current_data7  += final_x [m<13>(i, iax, ibx, icx, idx)] * final_yc[m<13>(i, iay, iby, icy, idy)] * final_z [m<13>(i, iaz, ibz, icz, idz)];
              *current_data8  += final_x [m<13>(i, iax, ibx, icx, idx)] * final_y [m<13>(i, iay, iby, icy, idy)] * final_zc[m<13>(i, iaz, ibz, icz, idz)];
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

          }}  
        }}  
      }}  
    }}  

  }

  stack_->release(b2*a2*c2*d2*13, final_zc);
  stack_->release(b2*a2*c2*d2*13, final_zb);
  stack_->release(b2*a2*c2*d2*13, final_za);
  stack_->release(b2*a2*c2*d2*13, final_yc);
  stack_->release(b2*a2*c2*d2*13, final_yb);
  stack_->release(b2*a2*c2*d2*13, final_ya);
  stack_->release(b2*a2*c2*d2*13, final_xc);
  stack_->release(b2*a2*c2*d2*13, final_xb);
  stack_->release(b2*a2*c2*d2*13, final_xa);
  stack_->release(b2*a2*c2*d2*13, final_z);
  stack_->release(b2*a2*c2*d2*13, final_y);
  stack_->release(b2*a2*c2*d2*13, final_x);

  stack_->release(b2*a2*(cmax_+1)*13, intermediate);

  stack_->release((cmax_+1)*c2*d2, trans2z);
  stack_->release((cmax_+1)*c2*d2, trans2y);
  stack_->release((cmax_+1)*c2*d2, trans2x);

  stack_->release((amax_+1)*a2*b2, transz);
  stack_->release((amax_+1)*a2*b2, transy);
  stack_->release((amax_+1)*a2*b2, transx);

  stack_->release(worksize*3, workx);
}

