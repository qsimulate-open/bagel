//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _gvrr_drv.cc
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

#ifndef __SRC_RYSINT__GVRR_DRV_H
#define __SRC_RYSINT__GVRR_DRV_H

#include <numeric>
#include <algorithm>
#include <array>
#include <src/integral/rys/int2d.h>
#include <src/integral/rys/scaledata.h>
#include <src/util/f77.h>

namespace bagel {

template<int a_, int b_, int c_, int d_, int rank_>
void gvrr_driver(double* out, const double* const roots, const double* const weights, const double& coeff,
                 const std::array<double,3>& a, const std::array<double,3>& b, const std::array<double,3>& c, const std::array<double,3>& d,
                 const double* const p, const double* const q, const double& xp, const double& xq, const size_t& size_block,
                 const double* const expo, const double* const transx, const double* const transy, const double* const transz,
                 const double* const trans2x, const double* const trans2y, const double* const trans2z, double* const intermediate,
                 double* const final_x, double* const final_y, double* const final_z,
                 double* const final_xa, double* const final_xb, double* const final_xc,
                 double* const final_ya, double* const final_yb, double* const final_yc,
                 double* const final_za, double* const final_zb, double* const final_zc,
                 double* const workx, double* const worky, double* const workz, const std::array<bool,4>& dummy) {

  constexpr int amax_ = a_+b_+1;
  constexpr int cmax_ = c_+d_+1;

  constexpr int a2 = a_+2;
  constexpr int b2 = b_+2;
  constexpr int c2 = c_+2;
  constexpr int d2 = d_+2;

  constexpr int isize = (amax_+1) * (cmax_+1);
  constexpr int worksize = rank_ * isize;

  const double oxp2 = 0.5 / xp;
  const double oxq2 = 0.5 / xq;
  const double opq = 1.0 / (xp + xq);
  int2d<amax_,cmax_,rank_>(p[0], q[0], a[0], b[0], c[0], d[0], xp, xq, oxp2, oxq2, opq, roots, workx);
  scaledata<rank_, worksize>(workx, weights, coeff , workx);

  // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)
  for (int i = 0; i <= cmax_; ++i)
    dgemm_("N", "N", rank_, b2*a2, amax_+1, 1.0, workx+i*rank_*(amax_+1), rank_, transx, amax_+1, 0.0, intermediate+i*rank_*b2*a2, rank_);
  dgemm_("N", "N", rank_*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, rank_*b2*a2, trans2x, cmax_+1, 0.0, final_x, rank_*b2*a2);

  int2d<amax_,cmax_,rank_>(p[1], q[1], a[1], b[1], c[1], d[1], xp, xq, oxp2, oxq2, opq, roots, worky);

  for (int i = 0; i <= cmax_; ++i)
    dgemm_("N", "N", rank_, b2*a2, amax_+1, 1.0, worky+i*rank_*(amax_+1), rank_, transy, amax_+1, 0.0, intermediate+i*rank_*b2*a2, rank_);
  dgemm_("N", "N", rank_*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, rank_*b2*a2, trans2y, cmax_+1, 0.0, final_y, rank_*b2*a2);

  int2d<amax_,cmax_,rank_>(p[2], q[2], a[2], b[2], c[2], d[2], xp, xq, oxp2, oxq2, opq, roots, workz);

  for (int i = 0; i <= cmax_; ++i)
    dgemm_("N", "N", rank_, b2*a2, amax_+1, 1.0, workz+i*rank_*(amax_+1), rank_, transz, amax_+1, 0.0, intermediate+i*rank_*b2*a2, rank_);
  dgemm_("N", "N", rank_*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, rank_*b2*a2, trans2z, cmax_+1, 0.0, final_z, rank_*b2*a2);

  // if 3 is dummy, 2 is skipped
  if (!(dummy[2] || dummy[3]))
  for (int id = 0; id <= d_; ++id)
    for (int ic = 0; ic <= c_; ++ic)
      for (int ib = 0; ib <= b_; ++ib)
        for (int ia = 0; ia <= a_; ++ia)
          for (int r = 0; r != rank_; ++r) {
            final_xc[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[2]*final_x[r+rank_*(ia+a2*(ib+b2*(ic+1+c2*id)))] - (ic == 0 ? 0.0 : ic*final_x[r+rank_*(ia+a2*(ib+b2*(ic-1+c2*id)))]);
            final_yc[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[2]*final_y[r+rank_*(ia+a2*(ib+b2*(ic+1+c2*id)))] - (ic == 0 ? 0.0 : ic*final_y[r+rank_*(ia+a2*(ib+b2*(ic-1+c2*id)))]);
            final_zc[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[2]*final_z[r+rank_*(ia+a2*(ib+b2*(ic+1+c2*id)))] - (ic == 0 ? 0.0 : ic*final_z[r+rank_*(ia+a2*(ib+b2*(ic-1+c2*id)))]);
          }

  assert(!dummy[3] || !dummy[2]);
  if (!dummy[1])
  for (int id = 0; id <= d_; ++id)
    for (int ic = 0; ic <= c_; ++ic)
      for (int ib = 0; ib <= b_; ++ib)
        for (int ia = 0; ia <= a_; ++ia)
          for (int r = 0; r != rank_; ++r) {
            final_xb[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[1]*final_x[r+rank_*(ia+a2*(ib+1+b2*(ic+c2*id)))] - (ib == 0 ? 0.0 : ib*final_x[r+rank_*(ia+a2*(ib-1+b2*(ic+c2*id)))]);
            final_yb[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[1]*final_y[r+rank_*(ia+a2*(ib+1+b2*(ic+c2*id)))] - (ib == 0 ? 0.0 : ib*final_y[r+rank_*(ia+a2*(ib-1+b2*(ic+c2*id)))]);
            final_zb[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[1]*final_z[r+rank_*(ia+a2*(ib+1+b2*(ic+c2*id)))] - (ib == 0 ? 0.0 : ib*final_z[r+rank_*(ia+a2*(ib-1+b2*(ic+c2*id)))]);
          }

  if (!dummy[0])
  for (int id = 0; id <= d_; ++id)
    for (int ic = 0; ic <= c_; ++ic)
      for (int ib = 0; ib <= b_; ++ib)
        for (int ia = 0; ia <= a_; ++ia)
          for (int r = 0; r != rank_; ++r) {
            final_xa[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[0]*final_x[r+rank_*(ia+1+a2*(ib+b2*(ic+c2*id)))] - (ia == 0 ? 0.0 : ia*final_x[r+rank_*(ia-1+a2*(ib+b2*(ic+c2*id)))]);
            final_ya[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[0]*final_y[r+rank_*(ia+1+a2*(ib+b2*(ic+c2*id)))] - (ia == 0 ? 0.0 : ia*final_y[r+rank_*(ia-1+a2*(ib+b2*(ic+c2*id)))]);
            final_za[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[0]*final_z[r+rank_*(ia+1+a2*(ib+b2*(ic+c2*id)))] - (ia == 0 ? 0.0 : ia*final_z[r+rank_*(ia-1+a2*(ib+b2*(ic+c2*id)))]);
          }


  double* current_data0  = out;
  double* current_data1  = out + size_block;
  double* current_data2  = out + size_block* 2;
  double* current_data3  = out + size_block* 3;
  double* current_data4  = out + size_block* 4;
  double* current_data5  = out + size_block* 5;
  double* current_data6  = out + size_block* 6;
  double* current_data7  = out + size_block* 7;
  double* current_data8  = out + size_block* 8;

  // CAUTION!
  // integrals in the 0(1(2(3(x2(x3(x0(x1))))))) order
  if (!(dummy[2] || dummy[3]))
  for (int icz = 0; icz <= c_; ++icz) {
  for (int icy = 0; icy <= c_ - icz; ++icy) {
  const int icx = c_ - icz - icy;
    for (int idz = 0; idz <= d_; ++idz) {
    for (int idy = 0; idy <= d_ - idz; ++idy) {
    const int idx = d_ - idz - idy;
      for (int iaz = 0; iaz <= a_; ++iaz) {
      for (int iay = 0; iay <= a_ - iaz; ++iay) {
      const int iax = a_ - iaz - iay;
        for (int ibz = 0; ibz <= b_; ++ibz) {
        for (int iby = 0; iby <= b_ - ibz; ++iby) {
        const int ibx = b_ - ibz - iby;
          for (int i = 0; i != rank_; ++i) {
            *current_data6  += final_xc[i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_y [i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_z [i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];
            *current_data7  += final_x [i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_yc[i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_z [i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];
            *current_data8  += final_x [i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_y [i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_zc[i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];
          }
          ++current_data6;
          ++current_data7;
          ++current_data8;
        }}
      }}
    }}
  }}

  if (!dummy[1])
  for (int icz = 0; icz <= c_; ++icz) {
  for (int icy = 0; icy <= c_ - icz; ++icy) {
  const int icx = c_ - icz - icy;
    for (int idz = 0; idz <= d_; ++idz) {
    for (int idy = 0; idy <= d_ - idz; ++idy) {
    const int idx = d_ - idz - idy;
      for (int iaz = 0; iaz <= a_; ++iaz) {
      for (int iay = 0; iay <= a_ - iaz; ++iay) {
      const int iax = a_ - iaz - iay;
        for (int ibz = 0; ibz <= b_; ++ibz) {
        for (int iby = 0; iby <= b_ - ibz; ++iby) {
        const int ibx = b_ - ibz - iby;
          for (int i = 0; i != rank_; ++i) {
            *current_data3  += final_xb[i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_y [i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_z [i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];
            *current_data4  += final_x [i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_yb[i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_z [i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];
            *current_data5  += final_x [i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_y [i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_zb[i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];
          }
          ++current_data3;
          ++current_data4;
          ++current_data5;
        }}
      }}
    }}
  }}

  if (!dummy[0])
  for (int icz = 0; icz <= c_; ++icz) {
  for (int icy = 0; icy <= c_ - icz; ++icy) {
  const int icx = c_ - icz - icy;
    for (int idz = 0; idz <= d_; ++idz) {
    for (int idy = 0; idy <= d_ - idz; ++idy) {
    const int idx = d_ - idz - idy;
      for (int iaz = 0; iaz <= a_; ++iaz) {
      for (int iay = 0; iay <= a_ - iaz; ++iay) {
      const int iax = a_ - iaz - iay;
        for (int ibz = 0; ibz <= b_; ++ibz) {
        for (int iby = 0; iby <= b_ - ibz; ++iby) {
        const int ibx = b_ - ibz - iby;
          for (int i = 0; i != rank_; ++i) {
            *current_data0  += final_xa[i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_y [i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_z [i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];
            *current_data1  += final_x [i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_ya[i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_z [i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];
            *current_data2  += final_x [i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_y [i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_za[i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];
          }
          ++current_data0;
          ++current_data1;
          ++current_data2;
        }}
      }}
    }}
  }}
}

}

#endif
