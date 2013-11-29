//
// BAGEL - Parallel electron correlation program.
// Filename: _svrr_drv.h
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


// replaces generated codes _svrr_xxxx.cc etc

#ifndef __SRC_INTEGRAL_RYS___SVRR_DRIVER_H
#define __SRC_INTEGRAL_RYS___SVRR_DRIVER_H

#include <numeric>
#include <algorithm>
#include <array>
#include <src/integral/rys/int2d.h>
#include <src/integral/rys/scaledata.h>

namespace bagel {

template<int a_, int b_, int c_, int d_, int rank_>
void svrr_driver(double* out, const double* const roots, const double* const weights, const double& coeff,
                 const std::array<double,3>& a, const std::array<double,3>& b, const std::array<double,3>& c, const std::array<double,3>& d,
                 const double* const p, const double* const q, const double& xp, const double& xq,
                 const int* const amap, const int* const cmap, const int& asize_, double* const workx, double* const worky, double* const workz) {

  // compile time
  constexpr int amax_ = a_+b_;
  constexpr int cmax_ = c_+d_;
  constexpr int amax1_ = a_+b_+1;
  constexpr int cmax1_ = c_+d_+1;
  constexpr int amin_ = a_;
  constexpr int cmin_ = c_;

  constexpr int isize = (amax_ + 1) * (cmax_ + 1);
  constexpr int worksize = rank_ * isize;

  alignas(32) double iyiz[rank_];
  alignas(32) double womt[rank_];
  alignas(32) double wt[rank_];

  const double oxp2 = 0.5 / xp;
  const double oxq2 = 0.5 / xq;
  const double opq = 1.0 / (xp + xq);

  int2d<amax_,cmax_,rank_>(p[0], q[0], a[0], b[0], c[0], d[0], xp, xq, oxp2, oxq2, opq, roots, workx);
  for (int i = 0; i != rank_; ++i) {
    wt[i] = weights[i] * roots[i];
    womt[i] = weights[i] - wt[i];
  }
  scaledata<rank_, worksize>(workx, womt, coeff, workx);

  int2d<amax_,cmax_,rank_>(p[1], q[1], a[1], b[1], c[1], d[1], xp, xq, oxp2, oxq2, opq, roots, worky);
  int2d<amax_,cmax_,rank_>(p[2], q[2], a[2], b[2], c[2], d[2], xp, xq, oxp2, oxq2, opq, roots, workz);

  for (int iz = 0; iz <= cmax_; ++iz) {
    for (int iy = 0; iy <= cmax_ - iz; ++iy) {
      const int iyz = cmax1_ * (iy + cmax1_ * iz);
      for (int jz = 0; jz <= amax_; ++jz) {
        const int offsetz = rank_ * (amax1_ * iz + jz);
        for (int jy = 0; jy <= amax_ - jz; ++jy) {
          const int offsety = rank_ * (amax1_ * iy + jy);
          const int jyz = amax1_ * (jy + amax1_ * jz);
          for (int i = 0; i != rank_; ++i)
            iyiz[i] = worky[offsety + i] * workz[offsetz + i];
          for (int ix = std::max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
            const int iposition = cmap[ix + iyz];
            const int ipos_asize = iposition * asize_;
            for (int jx =std::max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
              const int offsetx = rank_ * (amax1_ * ix + jx);
              const int jposition = amap[jx + jyz];
              const int ijposition = jposition + ipos_asize;
              out[ijposition] = std::inner_product(iyiz, iyiz+rank_, workx+offsetx, 0.0);
            }
          }
        }
      }
    }
  }

}


template<int a_, int b_, int c_, int d_, int rank_>
void usvrr_driver(double* out, double* out2, const double* const roots, const double* const weights, const double& coeff, const double& coeffy,
                  const std::array<double,3>& a, const std::array<double,3>& b, const std::array<double,3>& c, const std::array<double,3>& d,
                  const double* const p, const double* const q, const double& xp, const double& xq,
                  const int* const amap, const int* const cmap, const int& asize_, double* const workx, double* const worky, double* const workz, double* const workx2) {

  // compile time
  constexpr int amax_ = a_+b_;
  constexpr int cmax_ = c_+d_;
  constexpr int amax1_ = a_+b_+1;
  constexpr int cmax1_ = c_+d_+1;
  constexpr int amin_ = a_;
  constexpr int cmin_ = c_;

  constexpr int isize = (amax_ + 1) * (cmax_ + 1);
  constexpr int worksize = rank_ * isize;

  alignas(32) double iyiz[rank_];
  alignas(32) double womt[rank_];
  alignas(32) double wt[rank_];

  const double oxp2 = 0.5 / xp;
  const double oxq2 = 0.5 / xq;
  const double opq = 1.0 / (xp + xq);

  int2d<amax_,cmax_,rank_>(p[0], q[0], a[0], b[0], c[0], d[0], xp, xq, oxp2, oxq2, opq, roots, workx);
  for (int i = 0; i != rank_; ++i) {
    wt[i] = weights[i] * roots[i];
    womt[i] = weights[i] - wt[i];
  }
  scaledata<rank_, worksize>(workx2, wt, coeffy, workx);
  scaledata<rank_, worksize>(workx, womt, coeff, workx);

  int2d<amax_,cmax_,rank_>(p[1], q[1], a[1], b[1], c[1], d[1], xp, xq, oxp2, oxq2, opq, roots, worky);
  int2d<amax_,cmax_,rank_>(p[2], q[2], a[2], b[2], c[2], d[2], xp, xq, oxp2, oxq2, opq, roots, workz);

  for (int iz = 0; iz <= cmax_; ++iz) {
    for (int iy = 0; iy <= cmax_ - iz; ++iy) {
      const int iyz = cmax1_ * (iy + cmax1_ * iz);
      for (int jz = 0; jz <= amax_; ++jz) {
        const int offsetz = rank_ * (amax1_ * iz + jz);
        for (int jy = 0; jy <= amax_ - jz; ++jy) {
          const int offsety = rank_ * (amax1_ * iy + jy);
          const int jyz = amax1_ * (jy + amax1_ * jz);
          for (int i = 0; i != rank_; ++i)
            iyiz[i] = worky[offsety + i] * workz[offsetz + i];
          for (int ix = std::max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
            const int iposition = cmap[ix + iyz];
            const int ipos_asize = iposition * asize_;
            for (int jx =std::max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
              const int offsetx = rank_ * (amax1_ * ix + jx);
              const int jposition = amap[jx + jyz];
              const int ijposition = jposition + ipos_asize;
              out[ijposition] = std::inner_product(iyiz, iyiz+rank_, workx+offsetx, 0.0);
              out2[ijposition] = std::inner_product(iyiz, iyiz+rank_, workx2+offsetx, 0.0);
            }
          }
        }
      }
    }
  }

}


}
#endif
