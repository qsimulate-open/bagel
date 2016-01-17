//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: _bvrr_drv.h
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

#ifndef __SRC_INTEGRAL_RYS____BVRR_DRIVER_H
#define __SRC_INTEGRAL_RYS____BVRR_DRIVER_H

#include <numeric>
#include <algorithm>
#include <array>
#include <src/integral/rys/int2d.h>
#include <src/integral/rys/scaledata.h>

namespace bagel {

template<int a_, int b_, int c_, int d_, int rank_>
void bvrr_driver(double* out, const double* const roots, const double* const weights, const double& coeff,
                const std::array<double,3>& a, const std::array<double,3>& b, const std::array<double,3>& c, const std::array<double,3>& d,
                const double* const p, const double* const q, const double& xp, const double& xq, const size_t& size_block,
                const int* const amap, const int* const cmap, const int& asize_, double* const workx, double* const worky, double* const workz,
                double* const worktx, double* const workty, double* const worktz, double* const worksx, double* const worksy, double* const worksz) {

  constexpr int amin_ = a_;
  constexpr int cmin_ = c_;
  constexpr int amax_ = a_+b_;
  constexpr int cmax_ = c_+d_;
  constexpr int amax1_ = amax_+1;
  constexpr int cmax1_ = cmax_+1;
  constexpr int amax2 = amax_+2;
  constexpr int cmax2 = cmax_+2;
  constexpr int isize = amax2 * cmax2;
  constexpr int worksize = rank_ * isize;

  double iyiz_nn[rank_];
  double iyiz_tn[rank_];
  double iyiz_nt[rank_];
  double iyiz_tt[rank_];
  double iyiz_sn[rank_];
  double iyiz_ns[rank_];

  const double oxp2 = 0.5 / xp;
  const double oxq2 = 0.5 / xq;
  const double opq = 1.0 / (xp + xq);

  int2d<amax1_,cmax1_,rank_>(p[0], q[0], a[0], b[0], c[0], d[0], xp, xq, oxp2, oxq2, opq, roots, workx);
  scaledata<rank_, worksize>(workx, weights, coeff*xp*xq*2.0*opq, workx);

  int2d<amax1_,cmax1_,rank_>(p[1], q[1], a[1], b[1], c[1], d[1], xp, xq, oxp2, oxq2, opq, roots, worky);
  int2d<amax1_,cmax1_,rank_>(p[2], q[2], a[2], b[2], c[2], d[2], xp, xq, oxp2, oxq2, opq, roots, workz);


  const double pq[3] = {p[0]-q[0], p[1]-q[1], p[2]-q[2]};
  const double ac[3] = {a[0]-c[0], a[1]-c[1], a[2]-c[2]};

  // next compute \\tidle{I}_x,y,z up to amax_, cmax_
  for (int ic = 0; ic <= cmax1_; ++ic)
    for (int ia = 0; ia <= amax1_; ++ia)
      for (int i = 0; i != rank_; ++i) {
        worktx[i+rank_*(ia+amax2*ic)] = pq[0]*workx[i+rank_*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workx[i+rank_*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workx[i+rank_*(ia+amax2*(ic-1))]);
        workty[i+rank_*(ia+amax2*ic)] = pq[1]*worky[i+rank_*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*worky[i+rank_*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*worky[i+rank_*(ia+amax2*(ic-1))]);
        worktz[i+rank_*(ia+amax2*ic)] = pq[2]*workz[i+rank_*(ia+amax2*ic)] + (ia==0 ? 0.0 : ia*oxp2*workz[i+rank_*(ia-1+amax2*ic)]) - (ic==0 ? 0.0 : ic*oxq2*workz[i+rank_*(ia+amax2*(ic-1))]);
      }
  // then compute \\tilde{\\tilde{I_x,y,z up to amax_-1, cmax_-1
  for (int ic = 0; ic != cmax1_; ++ic)
    for (int ia = 0; ia != amax1_; ++ia)
      for (int i = 0; i != rank_; ++i) {
        worksx[i+rank_*(ia+amax2*ic)] = worktx[i+rank_*((ia+1)+amax2*ic)] - worktx[i+rank_*(ia+amax2*(ic+1))] + ac[0]*worktx[i+rank_*(ia+amax2*ic)];
        worksy[i+rank_*(ia+amax2*ic)] = workty[i+rank_*((ia+1)+amax2*ic)] - workty[i+rank_*(ia+amax2*(ic+1))] + ac[1]*workty[i+rank_*(ia+amax2*ic)];
        worksz[i+rank_*(ia+amax2*ic)] = worktz[i+rank_*((ia+1)+amax2*ic)] - worktz[i+rank_*(ia+amax2*(ic+1))] + ac[2]*worktz[i+rank_*(ia+amax2*ic)];
      }

  double* const dataxx = out;
  double* const dataxy = dataxx + size_block;
  double* const dataxz = dataxy + size_block;
  double* const datayy = dataxz + size_block;
  double* const datayz = datayy + size_block;
  double* const datazz = datayz + size_block;

  // assemble up to amax_, cmax_
  for (int iz = 0; iz <= cmax_; ++iz) {
    for (int iy = 0; iy <= cmax_ - iz; ++iy) {
      for (int jz = 0; jz <= amax_; ++jz) {
        for (int jy = 0; jy <= amax_ - jz; ++jy) {
          const int offsetz = rank_ * (amax2 * iz + jz);
          const int offsety = rank_ * (amax2 * iy + jy);
          const int iyz = cmax1_ * (iy + cmax1_ * iz);
          const int jyz = amax1_ * (jy + amax1_ * jz);
          for (int i = 0; i != rank_; ++i) {
            iyiz_nn[i] = worky [offsety + i] * workz [offsetz + i];
            iyiz_tn[i] = workty[offsety + i] * workz [offsetz + i] * (1.0-roots[i]);
            iyiz_nt[i] = worky [offsety + i] * worktz[offsetz + i] * (1.0-roots[i]);
            iyiz_tt[i] = workty[offsety + i] * worktz[offsetz + i] * (1.0-roots[i]);
            iyiz_sn[i] = worksy[offsety + i] * workz [offsetz + i];
            iyiz_ns[i] = worky [offsety + i] * worksz[offsetz + i];
          }
          for (int ix = std::max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {
            const int iposition = cmap[ix + iyz];
            const int ipos_asize = iposition * asize_;
            for (int jx = std::max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {
              const int offsetx = rank_ * (amax2 * ix + jx);
              const int jposition = amap[jx + jyz];
              const int ijposition = jposition + ipos_asize;
              dataxx[ijposition] = blas::dot_product_noconj(iyiz_nn, rank_, worksx+offsetx);
              dataxy[ijposition] = blas::dot_product_noconj(iyiz_tn, rank_, worktx+offsetx);
              dataxz[ijposition] = blas::dot_product_noconj(iyiz_nt, rank_, worktx+offsetx);
              datayy[ijposition] = blas::dot_product_noconj(iyiz_sn, rank_, workx +offsetx);
              datayz[ijposition] = blas::dot_product_noconj(iyiz_tt, rank_, workx +offsetx);
              datazz[ijposition] = blas::dot_product_noconj(iyiz_ns, rank_, workx +offsetx);
            }
          }
        }
      }
    }
  }

}

}

#endif
