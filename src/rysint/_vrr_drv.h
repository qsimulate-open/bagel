//
// BAGEL - Parallel electron correlation program.
// Filename: _vrr_drv.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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


// replaces generated codes _vrr_xxxx.cc etc

#ifndef __SRC_RYSINT____VRR_DRIVER_H
#define __SRC_RYSINT____VRR_DRIVER_H

#include <algorithm>
#include <array>
#include <src/rysint/eribatch.h>
#include <src/rysint/int2d.h>
#include <src/rysint/scaledata.h>

namespace bagel {

template<int amax_, int cmax_, int rank_>
void vrr_driver(double* out, const double* const roots, const double* const weights, const double& coeff,
                const std::array<double,3>& a, const std::array<double,3>& b, const std::array<double,3>& c, const std::array<double,3>& d,
                const double* const p, const double* const q, const double& xp, const double& xq,
                const int* const amap, const int* const cmap, const int& amin_, const int& cmin_, const int& asize_) {

  const int isize = (amax_ + 1) * (cmax_ + 1);
  const int worksize = rank_ * isize;

  double workx[worksize]__attribute__((aligned(32)));
  double worky[worksize]__attribute__((aligned(32)));
  double workz[worksize]__attribute__((aligned(32)));
  double iyiz[rank_]__attribute__((aligned(32)));

  const double oxp2 = 0.5 / xp;
  const double oxq2 = 0.5 / xq;
  const double opq = 1.0 / (xp + xq);

  int2d<amax_,cmax_,rank_>(std::array<double,11>{{p[0], q[0], a[0], b[0], c[0], d[0], xp, xq, oxp2, oxq2, opq}}, roots, workx);
  scaledata<rank_>(workx, weights, coeff, workx, worksize);

  int2d<amax_,cmax_,rank_>(std::array<double,11>{{p[1], q[1], a[1], b[1], c[1], d[1], xp, xq, oxp2, oxq2, opq}}, roots, worky);
  int2d<amax_,cmax_,rank_>(std::array<double,11>{{p[2], q[2], a[2], b[2], c[2], d[2], xp, xq, oxp2, oxq2, opq}}, roots, workz);

  const int amax1_ = amax_+1;
  const int cmax1_ = cmax_+1;

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
              out[ijposition] = 0.0;
              for (int k = 0; k != rank_; ++k)
                out[ijposition] += iyiz[k] * workx[offsetx+k];
            }
          }
        }
      }
    }
  }

}


}
#endif
