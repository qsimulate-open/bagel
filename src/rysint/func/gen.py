#!/usr/bin/python

import math

driver = "\
//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: _vrr_drv.h\n\
// Copyright (C) 2012 Toru Shiozaki\n\
//\n\
// Author: Toru Shiozaki <shiozaki@northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// The BAGEL package is free software; you can redistribute it and\/or modify\n\
// it under the terms of the GNU Library General Public License as published by\n\
// the Free Software Foundation; either version 2, or (at your option)\n\
// any later version.\n\
//\n\
// The BAGEL package is distributed in the hope that it will be useful,\n\
// but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
// GNU Library General Public License for more details.\n\
//\n\
// You should have received a copy of the GNU Library General Public License\n\
// along with the BAGEL package; see COPYING.  If not, write to\n\
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.\n\
//\n\
\n\
\n\
// replaces generated codes _vrr_xxxx.cc etc\n\
\n\
#ifndef __SRC_RYSINT____VRR_DRIVER_H\n\
#define __SRC_RYSINT____VRR_DRIVER_H\n\
\n\
#include <numeric>\n\
#include <algorithm>\n\
#include <array>\n\
#include <src/rysint/int2d.h>\n\
#include <src/rysint/scaledata.h>\n\
\n\
namespace bagel {\n\
\n\
template<int a_, int b_, int c_, int d_, int rank_>\n\
void vrr_driver(double* out, const double* const roots, const double* const weights, const double& coeff,\n\
                const std::array<double,3>& a, const std::array<double,3>& b, const std::array<double,3>& c, const std::array<double,3>& d,\n\
                const double* const p, const double* const q, const double& xp, const double& xq,\n\
                const int* const amap, const int* const cmap, const int& asize_, double* const workx, double* const worky, double* const workz) {\n\
\n\
  // compile time\n\
  const int amax_ = a_+b_;\n\
  const int cmax_ = c_+d_;\n\
  const int amax1_ = a_+b_+1;\n\
  const int cmax1_ = c_+d_+1;\n\
  const int amin_ = a_;\n\
  const int cmin_ = c_;\n\
\n\
  const int isize = (amax_ + 1) * (cmax_ + 1);\n\
  const int worksize = rank_ * isize;\n\
\n\
  double iyiz[rank_]__attribute__((aligned(32)));\n\
\n\
  const double oxp2 = 0.5 / xp;\n\
  const double oxq2 = 0.5 / xq;\n\
  const double opq = 1.0 / (xp + xq);\n\
\n\
  int2d<amax_,cmax_,rank_>(p[0], q[0], a[0], b[0], c[0], d[0], xp, xq, oxp2, oxq2, opq, roots, workx);\n\
  scaledata<rank_, worksize>(workx, weights, coeff, workx);\n\
\n\
  int2d<amax_,cmax_,rank_>(p[1], q[1], a[1], b[1], c[1], d[1], xp, xq, oxp2, oxq2, opq, roots, worky);\n\
  int2d<amax_,cmax_,rank_>(p[2], q[2], a[2], b[2], c[2], d[2], xp, xq, oxp2, oxq2, opq, roots, workz);\n\
\n\
  for (int iz = 0; iz <= cmax_; ++iz) {\n\
    for (int iy = 0; iy <= cmax_ - iz; ++iy) {\n\
      const int iyz = cmax1_ * (iy + cmax1_ * iz);\n\
      for (int jz = 0; jz <= amax_; ++jz) {\n\
        const int offsetz = rank_ * (amax1_ * iz + jz);\n\
        for (int jy = 0; jy <= amax_ - jz; ++jy) {\n\
          const int offsety = rank_ * (amax1_ * iy + jy);\n\
          const int jyz = amax1_ * (jy + amax1_ * jz);\n\
          for (int i = 0; i != rank_; ++i)\n\
            iyiz[i] = worky[offsety + i] * workz[offsetz + i];\n\
          for (int ix = std::max(0, cmin_ - iy - iz); ix <= cmax_ - iy - iz; ++ix) {\n\
            const int iposition = cmap[ix + iyz];\n\
            const int ipos_asize = iposition * asize_;\n\
            for (int jx =std::max(0, amin_ - jy - jz); jx <= amax_ - jy - jz; ++jx) {\n\
              const int offsetx = rank_ * (amax1_ * ix + jx);\n\
              const int jposition = amap[jx + jyz];\n\
              const int ijposition = jposition + ipos_asize;\n\
              out[ijposition] = std::inner_product(iyiz, iyiz+rank_, workx+offsetx, 0.0);\n\
            }\n\
          }\n\
        }\n\
      }\n\
    }\n\
  }\n\
\n\
}\n\
\n\
struct VRR_Driver {\n"

#in order to minimize the parsing time, there is no header
files = ""

for a in range(0,7):
 for b in range(0,7):
  if a < b: continue
  for c in range(0,7):
   for d in range(0,7):
    if c < d: continue
    rank = int(math.ceil((a+b+c+d+1)*0.5-0.001))

    filename = "vrr_driver_" + str(a) + "_" + str(b) + "_" + str(c) + "_" + str(d) + ".cc"
    fp = open(filename, "w")

    func_dec = "vrr_driver_" + str(a) + "_" + str(b) + "_" + str(c) + "_" + str(d) + "(double* out, const double* const roots, const double* const weights, const double& coeff,\n\
    const std::array<double,3>& a, const std::array<double,3>& b, const std::array<double,3>& c, const std::array<double,3>& d,\n\
    const double* const p, const double* const q, const double& xp, const double& xq, \n\
    const int* const amap, const int* const cmap, const int& asize_, double* const workx, double* const worky, double* const workz)"

    driver += "  static void " + func_dec + ";\n" 
    contents = "\
#include <src/rysint/_vrr_drv.h>\n\
void bagel::VRR_Driver::" + func_dec + " {\n\
  bagel::vrr_driver<" + str(a) + "," + str(b) + "," + str(c) + "," + str(d) + "," + str(rank) + ">(out, roots, weights, coeff, a, b, c, d, p, q, xp, xq, amap, cmap, asize_, workx, worky, workz);\n\
}\n"
    fp.write(contents)
    fp.close()
    files += filename + " "

driver += "};\n\
}\n\
#endif"

fp = open("_vrr_drv.h", "w")
fp.write(driver)

#print files
