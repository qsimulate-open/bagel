#!/usr/bin/python

import math

driver = "\
//\n\
// BAGEL - Parallel electron correlation program.\n\
// Filename: _gvrr_drv.cc\n\
// Copyright (C) 2012 Toru Shiozaki\n\
//\n\
// Author: Toru Shiozaki <shiozaki@northwestern.edu>\n\
// Maintainer: Shiozaki group\n\
//\n\
// This file is part of the BAGEL package.\n\
//\n\
// The BAGEL package is free software; you can redistribute it and/or modify\n\
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
#ifndef __SRC_RYSINT__GVRR_DRV_H\n\
#define __SRC_RYSINT__GVRR_DRV_H\n\
\n\
#include <numeric>\n\
#include <algorithm>\n\
#include <array>\n\
#include <src/rysint/int2d.h>\n\
#include <src/rysint/scaledata.h>\n\
#include <src/util/f77.h>\n\
\n\
namespace bagel {\n\
\n\
template<int a_, int b_, int c_, int d_, int rank_>\n\
void gvrr_driver(double* out, const double* const roots, const double* const weights, const double& coeff,\n\
                 const std::array<double,3>& a, const std::array<double,3>& b, const std::array<double,3>& c, const std::array<double,3>& d,\n\
                 const double* const p, const double* const q, const double& xp, const double& xq, const size_t& size_block,\n\
                 const double* const expo, const double* const transx, const double* const transy, const double* const transz,\n\
                 const double* const trans2x, const double* const trans2y, const double* const trans2z, double* const intermediate,\n\
                 double* const final_x, double* const final_y, double* const final_z,\n\
                 double* const final_xa, double* const final_xb, double* const final_xc,\n\
                 double* const final_ya, double* const final_yb, double* const final_yc,\n\
                 double* const final_za, double* const final_zb, double* const final_zc,\n\
                 double* const workx, double* const worky, double* const workz) {\n\
\n\
  const int amax_ = a_+b_+1;\n\
  const int cmax_ = c_+d_+1;\n\
\n\
  const int a2 = a_+2;\n\
  const int b2 = b_+2;\n\
  const int c2 = c_+2;\n\
  const int d2 = d_+2;\n\
\n\
  const int isize = (amax_+1) * (cmax_+1);\n\
  const int worksize = rank_ * isize;\n\
\n\
  const double oxp2 = 0.5 / xp;\n\
  const double oxq2 = 0.5 / xq;\n\
  const double opq = 1.0 / (xp + xq);\n\
  int2d<amax_,cmax_,rank_>(p[0], q[0], a[0], b[0], c[0], d[0], xp, xq, oxp2, oxq2, opq, roots, workx);\n\
  scaledata<rank_, worksize>(workx, weights, coeff , workx);\n\
\n\
  // first (0:a+b, 0, 0:c+d, 0)-> (0:a+1, 0:b+1, 0:c+1, 0:d+1)\n\
  for (int i = 0; i <= cmax_; ++i)\n\
    dgemm_(\"N\", \"N\", rank_, b2*a2, amax_+1, 1.0, workx+i*rank_*(amax_+1), rank_, transx, amax_+1, 0.0, intermediate+i*rank_*b2*a2, rank_);\n\
  dgemm_(\"N\", \"N\", rank_*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, rank_*b2*a2, trans2x, cmax_+1, 0.0, final_x, rank_*b2*a2);\n\
\n\
  int2d<amax_,cmax_,rank_>(p[1], q[1], a[1], b[1], c[1], d[1], xp, xq, oxp2, oxq2, opq, roots, worky);\n\
\n\
  for (int i = 0; i <= cmax_; ++i)\n\
    dgemm_(\"N\", \"N\", rank_, b2*a2, amax_+1, 1.0, worky+i*rank_*(amax_+1), rank_, transy, amax_+1, 0.0, intermediate+i*rank_*b2*a2, rank_);\n\
  dgemm_(\"N\", \"N\", rank_*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, rank_*b2*a2, trans2y, cmax_+1, 0.0, final_y, rank_*b2*a2);\n\
\n\
  int2d<amax_,cmax_,rank_>(p[2], q[2], a[2], b[2], c[2], d[2], xp, xq, oxp2, oxq2, opq, roots, workz);\n\
\n\
  for (int i = 0; i <= cmax_; ++i)\n\
    dgemm_(\"N\", \"N\", rank_, b2*a2, amax_+1, 1.0, workz+i*rank_*(amax_+1), rank_, transz, amax_+1, 0.0, intermediate+i*rank_*b2*a2, rank_);\n\
  dgemm_(\"N\", \"N\", rank_*b2*a2, c2*d2, cmax_+1, 1.0, intermediate, rank_*b2*a2, trans2z, cmax_+1, 0.0, final_z, rank_*b2*a2);\n\
\n\
  for (int id = 0; id <= d_; ++id) {\n\
    for (int ic = 0; ic <= c_; ++ic) {\n\
      for (int ib = 0; ib <= b_; ++ib) {\n\
        for (int ia = 0; ia <= a_; ++ia) {\n\
          for (int r = 0; r != rank_; ++r) {\n\
                                                                              // v- this is a little dangerous, but perhaps the best\n\
            final_xa[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[0]*final_x[r+rank_*(ia+1+a2*(ib+b2*(ic+c2*id)))] - (ia == 0 ? 0.0 : ia*final_x[r+rank_*(ia-1+a2*(ib+b2*(ic+c2*id)))]);\n\
            final_xb[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[1]*final_x[r+rank_*(ia+a2*(ib+1+b2*(ic+c2*id)))] - (ib == 0 ? 0.0 : ib*final_x[r+rank_*(ia+a2*(ib-1+b2*(ic+c2*id)))]);\n\
            final_xc[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[2]*final_x[r+rank_*(ia+a2*(ib+b2*(ic+1+c2*id)))] - (ic == 0 ? 0.0 : ic*final_x[r+rank_*(ia+a2*(ib+b2*(ic-1+c2*id)))]);\n\
\n\
            final_ya[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[0]*final_y[r+rank_*(ia+1+a2*(ib+b2*(ic+c2*id)))] - (ia == 0 ? 0.0 : ia*final_y[r+rank_*(ia-1+a2*(ib+b2*(ic+c2*id)))]);\n\
            final_yb[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[1]*final_y[r+rank_*(ia+a2*(ib+1+b2*(ic+c2*id)))] - (ib == 0 ? 0.0 : ib*final_y[r+rank_*(ia+a2*(ib-1+b2*(ic+c2*id)))]);\n\
            final_yc[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[2]*final_y[r+rank_*(ia+a2*(ib+b2*(ic+1+c2*id)))] - (ic == 0 ? 0.0 : ic*final_y[r+rank_*(ia+a2*(ib+b2*(ic-1+c2*id)))]);\n\
\n\
            final_za[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[0]*final_z[r+rank_*(ia+1+a2*(ib+b2*(ic+c2*id)))] - (ia == 0 ? 0.0 : ia*final_z[r+rank_*(ia-1+a2*(ib+b2*(ic+c2*id)))]);\n\
            final_zb[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[1]*final_z[r+rank_*(ia+a2*(ib+1+b2*(ic+c2*id)))] - (ib == 0 ? 0.0 : ib*final_z[r+rank_*(ia+a2*(ib-1+b2*(ic+c2*id)))]);\n\
            final_zc[r+rank_*(ia+a2*(ib+b2*(ic+c2*id)))] = 2.0*expo[2]*final_z[r+rank_*(ia+a2*(ib+b2*(ic+1+c2*id)))] - (ic == 0 ? 0.0 : ic*final_z[r+rank_*(ia+a2*(ib+b2*(ic-1+c2*id)))]);\n\
          }\n\
        }\n\
      }\n\
    }\n\
  }\n\
\n\
  double* current_data0  = out;\n\
  double* current_data1  = out + size_block;\n\
  double* current_data2  = out + size_block* 2;\n\
  double* current_data3  = out + size_block* 3;\n\
  double* current_data4  = out + size_block* 4;\n\
  double* current_data5  = out + size_block* 5;\n\
  double* current_data6  = out + size_block* 6;\n\
  double* current_data7  = out + size_block* 7;\n\
  double* current_data8  = out + size_block* 8;\n\
\n\
  // CAUTION!\n\
  // integrals in the 0(1(2(3(x2(x3(x0(x1))))))) order\n\
  for (int icz = 0; icz <= c_; ++icz) {\n\
  for (int icy = 0; icy <= c_ - icz; ++icy) {\n\
  const int icx = c_ - icz - icy;\n\
\n\
    for (int idz = 0; idz <= d_; ++idz) {\n\
    for (int idy = 0; idy <= d_ - idz; ++idy) {\n\
    const int idx = d_ - idz - idy;\n\
\n\
      for (int iaz = 0; iaz <= a_; ++iaz) {\n\
      for (int iay = 0; iay <= a_ - iaz; ++iay) {\n\
      const int iax = a_ - iaz - iay;\n\
\n\
        for (int ibz = 0; ibz <= b_; ++ibz) {\n\
        for (int iby = 0; iby <= b_ - ibz; ++iby) {\n\
        const int ibx = b_ - ibz - iby;\n\
\n\
          for (int i = 0; i != rank_; ++i) {\n\
            *current_data0  += final_xa[i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_y [i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_z [i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];\n\
            *current_data1  += final_x [i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_ya[i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_z [i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];\n\
            *current_data2  += final_x [i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_y [i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_za[i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];\n\
            *current_data3  += final_xb[i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_y [i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_z [i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];\n\
            *current_data4  += final_x [i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_yb[i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_z [i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];\n\
            *current_data5  += final_x [i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_y [i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_zb[i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];\n\
            *current_data6  += final_xc[i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_y [i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_z [i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];\n\
            *current_data7  += final_x [i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_yc[i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_z [i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];\n\
            *current_data8  += final_x [i+rank_*(iax+a2*(ibx+b2*(icx+c2*idx)))] * final_y [i+rank_*(iay+a2*(iby+b2*(icy+c2*idy)))] * final_zc[i+rank_*(iaz+a2*(ibz+b2*(icz+c2*idz)))];\n\
          }\n\
          ++current_data0;\n\
          ++current_data1;\n\
          ++current_data2;\n\
          ++current_data3;\n\
          ++current_data4;\n\
          ++current_data5;\n\
          ++current_data6;\n\
          ++current_data7;\n\
          ++current_data8;\n\
\n\
        }}\n\
      }}\n\
    }}\n\
  }}\n\
\n\
}\n\
\n\
struct GVRR_Driver {\n\
\n"
#in order to minimize the parsing time, there is no header
files = ""

for a in range(0,7):
 for b in range(0,7):
  if a < b: continue
  for c in range(0,7):
   for d in range(0,7):
    if c < d: continue
    rank = int(math.ceil((a+b+c+d+2)*0.5-0.001))

    filename = "gvrr_driver_" + str(a) + "_" + str(b) + "_" + str(c) + "_" + str(d) + ".cc"
    fp = open(filename, "w")

    func_dec = "\
gvrr_driver_" + str(a) + "_" + str(b) + "_" + str(c) + "_" + str(d) + "(double* out, const double* const roots, const double* const weights, const double& coeff,\n\
    const std::array<double,3>& a, const std::array<double,3>& b, const std::array<double,3>& c, const std::array<double,3>& d,\n\
    const double* const p, const double* const q, const double& xp, const double& xq, const size_t& size_block,\n\
    const double* const expo, const double* const transx, const double* const transy, const double* const transz,\n\
    const double* const trans2x, const double* const trans2y, const double* const trans2z, double* const intermediate,\n\
    double* const final_x, double* const final_y, double* const final_z,\n\
    double* const final_xa, double* const final_xb, double* const final_xc,\n\
    double* const final_ya, double* const final_yb, double* const final_yc,\n\
    double* const final_za, double* const final_zb, double* const final_zc,\n\
    double* const workx, double* const worky, double* const workz)"

    driver += "  static void " + func_dec + ";\n" 
    contents = "\
#include <src/rysint/_gvrr_drv.h>\n\
void bagel::GVRR_Driver::" + func_dec + " {\n\
  bagel::gvrr_driver<" + str(a) + "," + str(b) + "," + str(c) + "," + str(d) + "," + str(rank) + ">(out, roots, weights, coeff, a, b, c, d, p, q, xp, xq, size_block,\n\
    expo, transx, transy, transz, trans2x, trans2y, trans2z, intermediate, final_x, final_y, final_z, final_xa, final_xb, final_xc, final_ya, final_yb, final_yc,\n\
    final_za, final_zb, final_zc, workx, worky, workz);\n\
}\n"
    fp.write(contents)
    fp.close()
    files += filename + " "
driver += "\
};\n\
\n\
}\n\
\n\
#endif\n\
"

fp = open("_gvrr_drv.h", "w")
fp.write(driver)

print files

