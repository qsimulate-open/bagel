//
// Newint - Parallel electron correlation program.
// Filename: _gvrr_b0a0.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/grad/gvrrlist.h>

// returns double array of length 1452
void GVRRList::_gvrr_b0a0(double* data_, const double* C00, const double* D00, const double* B00, const double* B01, const double* B10) {
#ifdef __GNUC__
  const double C00_[11]__attribute__((aligned(16))) = {C00[0], C00[1], C00[2], C00[3], C00[4], C00[5], C00[6], C00[7], C00[8], C00[9], C00[10]};
  const double D00_[11]__attribute__((aligned(16))) = {D00[0], D00[1], D00[2], D00[3], D00[4], D00[5], D00[6], D00[7], D00[8], D00[9], D00[10]};
  const double B00_[11]__attribute__((aligned(16))) = {B00[0], B00[1], B00[2], B00[3], B00[4], B00[5], B00[6], B00[7], B00[8], B00[9], B00[10]};
  const double B01_[11]__attribute__((aligned(16))) = {B01[0], B01[1], B01[2], B01[3], B01[4], B01[5], B01[6], B01[7], B01[8], B01[9], B01[10]};
  const double B10_[11]__attribute__((aligned(16))) = {B10[0], B10[1], B10[2], B10[3], B10[4], B10[5], B10[6], B10[7], B10[8], B10[9], B10[10]};
#else
  const double* C00_ = C00;
  const double* D00_ = D00;
  const double* B00_ = B00;
  const double* B01_ = B01;
  const double* B10_ = B10;
#endif

  for (int t = 0; t != 11; ++t)
    data_[0+t] = 1.0;

  for (int t = 0; t != 11; ++t)
    data_[11+t] = C00_[t];

  double B10_current[11];
  for (int t = 0; t != 11; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[22+t] = C00_[t] * data_[11+t] + B10_current[t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[33+t] = C00_[t] * data_[22+t] + B10_current[t] * data_[11+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[44+t] = C00_[t] * data_[33+t] + B10_current[t] * data_[22+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[55+t] = C00_[t] * data_[44+t] + B10_current[t] * data_[33+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[66+t] = C00_[t] * data_[55+t] + B10_current[t] * data_[44+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[77+t] = C00_[t] * data_[66+t] + B10_current[t] * data_[55+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[88+t] = C00_[t] * data_[77+t] + B10_current[t] * data_[66+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[99+t] = C00_[t] * data_[88+t] + B10_current[t] * data_[77+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[110+t] = C00_[t] * data_[99+t] + B10_current[t] * data_[88+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[121+t] = C00_[t] * data_[110+t] + B10_current[t] * data_[99+t];

  for (int t = 0; t != 11; ++t)
    data_[132+t] = D00_[t];

  double cB00_current[11];
  for (int t = 0; t != 11; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 11; ++t)
    data_[143+t] = C00_[t] * data_[132+t] + cB00_current[t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[154+t] = C00_[t] * data_[143+t] + B10_current[t] * data_[132+t] + cB00_current[t] * data_[11+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[165+t] = C00_[t] * data_[154+t] + B10_current[t] * data_[143+t] + cB00_current[t] * data_[22+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[176+t] = C00_[t] * data_[165+t] + B10_current[t] * data_[154+t] + cB00_current[t] * data_[33+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[187+t] = C00_[t] * data_[176+t] + B10_current[t] * data_[165+t] + cB00_current[t] * data_[44+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[198+t] = C00_[t] * data_[187+t] + B10_current[t] * data_[176+t] + cB00_current[t] * data_[55+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[209+t] = C00_[t] * data_[198+t] + B10_current[t] * data_[187+t] + cB00_current[t] * data_[66+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[220+t] = C00_[t] * data_[209+t] + B10_current[t] * data_[198+t] + cB00_current[t] * data_[77+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[231+t] = C00_[t] * data_[220+t] + B10_current[t] * data_[209+t] + cB00_current[t] * data_[88+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[242+t] = C00_[t] * data_[231+t] + B10_current[t] * data_[220+t] + cB00_current[t] * data_[99+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[253+t] = C00_[t] * data_[242+t] + B10_current[t] * data_[231+t] + cB00_current[t] * data_[110+t];

  double B01_current[11];
  for (int t = 0; t != 11; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 11; ++t)
    data_[264+t] = D00_[t] * data_[132+t] + B01_current[t];

  for (int t = 0; t != 11; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 11; ++t)
    data_[275+t] = C00_[t] * data_[264+t] + cB00_current[t] * data_[132+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[286+t] = C00_[t] * data_[275+t] + B10_current[t] * data_[264+t] + cB00_current[t] * data_[143+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[297+t] = C00_[t] * data_[286+t] + B10_current[t] * data_[275+t] + cB00_current[t] * data_[154+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[308+t] = C00_[t] * data_[297+t] + B10_current[t] * data_[286+t] + cB00_current[t] * data_[165+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[319+t] = C00_[t] * data_[308+t] + B10_current[t] * data_[297+t] + cB00_current[t] * data_[176+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[330+t] = C00_[t] * data_[319+t] + B10_current[t] * data_[308+t] + cB00_current[t] * data_[187+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[341+t] = C00_[t] * data_[330+t] + B10_current[t] * data_[319+t] + cB00_current[t] * data_[198+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[352+t] = C00_[t] * data_[341+t] + B10_current[t] * data_[330+t] + cB00_current[t] * data_[209+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[363+t] = C00_[t] * data_[352+t] + B10_current[t] * data_[341+t] + cB00_current[t] * data_[220+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[374+t] = C00_[t] * data_[363+t] + B10_current[t] * data_[352+t] + cB00_current[t] * data_[231+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[385+t] = C00_[t] * data_[374+t] + B10_current[t] * data_[363+t] + cB00_current[t] * data_[242+t];

  for (int t = 0; t != 11; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 11; ++t)
    data_[396+t] = D00_[t] * data_[264+t] + B01_current[t] * data_[132+t];

  for (int t = 0; t != 11; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 11; ++t)
    data_[407+t] = C00_[t] * data_[396+t] + cB00_current[t] * data_[264+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[418+t] = C00_[t] * data_[407+t] + B10_current[t] * data_[396+t] + cB00_current[t] * data_[275+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[429+t] = C00_[t] * data_[418+t] + B10_current[t] * data_[407+t] + cB00_current[t] * data_[286+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[440+t] = C00_[t] * data_[429+t] + B10_current[t] * data_[418+t] + cB00_current[t] * data_[297+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[451+t] = C00_[t] * data_[440+t] + B10_current[t] * data_[429+t] + cB00_current[t] * data_[308+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[462+t] = C00_[t] * data_[451+t] + B10_current[t] * data_[440+t] + cB00_current[t] * data_[319+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[473+t] = C00_[t] * data_[462+t] + B10_current[t] * data_[451+t] + cB00_current[t] * data_[330+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[484+t] = C00_[t] * data_[473+t] + B10_current[t] * data_[462+t] + cB00_current[t] * data_[341+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[495+t] = C00_[t] * data_[484+t] + B10_current[t] * data_[473+t] + cB00_current[t] * data_[352+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[506+t] = C00_[t] * data_[495+t] + B10_current[t] * data_[484+t] + cB00_current[t] * data_[363+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[517+t] = C00_[t] * data_[506+t] + B10_current[t] * data_[495+t] + cB00_current[t] * data_[374+t];

  for (int t = 0; t != 11; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 11; ++t)
    data_[528+t] = D00_[t] * data_[396+t] + B01_current[t] * data_[264+t];

  for (int t = 0; t != 11; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 11; ++t)
    data_[539+t] = C00_[t] * data_[528+t] + cB00_current[t] * data_[396+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[550+t] = C00_[t] * data_[539+t] + B10_current[t] * data_[528+t] + cB00_current[t] * data_[407+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[561+t] = C00_[t] * data_[550+t] + B10_current[t] * data_[539+t] + cB00_current[t] * data_[418+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[572+t] = C00_[t] * data_[561+t] + B10_current[t] * data_[550+t] + cB00_current[t] * data_[429+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[583+t] = C00_[t] * data_[572+t] + B10_current[t] * data_[561+t] + cB00_current[t] * data_[440+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[594+t] = C00_[t] * data_[583+t] + B10_current[t] * data_[572+t] + cB00_current[t] * data_[451+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[605+t] = C00_[t] * data_[594+t] + B10_current[t] * data_[583+t] + cB00_current[t] * data_[462+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[616+t] = C00_[t] * data_[605+t] + B10_current[t] * data_[594+t] + cB00_current[t] * data_[473+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[627+t] = C00_[t] * data_[616+t] + B10_current[t] * data_[605+t] + cB00_current[t] * data_[484+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[638+t] = C00_[t] * data_[627+t] + B10_current[t] * data_[616+t] + cB00_current[t] * data_[495+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[649+t] = C00_[t] * data_[638+t] + B10_current[t] * data_[627+t] + cB00_current[t] * data_[506+t];

  for (int t = 0; t != 11; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 11; ++t)
    data_[660+t] = D00_[t] * data_[528+t] + B01_current[t] * data_[396+t];

  for (int t = 0; t != 11; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 11; ++t)
    data_[671+t] = C00_[t] * data_[660+t] + cB00_current[t] * data_[528+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[682+t] = C00_[t] * data_[671+t] + B10_current[t] * data_[660+t] + cB00_current[t] * data_[539+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[693+t] = C00_[t] * data_[682+t] + B10_current[t] * data_[671+t] + cB00_current[t] * data_[550+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[704+t] = C00_[t] * data_[693+t] + B10_current[t] * data_[682+t] + cB00_current[t] * data_[561+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[715+t] = C00_[t] * data_[704+t] + B10_current[t] * data_[693+t] + cB00_current[t] * data_[572+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[726+t] = C00_[t] * data_[715+t] + B10_current[t] * data_[704+t] + cB00_current[t] * data_[583+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[737+t] = C00_[t] * data_[726+t] + B10_current[t] * data_[715+t] + cB00_current[t] * data_[594+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[748+t] = C00_[t] * data_[737+t] + B10_current[t] * data_[726+t] + cB00_current[t] * data_[605+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[759+t] = C00_[t] * data_[748+t] + B10_current[t] * data_[737+t] + cB00_current[t] * data_[616+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[770+t] = C00_[t] * data_[759+t] + B10_current[t] * data_[748+t] + cB00_current[t] * data_[627+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[781+t] = C00_[t] * data_[770+t] + B10_current[t] * data_[759+t] + cB00_current[t] * data_[638+t];

  for (int t = 0; t != 11; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 11; ++t)
    data_[792+t] = D00_[t] * data_[660+t] + B01_current[t] * data_[528+t];

  for (int t = 0; t != 11; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 11; ++t)
    data_[803+t] = C00_[t] * data_[792+t] + cB00_current[t] * data_[660+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[814+t] = C00_[t] * data_[803+t] + B10_current[t] * data_[792+t] + cB00_current[t] * data_[671+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[825+t] = C00_[t] * data_[814+t] + B10_current[t] * data_[803+t] + cB00_current[t] * data_[682+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[836+t] = C00_[t] * data_[825+t] + B10_current[t] * data_[814+t] + cB00_current[t] * data_[693+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[847+t] = C00_[t] * data_[836+t] + B10_current[t] * data_[825+t] + cB00_current[t] * data_[704+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[858+t] = C00_[t] * data_[847+t] + B10_current[t] * data_[836+t] + cB00_current[t] * data_[715+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[869+t] = C00_[t] * data_[858+t] + B10_current[t] * data_[847+t] + cB00_current[t] * data_[726+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[880+t] = C00_[t] * data_[869+t] + B10_current[t] * data_[858+t] + cB00_current[t] * data_[737+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[891+t] = C00_[t] * data_[880+t] + B10_current[t] * data_[869+t] + cB00_current[t] * data_[748+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[902+t] = C00_[t] * data_[891+t] + B10_current[t] * data_[880+t] + cB00_current[t] * data_[759+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[913+t] = C00_[t] * data_[902+t] + B10_current[t] * data_[891+t] + cB00_current[t] * data_[770+t];

  for (int t = 0; t != 11; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 11; ++t)
    data_[924+t] = D00_[t] * data_[792+t] + B01_current[t] * data_[660+t];

  for (int t = 0; t != 11; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 11; ++t)
    data_[935+t] = C00_[t] * data_[924+t] + cB00_current[t] * data_[792+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[946+t] = C00_[t] * data_[935+t] + B10_current[t] * data_[924+t] + cB00_current[t] * data_[803+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[957+t] = C00_[t] * data_[946+t] + B10_current[t] * data_[935+t] + cB00_current[t] * data_[814+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[968+t] = C00_[t] * data_[957+t] + B10_current[t] * data_[946+t] + cB00_current[t] * data_[825+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[979+t] = C00_[t] * data_[968+t] + B10_current[t] * data_[957+t] + cB00_current[t] * data_[836+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[990+t] = C00_[t] * data_[979+t] + B10_current[t] * data_[968+t] + cB00_current[t] * data_[847+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1001+t] = C00_[t] * data_[990+t] + B10_current[t] * data_[979+t] + cB00_current[t] * data_[858+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1012+t] = C00_[t] * data_[1001+t] + B10_current[t] * data_[990+t] + cB00_current[t] * data_[869+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1023+t] = C00_[t] * data_[1012+t] + B10_current[t] * data_[1001+t] + cB00_current[t] * data_[880+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1034+t] = C00_[t] * data_[1023+t] + B10_current[t] * data_[1012+t] + cB00_current[t] * data_[891+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1045+t] = C00_[t] * data_[1034+t] + B10_current[t] * data_[1023+t] + cB00_current[t] * data_[902+t];

  for (int t = 0; t != 11; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 11; ++t)
    data_[1056+t] = D00_[t] * data_[924+t] + B01_current[t] * data_[792+t];

  for (int t = 0; t != 11; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 11; ++t)
    data_[1067+t] = C00_[t] * data_[1056+t] + cB00_current[t] * data_[924+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1078+t] = C00_[t] * data_[1067+t] + B10_current[t] * data_[1056+t] + cB00_current[t] * data_[935+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1089+t] = C00_[t] * data_[1078+t] + B10_current[t] * data_[1067+t] + cB00_current[t] * data_[946+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1100+t] = C00_[t] * data_[1089+t] + B10_current[t] * data_[1078+t] + cB00_current[t] * data_[957+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1111+t] = C00_[t] * data_[1100+t] + B10_current[t] * data_[1089+t] + cB00_current[t] * data_[968+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1122+t] = C00_[t] * data_[1111+t] + B10_current[t] * data_[1100+t] + cB00_current[t] * data_[979+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1133+t] = C00_[t] * data_[1122+t] + B10_current[t] * data_[1111+t] + cB00_current[t] * data_[990+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1144+t] = C00_[t] * data_[1133+t] + B10_current[t] * data_[1122+t] + cB00_current[t] * data_[1001+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1155+t] = C00_[t] * data_[1144+t] + B10_current[t] * data_[1133+t] + cB00_current[t] * data_[1012+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1166+t] = C00_[t] * data_[1155+t] + B10_current[t] * data_[1144+t] + cB00_current[t] * data_[1023+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1177+t] = C00_[t] * data_[1166+t] + B10_current[t] * data_[1155+t] + cB00_current[t] * data_[1034+t];

  for (int t = 0; t != 11; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 11; ++t)
    data_[1188+t] = D00_[t] * data_[1056+t] + B01_current[t] * data_[924+t];

  for (int t = 0; t != 11; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 11; ++t)
    data_[1199+t] = C00_[t] * data_[1188+t] + cB00_current[t] * data_[1056+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1210+t] = C00_[t] * data_[1199+t] + B10_current[t] * data_[1188+t] + cB00_current[t] * data_[1067+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1221+t] = C00_[t] * data_[1210+t] + B10_current[t] * data_[1199+t] + cB00_current[t] * data_[1078+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1232+t] = C00_[t] * data_[1221+t] + B10_current[t] * data_[1210+t] + cB00_current[t] * data_[1089+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1243+t] = C00_[t] * data_[1232+t] + B10_current[t] * data_[1221+t] + cB00_current[t] * data_[1100+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1254+t] = C00_[t] * data_[1243+t] + B10_current[t] * data_[1232+t] + cB00_current[t] * data_[1111+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1265+t] = C00_[t] * data_[1254+t] + B10_current[t] * data_[1243+t] + cB00_current[t] * data_[1122+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1276+t] = C00_[t] * data_[1265+t] + B10_current[t] * data_[1254+t] + cB00_current[t] * data_[1133+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1287+t] = C00_[t] * data_[1276+t] + B10_current[t] * data_[1265+t] + cB00_current[t] * data_[1144+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1298+t] = C00_[t] * data_[1287+t] + B10_current[t] * data_[1276+t] + cB00_current[t] * data_[1155+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1309+t] = C00_[t] * data_[1298+t] + B10_current[t] * data_[1287+t] + cB00_current[t] * data_[1166+t];

  for (int t = 0; t != 11; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 11; ++t)
    data_[1320+t] = D00_[t] * data_[1188+t] + B01_current[t] * data_[1056+t];

  for (int t = 0; t != 11; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 11; ++t)
    data_[1331+t] = C00_[t] * data_[1320+t] + cB00_current[t] * data_[1188+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1342+t] = C00_[t] * data_[1331+t] + B10_current[t] * data_[1320+t] + cB00_current[t] * data_[1199+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1353+t] = C00_[t] * data_[1342+t] + B10_current[t] * data_[1331+t] + cB00_current[t] * data_[1210+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1364+t] = C00_[t] * data_[1353+t] + B10_current[t] * data_[1342+t] + cB00_current[t] * data_[1221+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1375+t] = C00_[t] * data_[1364+t] + B10_current[t] * data_[1353+t] + cB00_current[t] * data_[1232+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1386+t] = C00_[t] * data_[1375+t] + B10_current[t] * data_[1364+t] + cB00_current[t] * data_[1243+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1397+t] = C00_[t] * data_[1386+t] + B10_current[t] * data_[1375+t] + cB00_current[t] * data_[1254+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1408+t] = C00_[t] * data_[1397+t] + B10_current[t] * data_[1386+t] + cB00_current[t] * data_[1265+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1419+t] = C00_[t] * data_[1408+t] + B10_current[t] * data_[1397+t] + cB00_current[t] * data_[1276+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1430+t] = C00_[t] * data_[1419+t] + B10_current[t] * data_[1408+t] + cB00_current[t] * data_[1287+t];

  for (int t = 0; t != 11; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 11; ++t)
    data_[1441+t] = C00_[t] * data_[1430+t] + B10_current[t] * data_[1419+t] + cB00_current[t] * data_[1298+t];
}

