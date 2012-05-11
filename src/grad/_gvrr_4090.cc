//
// Newint - Parallel electron correlation program.
// Filename: _gvrr_4090.cc
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

// returns double array of length 350
void GVRRList::_gvrr_4090(double* data_, const double* C00, const double* D00, const double* B00, const double* B01, const double* B10) {
#ifdef __GNUC__
  double C00_[7]__attribute__((aligned(16))) = {C00[0], C00[1], C00[2], C00[3], C00[4], C00[5], C00[6]};
  double D00_[7]__attribute__((aligned(16))) = {D00[0], D00[1], D00[2], D00[3], D00[4], D00[5], D00[6]};
  double B00_[7]__attribute__((aligned(16))) = {B00[0], B00[1], B00[2], B00[3], B00[4], B00[5], B00[6]};
  double B01_[7]__attribute__((aligned(16))) = {B01[0], B01[1], B01[2], B01[3], B01[4], B01[5], B01[6]};
  double B10_[7]__attribute__((aligned(16))) = {B10[0], B10[1], B10[2], B10[3], B10[4], B10[5], B10[6]};
#else
  double* C00_ = C00;
  double* D00_ = D00;
  double* B00_ = B00;
  double* B01_ = B01;
  double* B10_ = B10;
#endif

  for (int t = 0; t != 7; ++t)
    data_[0+t] = 1.0;

  for (int t = 0; t != 7; ++t)
    data_[7+t] = C00_[t];

  double B10_current[7];
  for (int t = 0; t != 7; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[14+t] = C00_[t] * data_[7+t] + B10_current[t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[21+t] = C00_[t] * data_[14+t] + B10_current[t] * data_[7+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[28+t] = C00_[t] * data_[21+t] + B10_current[t] * data_[14+t];

  for (int t = 0; t != 7; ++t)
    data_[35+t] = D00_[t];

  double cB00_current[7];
  for (int t = 0; t != 7; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 7; ++t)
    data_[42+t] = C00_[t] * data_[35+t] + cB00_current[t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[49+t] = C00_[t] * data_[42+t] + B10_current[t] * data_[35+t] + cB00_current[t] * data_[7+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[56+t] = C00_[t] * data_[49+t] + B10_current[t] * data_[42+t] + cB00_current[t] * data_[14+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[63+t] = C00_[t] * data_[56+t] + B10_current[t] * data_[49+t] + cB00_current[t] * data_[21+t];

  double B01_current[7];
  for (int t = 0; t != 7; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 7; ++t)
    data_[70+t] = D00_[t] * data_[35+t] + B01_current[t];

  for (int t = 0; t != 7; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 7; ++t)
    data_[77+t] = C00_[t] * data_[70+t] + cB00_current[t] * data_[35+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[84+t] = C00_[t] * data_[77+t] + B10_current[t] * data_[70+t] + cB00_current[t] * data_[42+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[91+t] = C00_[t] * data_[84+t] + B10_current[t] * data_[77+t] + cB00_current[t] * data_[49+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[98+t] = C00_[t] * data_[91+t] + B10_current[t] * data_[84+t] + cB00_current[t] * data_[56+t];

  for (int t = 0; t != 7; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 7; ++t)
    data_[105+t] = D00_[t] * data_[70+t] + B01_current[t] * data_[35+t];

  for (int t = 0; t != 7; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 7; ++t)
    data_[112+t] = C00_[t] * data_[105+t] + cB00_current[t] * data_[70+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[119+t] = C00_[t] * data_[112+t] + B10_current[t] * data_[105+t] + cB00_current[t] * data_[77+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[126+t] = C00_[t] * data_[119+t] + B10_current[t] * data_[112+t] + cB00_current[t] * data_[84+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[133+t] = C00_[t] * data_[126+t] + B10_current[t] * data_[119+t] + cB00_current[t] * data_[91+t];

  for (int t = 0; t != 7; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 7; ++t)
    data_[140+t] = D00_[t] * data_[105+t] + B01_current[t] * data_[70+t];

  for (int t = 0; t != 7; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 7; ++t)
    data_[147+t] = C00_[t] * data_[140+t] + cB00_current[t] * data_[105+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[154+t] = C00_[t] * data_[147+t] + B10_current[t] * data_[140+t] + cB00_current[t] * data_[112+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[161+t] = C00_[t] * data_[154+t] + B10_current[t] * data_[147+t] + cB00_current[t] * data_[119+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[168+t] = C00_[t] * data_[161+t] + B10_current[t] * data_[154+t] + cB00_current[t] * data_[126+t];

  for (int t = 0; t != 7; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 7; ++t)
    data_[175+t] = D00_[t] * data_[140+t] + B01_current[t] * data_[105+t];

  for (int t = 0; t != 7; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 7; ++t)
    data_[182+t] = C00_[t] * data_[175+t] + cB00_current[t] * data_[140+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[189+t] = C00_[t] * data_[182+t] + B10_current[t] * data_[175+t] + cB00_current[t] * data_[147+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[196+t] = C00_[t] * data_[189+t] + B10_current[t] * data_[182+t] + cB00_current[t] * data_[154+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[203+t] = C00_[t] * data_[196+t] + B10_current[t] * data_[189+t] + cB00_current[t] * data_[161+t];

  for (int t = 0; t != 7; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 7; ++t)
    data_[210+t] = D00_[t] * data_[175+t] + B01_current[t] * data_[140+t];

  for (int t = 0; t != 7; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 7; ++t)
    data_[217+t] = C00_[t] * data_[210+t] + cB00_current[t] * data_[175+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[224+t] = C00_[t] * data_[217+t] + B10_current[t] * data_[210+t] + cB00_current[t] * data_[182+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[231+t] = C00_[t] * data_[224+t] + B10_current[t] * data_[217+t] + cB00_current[t] * data_[189+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[238+t] = C00_[t] * data_[231+t] + B10_current[t] * data_[224+t] + cB00_current[t] * data_[196+t];

  for (int t = 0; t != 7; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 7; ++t)
    data_[245+t] = D00_[t] * data_[210+t] + B01_current[t] * data_[175+t];

  for (int t = 0; t != 7; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 7; ++t)
    data_[252+t] = C00_[t] * data_[245+t] + cB00_current[t] * data_[210+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[259+t] = C00_[t] * data_[252+t] + B10_current[t] * data_[245+t] + cB00_current[t] * data_[217+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[266+t] = C00_[t] * data_[259+t] + B10_current[t] * data_[252+t] + cB00_current[t] * data_[224+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[273+t] = C00_[t] * data_[266+t] + B10_current[t] * data_[259+t] + cB00_current[t] * data_[231+t];

  for (int t = 0; t != 7; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 7; ++t)
    data_[280+t] = D00_[t] * data_[245+t] + B01_current[t] * data_[210+t];

  for (int t = 0; t != 7; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 7; ++t)
    data_[287+t] = C00_[t] * data_[280+t] + cB00_current[t] * data_[245+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[294+t] = C00_[t] * data_[287+t] + B10_current[t] * data_[280+t] + cB00_current[t] * data_[252+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[301+t] = C00_[t] * data_[294+t] + B10_current[t] * data_[287+t] + cB00_current[t] * data_[259+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[308+t] = C00_[t] * data_[301+t] + B10_current[t] * data_[294+t] + cB00_current[t] * data_[266+t];

  for (int t = 0; t != 7; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 7; ++t)
    data_[315+t] = D00_[t] * data_[280+t] + B01_current[t] * data_[245+t];

  for (int t = 0; t != 7; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 7; ++t)
    data_[322+t] = C00_[t] * data_[315+t] + cB00_current[t] * data_[280+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[329+t] = C00_[t] * data_[322+t] + B10_current[t] * data_[315+t] + cB00_current[t] * data_[287+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[336+t] = C00_[t] * data_[329+t] + B10_current[t] * data_[322+t] + cB00_current[t] * data_[294+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[343+t] = C00_[t] * data_[336+t] + B10_current[t] * data_[329+t] + cB00_current[t] * data_[301+t];
}

