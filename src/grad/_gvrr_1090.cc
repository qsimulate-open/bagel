//
// Newint - Parallel electron correlation program.
// Filename: _gvrr_1090.cc
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

// returns double array of length 100
void GVRRList::_gvrr_1090(double* data_, const double* C00, const double* D00, const double* B00, const double* B01, const double* B10) {
#ifdef __GNUC__
  double C00_[5]__attribute__((aligned(16))) = {C00[0], C00[1], C00[2], C00[3], C00[4]};
  double D00_[5]__attribute__((aligned(16))) = {D00[0], D00[1], D00[2], D00[3], D00[4]};
  double B00_[5]__attribute__((aligned(16))) = {B00[0], B00[1], B00[2], B00[3], B00[4]};
  double B01_[5]__attribute__((aligned(16))) = {B01[0], B01[1], B01[2], B01[3], B01[4]};
  double B10_[5]__attribute__((aligned(16))) = {B10[0], B10[1], B10[2], B10[3], B10[4]};
#else
  double* C00_ = C00;
  double* D00_ = D00;
  double* B00_ = B00;
  double* B01_ = B01;
  double* B10_ = B10;
#endif

  for (int t = 0; t != 5; ++t)
    data_[0+t] = 1.0;

  for (int t = 0; t != 5; ++t)
    data_[5+t] = C00_[t];

  for (int t = 0; t != 5; ++t)
    data_[10+t] = D00_[t];

  double cB00_current[5];
  for (int t = 0; t != 5; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 5; ++t)
    data_[15+t] = C00_[t] * data_[10+t] + cB00_current[t];

  double B01_current[5];
  for (int t = 0; t != 5; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 5; ++t)
    data_[20+t] = D00_[t] * data_[10+t] + B01_current[t];

  for (int t = 0; t != 5; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 5; ++t)
    data_[25+t] = C00_[t] * data_[20+t] + cB00_current[t] * data_[10+t];

  for (int t = 0; t != 5; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 5; ++t)
    data_[30+t] = D00_[t] * data_[20+t] + B01_current[t] * data_[10+t];

  for (int t = 0; t != 5; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 5; ++t)
    data_[35+t] = C00_[t] * data_[30+t] + cB00_current[t] * data_[20+t];

  for (int t = 0; t != 5; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 5; ++t)
    data_[40+t] = D00_[t] * data_[30+t] + B01_current[t] * data_[20+t];

  for (int t = 0; t != 5; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 5; ++t)
    data_[45+t] = C00_[t] * data_[40+t] + cB00_current[t] * data_[30+t];

  for (int t = 0; t != 5; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 5; ++t)
    data_[50+t] = D00_[t] * data_[40+t] + B01_current[t] * data_[30+t];

  for (int t = 0; t != 5; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 5; ++t)
    data_[55+t] = C00_[t] * data_[50+t] + cB00_current[t] * data_[40+t];

  for (int t = 0; t != 5; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 5; ++t)
    data_[60+t] = D00_[t] * data_[50+t] + B01_current[t] * data_[40+t];

  for (int t = 0; t != 5; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 5; ++t)
    data_[65+t] = C00_[t] * data_[60+t] + cB00_current[t] * data_[50+t];

  for (int t = 0; t != 5; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 5; ++t)
    data_[70+t] = D00_[t] * data_[60+t] + B01_current[t] * data_[50+t];

  for (int t = 0; t != 5; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 5; ++t)
    data_[75+t] = C00_[t] * data_[70+t] + cB00_current[t] * data_[60+t];

  for (int t = 0; t != 5; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 5; ++t)
    data_[80+t] = D00_[t] * data_[70+t] + B01_current[t] * data_[60+t];

  for (int t = 0; t != 5; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 5; ++t)
    data_[85+t] = C00_[t] * data_[80+t] + cB00_current[t] * data_[70+t];

  for (int t = 0; t != 5; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 5; ++t)
    data_[90+t] = D00_[t] * data_[80+t] + B01_current[t] * data_[70+t];

  for (int t = 0; t != 5; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 5; ++t)
    data_[95+t] = C00_[t] * data_[90+t] + cB00_current[t] * data_[80+t];
}

