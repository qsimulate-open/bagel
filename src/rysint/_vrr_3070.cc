//
// Newint - Parallel electron correlation program.
// Filename: _vrr_3070.cc
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


#include <src/rysint/vrrlist.h>

// returns double array of length 192
void VRRList::_vrr_3070(double* data_, const double* C00, const double* D00, const double* B00, const double* B01, const double* B10) {
#ifdef __GNUC__
  const double C00_[6]__attribute__((aligned(32))) = {C00[0], C00[1], C00[2], C00[3], C00[4], C00[5]};
  const double D00_[6]__attribute__((aligned(32))) = {D00[0], D00[1], D00[2], D00[3], D00[4], D00[5]};
  const double B00_[6]__attribute__((aligned(32))) = {B00[0], B00[1], B00[2], B00[3], B00[4], B00[5]};
  const double B01_[6]__attribute__((aligned(32))) = {B01[0], B01[1], B01[2], B01[3], B01[4], B01[5]};
  const double B10_[6]__attribute__((aligned(32))) = {B10[0], B10[1], B10[2], B10[3], B10[4], B10[5]};
#else
  const double* C00_ = C00;
  const double* D00_ = D00;
  const double* B00_ = B00;
  const double* B01_ = B01;
  const double* B10_ = B10;
#endif

  for (int t = 0; t != 6; ++t)
    data_[0+t] = 1.0;

  for (int t = 0; t != 6; ++t)
    data_[6+t] = C00_[t];

  double B10_current[6];
  for (int t = 0; t != 6; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[12+t] = C00_[t] * data_[6+t] + B10_current[t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[18+t] = C00_[t] * data_[12+t] + B10_current[t] * data_[6+t];

  for (int t = 0; t != 6; ++t)
    data_[24+t] = D00_[t];

  double cB00_current[6];
  for (int t = 0; t != 6; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 6; ++t)
    data_[30+t] = C00_[t] * data_[24+t] + cB00_current[t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[36+t] = C00_[t] * data_[30+t] + B10_current[t] * data_[24+t] + cB00_current[t] * data_[6+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[42+t] = C00_[t] * data_[36+t] + B10_current[t] * data_[30+t] + cB00_current[t] * data_[12+t];

  double B01_current[6];
  for (int t = 0; t != 6; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 6; ++t)
    data_[48+t] = D00_[t] * data_[24+t] + B01_current[t];

  for (int t = 0; t != 6; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 6; ++t)
    data_[54+t] = C00_[t] * data_[48+t] + cB00_current[t] * data_[24+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[60+t] = C00_[t] * data_[54+t] + B10_current[t] * data_[48+t] + cB00_current[t] * data_[30+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[66+t] = C00_[t] * data_[60+t] + B10_current[t] * data_[54+t] + cB00_current[t] * data_[36+t];

  for (int t = 0; t != 6; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 6; ++t)
    data_[72+t] = D00_[t] * data_[48+t] + B01_current[t] * data_[24+t];

  for (int t = 0; t != 6; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 6; ++t)
    data_[78+t] = C00_[t] * data_[72+t] + cB00_current[t] * data_[48+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[84+t] = C00_[t] * data_[78+t] + B10_current[t] * data_[72+t] + cB00_current[t] * data_[54+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[90+t] = C00_[t] * data_[84+t] + B10_current[t] * data_[78+t] + cB00_current[t] * data_[60+t];

  for (int t = 0; t != 6; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 6; ++t)
    data_[96+t] = D00_[t] * data_[72+t] + B01_current[t] * data_[48+t];

  for (int t = 0; t != 6; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 6; ++t)
    data_[102+t] = C00_[t] * data_[96+t] + cB00_current[t] * data_[72+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[108+t] = C00_[t] * data_[102+t] + B10_current[t] * data_[96+t] + cB00_current[t] * data_[78+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[114+t] = C00_[t] * data_[108+t] + B10_current[t] * data_[102+t] + cB00_current[t] * data_[84+t];

  for (int t = 0; t != 6; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 6; ++t)
    data_[120+t] = D00_[t] * data_[96+t] + B01_current[t] * data_[72+t];

  for (int t = 0; t != 6; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 6; ++t)
    data_[126+t] = C00_[t] * data_[120+t] + cB00_current[t] * data_[96+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[132+t] = C00_[t] * data_[126+t] + B10_current[t] * data_[120+t] + cB00_current[t] * data_[102+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[138+t] = C00_[t] * data_[132+t] + B10_current[t] * data_[126+t] + cB00_current[t] * data_[108+t];

  for (int t = 0; t != 6; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 6; ++t)
    data_[144+t] = D00_[t] * data_[120+t] + B01_current[t] * data_[96+t];

  for (int t = 0; t != 6; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 6; ++t)
    data_[150+t] = C00_[t] * data_[144+t] + cB00_current[t] * data_[120+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[156+t] = C00_[t] * data_[150+t] + B10_current[t] * data_[144+t] + cB00_current[t] * data_[126+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[162+t] = C00_[t] * data_[156+t] + B10_current[t] * data_[150+t] + cB00_current[t] * data_[132+t];

  for (int t = 0; t != 6; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 6; ++t)
    data_[168+t] = D00_[t] * data_[144+t] + B01_current[t] * data_[120+t];

  for (int t = 0; t != 6; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 6; ++t)
    data_[174+t] = C00_[t] * data_[168+t] + cB00_current[t] * data_[144+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[180+t] = C00_[t] * data_[174+t] + B10_current[t] * data_[168+t] + cB00_current[t] * data_[150+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[186+t] = C00_[t] * data_[180+t] + B10_current[t] * data_[174+t] + cB00_current[t] * data_[156+t];
}

