//
// Newint - Parallel electron correlation program.
// Filename: _svrr_5040.cc
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


#include "svrrlist.h"

// returns double array of length 180
void SVRRList::_svrr_5040(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[24+t] = C00_[t] * data_[18+t] + B10_current[t] * data_[12+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[30+t] = C00_[t] * data_[24+t] + B10_current[t] * data_[18+t];

  for (int t = 0; t != 6; ++t)
    data_[36+t] = D00_[t];

  double cB00_current[6];
  for (int t = 0; t != 6; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 6; ++t)
    data_[42+t] = C00_[t] * data_[36+t] + cB00_current[t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[48+t] = C00_[t] * data_[42+t] + B10_current[t] * data_[36+t] + cB00_current[t] * data_[6+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[54+t] = C00_[t] * data_[48+t] + B10_current[t] * data_[42+t] + cB00_current[t] * data_[12+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[60+t] = C00_[t] * data_[54+t] + B10_current[t] * data_[48+t] + cB00_current[t] * data_[18+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[66+t] = C00_[t] * data_[60+t] + B10_current[t] * data_[54+t] + cB00_current[t] * data_[24+t];

  double B01_current[6];
  for (int t = 0; t != 6; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 6; ++t)
    data_[72+t] = D00_[t] * data_[36+t] + B01_current[t];

  for (int t = 0; t != 6; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 6; ++t)
    data_[78+t] = C00_[t] * data_[72+t] + cB00_current[t] * data_[36+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[84+t] = C00_[t] * data_[78+t] + B10_current[t] * data_[72+t] + cB00_current[t] * data_[42+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[90+t] = C00_[t] * data_[84+t] + B10_current[t] * data_[78+t] + cB00_current[t] * data_[48+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[96+t] = C00_[t] * data_[90+t] + B10_current[t] * data_[84+t] + cB00_current[t] * data_[54+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[102+t] = C00_[t] * data_[96+t] + B10_current[t] * data_[90+t] + cB00_current[t] * data_[60+t];

  for (int t = 0; t != 6; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 6; ++t)
    data_[108+t] = D00_[t] * data_[72+t] + B01_current[t] * data_[36+t];

  for (int t = 0; t != 6; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 6; ++t)
    data_[114+t] = C00_[t] * data_[108+t] + cB00_current[t] * data_[72+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[120+t] = C00_[t] * data_[114+t] + B10_current[t] * data_[108+t] + cB00_current[t] * data_[78+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[126+t] = C00_[t] * data_[120+t] + B10_current[t] * data_[114+t] + cB00_current[t] * data_[84+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[132+t] = C00_[t] * data_[126+t] + B10_current[t] * data_[120+t] + cB00_current[t] * data_[90+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[138+t] = C00_[t] * data_[132+t] + B10_current[t] * data_[126+t] + cB00_current[t] * data_[96+t];

  for (int t = 0; t != 6; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 6; ++t)
    data_[144+t] = D00_[t] * data_[108+t] + B01_current[t] * data_[72+t];

  for (int t = 0; t != 6; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 6; ++t)
    data_[150+t] = C00_[t] * data_[144+t] + cB00_current[t] * data_[108+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[156+t] = C00_[t] * data_[150+t] + B10_current[t] * data_[144+t] + cB00_current[t] * data_[114+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[162+t] = C00_[t] * data_[156+t] + B10_current[t] * data_[150+t] + cB00_current[t] * data_[120+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[168+t] = C00_[t] * data_[162+t] + B10_current[t] * data_[156+t] + cB00_current[t] * data_[126+t];

  for (int t = 0; t != 6; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 6; ++t)
    data_[174+t] = C00_[t] * data_[168+t] + B10_current[t] * data_[162+t] + cB00_current[t] * data_[132+t];
}

