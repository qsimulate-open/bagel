//
// BAGEL - Parallel electron correlation program.
// Filename: _svrr_1040.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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

#include "svrrlist.h"

using namespace bagel;

// returns double array of length 40
void SVRRList::_svrr_1040(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  for (int t = 0; t != 4; ++t)
    data_[0+t] = 1.0;

  for (int t = 0; t != 4; ++t)
    data_[4+t] = C00_[t];

  for (int t = 0; t != 4; ++t)
    data_[8+t] = D00_[t];

  double cB00_current[4];
  for (int t = 0; t != 4; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 4; ++t)
    data_[12+t] = C00_[t] * data_[8+t] + cB00_current[t];

  double B01_current[4];
  for (int t = 0; t != 4; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 4; ++t)
    data_[16+t] = D00_[t] * data_[8+t] + B01_current[t];

  for (int t = 0; t != 4; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 4; ++t)
    data_[20+t] = C00_[t] * data_[16+t] + cB00_current[t] * data_[8+t];

  for (int t = 0; t != 4; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 4; ++t)
    data_[24+t] = D00_[t] * data_[16+t] + B01_current[t] * data_[8+t];

  for (int t = 0; t != 4; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 4; ++t)
    data_[28+t] = C00_[t] * data_[24+t] + cB00_current[t] * data_[16+t];

  for (int t = 0; t != 4; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 4; ++t)
    data_[32+t] = D00_[t] * data_[24+t] + B01_current[t] * data_[16+t];

  for (int t = 0; t != 4; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 4; ++t)
    data_[36+t] = C00_[t] * data_[32+t] + cB00_current[t] * data_[24+t];
}

