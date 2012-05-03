//
// Newint - Parallel electron correlation program.
// Filename: _gvrr_3020.cc
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


#include "gvrrlist.h"

// returns double array of length 36
void GVRRList::_gvrr_3020(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  for (int t = 0; t != 3; ++t)
    data_[0+t] = 1.0;

  for (int t = 0; t != 3; ++t)
    data_[3+t] = C00_[t];

  double B10_current[3];
  for (int t = 0; t != 3; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 3; ++t)
    data_[6+t] = C00_[t] * data_[3+t] + B10_current[t];

  for (int t = 0; t != 3; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 3; ++t)
    data_[9+t] = C00_[t] * data_[6+t] + B10_current[t] * data_[3+t];

  for (int t = 0; t != 3; ++t)
    data_[12+t] = D00_[t];

  double cB00_current[3];
  for (int t = 0; t != 3; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 3; ++t)
    data_[15+t] = C00_[t] * data_[12+t] + cB00_current[t];

  for (int t = 0; t != 3; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 3; ++t)
    data_[18+t] = C00_[t] * data_[15+t] + B10_current[t] * data_[12+t] + cB00_current[t] * data_[3+t];

  for (int t = 0; t != 3; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 3; ++t)
    data_[21+t] = C00_[t] * data_[18+t] + B10_current[t] * data_[15+t] + cB00_current[t] * data_[6+t];


  for (int t = 0; t != 3; ++t)
    data_[24+t] = D00_[t] * data_[12+t] + B01_[t];

  for (int t = 0; t != 3; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 3; ++t)
    data_[27+t] = C00_[t] * data_[24+t] + cB00_current[t] * data_[12+t];

  for (int t = 0; t != 3; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 3; ++t)
    data_[30+t] = C00_[t] * data_[27+t] + B10_current[t] * data_[24+t] + cB00_current[t] * data_[15+t];

  for (int t = 0; t != 3; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 3; ++t)
    data_[33+t] = C00_[t] * data_[30+t] + B10_current[t] * data_[27+t] + cB00_current[t] * data_[18+t];
}

