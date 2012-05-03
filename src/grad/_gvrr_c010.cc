//
// Newint - Parallel electron correlation program.
// Filename: _gvrr_c010.cc
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

// returns double array of length 182
void GVRRList::_gvrr_c010(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[35+t] = C00_[t] * data_[28+t] + B10_current[t] * data_[21+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[42+t] = C00_[t] * data_[35+t] + B10_current[t] * data_[28+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[49+t] = C00_[t] * data_[42+t] + B10_current[t] * data_[35+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[56+t] = C00_[t] * data_[49+t] + B10_current[t] * data_[42+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[63+t] = C00_[t] * data_[56+t] + B10_current[t] * data_[49+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[70+t] = C00_[t] * data_[63+t] + B10_current[t] * data_[56+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[77+t] = C00_[t] * data_[70+t] + B10_current[t] * data_[63+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[84+t] = C00_[t] * data_[77+t] + B10_current[t] * data_[70+t];

  for (int t = 0; t != 7; ++t)
    data_[91+t] = D00_[t];

  for (int t = 0; t != 7; ++t)
    data_[98+t] = C00_[t] * data_[91+t] + B00_[t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[105+t] = C00_[t] * data_[98+t] + B10_current[t] * data_[91+t] + B00_[t] * data_[7+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[112+t] = C00_[t] * data_[105+t] + B10_current[t] * data_[98+t] + B00_[t] * data_[14+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[119+t] = C00_[t] * data_[112+t] + B10_current[t] * data_[105+t] + B00_[t] * data_[21+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[126+t] = C00_[t] * data_[119+t] + B10_current[t] * data_[112+t] + B00_[t] * data_[28+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[133+t] = C00_[t] * data_[126+t] + B10_current[t] * data_[119+t] + B00_[t] * data_[35+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[140+t] = C00_[t] * data_[133+t] + B10_current[t] * data_[126+t] + B00_[t] * data_[42+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[147+t] = C00_[t] * data_[140+t] + B10_current[t] * data_[133+t] + B00_[t] * data_[49+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[154+t] = C00_[t] * data_[147+t] + B10_current[t] * data_[140+t] + B00_[t] * data_[56+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[161+t] = C00_[t] * data_[154+t] + B10_current[t] * data_[147+t] + B00_[t] * data_[63+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[168+t] = C00_[t] * data_[161+t] + B10_current[t] * data_[154+t] + B00_[t] * data_[70+t];

  for (int t = 0; t != 7; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 7; ++t)
    data_[175+t] = C00_[t] * data_[168+t] + B10_current[t] * data_[161+t] + B00_[t] * data_[77+t];
}

