//
// Newint - Parallel electron correlation program.
// Filename: _vrr_0040.cc
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


#include "vrrlist.h"

// returns double array of length 15
void VRRList::_vrr_0040(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;

  data_[3] = D00_[0];
  data_[4] = D00_[1];
  data_[5] = D00_[2];

  double B01_current[3];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];

  data_[6] = D00_[0] * data_[3] + B01_current[0];
  data_[7] = D00_[1] * data_[4] + B01_current[1];
  data_[8] = D00_[2] * data_[5] + B01_current[2];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];

  data_[9] = D00_[0] * data_[6] + B01_current[0] * data_[3];
  data_[10] = D00_[1] * data_[7] + B01_current[1] * data_[4];
  data_[11] = D00_[2] * data_[8] + B01_current[2] * data_[5];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];

  data_[12] = D00_[0] * data_[9] + B01_current[0] * data_[6];
  data_[13] = D00_[1] * data_[10] + B01_current[1] * data_[7];
  data_[14] = D00_[2] * data_[11] + B01_current[2] * data_[8];
}

