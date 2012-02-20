//
// Newint - Parallel electron correlation program.
// Filename: _svrr_3010.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


#include <src/slater/svrrlist.h>

// returns double array of length 24
void SVRRList::_svrr_3010(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  data_[0] = 1.0;
  data_[1] = 1.0;
  data_[2] = 1.0;

  data_[3] = C00_[0];
  data_[4] = C00_[1];
  data_[5] = C00_[2];

  double B10_current[3];
  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];

  data_[6] = C00_[0] * data_[3] + B10_current[0];
  data_[7] = C00_[1] * data_[4] + B10_current[1];
  data_[8] = C00_[2] * data_[5] + B10_current[2];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];

  data_[9] = C00_[0] * data_[6] + B10_current[0] * data_[3];
  data_[10] = C00_[1] * data_[7] + B10_current[1] * data_[4];
  data_[11] = C00_[2] * data_[8] + B10_current[2] * data_[5];

  data_[12] = D00_[0];
  data_[13] = D00_[1];
  data_[14] = D00_[2];

  data_[15] = C00_[0] * data_[12] + B00_[0];
  data_[16] = C00_[1] * data_[13] + B00_[1];
  data_[17] = C00_[2] * data_[14] + B00_[2];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];

  data_[18] = C00_[0] * data_[15] + B10_current[0] * data_[12] + B00_[0] * data_[3];
  data_[19] = C00_[1] * data_[16] + B10_current[1] * data_[13] + B00_[1] * data_[4];
  data_[20] = C00_[2] * data_[17] + B10_current[2] * data_[14] + B00_[2] * data_[5];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];

  data_[21] = C00_[0] * data_[18] + B10_current[0] * data_[15] + B00_[0] * data_[6];
  data_[22] = C00_[1] * data_[19] + B10_current[1] * data_[16] + B00_[1] * data_[7];
  data_[23] = C00_[2] * data_[20] + B10_current[2] * data_[17] + B00_[2] * data_[8];
}

