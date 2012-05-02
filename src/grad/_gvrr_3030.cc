//
// Newint - Parallel electron correlation program.
// Filename: _gvrr_3030.cc
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

// returns double array of length 48
void GVRRList::_gvrr_3030(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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

  double cB00_current[3];
  cB00_current[0] = B00_[0];
  cB00_current[1] = B00_[1];
  cB00_current[2] = B00_[2];

  data_[15] = C00_[0] * data_[12] + cB00_current[0];
  data_[16] = C00_[1] * data_[13] + cB00_current[1];
  data_[17] = C00_[2] * data_[14] + cB00_current[2];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];

  data_[18] = C00_[0] * data_[15] + B10_current[0] * data_[12] + cB00_current[0] * data_[3];
  data_[19] = C00_[1] * data_[16] + B10_current[1] * data_[13] + cB00_current[1] * data_[4];
  data_[20] = C00_[2] * data_[17] + B10_current[2] * data_[14] + cB00_current[2] * data_[5];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];

  data_[21] = C00_[0] * data_[18] + B10_current[0] * data_[15] + cB00_current[0] * data_[6];
  data_[22] = C00_[1] * data_[19] + B10_current[1] * data_[16] + cB00_current[1] * data_[7];
  data_[23] = C00_[2] * data_[20] + B10_current[2] * data_[17] + cB00_current[2] * data_[8];

  double B01_current[3];
  B01_current[0] = B01_[0];
  B01_current[1] = B01_[1];
  B01_current[2] = B01_[2];

  data_[24] = D00_[0] * data_[12] + B01_current[0];
  data_[25] = D00_[1] * data_[13] + B01_current[1];
  data_[26] = D00_[2] * data_[14] + B01_current[2];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];

  data_[27] = C00_[0] * data_[24] + cB00_current[0] * data_[12];
  data_[28] = C00_[1] * data_[25] + cB00_current[1] * data_[13];
  data_[29] = C00_[2] * data_[26] + cB00_current[2] * data_[14];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];

  data_[30] = C00_[0] * data_[27] + B10_current[0] * data_[24] + cB00_current[0] * data_[15];
  data_[31] = C00_[1] * data_[28] + B10_current[1] * data_[25] + cB00_current[1] * data_[16];
  data_[32] = C00_[2] * data_[29] + B10_current[2] * data_[26] + cB00_current[2] * data_[17];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];

  data_[33] = C00_[0] * data_[30] + B10_current[0] * data_[27] + cB00_current[0] * data_[18];
  data_[34] = C00_[1] * data_[31] + B10_current[1] * data_[28] + cB00_current[1] * data_[19];
  data_[35] = C00_[2] * data_[32] + B10_current[2] * data_[29] + cB00_current[2] * data_[20];

  B01_current[0] += B01_[0];
  B01_current[1] += B01_[1];
  B01_current[2] += B01_[2];

  data_[36] = D00_[0] * data_[24] + B01_current[0] * data_[12];
  data_[37] = D00_[1] * data_[25] + B01_current[1] * data_[13];
  data_[38] = D00_[2] * data_[26] + B01_current[2] * data_[14];

  cB00_current[0] += B00_[0];
  cB00_current[1] += B00_[1];
  cB00_current[2] += B00_[2];

  data_[39] = C00_[0] * data_[36] + cB00_current[0] * data_[24];
  data_[40] = C00_[1] * data_[37] + cB00_current[1] * data_[25];
  data_[41] = C00_[2] * data_[38] + cB00_current[2] * data_[26];

  B10_current[0] = B10_[0];
  B10_current[1] = B10_[1];
  B10_current[2] = B10_[2];

  data_[42] = C00_[0] * data_[39] + B10_current[0] * data_[36] + cB00_current[0] * data_[27];
  data_[43] = C00_[1] * data_[40] + B10_current[1] * data_[37] + cB00_current[1] * data_[28];
  data_[44] = C00_[2] * data_[41] + B10_current[2] * data_[38] + cB00_current[2] * data_[29];

  B10_current[0] += B10_[0];
  B10_current[1] += B10_[1];
  B10_current[2] += B10_[2];

  data_[45] = C00_[0] * data_[42] + B10_current[0] * data_[39] + cB00_current[0] * data_[30];
  data_[46] = C00_[1] * data_[43] + B10_current[1] * data_[40] + cB00_current[1] * data_[31];
  data_[47] = C00_[2] * data_[44] + B10_current[2] * data_[41] + cB00_current[2] * data_[32];
}

