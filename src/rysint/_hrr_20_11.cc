//
// Newint - Parallel electron correlation program.
// Filename: _hrr_20_11.cc
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


#include "hrrlist.h"
#include <algorithm>

using namespace std;

// modified by hand

void HRRList::perform_HRR_20_11(const int nloop, const double* data_start, const double* AB, double* data_out) {
  const double* current_data = data_start;
  double* current_out  = data_out;
  for (int c = 0; c != nloop; ++c, current_data += 9, current_out += 9) {
    current_out[0] = current_data[3] + AB[0] * current_data[0]; // a0_x
    current_out[1] = current_data[4] + AB[1] * current_data[0]; // a0_y
    current_out[2] = current_data[6] + AB[2] * current_data[0]; // a0_z
    current_out[3] = current_data[4] + AB[0] * current_data[1]; // a0_x
    current_out[4] = current_data[5] + AB[1] * current_data[1]; // a0_y
    current_out[5] = current_data[7] + AB[2] * current_data[1]; // a0_z
    current_out[6] = current_data[6] + AB[0] * current_data[2]; // a0_x
    current_out[7] = current_data[7] + AB[1] * current_data[2]; // a0_y
    current_out[8] = current_data[8] + AB[2] * current_data[2]; // a0_z
  }
}

