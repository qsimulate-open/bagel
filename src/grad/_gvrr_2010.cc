//
// BAGEL - Parallel electron correlation program.
// Filename: _gvrr_2010.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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


#include <src/grad/gvrrlist.h>

// returns double array of length 12
void GVRRList::_gvrr_2010(double* data_, const double* C00, const double* D00, const double* B00, const double* B01, const double* B10) {
#ifdef __GNUC__
  const double C00_[2]__attribute__((aligned(32))) = {C00[0], C00[1]};
  const double D00_[2]__attribute__((aligned(32))) = {D00[0], D00[1]};
  const double B00_[2]__attribute__((aligned(32))) = {B00[0], B00[1]};
  const double B01_[2]__attribute__((aligned(32))) = {B01[0], B01[1]};
  const double B10_[2]__attribute__((aligned(32))) = {B10[0], B10[1]};
#else
  const double* C00_ = C00;
  const double* D00_ = D00;
  const double* B00_ = B00;
  const double* B01_ = B01;
  const double* B10_ = B10;
#endif

  for (int t = 0; t != 2; ++t)
    data_[0+t] = 1.0;

  for (int t = 0; t != 2; ++t)
    data_[2+t] = C00_[t];

  for (int t = 0; t != 2; ++t)
    data_[4+t] = C00_[t] * data_[2+t] + B10_[t];

  for (int t = 0; t != 2; ++t)
    data_[6+t] = D00_[t];

  for (int t = 0; t != 2; ++t)
    data_[8+t] = C00_[t] * data_[6+t] + B00_[t];

  for (int t = 0; t != 2; ++t)
    data_[10+t] = C00_[t] * data_[8+t] + B10_[t] * data_[6+t] + B00_[t] * data_[2+t];
}

