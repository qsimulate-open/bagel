//
// Newint - Parallel electron correlation program.
// Filename: _vrr_c020.cc
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

// returns double array of length 312
void VRRList::_vrr_c020(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  for (int t = 0; t != 8; ++t)
    data_[0+t] = 1.0;

  for (int t = 0; t != 8; ++t)
    data_[8+t] = C00_[t];

  double B10_current[8];
  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[16+t] = C00_[t] * data_[8+t] + B10_current[t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[24+t] = C00_[t] * data_[16+t] + B10_current[t] * data_[8+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[32+t] = C00_[t] * data_[24+t] + B10_current[t] * data_[16+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[40+t] = C00_[t] * data_[32+t] + B10_current[t] * data_[24+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[48+t] = C00_[t] * data_[40+t] + B10_current[t] * data_[32+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[56+t] = C00_[t] * data_[48+t] + B10_current[t] * data_[40+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[64+t] = C00_[t] * data_[56+t] + B10_current[t] * data_[48+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[72+t] = C00_[t] * data_[64+t] + B10_current[t] * data_[56+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[80+t] = C00_[t] * data_[72+t] + B10_current[t] * data_[64+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[88+t] = C00_[t] * data_[80+t] + B10_current[t] * data_[72+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[96+t] = C00_[t] * data_[88+t] + B10_current[t] * data_[80+t];

  for (int t = 0; t != 8; ++t)
    data_[104+t] = D00_[t];

  double cB00_current[8];
  for (int t = 0; t != 8; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[112+t] = C00_[t] * data_[104+t] + cB00_current[t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[120+t] = C00_[t] * data_[112+t] + B10_current[t] * data_[104+t] + cB00_current[t] * data_[8+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[128+t] = C00_[t] * data_[120+t] + B10_current[t] * data_[112+t] + cB00_current[t] * data_[16+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[136+t] = C00_[t] * data_[128+t] + B10_current[t] * data_[120+t] + cB00_current[t] * data_[24+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[144+t] = C00_[t] * data_[136+t] + B10_current[t] * data_[128+t] + cB00_current[t] * data_[32+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[152+t] = C00_[t] * data_[144+t] + B10_current[t] * data_[136+t] + cB00_current[t] * data_[40+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[160+t] = C00_[t] * data_[152+t] + B10_current[t] * data_[144+t] + cB00_current[t] * data_[48+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[168+t] = C00_[t] * data_[160+t] + B10_current[t] * data_[152+t] + cB00_current[t] * data_[56+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[176+t] = C00_[t] * data_[168+t] + B10_current[t] * data_[160+t] + cB00_current[t] * data_[64+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[184+t] = C00_[t] * data_[176+t] + B10_current[t] * data_[168+t] + cB00_current[t] * data_[72+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[192+t] = C00_[t] * data_[184+t] + B10_current[t] * data_[176+t] + cB00_current[t] * data_[80+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[200+t] = C00_[t] * data_[192+t] + B10_current[t] * data_[184+t] + cB00_current[t] * data_[88+t];


  for (int t = 0; t != 8; ++t)
    data_[208+t] = D00_[t] * data_[104+t] + B01_[t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[216+t] = C00_[t] * data_[208+t] + cB00_current[t] * data_[104+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[224+t] = C00_[t] * data_[216+t] + B10_current[t] * data_[208+t] + cB00_current[t] * data_[112+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[232+t] = C00_[t] * data_[224+t] + B10_current[t] * data_[216+t] + cB00_current[t] * data_[120+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[240+t] = C00_[t] * data_[232+t] + B10_current[t] * data_[224+t] + cB00_current[t] * data_[128+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[248+t] = C00_[t] * data_[240+t] + B10_current[t] * data_[232+t] + cB00_current[t] * data_[136+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[256+t] = C00_[t] * data_[248+t] + B10_current[t] * data_[240+t] + cB00_current[t] * data_[144+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[264+t] = C00_[t] * data_[256+t] + B10_current[t] * data_[248+t] + cB00_current[t] * data_[152+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[272+t] = C00_[t] * data_[264+t] + B10_current[t] * data_[256+t] + cB00_current[t] * data_[160+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[280+t] = C00_[t] * data_[272+t] + B10_current[t] * data_[264+t] + cB00_current[t] * data_[168+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[288+t] = C00_[t] * data_[280+t] + B10_current[t] * data_[272+t] + cB00_current[t] * data_[176+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[296+t] = C00_[t] * data_[288+t] + B10_current[t] * data_[280+t] + cB00_current[t] * data_[184+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[304+t] = C00_[t] * data_[296+t] + B10_current[t] * data_[288+t] + cB00_current[t] * data_[192+t];
}

