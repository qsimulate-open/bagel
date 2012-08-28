//
// BAGEL - Parallel electron correlation program.
// Filename: _svrr_20c0.cc
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


#include "svrrlist.h"

// returns double array of length 312
void SVRRList::_svrr_20c0(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  for (int t = 0; t != 8; ++t)
    data_[0+t] = 1.0;

  for (int t = 0; t != 8; ++t)
    data_[8+t] = C00_[t];

  for (int t = 0; t != 8; ++t)
    data_[16+t] = C00_[t] * data_[8+t] + B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[24+t] = D00_[t];

  double cB00_current[8];
  for (int t = 0; t != 8; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[32+t] = C00_[t] * data_[24+t] + cB00_current[t];

  for (int t = 0; t != 8; ++t)
    data_[40+t] = C00_[t] * data_[32+t] + B10_[t] * data_[24+t] + cB00_current[t] * data_[8+t];

  double B01_current[8];
  for (int t = 0; t != 8; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[48+t] = D00_[t] * data_[24+t] + B01_current[t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[56+t] = C00_[t] * data_[48+t] + cB00_current[t] * data_[24+t];

  for (int t = 0; t != 8; ++t)
    data_[64+t] = C00_[t] * data_[56+t] + B10_[t] * data_[48+t] + cB00_current[t] * data_[32+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[72+t] = D00_[t] * data_[48+t] + B01_current[t] * data_[24+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[80+t] = C00_[t] * data_[72+t] + cB00_current[t] * data_[48+t];

  for (int t = 0; t != 8; ++t)
    data_[88+t] = C00_[t] * data_[80+t] + B10_[t] * data_[72+t] + cB00_current[t] * data_[56+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[96+t] = D00_[t] * data_[72+t] + B01_current[t] * data_[48+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[104+t] = C00_[t] * data_[96+t] + cB00_current[t] * data_[72+t];

  for (int t = 0; t != 8; ++t)
    data_[112+t] = C00_[t] * data_[104+t] + B10_[t] * data_[96+t] + cB00_current[t] * data_[80+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[120+t] = D00_[t] * data_[96+t] + B01_current[t] * data_[72+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[128+t] = C00_[t] * data_[120+t] + cB00_current[t] * data_[96+t];

  for (int t = 0; t != 8; ++t)
    data_[136+t] = C00_[t] * data_[128+t] + B10_[t] * data_[120+t] + cB00_current[t] * data_[104+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[144+t] = D00_[t] * data_[120+t] + B01_current[t] * data_[96+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[152+t] = C00_[t] * data_[144+t] + cB00_current[t] * data_[120+t];

  for (int t = 0; t != 8; ++t)
    data_[160+t] = C00_[t] * data_[152+t] + B10_[t] * data_[144+t] + cB00_current[t] * data_[128+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[168+t] = D00_[t] * data_[144+t] + B01_current[t] * data_[120+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[176+t] = C00_[t] * data_[168+t] + cB00_current[t] * data_[144+t];

  for (int t = 0; t != 8; ++t)
    data_[184+t] = C00_[t] * data_[176+t] + B10_[t] * data_[168+t] + cB00_current[t] * data_[152+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[192+t] = D00_[t] * data_[168+t] + B01_current[t] * data_[144+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[200+t] = C00_[t] * data_[192+t] + cB00_current[t] * data_[168+t];

  for (int t = 0; t != 8; ++t)
    data_[208+t] = C00_[t] * data_[200+t] + B10_[t] * data_[192+t] + cB00_current[t] * data_[176+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[216+t] = D00_[t] * data_[192+t] + B01_current[t] * data_[168+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[224+t] = C00_[t] * data_[216+t] + cB00_current[t] * data_[192+t];

  for (int t = 0; t != 8; ++t)
    data_[232+t] = C00_[t] * data_[224+t] + B10_[t] * data_[216+t] + cB00_current[t] * data_[200+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[240+t] = D00_[t] * data_[216+t] + B01_current[t] * data_[192+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[248+t] = C00_[t] * data_[240+t] + cB00_current[t] * data_[216+t];

  for (int t = 0; t != 8; ++t)
    data_[256+t] = C00_[t] * data_[248+t] + B10_[t] * data_[240+t] + cB00_current[t] * data_[224+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[264+t] = D00_[t] * data_[240+t] + B01_current[t] * data_[216+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[272+t] = C00_[t] * data_[264+t] + cB00_current[t] * data_[240+t];

  for (int t = 0; t != 8; ++t)
    data_[280+t] = C00_[t] * data_[272+t] + B10_[t] * data_[264+t] + cB00_current[t] * data_[248+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[288+t] = D00_[t] * data_[264+t] + B01_current[t] * data_[240+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[296+t] = C00_[t] * data_[288+t] + cB00_current[t] * data_[264+t];

  for (int t = 0; t != 8; ++t)
    data_[304+t] = C00_[t] * data_[296+t] + B10_[t] * data_[288+t] + cB00_current[t] * data_[272+t];
}

