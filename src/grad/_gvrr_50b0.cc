//
// Newint - Parallel electron correlation program.
// Filename: _gvrr_50b0.cc
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

// returns double array of length 576
void GVRRList::_gvrr_50b0(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
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
    data_[48+t] = D00_[t];

  double cB00_current[8];
  for (int t = 0; t != 8; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[56+t] = C00_[t] * data_[48+t] + cB00_current[t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[64+t] = C00_[t] * data_[56+t] + B10_current[t] * data_[48+t] + cB00_current[t] * data_[8+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[72+t] = C00_[t] * data_[64+t] + B10_current[t] * data_[56+t] + cB00_current[t] * data_[16+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[80+t] = C00_[t] * data_[72+t] + B10_current[t] * data_[64+t] + cB00_current[t] * data_[24+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[88+t] = C00_[t] * data_[80+t] + B10_current[t] * data_[72+t] + cB00_current[t] * data_[32+t];

  double B01_current[8];
  for (int t = 0; t != 8; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[96+t] = D00_[t] * data_[48+t] + B01_current[t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[104+t] = C00_[t] * data_[96+t] + cB00_current[t] * data_[48+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[112+t] = C00_[t] * data_[104+t] + B10_current[t] * data_[96+t] + cB00_current[t] * data_[56+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[120+t] = C00_[t] * data_[112+t] + B10_current[t] * data_[104+t] + cB00_current[t] * data_[64+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[128+t] = C00_[t] * data_[120+t] + B10_current[t] * data_[112+t] + cB00_current[t] * data_[72+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[136+t] = C00_[t] * data_[128+t] + B10_current[t] * data_[120+t] + cB00_current[t] * data_[80+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[144+t] = D00_[t] * data_[96+t] + B01_current[t] * data_[48+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[152+t] = C00_[t] * data_[144+t] + cB00_current[t] * data_[96+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[160+t] = C00_[t] * data_[152+t] + B10_current[t] * data_[144+t] + cB00_current[t] * data_[104+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[168+t] = C00_[t] * data_[160+t] + B10_current[t] * data_[152+t] + cB00_current[t] * data_[112+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[176+t] = C00_[t] * data_[168+t] + B10_current[t] * data_[160+t] + cB00_current[t] * data_[120+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[184+t] = C00_[t] * data_[176+t] + B10_current[t] * data_[168+t] + cB00_current[t] * data_[128+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[192+t] = D00_[t] * data_[144+t] + B01_current[t] * data_[96+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[200+t] = C00_[t] * data_[192+t] + cB00_current[t] * data_[144+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[208+t] = C00_[t] * data_[200+t] + B10_current[t] * data_[192+t] + cB00_current[t] * data_[152+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[216+t] = C00_[t] * data_[208+t] + B10_current[t] * data_[200+t] + cB00_current[t] * data_[160+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[224+t] = C00_[t] * data_[216+t] + B10_current[t] * data_[208+t] + cB00_current[t] * data_[168+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[232+t] = C00_[t] * data_[224+t] + B10_current[t] * data_[216+t] + cB00_current[t] * data_[176+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[240+t] = D00_[t] * data_[192+t] + B01_current[t] * data_[144+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[248+t] = C00_[t] * data_[240+t] + cB00_current[t] * data_[192+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[256+t] = C00_[t] * data_[248+t] + B10_current[t] * data_[240+t] + cB00_current[t] * data_[200+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[264+t] = C00_[t] * data_[256+t] + B10_current[t] * data_[248+t] + cB00_current[t] * data_[208+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[272+t] = C00_[t] * data_[264+t] + B10_current[t] * data_[256+t] + cB00_current[t] * data_[216+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[280+t] = C00_[t] * data_[272+t] + B10_current[t] * data_[264+t] + cB00_current[t] * data_[224+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[288+t] = D00_[t] * data_[240+t] + B01_current[t] * data_[192+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[296+t] = C00_[t] * data_[288+t] + cB00_current[t] * data_[240+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[304+t] = C00_[t] * data_[296+t] + B10_current[t] * data_[288+t] + cB00_current[t] * data_[248+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[312+t] = C00_[t] * data_[304+t] + B10_current[t] * data_[296+t] + cB00_current[t] * data_[256+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[320+t] = C00_[t] * data_[312+t] + B10_current[t] * data_[304+t] + cB00_current[t] * data_[264+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[328+t] = C00_[t] * data_[320+t] + B10_current[t] * data_[312+t] + cB00_current[t] * data_[272+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[336+t] = D00_[t] * data_[288+t] + B01_current[t] * data_[240+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[344+t] = C00_[t] * data_[336+t] + cB00_current[t] * data_[288+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[352+t] = C00_[t] * data_[344+t] + B10_current[t] * data_[336+t] + cB00_current[t] * data_[296+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[360+t] = C00_[t] * data_[352+t] + B10_current[t] * data_[344+t] + cB00_current[t] * data_[304+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[368+t] = C00_[t] * data_[360+t] + B10_current[t] * data_[352+t] + cB00_current[t] * data_[312+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[376+t] = C00_[t] * data_[368+t] + B10_current[t] * data_[360+t] + cB00_current[t] * data_[320+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[384+t] = D00_[t] * data_[336+t] + B01_current[t] * data_[288+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[392+t] = C00_[t] * data_[384+t] + cB00_current[t] * data_[336+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[400+t] = C00_[t] * data_[392+t] + B10_current[t] * data_[384+t] + cB00_current[t] * data_[344+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[408+t] = C00_[t] * data_[400+t] + B10_current[t] * data_[392+t] + cB00_current[t] * data_[352+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[416+t] = C00_[t] * data_[408+t] + B10_current[t] * data_[400+t] + cB00_current[t] * data_[360+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[424+t] = C00_[t] * data_[416+t] + B10_current[t] * data_[408+t] + cB00_current[t] * data_[368+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[432+t] = D00_[t] * data_[384+t] + B01_current[t] * data_[336+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[440+t] = C00_[t] * data_[432+t] + cB00_current[t] * data_[384+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[448+t] = C00_[t] * data_[440+t] + B10_current[t] * data_[432+t] + cB00_current[t] * data_[392+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[456+t] = C00_[t] * data_[448+t] + B10_current[t] * data_[440+t] + cB00_current[t] * data_[400+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[464+t] = C00_[t] * data_[456+t] + B10_current[t] * data_[448+t] + cB00_current[t] * data_[408+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[472+t] = C00_[t] * data_[464+t] + B10_current[t] * data_[456+t] + cB00_current[t] * data_[416+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[480+t] = D00_[t] * data_[432+t] + B01_current[t] * data_[384+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[488+t] = C00_[t] * data_[480+t] + cB00_current[t] * data_[432+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[496+t] = C00_[t] * data_[488+t] + B10_current[t] * data_[480+t] + cB00_current[t] * data_[440+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[504+t] = C00_[t] * data_[496+t] + B10_current[t] * data_[488+t] + cB00_current[t] * data_[448+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[512+t] = C00_[t] * data_[504+t] + B10_current[t] * data_[496+t] + cB00_current[t] * data_[456+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[520+t] = C00_[t] * data_[512+t] + B10_current[t] * data_[504+t] + cB00_current[t] * data_[464+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[528+t] = D00_[t] * data_[480+t] + B01_current[t] * data_[432+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[536+t] = C00_[t] * data_[528+t] + cB00_current[t] * data_[480+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[544+t] = C00_[t] * data_[536+t] + B10_current[t] * data_[528+t] + cB00_current[t] * data_[488+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[552+t] = C00_[t] * data_[544+t] + B10_current[t] * data_[536+t] + cB00_current[t] * data_[496+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[560+t] = C00_[t] * data_[552+t] + B10_current[t] * data_[544+t] + cB00_current[t] * data_[504+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[568+t] = C00_[t] * data_[560+t] + B10_current[t] * data_[552+t] + cB00_current[t] * data_[512+t];
}

