//
// Newint - Parallel electron correlation program.
// Filename: _gvrr_7080.cc
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


#include <src/grad/gvrrlist.h>

// returns double array of length 576
void GVRRList::_gvrr_7080(double* data_, const double* C00, const double* D00, const double* B00, const double* B01, const double* B10) {
#ifdef __GNUC__
  double C00_[8]__attribute__((aligned(16))) = {C00[0], C00[1], C00[2], C00[3], C00[4], C00[5], C00[6], C00[7]};
  double D00_[8]__attribute__((aligned(16))) = {D00[0], D00[1], D00[2], D00[3], D00[4], D00[5], D00[6], D00[7]};
  double B00_[8]__attribute__((aligned(16))) = {B00[0], B00[1], B00[2], B00[3], B00[4], B00[5], B00[6], B00[7]};
  double B01_[8]__attribute__((aligned(16))) = {B01[0], B01[1], B01[2], B01[3], B01[4], B01[5], B01[6], B01[7]};
  double B10_[8]__attribute__((aligned(16))) = {B10[0], B10[1], B10[2], B10[3], B10[4], B10[5], B10[6], B10[7]};
#else
  double* C00_ = C00;
  double* D00_ = D00;
  double* B00_ = B00;
  double* B01_ = B01;
  double* B10_ = B10;
#endif

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
    data_[64+t] = D00_[t];

  double cB00_current[8];
  for (int t = 0; t != 8; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[72+t] = C00_[t] * data_[64+t] + cB00_current[t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[80+t] = C00_[t] * data_[72+t] + B10_current[t] * data_[64+t] + cB00_current[t] * data_[8+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[88+t] = C00_[t] * data_[80+t] + B10_current[t] * data_[72+t] + cB00_current[t] * data_[16+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[96+t] = C00_[t] * data_[88+t] + B10_current[t] * data_[80+t] + cB00_current[t] * data_[24+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[104+t] = C00_[t] * data_[96+t] + B10_current[t] * data_[88+t] + cB00_current[t] * data_[32+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[112+t] = C00_[t] * data_[104+t] + B10_current[t] * data_[96+t] + cB00_current[t] * data_[40+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[120+t] = C00_[t] * data_[112+t] + B10_current[t] * data_[104+t] + cB00_current[t] * data_[48+t];

  double B01_current[8];
  for (int t = 0; t != 8; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[128+t] = D00_[t] * data_[64+t] + B01_current[t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[136+t] = C00_[t] * data_[128+t] + cB00_current[t] * data_[64+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[144+t] = C00_[t] * data_[136+t] + B10_current[t] * data_[128+t] + cB00_current[t] * data_[72+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[152+t] = C00_[t] * data_[144+t] + B10_current[t] * data_[136+t] + cB00_current[t] * data_[80+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[160+t] = C00_[t] * data_[152+t] + B10_current[t] * data_[144+t] + cB00_current[t] * data_[88+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[168+t] = C00_[t] * data_[160+t] + B10_current[t] * data_[152+t] + cB00_current[t] * data_[96+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[176+t] = C00_[t] * data_[168+t] + B10_current[t] * data_[160+t] + cB00_current[t] * data_[104+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[184+t] = C00_[t] * data_[176+t] + B10_current[t] * data_[168+t] + cB00_current[t] * data_[112+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[192+t] = D00_[t] * data_[128+t] + B01_current[t] * data_[64+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[200+t] = C00_[t] * data_[192+t] + cB00_current[t] * data_[128+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[208+t] = C00_[t] * data_[200+t] + B10_current[t] * data_[192+t] + cB00_current[t] * data_[136+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[216+t] = C00_[t] * data_[208+t] + B10_current[t] * data_[200+t] + cB00_current[t] * data_[144+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[224+t] = C00_[t] * data_[216+t] + B10_current[t] * data_[208+t] + cB00_current[t] * data_[152+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[232+t] = C00_[t] * data_[224+t] + B10_current[t] * data_[216+t] + cB00_current[t] * data_[160+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[240+t] = C00_[t] * data_[232+t] + B10_current[t] * data_[224+t] + cB00_current[t] * data_[168+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[248+t] = C00_[t] * data_[240+t] + B10_current[t] * data_[232+t] + cB00_current[t] * data_[176+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[256+t] = D00_[t] * data_[192+t] + B01_current[t] * data_[128+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[264+t] = C00_[t] * data_[256+t] + cB00_current[t] * data_[192+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[272+t] = C00_[t] * data_[264+t] + B10_current[t] * data_[256+t] + cB00_current[t] * data_[200+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[280+t] = C00_[t] * data_[272+t] + B10_current[t] * data_[264+t] + cB00_current[t] * data_[208+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[288+t] = C00_[t] * data_[280+t] + B10_current[t] * data_[272+t] + cB00_current[t] * data_[216+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[296+t] = C00_[t] * data_[288+t] + B10_current[t] * data_[280+t] + cB00_current[t] * data_[224+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[304+t] = C00_[t] * data_[296+t] + B10_current[t] * data_[288+t] + cB00_current[t] * data_[232+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[312+t] = C00_[t] * data_[304+t] + B10_current[t] * data_[296+t] + cB00_current[t] * data_[240+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[320+t] = D00_[t] * data_[256+t] + B01_current[t] * data_[192+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[328+t] = C00_[t] * data_[320+t] + cB00_current[t] * data_[256+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[336+t] = C00_[t] * data_[328+t] + B10_current[t] * data_[320+t] + cB00_current[t] * data_[264+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[344+t] = C00_[t] * data_[336+t] + B10_current[t] * data_[328+t] + cB00_current[t] * data_[272+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[352+t] = C00_[t] * data_[344+t] + B10_current[t] * data_[336+t] + cB00_current[t] * data_[280+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[360+t] = C00_[t] * data_[352+t] + B10_current[t] * data_[344+t] + cB00_current[t] * data_[288+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[368+t] = C00_[t] * data_[360+t] + B10_current[t] * data_[352+t] + cB00_current[t] * data_[296+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[376+t] = C00_[t] * data_[368+t] + B10_current[t] * data_[360+t] + cB00_current[t] * data_[304+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[384+t] = D00_[t] * data_[320+t] + B01_current[t] * data_[256+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[392+t] = C00_[t] * data_[384+t] + cB00_current[t] * data_[320+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[400+t] = C00_[t] * data_[392+t] + B10_current[t] * data_[384+t] + cB00_current[t] * data_[328+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[408+t] = C00_[t] * data_[400+t] + B10_current[t] * data_[392+t] + cB00_current[t] * data_[336+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[416+t] = C00_[t] * data_[408+t] + B10_current[t] * data_[400+t] + cB00_current[t] * data_[344+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[424+t] = C00_[t] * data_[416+t] + B10_current[t] * data_[408+t] + cB00_current[t] * data_[352+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[432+t] = C00_[t] * data_[424+t] + B10_current[t] * data_[416+t] + cB00_current[t] * data_[360+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[440+t] = C00_[t] * data_[432+t] + B10_current[t] * data_[424+t] + cB00_current[t] * data_[368+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[448+t] = D00_[t] * data_[384+t] + B01_current[t] * data_[320+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[456+t] = C00_[t] * data_[448+t] + cB00_current[t] * data_[384+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[464+t] = C00_[t] * data_[456+t] + B10_current[t] * data_[448+t] + cB00_current[t] * data_[392+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[472+t] = C00_[t] * data_[464+t] + B10_current[t] * data_[456+t] + cB00_current[t] * data_[400+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[480+t] = C00_[t] * data_[472+t] + B10_current[t] * data_[464+t] + cB00_current[t] * data_[408+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[488+t] = C00_[t] * data_[480+t] + B10_current[t] * data_[472+t] + cB00_current[t] * data_[416+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[496+t] = C00_[t] * data_[488+t] + B10_current[t] * data_[480+t] + cB00_current[t] * data_[424+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[504+t] = C00_[t] * data_[496+t] + B10_current[t] * data_[488+t] + cB00_current[t] * data_[432+t];

  for (int t = 0; t != 8; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 8; ++t)
    data_[512+t] = D00_[t] * data_[448+t] + B01_current[t] * data_[384+t];

  for (int t = 0; t != 8; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 8; ++t)
    data_[520+t] = C00_[t] * data_[512+t] + cB00_current[t] * data_[448+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[528+t] = C00_[t] * data_[520+t] + B10_current[t] * data_[512+t] + cB00_current[t] * data_[456+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[536+t] = C00_[t] * data_[528+t] + B10_current[t] * data_[520+t] + cB00_current[t] * data_[464+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[544+t] = C00_[t] * data_[536+t] + B10_current[t] * data_[528+t] + cB00_current[t] * data_[472+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[552+t] = C00_[t] * data_[544+t] + B10_current[t] * data_[536+t] + cB00_current[t] * data_[480+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[560+t] = C00_[t] * data_[552+t] + B10_current[t] * data_[544+t] + cB00_current[t] * data_[488+t];

  for (int t = 0; t != 8; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 8; ++t)
    data_[568+t] = C00_[t] * data_[560+t] + B10_current[t] * data_[552+t] + cB00_current[t] * data_[496+t];
}

