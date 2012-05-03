//
// Newint - Parallel electron correlation program.
// Filename: _gvrr_d0c0.cc
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

// returns double array of length 2366
void GVRRList::_gvrr_d0c0(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  for (int t = 0; t != 13; ++t)
    data_[0+t] = 1.0;

  for (int t = 0; t != 13; ++t)
    data_[13+t] = C00_[t];

  double B10_current[13];
  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[26+t] = C00_[t] * data_[13+t] + B10_current[t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[39+t] = C00_[t] * data_[26+t] + B10_current[t] * data_[13+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[52+t] = C00_[t] * data_[39+t] + B10_current[t] * data_[26+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[65+t] = C00_[t] * data_[52+t] + B10_current[t] * data_[39+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[78+t] = C00_[t] * data_[65+t] + B10_current[t] * data_[52+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[91+t] = C00_[t] * data_[78+t] + B10_current[t] * data_[65+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[104+t] = C00_[t] * data_[91+t] + B10_current[t] * data_[78+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[117+t] = C00_[t] * data_[104+t] + B10_current[t] * data_[91+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[130+t] = C00_[t] * data_[117+t] + B10_current[t] * data_[104+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[143+t] = C00_[t] * data_[130+t] + B10_current[t] * data_[117+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[156+t] = C00_[t] * data_[143+t] + B10_current[t] * data_[130+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[169+t] = C00_[t] * data_[156+t] + B10_current[t] * data_[143+t];

  for (int t = 0; t != 13; ++t)
    data_[182+t] = D00_[t];

  double cB00_current[13];
  for (int t = 0; t != 13; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 13; ++t)
    data_[195+t] = C00_[t] * data_[182+t] + cB00_current[t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[208+t] = C00_[t] * data_[195+t] + B10_current[t] * data_[182+t] + cB00_current[t] * data_[13+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[221+t] = C00_[t] * data_[208+t] + B10_current[t] * data_[195+t] + cB00_current[t] * data_[26+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[234+t] = C00_[t] * data_[221+t] + B10_current[t] * data_[208+t] + cB00_current[t] * data_[39+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[247+t] = C00_[t] * data_[234+t] + B10_current[t] * data_[221+t] + cB00_current[t] * data_[52+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[260+t] = C00_[t] * data_[247+t] + B10_current[t] * data_[234+t] + cB00_current[t] * data_[65+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[273+t] = C00_[t] * data_[260+t] + B10_current[t] * data_[247+t] + cB00_current[t] * data_[78+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[286+t] = C00_[t] * data_[273+t] + B10_current[t] * data_[260+t] + cB00_current[t] * data_[91+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[299+t] = C00_[t] * data_[286+t] + B10_current[t] * data_[273+t] + cB00_current[t] * data_[104+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[312+t] = C00_[t] * data_[299+t] + B10_current[t] * data_[286+t] + cB00_current[t] * data_[117+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[325+t] = C00_[t] * data_[312+t] + B10_current[t] * data_[299+t] + cB00_current[t] * data_[130+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[338+t] = C00_[t] * data_[325+t] + B10_current[t] * data_[312+t] + cB00_current[t] * data_[143+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[351+t] = C00_[t] * data_[338+t] + B10_current[t] * data_[325+t] + cB00_current[t] * data_[156+t];

  double B01_current[13];
  for (int t = 0; t != 13; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 13; ++t)
    data_[364+t] = D00_[t] * data_[182+t] + B01_current[t];

  for (int t = 0; t != 13; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 13; ++t)
    data_[377+t] = C00_[t] * data_[364+t] + cB00_current[t] * data_[182+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[390+t] = C00_[t] * data_[377+t] + B10_current[t] * data_[364+t] + cB00_current[t] * data_[195+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[403+t] = C00_[t] * data_[390+t] + B10_current[t] * data_[377+t] + cB00_current[t] * data_[208+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[416+t] = C00_[t] * data_[403+t] + B10_current[t] * data_[390+t] + cB00_current[t] * data_[221+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[429+t] = C00_[t] * data_[416+t] + B10_current[t] * data_[403+t] + cB00_current[t] * data_[234+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[442+t] = C00_[t] * data_[429+t] + B10_current[t] * data_[416+t] + cB00_current[t] * data_[247+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[455+t] = C00_[t] * data_[442+t] + B10_current[t] * data_[429+t] + cB00_current[t] * data_[260+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[468+t] = C00_[t] * data_[455+t] + B10_current[t] * data_[442+t] + cB00_current[t] * data_[273+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[481+t] = C00_[t] * data_[468+t] + B10_current[t] * data_[455+t] + cB00_current[t] * data_[286+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[494+t] = C00_[t] * data_[481+t] + B10_current[t] * data_[468+t] + cB00_current[t] * data_[299+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[507+t] = C00_[t] * data_[494+t] + B10_current[t] * data_[481+t] + cB00_current[t] * data_[312+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[520+t] = C00_[t] * data_[507+t] + B10_current[t] * data_[494+t] + cB00_current[t] * data_[325+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[533+t] = C00_[t] * data_[520+t] + B10_current[t] * data_[507+t] + cB00_current[t] * data_[338+t];

  for (int t = 0; t != 13; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 13; ++t)
    data_[546+t] = D00_[t] * data_[364+t] + B01_current[t] * data_[182+t];

  for (int t = 0; t != 13; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 13; ++t)
    data_[559+t] = C00_[t] * data_[546+t] + cB00_current[t] * data_[364+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[572+t] = C00_[t] * data_[559+t] + B10_current[t] * data_[546+t] + cB00_current[t] * data_[377+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[585+t] = C00_[t] * data_[572+t] + B10_current[t] * data_[559+t] + cB00_current[t] * data_[390+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[598+t] = C00_[t] * data_[585+t] + B10_current[t] * data_[572+t] + cB00_current[t] * data_[403+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[611+t] = C00_[t] * data_[598+t] + B10_current[t] * data_[585+t] + cB00_current[t] * data_[416+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[624+t] = C00_[t] * data_[611+t] + B10_current[t] * data_[598+t] + cB00_current[t] * data_[429+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[637+t] = C00_[t] * data_[624+t] + B10_current[t] * data_[611+t] + cB00_current[t] * data_[442+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[650+t] = C00_[t] * data_[637+t] + B10_current[t] * data_[624+t] + cB00_current[t] * data_[455+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[663+t] = C00_[t] * data_[650+t] + B10_current[t] * data_[637+t] + cB00_current[t] * data_[468+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[676+t] = C00_[t] * data_[663+t] + B10_current[t] * data_[650+t] + cB00_current[t] * data_[481+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[689+t] = C00_[t] * data_[676+t] + B10_current[t] * data_[663+t] + cB00_current[t] * data_[494+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[702+t] = C00_[t] * data_[689+t] + B10_current[t] * data_[676+t] + cB00_current[t] * data_[507+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[715+t] = C00_[t] * data_[702+t] + B10_current[t] * data_[689+t] + cB00_current[t] * data_[520+t];

  for (int t = 0; t != 13; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 13; ++t)
    data_[728+t] = D00_[t] * data_[546+t] + B01_current[t] * data_[364+t];

  for (int t = 0; t != 13; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 13; ++t)
    data_[741+t] = C00_[t] * data_[728+t] + cB00_current[t] * data_[546+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[754+t] = C00_[t] * data_[741+t] + B10_current[t] * data_[728+t] + cB00_current[t] * data_[559+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[767+t] = C00_[t] * data_[754+t] + B10_current[t] * data_[741+t] + cB00_current[t] * data_[572+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[780+t] = C00_[t] * data_[767+t] + B10_current[t] * data_[754+t] + cB00_current[t] * data_[585+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[793+t] = C00_[t] * data_[780+t] + B10_current[t] * data_[767+t] + cB00_current[t] * data_[598+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[806+t] = C00_[t] * data_[793+t] + B10_current[t] * data_[780+t] + cB00_current[t] * data_[611+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[819+t] = C00_[t] * data_[806+t] + B10_current[t] * data_[793+t] + cB00_current[t] * data_[624+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[832+t] = C00_[t] * data_[819+t] + B10_current[t] * data_[806+t] + cB00_current[t] * data_[637+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[845+t] = C00_[t] * data_[832+t] + B10_current[t] * data_[819+t] + cB00_current[t] * data_[650+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[858+t] = C00_[t] * data_[845+t] + B10_current[t] * data_[832+t] + cB00_current[t] * data_[663+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[871+t] = C00_[t] * data_[858+t] + B10_current[t] * data_[845+t] + cB00_current[t] * data_[676+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[884+t] = C00_[t] * data_[871+t] + B10_current[t] * data_[858+t] + cB00_current[t] * data_[689+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[897+t] = C00_[t] * data_[884+t] + B10_current[t] * data_[871+t] + cB00_current[t] * data_[702+t];

  for (int t = 0; t != 13; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 13; ++t)
    data_[910+t] = D00_[t] * data_[728+t] + B01_current[t] * data_[546+t];

  for (int t = 0; t != 13; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 13; ++t)
    data_[923+t] = C00_[t] * data_[910+t] + cB00_current[t] * data_[728+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[936+t] = C00_[t] * data_[923+t] + B10_current[t] * data_[910+t] + cB00_current[t] * data_[741+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[949+t] = C00_[t] * data_[936+t] + B10_current[t] * data_[923+t] + cB00_current[t] * data_[754+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[962+t] = C00_[t] * data_[949+t] + B10_current[t] * data_[936+t] + cB00_current[t] * data_[767+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[975+t] = C00_[t] * data_[962+t] + B10_current[t] * data_[949+t] + cB00_current[t] * data_[780+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[988+t] = C00_[t] * data_[975+t] + B10_current[t] * data_[962+t] + cB00_current[t] * data_[793+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1001+t] = C00_[t] * data_[988+t] + B10_current[t] * data_[975+t] + cB00_current[t] * data_[806+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1014+t] = C00_[t] * data_[1001+t] + B10_current[t] * data_[988+t] + cB00_current[t] * data_[819+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1027+t] = C00_[t] * data_[1014+t] + B10_current[t] * data_[1001+t] + cB00_current[t] * data_[832+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1040+t] = C00_[t] * data_[1027+t] + B10_current[t] * data_[1014+t] + cB00_current[t] * data_[845+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1053+t] = C00_[t] * data_[1040+t] + B10_current[t] * data_[1027+t] + cB00_current[t] * data_[858+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1066+t] = C00_[t] * data_[1053+t] + B10_current[t] * data_[1040+t] + cB00_current[t] * data_[871+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1079+t] = C00_[t] * data_[1066+t] + B10_current[t] * data_[1053+t] + cB00_current[t] * data_[884+t];

  for (int t = 0; t != 13; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 13; ++t)
    data_[1092+t] = D00_[t] * data_[910+t] + B01_current[t] * data_[728+t];

  for (int t = 0; t != 13; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 13; ++t)
    data_[1105+t] = C00_[t] * data_[1092+t] + cB00_current[t] * data_[910+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1118+t] = C00_[t] * data_[1105+t] + B10_current[t] * data_[1092+t] + cB00_current[t] * data_[923+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1131+t] = C00_[t] * data_[1118+t] + B10_current[t] * data_[1105+t] + cB00_current[t] * data_[936+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1144+t] = C00_[t] * data_[1131+t] + B10_current[t] * data_[1118+t] + cB00_current[t] * data_[949+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1157+t] = C00_[t] * data_[1144+t] + B10_current[t] * data_[1131+t] + cB00_current[t] * data_[962+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1170+t] = C00_[t] * data_[1157+t] + B10_current[t] * data_[1144+t] + cB00_current[t] * data_[975+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1183+t] = C00_[t] * data_[1170+t] + B10_current[t] * data_[1157+t] + cB00_current[t] * data_[988+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1196+t] = C00_[t] * data_[1183+t] + B10_current[t] * data_[1170+t] + cB00_current[t] * data_[1001+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1209+t] = C00_[t] * data_[1196+t] + B10_current[t] * data_[1183+t] + cB00_current[t] * data_[1014+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1222+t] = C00_[t] * data_[1209+t] + B10_current[t] * data_[1196+t] + cB00_current[t] * data_[1027+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1235+t] = C00_[t] * data_[1222+t] + B10_current[t] * data_[1209+t] + cB00_current[t] * data_[1040+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1248+t] = C00_[t] * data_[1235+t] + B10_current[t] * data_[1222+t] + cB00_current[t] * data_[1053+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1261+t] = C00_[t] * data_[1248+t] + B10_current[t] * data_[1235+t] + cB00_current[t] * data_[1066+t];

  for (int t = 0; t != 13; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 13; ++t)
    data_[1274+t] = D00_[t] * data_[1092+t] + B01_current[t] * data_[910+t];

  for (int t = 0; t != 13; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 13; ++t)
    data_[1287+t] = C00_[t] * data_[1274+t] + cB00_current[t] * data_[1092+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1300+t] = C00_[t] * data_[1287+t] + B10_current[t] * data_[1274+t] + cB00_current[t] * data_[1105+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1313+t] = C00_[t] * data_[1300+t] + B10_current[t] * data_[1287+t] + cB00_current[t] * data_[1118+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1326+t] = C00_[t] * data_[1313+t] + B10_current[t] * data_[1300+t] + cB00_current[t] * data_[1131+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1339+t] = C00_[t] * data_[1326+t] + B10_current[t] * data_[1313+t] + cB00_current[t] * data_[1144+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1352+t] = C00_[t] * data_[1339+t] + B10_current[t] * data_[1326+t] + cB00_current[t] * data_[1157+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1365+t] = C00_[t] * data_[1352+t] + B10_current[t] * data_[1339+t] + cB00_current[t] * data_[1170+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1378+t] = C00_[t] * data_[1365+t] + B10_current[t] * data_[1352+t] + cB00_current[t] * data_[1183+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1391+t] = C00_[t] * data_[1378+t] + B10_current[t] * data_[1365+t] + cB00_current[t] * data_[1196+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1404+t] = C00_[t] * data_[1391+t] + B10_current[t] * data_[1378+t] + cB00_current[t] * data_[1209+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1417+t] = C00_[t] * data_[1404+t] + B10_current[t] * data_[1391+t] + cB00_current[t] * data_[1222+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1430+t] = C00_[t] * data_[1417+t] + B10_current[t] * data_[1404+t] + cB00_current[t] * data_[1235+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1443+t] = C00_[t] * data_[1430+t] + B10_current[t] * data_[1417+t] + cB00_current[t] * data_[1248+t];

  for (int t = 0; t != 13; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 13; ++t)
    data_[1456+t] = D00_[t] * data_[1274+t] + B01_current[t] * data_[1092+t];

  for (int t = 0; t != 13; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 13; ++t)
    data_[1469+t] = C00_[t] * data_[1456+t] + cB00_current[t] * data_[1274+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1482+t] = C00_[t] * data_[1469+t] + B10_current[t] * data_[1456+t] + cB00_current[t] * data_[1287+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1495+t] = C00_[t] * data_[1482+t] + B10_current[t] * data_[1469+t] + cB00_current[t] * data_[1300+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1508+t] = C00_[t] * data_[1495+t] + B10_current[t] * data_[1482+t] + cB00_current[t] * data_[1313+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1521+t] = C00_[t] * data_[1508+t] + B10_current[t] * data_[1495+t] + cB00_current[t] * data_[1326+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1534+t] = C00_[t] * data_[1521+t] + B10_current[t] * data_[1508+t] + cB00_current[t] * data_[1339+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1547+t] = C00_[t] * data_[1534+t] + B10_current[t] * data_[1521+t] + cB00_current[t] * data_[1352+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1560+t] = C00_[t] * data_[1547+t] + B10_current[t] * data_[1534+t] + cB00_current[t] * data_[1365+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1573+t] = C00_[t] * data_[1560+t] + B10_current[t] * data_[1547+t] + cB00_current[t] * data_[1378+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1586+t] = C00_[t] * data_[1573+t] + B10_current[t] * data_[1560+t] + cB00_current[t] * data_[1391+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1599+t] = C00_[t] * data_[1586+t] + B10_current[t] * data_[1573+t] + cB00_current[t] * data_[1404+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1612+t] = C00_[t] * data_[1599+t] + B10_current[t] * data_[1586+t] + cB00_current[t] * data_[1417+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1625+t] = C00_[t] * data_[1612+t] + B10_current[t] * data_[1599+t] + cB00_current[t] * data_[1430+t];

  for (int t = 0; t != 13; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 13; ++t)
    data_[1638+t] = D00_[t] * data_[1456+t] + B01_current[t] * data_[1274+t];

  for (int t = 0; t != 13; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 13; ++t)
    data_[1651+t] = C00_[t] * data_[1638+t] + cB00_current[t] * data_[1456+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1664+t] = C00_[t] * data_[1651+t] + B10_current[t] * data_[1638+t] + cB00_current[t] * data_[1469+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1677+t] = C00_[t] * data_[1664+t] + B10_current[t] * data_[1651+t] + cB00_current[t] * data_[1482+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1690+t] = C00_[t] * data_[1677+t] + B10_current[t] * data_[1664+t] + cB00_current[t] * data_[1495+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1703+t] = C00_[t] * data_[1690+t] + B10_current[t] * data_[1677+t] + cB00_current[t] * data_[1508+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1716+t] = C00_[t] * data_[1703+t] + B10_current[t] * data_[1690+t] + cB00_current[t] * data_[1521+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1729+t] = C00_[t] * data_[1716+t] + B10_current[t] * data_[1703+t] + cB00_current[t] * data_[1534+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1742+t] = C00_[t] * data_[1729+t] + B10_current[t] * data_[1716+t] + cB00_current[t] * data_[1547+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1755+t] = C00_[t] * data_[1742+t] + B10_current[t] * data_[1729+t] + cB00_current[t] * data_[1560+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1768+t] = C00_[t] * data_[1755+t] + B10_current[t] * data_[1742+t] + cB00_current[t] * data_[1573+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1781+t] = C00_[t] * data_[1768+t] + B10_current[t] * data_[1755+t] + cB00_current[t] * data_[1586+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1794+t] = C00_[t] * data_[1781+t] + B10_current[t] * data_[1768+t] + cB00_current[t] * data_[1599+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1807+t] = C00_[t] * data_[1794+t] + B10_current[t] * data_[1781+t] + cB00_current[t] * data_[1612+t];

  for (int t = 0; t != 13; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 13; ++t)
    data_[1820+t] = D00_[t] * data_[1638+t] + B01_current[t] * data_[1456+t];

  for (int t = 0; t != 13; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 13; ++t)
    data_[1833+t] = C00_[t] * data_[1820+t] + cB00_current[t] * data_[1638+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1846+t] = C00_[t] * data_[1833+t] + B10_current[t] * data_[1820+t] + cB00_current[t] * data_[1651+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1859+t] = C00_[t] * data_[1846+t] + B10_current[t] * data_[1833+t] + cB00_current[t] * data_[1664+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1872+t] = C00_[t] * data_[1859+t] + B10_current[t] * data_[1846+t] + cB00_current[t] * data_[1677+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1885+t] = C00_[t] * data_[1872+t] + B10_current[t] * data_[1859+t] + cB00_current[t] * data_[1690+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1898+t] = C00_[t] * data_[1885+t] + B10_current[t] * data_[1872+t] + cB00_current[t] * data_[1703+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1911+t] = C00_[t] * data_[1898+t] + B10_current[t] * data_[1885+t] + cB00_current[t] * data_[1716+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1924+t] = C00_[t] * data_[1911+t] + B10_current[t] * data_[1898+t] + cB00_current[t] * data_[1729+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1937+t] = C00_[t] * data_[1924+t] + B10_current[t] * data_[1911+t] + cB00_current[t] * data_[1742+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1950+t] = C00_[t] * data_[1937+t] + B10_current[t] * data_[1924+t] + cB00_current[t] * data_[1755+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1963+t] = C00_[t] * data_[1950+t] + B10_current[t] * data_[1937+t] + cB00_current[t] * data_[1768+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1976+t] = C00_[t] * data_[1963+t] + B10_current[t] * data_[1950+t] + cB00_current[t] * data_[1781+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[1989+t] = C00_[t] * data_[1976+t] + B10_current[t] * data_[1963+t] + cB00_current[t] * data_[1794+t];

  for (int t = 0; t != 13; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 13; ++t)
    data_[2002+t] = D00_[t] * data_[1820+t] + B01_current[t] * data_[1638+t];

  for (int t = 0; t != 13; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 13; ++t)
    data_[2015+t] = C00_[t] * data_[2002+t] + cB00_current[t] * data_[1820+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2028+t] = C00_[t] * data_[2015+t] + B10_current[t] * data_[2002+t] + cB00_current[t] * data_[1833+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2041+t] = C00_[t] * data_[2028+t] + B10_current[t] * data_[2015+t] + cB00_current[t] * data_[1846+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2054+t] = C00_[t] * data_[2041+t] + B10_current[t] * data_[2028+t] + cB00_current[t] * data_[1859+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2067+t] = C00_[t] * data_[2054+t] + B10_current[t] * data_[2041+t] + cB00_current[t] * data_[1872+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2080+t] = C00_[t] * data_[2067+t] + B10_current[t] * data_[2054+t] + cB00_current[t] * data_[1885+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2093+t] = C00_[t] * data_[2080+t] + B10_current[t] * data_[2067+t] + cB00_current[t] * data_[1898+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2106+t] = C00_[t] * data_[2093+t] + B10_current[t] * data_[2080+t] + cB00_current[t] * data_[1911+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2119+t] = C00_[t] * data_[2106+t] + B10_current[t] * data_[2093+t] + cB00_current[t] * data_[1924+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2132+t] = C00_[t] * data_[2119+t] + B10_current[t] * data_[2106+t] + cB00_current[t] * data_[1937+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2145+t] = C00_[t] * data_[2132+t] + B10_current[t] * data_[2119+t] + cB00_current[t] * data_[1950+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2158+t] = C00_[t] * data_[2145+t] + B10_current[t] * data_[2132+t] + cB00_current[t] * data_[1963+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2171+t] = C00_[t] * data_[2158+t] + B10_current[t] * data_[2145+t] + cB00_current[t] * data_[1976+t];

  for (int t = 0; t != 13; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 13; ++t)
    data_[2184+t] = D00_[t] * data_[2002+t] + B01_current[t] * data_[1820+t];

  for (int t = 0; t != 13; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 13; ++t)
    data_[2197+t] = C00_[t] * data_[2184+t] + cB00_current[t] * data_[2002+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2210+t] = C00_[t] * data_[2197+t] + B10_current[t] * data_[2184+t] + cB00_current[t] * data_[2015+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2223+t] = C00_[t] * data_[2210+t] + B10_current[t] * data_[2197+t] + cB00_current[t] * data_[2028+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2236+t] = C00_[t] * data_[2223+t] + B10_current[t] * data_[2210+t] + cB00_current[t] * data_[2041+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2249+t] = C00_[t] * data_[2236+t] + B10_current[t] * data_[2223+t] + cB00_current[t] * data_[2054+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2262+t] = C00_[t] * data_[2249+t] + B10_current[t] * data_[2236+t] + cB00_current[t] * data_[2067+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2275+t] = C00_[t] * data_[2262+t] + B10_current[t] * data_[2249+t] + cB00_current[t] * data_[2080+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2288+t] = C00_[t] * data_[2275+t] + B10_current[t] * data_[2262+t] + cB00_current[t] * data_[2093+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2301+t] = C00_[t] * data_[2288+t] + B10_current[t] * data_[2275+t] + cB00_current[t] * data_[2106+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2314+t] = C00_[t] * data_[2301+t] + B10_current[t] * data_[2288+t] + cB00_current[t] * data_[2119+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2327+t] = C00_[t] * data_[2314+t] + B10_current[t] * data_[2301+t] + cB00_current[t] * data_[2132+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2340+t] = C00_[t] * data_[2327+t] + B10_current[t] * data_[2314+t] + cB00_current[t] * data_[2145+t];

  for (int t = 0; t != 13; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 13; ++t)
    data_[2353+t] = C00_[t] * data_[2340+t] + B10_current[t] * data_[2327+t] + cB00_current[t] * data_[2158+t];
}

