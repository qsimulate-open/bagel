//
// BAGEL - Parallel electron correlation program.
// Filename: _svrr_c0a0.cc
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

// returns double array of length 1716
void SVRRList::_svrr_c0a0(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  for (int t = 0; t != 12; ++t)
    data_[0+t] = 1.0;

  for (int t = 0; t != 12; ++t)
    data_[12+t] = C00_[t];

  double B10_current[12];
  for (int t = 0; t != 12; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[24+t] = C00_[t] * data_[12+t] + B10_current[t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[36+t] = C00_[t] * data_[24+t] + B10_current[t] * data_[12+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[48+t] = C00_[t] * data_[36+t] + B10_current[t] * data_[24+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[60+t] = C00_[t] * data_[48+t] + B10_current[t] * data_[36+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[72+t] = C00_[t] * data_[60+t] + B10_current[t] * data_[48+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[84+t] = C00_[t] * data_[72+t] + B10_current[t] * data_[60+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[96+t] = C00_[t] * data_[84+t] + B10_current[t] * data_[72+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[108+t] = C00_[t] * data_[96+t] + B10_current[t] * data_[84+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[120+t] = C00_[t] * data_[108+t] + B10_current[t] * data_[96+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[132+t] = C00_[t] * data_[120+t] + B10_current[t] * data_[108+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[144+t] = C00_[t] * data_[132+t] + B10_current[t] * data_[120+t];

  for (int t = 0; t != 12; ++t)
    data_[156+t] = D00_[t];

  double cB00_current[12];
  for (int t = 0; t != 12; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 12; ++t)
    data_[168+t] = C00_[t] * data_[156+t] + cB00_current[t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[180+t] = C00_[t] * data_[168+t] + B10_current[t] * data_[156+t] + cB00_current[t] * data_[12+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[192+t] = C00_[t] * data_[180+t] + B10_current[t] * data_[168+t] + cB00_current[t] * data_[24+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[204+t] = C00_[t] * data_[192+t] + B10_current[t] * data_[180+t] + cB00_current[t] * data_[36+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[216+t] = C00_[t] * data_[204+t] + B10_current[t] * data_[192+t] + cB00_current[t] * data_[48+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[228+t] = C00_[t] * data_[216+t] + B10_current[t] * data_[204+t] + cB00_current[t] * data_[60+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[240+t] = C00_[t] * data_[228+t] + B10_current[t] * data_[216+t] + cB00_current[t] * data_[72+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[252+t] = C00_[t] * data_[240+t] + B10_current[t] * data_[228+t] + cB00_current[t] * data_[84+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[264+t] = C00_[t] * data_[252+t] + B10_current[t] * data_[240+t] + cB00_current[t] * data_[96+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[276+t] = C00_[t] * data_[264+t] + B10_current[t] * data_[252+t] + cB00_current[t] * data_[108+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[288+t] = C00_[t] * data_[276+t] + B10_current[t] * data_[264+t] + cB00_current[t] * data_[120+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[300+t] = C00_[t] * data_[288+t] + B10_current[t] * data_[276+t] + cB00_current[t] * data_[132+t];

  double B01_current[12];
  for (int t = 0; t != 12; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 12; ++t)
    data_[312+t] = D00_[t] * data_[156+t] + B01_current[t];

  for (int t = 0; t != 12; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 12; ++t)
    data_[324+t] = C00_[t] * data_[312+t] + cB00_current[t] * data_[156+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[336+t] = C00_[t] * data_[324+t] + B10_current[t] * data_[312+t] + cB00_current[t] * data_[168+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[348+t] = C00_[t] * data_[336+t] + B10_current[t] * data_[324+t] + cB00_current[t] * data_[180+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[360+t] = C00_[t] * data_[348+t] + B10_current[t] * data_[336+t] + cB00_current[t] * data_[192+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[372+t] = C00_[t] * data_[360+t] + B10_current[t] * data_[348+t] + cB00_current[t] * data_[204+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[384+t] = C00_[t] * data_[372+t] + B10_current[t] * data_[360+t] + cB00_current[t] * data_[216+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[396+t] = C00_[t] * data_[384+t] + B10_current[t] * data_[372+t] + cB00_current[t] * data_[228+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[408+t] = C00_[t] * data_[396+t] + B10_current[t] * data_[384+t] + cB00_current[t] * data_[240+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[420+t] = C00_[t] * data_[408+t] + B10_current[t] * data_[396+t] + cB00_current[t] * data_[252+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[432+t] = C00_[t] * data_[420+t] + B10_current[t] * data_[408+t] + cB00_current[t] * data_[264+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[444+t] = C00_[t] * data_[432+t] + B10_current[t] * data_[420+t] + cB00_current[t] * data_[276+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[456+t] = C00_[t] * data_[444+t] + B10_current[t] * data_[432+t] + cB00_current[t] * data_[288+t];

  for (int t = 0; t != 12; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 12; ++t)
    data_[468+t] = D00_[t] * data_[312+t] + B01_current[t] * data_[156+t];

  for (int t = 0; t != 12; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 12; ++t)
    data_[480+t] = C00_[t] * data_[468+t] + cB00_current[t] * data_[312+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[492+t] = C00_[t] * data_[480+t] + B10_current[t] * data_[468+t] + cB00_current[t] * data_[324+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[504+t] = C00_[t] * data_[492+t] + B10_current[t] * data_[480+t] + cB00_current[t] * data_[336+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[516+t] = C00_[t] * data_[504+t] + B10_current[t] * data_[492+t] + cB00_current[t] * data_[348+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[528+t] = C00_[t] * data_[516+t] + B10_current[t] * data_[504+t] + cB00_current[t] * data_[360+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[540+t] = C00_[t] * data_[528+t] + B10_current[t] * data_[516+t] + cB00_current[t] * data_[372+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[552+t] = C00_[t] * data_[540+t] + B10_current[t] * data_[528+t] + cB00_current[t] * data_[384+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[564+t] = C00_[t] * data_[552+t] + B10_current[t] * data_[540+t] + cB00_current[t] * data_[396+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[576+t] = C00_[t] * data_[564+t] + B10_current[t] * data_[552+t] + cB00_current[t] * data_[408+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[588+t] = C00_[t] * data_[576+t] + B10_current[t] * data_[564+t] + cB00_current[t] * data_[420+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[600+t] = C00_[t] * data_[588+t] + B10_current[t] * data_[576+t] + cB00_current[t] * data_[432+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[612+t] = C00_[t] * data_[600+t] + B10_current[t] * data_[588+t] + cB00_current[t] * data_[444+t];

  for (int t = 0; t != 12; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 12; ++t)
    data_[624+t] = D00_[t] * data_[468+t] + B01_current[t] * data_[312+t];

  for (int t = 0; t != 12; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 12; ++t)
    data_[636+t] = C00_[t] * data_[624+t] + cB00_current[t] * data_[468+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[648+t] = C00_[t] * data_[636+t] + B10_current[t] * data_[624+t] + cB00_current[t] * data_[480+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[660+t] = C00_[t] * data_[648+t] + B10_current[t] * data_[636+t] + cB00_current[t] * data_[492+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[672+t] = C00_[t] * data_[660+t] + B10_current[t] * data_[648+t] + cB00_current[t] * data_[504+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[684+t] = C00_[t] * data_[672+t] + B10_current[t] * data_[660+t] + cB00_current[t] * data_[516+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[696+t] = C00_[t] * data_[684+t] + B10_current[t] * data_[672+t] + cB00_current[t] * data_[528+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[708+t] = C00_[t] * data_[696+t] + B10_current[t] * data_[684+t] + cB00_current[t] * data_[540+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[720+t] = C00_[t] * data_[708+t] + B10_current[t] * data_[696+t] + cB00_current[t] * data_[552+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[732+t] = C00_[t] * data_[720+t] + B10_current[t] * data_[708+t] + cB00_current[t] * data_[564+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[744+t] = C00_[t] * data_[732+t] + B10_current[t] * data_[720+t] + cB00_current[t] * data_[576+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[756+t] = C00_[t] * data_[744+t] + B10_current[t] * data_[732+t] + cB00_current[t] * data_[588+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[768+t] = C00_[t] * data_[756+t] + B10_current[t] * data_[744+t] + cB00_current[t] * data_[600+t];

  for (int t = 0; t != 12; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 12; ++t)
    data_[780+t] = D00_[t] * data_[624+t] + B01_current[t] * data_[468+t];

  for (int t = 0; t != 12; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 12; ++t)
    data_[792+t] = C00_[t] * data_[780+t] + cB00_current[t] * data_[624+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[804+t] = C00_[t] * data_[792+t] + B10_current[t] * data_[780+t] + cB00_current[t] * data_[636+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[816+t] = C00_[t] * data_[804+t] + B10_current[t] * data_[792+t] + cB00_current[t] * data_[648+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[828+t] = C00_[t] * data_[816+t] + B10_current[t] * data_[804+t] + cB00_current[t] * data_[660+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[840+t] = C00_[t] * data_[828+t] + B10_current[t] * data_[816+t] + cB00_current[t] * data_[672+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[852+t] = C00_[t] * data_[840+t] + B10_current[t] * data_[828+t] + cB00_current[t] * data_[684+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[864+t] = C00_[t] * data_[852+t] + B10_current[t] * data_[840+t] + cB00_current[t] * data_[696+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[876+t] = C00_[t] * data_[864+t] + B10_current[t] * data_[852+t] + cB00_current[t] * data_[708+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[888+t] = C00_[t] * data_[876+t] + B10_current[t] * data_[864+t] + cB00_current[t] * data_[720+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[900+t] = C00_[t] * data_[888+t] + B10_current[t] * data_[876+t] + cB00_current[t] * data_[732+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[912+t] = C00_[t] * data_[900+t] + B10_current[t] * data_[888+t] + cB00_current[t] * data_[744+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[924+t] = C00_[t] * data_[912+t] + B10_current[t] * data_[900+t] + cB00_current[t] * data_[756+t];

  for (int t = 0; t != 12; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 12; ++t)
    data_[936+t] = D00_[t] * data_[780+t] + B01_current[t] * data_[624+t];

  for (int t = 0; t != 12; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 12; ++t)
    data_[948+t] = C00_[t] * data_[936+t] + cB00_current[t] * data_[780+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[960+t] = C00_[t] * data_[948+t] + B10_current[t] * data_[936+t] + cB00_current[t] * data_[792+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[972+t] = C00_[t] * data_[960+t] + B10_current[t] * data_[948+t] + cB00_current[t] * data_[804+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[984+t] = C00_[t] * data_[972+t] + B10_current[t] * data_[960+t] + cB00_current[t] * data_[816+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[996+t] = C00_[t] * data_[984+t] + B10_current[t] * data_[972+t] + cB00_current[t] * data_[828+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1008+t] = C00_[t] * data_[996+t] + B10_current[t] * data_[984+t] + cB00_current[t] * data_[840+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1020+t] = C00_[t] * data_[1008+t] + B10_current[t] * data_[996+t] + cB00_current[t] * data_[852+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1032+t] = C00_[t] * data_[1020+t] + B10_current[t] * data_[1008+t] + cB00_current[t] * data_[864+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1044+t] = C00_[t] * data_[1032+t] + B10_current[t] * data_[1020+t] + cB00_current[t] * data_[876+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1056+t] = C00_[t] * data_[1044+t] + B10_current[t] * data_[1032+t] + cB00_current[t] * data_[888+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1068+t] = C00_[t] * data_[1056+t] + B10_current[t] * data_[1044+t] + cB00_current[t] * data_[900+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1080+t] = C00_[t] * data_[1068+t] + B10_current[t] * data_[1056+t] + cB00_current[t] * data_[912+t];

  for (int t = 0; t != 12; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 12; ++t)
    data_[1092+t] = D00_[t] * data_[936+t] + B01_current[t] * data_[780+t];

  for (int t = 0; t != 12; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 12; ++t)
    data_[1104+t] = C00_[t] * data_[1092+t] + cB00_current[t] * data_[936+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1116+t] = C00_[t] * data_[1104+t] + B10_current[t] * data_[1092+t] + cB00_current[t] * data_[948+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1128+t] = C00_[t] * data_[1116+t] + B10_current[t] * data_[1104+t] + cB00_current[t] * data_[960+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1140+t] = C00_[t] * data_[1128+t] + B10_current[t] * data_[1116+t] + cB00_current[t] * data_[972+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1152+t] = C00_[t] * data_[1140+t] + B10_current[t] * data_[1128+t] + cB00_current[t] * data_[984+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1164+t] = C00_[t] * data_[1152+t] + B10_current[t] * data_[1140+t] + cB00_current[t] * data_[996+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1176+t] = C00_[t] * data_[1164+t] + B10_current[t] * data_[1152+t] + cB00_current[t] * data_[1008+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1188+t] = C00_[t] * data_[1176+t] + B10_current[t] * data_[1164+t] + cB00_current[t] * data_[1020+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1200+t] = C00_[t] * data_[1188+t] + B10_current[t] * data_[1176+t] + cB00_current[t] * data_[1032+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1212+t] = C00_[t] * data_[1200+t] + B10_current[t] * data_[1188+t] + cB00_current[t] * data_[1044+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1224+t] = C00_[t] * data_[1212+t] + B10_current[t] * data_[1200+t] + cB00_current[t] * data_[1056+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1236+t] = C00_[t] * data_[1224+t] + B10_current[t] * data_[1212+t] + cB00_current[t] * data_[1068+t];

  for (int t = 0; t != 12; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 12; ++t)
    data_[1248+t] = D00_[t] * data_[1092+t] + B01_current[t] * data_[936+t];

  for (int t = 0; t != 12; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 12; ++t)
    data_[1260+t] = C00_[t] * data_[1248+t] + cB00_current[t] * data_[1092+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1272+t] = C00_[t] * data_[1260+t] + B10_current[t] * data_[1248+t] + cB00_current[t] * data_[1104+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1284+t] = C00_[t] * data_[1272+t] + B10_current[t] * data_[1260+t] + cB00_current[t] * data_[1116+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1296+t] = C00_[t] * data_[1284+t] + B10_current[t] * data_[1272+t] + cB00_current[t] * data_[1128+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1308+t] = C00_[t] * data_[1296+t] + B10_current[t] * data_[1284+t] + cB00_current[t] * data_[1140+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1320+t] = C00_[t] * data_[1308+t] + B10_current[t] * data_[1296+t] + cB00_current[t] * data_[1152+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1332+t] = C00_[t] * data_[1320+t] + B10_current[t] * data_[1308+t] + cB00_current[t] * data_[1164+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1344+t] = C00_[t] * data_[1332+t] + B10_current[t] * data_[1320+t] + cB00_current[t] * data_[1176+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1356+t] = C00_[t] * data_[1344+t] + B10_current[t] * data_[1332+t] + cB00_current[t] * data_[1188+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1368+t] = C00_[t] * data_[1356+t] + B10_current[t] * data_[1344+t] + cB00_current[t] * data_[1200+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1380+t] = C00_[t] * data_[1368+t] + B10_current[t] * data_[1356+t] + cB00_current[t] * data_[1212+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1392+t] = C00_[t] * data_[1380+t] + B10_current[t] * data_[1368+t] + cB00_current[t] * data_[1224+t];

  for (int t = 0; t != 12; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 12; ++t)
    data_[1404+t] = D00_[t] * data_[1248+t] + B01_current[t] * data_[1092+t];

  for (int t = 0; t != 12; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 12; ++t)
    data_[1416+t] = C00_[t] * data_[1404+t] + cB00_current[t] * data_[1248+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1428+t] = C00_[t] * data_[1416+t] + B10_current[t] * data_[1404+t] + cB00_current[t] * data_[1260+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1440+t] = C00_[t] * data_[1428+t] + B10_current[t] * data_[1416+t] + cB00_current[t] * data_[1272+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1452+t] = C00_[t] * data_[1440+t] + B10_current[t] * data_[1428+t] + cB00_current[t] * data_[1284+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1464+t] = C00_[t] * data_[1452+t] + B10_current[t] * data_[1440+t] + cB00_current[t] * data_[1296+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1476+t] = C00_[t] * data_[1464+t] + B10_current[t] * data_[1452+t] + cB00_current[t] * data_[1308+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1488+t] = C00_[t] * data_[1476+t] + B10_current[t] * data_[1464+t] + cB00_current[t] * data_[1320+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1500+t] = C00_[t] * data_[1488+t] + B10_current[t] * data_[1476+t] + cB00_current[t] * data_[1332+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1512+t] = C00_[t] * data_[1500+t] + B10_current[t] * data_[1488+t] + cB00_current[t] * data_[1344+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1524+t] = C00_[t] * data_[1512+t] + B10_current[t] * data_[1500+t] + cB00_current[t] * data_[1356+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1536+t] = C00_[t] * data_[1524+t] + B10_current[t] * data_[1512+t] + cB00_current[t] * data_[1368+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1548+t] = C00_[t] * data_[1536+t] + B10_current[t] * data_[1524+t] + cB00_current[t] * data_[1380+t];

  for (int t = 0; t != 12; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 12; ++t)
    data_[1560+t] = D00_[t] * data_[1404+t] + B01_current[t] * data_[1248+t];

  for (int t = 0; t != 12; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 12; ++t)
    data_[1572+t] = C00_[t] * data_[1560+t] + cB00_current[t] * data_[1404+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1584+t] = C00_[t] * data_[1572+t] + B10_current[t] * data_[1560+t] + cB00_current[t] * data_[1416+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1596+t] = C00_[t] * data_[1584+t] + B10_current[t] * data_[1572+t] + cB00_current[t] * data_[1428+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1608+t] = C00_[t] * data_[1596+t] + B10_current[t] * data_[1584+t] + cB00_current[t] * data_[1440+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1620+t] = C00_[t] * data_[1608+t] + B10_current[t] * data_[1596+t] + cB00_current[t] * data_[1452+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1632+t] = C00_[t] * data_[1620+t] + B10_current[t] * data_[1608+t] + cB00_current[t] * data_[1464+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1644+t] = C00_[t] * data_[1632+t] + B10_current[t] * data_[1620+t] + cB00_current[t] * data_[1476+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1656+t] = C00_[t] * data_[1644+t] + B10_current[t] * data_[1632+t] + cB00_current[t] * data_[1488+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1668+t] = C00_[t] * data_[1656+t] + B10_current[t] * data_[1644+t] + cB00_current[t] * data_[1500+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1680+t] = C00_[t] * data_[1668+t] + B10_current[t] * data_[1656+t] + cB00_current[t] * data_[1512+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1692+t] = C00_[t] * data_[1680+t] + B10_current[t] * data_[1668+t] + cB00_current[t] * data_[1524+t];

  for (int t = 0; t != 12; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 12; ++t)
    data_[1704+t] = C00_[t] * data_[1692+t] + B10_current[t] * data_[1680+t] + cB00_current[t] * data_[1536+t];
}

