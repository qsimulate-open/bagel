//
// Newint - Parallel electron correlation program.
// Filename: _svrr_a070.cc
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


#include "svrrlist.h"

// returns double array of length 880
void SVRRList::_svrr_a070(double* data_, const double* C00_, const double* D00_, const double* B00_, const double* B01_, const double* B10_) {
  for (int t = 0; t != 10; ++t)
    data_[0+t] = 1.0;

  for (int t = 0; t != 10; ++t)
    data_[10+t] = C00_[t];

  double B10_current[10];
  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[20+t] = C00_[t] * data_[10+t] + B10_current[t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[30+t] = C00_[t] * data_[20+t] + B10_current[t] * data_[10+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[40+t] = C00_[t] * data_[30+t] + B10_current[t] * data_[20+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[50+t] = C00_[t] * data_[40+t] + B10_current[t] * data_[30+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[60+t] = C00_[t] * data_[50+t] + B10_current[t] * data_[40+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[70+t] = C00_[t] * data_[60+t] + B10_current[t] * data_[50+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[80+t] = C00_[t] * data_[70+t] + B10_current[t] * data_[60+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[90+t] = C00_[t] * data_[80+t] + B10_current[t] * data_[70+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[100+t] = C00_[t] * data_[90+t] + B10_current[t] * data_[80+t];

  for (int t = 0; t != 10; ++t)
    data_[110+t] = D00_[t];

  double cB00_current[10];
  for (int t = 0; t != 10; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[120+t] = C00_[t] * data_[110+t] + cB00_current[t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[130+t] = C00_[t] * data_[120+t] + B10_current[t] * data_[110+t] + cB00_current[t] * data_[10+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[140+t] = C00_[t] * data_[130+t] + B10_current[t] * data_[120+t] + cB00_current[t] * data_[20+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[150+t] = C00_[t] * data_[140+t] + B10_current[t] * data_[130+t] + cB00_current[t] * data_[30+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[160+t] = C00_[t] * data_[150+t] + B10_current[t] * data_[140+t] + cB00_current[t] * data_[40+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[170+t] = C00_[t] * data_[160+t] + B10_current[t] * data_[150+t] + cB00_current[t] * data_[50+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[180+t] = C00_[t] * data_[170+t] + B10_current[t] * data_[160+t] + cB00_current[t] * data_[60+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[190+t] = C00_[t] * data_[180+t] + B10_current[t] * data_[170+t] + cB00_current[t] * data_[70+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[200+t] = C00_[t] * data_[190+t] + B10_current[t] * data_[180+t] + cB00_current[t] * data_[80+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[210+t] = C00_[t] * data_[200+t] + B10_current[t] * data_[190+t] + cB00_current[t] * data_[90+t];

  double B01_current[10];
  for (int t = 0; t != 10; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[220+t] = D00_[t] * data_[110+t] + B01_current[t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[230+t] = C00_[t] * data_[220+t] + cB00_current[t] * data_[110+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[240+t] = C00_[t] * data_[230+t] + B10_current[t] * data_[220+t] + cB00_current[t] * data_[120+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[250+t] = C00_[t] * data_[240+t] + B10_current[t] * data_[230+t] + cB00_current[t] * data_[130+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[260+t] = C00_[t] * data_[250+t] + B10_current[t] * data_[240+t] + cB00_current[t] * data_[140+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[270+t] = C00_[t] * data_[260+t] + B10_current[t] * data_[250+t] + cB00_current[t] * data_[150+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[280+t] = C00_[t] * data_[270+t] + B10_current[t] * data_[260+t] + cB00_current[t] * data_[160+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[290+t] = C00_[t] * data_[280+t] + B10_current[t] * data_[270+t] + cB00_current[t] * data_[170+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[300+t] = C00_[t] * data_[290+t] + B10_current[t] * data_[280+t] + cB00_current[t] * data_[180+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[310+t] = C00_[t] * data_[300+t] + B10_current[t] * data_[290+t] + cB00_current[t] * data_[190+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[320+t] = C00_[t] * data_[310+t] + B10_current[t] * data_[300+t] + cB00_current[t] * data_[200+t];

  for (int t = 0; t != 10; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[330+t] = D00_[t] * data_[220+t] + B01_current[t] * data_[110+t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[340+t] = C00_[t] * data_[330+t] + cB00_current[t] * data_[220+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[350+t] = C00_[t] * data_[340+t] + B10_current[t] * data_[330+t] + cB00_current[t] * data_[230+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[360+t] = C00_[t] * data_[350+t] + B10_current[t] * data_[340+t] + cB00_current[t] * data_[240+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[370+t] = C00_[t] * data_[360+t] + B10_current[t] * data_[350+t] + cB00_current[t] * data_[250+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[380+t] = C00_[t] * data_[370+t] + B10_current[t] * data_[360+t] + cB00_current[t] * data_[260+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[390+t] = C00_[t] * data_[380+t] + B10_current[t] * data_[370+t] + cB00_current[t] * data_[270+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[400+t] = C00_[t] * data_[390+t] + B10_current[t] * data_[380+t] + cB00_current[t] * data_[280+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[410+t] = C00_[t] * data_[400+t] + B10_current[t] * data_[390+t] + cB00_current[t] * data_[290+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[420+t] = C00_[t] * data_[410+t] + B10_current[t] * data_[400+t] + cB00_current[t] * data_[300+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[430+t] = C00_[t] * data_[420+t] + B10_current[t] * data_[410+t] + cB00_current[t] * data_[310+t];

  for (int t = 0; t != 10; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[440+t] = D00_[t] * data_[330+t] + B01_current[t] * data_[220+t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[450+t] = C00_[t] * data_[440+t] + cB00_current[t] * data_[330+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[460+t] = C00_[t] * data_[450+t] + B10_current[t] * data_[440+t] + cB00_current[t] * data_[340+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[470+t] = C00_[t] * data_[460+t] + B10_current[t] * data_[450+t] + cB00_current[t] * data_[350+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[480+t] = C00_[t] * data_[470+t] + B10_current[t] * data_[460+t] + cB00_current[t] * data_[360+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[490+t] = C00_[t] * data_[480+t] + B10_current[t] * data_[470+t] + cB00_current[t] * data_[370+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[500+t] = C00_[t] * data_[490+t] + B10_current[t] * data_[480+t] + cB00_current[t] * data_[380+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[510+t] = C00_[t] * data_[500+t] + B10_current[t] * data_[490+t] + cB00_current[t] * data_[390+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[520+t] = C00_[t] * data_[510+t] + B10_current[t] * data_[500+t] + cB00_current[t] * data_[400+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[530+t] = C00_[t] * data_[520+t] + B10_current[t] * data_[510+t] + cB00_current[t] * data_[410+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[540+t] = C00_[t] * data_[530+t] + B10_current[t] * data_[520+t] + cB00_current[t] * data_[420+t];

  for (int t = 0; t != 10; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[550+t] = D00_[t] * data_[440+t] + B01_current[t] * data_[330+t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[560+t] = C00_[t] * data_[550+t] + cB00_current[t] * data_[440+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[570+t] = C00_[t] * data_[560+t] + B10_current[t] * data_[550+t] + cB00_current[t] * data_[450+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[580+t] = C00_[t] * data_[570+t] + B10_current[t] * data_[560+t] + cB00_current[t] * data_[460+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[590+t] = C00_[t] * data_[580+t] + B10_current[t] * data_[570+t] + cB00_current[t] * data_[470+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[600+t] = C00_[t] * data_[590+t] + B10_current[t] * data_[580+t] + cB00_current[t] * data_[480+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[610+t] = C00_[t] * data_[600+t] + B10_current[t] * data_[590+t] + cB00_current[t] * data_[490+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[620+t] = C00_[t] * data_[610+t] + B10_current[t] * data_[600+t] + cB00_current[t] * data_[500+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[630+t] = C00_[t] * data_[620+t] + B10_current[t] * data_[610+t] + cB00_current[t] * data_[510+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[640+t] = C00_[t] * data_[630+t] + B10_current[t] * data_[620+t] + cB00_current[t] * data_[520+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[650+t] = C00_[t] * data_[640+t] + B10_current[t] * data_[630+t] + cB00_current[t] * data_[530+t];

  for (int t = 0; t != 10; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[660+t] = D00_[t] * data_[550+t] + B01_current[t] * data_[440+t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[670+t] = C00_[t] * data_[660+t] + cB00_current[t] * data_[550+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[680+t] = C00_[t] * data_[670+t] + B10_current[t] * data_[660+t] + cB00_current[t] * data_[560+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[690+t] = C00_[t] * data_[680+t] + B10_current[t] * data_[670+t] + cB00_current[t] * data_[570+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[700+t] = C00_[t] * data_[690+t] + B10_current[t] * data_[680+t] + cB00_current[t] * data_[580+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[710+t] = C00_[t] * data_[700+t] + B10_current[t] * data_[690+t] + cB00_current[t] * data_[590+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[720+t] = C00_[t] * data_[710+t] + B10_current[t] * data_[700+t] + cB00_current[t] * data_[600+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[730+t] = C00_[t] * data_[720+t] + B10_current[t] * data_[710+t] + cB00_current[t] * data_[610+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[740+t] = C00_[t] * data_[730+t] + B10_current[t] * data_[720+t] + cB00_current[t] * data_[620+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[750+t] = C00_[t] * data_[740+t] + B10_current[t] * data_[730+t] + cB00_current[t] * data_[630+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[760+t] = C00_[t] * data_[750+t] + B10_current[t] * data_[740+t] + cB00_current[t] * data_[640+t];

  for (int t = 0; t != 10; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[770+t] = D00_[t] * data_[660+t] + B01_current[t] * data_[550+t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[780+t] = C00_[t] * data_[770+t] + cB00_current[t] * data_[660+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[790+t] = C00_[t] * data_[780+t] + B10_current[t] * data_[770+t] + cB00_current[t] * data_[670+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[800+t] = C00_[t] * data_[790+t] + B10_current[t] * data_[780+t] + cB00_current[t] * data_[680+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[810+t] = C00_[t] * data_[800+t] + B10_current[t] * data_[790+t] + cB00_current[t] * data_[690+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[820+t] = C00_[t] * data_[810+t] + B10_current[t] * data_[800+t] + cB00_current[t] * data_[700+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[830+t] = C00_[t] * data_[820+t] + B10_current[t] * data_[810+t] + cB00_current[t] * data_[710+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[840+t] = C00_[t] * data_[830+t] + B10_current[t] * data_[820+t] + cB00_current[t] * data_[720+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[850+t] = C00_[t] * data_[840+t] + B10_current[t] * data_[830+t] + cB00_current[t] * data_[730+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[860+t] = C00_[t] * data_[850+t] + B10_current[t] * data_[840+t] + cB00_current[t] * data_[740+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[870+t] = C00_[t] * data_[860+t] + B10_current[t] * data_[850+t] + cB00_current[t] * data_[750+t];
}

