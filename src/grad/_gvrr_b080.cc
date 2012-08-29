//
// BAGEL - Parallel electron correlation program.
// Filename: _gvrr_b080.cc
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

// returns double array of length 1080
void GVRRList::_gvrr_b080(double* data_, const double* C00, const double* D00, const double* B00, const double* B01, const double* B10) {
#ifdef __GNUC__
  const double C00_[10]__attribute__((aligned(32))) = {C00[0], C00[1], C00[2], C00[3], C00[4], C00[5], C00[6], C00[7], C00[8], C00[9]};
  const double D00_[10]__attribute__((aligned(32))) = {D00[0], D00[1], D00[2], D00[3], D00[4], D00[5], D00[6], D00[7], D00[8], D00[9]};
  const double B00_[10]__attribute__((aligned(32))) = {B00[0], B00[1], B00[2], B00[3], B00[4], B00[5], B00[6], B00[7], B00[8], B00[9]};
  const double B01_[10]__attribute__((aligned(32))) = {B01[0], B01[1], B01[2], B01[3], B01[4], B01[5], B01[6], B01[7], B01[8], B01[9]};
  const double B10_[10]__attribute__((aligned(32))) = {B10[0], B10[1], B10[2], B10[3], B10[4], B10[5], B10[6], B10[7], B10[8], B10[9]};
#else
  const double* C00_ = C00;
  const double* D00_ = D00;
  const double* B00_ = B00;
  const double* B01_ = B01;
  const double* B10_ = B10;
#endif

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
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[110+t] = C00_[t] * data_[100+t] + B10_current[t] * data_[90+t];

  for (int t = 0; t != 10; ++t)
    data_[120+t] = D00_[t];

  double cB00_current[10];
  for (int t = 0; t != 10; ++t)
    cB00_current[t] = B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[130+t] = C00_[t] * data_[120+t] + cB00_current[t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[140+t] = C00_[t] * data_[130+t] + B10_current[t] * data_[120+t] + cB00_current[t] * data_[10+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[150+t] = C00_[t] * data_[140+t] + B10_current[t] * data_[130+t] + cB00_current[t] * data_[20+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[160+t] = C00_[t] * data_[150+t] + B10_current[t] * data_[140+t] + cB00_current[t] * data_[30+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[170+t] = C00_[t] * data_[160+t] + B10_current[t] * data_[150+t] + cB00_current[t] * data_[40+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[180+t] = C00_[t] * data_[170+t] + B10_current[t] * data_[160+t] + cB00_current[t] * data_[50+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[190+t] = C00_[t] * data_[180+t] + B10_current[t] * data_[170+t] + cB00_current[t] * data_[60+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[200+t] = C00_[t] * data_[190+t] + B10_current[t] * data_[180+t] + cB00_current[t] * data_[70+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[210+t] = C00_[t] * data_[200+t] + B10_current[t] * data_[190+t] + cB00_current[t] * data_[80+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[220+t] = C00_[t] * data_[210+t] + B10_current[t] * data_[200+t] + cB00_current[t] * data_[90+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[230+t] = C00_[t] * data_[220+t] + B10_current[t] * data_[210+t] + cB00_current[t] * data_[100+t];

  double B01_current[10];
  for (int t = 0; t != 10; ++t)
    B01_current[t] = B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[240+t] = D00_[t] * data_[120+t] + B01_current[t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[250+t] = C00_[t] * data_[240+t] + cB00_current[t] * data_[120+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[260+t] = C00_[t] * data_[250+t] + B10_current[t] * data_[240+t] + cB00_current[t] * data_[130+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[270+t] = C00_[t] * data_[260+t] + B10_current[t] * data_[250+t] + cB00_current[t] * data_[140+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[280+t] = C00_[t] * data_[270+t] + B10_current[t] * data_[260+t] + cB00_current[t] * data_[150+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[290+t] = C00_[t] * data_[280+t] + B10_current[t] * data_[270+t] + cB00_current[t] * data_[160+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[300+t] = C00_[t] * data_[290+t] + B10_current[t] * data_[280+t] + cB00_current[t] * data_[170+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[310+t] = C00_[t] * data_[300+t] + B10_current[t] * data_[290+t] + cB00_current[t] * data_[180+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[320+t] = C00_[t] * data_[310+t] + B10_current[t] * data_[300+t] + cB00_current[t] * data_[190+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[330+t] = C00_[t] * data_[320+t] + B10_current[t] * data_[310+t] + cB00_current[t] * data_[200+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[340+t] = C00_[t] * data_[330+t] + B10_current[t] * data_[320+t] + cB00_current[t] * data_[210+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[350+t] = C00_[t] * data_[340+t] + B10_current[t] * data_[330+t] + cB00_current[t] * data_[220+t];

  for (int t = 0; t != 10; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[360+t] = D00_[t] * data_[240+t] + B01_current[t] * data_[120+t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[370+t] = C00_[t] * data_[360+t] + cB00_current[t] * data_[240+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[380+t] = C00_[t] * data_[370+t] + B10_current[t] * data_[360+t] + cB00_current[t] * data_[250+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[390+t] = C00_[t] * data_[380+t] + B10_current[t] * data_[370+t] + cB00_current[t] * data_[260+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[400+t] = C00_[t] * data_[390+t] + B10_current[t] * data_[380+t] + cB00_current[t] * data_[270+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[410+t] = C00_[t] * data_[400+t] + B10_current[t] * data_[390+t] + cB00_current[t] * data_[280+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[420+t] = C00_[t] * data_[410+t] + B10_current[t] * data_[400+t] + cB00_current[t] * data_[290+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[430+t] = C00_[t] * data_[420+t] + B10_current[t] * data_[410+t] + cB00_current[t] * data_[300+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[440+t] = C00_[t] * data_[430+t] + B10_current[t] * data_[420+t] + cB00_current[t] * data_[310+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[450+t] = C00_[t] * data_[440+t] + B10_current[t] * data_[430+t] + cB00_current[t] * data_[320+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[460+t] = C00_[t] * data_[450+t] + B10_current[t] * data_[440+t] + cB00_current[t] * data_[330+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[470+t] = C00_[t] * data_[460+t] + B10_current[t] * data_[450+t] + cB00_current[t] * data_[340+t];

  for (int t = 0; t != 10; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[480+t] = D00_[t] * data_[360+t] + B01_current[t] * data_[240+t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[490+t] = C00_[t] * data_[480+t] + cB00_current[t] * data_[360+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[500+t] = C00_[t] * data_[490+t] + B10_current[t] * data_[480+t] + cB00_current[t] * data_[370+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[510+t] = C00_[t] * data_[500+t] + B10_current[t] * data_[490+t] + cB00_current[t] * data_[380+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[520+t] = C00_[t] * data_[510+t] + B10_current[t] * data_[500+t] + cB00_current[t] * data_[390+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[530+t] = C00_[t] * data_[520+t] + B10_current[t] * data_[510+t] + cB00_current[t] * data_[400+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[540+t] = C00_[t] * data_[530+t] + B10_current[t] * data_[520+t] + cB00_current[t] * data_[410+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[550+t] = C00_[t] * data_[540+t] + B10_current[t] * data_[530+t] + cB00_current[t] * data_[420+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[560+t] = C00_[t] * data_[550+t] + B10_current[t] * data_[540+t] + cB00_current[t] * data_[430+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[570+t] = C00_[t] * data_[560+t] + B10_current[t] * data_[550+t] + cB00_current[t] * data_[440+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[580+t] = C00_[t] * data_[570+t] + B10_current[t] * data_[560+t] + cB00_current[t] * data_[450+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[590+t] = C00_[t] * data_[580+t] + B10_current[t] * data_[570+t] + cB00_current[t] * data_[460+t];

  for (int t = 0; t != 10; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[600+t] = D00_[t] * data_[480+t] + B01_current[t] * data_[360+t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[610+t] = C00_[t] * data_[600+t] + cB00_current[t] * data_[480+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[620+t] = C00_[t] * data_[610+t] + B10_current[t] * data_[600+t] + cB00_current[t] * data_[490+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[630+t] = C00_[t] * data_[620+t] + B10_current[t] * data_[610+t] + cB00_current[t] * data_[500+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[640+t] = C00_[t] * data_[630+t] + B10_current[t] * data_[620+t] + cB00_current[t] * data_[510+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[650+t] = C00_[t] * data_[640+t] + B10_current[t] * data_[630+t] + cB00_current[t] * data_[520+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[660+t] = C00_[t] * data_[650+t] + B10_current[t] * data_[640+t] + cB00_current[t] * data_[530+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[670+t] = C00_[t] * data_[660+t] + B10_current[t] * data_[650+t] + cB00_current[t] * data_[540+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[680+t] = C00_[t] * data_[670+t] + B10_current[t] * data_[660+t] + cB00_current[t] * data_[550+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[690+t] = C00_[t] * data_[680+t] + B10_current[t] * data_[670+t] + cB00_current[t] * data_[560+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[700+t] = C00_[t] * data_[690+t] + B10_current[t] * data_[680+t] + cB00_current[t] * data_[570+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[710+t] = C00_[t] * data_[700+t] + B10_current[t] * data_[690+t] + cB00_current[t] * data_[580+t];

  for (int t = 0; t != 10; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[720+t] = D00_[t] * data_[600+t] + B01_current[t] * data_[480+t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[730+t] = C00_[t] * data_[720+t] + cB00_current[t] * data_[600+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[740+t] = C00_[t] * data_[730+t] + B10_current[t] * data_[720+t] + cB00_current[t] * data_[610+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[750+t] = C00_[t] * data_[740+t] + B10_current[t] * data_[730+t] + cB00_current[t] * data_[620+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[760+t] = C00_[t] * data_[750+t] + B10_current[t] * data_[740+t] + cB00_current[t] * data_[630+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[770+t] = C00_[t] * data_[760+t] + B10_current[t] * data_[750+t] + cB00_current[t] * data_[640+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[780+t] = C00_[t] * data_[770+t] + B10_current[t] * data_[760+t] + cB00_current[t] * data_[650+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[790+t] = C00_[t] * data_[780+t] + B10_current[t] * data_[770+t] + cB00_current[t] * data_[660+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[800+t] = C00_[t] * data_[790+t] + B10_current[t] * data_[780+t] + cB00_current[t] * data_[670+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[810+t] = C00_[t] * data_[800+t] + B10_current[t] * data_[790+t] + cB00_current[t] * data_[680+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[820+t] = C00_[t] * data_[810+t] + B10_current[t] * data_[800+t] + cB00_current[t] * data_[690+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[830+t] = C00_[t] * data_[820+t] + B10_current[t] * data_[810+t] + cB00_current[t] * data_[700+t];

  for (int t = 0; t != 10; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[840+t] = D00_[t] * data_[720+t] + B01_current[t] * data_[600+t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[850+t] = C00_[t] * data_[840+t] + cB00_current[t] * data_[720+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[860+t] = C00_[t] * data_[850+t] + B10_current[t] * data_[840+t] + cB00_current[t] * data_[730+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[870+t] = C00_[t] * data_[860+t] + B10_current[t] * data_[850+t] + cB00_current[t] * data_[740+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[880+t] = C00_[t] * data_[870+t] + B10_current[t] * data_[860+t] + cB00_current[t] * data_[750+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[890+t] = C00_[t] * data_[880+t] + B10_current[t] * data_[870+t] + cB00_current[t] * data_[760+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[900+t] = C00_[t] * data_[890+t] + B10_current[t] * data_[880+t] + cB00_current[t] * data_[770+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[910+t] = C00_[t] * data_[900+t] + B10_current[t] * data_[890+t] + cB00_current[t] * data_[780+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[920+t] = C00_[t] * data_[910+t] + B10_current[t] * data_[900+t] + cB00_current[t] * data_[790+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[930+t] = C00_[t] * data_[920+t] + B10_current[t] * data_[910+t] + cB00_current[t] * data_[800+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[940+t] = C00_[t] * data_[930+t] + B10_current[t] * data_[920+t] + cB00_current[t] * data_[810+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[950+t] = C00_[t] * data_[940+t] + B10_current[t] * data_[930+t] + cB00_current[t] * data_[820+t];

  for (int t = 0; t != 10; ++t)
    B01_current[t] += B01_[t];

  for (int t = 0; t != 10; ++t)
    data_[960+t] = D00_[t] * data_[840+t] + B01_current[t] * data_[720+t];

  for (int t = 0; t != 10; ++t)
    cB00_current[t] += B00_[t];

  for (int t = 0; t != 10; ++t)
    data_[970+t] = C00_[t] * data_[960+t] + cB00_current[t] * data_[840+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] = B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[980+t] = C00_[t] * data_[970+t] + B10_current[t] * data_[960+t] + cB00_current[t] * data_[850+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[990+t] = C00_[t] * data_[980+t] + B10_current[t] * data_[970+t] + cB00_current[t] * data_[860+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[1000+t] = C00_[t] * data_[990+t] + B10_current[t] * data_[980+t] + cB00_current[t] * data_[870+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[1010+t] = C00_[t] * data_[1000+t] + B10_current[t] * data_[990+t] + cB00_current[t] * data_[880+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[1020+t] = C00_[t] * data_[1010+t] + B10_current[t] * data_[1000+t] + cB00_current[t] * data_[890+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[1030+t] = C00_[t] * data_[1020+t] + B10_current[t] * data_[1010+t] + cB00_current[t] * data_[900+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[1040+t] = C00_[t] * data_[1030+t] + B10_current[t] * data_[1020+t] + cB00_current[t] * data_[910+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[1050+t] = C00_[t] * data_[1040+t] + B10_current[t] * data_[1030+t] + cB00_current[t] * data_[920+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[1060+t] = C00_[t] * data_[1050+t] + B10_current[t] * data_[1040+t] + cB00_current[t] * data_[930+t];

  for (int t = 0; t != 10; ++t)
    B10_current[t] += B10_[t];

  for (int t = 0; t != 10; ++t)
    data_[1070+t] = C00_[t] * data_[1060+t] + B10_current[t] * data_[1050+t] + cB00_current[t] * data_[940+t];
}

