//
// Newint - Parallel electron correlation program.
// Filename: svrrlist.h
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

//
// stand alone class for VRR functions.
// can be called from outside through vrrfunc_call.
//
  
#ifndef __slater_svrrlist_h
#define __slater_svrrlist_h

#include <src/rysint/macros.h>
 
struct SVRRList {
  SVRRList();
  ~SVRRList();
//
// commented-out lines are hand-written, which is implimented in vrr_template.cc 
//
  static void _svrr_0000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_0010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_0020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_0030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_0040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_0050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_0060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_0070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_0080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_0090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_00a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_00b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_00c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _svrr_1000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_1010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_1020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_1030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_1040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_1050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_1060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_1070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_1080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_1090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_10a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_10b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_10c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _svrr_2000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_2010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_2020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_2030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_2040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_2050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_2060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_2070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_2080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_2090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_20a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_20b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_20c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _svrr_3000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_3010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_3020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_3030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_3040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_3050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_3060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_3070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_3080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_3090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_30a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_30b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_30c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _svrr_4000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_4010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_4020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_4030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_4040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_4050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_4060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_4070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_4080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_4090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_40a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_40b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_40c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _svrr_5000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_5010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_5020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_5030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_5040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_5050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_5060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_5070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_5080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_5090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_50a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_50b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_50c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _svrr_6000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_6010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_6020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_6030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_6040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_6050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_6060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_6070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_6080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_6090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_60a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_60b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_60c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _svrr_7000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_7010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_7020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_7030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_7040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_7050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_7060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_7070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_7080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_7090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_70a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_70b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_70c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _svrr_8000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_8010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_8020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_8030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_8040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_8050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_8060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_8070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_8080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_8090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_80a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_80b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_80c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _svrr_9000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_9010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_9020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_9030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_9040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_9050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_9060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_9070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_9080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_9090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_90a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_90b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_90c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _svrr_a000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_a010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_a020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_a030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_a040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_a050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_a060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_a070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_a080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_a090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_a0a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_a0b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_a0c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _svrr_b000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_b010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_b020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_b030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_b040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_b050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_b060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_b070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_b080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_b090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_b0a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_b0b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_b0c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _svrr_c000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_c010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_c020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_c030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_c040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_c050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_c060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_c070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_c080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_c090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_c0a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_c0b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _svrr_c0c0(double*, const double*, const double*, const double*, const double*, const double*);

  void svrrfunc_call(const unsigned int i, double* a, const double* b, const double* c, const double* d, const double* e, const double* f) {
    return (svrrfunc[i])(a, b, c, d, e, f);
  };

  void (*svrrfunc[ANG_VRR_END * ANG_VRR_END])(double*, const double*, const double*, const double*, const double*, const double*);
};

#endif
