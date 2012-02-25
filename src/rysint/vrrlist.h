//
// Newint - Parallel electron correlation program.
// Filename: vrrlist.h
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
  
#ifndef __rysint_vrrlist_h
#define __rysint_vrrlist_h

#include <src/rysint/macros.h>
 
struct VRRList {
  VRRList();
  ~VRRList();
//
// commented-out lines are hand-written, which is implimented in vrr_template.cc 
//
//static void _vrr_0000(double*, const double*, const double*, const double*, const double*, const double*);
//static void _vrr_0010(double*, const double*, const double*, const double*, const double*, const double*);
//static void _vrr_0020(double*, const double*, const double*, const double*, const double*, const double*);
//static void _vrr_0030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_0040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_0050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_0060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_0070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_0080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_0090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_00a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_00b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_00c0(double*, const double*, const double*, const double*, const double*, const double*);

//static void _vrr_1000(double*, const double*, const double*, const double*, const double*, const double*);
//static void _vrr_1010(double*, const double*, const double*, const double*, const double*, const double*);
//static void _vrr_1020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_1030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_1040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_1050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_1060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_1070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_1080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_1090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_10a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_10b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_10c0(double*, const double*, const double*, const double*, const double*, const double*);

//static void _vrr_2000(double*, const double*, const double*, const double*, const double*, const double*);
//static void _vrr_2010(double*, const double*, const double*, const double*, const double*, const double*);
//static void _vrr_2020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_2030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_2040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_2050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_2060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_2070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_2080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_2090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_20a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_20b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_20c0(double*, const double*, const double*, const double*, const double*, const double*);

//static void _vrr_3000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_3010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_3020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_3030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_3040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_3050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_3060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_3070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_3080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_3090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_30a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_30b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_30c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _vrr_4000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_4010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_4020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_4030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_4040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_4050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_4060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_4070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_4080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_4090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_40a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_40b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_40c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _vrr_5000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_5010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_5020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_5030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_5040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_5050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_5060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_5070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_5080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_5090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_50a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_50b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_50c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _vrr_6000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_6010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_6020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_6030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_6040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_6050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_6060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_6070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_6080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_6090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_60a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_60b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_60c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _vrr_7000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_7010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_7020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_7030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_7040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_7050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_7060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_7070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_7080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_7090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_70a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_70b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_70c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _vrr_8000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_8010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_8020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_8030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_8040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_8050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_8060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_8070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_8080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_8090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_80a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_80b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_80c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _vrr_9000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_9010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_9020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_9030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_9040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_9050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_9060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_9070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_9080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_9090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_90a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_90b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_90c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _vrr_a000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_a010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_a020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_a030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_a040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_a050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_a060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_a070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_a080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_a090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_a0a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_a0b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_a0c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _vrr_b000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_b010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_b020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_b030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_b040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_b050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_b060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_b070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_b080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_b090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_b0a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_b0b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_b0c0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _vrr_c000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_c010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_c020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_c030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_c040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_c050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_c060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_c070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_c080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_c090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_c0a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_c0b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _vrr_c0c0(double*, const double*, const double*, const double*, const double*, const double*);

  void vrrfunc_call(const unsigned int i, double* a, const double* b, const double* c, const double* d, const double* e, const double* f) {
    return (vrrfunc[i])(a, b, c, d, e, f);
  };

  void (*vrrfunc[ANG_VRR_END * ANG_VRR_END])(double*, const double*, const double*, const double*, const double*, const double*);
};

#endif
