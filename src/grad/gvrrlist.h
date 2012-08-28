//
// BAGEL - Parallel electron correlation program.
// Filename: gvrrlist.h
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

//
// stand alone class for VRR functions.
// can be called from outside through vrrfunc_call.
//
  
#ifndef __grad_gvrrlist_h
#define __grad_gvrrlist_h

#include <src/rysint/vrrlist.h>

struct GVRRList : public VRRListBase {
  GVRRList();
  ~GVRRList();
  
  static void _gvrr_0000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_0010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_0020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_0030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_0040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_0050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_0060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_0070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_0080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_0090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_00a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_00b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_00c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_00d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_1000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_1010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_1020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_1030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_1040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_1050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_1060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_1070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_1080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_1090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_10a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_10b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_10c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_10d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_2000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_2010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_2020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_2030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_2040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_2050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_2060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_2070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_2080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_2090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_20a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_20b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_20c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_20d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_3000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_3010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_3020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_3030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_3040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_3050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_3060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_3070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_3080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_3090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_30a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_30b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_30c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_30d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_4000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_4010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_4020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_4030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_4040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_4050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_4060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_4070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_4080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_4090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_40a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_40b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_40c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_40d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_5000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_5010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_5020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_5030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_5040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_5050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_5060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_5070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_5080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_5090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_50a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_50b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_50c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_50d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_6000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_6010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_6020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_6030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_6040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_6050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_6060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_6070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_6080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_6090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_60a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_60b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_60c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_60d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_7000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_7010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_7020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_7030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_7040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_7050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_7060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_7070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_7080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_7090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_70a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_70b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_70c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_70d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_8000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_8010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_8020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_8030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_8040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_8050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_8060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_8070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_8080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_8090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_80a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_80b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_80c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_80d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_9000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_9010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_9020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_9030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_9040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_9050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_9060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_9070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_9080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_9090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_90a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_90b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_90c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_90d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_a000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a0a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a0b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a0c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_a0d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_b000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b0a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b0b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b0c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_b0d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_c000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c0a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c0b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c0c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_c0d0(double*, const double*, const double*, const double*, const double*, const double*);

  static void _gvrr_d000(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d010(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d020(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d030(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d040(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d050(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d060(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d070(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d080(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d090(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d0a0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d0b0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d0c0(double*, const double*, const double*, const double*, const double*, const double*);
  static void _gvrr_d0d0(double*, const double*, const double*, const double*, const double*, const double*);
};

#endif
