//
// BAGEL - Parallel electron correlation program.
// Filename: gvrrlist.cc
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

GVRRList::GVRRList() {
      vrrfunc[                  0] = &_gvrr_0000;
      vrrfunc[                  1] = &_gvrr_0010;
      vrrfunc[                  2] = &_gvrr_0020;
      vrrfunc[                  3] = &_gvrr_0030;
      vrrfunc[                  4] = &_gvrr_0040;
      vrrfunc[                  5] = &_gvrr_0050;
      vrrfunc[                  6] = &_gvrr_0060;
      vrrfunc[                  7] = &_gvrr_0070;
      vrrfunc[                  8] = &_gvrr_0080;
      vrrfunc[                  9] = &_gvrr_0090;
      vrrfunc[                 10] = &_gvrr_00a0;
      vrrfunc[                 11] = &_gvrr_00b0;
      vrrfunc[                 12] = &_gvrr_00c0;
      vrrfunc[                 13] = &_gvrr_00d0;
      vrrfunc[ANG_VRR_END        ] = &_gvrr_1000;
      vrrfunc[ANG_VRR_END     + 1] = &_gvrr_1010;
      vrrfunc[ANG_VRR_END     + 2] = &_gvrr_1020;
      vrrfunc[ANG_VRR_END     + 3] = &_gvrr_1030;
      vrrfunc[ANG_VRR_END     + 4] = &_gvrr_1040;
      vrrfunc[ANG_VRR_END     + 5] = &_gvrr_1050;
      vrrfunc[ANG_VRR_END     + 6] = &_gvrr_1060;
      vrrfunc[ANG_VRR_END     + 7] = &_gvrr_1070;
      vrrfunc[ANG_VRR_END     + 8] = &_gvrr_1080;
      vrrfunc[ANG_VRR_END     + 9] = &_gvrr_1090;
      vrrfunc[ANG_VRR_END     +10] = &_gvrr_10a0;
      vrrfunc[ANG_VRR_END     +11] = &_gvrr_10b0;
      vrrfunc[ANG_VRR_END     +12] = &_gvrr_10c0;
      vrrfunc[ANG_VRR_END     +13] = &_gvrr_10d0;
      vrrfunc[ANG_VRR_END * 2    ] = &_gvrr_2000;
      vrrfunc[ANG_VRR_END * 2 + 1] = &_gvrr_2010;
      vrrfunc[ANG_VRR_END * 2 + 2] = &_gvrr_2020;
      vrrfunc[ANG_VRR_END * 2 + 3] = &_gvrr_2030;
      vrrfunc[ANG_VRR_END * 2 + 4] = &_gvrr_2040;
      vrrfunc[ANG_VRR_END * 2 + 5] = &_gvrr_2050;
      vrrfunc[ANG_VRR_END * 2 + 6] = &_gvrr_2060;
      vrrfunc[ANG_VRR_END * 2 + 7] = &_gvrr_2070;
      vrrfunc[ANG_VRR_END * 2 + 8] = &_gvrr_2080;
      vrrfunc[ANG_VRR_END * 2 + 9] = &_gvrr_2090;
      vrrfunc[ANG_VRR_END * 2 +10] = &_gvrr_20a0;
      vrrfunc[ANG_VRR_END * 2 +11] = &_gvrr_20b0;
      vrrfunc[ANG_VRR_END * 2 +12] = &_gvrr_20c0;
      vrrfunc[ANG_VRR_END * 2 +13] = &_gvrr_20d0;
      vrrfunc[ANG_VRR_END * 3    ] = &_gvrr_3000;
      vrrfunc[ANG_VRR_END * 3 + 1] = &_gvrr_3010;
      vrrfunc[ANG_VRR_END * 3 + 2] = &_gvrr_3020;
      vrrfunc[ANG_VRR_END * 3 + 3] = &_gvrr_3030;
      vrrfunc[ANG_VRR_END * 3 + 4] = &_gvrr_3040;
      vrrfunc[ANG_VRR_END * 3 + 5] = &_gvrr_3050;
      vrrfunc[ANG_VRR_END * 3 + 6] = &_gvrr_3060;
      vrrfunc[ANG_VRR_END * 3 + 7] = &_gvrr_3070;
      vrrfunc[ANG_VRR_END * 3 + 8] = &_gvrr_3080;
      vrrfunc[ANG_VRR_END * 3 + 9] = &_gvrr_3090;
      vrrfunc[ANG_VRR_END * 3 +10] = &_gvrr_30a0;
      vrrfunc[ANG_VRR_END * 3 +11] = &_gvrr_30b0;
      vrrfunc[ANG_VRR_END * 3 +12] = &_gvrr_30c0;
      vrrfunc[ANG_VRR_END * 3 +13] = &_gvrr_30d0;
      vrrfunc[ANG_VRR_END * 4    ] = &_gvrr_4000;
      vrrfunc[ANG_VRR_END * 4 + 1] = &_gvrr_4010;
      vrrfunc[ANG_VRR_END * 4 + 2] = &_gvrr_4020;
      vrrfunc[ANG_VRR_END * 4 + 3] = &_gvrr_4030;
      vrrfunc[ANG_VRR_END * 4 + 4] = &_gvrr_4040;
      vrrfunc[ANG_VRR_END * 4 + 5] = &_gvrr_4050;
      vrrfunc[ANG_VRR_END * 4 + 6] = &_gvrr_4060;
      vrrfunc[ANG_VRR_END * 4 + 7] = &_gvrr_4070;
      vrrfunc[ANG_VRR_END * 4 + 8] = &_gvrr_4080;
      vrrfunc[ANG_VRR_END * 4 + 9] = &_gvrr_4090;
      vrrfunc[ANG_VRR_END * 4 +10] = &_gvrr_40a0;
      vrrfunc[ANG_VRR_END * 4 +11] = &_gvrr_40b0;
      vrrfunc[ANG_VRR_END * 4 +12] = &_gvrr_40c0;
      vrrfunc[ANG_VRR_END * 4 +13] = &_gvrr_40d0;
      vrrfunc[ANG_VRR_END * 5    ] = &_gvrr_5000;
      vrrfunc[ANG_VRR_END * 5 + 1] = &_gvrr_5010;
      vrrfunc[ANG_VRR_END * 5 + 2] = &_gvrr_5020;
      vrrfunc[ANG_VRR_END * 5 + 3] = &_gvrr_5030;
      vrrfunc[ANG_VRR_END * 5 + 4] = &_gvrr_5040;
      vrrfunc[ANG_VRR_END * 5 + 5] = &_gvrr_5050;
      vrrfunc[ANG_VRR_END * 5 + 6] = &_gvrr_5060;
      vrrfunc[ANG_VRR_END * 5 + 7] = &_gvrr_5070;
      vrrfunc[ANG_VRR_END * 5 + 8] = &_gvrr_5080;
      vrrfunc[ANG_VRR_END * 5 + 9] = &_gvrr_5090;
      vrrfunc[ANG_VRR_END * 5 +10] = &_gvrr_50a0;
      vrrfunc[ANG_VRR_END * 5 +11] = &_gvrr_50b0;
      vrrfunc[ANG_VRR_END * 5 +12] = &_gvrr_50c0;
      vrrfunc[ANG_VRR_END * 5 +13] = &_gvrr_50d0;
      vrrfunc[ANG_VRR_END * 6    ] = &_gvrr_6000;
      vrrfunc[ANG_VRR_END * 6 + 1] = &_gvrr_6010;
      vrrfunc[ANG_VRR_END * 6 + 2] = &_gvrr_6020;
      vrrfunc[ANG_VRR_END * 6 + 3] = &_gvrr_6030;
      vrrfunc[ANG_VRR_END * 6 + 4] = &_gvrr_6040;
      vrrfunc[ANG_VRR_END * 6 + 5] = &_gvrr_6050;
      vrrfunc[ANG_VRR_END * 6 + 6] = &_gvrr_6060;
      vrrfunc[ANG_VRR_END * 6 + 7] = &_gvrr_6070;
      vrrfunc[ANG_VRR_END * 6 + 8] = &_gvrr_6080;
      vrrfunc[ANG_VRR_END * 6 + 9] = &_gvrr_6090;
      vrrfunc[ANG_VRR_END * 6 +10] = &_gvrr_60a0;
      vrrfunc[ANG_VRR_END * 6 +11] = &_gvrr_60b0;
      vrrfunc[ANG_VRR_END * 6 +12] = &_gvrr_60c0;
      vrrfunc[ANG_VRR_END * 6 +13] = &_gvrr_60d0;
      vrrfunc[ANG_VRR_END * 7    ] = &_gvrr_7000;
      vrrfunc[ANG_VRR_END * 7 + 1] = &_gvrr_7010;
      vrrfunc[ANG_VRR_END * 7 + 2] = &_gvrr_7020;
      vrrfunc[ANG_VRR_END * 7 + 3] = &_gvrr_7030;
      vrrfunc[ANG_VRR_END * 7 + 4] = &_gvrr_7040;
      vrrfunc[ANG_VRR_END * 7 + 5] = &_gvrr_7050;
      vrrfunc[ANG_VRR_END * 7 + 6] = &_gvrr_7060;
      vrrfunc[ANG_VRR_END * 7 + 7] = &_gvrr_7070;
      vrrfunc[ANG_VRR_END * 7 + 8] = &_gvrr_7080;
      vrrfunc[ANG_VRR_END * 7 + 9] = &_gvrr_7090;
      vrrfunc[ANG_VRR_END * 7 +10] = &_gvrr_70a0;
      vrrfunc[ANG_VRR_END * 7 +11] = &_gvrr_70b0;
      vrrfunc[ANG_VRR_END * 7 +12] = &_gvrr_70c0;
      vrrfunc[ANG_VRR_END * 7 +13] = &_gvrr_70d0;
      vrrfunc[ANG_VRR_END * 8    ] = &_gvrr_8000;
      vrrfunc[ANG_VRR_END * 8 + 1] = &_gvrr_8010;
      vrrfunc[ANG_VRR_END * 8 + 2] = &_gvrr_8020;
      vrrfunc[ANG_VRR_END * 8 + 3] = &_gvrr_8030;
      vrrfunc[ANG_VRR_END * 8 + 4] = &_gvrr_8040;
      vrrfunc[ANG_VRR_END * 8 + 5] = &_gvrr_8050;
      vrrfunc[ANG_VRR_END * 8 + 6] = &_gvrr_8060;
      vrrfunc[ANG_VRR_END * 8 + 7] = &_gvrr_8070;
      vrrfunc[ANG_VRR_END * 8 + 8] = &_gvrr_8080;
      vrrfunc[ANG_VRR_END * 8 + 9] = &_gvrr_8090;
      vrrfunc[ANG_VRR_END * 8 +10] = &_gvrr_80a0;
      vrrfunc[ANG_VRR_END * 8 +11] = &_gvrr_80b0;
      vrrfunc[ANG_VRR_END * 8 +12] = &_gvrr_80c0;
      vrrfunc[ANG_VRR_END * 8 +13] = &_gvrr_80d0;
      vrrfunc[ANG_VRR_END * 9    ] = &_gvrr_9000;
      vrrfunc[ANG_VRR_END * 9 + 1] = &_gvrr_9010;
      vrrfunc[ANG_VRR_END * 9 + 2] = &_gvrr_9020;
      vrrfunc[ANG_VRR_END * 9 + 3] = &_gvrr_9030;
      vrrfunc[ANG_VRR_END * 9 + 4] = &_gvrr_9040;
      vrrfunc[ANG_VRR_END * 9 + 5] = &_gvrr_9050;
      vrrfunc[ANG_VRR_END * 9 + 6] = &_gvrr_9060;
      vrrfunc[ANG_VRR_END * 9 + 7] = &_gvrr_9070;
      vrrfunc[ANG_VRR_END * 9 + 8] = &_gvrr_9080;
      vrrfunc[ANG_VRR_END * 9 + 9] = &_gvrr_9090;
      vrrfunc[ANG_VRR_END * 9 +10] = &_gvrr_90a0;
      vrrfunc[ANG_VRR_END * 9 +11] = &_gvrr_90b0;
      vrrfunc[ANG_VRR_END * 9 +12] = &_gvrr_90c0;
      vrrfunc[ANG_VRR_END * 9 +13] = &_gvrr_90d0;
      vrrfunc[ANG_VRR_END *10    ] = &_gvrr_a000;
      vrrfunc[ANG_VRR_END *10 + 1] = &_gvrr_a010;
      vrrfunc[ANG_VRR_END *10 + 2] = &_gvrr_a020;
      vrrfunc[ANG_VRR_END *10 + 3] = &_gvrr_a030;
      vrrfunc[ANG_VRR_END *10 + 4] = &_gvrr_a040;
      vrrfunc[ANG_VRR_END *10 + 5] = &_gvrr_a050;
      vrrfunc[ANG_VRR_END *10 + 6] = &_gvrr_a060;
      vrrfunc[ANG_VRR_END *10 + 7] = &_gvrr_a070;
      vrrfunc[ANG_VRR_END *10 + 8] = &_gvrr_a080;
      vrrfunc[ANG_VRR_END *10 + 9] = &_gvrr_a090;
      vrrfunc[ANG_VRR_END *10 +10] = &_gvrr_a0a0;
      vrrfunc[ANG_VRR_END *10 +11] = &_gvrr_a0b0;
      vrrfunc[ANG_VRR_END *10 +12] = &_gvrr_a0c0;
      vrrfunc[ANG_VRR_END *10 +13] = &_gvrr_a0d0;
      vrrfunc[ANG_VRR_END *11    ] = &_gvrr_b000;
      vrrfunc[ANG_VRR_END *11 + 1] = &_gvrr_b010;
      vrrfunc[ANG_VRR_END *11 + 2] = &_gvrr_b020;
      vrrfunc[ANG_VRR_END *11 + 3] = &_gvrr_b030;
      vrrfunc[ANG_VRR_END *11 + 4] = &_gvrr_b040;
      vrrfunc[ANG_VRR_END *11 + 5] = &_gvrr_b050;
      vrrfunc[ANG_VRR_END *11 + 6] = &_gvrr_b060;
      vrrfunc[ANG_VRR_END *11 + 7] = &_gvrr_b070;
      vrrfunc[ANG_VRR_END *11 + 8] = &_gvrr_b080;
      vrrfunc[ANG_VRR_END *11 + 9] = &_gvrr_b090;
      vrrfunc[ANG_VRR_END *11 +10] = &_gvrr_b0a0;
      vrrfunc[ANG_VRR_END *11 +11] = &_gvrr_b0b0;
      vrrfunc[ANG_VRR_END *11 +12] = &_gvrr_b0c0;
      vrrfunc[ANG_VRR_END *11 +13] = &_gvrr_b0d0;
      vrrfunc[ANG_VRR_END *12    ] = &_gvrr_c000;
      vrrfunc[ANG_VRR_END *12 + 1] = &_gvrr_c010;
      vrrfunc[ANG_VRR_END *12 + 2] = &_gvrr_c020;
      vrrfunc[ANG_VRR_END *12 + 3] = &_gvrr_c030;
      vrrfunc[ANG_VRR_END *12 + 4] = &_gvrr_c040;
      vrrfunc[ANG_VRR_END *12 + 5] = &_gvrr_c050;
      vrrfunc[ANG_VRR_END *12 + 6] = &_gvrr_c060;
      vrrfunc[ANG_VRR_END *12 + 7] = &_gvrr_c070;
      vrrfunc[ANG_VRR_END *12 + 8] = &_gvrr_c080;
      vrrfunc[ANG_VRR_END *12 + 9] = &_gvrr_c090;
      vrrfunc[ANG_VRR_END *12 +10] = &_gvrr_c0a0;
      vrrfunc[ANG_VRR_END *12 +11] = &_gvrr_c0b0;
      vrrfunc[ANG_VRR_END *12 +12] = &_gvrr_c0c0;
      vrrfunc[ANG_VRR_END *12 +13] = &_gvrr_c0d0;
      vrrfunc[ANG_VRR_END *13    ] = &_gvrr_d000;
      vrrfunc[ANG_VRR_END *13 + 1] = &_gvrr_d010;
      vrrfunc[ANG_VRR_END *13 + 2] = &_gvrr_d020;
      vrrfunc[ANG_VRR_END *13 + 3] = &_gvrr_d030;
      vrrfunc[ANG_VRR_END *13 + 4] = &_gvrr_d040;
      vrrfunc[ANG_VRR_END *13 + 5] = &_gvrr_d050;
      vrrfunc[ANG_VRR_END *13 + 6] = &_gvrr_d060;
      vrrfunc[ANG_VRR_END *13 + 7] = &_gvrr_d070;
      vrrfunc[ANG_VRR_END *13 + 8] = &_gvrr_d080;
      vrrfunc[ANG_VRR_END *13 + 9] = &_gvrr_d090;
      vrrfunc[ANG_VRR_END *13 +10] = &_gvrr_d0a0;
      vrrfunc[ANG_VRR_END *13 +11] = &_gvrr_d0b0;
      vrrfunc[ANG_VRR_END *13 +12] = &_gvrr_d0c0;
      vrrfunc[ANG_VRR_END *13 +13] = &_gvrr_d0d0;
}


GVRRList::~GVRRList() {

}

