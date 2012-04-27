//
// Newint - Parallel electron correlation program.
// Filename: svrrlist.cc
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


#include <src/slater/svrrlist.h>

SVRRList::SVRRList() {
      svrrfunc[                  0] = &_svrr_0000;
      svrrfunc[                  1] = &_svrr_0010;
      svrrfunc[                  2] = &_svrr_0020;
      svrrfunc[                  3] = &_svrr_0030;
      svrrfunc[                  4] = &_svrr_0040;
      svrrfunc[                  5] = &_svrr_0050;
      svrrfunc[                  6] = &_svrr_0060;
      svrrfunc[                  7] = &_svrr_0070;
      svrrfunc[                  8] = &_svrr_0080;
      svrrfunc[                  9] = &_svrr_0090;
      svrrfunc[                 10] = &_svrr_00a0;
      svrrfunc[                 11] = &_svrr_00b0;
      svrrfunc[                 12] = &_svrr_00c0;
      svrrfunc[ANG_VRR_END        ] = &_svrr_1000;
      svrrfunc[ANG_VRR_END     + 1] = &_svrr_1010;
      svrrfunc[ANG_VRR_END     + 2] = &_svrr_1020;
      svrrfunc[ANG_VRR_END     + 3] = &_svrr_1030;
      svrrfunc[ANG_VRR_END     + 4] = &_svrr_1040;
      svrrfunc[ANG_VRR_END     + 5] = &_svrr_1050;
      svrrfunc[ANG_VRR_END     + 6] = &_svrr_1060;
      svrrfunc[ANG_VRR_END     + 7] = &_svrr_1070;
      svrrfunc[ANG_VRR_END     + 8] = &_svrr_1080;
      svrrfunc[ANG_VRR_END     + 9] = &_svrr_1090;
      svrrfunc[ANG_VRR_END     +10] = &_svrr_10a0;
      svrrfunc[ANG_VRR_END     +11] = &_svrr_10b0;
      svrrfunc[ANG_VRR_END     +12] = &_svrr_10c0;
      svrrfunc[ANG_VRR_END * 2    ] = &_svrr_2000;
      svrrfunc[ANG_VRR_END * 2 + 1] = &_svrr_2010;
      svrrfunc[ANG_VRR_END * 2 + 2] = &_svrr_2020;
      svrrfunc[ANG_VRR_END * 2 + 3] = &_svrr_2030;
      svrrfunc[ANG_VRR_END * 2 + 4] = &_svrr_2040;
      svrrfunc[ANG_VRR_END * 2 + 5] = &_svrr_2050;
      svrrfunc[ANG_VRR_END * 2 + 6] = &_svrr_2060;
      svrrfunc[ANG_VRR_END * 2 + 7] = &_svrr_2070;
      svrrfunc[ANG_VRR_END * 2 + 8] = &_svrr_2080;
      svrrfunc[ANG_VRR_END * 2 + 9] = &_svrr_2090;
      svrrfunc[ANG_VRR_END * 2 +10] = &_svrr_20a0;
      svrrfunc[ANG_VRR_END * 2 +11] = &_svrr_20b0;
      svrrfunc[ANG_VRR_END * 2 +12] = &_svrr_20c0;
      svrrfunc[ANG_VRR_END * 3    ] = &_svrr_3000;
      svrrfunc[ANG_VRR_END * 3 + 1] = &_svrr_3010;
      svrrfunc[ANG_VRR_END * 3 + 2] = &_svrr_3020;
      svrrfunc[ANG_VRR_END * 3 + 3] = &_svrr_3030;
      svrrfunc[ANG_VRR_END * 3 + 4] = &_svrr_3040;
      svrrfunc[ANG_VRR_END * 3 + 5] = &_svrr_3050;
      svrrfunc[ANG_VRR_END * 3 + 6] = &_svrr_3060;
      svrrfunc[ANG_VRR_END * 3 + 7] = &_svrr_3070;
      svrrfunc[ANG_VRR_END * 3 + 8] = &_svrr_3080;
      svrrfunc[ANG_VRR_END * 3 + 9] = &_svrr_3090;
      svrrfunc[ANG_VRR_END * 3 +10] = &_svrr_30a0;
      svrrfunc[ANG_VRR_END * 3 +11] = &_svrr_30b0;
      svrrfunc[ANG_VRR_END * 3 +12] = &_svrr_30c0;
      svrrfunc[ANG_VRR_END * 4    ] = &_svrr_4000;
      svrrfunc[ANG_VRR_END * 4 + 1] = &_svrr_4010;
      svrrfunc[ANG_VRR_END * 4 + 2] = &_svrr_4020;
      svrrfunc[ANG_VRR_END * 4 + 3] = &_svrr_4030;
      svrrfunc[ANG_VRR_END * 4 + 4] = &_svrr_4040;
      svrrfunc[ANG_VRR_END * 4 + 5] = &_svrr_4050;
      svrrfunc[ANG_VRR_END * 4 + 6] = &_svrr_4060;
      svrrfunc[ANG_VRR_END * 4 + 7] = &_svrr_4070;
      svrrfunc[ANG_VRR_END * 4 + 8] = &_svrr_4080;
      svrrfunc[ANG_VRR_END * 4 + 9] = &_svrr_4090;
      svrrfunc[ANG_VRR_END * 4 +10] = &_svrr_40a0;
      svrrfunc[ANG_VRR_END * 4 +11] = &_svrr_40b0;
      svrrfunc[ANG_VRR_END * 4 +12] = &_svrr_40c0;
      svrrfunc[ANG_VRR_END * 5    ] = &_svrr_5000;
      svrrfunc[ANG_VRR_END * 5 + 1] = &_svrr_5010;
      svrrfunc[ANG_VRR_END * 5 + 2] = &_svrr_5020;
      svrrfunc[ANG_VRR_END * 5 + 3] = &_svrr_5030;
      svrrfunc[ANG_VRR_END * 5 + 4] = &_svrr_5040;
      svrrfunc[ANG_VRR_END * 5 + 5] = &_svrr_5050;
      svrrfunc[ANG_VRR_END * 5 + 6] = &_svrr_5060;
      svrrfunc[ANG_VRR_END * 5 + 7] = &_svrr_5070;
      svrrfunc[ANG_VRR_END * 5 + 8] = &_svrr_5080;
      svrrfunc[ANG_VRR_END * 5 + 9] = &_svrr_5090;
      svrrfunc[ANG_VRR_END * 5 +10] = &_svrr_50a0;
      svrrfunc[ANG_VRR_END * 5 +11] = &_svrr_50b0;
      svrrfunc[ANG_VRR_END * 5 +12] = &_svrr_50c0;
      svrrfunc[ANG_VRR_END * 6    ] = &_svrr_6000;
      svrrfunc[ANG_VRR_END * 6 + 1] = &_svrr_6010;
      svrrfunc[ANG_VRR_END * 6 + 2] = &_svrr_6020;
      svrrfunc[ANG_VRR_END * 6 + 3] = &_svrr_6030;
      svrrfunc[ANG_VRR_END * 6 + 4] = &_svrr_6040;
      svrrfunc[ANG_VRR_END * 6 + 5] = &_svrr_6050;
      svrrfunc[ANG_VRR_END * 6 + 6] = &_svrr_6060;
      svrrfunc[ANG_VRR_END * 6 + 7] = &_svrr_6070;
      svrrfunc[ANG_VRR_END * 6 + 8] = &_svrr_6080;
      svrrfunc[ANG_VRR_END * 6 + 9] = &_svrr_6090;
      svrrfunc[ANG_VRR_END * 6 +10] = &_svrr_60a0;
      svrrfunc[ANG_VRR_END * 6 +11] = &_svrr_60b0;
      svrrfunc[ANG_VRR_END * 6 +12] = &_svrr_60c0;
      svrrfunc[ANG_VRR_END * 7    ] = &_svrr_7000;
      svrrfunc[ANG_VRR_END * 7 + 1] = &_svrr_7010;
      svrrfunc[ANG_VRR_END * 7 + 2] = &_svrr_7020;
      svrrfunc[ANG_VRR_END * 7 + 3] = &_svrr_7030;
      svrrfunc[ANG_VRR_END * 7 + 4] = &_svrr_7040;
      svrrfunc[ANG_VRR_END * 7 + 5] = &_svrr_7050;
      svrrfunc[ANG_VRR_END * 7 + 6] = &_svrr_7060;
      svrrfunc[ANG_VRR_END * 7 + 7] = &_svrr_7070;
      svrrfunc[ANG_VRR_END * 7 + 8] = &_svrr_7080;
      svrrfunc[ANG_VRR_END * 7 + 9] = &_svrr_7090;
      svrrfunc[ANG_VRR_END * 7 +10] = &_svrr_70a0;
      svrrfunc[ANG_VRR_END * 7 +11] = &_svrr_70b0;
      svrrfunc[ANG_VRR_END * 7 +12] = &_svrr_70c0;
      svrrfunc[ANG_VRR_END * 8    ] = &_svrr_8000;
      svrrfunc[ANG_VRR_END * 8 + 1] = &_svrr_8010;
      svrrfunc[ANG_VRR_END * 8 + 2] = &_svrr_8020;
      svrrfunc[ANG_VRR_END * 8 + 3] = &_svrr_8030;
      svrrfunc[ANG_VRR_END * 8 + 4] = &_svrr_8040;
      svrrfunc[ANG_VRR_END * 8 + 5] = &_svrr_8050;
      svrrfunc[ANG_VRR_END * 8 + 6] = &_svrr_8060;
      svrrfunc[ANG_VRR_END * 8 + 7] = &_svrr_8070;
      svrrfunc[ANG_VRR_END * 8 + 8] = &_svrr_8080;
      svrrfunc[ANG_VRR_END * 8 + 9] = &_svrr_8090;
      svrrfunc[ANG_VRR_END * 8 +10] = &_svrr_80a0;
      svrrfunc[ANG_VRR_END * 8 +11] = &_svrr_80b0;
      svrrfunc[ANG_VRR_END * 8 +12] = &_svrr_80c0;
      svrrfunc[ANG_VRR_END * 9    ] = &_svrr_9000;
      svrrfunc[ANG_VRR_END * 9 + 1] = &_svrr_9010;
      svrrfunc[ANG_VRR_END * 9 + 2] = &_svrr_9020;
      svrrfunc[ANG_VRR_END * 9 + 3] = &_svrr_9030;
      svrrfunc[ANG_VRR_END * 9 + 4] = &_svrr_9040;
      svrrfunc[ANG_VRR_END * 9 + 5] = &_svrr_9050;
      svrrfunc[ANG_VRR_END * 9 + 6] = &_svrr_9060;
      svrrfunc[ANG_VRR_END * 9 + 7] = &_svrr_9070;
      svrrfunc[ANG_VRR_END * 9 + 8] = &_svrr_9080;
      svrrfunc[ANG_VRR_END * 9 + 9] = &_svrr_9090;
      svrrfunc[ANG_VRR_END * 9 +10] = &_svrr_90a0;
      svrrfunc[ANG_VRR_END * 9 +11] = &_svrr_90b0;
      svrrfunc[ANG_VRR_END * 9 +12] = &_svrr_90c0;
      svrrfunc[ANG_VRR_END *10    ] = &_svrr_a000;
      svrrfunc[ANG_VRR_END *10 + 1] = &_svrr_a010;
      svrrfunc[ANG_VRR_END *10 + 2] = &_svrr_a020;
      svrrfunc[ANG_VRR_END *10 + 3] = &_svrr_a030;
      svrrfunc[ANG_VRR_END *10 + 4] = &_svrr_a040;
      svrrfunc[ANG_VRR_END *10 + 5] = &_svrr_a050;
      svrrfunc[ANG_VRR_END *10 + 6] = &_svrr_a060;
      svrrfunc[ANG_VRR_END *10 + 7] = &_svrr_a070;
      svrrfunc[ANG_VRR_END *10 + 8] = &_svrr_a080;
      svrrfunc[ANG_VRR_END *10 + 9] = &_svrr_a090;
      svrrfunc[ANG_VRR_END *10 +10] = &_svrr_a0a0;
      svrrfunc[ANG_VRR_END *10 +11] = &_svrr_a0b0;
      svrrfunc[ANG_VRR_END *10 +12] = &_svrr_a0c0;
      svrrfunc[ANG_VRR_END *11    ] = &_svrr_b000;
      svrrfunc[ANG_VRR_END *11 + 1] = &_svrr_b010;
      svrrfunc[ANG_VRR_END *11 + 2] = &_svrr_b020;
      svrrfunc[ANG_VRR_END *11 + 3] = &_svrr_b030;
      svrrfunc[ANG_VRR_END *11 + 4] = &_svrr_b040;
      svrrfunc[ANG_VRR_END *11 + 5] = &_svrr_b050;
      svrrfunc[ANG_VRR_END *11 + 6] = &_svrr_b060;
      svrrfunc[ANG_VRR_END *11 + 7] = &_svrr_b070;
      svrrfunc[ANG_VRR_END *11 + 8] = &_svrr_b080;
      svrrfunc[ANG_VRR_END *11 + 9] = &_svrr_b090;
      svrrfunc[ANG_VRR_END *11 +10] = &_svrr_b0a0;
      svrrfunc[ANG_VRR_END *11 +11] = &_svrr_b0b0;
      svrrfunc[ANG_VRR_END *11 +12] = &_svrr_b0c0;
      svrrfunc[ANG_VRR_END *12    ] = &_svrr_c000;
      svrrfunc[ANG_VRR_END *12 + 1] = &_svrr_c010;
      svrrfunc[ANG_VRR_END *12 + 2] = &_svrr_c020;
      svrrfunc[ANG_VRR_END *12 + 3] = &_svrr_c030;
      svrrfunc[ANG_VRR_END *12 + 4] = &_svrr_c040;
      svrrfunc[ANG_VRR_END *12 + 5] = &_svrr_c050;
      svrrfunc[ANG_VRR_END *12 + 6] = &_svrr_c060;
      svrrfunc[ANG_VRR_END *12 + 7] = &_svrr_c070;
      svrrfunc[ANG_VRR_END *12 + 8] = &_svrr_c080;
      svrrfunc[ANG_VRR_END *12 + 9] = &_svrr_c090;
      svrrfunc[ANG_VRR_END *12 +10] = &_svrr_c0a0;
      svrrfunc[ANG_VRR_END *12 +11] = &_svrr_c0b0;
      svrrfunc[ANG_VRR_END *12 +12] = &_svrr_c0c0;
}

 
SVRRList::~SVRRList() {

}

