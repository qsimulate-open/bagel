//
// Newint - Parallel electron correlation program.
// Filename: svrrlist.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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
      svrrfunc[ANG_VRR_END        ] = &_svrr_1000;
      svrrfunc[ANG_VRR_END     + 1] = &_svrr_1010;
      svrrfunc[ANG_VRR_END     + 2] = &_svrr_1020;
      svrrfunc[ANG_VRR_END     + 3] = &_svrr_1030;
      svrrfunc[ANG_VRR_END     + 4] = &_svrr_1040;
      svrrfunc[ANG_VRR_END     + 5] = &_svrr_1050;
      svrrfunc[ANG_VRR_END     + 6] = &_svrr_1060;
      svrrfunc[ANG_VRR_END     + 7] = &_svrr_1070;
      svrrfunc[ANG_VRR_END     + 8] = &_svrr_1080;
      svrrfunc[ANG_VRR_END * 2    ] = &_svrr_2000;
      svrrfunc[ANG_VRR_END * 2 + 1] = &_svrr_2010;
      svrrfunc[ANG_VRR_END * 2 + 2] = &_svrr_2020;
      svrrfunc[ANG_VRR_END * 2 + 3] = &_svrr_2030;
      svrrfunc[ANG_VRR_END * 2 + 4] = &_svrr_2040;
      svrrfunc[ANG_VRR_END * 2 + 5] = &_svrr_2050;
      svrrfunc[ANG_VRR_END * 2 + 6] = &_svrr_2060;
      svrrfunc[ANG_VRR_END * 2 + 7] = &_svrr_2070;
      svrrfunc[ANG_VRR_END * 2 + 8] = &_svrr_2080;
      svrrfunc[ANG_VRR_END * 3    ] = &_svrr_3000;
      svrrfunc[ANG_VRR_END * 3 + 1] = &_svrr_3010;
      svrrfunc[ANG_VRR_END * 3 + 2] = &_svrr_3020;
      svrrfunc[ANG_VRR_END * 3 + 3] = &_svrr_3030;
      svrrfunc[ANG_VRR_END * 3 + 4] = &_svrr_3040;
      svrrfunc[ANG_VRR_END * 3 + 5] = &_svrr_3050;
      svrrfunc[ANG_VRR_END * 3 + 6] = &_svrr_3060;
      svrrfunc[ANG_VRR_END * 3 + 7] = &_svrr_3070;
      svrrfunc[ANG_VRR_END * 3 + 8] = &_svrr_3080;
      svrrfunc[ANG_VRR_END * 4    ] = &_svrr_4000;
      svrrfunc[ANG_VRR_END * 4 + 1] = &_svrr_4010;
      svrrfunc[ANG_VRR_END * 4 + 2] = &_svrr_4020;
      svrrfunc[ANG_VRR_END * 4 + 3] = &_svrr_4030;
      svrrfunc[ANG_VRR_END * 4 + 4] = &_svrr_4040;
      svrrfunc[ANG_VRR_END * 4 + 5] = &_svrr_4050;
      svrrfunc[ANG_VRR_END * 4 + 6] = &_svrr_4060;
      svrrfunc[ANG_VRR_END * 4 + 7] = &_svrr_4070;
      svrrfunc[ANG_VRR_END * 4 + 8] = &_svrr_4080;
      svrrfunc[ANG_VRR_END * 5    ] = &_svrr_5000;
      svrrfunc[ANG_VRR_END * 5 + 1] = &_svrr_5010;
      svrrfunc[ANG_VRR_END * 5 + 2] = &_svrr_5020;
      svrrfunc[ANG_VRR_END * 5 + 3] = &_svrr_5030;
      svrrfunc[ANG_VRR_END * 5 + 4] = &_svrr_5040;
      svrrfunc[ANG_VRR_END * 5 + 5] = &_svrr_5050;
      svrrfunc[ANG_VRR_END * 5 + 6] = &_svrr_5060;
      svrrfunc[ANG_VRR_END * 5 + 7] = &_svrr_5070;
      svrrfunc[ANG_VRR_END * 5 + 8] = &_svrr_5080;
      svrrfunc[ANG_VRR_END * 6    ] = &_svrr_6000;
      svrrfunc[ANG_VRR_END * 6 + 1] = &_svrr_6010;
      svrrfunc[ANG_VRR_END * 6 + 2] = &_svrr_6020;
      svrrfunc[ANG_VRR_END * 6 + 3] = &_svrr_6030;
      svrrfunc[ANG_VRR_END * 6 + 4] = &_svrr_6040;
      svrrfunc[ANG_VRR_END * 6 + 5] = &_svrr_6050;
      svrrfunc[ANG_VRR_END * 6 + 6] = &_svrr_6060;
      svrrfunc[ANG_VRR_END * 6 + 7] = &_svrr_6070;
      svrrfunc[ANG_VRR_END * 6 + 8] = &_svrr_6080;
      svrrfunc[ANG_VRR_END * 7    ] = &_svrr_7000;
      svrrfunc[ANG_VRR_END * 7 + 1] = &_svrr_7010;
      svrrfunc[ANG_VRR_END * 7 + 2] = &_svrr_7020;
      svrrfunc[ANG_VRR_END * 7 + 3] = &_svrr_7030;
      svrrfunc[ANG_VRR_END * 7 + 4] = &_svrr_7040;
      svrrfunc[ANG_VRR_END * 7 + 5] = &_svrr_7050;
      svrrfunc[ANG_VRR_END * 7 + 6] = &_svrr_7060;
      svrrfunc[ANG_VRR_END * 7 + 7] = &_svrr_7070;
      svrrfunc[ANG_VRR_END * 7 + 8] = &_svrr_7080;
      svrrfunc[ANG_VRR_END * 8    ] = &_svrr_8000;
      svrrfunc[ANG_VRR_END * 8 + 1] = &_svrr_8010;
      svrrfunc[ANG_VRR_END * 8 + 2] = &_svrr_8020;
      svrrfunc[ANG_VRR_END * 8 + 3] = &_svrr_8030;
      svrrfunc[ANG_VRR_END * 8 + 4] = &_svrr_8040;
      svrrfunc[ANG_VRR_END * 8 + 5] = &_svrr_8050;
      svrrfunc[ANG_VRR_END * 8 + 6] = &_svrr_8060;
      svrrfunc[ANG_VRR_END * 8 + 7] = &_svrr_8070;
      svrrfunc[ANG_VRR_END * 8 + 8] = &_svrr_8080;
}

 
SVRRList::~SVRRList() {

}

