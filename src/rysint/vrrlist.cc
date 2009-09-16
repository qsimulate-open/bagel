//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include <src/rysint/vrrlist.h>

VRRList::VRRList() {
//    vrrfunc[                  0] = &_vrr_0000;
//    vrrfunc[                  1] = &_vrr_0010;
//    vrrfunc[                  2] = &_vrr_0020;
//    vrrfunc[                  3] = &_vrr_0030;
      vrrfunc[                  4] = &_vrr_0040;
      vrrfunc[                  5] = &_vrr_0050;
      vrrfunc[                  6] = &_vrr_0060;
      vrrfunc[                  7] = &_vrr_0070;
      vrrfunc[                  8] = &_vrr_0080;
//    vrrfunc[ANG_VRR_END        ] = &_vrr_1000;
//    vrrfunc[ANG_VRR_END     + 1] = &_vrr_1010;
//    vrrfunc[ANG_VRR_END     + 2] = &_vrr_1020;
      vrrfunc[ANG_VRR_END     + 3] = &_vrr_1030;
      vrrfunc[ANG_VRR_END     + 4] = &_vrr_1040;
      vrrfunc[ANG_VRR_END     + 5] = &_vrr_1050;
      vrrfunc[ANG_VRR_END     + 6] = &_vrr_1060;
      vrrfunc[ANG_VRR_END     + 7] = &_vrr_1070;
      vrrfunc[ANG_VRR_END     + 8] = &_vrr_1080;
//    vrrfunc[ANG_VRR_END * 2    ] = &_vrr_2000;
//    vrrfunc[ANG_VRR_END * 2 + 1] = &_vrr_2010;
//    vrrfunc[ANG_VRR_END * 2 + 2] = &_vrr_2020;
      vrrfunc[ANG_VRR_END * 2 + 3] = &_vrr_2030;
      vrrfunc[ANG_VRR_END * 2 + 4] = &_vrr_2040;
      vrrfunc[ANG_VRR_END * 2 + 5] = &_vrr_2050;
      vrrfunc[ANG_VRR_END * 2 + 6] = &_vrr_2060;
      vrrfunc[ANG_VRR_END * 2 + 7] = &_vrr_2070;
      vrrfunc[ANG_VRR_END * 2 + 8] = &_vrr_2080;
//    vrrfunc[ANG_VRR_END * 3    ] = &_vrr_3000;
      vrrfunc[ANG_VRR_END * 3 + 1] = &_vrr_3010;
      vrrfunc[ANG_VRR_END * 3 + 2] = &_vrr_3020;
      vrrfunc[ANG_VRR_END * 3 + 3] = &_vrr_3030;
      vrrfunc[ANG_VRR_END * 3 + 4] = &_vrr_3040;
      vrrfunc[ANG_VRR_END * 3 + 5] = &_vrr_3050;
      vrrfunc[ANG_VRR_END * 3 + 6] = &_vrr_3060;
      vrrfunc[ANG_VRR_END * 3 + 7] = &_vrr_3070;
      vrrfunc[ANG_VRR_END * 3 + 8] = &_vrr_3080;
      vrrfunc[ANG_VRR_END * 4    ] = &_vrr_4000;
      vrrfunc[ANG_VRR_END * 4 + 1] = &_vrr_4010;
      vrrfunc[ANG_VRR_END * 4 + 2] = &_vrr_4020;
      vrrfunc[ANG_VRR_END * 4 + 3] = &_vrr_4030;
      vrrfunc[ANG_VRR_END * 4 + 4] = &_vrr_4040;
      vrrfunc[ANG_VRR_END * 4 + 5] = &_vrr_4050;
      vrrfunc[ANG_VRR_END * 4 + 6] = &_vrr_4060;
      vrrfunc[ANG_VRR_END * 4 + 7] = &_vrr_4070;
      vrrfunc[ANG_VRR_END * 4 + 8] = &_vrr_4080;
      vrrfunc[ANG_VRR_END * 5    ] = &_vrr_5000;
      vrrfunc[ANG_VRR_END * 5 + 1] = &_vrr_5010;
      vrrfunc[ANG_VRR_END * 5 + 2] = &_vrr_5020;
      vrrfunc[ANG_VRR_END * 5 + 3] = &_vrr_5030;
      vrrfunc[ANG_VRR_END * 5 + 4] = &_vrr_5040;
      vrrfunc[ANG_VRR_END * 5 + 5] = &_vrr_5050;
      vrrfunc[ANG_VRR_END * 5 + 6] = &_vrr_5060;
      vrrfunc[ANG_VRR_END * 5 + 7] = &_vrr_5070;
      vrrfunc[ANG_VRR_END * 5 + 8] = &_vrr_5080;
      vrrfunc[ANG_VRR_END * 6    ] = &_vrr_6000;
      vrrfunc[ANG_VRR_END * 6 + 1] = &_vrr_6010;
      vrrfunc[ANG_VRR_END * 6 + 2] = &_vrr_6020;
      vrrfunc[ANG_VRR_END * 6 + 3] = &_vrr_6030;
      vrrfunc[ANG_VRR_END * 6 + 4] = &_vrr_6040;
      vrrfunc[ANG_VRR_END * 6 + 5] = &_vrr_6050;
      vrrfunc[ANG_VRR_END * 6 + 6] = &_vrr_6060;
      vrrfunc[ANG_VRR_END * 6 + 7] = &_vrr_6070;
      vrrfunc[ANG_VRR_END * 6 + 8] = &_vrr_6080;
      vrrfunc[ANG_VRR_END * 7    ] = &_vrr_7000;
      vrrfunc[ANG_VRR_END * 7 + 1] = &_vrr_7010;
      vrrfunc[ANG_VRR_END * 7 + 2] = &_vrr_7020;
      vrrfunc[ANG_VRR_END * 7 + 3] = &_vrr_7030;
      vrrfunc[ANG_VRR_END * 7 + 4] = &_vrr_7040;
      vrrfunc[ANG_VRR_END * 7 + 5] = &_vrr_7050;
      vrrfunc[ANG_VRR_END * 7 + 6] = &_vrr_7060;
      vrrfunc[ANG_VRR_END * 7 + 7] = &_vrr_7070;
      vrrfunc[ANG_VRR_END * 7 + 8] = &_vrr_7080;
      vrrfunc[ANG_VRR_END * 8    ] = &_vrr_8000;
      vrrfunc[ANG_VRR_END * 8 + 1] = &_vrr_8010;
      vrrfunc[ANG_VRR_END * 8 + 2] = &_vrr_8020;
      vrrfunc[ANG_VRR_END * 8 + 3] = &_vrr_8030;
      vrrfunc[ANG_VRR_END * 8 + 4] = &_vrr_8040;
      vrrfunc[ANG_VRR_END * 8 + 5] = &_vrr_8050;
      vrrfunc[ANG_VRR_END * 8 + 6] = &_vrr_8060;
      vrrfunc[ANG_VRR_END * 8 + 7] = &_vrr_8070;
      vrrfunc[ANG_VRR_END * 8 + 8] = &_vrr_8080;
}

 
VRRList::~VRRList() {

}

