//
// Author: Toru Shiozaki
// Date  : April 2009
//

#include <src/rysint/hrrlist.h>

HRRList::HRRList() {
      hrrfunc[ANG_HRR_END * 1 + 1] = &perform_HRR_20_11;
      hrrfunc[ANG_HRR_END * 2 + 1] = &perform_HRR_30_21;
      hrrfunc[ANG_HRR_END * 2 + 2] = &perform_HRR_40_22;
      hrrfunc[ANG_HRR_END * 3 + 1] = &perform_HRR_40_31;
      hrrfunc[ANG_HRR_END * 3 + 2] = &perform_HRR_50_32;
      hrrfunc[ANG_HRR_END * 4 + 1] = &perform_HRR_50_41;
      hrrfunc[ANG_HRR_END * 3 + 3] = &perform_HRR_60_33;
      hrrfunc[ANG_HRR_END * 4 + 2] = &perform_HRR_60_42;
      hrrfunc[ANG_HRR_END * 4 + 3] = &perform_HRR_70_43;
      hrrfunc[ANG_HRR_END * 4 + 4] = &perform_HRR_80_44;
}

 
HRRList::~HRRList() {

}

