//
// Author: Toru Shiozaki
// Date  : May 2009
//

#include <src/rysint/carsphlist.h>

// same convention as HRR

CarSphList::CarSphList() {
      carsphfunc[ANG_HRR_END * 0 + 0] = &carsph_00;
      carsphfunc[ANG_HRR_END * 1 + 0] = &carsph_10;
      carsphfunc[ANG_HRR_END * 2 + 0] = &carsph_20;
      carsphfunc[ANG_HRR_END * 3 + 0] = &carsph_30;
      carsphfunc[ANG_HRR_END * 4 + 0] = &carsph_40;
      carsphfunc[ANG_HRR_END * 1 + 1] = &carsph_11;
      carsphfunc[ANG_HRR_END * 2 + 1] = &carsph_21;
      carsphfunc[ANG_HRR_END * 3 + 1] = &carsph_31;
      carsphfunc[ANG_HRR_END * 4 + 1] = &carsph_41;
      carsphfunc[ANG_HRR_END * 2 + 2] = &carsph_22;
      carsphfunc[ANG_HRR_END * 3 + 2] = &carsph_32;
      carsphfunc[ANG_HRR_END * 4 + 2] = &carsph_42;
      carsphfunc[ANG_HRR_END * 3 + 3] = &carsph_33;
      carsphfunc[ANG_HRR_END * 4 + 3] = &carsph_43;
      carsphfunc[ANG_HRR_END * 4 + 4] = &carsph_44;
}

 
CarSphList::~CarSphList() {

}

