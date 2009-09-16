//
// Author : Toru Shiozaki
// Date   : May 2009
//

#include "carsphlist.h"
#include "f77.h"


void CarSphList::carsph_00(const int nloop, const double* source, double* target) {
  const int size = nloop;
  const int unit = 1;
  dcopy_(&size, source, &unit, target, &unit);
}

