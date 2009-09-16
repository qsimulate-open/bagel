//
// Author : Toru Shiozaki
// Date   : May 2009
//

#include "carsphlist.h"
#include "f77.h"


void CarSphList::carsph_10(const int nloop, const double* source, double* target) {
  const int size = nloop * 3;
  const int unit = 1;
  dcopy_(&size, source, &unit, target, &unit);
}

