//
// Author : Toru Shiozaki
// Date   : May 2009
//

#include "carsphlist.h"
#include <cstring>

void CarSphList::carsph_10(const int nloop, const double* source, double* target) {
  const int size = nloop * 3;
  ::memcpy(target, source, (nloop * 3) * sizeof(double));
}

