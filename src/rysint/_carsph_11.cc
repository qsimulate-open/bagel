//
// Author : Toru Shiozaki
// Date   : May 2009
//

#include "carsphlist.h"
#include <cstring>

void CarSphList::carsph_11(const int nloop, const double* source, double* target) {
  ::memcpy(target, source, (nloop * 9) * sizeof(double));
}

