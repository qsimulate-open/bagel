//
// Author : Toru Shiozaki
// Date   : May 2009
//

// edited by hand a little

#include "carsphlist.h"
#include "f77.h"


void CarSphList::carsph_20(const int nloop, const double* source, double* target) {
  const double c1 = 1.7320508075688772;
  const double c0 = 0.8660254037844386;
  const double c2 = 0.5;
  for (int iloop = 0; iloop != nloop; ++iloop, target += 5, source += 6) {
    target[0] =  c0 * (source[0] - source[2]);
    target[1] =  c1 * source[1];
    target[2] =  c1 * source[3];
    target[3] =  c1 * source[4];
    target[4] =  source[5] - c2 * (source[0] + source[2]);
  }
}

